import heapq
import json
import multiprocessing as mp
import multiprocessing.dummy as mpd
import numpy as np
import parasail
import sys

from pysam import AlignmentFile
from sklearn.mixture import GaussianMixture
from typing import List, Optional, Tuple

from strkit.call.allele import get_n_alleles, call_alleles
from strkit.utils import apply_or_none

__all__ = [
    "call_sample",
]


def log(fd=sys.stderr, level: str = "ERROR"):
    def inner(message: str):
        fd.write(f"[strkit.call] {level}: {message}\n")
        fd.flush()

    return inner


log_error = log(level="ERROR")
log_warning = log(level="WARNING")
log_info = log(level="INFO")
log_debug_ = log(level="DEBUG")

debug = False


def log_debug(*args, **kwargs):
    if debug:
        log_debug_(*args, **kwargs)


match_score = 2
mismatch_penalty = 7
indel_penalty = 5
min_read_score = 0.9  # TODO: parametrize

dna_matrix = parasail.matrix_create("ACGT", match_score, -1 * mismatch_penalty)


def get_repeat_count(
        start_count: int,
        tr_seq: str,
        flank_left_seq: str,
        flank_right_seq: str,
        motif: str
) -> Tuple[int, int]:
    moving = 0
    to_explore = [(start_count - 1, -1), (start_count + 1, 1), (start_count, 0)]
    sizes_and_scores = {}

    flsub = flank_left_seq
    frsub = flank_right_seq

    db_seq = flank_left_seq + tr_seq + flank_right_seq

    while to_explore:
        size_to_explore, direction = to_explore.pop()
        szs = []

        for i in range(size_to_explore - (1 if moving < 1 else 0), size_to_explore + (2 if moving > -1 else 1)):
            if i < 0:
                continue

            rs = sizes_and_scores.get(i)
            if rs is None:
                mm = motif * i
                r_fwd = parasail.sg_de_stats_rowcol_scan_sat(
                        flsub + mm, db_seq, indel_penalty, indel_penalty, dna_matrix)
                r_rev = parasail.sg_db_stats_rowcol_scan_sat(
                        mm + frsub, db_seq, indel_penalty, indel_penalty, dna_matrix)
                sizes_and_scores[i] = rs = max(r_fwd.score, r_rev.score)

            szs.append((i, rs))

        mv: Tuple[int, int] = max(szs, key=lambda x: x[1])
        if mv[0] > size_to_explore and (new_rc := mv[0] + 1) not in sizes_and_scores:
            to_explore.append((new_rc, 1))
        if mv[0] < size_to_explore and (new_rc := mv[0] - 1) not in sizes_and_scores:
            to_explore.append((new_rc, -1))

    # noinspection PyTypeChecker
    res: Tuple[int, int] = max(sizes_and_scores.items(), key=lambda x: x[1])
    return res


def normalize_contig(contig: str, has_chr: bool):
    return ("chr" if has_chr else "") + contig.replace("chr", "")


def call_locus(t_idx: int, t: tuple, bfs: Tuple[AlignmentFile, ...], ref, min_reads: int, min_allele_reads: int,
               num_bootstrap: int, flank_size: int, sex_chroms: Optional[str] = None, targeted: bool = False,
               read_file_has_chr: bool = True, ref_file_has_chr: bool = True) -> Optional[dict]:
    # TODO: Figure out coords properly!!!

    contig: str = t[0]
    read_contig = normalize_contig(contig, read_file_has_chr)
    ref_contig = normalize_contig(contig, ref_file_has_chr)

    motif: str = t[-1]
    motif_size = len(motif)

    left_coord = int(t[1])
    right_coord = int(t[2])

    left_flank_coord = left_coord - flank_size - 1
    right_flank_coord = right_coord + flank_size

    ref_left_flank_seq = ""
    ref_right_flank_seq = ""
    ref_seq = ""
    raised = False

    try:
        ref_left_flank_seq = ref.fetch(ref_contig, left_flank_coord, left_coord)
        ref_right_flank_seq = ref.fetch(ref_contig, right_coord - 1, right_flank_coord)
        ref_seq = ref.fetch(ref_contig, left_coord, right_coord - 1)
    except IndexError:
        log_warning(
            f"Coordinates out of range in provided reference FASTA for region {ref_contig} with flank size "
            f"{flank_size}: [{left_flank_coord}, {right_flank_coord}] (skipping locus {t_idx})")
        raised = True
    except ValueError:
        log_error(f"Invalid region '{ref_contig}' for provided reference FASTA (skipping locus {t_idx})")
        raised = True

    if len(ref_left_flank_seq) < flank_size or len(ref_right_flank_seq) < flank_size:
        if not raised:  # flank sequence too small for another reason
            log_warning(f"Reference flank size too small for locus {t_idx} (skipping)")
            return None

    if raised:
        return None

    # Get reference repeat count by our method, so we can calculate offsets from reference
    ref_size = round(len(ref_seq) / motif_size)
    rc = get_repeat_count(ref_size, ref_seq, ref_left_flank_seq, ref_right_flank_seq, motif)

    read_cn_dict = {}
    read_weight_dict = {}

    overlapping_segments = [
        segment
        for bf in bfs
        for segment in bf.fetch(read_contig, left_flank_coord, right_flank_coord)
    ]
    read_lengths = [segment.query_alignment_length for segment in overlapping_segments]
    sorted_read_lengths = sorted(read_lengths)

    for segment, read_len in zip(overlapping_segments, read_lengths):
        qs = segment.query_sequence

        if qs is None:  # No aligned segment, I guess
            log_debug(f"Skipping read {segment.query_name} (no aligned segment)")
            continue

        left_flank_start_idx = -1
        left_flank_end_idx = -1
        right_flank_start_idx = -1
        right_flank_end_idx = -1

        for pair in segment.get_aligned_pairs(matches_only=True):
            # Skip gaps on either side to find mapped flank indices

            if pair[1] <= left_flank_coord:
                left_flank_start_idx = pair[0]
            elif pair[1] < left_coord:
                # Coordinate here is exclusive - we don't want to include a gap between the flanking region and
                # the STR; if we include the left-most base of the STR, we will have a giant flanking region which
                # will include part of the tandem repeat itself.
                left_flank_end_idx = pair[0] + 1  # Add 1 to make it exclusive
            elif pair[1] < right_coord:
                right_flank_start_idx = pair[0]
            elif pair[1] >= right_flank_coord:
                right_flank_end_idx = pair[0]
                break

        if any(v == -1 for v in (
                left_flank_start_idx,
                left_flank_end_idx,
                right_flank_start_idx,
                right_flank_end_idx,
        )):
            log_debug(
                f"Skipping read {segment.query_name} in locus {t_idx}: could not get sufficient flanking "
                f"sequence")
            continue

        tr_read_seq = qs[left_flank_end_idx:right_flank_start_idx]

        # Truncate to flank_size (plus some leeway for small indels in flanking region) to stop any expansion sequences
        # from accidentally being included in the flanking region; e.g. if the insert gets mapped onto bases outside
        # the definition coordinates.
        # The +10 here won't include any real TR region if the mapping is solid, since the flank coordinates will
        # contain a correctly-sized sequence.
        flank_left_seq = qs[left_flank_start_idx:left_flank_end_idx][:flank_size+10]
        flank_right_seq = qs[right_flank_start_idx:right_flank_end_idx][-(flank_size+10):]

        tr_len = len(tr_read_seq)
        flank_len = len(flank_left_seq) + len(flank_right_seq)
        tr_len_w_flank = tr_len + flank_len

        read_rc = get_repeat_count(
            start_count=round(tr_len / motif_size),
            tr_seq=tr_read_seq,
            flank_left_seq=flank_left_seq,
            flank_right_seq=flank_right_seq,
            motif=motif,
        )

        # TODO: need to rethink this; it should maybe quantify mismatches/indels in the flanking regions
        read_adj_score = match_score if tr_len == 0 else read_rc[1] / tr_len_w_flank
        if read_adj_score < min_read_score:
            log_debug(f"Skipping read {segment.query_name} (scored {read_adj_score} < {min_read_score})")
            continue

        read_cn_dict[segment.query_name] = read_rc[0]

        # When we don't have targeted sequencing, the probability of a read containing the TR region, given that it
        # overlaps the region, is P(read is large enough to contain) * P(
        partition_idx = next((i for i in range(len(read_lengths)) if sorted_read_lengths[i] >= tr_len_w_flank), None)
        if partition_idx is None:
            # Fatal
            # TODO: Just skip this locus
            log_error(
                f"Something strange happened; could not find an encompassing read where one should be guaranteed. "
                f"TRF row: {t}; TR length with flank: {tr_len_w_flank}; read lengths: {sorted_read_lengths}")
            exit(1)
        mean_containing_size = read_len if targeted else np.mean(sorted_read_lengths[partition_idx:])
        read_weight_dict[segment.query_name] = (
                (mean_containing_size + tr_len_w_flank - 2) / (mean_containing_size - tr_len_w_flank + 1))

    n_alleles = get_n_alleles(2, sex_chroms, contig)
    if n_alleles is None:
        return None

    # Dicts are ordered in Python; very nice :)
    read_sizes = np.array(list(read_cn_dict.values()))
    read_weights = np.array(list(read_weight_dict.values()))
    read_weights = read_weights / np.sum(read_weights)  # Normalize to probabilities

    call = call_alleles(
        read_sizes, (),
        read_weights, (),
        bootstrap_iterations=num_bootstrap,
        min_reads=min_reads,
        min_allele_reads=min_allele_reads,
        n_alleles=n_alleles,
        separate_strands=False,
        read_bias_corr_min=0,
        gm_filter_factor=3,
        force_int=True,
    ) or {}  # Still false-y

    peaks_data = {
        "means": list(call["peaks"]),  # from np.ndarray
        "weights": list(call["peak_weights"]),  # from np.ndarray
        "stdevs": list(call["peak_stdevs"]),  # from np.ndarray
        "modal_n": int(call["modal_n_peaks"]),  # from np.int64
    } if call else None

    read_peak_labels = None
    # We cannot call read-level cluster labels with >2 peaks;
    # don't know how re-sampling has occurred.
    if peaks_data and peaks_data["modal_n"] <= 2:
        mn = peaks_data["modal_n"]
        ws = peaks_data["weights"][:mn]
        final_model = GaussianMixture(
            n_components=mn,
            covariance_type="spherical",
            max_iter=1,  # Lowest iteration # to keep it close to predicted parameters
            weights_init=ws/np.sum(ws),
            means_init=np.array(peaks_data["means"][:mn]).reshape(-1, 1),
            precisions_init=1 / (np.array(peaks_data["stdevs"][:mn]) ** 2),  # TODO: Check, this looks wrong
        )
        rvs = np.array(list(read_cn_dict.values())).reshape(-1, 1)
        res = final_model.fit_predict(rvs)
        read_peak_labels = {k: int(v) for k, v in zip(read_cn_dict.keys(), res)}

    def _int_ndarray_serialize(x):
        return [int(y) for y in x]

    def _nested_int_ndarray_serialize(x):
        return [_int_ndarray_serialize(y) for y in x]

    return {
        "locus_index": t_idx,
        "contig": contig,
        "start": left_coord,
        "end": right_coord,
        "motif": motif,
        "ref_cn": rc[0],
        "call": apply_or_none(_int_ndarray_serialize, call.get("call")),
        "call_95_cis": apply_or_none(_nested_int_ndarray_serialize, call.get("call_95_cis")),
        "call_99_cis": apply_or_none(_nested_int_ndarray_serialize, call.get("call_99_cis")),
        "peaks": peaks_data,
        "read_cns": read_cn_dict,
        "read_weights": read_weight_dict,
        "read_peak_labels": read_peak_labels,
    }


def locus_worker(
        read_files: Tuple[str, ...],
        reference_file: str,
        min_reads: int,
        min_allele_reads: int,
        num_bootstrap: int,
        flank_size: int,
        sex_chroms: Optional[str],
        targeted: bool,
        locus_queue: mp.Queue) -> List[dict]:

    import pysam as p

    ref = p.FastaFile(reference_file)
    bfs = tuple(p.AlignmentFile(rf, reference_filename=reference_file) for rf in read_files)

    ref_file_has_chr = any(r.startswith("chr") for r in ref.references)
    read_file_has_chr = any(r.startswith("chr") for bf in bfs for r in bf.references)

    results: List[dict] = []

    while True:
        td = locus_queue.get()
        sys.stdout.flush()

        if td is None:  # Kill signal
            break

        t_idx, t = td
        res = call_locus(
            t_idx, t, bfs, ref,
            min_reads=min_reads,
            min_allele_reads=min_allele_reads,
            num_bootstrap=num_bootstrap,
            flank_size=flank_size,
            sex_chroms=sex_chroms,
            targeted=targeted,
            read_file_has_chr=read_file_has_chr,
            ref_file_has_chr=ref_file_has_chr,
        )

        if res is not None:
            results.append(res)

    # Sort worker results; we will merge them after
    return sorted(results, key=lambda x: x["locus_index"])


def parse_loci_bed(loci_file: str):
    with open(loci_file, "r") as tf:
        yield from (tuple(line.split("\t")) for line in (s.strip() for s in tf) if line)


def call_sample(
        read_files: Tuple[str, ...],
        reference_file: str,
        loci_file: str,
        min_reads: int = 4,
        min_allele_reads: int = 2,
        num_bootstrap: int = 100,
        flank_size: int = 70,
        sex_chroms: Optional[str] = None,
        targeted: bool = False,
        json_path: Optional[str] = None,
        output_tsv: bool = True,
        processes: int = 1):

    manager = mp.Manager()
    locus_queue = manager.Queue()

    job_args = (
        read_files,
        reference_file,
        min_reads,
        min_allele_reads,
        num_bootstrap,
        flank_size,
        sex_chroms,
        targeted,
        locus_queue,
    )
    result_lists = []

    pool_class = mp.Pool if processes > 1 else mpd.Pool
    with pool_class(processes) as p:
        # Spin up the jobs
        jobs = [p.apply_async(locus_worker, job_args) for _ in range(processes)]

        # Add all loci from the BED file to the queue, allowing each job
        # to pull from the queue as it becomes freed up to do so.
        for t_idx, t in enumerate(parse_loci_bed(loci_file), 1):
            locus_queue.put((t_idx, t))

        # At the end of the queue, add a None value (* the # of processes).
        # When a job encounters a None value, it will terminate.
        for _ in range(processes):
            locus_queue.put(None)

        # Gather the process-specific results for combining.
        for j in jobs:
            result_lists.append(j.get())

    # Merge sorted result lists into single sorted list.
    results: Tuple[dict, ...] = tuple(heapq.merge(*result_lists, key=lambda x: x["locus_index"]))

    if output_tsv:
        for res in results:
            has_call = res["call"] is not None
            # n_peaks = res["peaks"]["modal_n"]

            sys.stdout.write("\t".join((
                res["contig"],
                str(res["start"]),
                str(res["end"]),
                res["motif"],
                str(res["ref_cn"]),
                ",".join(map(str, sorted(res["read_cns"].values()))),
                "|".join(map(str, res["call"])) if has_call else ".",
                ("|".join("-".join(map(str, gc)) for gc in res["call_95_cis"]) if has_call else "."),

                # ("|".join(map(lambda x: f"{x:.5f}", res["peaks"]["means"][:n_peaks]))
                #  if has_call and n_peaks <= 2 else "."),
                # ("|".join(map(lambda x: f"{x:.5f}", res["peaks"]["weights"][:n_peaks]))
                #  if has_call and n_peaks <= 2 else "."),
                # ("|".join(map(lambda x: f"{x:.5f}", res["peaks"]["stdevs"][:n_peaks]))
                #  if has_call and n_peaks <= 2 else "."),
            )) + "\n")

    if json_path:
        if json_path == "stdout":
            json.dump(results, sys.stdout)
            sys.stdout.write("\n")
            sys.stdout.flush()
        else:
            with open(json_path, "w") as jf:
                json.dump(results, jf)
