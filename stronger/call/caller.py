import heapq
import json
import multiprocessing as mp
import multiprocessing.dummy as mpd
import parasail
import sys

from typing import Iterable, Optional

from stronger import constants as cc
from stronger.call.allele import call_alleles

__all__ = [
    "call_sample",
]

dna_matrix = parasail.matrix_create("ACGT", 2, -7)


def log(fd=sys.stderr, level: str = "ERROR"):
    def inner(message: str):
        fd.write(f"[stronger.call] {level}: {message}\n")
        fd.flush()

    return inner


log_error = log(level="ERROR")
log_warning = log(level="WARNING")
log_info = log(level="INFO")
log_debug = log(level="DEBUG")

debug = False

indel_penalty = 7


def get_repeat_count(start_count: int, tr_seq: str, flank_left_seq: str, flank_right_seq: str, motif: str) -> tuple:
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
                r_fwd = parasail.sg_de_stats_rowcol_striped_16(
                        flsub + mm, db_seq, indel_penalty, indel_penalty, dna_matrix)
                r_rev = parasail.sg_db_stats_rowcol_striped_16(
                        mm + frsub, db_seq, indel_penalty, indel_penalty, dna_matrix)
                sizes_and_scores[i] = rs = max(r_fwd.score, r_rev.score)

            szs.append((i, rs))

        mv: tuple[int, int] = max(szs, key=lambda x: x[1])
        if mv[0] > size_to_explore and (new_rc := mv[0] + 1) not in sizes_and_scores:
            to_explore.append((new_rc, 1))
        if mv[0] < size_to_explore and (new_rc := mv[0] - 1) not in sizes_and_scores:
            to_explore.append((new_rc, -1))

    # noinspection PyTypeChecker
    res: tuple[int, int] = max(sizes_and_scores.items(), key=lambda x: x[1])
    return res


def call_locus(t_idx: int, t: tuple, bf, ref, min_reads: int, min_allele_reads: int, num_bootstrap: int,
               flank_size: int, sex_chroms: Optional[str] = None,
               read_file_has_chr: bool = True, ref_file_has_chr: bool = True):
    # TODO: Figure out coords properly!!!

    contig: str = t[0]

    read_contig = ("chr" if read_file_has_chr else "") + contig.removeprefix("chr")
    ref_contig = ("chr" if ref_file_has_chr else "") + contig.removeprefix("chr")

    motif = t[-1]
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

    read_size_dict = {}

    for segment in bf.fetch(read_contig, left_flank_coord, right_flank_coord):
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
            if debug:
                log_debug(
                    f"Skipping read {segment.query_name} in locus {t_idx}: could not get sufficient flanking "
                    f"sequence")
            continue

        tr_read_seq = segment.query_sequence[left_flank_end_idx:right_flank_start_idx]

        # Truncate to flank_size (plus some leeway for small indels in flanking region) to stop any expansion sequences
        # from accidentally being included in the flanking region; e.g. if the insert gets mapped onto bases outside
        # the definition coordinates.
        # The +10 here won't include any real TR region if the mapping is solid, since the flank coordinates will
        # contain a correctly-sized sequence.
        flank_left_seq = segment.query_sequence[left_flank_start_idx:left_flank_end_idx][:flank_size+10]
        flank_right_seq = segment.query_sequence[right_flank_start_idx:right_flank_end_idx][-(flank_size+10):]

        read_rc = get_repeat_count(
            start_count=round(len(tr_read_seq) / motif_size),
            tr_seq=tr_read_seq,
            flank_left_seq=flank_left_seq,
            flank_right_seq=flank_right_seq,
            motif=motif,
            # lid=left_coord
        )
        read_size_dict[segment.query_name] = read_rc[0]

    n_alleles = 2
    if contig in ("chrM", "M"):
        n_alleles = 1
    if contig in cc.SEX_CHROMOSOMES:
        if sex_chroms is None:
            return None
        if contig in cc.X_CHROMOSOME_NAMES:
            n_alleles = sex_chroms.count("X")
        if contig in cc.Y_CHROMOSOME_NAMES:
            n_alleles = sex_chroms.count("Y")

    call = call_alleles(
        list(read_size_dict.values()), [],
        bootstrap_iterations=num_bootstrap,
        min_reads=min_reads,
        min_allele_reads=min_allele_reads,
        n_alleles=n_alleles,
        separate_strands=False,
        read_bias_corr_min=0,
        gm_filter_factor=3,
        force_int=True,
    )

    return {
        "locus_index": t_idx,
        "contig": contig,
        "start": left_coord,
        "end": right_coord,
        "motif": motif,
        "ref_cn": rc[0],
        "ref_cn_trf": t[5],
        "call": list(call[0]) if call[0] is not None else None,
        "call_95_cis": list(call[1]) if call[1] is not None else None,
        "call_99_cis": list(call[2]) if call[2] is not None else None,
        "read_cns": read_size_dict,
    }


def locus_worker(
        read_file: str,
        reference_file: str,
        min_reads: int,
        min_allele_reads: int,
        num_bootstrap: int,
        flank_size: int,
        sex_chroms: Optional[str],
        locus_queue: mp.Queue) -> list:
    import pysam as p

    ref = p.FastaFile(reference_file)
    bf = p.AlignmentFile(read_file, reference_filename=reference_file)

    ref_file_has_chr = any(r.startswith("chr") for r in ref.references)
    read_file_has_chr = any(r.startswith("chr") for r in bf.references)

    results = []
    while True:
        td = locus_queue.get()
        sys.stdout.flush()

        if td is None:  # Kill signal
            break

        # print(
        #     "read_file:", read_file,
        #     "reference_file:", reference_file,
        #     "flank_size:", flank_size,
        #     "min_reads:", min_reads,
        #     "min_allele_reads:", min_allele_reads,
        #     "num_bootstrap:", num_bootstrap,
        #     "sex_chroms", sex_chroms,
        # )

        t_idx, t = td
        res = call_locus(
            t_idx, t, bf, ref,
            min_reads=min_reads,
            min_allele_reads=min_allele_reads,
            num_bootstrap=num_bootstrap,
            flank_size=flank_size,
            sex_chroms=sex_chroms,
            read_file_has_chr=read_file_has_chr,
            ref_file_has_chr=ref_file_has_chr,
        )

        if res is not None:
            results.append(res)

    # Sort worker results; we will merge them after
    return sorted(results, key=lambda x: x["locus_index"])


def call_sample(
        read_file: str,
        reference_file: str,
        loci_file: str,
        min_reads: int = 4,
        min_allele_reads: int = 2,
        num_bootstrap: int = 100,
        flank_size: int = 70,
        sex_chroms: Optional[str] = None,
        output_format: str = "tsv",
        processes: int = 1):
    with open(loci_file, "r") as tf:
        trf_lines = [tuple(line.split("\t")) for line in (s.strip() for s in tf) if line]

    trf_lines = tuple(trf_lines)

    manager = mp.Manager()
    locus_queue = manager.Queue()

    result_lists = []

    pool_class = mp.Pool if processes > 1 else mpd.Pool

    with pool_class(processes) as p:
        jobs = []
        for _ in range(processes):
            jobs.append(p.apply_async(locus_worker, (
                read_file,
                reference_file,
                min_reads,
                min_allele_reads,
                num_bootstrap,
                flank_size,
                sex_chroms,
                locus_queue
            )))

        for t_idx, t in enumerate(trf_lines, 1):
            locus_queue.put((t_idx, t))

        for _ in range(processes):
            locus_queue.put(None)  # Kill the locus processors

        for j in jobs:
            result_lists.append(j.get())

    # Merge sorted result lists into single sorted list
    results: Iterable[dict] = heapq.merge(*result_lists, key=lambda x: x["locus_index"])

    # todo: allow simultaenous writing of tsv and json

    if output_format == "tsv":
        for res in results:
            sys.stdout.write("\t".join((
                res["contig"],
                str(res["start"]),
                str(res["end"]),
                res["motif"],
                str(res["ref_cn"]),
                str(res["ref_cn_trf"]),
                ",".join(map(str, sorted(res["read_cns"].values()))),
                str(res["call"][0]) if res["call"] is not None else ".",
                str(res["call"][1]) if res["call"] is not None else ".",
                str("-".join(map(str, res["call_95_cis"][0]))) if res["call_95_cis"] is not None else ".",
                str("-".join(map(str, res["call_95_cis"][1]))) if res["call_95_cis"] is not None else ".",
            )) + "\n")

    if output_format == "json":
        sys.stdout.write(json.dumps(results))
