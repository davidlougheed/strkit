# Disable OpenMP multithreading since it adds enormous overhead when multiprocessing
import os
os.environ["OMP_NUM_THREADS"] = "1"

# ----------------------------------------------------------------------------------------------------------------------

import heapq
import json
import math
import multiprocessing as mp
import multiprocessing.dummy as mpd
import numpy as np
import parasail
import sys

from datetime import datetime
from pysam import AlignmentFile
from typing import Dict, List, Iterable, Optional, Tuple, Union

from strkit import __version__

from strkit.call.allele import get_n_alleles, call_alleles
from strkit.utils import apply_or_none, sign

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
local_search_range = 3  # TODO: parametrize
base_wildcard_threshold = 3

roughly_equiv_stdev_dist = 1

# TODO: Customize matrix based on error chances
# Create a substitution matrix for alignment.
# Include IUPAC wildcard bases to allow for motifs with multiple possible motifs.
# Include a wildcard base 'X' for very low-confidence base calls, to prevent needlessly harsh penalties - this is
# inserted into a read in place of bases with low PHRED scores.
dna_bases_str = "ACGTRYSWKMBDHVNX"
dna_bases = {b: i for i, b in enumerate(dna_bases_str)}
dna_codes = {
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("G", "C"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "B": ("C", "G", "T"),
    "D": ("A", "C", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "N": ("A", "C", "G", "T"),

    "X": ("A", "C", "G", "T"),  # Special character for matching low-quality bases
}
dna_matrix = parasail.matrix_create(dna_bases_str, match_score, -1 * mismatch_penalty)

for code, code_matches in dna_codes.items():
    for cm in code_matches:
        dna_matrix[dna_bases[code], dna_bases[cm]] = 2 if code != "X" else 0
        dna_matrix[dna_bases[cm], dna_bases[code]] = 2 if code != "X" else 0


def score_candidate(db_seq: str, tr_candidate: str, flank_left_seq: str, flank_right_seq: str,
                    **_kwargs) -> int:
    # TODO: sub-flank again, to avoid more errors in flanking region contributing to score?
    # Always assign parasail results to variables due to funky memory allocation behaviour
    r = parasail.sg_stats_scan_sat(
        flank_left_seq + tr_candidate + flank_right_seq, db_seq, indel_penalty, indel_penalty, dna_matrix)
    return r.score


def score_ref_boundaries(db_seq: str, tr_candidate: str, flank_left_seq: str, flank_right_seq: str,
                         **kwargs) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    ref_size = kwargs.pop("ref_size")

    # Always assign parasail results to variables due to funky memory allocation behaviour
    ext_r_seq = flank_left_seq + tr_candidate
    r_fwd = parasail.sg_de_stats_rowcol_scan_sat(
        ext_r_seq, db_seq, indel_penalty, indel_penalty, dna_matrix)
    r_adj = r_fwd.end_ref + 1 - len(flank_left_seq) - ref_size  # Amount to tweak boundary on the right side by

    db_seq_rev = db_seq[::-1]
    ext_l_seq = (tr_candidate + flank_right_seq[max(r_adj, 0):])[::-1]  # reverse

    r_rev = parasail.sg_de_stats_rowcol_scan_sat(
        ext_l_seq, db_seq_rev, indel_penalty, indel_penalty, dna_matrix)
    l_adj = r_rev.end_ref + 1 - len(flank_right_seq) - ref_size  # Amount to tweak boundary on the left side by

    return (r_fwd.score, r_adj), (r_rev.score, l_adj)


def gen_frac_repeats(motif: str, base_tr: str, j: int):
    tr_s = base_tr[abs(j):] if j < 0 else motif[j:] + base_tr
    tr_e = base_tr[:j] if j < 0 else base_tr + motif[:j]
    return tr_s, tr_e


def get_fractional_rc(top_int_res: List[Tuple[int, int]], motif: str, flank_left_seq: str, flank_right_seq: str,
                      db_seq: str) -> Tuple[float, int]:
    motif_size = len(motif)
    p_szs = {float(int_res[0]): (int_res[1], motif * int_res[0]) for int_res in top_int_res}

    j_range = range(-1 * motif_size + 1, motif_size)
    for int_res in top_int_res:
        i_mm = motif * int_res[0]  # Best integer candidate
        for j in j_range:
            if j == 0:
                # Already done
                continue

            frac_cn = int_res[0] + j / motif_size

            if frac_cn in p_szs:
                continue

            tr_s, tr_e = gen_frac_repeats(motif, i_mm, j)
            p_szs[frac_cn] = max((
                score_candidate(
                    db_seq=db_seq,
                    # Add or subtract partial copies at start
                    tr_candidate=tr_s,
                    flank_left_seq=flank_left_seq,
                    flank_right_seq=flank_right_seq,
                ),
                score_candidate(
                    db_seq=db_seq,
                    # Add or subtract partial copies at end
                    tr_candidate=tr_e,
                    flank_left_seq=flank_left_seq,
                    flank_right_seq=flank_right_seq,
                ),
            ))

    res: Tuple[float, int] = max(p_szs.items(), key=lambda x: x[1])
    return res


def get_repeat_count(
    start_count: int,
    tr_seq: str,
    flank_left_seq: str,
    flank_right_seq: str,
    motif: str,
    fractional: bool = False,
) -> Tuple[Union[int, float], int]:
    to_explore = [(start_count - 1, -1), (start_count + 1, 1), (start_count, 0)]
    sizes_and_scores: Dict[int, int] = {}

    db_seq = flank_left_seq + tr_seq + flank_right_seq

    while to_explore:
        size_to_explore, direction = to_explore.pop()
        if size_to_explore < 0:
            continue

        szs = []

        start_size = max(size_to_explore - (local_search_range if direction < 1 else 0), 0)
        end_size = size_to_explore + (local_search_range if direction > -1 else 0)

        for i in range(start_size, end_size + 1):
            rs = sizes_and_scores.get(i)
            if rs is None:
                # Generate a candidate TR tract by copying the provided motif 'i' times & score it
                # Separate this from the .get() to postpone computation to until we need it
                sizes_and_scores[i] = rs = score_candidate(db_seq, motif * i, flank_left_seq, flank_right_seq)

            szs.append((i, rs))

        mv: Tuple[int, int] = max(szs, key=lambda x: x[1])
        if mv[0] > size_to_explore and (new_rc := mv[0] + 1) not in sizes_and_scores:
            if new_rc >= 0:
                to_explore.append((new_rc, 1))
        if mv[0] < size_to_explore and (new_rc := mv[0] - 1) not in sizes_and_scores:
            if new_rc >= 0:
                to_explore.append((new_rc, -1))

    if fractional:
        # Refine further using partial copy numbers, starting from the best couple of integer copy numbers
        # noinspection PyTypeChecker
        top_int_res: List[Tuple[int, int]] = sorted(sizes_and_scores.items(), reverse=True, key=lambda x: x[1])[:3]
        return get_fractional_rc(top_int_res, motif, flank_left_seq, flank_right_seq, db_seq)

    # noinspection PyTypeChecker
    res: Tuple[int, int] = max(sizes_and_scores.items(), key=lambda x: x[1])
    return res[0], res[1]


def get_ref_repeat_count(
    start_count: int,
    tr_seq: str,
    flank_left_seq: str,
    flank_right_seq: str,
    motif: str,
    ref_size: int,
    fractional: bool = False,
):
    to_explore = [(start_count - 1, -1), (start_count + 1, 1), (start_count, 0)]

    fwd_sizes_scores_adj: Dict[Union[int, float], Tuple[int, int]] = {}
    rev_sizes_scores_adj: Dict[Union[int, float], Tuple[int, int]] = {}

    db_seq = flank_left_seq + tr_seq + flank_right_seq

    motif_size = len(motif)
    j_range = range(-1 * motif_size + 1, motif_size)

    while to_explore:
        size_to_explore, direction = to_explore.pop()
        if size_to_explore < 0:
            continue

        fwd_scores = []  # For right-side adjustment
        rev_scores = []  # For left-side adjustment

        start_size = max(size_to_explore - (local_search_range if direction < 1 else 0), 0)
        end_size = size_to_explore + (local_search_range if direction > -1 else 0)

        for i in range(start_size, end_size + 1):
            i_mm = motif * i

            if fractional:
                for j in j_range:  # j_range:
                    frac_cn = i + j/motif_size

                    fwd_rs = fwd_sizes_scores_adj.get(frac_cn)
                    rev_rs = rev_sizes_scores_adj.get(frac_cn)

                    if fwd_rs is None or rev_rs is None:
                        # Generate a candidate TR tract by copying the provided motif 'i +/- frac j' times & score it
                        # Separate this from the .get() to postpone computation to until we need it

                        tr_s, tr_e = gen_frac_repeats(motif, i_mm, j)

                        res_s = score_ref_boundaries(db_seq, tr_s, flank_left_seq, flank_right_seq, ref_size=ref_size)
                        res_e = score_ref_boundaries(db_seq, tr_e, flank_left_seq, flank_right_seq, ref_size=ref_size)

                        res = res_s if (res_s[0][0] + res_s[1][0]) > (res_e[0][0] + res_e[1][0]) else res_e

                        fwd_sizes_scores_adj[frac_cn] = fwd_rs = res[0]
                        rev_sizes_scores_adj[frac_cn] = rev_rs = res[1]

                    fwd_scores.append((frac_cn, fwd_rs, i))
                    rev_scores.append((frac_cn, rev_rs, i))

            else:
                fwd_rs = fwd_sizes_scores_adj.get(i)
                rev_rs = rev_sizes_scores_adj.get(i)

                if fwd_rs is None or rev_rs is None:
                    res = score_ref_boundaries(db_seq, i_mm, flank_left_seq, flank_right_seq, ref_size=ref_size)

                    fwd_sizes_scores_adj[i] = fwd_rs = res[0]
                    rev_sizes_scores_adj[i] = rev_rs = res[1]

                fwd_scores.append((i, fwd_rs, i))
                rev_scores.append((i, rev_rs, i))

        mv: Tuple[int, int] = max((*fwd_scores, *rev_scores), key=lambda x: x[1])
        if mv[2] > size_to_explore and (
                (new_rc := mv[2] + 1) not in fwd_sizes_scores_adj or new_rc not in rev_sizes_scores_adj):
            if new_rc >= 0:
                to_explore.append((new_rc, 1))
        if mv[2] < size_to_explore and (
                (new_rc := mv[2] - 1) not in fwd_sizes_scores_adj or new_rc not in rev_sizes_scores_adj):
            if new_rc >= 0:
                to_explore.append((new_rc, -1))

    # noinspection PyTypeChecker
    fwd_top_res: Tuple[Union[int, float], tuple] = max(fwd_sizes_scores_adj.items(), key=lambda x: x[1][0])
    # noinspection PyTypeChecker
    rev_top_res: Tuple[Union[int, float], tuple] = max(rev_sizes_scores_adj.items(), key=lambda x: x[1][0])

    # Ignore negative differences (contractions vs TRF definition), but follow expansions
    # TODO: Should we incorporate contractions? How would that work?

    if fractional:
        l_offset = rev_top_res[1][1]
        r_offset = fwd_top_res[1][1]
    else:
        l_offset = sign(rev_top_res[1][1]) * math.floor(abs(rev_top_res[1][1]) / motif_size) * motif_size
        r_offset = sign(fwd_top_res[1][1]) * math.floor(abs(fwd_top_res[1][1]) / motif_size) * motif_size

    if l_offset > 0:
        flank_left_seq = flank_left_seq[:-1*l_offset]
    if r_offset > 0:
        flank_right_seq = flank_right_seq[r_offset:]

    final_res = get_repeat_count(
        round(start_count + (max(0, l_offset) + max(0, r_offset)) / motif_size),  # always start with int here
        db_seq,
        flank_left_seq,
        flank_right_seq,
        motif)

    return final_res, l_offset, r_offset


def normalize_contig(contig: str, has_chr: bool):
    return ("chr" if has_chr else "") + contig.replace("chr", "")


def call_locus(
    t_idx: int,
    t: tuple,
    bfs: Tuple[AlignmentFile, ...],
    ref,
    min_reads: int,
    min_allele_reads: int,
    min_avg_phred: int,
    num_bootstrap: int,
    flank_size: int,
    seed: int,
    sex_chroms: Optional[str] = None,
    targeted: bool = False,
    fractional: bool = False,
    read_file_has_chr: bool = True,
    ref_file_has_chr: bool = True,
) -> Optional[dict]:
    rng = np.random.default_rng(seed=seed)

    contig: str = t[0]
    read_contig = normalize_contig(contig, read_file_has_chr)
    ref_contig = normalize_contig(contig, ref_file_has_chr)

    motif: str = t[-1]
    motif_size = len(motif)

    left_coord = int(t[1])
    right_coord = int(t[2])

    left_flank_coord = left_coord - flank_size
    right_flank_coord = right_coord + flank_size

    ref_left_flank_seq = ""
    ref_right_flank_seq = ""
    ref_seq = ""
    raised = False

    try:
        ref_left_flank_seq = ref.fetch(ref_contig, left_flank_coord, left_coord)
        ref_right_flank_seq = ref.fetch(ref_contig, right_coord, right_flank_coord)
        ref_seq = ref.fetch(ref_contig, left_coord, right_coord)
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
    ref_cn: Union[int, float]
    (ref_cn, _), l_offset, r_offset = get_ref_repeat_count(
        round(len(ref_seq) / motif_size),  # Initial estimate of copy number based on coordinates + motif size
        ref_seq,
        ref_left_flank_seq,
        ref_right_flank_seq,
        motif,
        ref_size=right_coord-left_coord,  # reference size, in terms of coordinates (not TRF-recorded size)
        fractional=fractional,
    )

    # If our reference repeat count getter has altered the TR boundaries a bit (which is done to allow for
    # more spaces in which an indel could end up), adjust our coordinates to match.
    # Currently, contractions of the TR region are ignored.
    if l_offset > 0:
        left_coord -= max(0, l_offset)
    if r_offset > 0:
        right_coord += max(0, r_offset)

    read_dict = {}

    overlapping_segments = [
        segment
        for bf in bfs
        for segment in bf.fetch(read_contig, left_flank_coord, right_flank_coord)
    ]
    read_lengths = np.fromiter((segment.query_alignment_length for segment in overlapping_segments), dtype=np.int)
    sorted_read_lengths = np.sort(read_lengths)

    for segment, read_len in zip(overlapping_segments, read_lengths):
        qs = segment.query_sequence

        if qs is None:  # No aligned segment, I guess
            log_debug(f"Skipping read {segment.query_name} (no aligned segment)")
            continue

        qq = segment.query_qualities

        left_flank_start = -1
        left_flank_end = -1
        right_flank_start = -1
        right_flank_end = -1

        last_idx = -1

        for pair in segment.get_aligned_pairs(matches_only=True):
            # Skip gaps on either side to find mapped flank indices

            if pair[1] <= left_flank_coord:
                left_flank_start = pair[0]
            elif pair[1] < left_coord:
                # Coordinate here is exclusive - we don't want to include a gap between the flanking region and
                # the STR; if we include the left-most base of the STR, we will have a giant flanking region which
                # will include part of the tandem repeat itself.
                left_flank_end = pair[0] + 1  # Add 1 to make it exclusive
            elif pair[1] >= right_coord and (
                # Reached end of TR region and haven't set end of TR region yet, or there was an indel with the motif
                # in it right after we finished due to a subtle mis-alignment - this can be seen in the HTT alignments
                # in bc1018, which were used as test data and thus will not be included in the paper...
                # TODO: do the same thing for the left side
                right_flank_start == -1 or
                (pair[0] - last_idx >= motif_size and (pair[1] - right_coord <= motif_size * 2) and
                 qs[last_idx:pair[0]].count(motif)/((pair[0]-last_idx)/motif_size) >= 0.5)
            ):
                right_flank_start = pair[0]
            elif pair[1] >= right_flank_coord:
                right_flank_end = pair[0]
                break

            last_idx = pair[0]

        if any(v == -1 for v in (left_flank_start, left_flank_end, right_flank_start, right_flank_end)):
            log_debug(
                f"Skipping read {segment.query_name} in locus {t_idx}: could not get sufficient flanking "
                f"sequence")
            continue

        qqs = np.array(qq[left_flank_end:right_flank_start])
        if qqs.shape[0] and (m_qqs := np.mean(qqs)) < min_avg_phred:
            log_debug(
                f"Skipping read {segment.query_name} due to low average base quality ({m_qqs} < {min_avg_phred}")
            continue

        # -----

        tr_read_seq = qs[left_flank_end:right_flank_start]

        # Truncate to flank_size (plus some leeway for small indels in flanking region) to stop any expansion sequences
        # from accidentally being included in the flanking region; e.g. if the insert gets mapped onto bases outside
        # the definition coordinates.
        # The +10 here won't include any real TR region if the mapping is solid, since the flank coordinates will
        # contain a correctly-sized sequence.
        flank_left_seq = qs[left_flank_start:left_flank_end][:flank_size+10]
        flank_right_seq = qs[right_flank_start:right_flank_end][-(flank_size+10):]

        tr_len = len(tr_read_seq)
        flank_len = len(flank_left_seq) + len(flank_right_seq)
        tr_len_w_flank = tr_len + flank_len

        # TODO: wildcards in flanking region too?
        tr_read_seq_wc = "".join(tr_read_seq[i] if qqs[i] > base_wildcard_threshold else "X" for i in np.arange(tr_len))

        read_cn, read_cn_score = get_repeat_count(
            start_count=round(tr_len / motif_size),  # Set initial integer copy number based on aligned TR size
            tr_seq=tr_read_seq_wc,
            flank_left_seq=flank_left_seq,
            flank_right_seq=flank_right_seq,
            motif=motif,
            fractional=fractional,
        )

        # TODO: need to rethink this; it should maybe quantify mismatches/indels in the flanking regions
        read_adj_score = match_score if tr_len == 0 else read_cn_score / tr_len_w_flank
        if read_adj_score < min_read_score:
            log_debug(f"Skipping read {segment.query_name} (scored {read_adj_score} < {min_read_score})")
            continue

        # When we don't have targeted sequencing, the probability of a read containing the TR region, given that it
        # overlaps the region, is P(read is large enough to contain) * P(  # TODO: complete this..
        partition_idx = np.searchsorted(sorted_read_lengths, tr_len_w_flank, side="right")
        if partition_idx == sorted_read_lengths.shape[0]:  # tr_len_w_flank is longer than the longest read... :(
            # Fatal
            # TODO: Just skip this locus
            log_error(
                f"Something strange happened; could not find an encompassing read where one should be guaranteed. "
                f"TRF row: {t}; TR length with flank: {tr_len_w_flank}; read lengths: {sorted_read_lengths}")
            exit(1)

        mean_containing_size = read_len if targeted else np.mean(sorted_read_lengths[partition_idx:])
        # TODO: re-examine weighting to possibly incorporate chance of drawing read large enough
        read_weight = (mean_containing_size + tr_len_w_flank - 2) / (mean_containing_size - tr_len_w_flank + 1)

        read_dict[segment.query_name] = {
            "cn": read_cn,
            "weight": read_weight,
        }

    n_alleles = get_n_alleles(2, sex_chroms, contig)
    if n_alleles is None:
        return None

    # Dicts are ordered in Python; very nice :)
    read_cns = np.fromiter((r["cn"] for r in read_dict.values()), dtype=np.float if fractional else np.int)
    read_weights = np.fromiter((r["weight"] for r in read_dict.values()), dtype=np.float)
    read_weights = read_weights / np.sum(read_weights)  # Normalize to probabilities

    call = call_alleles(
        read_cns, (),
        read_weights, (),
        bootstrap_iterations=num_bootstrap,
        min_reads=min_reads,
        min_allele_reads=min_allele_reads,
        n_alleles=n_alleles,
        separate_strands=False,
        read_bias_corr_min=0,
        gm_filter_factor=3,
        force_int=not fractional,
        seed=rng.integers(0, 4096).item(),
    ) or {}  # Still false-y

    call_peaks = call.get("peaks")
    call_weights = call.get("peak_weights")
    call_stdevs = call.get("peak_stdevs")
    call_modal_n = call.get("modal_n_peaks")

    # We cannot call read-level cluster labels with >2 peaks;
    # don't know how re-sampling has occurred.
    read_peaks_called = call_modal_n and call_modal_n <= 2
    if read_peaks_called:
        peaks = call_peaks[:call_modal_n]
        stdevs = call_stdevs[:call_modal_n]
        weights = call_weights[:call_modal_n]

        allele_reads = []
        for _ in range(call_modal_n):
            allele_reads.append([])

        for r, rd in read_dict.items():
            cn = rd["cn"]

            sd_dist = np.abs((peaks - cn) / stdevs)
            weighted_dist = np.abs(((peaks - cn) / stdevs) * weights)

            # Hack: if both peaks are 1 stdev away, pretend we aren't sure and fill in whichever allele has less
            peak: int = (
                # bool to int conversion: 1 if we add to allele_reads[1]
                int(len(allele_reads[0]) > len(allele_reads[1]))
                if call_modal_n == 2 and np.all(sd_dist < roughly_equiv_stdev_dist)
                else np.argmin(weighted_dist).item()
            )

            allele_reads[peak].append(r)
            rd["peak"] = peak

    def _round_to_base_pos(x) -> float:
        return round(float(x) * motif_size) / motif_size

    def _ndarray_serialize(x: Iterable) -> List[Union[int, float, np.int, np.float]]:
        return [(round(y) if not fractional else _round_to_base_pos(y)) for y in x]

    def _nested_ndarray_serialize(x: Iterable) -> List[List[Union[int, float, np.int, np.float]]]:
        return [_ndarray_serialize(y) for y in x]

    return {
        "locus_index": t_idx,
        "contig": contig,
        "start": left_coord,
        "end": right_coord,
        "motif": motif,
        "ref_cn": ref_cn,
        "call": apply_or_none(_ndarray_serialize, call.get("call")),
        "call_95_cis": apply_or_none(_nested_ndarray_serialize, call.get("call_95_cis")),
        "call_99_cis": apply_or_none(_nested_ndarray_serialize, call.get("call_99_cis")),
        "peaks": {
            "means": call_peaks.tolist(),  # from np.ndarray
            "weights": call_weights.tolist(),  # from np.ndarray
            "stdevs": call_stdevs.tolist(),  # from np.ndarray
            "modal_n": call_modal_n,
        } if call else None,
        "reads": read_dict,
        "read_peaks_called": read_peaks_called,
    }


def locus_worker(
        read_files: Tuple[str, ...],
        reference_file: str,
        min_reads: int,
        min_allele_reads: int,
        min_avg_phred: int,
        num_bootstrap: int,
        flank_size: int,
        sex_chroms: Optional[str],
        targeted: bool,
        fractional: bool,
        locus_queue: mp.Queue) -> List[dict]:

    import pysam as p

    ref = p.FastaFile(reference_file)
    bfs = tuple(p.AlignmentFile(rf, reference_filename=reference_file) for rf in read_files)

    ref_file_has_chr = any(r.startswith("chr") for r in ref.references)
    read_file_has_chr = any(r.startswith("chr") for bf in bfs for r in bf.references)

    results: List[dict] = []

    while True:
        td = locus_queue.get()
        if td is None:  # Kill signal
            break

        t_idx, t, locus_seed = td
        res = call_locus(
            t_idx, t, bfs, ref,
            min_reads=min_reads,
            min_allele_reads=min_allele_reads,
            min_avg_phred=min_avg_phred,
            num_bootstrap=num_bootstrap,
            flank_size=flank_size,
            seed=locus_seed,
            sex_chroms=sex_chroms,
            targeted=targeted,
            fractional=fractional,
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


def _cn_to_str(cn: Union[int, float]) -> str:
    return f"{cn:.1f}" if isinstance(cn, float) else str(cn)


def call_sample(
    read_files: Tuple[str, ...],
    reference_file: str,
    loci_file: str,
    min_reads: int = 4,
    min_allele_reads: int = 2,
    min_avg_phred: int = 13,
    num_bootstrap: int = 100,
    flank_size: int = 70,
    sex_chroms: Optional[str] = None,
    targeted: bool = False,
    fractional: bool = False,
    json_path: Optional[str] = None,
    output_tsv: bool = True,
    processes: int = 1,
    seed: Optional[int] = None,
):
    # Start the call timer
    start_time = datetime.now()

    # Seed the random number generator if a seed is provided, for replicability
    rng = np.random.default_rng(seed=seed)

    manager = mp.Manager()
    locus_queue = manager.Queue()

    # Order matters here!!
    job_params = {
        "read_files": read_files,
        "reference_file": reference_file,
        "min_reads": min_reads,
        "min_allele_reads": min_allele_reads,
        "min_avg_phred": min_avg_phred,
        "num_bootstrap": num_bootstrap,
        "flank_size": flank_size,
        "sex_chroms": sex_chroms,
        "targeted": targeted,
        "fractional": fractional,
    }

    job_args = (*job_params.values(), locus_queue)
    result_lists = []

    pool_class = mp.Pool if processes > 1 else mpd.Pool
    with pool_class(processes) as p:
        # Spin up the jobs
        jobs = [p.apply_async(locus_worker, job_args) for _ in range(processes)]

        # Add all loci from the BED file to the queue, allowing each job
        # to pull from the queue as it becomes freed up to do so.
        for t_idx, t in enumerate(parse_loci_bed(loci_file), 1):
            # We use locus-specific random seeds for replicability, no matter which order
            # the loci are yanked out of the queue / how many processes we have.
            # Tuple of (1-indexed locus index, locus data, locus-specific random seed)
            locus_queue.put((t_idx, t, rng.integers(0, 4096).item()))

        # At the end of the queue, add a None value (* the # of processes).
        # When a job encounters a None value, it will terminate.
        for _ in range(processes):
            locus_queue.put(None)

        # Gather the process-specific results for combining.
        for j in jobs:
            result_lists.append(j.get())

    # Merge sorted result lists into single sorted list.
    results: Tuple[dict, ...] = tuple(heapq.merge(*result_lists, key=lambda x: x["locus_index"]))

    time_taken = datetime.now() - start_time

    if output_tsv:
        for res in results:
            has_call = res["call"] is not None
            # n_peaks = res["peaks"]["modal_n"]

            sys.stdout.write("\t".join((
                res["contig"],
                str(res["start"]),
                str(res["end"]),
                res["motif"],
                _cn_to_str(res["ref_cn"]),
                ",".join(map(_cn_to_str, sorted(r["cn"] for r in res["reads"].values()))),
                "|".join(map(_cn_to_str, res["call"])) if has_call else ".",
                ("|".join("-".join(map(_cn_to_str, gc)) for gc in res["call_95_cis"]) if has_call else "."),

                # ("|".join(map(lambda x: f"{x:.5f}", res["peaks"]["means"][:n_peaks]))
                #  if has_call and n_peaks <= 2 else "."),
                # ("|".join(map(lambda x: f"{x:.5f}", res["peaks"]["weights"][:n_peaks]))
                #  if has_call and n_peaks <= 2 else "."),
                # ("|".join(map(lambda x: f"{x:.5f}", res["peaks"]["stdevs"][:n_peaks]))
                #  if has_call and n_peaks <= 2 else "."),
            )) + "\n")

    if json_path:
        json_report = {
            "caller": {
                "name": "strkit",
                "version": __version__,
            },
            "parameters": {
                **job_params,
                "processes": processes,
            },
            "runtime": time_taken.total_seconds(),
            "results": results,
        }

        if json_path == "stdout":
            json.dump(json_report, sys.stdout, indent=2)
            sys.stdout.write("\n")
            sys.stdout.flush()
        else:
            with open(json_path, "w") as jf:
                json.dump(json_report, jf, indent=2)
