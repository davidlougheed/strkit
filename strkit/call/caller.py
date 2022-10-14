from __future__ import annotations

# Disable OpenMP multithreading since it adds enormous overhead when multiprocessing
import os
os.environ["OMP_NUM_THREADS"] = "1"

# ----------------------------------------------------------------------------------------------------------------------

import heapq
import itertools
import logging
import math
import multiprocessing as mp
import multiprocessing.dummy as mpd
import numpy as np
import parasail
import queue
import sys

from collections import Counter
from ctypes import c_uint32
from datetime import datetime
from pysam import AlignmentFile, FastaFile
from typing import Iterable, Generator, Optional, Union

from strkit import __version__

from strkit.call.allele import get_n_alleles, call_alleles
from strkit.json import dumps_indented, json
from strkit.logger import logger
from strkit.utils import apply_or_none, sign

__all__ = [
    "call_sample",
]

PROFILE_LOCUS_CALLS = False

# TODO: Parameterize
LOG_PROGRESS_INTERVAL = 1000

match_score = 2  # TODO: parametrize
mismatch_penalty = 7  # TODO: parametrize
indel_penalty = 5  # TODO: parametrize
realign_indel_open_penalty = 7  # TODO: parametrize
realign_timeout = 5
min_read_score = 0.9  # TODO: parametrize
min_realign_score_ratio = 0.95  # TODO: parametrize
local_search_range = 3  # TODO: parametrize
base_wildcard_threshold = 3

roughly_equiv_stdev_dist = 1

force_realign = False

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
                         **kwargs) -> tuple[tuple[int, int], tuple[int, int]]:
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


def gen_frac_repeats(motif: str, base_tr: str, j: int) -> tuple[str, str]:
    tr_s = base_tr[abs(j):] if j < 0 else motif[j:] + base_tr
    tr_e = base_tr[:j] if j < 0 else base_tr + motif[:j]
    return tr_s, tr_e


def get_fractional_rc(
    top_int_res: list[tuple[int, int]],
    motif: str,
    flank_left_seq: str,
    flank_right_seq: str,
    db_seq: str,
) -> tuple[float, int]:
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

    res: tuple[float, int] = max(p_szs.items(), key=lambda x: x[1])
    return res


def get_repeat_count(
    start_count: int,
    tr_seq: str,
    flank_left_seq: str,
    flank_right_seq: str,
    motif: str,
    fractional: bool = False,
) -> tuple[Union[int, float], int]:
    to_explore = [(start_count - 1, -1), (start_count + 1, 1), (start_count, 0)]
    sizes_and_scores: dict[int, int] = {}

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

        mv: tuple[int, int] = max(szs, key=lambda x: x[1])
        if mv[0] > size_to_explore and (new_rc := mv[0] + 1) not in sizes_and_scores:
            if new_rc >= 0:
                to_explore.append((new_rc, 1))
        if mv[0] < size_to_explore and (new_rc := mv[0] - 1) not in sizes_and_scores:
            if new_rc >= 0:
                to_explore.append((new_rc, -1))

    if fractional:
        # Refine further using partial copy numbers, starting from the best couple of integer copy numbers
        # noinspection PyTypeChecker
        top_int_res: list[tuple[int, int]] = sorted(sizes_and_scores.items(), reverse=True, key=lambda x: x[1])[:3]
        return get_fractional_rc(top_int_res, motif, flank_left_seq, flank_right_seq, db_seq)

    # noinspection PyTypeChecker
    return max(sizes_and_scores.items(), key=lambda x: x[1])


def get_ref_repeat_count(
    start_count: int,
    tr_seq: str,
    flank_left_seq: str,
    flank_right_seq: str,
    motif: str,
    ref_size: int,
    fractional: bool = False,
    respect_coords: bool = False,
) -> tuple[tuple[Union[int, float], int], int, int]:
    l_offset: int = 0
    r_offset: int = 0

    db_seq: str = flank_left_seq + tr_seq + flank_right_seq
    motif_size = len(motif)

    if not respect_coords:  # Extend out coordinates from initial definition
        to_explore = [(start_count - 1, -1), (start_count + 1, 1), (start_count, 0)]

        fwd_sizes_scores_adj: dict[Union[int, float], tuple[int, int]] = {}
        rev_sizes_scores_adj: dict[Union[int, float], tuple[int, int]] = {}

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
                            # Generate a candidate TR tract by copying the provided motif 'i +/- frac j' times & score
                            # it Separate this from the .get() to postpone computation to until we need it

                            tr_s, tr_e = gen_frac_repeats(motif, i_mm, j)

                            res_s = score_ref_boundaries(
                                db_seq, tr_s, flank_left_seq, flank_right_seq, ref_size=ref_size)
                            res_e = score_ref_boundaries(
                                db_seq, tr_e, flank_left_seq, flank_right_seq, ref_size=ref_size)

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

            mv: tuple[int, int] = max((*fwd_scores, *rev_scores), key=lambda x: x[1])
            if mv[2] > size_to_explore and (
                    (new_rc := mv[2] + 1) not in fwd_sizes_scores_adj or new_rc not in rev_sizes_scores_adj):
                if new_rc >= 0:
                    to_explore.append((new_rc, 1))
            if mv[2] < size_to_explore and (
                    (new_rc := mv[2] - 1) not in fwd_sizes_scores_adj or new_rc not in rev_sizes_scores_adj):
                if new_rc >= 0:
                    to_explore.append((new_rc, -1))

        # noinspection PyTypeChecker
        fwd_top_res: tuple[Union[int, float], tuple] = max(fwd_sizes_scores_adj.items(), key=lambda x: x[1][0])
        # noinspection PyTypeChecker
        rev_top_res: tuple[Union[int, float], tuple] = max(rev_sizes_scores_adj.items(), key=lambda x: x[1][0])

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

    # ------------------------------------------------------------------------------------------------------------------

    final_res = get_repeat_count(
        round(start_count + (max(0, l_offset) + max(0, r_offset)) / motif_size),  # always start with int here
        db_seq,
        flank_left_seq,
        flank_right_seq,
        motif)

    return final_res, l_offset, r_offset


def normalize_contig(contig: str, has_chr: bool):
    return ("chr" if has_chr else "") + contig.replace("chr", "")


def get_read_coords_from_matched_pairs(
    left_flank_coord: int,
    left_coord: int,
    right_coord: int,
    right_flank_coord: int,
    motif: str,
    motif_size: int,
    query_seq: str,
    matched_pairs
) -> tuple[int, int, int, int]:
    left_flank_start = -1
    left_flank_end = -1
    right_flank_start = -1
    right_flank_end = -1

    last_idx = -1

    for query_coord, ref_coord in matched_pairs:
        # Skip gaps on either side to find mapped flank indices
        if ref_coord <= left_flank_coord:
            left_flank_start = query_coord
        elif ref_coord < left_coord:
            # Coordinate here is exclusive - we don't want to include a gap between the flanking region and
            # the STR; if we include the left-most base of the STR, we will have a giant flanking region which
            # will include part of the tandem repeat itself.
            left_flank_end = query_coord + 1  # Add 1 to make it exclusive
        elif ref_coord >= right_coord and (
                # Reached end of TR region and haven't set end of TR region yet, or there was an indel with the motif
                # in it right after we finished due to a subtle mis-alignment - this can be seen in the HTT alignments
                # in bc1018, which were used as test data and thus will not be included in the paper...
                # TODO: do the same thing for the left side
                right_flank_start == -1 or
                (query_coord - last_idx >= motif_size and (ref_coord - right_coord <= motif_size * 2) and
                 query_seq[last_idx:query_coord].count(motif) / ((query_coord - last_idx) / motif_size) >= 0.5)
        ):
            right_flank_start = query_coord
        elif ref_coord >= right_flank_coord:
            right_flank_end = query_coord
            break

    return left_flank_start, left_flank_end, right_flank_start, right_flank_end


def decode_cigar(encoded_cigar: list[int]) -> Generator[tuple[int, int], None, None]:
    for item in encoded_cigar:
        yield item & 15, item >> 4


def get_aligned_pairs_from_cigar(
    cigar: Iterable[tuple[int, int]],
    query_start: int = 0,
    ref_start: int = 0,
) -> Generator[tuple[Union[int, None], Union[int, None]], None, None]:
    qi = itertools.count(start=query_start)
    di = itertools.count(start=ref_start)

    for c in cigar:
        op, cnt = c
        rc = range(cnt)

        # TODO: Probably a nicer way to do this:

        cq = ((next(qi), None) for _ in rc)
        cr = ((None, next(di)) for _ in rc)
        cb = ((next(qi), next(di)) for _ in rc)

        yield from {
            0: cb,  # M
            1: cq,  # I
            2: cr,  # D
            3: cr,  # N
            4: cq,  # S
            5: (),  # H
            6: (),  # P
            7: cb,  # =
            8: cb,  # X
        }[op]


def realign_read(
    ref_seq: str,
    query_seq: str,
    left_flank_coord: int,
    flank_size: int,
    rn: str,
    t_idx: int,
    always_realign: bool,
    q: Optional[mp.Queue] = None,
    log_level: int = logging.WARNING,
) -> Optional[list[tuple[Optional[int], Optional[int]]]]:
    # Have to re-attach logger in separate process I guess

    from strkit.logger import logger as lg, attach_stream_handler
    attach_stream_handler(log_level)

    # flipped: 'ref sequence' as query here, since it should in general be shorter
    pr = parasail.sg_dx_trace_scan_sat(
        # fetch an extra base for the right flank coordinate check later (needs to be >= the exclusive coord)
        ref_seq, query_seq, realign_indel_open_penalty, 0, dna_matrix)

    if pr.score < (th := min_realign_score_ratio * (flank_size * 2 * match_score - realign_indel_open_penalty)):
        lg.debug(f"Realignment for {rn} scored below threshold ({pr.score} < {th:.2f})")
        q.put(None)
        return None

    lg.debug(
        f"Realigned {rn} in locus {t_idx}{' (due to soft clipping)' if not always_realign else ''}: scored {pr.score}; "
        f"CIGAR: {pr.cigar.decode.decode('ascii')}")
    # noinspection PyTypeChecker
    res: list[tuple[Optional[int], Optional[int]]] = [
        # reverse to get (ref, query) instead of (query, ref) due to the flip
        tuple(reversed(p))
        # query_start instead of ref_start, due to the flip
        for p in get_aligned_pairs_from_cigar(decode_cigar(pr.cigar.seq), query_start=left_flank_coord)
        if p[0] is not None and p[1] is not None
    ]
    if q:
        q.put(res)
    return res


def calculate_seq_with_wildcards(qs: str, quals: list[int]):
    return "".join(qs[i] if quals[i] > base_wildcard_threshold else "X" for i in np.arange(len(qs)))


def round_to_base_pos(x, motif_size) -> float:
    return round(float(x) * motif_size) / motif_size


def call_locus(
    t_idx: int,
    t: tuple,
    bfs: tuple[AlignmentFile, ...],
    ref: FastaFile,
    min_reads: int,
    min_allele_reads: int,
    min_avg_phred: int,
    num_bootstrap: int,
    flank_size: int,
    seed: int,
    sex_chroms: Optional[str] = None,
    realign: bool = False,
    hq: bool = False,
    targeted: bool = False,
    fractional: bool = False,
    respect_ref: bool = False,
    count_kmers: str = "none",  # "none" | "peak" | "read"
    log_level: int = logging.WARNING,
    read_file_has_chr: bool = True,
    ref_file_has_chr: bool = True,
) -> Optional[dict]:
    call_timer = datetime.now()

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
        logger.warning(
            f"Coordinates out of range in provided reference FASTA for region {ref_contig} with flank size "
            f"{flank_size}: [{left_flank_coord}, {right_flank_coord}] (skipping locus {t_idx})")
        raised = True
    except ValueError:
        logger.error(f"Invalid region '{ref_contig}' for provided reference FASTA (skipping locus {t_idx})")
        raised = True

    if len(ref_left_flank_seq) < flank_size or len(ref_right_flank_seq) < flank_size:
        if not raised:  # flank sequence too small for another reason
            logger.warning(f"Reference flank size too small for locus {t_idx} (skipping)")
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
        respect_coords=respect_ref,
    )

    # If our reference repeat count getter has altered the TR boundaries a bit (which is done to allow for
    # more spaces in which an indel could end up), adjust our coordinates to match.
    # Currently, contractions of the TR region are ignored.
    left_coord_adj = left_coord if respect_ref else left_coord - max(0, l_offset)
    right_coord_adj = right_coord if respect_ref else right_coord + max(0, r_offset)

    read_dict: dict[str, dict] = {}
    chimeric_read_status: dict[str, int] = {}

    overlapping_segments = []
    seen_reads: set[str] = set()
    read_lengths = []

    for segment in (s for bf in bfs for s in bf.fetch(read_contig, left_flank_coord, right_flank_coord)):
        rn = segment.query_name
        supp = segment.flag & 2048

        # If we have two overlapping alignments for the same read, we have a chimeric read within the TR
        # (so probably a large expansion...)
        qn = segment.query_name
        chimeric_read_status[qn] = chimeric_read_status.get(qn, 0) | (2 if supp else 1)

        if supp:  # Skip supplemental alignments
            logger.debug(f"Skipping entry for read {rn} (supplemental)")
            continue

        if rn in seen_reads:
            logger.debug(f"Skipping entry for read {rn} (already seen)")
            continue

        if segment.query_sequence is None:
            logger.debug(f"Skipping entry for read {rn} (no aligned segment)")
            continue

        seen_reads.add(rn)
        overlapping_segments.append(segment)
        read_lengths.append(segment.query_alignment_length)

    sorted_read_lengths = np.sort(read_lengths)

    for segment, read_len in zip(overlapping_segments, read_lengths):
        rn = segment.query_name
        qs = segment.query_sequence

        c1: tuple[int, int] = segment.cigar[0]
        c2: tuple[int, int] = segment.cigar[-1]

        fqqs = segment.query_qualities

        realigned = False
        pairs = None

        # Soft-clipping in large insertions can result from mapping difficulties.
        # If we have a soft clip which overlaps with our TR region (+ flank), we can try to recover it
        # via realignment with parasail.
        # 4: BAM code for soft clip
        # TODO: if some of the BAM alignment is present, use it to reduce realignment overhead?
        #  - use start point + flank*3 or end point - flank*3 or something like that
        if realign and (force_realign or (
            (c1[0] == 4 and segment.reference_start > left_flank_coord >= segment.reference_start - c1[1]) or
            (c2[0] == 4 and segment.reference_end < right_flank_coord <= segment.reference_end + c2[1])
        )):
            # Run the realignment in a separate process, to give us a timeout mechanism.
            # This means we're spawning a second process for this job, just momentarily, beyond the pool size.

            q = mp.Queue()
            proc = mp.Process(target=realign_read, kwargs=dict(
                # fetch an extra base for the right flank coordinate check later (needs to be >= the exclusive coord)
                ref_seq=ref_left_flank_seq + ref_seq + ref.fetch(ref_contig, right_coord, right_flank_coord + 1),
                query_seq=calculate_seq_with_wildcards(qs, fqqs),
                left_flank_coord=left_flank_coord,
                flank_size=flank_size,
                rn=rn,
                t_idx=t_idx,
                always_realign=force_realign,
                q=q,
                log_level=log_level,
            ))
            proc.start()
            pairs_new = None
            try:
                pairs_new = q.get(timeout=realign_timeout)
                proc.join()
            except queue.Empty:
                logger.warning(
                    f"Experienced timeout while re-aligning read {rn} in locus {t_idx}. Reverting to BAM alignment.")
            finally:
                proc.close()

            if pairs_new is not None:
                pairs = pairs_new
                realigned = True

        if pairs is None:
            pairs = segment.get_aligned_pairs(matches_only=True)

        left_flank_start, left_flank_end, right_flank_start, right_flank_end = get_read_coords_from_matched_pairs(
            left_flank_coord,
            left_coord_adj,
            right_coord_adj,
            right_flank_coord,
            motif,
            motif_size,
            query_seq=qs,
            matched_pairs=pairs
        )

        if any(v == -1 for v in (left_flank_start, left_flank_end, right_flank_start, right_flank_end)):
            logger.debug(
                f"Skipping read {rn} in locus {t_idx}: could not get sufficient flanking sequence"
                f"{' (post-realignment)' if realigned else ''}")
            continue

        qqs = np.array(fqqs[left_flank_end:right_flank_start])
        if qqs.shape[0] and (m_qqs := np.mean(qqs)) < min_avg_phred:  # TODO: check flank?
            logger.debug(
                f"Skipping read {rn} in locus {t_idx} due to low average base quality ({m_qqs:.2f} < {min_avg_phred})")
            continue

        # -----

        tr_read_seq = qs[left_flank_end:right_flank_start]

        # Truncate to flank_size (plus some leeway for small indels in flanking region) to stop any expansion sequences
        # from accidentally being included in the flanking region; e.g. if the insert gets mapped onto bases outside
        # the definition coordinates.
        # The +10 here won't include any real TR region if the mapping is solid, since the flank coordinates will
        # contain a correctly-sized sequence.

        # TODO: wildcards in flanking region too?
        flank_left_seq = qs[left_flank_start:left_flank_end][:flank_size+10]
        flank_right_seq = qs[right_flank_start:right_flank_end][-(flank_size+10):]

        tr_len = len(tr_read_seq)
        flank_len = len(flank_left_seq) + len(flank_right_seq)
        tr_len_w_flank = tr_len + flank_len

        tr_read_seq_wc = calculate_seq_with_wildcards(qs[left_flank_end:right_flank_start], qqs)

        read_kmers = Counter()
        if count_kmers != "none":
            for i in range(0, tr_len - motif_size + 1):
                read_kmers.update((tr_read_seq_wc[i:i+motif_size],))

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
            logger.debug(f"Skipping read {segment.query_name} (scored {read_adj_score} < {min_read_score})")
            continue

        # When we don't have targeted sequencing, the probability of a read containing the TR region, given that it
        # overlaps the region, is P(read is large enough to contain) * P(  # TODO: complete this..
        partition_idx = np.searchsorted(sorted_read_lengths, tr_len_w_flank, side="right")
        if partition_idx == sorted_read_lengths.shape[0]:  # tr_len_w_flank is longer than the longest read... :(
            # Fatal
            # TODO: Just skip this locus
            logger.error(
                f"Something strange happened; could not find an encompassing read where one should be guaranteed. "
                f"TRF row: {t}; TR length with flank: {tr_len_w_flank}; read lengths: {sorted_read_lengths}")
            exit(1)

        mean_containing_size = read_len if targeted else np.mean(sorted_read_lengths[partition_idx:])
        # TODO: re-examine weighting to possibly incorporate chance of drawing read large enough
        read_weight = (mean_containing_size + tr_len_w_flank - 2) / (mean_containing_size - tr_len_w_flank + 1)

        crs_cir = chimeric_read_status[rn] == 3  # Chimera within the TR region, indicating a potential large expansion
        read_dict[rn] = {
            "s": "+" if segment.is_forward else "-",
            "cn": read_cn,
            "w": read_weight.item(),
            **({"realn": realigned} if realign and realigned else {}),
            **({"chimeric_in_region": crs_cir} if crs_cir else {}),
            **({"kmers": dict(read_kmers)} if count_kmers != "none" else {}),
        }

    n_alleles = get_n_alleles(2, sex_chroms, contig)
    if n_alleles is None:
        return None

    # Dicts are ordered in Python; very nice :)
    rdvs = tuple(read_dict.values())
    rcns = tuple(r["cn"] for r in rdvs)
    read_cns = np.fromiter(rcns, dtype=np.float if fractional else np.int)
    read_weights = np.fromiter((r["w"] for r in rdvs), dtype=np.float)
    read_weights = read_weights / np.sum(read_weights)  # Normalize to probabilities

    call_dict_base = {
        "locus_index": t_idx,
        "contig": contig,
        "start": left_coord,
        "end": right_coord,
        **({} if respect_ref else {
            "start_adj": left_coord_adj,
            "end_adj": right_coord_adj,
        }),
        "motif": motif,
        "ref_cn": ref_cn,
        "reads": read_dict,
    }

    # Check now if we definitely don't have enough reads to make a call
    # We also check again later when we calculate all the flanking stuff
    if len(overlapping_segments) < min_reads:
        return {
            **call_dict_base,
            "call": None,
            "call_95_cis": None,
            "call_99_cis": None,
            "peaks": None,
            "read_peaks_called": False,
            "time": (datetime.now() - call_timer).total_seconds(),
        }

    # If the locus only has one value, don't bother bootstrapping
    if hq and len(set(rcns)) == 1:
        logger.debug(f"Skipping bootstrap for locus at {contig}:{left_coord}-{right_coord} (single value)")
        num_bootstrap = 1

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
        hq=hq,
        force_int=not fractional,
        seed=rng.integers(0, 4096).item(),
    ) or {}  # Still false-y

    call_peaks = call.get("peaks")
    call_weights = call.get("peak_weights")
    call_stdevs = call.get("peak_stdevs")
    call_modal_n = call.get("modal_n_peaks")
    call_peak_n_reads = []

    # We cannot call read-level cluster labels with >2 peaks;
    # don't know how re-sampling has occurred.
    peak_kmers: list[Counter] = [Counter() for _ in range(call_modal_n or 0)]
    if read_peaks_called := call_modal_n and call_modal_n <= 2:
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
            rd["p"] = peak

            if count_kmers in ("peak", "both"):
                peak_kmers[peak] += rd["kmers"]

                # If we aren't reporting read-level k-mers, we have to delete them (space-saving!)
                if count_kmers == "peak":
                    del rd["kmers"]

        call_peak_n_reads = list(map(len, allele_reads))

    def _ndarray_serialize(x: Iterable) -> list[Union[int, float, np.int, np.float]]:
        return [(round(y) if not fractional else round_to_base_pos(y, motif_size)) for y in x]

    def _nested_ndarray_serialize(x: Iterable) -> list[list[Union[int, float, np.int, np.float]]]:
        return [_ndarray_serialize(y) for y in x]

    return {
        **call_dict_base,
        "call": apply_or_none(_ndarray_serialize, call.get("call")),
        "call_95_cis": apply_or_none(_nested_ndarray_serialize, call.get("call_95_cis")),
        "call_99_cis": apply_or_none(_nested_ndarray_serialize, call.get("call_99_cis")),
        "peaks": {
            "means": call_peaks.tolist(),  # from np.ndarray
            "weights": call_weights.tolist(),  # from np.ndarray
            "stdevs": call_stdevs.tolist(),  # from np.ndarray
            "modal_n": call_modal_n,
            "n_reads": call_peak_n_reads,
            **({"kmers": [dict(c) for c in peak_kmers]} if count_kmers in ("peak", "both") else {}),
        } if call else None,
        "read_peaks_called": read_peaks_called,
        "time": (datetime.now() - call_timer).total_seconds(),
    }


def locus_worker(
    read_files: tuple[str, ...],
    reference_file: str,
    min_reads: int,
    min_allele_reads: int,
    min_avg_phred: int,
    num_bootstrap: int,
    flank_size: int,
    sex_chroms: Optional[str],
    realign: bool,
    hq: bool,
    targeted: bool,
    fractional: bool,
    respect_ref: bool,
    count_kmers: str,
    log_level: int,
    start_time: datetime,
    locus_queue: mp.Queue,
    num_processed: mp.Value,
    is_single_processed: bool,
) -> list[dict]:
    if PROFILE_LOCUS_CALLS:
        import cProfile
        pr = cProfile.Profile()
        pr.enable()
    else:
        pr = None

    import pysam as p

    ref = p.FastaFile(reference_file)
    bfs = tuple(p.AlignmentFile(rf, reference_filename=reference_file) for rf in read_files)

    ref_file_has_chr = any(r.startswith("chr") for r in ref.references)
    read_file_has_chr = any(r.startswith("chr") for bf in bfs for r in bf.references)

    results: list[dict] = []

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
            realign=realign,
            hq=hq,
            targeted=targeted,
            fractional=fractional,
            respect_ref=respect_ref,
            count_kmers=count_kmers,
            log_level=log_level,
            read_file_has_chr=read_file_has_chr,
            ref_file_has_chr=ref_file_has_chr,
        )

        with num_processed.get_lock():
            n_proc = num_processed.value + 1
            num_processed.value = n_proc
            # Release the lock early (before logging)

        if n_proc % LOG_PROGRESS_INTERVAL == 0:
            logger.info(
                f"Processed {n_proc} loci in "
                f"{(datetime.now() - start_time).total_seconds():.1f} seconds")

        if res is not None:
            results.append(res)

    if PROFILE_LOCUS_CALLS:
        pr.disable()
        pr.print_stats("tottime")

    # Sort worker results; we will merge them after
    return results if is_single_processed else sorted(results, key=lambda x: x["locus_index"])


def parse_loci_bed(loci_file: str):
    with open(loci_file, "r") as tf:
        yield from (tuple(line.split("\t")) for line in (s.strip() for s in tf) if line and not line.startswith("#"))


def _cn_to_str(cn: Union[int, float]) -> str:
    return f"{cn:.1f}" if isinstance(cn, float) else str(cn)


def call_sample(
    read_files: tuple[str, ...],
    reference_file: str,
    loci_file: str,
    sample_id: Optional[str],
    min_reads: int = 4,
    min_allele_reads: int = 2,
    min_avg_phred: int = 13,
    num_bootstrap: int = 100,
    flank_size: int = 70,
    sex_chroms: Optional[str] = None,
    realign: bool = False,
    hq: bool = False,
    targeted: bool = False,
    fractional: bool = False,
    respect_ref: bool = False,
    count_kmers: str = "none",  # "none" | "peak" | "read"
    log_level: int = logging.WARNING,
    json_path: Optional[str] = None,
    indent_json: bool = False,
    output_tsv: bool = True,
    processes: int = 1,
    seed: Optional[int] = None,
) -> None:
    # Start the call timer
    start_time = datetime.now()

    # Get sample ID from read file(s)

    bfs = tuple(AlignmentFile(rf, reference_filename=reference_file) for rf in read_files)

    # noinspection PyTypeChecker
    bfhs = [dict(bf.header) for bf in bfs]

    sns: set[str] = {e.get("SM") for bfh in bfhs for e in bfh.get("RG", ()) if e.get("SM")}
    bam_sample_id: Optional[str] = None

    if len(sns) > 1:
        # Error or warning or what?
        sns_str = "', '".join(sns)
        logger.warning(f"Found more than one sample ID in BAM file(s): '{sns_str}'")
    elif not sns:
        logger.warning("Could not find sample ID in BAM file(s)")
    else:
        bam_sample_id = sns.pop()

    sample_id: Optional[str] = sample_id or bam_sample_id

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
        "realign": realign,
        "hq": hq,
        "targeted": targeted,
        "fractional": fractional,
        "respect_ref": respect_ref,
        "count_kmers": count_kmers,
        "log_level": log_level,
    }

    num_processed = mp.Value(c_uint32)
    is_single_processed = processes == 1

    job_args = (*job_params.values(), start_time, locus_queue, num_processed, is_single_processed)
    result_lists = []

    pool_class = mpd.Pool if is_single_processed else mp.Pool
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
    results: tuple[dict, ...] = tuple(heapq.merge(*result_lists, key=lambda x: x["locus_index"]))

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
            "sample_id": sample_id,
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

        dfn = dumps_indented if indent_json else json.dumps
        if json_path == "stdout":
            sys.stdout.buffer.write(dfn(json_report))
            sys.stdout.write("\n")
            sys.stdout.flush()
        else:
            with open(json_path, "wb") as jf:
                jf.write(dfn(json_report))
                jf.write(b"\n")
