import math
import parasail

from typing import Union

from strkit.utils import sign
from .align_matrix import dna_matrix, indel_penalty

__all__ = [
    "get_repeat_count",
    "get_ref_repeat_count",
]

local_search_range = 3  # TODO: parametrize


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

    res: tuple[float, int] = max(p_szs.items())
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
            if i not in sizes_and_scores:
                # Generate a candidate TR tract by copying the provided motif 'i' times & score it
                # Separate this from the .get() to postpone computation to until we need it
                sizes_and_scores[i] = score_candidate(db_seq, motif * i, flank_left_seq, flank_right_seq)

            szs.append((i, sizes_and_scores[i]))

        mv: tuple[int, int] = max(szs, key=lambda x: x[1])
        if mv[0] > size_to_explore and (new_rc := mv[0] + 1) not in sizes_and_scores:
            if new_rc >= 0:
                to_explore.append((new_rc, 1))
        elif mv[0] < size_to_explore and (new_rc := mv[0] - 1) not in sizes_and_scores:
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

            fwd_scores: list[tuple[Union[float, int], tuple[int, int], int]] = []  # For right-side adjustment
            rev_scores: list[tuple[Union[float, int], tuple[int, int], int]] = []  # For left-side adjustment

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

            mv: tuple[Union[float, int], tuple[int, int], int] = max((*fwd_scores, *rev_scores), key=lambda x: x[1])
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
