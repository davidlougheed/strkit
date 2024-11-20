import parasail

from functools import lru_cache
from typing import Literal, Union

from strkit_rust_ext import get_repeat_count as _get_repeat_count

from .align_matrix import dna_matrix, indel_penalty
from .utils import idx_1_getter

__all__ = [
    "get_repeat_count",
    "get_ref_repeat_count",
]


DEFAULT_LOCAL_SEARCH_RANGE = 3


def score_candidate_with_string(db_seq_profile: parasail.Profile, tr_seq: str) -> int:
    # TODO: sub-flank again, to avoid more errors in flanking region contributing to score?
    # Always assign parasail results to variables due to funky memory allocation behaviour
    #  - switch 'db' and 'query' here so we can use the db sequence as the profile for a "database" search against
    #    candidate sequences. order doesn't end up mattering, since we're using semi-global alignment.
    r = parasail.sg_striped_profile_sat(db_seq_profile, tr_seq, indel_penalty, indel_penalty)
    return r.score


def score_candidate(
    db_seq_profile: parasail.Profile,
    motif: str,
    motif_count: int,
    flank_left_seq: str,
    flank_right_seq: str,
) -> int:
    return score_candidate_with_string(db_seq_profile, f"{flank_left_seq}{motif * motif_count}{flank_right_seq}")


def score_ref_boundaries(
    db_seq_profile: parasail.Profile,
    db_seq_rev_profile: parasail.Profile,
    tr_candidate: str,
    flank_left_seq: str,
    flank_right_seq: str,
    ref_size: int,
) -> tuple[tuple[int, int], tuple[int, int]]:
    # Always assign parasail results to variables due to funky memory allocation behaviour
    ext_r_seq = f"{flank_left_seq}{tr_candidate}"
    r_fwd = parasail.sg_qe_scan_profile_sat(db_seq_profile, ext_r_seq, indel_penalty, indel_penalty)
    r_adj = r_fwd.end_query + 1 - len(flank_left_seq) - ref_size  # Amount to tweak boundary on the right side by

    ext_l_seq = f"{tr_candidate}{flank_right_seq[max(r_adj, 0):]}"[::-1]  # reverse

    r_rev = parasail.sg_qe_scan_profile_sat(db_seq_rev_profile, ext_l_seq, indel_penalty, indel_penalty)
    l_adj = r_rev.end_query + 1 - len(flank_right_seq) - ref_size  # Amount to tweak boundary on the left side by

    return (r_fwd.score, r_adj), (r_rev.score, l_adj)


# TODO: instead of lru_cache, some more custom mechanism for sharing?
@lru_cache(maxsize=512)
def get_repeat_count(
    start_count: int,
    tr_seq: str,
    flank_left_seq: str,
    flank_right_seq: str,
    motif: str,
    max_iters: int,
    local_search_range: int = DEFAULT_LOCAL_SEARCH_RANGE,  # TODO: Parametrize for user
    step_size: int = 1,
) -> tuple[tuple[int, int], int, int]:
    return _get_repeat_count(
        start_count, tr_seq, flank_left_seq, flank_right_seq, motif, max_iters, local_search_range, step_size
    )


def get_ref_repeat_count(
    start_count: int,
    tr_seq: str,
    flank_left_seq: str,
    flank_right_seq: str,
    motif: str,
    ref_size: int,
    vcf_anchor_size: int,
    max_iters: int,
    respect_coords: bool = False,
    local_search_range: int = DEFAULT_LOCAL_SEARCH_RANGE,  # TODO: Parametrize for user
    step_size: int = 1,
) -> tuple[tuple[Union[int, float], int], int, int, tuple[int, int], tuple[str, str, str]]:
    l_offset: int = 0
    r_offset: int = 0

    db_seq: str = f"{flank_left_seq}{tr_seq}{flank_right_seq}"
    db_seq_profile: parasail.Profile = parasail.profile_create_sat(db_seq, dna_matrix)
    db_seq_rev_profile: parasail.Profile = parasail.profile_create_sat(db_seq[::-1], dna_matrix)

    motif_size = len(motif)

    n_offset_scores: int = 0

    if not respect_coords:  # Extend out coordinates from initial definition
        to_explore: list[tuple[int, Literal[-1, 0, 1]]] = [
            (start_count - step_size, -1), (start_count + step_size, 1), (start_count, 0)]

        fwd_sizes_scores_adj: dict[Union[int, float], tuple[int, int]] = {}
        rev_sizes_scores_adj: dict[Union[int, float], tuple[int, int]] = {}

        while to_explore and n_offset_scores < max_iters:
            size_to_explore, direction = to_explore.pop()
            if size_to_explore < 0:
                continue

            fwd_scores: list[tuple[Union[float, int], tuple[int, int], int]] = []  # For right-side adjustment
            rev_scores: list[tuple[Union[float, int], tuple[int, int], int]] = []  # For left-side adjustment

            start_size = max(
                size_to_explore - (local_search_range if (direction < 1 or step_size > local_search_range) else 0), 0)
            end_size = size_to_explore + (local_search_range if (direction > -1 or step_size > local_search_range)
                                          else 0)

            for i in range(start_size, end_size + 1):
                fwd_rs = fwd_sizes_scores_adj.get(i)
                rev_rs = rev_sizes_scores_adj.get(i)

                if fwd_rs is None or rev_rs is None:
                    res = score_ref_boundaries(
                        db_seq_profile, db_seq_rev_profile, motif * i, flank_left_seq, flank_right_seq, ref_size)

                    fwd_sizes_scores_adj[i] = fwd_rs = res[0]
                    rev_sizes_scores_adj[i] = rev_rs = res[1]

                    n_offset_scores += 1

                fwd_scores.append((i, fwd_rs, i))
                rev_scores.append((i, rev_rs, i))

            mv: tuple[Union[float, int], tuple[int, int], int] = max((*fwd_scores, *rev_scores), key=idx_1_getter)
            if mv[2] > size_to_explore and (
                    (new_rc := mv[2] + step_size) not in fwd_sizes_scores_adj or new_rc not in rev_sizes_scores_adj):
                if new_rc >= 0:
                    to_explore.append((new_rc, 1))
            if mv[2] < size_to_explore and (
                    (new_rc := mv[2] - step_size) not in fwd_sizes_scores_adj or new_rc not in rev_sizes_scores_adj):
                if new_rc >= 0:
                    to_explore.append((new_rc, -1))

        # noinspection PyTypeChecker
        fwd_top_res: tuple[Union[int, float], tuple] = max(fwd_sizes_scores_adj.items(), key=lambda x: x[1][0])
        # noinspection PyTypeChecker
        rev_top_res: tuple[Union[int, float], tuple] = max(rev_sizes_scores_adj.items(), key=lambda x: x[1][0])

        # Ignore negative differences (contractions vs TRF definition), but follow expansions
        # TODO: Should we incorporate contractions? How would that work?

        l_offset = rev_top_res[1][1]
        r_offset = fwd_top_res[1][1]

        if l_offset >= len(flank_left_seq) - vcf_anchor_size:
            # don't do anything weird if we're removing the entire flank sequence
            # TODO: this can be caused by NNNNNNN - see chr5:139453668-139454525 in GRCh38
            l_offset = 0
        if r_offset >= len(flank_right_seq):
            r_offset = 0  # same here

        if l_offset > 0:
            tr_seq = flank_left_seq[-1*l_offset:] + tr_seq  # first, move a chunk of the left flank to the TR seq
            flank_left_seq = flank_left_seq[:-1*l_offset]  # then, remove that chunk from the left flank
        if r_offset > 0:
            tr_seq = tr_seq + flank_right_seq[:r_offset]  # same, but for the right flank
            flank_right_seq = flank_right_seq[r_offset:]

    # ------------------------------------------------------------------------------------------------------------------

    final_res, n_iters_final_count, _ = get_repeat_count(
        # always start with int here:
        round(((start_count * motif_size) + (max(0, l_offset) + max(0, r_offset))) / motif_size),
        tr_seq,
        flank_left_seq,
        flank_right_seq,
        motif,
        max_iters=max_iters,
        step_size=step_size,
    )

    return (
        final_res, l_offset, r_offset, (n_offset_scores, n_iters_final_count), (flank_left_seq, tr_seq, flank_right_seq)
    )
