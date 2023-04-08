import bisect
import numpy as np
from collections import Counter

from numpy.typing import NDArray
from typing import Callable, Optional

from strkit.logger import logger

from .types import ReadDict


__all__ = [
    "SNV_OUT_OF_RANGE_CHAR",
    "get_read_snvs",
    "calculate_useful_snvs",
    "call_and_filter_useful_snvs",
]

SNV_OUT_OF_RANGE_CHAR = "-"

# TODO: annotate with rsID if file provided


def _get_read_snvs_meticulous(
    query_sequence: str,
    pairs: list[tuple[int, int]],
    ref_seq: str,
    ref_coord_start: int,
    tr_start_pos: int,
    tr_end_pos: int,
    contiguous_threshold: int = 5,
    max_snv_group_size: int = 5,
) -> dict[int, str]:
    """
    Given a list of tuples of aligned (read pos, ref pos) pairs, this function finds non-reference SNVs which are
    surrounded by a stretch of aligned bases of a specified size on either side.
    :return: Dictionary of {position: base}
    """

    snvs: dict[int, str] = {}

    lhs_contiguous: int = 0
    rhs_contiguous: int = 0
    last_rp: int = -1

    snv_group: list[tuple[int, str]] = []

    for read_pos, ref_pos in pairs:
        if tr_start_pos <= ref_pos < tr_end_pos:  # base is in the tandem repeat itself; skip it
            continue

        read_base = query_sequence[read_pos]
        ref_base = ref_seq[ref_pos - ref_coord_start]

        if read_base == ref_base and (ref_pos - last_rp == 1 or last_rp == -1):
            if snv_group:
                rhs_contiguous += 1
            else:
                lhs_contiguous += 1

            if lhs_contiguous > contiguous_threshold and rhs_contiguous > contiguous_threshold:
                if len(snv_group) <= max_snv_group_size:
                    snvs.update(snv_group)
                # Otherwise, it might be a little mismapped area or a longer deletion vs reference, so ignore it.
                lhs_contiguous = 0
                rhs_contiguous = 0
                snv_group.clear()

            last_rp = ref_pos
            continue

        if ref_pos - last_rp > 1:
            lhs_contiguous = 0
            last_rp = ref_pos
            continue

        if read_base != ref_base:
            snv_group.append((ref_pos, read_base))
            # Don't reset either contiguous variable; instead, take this as part of a SNP group
            last_rp = ref_pos

    return snvs


def _get_read_snvs_simple_py(
    query_sequence: str,
    pairs: list[tuple[int, int]],
    ref_seq: str,
    ref_coord_start: int,
    tr_start_pos: int,
    tr_end_pos: int,
) -> dict[int, str]:
    snvs: dict[int, str] = {}
    for read_pos, ref_pos in pairs:
        if tr_start_pos <= ref_pos < tr_end_pos:  # base is in the tandem repeat itself; skip it
            continue
        if (read_base := query_sequence[read_pos]) != ref_seq[ref_pos - ref_coord_start]:
            snvs[ref_pos] = read_base
    return snvs


get_read_snvs_simple: Callable[[str, list[tuple[int, int]], str, int, int, int], dict[int, str]]
try:
    from strkit_rust_ext import mk_snvs_dict as get_read_snvs_simple
    logger.debug("Found STRkit Rust component, importing get_read_snvs_simple")
except ImportError:
    get_read_snvs_simple = _get_read_snvs_simple_py


def get_read_snvs(
    query_sequence: str,
    pairs: list[tuple[int, int]],
    ref_seq: str,
    ref_coord_start: int,
    tr_start_pos: int,
    tr_end_pos: int,
    contiguous_threshold: int = 5,
    max_snv_group_size: int = 5,
    too_many_snvs_threshold: int = 150,
) -> dict[int, str]:
    """
    Given a list of tuples of aligned (read pos, ref pos) pairs, this function finds non-reference SNVs which are
    surrounded by a stretch of aligned bases of a specified size on either side.
    :return: Dictionary of {position: base}
    """

    # Tried to vectorize this with numpy, and it ended up slower... oh well

    snvs: dict[int, str] = get_read_snvs_simple(
        query_sequence, pairs, ref_seq, ref_coord_start, tr_start_pos, tr_end_pos)

    if len(snvs) >= too_many_snvs_threshold:  # TOO MANY, some kind of mismapping going on?
        return _get_read_snvs_meticulous(
            query_sequence,
            pairs,
            ref_seq,
            ref_coord_start,
            tr_start_pos,
            tr_end_pos,
            contiguous_threshold,
            max_snv_group_size,
        )

    return snvs


def _find_base_at_pos(query_sequence: str, pairs_for_read: list[tuple[int, int]], n_pairs: int, t: int) -> str:
    idx = bisect.bisect_left(pairs_for_read, t, key=lambda p: p[1])

    if idx != n_pairs and (pair := pairs_for_read[idx])[1] == t:
        # Even if not in SNV set, it is not guaranteed to be a reference base, since
        # it's possible it was surrounded by too much other variation during the original
        # SNV getter algorithm.
        return query_sequence[pair[0]]

    # Nothing found, so must have been a gap
    return "_"


def calculate_useful_snvs(
    n_reads: int,
    read_dict_items: tuple[tuple[str, ReadDict], ...],
    read_dict_extra: dict[str, dict],
    read_match_pairs: dict[str, list[tuple[int, int]]],
    locus_snvs: set[int],
) -> list[tuple[int, int]]:
    sorted_snvs: list[int] = sorted(locus_snvs)
    snv_counters: dict[int, Counter] = {sp: Counter() for sp in sorted_snvs}

    for rn, read in read_dict_items:
        pairs_for_read = read_match_pairs[rn]
        n_pairs_for_read: int = len(pairs_for_read)
        extra_data = read_dict_extra[rn]

        snvs: dict[int, str] = extra_data["snv"]

        # Know this to not be None since we were passed only segments with non-None strings earlier
        qs: str = extra_data["_qs"]

        segment_start: int = extra_data["_ref_start"]
        segment_end: int = extra_data["_ref_end"]

        snv_list: list[str] = []

        for snv_pos in sorted_snvs:
            base: str = SNV_OUT_OF_RANGE_CHAR
            if segment_start <= snv_pos <= segment_end:
                if bb := snvs.get(snv_pos):
                    base = bb
                else:
                    # Binary search for base from correct pair
                    base = _find_base_at_pos(qs, pairs_for_read, n_pairs_for_read, snv_pos)

            # Otherwise, leave as out-of-range

            snv_list.append(base)
            snv_counters[snv_pos][base] += 1

        extra_data["snv_bases"] = tuple(snv_list)

    # Enough reads to try for SNV based separation
    useful_snvs: list[int] = []
    for si, (snv_counted, snv_counter) in enumerate(snv_counters.items()):
        read_threshold = max(round(n_reads / 5), 2)  # TODO: parametrize
        n_alleles_meeting_threshold = 0
        for k in snv_counter:
            if k == SNV_OUT_OF_RANGE_CHAR:
                continue
            if snv_counter[k] >= read_threshold:
                n_alleles_meeting_threshold += 1
        if n_alleles_meeting_threshold >= 2:
            useful_snvs.append(si)

    return [(si, sorted_snvs[si]) for si in useful_snvs]  # Tuples of (index in STR list, ref position)


def call_and_filter_useful_snvs(
    n_alleles: int,
    read_dict: dict[str, ReadDict],
    useful_snvs: list[tuple[int, int]],
    peak_order: NDArray[np.int_],
    locus_log_str: str,
    logger,
) -> list[dict]:
    """
    Call useful SNVs at a locus level from read-level SNV data.
    :param n_alleles: The number of alleles called for this locus.
    :param read_dict: Dictionary of read data. Must already have peaks assigned.
    :param useful_snvs: List of tuples representing useful SNVs: (SNV index, reference position)
    :param peak_order: Indices for rearranging the call arrays into the final order, to match the sorted copy numbers.
    :param locus_log_str: Locus string representation for logging purposes.
    :param logger: Python logger object.
    :return: List of called SNVs for the locus.
    """

    # Since these have already been classified as 'useful' earlier in the pipeline,
    # we have some guarantees that these values should be fairly internally consistent
    # for a given peak... most of the time.

    allele_range = tuple(range(n_alleles))
    peak_base_counts: dict[int, dict[int, Counter]] = {}

    for _, u_ref in useful_snvs:
        peak_base_counts[u_ref] = {p: Counter() for p in allele_range}

    for rn, read in read_dict.items():
        p: Optional[int] = read.get("p")
        if p is None:
            logger.debug(f"call_useful_snvs: no peak found for read {rn}")
            continue
        for u_idx, (_, u_ref) in enumerate(useful_snvs):
            peak_base_counts[u_ref][p].update((read["snvu"][u_idx],))

    called_snvs: list[dict] = []
    skipped_snvs: set[int] = set()

    for u_idx, (u_ref, peak_counts) in enumerate(peak_base_counts.items()):
        call: list[str] = []
        rs: list[int] = []

        skipped: bool = False

        for a in allele_range:
            mc = peak_counts[a].most_common(2)
            mcc = mc[0]
            try:
                if mcc[0] == SNV_OUT_OF_RANGE_CHAR:  # Chose most common non-uncalled value
                    mcc = mc[1]
            except IndexError:  # '-' is the only value, somehow
                logger.warn(
                    f"{locus_log_str} - for SNV {u_ref}, found only '{SNV_OUT_OF_RANGE_CHAR}' with {mcc[1]} reads")
                logger.debug(f"{locus_log_str} - for SNV {u_ref}: {mc=}, {peak_counts[a]=}")
                skipped = True
                break

            call.append(mcc[0])
            rs.append(mcc[1])

        if skipped:
            skipped_snvs.add(u_idx)  # Skip this useful SNV, since it isn't actually useful
            continue

        called_snvs.append({
            "pos": u_ref,
            "call": np.array(call)[peak_order].tolist(),
            "rs": np.array(rs)[peak_order].tolist(),
        })

    # If we've skipped any SNVs, filter them out of the read dict - MUTATION
    if skipped_snvs:
        for read in read_dict.values():
            read["snvu"] = tuple(b for i, b in enumerate(read["snvu"]) if i not in skipped_snvs)
        logger.debug(f"{locus_log_str} - filtered out {len(skipped_snvs)} not-actually-useful SNVs")

    return called_snvs
