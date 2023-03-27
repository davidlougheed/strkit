import numpy as np
import pysam
from collections import Counter

from numpy.typing import NDArray
from typing import Optional

from .types import ReadDict


__all__ = [
    "SNV_OUT_OF_RANGE_CHAR",
    "get_read_snvs",
    "calculate_useful_snvs",
    "call_useful_snvs",
]

SNV_OUT_OF_RANGE_CHAR = "-"

# TODO: annotate with rsID if file provided


def _get_read_snvs_meticulous(
    query_sequence: str,
    pairs: list[tuple[int, int], ...],
    ref_seq: str,
    ref_coord_start: int,
    tr_start_pos: int,
    tr_end_pos: int,
    contiguous_threshold: int = 5,
    max_snv_group_size: int = 5,
):
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


def get_read_snvs(
    query_sequence: str,
    pairs: list[tuple[int, int], ...],
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

    snvs: dict[int, str] = {}

    for read_pos, ref_pos in pairs:
        if tr_start_pos <= ref_pos < tr_end_pos:  # base is in the tandem repeat itself; skip it
            continue
        if (read_base := query_sequence[read_pos]) != ref_seq[ref_pos - ref_coord_start]:
            snvs[ref_pos] = read_base

    if len(snvs) >= too_many_snvs_threshold:  # TOO MANY, some kind of mismapping going on
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


def calculate_useful_snvs(
    n_reads: int,
    read_dict: dict[str, ReadDict],
    read_dict_items: tuple[tuple[str, ReadDict], ...],
    read_match_pairs: dict[str, list[tuple[int, int]]],
    locus_snvs: set[int],
) -> list[tuple[int, int]]:
    sorted_snvs: list[int] = sorted(locus_snvs)
    snv_counters: dict[int, Counter] = {sp: Counter() for sp in sorted_snvs}

    for rn, read in read_dict_items:
        snvs: dict[int, str] = read_dict[rn]["snv"]

        # Know this to not be None since we were passed only segments with non-None strings earlier
        qs: str = read["_qs"]

        segment_start: int = read["_ref_start"]
        segment_end: int = read["_ref_end"]

        snv_list: list[str] = []

        for snv_pos in sorted_snvs:
            base: str = SNV_OUT_OF_RANGE_CHAR
            if segment_start <= snv_pos <= segment_end:
                if bb := snvs.get(snv_pos):
                    base = bb
                else:
                    # Binary search for pair set
                    pairs_for_read = read_match_pairs[rn]

                    def _bin_search() -> str:
                        lhs: int = 0
                        rhs: int = len(pairs_for_read) - 1

                        while lhs <= rhs:
                            pivot: int = (lhs + rhs) // 2
                            pair: tuple[int, int] = pairs_for_read[pivot]
                            if pair[1] < snv_pos:
                                lhs = pivot + 1
                            elif pair[1] > snv_pos:  # pair[1] > snv_pos
                                rhs = pivot - 1
                            else:
                                # Even if not in SNV set, it is not guaranteed to be a reference base, since
                                # it's possible it was surrounded by too much other variation during the original
                                # SNV getter algorithm.
                                return qs[pair[0]]

                        # Nothing found, so must have been a gap
                        return "_"

                    base = _bin_search()

            # Otherwise, leave as gap

            snv_list.append(base)
            snv_counters[snv_pos][base] += 1

        read_dict[rn]["snv_bases"] = tuple(snv_list)

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


def call_useful_snvs(
    n_alleles: int,
    read_dict: dict[str, ReadDict],
    useful_snvs: list[tuple[int, int]],
    peak_order: NDArray[int],
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
            continue
        for u_idx, (_, u_ref) in enumerate(useful_snvs):
            peak_base_counts[u_ref][p].update((read["snvu"][u_idx],))

    called_snvs: list[dict] = []

    for u_ref, peak_counts in peak_base_counts.items():
        call: list[str] = []
        rs: list[int] = []

        for a in allele_range:
            mc = peak_counts[a].most_common(2)
            mcc = mc[0]
            try:
                if mcc[0] == SNV_OUT_OF_RANGE_CHAR:  # Chose most common non-uncalled value
                    mcc = mc[1]
            except IndexError:  # - is the only value, somehow
                logger.warn(
                    f"{locus_log_str} - for SNV {u_ref}, found only '{SNV_OUT_OF_RANGE_CHAR}' with {mcc[1]} reads")
                logger.debug(f"{locus_log_str} - for SNV {u_ref}: {mc=}, {peak_counts[a]=}")
                pass  # TODO: should we set mcc[1] to 0 here?
            call.append(mcc[0])
            rs.append(mcc[1])

        called_snvs.append({
            "pos": u_ref,
            "call": np.array(call)[peak_order].tolist(),
            "rs": np.array(rs)[peak_order].tolist(),
        })

    return called_snvs
