import logging
import math
import numpy as np
from collections import Counter

from typing import Callable, Optional

from strkit.logger import logger

from .types import ReadDict, CandidateSNV
from .utils import find_pair_by_ref_pos


__all__ = [
    "SNV_OUT_OF_RANGE_CHAR",
    "get_read_snvs_dbsnp",
    "get_read_snvs",
    "calculate_useful_snvs",
    "call_and_filter_useful_snvs",
]

SNV_OUT_OF_RANGE_CHAR = "-"
SNV_GAP_CHAR = "_"

# TODO: annotate with rsID if file provided


# We check entropy against a threshold in order to make sure the SNVs we find are somewhat
# useful, i.e., surrounded by a nice mixture of different bases rather than some possibly
# mis-mappable base inside a homopolymer (I used to get results like AAAAAAAAAACAAAAAAAA
# where one of those As could be a gap instead... really useless stuff that made it
# through the filters.)
#
# Below is a Python implementation of https://en.wikipedia.org/wiki/Entropy_(information_theory)
# but a Rust one exists too in strkit_rust_ext, which gets used by the Rust SNV-finding implementations.

def shannon_entropy(seq: str) -> float:
    seq_len = len(seq)
    return -1.0 * sum(p * math.log2(p) for p in (c / seq_len for c in Counter(seq).values()))


def _get_read_snvs_meticulous_py(
    query_sequence: str,
    pairs: list[tuple[int, int]],
    ref_seq: str,
    ref_coord_start: int,
    tr_start_pos: int,
    tr_end_pos: int,
    contiguous_threshold: int,
    max_snv_group_size: int,
    entropy_flank_size: int,
    entropy_threshold: float,
) -> dict[int, str]:
    """
    Given a list of tuples of aligned (read pos, ref pos) pairs, this function finds non-reference SNVs which are
    surrounded by a stretch of aligned bases of a specified size on either side.
    :return: Dictionary of {position: base}
    """

    query_sequence_len = len(query_sequence)

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

        if read_base == ref_base and (contiguous_threshold == 0 or ref_pos - last_rp == 1 or last_rp == -1):
            if snv_group:
                rhs_contiguous += 1
            else:
                lhs_contiguous += 1

            if lhs_contiguous >= contiguous_threshold and rhs_contiguous >= contiguous_threshold:
                if len(snv_group) <= max_snv_group_size:
                    snvs.update(snv_group)
                # Otherwise, it might be a little mismapped area or a longer deletion vs reference, so ignore it.
                lhs_contiguous = 0
                rhs_contiguous = 0
                snv_group.clear()

            last_rp = ref_pos
            continue

        if ref_pos - last_rp > 1:  # Non-contiguous jump; insertion in query
            lhs_contiguous = 0
            last_rp = ref_pos
            continue

        if read_base != ref_base:
            seq = query_sequence[max(read_pos-entropy_flank_size, 0):min(
                read_pos+entropy_flank_size, query_sequence_len)]
            if shannon_entropy(seq) >= entropy_threshold:
                snv_group.append((ref_pos, read_base))
            # Don't reset either contiguous variable; instead, take this as part of a SNP group
            last_rp = ref_pos

    # Special case: if we have stuff in the SNV group with no contiguous requirements,
    # add it to the SNV dict.
    if contiguous_threshold == 0 and (0 < len(snv_group) <= max_snv_group_size):
        snvs.update(snv_group)

    return snvs


def _get_read_snvs_simple_py(
    query_sequence: str,
    pairs: list[tuple[int, int]],
    ref_seq: str,
    ref_coord_start: int,
    tr_start_pos: int,
    tr_end_pos: int,
    entropy_flank_size: int,
    entropy_threshold: float,
) -> dict[int, str]:
    query_sequence_len = len(query_sequence)
    snvs: dict[int, str] = {}
    for read_pos, ref_pos in pairs:
        if tr_start_pos <= ref_pos < tr_end_pos:  # base is in the tandem repeat itself; skip it
            continue
        if (read_base := query_sequence[read_pos]) != ref_seq[ref_pos - ref_coord_start]:
            seq = query_sequence[max(read_pos - entropy_flank_size, 0):min(
                read_pos + entropy_flank_size, query_sequence_len)]
            if shannon_entropy(seq) >= entropy_threshold:
                snvs[ref_pos] = read_base
    return snvs


get_read_snvs_meticulous: Callable[
    [str, list[tuple[int, int]], str, int, int, int, int, int, int, float], dict[int, str]]
try:
    from strkit_rust_ext import get_snvs_meticulous as get_read_snvs_meticulous
    logger.debug("Found STRkit Rust component, importing get_read_snvs_meticulous")
except ImportError:
    get_read_snvs_meticulous = _get_read_snvs_meticulous_py


get_read_snvs_simple: Callable[[str, list[tuple[int, int]], str, int, int, int, int, float], dict[int, str]]
try:
    from strkit_rust_ext import get_snvs_simple as get_read_snvs_simple
    logger.debug("Found STRkit Rust component, importing get_read_snvs_simple")
except ImportError:
    get_read_snvs_simple = _get_read_snvs_simple_py


def get_read_snvs_dbsnp(
    candidate_snvs_dict: dict[int, CandidateSNV],
    query_sequence: str,
    pairs: list[tuple[int, int]],
    tr_start_pos: int,
    tr_end_pos: int,
):
    pairs_l = pairs[0][1]
    pairs_r = pairs[-1][1]

    snvs: dict[int, str] = {}

    for pos, c_snv in candidate_snvs_dict.items():
        if pos < pairs_l:
            continue
        if pos > pairs_r:
            break
        if tr_start_pos <= pos <= tr_end_pos:
            continue
        ref = c_snv["ref"]
        alts = c_snv["alts"]
        if ref is not None and len(ref) == 1 and alts and all(len(a) == 1 for a in alts):
            read_base = _find_base_at_pos(query_sequence, pairs, pos)
            if read_base == ref or read_base in alts:
                snvs[pos] = read_base

    return snvs


def get_read_snvs(
    query_sequence: str,
    pairs: list[tuple[int, int]],
    ref_seq: str,
    ref_coord_start: int,
    tr_start_pos: int,
    tr_end_pos: int,
    contiguous_threshold: int = 5,
    max_snv_group_size: int = 5,
    too_many_snvs_threshold: int = 20,
    entropy_flank_size: int = 10,
    entropy_threshold: float = 1.8,
) -> dict[int, str]:
    """
    Given a list of tuples of aligned (read pos, ref pos) pairs, this function finds non-reference SNVs which are
    surrounded by a stretch of aligned bases of a specified size on either side.
    :return: Dictionary of {position: base}
    """

    # Tried to vectorize this with numpy, and it ended up slower... oh well

    snvs: dict[int, str] = get_read_snvs_simple(
        query_sequence,
        pairs,
        ref_seq,
        ref_coord_start,
        tr_start_pos,
        tr_end_pos,
        entropy_flank_size,
        entropy_threshold,
    )

    if len(snvs) >= too_many_snvs_threshold:  # TOO MANY, some kind of mismapping going on?
        return get_read_snvs_meticulous(
            query_sequence,
            pairs,
            ref_seq,
            ref_coord_start,
            tr_start_pos,
            tr_end_pos,
            contiguous_threshold,
            max_snv_group_size,
            entropy_flank_size,
            entropy_threshold,
        )

    return snvs


def _find_base_at_pos(query_sequence: str, pairs_for_read: list[tuple[int, int]], t: int) -> str:
    idx, found = find_pair_by_ref_pos(pairs_for_read, t)

    if found:
        # Even if not in SNV set, it is not guaranteed to be a reference base, since
        # it's possible it was surrounded by too much other variation during the original
        # SNV getter algorithm.
        return query_sequence[pairs_for_read[idx][0]]

    # Nothing found, so must have been a gap
    return SNV_GAP_CHAR


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
        extra_data = read_dict_extra[rn]

        if "snv" not in extra_data:
            continue

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
                    base = _find_base_at_pos(qs, pairs_for_read, snv_pos)

            # Otherwise, leave as out-of-range

            snv_list.append(base)
            snv_counters[snv_pos][base] += 1

        extra_data["snv_bases"] = tuple(snv_list)

    # Enough reads to try for SNV based separation
    useful_snvs: list[int] = []
    allele_read_threshold: int = max(round(n_reads / 5), 2)  # TODO: parametrize
    total_read_threshold: int = max(round(n_reads * 0.55), 5)
    for si, (snv_counted, snv_counter) in enumerate(snv_counters.items()):
        n_alleles_meeting_threshold: int = 0
        n_reads_with_this_snv_called: int = 0
        for k in snv_counter:
            if k in (SNV_OUT_OF_RANGE_CHAR, SNV_GAP_CHAR):
                continue
            n_reads_for_allele = snv_counter[k]
            n_reads_with_this_snv_called += n_reads_for_allele
            if n_reads_for_allele >= allele_read_threshold:
                n_alleles_meeting_threshold += 1
        # require 2 alleles for the SNV, both with at least 1/5 of the reads, in order to differentiate alleles.
        # require over ~55% of the reads to have the SNV; otherwise it becomes too easy to occasionally get cases with
        # disjoint sets of SNVs.
        if n_alleles_meeting_threshold >= 2 and n_reads_with_this_snv_called >= total_read_threshold:
            useful_snvs.append(si)

    return [(si, sorted_snvs[si]) for si in useful_snvs]  # Tuples of (index in STR list, ref position)


def call_and_filter_useful_snvs(
    n_alleles: int,
    read_dict: dict[str, ReadDict],
    useful_snvs: list[tuple[int, int]],
    candidate_snvs_dict: dict[int, CandidateSNV],
    locus_log_str: str,
    logger_: logging.Logger,
) -> list[dict]:
    """
    Call useful SNVs at a locus level from read-level SNV data.
    :param n_alleles: The number of alleles called for this locus.
    :param read_dict: Dictionary of read data. Must already have peaks assigned.
    :param useful_snvs: List of tuples representing useful SNVs: (SNV index, reference position)
    :param candidate_snvs_dict: A dictionary of useful SNVs, indexed by reference position. Used to look up IDs.
    :param locus_log_str: Locus string representation for logging purposes.
    :param logger_: Python logger object.
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
        if p is None:  # No peak; read wasn't used to call peaks
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
                logger_.warning(
                    f"{locus_log_str} - for SNV {u_ref}, found only '{SNV_OUT_OF_RANGE_CHAR}' with {mcc[1]} reads")
                logger_.debug(f"{locus_log_str} - for SNV {u_ref}: {mc=}, {peak_counts[a]=}")
                skipped = True
                break

            call.append(mcc[0])
            rs.append(mcc[1])

        if skipped:
            skipped_snvs.add(u_idx)  # Skip this useful SNV, since it isn't actually useful
            continue

        snv_rec = candidate_snvs_dict.get(u_ref)
        called_snvs.append({
            **({"id": snv_rec["id"], "ref": snv_rec["ref"]} if snv_rec is not None else {}),
            "pos": u_ref,
            "call": np.array(call).tolist(),
            "rcs": np.array(rs).tolist(),
        })

    # If we've skipped any SNVs, filter them out of the read dict - MUTATION
    if skipped_snvs:
        for read in read_dict.values():
            read["snvu"] = tuple(b for i, b in enumerate(read["snvu"]) if i not in skipped_snvs)
        logger.debug(f"{locus_log_str} - filtered out {len(skipped_snvs)} not-actually-useful SNVs")

    return called_snvs
