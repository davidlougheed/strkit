import logging
import numpy as np
import pysam
from collections import Counter

from typing import Literal, Optional

from strkit.logger import logger

from .types import ReadDict, CandidateSNV
from .utils import find_pair_by_ref_pos


__all__ = [
    "SNV_OUT_OF_RANGE_CHAR",
    "SNV_GAP_CHAR",
    "get_candidate_snvs",
    "get_read_snvs",
    "calculate_useful_snvs",
    "call_and_filter_useful_snvs",
    "process_read_snvs_for_locus_and_calculate_useful_snvs",
]

SNV_OUT_OF_RANGE_CHAR = "-"
SNV_GAP_CHAR = "_"


def _human_chrom_to_refseq_accession(contig: str, snv_vcf_contigs: tuple[str, ...]) -> Optional[str]:
    contig = contig.removeprefix("chr")
    if contig == "X":
        contig = "23"
    if contig == "Y":
        contig = "24"
    if contig == "M":
        contig = "12920"
    contig = f"NC_{contig.zfill(6)}"

    for vcf_contig in snv_vcf_contigs:
        if vcf_contig.startswith(contig):
            return vcf_contig

    return None


def get_candidate_snvs(
    snv_vcf_file: pysam.VariantFile,
    snv_vcf_contigs: tuple[str, ...],
    snv_vcf_file_format: Literal["chr", "num", "acc", ""],
    contig: str,
    left_most_coord: int,
    right_most_coord: int,
) -> dict[int, CandidateSNV]:
    candidate_snvs_dict: dict[int, CandidateSNV] = {}  # Lookup dictionary for candidate SNVs by position

    snv_contig: str = contig
    if snv_contig not in snv_vcf_contigs:
        if snv_vcf_file_format == "num":
            snv_contig = snv_contig.removeprefix("chr")
        elif snv_vcf_file_format == "acc":
            snv_contig = _human_chrom_to_refseq_accession(snv_contig, snv_vcf_contigs)
        # Otherwise, leave as-is

    for snv in snv_vcf_file.fetch(snv_contig, left_most_coord, right_most_coord + 1):
        snv_ref = snv.ref
        snv_alts = snv.alts
        # check actually is SNV
        if snv_ref is not None and len(snv_ref) == 1 and snv_alts and any(len(a) == 1 for a in snv_alts):
            # Convert from 1-based indexing to 0-based indexing!!!
            # See https://pysam.readthedocs.io/en/latest/api.html#pysam.VariantRecord.pos
            candidate_snvs_dict[snv.pos - 1] = CandidateSNV(id=snv.id, ref=snv.ref, alts=snv_alts)

    # candidate_snvs_dict_items: list[tuple[int, CandidateSNV]] = list(candidate_snvs_dict.items())
    # This flattened version is useful for passing to the Rust extension
    # candidate_snvs_dict_items_flat: list[tuple[int, str, str, list[str]]] = [
    #     (k, v["id"], v["ref"], list(v["alts"])) for k, v in candidate_snvs_dict_items]

    return candidate_snvs_dict


from strkit_rust_ext import get_read_snvs, process_read_snvs_for_locus_and_calculate_useful_snvs


def _find_base_at_pos(query_sequence: str, pairs_for_read: tuple[list[int], list[int]], t: int,
                      start_left: int = 0) -> tuple[str, int]:
    idx, found = find_pair_by_ref_pos(pairs_for_read[1], t, start_left=start_left)

    if found:
        # Even if not in SNV set, it is not guaranteed to be a reference base, since
        # it's possible it was surrounded by too much other variation during the original
        # SNV getter algorithm.
        return query_sequence[pairs_for_read[0][idx]], idx

    # Nothing found, so must have been a gap
    return SNV_GAP_CHAR, idx


def calculate_useful_snvs(
    n_reads: int,
    read_dict_items: tuple[tuple[str, ReadDict], ...],
    read_dict_extra: dict[str, dict],
    read_match_pairs: dict[str, tuple[list[int], list[int]]],
    locus_snvs: set[int],
    min_allele_reads: int,
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
        last_pair_idx: int = 0

        for snv_pos in sorted_snvs:
            base: str = SNV_OUT_OF_RANGE_CHAR
            if segment_start <= snv_pos <= segment_end:
                if bb := snvs.get(snv_pos):
                    base = bb
                else:
                    # Binary search for base from correct pair
                    # - We go in order, so we don't have to search left of the last pair index we tried.
                    base, pair_idx = _find_base_at_pos(qs, pairs_for_read, snv_pos, start_left=last_pair_idx)
                    last_pair_idx = pair_idx

            # Otherwise, leave as out-of-range

            snv_list.append(base)

            if base != SNV_OUT_OF_RANGE_CHAR and base != SNV_GAP_CHAR:
                # Only count SNV bases where we've actually found the base for the read.
                snv_counters[snv_pos][base] += 1

        extra_data["snv_bases"] = tuple(snv_list)

    # Enough reads to try for SNV based separation

    # require 2 alleles for the SNV, both with at least 1/5 of the reads, in order to differentiate alleles.
    # require over ~55% of the reads to have the SNV; otherwise it becomes too easy to occasionally get cases with
    # disjoint sets of SNVs.

    useful_snvs: list[int] = []
    allele_read_threshold: int = max(round(n_reads / 5), min_allele_reads)  # TODO: parametrize proportion
    total_read_threshold: int = max(round(n_reads * 0.55), 5)  # TODO: parametrize

    # snv_counters is guaranteed by the previous inner loop to not have SNV_OUT_OF_RANGE_CHAR or SNV_GAP_CHAR

    for s_idx, (_, snv_counter) in enumerate(snv_counters.items()):
        n_alleles_meeting_threshold: int = sum(
            1 for n_reads_for_allele in snv_counter.values() if n_reads_for_allele >= allele_read_threshold)
        n_reads_with_this_snv_called: int = snv_counter.total()
        if n_alleles_meeting_threshold >= 2 and n_reads_with_this_snv_called >= total_read_threshold:
            useful_snvs.append(s_idx)

    return [(s_idx, sorted_snvs[s_idx]) for s_idx in useful_snvs]  # Tuples of (index in STR list, ref position)


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
