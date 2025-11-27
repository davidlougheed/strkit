from __future__ import annotations

from collections import Counter
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from logging import Logger
    # noinspection PyProtectedMember
    from multiprocessing.managers import DictProxy

from strkit_rust_ext import process_read_snvs_for_locus_and_calculate_useful_snvs, CandidateSNVs
from strkit.utils import idx_1_getter

from .types import ReadDict, CalledSNV


__all__ = [
    "SNV_OUT_OF_RANGE_CHAR",
    "SNV_GAP_CHAR",
    "SNV_NA_CHARS",
    "call_and_filter_useful_snvs",
    "process_read_snvs_for_locus_and_calculate_useful_snvs",
]

SNV_OUT_OF_RANGE_CHAR = "-"
SNV_GAP_CHAR = "_"
SNV_NA_CHARS = (SNV_OUT_OF_RANGE_CHAR, SNV_GAP_CHAR)


def call_and_filter_useful_snvs(
    contig: str,
    n_alleles: int,
    read_dict: dict[str, ReadDict],
    useful_snvs: list[tuple[int, int]],
    candidate_snvs: CandidateSNVs,
    # ---
    snv_quality_threshold: int,
    # ---
    snv_genotype_cache: DictProxy,
    # ---
    locus_log_str: str,
    logger_: Logger,
) -> list[CalledSNV]:
    """
    Call useful SNVs at a locus level from read-level SNV data.
    :param contig: The contig of the SNVs. Used for generating an ID if one does not exist.
    :param n_alleles: The number of alleles called for this locus.
    :param read_dict: Dictionary of read data. Must already have peaks assigned.
    :param useful_snvs: List of tuples representing useful SNVs: (SNV index, reference position)
    :param candidate_snvs: A dictionary of useful SNVs, indexed by reference position. Used to look up IDs.
    :param snv_quality_threshold: Minimum PHRED score needed to incorporate a read base into the genotype.
    :param snv_genotype_cache: Cache for SNV genotype/phase set information.
    :param locus_log_str: Locus string representation for logging purposes.
    :param logger_: Python logger object.
    :return: List of called SNVs for the locus.
    """

    # Since these have already been classified as 'useful' earlier in the pipeline,
    # we have some guarantees that these values should be fairly internally consistent
    # for a given peak... most of the time.

    allele_range: tuple[int, ...] = tuple(range(n_alleles))
    peak_base_counts: dict[int, dict[int, Counter]] = {
        u_ref: {p: Counter() for p in allele_range}
        for _, u_ref in useful_snvs
    }

    for rn, read in read_dict.items():
        p: int | None = read.get("p")
        if p is None:  # No peak; read wasn't used to call peaks
            continue
        for u_idx, (_, u_ref) in enumerate(useful_snvs):
            su, su_q = read["snvu"][u_idx]

            if su == SNV_GAP_CHAR or su_q >= snv_quality_threshold:
                peak_base_counts[u_ref][p].update((su,))

    called_snvs: list[CalledSNV] = []
    skipped_snvs: set[int] = set()

    for u_idx, (u_ref, peak_counts) in enumerate(peak_base_counts.items()):
        call: list[str] = []
        rs: list[int] = []

        skipped: bool = False

        for a in allele_range:
            if skipped:
                break

            peak_counts_a = peak_counts[a]
            a_total = peak_counts[a].total()

            if a_total == 0:  # probably due to quality filtering
                skipped = True
                logger_.warning(f"{locus_log_str} - for SNV {u_ref}, found a 0-total for allele {a} (a)")
                break

            mc = peak_counts_a.most_common(2)
            mcc = mc[0]

            try:
                if mcc[0] == SNV_OUT_OF_RANGE_CHAR:  # Chose most common non-uncalled value
                    mcc = mc[1]

                for b in allele_range:
                    if b == a:
                        continue

                    peak_counts_b = peak_counts[b]
                    b_total = peak_counts_b.total()

                    if b_total == 0:  # probably due to quality filtering
                        skipped = True
                        logger_.warning("%s - for SNV %d, found a 0-total for allele %d (b)", locus_log_str, u_ref, b)
                        break

                    if (peak_counts_b[mcc[0]] / b_total) > (peak_counts_a[mcc[0]] / a_total / 2):  # TODO: parametrize
                        logger_.debug(
                            "%s - for SNV position %d: got uninformative peak counts (cross-talk) - peak_counts=%s",
                            locus_log_str, u_ref, str(peak_counts))
                        skipped = True
                        break

            except IndexError:  # '-' is the only value, somehow
                logger_.debug(
                    f"%s - for SNV %d, found only '{SNV_OUT_OF_RANGE_CHAR}' with %d reads",
                    locus_log_str, u_ref, mcc[1])
                logger_.debug(
                    "%s - for SNV position %d: mc=%s, peak_counts[a]=%s",
                    locus_log_str, u_ref, str(mc), str(peak_counts[a]))
                skipped = True
                break

            if not skipped:
                call.append(mcc[0])
                rs.append(mcc[1])

        snv_call_set = set(call)

        if not skipped and len(snv_call_set) == 1:
            logger_.warning(
                "%s - for SNV position %d: got degenerate call %s from peak_counts=%s",
                locus_log_str, u_ref, call, str(peak_counts))
            skipped = True

        snv_rec = candidate_snvs.get(u_ref)
        if snv_rec is not None:
            snv_id = snv_rec["id"]
            if snv_id == ".":
                snv_id = f"{contig}_{u_ref}"
        else:
            snv_id = f"{contig}_{u_ref}"

        if not skipped:
            cached_snv_genotype = snv_genotype_cache.get(snv_id)
            if cached_snv_genotype is not None and (cgt := set(cached_snv_genotype[0])) != snv_call_set:
                logger_.warning(
                    f"%s - got mismatch for SNV %s (position %d); cache genotype set {cgt} != "
                    f"current genotype set {snv_call_set}", locus_log_str, snv_id, u_ref)
                skipped = True

        if skipped:
            skipped_snvs.add(u_idx)  # Skip this useful SNV, since it isn't actually useful
            continue

        called_snvs.append(
            CalledSNV(
                id=snv_id,
                pos=u_ref,
                call=tuple(call),
                rcs=rs,
                **({"ref": snv_rec["ref_base"]} if snv_rec is not None else {}),
            )
        )

    # If we've skipped any SNVs, filter them out of the read dict - MUTATION
    if skipped_snvs:
        for read in read_dict.values():
            if "snvu" not in read:
                continue
            read["snvu"] = tuple(map(idx_1_getter, filter(lambda e: e[0] not in skipped_snvs, enumerate(read["snvu"]))))
        logger_.debug("%s - filtered out %d not-actually-useful SNVs", locus_log_str, len(skipped_snvs))

    return called_snvs
