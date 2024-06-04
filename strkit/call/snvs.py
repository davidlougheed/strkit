import logging
import multiprocessing.managers as mmg
import threading

from collections import Counter
from typing import Optional

from strkit_rust_ext import get_read_snvs, process_read_snvs_for_locus_and_calculate_useful_snvs, CandidateSNVs

from strkit.logger import logger
from .types import ReadDict, CalledSNV
from .utils import idx_1_getter


__all__ = [
    "SNV_OUT_OF_RANGE_CHAR",
    "SNV_GAP_CHAR",
    "get_read_snvs",
    "call_and_filter_useful_snvs",
    "process_read_snvs_for_locus_and_calculate_useful_snvs",
]

SNV_OUT_OF_RANGE_CHAR = "-"
SNV_GAP_CHAR = "_"


def call_and_filter_useful_snvs(
    contig: str,
    n_alleles: int,
    read_dict: dict[str, ReadDict],
    useful_snvs: list[tuple[int, int]],
    candidate_snvs_dict: CandidateSNVs,
    # ---
    snv_genotype_update_lock: threading.Lock,
    snv_genotype_cache: mmg.DictProxy,
    # ---
    locus_log_str: str,
    logger_: logging.Logger,
) -> list[CalledSNV]:
    """
    Call useful SNVs at a locus level from read-level SNV data.
    :param contig: The contig of the SNVs. Used for generating an ID if one does not exist.
    :param n_alleles: The number of alleles called for this locus.
    :param read_dict: Dictionary of read data. Must already have peaks assigned.
    :param useful_snvs: List of tuples representing useful SNVs: (SNV index, reference position)
    :param candidate_snvs_dict: A dictionary of useful SNVs, indexed by reference position. Used to look up IDs.
    :param snv_genotype_update_lock: Lock for updating/querying the SNV genotype/phase set cache for a long time.
    :param snv_genotype_cache: Cache for SNV genotype/phase set information.
    :param locus_log_str: Locus string representation for logging purposes.
    :param logger_: Python logger object.
    :return: List of called SNVs for the locus.
    """

    # Since these have already been classified as 'useful' earlier in the pipeline,
    # we have some guarantees that these values should be fairly internally consistent
    # for a given peak... most of the time.

    allele_range = tuple(range(n_alleles))
    peak_base_counts: dict[int, dict[int, Counter]] = {
        u_ref: {p: Counter() for p in allele_range}
        for _, u_ref in useful_snvs
    }

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
            if skipped:
                break

            mc = peak_counts[a].most_common(2)
            mcc = mc[0]

            try:
                if mcc[0] == SNV_OUT_OF_RANGE_CHAR:  # Chose most common non-uncalled value
                    mcc = mc[1]

                for b in allele_range:
                    if b == a:
                        continue
                    if (peak_counts[b][mcc[0]] / peak_counts[b].total()) > (
                            peak_counts[a][mcc[0]] / peak_counts[a].total() / 2):  # TODO: parametrize
                        logger_.debug(
                            f"{locus_log_str} - for SNV position {u_ref}: got uninformative peak counts (cross-talk) - "
                            f"{peak_counts=}")
                        skipped = True
                        break

            except IndexError:  # '-' is the only value, somehow
                logger_.debug(
                    f"{locus_log_str} - for SNV {u_ref}, found only '{SNV_OUT_OF_RANGE_CHAR}' with {mcc[1]} reads")
                logger_.debug(f"{locus_log_str} - for SNV position {u_ref}: {mc=}, {peak_counts[a]=}")
                skipped = True
                break

            call.append(mcc[0])
            rs.append(mcc[1])

        snv_call_set = set(call)

        if not skipped and len(snv_call_set) == 1:
            # print(u_idx, u_ref, peak_counts, call, rs)
            logger_.warning(f"{locus_log_str} - for SNV position {u_ref}: got degenerate call {call} from {peak_counts=}")
            skipped = True

        snv_rec = candidate_snvs_dict.get(u_ref)
        if snv_rec is not None:
            snv_id = snv_rec["id"]
            if snv_id == ".":
                snv_id = f"{contig}_{u_ref}"
        else:
            snv_id = f"{contig}_{u_ref}"

        if not skipped:
            snv_genotype_update_lock.acquire(timeout=600)
            if snv_id in snv_genotype_cache and (cgt := set(snv_genotype_cache[snv_id][0])) != snv_call_set:
                logger_.warning(
                    f"{locus_log_str} - got mismatch for SNV {snv_id} (position {u_ref}); cache genotype set {cgt} != "
                    f"current genotype set {snv_call_set}")
                skipped = True
            snv_genotype_update_lock.release()

        if skipped:
            skipped_snvs.add(u_idx)  # Skip this useful SNV, since it isn't actually useful
            continue

        called_snvs.append({
            "id": snv_id,
            **({"ref": snv_rec["ref_base"]} if snv_rec is not None else {}),
            "pos": u_ref,
            "call": call,
            "rcs": rs,
        })

    # If we've skipped any SNVs, filter them out of the read dict - MUTATION
    if skipped_snvs:
        for read in read_dict.values():
            if "snvu" not in read:
                continue
            read["snvu"] = tuple(map(idx_1_getter, filter(lambda e: e[0] not in skipped_snvs, enumerate(read["snvu"]))))
        logger.debug(f"{locus_log_str} - filtered out {len(skipped_snvs)} not-actually-useful SNVs")

    return called_snvs
