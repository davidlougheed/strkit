import logging
import multiprocessing.managers as mmg
import numpy as np
import pysam
import threading

from collections import Counter
from typing import Literal, Optional

from strkit_rust_ext import get_read_snvs, process_read_snvs_for_locus_and_calculate_useful_snvs

from strkit.logger import logger
from .types import ReadDict, CandidateSNV


__all__ = [
    "SNV_OUT_OF_RANGE_CHAR",
    "SNV_GAP_CHAR",
    "get_candidate_snvs",
    "get_read_snvs",
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


def call_and_filter_useful_snvs(
    contig: str,
    n_alleles: int,
    read_dict: dict[str, ReadDict],
    useful_snvs: list[tuple[int, int]],
    candidate_snvs_dict: dict[int, CandidateSNV],
    # ---
    snv_genotype_update_lock: threading.Lock,
    snv_genotype_cache: mmg.DictProxy,
    # ---
    locus_log_str: str,
    logger_: logging.Logger,
) -> list[dict]:
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

        if not skipped and len(set(call)) == 1:
            # print(u_idx, u_ref, peak_counts, call, rs)
            logger_.warning(f"{locus_log_str} - for SNV position {u_ref}: got degenerate call {call} from {peak_counts=}")
            skipped = True

        snv_rec = candidate_snvs_dict.get(u_ref)
        snv_id = snv_rec["id"] if snv_rec is not None else f"{contig}_{u_ref}"
        snv_call = np.array(call).tolist()

        if not skipped:
            snv_genotype_update_lock.acquire(timeout=30)
            if snv_id in snv_genotype_cache and (cgt := set(snv_genotype_cache[snv_id][0])) != (sgt := set(snv_call)):
                logger_.warning(
                    f"{locus_log_str} - got mismatch for SNV {snv_id} (position {u_ref}); cache genotype set {cgt} != "
                    f"current genotype set {sgt}")
                skipped = True
            snv_genotype_update_lock.release()

        if skipped:
            skipped_snvs.add(u_idx)  # Skip this useful SNV, since it isn't actually useful
            continue

        called_snvs.append({
            "id": snv_id,
            **({"ref": snv_rec["ref"]} if snv_rec is not None else {}),
            "pos": u_ref,
            "call": snv_call,
            "rcs": np.array(rs).tolist(),
        })

    # If we've skipped any SNVs, filter them out of the read dict - MUTATION
    if skipped_snvs:
        for read in read_dict.values():
            read["snvu"] = tuple(b for i, b in enumerate(read["snvu"]) if i not in skipped_snvs)
        logger.debug(f"{locus_log_str} - filtered out {len(skipped_snvs)} not-actually-useful SNVs")

    return called_snvs
