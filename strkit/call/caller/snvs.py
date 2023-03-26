import numpy as np
import pysam
from collections import Counter

from numpy.typing import NDArray
from typing import Optional

from .types import ReadDict


__all__ = [
    "get_read_snvs",
    "call_useful_snvs",
]

# TODO: annotate with rsID if file provided


def get_read_snvs(
    query_sequence: str,
    pairs: list[tuple[int, int], ...],
    contig: str,
    ref: pysam.FastaFile,
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

    fm_qp, fm_rp = pairs[0]
    lm_qp, lm_rp = pairs[-1]

    query_sequence = query_sequence.upper()
    ref_sequence: str = ref.fetch(contig, fm_rp, lm_rp + 1).upper()

    lhs_contiguous: int = 0
    rhs_contiguous: int = 0
    last_rp: int = -1

    snv_group: list[tuple[int, str]] = []

    for read_pos, ref_pos in pairs:
        if tr_start_pos <= ref_pos < tr_end_pos:  # base is in the tandem repeat itself; skip it
            continue

        read_base = query_sequence[read_pos]
        ref_base = ref_sequence[ref_pos - fm_rp]

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
                if mcc[0] == "-":  # Chose most common non-uncalled value
                    mcc = mc[1]
            except IndexError:  # - is the only value, somehow
                logger.warn(f"{locus_log_str} - for SNV {u_ref}, found only '-' with {mcc[1]} reads")
                pass  # TODO: should we set mcc[1] to 0 here?
            call.append(mcc[0])
            rs.append(mcc[1])

        called_snvs.append({
            "pos": u_ref,
            "call": np.array(call)[peak_order].tolist(),
            "rs": np.array(rs)[peak_order].tolist(),
        })

    return called_snvs
