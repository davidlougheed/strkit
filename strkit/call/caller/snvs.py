import pysam


__all__ = [
    "get_read_snvs",
]

# vcf = pysam.VariantFile("./00-common_all.vcf.gz")


def get_read_snvs(
    query_sequence: str,
    pairs: list[tuple[int, int], ...],
    contig: str,
    ref: pysam.FastaFile,
    tr_start_pos: int,
    tr_end_pos: int,
    contiguous_threshold: int = 5,
    max_snv_group_size: int = 5,
):
    """
    Given a list of tuples of aligned (read pos, ref pos) pairs, this function finds non-reference SNVs which are
    surrounded by a stretch of aligned bases of a specified size on either side.
    :return:
    """

    snvs = {}

    fm_qp, fm_rp = pairs[0]
    lm_qp, lm_rp = pairs[-1]

    ref_sequence = ref.fetch(contig, fm_rp, lm_rp + 1)

    lhs_contiguous = 0
    rhs_contiguous = 0
    last_rp = -1

    snv_group = []

    for qp, rp in pairs:
        if rp is None:
            continue  # insertion in read; skip this position

        if tr_start_pos <= rp < tr_end_pos:  # base is in the tandem repeat itself; skip it
            continue

        qb = query_sequence[qp] if qp else "_"  # If none, turn into a blank
        rb = ref_sequence[rp - fm_rp].upper()

        if qb == rb and (rp - last_rp == 1 or last_rp == -1):
            if snv_group:
                rhs_contiguous += 1
            else:
                lhs_contiguous += 1

            if lhs_contiguous > contiguous_threshold and rhs_contiguous > contiguous_threshold:
                if len(snv_group) <= max_snv_group_size:
                    for snv in snv_group:
                        snvs[snv[:2]] = snv[-1]
                # Otherwise, it might be a little mismapped area or a longer deletion vs reference, so ignore it.
                lhs_contiguous = 0
                rhs_contiguous = 0
                snv_group.clear()

            last_rp = rp
            continue

        if rp - last_rp > 1:
            lhs_contiguous = 0
            last_rp = rp
            continue

        if qb != rb:
            snv_group.append((rp, rb, qb))
            # Don't reset either contiguous variable; instead, take this as part of a SNP group
            last_rp = rp

    return snvs
