class SkipWritingLocus(Exception):
    pass


def blank_vcf_genotype(n_alleles: int) -> tuple[None, ...]:
    """
    Generates a blank VCF genotype for a sample locus.
    :param n_alleles: Number of alleles to use for the blank genotype.
    :return: Tuple of [None] * <number of alleles>
    """
    return tuple([None] * n_alleles)


def genotype_indices(alleles: tuple[str, ...], call: tuple[str, ...] | None, n_alleles: int) -> tuple[int | None, ...]:
    """
    Given possible variant alleles for a particular locus and a sample's call (the sample's alleles), plus the number of
    possible alleles for the sample at the locus, return a VCF-compatible genotype index tuple.
    :param alleles: Possible variant alleles: (ref, alt1, alt2, alt3, ...)
    :param call: Call for the sample: (alt3, alt1) or None if not called.
    :param n_alleles: Number of alleles for the sample locus (used for computing blank genotypes if not called.)
    :return: Genotype index tuple, e.g., (3, 1).
    """
    return tuple(map(alleles.index, call)) if call is not None else blank_vcf_genotype(n_alleles)
