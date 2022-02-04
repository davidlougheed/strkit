from typing import Optional, Tuple

from .allele import call_allele
from ..constants import SEX_CHROMOSOMES

__all__ = [
    "call_tandem_genotypes",
]


def call_tandem_genotypes(args: Tuple[Optional[str], str, int, int, int, int, str]) -> str:
    contig: Optional[str] = args[0]
    sex_chr: str = args[1]
    bootstrap_iterations: int = args[2]
    min_reads: int = args[3]
    min_allele_reads: int = args[4]
    read_bias_corr_min: int = args[5]
    line: str = args[6]

    n_alleles: int = 2
    gm_filter_factor: int = 3

    if line.startswith("#"):  # comment, strip it out
        return ""

    data: list = line.strip().split()
    null_response = line.strip() + "\t" + "\t".join(["."] * (n_alleles * 2)) + "\n"

    locus_chr = data[0]

    if contig is not None and locus_chr != contig:
        return ""

    if sex_chr == "NONE" and locus_chr in SEX_CHROMOSOMES:
        return ""  # No calling of sex chromosomes if we're not given the sex chromosomes

    if locus_chr in SEX_CHROMOSOMES and sex_chr == "XY":
        n_alleles = 1

    allele_estimates, allele_cis_95, allele_cis_99 = call_allele(
        tuple(map(int, data[6].split(","))) if data[6] != "." else (),
        tuple(map(int, data[7].split(","))) if data[7] != "." else (),
        bootstrap_iterations=bootstrap_iterations,
        min_reads=min_reads,
        min_allele_reads=min_allele_reads,
        n_alleles=n_alleles,
        separate_strands=True,
        read_bias_corr_min=read_bias_corr_min,
        gm_filter_factor=gm_filter_factor,
        force_int=True,
    )

    if allele_estimates is None:
        return null_response

    return (
        line.strip() + "\t" +
        "\t".join((
            *map(str, allele_estimates),
            *(",".join(map(str, ci)) for ci in allele_cis_95),
            *(",".join(map(str, ci)) for ci in allele_cis_99),
        )) + "\n")
