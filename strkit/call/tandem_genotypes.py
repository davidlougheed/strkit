from typing import Optional, Tuple

from .allele import get_n_alleles, call_alleles

__all__ = [
    "call_tandem_genotypes",
]


def call_tandem_genotypes(args: Tuple[Optional[str], Optional[str], int, int, int, int, str]) -> str:
    contig: Optional[str] = args[0]
    sex_chr: Optional[str] = args[1]
    bootstrap_iterations: int = args[2]
    min_reads: int = args[3]
    min_allele_reads: int = args[4]
    read_bias_corr_min: int = args[5]
    line: str = args[6]

    gm_filter_factor: int = 3

    if line.startswith("#"):  # comment, strip it out
        return ""

    data: list = line.strip().split()

    locus_chr = data[0]

    if contig is not None and locus_chr != contig:
        return ""

    n_alleles: Optional[int] = get_n_alleles(2, sex_chr, locus_chr)
    if n_alleles is None:
        return ""  # No calling of sex chromosomes if we're not given the sex chromosome configuration

    call = call_alleles(
        tuple(map(int, data[-2].split(","))) if data[-2] != "." else (),
        tuple(map(int, data[-1].split(","))) if data[-1] != "." else (),
        None, None,
        bootstrap_iterations=bootstrap_iterations,
        min_reads=min_reads,
        min_allele_reads=min_allele_reads,
        n_alleles=n_alleles,
        separate_strands=True,
        read_bias_corr_min=read_bias_corr_min,
        gm_filter_factor=gm_filter_factor,
        force_int=True,
    )

    if call is None:
        # No call response
        return line.strip() + "\t" + "\t".join(["."] * (n_alleles * 3)) + "\n"  # 3: 1 exact + 2 CIs (95, 99)

    return (
        line.strip() + "\t" +
        "\t".join((
            *map(str, call["call"]),
            *(",".join(map(str, ci)) for ci in call["call_95_cis"]),
            *(",".join(map(str, ci)) for ci in call["call_99_cis"]),
        )) + "\n")
