from typing import Optional, Tuple

from .allele import call_alleles
from ..constants import SEX_CHROMOSOMES

__all__ = [
    "preprocess_lines_straglr",
    "call_straglr",
]


def preprocess_lines_straglr(lines: list) -> list:
    # Need to group lines by locus since straglr does read-level stuff
    new_lines = []
    acc = None

    for line in lines:
        if line.startswith("#"):
            continue

        data = line.strip().split("\t")

        if acc is not None and tuple(data[:4]) != acc[0]:
            new_lines.append(acc)
            acc = None

        if acc is None:
            # (chrom start end repeat_unit) | genotypes read_ids copy_numbers sizes read_starts | allele
            acc = (tuple(data[:4]), data[4], [], [], [], [], data[9])

        acc[2].append(data[5])
        acc[3].append(data[6])
        acc[4].append(data[7])
        acc[5].append(data[8])

    new_lines.append(acc)

    return new_lines


def _tenths_str(s: float):
    return f"{round(s * 10) / 10:.1f}"


def call_straglr(args: Tuple[Optional[str], str, int, int, int, int, tuple]) -> str:
    contig: Optional[str] = args[0]
    sex_chr: str = args[1]
    bootstrap_iterations: int = args[2]
    min_reads: int = args[3]
    min_allele_reads: int = args[4]
    # read_bias_corr_min isn't used here (args[5])
    data: tuple = args[6]

    n_alleles: int = 2
    gm_filter_factor: int = 3

    line: str = "\t".join((*data[0], data[1])) + "\t" + "/".join(data[3])
    null_response: str = line + "\t" + "\t".join(["."] * (n_alleles * 2)) + "\n"

    locus_chr = data[0][0]

    if contig is not None and locus_chr != contig:
        return ""

    if sex_chr == "NONE" and locus_chr in SEX_CHROMOSOMES:
        return ""  # No calling of sex chromosomes if we're not given the sex chromosomes

    if locus_chr in SEX_CHROMOSOMES and sex_chr == "XY":
        n_alleles = 1

    allele_estimates, allele_cis_95, allele_cis_99 = call_alleles(
        tuple(map(float, data[3])),
        (),
        bootstrap_iterations=bootstrap_iterations,
        min_reads=min_reads,
        min_allele_reads=min_allele_reads,
        n_alleles=n_alleles,
        separate_strands=False,
        read_bias_corr_min=0,  # Not relevant
        gm_filter_factor=gm_filter_factor,
        force_int=False,
    )

    if allele_estimates is None:
        return null_response

    return (
        line + "\t" +
        "\t".join((
            *map(_tenths_str, allele_estimates),
            *(",".join(map(_tenths_str, ci)) for ci in allele_cis_95),
            *(",".join(map(_tenths_str, ci)) for ci in allele_cis_99),
        )) + "\n")
