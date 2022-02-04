import ast
import re

from typing import Optional, Tuple

from .allele import call_alleles
from ..constants import SEX_CHROMOSOMES

__all__ = [
    "call_repeathmm",
]

WSP_SPLIT = re.compile(r"\s{3,7}")


def call_repeathmm(args: Tuple[Optional[str], str, int, int, int, int, str]) -> str:
    contig: Optional[str] = args[0]
    sex_chr = args[1]
    bootstrap_iterations: int = args[2]
    min_reads: int = args[3]
    min_allele_reads: int = args[4]
    # skip read_bias_corr_min -- not relevant to repeatHMM
    line: str = args[6]

    n_alleles: int = 2
    gm_filter_factor: int = 3

    data: list = WSP_SPLIT.split(line.strip())

    no_call_response = "\t".join((
        data[0],
        *(["."] * (n_alleles * 2)),
        data[1],
    )) + "\n"

    if "allocr:" not in data[1]:
        # No call for this locus
        return no_call_response

    locus_chr = data[0].split(":")[0]

    if contig is not None and locus_chr != contig:
        return ""  # No output if we're not processing this contig

    if sex_chr == "NONE" and locus_chr in SEX_CHROMOSOMES:
        return ""  # No calling of sex chromosomes if we're not given the sex chromosomes

    if locus_chr in SEX_CHROMOSOMES and sex_chr == "XY":
        n_alleles = 1

    # Cut out 'p2sp'/'p2bamhmm' and ending '><'
    call_and_read_data = ast.literal_eval(
        data[1].lstrip("p2sp ").lstrip("p2bamhmm ").rstrip("><"))

    if call_and_read_data[4] < min_reads:
        return no_call_response  # Not enough supporting reads; don't call anything

    read_counts = [
        tuple(map(int, v.lstrip("allocr:").split(":")))
        for v in (w.strip() for w in call_and_read_data[3].split(","))
        if v
    ]

    read_counts_expanded = []
    for rv, rc in read_counts:
        if rv == 0:
            # Noisy read from RepeatHMM, skip it.
            continue
        read_counts_expanded.extend([rv] * rc)

    allele_estimates, allele_cis_95, allele_cis_99 = call_alleles(
        read_counts_expanded,
        (),  # RepeatHMM doesn't separate read counts by strand
        bootstrap_iterations=bootstrap_iterations,
        min_reads=min_reads,
        min_allele_reads=min_allele_reads,
        n_alleles=n_alleles,
        separate_strands=False,
        read_bias_corr_min=0,  # Not relevant since we don't have separate strand data
        gm_filter_factor=gm_filter_factor,
        force_int=True,
    )

    if allele_estimates is None:
        return no_call_response

    return (
        "\t".join((
            data[0],
            *map(str, allele_estimates),
            *(",".join(map(str, ci)) for ci in allele_cis_95),
            *(",".join(map(str, ci)) for ci in allele_cis_99),
            data[1],
        )) + "\n")
