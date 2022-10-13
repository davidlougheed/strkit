from __future__ import annotations

import ast
import re

from typing import Optional

from .allele import get_n_alleles, call_alleles

__all__ = [
    "call_repeathmm",
]

WSP_SPLIT = re.compile(r"\s{3,7}")


def call_repeathmm(args: tuple[Optional[str], Optional[str], int, int, int, int, int, str]) -> str:
    contig: Optional[str] = args[0]
    sex_chr: Optional[str] = args[1]
    bootstrap_iterations: int = args[2]
    min_reads: int = args[3]
    min_allele_reads: int = args[4]
    # skip read_bias_corr_min -- not relevant to repeatHMM
    seed: int = args[6]
    line: str = args[7]

    gm_filter_factor: int = 3

    data: list = WSP_SPLIT.split(line.strip())

    locus_chr = data[0].split(":")[0]

    if contig is not None and locus_chr != contig:
        return ""  # No output if we're not processing this contig

    n_alleles: Optional[int] = get_n_alleles(2, sex_chr, locus_chr)
    if n_alleles is None:
        return ""  # No calling of sex chromosomes if we're not given the sex chromosome configuration

    no_call_response = "\t".join((
        data[0],
        *(["."] * (n_alleles * 3)),  # 3: 1 exact + 2 CIs (95, 99)
        data[1],
    )) + "\n"

    if "allocr:" not in data[1]:
        # No call for this locus
        return no_call_response

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

    call = call_alleles(
        read_counts_expanded, (),  # RepeatHMM doesn't separate read counts by strand
        None, None,
        bootstrap_iterations=bootstrap_iterations,
        min_reads=min_reads,
        min_allele_reads=min_allele_reads,
        n_alleles=n_alleles,
        separate_strands=False,
        read_bias_corr_min=0,  # Not relevant since we don't have separate strand data
        gm_filter_factor=gm_filter_factor,
        hq=False,
        force_int=True,
        seed=seed,
    )

    if call is None:
        return no_call_response

    return (
        "\t".join((
            data[0],
            *map(str, call["call"]),
            *(",".join(map(str, ci)) for ci in call["call_95_cis"]),
            *(",".join(map(str, ci)) for ci in call["call_99_cis"]),
            data[1],
        )) + "\n")
