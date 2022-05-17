import pysam
import sys

from typing import List, Literal, Optional, Tuple
from strkit.call.allele import call_alleles

__all__ = [
    "aligned_lengths",
    "aligned_lengths_cmd",
]


def aligned_lengths(file: str, region: str) -> Tuple[Optional[List[int]], Optional[List[int]]]:
    try:
        chrom, coords = region.split(":")
        pos1, pos2 = map(int, coords.split("-"))
    except ValueError:
        sys.stderr.write("Error: please provide region in chr#:pos1-pos2 format\n")
        return None, None

    mode: Literal["r", "rb", "rc"]
    if file.endswith(".sam"):
        mode = "r"
    elif file.endswith(".bam"):
        mode = "rb"
    elif file.endswith(".cram"):
        mode = "rc"
    else:
        sys.stderr.write("Error: please provide a .(b|cr|s)am file\n")
        return None, None

    align = pysam.AlignmentFile(file, mode)
    fwd_region_lengths = []
    rev_region_lengths = []

    for read in align.fetch(chrom, pos1, pos2):
        region_length = 0
        for read_pos, ref_pos in read.get_aligned_pairs():
            if ref_pos is not None and ref_pos < pos1:
                continue
            if ref_pos is not None and ref_pos > pos2:
                break
            region_length += 1

        # print(read.get_aligned_pairs())
        if read.is_reverse:
            rev_region_lengths.append(region_length)
        else:
            fwd_region_lengths.append(region_length)

    return fwd_region_lengths, rev_region_lengths


def aligned_lengths_cmd(file: str, region: str) -> int:
    r1, r2 = aligned_lengths(file, region)
    call = call_alleles(
        repeats_fwd=r1,
        repeats_rev=r2,
        read_weights_fwd=None,
        read_weights_rev=None,
        bootstrap_iterations=100,
        min_reads=2,
        min_allele_reads=1,
        n_alleles=2,
        separate_strands=False,
        read_bias_corr_min=4,
        gm_filter_factor=3,
        force_int=True,
    )
    fwd_str = ',\t'.join(map(str, r1))
    rev_str = ',\t'.join(map(str, r2))
    print(f"Aligned lengths:")
    print(f" Forward strand: {fwd_str}")
    print(f" Reverse strand: {rev_str}")
    print(f"Best guess for allele lengths: {call['call'][0]}, {call['call'][1]}")
    print(f"                      95% CIs: {call['call_95_cis'][0]}, {call['call_95_cis'][1]}")
    print(f"                      99% CIs: {call['call_99_cis'][0]}, {call['call_99_cis'][1]}")
    return 1 if r1 is None else 0
