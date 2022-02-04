import pysam
import sys

__all__ = [
    "aligned_lengths",
]


def aligned_lengths(file: str, region: str) -> int:
    try:
        chrom, coords = region.split(":")
        pos1, pos2 = map(int, coords.split("-"))
    except ValueError:
        sys.stderr.write("Error: please provide region in chr#:pos1-pos2 format\n")
        return 1

    if file.endswith(".sam"):
        mode = "r"
    elif file.endswith(".bam"):
        mode = "rb"
    elif file.endswith(".cram"):
        mode = "rc"
    else:
        sys.stderr.write("Error: please provide a .(b|cr|s)am file\n")
        return 1

    align = pysam.AlignmentFile(file, mode)
    region_lengths = []

    for read in align.fetch(chrom, pos1, pos2):
        region_lengths.append(read.query_alignment_length)

    print(f"Aligned lengths: {', '.join(map(str, region_lengths))}")
