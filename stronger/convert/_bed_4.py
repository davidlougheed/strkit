import sys

__all__ = [
    "to_bed_4",
]


def to_bed_4(trf_data: list):
    for i, item in enumerate(trf_data, 1):
        sys.stdout.write("\t".join((*item[:3], item[-1])) + "\n")
