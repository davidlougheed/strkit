import sys
from logging import Logger

__all__ = [
    "trf_to_bed_4",
]


def trf_to_bed_4(trf_data: list, _logger: Logger):
    for item in trf_data:
        sys.stdout.write("\t".join((*item[:3], item[-1])) + "\n")
