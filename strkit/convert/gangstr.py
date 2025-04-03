import sys
from logging import Logger

__all__ = [
    "trf_bed_to_gangstr",
]


def trf_bed_to_gangstr(trf_data: list, _logger: Logger):
    for i, item in enumerate(trf_data, 1):
        sys.stdout.write("\t".join((*item[:3], str(len(item[-1])), item[-1])) + "\n")
