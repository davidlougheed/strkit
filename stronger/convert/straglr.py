__all__ = [
    "trf_bed_to_straglr",
]

from ._bed_4 import to_bed_4


def trf_bed_to_straglr(trf_data: list):
    to_bed_4(trf_data)
