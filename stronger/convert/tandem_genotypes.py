__all__ = [
    "trf_bed_to_tandem_genotypes",
]

from ._bed_4 import to_bed_4


def trf_bed_to_tandem_genotypes(trf_data: list):
    to_bed_4(trf_data)
