import sys

__all__ = [
    "trf_bed_to_hipstr",
]


def trf_bed_to_hipstr(trf_data: list):
    for i, item in enumerate(trf_data, 1):
        sys.stdout.write("\t".join((*item[:3], str(len(item[-1])), str(round(float(item[5]))))) + "\n")
