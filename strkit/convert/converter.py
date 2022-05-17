import sys

from .expansionhunter import trf_bed_to_eh
from .hipstr import trf_bed_to_hipstr
from .gangstr import trf_bed_to_gangstr
from .straglr import trf_bed_to_straglr
from .tandem_genotypes import trf_bed_to_tandem_genotypes

import strkit.constants as c

__all__ = [
    "convert",
]

convert_formats = {
    c.CALLER_EXPANSIONHUNTER: trf_bed_to_eh,
    c.CALLER_HIPSTR: trf_bed_to_hipstr,
    c.CALLER_GANGSTR: trf_bed_to_gangstr,
    c.CALLER_REPEATHMM: lambda x: x,
    c.CALLER_STRAGLR: trf_bed_to_straglr,
    c.CALLER_TANDEM_GENOTYPES: trf_bed_to_tandem_genotypes,
}


def convert(trf_file: str, out_format: str) -> int:
    out_format = out_format.lower()

    if out_format == c.CALLER_REPEATHMM:
        sys.stderr.write(f"No need to convert for '{out_format}'; TRF BED files are accepted as input")
        return 1

    if out_format not in convert_formats:
        sys.stderr.write(f"Unsupported outout format: '{out_format}'")
        return 1

    with open(trf_file, "r") as tf:
        data = [line.strip().split("\t") for line in tf]

    convert_formats[out_format](data)
    return 0
