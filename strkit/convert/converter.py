from logging import Logger
from typing import Callable

from ._bed_4 import trf_to_bed_4
from .constants import IN_FORMAT_TRF, IN_FORMAT_TRGT, CONVERTER_IN_FORMATS
from .expansionhunter import trf_bed_to_eh
from .hipstr import trf_bed_to_hipstr
from .gangstr import trf_bed_to_gangstr
from .trgt import trgt_bed_to_bed4, trf_or_strkit_bed_to_trgt

import strkit.constants as c

__all__ = [
    "CONVERTER_OUTPUT_FORMATS",
    "convert",
]

convert_formats: dict[tuple[str, str], Callable[[list, Logger], None]] = {
    # TRF converters:
    (IN_FORMAT_TRF, c.CALLER_EXPANSIONHUNTER): trf_bed_to_eh,
    (IN_FORMAT_TRF, c.CALLER_HIPSTR): trf_bed_to_hipstr,
    (IN_FORMAT_TRF, c.CALLER_GANGSTR): trf_bed_to_gangstr,
    (IN_FORMAT_TRF, c.CALLER_REPEATHMM): lambda x: x,
    (IN_FORMAT_TRF, c.CALLER_STRAGLR): trf_to_bed_4,
    (IN_FORMAT_TRF, c.CALLER_STRKIT): trf_to_bed_4,  # or can just leave -asis
    (IN_FORMAT_TRF, c.CALLER_TANDEM_GENOTYPES): trf_to_bed_4,
    (IN_FORMAT_TRF, c.CALLER_TRGT): trf_or_strkit_bed_to_trgt,
    # TRGT converters:
    (IN_FORMAT_TRGT, c.CALLER_STRAGLR): trgt_bed_to_bed4,
    (IN_FORMAT_TRGT, c.CALLER_STRKIT): trgt_bed_to_bed4,
    (IN_FORMAT_TRGT, c.CALLER_TANDEM_GENOTYPES): trgt_bed_to_bed4,
}

CONVERTER_OUTPUT_FORMATS: tuple[str, ...] = tuple(sorted(set(k[1] for k in convert_formats)))


def convert(in_file: str, in_format: str, out_format: str, logger: Logger) -> int:
    out_format = out_format.lower()

    if in_format == IN_FORMAT_TRF:
        if out_format == c.CALLER_REPEATHMM:
            logger.critical(f"No need to convert for '{out_format}'; TRF BED files are accepted as input")
            return 1
        elif out_format == c.CALLER_STRKIT:
            logger.info("STRkit can use TRF BED files as-is; will convert to a BED4 file")

    if in_format not in CONVERTER_IN_FORMATS:
        logger.critical(f"Unsupported input format: {in_format}")

    if (in_format, out_format) not in convert_formats:
        logger.critical(f"Unsupported conversion: {in_format} -> {out_format} (no converter defined)")
        return 1

    with open(in_file, "r") as tf:
        data = [line.strip().split("\t") for line in tf]

    convert_formats[(in_format, out_format)](data, logger)
    return 0
