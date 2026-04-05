import sys

from logging import Logger
from typing import Callable, Generator, Iterable, TextIO, TypeAlias

from ._bed_4 import trf_to_bed_4
from .constants import FORMAT_BED4, FORMAT_TRF, FORMAT_TRGT, CONVERTER_IN_FORMATS
from .expansionhunter import trf_bed_to_eh
from .hipstr import trf_bed_to_hipstr
from .gangstr import trf_bed_to_gangstr
from .trf import trf_dat_to_bed, trf_passthrough
from .trgt import trgt_bed_to_bed4, trgt_bed_to_strkit_bed, trf_or_strkit_bed_to_trgt

import strkit.constants as c

__all__ = [
    "CONVERTER_OUTPUT_FORMATS",
    "convert",
]

ConverterFn: TypeAlias = Callable[[Iterable, Logger], Generator[str, None, None]]

convert_formats: dict[tuple[str | tuple[str, ...], str], ConverterFn] = {
    # TRF DAT to BED converter:
    (FORMAT_TRF, FORMAT_TRF): trf_passthrough,
    # BED4/TRF BED converters:
    ((FORMAT_BED4, FORMAT_TRF), c.CALLER_EXPANSIONHUNTER): trf_bed_to_eh,
    ((FORMAT_BED4, FORMAT_TRF), c.CALLER_HIPSTR): trf_bed_to_hipstr,
    ((FORMAT_BED4, FORMAT_TRF), c.CALLER_GANGSTR): trf_bed_to_gangstr,
    ((FORMAT_BED4, FORMAT_TRF), c.CALLER_REPEATHMM): lambda x: x,
    ((FORMAT_BED4, FORMAT_TRF), FORMAT_BED4): trf_to_bed_4,
    ((FORMAT_BED4, FORMAT_TRF), c.CALLER_STRAGLR): trf_to_bed_4,
    ((FORMAT_BED4, FORMAT_TRF), c.CALLER_STRKIT): trf_to_bed_4,  # or can just leave as-is
    ((FORMAT_BED4, FORMAT_TRF), c.CALLER_TANDEM_GENOTYPES): trf_to_bed_4,
    ((FORMAT_BED4, FORMAT_TRF), c.CALLER_TRGT): trf_or_strkit_bed_to_trgt,
    # TRGT converters:
    (FORMAT_TRGT, c.CALLER_STRAGLR): trgt_bed_to_bed4,
    (FORMAT_TRGT, c.CALLER_STRKIT): trgt_bed_to_strkit_bed,
    (FORMAT_TRGT, c.CALLER_TANDEM_GENOTYPES): trgt_bed_to_bed4,
}

CONVERTER_OUTPUT_FORMATS: tuple[str, ...] = tuple(sorted(set(k[1] for k in convert_formats)))


def _load_bed_like(fh: TextIO, sort: bool) -> list:
    data = []
    contig_order = []

    for line in map(lambda s: s.split("\t"), map(str.strip, fh)):
        if not contig_order or line[0] != contig_order[-1]:
            contig_order.append(line[0])
        data.append([line[0], int(line[1]), int(line[2]), *line[3:]] if sort else line)

    if sort:
        data.sort(key=lambda dd: (contig_order.index(dd[0]), dd[1], dd[2]))
        for d in data:
            d[1] = str(d[1])
            d[2] = str(d[2])

    return data


def convert(in_file: str, in_format: str, out_format: str, sort: bool, logger: Logger) -> int:
    out_format = out_format.lower()

    is_trf_dat: bool = False

    if in_format == FORMAT_TRF:
        if in_file.endswith(".dat"):
            is_trf_dat = True
        elif out_format == FORMAT_TRF:
            logger.critical(f"No need to convert from TRF BED to TRF BED")
            exit(1)
        elif out_format == c.CALLER_REPEATHMM:
            logger.critical(f"No need to convert for '{out_format}'; TRF BED files are accepted as input")
            return 1
        elif out_format == c.CALLER_STRKIT:
            logger.info("STRkit can use TRF BED files as-is; will convert to a BED4 file")

    if in_format not in CONVERTER_IN_FORMATS:
        logger.critical(f"Unsupported input format: {in_format}")

    converter: ConverterFn | None = None
    for (i, o), cc in convert_formats.items():
        if in_format in i and out_format == o:
            converter = cc

    if converter is None:
        logger.critical(f"Unsupported conversion: {in_format} -> {out_format} (no converter defined)")
        return 1

    with open(in_file, "r") as tf:
        if is_trf_dat:
            data = trf_dat_to_bed(tf, sort)
        else:
            # in all other cases, we have various forms of BED files
            data = _load_bed_like(tf, sort)
        # noinspection PyCallingNonCallable
        for line in converter(data, logger):
            sys.stdout.write(line)

    return 0
