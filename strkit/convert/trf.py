import sys

from logging import Logger
from typing import Generator, Iterable, TextIO

__all__ = ["trf_dat_to_bed", "trf_passthrough"]

TRF_SEQUENCE_PREFIX = "Sequence: "

# Input TRF dat file columns: see https://tandem.bu.edu/trf/definitions#table
# start
# end
# period
# copy number
# consensus pattern size
# % matches
# % indels
# align score
# % A
# % C
# % G
# % T
# entropy
# motif
# sequence

# Output BED columns: same as above, but with contig, with "trf" as the 4th column, and without sequence


def trf_dat_to_bed(file: TextIO) -> Generator[list, None, None]:
    current_contig = ""
    for line in map(str.strip, file):
        if not line:
            continue
        if line.startswith(TRF_SEQUENCE_PREFIX):
            current_contig = line.removeprefix(TRF_SEQUENCE_PREFIX)
        data = line.split(" ")
        if data[0].isdigit() and len(data) == 15:  # data column
            yield [current_contig, *data[:2], "trf", *data[2:]]


def trf_passthrough(trf_data: Iterable[list], _logger: Logger):
    for data in trf_data:
        sys.stdout.write("\t".join(data) + "\n")
