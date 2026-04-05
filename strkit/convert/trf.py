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


def trf_dat_to_bed(file: TextIO, sort: bool) -> Generator[list, None, None]:
    res = []  # only used if sort=True

    def _output_sorted_res():
        res.sort(key=lambda dd: (dd[1], dd[2]))
        for r in res:
            r[1] = str(r[1])
            r[2] = str(r[2])
            yield r
        res.clear()

    current_contig = ""
    for line in map(str.strip, file):
        if not line:
            continue

        if line.startswith(TRF_SEQUENCE_PREFIX):
            current_contig = line.removeprefix(TRF_SEQUENCE_PREFIX)
            if sort:
                yield from _output_sorted_res()

        data = line.split(" ")
        if data[0].isdigit() and len(data) == 15:  # data column
            # TRF start coord is 1-indexed
            start = int(data[0]) - 1
            end = int(data[1])
            if sort:
                res.append([current_contig, start, end, "trf", *data[2:14]])
            else:
                yield [current_contig, str(start), str(end), "trf", *data[2:14]]

    yield from _output_sorted_res()


def trf_passthrough(trf_data: Iterable[list], _logger: Logger) -> Generator[str, None, None]:
    for data in trf_data:
        yield "\t".join(data) + "\n"
