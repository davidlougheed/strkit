from collections.abc import Generator, Iterable
from logging import Logger

__all__ = [
    "trf_bed_to_gangstr",
]


def trf_bed_to_gangstr(trf_data: Iterable[list], _logger: Logger) -> Generator[str, None, None]:
    for i, item in enumerate(trf_data, 1):
        yield "\t".join((*item[:3], str(len(item[-1])), item[-1])) + "\n"
