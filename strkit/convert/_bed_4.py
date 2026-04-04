from logging import Logger
from typing import Generator, Iterable

__all__ = [
    "trf_to_bed_4",
]


def trf_to_bed_4(trf_data: Iterable[list], _logger: Logger) -> Generator[str, None, None]:
    for item in trf_data:
        yield "\t".join((*item[:3], item[-1])) + "\n"
