from typing import Iterable, Union

from strkit_rust_ext import get_aligned_pair_matches

__all__ = [
    "CoordPair",
    "decode_cigar",
    "get_aligned_pair_matches",
]


CoordPair = tuple[Union[int, None], Union[int, None]]


def _decode_cigar_item(item: int) -> tuple[int, int]:
    return item & 15, item >> 4


def decode_cigar(encoded_cigar: list[int]) -> Iterable[tuple[int, int]]:
    return map(_decode_cigar_item, encoded_cigar)
