import numpy as np
from numpy.typing import NDArray
from typing import Iterable, Union

from strkit_rust_ext import get_aligned_pair_matches

__all__ = [
    "CoordPair",
    "decode_cigar",
    "decode_cigar_np",
    "get_aligned_pair_matches",
]


CoordPair = tuple[Union[int, None], Union[int, None]]


def _decode_cigar_item(item: int) -> tuple[int, int]:
    return item & 15, item >> 4


def decode_cigar(encoded_cigar: list[int]) -> Iterable[tuple[int, int]]:
    return map(_decode_cigar_item, encoded_cigar)


def decode_cigar_np(encoded_cigar: NDArray[np.uint32]) -> Iterable[tuple[int, int]]:
    return zip(np.bitwise_and(encoded_cigar, 15, dtype=int), np.right_shift(encoded_cigar, 4, dtype=int))
