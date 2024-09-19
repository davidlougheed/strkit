import numpy as np
from numpy.typing import NDArray
from typing import Iterable, Union

from strkit_rust_ext import get_aligned_pair_matches

__all__ = [
    "decode_cigar_np",
    "get_aligned_pair_matches",
]


def decode_cigar_np(encoded_cigar: NDArray[np.uint32]) -> Iterable[tuple[int, int]]:
    return zip(np.bitwise_and(encoded_cigar, 15, dtype=int), np.right_shift(encoded_cigar, 4, dtype=int))
