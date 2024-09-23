import numpy as np
from numpy.typing import NDArray

from strkit_rust_ext import get_aligned_pair_matches

__all__ = [
    "decode_cigar_np",
    "get_aligned_pair_matches",
]


def decode_cigar_np(encoded_cigar: NDArray[np.uint32]) -> NDArray[np.uint32]:
    return np.stack((np.bitwise_and(encoded_cigar, 15), np.right_shift(encoded_cigar, 4)), axis=1)
