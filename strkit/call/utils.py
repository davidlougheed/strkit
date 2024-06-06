import bisect
import numpy as np
import operator

from functools import partial
from numpy.typing import NDArray

__all__ = [
    "idx_0_getter",
    "idx_1_getter",
    "neq_blank",
    "find_pair_by_ref_pos",
    "normalize_contig",
    "round_to_base_pos",
    "get_new_seed",
]


idx_0_getter = operator.itemgetter(0)
idx_1_getter = operator.itemgetter(1)
neq_blank = partial(operator.ne, "")


def find_pair_by_ref_pos(r_coords: NDArray[np.uint64], target: int, start_left: int = 0) -> tuple[int, bool]:
    n_pairs: int = len(r_coords)
    idx = bisect.bisect_left(r_coords, target, start_left, n_pairs)
    return idx, idx < n_pairs and r_coords[idx] == target


def normalize_contig(contig: str, has_chr: bool) -> str:
    return ("chr" if has_chr else "") + contig.replace("chr", "")


def round_to_base_pos(x, motif_size: int) -> float:
    return round(float(x) * motif_size) / motif_size


def get_new_seed(rng: np.random.Generator) -> int:
    return rng.integers(0, 4096, dtype=int)
