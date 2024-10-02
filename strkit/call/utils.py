import numpy as np
import operator

from functools import cache, partial
from numpy.typing import NDArray
from typing import Optional

from ..utils import cat_strs

__all__ = [
    "idx_0_getter",
    "idx_1_getter",
    "cn_getter",
    "neq_blank",
    "find_pair_by_ref_pos",
    "normalize_contig",
    "round_to_base_pos",
    "get_new_seed",
    "calculate_seq_with_wildcards",
]


# index/property getters and other partials
idx_0_getter = operator.itemgetter(0)
idx_1_getter = operator.itemgetter(1)
cn_getter = operator.itemgetter("cn")
neq_blank = partial(operator.ne, "")


def find_pair_by_ref_pos(r_coords: NDArray[np.uint64], target: int, start_left: int = 0) -> tuple[int, bool]:
    n_pairs: int = len(r_coords)
    idx = start_left + np.searchsorted(r_coords[start_left:], target)
    return idx, idx < n_pairs and r_coords[idx] == target


def normalize_contig(contig: str, has_chr: bool) -> str:
    return ("chr" if has_chr else "") + contig.replace("chr", "")


def round_to_base_pos(x, motif_size: int) -> float:
    return round(float(x) * motif_size) / motif_size


def get_new_seed(rng: np.random.Generator) -> int:
    return rng.integers(0, 4096, dtype=int)


@cache  # TODO: parametrize base_wildcard_threshold
def _mask_low_q_base(base_and_qual: tuple[str, int], base_wildcard_threshold: int = 3) -> str:
    return base_and_qual[0] if base_and_qual[1] > base_wildcard_threshold else "X"


def calculate_seq_with_wildcards(qs: str, quals: Optional[NDArray[np.uint8]]) -> str:
    if quals is None:
        return qs  # No quality information, so don't do anything
    return cat_strs(map(_mask_low_q_base, zip(qs, quals)))
