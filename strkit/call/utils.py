import numpy as np
import sys
from typing import Callable

__all__ = [
    "cat_strs",
    "find_pair_by_ref_pos_py",
    "find_pair_by_ref_pos",
    "normalize_contig",
    "round_to_base_pos",
    "get_new_seed",
]


cat_strs = "".join


def find_pair_by_ref_pos_py(pairs: list[tuple[int, int]], target: int, start_left: int = 0) -> tuple[int, bool]:
    lhs: int = start_left
    rhs: int = len(pairs) - 1

    while lhs <= rhs:
        pivot: int = (lhs + rhs) // 2
        pair: tuple[int, int] = pairs[pivot]
        if pair[1] < target:
            lhs = pivot + 1
        elif pair[1] > target:  # pair[1] > snv_pos
            rhs = pivot - 1
        else:
            return pivot, True  # Found!

    # Nothing found, so must have been a gap
    # LHS should be the insertion point
    return lhs, False


find_pair_by_ref_pos: Callable[[list[tuple[int, int]], int], tuple[int, bool]]


if sys.version_info[0] >= 3 and (sys.version_info[0] > 3 or sys.version_info[1] >= 10):
    import bisect

    def find_pair_by_ref_pos_bisect(pairs: list[tuple[int, int]], target: int, start_left: int = 0) -> tuple[int, bool]:
        n_pairs: int = len(pairs)
        idx = bisect.bisect_left(pairs, target, start_left, n_pairs, key=lambda x: x[1])
        return idx, idx < n_pairs and pairs[idx][1] == target

    find_pair_by_ref_pos = find_pair_by_ref_pos_bisect

else:
    find_pair_by_ref_pos = find_pair_by_ref_pos_py


def normalize_contig(contig: str, has_chr: bool) -> str:
    return ("chr" if has_chr else "") + contig.replace("chr", "")


def round_to_base_pos(x, motif_size: int) -> float:
    return round(float(x) * motif_size) / motif_size


def get_new_seed(rng: np.random.Generator) -> int:
    return rng.integers(0, 4096, dtype=int)
