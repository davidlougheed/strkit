import numpy as np
import operator

from functools import cache, partial
from numpy.typing import NDArray

from ..utils import cat_strs

__all__ = [
    "cn_getter",
    "neq_blank",
    "normalize_contig",
    "get_new_seed",
    "calculate_seq_with_wildcards",
]


# index/property getters and other partials
cn_getter = operator.itemgetter("cn")
neq_blank = partial(operator.ne, "")


def normalize_contig(contig: str, has_chr: bool) -> str:
    return ("chr" if has_chr else "") + contig.replace("chr", "")


def get_new_seed(rng: np.random.Generator) -> int:
    return rng.integers(0, 4096, dtype=int)


@cache  # TODO: parametrize base_wildcard_threshold
def _mask_low_q_base(base_and_qual: tuple[str, int], base_wildcard_threshold: int = 3) -> str:
    return base_and_qual[0] if base_and_qual[1] > base_wildcard_threshold else "X"


def calculate_seq_with_wildcards(qs: str, quals: NDArray[np.uint8] | None) -> str:
    if quals is None:
        return qs  # No quality information, so don't do anything
    return cat_strs(map(_mask_low_q_base, zip(qs, quals)))
