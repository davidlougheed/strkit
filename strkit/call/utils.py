import numpy as np
import operator

from collections import deque
from functools import cache, partial
from numpy.typing import NDArray
from typing import Iterable

from ..utils import cat_strs

__all__ = [
    "cn_getter",
    "neq_blank",
    "normalize_contig",
    "get_new_seed",
    "calculate_seq_with_wildcards",
    "motif_rotations",
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


def motif_rotations(motif: str) -> Iterable[str]:
    yield motif

    if not motif:
        return

    already_sent: set[str] = {motif}

    motif_chars = deque(motif)
    for _ in range(len(motif) - 1):
        motif_chars.rotate(1)
        if (m := cat_strs(motif_chars)) not in already_sent:
            already_sent.add(m)
            yield m
