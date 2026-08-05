from __future__ import annotations

from functools import partial
from operator import eq, itemgetter, ne
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from numpy.random import Generator

__all__ = [
    "cn_getter",
    "sl_getter",
    "weight_getter",
    "eq_0",
    "neq_blank",
    "get_new_seed",
]


# index/property getters and other partials
cn_getter = itemgetter("cn")
sl_getter = itemgetter("sl")
weight_getter = itemgetter("w")
eq_0 = partial(eq, 0)
neq_blank = partial(ne, "")


def get_new_seed(rng: Generator) -> int:
    return rng.integers(0, 4096, dtype=int)
