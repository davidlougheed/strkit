from __future__ import annotations

from functools import partial
from operator import itemgetter, ne as operator_ne
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from numpy.random import Generator

__all__ = [
    "cn_getter",
    "sl_getter",
    "neq_blank",
    "get_new_seed",
]


# index/property getters and other partials
cn_getter = itemgetter("cn")
sl_getter = itemgetter("sl")
neq_blank = partial(operator_ne, "")


def get_new_seed(rng: Generator) -> int:
    return rng.integers(0, 4096, dtype=int)
