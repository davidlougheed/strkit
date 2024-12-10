from __future__ import annotations

import math
import operator
from functools import partial
from typing import Any, Callable, Iterable

__all__ = [
    "cat_strs",
    "is_none",
    "idx_0_getter",
    "idx_1_getter",
    "apply_or_none",
    "int_tuple",
    "float_tuple",
    "parse_ci",
    "parse_cis",
    "cis_overlap",
    "sign",
]


# index/property getters and other partials
cat_strs = "".join
is_none = partial(operator.is_, None)
idx_0_getter = operator.itemgetter(0)
idx_1_getter = operator.itemgetter(1)


def apply_or_none(fn: Callable, x: Any) -> Any:
    # Python: add any type of monad functionality challenge [IMPOSSIBLE]
    return fn(x) if x is not None else None


def int_tuple(x: Iterable) -> tuple[int, ...]:
    return tuple(map(int, x))


def float_tuple(x: Iterable) -> tuple[float, ...]:
    return tuple(map(float, x))


def parse_ci(ci: str, commas=False, dtype=int) -> tuple[int, int] | tuple[float, float]:
    ci_s = ci.split("," if commas else "-")
    return dtype(ci_s[0]), dtype(ci_s[1])


def parse_cis(
    cis: Iterable[str], commas=False, dtype=int
) -> tuple[tuple[int, ...], ...] | tuple[tuple[float, ...], ...]:
    return tuple(map(lambda ci: parse_ci(ci, commas, dtype), cis))


def cis_overlap(ci1, ci2) -> bool:
    epsilon = -0.0001

    # []: ci1
    # (): ci2
    # [   (    ]   )  or  [   (    )   ]  or  (   [    )   ]  or  (   [    ]   )
    # int logic: ci1[0] <= ci2[1] and ci2[0] <= ci1[1]
    # float logic: lets add some epsilon to prevent little issues
    return (ci2[1] - ci1[0]) > epsilon and (ci1[1] - ci2[0]) > epsilon


def sign(x: int | float) -> int:
    return round(math.copysign(1, x))
