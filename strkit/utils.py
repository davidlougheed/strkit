from __future__ import annotations

import math
from typing import Any, Callable, Iterable, Union

__all__ = [
    "apply_or_none",
    "int_tuple",
    "float_tuple",
    "parse_cis",
    "cis_overlap",
    "sign",
]


def apply_or_none(fn: Callable, x: Any) -> Any:
    # Python: add any type of monad functionality challenge [IMPOSSIBLE]
    return fn(x) if x is not None else None


def int_tuple(x: Iterable) -> tuple[int, ...]:
    return tuple(map(int, x))


def float_tuple(x: Iterable) -> tuple[float, ...]:
    return tuple(map(float, x))


def parse_cis(cis, commas=False, dtype=int):
    return tuple(map(lambda ci: tuple(map(dtype, ci.split("," if commas else "-"))), cis))


def cis_overlap(ci1, ci2) -> bool:
    epsilon = -0.0001

    # []: ci1
    # (): ci2
    # [   (    ]   )  or  [   (    )   ]  or  (   [    )   ]  or  (   [    ]   )
    # int logic: ci1[0] <= ci2[1] and ci2[0] <= ci1[1]
    # float logic: lets add some epsilon to prevent little issues
    return (ci2[1] - ci1[0]) > epsilon and (ci1[1] - ci2[0]) > epsilon


def sign(x: Union[int, float]) -> int:
    return round(math.copysign(1, x))
