from typing import Iterable, Tuple

__all__ = [
    "int_tuple",
    "parse_cis",
    "cis_overlap",
]


def int_tuple(x: Iterable) -> Tuple[int, ...]:
    return tuple(map(int, x))


def parse_cis(cis, commas=False):
    return tuple(map(lambda ci: tuple(map(int, ci.split("," if commas else "-"))), cis))


def cis_overlap(ci1, ci2) -> bool:
    epsilon = -0.0001

    # []: ci1
    # (): ci2
    # [   (    ]   )  or  [   (    )   ]  or  (   [    )   ]  or  (   [    ]   )
    # int logic: ci1[0] <= ci2[1] and ci2[0] <= ci1[1]
    # float logic: lets add some epsilon to prevent little issues
    return (ci2[1] - ci1[0]) > epsilon and (ci1[1] - ci2[0]) > epsilon
