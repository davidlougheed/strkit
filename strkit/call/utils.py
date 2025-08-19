import numpy as np
import operator

from functools import partial

__all__ = [
    "cn_getter",
    "neq_blank",
    "normalize_contig",
    "get_new_seed",
]


# index/property getters and other partials
cn_getter = operator.itemgetter("cn")
neq_blank = partial(operator.ne, "")


def normalize_contig(contig: str, has_chr: bool) -> str:
    return ("chr" if has_chr else "") + contig.replace("chr", "")


def get_new_seed(rng: np.random.Generator) -> int:
    return rng.integers(0, 4096, dtype=int)
