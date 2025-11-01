import numpy as np
import orjson as json


__all__ = [
    "Serializable",
    "json",
    "dumps",
    "dumps_indented",
]


Serializable = dict | list | tuple | str | int | float


def _dumps_default(x):
    # Don't know why this happens
    if isinstance(x, np.ndarray) and not x.flags.c_contiguous:
        return x.tolist()
    raise TypeError


def dumps(v: Serializable) -> bytes:
    return json.dumps(v, option=json.OPT_NON_STR_KEYS | json.OPT_SERIALIZE_NUMPY, default=_dumps_default)


def dumps_indented(v: Serializable) -> bytes:
    return json.dumps(
        v, option=json.OPT_NON_STR_KEYS | json.OPT_INDENT_2 | json.OPT_SERIALIZE_NUMPY, default=_dumps_default
    )
