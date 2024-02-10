import orjson as json
from typing import Union

Serializable = Union[dict, list, tuple, str, int, float]


def dumps_indented(v: Serializable) -> bytes:
    return json.dumps(v, option=json.OPT_NON_STR_KEYS | json.OPT_INDENT_2)


__all__ = [
    "json",
    "dumps_indented",
]
