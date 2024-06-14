import orjson as json
from typing import Union


__all__ = [
    "Serializable",
    "json",
    "dumps",
    "dumps_indented",
]


Serializable = Union[dict, list, tuple, str, int, float]


def dumps(v: Serializable) -> bytes:
    return json.dumps(v, option=json.OPT_NON_STR_KEYS | json.OPT_SERIALIZE_NUMPY)


def dumps_indented(v: Serializable) -> bytes:
    return json.dumps(v, option=json.OPT_NON_STR_KEYS | json.OPT_INDENT_2 | json.OPT_SERIALIZE_NUMPY)
