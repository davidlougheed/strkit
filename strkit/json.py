from typing import Union

Serializable = Union[dict, list, tuple, str, int, float]

try:
    import orjson as json

    def dumps_indented(v: Serializable) -> bytes:
        return json.dumps(v, option=json.OPT_INDENT_2)

except ImportError:
    import json

    def dumps_indented(v: Serializable) -> bytes:
        return json.dumps(v, indent=2).encode("utf-8")

__all__ = [
    "json",
    "dumps_indented",
]
