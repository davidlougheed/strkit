import sys
from datetime import timedelta
from typing import Callable, Literal

from strkit import __version__
from strkit.json import Serializable, dumps, dumps_indented

from ..params import CallParams
from ..types import LocusResult

__all__ = [
    "output_json_report_header",
    "output_json_report_results",
    "output_json_report_footer",
]


def _get_dfn(indent_json: bool) -> Callable[[Serializable], bytes]:
    return dumps_indented if indent_json else dumps


def _write_bytes(b: bytes, json_path: str, mode: Literal["wb", "ab"]):
    if json_path == "stdout":
        sys.stdout.buffer.write(b)
        sys.stdout.flush()
    else:
        with open(json_path, mode) as jf:
            # noinspection PyTypeChecker
            jf.write(b)


def output_json_report_header(params: CallParams, contig_set: set[str], json_path: str, indent_json: bool):
    json_report_header = {
        "sample_id": params.sample_id,
        "caller": {
            "name": "strkit",
            "version": __version__,
        },
        "parameters": params.to_dict(as_inputted=True),
        "contigs": tuple(contig_set),
    }

    dfn = _get_dfn(indent_json)
    header_serialized: bytes = dfn(json_report_header)[:(-2 if indent_json else -1)]  # remove trailing ending brace

    # kludge: build up a portion of the JSON file, so we can output contig results as they come instead of storing them
    # in memory until the end of the run.
    header_serialized += b","
    if indent_json:
        header_serialized += b'\n  "results": [\n'
    else:
        header_serialized += b'"results":['

    # write partial JSON
    _write_bytes(header_serialized, json_path, "wb")


def output_json_report_results(results: tuple[LocusResult, ...], is_last: bool, json_path: str, indent_json: bool):
    dfn = _get_dfn(indent_json)
    results_bytes: bytes = dfn(results)

    if indent_json:
        results_bytes = results_bytes[2:-2]  # remove opening and closing "[]" + trailing newline
        if not is_last:
            results_bytes += b",\n"
    else:
        results_bytes = results_bytes[1:-1]  # remove opening and closing "[]"
        if not is_last:
            results_bytes += b","

    # write results "rows"
    _write_bytes(results_bytes, json_path, "ab")


def output_json_report_footer(time_taken: timedelta, json_path: str, indent_json: bool):
    runtime_bytes = dumps(time_taken.total_seconds())
    if indent_json:
        footer_bytes = b'\n  ],\n  "runtime": ' + runtime_bytes + b'\n}\n'
    else:
        footer_bytes = b'],"runtime":' + runtime_bytes + b'}\n'

    # write partial JSON
    _write_bytes(footer_bytes, json_path, "ab")
