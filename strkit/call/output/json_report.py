import sys
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


def _indent_lines(json_bytes: bytes, indent: int) -> bytes:
    """
    Given serialized, indented JSON bytes, indent them by extra spaces.
    :param json_bytes: Serialized, indented JSON bytes.
    :param indent: Number of spaces to indent the JSON by.
    :return: Extra-indented JSON bytes.
    """
    return b"\n".join(
        b" " * indent + line if line and i > 0 else line  # skip first indent since this will be attached to a key
        for i, line in enumerate(json_bytes.split(b"\n"))
    )


def output_json_report_footer(
    avg_read_depth: float,
    avg_read_depths_by_contig: dict[str, float],
    time_taken: float,
    json_path: str,
    indent_json: bool,
):
    avg_read_depth_bytes = dumps(avg_read_depth)  # ensure JSON-formatted float (probably needless)
    avg_read_depths_by_contig_bytes = (dumps_indented if indent_json else dumps)(avg_read_depths_by_contig)
    runtime_bytes = dumps(time_taken)
    if indent_json:
        footer_bytes = (
            b'\n  ],\n  "avg_read_depth": ' + avg_read_depth_bytes + b',\n  "avg_read_depths_by_contig": '
            + _indent_lines(avg_read_depths_by_contig_bytes, indent=2)
            + b',\n  "runtime": '
            + runtime_bytes + b'\n}\n'
        )
    else:
        footer_bytes = (
            b'],"avg_read_depth":' + avg_read_depth_bytes + b',"avg_read_depths_by_contig":'
            + avg_read_depths_by_contig_bytes
            + b',"runtime":'
            + runtime_bytes
            + b'}\n'
        )

    # write partial JSON
    _write_bytes(footer_bytes, json_path, "ab")
