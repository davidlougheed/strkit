import sys
from datetime import timedelta

from strkit import __version__
from strkit.json import dumps_indented, json

from ..params import CallParams

__all__ = ["output_json_report"]


def output_json_report(
    params: CallParams,
    time_taken: timedelta,
    contig_set: set[str],
    results: list[dict],
    json_path: str,
    indent_json: bool,
):
    json_report = {
        "sample_id": params.sample_id,
        "caller": {
            "name": "strkit",
            "version": __version__,
        },
        "parameters": params.to_dict(as_inputted=True),
        "runtime": time_taken.total_seconds(),
        "contigs": tuple(contig_set),
        "results": results,
    }

    dfn = dumps_indented if indent_json else json.dumps
    report_data = dfn(json_report)

    if json_path == "stdout":
        sys.stdout.buffer.write(report_data)
        sys.stdout.write("\n")
        sys.stdout.flush()
    else:
        with open(json_path, "wb") as jf:
            jf.write(report_data)
            jf.write(b"\n")
