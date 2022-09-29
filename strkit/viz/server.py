from flask import Flask, render_template, request, send_file
from werkzeug.exceptions import NotFound

__all__ = [
    "run_server",
]

app = Flask(__name__)


@app.route("/")
def browser():
    return render_template(
        "browser.html",
        **app.config["PARAMS"])


@app.route("/report-metadata")
def get_report_metadata():
    return {k: v for k, v in app.config["CALL_REPORT"].items() if k != "results"}


@app.route("/params")
def get_params():
    return {
        "cmd": app.config["PARAMS"],
        "report": app.config["CALL_REPORT"]["parameters"],
    }


@app.route("/loci")
def get_loci():
    cr = app.config["CALL_REPORT"]
    ecd = list(enumerate(cr["results"]))  # TODO: cache

    q = request.args.get("q", "").strip()
    if q:
        res = list(filter(lambda x: q.lower() in f"{x[1]['contig']}:{x[1]['start']}-{x[1]['end']}", ecd))  # TODO
    else:
        # TODO: nicer priority
        res = ecd[:10]

    return {
        "results": list(map(
            lambda x: {"i": x[0], "contig": x[1]["contig"], "start": x[1]["start"], "end": x[1]["end"]},
            res)),
    }


@app.route("/call_data/<int:i>")
def get_call_data(i: int):
    cr = app.config["CALL_REPORT"]
    cr_res = cr["results"]
    if i < 0 or i > len(cr_res) - 1:
        raise NotFound()
    return cr_res[i]


# @app.route("/ref")
# def get_ref_file():
#     return send_file(app.config["PARAMS"]["ref"], conditional=True)
#
#
# @app.route("/ref_index")
# def get_ref_index_file():
#     return send_file(app.config["PARAMS"]["ref_index"], conditional=True)


@app.route("/align_files/<int:i>")
def get_align_file(i: int):
    return send_file(app.config["PARAMS"]["align_files"][i], conditional=True)


@app.route("/align_indices/<int:i>")
def get_align_index_file(i: int):
    return send_file(app.config["PARAMS"]["align_indices"][i], conditional=True)


def run_server(call_report, **kwargs):
    app.config.from_mapping(dict(CALL_REPORT=call_report, PARAMS=kwargs))
    app.run(host="localhost", port=5011, debug=True)
