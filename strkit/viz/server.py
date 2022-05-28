from flask import Flask, send_file, render_template
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


@app.route("/call_data/<int:i>")
def get_call_data(i: int):
    cd = app.config["CALL_DATA"]
    if i < 0 or i > len(cd) - 1:
        raise NotFound()
    return app.config["CALL_DATA"][i]


# @app.route("/ref")
# def get_ref_file():
#     return send_file(app.config["PARAMS"]["ref"], conditional=True)
#
#
# @app.route("/ref_index")
# def get_ref_index_file():
#     return send_file(app.config["PARAMS"]["ref_index"], conditional=True)


@app.route("/align")
def get_align_file():
    return send_file(app.config["PARAMS"]["align"], conditional=True)


@app.route("/align_index")
def get_align_index_file():
    return send_file(app.config["PARAMS"]["align_index"], conditional=True)


def run_server(call_data, **kwargs):
    app.config.from_mapping(dict(CALL_DATA=call_data, PARAMS=kwargs))
    app.run(host="localhost", port=5011, debug=True)
