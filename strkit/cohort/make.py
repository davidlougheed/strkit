from __future__ import annotations

import os
import pathlib
import sqlite3

from strkit.exceptions import ParamError
from strkit.json import json

SCHEMA_PATH = pathlib.Path(__file__).parent / "schema.sql"


def run_schema(conn: sqlite3.Connection):
    with conn, open(SCHEMA_PATH, "r") as sf:
        conn.executescript(sf.read())


def ingest_json(sample_files: list[str]):
    # TODO: Initialize tr_ids and sample_ids from database if it already exists

    tr_id_it = 1
    tr_ids = {}

    sample_id_it = 1
    sample_ids = {}

    for sf_path in sample_files:
        if sf_path in sample_ids:
            s_id = sample_ids[sf_path]
        else:
            sample_ids[sf_path] = (s_id := sample_id_it)
            sample_id_it += 1

        with open(sf_path, "r") as sf:
            data = json.loads(sf.read())

        for res in data["results"]:
            k = (res["contig"], res["start"], res["end"], res["motif"])
            if k in tr_ids:
                res_id = tr_ids[k]
            else:
                tr_ids[k] = (res_id := tr_id_it)
                tr_id_it += 1


def cohort_make(sample_format: str, sample_files: list[str], out_file: str):
    if not sample_files:
        raise ParamError("No sample files specified")

    if sample_format not in ("tsv", "json"):
        raise ParamError(f"Invalid sample format: '{sample_format}'")

    if os.path.exists(out_file):
        raise ParamError(f"File already exists: '{out_file}'")

    conn = sqlite3.connect(out_file)

    try:
        run_schema(conn)

        if sample_format == "tsv":
            pass
        elif sample_files == "json":
            ingest_json(sample_files)
        else:
            pass

        # TODO: add data to tables

    finally:
        conn.close()

    # TODO: want to enable
    #  - small/large allele (or single allele) distribution comparisons
    #     - somehow testing both AD and AR type disorders ? maybe combined is enough here
    #  - small/large allele motif composition differences
    #  - maybe calculate distributions and histograms for each one and allow this to be visualized?
