from __future__ import annotations

import heapq
import logging
import multiprocessing as mp
import multiprocessing.dummy as mpd
import numpy as np
import os
import sys
import time

from datetime import datetime
from pysam import AlignmentFile

from typing import Optional, Union

from strkit import __version__

from strkit.exceptions import ParamError
from strkit.json import dumps_indented, json
from strkit.logger import logger

from .call_locus import call_locus

__all__ = [
    "call_sample",
]


# TODO: Parameterize
LOG_PROGRESS_INTERVAL = 120  # seconds
PROFILE_LOCUS_CALLS = False


def locus_worker(
    read_files: tuple[str, ...],
    reference_file: str,
    min_reads: int,
    min_allele_reads: int,
    min_avg_phred: int,
    num_bootstrap: int,
    flank_size: int,
    sex_chroms: Optional[str],
    realign: bool,
    hq: bool,
    incorporate_snvs: bool,
    targeted: bool,
    fractional: bool,
    respect_ref: bool,
    count_kmers: str,
    log_level: int,
    locus_queue: mp.Queue,
    is_single_processed: bool,
) -> list[dict]:
    if PROFILE_LOCUS_CALLS:
        import cProfile
        pr = cProfile.Profile()
        pr.enable()
    else:
        pr = None

    import pysam as p

    if is_single_processed:
        lg = logger
    else:
        from strkit.logger import create_process_logger
        lg = create_process_logger(os.getpid(), log_level)

    ref = p.FastaFile(reference_file)
    bfs = tuple(p.AlignmentFile(rf, reference_filename=reference_file) for rf in read_files)

    ref_file_has_chr = any(r.startswith("chr") for r in ref.references)
    read_file_has_chr = any(r.startswith("chr") for bf in bfs for r in bf.references)

    results: list[dict] = []

    while True:
        td = locus_queue.get()
        if td is None:  # Kill signal
            break

        t_idx, t, locus_seed = td
        res = call_locus(
            t_idx, t, bfs, ref,
            min_reads=min_reads,
            min_allele_reads=min_allele_reads,
            min_avg_phred=min_avg_phred,
            num_bootstrap=num_bootstrap,
            flank_size=flank_size,
            seed=locus_seed,
            logger_=lg,
            sex_chroms=sex_chroms,
            realign=realign,
            hq=hq,
            incorporate_snvs=incorporate_snvs,
            targeted=targeted,
            fractional=fractional,
            respect_ref=respect_ref,
            count_kmers=count_kmers,
            log_level=log_level,
            read_file_has_chr=read_file_has_chr,
            ref_file_has_chr=ref_file_has_chr,
        )

        if res is not None:
            results.append(res)

    if PROFILE_LOCUS_CALLS:
        pr.disable()
        pr.print_stats("tottime")

    # Sort worker results; we will merge them after
    return results if is_single_processed else sorted(results, key=lambda x: x["locus_index"])


def progress_worker(
    sample_id: Optional[str],
    start_time: datetime,
    log_level: int,
    locus_queue: mp.Queue,
    num_loci: int,
    num_workers: int,
    event: mp.Event,
):
    import os
    try:
        os.nice(20)
    except (AttributeError, OSError):
        pass

    from strkit.logger import create_process_logger
    lg = create_process_logger(os.getpid(), log_level)

    def _log():
        try:
            lg.info(
                f"{sample_id}: processed {num_loci - max(locus_queue.qsize(), num_workers) + num_workers} loci in "
                f"{round((datetime.now() - start_time).total_seconds())} seconds")
        except NotImplementedError:
            pass

    timer = 0
    while not event.is_set():
        time.sleep(1)
        timer += 1
        if timer >= LOG_PROGRESS_INTERVAL:
            _log()
            timer = 0

    _log()


def parse_loci_bed(loci_file: str):
    with open(loci_file, "r") as tf:
        yield from (tuple(line.split("\t")) for line in (s.strip() for s in tf) if line and not line.startswith("#"))


def _cn_to_str(cn: Union[int, float]) -> str:
    return f"{cn:.1f}" if isinstance(cn, float) else str(cn)


def call_sample(
    read_files: tuple[str, ...],
    reference_file: str,
    loci_file: str,
    sample_id: Optional[str],
    min_reads: int = 4,
    min_allele_reads: int = 2,
    min_avg_phred: int = 13,
    num_bootstrap: int = 100,
    flank_size: int = 70,
    sex_chroms: Optional[str] = None,
    realign: bool = False,
    hq: bool = False,
    incorporate_snvs: bool = False,
    targeted: bool = False,
    fractional: bool = False,
    respect_ref: bool = False,
    count_kmers: str = "none",  # "none" | "peak" | "read"
    log_level: int = logging.WARNING,
    json_path: Optional[str] = None,
    indent_json: bool = False,
    output_tsv: bool = True,
    processes: int = 1,
    seed: Optional[int] = None,
) -> None:
    # Check parameter validity

    if incorporate_snvs and not hq:
        raise ParamError("Cannot use --incorporate-snvs without --hq.")

    # Start the call timer
    start_time = datetime.now()

    # Get sample ID from read file(s)

    bfs = tuple(AlignmentFile(rf, reference_filename=reference_file) for rf in read_files)

    # noinspection PyTypeChecker
    bfhs = [dict(bf.header) for bf in bfs]

    sns: set[str] = {e.get("SM") for bfh in bfhs for e in bfh.get("RG", ()) if e.get("SM")}
    bam_sample_id: Optional[str] = None

    if len(sns) > 1:
        # Error or warning or what?
        sns_str = "', '".join(sns)
        logger.warning(f"Found more than one sample ID in BAM file(s): '{sns_str}'")
    elif not sns:
        logger.warning("Could not find sample ID in BAM file(s)")
    else:
        bam_sample_id = sns.pop()

    sample_id_final: Optional[str] = sample_id or bam_sample_id

    # Seed the random number generator if a seed is provided, for replicability
    rng = np.random.default_rng(seed=seed)

    manager = mp.Manager()
    locus_queue = manager.Queue()

    # Add all loci from the BED file to the queue, allowing each job
    # to pull from the queue as it becomes freed up to do so.
    num_loci = 0
    for t_idx, t in enumerate(parse_loci_bed(loci_file), 1):
        # We use locus-specific random seeds for replicability, no matter which order
        # the loci are yanked out of the queue / how many processes we have.
        # Tuple of (1-indexed locus index, locus data, locus-specific random seed)
        locus_queue.put((t_idx, t, rng.integers(0, 4096).item()))
        num_loci += 1

    # At the end of the queue, add a None value (* the # of processes).
    # When a job encounters a None value, it will terminate.
    for _ in range(processes):
        locus_queue.put(None)

    # Order matters here!!
    job_params = {
        "read_files": read_files,
        "reference_file": reference_file,
        "min_reads": min_reads,
        "min_allele_reads": min_allele_reads,
        "min_avg_phred": min_avg_phred,
        "num_bootstrap": num_bootstrap,
        "flank_size": flank_size,
        "sex_chroms": sex_chroms,
        "realign": realign,
        "hq": hq,
        "incorporate_snvs": incorporate_snvs,
        "targeted": targeted,
        "fractional": fractional,
        "respect_ref": respect_ref,
        "count_kmers": count_kmers,
        "log_level": log_level,
    }

    is_single_processed = processes == 1
    job_args = (*job_params.values(), locus_queue, is_single_processed)
    result_lists = []

    pool_class = mpd.Pool if is_single_processed else mp.Pool
    finish_event = mp.Event()
    with pool_class(processes) as p:
        # Start the progress tracking process
        progress_job = mp.Process(
            target=progress_worker,
            args=(sample_id_final, start_time, log_level, locus_queue, num_loci, processes, finish_event))
        progress_job.start()

        # Spin up the jobs
        jobs = [p.apply_async(locus_worker, job_args) for _ in range(processes)]

        # Gather the process-specific results for combining.
        for j in jobs:
            result_lists.append(j.get())

        finish_event.set()
        progress_job.join()

    # Merge sorted result lists into single sorted list.
    results: tuple[dict, ...] = tuple(heapq.merge(*result_lists, key=lambda x: x["locus_index"]))

    time_taken = datetime.now() - start_time

    if output_tsv:
        for res in results:
            has_call = res["call"] is not None
            # n_peaks = res["peaks"]["modal_n"]

            sys.stdout.write("\t".join((
                res["contig"],
                str(res["start"]),
                str(res["end"]),
                res["motif"],
                _cn_to_str(res["ref_cn"]),
                ",".join(map(_cn_to_str, sorted(r["cn"] for r in res["reads"].values()))),
                "|".join(map(_cn_to_str, res["call"])) if has_call else ".",
                ("|".join("-".join(map(_cn_to_str, gc)) for gc in res["call_95_cis"]) if has_call else "."),

                # ("|".join(map(lambda x: f"{x:.5f}", res["peaks"]["means"][:n_peaks]))
                #  if has_call and n_peaks <= 2 else "."),
                # ("|".join(map(lambda x: f"{x:.5f}", res["peaks"]["weights"][:n_peaks]))
                #  if has_call and n_peaks <= 2 else "."),
                # ("|".join(map(lambda x: f"{x:.5f}", res["peaks"]["stdevs"][:n_peaks]))
                #  if has_call and n_peaks <= 2 else "."),
            )) + "\n")

    if json_path:
        json_report = {
            "sample_id": sample_id_final,
            "caller": {
                "name": "strkit",
                "version": __version__,
            },
            "parameters": {
                **job_params,
                "processes": processes,
            },
            "runtime": time_taken.total_seconds(),
            "results": results,
        }

        dfn = dumps_indented if indent_json else json.dumps
        if json_path == "stdout":
            sys.stdout.buffer.write(dfn(json_report))
            sys.stdout.write("\n")
            sys.stdout.flush()
        else:
            with open(json_path, "wb") as jf:
                jf.write(dfn(json_report))
                jf.write(b"\n")
