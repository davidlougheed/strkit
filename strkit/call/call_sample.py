from __future__ import annotations

import heapq
import logging
import multiprocessing as mp
import multiprocessing.dummy as mpd
import multiprocessing.managers as mmg
import numpy as np
import os
import re
import threading
import time

from datetime import datetime
from multiprocessing.synchronize import Event as EventClass  # For type hinting

from typing import Literal, Optional
from strkit.logger import logger

from .allele import get_n_alleles
from .call_locus import call_locus
from .non_daemonic_pool import NonDaemonicPool
from .params import CallParams
from .output import output_json_report, output_tsv as output_tsv_fn, output_vcf
from .utils import get_new_seed

__all__ = [
    "call_sample",
]


# TODO: Parameterize
LOG_PROGRESS_INTERVAL: int = 120  # seconds
PROFILE_LOCUS_CALLS: bool = False

NUMERAL_CONTIG_PATTERN = re.compile(r"^(\d{1,2}|X|Y)$")


def locus_worker(
    params: CallParams,
    locus_queue: mp.Queue,
    phase_set_lock: threading.Lock,
    phase_set_counter: mmg.ValueProxy,
    phase_set_remap: mmg.DictProxy,
    snv_genotype_update_lock: threading.Lock,
    snv_genotype_cache: mmg.DictProxy,
    is_single_processed: bool,
) -> list[dict]:
    if PROFILE_LOCUS_CALLS:
        import cProfile
        pr = cProfile.Profile()
        pr.enable()
    else:
        pr = None

    import pysam as p

    lg: logging.Logger
    if is_single_processed:
        lg = logger
    else:
        from strkit.logger import create_process_logger
        lg = create_process_logger(os.getpid(), params.log_level)

    ref = p.FastaFile(params.reference_file)
    bfs = tuple(p.AlignmentFile(rf, reference_filename=params.reference_file) for rf in params.read_files)

    snv_vcf_file = p.VariantFile(params.snv_vcf) if params.snv_vcf else None
    snv_vcf_contigs: list[str] = []
    vcf_file_format: Literal["chr", "num", "acc", ""] = ""
    if snv_vcf_file is not None:
        snv_vcf_contigs = list(map(lambda c: c.name, snv_vcf_file.header.contigs.values()))
        if not snv_vcf_contigs or snv_vcf_contigs[0].startswith("chr"):
            vcf_file_format = "chr"
        elif NUMERAL_CONTIG_PATTERN.match(snv_vcf_contigs[0]):
            vcf_file_format = "num"
        else:  # assume stupid accession format
            vcf_file_format = "acc"

    ref_file_has_chr = any(r.startswith("chr") for r in ref.references)
    read_file_has_chr = any(r.startswith("chr") for bf in bfs for r in bf.references)

    results: list[dict] = []

    while True:
        td = locus_queue.get()
        if td is None:  # Kill signal
            break

        t_idx, t, n_alleles, locus_seed = td
        res = call_locus(
            t_idx, t, n_alleles, bfs, ref, params,
            phase_set_lock,
            phase_set_counter,
            phase_set_remap,
            snv_genotype_update_lock,
            snv_genotype_cache,
            seed=locus_seed,
            logger_=lg,
            snv_vcf_file=snv_vcf_file,
            snv_vcf_contigs=tuple(snv_vcf_contigs),
            snv_vcf_file_format=vcf_file_format,
            read_file_has_chr=read_file_has_chr,
            ref_file_has_chr=ref_file_has_chr,
        )

        if res is not None:
            results.append(res)

    if pr:
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
    event: EventClass,
):
    import os
    try:
        os.nice(20)
    except (AttributeError, OSError):
        pass

    from strkit.logger import create_process_logger
    lg = create_process_logger(os.getpid(), log_level)

    def _log(q_size: int):
        try:
            lg.info(
                f"{sample_id}: processed {num_loci - max(q_size, num_workers) + num_workers} loci in "
                f"{(datetime.now() - start_time).total_seconds():.1f} seconds")
        except NotImplementedError:
            pass

    last_qsize = locus_queue.qsize()
    last_qsize_n_stuck = 0
    timer = 0
    while not event.is_set():
        time.sleep(1)
        timer += 1
        if timer >= LOG_PROGRESS_INTERVAL:
            qsize = locus_queue.qsize()
            _log(qsize)
            if last_qsize == locus_queue.qsize():  # stuck
                last_qsize_n_stuck += 1
            if last_qsize_n_stuck >= 3:
                # zombie worker, exit
                return
            timer = 0

    _log(locus_queue.qsize())


def parse_loci_bed(loci_file: str):
    with open(loci_file, "r") as tf:
        yield from (tuple(line.split("\t")) for line in (s.strip() for s in tf) if line and not line.startswith("#"))


def call_sample(
    params: CallParams,
    json_path: Optional[str] = None,
    vcf_path: Optional[str] = None,
    indent_json: bool = False,
    output_tsv: bool = True,
) -> None:
    # Start the call timer
    start_time = datetime.now()

    logger.info(
        f"Starting STR genotyping; sample={params.sample_id}, hq={params.hq}, targeted={params.targeted}, "
        f"HP={params.use_hp}, SNVs={params.snv_vcf is not None}; seed={params.seed}")

    # Seed the random number generator if a seed is provided, for replicability
    rng: np.random.Generator = np.random.default_rng(seed=params.seed)

    manager: mmg.SyncManager = mp.Manager()
    locus_queue = manager.Queue()

    # Add all loci from the BED file to the queue, allowing each job
    # to pull from the queue as it becomes freed up to do so.
    num_loci: int = 0
    # Keep track of all contigs we are processing to speed up downstream Mendelian inheritance analysis.
    contig_set: set[str] = set()
    for t_idx, t in enumerate(parse_loci_bed(params.loci_file), 1):
        contig = t[0]

        n_alleles: Optional[int] = get_n_alleles(2, params.sex_chroms, contig)
        if n_alleles is None:  # Sex chromosome, but we don't have a specified sex chromosome karyotype
            continue  # --> skip the locus

        # We use locus-specific random seeds for replicability, no matter which order
        # the loci are yanked out of the queue / how many processes we have.
        # Tuple of (1-indexed locus index, locus data, locus-specific random seed)
        locus_queue.put((t_idx, t, n_alleles, get_new_seed(rng)))
        contig_set.add(contig)
        num_loci += 1

    # At the end of the queue, add a None value (* the # of processes).
    # When a job encounters a None value, it will terminate.
    for _ in range(params.processes):
        locus_queue.put(None)

    is_single_processed = params.processes == 1
    result_lists = []

    pool_class = mpd.Pool if is_single_processed else NonDaemonicPool
    finish_event = mp.Event()
    with pool_class(params.processes) as p:
        # Start the progress tracking process
        progress_job = mp.Process(
            target=progress_worker,
            daemon=True,
            args=(params.sample_id, start_time, params.log_level, locus_queue, num_loci, params.processes,
                  finish_event))
        progress_job.start()

        # Set up shared memory for cross-locus phasing tasks
        phase_set_lock = manager.Lock()  # Lock for generating phase set integers
        phase_set_counter = manager.Value('I', 1)  # Generates unique phase set IDs
        #  - We need to remap haplotagged BAM phase sets to unique phase set IDs generated from the above counter:
        phase_set_remap = manager.dict()
        #  - When updating SNV genotype cache, we need a lock to make sure we don't get weird phasing artifacts where
        #    the cache is updated simultaneously:
        snv_genotype_update_lock = manager.Lock()
        snv_genotype_cache = manager.dict()  # SNV genotype cache for cross-locus phasing

        # Set up the argument tuple and spin up the jobs
        job_args = (
            params,
            locus_queue,
            phase_set_lock,
            phase_set_counter,
            phase_set_remap,
            snv_genotype_update_lock,
            snv_genotype_cache,
            is_single_processed,
        )
        jobs = [p.apply_async(locus_worker, job_args) for _ in range(params.processes)]

        # Gather the process-specific results for combining.
        for j in jobs:
            result_lists.append(j.get())

        finish_event.set()
        progress_job.join()

    # Merge sorted result lists into single sorted list.
    results: tuple[dict, ...] = tuple(heapq.merge(*result_lists, key=lambda x: x["locus_index"]))

    time_taken = datetime.now() - start_time

    if output_tsv:
        output_tsv_fn(results, has_snv_vcf=params.snv_vcf is not None)

    if json_path:
        output_json_report(
            params,
            time_taken,
            contig_set,
            results,
            json_path,
            indent_json,
        )

    if vcf_path:
        output_vcf(params, results, vcf_path)
