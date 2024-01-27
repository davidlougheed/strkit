from __future__ import annotations

import heapq
import logging
import multiprocessing as mp
import multiprocessing.dummy as mpd
import pathlib

import numpy as np
import os
import re
import time

from datetime import datetime
from multiprocessing.synchronize import Event as EventClass  # For type hinting
from pysam import AlignmentFile

from typing import Literal, Optional
from strkit.logger import logger

from .allele import get_n_alleles
from .call_locus import call_locus
from .non_daemonic_pool import NonDaemonicPool
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
    read_files: tuple[str, ...],
    reference_file: str,
    min_reads: int,
    min_allele_reads: int,
    min_avg_phred: int,
    num_bootstrap: int,
    flank_size: int,
    sample_id: Optional[str],
    realign: bool,
    hq: bool,
    # incorporate_snvs: bool,
    snv_vcf: Optional[str],
    targeted: bool,
    fractional: bool,
    respect_ref: bool,
    count_kmers: str,
    consensus: bool,
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

    lg: logging.Logger
    if is_single_processed:
        lg = logger
    else:
        from strkit.logger import create_process_logger
        lg = create_process_logger(os.getpid(), log_level)

    ref = p.FastaFile(reference_file)
    bfs = tuple(p.AlignmentFile(rf, reference_filename=reference_file) for rf in read_files)

    snv_vcf_file = p.VariantFile(snv_vcf) if snv_vcf else None
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
            t_idx, t, n_alleles, bfs, ref,
            min_reads=min_reads,
            min_allele_reads=min_allele_reads,
            min_avg_phred=min_avg_phred,
            num_bootstrap=num_bootstrap,
            flank_size=flank_size,
            seed=locus_seed,
            logger_=lg,
            sample_id=sample_id,
            realign=realign,
            hq=hq,
            # incorporate_snvs=incorporate_snvs,
            snv_vcf_file=snv_vcf_file,
            snv_vcf_contigs=tuple(snv_vcf_contigs),
            snv_vcf_file_format=vcf_file_format,
            targeted=targeted,
            fractional=fractional,
            respect_ref=respect_ref,
            count_kmers=count_kmers,
            consensus=consensus,
            log_level=log_level,
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
    # incorporate_snvs: bool = False,
    snv_vcf: Optional[pathlib.Path] = None,
    targeted: bool = False,
    fractional: bool = False,
    respect_ref: bool = False,
    count_kmers: str = "none",  # "none" | "peak" | "read"
    consensus: bool = False,
    log_level: int = logging.WARNING,
    json_path: Optional[str] = None,
    vcf_path: Optional[str] = None,
    indent_json: bool = False,
    output_tsv: bool = True,
    processes: int = 1,
    seed: Optional[int] = None,
) -> None:
    # Start the call timer
    start_time = datetime.now()

    # Get sample ID from read file(s)

    logger.debug(f"Using read file(s): {'; '.join(read_files)}")
    bfs = tuple(AlignmentFile(rf, reference_filename=reference_file) for rf in read_files)

    # noinspection PyTypeChecker
    bfhs = [bf.header.to_dict() for bf in bfs]

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

    logger.info(
        f"Starting STR genotyping; sample={sample_id_final}, hq={hq}, targeted={targeted}, SNVs={snv_vcf is not None}; "
        f"seed={seed}")

    # Seed the random number generator if a seed is provided, for replicability
    rng: np.random.Generator = np.random.default_rng(seed=seed)

    manager = mp.Manager()
    locus_queue = manager.Queue()

    # Add all loci from the BED file to the queue, allowing each job
    # to pull from the queue as it becomes freed up to do so.
    num_loci: int = 0
    # Keep track of all contigs we are processing to speed up downstream Mendelian inheritance analysis.
    contig_set: set[str] = set()
    for t_idx, t in enumerate(parse_loci_bed(loci_file), 1):
        contig = t[0]

        n_alleles: Optional[int] = get_n_alleles(2, sex_chroms, contig)
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
        "sample_id": sample_id_final,
        "realign": realign,
        "hq": hq,
        # "incorporate_snvs": incorporate_snvs,
        "snv_vcf": str(snv_vcf) if snv_vcf else None,
        "targeted": targeted,
        "fractional": fractional,
        "respect_ref": respect_ref,
        "count_kmers": count_kmers,
        "consensus": consensus,
        "log_level": log_level,
    }

    is_single_processed = processes == 1
    job_args = (*job_params.values(), locus_queue, is_single_processed)
    result_lists = []

    pool_class = mpd.Pool if is_single_processed else NonDaemonicPool
    finish_event = mp.Event()
    with pool_class(processes) as p:
        # Start the progress tracking process
        progress_job = mp.Process(
            target=progress_worker,
            daemon=True,
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
        output_tsv_fn(results, has_snv_vcf=snv_vcf is not None)

    if json_path:
        output_json_report(
            sample_id_final,
            job_params,
            processes,
            time_taken,
            contig_set,
            results,
            json_path,
            indent_json,
        )

    if vcf_path:
        output_vcf(sample_id_final, reference_file, results, vcf_path)
