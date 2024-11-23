from __future__ import annotations

import importlib.metadata
import logging
import multiprocessing as mp
import multiprocessing.dummy as mpd
import multiprocessing.managers as mmg
import re
import time

from heapq import merge as heapq_merge
from multiprocessing.synchronize import Event as EventClass  # For type hinting
from numpy.random import Generator as NPRandomGenerator, default_rng as np_default_rng
from operator import itemgetter
from pysam import VariantFile as PySamVariantFile
from queue import Empty as QueueEmpty
from threading import Lock
from typing import Iterable, Literal, Optional

from .allele import get_n_alleles
from .call_locus import call_locus
from .non_daemonic_pool import NonDaemonicPool
from .params import CallParams
from .output import (
    output_json_report_header,
    output_json_report_results,
    output_json_report_footer,
    output_tsv as output_tsv_fn,
    build_vcf_header,
    output_contig_vcf_lines,
)
from .types import VCFContigFormat, LocusResult
from .utils import get_new_seed

__all__ = [
    "call_sample",
]


# TODO: Parameterize
LOG_PROGRESS_INTERVAL: int = 120  # seconds
PROFILE_LOCUS_CALLS: bool = False

NUMERAL_CONTIG_PATTERN = re.compile(r"^(\d{1,2}|X|Y)$")
ACCESSION_PATTERN = re.compile(r"^NC_\d+")

get_locus_index = itemgetter("locus_index")


def get_vcf_contig_format(snv_vcf_contigs: list[str]) -> VCFContigFormat:
    if not snv_vcf_contigs or snv_vcf_contigs[0].startswith("chr"):
        return "chr"
    elif NUMERAL_CONTIG_PATTERN.match(snv_vcf_contigs[0]):
        return "num"
    elif ACCESSION_PATTERN.match(snv_vcf_contigs[0]):
        return "acc"
    else:
        return ""


def locus_worker(
    worker_id: int,
    params: CallParams,
    locus_queue: mp.Queue,
    locus_counter_lock: Lock,
    locus_counter: mmg.ValueProxy,
    phase_set_lock: Lock,
    phase_set_counter: mmg.ValueProxy,
    phase_set_remap: mmg.DictProxy,
    phase_set_synonymous: mmg.DictProxy,
    snv_genotype_update_lock: Lock,
    snv_genotype_cache: mmg.DictProxy,
    is_single_processed: bool,
) -> list[LocusResult]:
    if PROFILE_LOCUS_CALLS:
        import cProfile
        pr = cProfile.Profile()
        pr.enable()
    else:
        pr = None

    from os import getpid
    from pysam import FastaFile, VariantFile
    from strkit_rust_ext import STRkitBAMReader, STRkitVCFReader

    lg: logging.Logger
    if is_single_processed:
        from strkit.logger import get_main_logger
        lg = get_main_logger()
    else:
        from strkit.logger import create_process_logger
        lg = create_process_logger(getpid(), params.log_level)

    sample_id = params.sample_id

    ref = FastaFile(params.reference_file)
    bf = STRkitBAMReader(params.read_file, params.reference_file)

    snv_vcf_contigs: list[str] = []
    if params.snv_vcf:
        with VariantFile(params.snv_vcf) as snv_vcf_file:
            snv_vcf_contigs.extend(map(lambda c: c.name, snv_vcf_file.header.contigs.values()))

    vcf_file_format: Literal["chr", "num", "acc", ""] = get_vcf_contig_format(snv_vcf_contigs)

    ref_file_has_chr = any(r.startswith("chr") for r in ref.references)
    read_file_has_chr = any(r.startswith("chr") for r in bf.references)

    snv_vcf_reader = STRkitVCFReader(str(params.snv_vcf)) if params.snv_vcf else None

    current_contig: Optional[str] = None
    results: list[LocusResult] = []

    while True:
        try:
            td = locus_queue.get_nowait()
            if td is None:  # Kill signal
                lg.debug(f"worker %d finished current contig: %s", worker_id, current_contig)
                break
        except QueueEmpty:
            lg.debug(f"worker %d encountered queue.Empty", worker_id)
            break

        t_idx, contig, left_coord, right_coord, motif, n_alleles, locus_seed = td

        if current_contig is None:
            current_contig = contig

        # String representation of locus for logging purposes
        locus_log_str: str = (
            f"[w{worker_id}] {sample_id or ''}{' ' if sample_id else ''}locus {t_idx}: "
            f"{contig}:{left_coord}-{right_coord} [{motif}]"
        )

        lg.debug("%s - working on locus", locus_log_str)

        try:
            res = call_locus(
                t_idx,
                contig,
                left_coord,
                right_coord,
                motif,
                # ---
                n_alleles,
                bf,
                ref,
                params,
                # ---
                phase_set_lock,
                phase_set_counter,
                phase_set_remap,
                phase_set_synonymous,
                snv_genotype_update_lock,
                snv_genotype_cache,
                # ---
                seed=locus_seed,
                logger_=lg,
                locus_log_str=locus_log_str,
                snv_vcf_file=snv_vcf_reader,
                snv_vcf_contigs=tuple(snv_vcf_contigs),
                snv_vcf_file_format=vcf_file_format,
                read_file_has_chr=read_file_has_chr,
                ref_file_has_chr=ref_file_has_chr,
            )

        except Exception as e:
            import traceback
            res = None
            lg.error(f"{locus_log_str} - encountered exception while genotyping ({t_idx=}, {n_alleles=}): {repr(e)}")
            lg.error(f"{locus_log_str} - {traceback.format_exc()}")

        locus_counter_lock.acquire(timeout=30)
        locus_counter.set(locus_counter.get() + 1)
        locus_counter_lock.release()

        if res is not None:
            results.append(res)

    ref.close()
    del bf
    del snv_vcf_reader

    if pr:
        pr.disable()
        pr.print_stats("tottime")

    lg.debug(f"worker {worker_id} - returning batch of {len(results)} locus results")

    # Sort worker results; we will merge them after
    return results if is_single_processed else sorted(results, key=get_locus_index)


def progress_worker(
    sample_id: Optional[str],
    start_time: float,
    log_level: int,
    locus_queue: mp.Queue,
    locus_counter: mmg.ValueProxy,
    num_loci: int,
    event: EventClass,
):
    from os import nice as os_nice, getpid
    try:
        os_nice(20)
    except (AttributeError, OSError):
        pass

    from strkit.logger import create_process_logger
    lg = create_process_logger(getpid(), log_level)

    def _log():
        try:
            processed_loci = int(locus_counter.get())
            n_seconds = time.perf_counter() - start_time
            loci_per_second = processed_loci / n_seconds
            est_time_remaining = ((num_loci - processed_loci) / loci_per_second) if loci_per_second else float("inf")
            lg.info(
                f"{sample_id}: processed {processed_loci}/{num_loci} loci ({processed_loci / num_loci * 100:.1f}%) in "
                f"{n_seconds:.1f} seconds (~{processed_loci / n_seconds:.0f} l/s; est. time remaining: "
                f"{est_time_remaining:.0f}s)")
        except NotImplementedError:
            pass

    qsize = locus_queue.qsize()
    last_qsize_n_stuck: int = 0
    timer: int = 0
    while not event.is_set():
        time.sleep(1)
        timer += 1  # one second has elapsed
        if timer >= LOG_PROGRESS_INTERVAL:
            last_qsize = qsize  # qsize from {LOG_PROGRESS_INTERVAL} seconds ago
            qsize = locus_queue.qsize()

            _log()  # log every {LOG_PROGRESS_INTERVAL} seconds

            if last_qsize == qsize:  # stuck
                last_qsize_n_stuck += 1
            else:
                last_qsize_n_stuck = 0

            if last_qsize_n_stuck >= 30:
                # zombie worker or stuck, exit
                lg.error(f"Terminating progress worker; seems to be stuck with {qsize=}")
                return

            timer = 0  # reset timer

    _log()


def parse_loci_bed(loci_file: str) -> Iterable[tuple[str, ...]]:
    with open(loci_file, "r") as tf:
        yield from (
            tuple(line.split("\t"))
            for line in (s.strip() for s in tf)
            if line and not line.startswith("#")  # skip blank lines and comment lines
        )


def call_sample(
    params: CallParams,
    json_path: Optional[str] = None,
    vcf_path: Optional[str] = None,
    indent_json: bool = False,
    output_tsv: bool = True,
) -> None:
    from strkit.logger import get_main_logger
    logger = get_main_logger()

    # Start the call timer
    start_time = time.perf_counter()

    logger.info(f"strkit_rust_ext version: {importlib.metadata.version('strkit_rust_ext')}")
    logger.info(
        f"Starting STR genotyping; sample={params.sample_id}, hq={params.hq}, targeted={params.targeted}, "
        f"HP={params.use_hp}, SNVs={params.snv_vcf is not None}; seed={params.seed}")

    # Seed the random number generator if a seed is provided, for replicability
    rng: NPRandomGenerator = np_default_rng(seed=params.seed)

    manager: mmg.SyncManager = mp.Manager()
    locus_queue = manager.Queue()  # TODO: one queue per contig?

    # Add all loci from the BED file to the queue, allowing each job
    # to pull from the queue as it becomes freed up to do so.
    num_loci: int = 0
    # Keep track of all contigs we are processing to speed up downstream Mendelian inheritance analysis.
    contig_set: set[str] = set()
    last_contig: Optional[str] = None
    last_none_append_n_loci: int = 0
    for t_idx, t in enumerate(parse_loci_bed(params.loci_file), 1):
        contig = t[0]

        n_alleles: Optional[int] = get_n_alleles(2, params.sex_chroms, contig)
        if n_alleles is None:  # Sex chromosome, but we don't have a specified sex chromosome karyotype
            last_contig = contig
            continue  # --> skip the locus

        # For each contig, add None breaks to make the jobs terminate. Then, as long as we have entries in the locus
        # queue, we can get contig-chunks to order and write to disk instead of keeping everything in RAM.
        if last_contig != contig:
            if last_contig is not None:
                for _ in range(params.processes):
                    locus_queue.put(None)
                last_none_append_n_loci = num_loci

            contig_set.add(contig)
            last_contig = contig

        # We use locus-specific random seeds for replicability, no matter which order
        # the loci are yanked out of the queue / how many processes we have.
        # Tuple of (1-indexed locus index, contig, left coord, right coord, motif, locus-specific random seed)
        locus_queue.put((t_idx, contig, int(t[1]), int(t[2]), t[-1], n_alleles, get_new_seed(rng)))
        num_loci += 1

    if num_loci > last_none_append_n_loci:  # have progressed since the last None-append
        # At the end of the queue, add a None value (* the # of processes).
        # When a job encounters a None value, it will terminate.
        for _ in range(params.processes):
            locus_queue.put(None)

    logger.info(f"Loaded {num_loci} loci")

    # ---

    # If we're outputting JSON, write the header
    if json_path is not None:
        output_json_report_header(params, contig_set, json_path, indent_json)

    # If we're outputting a VCF, open the file and write the header
    sample_id_str = params.sample_id or "sample"
    vf: Optional[PySamVariantFile] = None
    if vcf_path is not None:
        vh = build_vcf_header(sample_id_str, params.reference_file)
        vf = PySamVariantFile(vcf_path if vcf_path != "stdout" else "-", "w", header=vh)

    # ---

    is_single_processed = params.processes == 1

    if not is_single_processed:
        logger.info(f"Using {params.processes} workers")

    pool_class = mpd.Pool if is_single_processed else NonDaemonicPool
    finish_event = mp.Event()
    with pool_class(params.processes) as p:
        locus_counter_lock = manager.Lock()
        locus_counter = manager.Value("I", 0)

        # Start the progress tracking process
        progress_job = mp.Process(
            target=progress_worker,
            daemon=True,
            args=(
                params.sample_id,
                start_time,
                params.log_level,
                locus_queue,
                locus_counter,
                num_loci,
                finish_event,
            ))
        progress_job.start()

        # Set up shared memory for cross-locus phasing tasks
        phase_set_lock = manager.Lock()  # Lock for generating phase set integers
        phase_set_counter = manager.Value("I", 1)  # Generates unique phase set IDs
        #  - We need to remap haplotagged BAM phase sets to unique phase set IDs generated from the above counter:
        phase_set_remap = manager.dict()
        #  - Some phase sets end up being in truth the same when parallel processing; we need to reconcile this at the
        #    end. We can do this by keeping track of synonymous phase sets.
        phase_set_synonymous = manager.dict()
        #  - When updating SNV genotype cache, we need a lock to make sure we don't get weird phasing artifacts where
        #    the cache is updated simultaneously:
        snv_genotype_update_lock = manager.Lock()
        snv_genotype_cache = manager.dict()  # SNV genotype cache for cross-locus phasing

        # Set up the argument tuple and spin up the jobs
        job_args = (
            params,
            locus_queue,
            locus_counter_lock,
            locus_counter,
            phase_set_lock,
            phase_set_counter,
            phase_set_remap,
            phase_set_synonymous,
            snv_genotype_update_lock,
            snv_genotype_cache,
            is_single_processed,
        )

        n_contigs_processed: int = 0
        qsize: int = locus_queue.qsize()
        while qsize > 0:  # should have one loop iteration per contig
            jobs = [p.apply_async(locus_worker, (i + 1, *job_args)) for i in range(params.processes)]
            # jobs = [locus_worker(i + 1, *job_args) for i in range(params.processes)]

            # Write results
            #  - gather the process-specific results for combining
            #  - merge sorted result lists into single sorted list for the current contig
            results: tuple[LocusResult, ...] = tuple(heapq_merge(*(j.get() for j in jobs), key=get_locus_index))
            # results: tuple[LocusResult, ...] = tuple(heapq_merge(*jobs, key=get_locus_index))
            del jobs

            #  - fix-up phase sets based on phase_set_synonymous
            #     - if determined necessary while traversing the phase set graph, flip all peak data (i.e., if phase
            #       sets are synonymous but for phasing to work we need to flip some of the loci.)
            if phase_set_synonymous:
                r: LocusResult
                for r in filter(lambda rr: rr.get("ps") is not None, results):
                    should_flip = False

                    while r["ps"] in phase_set_synonymous:
                        r_ps, r_sf = phase_set_synonymous[r["ps"]]
                        r["ps"] = r_ps
                        should_flip = (not should_flip) if r_sf else should_flip

                    if should_flip and r["call"]:
                        # Need to flip EVERYTHING which is phased (calls, reads, SNVs, peaks):
                        r["call"].reverse()
                        r["call_95_cis"].reverse()
                        if r["call_99_cis"]:
                            r["call_99_cis"].reverse()
                        if "reads" in r:
                            for k in r["reads"]:
                                if "p" in r["reads"][k]:
                                    r["reads"][k]["p"] = int(not bool(r["reads"][k]["p"]))
                        if "snvs" in r:
                            for s in r["snvs"]:
                                s["call"] = tuple(reversed(s["call"]))
                                s["rcs"].reverse()
                        if r.get("peaks"):
                            for k in r["peaks"]:
                                if k in {"means", "weights", "stdevs", "n_reads", "kmers", "seqs"}:  # peak/list keys
                                    r["peaks"][k] = r["peaks"][k][::-1]

            #  - write partial results to stdout if we're writing a stdout TSV
            if output_tsv:
                output_tsv_fn(results, has_snv_vcf=params.snv_vcf is not None)

            if json_path is not None:
                output_json_report_results(results, n_contigs_processed == len(contig_set) - 1, json_path, indent_json)

            #  - write partial results to VCF if we're writing a VCF - this must be a whole contig's worth of calls, or
            #    at least the catalog's perception of a whole contig (i.e., all regions from one contig)
            if vf is not None:
                output_contig_vcf_lines(params, sample_id_str, vf, results, logger)

            #  - we're done with this tuple, so delete it as early as possible
            del results

            #  - clean up caches, since we're changing contigs. We don't need the locks, since at this point only the
            #    main process is operational.

            phase_set_synonymous.clear()
            phase_set_remap.clear()
            snv_genotype_cache.clear()

            #  - mark that we've processed another contig (equivalent to a batch):
            n_contigs_processed += 1

            #  - check that we're progressing:

            last_qsize = qsize
            qsize = locus_queue.qsize()
            if last_qsize == qsize:
                # If this happens, we've stalled out on completing work while having a
                # positive qsize, so we need to terminate (but something went wrong.)
                logger.warning(f"Terminating with non-zero queue size: {qsize}")
                break

        finish_event.set()
        progress_job.join()

    time_taken = time.perf_counter() - start_time

    logger.info(f"Finished STR genotyping in {time_taken:.1f}s")

    if json_path:
        output_json_report_footer(time_taken, json_path, indent_json)
