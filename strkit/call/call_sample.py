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
from strkit_rust_ext import STRkitLocus
from queue import Empty as QueueEmpty
from threading import Lock
from typing import Literal

from .call_locus import call_locus
from .loci import load_loci
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
from .utils import get_new_seed, normalize_contig

__all__ = [
    "call_sample",
]


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
    seed: int,
) -> list[LocusResult]:
    if params.profile:
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

    rng = np_default_rng(seed=seed)

    sample_id = params.sample_id

    ref = FastaFile(params.reference_file)
    bf = STRkitBAMReader(
        params.read_file,
        params.reference_file,
        params.max_reads,
        params.skip_supplementary,
        params.skip_secondary,
        params.use_hp,
        lg,
        params.log_level == logging.DEBUG,
    )

    snv_vcf_contigs: list[str] = []
    if params.snv_vcf:
        with VariantFile(str(params.snv_vcf)) as snv_vcf_file:
            snv_vcf_contigs.extend(map(lambda c: c.name, snv_vcf_file.header.contigs.values()))

    vcf_file_format: Literal["chr", "num", "acc", ""] = get_vcf_contig_format(snv_vcf_contigs)

    ref_file_has_chr = any(r.startswith("chr") for r in ref.references)
    read_file_has_chr = any(r.startswith("chr") for r in bf.references)

    snv_vcf_reader = STRkitVCFReader(str(params.snv_vcf)) if params.snv_vcf else None

    current_contig: str | None = None
    results: list[LocusResult] = []

    while True:
        try:
            lb = locus_queue.get_nowait()
            if lb is None:  # Kill signal
                lg.debug(f"worker %d finished current contig: %s", worker_id, current_contig)
                break
            locus_block: tuple[STRkitLocus, ...]
            locus_block_left: int
            locus_block_right: int
            locus_block_left, locus_block_right, locus_block = lb
        except QueueEmpty:
            lg.debug(f"worker %d encountered queue.Empty", worker_id)
            break

        if current_contig is None:
            current_contig = locus_block[0].contig

        # Get all aligned segments for this locus block, reducing the number of round-trip BAM IO ops we have to do!
        block_segments = bf.get_overlapping_segments_and_related_data_for_block(
            normalize_contig(current_contig, read_file_has_chr),
            locus_block_left,
            locus_block_right,
            f"[block {current_contig}:{locus_block_left}-{locus_block_right}]",
        )

        # Find candidate SNVs for this locus block, if we're using SNV data
        block_candidate_snvs = snv_vcf_reader.get_candidate_snvs(
            tuple(snv_vcf_contigs), vcf_file_format, current_contig, locus_block_left, locus_block_right
        ) if snv_vcf_reader is not None else None

        for locus in locus_block:
            # String representation of locus for logging purposes
            locus_log_str: str = f"[w{worker_id}] {sample_id or ''}{' ' if sample_id else ''}{locus.log_str()}"

            lg.debug("%s - working on locus", locus_log_str)

            try:
                res = call_locus(
                    locus,
                    # ---
                    block_segments,
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
                    rng=rng,
                    logger_=lg,
                    locus_log_str=locus_log_str,
                    candidate_snvs=block_candidate_snvs,
                    ref_file_has_chr=ref_file_has_chr,
                )

            except Exception as e:
                import traceback
                res = None
                lg.error(f"{locus_log_str} - encountered exception while genotyping ({locus.n_alleles=}): {repr(e)}")
                lg.error(f"{locus_log_str} - {traceback.format_exc()}")

            if res is not None:
                results.append(res)

        # Instead of updating the locus counter every time a single locus finishes, do it in blocks. This is less
        # accurate (which isn't that important) and saves time with locking/unlocking.
        locus_counter_lock.acquire(timeout=30)
        locus_counter.set(locus_counter.get() + len(locus_block))
        locus_counter_lock.release()

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
    sample_id: str | None,
    start_time: float,
    log_level: int,
    log_progress_interval: int,
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
        if timer >= log_progress_interval:
            last_qsize = qsize  # qsize from {log_progress_interval} seconds ago
            qsize = locus_queue.qsize()

            _log()  # log every {log_progress_interval} seconds

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


def call_sample(
    params: CallParams,
    json_path: str | None = None,
    vcf_path: str | None = None,
    indent_json: bool = False,
    output_tsv: bool = True,
) -> None:
    from strkit.logger import get_main_logger
    logger = get_main_logger()

    # Start the call timer
    start_time = time.perf_counter()

    logger.info(f"strkit_rust_ext version: {importlib.metadata.version('strkit_rust_ext')}")

    logger.info("Starting STR genotyping with parameters:")
    for k, v in params.to_dict().items():
        v_fmt = f'"{v}"' if isinstance(v, str) else v
        logger.info(f"    {k.rjust(22)} = {v_fmt}")

    # Seed the random number generator if a seed is provided, for replicability
    rng: NPRandomGenerator = np_default_rng(seed=params.seed)

    manager: mmg.SyncManager = mp.Manager()

    locus_queue = manager.Queue()
    # Loads actual STRkitLocus definitions into locus_queue, but returns # of loci + the contig set:
    num_loci, contig_set = load_loci(params, locus_queue, logger)

    # ---

    # If we're outputting JSON, write the header
    if json_path is not None:
        output_json_report_header(params, contig_set, json_path, indent_json)

    # If we're outputting a VCF, open the file and write the header
    sample_id_str = params.sample_id or "sample"
    vf: PySamVariantFile | None = None
    if vcf_path is not None:
        vh = build_vcf_header(
            sample_id_str, params.reference_file, partial_phasing=params.snv_vcf is not None or params.use_hp
        )
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
                params.log_progress_interval,
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
            get_new_seed(rng),
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
