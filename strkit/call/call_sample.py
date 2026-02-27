from __future__ import annotations

from operator import itemgetter
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # For type hinting only
    from logging import Logger
    from multiprocessing import Queue as MpQueue
    # noinspection PyProtectedMember
    from multiprocessing.managers import DictProxy, SyncManager, ValueProxy
    from multiprocessing.synchronize import Event as EventClass
    from numpy.random import Generator as NPRandomGenerator
    from strkit_rust_ext import STRkitLocus, STRkitLocusBlock
    from threading import Lock

    from .params import CallParams
    from .types import LocusResult


__all__ = [
    "call_sample",
]

get_locus_index = itemgetter("locus_index")


def locus_worker(
    worker_id: int,
    params: CallParams,
    locus_queue: MpQueue,
    locus_counter_lock: Lock,
    locus_counter: ValueProxy,
    phase_set_lock: Lock,
    phase_set_counter: ValueProxy,
    phase_set_remap: DictProxy,
    phase_set_synonymous: DictProxy,
    snv_genotype_update_lock: Lock,
    snv_genotype_cache: DictProxy,
    is_single_processed: bool,
    seed: int,
) -> list[LocusResult]:
    from logging import DEBUG
    from numpy.random import default_rng as np_default_rng
    from os import getpid
    from pysam import FastaFile
    from queue import Empty as QueueEmpty
    from strkit_rust_ext import STRkitBAMReader, STRkitVCFReader

    lg: Logger
    if is_single_processed:
        from strkit.logger import get_main_logger
        lg = get_main_logger()
    else:
        from strkit.logger import create_process_logger
        lg = create_process_logger(getpid(), params.log_level)

    if params.profile:
        import cProfile
        pr = cProfile.Profile()
        pr.enable()
        lg.info("Profiling is enabled")
    else:
        pr = None

    from .call_locus import call_locus

    # ------------------------------------------------------------------------------------------------------------------

    rng = np_default_rng(seed=seed)

    sample_id = params.sample_id
    log_str_prefix: str = f"[w{worker_id}] {sample_id or ''}{' ' if sample_id else ''}"

    ref = FastaFile(params.reference_file)
    bf = STRkitBAMReader(
        params.read_file,
        params.reference_file,
        params.max_reads,
        params.skip_supplementary,
        params.skip_secondary,
        params.use_hp,
        params.significant_clip_threshold,
        lg,
        params.log_level == DEBUG,
    )

    ref_file_has_chr = any(r.startswith("chr") for r in ref.references)

    snv_vcf_reader = STRkitVCFReader(
        str(params.snv_vcf), is_sample_vcf=False, sample_id=None
    ) if params.snv_vcf else None

    current_contig: str | None = None
    block_results: list[LocusResult] = []  # in the future, maybe we can use this for some kind of extra phasing info
    results: list[LocusResult] = []

    while True:
        try:
            locus_block: STRkitLocusBlock = locus_queue.get_nowait()
            if locus_block is None:  # Kill signal
                lg.debug("worker %d finished current contig: %s", worker_id, current_contig)
                break
        except QueueEmpty:
            lg.debug("worker %d encountered queue.Empty", worker_id)
            break

        # TODO: locus block ID
        locus_block_size = len(locus_block)
        lg.debug("working on locus block (len=%d)", locus_block_size)

        if current_contig is None:
            current_contig = locus_block.contig

        # Get all aligned segments for this locus block, reducing the number of round-trip BAM IO ops we have to do!
        block_segments = bf.get_overlapping_segments_and_related_data_for_block(locus_block)

        # Find candidate SNVs for this locus block, if we're using SNV data
        block_candidate_snvs = snv_vcf_reader.get_candidate_snvs(locus_block) if snv_vcf_reader is not None else None

        # TODO: could do SNV pileup for this block in the future.
        #  - this should improve SNV genotyping, as we can get the true genotypes of the candidate SNVs in the sample
        #    before we extract the values for the reads.

        locus: STRkitLocus
        for locus in locus_block:
            # String representation of locus for logging purposes
            locus_log_str: str = log_str_prefix + locus.log_str()

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
                res = None
                lg.exception(
                    "%s - encountered exception while genotyping (n_alleles=%d)",
                    locus_log_str, locus.n_alleles, exc_info=e)

            if res is not None:
                block_results.append(res)

        results.extend(block_results)

        # Done with these
        block_results.clear()
        del locus_block
        del block_segments
        del block_candidate_snvs

        # Instead of updating the locus counter every time a single locus finishes, do it in blocks. This is less
        # accurate (which isn't that important) and saves time with locking/unlocking.
        locus_counter_lock.acquire(timeout=30)
        locus_counter.set(locus_counter.get() + locus_block_size)
        locus_counter_lock.release()

    ref.close()
    del bf
    del snv_vcf_reader

    if pr:
        pr.disable()
        pr.print_stats("tottime")

    lg.debug("worker %d - returning batch of %d locus results", worker_id, len(results))

    if not is_single_processed:
        # If we're multiprocessing, sort worker results in-place; we will merge them after
        results.sort(key=get_locus_index)

    return results


def progress_worker(
    sample_id: str | None,
    start_time: float,
    log_level: int,
    log_progress_interval: int,
    locus_queue: MpQueue,
    locus_counter: ValueProxy,
    num_loci: int,
    event: EventClass,
):
    import time

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
    import importlib.metadata
    import multiprocessing as mp
    import time

    from heapq import merge as heapq_merge
    from numpy.random import default_rng as np_default_rng
    from pathlib import Path
    from pysam import VariantFile

    from strkit.logger import get_main_logger
    from strkit.vcf_utils.header import build_vcf_header

    from .loci import load_loci
    from .output import (
        output_json_report_header,
        output_json_report_results,
        output_json_report_footer,
        output_tsv as output_tsv_fn,
        output_contig_vcf_lines,
    )
    from .utils import get_new_seed

    # ------------------------------------------------------------------------------------------------------------------

    logger = get_main_logger()

    # Start the call timer
    start_time = time.perf_counter()

    logger.info("strkit_rust_ext version: %s", importlib.metadata.version("strkit_rust_ext"))

    logger.info("Starting STR genotyping with parameters:")
    for k, v in params.to_dict().items():
        logger.info("    %s = %s", k.rjust(27), f'"{v}"' if isinstance(v, str) else v)

    # Seed the random number generator if a seed is provided, for replicability
    rng: NPRandomGenerator = np_default_rng(seed=params.seed)

    manager: SyncManager = mp.Manager()

    locus_queue = manager.Queue()
    # Loads actual STRkitLocus definitions into locus_queue, but returns # of loci + the contig set:
    num_loci, contig_set, loci_hash = load_loci(params, locus_queue, logger)

    # ---

    # If we're outputting JSON, write the header
    if json_path is not None:
        output_json_report_header(params, contig_set, json_path, indent_json, num_loci, loci_hash)

    # If we're outputting a VCF, open the file and write the header
    sample_id_str = params.sample_id or "sample"
    vf: VariantFile | None = None
    if vcf_path is not None:
        vh = build_vcf_header(
            "call",
            (sample_id_str,),
            Path(params.reference_file),
            partial_phasing=params.snv_vcf is not None or params.use_hp,
            num_loci=num_loci,
            loci_hash=loci_hash,
        )
        vf = VariantFile(vcf_path if vcf_path != "stdout" else "-", "w", header=vh)

    # ---

    is_single_processed = params.processes == 1

    if is_single_processed:
        from multiprocessing.dummy import Pool
        pool_class = Pool
    else:
        logger.info("Using %d workers", params.processes)
        from .non_daemonic_pool import NonDaemonicPool
        pool_class = NonDaemonicPool

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
        avg_read_depths_by_contig: dict[str, float] = {}
        total_read_depth: int = 0
        total_loci_with_reads: int = 0

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

            #  - calculate average read depth for contig for a little logging report at the end
            if results:
                rd: int = 0
                nwr: int = 0
                for r in results:
                    if "reads" not in r:
                        continue
                    rd += len(r["reads"])
                    nwr += 1
                avg_read_depths_by_contig[results[0]["contig"]] = rd / nwr
                total_read_depth += rd
                total_loci_with_reads += nwr
            else:
                avg_read_depths_by_contig[f"unknown{n_contigs_processed}"] = 0.0

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

    avg_read_depth = total_read_depth / total_loci_with_reads

    logger.info("    Average read depth (called loci only): %.1f", avg_read_depth)
    for rc, rr in avg_read_depths_by_contig.items():
        logger.info("    Average read depth (called loci only) for contig %s: %.1f", rc, rr)

    if json_path:
        output_json_report_footer(avg_read_depth, avg_read_depths_by_contig, time_taken, json_path, indent_json)
