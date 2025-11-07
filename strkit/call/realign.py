from __future__ import annotations

import os
import time

from multiprocessing import Queue as MpQueue, Process
from parasail import sg_dx_trace_scan_16
from queue import Empty as QueueEmpty
from strkit_rust_ext import calculate_seq_with_wildcards
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from logging import Logger
    from strkit_rust_ext import STRkitAlignedCoords, STRkitAlignedSegment, STRkitLocusWithRefData
    from .params import CallParams

from .align_matrix import match_score, dna_matrix
from .cigar import decode_cigar_np, get_aligned_pair_matches

__all__ = [
    "realign_read",
    "perform_realign",
]


min_realign_score_ratio: float = 0.95  # TODO: parametrize
realign_indel_open_penalty: int = 7  # TODO: parametrize
max_ref_len_for_same_proc: int = 2000  # TODO: parametrize
max_read_len_for_same_proc: int = 25000  # TODO: parametrize


def realign_read(
    ref_seq: str,
    query_seq: str,
    left_flank_coord: int,
    flank_size: int,
    q,  # mp.Queue | None
    read_log_str: str,
    log_level: int,
) -> STRkitAlignedCoords | None:
    # Have to re-attach logger in separate process I guess

    def ret_q(v: STRkitAlignedCoords | None) -> STRkitAlignedCoords | None:
        if q:
            q.put(v)
            q.close()
        return v

    from strkit.logger import create_process_logger
    lg = create_process_logger(os.getpid(), log_level)

    # flipped: 'ref sequence' as query here, since it should in general be shorter (!)
    pr = sg_dx_trace_scan_16(
        # fetch an extra base for the right flank coordinate check later (needs to be >= the exclusive coord)
        ref_seq, query_seq, realign_indel_open_penalty, 0, dna_matrix)

    if pr.score < (th := min_realign_score_ratio * (flank_size * 2 * match_score - realign_indel_open_penalty)):
        lg.debug("Realignment for %s scored below threshold (%d < %.2f)", read_log_str, pr.score, th)
        return ret_q(None)

    lg.debug("Realigned %s: scored %d; Flipped CIGAR: %s", read_log_str, pr.score, pr.cigar.decode.decode("ascii"))

    res: STRkitAlignedCoords = get_aligned_pair_matches(decode_cigar_np(pr.cigar.seq), left_flank_coord, 0)
    return ret_q(res)


def perform_realign(
    locus_with_ref_data: STRkitLocusWithRefData,
    segment: STRkitAlignedSegment,
    # ---
    params: CallParams,
    # ---
    logger_: Logger,
    locus_log_str: str,
) -> STRkitAlignedCoords | None:
    rn = segment.name

    # TODO: add a segment method for doing this on a specified slice of query sequence or something...
    qs_wc = calculate_seq_with_wildcards(segment.query_sequence, segment.query_qualities, 3)

    ref_total_seq = locus_with_ref_data.ref_total_seq
    ref_seq_len = len(ref_total_seq)
    qs_len = len(qs_wc)
    left_flank_coord = locus_with_ref_data.locus_def.left_flank_coord

    read_log_str = f"{rn} in locus {locus_with_ref_data.locus_def.t_idx}"

    if ref_seq_len <= max_ref_len_for_same_proc and qs_len <= max_read_len_for_same_proc:
        # Don't start process for short realigns, since then process startup dominates the total time taken
        # TODO: more robust solution; realign worker somehow? How to do timeout?
        return realign_read(
            ref_total_seq, qs_wc, left_flank_coord, params.flank_size, None, locus_log_str, params.log_level
        )

    t = time.time()

    q: MpQueue = MpQueue()
    proc = Process(target=realign_read, daemon=False, kwargs=dict(
        # fetch an extra base for the right flank coordinate check later (needs to be >= the exclusive coord)
        ref_seq=ref_total_seq,  # TODO: with the plus 1, really?
        query_seq=qs_wc,
        left_flank_coord=left_flank_coord,
        flank_size=params.flank_size,
        q=q,
        read_log_str=read_log_str,
        log_level=params.log_level,
    ))
    proc.start()

    pairs_new: STRkitAlignedCoords | None = None
    try:
        pairs_new = q.get(timeout=params.realign_timeout)
        proc.join()
    except QueueEmpty:
        logger_.warning(
            "%s - experienced timeout while re-aligning read %s. Reverting to initial alignment.", locus_log_str, rn)
        proc.terminate()
        time.sleep(0.1)  # wait a little for the process to terminate
    finally:
        wait_count: int = 0
        while proc.is_alive():
            logger_.warning("%s - realign job has still not exited. Waiting 0.5 seconds...", locus_log_str)
            time.sleep(0.5)
            wait_count += 1
            if wait_count > 30:
                logger_.fatal("%s - realign job never exited. Terminating...", locus_log_str)
                exit(1)
        proc.close()

    logger_.debug(
        "%s - %s: long realign job completed in %.4fs (ref_seq_len=%d, qs_len=%d)",
        locus_log_str, rn, time.time() - t, ref_seq_len, qs_len)

    return pairs_new
