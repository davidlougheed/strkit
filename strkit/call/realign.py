import logging
import multiprocessing as mp
import numpy as np
import os
import parasail
import queue
import time

from numpy.typing import NDArray
from strkit_rust_ext import STRkitAlignedCoords

from .align_matrix import match_score, dna_matrix
from .cigar import decode_cigar_np, get_aligned_pair_matches
from .params import CallParams
from .utils import calculate_seq_with_wildcards

__all__ = [
    "MatchedCoordPairListOrNone",
    "realign_read",
    "perform_realign",
]


min_realign_score_ratio: float = 0.95  # TODO: parametrize
realign_indel_open_penalty: int = 7  # TODO: parametrize
max_ref_len_for_same_proc: int = 2000  # TODO: parametrize
max_read_len_for_same_proc: int = 25000  # TODO: parametrize


MatchedCoordPairList = tuple[NDArray[np.uint64], NDArray[np.uint64]]
MatchedCoordPairListOrNone = MatchedCoordPairList | None


def realign_read(
    ref_seq: str,
    query_seq: str,
    left_flank_coord: int,
    flank_size: int,
    q,  # mp.Queue | None
    read_log_str: str,
    log_level: int = logging.WARNING,
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
    pr = parasail.sg_dx_trace_scan_16(
        # fetch an extra base for the right flank coordinate check later (needs to be >= the exclusive coord)
        ref_seq, query_seq, realign_indel_open_penalty, 0, dna_matrix)

    if pr.score < (th := min_realign_score_ratio * (flank_size * 2 * match_score - realign_indel_open_penalty)):
        lg.debug(f"Realignment for {read_log_str} scored below threshold ({pr.score} < {th:.2f})")
        return ret_q(None)

    lg.debug(f"Realigned {read_log_str}: scored {pr.score}; Flipped CIGAR: {pr.cigar.decode.decode('ascii')}")

    res: STRkitAlignedCoords = get_aligned_pair_matches(decode_cigar_np(pr.cigar.seq), left_flank_coord, 0)
    return ret_q(res)


def perform_realign(
    t_idx: int,
    left_flank_coord: int,
    ref_total_seq: str,
    rn: str,
    qs: str,
    fqqs: NDArray[np.uint8],
    # ---
    params: CallParams,
    realign_timeout: int,
    # ---
    logger_: logging.Logger,
    locus_log_str: str,
) -> STRkitAlignedCoords | None:
    qs_wc = calculate_seq_with_wildcards(qs, fqqs)

    ref_seq_len = len(ref_total_seq)
    qs_len = len(qs_wc)

    read_log_str = f"{rn} in locus {t_idx}"

    if ref_seq_len <= max_ref_len_for_same_proc and qs_len <= max_read_len_for_same_proc:
        # Don't start process for short realigns, since then process startup dominates the total time taken
        # TODO: more robust solution; realign worker somehow? How to do timeout?
        return realign_read(
            ref_total_seq, qs_wc, left_flank_coord, params.flank_size, None, locus_log_str, params.log_level
        )

    t = time.time()

    q: mp.Queue = mp.Queue()
    proc = mp.Process(target=realign_read, daemon=False, kwargs=dict(
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
        pairs_new = q.get(timeout=realign_timeout)
        proc.join()
    except queue.Empty:
        logger_.warning(
            f"{locus_log_str} - experienced timeout while re-aligning read {rn}. Reverting to initial "
            f"alignment.")
        proc.terminate()
        time.sleep(0.1)  # wait a little for the process to terminate
    finally:
        wait_count: int = 0
        while proc.is_alive():
            logger_.warning(f"{locus_log_str} - realign job has still not exited. Waiting 0.5 seconds...")
            time.sleep(0.5)
            wait_count += 1
            if wait_count > 30:
                logger_.fatal(f"{locus_log_str} - realign job never exited. Terminating...")
                exit(1)
        proc.close()

    logger_.debug(
        f"{locus_log_str} - {rn}: long realign job completed in {time.time() - t:.4f}s ({ref_seq_len=}, {qs_len=})")

    return pairs_new
