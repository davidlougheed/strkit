import logging
import multiprocessing as mp
import os
import parasail

from typing import Optional

from .align_matrix import match_score, dna_matrix
from .cigar import get_aligned_pairs_from_cigar, decode_cigar

min_realign_score_ratio = 0.95  # TODO: parametrize
realign_indel_open_penalty = 7  # TODO: parametrize


def realign_read(
    ref_seq: str,
    query_seq: str,
    left_flank_coord: int,
    flank_size: int,
    rn: str,
    t_idx: int,
    always_realign: bool,
    q: Optional[mp.Queue] = None,
    log_level: int = logging.WARNING,
) -> Optional[list[tuple[Optional[int], Optional[int]]]]:
    # Have to re-attach logger in separate process I guess

    from strkit.logger import create_process_logger
    lg = create_process_logger(os.getpid(), log_level)

    # flipped: 'ref sequence' as query here, since it should in general be shorter
    pr = parasail.sg_dx_trace_scan_sat(
        # fetch an extra base for the right flank coordinate check later (needs to be >= the exclusive coord)
        ref_seq, query_seq, realign_indel_open_penalty, 0, dna_matrix)

    if pr.score < (th := min_realign_score_ratio * (flank_size * 2 * match_score - realign_indel_open_penalty)):
        lg.debug(f"Realignment for {rn} scored below threshold ({pr.score} < {th:.2f})")
        q.put(None)
        return None

    lg.debug(
        f"Realigned {rn} in locus {t_idx}{' (due to soft clipping)' if not always_realign else ''}: scored {pr.score}; "
        f"CIGAR: {pr.cigar.decode.decode('ascii')}")
    # noinspection PyTypeChecker
    res: list[tuple[Optional[int], Optional[int]]] = [
        # reverse to get (ref, query) instead of (query, ref) due to the flip
        tuple(reversed(p))
        # query_start instead of ref_start, due to the flip
        for p in get_aligned_pairs_from_cigar(decode_cigar(pr.cigar.seq), query_start=left_flank_coord)
        if p[0] is not None and p[1] is not None
    ]
    if q:
        q.put(res)
    return res
