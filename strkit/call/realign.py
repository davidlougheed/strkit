import logging
import multiprocessing as mp
import numpy as np
import os
import parasail

from numpy.typing import NDArray
from typing import Optional

from .align_matrix import match_score, dna_matrix
from .cigar import decode_cigar, get_aligned_pair_matches

min_realign_score_ratio: float = 0.95  # TODO: parametrize
realign_indel_open_penalty: int = 7  # TODO: parametrize


MatchedCoordPairList = tuple[NDArray[np.uint64], NDArray[np.uint64]]
MatchedCoordPairListOrNone = Optional[MatchedCoordPairList]


def realign_read(
    ref_seq: str,
    query_seq: str,
    left_flank_coord: int,
    flank_size: int,
    rn: str,
    t_idx: int,
    always_realign: bool,
    q: Optional[mp.Queue] = None,  # TODO: why was this optional, again...
    log_level: int = logging.WARNING,
) -> MatchedCoordPairListOrNone:
    # Have to re-attach logger in separate process I guess

    def ret_q(v: MatchedCoordPairListOrNone) -> MatchedCoordPairListOrNone:
        if q:
            q.put(v)
            q.close()
        return v

    from strkit.logger import create_process_logger
    lg = create_process_logger(os.getpid(), log_level)

    # flipped: 'ref sequence' as query here, since it should in general be shorter (!)
    pr = parasail.sg_dx_trace_scan_sat(
        # fetch an extra base for the right flank coordinate check later (needs to be >= the exclusive coord)
        ref_seq, query_seq, realign_indel_open_penalty, 0, dna_matrix)

    if pr.score < (th := min_realign_score_ratio * (flank_size * 2 * match_score - realign_indel_open_penalty)):
        lg.debug(f"Realignment for {rn} scored below threshold ({pr.score} < {th:.2f})")
        return ret_q(None)

    lg.debug(
        f"Realigned {rn} in locus {t_idx}{' (due to soft clipping)' if not always_realign else ''}: scored {pr.score}; "
        f"Flipped CIGAR: {pr.cigar.decode.decode('ascii')}")

    matches = get_aligned_pair_matches(list(decode_cigar(pr.cigar.seq)), left_flank_coord, 0)
    res: MatchedCoordPairList = (np.array(matches[1], dtype=np.uint64), np.array(matches[0], dtype=np.uint64))
    return ret_q(res)
