import logging
from strkit_rust_ext import best_representative, consensus_seq as _consensus_seq
from typing import Iterable, Optional

from .types import ConsensusMethod

__all__ = [
    "best_representative",
    "consensus_seq",
]


def _run_best_representative(seqs: list[str], logger: logging.Logger) -> Optional[tuple[str, ConsensusMethod]]:
    res = best_representative(seqs)
    method: ConsensusMethod = "best_rep"
    if res is None:
        logger.debug(f"Got no best representative from sequences {seqs}")
        return None
    return res, method


def consensus_seq(
    seqs: Iterable[str], logger: logging.Logger, max_mdn_poa_length: int
) -> Optional[tuple[str, ConsensusMethod]]:
    # Return a stringified, gapless version of the column-wise mode for the MSA
    # If the consensus fails, try a best-representative strategy instead. If that fails, something's gone wrong...

    seqs_l = list(seqs)
    n_seqs = len(seqs_l)

    if n_seqs == 0:
        return None
    elif len(set(seqs_l)) == 1:
        return seqs_l[0], "single"

    n_blanks = seqs_l.count("")

    if n_blanks >= round(n_seqs / 2):
        # blanks make up majority, so blank is the consensus
        return "", "best_rep"
    elif n_blanks > 0:
        # blanks make up minority, so filter them out for consensus
        seqs_l = [s for s in seqs_l if s != ""]
        n_seqs -= n_blanks

    seqs_l.sort()
    if len(seqs_l[len(seqs_l) // 2]) > max_mdn_poa_length:
        return _run_best_representative(seqs_l, logger)

    res: Optional[str] = _consensus_seq(seqs_l)
    method: ConsensusMethod = "poa"

    if res is None:
        logger.error(f"Got no consensus sequence from sequences {seqs_l}; trying best representative strategy")
        return _run_best_representative(seqs_l, logger)

    return res, method
