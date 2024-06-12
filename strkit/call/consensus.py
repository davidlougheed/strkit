import logging
from random import choice
from strkit_rust_ext import best_representatives as _best_representatives, consensus_seq as _consensus_seq
from typing import Iterable, Literal, Optional, Sequence

__all__ = [
    "ConsensusMethod",
    "best_representative",
    "consensus_seq",
]


ConsensusMethod = Literal["single", "poa", "best_rep"]


def best_representative(seqs: Sequence[str]) -> Optional[str]:
    """
    Slightly different from a true consensus - returns the string with the minimum Levenshtein distance to all other
    strings for a particular allele. This roughly approximates a true consensus when |seqs| is large. If more than one
    best representative exist, one is chosen at random. If |best| == |seqs| or |best| == 0, None is returned since there
    is effectively no true consensus.
    :param seqs: An iterable of sequences to find the best representative of.
    :return: One of the best representative sequences from the passed sequences.
    """
    res = _best_representatives(seqs)
    if len(res) == 1:
        return res.pop()
    elif len(res) == 0:
        return None
    else:
        return choice(tuple(res))


def _run_best_representative(seqs: list[str], logger: logging.Logger) -> Optional[tuple[str, ConsensusMethod]]:
    res = best_representative(seqs)
    method: ConsensusMethod = "best_rep"
    if res is None:
        logger.debug(f"Got no best representative from sequences {seqs}")
        return None
    return res, method


def consensus_seq(
    seqs: Iterable[str], logger: logging.Logger, max_poa_length: int
) -> Optional[tuple[str, ConsensusMethod]]:
    # Return a stringified, gapless version of the column-wise mode for the MSA
    # If the consensus fails, try a best-representative strategy instead. If that fails, something's gone wrong...

    seqs_l = list(seqs)

    if len(seqs_l) == 0:
        return None
    elif len(set(seqs_l)) == 1:
        return seqs_l[0], "single"

    seqs_l.sort()
    if any(s > max_poa_length for s in map(len, seqs_l)):
        return _run_best_representative(seqs_l, logger)

    res: Optional[str] = _consensus_seq(seqs_l)
    method: ConsensusMethod = "poa"

    if res is None:
        logger.error(f"Got no consensus sequence from sequences {seqs_l}; trying best representative strategy")
        return _run_best_representative(seqs_l, logger)

    return res, method
