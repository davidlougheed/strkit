import logging
from random import choice
from strkit_rust_ext import best_representatives as _best_representatives, consensus_seq as _consensus_seq
from typing import Iterable, Optional, Sequence

__all__ = [
    "best_representative",
    "consensus_seq",
]


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


def consensus_seq(seqs: Iterable[str], logger: logging.Logger) -> Optional[str]:
    # Return a stringified, gapless version of the column-wise mode for the MSA
    # If the consensus fails, try a best-representative strategy instead. If that fails, something's gone wrong...

    seqs_l = list(seqs)
    if len(set(seqs_l)) == 1:
        return seqs_l[0]

    seqs_l.sort()
    res: Optional[str] = _consensus_seq(seqs_l)
    if res is None:
        logger.error(f"Got no consensus sequence from sequences {seqs_l}; trying best representative strategy")
        res = best_representative(seqs_l)
        if res is None:
            logger.debug(f"Got no best representative from sequences {seqs_l}")
    return res
