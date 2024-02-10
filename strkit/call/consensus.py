from random import choice
from strkit_rust_ext import best_representatives, consensus_seq as _consensus_seq
from typing import Iterable, Optional

__all__ = [
    "best_representative",
    "consensus_seq",
]


def best_representative(seqs: Iterable[str]) -> Optional[str]:
    """
    Slightly different from a true consensus - returns the string with the minimum Levenshtein distance to all other
    strings for a particular allele. This roughly approximates a true consensus when |seqs| is large. If more than one
    best representative exist, one is chosen at random. If |best| == |seqs| or |best| == 0, None is returned since there
    is effectively no true consensus.
    :param seqs: An iterable of sequences to find the best representative of.
    :return: One of the best representative sequences from the passed sequences.
    """
    seqs_t = tuple(seqs)
    res = best_representatives(seqs_t)
    if len(res) == 1:
        return res.pop()
    elif len(res) == 0 or len(res) == len(seqs_t):
        return None
    else:
        res_t = tuple(res)
        print("\n".join(res_t))
        return choice(res_t)


def consensus_seq(seqs: Iterable[str]) -> Optional[str]:
    # Return a stringified, gapless version of the column-wise mode for the MSA
    # - Filter out blanks and if the consensus fails, try eliminating the outlier VERY naively
    #   via just comparing sorted values
    seqs_t = tuple(sorted(filter(lambda s: s != "", seqs)))
    res = _consensus_seq(seqs_t)
    if res is None:
        print(seqs_t)
    return res
