import numpy as np
import pyfamsa
from scipy.stats import mode
from typing import Iterable

__all__ = [
    "consensus_seq",
]

GAP_ORD = ord("-")

aligner = pyfamsa.Aligner(guide_tree="upgma")
join_str = "".join


def _not_gap(x: int) -> bool:
    return x != GAP_ORD


def consensus_seq(seqs: Iterable[str]) -> str:
    msa = aligner.align(pyfamsa.Sequence(f"s{i}".encode("ascii"), seq.encode("ascii")) for i, seq in enumerate(seqs))
    msa_grid = np.array(list(map(lambda aseq: list(map(int, aseq.sequence)), msa)))
    return join_str(map(chr, filter(_not_gap, getattr(mode(msa_grid), "mode", ()))))
