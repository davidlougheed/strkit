import numpy as np
import pyfamsa
from scipy.stats import mode
from typing import Iterable
from .utils import cat_strs

__all__ = [
    "consensus_seq",
]

GAP_ORD = ord("-")

aligner = pyfamsa.Aligner(threads=1, guide_tree="upgma", keep_duplicates=True)


def _not_gap(x: int) -> bool:
    return x != GAP_ORD


def consensus_seq(seqs: Iterable[str]) -> str:
    seqs_t = tuple(seqs)

    if seqs_t.count("") > len(seqs_t) // 2:
        return ""  # FAMSA crashes in various ways with blanks, so if we have a majority blank just return ""

    # Align all the allele tandem repeat sequences together into a multiple sequence alignment
    #  - Skip blanks, since FAMSA crashes in various ways with blank strings  TODO
    msa = aligner.align(
        pyfamsa.Sequence(f"s{i}".encode("ascii"), seq.encode("ascii"))
        for i, seq in enumerate(seqs_t) if seq)

    # Turn the MSA into a numpy matrix as setup for calculating the modal value at each position
    msa_grid = np.array(list(map(lambda aseq: list(map(int, aseq.sequence)), msa)))

    # Return a stringified, gapless version of the column-wise mode for the MSA
    return cat_strs(map(chr, filter(_not_gap, getattr(mode(msa_grid), "mode", ()))))
