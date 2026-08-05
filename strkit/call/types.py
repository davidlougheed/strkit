from typing import Literal, NotRequired, Required, TypedDict

import numpy as np
from numpy.typing import NDArray
from strkit_rust_ext import CallData

__all__ = [
    "AssignMethod",
    "AssignMethodWithHP",
    "ConsensusMethod",
    # ---
    "ReadDict",
    "SnvBases",
    "CalledSNV",
    "LocusResult",
]

AssignMethod = Literal["dist", "snv", "snv+dist", "single"]
AssignMethodWithHP = AssignMethod | Literal["hp"]

ConsensusMethod = Literal["single", "poa", "best_rep"]


class ReadDict(TypedDict, total=False):
    s: Required[Literal["-", "+"]]  # DNA strand alignment
    cn: Required[int]  # Copy number
    w: Required[float]  # Weight
    sc: Required[float | None]  # Adjusted read model align score (None if TR is missing)

    # ------------------------------------------------------------------------------------------------------------------

    # Whether the read was realigned by hand using a local alignment algorithm.
    realn: bool

    # Whether the read appears to be chimeric within the locus region,
    # i.e. aligned twice with different soft-clipping.
    chimeric_in_region: bool

    p: int  # Peak (allele)

    kmers: dict[str, int]  # Dictionary of {kmer: count}

    # Only added if HP tags from a haplotagged alignment file are being incorporated:
    hp: int
    ps: int

    # Only added if SNVs are being incorporated:
    #  - After including only useful SNVs, this contains a tuple of bases for just those + corresponding qualities
    snvu: tuple[tuple[str, int], ...]

    # BEGIN: only added if consensus sequences are being calculated:
    sl: int
    #  - the below two are only output if --json-read-seq is passed
    start_anchor_seq: str  # Left anchor for calculated allele sequence (usually 1 base)
    seq: str  # Tandem repeat sequence
    # END: only added if consensus sequences are being calculated

    # BEGIN: only added if methylation is being incorporated:
    m: float  # methylation proportion
    mc: int  # methylation site count
    # END: only added if methylation is being incorporated


SnvBases = tuple[tuple[str, int], ...]


class CalledSNV(TypedDict):
    id: str
    pos: int
    call: tuple[str, ...]
    rcs: list[int]
    ref: NotRequired[str]


class PeakData(TypedDict):
    means: Required[NDArray[np.float64]]
    weights: Required[NDArray[np.float64]]
    stdevs: Required[NDArray[np.float64]]
    modal_int: Required[int]
    n_reads: Required[list[int]]

    # --------------------------------------------------

    kmers: dict[str, int]
    seqs: list[tuple[str, ConsensusMethod]]
    start_anchor_seqs: list[tuple[str, ConsensusMethod]]
    am: list[float]
    amc: list[float]


class LocusResult(TypedDict, total=False):
    # BEGIN from STRkitLocus
    locus_index: Required[int]
    locus_id: Required[str]
    contig: Required[str]
    start: Required[int]
    end: Required[int]

    motif: Required[str]

    annotations: Required[list[str]]
    # END from STRkitLocus

    assign_method: Required[AssignMethodWithHP | None]
    call_data: Required[CallData | None]
    call: Required[list[int] | None]
    call_95_cis: Required[list[list[int]] | None]
    call_99_cis: Required[list[list[int]] | None]

    # -------------------------------------------------------
    start_adj: int
    end_adj: int

    ref_cn: int

    ps: int | None
    peaks: PeakData | None
    read_peaks_called: bool
    time: float

    # if we're in consensus mode: ---
    ref_start_anchor: str
    ref_seq: str
    # ---

    reads: dict[str, ReadDict]
    snvs: list[CalledSNV]

    # Mean model (candidate TR sequence) alignment score across reads.
    mean_model_align_score: float | None
