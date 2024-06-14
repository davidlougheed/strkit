# import pysam
import numpy as np
from typing import Literal, Optional, TypedDict, Union
from numpy.typing import NDArray


__all__ = [
    "VCFContigFormat",
    "AssignMethod",
    "AssignMethodWithHP",
    "ConsensusMethod",
    # ---
    "ReadDict",
    "ReadDictExtra",
    "CalledSNV",
    "LocusResult",
]

# TODO: py3.10: new Required[] TypedDict structuring


VCFContigFormat = Literal["chr", "num", "acc", ""]

AssignMethod = Literal["dist", "snv", "snv+dist", "single"]
AssignMethodWithHP = Union[AssignMethod, Literal["hp"]]

ConsensusMethod = Literal["single", "poa", "best_rep"]


class _ReadDictBase(TypedDict):
    s: Literal["-", "+"]  # DNA strand alignment
    cn: Union[int, float]  # Copy number
    w: float  # Weight


class ReadDict(_ReadDictBase, total=False):
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


class ReadDictExtra(TypedDict, total=False):
    _ref_start: int  # Read start in ref coordinates
    _ref_end: int  # Read end in ref coordinates

    _tr_seq: str  # Tandem repeat sequence... only added if consensus is being calculated

    # Below are only added if SNVs are being incorporated:

    _qs: str  # Query (read) sequence
    _fqqs: NDArray[np.uint8]  # Query (read) base qualities

    sig_clip_left: bool  # Significant amounts of clipping (5' of read)
    sig_clip_right: bool  # Significant amounts of clipping (3' of read)

    snv: dict[int, str]  # Intermediate result: dictionary of a bunch of SNVs for this read {position: base}
    # Intermediate result: tuple of bases/qualities for the set of SNVs across all reads
    snv_bases: tuple[tuple[str, int], ...]


class _CalledSNVBase(TypedDict):
    id: str
    pos: int
    call: tuple[str, ...]
    rcs: list[int]


class CalledSNV(_CalledSNVBase, total=False):
    ref: str


class BasePeakData(TypedDict):
    means: NDArray[np.float32]
    weights: NDArray[np.float32]
    stdevs: NDArray[np.float32]
    modal_int: int
    n_reads: list[int]


class PeakData(BasePeakData):
    kmers: dict[str, int]
    seqs: list[tuple[str, ConsensusMethod]]  # really "list-tyup


class BaseLocusResult(TypedDict):
    locus_index: int
    contig: str
    start: int
    end: int

    motif: str

    assign_method: Optional[AssignMethodWithHP]
    call: Optional[list[int]]
    call_95_cis: Optional[list[list[int]]]
    call_99_cis: Optional[list[list[int]]]


class LocusResult(BaseLocusResult, total=False):
    start_adj: int
    end_adj: int

    ref_cn: int

    ps: Optional[int]
    peaks: Optional[PeakData]
    read_peaks_called: bool
    time: float

    # if we're in consensus mode: ---
    ref_start_anchor: str
    ref_seq: str
    # ---

    reads: dict[str, ReadDict]
    snvs: list[CalledSNV]
