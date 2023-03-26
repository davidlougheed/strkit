from typing import Literal, TypedDict, Union


__all__ = [
    "ReadDict",
]

# TODO: py3.10: new Required[] TypedDict structuring


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

    snv: dict[int, str]  # Intermediate result: dictionary of a bunch of SNVs for this read {position: base}
    snv_bases: tuple[str, ...]  # Intermediate result: tuple of bases for the set of SNVs across all reads
    snvu: tuple[str, ...]  # After including only useful SNVs, this contains a tuple of bases for just those
