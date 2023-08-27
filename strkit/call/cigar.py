import itertools

from typing import Generator, Iterable, Union

__all__ = [
    "CIGAR_OP_MATCH",
    "CIGAR_OP_INSERTION",
    "CIGAR_OP_DELETION",
    "CIGAR_OP_SKIPPED",
    "CIGAR_OP_SOFT_CLIPPED",
    "CIGAR_OP_HARD_CLIPPED",
    "CIGAR_OP_PADDING",
    "CIGAR_OP_SEQ_MATCH",
    "CIGAR_OP_SEQ_MISMATCH",

    "get_aligned_pairs_from_cigar",
    "decode_cigar",
]


CIGAR_OP_MATCH = 0  # M
CIGAR_OP_INSERTION = 1  # I
CIGAR_OP_DELETION = 2  # D
CIGAR_OP_SKIPPED = 3  # N
CIGAR_OP_SOFT_CLIPPED = 4  # S
CIGAR_OP_HARD_CLIPPED = 5  # H
CIGAR_OP_PADDING = 6  # P
CIGAR_OP_SEQ_MATCH = 7  # =
CIGAR_OP_SEQ_MISMATCH = 8  # X

NONE_GENERATOR = itertools.repeat(None)


def _cq(count: int, qi: Iterable, _di: Iterable):
    return zip(itertools.islice(qi, count), NONE_GENERATOR)


def _cr(count: int, _qi: Iterable, di: Iterable):
    return zip(NONE_GENERATOR, itertools.islice(di, count))


def _cb(count: int, qi: Iterable, di: Iterable):
    return zip(itertools.islice(qi, count), itertools.islice(di, count))


def _no(_count: int, _qi: Iterable, _di: Iterable):
    return ()


CIGAR_OPS = (
    _cb,  # MATCH        | 0 | M
    _cq,  # INSERTION    | 1 | I
    _cr,  # DELETION     | 2 | D
    _cr,  # SKIPPED      | 3 | N
    _cq,  # SOFT_CLIPPED | 4 | S
    _no,  # HARD_CLIPPED | 5 | H
    _no,  # PADDING      | 6 | P
    _cb,  # SEQ_MATCH    | 7 | =
    _cb,  # SEQ_MISMATCH | 8 | X
)


def get_aligned_pairs_from_cigar(
    cigar: Iterable[tuple[int, int]],
    query_start: int = 0,
    ref_start: int = 0,
) -> Generator[tuple[Union[int, None], Union[int, None]], None, None]:
    """
    Given an iterable of CIGAR operations (op, count), yield aligned pairs of (query, ref).
    :param cigar: Iterable of CIGAR operations
    :param query_start: The starting query coordinate for the alignment.
    :param ref_start: The starting reference coordinate for the alignment.
    :return: Generator of aligned pairs of (query_coord, ref_coord) (None if a gap is aligned.)
    """

    qi = itertools.count(start=query_start)
    di = itertools.count(start=ref_start)

    for op, cnt in cigar:
        yield from CIGAR_OPS[op](cnt, qi, di)


def decode_cigar(encoded_cigar: list[int]) -> Generator[tuple[int, int], None, None]:
    for item in encoded_cigar:
        yield item & 15, item >> 4
