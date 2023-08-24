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

    for c in cigar:
        op, cnt = c

        cq = zip(itertools.islice(qi, cnt), NONE_GENERATOR)
        cr = zip(NONE_GENERATOR, itertools.islice(di, cnt))
        cb = zip(itertools.islice(qi, cnt), itertools.islice(di, cnt))

        yield from {
            CIGAR_OP_MATCH: cb,
            CIGAR_OP_INSERTION: cq,
            CIGAR_OP_DELETION: cr,
            CIGAR_OP_SKIPPED: cr,
            CIGAR_OP_SOFT_CLIPPED: cq,
            CIGAR_OP_HARD_CLIPPED: (),
            CIGAR_OP_PADDING: (),
            CIGAR_OP_SEQ_MATCH: cb,
            CIGAR_OP_SEQ_MISMATCH: cb,
        }[op]


def decode_cigar(encoded_cigar: list[int]) -> Generator[tuple[int, int], None, None]:
    for item in encoded_cigar:
        yield item & 15, item >> 4
