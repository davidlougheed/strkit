import itertools

from typing import Callable, Iterable, Union

from ..logger import logger

__all__ = [
    "CoordPair",
    "get_aligned_pairs_from_cigar",
    "decode_cigar",
    "get_aligned_pair_matches",
]


NONE_GENERATOR = itertools.repeat(None)


def _cq(count: int, qi: Iterable, _di: Iterable):
    return zip(itertools.islice(qi, count), NONE_GENERATOR)


def _cqc(count: int, qi: Iterable, _di: Iterable):
    tuple(itertools.islice(qi, count))  # consume w/o returning
    return ()


def _cr(count: int, _qi: Iterable, di: Iterable):
    return zip(NONE_GENERATOR, itertools.islice(di, count))


def _crc(count: int, _qi: Iterable, di: Iterable):
    tuple(itertools.islice(di, count))  # consume w/o returning
    return ()


def _cb(count: int, qi: Iterable, di: Iterable):
    # only need to .islice(...) one, since zip will auto-terminate when the shortest one is consumed:
    return zip(itertools.islice(qi, count), di)


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

CIGAR_OPS_MATCHES_ONLY = (
    _cb,   # MATCH        | 0 | M
    _cqc,  # INSERTION    | 1 | I
    _crc,  # DELETION     | 2 | D
    _crc,  # SKIPPED      | 3 | N
    _cqc,  # SOFT_CLIPPED | 4 | S
    _no,   # HARD_CLIPPED | 5 | H
    _no,   # PADDING      | 6 | P
    _cb,   # SEQ_MATCH    | 7 | =
    _cb,   # SEQ_MISMATCH | 8 | X
)


CoordPair = tuple[Union[int, None], Union[int, None]]


def get_aligned_pairs_from_cigar(
    cigar: Iterable[tuple[int, int]],
    query_start: int = 0,
    ref_start: int = 0,
    matches_only: bool = False,
) -> Iterable[CoordPair]:
    """
    Given an iterable of CIGAR operations (op, count), yield aligned pairs of (query, ref).
    :param cigar: Iterable of CIGAR operations
    :param query_start: The starting query coordinate for the alignment.
    :param ref_start: The starting reference coordinate for the alignment.
    :param matches_only: Whether to only return matches, i.e., no Nones.
    :return: Generator of aligned pairs of (query_coord, ref_coord) (None if a gap is aligned.)
    """

    qi = itertools.count(start=query_start)
    di = itertools.count(start=ref_start)
    ops = CIGAR_OPS_MATCHES_ONLY if matches_only else CIGAR_OPS

    return itertools.chain.from_iterable(map(lambda c: ops[c[0]](c[1], qi, di), cigar))


def _decode_cigar_item(item: int) -> tuple[int, int]:
    return item & 15, item >> 4


def decode_cigar(encoded_cigar: list[int]) -> Iterable[tuple[int, int]]:
    return map(_decode_cigar_item, encoded_cigar)


def _get_aligned_pair_matches_py(
    cigar_tuples: list[tuple[int, int]], query_start: int, ref_start: int
) -> tuple[list[int], list[int]]:
    # noinspection PyTypeChecker
    return tuple(map(list, zip(*get_aligned_pairs_from_cigar(cigar_tuples, query_start, ref_start, matches_only=True))))


get_aligned_pair_matches: Callable[[list[tuple[int, int]], int, int], tuple[list[int], list[int]]]
try:
    from strkit_rust_ext import get_aligned_pair_matches
    logger.debug("Found STRkit Rust component, importing get_aligned_pair_matches")
except ImportError:
    get_aligned_pair_matches = _get_aligned_pair_matches_py
