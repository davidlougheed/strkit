import itertools

from typing import Generator, Iterable, Union

__all__ = [
    "get_aligned_pairs_from_cigar",
    "decode_cigar",
]


def get_aligned_pairs_from_cigar(
    cigar: Iterable[tuple[int, int]],
    query_start: int = 0,
    ref_start: int = 0,
) -> Generator[tuple[Union[int, None], Union[int, None]], None, None]:
    qi = itertools.count(start=query_start)
    di = itertools.count(start=ref_start)

    for c in cigar:
        op, cnt = c
        rc = range(cnt)

        # TODO: Probably a nicer way to do this:

        cq = ((next(qi), None) for _ in rc)
        cr = ((None, next(di)) for _ in rc)
        cb = ((next(qi), next(di)) for _ in rc)

        yield from {
            0: cb,  # M
            1: cq,  # I
            2: cr,  # D
            3: cr,  # N
            4: cq,  # S
            5: (),  # H
            6: (),  # P
            7: cb,  # =
            8: cb,  # X
        }[op]


def decode_cigar(encoded_cigar: list[int]) -> Generator[tuple[int, int], None, None]:
    for item in encoded_cigar:
        yield item & 15, item >> 4
