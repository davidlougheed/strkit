import bisect
from pathlib import Path
from typing import Iterable, Optional, Union


def _line_filter_fn(s: str) -> bool:
    """
    Filter function to skip blank lines and comments
    :param s: line of a file
    :return: whether the line is not blank and is not a comment
    """
    return s and not s.startswith("#")


# key: contig, value: dict of (key: coordinate interval, value: list of extra values)
LociDictOfDict = dict[str, dict[tuple[int, int], list[str]]]

# key: contig, value: list of coordinate intervals
LociDictOfList = dict[str, list[tuple[int, int]]]


def build_loci_dict_of_dict_from_file(loci_path: Optional[Union[str, Path]]) -> LociDictOfDict:
    if not loci_path:
        return {}

    res: LociDictOfDict = {}

    with open(loci_path, "r") as lf:
        for line in filter(_line_filter_fn, map(str.strip, lf)):
            ls = line.split("\t")

            contig, ss, es = ls[:3]

            if contig not in res:
                res[contig] = {}

            res[contig][int(ss), int(es)] = ls[3:]

    return res


def build_loci_dict_of_list_from_file(loci_path: Optional[Union[str, Path]]) -> LociDictOfList:
    if not loci_path:
        return {}

    res: dict[str, list[tuple[int, int]]] = {}

    with open(loci_path, "r") as lf:
        for line in filter(_line_filter_fn, map(str.strip, lf)):
            ls = line.split("\t")
            contig = ls[0]
            if contig not in res:
                res[contig] = []
            res[contig].append((int(ls[1]), int(ls[2])))

    return res


def overlapping_loci_dict_of_dict(
    contig: str, start: int, end: int, d: LociDictOfDict
) -> list[tuple[int, int, list[str]]]:
    if contig not in d:
        return []

    c_dict = d[contig]
    c_keys = tuple(c_dict.keys())
    c_lhs = tuple(map(lambda k: k[0], c_keys))
    i = bisect.bisect_left(c_lhs, end)  # use _left since end is exclusive

    possible_overlaps = c_keys[:i]

    acc: list[tuple[int, int, list[str]]] = []

    for ov in possible_overlaps:
        if start < ov[1]:  # end is exclusive
            acc.append((ov[0], ov[1], c_dict[ov]))

    return acc


def overlapping_loci_dict_of_list(contig: str, start: int, end: int, d: LociDictOfList) -> Iterable[tuple[int, int]]:
    if contig not in d:
        yield from ()
        return

    c_ints = d[contig]
    c_lhs = tuple(map(lambda k: k[0], c_ints))
    i = bisect.bisect_left(c_lhs, end)  # use _left since end is exclusive

    yield from ((ov[0], ov[1]) for ov in c_ints[:i] if start < ov[1])
