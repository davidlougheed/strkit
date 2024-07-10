import pathlib
import pytest

from strkit.mi.intervals import (
    build_loci_dict_of_dict_from_file,
    overlapping_loci_dict_of_dict,
    build_loci_dict_of_list_from_file,
    overlapping_loci_dict_of_list,
)

TEST_LOCI = pathlib.Path(__file__).parent / "data" / "test_loci.bed"

BED_CASES = [
    ("chr1", 50, 70, 0),
    ("chr1", 205, 210, 1),
    ("chr1", 50, 1000, 3),
    ("chr1", 320, 500, 2),
    ("chr1", 400, 450, 1),
    ("chr1", 1000, 1001, 0),
    ("chr2", 100, 101, 1),
    ("chr2", 100, 200, 1),
    ("asdf", 50, 1000, 0),
]


@pytest.mark.parametrize("contig,start,end,nr", BED_CASES)
def test_loci_dict_of_dict(contig: str, start: int, end: int, nr: int):
    d = build_loci_dict_of_dict_from_file(TEST_LOCI)
    assert len(overlapping_loci_dict_of_dict(contig, start, end, d)) == nr


@pytest.mark.parametrize("contig,start,end,nr", BED_CASES)
def test_loci_dict_of_list(contig: str, start: int, end: int, nr: int):
    d = build_loci_dict_of_list_from_file(TEST_LOCI)
    assert len(tuple(overlapping_loci_dict_of_list(contig, start, end, d))) == nr
