from strkit.call.utils import find_pair_by_ref_pos, normalize_contig

#         A  A       T  T       C  G       C  C       C  C       A  A       A  A       A  C
PAIRS = [(0, 1000), (1, 1001), (2, 1003), (3, 1004), (4, 1005), (5, 1006), (6, 1008), (7, 1009)]
SNVS = ((1003, "C"), (1009, "A"))
PAIRS_Q = list(p[0] for p in PAIRS)
PAIRS_R = list(p[1] for p in PAIRS)


def test_find_pair_by_ref_pos():
    assert find_pair_by_ref_pos(PAIRS_R, 1004) == (3, True)
    assert find_pair_by_ref_pos(PAIRS_R, 1007) == (6, False)


def test_normalize_contig():
    assert normalize_contig("chr5", True) == "chr5"
    assert normalize_contig("5", True) == "chr5"
    assert normalize_contig("X", True) == "chrX"
    assert normalize_contig("chr5", False) == "5"
    assert normalize_contig("chrX", False) == "X"
