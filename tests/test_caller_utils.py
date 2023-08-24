from strkit.call.utils import find_pair_by_ref_pos_py, find_pair_by_ref_pos, normalize_contig, round_to_base_pos

#         A  A       T  T       C  G       C  C       C  C       A  A       A  A       A  C
PAIRS = [(0, 1000), (1, 1001), (2, 1003), (3, 1004), (4, 1005), (5, 1006), (6, 1008), (7, 1009)]
SNVS = ((1003, "C"), (1009, "A"))


def test_find_pair_by_ref_pos_py():
    assert find_pair_by_ref_pos_py(PAIRS, 1004) == (3, True)
    assert find_pair_by_ref_pos_py(PAIRS, 1007) == (6, False)


def test_find_pair_by_ref_pos():
    assert find_pair_by_ref_pos(PAIRS, 1004) == (3, True)
    assert find_pair_by_ref_pos(PAIRS, 1007) == (6, False)

    assert find_pair_by_ref_pos(PAIRS, 1004) == find_pair_by_ref_pos_py(PAIRS, 1004)
    assert find_pair_by_ref_pos(PAIRS, 1007) == find_pair_by_ref_pos_py(PAIRS, 1007)


def test_normalize_contig():
    assert normalize_contig("chr5", True) == "chr5"
    assert normalize_contig("5", True) == "chr5"
    assert normalize_contig("X", True) == "chrX"
    assert normalize_contig("chr5", False) == "5"
    assert normalize_contig("chrX", False) == "X"


def test_round_to_base_pos():
    assert round_to_base_pos(2.5, 3) == 8/3  # 2.666666
    assert round_to_base_pos(1.5, 2) == 1.5
    assert round_to_base_pos(1.3, 4) == 1.25
    assert round_to_base_pos(5.9, 4) == 6
