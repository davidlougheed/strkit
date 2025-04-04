from strkit.iupac import get_iupac_code_for_nt_set


def test_get_iupac_code():
    assert get_iupac_code_for_nt_set({"A", "T"}) == "W"
    assert get_iupac_code_for_nt_set({"A", "C", "G", "T"}) == "N"
    assert get_iupac_code_for_nt_set({"A", "T", "C", "G"}) == "N"
    assert get_iupac_code_for_nt_set({"A", "T", "C"}) == "H"
    assert get_iupac_code_for_nt_set({"A", "T", "C", "Z"}) is None
    assert get_iupac_code_for_nt_set({"A", "T", "C", ":)"}) is None
    assert get_iupac_code_for_nt_set({"A", "T", "C", ""}) is None
