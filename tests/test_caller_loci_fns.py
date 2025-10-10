import pytest
from strkit_rust_ext import STRkitLocus
from strkit.call.loci import LocusValidationError, valid_motif, validate_locus, parse_last_column


@pytest.mark.parametrize("motif,valid", [
    ("CAG", True),
    ("CAGN", True),
    ("CAGX", False),
    ("(CAG)n", False),
    ("XX", False),
])
def test_valid_motif(motif, valid):
    assert valid_motif(motif) == valid


def test_validate_locus():
    with pytest.raises(LocusValidationError):
        # start > end, invalid
        validate_locus(STRkitLocus(1, "locus1", "1", 1000, 500, "CAG", 2, 70))

    with pytest.raises(LocusValidationError):
        # start == end, invalid
        validate_locus(STRkitLocus(1, "locus1", "1", 1000, 1000, "CAG", 2, 70))

    with pytest.raises(LocusValidationError):
        # invalid motif
        validate_locus(STRkitLocus(1, "locus1", "1", 1000, 1200, "(CAG)n", 2, 70))


@pytest.mark.parametrize("t_idx,val,parsed", [
    (0, "CAG", {"id": "locus0", "motif": "CAG"}),
    (0, "MOTIF=CAG", {"id": "locus0", "motif": "CAG"}),
    (0, "Motif = CAG", {"id": "locus0", "motif": "CAG"}),
    (0, "ID=HTT; MOTIF = CAG", {"id": "HTT", "motif": "CAG"}),
    (0, "ID = HTT ;motif = CAG", {"id": "HTT", "motif": "CAG"}),
    (0, "id=HTT ; motif=cag", {"id": "HTT", "motif": "CAG"}),
])
def test_parse_last_column(t_idx: int, val: str, parsed: dict):
    assert parse_last_column(t_idx, val) == parsed
