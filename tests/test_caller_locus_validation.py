import pytest
from strkit.call.validation import LocusValidationError, valid_motif, validate_locus


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
        validate_locus(1, 1000, 500, "CAG")

    with pytest.raises(LocusValidationError):
        # start == end, invalid
        validate_locus(1, 1000, 1000, "CAG")

    with pytest.raises(LocusValidationError):
        # invalid motif
        validate_locus(1, 1000, 1200, "(CAG)n")
