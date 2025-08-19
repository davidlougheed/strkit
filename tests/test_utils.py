from strkit.utils import cis_overlap


def test_cis_overlap():
    assert cis_overlap((1.0, 2.0), (1.5, 2.5))
    assert cis_overlap((1.0, 2.0), (1.999, 2.5))
    assert not cis_overlap((1.0, 2.0), (2.001, 2.5))
    assert not cis_overlap((2.0, 3.0), (1.0, 1.999))
