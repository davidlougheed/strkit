import pytest
from strkit.ploidy import PloidyConfig, load_ploidy_config


def test_ploidy_validation():
    c = PloidyConfig.model_validate({
        "default": 2,
    })
    assert c.n_of("NC_000001") == 2

    c = PloidyConfig.model_validate({
        "default": 2,
        "overrides": {
            "chr(M|Y)": 1,
        },
        "ignore": ["chrZ"],
    })
    assert c.n_of("chr3") == 2
    assert c.n_of("chrY") == 1
    assert c.n_of("chrZ") is None

    with pytest.raises(ValueError):
        PloidyConfig.model_validate({"overrides": {"chrM": 1}})

    with pytest.raises(ValueError):
        PloidyConfig.model_validate({"default": 2, "overrides": ["chrY"]})


def test_ploidy_load_and_use():
    haploid = load_ploidy_config("haploid")
    assert haploid.n_of("chr1") == 1
    assert haploid.n_of("chrX") == 1
    assert haploid.n_of("chrM") == 1

    diploid_xx = load_ploidy_config("XX")
    assert diploid_xx.n_of("chr1") == 2
    assert diploid_xx.n_of("chr22") == 2
    assert diploid_xx.n_of("chrX") == 2
    assert diploid_xx.n_of("chrY") == 0
    assert diploid_xx.n_of("chrM") == 1

    diploid_xy = load_ploidy_config("XY")
    assert diploid_xy.n_of("chr1") == 2
    assert diploid_xy.n_of("chr22") == 2
    assert diploid_xy.n_of("chrX") == 1
    assert diploid_xy.n_of("chrM") == 1

    diploid_autosomes = load_ploidy_config("diploid_autosomes")
    assert diploid_autosomes.n_of("chr1") == 2
    assert diploid_autosomes.n_of("chr22") == 2
    assert diploid_autosomes.n_of("chrX") is None
    assert diploid_autosomes.n_of("chrM") == 1
