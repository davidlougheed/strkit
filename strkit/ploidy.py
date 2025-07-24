import re
from pathlib import Path
from pydantic import BaseModel, Field

__all__ = [
    "BUNDLED_PLOIDY_CONFIGS",
    "PLOIDY_OPTIONS_HELP_TEXT",
    "PloidyConfig",
    "load_ploidy_config",
]

PLOIDY_CONFIG_BASE_PATH = Path(__file__).parent / "data" / "ploidy_configs"
DIPLOID_XX_PATH = PLOIDY_CONFIG_BASE_PATH / "diploid_xx.json"
DIPLOID_XY_PATH = PLOIDY_CONFIG_BASE_PATH / "diploid_xy.json"

BUNDLED_PLOIDY_CONFIGS: dict[str, Path] = {
    "haploid": PLOIDY_CONFIG_BASE_PATH / "haploid.json",
    "diploid_autosomes": PLOIDY_CONFIG_BASE_PATH / "diploid_autosomes.json",
    # aliases for diploid_xx.json:
    "diploid_xx": DIPLOID_XX_PATH,
    "XX": DIPLOID_XX_PATH,
    # aliases for diploid_xy.json:
    "diploid_xy": DIPLOID_XY_PATH,
    "XY": DIPLOID_XY_PATH,
}
BUNDLED_PLOIDY_CONFIGS_KEYS = tuple(BUNDLED_PLOIDY_CONFIGS.keys())


def _build_ploidy_options_help():
    ht = ""

    for i, (k, v) in enumerate(BUNDLED_PLOIDY_CONFIGS.items()):
        if i > 0:
            if BUNDLED_PLOIDY_CONFIGS[BUNDLED_PLOIDY_CONFIGS_KEYS[i - 1]] == v:
                ht += f"/{k}"
            else:
                ht += f", {k}"
        else:
            ht += k

    return ht


PLOIDY_OPTIONS_HELP_TEXT = _build_ploidy_options_help()

VALID_KEY_SETS: frozenset[frozenset[str]] = frozenset([
    frozenset({"default"}),
    frozenset({"default", "overrides"}),
    frozenset({"default", "ignore"}),
    frozenset({"default", "ignore", "overrides"}),
])


class PloidyConfig(BaseModel):
    default: int
    overrides: dict[str, int] = Field(default_factory=lambda: {})
    ignore: frozenset[str] = frozenset([])

    def n_of(self, contig: str) -> int | None:
        for k, v in self.overrides.items():
            if re.fullmatch(k, contig):
                return v
        return self.default


def load_ploidy_config(id_or_path: Path | str) -> PloidyConfig:
    """
    Load a ploidy configuration dictionary from a specified location - either a string ID of a built-in ploidy
    configuration or a path to a JSON file.
    :param id_or_path: String ID of a preconfigured ploidy, or a string/Path object pointing to a JSON file.
    :return: Validated, typed version of the loaded ploidy configuration dictionary.
    """

    ploidy_config_path: Path
    if isinstance(id_or_path, str):
        if id_or_path in BUNDLED_PLOIDY_CONFIGS:
            ploidy_config_path = BUNDLED_PLOIDY_CONFIGS[id_or_path]
        else:
            ploidy_config_path = Path(id_or_path)
    else:
        ploidy_config_path = id_or_path

    with open(ploidy_config_path, "r") as fh:
        return PloidyConfig.model_validate_json(fh.read())
