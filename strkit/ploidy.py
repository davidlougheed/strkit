import re
from dataclasses import dataclass
from pathlib import Path
from strkit.json import json

__all__ = [
    "BUNDLED_PLOIDY_CONFIGS",
    "PloidyConfig",
    "validate_ploidy_config",
    "load_ploidy_config",
]


BUNDLED_PLOIDY_CONFIGS: dict[str, Path] = {
    "haploid": Path(__file__).parent / "data" / "ploidy_configs" / "haploid.json",
    "XX": Path(__file__).parent / "data" / "ploidy_configs" / "diploid_xx.json",
    "XY": Path(__file__).parent / "data" / "ploidy_configs" / "diploid_xy.json",
}


@dataclass
class PloidyConfig:
    default: int
    overrides: dict[str, int]

    def n_of(self, contig: str):
        for k, v in self.overrides.items():
            if re.fullmatch(k, contig):
                return v
        return self.default


def validate_ploidy_config(p: dict) -> PloidyConfig:
    """
    Validate the structure of a ploidy configuration dictionary and return the validated, typed version of the ploidy
    configuration dictionary.
    :param p: Unvalidated ploidy configuration dictionary.
    :return: Validated, typed version of the ploidy configuration dictionary.
    """

    p_keys = set(p.keys())

    if p_keys != {"default"} and p_keys != {"default", "overrides"}:
        raise ValueError(f"Invalid ploidy configuration: invalid set of keys {p_keys}")

    default_ploidy = p["default"]
    if not isinstance(default_ploidy, int):
        raise ValueError("Ploidy configuration default must be integer")

    overrides = p.get("overrides", {})

    for k, v in overrides.items():
        if not isinstance(v, int):
            raise ValueError("Ploidy configuration override value must be integer")

    return PloidyConfig(default=default_ploidy, overrides=overrides)


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
        return validate_ploidy_config(json.loads(fh.read()))
