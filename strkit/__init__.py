from pathlib import Path

__all__ = [
    "__version__",
]

with open(Path(__file__).parent / "VERSION", "r") as vf:
    __version__ = vf.read().strip()
