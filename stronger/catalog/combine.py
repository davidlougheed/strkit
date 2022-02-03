import sys
from typing import List
from ..constants import CALLER_STRAGLR, CHROMOSOMES

__all__ = [
    "combine_catalogs",
]


def combine_catalogs(caller: str, paths: List[str]) -> int:
    if caller != CALLER_STRAGLR:
        sys.stderr.write(f"Error: This command only supports caller '{CALLER_STRAGLR}'\n")
        return 1

    lines = set()

    for path in paths:
        if not path.endswith(".bed"):
            sys.stderr.write(f"Error: Please supply only .bed files from '{CALLER_STRAGLR}'\n")
            return 1

        with open(path, "r") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue

                raw_data = line.strip().split("\t")
                lines.add((raw_data[0], int(raw_data[1]), int(raw_data[2]), raw_data[3]))

    for line in sorted(lines, key=lambda x: (CHROMOSOMES.index(x[0]), x[1])):
        sys.stdout.write("\t".join(map(str, line)) + "\n")

    return 0
