import re
from logging import Logger

__all__ = [
    "LocusValidationError",
    "valid_motif",
    "validate_locus",
]

# patterns
RE_VALID_MOTIF = re.compile(r"^[ACGTRYSWKMBDHVN]+$")


# exceptions
class LocusValidationError(ValueError):
    def __init__(self, error_str: str, hint_msg: str):
        self._error_str = error_str
        self._hint_msg = hint_msg
        super().__init__(error_str)

    def log_error(self, logger: Logger) -> None:
        logger.critical(self._error_str)
        logger.critical(self._hint_msg)


# functions
def valid_motif(motif: str) -> bool:
    return RE_VALID_MOTIF.match(motif) is not None


def validate_locus(line: int, start: int, end: int, motif: str) -> None:
    if start >= end:
        raise LocusValidationError(
            f"BED catalog format error: invalid coordinates on line {line}: start ({start}) >= end ({end})",
            "BED catalog: coordinates must be 0-based, half-open - [start, end)",
        )

    if not valid_motif(motif):
        raise LocusValidationError(
            f"BED catalog format error: invalid motif on line {line}: {motif}",
            "BED catalog: motifs must contain only valid IUPAC nucleotide codes.",
        )
