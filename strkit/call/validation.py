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
    """
    Determines whether a motif is valid, i.e., can be used by `strkit call`. Here, valid means "composed of IUPAC
    nucleotide codes and no other characters."
    :param motif: The motif to assess the validity of.
    :return: Whether the motif is valid or not.
    """
    return RE_VALID_MOTIF.match(motif) is not None


def validate_locus(line: int, start: int, end: int, motif: str) -> None:
    """
    Validate a locus definition for use by STRkit.
    :param line: Line number, for logging errors in a catalog BED file.
    :param start: Start coordinate; 0-based, inclusive.
    :param end: End coordinate; 0-based, exclusive.
    :param motif: Motif sequence (to be validated).
    """

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
