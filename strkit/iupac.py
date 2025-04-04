__all__ = [
    "IUPAC_NUCLEOTIDE_CODES",
    "IUPAC_NUCLEOTIDE_CODES_REVERSE",
    "get_iupac_code_for_nt_set",
]

# IUPAC nucleotide codes representing >1 nucleotide (quasi-"wildcards"):
#  - It's important that the values remain sorted, so we can do a reverse-lookup (see below)
IUPAC_NUCLEOTIDE_CODES: dict[str, tuple[str, ...]] = {
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("C", "G"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "B": ("C", "G", "T"),
    "D": ("A", "C", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "N": ("A", "C", "G", "T"),
}

# Lookup table of {(sorted nucleotides): "<IUPAC code>"}
IUPAC_NUCLEOTIDE_CODES_REVERSE: dict[tuple[str, ...], str] = {
    v: k for k, v in IUPAC_NUCLEOTIDE_CODES.items()
}


def get_iupac_code_for_nt_set(nt_set: set[str]) -> str | None:
    """
    Given a set of standard nucleotides (ATGC), return an IUPAC code which represents the set.
    :param nt_set: A set of nucleotides (A, T, G, or C). Any other base will result in a None return.
    :return: An IUPAC nucleotide code representing the set of nucleotides, or None given an invalid nucleotide set.
    """
    return IUPAC_NUCLEOTIDE_CODES_REVERSE.get(tuple(sorted(nt_set)))
