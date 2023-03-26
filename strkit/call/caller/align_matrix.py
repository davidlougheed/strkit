import parasail

__all__ = [
    "match_score",
    "mismatch_penalty",
    "indel_penalty",
    "dna_bases",
    "dna_matrix",
]


match_score: int = 2  # TODO: parametrize
mismatch_penalty: int = 7  # TODO: parametrize
indel_penalty: int = 5  # TODO: parametrize


# TODO: Customize matrix based on error chances
# Create a substitution matrix for alignment.
# Include IUPAC wildcard bases to allow for motifs with multiple possible motifs.
# Include a wildcard base 'X' for very low-confidence base calls, to prevent needlessly harsh penalties - this is
# inserted into a read in place of bases with low PHRED scores.
dna_bases_str: str = "ACGTRYSWKMBDHVNX"
dna_bases: dict[str, int] = {b: i for i, b in enumerate(dna_bases_str)}
dna_codes: dict[str, tuple[str, ...]] = {
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("G", "C"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "B": ("C", "G", "T"),
    "D": ("A", "C", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "N": ("A", "C", "G", "T"),

    "X": ("A", "C", "G", "T"),  # Special character for matching low-quality bases
}
dna_matrix = parasail.matrix_create(dna_bases_str, match_score, -1 * mismatch_penalty)

for code, code_matches in dna_codes.items():
    for cm in code_matches:
        dna_matrix[dna_bases[code], dna_bases[cm]] = 2 if code != "X" else 0
        dna_matrix[dna_bases[cm], dna_bases[code]] = 2 if code != "X" else 0
