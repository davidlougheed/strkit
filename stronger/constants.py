__all__ = [
    "CALLER_EXPANSIONHUNTER",
    "CALLER_HIPSTR",
    "CALLER_GANGSTR",
    "CALLER_REPEATHMM",
    "CALLER_REPEATHMM_RECALL",
    "CALLER_STRAGLR",
    "CALLER_STRAGLR_RECALL",
    "CALLER_TANDEM_GENOTYPES",
    "CALLER_TANDEM_GENOTYPES_RECALL",
    "CALL_SUPPORTED_CALLERS",

    "SEX_CHROMOSOMES",
    "AUTOSOMES",
    "CHROMOSOMES",

    "MI_CALLERS",
]

CALLER_EXPANSIONHUNTER = "expansionhunter"
CALLER_HIPSTR = "hipstr"
CALLER_GANGSTR = "gangstr"
CALLER_REPEATHMM = "repeathmm"
CALLER_REPEATHMM_RECALL = "repeathmm-recall"
CALLER_STRAGLR = "straglr"
CALLER_STRAGLR_RECALL = "straglr-recall"
CALLER_TANDEM_GENOTYPES = "tandem-genotypes"
CALLER_TANDEM_GENOTYPES_RECALL = "tandem-genotypes-recall"

CALL_SUPPORTED_CALLERS = (
    CALLER_REPEATHMM,
    CALLER_STRAGLR,
    CALLER_TANDEM_GENOTYPES,
)

SEX_CHROMOSOMES = ("chrX", "chrY", "X", "Y")

AUTOSOMES = (
    *map(str, range(1, 23)),
    *(f"chr{i}" for i in range(1, 23)),
)

CHROMOSOMES = (
    *AUTOSOMES,
    *SEX_CHROMOSOMES,
)


MI_CALLERS = (
    CALLER_EXPANSIONHUNTER,
    CALLER_GANGSTR,
    "generic-vcf",
    "hipstr",
    "lobstr",
    "pacmonstr",
    CALLER_REPEATHMM,
    CALLER_REPEATHMM_RECALL,
    CALLER_STRAGLR,
    CALLER_STRAGLR_RECALL,
    CALLER_TANDEM_GENOTYPES,
    CALLER_TANDEM_GENOTYPES_RECALL,
)
