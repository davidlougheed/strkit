import sys
from logging import Logger

__all__ = [
    "trgt_bed_to_bed4",
    "trf_or_strkit_bed_to_trgt",
]


def trgt_bed_to_bed4(trgt_data: list, logger: Logger):
    """
    Converts a TRGT repeat catalog to the STRkit/BED4 catalog format.
    :param trgt_data: The loaded TRGT catalog (split by tab).
    :param logger: A logger instance for issuing conversion failure warnings.
    """

    for line, data in enumerate(trgt_data, 1):
        structure_data = {j[0]: j[1] for j in (i.split("=") for i in data[3].split(";"))}
        motifs = structure_data["MOTIFS"].split(",")

        if len(motifs) > 1:
            data_str = "\t".join(data)
            logger.warning(f"Could not convert complex locus at line {line}: {data_str}")
            continue

        sys.stdout.write("\t".join((*data[:3], motifs[0])) + "\n")


def trf_or_strkit_bed_to_trgt(trf_data: list, _logger: Logger):
    """
    Convets a TRF- or STRkit-formatted BED (motif-last) to a basic version of a TRGT catalog.
    :param trf_data: The loaded BED catalog data.
    :param _logger: Logger instance (unused).
    """

    for i, item in enumerate(trf_data):
        motif = trf_data[-1]
        sys.stdout.write("\t".join((*trf_data[:3], f"ID=locus{i};MOTIFS={motif};STRUC=({motif})n")) + "\n")
