import json
from logging import Logger
from typing import Generator, Iterable

__all__ = [
    "trf_bed_to_eh",
]


def trf_bed_to_eh(trf_data: Iterable[list], _logger: Logger) -> Generator[str, None, None]:
    eh_formatted_loci = []

    for i, item in enumerate(trf_data, 1):
        eh_formatted_loci.append({
            "LocusId": f"Locus{i}",
            "LocusStructure": f"({item[-1]})*",
            "ReferenceRegion": f"{item[0]}:{item[1]}-{item[2]}",
            "VariantType": "Repeat",
        })

    yield json.dumps(eh_formatted_loci, indent=2)
