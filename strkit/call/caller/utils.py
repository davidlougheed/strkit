__all__ = [
    "normalize_contig",
    "round_to_base_pos",
]


def normalize_contig(contig: str, has_chr: bool) -> str:
    return ("chr" if has_chr else "") + contig.replace("chr", "")


def round_to_base_pos(x, motif_size: int) -> float:
    return round(float(x) * motif_size) / motif_size
