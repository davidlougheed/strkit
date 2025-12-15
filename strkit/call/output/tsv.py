from __future__ import annotations
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..types import LocusResult

__all__ = ["output_tsv"]


def output_tsv(results: tuple[LocusResult, ...], has_snv_vcf: bool):
    for res in results:
        has_call = res["call"] is not None
        # n_peaks = res["peaks"]["modal_n"]

        ref_cn = res.get("ref_cn")
        reads = res.get("reads")

        sys.stdout.write("\t".join((
            res["contig"],
            str(res["start"]),
            str(res["end"]),
            res["motif"],
            str(ref_cn) if ref_cn is not None else ".",
            ",".join(map(str, sorted(r["cn"] for r in reads.values()))) if reads else ".",
            "|".join(map(str, res["call"])) if has_call else ".",
            ("|".join("-".join(map(str, gc)) for gc in res["call_95_cis"]) if has_call else "."),
            # *((res["assign_method"] if has_call else ".",) if incorporate_snvs else ()),
            *((res["assign_method"] if has_call else ".",) if has_snv_vcf else ()),

            # ("|".join(map(lambda x: f"{x:.5f}", res["peaks"]["means"][:n_peaks]))
            #  if has_call and n_peaks <= 2 else "."),
            # ("|".join(map(lambda x: f"{x:.5f}", res["peaks"]["weights"][:n_peaks]))
            #  if has_call and n_peaks <= 2 else "."),
            # ("|".join(map(lambda x: f"{x:.5f}", res["peaks"]["stdevs"][:n_peaks]))
            #  if has_call and n_peaks <= 2 else "."),
        )) + "\n")
