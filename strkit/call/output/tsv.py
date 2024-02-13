import sys
from typing import Union

__all__ = ["output_tsv"]


def _cn_to_str(cn: Union[int, float]) -> str:
    return f"{cn:.1f}" if isinstance(cn, float) else str(cn)


def output_tsv(results: tuple[dict, ...], has_snv_vcf: bool):
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
            _cn_to_str(ref_cn) if ref_cn is not None else ".",
            ",".join(map(_cn_to_str, sorted(r["cn"] for r in reads.values()))) if reads else ".",
            "|".join(map(_cn_to_str, res["call"])) if has_call else ".",
            ("|".join("-".join(map(_cn_to_str, gc)) for gc in res["call_95_cis"]) if has_call else "."),
            # *((res["assign_method"] if has_call else ".",) if incorporate_snvs else ()),
            *((res["assign_method"] if has_call else ".",) if has_snv_vcf else ()),

            # ("|".join(map(lambda x: f"{x:.5f}", res["peaks"]["means"][:n_peaks]))
            #  if has_call and n_peaks <= 2 else "."),
            # ("|".join(map(lambda x: f"{x:.5f}", res["peaks"]["weights"][:n_peaks]))
            #  if has_call and n_peaks <= 2 else "."),
            # ("|".join(map(lambda x: f"{x:.5f}", res["peaks"]["stdevs"][:n_peaks]))
            #  if has_call and n_peaks <= 2 else "."),
        )) + "\n")
