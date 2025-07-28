import functools
import logging

from collections import Counter
from os.path import commonprefix
from pathlib import Path
from pysam import FastaFile, VariantFile, VariantHeader, VariantRecord
from typing import Iterable

from strkit.utils import cat_strs, is_none, idx_0_getter
from ..params import CallParams
from ..types import LocusResult
from ..utils import cn_getter

__all__ = [
    "build_vcf_header",
    "output_contig_vcf_lines",
]


# VCF_ALLELE_CNV_TR = "<CNV:TR>"

# VCF_TR_INFO_RECORDS: tuple[tuple[str, str, str, str], ...] = (
#     ("SVLEN", "A", "Integer", "Length of the structural variant"),
#     ("CN", "A", "Float", "Copy number of allele"),
#     ("RN", "A", "Integer", "Total number of repeat sequences in this allele"),
#     ("RUS", ".", "String", "Repeat unit sequence of the corresponding repeat sequence"),
#     ("RUL", ".", "Integer", "Repeat unit length of the corresponding repeat sequence"),
#     ("RB", ".", "Integer", "Total number of bases in the corresponding repeat sequence"),
#     ("CIRUC", ".", "Float", "Confidence interval around RUC"),
#     ("CIRB", ".", "Integer", "Confidence interval around RB"),
# )

VCF_INFO_VT = "VT"
VCF_INFO_MOTIF = "MOTIF"
VCF_INFO_REFMC = "REFMC"
VCF_INFO_ANCH = "ANCH"

VT_STR = "str"
VT_SNV = "snv"


def iter_to_upper(x: Iterable[str]) -> Iterable[str]:
    # noinspection PyTypeChecker
    return map(str.upper, x)


def build_vcf_header(sample_id: str, reference_file: str) -> VariantHeader:
    vh = VariantHeader()  # automatically sets VCF version to 4.2

    # Add an absolute path to the reference genome
    vh.add_meta("reference", f"file://{str(Path(reference_file).resolve().absolute())}")

    # Add all contigs from the reference genome file + lengths
    with FastaFile(reference_file) as rf:
        for contig in rf.references:
            vh.contigs.add(contig, length=rf.get_reference_length(contig))

    # Add CNV:TR alt type (symbolic allele: tandem repeat)
    # vh.add_meta("ALT", "<ID=CNV:TR,Description=\"Tandem repeat\">")

    # Set up basic VCF formats
    vh.formats.add("AD", ".", "Integer", "Read depth for each allele")
    vh.formats.add("ANCL", ".", "Integer", "Anchor length for the ref and each alt, five-prime of TR sequence")
    vh.formats.add("CONS", ".", "String", "Consensus methods used for each alt (single/poa/best_rep)")
    vh.formats.add("DP", 1, "Integer", "Read depth")
    vh.formats.add("DPS", 1, "Integer", "Read depth (supporting reads only)")
    vh.formats.add("GT", 1, "String", "Genotype")
    vh.formats.add("MC", ".", "Integer", "Motif copy number for each allele")
    vh.formats.add("MCCI", ".", "String", "Motif copy number 95% confidence interval for each allele")
    vh.formats.add("MCRL", ".", "String", "Read-level motif copy numbers for each allele")
    vh.formats.add("MMAS", 1, "Float", "Mean model (candidate TR sequence) alignment score across reads.")
    vh.formats.add("NSNV", 1, "Integer", "Number of supporting SNVs for the STR peak-call")
    vh.formats.add("PS", 1, "Integer", "Phase set")
    vh.formats.add("PM", 1, "String", "Peak-calling method (dist/snv+dist/snv/hp)")

    # Set up VCF info fields
    vh.info.add(VCF_INFO_VT, 1, "String", "Variant record type (str/snv)")
    vh.info.add(VCF_INFO_MOTIF, 1, "String", "Motif string")
    vh.info.add(VCF_INFO_REFMC, 1, "Integer", "Motif copy number in the reference genome")
    vh.info.add(VCF_INFO_ANCH, 1, "Integer", "Five-prime anchor size")

    # Add INFO records for tandem repeat copies - these are new to VCF4.4!  TODO
    # for iv in VCF_TR_INFO_RECORDS:
    #     vh.info.add(*iv)

    # Add the sample
    vh.add_sample(sample_id)

    return vh


def _vr_pos_key(vr: VariantRecord) -> int:
    return vr.pos


def _reversed_str(s: str) -> str:
    return cat_strs(reversed(s))


@functools.cache
def _blank_entry(n_alleles: int) -> tuple[None, ...]:
    return tuple([None] * n_alleles)


def output_contig_vcf_lines(
    params: CallParams,
    sample_id: str,
    variant_file: VariantFile,
    results: tuple[LocusResult, ...],
    logger: logging.Logger,
) -> None:
    variant_records: list[VariantRecord] = []

    # has_at_least_one_snv_set = next((r.get("snvs") is not None for r in results), None) is not None
    snvs_written: set[str] = set()

    for result_idx, result in enumerate(results, 1):
        contig = result["contig"]
        start = result["start"]

        if "ref_start_anchor" not in result:
            logger.debug(f"No ref anchor for {contig}:{start}; skipping VCF output for locus")
            continue

        ref_start_anchor = result["ref_start_anchor"].upper()
        ref_seq = result["ref_seq"].upper()

        n_alleles: int = params.ploidy_config.n_of(contig)

        res_reads = result["reads"]
        res_peaks = result["peaks"] or {}

        peak_seqs_and_methods = {(seq.upper() if seq else seq): method for seq, method in res_peaks.get("seqs", [])}
        peak_seqs: tuple[str, ...] = tuple(peak_seqs_and_methods.keys())
        peak_start_anchor_seqs: list[str] = list(map(idx_0_getter, res_peaks.get("start_anchor_seqs", [])))

        if any(map(is_none, peak_seqs)):  # Occurs when no consensus for one of the peaks
            logger.error(f"Encountered None in results[{result_idx}].peaks.seqs: {peak_seqs}")
            continue

        if any(map(is_none, peak_start_anchor_seqs)):  # Occurs when no consensus for one of the peaks
            logger.error(f"Encountered None in results[{result_idx}].peaks.start_anchor_seqs: {peak_start_anchor_seqs}")
            continue

        peak_start_anchor_seqs_upper = tuple(iter_to_upper(peak_start_anchor_seqs))
        common_anchor_prefix = commonprefix([ref_start_anchor, *peak_start_anchor_seqs_upper])
        # anchor_offset = how many bases we can cut off from the front of the anchor
        # since they're shared between all alleles - yields a more compact representation.
        #  - we need to leave one base as an anchor for VCF compliance though, thus the min(...)
        anchor_offset = min(len(common_anchor_prefix), params.vcf_anchor_size - 1)

        ref_start_anchor = ref_start_anchor[anchor_offset:]
        ref_seq_with_anchor = ref_start_anchor + ref_seq

        seqs_with_anchors: list[tuple[str, str]] = list(
            zip(peak_seqs, map(lambda a: a[anchor_offset:], peak_start_anchor_seqs_upper))
        )

        if 0 < len(peak_seqs) < n_alleles:
            peak_seqs = tuple([peak_seqs[0]] * n_alleles)
            seqs_with_anchors = [seqs_with_anchors[0]] * n_alleles

        seq_alts = sorted(
            set(filter(lambda c: not (c[1] + c[0] == ref_seq_with_anchor), seqs_with_anchors)),
            key=lambda c: c[1] + c[0]
        )

        call = result["call"]
        call_95_cis = result["call_95_cis"]

        seq_alleles_raw: tuple[str | None, ...] = (
            ((ref_seq, ref_start_anchor), *(seq_alts or (None,)))
            if call is not None
            else ()
        )

        seq_alleles: list[str] = [ref_seq_with_anchor]

        if call is not None and seq_alts:
            # If we have a complete deletion, including the anchor, use a symbolic allele meaning "upstream deletion"
            for alt_tr_seq, alt_anchor in seq_alts:
                if not alt_tr_seq and not alt_anchor:
                    seq_alleles.append("*")
                    continue
                seq_alleles.append(alt_anchor + alt_tr_seq)
        else:
            seq_alleles.append(".")

        start = result.get("start_adj", start) - len(ref_start_anchor)

        vr: VariantRecord = variant_file.new_record(
            contig=contig,
            start=start,
            alleles=seq_alleles,
        )

        vr.info[VCF_INFO_VT] = VT_STR
        vr.info[VCF_INFO_MOTIF] = result["motif"]
        vr.info[VCF_INFO_REFMC] = result["ref_cn"]
        vr.info[VCF_INFO_ANCH] = params.vcf_anchor_size - anchor_offset

        vr.samples[sample_id]["GT"] = (
            tuple(map(seq_alleles_raw.index, seqs_with_anchors))
            if call is not None and peak_seqs
            else _blank_entry(n_alleles)
        )

        if am := result.get("assign_method"):
            vr.samples[sample_id]["PM"] = am

        str_snvs = result.get("snvs", ())
        if str_snvs:
            # Record number of support SNVs for the locus
            vr.samples[sample_id]["NSNV"] = len(str_snvs)

        vr.samples[sample_id]["DP"] = len(res_reads)
        vr.samples[sample_id]["MMAS"] = result.get("mean_model_align_score")

        if call is not None and res_peaks:
            vr.samples[sample_id]["DPS"] = sum(res_peaks["n_reads"])
            vr.samples[sample_id]["AD"] = tuple(res_peaks["n_reads"])
            vr.samples[sample_id]["MC"] = tuple(map(int, call))
            vr.samples[sample_id]["MCCI"] = tuple(f"{x[0]}-{x[1]}" for x in call_95_cis)

            vr.samples[sample_id]["ANCL"] = tuple(len(ar[1]) for ar in seq_alleles_raw if ar is not None)

            # For each alt, mention which consensus method was used to obtain the sequence.
            cons = tuple(
                peak_seqs_and_methods[ar[0]] for ar in seq_alleles_raw[1:] if ar is not None
            )
            vr.samples[sample_id]["CONS"] = cons if cons else (".",)

            # Produces a histogram-like format for read-level copy numbers
            # e.g., for two alleles with 8 and 9 copy-number respectively, we may get: 7x1|8x10|9x1,8x2|9x12
            vr.samples[sample_id]["MCRL"] = tuple(
                "|".join(
                    map(
                        lambda pair: "x".join(map(str, pair)),
                        sorted(
                            Counter(
                                map(cn_getter, filter(lambda r: r.get("p") == pi, res_reads.values()))
                            ).items()
                        )
                    )
                )
                for pi in range(res_peaks["modal_n"])
            )

            ps = result["ps"]

            try:
                if ps is not None:  # have phase set on call, so mark as phased
                    vr.samples[sample_id].phased = True
                    vr.samples[sample_id]["PS"] = ps
            except TypeError:
                vr.samples[sample_id].phased = False
                logger.error(f"Received bad PS value while writing VCF record at {contig}:{start} - {ps}")
                ps = None

            for snv in str_snvs:
                snv_id = snv["id"]
                if snv_id in snvs_written:
                    continue
                snvs_written.add(snv_id)

                ref = snv["ref"]
                snv_alts = tuple(sorted(set(filter(lambda v: v != ref, snv["call"]))))
                snv_alleles = (ref, *snv_alts)
                snv_pos = snv["pos"]

                if len(snv_alleles) < 2:
                    logger.error(f"Error while writing VCF: SNV ({snv_id}) at {contig}:{snv_pos+1} has no alts")
                    continue

                snv_vr: VariantRecord = variant_file.new_record(
                    contig=contig,
                    id=snv_id,
                    start=snv_pos,
                    stop=snv_pos + 1,
                    alleles=snv_alleles,
                )

                snv_vr.info[VCF_INFO_VT] = VT_SNV

                snv_vr.samples[sample_id]["GT"] = tuple(map(snv_alleles.index, snv["call"]))
                snv_vr.samples[sample_id]["DP"] = sum(snv["rcs"])
                snv_vr.samples[sample_id]["AD"] = snv["rcs"]

                if ps is not None:
                    snv_vr.samples[sample_id].phased = True
                    snv_vr.samples[sample_id]["PS"] = ps

                variant_records.append(snv_vr)

        variant_records.append(vr)

    # sort the variant records by position
    variant_records.sort(key=_vr_pos_key)

    # write them to the VCF
    for vrr in variant_records:
        variant_file.write(vrr)
