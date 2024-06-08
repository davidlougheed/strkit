import functools
import logging
import pathlib
import pysam
from os.path import commonprefix
from typing import Optional

from strkit.utils import cat_strs, is_none
from ..allele import get_n_alleles
from ..params import CallParams

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


def build_vcf_header(sample_id: str, reference_file: str) -> pysam.VariantHeader:
    vh = pysam.VariantHeader()  # automatically sets VCF version to 4.2

    # Add an absolute path to the reference genome
    vh.add_meta("reference", f"file://{str(pathlib.Path(reference_file).resolve().absolute())}")

    # Add all contigs from the reference genome file + lengths
    rf = pysam.FastaFile(reference_file)
    try:
        for contig in rf.references:
            vh.contigs.add(contig, length=rf.get_reference_length(contig))
    finally:
        rf.close()

    # Add CNV:TR alt type (symbolic allele: tandem repeat)
    vh.add_meta("ALT", "<ID=CNV:TR,Description=\"Tandem repeat\">")

    # Set up basic VCF formats
    vh.formats.add("AD", ".", "Integer", "Read depth for each allele")
    vh.formats.add("DP", 1, "Integer", "Read depth")
    vh.formats.add("GT", 1, "String", "Genotype")
    vh.formats.add("MC", ".", "Integer", "Motif copy number for each allele")
    vh.formats.add("PS", 1, "Integer", "Phase set")
    vh.formats.add("PM", 1, "String", "Peak-calling method (dist/snv+dist/snv/hp)")

    # Set up VCF info fields
    vh.info.add("MOTIF", 1, "String", "Motif string")
    vh.info.add("REFMC", 1, "Integer", "Motif copy number in the reference genome")

    # Add INFO records for tandem repeat copies - these are new to VCF4.4!  TODO
    # for iv in VCF_TR_INFO_RECORDS:
    #     vh.info.add(*iv)

    # Add the sample
    vh.add_sample(sample_id)

    return vh


def _vr_pos_key(vr: pysam.VariantRecord) -> int:
    return vr.pos


def _reversed_str(s: str) -> str:
    return cat_strs(reversed(s))


@functools.cache
def _blank_entry(n_alleles: int) -> tuple[None, ...]:
    return tuple([None] * n_alleles)


def output_contig_vcf_lines(
    params: CallParams,
    sample_id: str,
    variant_file: pysam.VariantFile,
    results: tuple[dict, ...],
    logger: logging.Logger,
):
    variant_records: list[pysam.VariantRecord] = []

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

        n_alleles: int = get_n_alleles(2, params.sex_chroms, contig) or 2

        res_peaks = result["peaks"] or {}
        peak_seqs = res_peaks.get("seqs", ())
        if any(map(is_none, peak_seqs)):  # Occurs when no consensus for one of the peaks
            logger.error(f"Encountered None in results[{result_idx}].peaks.seqs: {peak_seqs}")
            continue

        seqs = tuple(map(str.upper, peak_seqs))
        if 0 < len(seqs) < n_alleles:
            seqs = tuple([seqs[0]] * n_alleles)

        seq_alts = sorted(set(filter(lambda c: c != ref_seq, seqs)))
        common_suffix_idx = -1 * len(commonprefix(tuple(map(_reversed_str, (ref_seq, *seqs)))))

        call = result["call"]

        seq_alleles_raw: tuple[Optional[str], ...] = (ref_seq, *(seq_alts or (None,))) if call is not None else (".",)
        seq_alleles: list[str] = [ref_start_anchor + ref_seq[:common_suffix_idx]]

        if call is not None and seq_alts:
            seq_alleles.extend(ref_start_anchor + a[:common_suffix_idx] for a in seq_alts)
        else:
            seq_alleles.append(".")

        start = result.get("start_adj", start) - len(ref_start_anchor)

        vr: pysam.VariantRecord = variant_file.new_record(
            contig=contig,
            start=start,
            alleles=seq_alleles,
        )

        vr.info["MOTIF"] = result["motif"]
        vr.info["REFMC"] = result["ref_cn"]

        vr.samples[sample_id]["GT"] = tuple(map(seq_alleles_raw.index, seqs)) if call is not None and seqs \
            else _blank_entry(n_alleles)

        if am := result.get("assign_method"):
            vr.samples[sample_id]["PM"] = am

        if call is not None and res_peaks:
            vr.samples[sample_id]["DP"] = sum(res_peaks["n_reads"])
            vr.samples[sample_id]["AD"] = tuple(res_peaks["n_reads"])
            vr.samples[sample_id]["MC"] = tuple(map(int, call))

            ps = result["ps"]

            try:
                if ps is not None:  # have phase set on call, so mark as phased
                    vr.samples[sample_id].phased = True
                    vr.samples[sample_id]["PS"] = ps
            except TypeError:
                vr.samples[sample_id].phased = False
                logger.error(f"Received bad PS value while writing VCF record at {contig}:{start} - {ps}")
                ps = None

            for snv in result.get("snvs", ()):
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

                snv_vr = variant_file.new_record(
                    contig=contig,
                    id=snv_id,
                    start=snv_pos,
                    stop=snv_pos + 1,
                    alleles=snv_alleles,
                )

                # TODO: write "rcs" for sample SNV genotypes - list of #reads per allele

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
