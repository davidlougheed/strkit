import logging
import pathlib
import pysam
from os.path import commonprefix

from ..utils import cat_strs

__all__ = [
    "build_vcf_header",
    "output_vcf_lines",
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
    vh.add_meta("reference", f"file://{str(pathlib.Path(reference_file).absolute())}")

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


def output_vcf_lines(
    sample_id: str,
    variant_file: pysam.VariantFile,
    results: tuple[dict, ...],
    logger: logging.Logger,
):
    contig_vrs: list[pysam.VariantRecord] = []

    def _write_contig_vrs():
        # sort the variant records by position
        contig_vrs.sort(key=_vr_pos_key)

        # write them to the VCF
        for contig_vr in contig_vrs:
            variant_file.write(contig_vr)

        # clear the contig variant record list for the new contig
        contig_vrs.clear()

    last_contig = results[0]["contig"] if results else ""

    # has_at_least_one_snv_set = next((r.get("snvs") is not None for r in results), None) is not None
    snvs_written: set[str] = set()

    for result_idx, result in enumerate(results, 1):
        contig = result["contig"]

        if contig != last_contig:
            # we moved on from the last contig, so write the last batch of variant records to the VCF
            _write_contig_vrs()

        ref_start_anchor = result["ref_start_anchor"].upper()

        ref_seq = result["ref_seq"].upper()

        seqs = tuple(map(str.upper, (result["peaks"] or {}).get("seqs", ())))

        seq_alts = sorted(set(filter(lambda c: c != ref_seq, seqs)))
        common_suffix_idx = -1 * len(commonprefix(tuple(map(_reversed_str, (ref_seq, *seqs)))))

        call = result["call"]
        seq_alleles_raw: tuple[str, ...] = (ref_seq, *(seq_alts or (".",))) if call is not None else (".",)
        seq_alleles: list[str] = []

        if call is not None:
            seq_alleles.append(ref_start_anchor + ref_seq[:common_suffix_idx])
            if seq_alts:
                seq_alleles.extend(ref_start_anchor + a[:common_suffix_idx] for a in seq_alts)
            else:
                seq_alleles.append(".")

        start = result.get("start_adj", result["start"]) - len(ref_start_anchor)

        if len(seq_alleles) < 2:
            logger.error(f"{contig}:{start} - Don't have enough alleles to write record; have {seq_alleles}")
            continue

        vr: pysam.VariantRecord = variant_file.new_record(
            contig=contig,
            start=start,
            alleles=seq_alleles,
        )

        vr.samples[sample_id]["GT"] = tuple(map(seq_alleles_raw.index, seqs)) if seqs else (".",)
        vr.samples[sample_id]["DP"] = sum(result["peaks"]["n_reads"])
        vr.samples[sample_id]["AD"] = tuple(result["peaks"]["n_reads"])
        vr.samples[sample_id]["MC"] = tuple(map(int, result["call"]))  # TODO: support fractional

        ps = result["ps"]

        if ps is not None:  # have phase set on call, so mark as phased
            vr.samples[sample_id].phased = True
            vr.samples[sample_id]["PS"] = ps

        for snv in result.get("snvs", ()):
            snv_id = snv["id"]
            if snv_id in snvs_written:
                continue
            snvs_written.add(snv_id)

            ref = snv["ref"]
            snv_alts = tuple(filter(lambda v: v != ref, snv["call"]))
            snv_alleles = (ref, *snv_alts)
            snv_pos = snv["pos"]

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

            contig_vrs.append(snv_vr)

        contig_vrs.append(vr)

    _write_contig_vrs()  # write the final contig's worth of variant records to the VCF at the end
