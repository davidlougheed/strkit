import pathlib
import pysam
from os.path import commonprefix

from typing import Optional

from ..utils import cat_strs

__all__ = ["output_vcf"]


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


def _build_variant_header(sample_id: str, reference_file: str) -> pysam.VariantHeader:
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
    # vh.formats.add("PS", 1, "String", "Phase set")  TODO

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


def output_vcf(
    sample_id: Optional[str],
    reference_file: str,
    results: tuple[dict, ...],
    vcf_path: str,
):
    sample_id_str: str = sample_id or "sample"

    vh = _build_variant_header(sample_id_str, reference_file)
    vf = pysam.VariantFile(vcf_path if vcf_path != "stdout" else "-", "w", header=vh)

    contig_vrs: list[pysam.VariantRecord] = []

    def _write_contig_vrs():
        # sort the variant records by position
        contig_vrs.sort(key=_vr_pos_key)

        # write them to the VCF
        for contig_vr in contig_vrs:
            vf.write(contig_vr)

        # clear the contig variant record list for the new contig
        contig_vrs.clear()

    try:
        last_contig = results[0]["contig"] if results else ""

        # has_at_least_one_snv_set = next((r.get("snvs") is not None for r in results), None) is not None
        snvs_written: set[str] = set()

        for result_idx, result in enumerate(results, 1):
            contig = result["contig"]

            if contig != last_contig:
                # we moved on from the last contig, so write the last batch of variant records to the VCF
                _write_contig_vrs()

            ref_start_anchor = result["ref_start_anchor"].upper()

            # ref_cn = result["ref_cn"]
            ref_seq = result["ref_seq"].upper()

            seqs = tuple(map(str.upper, (result["peaks"] or {}).get("seqs", ())))

            # cn_alts = sorted(set(c for c in call if c != ref_cn)) if call is not None else ()
            seq_alts = sorted(set(filter(lambda c: c != ref_seq, seqs)))
            common_suffix_idx = -1 * len(commonprefix(tuple(map(_reversed_str, (ref_seq, *seqs)))))

            call = result["call"]
            # cn_alleles = (ref_cn, *cn_alts) if call is not None else (".",)
            seq_alleles_raw: tuple[str, ...] = (ref_seq, *(seq_alts or (".",))) if call is not None else (".",)
            seq_alleles: list[str] = []

            if call is not None:
                seq_alleles.append(ref_start_anchor + ref_seq[:common_suffix_idx])
                if seq_alts:
                    seq_alleles.extend(ref_start_anchor + a[:common_suffix_idx] for a in seq_alts)
                else:
                    seq_alleles.append(".")

            # seq_alleles = (ref_start_anchor + ref_seq, *(seq_alts or (".",))) if call is not None else (".",)

            start = result.get("start_adj", result["start"]) - len(ref_start_anchor)
            vr: pysam.VariantRecord = vf.new_record(
                contig=contig,
                start=start,
                alleles=seq_alleles,
            )

            # TODO: right now we only call one repeat in one allele and use the canonical motif sequence, even though
            #  there may be multiple motifs present.
            # motif = result["motif"]
            # vr.info["SVLEN"] = math.ceil(end - start)
            # vr.info["CN"] = tuple(a / ref_cn for a in alts) if ref_cn else tuple([1] * len(alts))
            # vr.info["RN"] = tuple([1] * len(alts))
            # vr.info["RUS"] = tuple([motif] * len(alts))
            # vr.info["RUL"] = tuple([len(motif)] * len(alts))

            # TODO: set up tandem repeat info fields

            vr.samples[sample_id_str]["GT"] = tuple(map(seq_alleles_raw.index, seqs)) if seqs else (".",)
            vr.samples[sample_id_str]["DP"] = sum(result["peaks"]["n_reads"])
            vr.samples[sample_id_str]["AD"] = tuple(result["peaks"]["n_reads"])
            vr.samples[sample_id_str]["MC"] = tuple(map(int, result["call"]))  # TODO: support fractional

            # TODO: output phased SNVs
            # if has_at_least_one_snv_set:
            #     vr.samples[sample_id_str].phased = True
            #     vr.samples[sample_id_str]["PS"] = str(result_idx)

            for snv in result.get("snvs", ()):
                snv_id = snv["id"]
                if snv_id in snvs_written:
                    continue
                snvs_written.add(snv_id)

                ref = snv["ref"]
                snv_alts = tuple(filter(lambda v: v != ref, snv["call"]))
                snv_alleles = (ref, *snv_alts)
                snv_pos = snv["pos"]

                snv_vr = vf.new_record(
                    contig=contig,
                    id=snv_id,
                    start=snv_pos,
                    stop=snv_pos + 1,
                    alleles=snv_alleles,
                )

                # TODO: write "rcs" for sample SNV genotypes - list of #reads per allele

                snv_vr.samples[sample_id_str]["GT"] = tuple(map(snv_alleles.index, snv["call"]))
                snv_vr.samples[sample_id_str]["DP"] = sum(snv["rcs"])
                snv_vr.samples[sample_id_str]["AD"] = snv["rcs"]
                # snv_vr.samples[sample_id_str].phased = True
                # snv_vr.samples[sample_id_str]["PS"] = str(result_idx)

                contig_vrs.append(snv_vr)

                # TODO: now is the time to figure out phasing SNVs across different TRs

            contig_vrs.append(vr)

        _write_contig_vrs()  # write the final contig's worth of variant records to the VCF at the end

    finally:
        vf.close()
