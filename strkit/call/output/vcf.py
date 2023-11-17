import pathlib
import pysam

from typing import Optional

__all__ = ["output_vcf"]


VCF_ALLELE_CNV_TR = "<CNV:TR>"

VCF_TR_INFO_RECORDS: tuple[tuple[str, str, str, str], ...] = (
    ("SVLEN", "A", "Integer", "Length of the structural variant"),
    ("CN", "A", "Float", "Copy number of allele"),
    ("RN", "A", "Integer", "Total number of repeat sequences in this allele"),
    ("RUS", ".", "String", "Repeat unit sequence of the corresponding repeat sequence"),
    ("RUL", ".", "Integer", "Repeat unit length of the corresponding repeat sequence"),
    ("RB", ".", "Integer", "Total number of bases in the corresponding repeat sequence"),
    ("CIRUC", ".", "Float", "Confidence interval around RUC"),
    ("CIRB", ".", "Integer", "Confidence interval around RB"),
)


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
    vh.formats.add("AD", "R", "String", "Read depth for each allele")
    vh.formats.add("DP", 1, "String", "Read depth")
    vh.formats.add("GT", 1, "String", "Genotype")
    vh.formats.add("PS", 1, "String", "Phase set")

    # Add INFO records for tandem repeat copies - these are new to VCF4.4!
    for iv in VCF_TR_INFO_RECORDS:
        vh.info.add(*iv)

    # Add the sample
    vh.add_sample(sample_id)

    return vh


def _vr_pos_key(vr: pysam.VariantRecord) -> int:
    return vr.pos


def output_vcf(
    sample_id: Optional[str],
    reference_file: str,
    results: tuple[dict, ...],
    vcf_path: str,
):
    sample_id_str: str = sample_id or "sample"

    vh = _build_variant_header(sample_id_str, reference_file)
    vf = pysam.VariantFile(vcf_path, "w", header=vh)

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

        has_at_least_one_snv_set = next((r.get("snvs") is not None for r in results), None) is not None

        for result_idx, result in enumerate(results, 1):
            contig = result["contig"]

            if contig != last_contig:
                # we moved on from the last contig, so write the last batch of variant records to the VCF
                _write_contig_vrs()

            # TODO: fix-up coordinates / note positioning
            start = result["start_adj"]
            end = result["end_adj"]
            call = result["call"]

            # ref_cn = result["ref_cn"]
            ref_seq = result["ref_seq"].lower()

            seqs = tuple(map(str.lower, (result["peaks"] or {}).get("seqs", ())))

            # cn_alts = sorted(set(c for c in call if c != ref_cn)) if call is not None else ()
            seq_alts = sorted(set(filter(lambda c: c != ref_seq, seqs)))

            # cn_alleles = (ref_cn, *cn_alts) if call is not None else (".",)
            seq_alleles = (ref_seq, *(seq_alts or (".",))) if call is not None else (".",)

            vr: pysam.VariantRecord = vf.new_record(
                contig=contig,
                start=start,
                stop=end,
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

            vr.samples[sample_id_str]["GT"] = tuple(map(seq_alleles.index, seqs)) if seqs else (".",)

            if has_at_least_one_snv_set:
                vr.samples[sample_id_str].phased = True
                vr.samples[sample_id_str]["PS"] = str(result_idx)

            for snv in result.get("snvs", ()):
                print(snv)

                ref = snv["ref"]
                snv_alts = tuple(filter(lambda v: v != ref, snv["call"]))
                snv_alleles = (ref, *snv_alts)

                snv_vr = vf.new_record(
                    contig=contig,
                    id=snv["id"],
                    # TODO: fix-up coordinates
                    start=snv["pos"],
                    stop=snv["pos"] + 1,
                    alleles=snv_alleles,
                )

                # TODO: write "rcs" for sample SNV genotypes - list of #reads per allele

                snv_vr.samples[sample_id_str]["GT"] = tuple(map(snv_alleles.index, snv["call"]))
                # snv_vr.samples[sample_id_str].phased = True
                snv_vr.samples[sample_id_str]["PS"] = str(result_idx)

                contig_vrs.append(snv_vr)

                # TODO: now is the time to figure out phasing SNVs across different TRs

            contig_vrs.append(vr)

        _write_contig_vrs()  # write the final contig's worth of variant records to the VCF at the end

    finally:
        vf.close()
