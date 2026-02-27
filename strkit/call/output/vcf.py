from __future__ import annotations

from collections import Counter
from os.path import commonprefix
from typing import Iterable, TYPE_CHECKING

from strkit import vcf_utils as vu
from strkit.utils import is_none, idx_0_getter
from ..utils import cn_getter

if TYPE_CHECKING:
    from logging import Logger
    from pysam import VariantFile, VariantRecord
    from ..params import CallParams
    from ..types import LocusResult

__all__ = [
    "output_contig_vcf_lines",
]


def iter_to_upper(x: Iterable[str]) -> Iterable[str]:
    # noinspection PyTypeChecker
    return map(str.upper, x)


def _vr_pos_key(vr: VariantRecord) -> int:
    return vr.pos


def create_result_vcf_records(
    params: CallParams,
    variant_file: VariantFile,
    sample_id: str,
    result_idx: int,
    result: LocusResult,
    snvs_written: set[str],  # to avoid writing the same SNV twice across records
    logger: Logger,
):
    variant_records: list[VariantRecord] = []

    contig = result["contig"]
    start = result["start"]

    if "ref_start_anchor" not in result:
        logger.debug("No ref anchor for %s:%d; skipping VCF output for locus", contig, start)
        raise vu.record.SkipWritingLocus()

    ref_start_anchor = result["ref_start_anchor"].upper()
    ref_seq = result["ref_seq"].upper()

    n_alleles: int = params.ploidy_config.n_of(contig)

    res_reads = result["reads"]
    res_peaks = result["peaks"] or {}

    peak_seqs_and_methods = {(seq.upper() if seq else seq): method for seq, method in res_peaks.get("seqs", [])}
    peak_seqs: tuple[str, ...] = tuple(peak_seqs_and_methods.keys())
    peak_start_anchor_seqs: list[str] = list(map(idx_0_getter, res_peaks.get("start_anchor_seqs", [])))

    if any(map(is_none, peak_seqs)):  # Occurs when no consensus for one of the peaks
        logger.error("Encountered None in results[%d].peaks.seqs: %s", result_idx, peak_seqs)
        raise vu.record.SkipWritingLocus()

    if any(map(is_none, peak_start_anchor_seqs)):  # Occurs when no consensus for one of the peaks
        logger.error("Encountered None in results[%d].peaks.start_anchor_seqs: %s", result_idx, peak_start_anchor_seqs)
        raise vu.record.SkipWritingLocus()

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
        id=result["locus_id"],
        contig=contig,
        start=start,
        alleles=seq_alleles,
    )

    vr.info[vu.header.VCF_INFO_VT.key] = vu.header.VT_STR
    vr.info[vu.header.VCF_INFO_MOTIF.key] = result["motif"]
    vr.info[vu.header.VCF_INFO_REFMC.key] = result["ref_cn"]
    vr.info[vu.header.VCF_INFO_BED_START.key] = result["start"]
    vr.info[vu.header.VCF_INFO_BED_END.key] = result["end"]
    vr.info[vu.header.VCF_INFO_ANCH.key] = params.vcf_anchor_size - anchor_offset

    try:
        vr.samples[sample_id]["GT"] = vu.record.genotype_indices(
            alleles=seq_alleles_raw,
            call=None if call is None or not peak_seqs else seqs_with_anchors,
            n_alleles=n_alleles,
        )
    except ValueError:
        logger.error(
            "results[%d] (locus_id=%s): one of %s not in %s",
            result_idx, result["locus_id"], seqs_with_anchors, seq_alleles_raw)
        raise vu.record.SkipWritingLocus()

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
            logger.error("Received bad PS value while writing VCF record at %s:%d - %s", contig, start, ps)
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
                logger.error("Error while writing VCF: SNV (%s) at %s:%d has no alts", snv_id, contig, snv_pos + 1)
                continue

            snv_vr: VariantRecord = variant_file.new_record(
                contig=contig,
                id=snv_id,
                start=snv_pos,
                stop=snv_pos + 1,
                alleles=snv_alleles,
            )

            snv_vr.info[vu.header.VCF_INFO_VT.key] = vu.header.VT_SNV

            snv_vr.samples[sample_id]["GT"] = vu.record.genotype_indices(snv_alleles, snv["call"], n_alleles)
            snv_vr.samples[sample_id]["DP"] = sum(snv["rcs"])
            snv_vr.samples[sample_id]["AD"] = snv["rcs"]

            if ps is not None:
                snv_vr.samples[sample_id].phased = True
                snv_vr.samples[sample_id]["PS"] = ps

            variant_records.append(snv_vr)

        variant_records.append(vr)

    return variant_records


def output_contig_vcf_lines(
    params: CallParams,
    sample_id: str,
    variant_file: VariantFile,
    results: tuple[LocusResult, ...],
    logger: Logger,
) -> None:
    variant_records: list[VariantRecord] = []
    snvs_written: set[str] = set()

    for result_idx, result in enumerate(results, 1):
        try:
            variant_records.extend(
                create_result_vcf_records(
                    params,
                    variant_file,
                    sample_id,
                    result_idx,
                    result,
                    snvs_written,
                    logger,
                )
            )
        except vu.record.SkipWritingLocus:
            pass  # just skipping the locus, nothing to do here as we've already logged
        except Exception as e:  # fallback if we didn't handle a case properly in create_result_vcf_records
            logger.exception("Error while writing VCF: unhandled exception at results[%d]", result_idx, exc_info=e)

    # sort the variant records by position to ensure we write a properly bgzippable/tabix-indexable VCF
    variant_records.sort(key=_vr_pos_key)

    # write them to the VCF
    for vrr in variant_records:
        variant_file.write(vrr)
