from __future__ import annotations

from pathlib import Path
from pysam import VariantFile
from typing import TYPE_CHECKING, Callable

from ..call.output import vcf
from ..logger import get_main_logger

if TYPE_CHECKING:
    from logging import Logger
    from pysam import VariantHeader, VariantRecord

__all__ = ["merge_vcfs"]


COMMON_STR_INFO_KEYS = (
    vcf.VCF_INFO_VT,
    vcf.VCF_INFO_MOTIF,
    vcf.VCF_INFO_REFMC,
    vcf.VCF_INFO_BED_START,
    vcf.VCF_INFO_BED_END,
)
COMMON_SNV_INFO_KEYS = (
    vcf.VCF_INFO_VT,
)


def check_headers(headers: tuple[VariantHeader, ...], logger: Logger):
    header_contigs: set[tuple[str, ...]] = set()
    strkit_versions: set[str] = set()
    strkit_catalog_hashes: set[str] = set()

    for idx, hdr in enumerate(headers, 1):
        header_contigs.add(tuple(sorted(hdr.contigs)))
        strkit_version_found = False
        strkit_catalog_hash_found = False
        for line in hdr.records:
            if line.key == "strkitVersion":
                strkit_versions.add(line.value)
                strkit_version_found = True
            elif line.key == "strkitCatalogLociHash":
                strkit_catalog_hashes.add(line.value)
                strkit_catalog_hash_found = True
        if not strkit_version_found:
            logger.critical("did not find strkitVersion in VCF #%d", idx)
            exit(1)
        if not strkit_catalog_hash_found:
            logger.critical("did not find strkitCatalogLociHash in VCF #%d", idx)
            exit(1)

    # check references
    # check info field overlaps

    passed_validation = True

    # check contigs
    if len(header_contigs) > 1:
        logger.warning("more than one header contig set found - are all VCFs using the same reference genome?")
        passed_validation = False
    # check STRkit versions
    if len(strkit_versions) > 1:
        logger.warning(
            "more than one STRkit version found (%s) - output may be different across versions",
            ", ".join(sorted(strkit_versions)),
        )
        passed_validation = False
    # check catalog hashes
    if len(strkit_catalog_hashes) > 1:
        logger.warning(
            "more than one catalog hash found (%s) - are all locus catalogs the same? merging may not work with "
            "different catalogs"
        )
        passed_validation = False

    if passed_validation:
        logger.info("VCF headers passed validation")
    else:
        logger.warning("VCF headers failed validation")

    # TODO: return some info to help?


def merge_str_records(header: VariantHeader, records: list[VariantRecord], logger: Logger) -> VariantRecord:
    contig = records[0].contig

    refs = tuple(p.ref for p in records)
    if len(set(refs)) > 1:
        logger.debug("merge_str_records: encountered >1 STR reference sequence")

    infos = set(tuple(p.info[k] for k in COMMON_STR_INFO_KEYS) for p in records)

    if len(infos) > 1:
        logger.warning("merge_str_records: encountered >1 set of INFO field values for STR")
        # TODO: more reporting

    if len(set(refs)) > 1:
        print("TODO ref consolidation!")
        record_max_anchor = max(records, key=lambda p: p.info["ANCH"])
        max_anchor = record_max_anchor.info["ANCH"]
        ref_max_anchor = record_max_anchor.ref
        start = record_max_anchor.start
        added_prefix_lengths = tuple(max_anchor - p.info["ANCH"] for p in records)

        # TODO: also consider if there's extra on the end or something?
        ref = ref_max_anchor
        # TODO: this misses SNVs that are in the VCF anchor. Maybe STRkit should genotype these.
        # TODO: filter out ref allele
        alts = sorted(set(ref_max_anchor[:added_prefix_lengths[i]] + a for i, p in enumerate(records) for a in p.alts))
    else:
        ref = refs[0]
        alts = sorted(set(a for p in records for a in (p.alts or ()) if a != ref), key=lambda a: (len(a), a))
        start = records[0].start

    # TODO: merge info/format fields

    # all same locus, can merge
    alleles = (ref, *(alts or (".",)))
    print(alleles)
    rec = header.new_record(contig, start=start, alleles=alleles, id=records[0].id)

    # -- add INFO (common across samples) to merged STR record ---------------------------------------------------------
    common_info = next(iter(infos))
    for ki, k in enumerate(COMMON_STR_INFO_KEYS):
        rec.info[k] = common_info[ki]
    # ------------------------------------------------------------------------------------------------------------------

    # Done common stuff, now need to add/transform sample data
    for sample in header.samples:
        pass

    return rec


def merge_snv_records(header: VariantHeader, records: list[VariantRecord], logger: Logger) -> VariantRecord:
    print("TODO: SNV MERGE")
    print([str(r) for r in records])

    # confirm we have a single base ref
    refs = tuple(p.ref for p in records)
    if len(ref_set := set(refs)) > 1:
        logger.critical("merge_snv_records: encountered >1 SNP reference sequence (%s)", ref_set)
        exit(1)
    ref = refs[0]

    infos = set(tuple(p.info[k] for k in COMMON_SNV_INFO_KEYS) for p in records)

    if len(infos) > 1:
        logger.warning("merge_str_records: encountered >1 set of INFO field values for STR")

    alts = sorted(set(a for p in records for a in (p.alts or ()) if a != ref))
    alleles = (ref, *alts)
    rec = header.new_record(contig=records[0].contig, start=records[0].start, alleles=alleles)  # TODO

    common_info = next(iter(infos))
    for ki, k in enumerate(COMMON_SNV_INFO_KEYS):
        rec.info[k] = common_info[ki]

    # Done common stuff, now need to add/transform sample data
    for sample in header.samples:
        pass

    return rec


MERGE_FUNCTION: dict[str, Callable[[VariantHeader, list[VariantRecord], Logger], VariantRecord]] = {
    vcf.VT_STR: merge_str_records,
    vcf.VT_SNV: merge_snv_records,
}


def merge_vcfs(paths: tuple[Path, ...]):
    logger = get_main_logger()  # TODO: inject

    fhs: tuple[VariantFile, ...] = tuple(VariantFile(str(p)) for p in paths)

    try:
        headers = tuple(f.header for f in fhs)

        # Check that headers indicate these files can be merged (they are STRkit headers and are compatible)
        check_headers(headers, logger)  # TODO

        # --------------------------------------------------------------------------------------------------------------

        # Try merging headers
        merged_header = headers[0]
        for hdr in headers[1:]:
            merged_header = merged_header.merge(hdr)

        print(merged_header, end="")

        # --------------------------------------------------------------------------------------------------------------

        iters = tuple(f.fetch() for f in fhs)
        ptrs = [next(i, None) for i in iters]

        while any(p is not None for p in ptrs):
            # find compatible records for merging based on ID (requires STRkit to output variant IDs for everything)
            if all(p is not None for p in ptrs) and len(set(p.id for p in ptrs)) == 1:
                print(MERGE_FUNCTION[ptrs[0].info["VT"]](merged_header, ptrs, logger), end="")
                ptrs = [next(i, None) for i in iters]
            else:
                # we have a record which doesn't apply to all samples - either some VCFs are finished, or the pointers
                # are at different points in the file.
                print("TODO some samples!")
                # TODO: selectively advance

        # TODO: for each STR, need to:
        #  - ensure they're compatible (same ID with same # of loci perhaps?)
        #  - resolve ref (which can vary between) - probably selecting the longest one and
        #     pasting any needed prefix (diff. between sample's ref and current ref) onto it.
        # TODO: for SNVs, we can only deliver heterozygous calls. need to note this for the user.

    finally:
        for fh in fhs:
            fh.close()
