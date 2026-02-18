from __future__ import annotations

from pathlib import Path
from pysam import VariantFile, VariantHeader
from typing import TYPE_CHECKING

from .. import __version__
from ..logger import get_main_logger

if TYPE_CHECKING:
    from logging import Logger

__all__ = ["merge_vcfs"]


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
            if len(set(p.id for p in ptrs)) == 1:
                contig = ptrs[0].contig

                refs = tuple(p.ref for p in ptrs)

                if len(set(refs)) > 1:
                    print("TODO ref consolidation!")
                    ref = refs[0]
                    # TODO: if this happens, we need to also
                else:
                    ref = refs[0]

                print(refs)

                # all same locus, can merge
                rec = merged_header.new_record(contig, alleles=(ref,))
                pass
                ptrs = [next(i, None) for i in iters]
            else:
                # we have a record which doesn't apply to all samples
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
