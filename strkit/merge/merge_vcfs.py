from pathlib import Path
from pysam import VariantFile, VariantHeader

from .. import __version__

__all__ = ["merge_vcfs"]


def check_headers(hdrs: tuple[VariantHeader, ...]):
    # check contigs
    # check STRkit versions
    # check references
    # check info field overlaps
    # TODO: return some info to help?
    pass


def merge_vcfs(paths: tuple[Path, ...], out_file: Path):
    fhs: tuple[VariantFile, ...] = tuple(VariantFile(str(p)) for p in paths)

    try:
        hdrs = tuple(f.header for f in fhs)

        # Check that headers indicate these files can be merged (they are STRkit headers and are compatible)
        check_headers(hdrs)  # TODO

        # TODO: for each STR, need to:
        #  - ensure they're compatible (same ID with same # of loci perhaps?)
        #  - resolve ref (which can vary between) - probably selecting the longest one and
        #     pasting any needed prefix (diff. between sample's ref and current ref) onto it.
        # TODO: for SNVs, we can only deliver heterozygous calls. need to note this for the user.

    finally:
        for fh in fhs:
            fh.close()
