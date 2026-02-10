from __future__ import annotations

import re
import time

from functools import cache
from strkit_rust_ext import STRkitLocus, STRkitLocusBlock
from typing import Iterable, TypedDict, TYPE_CHECKING

if TYPE_CHECKING:
    from logging import Logger
    from queue import Queue
    from .params import CallParams


__all__ = [
    "LocusValidationError",
    "valid_motif",
    "validate_locus",
    "parse_loci_bed",
    "parse_last_column",
    "load_loci",
]

# patterns
RE_VALID_MOTIF = re.compile(r"^[ACGTRYSWKMBDHVN]+$")
RE_LOCUS_INFO_DELIM = re.compile(r"; *")
RE_LOCUS_INFO_ASSIGN = re.compile(r" *= *")


# exceptions

class LocusValidationError(ValueError):
    def __init__(self, error_str: str, hint_msg: str):
        self._error_str = error_str
        self._hint_msg = hint_msg
        super().__init__(error_str)

    def log_error(self, logger: Logger) -> None:
        logger.critical(self._error_str)
        logger.critical(self._hint_msg)


# types

class LastColumnData(TypedDict, total=False):
    id: str
    motif: str


# functions

def valid_motif(motif: str) -> bool:
    """
    Determines whether a motif is valid, i.e., can be used by `strkit call`. Here, valid means "composed of IUPAC
    nucleotide codes and no other characters."
    :param motif: The motif to assess the validity of.
    :return: Whether the motif is valid or not.
    """
    return RE_VALID_MOTIF.match(motif) is not None


def validate_locus(locus: STRkitLocus) -> None:
    """
    Validate a locus definition for use by STRkit.
    :param locus: Locus definition data class (implemented in strkit_rust_ext.)
    """

    line = locus.t_idx

    if locus.left_coord >= locus.right_coord:
        raise LocusValidationError(
            f"BED catalog format error: invalid coordinates on line {line}: start ({locus.left_coord}) >= end "
            f"({locus.right_coord})",
            "BED catalog: coordinates must be 0-based, half-open - [start, end)",
        )

    if not valid_motif(locus.motif):
        raise LocusValidationError(
            f"BED catalog format error: invalid motif on line {line}: {locus.motif}",
            "BED catalog: motifs must contain only valid IUPAC nucleotide codes.",
        )


def parse_loci_bed(loci_file: str) -> Iterable[tuple[str, ...]]:
    with open(loci_file, "r") as tf:
        yield from (
            tuple(line.split("\t"))
            for line in (s.strip() for s in tf)
            if line and not line.startswith("#")  # skip blank lines and comment lines
        )


def parse_last_column(t_idx: int, val: str) -> LastColumnData:
    """
    Given the last line of a 4-column BED file, parse it into a LastColumnData dictionary. For STRkit, the last column
    of a catalog BED file line can either be a motif or a list of semicolon-delimited key-value pairs with = used for
    assignment. These key-value pairs specify the locus motif + other locus information, like a custom locus ID.
    For example, possible valid values are: CAG    MOTIF=CAG    ID=HTT;MOTIF=CAG
    :param t_idx: Locus index in BED file.
    :param val: Value of the last column of the locus BED file entry.
    :return: Dictionary of {"id": "<...>", "motif": "<...>"}
    """

    # We can have a straightforward motif rather than key/value pairs:
    if ";" not in val and "=" not in val:
        return {"id": f"locus{t_idx}", "motif": val}

    # Otherwise, we need to parse out the key/value pairs.

    hint = "BED catalog: last column must either be motif or ID=locusID;MOTIF=motif"
    err = LocusValidationError(f"BED catalog format error: could not parse last column on line {t_idx}: {val}", hint)

    res: LastColumnData = {"id": f"locus{t_idx}"}

    vals: list[str] = RE_LOCUS_INFO_DELIM.split(val)
    for part in vals:
        assign_parts = RE_LOCUS_INFO_ASSIGN.split(part)
        if len(assign_parts) != 2:
            raise err

        k = assign_parts[0].lower()
        v = assign_parts[1].strip()

        # validate k/v
        if k not in ("id", "motif"):  # key needs to be known
            raise err
        if not v:  # we need to have a non-blank value for the key
            raise LocusValidationError(
                f"BED catalog format error: cannot have empty value in last column on line {t_idx} for key {k}", hint
            )

        # post-processing of motif
        if k == "motif":
            v = v.upper()

        res[k] = v

    # to be valid, we need to have both a non-blank locus ID and a non-blank motif.
    #  - the former is guaranteed by the combination of starting with a default value and checking value blankness above
    #  - we need to check for the latter in case a motif did not get specified at all
    if res.get("motif"):
        return res

    raise err


def load_loci(params: CallParams, locus_queue: Queue, logger: Logger) -> tuple[int, set[str]]:
    @cache  # Cache get_n_alleles calls for contigs
    def _get_contig_n_alleles(ctg: str):
        return params.ploidy_config.n_of(ctg)

    load_start_time = time.perf_counter()

    with open(params.loci_file) as tfh:
        n_loci: int = sum(1 for line in tfh if not line.startswith("#"))

    # Add all loci from the BED file to the queue, allowing each job
    # to pull from the queue as it becomes freed up to do so.
    num_loci: int = 0
    # Keep track of all contigs we are processing to speed up downstream Mendelian inheritance analysis.
    contig_set: set[str] = set()
    last_contig: str | None = None
    last_none_append_n_loci: int = 0

    # We will pass blocks of loci for the queue, so that we can fetch all overlapping segments at once for a block and
    # re-use them via an interval tree during the calling process
    max_block_size: int = min(100, max(n_loci // params.processes, 1))
    max_block_inter_read_dist: int = 20000
    current_block: list[STRkitLocus] = []
    current_block_left: int = 9999999999999999
    current_block_right: int = -1

    def _put_current_block():
        nonlocal current_block_left
        nonlocal current_block_right

        locus_queue.put((current_block_left, current_block_right, tuple(current_block)))
        current_block.clear()
        current_block_left = 9999999999999999
        current_block_right = -1

    for t_idx, t in enumerate(parse_loci_bed(params.loci_file), 1):
        contig = t[0]

        n_alleles: int | None = _get_contig_n_alleles(contig)
        if (
            n_alleles is None  # Sex chromosome, but we don't have a specified sex chromosome karyotype
            or n_alleles == 0  # Don't have this chromosome, e.g., Y chromosome for an XX individual
        ):
            last_contig = contig
            continue  # --> skip the locus

        # For each contig, add None breaks to make the jobs terminate. Then, as long as we have entries in the locus
        # queue, we can get contig-chunks to order and write to disk instead of keeping everything in RAM.
        if last_contig != contig:
            if last_contig is not None:
                if current_block:
                    _put_current_block()

                for _ in range(params.processes):
                    locus_queue.put(None)

                last_none_append_n_loci = num_loci

            contig_set.add(contig)
            last_contig = contig

        # Extract / parse values from the catalog BED file and validate the locus
        try:
            last = parse_last_column(t_idx, t[-1])
            locus = STRkitLocus(
                t_idx, last["id"], contig, int(t[1]), int(t[2]), last["motif"], n_alleles, flank_size=params.flank_size
            )
            validate_locus(locus)
        except LocusValidationError as e:
            e.log_error(logger)
            exit(1)

        if len(current_block) == max_block_size or (
            current_block and locus.left_coord - current_block[-1].right_coord > max_block_inter_read_dist
        ):
            _put_current_block()

        current_block.append(locus)
        current_block_left = min(current_block_left, locus.left_flank_coord)
        current_block_right = max(current_block_right, locus.right_flank_coord)

        num_loci += 1

    if current_block:
        _put_current_block()

    del last_contig  # more as a marker for when we're finished with this, for later refactoring

    if num_loci > last_none_append_n_loci:  # have progressed since the last None-append
        # At the end of the queue, add a None value (* the # of processes).
        # When a job encounters a None value, it will terminate.
        for _ in range(params.processes):
            locus_queue.put(None)

    del last_none_append_n_loci  # more as a marker for when we're finished with this, for later refactoring

    logger.info(f"Loaded {num_loci} loci in {(time.perf_counter() - load_start_time):.2f}s")
    del load_start_time

    return num_loci, contig_set
