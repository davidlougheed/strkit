from __future__ import annotations

import functools
import itertools
import logging
import multiprocessing as mp
import multiprocessing.managers as mmg
import time

import numpy as np
import pysam
import operator
import queue
import threading

from collections import Counter
from collections.abc import Sequence
from datetime import datetime
from pysam import AlignmentFile, FastaFile
from sklearn.cluster import AgglomerativeClustering

from numpy.typing import NDArray
from typing import Iterable, Literal, Optional, Union

from strkit.call.allele import CallDict, call_alleles
from strkit.utils import apply_or_none

from .align_matrix import match_score
from .cigar import get_aligned_pair_matches
from .consensus import consensus_seq
from .params import CallParams
from .realign import realign_read
from .repeats import get_repeat_count, get_ref_repeat_count
from .snvs import (
    SNV_OUT_OF_RANGE_CHAR,
    get_candidate_snvs,
    get_read_snvs,
    calculate_useful_snvs,
    call_and_filter_useful_snvs,
)
from .types import ReadDict, ReadDictExtra, CandidateSNV
from .utils import cat_strs, find_pair_by_ref_pos, normalize_contig, round_to_base_pos, get_new_seed


__all__ = [
    "call_locus",
]


# TODO: Parameterize
CALL_WARN_TIME = 3  # seconds

realign_timeout = 5
min_read_score = 0.9  # TODO: parametrize
base_wildcard_threshold = 3

roughly_equiv_stdev_dist = 1

force_realign = False
many_realigns_threshold = 2

significant_clip_threshold = 100
significant_clip_snv_take_in = 250


# property getters & other partials
cn_getter = operator.itemgetter("cn")
weight_getter = operator.itemgetter("w")


@functools.cache
def _mask_low_q_base(base_and_qual: tuple[str, int]) -> str:
    return base_and_qual[0] if base_and_qual[1] > base_wildcard_threshold else "X"


def calculate_seq_with_wildcards(qs: str, quals: Optional[list[int]]) -> str:
    if quals is None:
        return qs  # No quality information, so don't do anything
    return cat_strs(map(_mask_low_q_base, zip(qs, quals)))


def get_read_coords_from_matched_pairs(
    left_flank_coord: int,
    left_coord: int,
    right_coord: int,
    right_flank_coord: int,
    motif: str,
    motif_size: int,
    query_seq: str,
    q_coords: list[int],
    r_coords: list[int],
    # matched_pairs: list[tuple[int, int]],
) -> tuple[int, int, int, int]:
    left_flank_end = -1
    right_flank_start = -1
    right_flank_end = -1

    last_idx = -1

    # Skip gaps on either side to find mapped flank indices

    # Binary search for left flank start -------------------------------------------------------------------------------
    lhs, found = find_pair_by_ref_pos(r_coords, left_flank_coord)

    # lhs now contains the index for the closest starting coordinate to left_flank_coord

    if not found and (lhs == 0 or lhs == len(r_coords)):
        # Completely out of bounds; either right at the start or inserting after the end
        return -1, -1, -1, -1

    if not found:
        # Choose pair to the left of where we'd insert the pair to maintain sorted order, since we want the closest
        # starting coordinate to left_flank_coord which gives us enough flanking material.
        lhs -= 1

    left_flank_start = q_coords[lhs]
    # ------------------------------------------------------------------------------------------------------------------

    # Binary search for left flank end ---------------------------------------------------------------------------------
    # TODO
    # ------------------------------------------------------------------------------------------------------------------

    for i in range(lhs + 1, len(q_coords)):
        query_coord = q_coords[i]
        ref_coord = r_coords[i]

        # Skip gaps on either side to find mapped flank indices
        if ref_coord < left_coord:
            # Coordinate here is exclusive - we don't want to include a gap between the flanking region and
            # the STR; if we include the left-most base of the STR, we will have a giant flanking region which
            # will include part of the tandem repeat itself.
            left_flank_end = query_coord + 1  # Add 1 to make it exclusive
        elif ref_coord >= right_coord and (
                # Reached end of TR region and haven't set end of TR region yet, or there was an indel with the motif
                # in it right after we finished due to a subtle mis-alignment - this can be seen in the HTT alignments
                # in bc1018
                # TODO: do the same thing for the left side
                right_flank_start == -1 or
                (query_coord - last_idx >= motif_size and (ref_coord - right_coord <= motif_size * 2) and
                 query_seq[last_idx:query_coord].count(motif) / ((query_coord - last_idx) / motif_size) >= 0.5)
        ):
            right_flank_start = query_coord
        elif ref_coord >= right_flank_coord:
            right_flank_end = query_coord
            break

        last_idx = query_coord

    return left_flank_start, left_flank_end, right_flank_start, right_flank_end


def get_overlapping_segments_and_related_data(
    bfs: tuple[pysam.AlignmentFile, ...],
    read_contig: str,
    left_flank_coord: int,
    right_flank_coord: int,
    logger_: logging.Logger,
    locus_log_str: str,
) -> tuple[list[pysam.AlignedSegment], list[int], dict[str, int], int, int]:

    left_most_coord: int = 999999999999999
    right_most_coord: int = 0

    overlapping_segments: list[pysam.AlignedSegment] = []
    seen_reads: set[str] = set()
    read_lengths: list[int] = []

    chimeric_read_status: dict[str, int] = {}

    segment: pysam.AlignedSegment

    for segment in itertools.chain.from_iterable(
        map(lambda bf: bf.fetch(read_contig, left_flank_coord, right_flank_coord), bfs)
    ):
        rn = segment.query_name

        if rn is None:  # Skip reads with no name
            logger_.debug(f"{locus_log_str} - skipping entry for read with no name")
            continue

        supp = segment.flag & 2048

        # If we have two overlapping alignments for the same read, we have a chimeric read within the TR
        # (so probably a large expansion...)

        chimeric_read_status[rn] = chimeric_read_status.get(rn, 0) | (2 if supp else 1)

        if supp:  # Skip supplemental alignments
            logger_.debug(f"{locus_log_str} - skipping entry for read {rn} (supplemental)")
            continue

        if rn in seen_reads:
            logger_.debug(f"{locus_log_str} - skipping entry for read {rn} (already seen)")
            continue

        if segment.query_sequence is None:
            logger_.debug(f"{locus_log_str} - skipping entry for read {rn} (no aligned segment)")
            continue

        if segment.reference_end is None:
            logger_.debug(f"{locus_log_str} - skipping entry for read {rn} (reference_end is None, unmapped?)")
            continue

        seen_reads.add(rn)
        overlapping_segments.append(segment)
        read_lengths.append(segment.query_alignment_length)

        left_most_coord = min(left_most_coord, segment.reference_start)
        right_most_coord = max(right_most_coord, segment.reference_end)

    return overlapping_segments, read_lengths, chimeric_read_status, left_most_coord, right_most_coord


def calculate_read_distance(
    n_reads: int,
    read_dict_items: Sequence[tuple[str, ReadDict]],
    pure_snv_peak_assignment: bool,
    n_useful_snvs: int,
    relative_cn_distance_weight_scaling_few: float = 0.2,
    relative_cn_distance_weight_scaling_many: float = 0.1,
    many_snvs_quantity: int = 3,
) -> NDArray[np.float_]:
    """
    Calculate pairwise distance for all reads using either SNVs ONLY or a mixture of SNVs and copy number.
    :param n_reads: Number of reads.
    :param read_dict_items: Itemized read dictionary entries: (read name, read data)
    :param pure_snv_peak_assignment: Whether to use just SNVs for peak assignment
    :param n_useful_snvs: Number of useful SNVs; length of snvu lists.
    :param relative_cn_distance_weight_scaling_few: How much to weight 1-difference in CN vs SNVs with few SNVs
            available (indels are more erroneous in CCS)
    :param relative_cn_distance_weight_scaling_many: How much to weight 1-difference in CN vs SNVs with many SNVs
            available (indels are more erroneous in CCS)
    :param many_snvs_quantity: How many SNVs constitutes "many" for the different weights.
    :return: The distance matrix.
    """

    useful_snvs_range: tuple[int, ...] = tuple(range(n_useful_snvs))

    # Initialize a distance matrix for all reads
    distance_matrix = np.zeros((n_reads, n_reads), dtype=np.float_)

    # Loop through and compare all vs. all reads. We can skip a few indices since the distance will be symmetrical.
    for i in range(n_reads - 1):
        r1 = read_dict_items[i][1]
        r1_snv_u = r1["snvu"]

        r1_out_of_range: set[int] = set(filter(lambda y: r1_snv_u[y] == SNV_OUT_OF_RANGE_CHAR, useful_snvs_range))

        for j in range(i + 1, n_reads):
            r2 = read_dict_items[j][1]
            r2_snv_u = r2["snvu"]

            n_not_equal: int = 0
            n_comparable: int = 0

            for z in useful_snvs_range:
                if z in r1_out_of_range:
                    continue
                r2_b = r2_snv_u[z]
                if r2_b == SNV_OUT_OF_RANGE_CHAR:
                    continue
                r1_b = r1_snv_u[z]
                if r1_b != r2_b:
                    n_not_equal += 1
                n_comparable += 1

            d: float = float(n_not_equal)
            if not pure_snv_peak_assignment:  # Add in copy number distance
                d += abs(r1["cn"] - r2["cn"]) * (
                    relative_cn_distance_weight_scaling_many if n_comparable >= many_snvs_quantity
                    else relative_cn_distance_weight_scaling_few)

            distance_matrix[i, j] = d
            distance_matrix[j, i] = d

    return distance_matrix


def process_read_snvs_for_locus(
    contig: str,
    left_coord_adj: int,
    right_coord_adj: int,
    left_most_coord: int,
    right_most_coord: int,
    ref: FastaFile,
    read_dict_items: tuple[tuple[str, ReadDict], ...],
    read_dict_extra: dict[str, ReadDictExtra],
    read_pairs: dict[str, tuple[list[int], list[int]]],
    candidate_snvs_dict: dict[int, CandidateSNV],
    only_known_snvs: bool,
    logger_: logging.Logger,
    locus_log_str: str,
) -> set[int]:
    # Loop through a second time if we are using SNVs. We do a second loop rather than just using the first loop
    # in order to have collected the edges of the reference sequence we can cache for faster SNV calculation.

    # Mutates: read_dict_extra

    locus_snvs: set[int] = set()

    ref_cache = ref.fetch(contig, left_most_coord, right_most_coord + 1).upper()

    for rn, read in read_dict_items:
        #   --> RE-ENABLE FOR DE NOVO SNV FINDER <--
        scl = read_dict_extra[rn]["sig_clip_left"]
        scr = read_dict_extra[rn]["sig_clip_right"]
        if scl or scr:
            logger_.debug(
                f"{locus_log_str} - {rn} has significant clipping; trimming pairs by "
                f"{significant_clip_snv_take_in} bp per side for SNV-finding")

        # snv_pairs = read_pairs[rn]
        q_coords, r_coords = read_pairs[rn]

        if len(q_coords) < (twox_takein := significant_clip_snv_take_in * 2):
            logger_.warning(f"{locus_log_str} - skipping SNV calculation for '{rn}' (<{twox_takein} pairs)")
            continue

        if scl:
            # snv_pairs = snv_pairs[significant_clip_snv_take_in:]
            q_coords = q_coords[significant_clip_snv_take_in:]
            r_coords = r_coords[significant_clip_snv_take_in:]
        if scr:
            # snv_pairs = snv_pairs[:-1 * significant_clip_snv_take_in]
            q_coords = q_coords[:-1 * significant_clip_snv_take_in]
            r_coords = r_coords[:-1 * significant_clip_snv_take_in]

        snvs = get_read_snvs(
            read_dict_extra[rn]["_qs"],
            ref_cache,
            q_coords,
            r_coords,
            left_most_coord,
            left_coord_adj,
            right_coord_adj,
            # below: magic values for skipping false positives / weird 'SNVs' that aren't helpful
            contiguous_threshold=5,
            max_snv_group_size=5,
            too_many_snvs_threshold=20,
            entropy_flank_size=10,
            entropy_threshold=1.8,
        )
        #   --> RE-ENABLE FOR dbSNP-BASED SNV FINDER <--
        # snvs = get_read_snvs_dbsnp(
        #     candidate_snvs_dict_items_flat,
        #     read_dict_extra[rn]["_qs"],
        #     read_pairs[rn],
        #     left_coord_adj,
        #     right_coord_adj,
        # )
        locus_snvs.update(filter(lambda p: (not only_known_snvs) or p in candidate_snvs_dict, snvs.keys()))
        read_dict_extra[rn]["snv"] = snvs

    return locus_snvs


def call_alleles_with_haplotags(
    params: CallParams,
    haplotags: list[str],
    ps_id: int,
    read_dict_items: tuple[tuple[str, ReadDict], ...],  # We could derive this again, but we already have before...
    rng: np.random.Generator,
    logger_: logging.Logger,
    locus_log_str: str,
) -> Optional[dict]:
    n_alleles: int = len(haplotags)

    hp_reads: list[tuple[ReadDict, ...]] = []
    cns: Union[list[list[int]], list[list[float]]] = []
    c_ws: list[Union[NDArray[np.int_], NDArray[np.float_]]] = []

    for hi, hp in enumerate(haplotags):
        # Find reads for cluster
        crs: tuple[ReadDict, ...] = tuple(
            r for i, (_, r) in enumerate(read_dict_items)
            if r.get("hp") == hp and r.get("ps") == ps_id)

        # Calculate copy number set
        cns.append(list(map(cn_getter, crs)))

        # Calculate weights array
        ws = np.fromiter(map(weight_getter, crs), dtype=np.float_)
        c_ws.append(ws / np.sum(ws))

        hp_reads.append(crs)

    cdd: list[CallDict] = []

    for hi, hp in enumerate(haplotags):
        cc: Optional[CallDict] = call_alleles(
            cns[hi], (),  # Don't bother separating by strand for now...
            c_ws[hi], (),
            params=params,
            min_reads=params.min_allele_reads,  # Calling alleles separately, so set min_reads=min_allele_reads
            n_alleles=1,  # Calling alleles separately: they were pre-separated by agglom. clustering
            separate_strands=False,
            read_bias_corr_min=0,  # separate_strands is false, so this is ignored
            gm_filter_factor=1,  # n_alleles=1, so this is ignored
            seed=get_new_seed(rng),
            logger_=logger_,
            debug_str=f"{locus_log_str} a{hi}"
        )

        if cc is None:  # Early escape
            return None

        # TODO: set peak weight [0] to the sum of read weights - we normalize this later, but this way
        #  call dicts with more reads will GET MORE WEIGHT! as it should be, instead of 50/50 for the peak.

        cdd.append(cc)

    # TODO: Multi-allele phasing across STRs

    for i in range(len(haplotags)):  # Cluster indices now referring to ordered ones
        for rd in hp_reads[i]:
            rd["p"] = i

    peak_weights_pre_adj = np.concatenate(tuple(cc["peak_weights"] for cc in cdd), axis=0)

    # All call_datas are truth-y; all arrays should be ordered by peak_order
    call_data = {
        "call": np.concatenate(tuple(cc["call"] for cc in cdd), axis=0),
        "call_95_cis": np.concatenate(tuple(cc["call_95_cis"] for cc in cdd), axis=0),
        "call_99_cis": np.concatenate(tuple(cc["call_99_cis"] for cc in cdd), axis=0),
        "peaks": np.concatenate(tuple(cc["peaks"] for cc in cdd), axis=None),

        # TODO: Readjust peak weights when combining or don't include
        # Make peak weights sum to 1
        "peak_weights": peak_weights_pre_adj / np.sum(peak_weights_pre_adj),

        "peak_stdevs": np.concatenate(tuple(cc["peak_stdevs"] for cc in cdd), axis=0),
        "modal_n_peaks": n_alleles,  # n. of alleles = n. of peaks always -- if we phased using SNVs

        "ps": ps_id,
    }

    return call_data


def call_alleles_with_incorporated_snvs(
    n_alleles: int,
    params: CallParams,
    read_dict: dict[str, ReadDict],
    read_dict_items: tuple[tuple[str, ReadDict], ...],  # We could derive this again, but we already have before...
    read_dict_extra: dict[str, dict],
    n_reads_in_dict: int,  # We could derive this again, but we already have before...
    useful_snvs: list[tuple[int, int]],
    candidate_snvs_dict: dict[int, CandidateSNV],
    # ---
    phase_set_lock: threading.Lock,
    phase_set_counter: mmg.ValueProxy,
    snv_genotype_update_lock: threading.Lock,
    snv_genotype_cache: mmg.DictProxy,
    # ---
    rng: np.random.Generator,
    logger_: logging.Logger,
    locus_log_str: str,
) -> tuple[Literal["dist", "snv", "snv+dist", "single"], Optional[tuple[dict, list[dict]]]]:
    assign_method: Literal["dist", "snv", "snv+dist", "single"] = "dist"

    # TODO: parametrize min 'enough to do pure SNV haplotyping' thresholds

    n_useful_snvs: int = len(useful_snvs)
    read_dict_items_with_many_snvs: list[tuple[str, ReadDict]] = []
    read_dict_items_with_at_least_one_snv: list[tuple[str, ReadDict]] = []
    read_dict_items_with_no_snvs: list[tuple[str, ReadDict]] = []

    print_snvs = False

    # TODO: Check that we don't have disjoint groups - how would we actually ascertain this?
    #  - each read as a set; merge sets when there is overlap, make sure we have 1 set at the end

    for read_item in read_dict_items:
        rn, read = read_item
        snv_bases = read_dict_extra[rn].get("snv_bases")

        if snv_bases is None:
            read_dict_items_with_no_snvs.append(read_item)
            continue

        read_useful_snv_bases = tuple(snv_bases[bi] for bi, _pos in useful_snvs)
        n_non_blank_read_useful_snv_bases = len(
            tuple(filter(lambda bb: bb != SNV_OUT_OF_RANGE_CHAR, read_useful_snv_bases))
        )

        if n_non_blank_read_useful_snv_bases:  # TODO: parametrize
            read_dict_items_with_at_least_one_snv.append(read_item)
            if n_non_blank_read_useful_snv_bases >= 2:  # TODO: parametrize
                read_dict_items_with_many_snvs.append(read_item)
        else:
            read_dict_items_with_no_snvs.append(read_item)

        read["snvu"] = read_useful_snv_bases  # Store read-level 'useful' SNVs

        if print_snvs:
            print(rn, f"\t{read['cn']:.0f}", "\t", cat_strs(read_useful_snv_bases), n_non_blank_read_useful_snv_bases)

    n_reads_with_many_snvs: int = len(read_dict_items_with_many_snvs)
    n_reads_with_at_least_one_snv: int = len(read_dict_items_with_at_least_one_snv)
    n_reads_with_no_snvs: int = len(read_dict_items_with_no_snvs)

    # Use purely SNVs if all reads which won't get discarded have many SNVs
    pure_snv_peak_assignment: bool = n_reads_with_many_snvs + n_reads_with_no_snvs == n_reads_in_dict

    # TODO: parametrize: how many reads with SNV information
    min_snv_incorporation_read_portion = 0.8  # At least 80% with 1+ SNV called
    min_snv_incorporation_read_abs = 16  # Or at least 16 reads with 1+ SNV called

    can_incorporate_snvs: bool = pure_snv_peak_assignment or (
        n_reads_with_at_least_one_snv >=
        min(n_reads_in_dict * min_snv_incorporation_read_portion, min_snv_incorporation_read_abs)
    )

    if not can_incorporate_snvs:
        # TODO: How to use partial data?
        return assign_method, None

    if n_reads_with_no_snvs:
        logger_.debug(f"{locus_log_str} - will discard {n_reads_with_no_snvs} reads with no SNV data")

    # Otherwise, we can use the SNV data --------------------------------------

    if pure_snv_peak_assignment:
        # We have enough SNVs in ALL reads, so we can phase purely based on SNVs
        logger_.debug(f"{locus_log_str} - haplotyping purely using SNVs ({n_useful_snvs=}, {n_reads_with_many_snvs=})")
        assign_method = "snv"
    else:
        # We have enough SNVs in lots of reads, so we can phase using a combined metric
        logger_.debug(
            f"{locus_log_str} - haplotyping using combined STR-SNV metric ("
            f"{n_useful_snvs=}, {n_reads_with_at_least_one_snv=}, {n_reads_in_dict=})")
        # TODO: Handle reads we didn't have SNVs for by retroactively assigning to groups
        assign_method = "snv+dist"

    # Calculate pairwise distance for all reads using either SNVs ONLY or
    # a mixture of SNVs and copy number:
    # dm = calculate_read_distance(n_reads_in_dict, read_dict_items, pure_snv_peak_assignment, n_useful_snvs)
    dm = calculate_read_distance(
        n_reads_with_at_least_one_snv, read_dict_items_with_at_least_one_snv, pure_snv_peak_assignment, n_useful_snvs)

    # Cluster reads together using the distance matrix, which incorporates
    # SNV and possibly copy number information.
    c = AgglomerativeClustering(n_clusters=n_alleles, metric="precomputed", linkage="average").fit(dm)

    # noinspection PyUnresolvedReferences
    cluster_labels = c.labels_
    cluster_indices = tuple(range(n_alleles))

    cluster_reads: list[tuple[ReadDict, ...]] = []
    cns: Union[list[list[int]], list[list[float]]] = []
    c_ws: list[Union[NDArray[np.int_], NDArray[np.float_]]] = []

    for ci in cluster_indices:
        # Find reads for cluster
        crs: tuple[ReadDict, ...] = tuple(
            r for i, (_, r) in enumerate(read_dict_items_with_at_least_one_snv)
            if cluster_labels[i] == ci
        )

        # Calculate copy number set
        cns.append(list(map(cn_getter, crs)))

        # Calculate weights array
        ws = np.fromiter(map(weight_getter, crs), dtype=np.float_)
        c_ws.append(ws / np.sum(ws))

        cluster_reads.append(crs)

    logger_.debug(f"{locus_log_str} - using {assign_method=}, got {cns=}")

    cdd: list[CallDict] = []

    for ci in cluster_indices:
        cc: Optional[CallDict] = call_alleles(
            cns[ci], (),  # Don't bother separating by strand for now...
            c_ws[ci], (),
            params,
            min_reads=params.min_allele_reads,  # Calling alleles separately, so set min_reads=min_allele_reads
            n_alleles=1,  # Calling alleles separately: they were pre-separated by agglom. clustering
            separate_strands=False,
            read_bias_corr_min=0,  # separate_strands is false, so this is ignored
            gm_filter_factor=1,  # n_alleles=1, so this is ignored
            seed=get_new_seed(rng),
            logger_=logger_,
            debug_str=f"{locus_log_str} a{ci}"
        )

        if cc is None:  # Early escape
            # One of the calls could not be made... what to do?
            # TODO: !!!!
            #  For now, revert to dist
            return "dist", None

        # TODO: set peak weight [0] to the sum of read weights - we normalize this later, but this way
        #  call dicts with more reads will GET MORE WEIGHT! as it should be, instead of 50/50 for the peak.

        cdd.append(cc)

    # We called these as single-allele (1 peak) loci as a sort of hack, so the return "call" key
    # is an array of length 1.

    cdd_sort_order_determiner = np.fromiter(
        map(lambda x: (x["peaks"][0], x["call_95_cis"][0][0]), cdd),
        dtype=[("p", np.float_), ("i", np.float_ if params.fractional else np.int_)])
    # To reorder call arrays in least-to-greatest by raw peak mean, and then by 95% CI left boundary:
    peak_order: NDArray[np.int_] = np.argsort(cdd_sort_order_determiner, order=("p", "i"))

    # Re-order the call data dictionary, now that we've established an ordering
    # noinspection PyTypeChecker
    cdd_ordered = [cdd[i] for i in peak_order]
    # noinspection PyTypeChecker
    cluster_reads_ordered = [cluster_reads[i] for i in peak_order]

    # Assign peak labels now that we've reordered
    for i in cluster_indices:  # Cluster indices now referring to ordered ones
        for rd in cluster_reads_ordered[i]:
            rd["p"] = i

    # Call useful SNVs to add to the final return dictionary ----------------------------------------------------------
    #  - This method needs the read_dict[p] value, so we need to run this after initial peaks have been calculated!
    called_useful_snvs: list[dict] = call_and_filter_useful_snvs(
        n_alleles, read_dict, useful_snvs, candidate_snvs_dict, locus_log_str, logger_)

    if not called_useful_snvs:  # No useful SNVs left, so revert to "dist" assignment method
        return "dist", None
    # ------------------------------------------------------------------------------------------------------------------

    # We may need to re-order (flip) calls based on SNVs. Check each SNV to see if it's in the SNV genotype/phase-set
    # dictionary; otherwise, assign a phase set to all reads which have been used for peak calling here.

    call_phase_set: Optional[int] = None

    found_snvs: list[str] = []
    snv_pss: list[int] = []
    should_flip: list[bool] = []

    snv_genotype_update_lock.acquire()
    try:
        for snv in called_useful_snvs:
            if snv["id"] in snv_genotype_cache:
                t_snv_genotype, snv_ps = snv_genotype_cache[snv["id"]]
                found_snvs.append(snv["id"])
                snv_pss.append(snv_ps)
                should_flip.append(len(t_snv_genotype) > 1 and tuple(t_snv_genotype) == tuple(reversed(snv["call"])))
        if not found_snvs:
            phase_set_lock.acquire()
            call_phase_set = int(phase_set_counter.value)
            phase_set_counter.set(call_phase_set + 1)
            phase_set_lock.release()

            for snv in called_useful_snvs:
                snv_genotype_cache[snv["id"]] = (tuple(snv["call"]), call_phase_set)
    finally:
        snv_genotype_update_lock.release()

    if found_snvs:  # else from the above `if not found_snvs`, but we can release the lock earlier
        # Have found SNVs, should flip/not flip and assign existing phase set
        all_should_flip = all(should_flip)
        flip_consensus = all_should_flip or all(map(operator.not_, should_flip))
        phase_set_consensus_set = tuple(sorted(set(snv_pss)))
        phase_set_consensus = len(phase_set_consensus_set) == 1
        if flip_consensus and phase_set_consensus:
            call_phase_set = phase_set_consensus_set[0]
            if all_should_flip:
                for r in read_dict.values():
                    if (p := r.get("p")) is not None:
                        # Flip peak - in SNV mode, we're not polyploid, so we can just do int(not bool(p))
                        r["p"] = int(not bool(p))
                        r["ps"] = call_phase_set
                # Then, reverse ordered call data list
                cdd_ordered.reverse()
            else:
                for r in read_dict.values():
                    if r.get("p") is not None:
                        # We're good as-is, so assign the phase set
                        r["ps"] = call_phase_set
        else:
            # TODO: what to do? can we recover something here? yes - we could join phase sets post-facto
            #  - sort snv_pss so we keep the smaller one

            logger_.warning(
                f"{locus_log_str} - no consensus for SNV phase set/should-flip: got {snv_pss}/{should_flip}; "
                f"skipping phasing")

    peak_weights_pre_adj = np.concatenate(tuple(cc["peak_weights"] for cc in cdd), axis=0)

    # All call_datas are truth-y; all arrays should be ordered by peak_order
    call_data = {
        "call": np.concatenate(tuple(cc["call"] for cc in cdd_ordered), axis=0),
        "call_95_cis": np.concatenate(tuple(cc["call_95_cis"] for cc in cdd_ordered), axis=0),
        "call_99_cis": np.concatenate(tuple(cc["call_99_cis"] for cc in cdd_ordered), axis=0),
        "peaks": np.concatenate(tuple(cc["peaks"] for cc in cdd_ordered), axis=None),

        # TODO: Readjust peak weights when combining or don't include
        # Make peak weights sum to 1
        "peak_weights": peak_weights_pre_adj / np.sum(peak_weights_pre_adj),

        "peak_stdevs": np.concatenate(tuple(cc["peak_stdevs"] for cc in cdd_ordered), axis=0),
        "modal_n_peaks": n_alleles,  # n. of alleles = n. of peaks always -- if we phased using SNVs

        **({"ps": call_phase_set} if call_phase_set is not None else {}),
    }

    # Return the SNV/SNV+dist method + our called peaks, and update the phase set information.

    return assign_method, (call_data, called_useful_snvs)


def debug_log_flanking_seq(logger_: logging.Logger, locus_log_str: str, rn: str, realigned: bool):
    logger_.debug(
        f"{locus_log_str} - skipping read {rn}: could not get sufficient flanking sequence"
        f"{' (post-realignment)' if realigned else ''}")


def call_locus(
    t_idx: int,
    t: tuple,
    n_alleles: int,
    bfs: tuple[AlignmentFile, ...],
    ref: FastaFile,
    params: CallParams,
    # ---
    phase_set_lock: threading.Lock,
    phase_set_counter: mmg.ValueProxy,
    phase_set_remap: mmg.DictProxy,
    snv_genotype_update_lock: threading.Lock,
    snv_genotype_cache: mmg.DictProxy,
    # ---
    seed: int,
    logger_: logging.Logger,
    # ---
    snv_vcf_file: Optional[pysam.VariantFile] = None,
    snv_vcf_contigs: tuple[str, ...] = (),
    snv_vcf_file_format: Literal["chr", "num", "acc", ""] = "",
    # ---
    read_file_has_chr: bool = True,
    ref_file_has_chr: bool = True,
) -> Optional[dict]:
    call_timer = datetime.now()

    # params de-structuring ------------
    consensus = params.consensus
    count_kmers = params.count_kmers
    flank_size = params.flank_size
    fractional = params.fractional
    log_level = params.log_level
    realign = params.realign
    respect_ref = params.respect_ref
    sample_id = params.sample_id
    # ----------------------------------

    rng = np.random.default_rng(seed=seed)

    contig: str = t[0]
    read_contig = normalize_contig(contig, read_file_has_chr)
    ref_contig = normalize_contig(contig, ref_file_has_chr)

    motif: str = t[-1]
    motif_size = len(motif)

    left_coord = int(t[1])
    right_coord = int(t[2])

    left_flank_coord = left_coord - flank_size
    right_flank_coord = right_coord + flank_size

    ref_left_flank_seq: str = ""
    ref_right_flank_seq: str = ""
    ref_right_flank_seq_plus_1: str = ""
    ref_seq: str = ""
    raised: bool = False

    # Currently, only support diploid use of SNVs. There's not much of a point with haploid loci,
    # and polyploidy is hard.
    # should_incorporate_snvs: bool = incorporate_snvs and n_alleles == 2
    should_incorporate_snvs: bool = snv_vcf_file is not None and n_alleles == 2
    only_known_snvs: bool = True  # TODO: parametrize

    try:
        ref_left_flank_seq = ref.fetch(ref_contig, left_flank_coord, left_coord)
        ref_right_flank_seq_plus_1 = ref.fetch(ref_contig, right_coord, right_flank_coord + 1)
        # ref_right_flank_seq = ref.fetch(ref_contig, right_coord, right_flank_coord)
        ref_right_flank_seq = ref_right_flank_seq_plus_1[:-1]
        ref_seq = ref.fetch(ref_contig, left_coord, right_coord)
    except IndexError:
        logger_.warning(
            f"Coordinates out of range in provided reference FASTA for region {ref_contig} with flank size "
            f"{params.flank_size}: [{left_flank_coord}, {right_flank_coord}] (skipping locus {t_idx})")
        raised = True
    except ValueError:
        logger_.error(f"Invalid region '{ref_contig}' for provided reference FASTA (skipping locus {t_idx})")
        raised = True

    if len(ref_left_flank_seq) < params.flank_size or len(ref_right_flank_seq) < params.flank_size:
        if not raised:  # flank sequence too small for another reason
            logger_.warning(f"Reference flank size too small for locus {t_idx} (skipping)")
            return None

    if raised:
        return None

    # String representation of locus for logging purposes
    locus_log_str: str = \
        f"{sample_id or ''}{' ' if sample_id else ''}locus {t_idx}: {contig}:{left_coord}-{right_coord}"

    # Get reference repeat count by our method, so we can calculate offsets from reference
    ref_cn: Union[int, float]
    (ref_cn, _), l_offset, r_offset = get_ref_repeat_count(
        round(len(ref_seq) / motif_size),  # Initial estimate of copy number based on coordinates + motif size
        ref_seq,
        ref_left_flank_seq,
        ref_right_flank_seq,
        motif,
        ref_size=right_coord-left_coord,  # reference size, in terms of coordinates (not TRF-recorded size)
        fractional=fractional,
        respect_coords=respect_ref,
    )

    # If our reference repeat count getter has altered the TR boundaries a bit (which is done to allow for
    # more spaces in which an indel could end up), adjust our coordinates to match.
    # Currently, contractions of the TR region are ignored.
    left_coord_adj = left_coord if respect_ref else left_coord - max(0, l_offset)
    right_coord_adj = right_coord if respect_ref else right_coord + max(0, r_offset)

    # Find the initial set of overlapping aligned segments with associated read lengths + whether we have in-locus
    # chimera reads (i.e., reads which aligned twice with different soft-clipping, likely due to a large indel.) -------

    overlapping_segments: list[pysam.AlignedSegment]
    read_lengths: list[int]
    chimeric_read_status: dict[str, int]

    overlapping_segments, read_lengths, chimeric_read_status, left_most_coord, right_most_coord = \
        get_overlapping_segments_and_related_data(
            bfs, read_contig, left_flank_coord, right_flank_coord, logger_, locus_log_str)

    n_overlapping_reads = len(overlapping_segments)

    sorted_read_lengths = np.sort(read_lengths)

    # Find candidate SNVs, if we're using SNV data

    candidate_snvs_dict: dict[int, CandidateSNV] = {}  # Lookup dictionary for candidate SNVs by position
    if n_overlapping_reads and should_incorporate_snvs and snv_vcf_file:
        # ^^ n_overlapping_reads check since otherwise we will have invalid left/right_most_coord
        candidate_snvs_dict = get_candidate_snvs(
            snv_vcf_file, snv_vcf_contigs, snv_vcf_file_format, contig, left_most_coord, right_most_coord)

    # Build the read dictionary with segment information, copy number, weight, & more. ---------------------------------

    read_dict: dict[str, ReadDict] = {}
    read_dict_extra: dict[str, ReadDictExtra] = {}
    realign_count: int = 0  # Number of realigned reads

    # Various aggregators for if we have a phased alignment file:
    haplotagged_reads_count: int = 0  # Number of reads with HP tags
    haplotags: set[str] = set()
    phase_sets: Counter[int] = Counter()

    # Aggregations for additional read-level data
    read_kmers: Counter[str] = Counter()
    read_pairs: dict[str, tuple[list[int], list[int]]] = {}

    # Keep track of left-most and right-most coordinates
    # If SNV-based peak calling is enabled, we can use this to pre-fetch reference data for all reads to reduce the
    # fairly significant overhead involved in reading from the reference genome for each read to identifify SNVs.
    # left_most_coord: int = 99999999999999
    # right_most_coord: int = 0

    segment: pysam.AlignedSegment
    for segment, read_len in zip(overlapping_segments, read_lengths):
        rn: str = segment.query_name  # Know this is not None from overlapping_segments calculation
        segment_start: int = segment.reference_start
        segment_end: int = segment.reference_end  # Optional[int], but if it's here we know it isn't None

        # left_most_coord = min(left_most_coord, segment_start)
        # right_most_coord = max(right_most_coord, segment_end)

        # While .query_sequence is Optional[str], we know (because we skipped all segments with query_sequence is None
        # above) that this is guaranteed to be, in fact, not None.
        qs: str = segment.query_sequence

        fqqs: Optional[list[int]] = segment.query_qualities
        cigar_tuples: list[tuple[int, int]] = segment.cigartuples

        realigned: bool = False
        # pairs = None
        q_coords: Optional[list[int]] = None
        r_coords: Optional[list[int]] = None

        # Soft-clipping in large insertions can result from mapping difficulties.
        # If we have a soft clip which overlaps with our TR region (+ flank), we can try to recover it
        # via realignment with parasail.
        # 4: BAM code for soft clip CIGAR operation
        # TODO: if some alignment is present, use it to reduce realignment overhead?
        #  - use start point + flank*3 or end point - flank*3 or something like that
        if realign and (force_realign or (
            ((c1 := cigar_tuples[0])[0] == 4 and segment_start > left_flank_coord >= segment_start - c1[1]) or
            ((c2 := cigar_tuples[-1])[0] == 4 and segment_end < right_flank_coord <= segment_end + c2[1])
        )):
            # Run the realignment in a separate process, to give us a timeout mechanism.
            # This means we're spawning a second process for this job, just momentarily, beyond the pool size.

            q: mp.Queue = mp.Queue()
            proc = mp.Process(target=realign_read, daemon=False, kwargs=dict(
                # fetch an extra base for the right flank coordinate check later (needs to be >= the exclusive coord)
                ref_seq=f"{ref_left_flank_seq}{ref_seq}{ref_right_flank_seq_plus_1}",  # TODO: plus 1, really?
                query_seq=calculate_seq_with_wildcards(qs, fqqs),
                left_flank_coord=left_flank_coord,
                flank_size=flank_size,
                rn=rn,
                t_idx=t_idx,
                always_realign=force_realign,
                q=q,
                log_level=log_level,
            ))
            proc.start()
            pairs_new = None
            try:
                pairs_new = q.get(timeout=realign_timeout)
                proc.join()
            except queue.Empty:
                logger_.warning(
                    f"{locus_log_str} - experienced timeout while re-aligning read {rn}. Reverting to BAM alignment.")
                proc.terminate()
                time.sleep(0.1)  # wait a little for the process to terminate
            finally:
                wait_count: int = 0
                while proc.is_alive():
                    logger_.warning(f"{locus_log_str} - realign job has still not exited. Waiting 0.5 seconds...")
                    time.sleep(0.5)
                    wait_count += 1
                    if wait_count > 5:
                        logger_.fatal(f"{locus_log_str} - realign job never exited. Terminating...")
                        exit(1)
                proc.close()

            if pairs_new is not None:
                q_coords, r_coords = pairs_new
                realigned = True
                realign_count += 1

        if q_coords is None:
            if left_flank_coord < segment_start or right_flank_coord > segment_end:
                # Cannot find pair for LHS flank start or RHS flank end;
                # early-continue before we load pairs since that step is slow
                debug_log_flanking_seq(logger_, locus_log_str, rn, realigned)
                continue

            q_coords, r_coords = get_aligned_pair_matches(cigar_tuples, 0, segment_start)

        left_flank_start, left_flank_end, right_flank_start, right_flank_end = get_read_coords_from_matched_pairs(
            left_flank_coord,
            left_coord_adj,
            right_coord_adj,
            right_flank_coord,
            motif,
            motif_size,
            query_seq=qs,
            q_coords=q_coords,
            r_coords=r_coords,
        )

        if any(v == -1 for v in (left_flank_start, left_flank_end, right_flank_start, right_flank_end)):
            debug_log_flanking_seq(logger_, locus_log_str, rn, realigned)
            continue

        # we can fit PHRED scores in uint8
        qqs = np.fromiter(fqqs[left_flank_end:right_flank_start], dtype=np.uint8)
        if qqs.shape[0] and (m_qqs := np.mean(qqs)) < (min_avg_phred := params.min_avg_phred):  # TODO: check flank?
            logger_.debug(
                f"{locus_log_str} - skipping read {rn} due to low average base quality ({m_qqs:.2f} < {min_avg_phred})")
            continue

        # -----

        # Truncate to flank_size (plus some leeway for small indels in flanking region) to stop any expansion sequences
        # from accidentally being included in the flanking region; e.g. if the insert gets mapped onto bases outside
        # the definition coordinates.
        # The +10 here won't include any real TR region if the mapping is solid, since the flank coordinates will
        # contain a correctly-sized sequence.

        # TODO: wildcards in flanking region too?
        flank_left_seq: str = qs[left_flank_start:left_flank_end][:flank_size+10]
        flank_right_seq: str = qs[right_flank_start:right_flank_end][-(flank_size+10):]

        tr_len: int = right_flank_start - left_flank_end  # i.e., len(tr_read_seq)
        flank_len: int = len(flank_left_seq) + len(flank_right_seq)
        tr_len_w_flank: int = tr_len + flank_len

        tr_read_seq = qs[left_flank_end:right_flank_start]
        tr_read_seq_wc = calculate_seq_with_wildcards(qs[left_flank_end:right_flank_start], qqs)

        if count_kmers != "none":
            read_kmers.clear()
            for i in range(0, tr_len - motif_size + 1):
                read_kmers.update((tr_read_seq_wc[i:i+motif_size],))

        read_cn, read_cn_score = get_repeat_count(
            start_count=round(tr_len / motif_size),  # Set initial integer copy number based on aligned TR size
            tr_seq=tr_read_seq_wc,
            flank_left_seq=flank_left_seq,
            flank_right_seq=flank_right_seq,
            motif=motif,
            fractional=fractional,
        )

        # TODO: need to rethink this; it should maybe quantify mismatches/indels in the flanking regions
        read_adj_score: float = match_score if tr_len == 0 else read_cn_score / tr_len_w_flank
        if read_adj_score < min_read_score:
            logger_.debug(f"{locus_log_str} - skipping read {rn} (scored {read_adj_score} < {min_read_score})")
            continue

        # When we don't have targeted sequencing, the probability of a read containing the TR region, given that it
        # overlaps the region, is P(read is large enough to contain) * P(  # TODO: complete this..
        partition_idx = np.searchsorted(sorted_read_lengths, tr_len_w_flank, side="right")
        if partition_idx == n_overlapping_reads:  # tr_len_w_flank is longer than the longest read... :(
            # Fatal
            # TODO: Just skip this locus
            logger_.error(
                f"{locus_log_str} - something strange happened; could not find an encompassing read where one should "
                f"be guaranteed. TR length with flank: {tr_len_w_flank}; read lengths: {sorted_read_lengths}")
            exit(1)

        mean_containing_size = read_len if params.targeted else np.mean(sorted_read_lengths[partition_idx:]).item()
        # TODO: re-examine weighting to possibly incorporate chance of drawing read large enough
        read_weight = (mean_containing_size + tr_len_w_flank - 2) / (mean_containing_size - tr_len_w_flank + 1)

        crs_cir = chimeric_read_status[rn] == 3  # Chimera within the TR region, indicating a potential large expansion
        read_dict[rn] = {
            "s": "-" if segment.is_reverse else "+",
            "cn": read_cn,
            "w": read_weight,
            **({"realn": realigned} if realign and realigned else {}),
            **({"chimeric_in_region": crs_cir} if crs_cir else {}),
            **({"kmers": dict(read_kmers)} if count_kmers != "none" else {}),
        }
        read_dict_extra[rn] = {
            "_ref_start": segment_start,
            "_ref_end": segment_end,
            **({"_tr_seq": tr_read_seq} if consensus else {}),
        }

        # Reads can show up more than once - TODO - cache this information across loci

        if params.use_hp:
            tags = dict(segment.get_tags())
            if (hp := tags.get("HP")) is not None and (ps := tags.get("PS")) is not None:
                orig_ps = int(ps)
                ps_remapped: int
                if orig_ps in phase_set_remap:
                    ps_remapped = phase_set_remap[orig_ps]
                else:
                    phase_set_lock.acquire()
                    psc = int(phase_set_counter.value)
                    phase_set_counter.set(psc + 1)
                    phase_set_lock.release()

                    phase_set_remap[orig_ps] = psc
                    ps_remapped = psc

                read_dict[rn]["hp"] = hp
                read_dict[rn]["ps"] = ps_remapped
                haplotags.add(hp)
                haplotagged_reads_count += 1
                phase_sets[ps_remapped] += 1

        if should_incorporate_snvs:
            # Store the segment sequence in the read dict for the next go-around if we've enabled SNV incorporation,
            # in order to pass the query sequence to the get_read_snvs function with the cached ref string.
            read_dict_extra[rn]["_qs"] = qs

            # Observed significant increase in annoying, probably false SNVs near the edges of significantly
            # clipped reads in CCS data. Figure out if we have large clipping for later use here in the SNV finder.
            #   --> RE-ENABLE FOR DE NOVO SNV FINDER <--
            cigar_first_op = segment.cigartuples[0]
            cigar_last_op = segment.cigartuples[-1]

            read_dict_extra[rn]["sig_clip_left"] = (
                cigar_first_op[0] in (4, 5) and cigar_first_op[1] >= significant_clip_threshold)
            read_dict_extra[rn]["sig_clip_right"] = (
                    cigar_last_op[0] in (4, 5) and cigar_last_op[1] >= significant_clip_threshold)

            # Cache aligned pairs, since it takes a lot of time to extract, and we use it for calculate_useful_snvs
            read_pairs[rn] = q_coords, r_coords

    # End of first read loop -------------------------------------------------------------------------------------------

    n_reads_in_dict: int = len(read_dict)

    call_dict_base = {
        "locus_index": t_idx,
        "contig": contig,
        "start": left_coord,
        "end": right_coord,
        **({} if respect_ref else {
            "start_adj": left_coord_adj,
            "end_adj": right_coord_adj,
        }),
        "motif": motif,
        "ref_cn": ref_cn,
        **({"ref_start_anchor": ref_left_flank_seq[-1], "ref_seq": ref_seq} if consensus else {}),
        "reads": read_dict,
    }

    # Check now if we don't have enough reads to make a call. We can still return some read-level information!
    if n_reads_in_dict < params.min_reads:
        return {
            **call_dict_base,
            "call": None,
            "call_95_cis": None,
            "call_99_cis": None,
            "peaks": None,
            "read_peaks_called": False,
            "time": (datetime.now() - call_timer).total_seconds(),
        }

    # Now, we know we have enough reads to maybe make a call -----------------------------------------------------------

    call_data = {}
    # noinspection PyTypeChecker
    read_dict_items: tuple[tuple[str, ReadDict], ...] = tuple(read_dict.items())

    assign_method: Literal["dist", "snv", "snv+dist", "single", "hp"] = "dist"
    if n_alleles < 2:
        assign_method = "single"

    min_hp_read_coverage: int = 8  # TODO: parametrize
    min_snv_read_coverage: int = 8  # TODO: parametrize

    # Realigns are missing significant amounts of flanking information since the realignment only uses a portion of the
    # reference genome. If we have a rare realignment (e.g., a large expansion), we cannot use SNVs.
    have_rare_realigns: bool = False
    for rn, read in read_dict_items:
        read_cn = read["cn"]
        n_same_cn_no_realign = sum(1 for _, r2 in read_dict_items if not r2.get("realn") and r2["cn"] == read_cn)
        if read.get("realn") and n_same_cn_no_realign == 0:
            have_rare_realigns = True
            break

    if params.use_hp:
        top_ps = phase_sets.most_common(1)
        if (haplotagged_reads_count >= min_hp_read_coverage and len(haplotags) == n_alleles and top_ps and
                top_ps[0][1] >= min_hp_read_coverage):
            hp_sorted = sorted(haplotags)
            ps_id = top_ps[0][0]
            call_res = call_alleles_with_haplotags(
                params,
                haplotags=hp_sorted,
                ps_id=ps_id,
                read_dict_items=read_dict_items,
                rng=rng,
                logger_=logger_,
                locus_log_str=locus_log_str,
            )
            if call_res is not None:
                assign_method = "hp"
                call_data = call_res
        else:
            logger_.debug(
                f"{locus_log_str} - Not enough HP/PS tags for incorporation; one of {haplotagged_reads_count} < "
                f"{min_hp_read_coverage}, top PS {phase_sets.most_common(1)[0][1]} < {min_hp_read_coverage}, or "
                f"{len(haplotags)} != {n_alleles}")

    if should_incorporate_snvs and assign_method != "hp":
        if realign_count >= many_realigns_threshold or have_rare_realigns:
            logger_.warning(
                f"{locus_log_str} - cannot use SNVs; one of {realign_count=} >= {many_realigns_threshold} or "
                f"{have_rare_realigns=}")

        elif n_reads_in_dict >= min_snv_read_coverage and not have_rare_realigns:
            # LIMITATION: Currently can only use SNVs for haplotyping with haploid/diploid

            # Second read loop occurs in this function
            locus_snvs: set[int] = process_read_snvs_for_locus(
                contig, left_coord_adj, right_coord_adj, left_most_coord, right_most_coord, ref, read_dict_items,
                read_dict_extra, read_pairs, candidate_snvs_dict, only_known_snvs, logger_, locus_log_str)

            useful_snvs: list[tuple[int, int]] = calculate_useful_snvs(
                n_reads_in_dict, read_dict_items, read_dict_extra, read_pairs, locus_snvs, params.min_allele_reads)
            n_useful_snvs: int = len(useful_snvs)

            if not n_useful_snvs:
                logger_.debug(f"{locus_log_str} - no useful SNVs")
            else:
                am, call_res = call_alleles_with_incorporated_snvs(
                    n_alleles=n_alleles,
                    params=params,
                    read_dict=read_dict,
                    read_dict_items=read_dict_items,
                    read_dict_extra=read_dict_extra,
                    n_reads_in_dict=n_reads_in_dict,
                    useful_snvs=useful_snvs,
                    candidate_snvs_dict=candidate_snvs_dict,
                    phase_set_lock=phase_set_lock,
                    phase_set_counter=phase_set_counter,
                    snv_genotype_update_lock=snv_genotype_update_lock,
                    snv_genotype_cache=snv_genotype_cache,
                    rng=rng,
                    logger_=logger_,
                    locus_log_str=locus_log_str,
                )
                assign_method = am
                if call_res is not None:
                    call_data = call_res[0]  # Call data dictionary
                    call_dict_base["snvs"] = call_res[1]  # Called useful SNVs

        elif n_reads_in_dict < min_snv_read_coverage:
            logger_.debug(
                f"{locus_log_str} - not enough coverage for SNV incorporation "
                f"({n_reads_in_dict} < {min_snv_read_coverage})")

    single_or_dist_assign: bool = assign_method in ("single", "dist")

    if single_or_dist_assign:  # Didn't use SNVs, so call the 'old-fashioned' way - using only copy number
        # Dicts are ordered in Python; very nice :)
        rdvs = tuple(read_dict.values())
        rcns = tuple(map(cn_getter, rdvs))
        read_cns = np.fromiter(rcns, dtype=np.float_ if fractional else np.int_)
        read_weights = np.fromiter(map(weight_getter, rdvs), dtype=np.float_)
        read_weights = read_weights / np.sum(read_weights)  # Normalize to probabilities

        call_data = call_alleles(
            read_cns, (),
            read_weights, (),
            params=params,
            min_reads=params.min_reads,
            n_alleles=n_alleles,
            separate_strands=False,
            read_bias_corr_min=0,  # TODO: parametrize
            gm_filter_factor=3,  # TODO: parametrize
            seed=get_new_seed(rng),
            logger_=logger_,
            debug_str=locus_log_str,
        ) or {}  # Still false-y

    # Extract data from call_data --------------------------------------------------------------------------------------

    call = call_data.get("call")

    call_95_cis = call_data.get("call_95_cis")
    call_99_cis = call_data.get("call_99_cis")

    call_peaks = call_data.get("peaks")
    call_weights = call_data.get("peak_weights")
    call_stdevs = call_data.get("peak_stdevs")
    call_modal_n = call_data.get("modal_n_peaks")

    call_ps = call_data.get("ps")

    # Assign reads to peaks and compute peak k-mers (and optionally consensus sequences) -------------------------------

    # We cannot call read-level cluster labels with >2 peaks using distance alone;
    # don't know how re-sampling has occurred.
    call_peak_n_reads: list[int] = []
    peak_kmers: list[Counter] = [Counter() for _ in range(call_modal_n or 0)]
    call_seqs: list[str] = []
    if read_peaks_called := call_modal_n and call_modal_n <= 2:
        peaks: NDArray[np.float_] = call_peaks[:call_modal_n]
        stdevs: NDArray[np.float_] = call_stdevs[:call_modal_n]
        weights: NDArray[np.float_] = call_weights[:call_modal_n]

        allele_reads: list[list[str]] = [list() for _ in range(call_modal_n)]

        for r, rd in read_dict_items:
            # Need latter term for peaks that we overwrite if we revert to "dist" assignment:
            if not single_or_dist_assign:
                if (rp := rd.get("p")) is not None:
                    # Already has a peak from using SNV data; add it to the right allele_reads list.
                    allele_reads[rp].append(r)
                # Skip the rest; no peak assigned.
                continue

            cn = rd["cn"]

            if 0.0 in stdevs:
                # Hack: add small value to stdevs if we are exactly sure to make the below code work okay
                stdevs += 0.00001

            sd_dist = np.abs((peaks - cn) / stdevs)
            weighted_dist = np.abs(((peaks - cn) / stdevs) * weights)

            # Hack: if both peaks are 1 stdev away, pretend we aren't sure and fill in whichever allele has less
            peak: int = (
                # bool to int conversion: 1 if we add to allele_reads[1]
                int(len(allele_reads[0]) > len(allele_reads[1]))
                if call_modal_n == 2 and np.all(sd_dist < roughly_equiv_stdev_dist)
                else np.argmin(weighted_dist).item()
            )

            allele_reads[peak].append(r)
            rd["p"] = peak

            if count_kmers in ("peak", "both"):
                peak_kmers[peak] += Counter(rd["kmers"])

                # If we aren't reporting read-level k-mers, we have to delete them (space-saving!)
                if count_kmers == "peak":
                    del rd["kmers"]

        call_peak_n_reads = list(map(len, allele_reads))

        if consensus:
            call_seqs = list(
                map(lambda a: consensus_seq(map(lambda rr: read_dict_extra[rr]["_tr_seq"], a)), allele_reads)
            )

    peak_data = {
        "means": call_peaks.tolist(),  # from np.ndarray
        "weights": call_weights.tolist(),  # from np.ndarray
        "stdevs": call_stdevs.tolist(),  # from np.ndarray
        "modal_n": call_modal_n,
        "n_reads": call_peak_n_reads,
        **({"kmers": list(map(dict, peak_kmers))} if count_kmers in ("peak", "both") else {}),
        **({"seqs": call_seqs} if consensus else {}),
    } if call_data else None

    # Calculate call time ----------------------------------------------------------------------------------------------

    call_time = (datetime.now() - call_timer).total_seconds()

    if call_time > CALL_WARN_TIME:
        logger_.warning(
            f"{locus_log_str} - locus call time exceeded {CALL_WARN_TIME}s; {n_reads_in_dict} reads took {call_time}s")

    # Finally, compile the call into a dictionary with all information to return ---------------------------------------

    if fractional:
        def _ndarray_serialize(x: Iterable) -> list[Union[float, np.float_]]:
            return [round_to_base_pos(y, motif_size) for y in x]
    else:
        def _ndarray_serialize(x: Iterable) -> list[Union[int, float, np.int_, np.float_]]:
            return list(map(round, x))

    def _nested_ndarray_serialize(x: Iterable) -> list[list[Union[int, float, np.int_, np.float_]]]:
        return list(map(_ndarray_serialize, x))

    call_val = apply_or_none(_ndarray_serialize, call)
    call_95_cis_val = apply_or_none(_nested_ndarray_serialize, call_95_cis)

    logger_.debug(
        f"{locus_log_str} - got call: {call_val} (95% CIs: {call_95_cis_val}); peak assign method={assign_method}")

    return {
        **call_dict_base,
        "assign_method": assign_method,
        "call": call_val,
        "call_95_cis": call_95_cis_val,
        "call_99_cis": apply_or_none(_nested_ndarray_serialize, call_99_cis),
        "peaks": peak_data,
        "ps": call_ps,
        # make typecheck happy above by checking all of these are not None (even though if call is false-y, all of them
        # should be None and otherwise none of them should).
        "read_peaks_called": read_peaks_called,
        "time": call_time,
    }
