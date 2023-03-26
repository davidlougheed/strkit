from __future__ import annotations

import logging
import multiprocessing as mp
import numpy as np
import pysam
import queue

from collections import Counter
from collections.abc import Sequence
from datetime import datetime
from pysam import AlignmentFile, FastaFile
from sklearn.cluster import AgglomerativeClustering

from numpy.typing import NDArray
from typing import Iterable, Literal, Optional, Union

from strkit.call.allele import CallDict, get_n_alleles, call_alleles
from strkit.utils import apply_or_none

from .align_matrix import match_score
from .realign import realign_read
from .repeats import get_repeat_count, get_ref_repeat_count
from .snvs import get_read_snvs, call_useful_snvs
from .types import ReadDict
from .utils import normalize_contig, round_to_base_pos


__all__ = [
    "call_locus",
]


# TODO: Parameterize
CALL_WARN_TIME = 3  # seconds

SNV_OUT_OF_RANGE_CHAR = "-"


realign_timeout = 5
min_read_score = 0.9  # TODO: parametrize
base_wildcard_threshold = 3

roughly_equiv_stdev_dist = 1

force_realign = False


def calculate_seq_with_wildcards(qs: str, quals: Optional[list[int]]) -> str:
    if quals is None:
        return qs  # No quality information, so don't do anything
    return "".join(qs[i] if quals[i] > base_wildcard_threshold else "X" for i in np.arange(len(qs)))


def get_read_coords_from_matched_pairs(
    left_flank_coord: int,
    left_coord: int,
    right_coord: int,
    right_flank_coord: int,
    motif: str,
    motif_size: int,
    query_seq: str,
    matched_pairs
) -> tuple[int, int, int, int]:
    left_flank_end = -1
    right_flank_start = -1
    right_flank_end = -1

    last_idx = -1

    # Skip gaps on either side to find mapped flank indices

    # Binary search for left flank start
    # TODO: python >= 3.10: use bisect_left
    lhs = 0
    rhs = len(matched_pairs) - 1

    while rhs - lhs > 1:
        midpoint = (lhs + rhs) // 2
        rc = matched_pairs[midpoint][1]
        if rc > left_flank_coord:
            rhs = midpoint - 1
        elif rc < left_flank_coord:
            lhs = midpoint  # change from normal algorithm: inclusive, still may be closest
        else:  # equal
            lhs = midpoint
            break

    # lhs now contains the index for the closest starting coordinate to left_flank_coord
    qcc, rcc = matched_pairs[lhs]
    left_flank_start = qcc
    if rcc > left_flank_coord:
        # Error state from binary search
        return -1, -1, -1, -1

    for query_coord, ref_coord in matched_pairs[lhs+1:]:
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


def calculate_useful_snvs(
    n_reads: int,
    overlapping_segments: list[pysam.AlignedSegment],
    read_dict: dict[str, ReadDict],
    read_match_pairs: dict[str, list[tuple[int, int]]],
    locus_snvs: set[int],
) -> list[tuple[int, int]]:
    sorted_snvs: list[int] = sorted(locus_snvs)
    snv_counters: dict[int, Counter] = {sp: Counter() for sp in sorted_snvs}

    for segment in overlapping_segments:
        rn = segment.query_name
        if rn not in read_dict:
            continue
        snvs: dict[int, str] = read_dict[rn]["snv"]

        # Know this to not be None since we were passed only segments with non-None strings earlier
        qs: str = segment.query_sequence

        segment_start: int = segment.reference_start
        segment_end: int = segment.reference_end

        snv_list: list[str] = []

        for snv_pos in sorted_snvs:
            base: str = SNV_OUT_OF_RANGE_CHAR
            if snv_pos < segment_start or snv_pos > segment_end:
                # leave as gap
                pass
            else:
                if bb := snvs.get(snv_pos):
                    base = bb
                else:
                    # Binary search for pair set
                    pairs_for_read = read_match_pairs[rn]

                    def _bin_search() -> str:
                        lhs = 0
                        rhs = len(pairs_for_read) - 1

                        while lhs <= rhs:
                            pivot: int = (lhs + rhs) // 2
                            pair: tuple[int, int] = pairs_for_read[pivot]
                            if pair[1] < snv_pos:
                                lhs = pivot + 1
                            elif pair[1] > snv_pos:  # pair[1] > snv_pos
                                rhs = pivot - 1
                            else:
                                # Even if not in SNV set, it is not guaranteed to be a reference base, since
                                # it's possible it was surrounded by too much other variation during the original
                                # SNV getter algorithm.
                                return qs[pair[0]]

                        # Nothing found, so must have been a gap
                        return "_"

                    base = _bin_search()

            snv_list.append(base)
            snv_counters[snv_pos][base] += 1

        read_dict[rn]["snv_bases"] = tuple(snv_list)

    # Enough reads to try for SNV based separation
    good_snvs: list[int] = []
    for si, (snv_counted, snv_counter) in enumerate(snv_counters.items()):
        read_threshold = max(round(n_reads / 5), 2)  # TODO: parametrize
        n_alleles_meeting_threshold = 0
        for k in snv_counter:
            if k == SNV_OUT_OF_RANGE_CHAR:
                continue
            if snv_counter[k] >= read_threshold:
                n_alleles_meeting_threshold += 1
        if n_alleles_meeting_threshold >= 2:
            good_snvs.append(si)

    return [(si, sorted_snvs[si]) for si in good_snvs]  # Tuples of (index in STR list, ref position)


def calculate_read_distance(
    n_reads: int,
    read_dict_items: Sequence[tuple[str, ReadDict]],
    pure_snv_peak_assignment: bool,
    relative_cn_distance_weight_scaling: float = 0.5,  # TODO: CLI specifyable param
) -> NDArray[np.float_, np.float_]:
    """
    Calculate pairwise distance for all reads using either SNVs ONLY or a mixture of SNVs and copy number.
    :param n_reads: Number of reads.
    :param read_dict_items: Itemized read dictionary entries: (read name, read data)
    :param pure_snv_peak_assignment: Whether to use just SNVs for peak assignment
    :param relative_cn_distance_weight_scaling: How much to weight 1-difference in CN vs SNVs
            (indels are more erroneous in CCS)
    :return: The distance matrix.
    """

    # Initialize a distance matrix for all reads
    distance_matrix = np.zeros((n_reads, n_reads))

    # Loop through and compare all vs. all reads. We can skip a few indices since the distance will be symmetrical.
    for i in range(n_reads - 1):
        for j in range(i + 1, n_reads):
            r1 = read_dict_items[i][1]
            r2 = read_dict_items[j][1]
            r1_snv_u = r1["snvu"]
            r2_snv_u = r2["snvu"]

            d = 0 if pure_snv_peak_assignment else abs(r1["cn"] - r2["cn"]) * relative_cn_distance_weight_scaling
            for b1, b2 in zip(r1_snv_u, r2_snv_u):
                if b1 == SNV_OUT_OF_RANGE_CHAR or b2 == SNV_OUT_OF_RANGE_CHAR:
                    continue
                if b1 != b2:
                    d += 1

            distance_matrix[i, j] = d
            distance_matrix[j, i] = d

    return distance_matrix


def call_locus(
    t_idx: int,
    t: tuple,
    bfs: tuple[AlignmentFile, ...],
    ref: FastaFile,
    min_reads: int,
    min_allele_reads: int,
    min_avg_phred: int,
    num_bootstrap: int,
    flank_size: int,
    seed: int,
    logger_,
    sex_chroms: Optional[str] = None,
    realign: bool = False,
    hq: bool = False,
    incorporate_snvs: bool = False,
    targeted: bool = False,
    fractional: bool = False,
    respect_ref: bool = False,
    count_kmers: str = "none",  # "none" | "peak" | "read"
    log_level: int = logging.WARNING,
    read_file_has_chr: bool = True,
    ref_file_has_chr: bool = True,
) -> Optional[dict]:
    call_timer = datetime.now()

    rng = np.random.default_rng(seed=seed)

    def _get_new_seed() -> int:
        return rng.integers(0, 4096).item()

    contig: str = t[0]
    read_contig = normalize_contig(contig, read_file_has_chr)
    ref_contig = normalize_contig(contig, ref_file_has_chr)

    motif: str = t[-1]
    motif_size = len(motif)

    left_coord = int(t[1])
    right_coord = int(t[2])

    left_flank_coord = left_coord - flank_size
    right_flank_coord = right_coord + flank_size

    ref_left_flank_seq = ""
    ref_right_flank_seq = ""
    ref_seq = ""
    raised = False

    n_alleles: Optional[int] = get_n_alleles(2, sex_chroms, contig)
    if n_alleles is None:  # Sex chromosome, but we don't have a specified sex chromosome karyotype
        return None

    try:
        ref_left_flank_seq = ref.fetch(ref_contig, left_flank_coord, left_coord)
        ref_right_flank_seq = ref.fetch(ref_contig, right_coord, right_flank_coord)
        ref_seq = ref.fetch(ref_contig, left_coord, right_coord)
    except IndexError:
        logger_.warning(
            f"Coordinates out of range in provided reference FASTA for region {ref_contig} with flank size "
            f"{flank_size}: [{left_flank_coord}, {right_flank_coord}] (skipping locus {t_idx})")
        raised = True
    except ValueError:
        logger_.error(f"Invalid region '{ref_contig}' for provided reference FASTA (skipping locus {t_idx})")
        raised = True

    if len(ref_left_flank_seq) < flank_size or len(ref_right_flank_seq) < flank_size:
        if not raised:  # flank sequence too small for another reason
            logger_.warning(f"Reference flank size too small for locus {t_idx} (skipping)")
            return None

    if raised:
        return None

    # String representation of locus for logging purposes
    locus_log_str: str = f"Locus {t_idx}: {contig}:{left_coord}-{right_coord}"

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

    read_dict: dict[str, ReadDict] = {}
    chimeric_read_status: dict[str, int] = {}

    overlapping_segments: list[pysam.AlignedSegment] = []
    seen_reads: set[str] = set()
    read_lengths: list[int] = []

    segment: pysam.AlignedSegment

    for segment in (s for bf in bfs for s in bf.fetch(read_contig, left_flank_coord, right_flank_coord)):
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

        seen_reads.add(rn)
        overlapping_segments.append(segment)
        read_lengths.append(segment.query_alignment_length)

    n_overlapping_reads = len(overlapping_segments)
    sorted_read_lengths = np.sort(read_lengths)

    # Aggregations for additional read-level data
    read_kmers = Counter()
    locus_snvs: set[int] = set()

    read_pairs: dict[str, list] = {}

    for segment, read_len in zip(overlapping_segments, read_lengths):
        rn = segment.query_name
        segment_start = segment.reference_start
        segment_end = segment.reference_end

        # While .query_sequence is Optional[str], we know (because we skipped all segments with query_sequence is None
        # above) that this is guaranteed to be, in fact, not None.
        qs: str = segment.query_sequence

        c1: tuple[int, int] = segment.cigartuples[0]
        c2: tuple[int, int] = segment.cigartuples[-1]

        fqqs: Optional[list[int]] = segment.query_qualities

        realigned = False
        pairs = None

        # Soft-clipping in large insertions can result from mapping difficulties.
        # If we have a soft clip which overlaps with our TR region (+ flank), we can try to recover it
        # via realignment with parasail.
        # 4: BAM code for soft clip
        # TODO: if some of the BAM alignment is present, use it to reduce realignment overhead?
        #  - use start point + flank*3 or end point - flank*3 or something like that
        if realign and (force_realign or (
            (c1[0] == 4 and segment_start > left_flank_coord >= segment_start - c1[1]) or
            (c2[0] == 4 and segment_end < right_flank_coord <= segment_end + c2[1])
        )):
            # Run the realignment in a separate process, to give us a timeout mechanism.
            # This means we're spawning a second process for this job, just momentarily, beyond the pool size.

            q = mp.Queue()
            proc = mp.Process(target=realign_read, kwargs=dict(
                # fetch an extra base for the right flank coordinate check later (needs to be >= the exclusive coord)
                ref_seq=ref_left_flank_seq + ref_seq + ref.fetch(ref_contig, right_coord, right_flank_coord + 1),
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
            finally:
                proc.close()

            if pairs_new is not None:
                pairs = pairs_new
                realigned = True

        if pairs is None:
            pairs = segment.get_aligned_pairs(matches_only=True)

        # Cache aligned pairs, since it takes a lot of time to extract, and we use it for calculate_useful_snvs
        read_pairs[rn] = pairs

        left_flank_start, left_flank_end, right_flank_start, right_flank_end = get_read_coords_from_matched_pairs(
            left_flank_coord,
            left_coord_adj,
            right_coord_adj,
            right_flank_coord,
            motif,
            motif_size,
            query_seq=qs,
            matched_pairs=pairs
        )

        if any(v == -1 for v in (left_flank_start, left_flank_end, right_flank_start, right_flank_end)):
            logger_.debug(
                f"{locus_log_str} - skipping read {rn}: could not get sufficient flanking sequence"
                f"{' (post-realignment)' if realigned else ''}")
            continue

        qqs = np.array(fqqs[left_flank_end:right_flank_start])
        if qqs.shape[0] and (m_qqs := np.mean(qqs)) < min_avg_phred:  # TODO: check flank?
            logger_.debug(
                f"{locus_log_str} - skipping read {rn} due to low average base quality ({m_qqs:.2f} < {min_avg_phred})")
            continue

        # -----

        tr_read_seq = qs[left_flank_end:right_flank_start]

        # Truncate to flank_size (plus some leeway for small indels in flanking region) to stop any expansion sequences
        # from accidentally being included in the flanking region; e.g. if the insert gets mapped onto bases outside
        # the definition coordinates.
        # The +10 here won't include any real TR region if the mapping is solid, since the flank coordinates will
        # contain a correctly-sized sequence.

        # TODO: wildcards in flanking region too?
        flank_left_seq = qs[left_flank_start:left_flank_end][:flank_size+10]
        flank_right_seq = qs[right_flank_start:right_flank_end][-(flank_size+10):]

        tr_len = len(tr_read_seq)
        flank_len = len(flank_left_seq) + len(flank_right_seq)
        tr_len_w_flank = tr_len + flank_len

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
            logger_.debug(
                f"{locus_log_str} - skipping read {segment.query_name} (scored {read_adj_score} < {min_read_score})")
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

        mean_containing_size = read_len if targeted else np.mean(sorted_read_lengths[partition_idx:]).item()
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

        if incorporate_snvs:
            snvs = get_read_snvs(qs, pairs, contig, ref, left_coord_adj, right_coord_adj)
            locus_snvs |= set(snvs.keys())
            read_dict[rn]["snv"] = snvs

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
        "reads": read_dict,
    }

    # Check now if we definitely don't have enough reads to make a call
    # We also check again later when we calculate all the flanking stuff
    if len(overlapping_segments) < min_reads:
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
    n_reads_in_dict: int = len(read_dict)
    # noinspection PyTypeChecker
    read_dict_items: tuple[tuple[str, ReadDict], ...] = tuple(read_dict.items())

    assign_method: Literal["dist", "snv", "snv+dist", "single"] = "dist"
    if n_alleles < 2:
        assign_method = "single"
    min_snv_read_coverage: int = 10  # TODO: parametrize

    # LIMITATION: Currently can only use SNVs for haplotyping with haploid/diploid
    if n_alleles == 2 and incorporate_snvs:
        useful_snvs: list[tuple[int, int]] = (
            calculate_useful_snvs(n_reads_in_dict, overlapping_segments, read_dict, read_pairs, locus_snvs)
            if n_reads_in_dict >= min_snv_read_coverage else []
        )
        n_useful_snvs: int = len(useful_snvs)

        if not n_useful_snvs:
            logger_.debug(f"{locus_log_str} - no useful SNVs")
        else:
            # TODO: parametrize min 'enough to do pure SNV haplotyping' thresholds

            read_dict_items_with_many_snvs: list[tuple[str, ReadDict]] = []
            read_dict_items_with_some_snvs: list[tuple[str, ReadDict]] = []
            read_dict_items_with_few_or_no_snvs: list[tuple[str, ReadDict]] = []

            print_snvs = False

            for rn, read in read_dict_items:
                read_useful_snv_bases = tuple(read["snv_bases"][bi] for bi, _pos in useful_snvs)
                non_blank_read_useful_snv_bases = [bb for bb in read_useful_snv_bases if bb != SNV_OUT_OF_RANGE_CHAR]

                if (nbr := len(non_blank_read_useful_snv_bases)) >= 2:  # TODO: parametrize
                    read_dict_items_with_some_snvs.append((rn, read))
                    if nbr >= 3:  # TODO: parametrize
                        read_dict_items_with_many_snvs.append((rn, read))
                else:
                    read_dict_items_with_few_or_no_snvs.append((rn, read))

                read["snvu"] = read_useful_snv_bases  # Store read-level 'useful' SNVs

                if print_snvs:
                    print(rn, f"\t{read['cn']:.0f}", "\t", "".join(read_useful_snv_bases),
                          len(non_blank_read_useful_snv_bases))

            n_reads_with_many_snvs: int = len(read_dict_items_with_many_snvs)
            pure_snv_peak_assignment: bool = n_reads_with_many_snvs == n_reads_in_dict

            # TODO: parametrize: how many reads with SNV information
            min_snv_incorporation_read_portion = 0.5
            min_snv_incorporation_read_abs = 16
            if pure_snv_peak_assignment or (
                    len(read_dict_items_with_some_snvs) >=
                    min(n_reads_in_dict * min_snv_incorporation_read_portion, min_snv_incorporation_read_abs)):
                if pure_snv_peak_assignment:
                    # We have enough SNVs in ALL reads, so we can phase purely based on SNVs
                    logger_.debug(f"{locus_log_str} - haplotyping purely using SNVs")
                    assign_method = "snv"
                else:
                    # We have enough SNVs in lots of reads, so we can phase using a combined metric
                    logger_.debug(f"{locus_log_str} - haplotyping using combined STR-SNV metric")
                    # TODO: Handle reads we didn't have SNVs for by retroactively assigning to groups
                    assign_method = "snv+dist"

                # Calculate pairwise distance for all reads using either SNVs ONLY or
                # a mixture of SNVs and copy number:
                dm = calculate_read_distance(n_reads_in_dict, read_dict_items, pure_snv_peak_assignment)

                # Cluster reads together using the distance matrix, which incorporates
                # SNV and possibly copy number information.
                c = AgglomerativeClustering(n_clusters=n_alleles, metric="precomputed", linkage="average").fit(dm)

                # noinspection PyUnresolvedReferences
                cluster_labels = c.labels_
                cluster_indices = tuple(range(n_alleles))

                cns: Union[list[list[int]], list[list[float]]] = []
                c_ws: list[Union[NDArray[np.int_], NDArray[np.float_]]] = []

                for ci in cluster_indices:
                    crs = []
                    # Find reads and assign peaks
                    for i, (_, r) in enumerate(read_dict_items):
                        if cluster_labels[i] == ci:
                            crs.append(r)
                            r["p"] = ci  # Mutation: assign peak index to read data dictionary

                    # Calculate copy number set
                    cns.append([r["cn"] for r in crs])

                    # Calculate weights array
                    ws = np.array([r["w"] for r in crs])
                    ws = ws / np.sum(ws)
                    c_ws.append(ws)

                # noinspection PyUnresolvedReferences

                cdd: list[Optional[CallDict]] = [
                    call_alleles(
                        cns[ci], (),  # Don't bother separating by strand for now...
                        c_ws[ci], (),
                        bootstrap_iterations=num_bootstrap // n_alleles,  # Apportion the bootstrap iters across alleles
                        min_reads=min_allele_reads,  # Calling alleles separately, so set min_reads=min_allele_reads
                        min_allele_reads=min_allele_reads,
                        n_alleles=1,  # Calling alleles separately: they were pre-separated by agglom. clustering
                        separate_strands=False,
                        read_bias_corr_min=0,  # separate_strands is false, so this is ignored
                        gm_filter_factor=1,  # n_alleles=1, so this is ignored
                        hq=hq,
                        force_int=not fractional,
                        seed=_get_new_seed(),
                        logger_=logger_,
                        debug_str=f"{contig}:{left_coord}-{right_coord} a0"
                    )
                    for ci in cluster_indices
                ]

                if not any((not cd) for cd in cdd):  # cdd is confirmed list[CallDict], no Nones
                    # TODO: Multi-allele phasing across STRs

                    # We called these as single-allele (1 peak) loci as a sort of hack, so the return "call" key
                    # is an array of length 1.

                    cdd_calls = np.array([x["call"][0] for x in cdd])
                    peak_order = np.argsort(cdd_calls)  # To reorder call arrays in least-to-greatest by copy number

                    # All call_datas are truth-y; all arrays should be ordered by peak_order
                    call_data = {
                        "call": cdd_calls[peak_order],
                        "call_95_cis": np.concatenate(tuple(cc["call_95_cis"] for cc in cdd), axis=0)[peak_order],
                        "call_99_cis": np.concatenate(tuple(cc["call_99_cis"] for cc in cdd), axis=0)[peak_order],
                        "peaks": np.concatenate(tuple(cc["peaks"] for cc in cdd), axis=None)[peak_order],

                        # TODO: Readjust peak weights when combining or don't include
                        "peak_weights": np.concatenate(tuple(cc["peak_weights"] for cc in cdd), axis=0)[peak_order],

                        "peak_stdevs": np.concatenate(tuple(cc["peak_stdevs"] for cc in cdd), axis=0)[peak_order],
                        "modal_n_peaks": n_alleles,  # # alleles = # peaks always -- if we phased using SNVs
                    }

                    # Add SNV data to final return dictionary
                    call_dict_base["snvs"] = call_useful_snvs(
                        n_alleles, read_dict, useful_snvs, peak_order, locus_log_str, logger_)
                else:
                    # One of the calls could not be made... what to do?
                    # TODO: !!!!
                    #  For now, revert to dist
                    assign_method = "dist"
            else:
                # TODO: How to use partial data?
                pass

        # TODO:
        #  - use reads with at least 1 SNPs called to separate 'training reads'
        #  - eliminate 1-SNV reads with a unique haplotype vs any non-1-SNP reads... or something like that.
        #    (see m64012_190920_173625/31917565/ccs)
        #  - then, assign all reads - either based on what group they 'trained' or (for no SNP call ones)
        #    based on the GMMs (maybe some kind of certainty of assignment???)

        # TODO:
        #  - maybe 3 approaches:
        #    - if not enough SNV info / almost no reads have it, just do old method.
        #    - if we have some SNV info for all reads AND it's not a single read count (perfectly homozygous)
        #      (how to quantify?), do distance-based with read copy numbers AND SNV data - e.g. if we have just 1 SNP
        #    - if we have LOTS of SNV data for the majority of reads, do assignment just using SNV data and assign
        #      reads to peaks after (random assignment if equal chance of both groups)
        #    - if we have COMPLETE SNV data (4+ for every read) we can do phasing + peak assignment just with SNVs
        #      and just call the Gaussians from the separated reads (even just calculate stdev + mean for each bootstrap
        #      iteration, which might save us some time...).
        #    - keep track of which option as 'peak_calling_method' (pcm) or something

    single_or_dist_assign: bool = assign_method in ("single", "dist")

    if single_or_dist_assign:  # Didn't use SNVs, so call the 'old-fashioned' way - using only copy number
        # Dicts are ordered in Python; very nice :)
        rdvs = tuple(read_dict.values())
        rcns = tuple(r["cn"] for r in rdvs)
        read_cns = np.fromiter(rcns, dtype=np.float_ if fractional else np.int_)
        read_weights = np.fromiter((r["w"] for r in rdvs), dtype=np.float_)
        read_weights = read_weights / np.sum(read_weights)  # Normalize to probabilities

        call_data = call_alleles(
            read_cns, (),
            read_weights, (),
            bootstrap_iterations=num_bootstrap,
            min_reads=min_reads,
            min_allele_reads=min_allele_reads,
            n_alleles=n_alleles,
            separate_strands=False,
            read_bias_corr_min=0,  # TODO: parametrize
            gm_filter_factor=3,  # TODO: parametrize
            hq=hq,
            force_int=not fractional,
            seed=_get_new_seed(),
            logger_=logger_,
            debug_str=f"{contig}:{left_coord}-{right_coord}",
        ) or {}  # Still false-y

    # Clean up read dict  TODO: nicer way to do this -------------------------------------------------------------------

    for rn, rd in read_dict.items():
        if "snv" in rd:
            del rd["snv"]
        if "snv_bases" in rd:
            del rd["snv_bases"]

    # Extract data from call_data --------------------------------------------------------------------------------------

    call = call_data.get("call")

    call_95_cis = call_data.get("call_95_cis")
    call_99_cis = call_data.get("call_99_cis")

    call_peaks = call_data.get("peaks")
    call_weights = call_data.get("peak_weights")
    call_stdevs = call_data.get("peak_stdevs")
    call_modal_n = call_data.get("modal_n_peaks")

    # Assign reads to peaks and compute peak k-mers --------------------------------------------------------------------

    # We cannot call read-level cluster labels with >2 peaks using distance alone;
    # don't know how re-sampling has occurred.
    call_peak_n_reads = []
    peak_kmers: list[Counter] = [Counter() for _ in range(call_modal_n or 0)]
    if read_peaks_called := call_modal_n and call_modal_n <= 2:
        peaks: NDArray = call_peaks[:call_modal_n]
        stdevs: NDArray[np.float_] = call_stdevs[:call_modal_n]
        weights: NDArray[np.float_] = call_weights[:call_modal_n]

        allele_reads: list[list[str]] = []
        for _ in range(call_modal_n):
            allele_reads.append([])

        for r, rd in read_dict.items():
            # Need latter term for peaks that we overwrite if we revert to "dist" assignment:
            if (rp := rd.get("p")) is not None and not single_or_dist_assign:
                # Already has a peak from using SNV data; add it to the right allele_reads list and skip the rest.
                allele_reads[rp].append(r)
                continue

            cn = rd["cn"]

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
                peak_kmers[peak] += rd["kmers"]

                # If we aren't reporting read-level k-mers, we have to delete them (space-saving!)
                if count_kmers == "peak":
                    del rd["kmers"]

        call_peak_n_reads = list(map(len, allele_reads))

    peak_data = {
        "means": call_peaks.tolist(),  # from np.ndarray
        "weights": call_weights.tolist(),  # from np.ndarray
        "stdevs": call_stdevs.tolist(),  # from np.ndarray
        "modal_n": call_modal_n,
        "n_reads": call_peak_n_reads,
        **({"kmers": [dict(c) for c in peak_kmers]} if count_kmers in ("peak", "both") else {}),
    } if call_data else None

    # Calculate call time ----------------------------------------------------------------------------------------------

    call_time = (datetime.now() - call_timer).total_seconds()

    if call_time > CALL_WARN_TIME:
        logger_.warning(
            f"{locus_log_str} - locus call time exceeded {CALL_WARN_TIME}s; {n_reads_in_dict} reads took {call_time}s")

    # Finally, compile the call into a dictionary with all information to return ---------------------------------------

    def _ndarray_serialize(x: Iterable) -> list[Union[int, float, np.int_, np.float_]]:
        return [(round(y) if not fractional else round_to_base_pos(y, motif_size)) for y in x]

    def _nested_ndarray_serialize(x: Iterable) -> list[list[Union[int, float, np.int_, np.float_]]]:
        return [_ndarray_serialize(y) for y in x]

    return {
        **call_dict_base,
        "assign_method": assign_method,
        "call": apply_or_none(_ndarray_serialize, call),
        "call_95_cis": apply_or_none(_nested_ndarray_serialize, call_95_cis),
        "call_99_cis": apply_or_none(_nested_ndarray_serialize, call_99_cis),
        "peaks": peak_data,
        # make typecheck happy above by checking all of these are not None (even though if call is false-y, all of them
        # should be None and otherwise none of them should).
        "read_peaks_called": read_peaks_called,
        "time": call_time,
    }
