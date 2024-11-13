from __future__ import annotations

import pysam

from typing import Optional, Union

from .base import BaseCalculator
from .result import MIContigResult, MILocusData
from .vcf_utils import VCFCalculatorMixin
from ..utils import parse_ci

__all__ = ["TRGTCalculator"]


def _parse_allele(a: Union[int, str, None]) -> Optional[int]:
    if isinstance(a, str):
        if a == ".":
            return None
        return int(a)
    return a


def _unzip_gt(
    vals, motif_len: int
) -> Union[
    tuple[tuple[int, ...], tuple[tuple[int, ...], tuple[int, ...]]],
    tuple[tuple[None, None], tuple[None, None]],
]:
    try:
        return (
            (
                round(_parse_allele(vals[0][0]) / motif_len),
                round(_parse_allele(vals[1][0]) / motif_len),
            ),
            (
                tuple(map(lambda x: round(x / motif_len), parse_ci(vals[0][1]))),
                tuple(map(lambda x: round(x / motif_len), parse_ci(vals[1][1]))),
            ),
        )
    except (ValueError, TypeError):
        return (None, None), (None, None)


class TRGTCalculator(BaseCalculator, VCFCalculatorMixin):
    def _get_sample_contigs(self) -> tuple[set, set, set]:
        return self.get_contigs_from_files(self._mother_call_file, self._father_call_file, self._child_call_file)

    def calculate_contig(self, contig: str) -> MIContigResult:
        cr = MIContigResult(contig, includes_95_ci=True)

        mvf = pysam.VariantFile(str(self._mother_call_file))
        fvf = pysam.VariantFile(str(self._father_call_file))
        cvf = pysam.VariantFile(str(self._child_call_file))

        # We want all common loci, so loop through the child and then look for the loci in the parent calls

        for cv in cvf.fetch(contig):
            mv = next(mvf.fetch(contig, cv.start, cv.stop), None)
            fv = next(fvf.fetch(contig, cv.start, cv.stop), None)

            # TODO: Handle sex chromosomes

            k = (contig, cv.start, cv.stop)

            if self.should_skip_locus(*k):
                continue

            cr.seen_locus(*k)

            if mv is None or fv is None:
                # Variant isn't found in at least one of the parents, so we can't do anything with it.
                # TODO: We need to actually check calls, and check with sample ID, not just assume
                continue

            # TODO: Handle missing samples gracefully
            # TODO: Handle wrong formatted VCFs gracefully

            motif = cv.info["MOTIFS"][0]

            cs = cv.samples[self._child_id or 0]
            ms = mv.samples[self._mother_id or 0]
            fs = fv.samples[self._father_id or 0]

            cs_reps = tuple(sorted(zip(cs["AL"], cs["ALLR"]), key=lambda x: x[0]))
            ms_reps = tuple(sorted(zip(ms["AL"], ms["ALLR"]), key=lambda x: x[0]))
            fs_reps = tuple(sorted(zip(fs["AL"], fs["ALLR"]), key=lambda x: x[0]))

            c_gt, c_gt_95_ci = _unzip_gt(cs_reps, len(motif))
            m_gt, m_gt_95_ci = _unzip_gt(ms_reps, len(motif))
            f_gt, f_gt_95_ci = _unzip_gt(fs_reps, len(motif))

            if c_gt[0] is None or m_gt[0] is None or f_gt[0] is None:
                # None call in VCF, skip this call
                continue

            cr.append(MILocusData(
                contig=contig,
                start=cv.pos,
                end=cv.stop,
                motif=motif,

                child_gt=c_gt, mother_gt=m_gt, father_gt=f_gt,
                child_gt_95_ci=c_gt_95_ci, mother_gt_95_ci=m_gt_95_ci, father_gt_95_ci=f_gt_95_ci,
            ))

        return cr
