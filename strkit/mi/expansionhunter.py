from __future__ import annotations

import pysam

from .base import BaseCalculator
from .result import MIContigResult, MILocusData
from .vcf_utils import VCFCalculatorMixin
from ..utils import parse_cis

__all__ = ["ExpansionHunterCalculator"]


def _parse_allele(a: int | str | None) -> int | None:
    if isinstance(a, str):
        if a == ".":
            return None
        return int(a)
    return a


def _unzip_gt(vals) -> tuple[tuple[int | float | None, ...], tuple[int | float | None, ...]]:
    try:
        return (_parse_allele(vals[0][0]), _parse_allele(vals[1][0])), parse_cis((vals[0][1], vals[1][1]))
    except ValueError:
        return (None, None), (None, None)


class ExpansionHunterCalculator(BaseCalculator, VCFCalculatorMixin):
    def _get_sample_contigs(self) -> tuple[set, set, set]:
        return self.get_contigs_from_files(self._mother_call_file, self._father_call_file, self._child_call_file)

    def calculate_contig(self, contig: str) -> MIContigResult:
        cr = MIContigResult(contig, includes_95_ci=True)

        mvf = pysam.VariantFile(str(self._mother_call_file))
        fvf = pysam.VariantFile(str(self._father_call_file))
        cvf = pysam.VariantFile(str(self._child_call_file))

        # We want all common loci, so loop through the child and then look for the loci in the parent calls
        # TODO: What to do about filtering etc? !!!!!!!!!!!!!!!!!!!!!!!!
        #  !!!!!!!!!!!!!!!!
        #  - Q score
        #  - CIs are "proper" - not inverted or weird

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

            cs = cv.samples[self._child_id or 0]
            ms = mv.samples[self._mother_id or 0]
            fs = fv.samples[self._father_id or 0]

            cs_reps = tuple(sorted(zip(cs["REPCN"].split("/"), cs["REPCI"].split("/")), key=lambda x: x[0]))
            ms_reps = tuple(sorted(zip(ms["REPCN"].split("/"), ms["REPCI"].split("/")), key=lambda x: x[0]))
            fs_reps = tuple(sorted(zip(fs["REPCN"].split("/"), fs["REPCI"].split("/")), key=lambda x: x[0]))

            c_gt, c_gt_95_ci = _unzip_gt(cs_reps)
            m_gt, m_gt_95_ci = _unzip_gt(ms_reps)
            f_gt, f_gt_95_ci = _unzip_gt(fs_reps)

            if c_gt[0] is None or m_gt[0] is None or f_gt[0] is None:
                # None call in VCF, skip this call
                continue

            cr.append(MILocusData(
                contig=contig,
                start=cv.start,
                end=cv.stop,
                motif=cv.info["RU"],

                child_gt=c_gt, mother_gt=m_gt, father_gt=f_gt,
                child_gt_95_ci=c_gt_95_ci, mother_gt_95_ci=m_gt_95_ci, father_gt_95_ci=f_gt_95_ci,

                reference_copies=cv.info["REF"],
            ))

        return cr
