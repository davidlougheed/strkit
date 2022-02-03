import pysam

from typing import List, Tuple

from .base import BaseCalculator
from .vcf_utils import VCFCalculatorMixin
from ..utils import parse_cis

__all__ = ["GangSTRCalculator"]


class GangSTRCalculator(BaseCalculator, VCFCalculatorMixin):
    def _get_sample_contigs(self, include_sex_chromosomes: bool = False) -> Tuple[set, set, set]:
        return self.get_contigs_from_files(self._mother_call_file, self._father_call_file, self._child_call_file)

    def calculate_contig(self, contig: str) -> Tuple[int, int, int, List[Tuple]]:
        mvf = pysam.VariantFile(str(self._mother_call_file))
        fvf = pysam.VariantFile(str(self._father_call_file))
        cvf = pysam.VariantFile(str(self._child_call_file))

        value = 0  # Sum of 1s for the eventual MI % calculation
        value_ci = 0  # Sum of 1s for the eventual MI % calculation (including GT CI)
        n_loci = 0

        non_matching = []

        # We want all common loci, so loop through the child and then look for the loci in the parent calls
        # TODO: What to do about filtering etc? !!!!!!!!!!!!!!!!!!!!!!!!
        #  !!!!!!!!!!!!!!!!
        #  - Q score
        #  - CIs are "proper" - not inverted or weird

        for cv in cvf.fetch(contig):
            mv = next(mvf.fetch(contig, cv.pos, cv.pos + 1), None)  # TODO: Do we need to add one
            fv = next(fvf.fetch(contig, cv.pos, cv.pos + 1), None)  # TODO: Do we need to add one

            # TODO: Handle sex chromosomes

            if mv is None or fv is None:
                # Variant isn't found in at least one of the parents, so we can't do anything with it.
                # TODO: We need to actually check calls, and check with sample ID, not just assume
                continue

            # TODO: Handle missing samples gracefully
            # TODO: Handle wrong formatted VCFs gracefully

            cs = cv.samples[self._child_id or 0]
            ms = mv.samples[self._mother_id or 0]
            fs = fv.samples[self._father_id or 0]

            c_gt = cs["REPCN"]
            m_gt = ms["REPCN"]
            f_gt = fs["REPCN"]

            c_gt_ci = parse_cis(cs["REPCI"])
            m_gt_ci = parse_cis(ms["REPCI"])
            f_gt_ci = parse_cis(fs["REPCI"])

            if c_gt[0] is None or m_gt[0] is None or f_gt[0] is None:
                # None call in VCF, skip this call
                continue

            n_loci += 1

            respects_mi_strict, respects_mi_ci = self.gts_respect_mi(
                c_gt=c_gt, m_gt=m_gt, f_gt=f_gt,
                c_gt_ci=c_gt_ci, m_gt_ci=m_gt_ci, f_gt_ci=f_gt_ci
            )

            if respects_mi_strict:
                # Mendelian inheritance upheld for this locus - strict (MLE of GT)
                value += 1

            if respects_mi_ci:
                # Mendelian inheritance upheld for this locus - within 95% CI from GangSTR
                value_ci += 1
            else:
                non_matching.append((
                    contig,
                    cv.pos,
                    cv.info["END"],
                    cv.info["RU"],

                    c_gt,
                    c_gt_ci,

                    m_gt,
                    m_gt_ci,

                    f_gt,
                    f_gt_ci,

                    cv.info["REF"],
                ))

        return value, value_ci, n_loci, non_matching
