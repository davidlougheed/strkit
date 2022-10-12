from __future__ import annotations

import pysam

from .base import BaseCalculator
from .result import MIContigResult, MILocusData
from .vcf_utils import VCFCalculatorMixin
from ..utils import parse_cis

__all__ = ["GangSTRCalculator"]


class GangSTRCalculator(BaseCalculator, VCFCalculatorMixin):
    def _get_sample_contigs(self, include_sex_chromosomes: bool = False) -> tuple[set, set, set]:
        return self.get_contigs_from_files(self._mother_call_file, self._father_call_file, self._child_call_file)

    def calculate_contig(self, contig: str) -> MIContigResult:
        cr = MIContigResult(includes_95_ci=True)

        mvf = pysam.VariantFile(str(self._mother_call_file))
        fvf = pysam.VariantFile(str(self._father_call_file))
        cvf = pysam.VariantFile(str(self._child_call_file))

        # We want all common loci, so loop through the child and then look for the loci in the parent calls
        # TODO: What to do about filtering etc? !!!!!!!!!!!!!!!!!!!!!!!!
        #  !!!!!!!!!!!!!!!!
        #  - Q score
        #  - CIs are "proper" - not inverted or weird

        for cv in cvf.fetch(contig):
            mv = next(mvf.fetch(contig, cv.pos, cv.pos + 1), None)  # TODO: Do we need to add one
            fv = next(fvf.fetch(contig, cv.pos, cv.pos + 1), None)  # TODO: Do we need to add one

            # TODO: Handle sex chromosomes

            # Check to make sure call is present in TRF BED file, if it is specified
            k1 = (contig, str(cv.pos), str(cv.stop))
            k2 = (contig, str(cv.pos - 1), str(cv.stop))
            if self._loci_file and self._loci_dict and k1 not in self._loci_dict and k2 not in self._loci_dict:
                continue

            if self.should_exclude_locus(k1) or self.should_exclude_locus(k2):
                continue

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

            try:
                c_gt_95_ci = parse_cis(cs["REPCI"])
                m_gt_95_ci = parse_cis(ms["REPCI"])
                f_gt_95_ci = parse_cis(fs["REPCI"])
            except (ValueError, TypeError):
                # None call in VCF, skip this call
                continue

            if c_gt[0] is None or m_gt[0] is None or f_gt[0] is None:
                # None call in VCF, skip this call
                continue

            cr.append(MILocusData(
                contig=contig,
                start=cv.pos,
                end=cv.stop,
                motif=cv.info["RU"],

                child_gt=c_gt, mother_gt=m_gt, father_gt=f_gt,
                child_gt_95_ci=c_gt_95_ci, mother_gt_95_ci=m_gt_95_ci, father_gt_95_ci=f_gt_95_ci,

                reference_copies=cv.info["REF"],
            ))

        return cr
