from __future__ import annotations

import pysam

from .base import BaseCalculator
from .result import MIContigResult, MILocusData
from .vcf_utils import VCFCalculatorMixin

__all__ = ["LongTRCalculator"]


class LongTRCalculator(BaseCalculator, VCFCalculatorMixin):
    def _get_sample_contigs(self) -> tuple[set, set, set]:
        contigs = self.get_contigs_from_files(self._mother_call_file, self._father_call_file, self._child_call_file)
        self._logger.debug(
            "Got trio contigs - child: %d, mother: %d, father: %d",
            len(contigs[2]), len(contigs[0]), len(contigs[1]),
        )
        return contigs

    def calculate_contig(self, contig: str) -> MIContigResult:
        cr = MIContigResult(contig)

        mvf = pysam.VariantFile(str(self._mother_call_file))
        fvf = pysam.VariantFile(str(self._father_call_file))
        cvf = pysam.VariantFile(str(self._child_call_file))

        # We want all common loci, so loop through the child and then look for the loci in the parent calls

        for cv in cvf.fetch(contig):
            mv = next(mvf.fetch(contig, cv.start, cv.stop), None)
            fv = next(fvf.fetch(contig, cv.start, cv.stop), None)

            # TODO: Handle sex chromosomes

            k = (contig, cv.start, cv.stop)

            overlapping = self.get_loci_overlapping(k[0], k[1], k[2], True)

            if r := self.should_skip_locus(k[0], k[1], k[2], cached_overlapping=overlapping):
                self._logger.debug(f"Skipping locus {k}: {r}")
                continue

            cr.seen_locus(*k)

            if mv is None or fv is None:
                # Variant isn't found in at least one of the parents, so we can't do anything with it.
                # TODO: We need to actually check calls, and check with sample ID, not just assume
                self._logger.debug(f"Skipping locus {k}: mv or fv is None")
                continue

            # TODO: Handle missing samples gracefully
            # TODO: Handle wrong formatted VCFs gracefully

            # Need to dig up original motif from the locus file - thus, the original locus file is required.
            motif: str = overlapping[0][-1][0]
            if not motif:
                self._logger.debug(f"Skipping locus {k}: motif is false-y")
                continue

            motif_len = len(motif)

            cs = cv.samples[self._child_id or 0]
            ms = mv.samples[self._mother_id or 0]
            fs = fv.samples[self._father_id or 0]

            c_gt = (
                tuple(sorted(round(len(cv.alleles[g]) / motif_len) for g in cs["GT"])) if None not in cs["GT"] else None
            )
            m_gt = (
                tuple(sorted(round(len(mv.alleles[g]) / motif_len) for g in ms["GT"])) if None not in ms["GT"] else None
            )
            f_gt = (
                tuple(sorted(round(len(fv.alleles[g]) / motif_len) for g in fs["GT"])) if None not in fs["GT"] else None
            )

            if c_gt is None or m_gt is None or f_gt is None:
                # None call in VCF, skip this call
                continue

            cr.append(MILocusData(
                contig=contig,
                start=cv.pos,
                end=cv.stop,
                motif=motif,

                child_gt=c_gt, mother_gt=m_gt, father_gt=f_gt,
            ))

        return cr
