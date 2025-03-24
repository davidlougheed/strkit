from __future__ import annotations

import pysam

from .base import BaseCalculator
from .result import MIContigResult, MILocusData
from .vcf_utils import VCFCalculatorMixin

__all__ = ["GenericVCFLengthCalculator"]


class GenericVCFLengthCalculator(BaseCalculator, VCFCalculatorMixin):
    def _get_sample_contigs(self) -> tuple[set, set, set]:
        contigs = self.get_contigs_from_files(self._mother_call_file, self._father_call_file, self._child_call_file)
        self._logger.debug(
            "Got trio contigs - child: %d, mother: %d, father: %d",
            len(contigs[2]), len(contigs[0]), len(contigs[1]),
        )
        return contigs

    def calculate_contig(self, contig: str) -> MIContigResult:
        cr = MIContigResult(contig, includes_seq=True)

        mvf = pysam.VariantFile(str(self._mother_call_file))
        fvf = pysam.VariantFile(str(self._father_call_file))
        cvf = pysam.VariantFile(str(self._child_call_file))

        # We want all common loci, so loop through the child and then look for the loci in the parent calls

        for cv in cvf.fetch(contig):
            # child variant start/end, as determined by the reference allele sequence
            cv_start = cv.start
            cv_stop = cv.stop

            # hack for LongTR: if we override start/end in INFO, use those values as the true start/end in the context
            # of the locus boundaries
            if "START" in cv.info:
                cv_start = int(cv.info["START"]) - 1
                if "END" in cv.info:
                    cv_stop = int(cv.info["END"])

            mv = next(mvf.fetch(contig, cv_start, cv_stop), None)
            fv = next(fvf.fetch(contig, cv_start, cv_stop), None)

            # TODO: Handle sex chromosomes

            k = (contig, cv_start, cv_stop)

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

            c_seq_gt = tuple(sorted((cv.alleles[g] for g in cs["GT"]), key=len)) if None not in cs["GT"] else None
            c_gt = tuple(round(len(a) / motif_len) for a in c_seq_gt) if c_seq_gt is not None else None
            m_seq_gt = tuple(sorted((mv.alleles[g] for g in ms["GT"]), key=len)) if None not in ms["GT"] else None
            m_gt = tuple(round(len(a) / motif_len) for a in m_seq_gt) if m_seq_gt is not None else None
            f_seq_gt = tuple(sorted((fv.alleles[g] for g in fs["GT"]), key=len)) if None not in fs["GT"] else None
            f_gt = tuple(round(len(a) / motif_len) for a in f_seq_gt) if f_seq_gt is not None else None

            if c_gt is None or m_gt is None or f_gt is None:
                # None call in VCF, skip this call
                continue

            cr.append(MILocusData(
                contig=contig,
                start=cv_start,
                end=cv_stop,
                motif=motif,

                child_gt=c_gt, mother_gt=m_gt, father_gt=f_gt,

                # sequence may not line up with start/end if VCF record INFO START/END entries are used
                child_seq_gt=c_seq_gt, mother_seq_gt=m_seq_gt, father_seq_gt=f_seq_gt,
            ))

        return cr
