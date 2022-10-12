from __future__ import annotations

from .base import BaseCalculator
from .result import MIContigResult, MILocusData
from ..utils import int_tuple, parse_cis

__all__ = [
    "TandemGenotypesCalculator",
    "TandemGenotypesReCallCalculator",
]


class TandemGenotypesCalculator(BaseCalculator):
    @staticmethod
    def get_contigs_from_fh(fh) -> set[str]:
        return {ls[0] for ls in (line.split("\t") for line in fh if not line.startswith("#"))}

    @staticmethod
    def make_calls_dict(ph, contig):
        return {
            tuple(line[:4]): int_tuple(line[6:8]) if "." not in line[6:8] else (None, None)
            for line in (pv.strip().split("\t") for pv in ph if not pv.startswith("#"))
            if line[0] == contig
        }

    def _get_sample_contigs(self, include_sex_chromosomes: bool = False) -> tuple[set, set, set]:
        with open(self._mother_call_file, "r") as mvf, open(self._father_call_file, "r") as fvf, \
                open(self._child_call_file, "r") as cvf:

            mc = self.get_contigs_from_fh(mvf)
            fc = self.get_contigs_from_fh(fvf)
            cc = self.get_contigs_from_fh(cvf)

            return mc, fc, cc

    def calculate_contig(self, contig: str) -> MIContigResult:
        cr = MIContigResult()

        with open(self._mother_call_file) as mh:
            mother_calls = self.make_calls_dict(mh, contig)

        with open(self._father_call_file) as fh:
            father_calls = self.make_calls_dict(fh, contig)

        with open(self._child_call_file) as ch:
            for cv in ch:
                locus_data = cv.strip().split("\t")
                lookup = tuple(locus_data[:4])

                if locus_data[0] != contig:
                    continue

                bed_k = lookup[:3]

                # Check to make sure call is present in TRF BED file, if it is specified
                if self._loci_file and self._loci_dict and bed_k not in self._loci_dict:
                    continue

                if self.should_exclude_locus(bed_k):
                    continue

                # Check to make sure call is present in all trio individuals
                if lookup not in mother_calls or lookup not in father_calls:
                    continue

                child_calls = locus_data[6:8]

                if "." in child_calls:
                    continue

                cr.append(MILocusData(
                    contig=lookup[0],
                    start=int(lookup[1]),
                    end=int(lookup[2]),
                    motif=lookup[3],

                    child_gt=int_tuple(child_calls),
                    mother_gt=mother_calls[lookup],
                    father_gt=father_calls[lookup],
                ))

        return cr


class TandemGenotypesReCallCalculator(TandemGenotypesCalculator):
    @staticmethod
    def make_calls_dict(ph, contig):
        return {
            # Use relative indices because we may have the original calls lurking
            tuple(line[:4]): (
                int_tuple(line[-6:-4]),
                parse_cis(line[-4:-2], commas=True),
                parse_cis(line[-2:], commas=True),
            )
            for line in (pv.strip().split("\t") for pv in ph)
            if line[0] == contig and "." not in line[-6:-4]
        }

    def calculate_contig(self, contig: str) -> MIContigResult:
        cr = MIContigResult(includes_95_ci=True, includes_99_ci=True)

        with open(self._mother_call_file) as mh:
            mother_calls = self.make_calls_dict(mh, contig)

        with open(self._father_call_file) as fh:
            father_calls = self.make_calls_dict(fh, contig)

        with open(self._child_call_file) as ch:
            for cv in ch:
                locus_data = cv.strip().split("\t")
                lookup = tuple(locus_data[:4])

                if locus_data[0] != contig:
                    continue

                # Check to make sure call is present in all trio individuals
                if lookup not in mother_calls or lookup not in father_calls:
                    continue

                m_gt, m_gt_95_ci, m_gt_99_ci = mother_calls[lookup]
                f_gt, f_gt_95_ci, f_gt_99_ci = father_calls[lookup]

                calls = locus_data[-6:-4]

                if "." in calls:
                    # Failed call
                    continue

                cr.append(MILocusData(
                    contig=lookup[0],
                    start=int(lookup[1]),
                    end=int(lookup[2]),
                    motif=lookup[3],

                    child_gt=int_tuple(calls),
                    mother_gt=m_gt,
                    father_gt=f_gt,

                    child_gt_95_ci=parse_cis(locus_data[-4:-2], commas=True),
                    mother_gt_95_ci=m_gt_95_ci,
                    father_gt_95_ci=f_gt_95_ci,

                    child_gt_99_ci=parse_cis(locus_data[-2:], commas=True),
                    mother_gt_99_ci=m_gt_99_ci,
                    father_gt_99_ci=f_gt_99_ci,

                    # TODO: ref count
                ))

        return cr
