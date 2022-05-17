from typing import Tuple

from .base import BaseCalculator
from .result import MIContigResult, MILocusData
from ..utils import int_tuple, parse_cis

__all__ = [
    "StrKitCalculator",
]


class StrKitCalculator(BaseCalculator):
    @staticmethod
    def get_contigs_from_fh(fh) -> set:
        return {ls[0] for ls in (line.split("\t") for line in fh if not line.startswith("#"))}

    def _get_sample_contigs(self, include_sex_chromosomes: bool = False) -> Tuple[set, set, set]:
        with open(self._mother_call_file, "r") as mvf, open(self._father_call_file, "r") as fvf, \
                open(self._child_call_file, "r") as cvf:
            mc = self.get_contigs_from_fh(mvf)
            fc = self.get_contigs_from_fh(fvf)
            cc = self.get_contigs_from_fh(cvf)

            return mc, fc, cc

    @staticmethod
    def make_calls_dict(ph, contig):
        return {
            tuple(line[:4]): (
                int_tuple(line[-2].split("|")),
                parse_cis(line[-1].split("|")),
                None  # parse_cis(line[-1:].split("|")),
            )
            for line in (pv.strip().split("\t") for pv in ph)
            if line[0] == contig and "." not in line[-2].split("|")
        }

    def calculate_contig(self, contig: str) -> MIContigResult:
        cr = MIContigResult(includes_95_ci=True)

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

                m_gt, m_gt_95_ci, _ = mother_calls[lookup]
                f_gt, f_gt_95_ci, _ = father_calls[lookup]

                calls = locus_data[-2].split("|")

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

                    child_gt_95_ci=parse_cis(locus_data[-1].split("|")),
                    mother_gt_95_ci=m_gt_95_ci,
                    father_gt_95_ci=f_gt_95_ci,

                    # child_gt_99_ci=parse_cis(locus_data[-1:].split("|")),
                    # mother_gt_99_ci=m_gt_99_ci,
                    # father_gt_99_ci=f_gt_99_ci,

                    reference_copies=int(locus_data[4]),
                ))

        return cr
