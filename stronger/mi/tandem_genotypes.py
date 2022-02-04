from typing import Tuple

from .base import BaseCalculator
from ..utils import int_tuple, parse_cis

__all__ = [
    "TandemGenotypesCalculator",
    "TandemGenotypesReCallCalculator",
]


class TandemGenotypesCalculator(BaseCalculator):
    @staticmethod
    def get_contigs_from_fh(fh) -> set:
        return {ls[0] for ls in (line.split("\t") for line in fh if not line.startswith("#"))}

    @staticmethod
    def make_calls_dict(ph, contig):
        return {
            tuple(line[:4]): int_tuple(line[8:10])
            for line in (pv.strip().split("\t") for pv in ph if not pv.startswith("#"))
            if line[0] == contig
        }

    def _get_sample_contigs(self, include_sex_chromosomes: bool = False) -> Tuple[set, set, set]:
        with open(self._mother_call_file, "r") as mvf, open(self._father_call_file, "r") as fvf, \
                open(self._child_call_file, "r") as cvf:

            mc = self.get_contigs_from_fh(mvf)
            fc = self.get_contigs_from_fh(fvf)
            cc = self.get_contigs_from_fh(cvf)

            return mc, fc, cc

    def calculate_contig(self, contig: str):
        value = 0  # Sum of 1s for the eventual MI % calculation
        n_loci = 0

        non_matching = []

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

                # TODO: What do failed calls look like here?

                n_loci += 1

                c_gt = int_tuple(locus_data[8:10])
                m_gt = mother_calls[lookup]
                f_gt = father_calls[lookup]

                respects_mi_strict, _ = self.gts_respect_mi(c_gt, m_gt, f_gt)
                if respects_mi_strict:
                    # Mendelian inheritance upheld for this locus - strict
                    value += 1
                else:
                    non_matching.append((
                        *lookup,

                        c_gt, "",
                        m_gt, "",
                        f_gt, "",

                        "",
                    ))

        return value, None, n_loci, non_matching


class TandemGenotypesReCallCalculator(TandemGenotypesCalculator):
    @staticmethod
    def make_calls_dict(ph, contig):
        return {
            # Use relative indices because we may have the original calls lurking
            tuple(line[:4]): (
                tuple(map(int, line[-6:-4])),
                parse_cis(line[-4:-2], commas=True),
                parse_cis(line[-2:], commas=True),
            )
            for line in (pv.strip().split("\t") for pv in ph)
            if line[0] == contig and "." not in line[-6:-4]
        }

    def calculate_contig(self, contig: str):
        value = 0  # Sum of 1s for the eventual MI % calculation
        value_95_ci = 0
        n_loci = 0

        non_matching = []

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

                # TODO: What do failed calls look like here?

                m_gt, m_gt_95_ci, _ = mother_calls[lookup]
                f_gt, f_gt_95_ci, _ = father_calls[lookup]

                calls = locus_data[-6:-4]

                if "." in calls:
                    # Failed call
                    continue

                c_gt = int_tuple(calls)
                c_gt_95_ci = parse_cis(locus_data[-4:-2], commas=True)
                # c_gt_99_ci = parse_cis(locus_data[-2:], commas=True)

                n_loci += 1

                respects_mi_strict, respects_mi_95_ci = self.gts_respect_mi(
                    c_gt=c_gt, m_gt=m_gt, f_gt=f_gt,
                    c_gt_ci=c_gt_95_ci, m_gt_ci=m_gt_95_ci, f_gt_ci=f_gt_95_ci
                )

                if respects_mi_strict:
                    # Mendelian inheritance upheld for this locus - strict
                    value += 1

                if respects_mi_95_ci:
                    # Mendelian inheritance upheld for this locus - within 95% CI from TG2MM
                    value_95_ci += 1
                else:
                    non_matching.append((
                        *lookup,

                        c_gt, c_gt_95_ci,
                        m_gt, m_gt_95_ci,
                        f_gt, f_gt_95_ci,

                        "",
                    ))

        return value, value_95_ci, n_loci, non_matching
