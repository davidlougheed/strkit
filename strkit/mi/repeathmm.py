from __future__ import annotations

from .base import BaseCalculator
from .result import MIContigResult, MILocusData
from ..utils import int_tuple, parse_cis

__all__ = [
    "RepeatHMMCalculator",
    "RepeatHMMReCallCalculator",
]


class RepeatHMMCalculator(BaseCalculator):
    @staticmethod
    def get_contigs_from_fh(fh) -> set:
        return {ls[0] for ls in (line.split(":") for line in fh)}

    @staticmethod
    def make_calls_dict(ph, contig):
        return {
            tuple(k.split(":")): int_tuple(v.split("/"))
            for k, v in (pv.split() for pv in ph)
            if k.split(":")[0] == contig
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
                locus_data, call = cv.strip().split(" ")
                lookup = tuple(locus_data.split(":"))

                if lookup[0] != contig:
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

                c_gt = int_tuple(call.split("/"))
                m_gt = mother_calls[lookup]
                f_gt = father_calls[lookup]

                # Failed calls from RepeatHMM seem to be represented as 0/0, so skip this
                # TODO… Need to decide if we actually want to include these?
                #  or at least somehow record them
                if (0, 0) in (c_gt, m_gt, f_gt):
                    continue

                # TODO: Include ref copies... should be in file somewhere?
                cr.append(MILocusData(
                    lookup[0],
                    int(lookup[1]),
                    int(lookup[2]),
                    lookup[3],

                    child_gt=int_tuple(call.split("/")),
                    mother_gt=mother_calls[lookup],
                    father_gt=father_calls[lookup],
                ))

        return cr


class RepeatHMMReCallCalculator(RepeatHMMCalculator):
    @staticmethod
    def make_calls_dict(ph, contig):
        return {
            tuple(v[0].split(":")): (
                int_tuple(v[1:3]),
                parse_cis(v[3:5], commas=True),
                parse_cis(v[5:7], commas=True),
            )
            for v in (pv.split("\t") for pv in ph)
            if v[0].split(":")[0] == contig and "." not in v[1:3]
        }

    # TODO: Deduplicate with above
    def calculate_contig(self, contig: str) -> MIContigResult:
        cr = MIContigResult(includes_95_ci=True, includes_99_ci=True)

        with open(self._mother_call_file) as mh:
            mother_calls = self.make_calls_dict(mh, contig)

        with open(self._father_call_file) as fh:
            father_calls = self.make_calls_dict(fh, contig)

        with open(self._child_call_file) as ch:
            for cv in ch:
                locus_data = cv.strip().split("\t")
                lookup = tuple(locus_data[0].split(":"))

                if lookup[0] != contig:
                    continue

                # Check to make sure call is present in all trio individuals
                if lookup not in mother_calls or lookup not in father_calls:
                    continue

                # TODO: What will failed calls look like here? Do we also check for 0?

                m_gt, m_gt_95_ci, m_gt_99_ci = mother_calls[lookup]
                f_gt, f_gt_95_ci, f_gt_99_ci = father_calls[lookup]

                # print(m_gt, m_gt_95_ci, m_gt_99_ci)

                calls = locus_data[1:3]

                if "." in calls:
                    # Failed call
                    continue

                c_gt = int_tuple(calls)
                c_gt_95_ci = parse_cis(locus_data[3:5], commas=True)
                c_gt_99_ci = parse_cis(locus_data[5:7], commas=True)

                if (0, 0) in (c_gt, m_gt, f_gt):  # TODO
                    # Failed call
                    continue

                if self._debug:  # TODO: Real logging
                    print(f"c_gt={c_gt} c_gt_95_ci={c_gt_95_ci}")
                    print(f"m_gt={m_gt} m_gt_95_ci={m_gt_95_ci}")
                    print(f"f_gt={f_gt} f_gt_95_ci={f_gt_95_ci}")

                # TODO: Put ref # here, since we have it with the detail thing
                cr.append(MILocusData(
                    lookup[0],
                    int(lookup[1]),
                    int(lookup[2]),
                    lookup[3],

                    child_gt=int_tuple(calls),
                    mother_gt=m_gt,
                    father_gt=f_gt,

                    child_gt_95_ci=c_gt_95_ci,
                    mother_gt_95_ci=m_gt_95_ci,
                    father_gt_95_ci=f_gt_95_ci,

                    child_gt_99_ci=c_gt_99_ci,
                    mother_gt_99_ci=m_gt_99_ci,
                    father_gt_99_ci=f_gt_99_ci,
                ))

        return cr
