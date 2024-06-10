from __future__ import annotations

from .base import BaseCalculator
from .result import MIContigResult, MILocusData
from ..utils import int_tuple

__all__ = [
    "RepeatHMMCalculator",
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

                locus_start: int = int(lookup[1])
                locus_end: int = int(lookup[2])

                cr.seen_locus(contig, locus_start, locus_end)

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
                    locus_start,
                    locus_end,
                    lookup[3],

                    child_gt=int_tuple(call.split("/")),
                    mother_gt=mother_calls[lookup],
                    father_gt=father_calls[lookup],

                    logger=self._logger,
                ))

        return cr
