from __future__ import annotations

from .base import BaseCalculator
from .result import MIContigResult, MILocusData
from ..utils import int_tuple

__all__ = [
    "TandemGenotypesCalculator",
]


class TandemGenotypesCalculator(BaseCalculator):
    @staticmethod
    def get_contigs_from_fh(fh) -> set[str]:
        return {ls[0] for ls in (line.split("\t") for line in fh if not line.startswith("#"))}

    @staticmethod
    def make_calls_dict(ph, contig):
        return {
            tuple(line[:4]): int_tuple(line[6:8])
            for line in (pv.strip().split("\t") for pv in ph if not pv.startswith("#"))
            if line[0] == contig and "." not in line[6:8]
        }

    def _get_sample_contigs(self, include_sex_chromosomes: bool = False) -> tuple[set, set, set]:
        with open(self._mother_call_file, "r") as mvf, open(self._father_call_file, "r") as fvf, \
                open(self._child_call_file, "r") as cvf:

            mc = self.get_contigs_from_fh(mvf)
            fc = self.get_contigs_from_fh(fvf)
            cc = self.get_contigs_from_fh(cvf)

            return mc, fc, cc

    def calculate_contig(self, contig: str) -> MIContigResult:
        cr = MIContigResult(contig)

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

                k = (contig, int(lookup[1]), int(lookup[2]))

                # Check to make sure call is present in TRF BED file, if it is specified
                if self._loci_file and self._loci_dict and not self.get_loci_overlapping(*k):
                    continue

                # noinspection PyTypeChecker
                if self.should_exclude_locus(*k):
                    continue

                cr.seen_locus(*k)

                # Check to make sure call is present in all trio individuals
                if lookup not in mother_calls or lookup not in father_calls:
                    continue

                child_calls = locus_data[6:8]

                if "." in child_calls:
                    # Failed call
                    continue

                cr.append(MILocusData(
                    contig=contig,
                    start=k[1],
                    end=k[2],
                    motif=lookup[3],

                    child_gt=int_tuple(child_calls),
                    mother_gt=mother_calls[lookup],
                    father_gt=father_calls[lookup],
                ))

        return cr
