from __future__ import annotations

from typing import Optional

from .base import BaseCalculator
from .result import MILocusData, MIContigResult

__all__ = [
    "StraglrCalculator",
]


class StraglrCalculator(BaseCalculator):
    @staticmethod
    def get_contigs_from_fh(fh) -> set:
        return {ls[0] for ls in (line.split("\t") for line in fh if not line.startswith("#"))}

    def make_calls_dict(self, ph, contig, cr: Optional[MIContigResult] = None):
        # For reference, dicts are ordered in Python 3.7+ (guaranteed)

        calls = {}

        for pv in ph:
            if pv.startswith("#"):
                continue

            line = pv.strip().split("\t")

            if line[0] != contig:
                continue

            locus = tuple(line[:3])

            orig_motif = self._loci_dict.get(locus)
            orig_motif = orig_motif[-1] if orig_motif else None

            if not orig_motif:
                continue

            if self.should_exclude_locus(locus):
                continue

            if cr:
                cr.seen_locus(locus[0], int(locus[1]), int(locus[2]))

            # Transform the genotypes into something that is consistent across individuals,
            # using the file with the list of loci.
            gt_fact = len(line[3]) / len(orig_motif)

            gt = tuple(float(g.split("(")[0]) * gt_fact for g in line[4].split(";"))
            if len(gt) == 1:  # If it's homozygous, expand it out to length 2
                gt = gt + gt

            calls[locus + (orig_motif,)] = gt

        return calls

    def _get_sample_contigs(self, include_sex_chromosomes: bool = False) -> tuple[set, set, set]:
        with open(self._mother_call_file, "r") as mvf, open(self._father_call_file, "r") as fvf, \
                open(self._child_call_file, "r") as cvf:

            mc = self.get_contigs_from_fh(mvf)
            fc = self.get_contigs_from_fh(fvf)
            cc = self.get_contigs_from_fh(cvf)

            return mc, fc, cc

    def calculate_contig(self, contig: str):
        cr = MIContigResult()

        with open(self._mother_call_file, "r") as mh:
            mother_calls = self.make_calls_dict(mh, contig)

        with open(self._father_call_file, "r") as fh:
            father_calls = self.make_calls_dict(fh, contig)

        with open(self._child_call_file, "r") as ch:
            child_calls = self.make_calls_dict(ch, contig, cr)

        for locus_data, c_gt in child_calls.items():
            # Check to make sure call is present in all trio individuals
            if locus_data not in mother_calls or locus_data not in father_calls:
                continue

            cr.append(MILocusData(
                contig=locus_data[0],
                start=int(locus_data[1]),
                end=int(locus_data[2]),
                motif=locus_data[3],

                child_gt=c_gt,
                mother_gt=mother_calls[locus_data],
                father_gt=father_calls[locus_data],

                decimal=True,
            ))

        return cr
