from typing import Tuple

from .base import BaseCalculator

__all__ = [
    "StraglrCalculator",
    "StraglrReCallCalculator",
]


class StraglrCalculator(BaseCalculator):
    @staticmethod
    def get_contigs_from_fh(fh) -> set:
        return {ls[0] for ls in (line.split("\t") for line in fh if not line.startswith("#"))}

    def make_calls_dict(self, ph, contig):
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
            if not orig_motif:
                continue

            # Transform the genotypes into something that is consistent across individuals,
            # using the file with the list of loci.
            gt_fact = len(line[3]) / len(orig_motif)

            gt = tuple(float(g.split("(")[0]) * gt_fact for g in line[4].split(";"))
            if len(gt) == 1:  # If it's homozygous, expand it out to length 2
                gt = gt + gt

            calls[locus + (orig_motif,)] = gt

        return calls

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
            child_calls = self.make_calls_dict(ch, contig)

        for locus_data, c_gt in child_calls.items():
            if locus_data[0] != contig:
                continue

            # Check to make sure call is present in all trio individuals
            if locus_data not in mother_calls or locus_data not in father_calls:
                continue

            n_loci += 1
            # TODO: Divide n_loci by number of reads otherwise it'll be too large

            m_gt = mother_calls[locus_data]
            f_gt = father_calls[locus_data]

            respects_mi_strict, _ = self.gts_respect_mi(c_gt, m_gt, f_gt, decimal=True)
            if respects_mi_strict:
                # Mendelian inheritance upheld for this locus - strict
                value += 1
            else:
                non_matching.append((
                    *locus_data,

                    c_gt, "",
                    m_gt, "",
                    f_gt, "",

                    "",
                ))

        return value, None, n_loci, non_matching


class StraglrReCallCalculator(BaseCalculator):
    @staticmethod
    def get_contigs_from_fh(fh) -> set:
        return {ls[0] for ls in (line.split("\t") for line in fh if not line.startswith("#"))}

    def make_calls_dict(self, ph, contig):
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
            if not orig_motif:
                continue

            if "." in line[6:8]:
                continue

            # Transform the genotypes into something that is consistent across individuals,
            # using the file with the list of loci.
            # Round it to the nearest first decimal place, since that is what Straglr calls.

            gt_fact = len(line[3]) / len(orig_motif)

            def _to_tenth(x: str):
                return round(float(x) * gt_fact * 10) / 10

            gt = tuple(map(_to_tenth, line[6:8]))
            gt_95_ci = tuple(tuple(map(_to_tenth, ci.split(","))) for ci in line[8:10])
            gt_99_ci = tuple(tuple(map(_to_tenth, ci.split(","))) for ci in line[10:12])

            calls[locus + (orig_motif,)] = (gt, gt_95_ci, gt_99_ci)

        return calls

    def _get_sample_contigs(self, include_sex_chromosomes: bool = False) -> Tuple[set, set, set]:
        with open(self._mother_call_file, "r") as mvf, open(self._father_call_file, "r") as fvf, \
                open(self._child_call_file, "r") as cvf:

            mc = self.get_contigs_from_fh(mvf)
            fc = self.get_contigs_from_fh(fvf)
            cc = self.get_contigs_from_fh(cvf)

            return mc, fc, cc

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
            child_calls = self.make_calls_dict(ch, contig)

        for locus_data, c_gt_and_cis in child_calls.items():
            if locus_data[0] != contig:
                continue

            # Check to make sure call is present in all trio individuals
            if locus_data not in mother_calls or locus_data not in father_calls:
                continue

            c_gt, c_gt_95_ci, _ = c_gt_and_cis

            n_loci += 1

            m_gt, m_gt_95_ci, _ = mother_calls[locus_data]
            f_gt, f_gt_95_ci, _ = father_calls[locus_data]

            respects_mi_strict, respects_mi_95_ci = self.gts_respect_mi(
                c_gt=c_gt, m_gt=m_gt, f_gt=f_gt,
                c_gt_ci=c_gt_95_ci, m_gt_ci=m_gt_95_ci, f_gt_ci=f_gt_95_ci,
                decimal=True)

            if respects_mi_strict:
                # Mendelian inheritance upheld for this locus - strict
                value += 1

            if respects_mi_95_ci:
                # Mendelian inheritance upheld for this locus - within 95% CI from TG2MM
                value_95_ci += 1
            else:
                non_matching.append((
                    *locus_data,

                    c_gt, c_gt_95_ci,
                    m_gt, m_gt_95_ci,
                    f_gt, f_gt_95_ci,

                    "",
                ))

        return value, value_95_ci, n_loci, non_matching
