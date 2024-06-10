from __future__ import annotations

import itertools
import logging
import math
import numpy as np
import sys
import scipy.stats as sst

from statistics import mean
from statsmodels.stats.multitest import multipletests

from strkit.constants import CHROMOSOMES
from strkit.json import json, dumps_indented
from strkit.logger import get_main_logger
from strkit.utils import cat_strs, cis_overlap

from typing import Generator, Iterable, Optional, Union

__all__ = [
    "MILocusData",
    "MIContigResult",
    "MIResult",
]


OptionalReadCounts = Optional[tuple[tuple[int, ...], tuple[int, ...]]]

INHERITANCE_CONFIGS = (
    (0, 1, 0, 0),  # allele 0 maternal, allele 1 paternal, maternal 0 allele, paternal 0 allele
    (0, 1, 0, 1),  # allele 0 maternal, allele 1 paternal, maternal 0 allele, paternal 1 allele
    (0, 1, 1, 0),  # allele 0 maternal, allele 1 paternal, maternal 1 allele, paternal 0 allele
    (0, 1, 1, 1),  # allele 0 maternal, allele 1 paternal, maternal 1 allele, paternal 1 allele

    (1, 0, 0, 0),  # allele 1 maternal, allele 0 paternal, maternal 0 allele, paternal 0 allele
    (1, 0, 0, 1),  # allele 1 maternal, allele 0 paternal, maternal 0 allele, paternal 1 allele
    (1, 0, 1, 0),  # allele 1 maternal, allele 0 paternal, maternal 1 allele, paternal 0 allele
    (1, 0, 1, 1),  # allele 1 maternal, allele 0 paternal, maternal 1 allele, paternal 1 allele
)


class MILocusData:
    def __init__(
        self,
        contig: str,
        start: int,
        end: int,
        motif: str,

        child_gt, mother_gt, father_gt,
        child_gt_95_ci=None, mother_gt_95_ci=None, father_gt_95_ci=None,
        child_gt_99_ci=None, mother_gt_99_ci=None, father_gt_99_ci=None,

        child_read_counts: OptionalReadCounts = None,
        mother_read_counts: OptionalReadCounts = None,
        father_read_counts: OptionalReadCounts = None,

        reference_copies=None,
        decimal: bool = False,
        widen: float = 0,

        test_to_perform: str = "none",

        logger: Optional[logging.Logger] = None,
    ):
        self._contig = contig
        self._start = start
        self._end = end
        self._motif = motif

        self._child_gt = child_gt
        self._mother_gt = mother_gt
        self._father_gt = father_gt
        self._child_gt_95_ci = child_gt_95_ci
        self._mother_gt_95_ci = mother_gt_95_ci
        self._father_gt_95_ci = father_gt_95_ci
        self._child_gt_99_ci = child_gt_99_ci
        self._mother_gt_99_ci = mother_gt_99_ci
        self._father_gt_99_ci = father_gt_99_ci

        self._child_read_counts: OptionalReadCounts = child_read_counts
        self._child_read_counts_flattened: Optional[np.ndarray] = (
            np.array([z for y in child_read_counts for z in y]) if child_read_counts else None)
        self._mother_read_counts: OptionalReadCounts = mother_read_counts
        self._father_read_counts: OptionalReadCounts = father_read_counts

        self._reference_copies = reference_copies

        self._decimal = decimal
        self._decimal_threshold: float = 0.5
        self._widen: float = widen

        self._p_value = None
        self._adj_p_value = None
        if test_to_perform != "none":
            self._p_value = self.de_novo_test(test_to_perform)

        self._logger = logger or get_main_logger()

    @property
    def contig(self) -> str:
        return self._contig

    @property
    def start(self) -> int:
        return self._start

    @property
    def end(self) -> int:
        return self._end

    @property
    def motif(self) -> str:
        return self._motif

    @property
    def locus_str_data(self):
        return self._contig, str(self._start), str(self._end), self._motif

    @staticmethod
    def _res_str(res: Union[float, int]) -> str:
        return f"{res:.1f}" if isinstance(res, float) else str(res)

    @staticmethod
    def _gt_str(gt):
        return "|".join(map(MILocusData._res_str, gt)) if gt is not None else ""

    @staticmethod
    def _gt_ci_str(gt_ci):
        return "|".join(",".join(map(MILocusData._res_str, ci)) for ci in gt_ci) if gt_ci else "",

    @property
    def child_gt_str(self):
        return MILocusData._gt_str(self._child_gt)

    @property
    def child_gt_95_ci_str(self):
        return MILocusData._gt_str(self._child_gt_95_ci)

    @property
    def child_gt_99_ci_str(self):
        return MILocusData._gt_str(self._child_gt_99_ci)

    @property
    def mother_gt_str(self):
        return MILocusData._gt_str(self._mother_gt)

    @property
    def mother_gt_95_ci_str(self):
        return MILocusData._gt_str(self._mother_gt_95_ci)

    @property
    def mother_gt_99_ci_str(self):
        return MILocusData._gt_str(self._mother_gt_99_ci)

    @property
    def father_gt_str(self):
        return MILocusData._gt_str(self._father_gt)

    @property
    def father_gt_95_ci_str(self):
        return MILocusData._gt_str(self._father_gt_95_ci)

    @property
    def father_gt_99_ci_str(self):
        return MILocusData._gt_str(self._father_gt_99_ci)

    @property
    def reference_copies(self) -> Optional[int]:
        return self._reference_copies

    @property
    def p_value(self) -> Optional[float]:
        return self._p_value

    @property
    def adj_p_value(self) -> Optional[float]:
        return self._adj_p_value

    @adj_p_value.setter
    def adj_p_value(self, value: Optional[float]):
        try:
            assert value is None or 0 <= value <= 1
        except AssertionError as e:
            self._logger.error(f"Encountered unexpected value: {value}")
            raise e
        self._adj_p_value = value

    @staticmethod
    def _respects_strict_ci(c_gt, m_gt, f_gt) -> bool:
        # First hypothesis: first allele from mother, second from father
        # Second hypothesis: first allele from father, first from mother
        return any((
            (c_gt[0] in m_gt and c_gt[1] in f_gt),
            (c_gt[0] in f_gt and c_gt[1] in m_gt),
        ))

    def _respects_decimal_ci(self, c_gt, m_gt, f_gt) -> bool:
        t = self._decimal_threshold

        return any((
            # First hypothesis: first allele from mother, second from father
            abs(c_gt[0] - m_gt[0]) < t and abs(c_gt[1] - f_gt[0]) < t,
            abs(c_gt[0] - m_gt[0]) < t and abs(c_gt[1] - f_gt[1]) < t,
            abs(c_gt[0] - m_gt[1]) < t and abs(c_gt[1] - f_gt[0]) < t,
            abs(c_gt[0] - m_gt[1]) < t and abs(c_gt[1] - f_gt[1]) < t,

            abs(c_gt[1] - m_gt[0]) < t and abs(c_gt[0] - f_gt[0]) < t,
            abs(c_gt[1] - m_gt[0]) < t and abs(c_gt[0] - f_gt[1]) < t,
            abs(c_gt[1] - m_gt[1]) < t and abs(c_gt[0] - f_gt[0]) < t,
            abs(c_gt[1] - m_gt[1]) < t and abs(c_gt[0] - f_gt[1]) < t,
        ))

    def _respects_mi_ci(self, c_gt_ci, m_gt_ci, f_gt_ci, widen: float) -> Optional[bool]:
        for x in (c_gt_ci, m_gt_ci, f_gt_ci):
            if x is None:
                return None
            for y in x:
                if y is None:
                    return None

        try:
            m_gt_ci_0 = (m_gt_ci[0][0] - (m_gt_ci[0][0] * widen), m_gt_ci[0][1] + (m_gt_ci[0][1] * widen))
            m_gt_ci_1 = (m_gt_ci[1][0] - (m_gt_ci[1][0] * widen), m_gt_ci[1][1] + (m_gt_ci[1][1] * widen))
        except (IndexError, TypeError) as e:
            self._logger.error(f"Encountered invalid maternal confidence intervals: {m_gt_ci} ({e})")
            return None

        try:
            f_gt_ci_0 = (f_gt_ci[0][0] - (f_gt_ci[0][0] * widen), f_gt_ci[0][1] + (f_gt_ci[0][1] * widen))
            f_gt_ci_1 = (f_gt_ci[1][0] - (f_gt_ci[1][0] * widen), f_gt_ci[1][1] + (f_gt_ci[1][1] * widen))
        except (IndexError, TypeError) as e:
            self._logger.error(f"Encountered invalid paternal confidence intervals: {f_gt_ci} ({e})")
            return None

        try:
            return any((
                # First hypothesis: first allele from mother, second from father
                (cis_overlap(c_gt_ci[0], m_gt_ci_0) and cis_overlap(c_gt_ci[1], f_gt_ci_0)),
                (cis_overlap(c_gt_ci[0], m_gt_ci_0) and cis_overlap(c_gt_ci[1], f_gt_ci_1)),
                (cis_overlap(c_gt_ci[0], m_gt_ci_1) and cis_overlap(c_gt_ci[1], f_gt_ci_0)),
                (cis_overlap(c_gt_ci[0], m_gt_ci_1) and cis_overlap(c_gt_ci[1], f_gt_ci_1)),

                # Second hypothesis: first allele from father, first from mother
                (cis_overlap(c_gt_ci[1], m_gt_ci_0) and cis_overlap(c_gt_ci[0], f_gt_ci_0)),
                (cis_overlap(c_gt_ci[1], m_gt_ci_0) and cis_overlap(c_gt_ci[0], f_gt_ci_1)),
                (cis_overlap(c_gt_ci[1], m_gt_ci_1) and cis_overlap(c_gt_ci[0], f_gt_ci_0)),
                (cis_overlap(c_gt_ci[1], m_gt_ci_1) and cis_overlap(c_gt_ci[0], f_gt_ci_1)),
            ))
        except (IndexError, TypeError) as e:
            self._logger.error(f"Encountered invalid child confidence intervals: {c_gt_ci} ({e})")
            return None

    def respects_mi(self, widen: Optional[float] = None) -> tuple[bool, bool, Optional[bool], Optional[bool]]:
        fn = self._respects_decimal_ci if self._decimal else MILocusData._respects_strict_ci
        respects_mi_strict = fn(self._child_gt, self._mother_gt, self._father_gt)

        respects_mi_pm1 = self._respects_mi_ci(
            ((self._child_gt[0] - 1, self._child_gt[0] + 1), (self._child_gt[1] - 1, self._child_gt[1] + 1)),
            ((self._mother_gt[0], self._mother_gt[0]), (self._mother_gt[1], self._mother_gt[1])),
            ((self._father_gt[0], self._father_gt[0]), (self._father_gt[1], self._father_gt[1])),
            widen=0.0,
        )

        respects_mi_95_ci = self._respects_mi_ci(
            self._child_gt_95_ci, self._mother_gt_95_ci, self._father_gt_95_ci,
            widen=self._widen if widen is None else widen)

        respects_mi_99_ci = self._respects_mi_ci(
            self._child_gt_99_ci, self._mother_gt_99_ci, self._father_gt_99_ci,
            widen=self._widen if widen is None else widen)

        return respects_mi_strict, respects_mi_pm1, respects_mi_95_ci, respects_mi_99_ci

    def de_novo_test(self, test: str) -> Optional[float]:
        test_res = []

        cr_all = self._child_read_counts_flattened
        cr_all_set = set(cr_all)

        for config in INHERITANCE_CONFIGS[:4]:  # Don't need child allele configs here
            # TODO: I'm pretty sure this can be more numpy-ified / sped up

            mat_reads = self._mother_read_counts[config[2]]
            pat_reads = self._father_read_counts[config[3]]

            parent_all = np.array((*mat_reads, *pat_reads))

            if len(cr_all_set) == 1 and cr_all_set == set(parent_all) and test.startswith("wmw"):
                return None  # Complete overlap, NaN p-value for WMW test

            if not mat_reads or not pat_reads:
                w = []
                if not mat_reads:
                    w.append("maternal")
                if not pat_reads:
                    w.append("paternal")
                self._logger.warning(
                    f"{self.contig}:{self.start}-{self.end} - encountered empty readset for {' and '.join(w)} allele")
                sys.stderr.flush()
                return None

            if test == "x2":
                min_cn = min((*cr_all, *parent_all))
                max_cn = max((*cr_all, *parent_all))

                cr_all_resamp = cr_all

                n_bins = max_cn + 1 - min_cn
                child_freqs: list[int] = [0] * n_bins
                parent_freqs: list[int] = [0] * n_bins

                for c in cr_all_resamp:
                    child_freqs[c - min_cn] += 1

                for c in parent_all:
                    parent_freqs[c - min_cn] += 1

                empty_bins: set = {i for i in range(n_bins) if not child_freqs[i] and not parent_freqs[i]}
                child_freqs = [v for i, v in enumerate(child_freqs) if i not in empty_bins]
                parent_freqs = [v for i, v in enumerate(parent_freqs) if i not in empty_bins]

                # X2 test to check difference in distribution - don't apply Yates' correction
                p_val = sst.chi2_contingency(np.array([child_freqs, parent_freqs]), correction=False)[1]

            else:  # test == wmw or wmw-gt
                # Mann-Whitney U test; where we have H0 as the awkward 'equally likely to observe values
                # from the child that are less than observations from the parent as greater than' (something like that)
                # See https://www.tandfonline.com/doi/full/10.1080/00031305.2017.1305291
                # Same-variance requirements should hold almost all the time, unless there is suspected mosaicism
                # (and then X2 should be used instead)
                p_val = sst.mannwhitneyu(
                    cr_all,
                    parent_all,
                    alternative="greater" if test == "wmw-gt" else "two-sided",
                    use_continuity=False)[1]

            test_res.append(float(p_val))  # cast to make sure we're in native float type

        return max(test_res)

    def __iter__(self):  # for dict casting
        yield "contig", self._contig
        yield "start", self._start
        yield "end", self._end
        yield "motif", self._motif
        yield "child_gt", self._child_gt
        yield "child_gt_95_ci", self._child_gt_95_ci
        yield "child_read_counts", self._child_read_counts
        yield "mother_gt", self._mother_gt
        yield "mother_gt_95_ci", self._mother_gt_95_ci
        yield "mother_read_counts", self._mother_read_counts
        yield "father_gt", self._father_gt
        yield "father_gt_95_ci", self._father_gt_95_ci
        yield "father_read_counts", self._father_read_counts

        if self._p_value:
            yield "p", self._p_value
        if self._adj_p_value:
            yield "p_adj", self._adj_p_value

    def __str__(self):
        def _opt_ci_str(ci_str, level="95"):
            return "" if not ci_str else f" [{level}%: {ci_str}]"

        return (
            f"{self._contig}:{self._start}-{self._end} ("
            f"child: {self.child_gt_str}{_opt_ci_str(self.child_gt_95_ci_str)}, "
            f"mother: {self.mother_gt_str}{_opt_ci_str(self.mother_gt_95_ci_str)}, "
            f"father: {self.father_gt_str}{_opt_ci_str(self.father_gt_95_ci_str)}"
            f")"
        )


class MIContigResult:
    def __init__(self, includes_95_ci: bool = False, includes_99_ci: bool = False):
        self._loci_data: list[MILocusData] = []
        self._includes_95_ci: bool = includes_95_ci
        self._includes_99_ci: bool = includes_99_ci

        self._seen_loci: set[tuple[str, int, int]] = set()
        self._seen_loci_lengths: list[int] = []

    def seen_locus(self, contig: str, start: int, end: int):
        self._seen_loci.add((contig, start, end))

    @property
    def n_loci_seen(self) -> int:
        # Counter for number of loci which are in callsets, although call may have failed
        return len(self._seen_loci)

    @property
    def seen_loci_lengths(self) -> Generator[int, None, None]:
        # Lengths for binned counts
        yield from (t[2] - t[1] for t in self._seen_loci)

    def append(self, item: MILocusData):
        self._loci_data.append(item)

    def process_loci(
        self,
        calculate_non_matching: bool = True,
    ) -> tuple[tuple[int, int, Optional[int], Optional[int]], list[MILocusData]]:
        value = 0
        value_pm1 = 0
        value_95_ci = 0 if self._includes_95_ci else None
        value_99_ci = 0 if self._includes_99_ci else None
        non_matching = []

        for locus in self._loci_data:
            r = locus.respects_mi()

            if calculate_non_matching and not any(r[:2]):
                # TODO: Custom ability to choose level...
                non_matching.append(locus)

            r_strict, r_pm1, r_95ci, r_99ci = r
            value += r_strict
            value_pm1 += r_pm1
            if self._includes_95_ci and r_95ci is not None:
                value_95_ci += r_95ci
            if self._includes_99_ci and r_99ci is not None:
                value_99_ci += r_99ci

        return (value, value_pm1, value_95_ci, value_99_ci), non_matching

    def __bool__(self):
        return True  # Otherwise, it goes to len() which gives False if it's empty

    def __len__(self):
        return len(self._loci_data)

    def __iter__(self):
        yield from self._loci_data

    def __str__(self):
        return f"<MIContigResult #loci={len(self._loci_data)}>"


class MIResult:
    def __init__(
        self,
        mi_value: float,
        mi_value_pm1: float,
        mi_value_95_ci: Optional[float],
        mi_value_99_ci: Optional[float],
        contig_results: Iterable[MIContigResult],
        output_loci: list[MILocusData],
        widen: float = 0,
        test_to_perform: str = "none",
        sig_level: float = 0.05,
        mt_corr: str = "none",
        logger: Optional[logging.Logger] = None,
    ):
        self.mi_value: float = mi_value
        self.mi_value_pm1: float = mi_value_pm1
        self.mi_value_95_ci: Optional[float] = mi_value_95_ci
        self.mi_value_99_ci: Optional[float] = mi_value_99_ci
        self._contig_results: tuple[MIContigResult, ...] = tuple(contig_results)
        self._output_loci: list[MILocusData] = output_loci
        self.widen: float = widen
        self._test_to_perform: str = test_to_perform
        self._sig_level: float = sig_level
        self._mt_corr: str = mt_corr  # Method to use when correcting for multiple testing

        self._n_loci_seen: int = sum(cr.n_loci_seen for cr in self.contig_results)
        self._seen_loci_lengths: tuple[int, ...] = tuple(
            itertools.chain.from_iterable(cr.seen_loci_lengths for cr in self.contig_results))

        self._logger: logging.Logger = logger or get_main_logger()

    @property
    def contig_results(self) -> tuple[MIContigResult, ...]:
        return self._contig_results

    @property
    def output_loci(self) -> list[MILocusData]:
        return self._output_loci

    @property
    def test_to_perform(self) -> str:
        return self.test_to_perform

    @property
    def n_loci_seen(self) -> int:
        return self._n_loci_seen

    def correct_for_multiple_testing(self):
        if self._test_to_perform == "none":
            self._logger.warning("Cannot correct for multiple testing when test is not enabled")
            return

        p_values = [locus.p_value for cr in self.contig_results for locus in cr]
        not_tested = {i for i, p in enumerate(p_values) if p is None}

        loci_filt = [
            locus for i, locus in enumerate(locus for cr in self.contig_results for locus in cr)
            if i not in not_tested
        ]
        p_values_filt = np.fromiter((p for i, p in enumerate(p_values) if i not in not_tested), dtype=np.float64)

        if (mtm := self._mt_corr) != "none":
            rejected_hs, p_corr, _, _ = multipletests(
                p_values_filt,
                alpha=self._sig_level,
                method=mtm)  # Correct for multiple testing effects using the specified method
        else:  # No MT correction; just use p-values as-is
            p_corr = p_values_filt
            rejected_hs = (p_corr < self._sig_level)

        new_output_loci = []

        for pc, rej_h, locus in filter(lambda x: x[1], zip(p_corr, rejected_hs, loci_filt)):
            locus.adj_p_value = pc
            new_output_loci.append(locus)

        self._output_loci = new_output_loci

    def as_csv_row(self, sep=",") -> str:
        return f"{self.mi_value}{sep}{self.mi_value_95_ci}{sep}{self.mi_value_99_ci}\n"

    @staticmethod
    def _res_str(res: Union[float, int]) -> str:
        return f"{res:.1f}" if isinstance(res, float) else str(res)

    def _std_err_for_bernoulli(self, val: float) -> float:
        return math.sqrt((val * (1 - val)) / self._n_loci_seen)

    def _estimate_ci_for_bernoulli(self, val: float) -> tuple[float, float]:
        std_err = self._std_err_for_bernoulli(val)
        # 1.96 standard deviations for 95% CI, by CLM
        return val - 1.96*std_err, val + 1.96*std_err

    def write_report_json(self, json_path: str, bin_width: int = 10):
        hist, _bins = self.calculate_histogram(bin_width=bin_width)

        # get 95% CI from standard error, which we estimate from Bernoulli distribution of resulting percentage

        obj = {
            "mi": {
                "val": self.mi_value,
                "est_std_err": self._std_err_for_bernoulli(self.mi_value),
                "est_95_ci": self._estimate_ci_for_bernoulli(self.mi_value),
            },
            "mi_pm1": {
                "val": self.mi_value_pm1,
                "est_std_err": self._std_err_for_bernoulli(self.mi_value_pm1),
                "est_95_ci": self._estimate_ci_for_bernoulli(self.mi_value_pm1),
            },
            "mi_95": {
                "val": self.mi_value_95_ci,
                "est_std_err": self._std_err_for_bernoulli(self.mi_value_95_ci),
                "est_95_ci": self._estimate_ci_for_bernoulli(self.mi_value_95_ci),
            } if self.mi_value_95_ci else None,
            "mi_99": {
                "val": self.mi_value_99_ci,
                "est_std_err": self._std_err_for_bernoulli(self.mi_value_99_ci),
                "est_95_ci": self._estimate_ci_for_bernoulli(self.mi_value_99_ci),
            } if self.mi_value_99_ci else None,

            "n_loci_trio_called": sum((len(lr) for lr in self.contig_results)),
            "n_loci_total": self._n_loci_seen,

            "test": self._test_to_perform,
            "hist": hist,
            "significant_loci": [dict(nm) for nm in self.sorted_output_loci],
        }

        if json_path == "stdout":
            sys.stdout.buffer.write(json.dumps(obj))
            sys.stdout.write("\n")
            sys.stdout.flush()
            return

        with open(json_path, "wb") as jf:
            jf.write(dumps_indented(obj))

    @staticmethod
    def _nm_sort_key_pos(locus: MILocusData):
        return CHROMOSOMES.index(locus.contig), locus.start, locus.motif

    def _nm_sort_key_test(self, locus: MILocusData):
        if self._mt_corr != "none":
            return locus.adj_p_value
        return locus.p_value

    @property
    def sorted_output_loci(self):
        return sorted(
            self._output_loci, key=self._nm_sort_key_test if self._test_to_perform != "none" else self._nm_sort_key_pos)

    def locus_tsv(self, sep="\t") -> str:
        res = ""

        for locus in self.sorted_output_loci:
            res += sep.join((
                *locus.locus_str_data,

                # Child genotype + CIs if available
                locus.child_gt_str,
                locus.child_gt_95_ci_str,
                locus.child_gt_99_ci_str,

                # Mother genotype + CIs if available
                locus.mother_gt_str,
                locus.mother_gt_95_ci_str,
                locus.mother_gt_99_ci_str,

                # Father genotype + CIs if available
                locus.father_gt_str,
                locus.father_gt_95_ci_str,
                locus.father_gt_99_ci_str,

                str(locus.reference_copies or ""),  # Reference number of copies (or blank depending on caller)

                *((str(locus.p_value),) if self._test_to_perform != "none" else ()),
                *((str(locus.adj_p_value),) if self._test_to_perform != "none" and self._mt_corr != "none" else ()),
            )) + "\n"
        return res

    def __str__(self):
        widen_str = "" if self.widen < 0.00001 else f"; widened {self.widen * 100:.1f}%"
        header = ["MI%", "MI% (±1)"]
        mi_vals = [self.mi_value, self.mi_value_pm1]

        if self.mi_value_95_ci:
            header.append(f"MI% (95% CI{widen_str})")
            mi_vals.append(self.mi_value_95_ci)

        if self.mi_value_99_ci:
            header.append(f"MI% (99% CI{widen_str})")
            mi_vals.append(self.mi_value_99_ci)

        header_str = cat_strs(h.ljust(14) for h in header)
        mi_vals_str = cat_strs(map(lambda m: f"{m*100:.2f}".ljust(14), mi_vals))

        return f"{header_str}\n{mi_vals_str}"

    def calculate_histogram(self, bin_width: int = 10) -> tuple[list[dict], list]:
        # TODO: Don't duplicate with calculate()

        loci: list[MILocusData] = []
        for cr in self._contig_results:
            loci.extend(list(cr))

        bins = np.arange(0, max((locus.end - locus.start) for locus in loci) + bin_width, bin_width).tolist()
        hist = []

        vals_strict_by_bin = [[] for _ in bins]
        vals_pm1_by_bin = [[] for _ in bins]
        vals_95_ci_by_bin = [[] for _ in bins]
        vals_99_ci_by_bin = [[] for _ in bins]
        bin_totals = [0] * len(bins)

        for locus_len in self._seen_loci_lengths:
            if (bin_idx := locus_len // bin_width) < len(bin_totals):
                bin_totals[bin_idx] += 1

        for locus in loci:
            r_strict, r_pm1, r_95ci, r_99ci = locus.respects_mi()

            locus_len = locus.end - locus.start
            locus_bin_idx = locus_len // bin_width

            vals_strict_by_bin[locus_bin_idx].append(int(r_strict))
            vals_pm1_by_bin[locus_bin_idx].append(int(r_pm1))
            if r_95ci is not None:
                vals_95_ci_by_bin[locus_bin_idx].append(int(r_95ci))
            if r_99ci is not None:
                vals_99_ci_by_bin[locus_bin_idx].append(int(r_99ci))

        for i in range(len(bins)):
            vsb = vals_strict_by_bin[i]
            vpm1b = vals_pm1_by_bin[i]
            v95b = vals_95_ci_by_bin[i] if vals_95_ci_by_bin else None
            v99b = vals_99_ci_by_bin[i] if vals_95_ci_by_bin else None

            hist.append({
                "bin": bins[i],
                "bin_count": len(vsb),
                "bin_total": bin_totals[i],
                "mi": mean(vsb) if vsb else None,
                "mi_pm1": mean(vpm1b) if vpm1b else None,
                "mi_95": mean(v95b) if v95b else None,
                "mi_99": mean(v99b) if v99b else None,
            })

        return hist, bins

    def histogram_text(self, bin_width: int = 10) -> str:
        hist, bins = self.calculate_histogram(bin_width)

        bins_str = "\t".join(f"{b}-{b+bin_width-1}" for b in bins)
        bin_count_str = "\t".join(str(b["bin_count"]) for b in hist)

        def _format_means(bin_vals):
            return "\t".join(f"{b*100:.2f}" if b else "-" for b in bin_vals)

        vals_95_ci_by_bin = [b["mi_95"] for b in hist]
        vals_99_ci_by_bin = [b["mi_99"] for b in hist]

        mi_strict_by_bin_vals_str = _format_means([b["mi"] for b in hist])
        mi_pm1_by_bin_vals_str = "\n" + _format_means([b["mi_pm1"] for b in hist])
        mi_95_by_bin_vals_str = (("\n" + _format_means(vals_95_ci_by_bin)) if any(vals_95_ci_by_bin) else "")
        mi_99_by_bin_vals_str = (("\n" + _format_means(vals_99_ci_by_bin)) if any(vals_99_ci_by_bin) else "")

        return (
            f"{bins_str}\n{bin_count_str}\n{mi_strict_by_bin_vals_str}{mi_pm1_by_bin_vals_str}{mi_95_by_bin_vals_str}"
            f"{mi_99_by_bin_vals_str}")
