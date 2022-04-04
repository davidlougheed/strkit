import sys

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional, Tuple, Union

from .result import MIContigResult, MIResult
from ..utils import cis_overlap

__all__ = [
    "BaseCalculator",
]


# noinspection PyUnusedLocal
class BaseCalculator(ABC):
    def __init__(
            self,
            child_call_file: Path,
            mother_call_file: Path,
            father_call_file: Path,

            child_id: Optional[str] = None,
            mother_id: Optional[str] = None,
            father_id: Optional[str] = None,

            loci_file: Optional[str] = None,

            widen: float = 0,

            debug: bool = False,
    ):
        self._child_call_file: Path = child_call_file
        self._mother_call_file: Path = mother_call_file
        self._father_call_file: Path = father_call_file

        self._child_id: Optional[str] = child_id
        self._mother_id: Optional[str] = mother_id
        self._father_id: Optional[str] = father_id

        self._loci_file: Optional[str] = loci_file
        self._loci_dict = self._make_loci_dict()

        self._decimal_threshold: float = 0.5
        self._widen: float = widen

        self._debug: bool = debug

    def _make_loci_dict(self) -> dict:
        if not self._loci_file:
            return {}

        with open(self._loci_file, "r") as lf:
            return {
                tuple(d[:3]): d[3:]
                for d in (line.strip().split("\t") for line in lf)
            }

    @abstractmethod
    def _get_sample_contigs(self, include_sex_chromosomes: bool = False) -> Tuple[set, set, set]:
        return set(), set(), set()

    def get_trio_contigs(self, include_sex_chromosomes: bool = False) -> set:
        mc, fc, cc = self._get_sample_contigs()

        contig_set = mc.intersection(fc).intersection(cc)

        if include_sex_chromosomes:
            if "Y" in cc:
                contig_set = contig_set.union({"X", "Y"})
            elif "chrY" in cc:
                contig_set = contig_set.union({"chrX", "chrY"})
            elif "X" in cc:
                contig_set = contig_set.union({"X"})
            elif "chrX" in cc:
                contig_set = contig_set.union({"chrX"})

        return contig_set

    def gts_respect_mi(
            self,

            c_gt: Union[Tuple[int, ...], Tuple[float, ...]],
            m_gt: Union[Tuple[int, ...], Tuple[float, ...]],
            f_gt: Union[Tuple[int, ...], Tuple[float, ...]],

            c_gt_ci: Optional[Tuple] = None,
            m_gt_ci: Optional[Tuple] = None,
            f_gt_ci: Optional[Tuple] = None,

            decimal: bool = False,

            widen: Optional[float] = None,
    ) -> Tuple[bool, Optional[bool]]:
        # First hypothesis: first allele from mother, second from father
        # Second hypothesis: first allele from father, first from mother
        respects_mi_strict = any((
            (c_gt[0] in m_gt and c_gt[1] in f_gt),
            (c_gt[0] in f_gt and c_gt[1] in m_gt),
        ))

        if decimal:
            t = self._decimal_threshold

            respects_mi_strict = any((
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

        respects_mi_ci = None

        if c_gt_ci is not None and m_gt_ci is not None and f_gt_ci is not None:
            _widen = self._widen if widen is None else widen

            m_gt_ci_0 = (m_gt_ci[0][0] - (m_gt_ci[0][0] * _widen), m_gt_ci[0][1] + (m_gt_ci[0][1] * _widen))
            m_gt_ci_1 = (m_gt_ci[1][0] - (m_gt_ci[1][0] * _widen), m_gt_ci[1][1] + (m_gt_ci[1][1] * _widen))
            f_gt_ci_0 = (f_gt_ci[0][0] - (f_gt_ci[0][0] * _widen), f_gt_ci[0][1] + (f_gt_ci[0][1] * _widen))
            f_gt_ci_1 = (f_gt_ci[1][0] - (f_gt_ci[1][0] * _widen), f_gt_ci[1][1] + (f_gt_ci[1][1] * _widen))

            respects_mi_ci = any((
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

        return respects_mi_strict, respects_mi_ci

    @abstractmethod
    def calculate_contig(self, contig: str) -> MIContigResult:
        return MIContigResult()

    def calculate(self, included_contigs: set) -> Optional[MIResult]:
        res = 0
        res_95_ci = 0
        res_99_ci = 0
        n_total = 0

        contig_results = []
        non_matching = []

        for contig_result in map(self.calculate_contig, included_contigs):
            contig_results.append(contig_result)
            r, nm = contig_result.get_sums_and_non_matching()
            value, value_95_ci, value_99_ci = r
            res += value
            res_95_ci = None if value_95_ci is None else (res_95_ci + value_95_ci)
            res_99_ci = None if value_99_ci is None else (res_99_ci + value_99_ci)
            n_total += len(contig_result)
            non_matching.extend(nm)

        if n_total == 0:
            sys.stderr.write("Warning: no common loci found\n")
            return None

        res /= n_total
        res_95_ci = None if res_95_ci is None else (res_95_ci / n_total)
        res_99_ci = None if res_99_ci is None else (res_99_ci / n_total)

        return MIResult(res, res_95_ci, res_99_ci, contig_results, non_matching, self._widen)
