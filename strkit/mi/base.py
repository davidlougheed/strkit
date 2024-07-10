from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Optional

from strkit.logger import get_main_logger
from .intervals import (
    LociDictOfDict,
    LociDictOfList,
    build_loci_dict_of_dict_from_file,
    build_loci_dict_of_list_from_file,
    overlapping_loci_dict_of_dict,
    overlapping_loci_dict_of_list,
)
from .result import MIContigResult, MIResult

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
        exclude_file: Optional[str] = None,

        widen: float = 0,

        test_to_perform: str = "none",
        sig_level: float = 0.05,
        mt_corr: str = "none",

        debug: bool = False,
        logger: Optional[logging.Logger] = None,
    ):
        self._child_call_file: Path = child_call_file
        self._mother_call_file: Path = mother_call_file
        self._father_call_file: Path = father_call_file

        self._child_id: Optional[str] = child_id
        self._mother_id: Optional[str] = mother_id
        self._father_id: Optional[str] = father_id

        self._loci_file: Optional[str] = loci_file
        self._loci_dict: LociDictOfDict = build_loci_dict_of_dict_from_file(loci_file)

        self._exclude_file: Optional[str] = exclude_file
        self._exclude_dict: LociDictOfList = build_loci_dict_of_list_from_file(exclude_file)

        self._decimal_threshold: float = 0.5
        self._widen: float = widen

        self._test_to_perform: str = test_to_perform
        self._sig_level: float = sig_level
        self._mt_corr: str = mt_corr

        self._debug: bool = debug
        self._logger: logging.Logger = logger or get_main_logger()

        self._cache: dict[str, Any] = {}

    @property
    def test_to_perform(self) -> str:
        return self._test_to_perform

    @property
    def sig_level(self) -> float:
        return self._sig_level

    @property
    def mt_corr(self) -> str:
        return self._mt_corr

    def get_loci_overlapping(self, contig: str, start: int, end: int) -> list[tuple[int, int, list[str]]]:
        return overlapping_loci_dict_of_dict(contig, start, end, self._loci_dict)

    def should_exclude_locus(self, contig: str, start: int, end: int) -> bool:
        return any(True for _ in overlapping_loci_dict_of_list(contig, start, end, self._exclude_dict))

    @abstractmethod
    def _get_sample_contigs(self, include_sex_chromosomes: bool = False) -> tuple[set, set, set]:
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

        if self._loci_dict:
            # Limit contig set to only contigs which are in the locus dictionary if one is specified.
            contig_set = contig_set.intersection({k[0] for k in self._loci_dict})

        return contig_set

    @abstractmethod
    def calculate_contig(self, contig: str) -> MIContigResult:
        return MIContigResult()

    def calculate(self, included_contigs: set) -> Optional[MIResult]:
        res: float = 0
        res_pm1: float = 0
        res_95_ci: Optional[float] = None
        res_99_ci: Optional[float] = None
        n_total: int = 0

        contig_results = []
        output_loci = []

        for contig_result in map(self.calculate_contig, included_contigs):
            contig_results.append(contig_result)
            r, nm = contig_result.process_loci(calculate_non_matching=self.test_to_perform == "none")
            value, value_pm1, value_95_ci, value_99_ci = r
            res += value
            res_pm1 += value_pm1
            res_95_ci = None if value_95_ci is None else ((res_95_ci or 0) + value_95_ci)
            res_99_ci = None if value_99_ci is None else ((res_99_ci or 0) + value_99_ci)
            n_total += len(contig_result)
            output_loci.extend(nm)

        if n_total == 0:
            self._logger.warning("No common loci found")
            return None

        res /= n_total
        res_pm1 /= n_total
        res_95_ci = None if res_95_ci is None else (res_95_ci / n_total)
        res_99_ci = None if res_99_ci is None else (res_99_ci / n_total)

        mi_res = MIResult(
            res,
            res_pm1,
            res_95_ci,
            res_99_ci,
            contig_results,
            output_loci,
            self._widen,
            self.test_to_perform,
            self.sig_level,
            self.mt_corr,
            logger=self._logger,
        )

        if self.test_to_perform != "none":
            mi_res.correct_for_multiple_testing()  # Also calculates new output loci

        return mi_res
