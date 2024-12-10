from __future__ import annotations

import logging
import uuid
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

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
    "SEX_CHROMOSOMES",
    "BaseCalculator",
]


SEX_CHROMOSOMES = {"chrX", "X", "chrY", "Y"}  # TODO: proper parametrization


# noinspection PyUnusedLocal
class BaseCalculator(ABC):
    def __init__(
        self,
        child_call_file: Path,
        mother_call_file: Path,
        father_call_file: Path,

        child_id: str | None = None,
        mother_id: str | None = None,
        father_id: str | None = None,

        loci_file: str | None = None,
        exclude_file: str | None = None,

        widen: float = 0,

        test_to_perform: str = "none",
        sig_level: float = 0.05,
        mt_corr: str = "none",
        only_phased: bool = False,

        debug: bool = False,
        logger: logging.Logger | None = None,
    ):
        self._debug: bool = debug
        self._logger: logging.Logger = logger or get_main_logger()

        self._child_call_file: Path = child_call_file
        self._mother_call_file: Path = mother_call_file
        self._father_call_file: Path = father_call_file

        self._child_id: str | None = child_id
        self._mother_id: str | None = mother_id
        self._father_id: str | None = father_id

        self._loci_file: str | None = loci_file
        self._loci_dict: LociDictOfDict = build_loci_dict_of_dict_from_file(loci_file)
        self._loci_dict_cache_key: str = str(uuid.uuid4())
        if self._loci_file is not None:
            self._logger.debug(
                "Built loci dict of size %d with contigs %s",
                sum(len(loc) for loc in self._loci_dict.values()),
                tuple(self._loci_dict.keys()),
            )

        self._exclude_file: str | None = exclude_file
        self._exclude_dict: LociDictOfList = build_loci_dict_of_list_from_file(exclude_file)
        if self._exclude_file is not None:
            self._logger.debug(
                "Built exclude dict of size %d with contigs %s",
                len(self._loci_dict),
                tuple(self._exclude_dict.keys()),
            )

        self._decimal_threshold: float = 0.5
        self._widen: float = widen

        self._test_to_perform: str = test_to_perform
        self._sig_level: float = sig_level
        self._mt_corr: str = mt_corr
        self._only_phased: bool = only_phased

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

    def get_loci_overlapping(
        self, contig: str, start: int, end: int, first_only: bool
    ) -> list[tuple[int, int, list[str]]]:
        return overlapping_loci_dict_of_dict(
            contig, start, end, self._loci_dict, first_only, dict_cache_key=self._loci_dict_cache_key
        )

    def should_exclude_locus(self, contig: str, start: int, end: int) -> bool:
        return any(True for _ in overlapping_loci_dict_of_list(contig, start, end, self._exclude_dict, True))

    def should_skip_locus(
        self, contig: str, start: int, end: int, cached_overlapping: list | None = None
    ) -> str | None:
        # Returns either a reason string (if yes) or None (=== no)

        # Check to make sure call is present in TRF BED file, if it is specified
        # Check to make sure the locus is not excluded via overlap with exclude BED

        if not self._loci_file or not self._loci_dict:
            return None

        if not (cached_overlapping or self.get_loci_overlapping(contig, start, end, True)):
            return "no overlapping loci"

        if self.should_exclude_locus(contig, start, end):
            return "should_exclude_locus returned True"

        return None

    @abstractmethod
    def _get_sample_contigs(self) -> tuple[set, set, set]:
        return set(), set(), set()

    def get_trio_contigs(self, include_sex_chromosomes: bool = False) -> set:
        mc, fc, cc = self._get_sample_contigs()

        contig_set = mc.intersection(fc).intersection(cc)

        if include_sex_chromosomes:  # TODO: proper parametrization
            if "Y" in cc:
                contig_set = contig_set.union({"X", "Y"})
            elif "chrY" in cc:
                contig_set = contig_set.union({"chrX", "chrY"})
            elif "X" in cc:
                contig_set = contig_set.union({"X"})
            elif "chrX" in cc:
                contig_set = contig_set.union({"chrX"})
        else:
            contig_set = contig_set.difference(SEX_CHROMOSOMES)

        if self._loci_dict:
            # Limit contig set to only contigs which are in the locus dictionary if one is specified.
            contig_set = contig_set.intersection(self._loci_dict.keys())

        self._logger.debug("Got %d intersection trio contigs", len(contig_set))

        return contig_set

    @abstractmethod
    def calculate_contig(self, contig: str) -> MIContigResult:
        return MIContigResult(contig)

    @staticmethod
    def _updated_mi_res(res: float | None, v: int | float | None) -> float | None:
        return None if v is None else ((res or 0) + v)

    def calculate(self, included_contigs: set) -> MIResult | None:
        # copy number
        res: float = 0
        res_pm1: float = 0
        res_95_ci: float | None = None
        res_99_ci: float | None = None
        # sequence
        res_seq: float | None = None
        res_sl: float | None = None
        res_sl_pm1: float | None = None

        n_total: int = 0

        contig_results = []
        output_loci = []

        for contig in sorted(included_contigs):
            self._logger.info("Processing contig %s", contig)

            contig_result = self.calculate_contig(contig)
            contig_results.append(contig_result)

            r, nm = contig_result.process_loci(calculate_non_matching=self.test_to_perform == "none")

            value_95_ci = r["ci_95"]
            value_99_ci = r["ci_99"]
            value_seq = r["seq"]
            value_sl = r["sl"]
            value_sl_pm1 = r["sl_pm1"]

            res += r["strict"]
            res_pm1 += r["pm1"]
            res_95_ci = self._updated_mi_res(res_95_ci, value_95_ci)
            res_99_ci = self._updated_mi_res(res_99_ci, value_99_ci)
            res_seq = self._updated_mi_res(res_seq, value_seq)
            res_sl = self._updated_mi_res(res_sl, value_sl)
            res_sl_pm1 = self._updated_mi_res(res_sl_pm1, value_sl_pm1)

            n_total += len(contig_result)
            output_loci.extend(nm)

            logger_fmt = "Finished processing contig %s; n_total=%d. Current value: %.2f%%, ±1: %.2f%%"
            logger_args = [contig_result.contig, n_total, res / n_total * 100, res_pm1 / n_total * 100]

            extras = (
                (res_95_ci, "95% CI"),
                (res_99_ci, "99% CI"),
                (res_seq, "seq"),
                (res_sl, "s.l."),
                (res_sl_pm1, "s.l.±1"),
            )

            for val, fmt_txt in extras:
                if val is not None:
                    logger_fmt += f", {fmt_txt}: %.2f%%"
                    logger_args.append(val / n_total * 100)

            self._logger.info(logger_fmt, *logger_args)

        if n_total == 0:
            self._logger.warning("No common loci found")
            return None

        res /= n_total
        res_pm1 /= n_total
        res_95_ci = None if res_95_ci is None else (res_95_ci / n_total)
        res_99_ci = None if res_99_ci is None else (res_99_ci / n_total)
        res_seq = None if res_seq is None else (res_seq / n_total)
        res_sl = None if res_sl is None else (res_sl / n_total)
        res_sl_pm1 = None if res_sl is None else (res_sl_pm1 / n_total)

        mi_res = MIResult(
            {
                "strict": res,
                "pm1": res_pm1,
                "ci_95": res_95_ci,
                "ci_99": res_99_ci,
                "seq": res_seq,
                "sl": res_sl,
                "sl_pm1": res_sl_pm1,
            },
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
