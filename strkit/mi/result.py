import json
import numpy as np
import sys

from typing import List, Iterable, Optional, Tuple, Union

from statistics import mean
from strkit.constants import CHROMOSOMES
from strkit.utils import cis_overlap

__all__ = [
    "MILocusData",
    "MIContigResult",
    "MIResult",
]


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

            reference_copies=None,
            decimal: bool = False,
            widen: float = 0):
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

        self._reference_copies = reference_copies

        self._decimal = decimal
        self._decimal_threshold: float = 0.5
        self._widen: float = widen

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

    @staticmethod
    def _respects_mi_ci(c_gt_ci, m_gt_ci, f_gt_ci, widen) -> Optional[bool]:
        if any(x is None for x in (c_gt_ci, m_gt_ci, f_gt_ci)):
            return None

        m_gt_ci_0 = (m_gt_ci[0][0] - (m_gt_ci[0][0] * widen), m_gt_ci[0][1] + (m_gt_ci[0][1] * widen))
        m_gt_ci_1 = (m_gt_ci[1][0] - (m_gt_ci[1][0] * widen), m_gt_ci[1][1] + (m_gt_ci[1][1] * widen))
        f_gt_ci_0 = (f_gt_ci[0][0] - (f_gt_ci[0][0] * widen), f_gt_ci[0][1] + (f_gt_ci[0][1] * widen))
        f_gt_ci_1 = (f_gt_ci[1][0] - (f_gt_ci[1][0] * widen), f_gt_ci[1][1] + (f_gt_ci[1][1] * widen))

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

    def respects_mi(self, widen: Optional[float] = None) -> Tuple[bool, Optional[bool], Optional[bool]]:
        fn = self._respects_decimal_ci if self._decimal else MILocusData._respects_strict_ci
        respects_mi_strict = fn(self._child_gt, self._mother_gt, self._father_gt)

        respects_mi_95_ci = MILocusData._respects_mi_ci(
            self._child_gt_95_ci, self._mother_gt_95_ci, self._father_gt_95_ci,
            widen=self._widen if widen is None else widen)

        respects_mi_99_ci = MILocusData._respects_mi_ci(
            self._child_gt_99_ci, self._mother_gt_99_ci, self._father_gt_99_ci,
            widen=self._widen if widen is None else widen)

        return respects_mi_strict, respects_mi_95_ci, respects_mi_99_ci

    def __iter__(self):  # for dict casting
        yield "contig", self._contig
        yield "start", self._start
        yield "end", self._end
        yield "motif", self._motif
        yield "child_gt", self._child_gt
        yield "child_gt_95_ci", self._child_gt_95_ci
        yield "mother_gt", self._mother_gt
        yield "mother_gt_95_ci", self._mother_gt_95_ci
        yield "father_gt", self._father_gt
        yield "father_gt_95_ci", self._father_gt_95_ci

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
        self._loci_data: List[MILocusData] = []
        self._includes_95_ci = includes_95_ci
        self._includes_99_ci = includes_99_ci

    def append(self, item: MILocusData):
        self._loci_data.append(item)

    def get_sums_and_non_matching(self) -> Tuple[Tuple[int, Optional[int],  Optional[int]], List[MILocusData]]:
        value = 0
        value_95_ci = 0 if self._includes_95_ci else None
        value_99_ci = 0 if self._includes_99_ci else None
        non_matching = []

        for locus in self._loci_data:
            r = locus.respects_mi()
            if not any(r[:2]):
                # TODO: Custom ability to choose level...
                non_matching.append(locus)
            value += r[0]
            if self._includes_95_ci and r[1] is not None:
                value_95_ci += r[1]
            if self._includes_99_ci and r[2] is not None:
                value_99_ci += r[2]

        return (value, value_95_ci, value_99_ci), non_matching

    def __len__(self):
        return len(self._loci_data)

    def __iter__(self):
        yield from self._loci_data

    def __str__(self):
        return f"<MIContigResult #loci={len(self._loci_data)}>"


class MIResult:
    def __init__(self,
                 mi_value: float,
                 mi_value_95_ci: Optional[float],
                 mi_value_99_ci: Optional[float],
                 contig_results: Iterable[MIContigResult],
                 non_matching: List[MILocusData],
                 widen: float = 0):
        self.mi_value = mi_value
        self.mi_value_95_ci = mi_value_95_ci
        self.mi_value_99_ci = mi_value_99_ci
        self._contig_results: Tuple[MIContigResult] = tuple(contig_results)
        self._non_matching: List[MILocusData] = non_matching
        self.widen = widen

    @property
    def contig_results(self):
        return self._contig_results

    @property
    def non_matching(self):
        return self._non_matching

    def as_csv_row(self, sep=","):
        return f"{self.mi_value}{sep}{self.mi_value_95_ci}{sep}{self.mi_value_99_ci}\n"

    @staticmethod
    def _res_str(res: Union[float, int]) -> str:
        return f"{res:.1f}" if isinstance(res, float) else str(res)

    def write_report_json(self, json_path: str, bin_width: int = 10):
        hist, _bins = self.calculate_histogram(bin_width=bin_width)

        obj = {
            "mi": self.mi_value,
            "mi_95": self.mi_value_95_ci,
            "mi_99": self.mi_value_99_ci,
            "non_matching": [dict(nm) for nm in self._non_matching],
            "hist": hist,
        }

        if json_path == "stdout":
            json.dump(obj, sys.stdout)
            sys.stdout.write("\n")
            sys.stdout.flush()
            return

        with open(json_path, "w") as jf:
            json.dump(obj, jf, indent=2)

    def non_matching_tsv(self, sep="\t") -> str:
        res = ""
        for nm in sorted(self._non_matching, key=lambda x: (CHROMOSOMES.index(x.contig), x.start, x.motif)):
            res += sep.join((
                *nm.locus_str_data,

                # Child genotype + CIs if available
                nm.child_gt_str,
                nm.child_gt_95_ci_str,
                nm.child_gt_99_ci_str,

                # Mother genotype + CIs if available
                nm.mother_gt_str,
                nm.mother_gt_95_ci_str,
                nm.mother_gt_99_ci_str,

                # Father genotype + CIs if available
                nm.father_gt_str,
                nm.father_gt_95_ci_str,
                nm.father_gt_99_ci_str,

                str(nm.reference_copies or ""),  # Reference number of copies (or blank depending on caller)
            )) + "\n"
        return res

    def __str__(self):
        widen_str = "" if self.widen < 0.00001 else f"; widened {self.widen * 100:.1f}%"
        header = ["MI %"]
        mi_vals = [self.mi_value]

        if self.mi_value_95_ci:
            header.append(f"MI % (95% CI{widen_str})")
            mi_vals.append(self.mi_value_95_ci)

        if self.mi_value_99_ci:
            header.append(f"MI % (99% CI{widen_str})")
            mi_vals.append(self.mi_value_99_ci)

        header_str = "\t".join(header)
        mi_vals_str = "\t".join(map(lambda m: f"{m*100:.2f}", mi_vals))

        return f"{header_str}\n{mi_vals_str}"

    def calculate_histogram(self, bin_width: int = 10) -> Tuple[List[dict], list]:
        # TODO: Don't duplicate with calculate()

        loci: List[MILocusData] = []
        for cr in self._contig_results:
            loci.extend(list(cr))

        bins = np.arange(0, max((locus.end - locus.start) for locus in loci) + bin_width, bin_width).tolist()
        hist = []

        vals_strict_by_bin = [[] for _ in bins]
        vals_95_ci_by_bin = [[] for _ in bins]
        vals_99_ci_by_bin = [[] for _ in bins]

        for locus in loci:
            r = locus.respects_mi()

            locus_len = locus.end - locus.start
            locus_bin_idx = locus_len // bin_width

            vals_strict_by_bin[locus_bin_idx].append(int(r[0]))
            if r[1] is not None:
                vals_95_ci_by_bin[locus_bin_idx].append(int(r[1]))
            if r[2] is not None:
                vals_99_ci_by_bin[locus_bin_idx].append(int(r[2]))

        for i in range(len(bins)):
            vsb = vals_strict_by_bin[i]
            v95b = vals_95_ci_by_bin[i] if vals_95_ci_by_bin else None
            v99b = vals_99_ci_by_bin[i] if vals_95_ci_by_bin else None

            hist.append({
                "bin": bins[i],
                "bin_count": len(vsb),
                "mi": mean(vsb) if vsb else None,
                "mi_95": mean(v95b) if v95b else None,
                "mi_99": mean(v99b) if v99b else None,
            })

        return hist, bins

    def histogram_text(self, bin_width: int = 10) -> str:
        hist, bins = self.calculate_histogram(bin_width)

        bins_str = "\t".join(f"{b}-{b+bin_width-1}" for b in bins)
        bin_count_str = "\t".join(str(b["bin_count"]) for b in hist)

        def _format_means(bin_vals):
            return "\t".join(f"{mean(b)*100:.2f}" if b else "-" for b in bin_vals)

        vals_95_ci_by_bin = [b["mi_95"] for b in hist]
        vals_99_ci_by_bin = [b["mi_99"] for b in hist]

        mi_strict_by_bin_vals_str = _format_means([b["mi"] for b in hist])
        mi_95_by_bin_vals_str = (("\n" + _format_means(vals_95_ci_by_bin)) if any(vals_95_ci_by_bin) else "")
        mi_99_by_bin_vals_str = (("\n" + _format_means(vals_99_ci_by_bin)) if any(vals_99_ci_by_bin) else "")

        return (
            f"{bins_str}\n{bin_count_str}\n{mi_strict_by_bin_vals_str}{mi_95_by_bin_vals_str}{mi_99_by_bin_vals_str}")
