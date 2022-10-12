from __future__ import annotations

import numpy as np

from strkit.json import json

from .base import BaseCalculator
from .result import MIContigResult, MILocusData
from ..utils import int_tuple, float_tuple, parse_cis

__all__ = [
    "StrKitCalculator",
    "StrKitJSONCalculator",
]


class StrKitCalculator(BaseCalculator):
    fractional = False

    @staticmethod
    def get_contigs_from_fh(fh) -> set[str]:
        return {ls[0] for ls in (line.split("\t") for line in fh if not line.startswith("#"))}

    def _get_sample_contigs(self, include_sex_chromosomes: bool = False) -> tuple[set, set, set]:
        with open(self._mother_call_file, "r") as mvf:
            mc = self.get_contigs_from_fh(mvf)
        with open(self._father_call_file, "r") as fvf:
            fc = self.get_contigs_from_fh(fvf)
        with open(self._child_call_file, "r") as cvf:
            cc = self.get_contigs_from_fh(cvf)
        return mc, fc, cc

    @staticmethod
    def make_calls_dict(ph, contig):
        tuple_conv = float_tuple if StrKitCalculator.fractional else int_tuple
        dtype = float if StrKitCalculator.fractional else int

        return {
            tuple(line[:4]): (
                tuple_conv(line[-2].split("|")),
                parse_cis(line[-1].split("|"), dtype=dtype),
                None  # parse_cis(line[-1:].split("|")),
            )
            for line in (pv.strip().split("\t") for pv in ph)
            if line[0] == contig and "." not in line[-2].split("|")
        }

    def calculate_contig(self, contig: str) -> MIContigResult:
        tuple_conv = float_tuple if StrKitCalculator.fractional else int_tuple
        dtype = float if StrKitCalculator.fractional else int

        cr = MIContigResult(includes_95_ci=True)

        with open(self._mother_call_file) as mh:
            mother_calls = self.make_calls_dict(mh, contig)

        with open(self._father_call_file) as fh:
            father_calls = self.make_calls_dict(fh, contig)

        with open(self._child_call_file) as ch:
            for cv in ch:
                locus_data = cv.strip().split("\t")

                if locus_data[0] != contig:
                    continue

                lookup = tuple(locus_data[:4])

                # Check to make sure call is present in TRF BED file, if it is specified
                if self._loci_file and self._loci_dict and lookup[:3] not in self._loci_dict:
                    continue

                # Check to make sure call is present in all trio individuals
                if lookup not in mother_calls or lookup not in father_calls:
                    continue

                m_gt, m_gt_95_ci, _ = mother_calls[lookup]
                f_gt, f_gt_95_ci, _ = father_calls[lookup]

                calls = locus_data[-2].split("|")

                if "." in calls:
                    # Failed call
                    continue

                cr.append(MILocusData(
                    contig=lookup[0],
                    start=int(lookup[1]),
                    end=int(lookup[2]),
                    motif=lookup[3],

                    child_gt=tuple_conv(calls),
                    mother_gt=m_gt,
                    father_gt=f_gt,

                    child_gt_95_ci=parse_cis(locus_data[-1].split("|"), dtype=dtype),
                    mother_gt_95_ci=m_gt_95_ci,
                    father_gt_95_ci=f_gt_95_ci,

                    # child_gt_99_ci=parse_cis(locus_data[-1:].split("|")),
                    # mother_gt_99_ci=m_gt_99_ci,
                    # father_gt_99_ci=f_gt_99_ci,

                    reference_copies=dtype(locus_data[4]),

                    decimal=StrKitCalculator.fractional,
                ))

        return cr


class StrKitJSONCalculator(BaseCalculator):
    @staticmethod
    def get_contigs_from_fh(fh) -> set:
        report = json.loads(fh.read())
        return {res["contig"] for res in report["results"]}

    def _get_sample_contigs(self, include_sex_chromosomes: bool = False) -> tuple[set, set, set]:
        with open(self._mother_call_file, "r") as mvf:
            mc = self.get_contigs_from_fh(mvf)
        with open(self._father_call_file, "r") as fvf:
            fc = self.get_contigs_from_fh(fvf)
        with open(self._child_call_file, "r") as cvf:
            cc = self.get_contigs_from_fh(cvf)
        return mc, fc, cc

    @staticmethod
    def get_read_counts(res: dict, dtype=int):
        # TODO: This only works with diploids...

        read_cns = [r["cn"] for r in res["reads"].values()]
        read_peaks = [r["p"] for r in res["reads"].values()]

        n = res["peaks"]["modal_n"]

        if n < 2 or len(set(res["call"])) == 1:
            rcs = np.array(read_cns, dtype=dtype)
            np.random.shuffle(rcs)  # TODO: seed shuffle
            part = rcs.shape[0] // 2
            return tuple(rcs[:part].tolist()), tuple(rcs[part:].tolist())

        rc = []
        for _ in range(n):
            rc.append([])
        for r, cn in enumerate(read_cns):
            rc[read_peaks[r]].append(cn)
        return tuple(map(tuple, rc))

    @staticmethod
    def make_calls_dict(report: dict, contig):
        fractional = report["parameters"]["fractional"]
        tuple_conv = float_tuple if fractional else int_tuple
        dtype = float if fractional else int

        return {
            (res["contig"], res["start"], res["end"], res["motif"]): (
                tuple_conv(res["call"]),
                tuple(map(lambda x: tuple(map(dtype, x)), res["call_95_cis"])),
                None,  # Placeholder for 99% CI
                StrKitJSONCalculator.get_read_counts(res, dtype=dtype),
            )
            for res in report["results"]
            if res["contig"] == contig and res["call"] is not None
        }

    def calculate_contig(self, contig: str) -> MIContigResult:
        with open(self._child_call_file, "r") as ch:
            c_report = json.loads(ch.read())

        fractional = c_report["parameters"]["fractional"]

        tuple_conv = float_tuple if fractional else int_tuple
        dtype = float if fractional else int

        cr = MIContigResult(includes_95_ci=True)

        with open(self._mother_call_file) as mh:
            mother_data = self.make_calls_dict(json.loads(mh.read()), contig)

        with open(self._father_call_file) as fh:
            father_data = self.make_calls_dict(json.loads(fh.read()), contig)

        for res in c_report["results"]:
            if res["contig"] != contig:
                continue

            lookup = (res["contig"], res["start"], res["end"], res["motif"])
            # noinspection PyTypeChecker
            loci_lookup: tuple[str, str, str] = tuple(map(str, lookup[:3]))

            # Check to make sure call is present in TRF BED file, if it is specified
            if self._loci_file and self._loci_dict and loci_lookup not in self._loci_dict:
                continue

            if self.should_exclude_locus(loci_lookup):
                continue

            # Check to make sure call is present in all trio individuals
            if lookup not in mother_data or lookup not in father_data:
                continue

            m_gt, m_gt_95_ci, _, m_rcs = mother_data[lookup]
            f_gt, f_gt_95_ci, _, f_rcs = father_data[lookup]

            if res["call"] is None:
                # Failed call
                continue

            call = tuple_conv(res["call"])

            cr.append(MILocusData(
                contig=lookup[0],
                start=lookup[1],
                end=lookup[2],
                motif=lookup[3],

                child_gt=tuple_conv(call),
                mother_gt=m_gt,
                father_gt=f_gt,

                child_gt_95_ci=tuple(map(lambda x: tuple(map(dtype, x)), res["call_95_cis"])),
                mother_gt_95_ci=m_gt_95_ci,
                father_gt_95_ci=f_gt_95_ci,

                # child_gt_99_ci=parse_cis(locus_data[-1:].split("|")),
                # mother_gt_99_ci=m_gt_99_ci,
                # father_gt_99_ci=f_gt_99_ci,

                child_read_counts=StrKitJSONCalculator.get_read_counts(res, dtype=dtype),
                mother_read_counts=m_rcs,
                father_read_counts=f_rcs,

                reference_copies=dtype(res["ref_cn"]),

                decimal=fractional,

                test_to_perform=self.test_to_perform,
            ))

        return cr
