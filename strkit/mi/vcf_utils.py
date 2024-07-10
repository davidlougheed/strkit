from __future__ import annotations

import pysam

__all__ = ["VCFCalculatorMixin"]


class VCFCalculatorMixin:
    @staticmethod
    def get_contigs_from_files(mother_call_file, father_call_file, child_call_file) -> tuple[set, set, set]:
        with pysam.VariantFile(str(mother_call_file)) as mvf:
            mc = set(mvf.header.contigs)

        with pysam.VariantFile(str(father_call_file)) as fvf:
            fc = set(fvf.header.contigs)

        with pysam.VariantFile(str(child_call_file)) as cvf:
            cc = set(cvf.header.contigs)

        return mc, fc, cc
