import pysam

__all__ = ["VCFCalculatorMixin"]


class VCFCalculatorMixin:
    @staticmethod
    def get_contigs_from_files(mother_call_file, father_call_file, child_call_file):
        mvf = pysam.VariantFile(str(mother_call_file))
        mc = set(mvf.header.contigs)

        fvf = pysam.VariantFile(str(father_call_file))
        fc = set(fvf.header.contigs)

        cvf = pysam.VariantFile(str(child_call_file))
        cc = set(cvf.header.contigs)

        return mc, fc, cc
