from dataclasses import dataclass
from datetime import datetime
from pysam import VariantHeader, VariantHeaderMetadata

__all__ = [
    "VCF_INFO_VT",
    "VCF_INFO_MOTIF",
    "VCF_INFO_REFMC",
    "VCF_INFO_BED_START",
    "VCF_INFO_BED_END",
    "VCF_INFO_ANCH",
    "VT_STR",
    "VT_SNV",
    "SkipWritingLocus",
    "blank_vcf_genotype",
    "genotype_indices",
]


@dataclass
class VcfInfo:
    key: str
    n: int | str  # number of values in record row
    type: str  # VCF type
    description: str

    def add_to_info(self, info: VariantHeaderMetadata):
        """
        Add info field to variant file header metadata via mutation.
        :param info: VariantHeaderMetadata instance to add the record to. This object is mutated.
        """
        info.add(self.key, self.n, self.type, self.description)


VCF_INFO_VT = VcfInfo("VT", 1, "String", "Variant record type (str/snv)")
VCF_INFO_MOTIF = VcfInfo("MOTIF", 1, "String", "Motif string")
VCF_INFO_REFMC = VcfInfo("REFMC", 1, "Integer", "Motif copy number in the reference genome")
VCF_INFO_BED_START = VcfInfo(
    "BED_START",
    1,
    "Integer",
    "Original start position of the locus as defined in the catalog (0-based inclusive)",
)
VCF_INFO_BED_END = VcfInfo(
    "BED_END",
    1,
    "Integer",
    "Original end position of the locus as defined in the catalog (0-based exclusive, i.e., 1-based)",
)
VCF_INFO_ANCH = VcfInfo("ANCH", 1, "Integer", "Five-prime anchor size")

VT_STR = "str"
VT_SNV = "snv"


class SkipWritingLocus(Exception):
    pass


def add_vcf_file_date(vh: VariantHeader):
    """
    Adds a fileDate header metadata line to a variant header via mutation.
    :param vh: The variant header to add to.
    """
    now = datetime.now()
    vh.add_meta("fileDate", f"{now.year}{now.month:02d}{now.day:02d}")


def blank_vcf_genotype(n_alleles: int) -> tuple[None, ...]:
    """
    Generates a blank VCF genotype for a sample locus.
    :param n_alleles: Number of alleles to use for the blank genotype.
    :return: Tuple of [None] * <number of alleles>
    """
    return tuple([None] * n_alleles)


def genotype_indices(alleles: tuple[str, ...], call: tuple[str, ...] | None, n_alleles: int) -> tuple[int | None, ...]:
    """
    Given possible variant alleles for a particular locus and a sample's call (the sample's alleles), plus the number of
    possible alleles for the sample at the locus, return a VCF-compatible genotype index tuple.
    :param alleles: Possible variant alleles: (ref, alt1, alt2, alt3, ...)
    :param call: Call for the sample: (alt3, alt1) or None if not called.
    :param n_alleles: Number of alleles for the sample locus (used for computing blank genotypes if not called.)
    :return: Genotype index tuple, e.g., (3, 1).
    """
    return tuple(map(alleles.index, call)) if call is not None else blank_vcf_genotype(n_alleles)
