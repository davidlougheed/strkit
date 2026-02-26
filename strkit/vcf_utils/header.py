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


def add_vcf_file_date(vh: VariantHeader):
    """
    Adds a fileDate header metadata line to a variant header via mutation.
    :param vh: The variant header to add to.
    """
    now = datetime.now()
    vh.add_meta("fileDate", f"{now.year}{now.month:02d}{now.day:02d}")
