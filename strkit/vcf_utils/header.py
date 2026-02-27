from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pysam import VariantHeader, VariantHeaderMetadata
from typing import TYPE_CHECKING
from .. import __version__

if TYPE_CHECKING:
    from pathlib import Path

__all__ = [
    "VCF_INFO_VT",
    "VCF_INFO_MOTIF",
    "VCF_INFO_REFMC",
    "VCF_INFO_BED_START",
    "VCF_INFO_BED_END",
    "VCF_INFO_ANCH",
    "VT_STR",
    "VT_SNV",
    "build_vcf_header",
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

# VCF_ALLELE_CNV_TR = "<CNV:TR>"

# VCF_TR_INFO_RECORDS: tuple[tuple[str, str, str, str], ...] = (
#     ("SVLEN", "A", "Integer", "Length of the structural variant"),
#     ("CN", "A", "Float", "Copy number of allele"),
#     ("RN", "A", "Integer", "Total number of repeat sequences in this allele"),
#     ("RUS", ".", "String", "Repeat unit sequence of the corresponding repeat sequence"),
#     ("RUL", ".", "Integer", "Repeat unit length of the corresponding repeat sequence"),
#     ("RB", ".", "Integer", "Total number of bases in the corresponding repeat sequence"),
#     ("CIRUC", ".", "Float", "Confidence interval around RUC"),
#     ("CIRB", ".", "Integer", "Confidence interval around RB"),
# )


def add_vcf_file_date(vh: VariantHeader):
    """
    Adds a fileDate header metadata line to a variant header via mutation.
    :param vh: The variant header to add to.
    """
    now = datetime.now()
    vh.add_meta("fileDate", f"{now.year}{now.month:02d}{now.day:02d}")


def build_vcf_header(
    command: str,
    sample_ids: tuple[str, ...],
    reference_file: Path | str,  # if Path, turn into file URI. Otherwise, treat "as-is".
    partial_phasing: bool,
    num_loci: int,
    loci_hash: str,
    contigs: tuple[tuple[str, int], ...] | None = None,
) -> VariantHeader:
    vh = VariantHeader()  # automatically sets VCF version to 4.2

    # Add file date
    add_vcf_file_date(vh)

    # Add source
    vh.add_meta("source", "strkit")

    # Mark that we have partial phasing if we're using HP/SNVs
    if partial_phasing:
        vh.add_meta("phasing", "partial")

    # Add an absolute path to the reference genome. If we're passed a string (from merge), assume this is already a file
    # URI.
    vh.add_meta(
        "reference",
        f"file://{str(reference_file.resolve().absolute())}" if isinstance(reference_file, Path) else reference_file
    )

    if contigs is not None:
        # If we're passed contig records directly (from merge), add them rather than pulling them from the reference.
        for contig, contig_length in contigs:
            vh.contigs.add(contig, length=contig_length)
    else:
        # Add all contigs from the reference genome file + lengths
        from pysam import FastaFile
        with FastaFile(str(reference_file)) as rf:
            for contig in rf.references:
                vh.contigs.add(contig, length=rf.get_reference_length(contig))

    # Add STRkit-specific fields:
    #  - marking version
    vh.add_meta("strkitVersion", str(__version__))
    #  - the subcommand being used to generate this VCF (call|merge)
    vh.add_meta("strkitCommand", command)
    #  - indicating number of loci provided (i.e., catalogue size)
    vh.add_meta("strkitCatalogNumLoci", str(num_loci))
    #  - indicating hash of STRkitLocus objects (for checking catalogue sameness)
    vh.add_meta("strkitCatalogLociHash", loci_hash)

    # Add CNV:TR alt type (symbolic allele: tandem repeat)
    # vh.add_meta("ALT", "<ID=CNV:TR,Description=\"Tandem repeat\">")

    # Set up basic VCF formats
    vh.formats.add("AD", ".", "Integer", "Read depth for each allele")
    vh.formats.add("ANCL", ".", "Integer", "Anchor length for the ref and each alt, five-prime of TR sequence")
    vh.formats.add("CONS", ".", "String", "Consensus methods used for each alt (single/poa/best_rep)")
    vh.formats.add("DP", 1, "Integer", "Read depth")
    vh.formats.add("DPS", 1, "Integer", "Read depth (supporting reads only)")
    vh.formats.add("GT", 1, "String", "Genotype")
    vh.formats.add("MC", ".", "Integer", "Motif copy number for each allele")
    vh.formats.add("MCCI", ".", "String", "Motif copy number 95% confidence interval for each allele")
    vh.formats.add("MCRL", ".", "String", "Read-level motif copy numbers for each allele")
    vh.formats.add("MMAS", 1, "Float", "Mean model (candidate TR sequence) alignment score across reads.")
    vh.formats.add("NSNV", 1, "Integer", "Number of supporting SNVs for the STR peak-call")
    vh.formats.add("PS", 1, "Integer", "Phase set")
    vh.formats.add("PM", 1, "String", "Peak-calling method (dist/snv+dist/snv/hp)")

    # Set up VCF info fields (mutates vh.info)
    VCF_INFO_VT.add_to_info(vh.info)
    VCF_INFO_MOTIF.add_to_info(vh.info)
    VCF_INFO_REFMC.add_to_info(vh.info)
    VCF_INFO_BED_START.add_to_info(vh.info)
    VCF_INFO_BED_END.add_to_info(vh.info)
    VCF_INFO_ANCH.add_to_info(vh.info)

    # Add INFO records for tandem repeat copies - these are new to VCF4.4!  TODO
    # for iv in VCF_TR_INFO_RECORDS:
    #     vh.info.add(*iv)

    # Add the sample(s)
    for sample_id in sample_ids:
        vh.add_sample(sample_id)

    return vh
