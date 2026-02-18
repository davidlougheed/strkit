from __future__ import annotations

import argparse
import os
import sys

from pathlib import Path
from typing import Callable, Type, TYPE_CHECKING

import strkit.constants as c
from strkit import __version__
from strkit.exceptions import ParamError, InputError
from strkit.logger import get_main_logger, attach_stream_handler, log_levels
from strkit.ploidy import PLOIDY_OPTIONS_HELP_TEXT

if TYPE_CHECKING:
    from logging import Logger


def add_call_parser_args(call_parser):
    call_parser.add_argument(
        "read_file",
        type=str,
        help="Indexed BAM/CRAM file(s) with reads to call from.")

    call_parser.add_argument(
        "--ref", "-r",
        type=str,
        required=True,
        help="Path to a reference genome, FASTA-formatted and indexed.")

    call_parser.add_argument(
        "--loci", "-l",
        type=str,
        required=True,
        help="Specifies a BED file with all loci to call. The last column of the BED must contain the motif sequences "
             "of the loci.")

    call_parser.add_argument(
        "--sample-id", "-s",
        type=str,
        help="Set a sample ID, or override the alignment file sample ID.")

    call_parser.add_argument(
        "--realign", "-a",
        action="store_true",
        help="Whether to perform local re-alignment to attempt recovery of soft-clipped reads. Some aligners may "
             "soft-clip around large insertions, e.g. with an expansion. Recommended for CCS ONLY!")

    call_parser.add_argument(
        "--hq",
        action="store_true",
        help="Whether to treat provided reads as 'fairly' accurate, i.e. not liable to be extremely far away from the "
             "DNA truth. Recommended for CCS, and CCS ONLY!")

    call_parser.add_argument(
        "--use-hp",
        action="store_true",
        help="Whether to use HP tags from the alignment file (i.e., a haplotagged alignment file), if available.")

    call_parser.add_argument(
        "--incorporate-snvs", "--snv", "-v",
        type=Path,
        help="A path to a dbSNP VCF file with a list of validated SNVs to help phase with. Specifying this enables the "
             "use of read-phased SNVs to help properly call genotypes. This option can only be used if --hq "
             "is enabled. Use with CCS ONLY! This option is currently EXPERIMENTAL!")

    call_parser.add_argument(
        "--snv-min-base-qual", "--snv-mbq",
        type=int,
        default=20,
        help="Minimum PHRED quality score for bases of SNVs to use for phasing.")

    call_parser.add_argument(
        "--targeted", "-t",
        action="store_true",
        help="Whether these reads come from targeted (e.g. gene-of-interest or WES) or genome-wide sequencing (WGS).")

    call_parser.add_argument(
        "--respect-ref", "-e",
        action="store_true",
        help="Do not extend reference TR region coordinates out from what is specified in the catalog. This should "
             "bring results closer to other callers, but may reduce accuracy.")

    call_parser.add_argument(
        "--count-kmers", "-k",
        nargs="?",
        type=str,
        default="none",
        const="peak",
        choices=("none", "peak", "read", "both"),
        help="Count motif-sized k-mers in each read's tandem repeat region. Useful for detecting heterogeneous "
             "expansion patterns. If set to 'peak', read-level counts will be aggregated into peak-level counts for "
             "more efficient space and processing. This is the default if the flag is passed by itself.")

    call_parser.add_argument(
        "--consensus", "-c",
        action="store_true",
        help="Calculate consensus sequences for alleles. This flag increases runtime.")

    call_parser.add_argument(
        "--large-consensus-length",
        type=int,
        default=2000,
        help="If the smallest tandem repeat sequence is longer than this, the reads will be downsampled to "
             "--max-n-large-consensus-reads.")

    call_parser.add_argument(
        "--max-n-large-consensus-reads",
        type=int,
        default=20,
        help="When the smallest tandem repeat sequence is longer than --large-consensus-length, the number of reads "
             "will be downsampled to this value.")

    call_parser.add_argument(
        "--max-mdn-poa-length",
        type=int,
        default=5000,
        help="Maximum median STR sequence length before we use a best-representative-read method for the 'consensus' "
             "STR sequence rather than partial order alignment (POA) due to computational performance.")

    call_parser.add_argument(
        "--min-reads",
        type=int,
        default=4,
        help="Minimum number of supporting reads needed to call a locus.")

    call_parser.add_argument(
        "--min-allele-reads",
        type=int,
        default=2,
        help="Minimum number of supporting reads needed to call a specific allele peak.")

    call_parser.add_argument(
        "--max-reads",
        type=int,
        default=250,
        help="Maximum number of supporting reads to use for a locus. Additional reads will be ignored.")

    # Min PHRED score of 13 for 95% confidence avg. in bases; although this may be over-estimated in ONT.
    # TODO: cite overest. and put a note in README.
    call_parser.add_argument(
        "--min-avg-phred",
        type=int,
        default=13,
        help="Minimum average PHRED score for relevant (flanking + TR) bases required for a read to be used.")

    call_parser.add_argument(
        "--significant-clip-threshold",
        type=int,
        default=100,
        help=(
            "Minimum number of bases to consider the end of a read 'significantly' soft-clipped, meaning we should do "
            "some read trimming in the SNV-finding process to avoid SNV false positives."
        ))

    # Minimum read alignment score (fractional) to include a read in calling.
    call_parser.add_argument(
        "--min-read-align-score", "--mras",
        type=float,
        default=0.9,
        help=(
            "Minimum read alignment score (fractional) to include a read in calling. A good value for pure tandem "
            "repeats is 0.9. A good value for more lenient genotyping is ~0.4."
        ))

    # Maximum read copy number counting iterations before terminating copy number counting.
    call_parser.add_argument(
        "--max-rcn-iters", "--mrci",
        type=int,
        default=50,
        help=(
            "Maximum number of read copy-number counting iterations to perform. Loci which require a lot of "
            "iterations are probably impure tandem repeats, for which the resulting copy number will not be very "
            "accurate anyway."
        ))

    call_parser.add_argument(
        "--flank-size",
        type=int,
        default=70,
        help="Number of bases around the locus to use for context.")

    call_parser.add_argument(
        "--skip-supplementary", "--skip-supp",
        action="store_true",
        help="Whether to skip supplementary aligned reads.")

    call_parser.add_argument(
        "--skip-secondary", "--skip-sec",
        action="store_true",
        help="Whether to skip secondary aligned reads.")

    call_parser.add_argument(
        "--processes", "-p",
        type=int,
        default=1,
        help="Number of processes to use when calling.")

    call_parser.add_argument(
        "--num-bootstrap", "-b",
        type=int,
        default=100,
        help="Specifies the number of bootstrap samples to use to calculate confidence intervals and median genotype "
             "estimate.")

    call_parser.add_argument(
        "--gm-filter-factor",
        type=int,
        default=3,
        help=(
            "Filter factor for copy number Gaussian mixture models. With two alleles, a value of ex. 3 would mean "
            "peaks of weight 1/(3*[n_alleles=2]) would be filtered out."
        ))

    call_parser.add_argument(
        "--ploidy", "--sex-chr", "-x",
        type=str,
        default="diploid_autosomes",
        help=f"Sex chromosome configuration to use for this sample (XX, XY, etc.) If left out, sex chromosomes will "
             f"not be genotyped. Bundled options are: {PLOIDY_OPTIONS_HELP_TEXT}. See documentation for more "
             f"information on how these JSON configuration files can be formatted.")

    # BEGIN FILE OUTPUT ARGUMENTS ======================================================================================

    # Begin JSON output arguments --------------------------------------------------------------------------------------
    call_parser.add_argument(
        "--json", "-j",
        type=str,
        help="Path to write JSON-formatted calls to. If left blank, no JSON file will be written. If the value is set "
             "to 'stdout', JSON will be written to stdout, after the TSV unless TSV output is disabled.")

    call_parser.add_argument(
        "--indent-json", "-i",
        action="store_true",
        help="If passed alongside --json [x], the JSON output will be indented to be more human readable but "
             "less compact.")
    # End JSON output arguments ----------------------------------------------------------------------------------------

    # Begin VCF output arguments ---------------------------------------------------------------------------------------
    call_parser.add_argument(
        "--vcf",
        type=str,
        help="Path to write VCF-formatted calls to. STR alleles are currently formatted symbolically rather than as "
             "consensus sequences.")
    call_parser.add_argument(
        "--vcf-anchor-size",
        type=int,
        default=5,
        help="Size of upstream sequence (including potential small upstream variants) to include in VCF output "
             "consensus sequences.")
    # End VCF output arguments -----------------------------------------------------------------------------------------

    call_parser.add_argument(
        "--no-tsv",
        action="store_true",
        help="If passed, no TSV call output will be written to stdout.")

    # END FILE OUTPUT ARGUMENTS ========================================================================================

    call_parser.add_argument(
        "--seed",
        type=int,
        help="Random seed to pass to random number generators for result replicability.")

    call_parser.add_argument(
        "--log-progress-interval",
        type=int,
        default=120,
        help="How often (in seconds) to log calling progress INFO messages.")

    call_parser.add_argument("--profile", action="store_true", help="Profile function calls (for development.)")


def add_mi_parser_args(mi_parser):
    mi_parser.add_argument(
        "--debug",
        action="store_true",
        help="Whether to enable debug mode (with additional logging.)")

    mi_parser.add_argument(
        "--widen",
        type=float,
        help="Widen the CIs when calculating MI by a given proportion to account for small errors.")

    mi_parser.add_argument(
        "--caller",
        type=str,
        choices=c.MI_CALLERS,
        required=True,
        help="The program which called the TR genotypes.")

    mi_parser.add_argument("--contig", type=str, help="Specifies a specific contig to process (optional).")

    mi_parser.add_argument(
        "--motif-bed", "--trf-bed",
        type=str,
        help="Specifies a BED file with all called loci and motifs (required for LongTR and Straglr, since they don't "
             "contain information about the original locus motif.)")

    mi_parser.add_argument(
        "--exclude-loci-bed", "-e",
        type=str,
        help="Specifies a BED file with a list of loci to exclude calls from. Note that this does not handle regions "
             "right now (i.e., range calculations), just specific loci coordinates.")

    mi_parser.add_argument(
        "--one-based-bed",
        action="store_true",
        help="If passed, the BED(s) will be treated as non-standard, with one-based coordinates.")

    mi_parser.add_argument(
        "--json", "-j",
        type=str,
        help="Path to write a JSON-formatted Mendelian inheritance/error report to. If left blank, no JSON file will "
             "be written. If the value is set to 'stdout', JSON will be written to stdout, after the TSV unless TSV "
             "output is disabled.")

    mi_parser.add_argument(
        "--no-tsv",
        action="store_true",
        help="If passed, no TSV call output will be written to stdout.")

    mi_parser.add_argument(
        "--mismatch-out-mi",
        type=str,
        choices=("strict", "pm1", "ci_95", "ci_99", "seq", "sl", "sl_pm1"),
        default="pm1",
        help="This parameter sets which type of Mendelian inheritance will be used to determine which "
             "non-MI-respecting loci are output."
    )

    # Histogram-related arguments -------------------------------------------------------
    mi_parser.add_argument(
        "--hist",
        action="store_true",
        help="Generate a textual histogram of binned Mendelian inheritance (MI) rates.")
    mi_parser.add_argument(
        "--bin-width",
        type=int,
        default=10,
        help="Bin width for generated MI histograms.")
    # -----------------------------------------------------------------------------------

    # Test-related arguments ------------------------------------------------------------

    mi_parser.add_argument(
        "--test",
        nargs="?",
        type=str,
        default="none",
        const="x2",
        choices=("none", "x2", "wmw", "wmw-gt"),
        help="Perform statistical chi-squared (by default) or Wilcoxon-Mann-Whitney tests to detect de novo mutation "
             "in trios. Available for strkit-json and strkit-vcf only.")

    mi_parser.add_argument(
        "--sig-level",
        type=float,
        default=0.05,
        help=(
            "Significance level for de novo mutation test results. Without correction, this parameter functions as a "
            "simple threshold. When a multiple testing correction method is applied, this parameter specifies the "
            "desired significance level, alpha."
        ))

    mi_parser.add_argument(
        "--only-phased",
        action="store_true",
        help="Whether to only compare phasing-supported (SNVs, haplotags) STR loci. Available for strkit-vcf only.",
    )

    mi_parser.add_argument(
        "--mt-corr",
        type=str,
        choices=(
            "none",
            "bonferroni",
            "sidak",
            "holm-sidak",
            "holm",
            "simes-hochberg",
            "fdr_bh",
            "fdr_by",
            "fdr_tsbh",
            "fdr_tsbky",
        ),
        default="none",
        help=(
            "Method to use for multiple testing correction. Under the hood, this uses the statsmodels library, and"
            "so the options are (mostly) identical to theirs: https://www.statsmodels.org/v0.13.0/generated/statsmodels"
            ".stats.multitest.multipletests.html"
        ),
    )

    # -----------------------------------------------------------------------------------

    mi_parser.add_argument(
        "--child-id",
        type=str,
        help="Specifies the child's sample ID (if a combined VCF is used for the callset.)")
    mi_parser.add_argument(
        "--mother-id",
        type=str,
        help="Specifies the mother's sample ID (if a combined VCF is used for the callset.)")
    mi_parser.add_argument(
        "--father-id",
        type=str,
        help="Specifies the father's sample ID (if a combined VCF is used for the callset.)")

    mi_parser.add_argument(
        "child-calls",
        type=str,
        help="File (VCF, TSV, etc.) containing the genotype calls for the child (or all three individuals.)")
    mi_parser.add_argument(
        "mother-calls",
        nargs="?",
        type=str,
        help="File (VCF, TSV, etc.) containing the genotype calls for the mother. "
             "If unspecified, defaults to the first genotype file given.")
    mi_parser.add_argument(
        "father-calls",
        nargs="?",
        type=str,
        help="File (VCF, TSV, etc.) containing the genotype calls for the father. "
             "If unspecified, defaults to the first genotype file given.")


def add_cc_parser_args(cc_parser):
    cc_parser.add_argument(
        "--caller",
        type=str,
        choices=(c.CALLER_STRAGLR,),
        required=True,
        help="The program which called the TR genotypes.")

    cc_parser.add_argument("paths", type=str, nargs="+", help="Paths to the BED catalogs to combine.")


def add_cv_parser_args(al_parser):
    from strkit.convert.converter import CONVERTER_OUTPUT_FORMATS
    al_parser.add_argument("in_file", type=str, help="Input file to convert from.")
    al_parser.add_argument(
        "--in-format", type=str, choices=("trf", "trgt"), default="trf", help="Format to convert from."
    )
    al_parser.add_argument("--out-format", type=str, choices=CONVERTER_OUTPUT_FORMATS, help="Format to convert to.")


def add_vs_parser_args(vs_parser):
    vs_parser.add_argument("align_file", type=str, help="Alignment file to visualize.")
    vs_parser.add_argument(
        "--align-index",
        nargs="*",
        type=str,
        help="Index for alignment file to visualize.")
    vs_parser.add_argument(
        "--ref", type=str, default="hg38",
        help="Reference genome code used for visualization and calls. Default: hg38")
    vs_parser.add_argument("--json", type=str, help="JSON file with STRkit calls to read from.")
    vs_parser.add_argument(
        "-i", "--index", type=int, default=1,
        help="1-based index of the locus to visualize in the JSON file. Default: 0")
    vs_parser.add_argument("--port", type=int, default=5011, help="Port to start the visualization server on.")


def _main_logger(p_args) -> Logger:
    return get_main_logger(log_levels[p_args.log_level])


def _exec_call(p_args) -> None:
    from strkit.call.call_sample import call_sample
    from strkit.call.params import CallParams
    call_sample(
        CallParams.from_args(_main_logger(p_args), p_args),
        json_path=p_args.json,
        indent_json=p_args.indent_json,
        vcf_path=p_args.vcf,
        output_tsv=not p_args.no_tsv,
    )


def _exec_mi(p_args) -> None:
    from strkit.mi.base import BaseCalculator
    from strkit.mi.expansionhunter import ExpansionHunterCalculator
    from strkit.mi.gangstr import GangSTRCalculator
    from strkit.mi.generic_vcf import GenericVCFLengthCalculator
    from strkit.mi.repeathmm import RepeatHMMCalculator
    from strkit.mi.straglr import StraglrCalculator
    from strkit.mi.strkit import StrKitCalculator, StrKitJSONCalculator, StrKitVCFCalculator
    from strkit.mi.tandem_genotypes import TandemGenotypesCalculator
    from strkit.mi.trgt import TRGTCalculator

    calc_classes: dict[str, Type[BaseCalculator]] = {
        c.CALLER_EXPANSIONHUNTER: ExpansionHunterCalculator,
        c.CALLER_GANGSTR: GangSTRCalculator,
        c.CALLER_GENERIC_VCF: GenericVCFLengthCalculator,
        c.CALLER_LONGTR: GenericVCFLengthCalculator,  # Alias
        c.CALLER_REPEATHMM: RepeatHMMCalculator,
        c.CALLER_STRDUST: GenericVCFLengthCalculator,  # Alias
        c.CALLER_STRAGLR: StraglrCalculator,
        c.CALLER_STRKIT: StrKitCalculator,
        c.CALLER_STRKIT_JSON: StrKitJSONCalculator,
        c.CALLER_STRKIT_VCF: StrKitVCFCalculator,
        c.CALLER_TANDEM_GENOTYPES: TandemGenotypesCalculator,
        c.CALLER_TRGT: TRGTCalculator,
    }

    caller = p_args.caller.lower()

    motif_bed_file = getattr(p_args, "motif_bed") or None
    if motif_bed_file is None and caller in (c.CALLER_LONGTR, c.CALLER_STRAGLR):
        raise ParamError(f"Using `strkit mi` with the '{caller}' caller requires that the --motif-bed flag is used.")

    exclude_bed_file = p_args.exclude_loci_bed or None

    calc_class: Type[BaseCalculator] | None = calc_classes.get(caller)
    if not calc_class:
        raise ParamError(f"Unknown or unimplemented caller '{caller}'")

    mismatch_out_mi = p_args.mismatch_out_mi

    test_to_perform = p_args.test
    if caller not in (c.CALLER_STRKIT_JSON, c.CALLER_STRKIT_VCF) and test_to_perform != "none":
        raise ParamError(f"Caller '{caller}' does not support inheritance tests.")

    only_phased = p_args.only_phased
    if caller != c.CALLER_STRKIT_VCF and only_phased:
        raise ParamError(f"Caller '{caller}' does not support --only-phased.")

    child_file = getattr(p_args, "child-calls")  # First call file is not optional
    mother_file = getattr(p_args, "mother-calls", None)
    father_file = getattr(p_args, "father-calls", None)

    mother_id = getattr(p_args, "mother_id", None)
    father_id = getattr(p_args, "father_id", None)
    child_id = getattr(p_args, "child_id", None)

    # TODO: Check that caller supports "combined" call files

    if (father_file is None or mother_file is None) and (mother_id is None or father_id is None or child_id is None):
        raise ParamError("If less than 3 genotype files are specified, sample IDs must be given")

    logger = _main_logger(p_args)

    calc_inst = calc_class(
        child_call_file=child_file,
        mother_call_file=mother_file,
        father_call_file=father_file,

        child_id=child_id,
        mother_id=mother_id,
        father_id=father_id,

        loci_file=motif_bed_file,
        exclude_file=exclude_bed_file,
        one_based_loci=p_args.one_based_bed,

        widen=getattr(p_args, "widen", 0) or 0,

        mismatch_out_mi=mismatch_out_mi,
        test_to_perform=test_to_perform,
        sig_level=p_args.sig_level,
        mt_corr=p_args.mt_corr,
        only_phased=only_phased,

        debug=p_args.debug,
        logger=logger,
    )

    # TODO: Figure out contigs properly... sex chromosome stuff
    contigs = calc_inst.get_trio_contigs(include_sex_chromosomes=False)

    contig = getattr(p_args, "contig", None)

    if contig is not None and contig not in contigs:
        raise InputError(f"Could not find specified contig {p_args.contig} in trio contigs {contigs}")

    if contig is not None:
        contigs = {contig}

    logger.info("Executing MI calculator: %s", caller)
    logger.info("MI included contigs: %s", contigs)
    res = calc_inst.calculate(included_contigs=contigs)

    if not res:
        return

    if not p_args.no_tsv:
        print(str(res))
        if p_args.hist:
            print(res.histogram_text(bin_width=p_args.bin_width))
        print("---")
        sys.stdout.write(res.locus_tsv())
        sys.stdout.flush()

    if p_args.json:
        res.write_report_json(json_path=p_args.json, bin_width=p_args.bin_width)


def _exec_combine_catalogs(p_args):
    from strkit.catalog.combine import combine_catalogs
    return combine_catalogs(p_args.caller, p_args.paths)


def _exec_convert(p_args):
    from strkit.convert.converter import convert
    return convert(p_args.in_file, p_args.in_format, p_args.out_format, _main_logger(p_args))


def _exec_viz_server(p_args):
    from strkit.json import json
    from strkit.viz.server import run_server as viz_run_server

    align_file = str(Path(p_args.align_file).resolve())
    align_index = str(Path(p_args.align_index).resolve()) if p_args.align_index else None

    align_format = os.path.splitext(align_file)[-1].lstrip(".")
    if align_format not in ("bam", "cram"):
        raise ParamError(f"File type '{align_format}' not supported")

    if not align_index:
        align_index = f"{align_file}.{'crai' if align_format == 'cram' else 'bai'}"
        if not os.path.exists(align_index):
            raise ParamError(f"Missing index at '{align_index}'")

    # TODO: Conditionally use this code if ref looks like a path
    # ref = pathlib.Path(p_args.ref).resolve()
    #
    # if not ref.exists():
    #     print(f"Error: could not find reference genome at '{ref}'", file=sys.stderr)
    #     return 1
    #
    # ref_index = str(pathlib.Path(p_args.ref_index).resolve()) if p_args.ref_index else str(ref) + ".fai"
    # if not os.path.exists(ref_index):
    #     print(f"Error: missing .fai for reference genome at '{ref}'", file=sys.stderr)
    #     return 1

    if not (json_file := Path(p_args.json)).exists():
        raise ParamError(f"Could not find JSON call report file at '{json_file}'")

    with open(json_file, "r") as jf:
        call_report = json.loads(jf.read())

    idx = p_args.index
    if idx < 1 or idx > len(call_report["results"]):
        raise InputError(f"JSON offset out of bounds: '{idx}'")

    viz_run_server(
        call_report=call_report,
        port=p_args.port,
        # ---------------------------------------------
        initial_i=idx-1,
        ref=p_args.ref,
        # ref_index=ref_index,
        align_file=align_file,
        align_index=align_index,
        align_name=os.path.basename(align_file),
        align_index_name=os.path.basename(align_index),
        align_format=align_format,
    )

    return 0


def main(args: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="A toolkit for analyzing variation in short(ish) tandem repeats.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--version", "-v", action="version", version=__version__)

    subparsers = parser.add_subparsers()

    def _make_subparser(*names: str, help_text: str, exec_func: Callable, arg_func: Callable):
        sp = subparsers.add_parser(names[0], aliases=names[1:], help=help_text)
        sp.add_argument("--log-level", type=str, default="info", choices=("error", "warning", "info", "debug"))
        sp.add_argument(
            "--verbose",
            action="store_true",
            help="If --log-level is debug, this will yield many more read-level debug messages.",
        )
        sp.set_defaults(func=exec_func)
        arg_func(sp)

    _make_subparser(
        "call",
        help_text="A tandem repeat (TR) caller designed for high-fidelity long reads.",
        exec_func=_exec_call,
        arg_func=add_call_parser_args)

    _make_subparser(
        "mi",
        help_text="A Mendelian inheritance calculator for different TR genotyping callers.",
        exec_func=_exec_mi,
        arg_func=add_mi_parser_args)

    _make_subparser(
        "combine-catalogs",
        help_text="Combine Straglr result catalogs for use in calling with a consistent motif set.",
        exec_func=_exec_combine_catalogs,
        arg_func=add_cc_parser_args)

    _make_subparser(
        "convert",
        help_text="Convert a repeat catalog to other formats.",
        exec_func=_exec_convert,
        arg_func=add_cv_parser_args)

    _make_subparser(
        "visualize", "vis", "viz",
        help_text="Start a web server to visualize results from an STR genotyping report.",
        exec_func=_exec_viz_server,
        arg_func=add_vs_parser_args)

    args = args or sys.argv[1:]
    p_args = parser.parse_args(args)

    if hasattr(p_args, "log_level"):
        logger = _main_logger(p_args)
        attach_stream_handler(log_levels[p_args.log_level], logger)
    else:
        logger = get_main_logger()

    if not getattr(p_args, "func", None):
        p_args = parser.parse_args(("--help",))

    try:
        logger.info("STRkit version %s", __version__)
        p_args.func(p_args)
        return 0
    except ParamError as e:
        logger.critical(f"Parameter error: {e}")
        return 1
    except InputError as e:
        logger.critical(f"Input error: {e}")
        return 1


if __name__ == "__main__":
    exit(main())
