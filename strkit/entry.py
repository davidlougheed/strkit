from __future__ import annotations

import argparse
import pathlib
import os
import sys

from typing import Callable, Optional, Type

import strkit.constants as c
from strkit import __version__
from strkit.exceptions import ParamError, InputError
from strkit.logger import get_main_logger, attach_stream_handler, log_levels


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
        type=pathlib.Path,
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
        help="Maximum number of supporting reads before a locus is skipped.")

    # Min PHRED score of 13 for 95% confidence avg. in bases; although this may be over-estimated in ONT.
    # TODO: cite overest. and put a note in README.
    call_parser.add_argument(
        "--min-avg-phred",
        type=int,
        default=13,
        help="Minimum average PHRED score for relevant (flanking + TR) bases required for a read to be used.")

    call_parser.add_argument(
        "--flank-size",
        type=int,
        default=70,
        help="Number of bases around the locus to use for context.")

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
        "--sex-chr", "-x",
        type=str,
        help="Sex chromosome configuration to use for this sample (XX, XY, etc.) If left out, sex chromosomes will not "
             "be genotyped.")

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
        "--trf-bed",
        type=str,
        help="Specifies a TRF-derived BED file with all called loci and motifs (required for Straglr, since it doesn't "
             "respect the original motifs given to it.)")

    mi_parser.add_argument(
        "--exclude-loci-bed", "-e",
        type=str,
        help="Specifies a BED file with a list of loci to exclude calls from. Note that this does not handle regions "
             "right now (i.e., range calculations), just specific loci coordinates.")

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

    mi_parser.add_argument(
        "--test",
        nargs="?",
        type=str,
        default="none",
        const="x2",
        choices=("none", "x2", "wmw", "wmw-gt"),
        help="Perform statistical chi-squared (by default) or Wilcoxon-Mann-Whitney tests to detect de novo mutation "
             "in trios. Available for strkit-json only.")

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
    al_parser.add_argument("trf-bed", type=str, help="TRF BED file to convert from.")
    al_parser.add_argument("--caller", type=str, choices=(
        c.CALLER_EXPANSIONHUNTER,
        c.CALLER_GANGSTR,
        c.CALLER_STRAGLR,
        c.CALLER_TANDEM_GENOTYPES,
    ), help="Caller format to convert to.")


def add_vs_parser_args(vs_parser):
    vs_parser.add_argument("align_file", type=str, help="Alignment file to visualize.")
    vs_parser.add_argument(
        "--align-indices",
        nargs="*",
        type=str,
        help="Index/indices for alignment file to visualize. If specified, the number of arguments must match the "
             "number of given alignment files.")
    vs_parser.add_argument(
        "--ref", type=str, default="hg38",
        help="Reference genome code used for visualization and calls. Default: hg38")
    vs_parser.add_argument("--json", type=str, help="JSON file with STRkit calls to read from.")
    vs_parser.add_argument(
        "-i", type=int, default=1,
        help="1-based index of the locus to visualize in the JSON file. Default: 0")


def _exec_call(p_args) -> None:
    from strkit.call import call_sample, CallParams
    logger = get_main_logger(log_levels[p_args.log_level])
    call_sample(
        CallParams.from_args(logger, p_args),
        json_path=p_args.json,
        indent_json=p_args.indent_json,
        vcf_path=p_args.vcf,
        output_tsv=not p_args.no_tsv,
    )


def _exec_mi(p_args) -> None:
    from strkit.mi.base import BaseCalculator
    from strkit.mi.expansionhunter import ExpansionHunterCalculator
    from strkit.mi.gangstr import GangSTRCalculator
    from strkit.mi.repeathmm import RepeatHMMCalculator
    from strkit.mi.straglr import StraglrCalculator
    from strkit.mi.strkit import StrKitCalculator, StrKitJSONCalculator
    from strkit.mi.tandem_genotypes import TandemGenotypesCalculator

    calc_classes: dict[str, Type[BaseCalculator]] = {
        c.CALLER_EXPANSIONHUNTER: ExpansionHunterCalculator,
        c.CALLER_GANGSTR: GangSTRCalculator,
        c.CALLER_REPEATHMM: RepeatHMMCalculator,
        c.CALLER_STRAGLR: StraglrCalculator,
        c.CALLER_STRKIT: StrKitCalculator,
        c.CALLER_STRKIT_JSON: StrKitJSONCalculator,
        c.CALLER_TANDEM_GENOTYPES: TandemGenotypesCalculator,
    }

    caller = p_args.caller.lower()

    trf_bed_file = getattr(p_args, "trf_bed") or None
    if trf_bed_file is None and caller == c.CALLER_STRAGLR:
        raise ParamError("Using `strkit mi` with Straglr requires that the --trf-bed flag is used.")

    exclude_bed_file = p_args.exclude_loci_bed or None

    calc_class: Optional[Type[BaseCalculator]] = calc_classes.get(caller)
    if not calc_class:
        raise ParamError(f"Unknown or unimplemented caller '{caller}'")

    test_to_perform = p_args.test
    if caller != c.CALLER_STRKIT_JSON and test_to_perform != "none":
        raise ParamError(f"Caller '{caller}' does not support inheritance tests.")

    child_file = getattr(p_args, "child-calls")  # First call file is not optional
    mother_file = getattr(p_args, "mother-calls", None)
    father_file = getattr(p_args, "father-calls", None)

    mother_id = getattr(p_args, "mother_id", None)
    father_id = getattr(p_args, "father_id", None)
    child_id = getattr(p_args, "child_id", None)

    # TODO: Check that caller supports "combined" call files

    if (father_file is None or mother_file is None) and (mother_id is None or father_id is None or child_id is None):
        raise ParamError("If less than 3 genotype files are specified, sample IDs must be given")

    calc_inst = calc_class(
        child_call_file=child_file,
        mother_call_file=mother_file,
        father_call_file=father_file,

        child_id=child_id,
        mother_id=mother_id,
        father_id=father_id,

        loci_file=trf_bed_file,
        exclude_file=exclude_bed_file,

        widen=getattr(p_args, "widen", 0) or 0,

        test_to_perform=test_to_perform,
        sig_level=p_args.sig_level,
        mt_corr=p_args.mt_corr,

        debug=p_args.debug,
    )

    # TODO: Figure out contigs properly... sex chromosome stuff
    contigs = calc_inst.get_trio_contigs(include_sex_chromosomes=False)

    contig = getattr(p_args, "contig", None)

    if contig is not None and contig not in contigs:
        raise InputError(f"Could not find specified contig {p_args.contig} in trio contigs {contigs}")

    if contig is not None:
        contigs = {contig}

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
    return convert(getattr(p_args, "trf-file"), p_args.caller)


def _exec_viz_server(p_args):
    from strkit.json import json
    from strkit.viz.server import run_server as viz_run_server

    align_file = str(pathlib.Path(p_args.align_file).resolve())
    align_index = str(pathlib.Path(p_args.align_index).resolve()) if p_args.align_index else None

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

    if not (json_file := pathlib.Path(p_args.json)).exists():
        raise ParamError(f"Could not find JSON call report file at '{json_file}'")

    with open(json_file, "r") as jf:
        call_report = json.loads(jf.read())

    idx = p_args.i
    if idx < 1 or idx > len(call_report["results"]):
        raise InputError(f"JSON offset out of bounds: '{idx}'")

    viz_run_server(
        call_report=call_report,
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


def main(args: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="A toolkit for analyzing variation in short(ish) tandem repeats.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--version", "-v", action="version", version=__version__)

    subparsers = parser.add_subparsers()

    def _make_subparser(arg: str, help_text: str, exec_func: Callable, arg_func: Callable):
        sp = subparsers.add_parser(arg, help=help_text)
        sp.add_argument("--log-level", type=str, default="info", choices=("error", "warning", "info", "debug"))
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
        help_text="Convert TRF BED file to other caller formats.",
        exec_func=_exec_convert,
        arg_func=add_cv_parser_args)

    _make_subparser(
        "visualize",
        help_text="Start a web server to visualize results from an STR genotyping report.",
        exec_func=_exec_viz_server,
        arg_func=add_vs_parser_args)

    args = args or sys.argv[1:]
    p_args = parser.parse_args(args)

    ll = log_levels[p_args.log_level]


    if hasattr(p_args, "log_level"):
        logger = get_main_logger(ll)
        attach_stream_handler(ll, logger)
    else:
        logger = get_main_logger()

    if not getattr(p_args, "func", None):
        p_args = parser.parse_args(("--help",))

    try:
        logger.info(f"STRkit version {__version__}")
        p_args.func(p_args)
        return 0
    except ParamError as e:
        logger.critical(f"Paramter error: {e}")
        return 1
    except InputError as e:
        logger.critical(f"Input error: {e}")
        return 1


if __name__ == "__main__":
    exit(main())
