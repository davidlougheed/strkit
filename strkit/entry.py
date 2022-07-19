import argparse
import json
import pathlib
import os
import sys

from typing import Dict, List, Optional, Type

import strkit.constants as c
from strkit.aligned.lengths import aligned_lengths_cmd
from strkit.call import call_sample, re_call_all_alleles
from strkit.catalog.combine import combine_catalogs
from strkit.convert.converter import convert
from strkit.mi.base import BaseCalculator
from strkit.mi.expansionhunter import ExpansionHunterCalculator
from strkit.mi.gangstr import GangSTRCalculator
from strkit.mi.repeathmm import RepeatHMMCalculator, RepeatHMMReCallCalculator
from strkit.mi.straglr import StraglrCalculator, StraglrReCallCalculator
from strkit.mi.strkit import StrKitCalculator
from strkit.mi.tandem_genotypes import TandemGenotypesCalculator, TandemGenotypesReCallCalculator
from strkit.viz.server import run_server as viz_run_server


def add_call_parser_args(call_parser):
    call_parser.add_argument(
        "read_files",
        nargs="+",
        type=str,
        help="BAM file(s) with reads to call from. If multiple files are specified, the reads will be pooled.")

    call_parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="Path to a reference genome, FASTA-formatted and indexed.")

    call_parser.add_argument(
        "--loci",
        type=str,
        required=True,
        help="Specifies a BED file with all loci to call. The last column of the BED must contain the motif sequences "
             "of the loci.")

    call_parser.add_argument(
        "--targeted",
        action="store_true",
        help="Whether these reads come from targeted (e.g. gene-of-interest or WES) or genome-wide sequencing (WGS).")

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
        "--flank-size",
        type=int,
        default=70,
        help="Number of bases around the locus to use for context.")

    call_parser.add_argument(
        "--processes",
        type=int,
        default=1,
        help="Number of processes to use when calling.")

    call_parser.add_argument(
        "--num-bootstrap",
        type=int,
        default=100,
        help="Specifies the number of bootstrap samples to use to calculate confidence intervals and median genotype "
             "estimate.")

    call_parser.add_argument(
        "--sex-chr",
        type=str,
        help="Sex chromosome configuration to use for this sample (XX, XY, etc.) If left out, sex chromosomes will not "
             "be genotyped.")

    call_parser.add_argument(
        "--json",
        type=str,
        help="Path to write JSON-formatted calls to. If left blank, no JSON file will be written. If the value is set "
             "to 'stdout', JSON will be written to stdout, after the TSV unless TSV output is disabled.")

    call_parser.add_argument(
        "--no-tsv",
        action="store_true",
        help="If passed, no TSV call output will be written to stdout.")


def add_re_call_parser_args(re_call_parser):
    re_call_parser.add_argument(
        "--caller",
        type=str,
        choices=c.CALL_SUPPORTED_CALLERS,
        required=True,
        default=argparse.SUPPRESS,
        help="The program which called the TR genotypes.")

    re_call_parser.add_argument(
        "--num-bootstrap",
        type=int,
        default=100,
        help="Specifies the number of bootstrap samples to use to calculate confidence intervals and median genotype "
             "estimate.")

    re_call_parser.add_argument(
        "--min-reads",
        type=int,
        default=4,
        help="Minimum number of supporting reads needed to call a locus.")

    re_call_parser.add_argument(
        "--min-allele-reads",
        type=int,
        default=2,
        help="Minimum number of supporting reads needed to call a specific allele peak.")

    re_call_parser.add_argument(
        "--read-bias-corr-min",
        type=int,
        default=4,
        help="Minimum strand coverage (required on both strands) needed to attempt strand bias correction via "
             "resampling. If this is set too low, may lead to calling bias or dropout.")

    re_call_parser.add_argument(
        "--processes",
        type=int,
        default=1,
        help="Number of processes to use when calling.")

    re_call_parser.add_argument(
        "--contig",
        type=str,
        default=argparse.SUPPRESS,
        help="Specifies a specific contig to process (optional).")

    re_call_parser.add_argument(
        "--sex-chr",
        type=str,
        help="Sex chromosome configuration to use for this sample (XX, XY, etc.) If left out, sex chromosomes will not "
             "be genotyped.")


CALC_CLASSES: Dict[str, Type[BaseCalculator]] = {
    c.CALLER_EXPANSIONHUNTER: ExpansionHunterCalculator,
    c.CALLER_GANGSTR: GangSTRCalculator,
    c.CALLER_REPEATHMM: RepeatHMMCalculator,
    c.CALLER_REPEATHMM_RECALL: RepeatHMMReCallCalculator,
    c.CALLER_STRAGLR: StraglrCalculator,
    c.CALLER_STRAGLR_RECALL: StraglrReCallCalculator,
    c.CALLER_STRKIT: StrKitCalculator,
    c.CALLER_TANDEM_GENOTYPES: TandemGenotypesCalculator,
    c.CALLER_TANDEM_GENOTYPES_RECALL: TandemGenotypesReCallCalculator,
}


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
        "--json",
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


def add_al_parser_args(al_parser):
    al_parser.add_argument("file", type=str, help="Alignment file to query.")
    al_parser.add_argument("region", type=str, help="Region to query aligned lengths for.")


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
    vs_parser.add_argument("--align-index", type=str, default="", help="Index for alignment file to visualize.")
    vs_parser.add_argument(
        "--ref", type=str, default="hg38",
        help="Reference genome code used for visualization and calls. Default: hg38")
    vs_parser.add_argument("--json", type=str, help="JSON file with STRkit calls to read from.")
    vs_parser.add_argument(
        "-i", type=int, default=1,
        help="1-based index of the locus to visualize in the JSON file. Default: 0")


def _exec_call(p_args) -> int:
    call_sample(
        tuple(p_args.read_files),
        p_args.ref,
        p_args.loci,
        min_reads=p_args.min_reads,
        min_allele_reads=p_args.min_allele_reads,
        num_bootstrap=p_args.num_bootstrap,
        flank_size=p_args.flank_size,
        sex_chroms=p_args.sex_chr,
        targeted=p_args.targeted,
        json_path=p_args.json,
        output_tsv=not p_args.no_tsv,
        processes=p_args.processes,
    )

    return 0


def _exec_re_call(p_args) -> int:
    contig: Optional[str] = getattr(p_args, "contig", None)

    return re_call_all_alleles(
        contig=contig,
        sex_chr=p_args.sex_chr,
        bootstrap_iterations=p_args.num_bootstrap,
        min_reads=p_args.min_reads,
        min_allele_reads=p_args.min_allele_reads,
        read_bias_corr_min=p_args.read_bias_corr_min,
        caller=p_args.caller,
        processes=p_args.processes,
    )


def _exec_mi(p_args) -> int:
    caller = p_args.caller.lower()

    trf_bed_file = getattr(p_args, "trf_bed") or None
    if trf_bed_file is None and caller in (c.CALLER_STRAGLR, c.CALLER_STRAGLR_RECALL):
        sys.stderr.write(f"Error: Using mistr with Straglr requires that the --trf-bed flag is used.\n")
        exit(1)

    calc_class: Optional[Type[BaseCalculator]] = CALC_CLASSES.get(caller)
    if not calc_class:
        sys.stderr.write(f"Error: Unknown or unimplemented caller '{caller}'\n")
        exit(1)

    child_file = getattr(p_args, "child-calls")  # First call file is not optional
    mother_file = getattr(p_args, "mother-calls", None)
    father_file = getattr(p_args, "father-calls", None)

    mother_id = getattr(p_args, "mother_id", None)
    father_id = getattr(p_args, "father_id", None)
    child_id = getattr(p_args, "child_id", None)

    # TODO: Check that caller supports "combined" call files

    if (father_file is None or mother_file is None) and (mother_id is None or father_id is None or child_id is None):
        print("Error: if less than 3 genotype files are specified, sample IDs must be given")
        exit(1)

    calc_inst = calc_class(
        child_call_file=child_file,
        mother_call_file=mother_file,
        father_call_file=father_file,

        child_id=child_id,
        mother_id=mother_id,
        father_id=father_id,

        loci_file=trf_bed_file,

        widen=getattr(p_args, "widen", 0) or 0,

        debug=p_args.debug,
    )

    # TODO: Figure out contigs properly... sex chromosome stuff
    contigs = calc_inst.get_trio_contigs(include_sex_chromosomes=False)

    contig = getattr(p_args, "contig", None)

    if contig is not None and contig not in contigs:
        print(f"Error: Could not find specified contig {p_args.contig} in trio contigs {contigs}")
        return 1

    if contig is not None:
        contigs = {contig}

    res = calc_inst.calculate(included_contigs=contigs)

    if not res:
        return 0

    if not p_args.no_tsv:
        print(str(res))
        if p_args.hist:
            print(res.histogram_text(bin_width=p_args.bin_width))
        print("---")
        sys.stdout.write(res.non_matching_tsv())
        sys.stdout.flush()

    if p_args.json:
        res.write_report_json(json_path=p_args.json, bin_width=p_args.bin_width)

    return 0


def _exec_combine_catalogs(p_args):
    return combine_catalogs(p_args.caller, p_args.paths)


def _exec_aligned_lengths(p_args):
    return aligned_lengths_cmd(p_args.file, p_args.region)


def _exec_convert(p_args):
    return convert(getattr(p_args, "trf-file"), p_args.caller)


def _exec_viz_server(p_args):
    align_file = str(pathlib.Path(p_args.align_file).resolve())
    align_type = os.path.splitext(align_file)[-1].lstrip(".")

    if align_type not in ("bam", "cram"):
        print(f"Error: file type '{align_type}' not supported", file=sys.stderr)
        return 1

    if p_args.align_index:
        align_index = str(pathlib.Path(p_args.align_index).resolve())
    else:
        align_index = f"{align_file}.{'crai' if align_type == 'cram' else 'bai'}"

    if not os.path.exists(align_index):
        print(f"Error: missing index at '{align_index}'", file=sys.stderr)
        return 1

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

    json_file = pathlib.Path(p_args.json)

    if not json_file.exists():
        print(f"Error: could not find JSON call file at '{json_file}'", file=sys.stderr)
        return 1

    with open(json_file, "r") as jf:
        call_data = json.load(jf)

    idx = p_args.i
    if idx < 1 or idx > len(call_data):
        print(f"Error: JSON offset out of bounds: '{idx}'", file=sys.stderr)
        return 1

    viz_run_server(
        call_data=call_data,
        initial_i=idx-1,
        ref=p_args.ref,
        # ref_index=ref_index,
        align=align_file,
        align_index=align_index,
        align_name=os.path.basename(align_file),
        align_index_name=os.path.basename(align_index),
        align_format=align_type,
    )

    return 0


def main(args: Optional[List[str]] = None):
    parser = argparse.ArgumentParser(
        description="A toolkit for analyzing variation in short(ish) tandem repeats.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers()

    call_parser = subparsers.add_parser(
        "call",
        help="A tandem repeat (TR) caller designed for high-fidelity long reads.")
    call_parser.set_defaults(func=_exec_call)
    add_call_parser_args(call_parser)

    re_call_parser = subparsers.add_parser(
        "re-call",
        help="A long-read tandem repeat (TR) re-caller, designed to build upon existing TR genotyping "
             "methods to yield calls with confidence intervals.")
    re_call_parser.set_defaults(func=_exec_re_call)
    add_re_call_parser_args(re_call_parser)

    mi_parser = subparsers.add_parser(
        "mi",
        help="A Mendelian inheritance calculator for different TR genotyping callers.")
    mi_parser.set_defaults(func=_exec_mi)
    add_mi_parser_args(mi_parser)

    cc_parser = subparsers.add_parser(
        "combine-catalogs",
        help="Combine Straglr result catalogs for use in re-calling with a consistent motif set.")
    cc_parser.set_defaults(func=_exec_combine_catalogs)
    add_cc_parser_args(cc_parser)

    al_parser = subparsers.add_parser(
        "aligned-lengths",
        help="See aligned lengths of (repeat) regions in a file to validate calls.")
    al_parser.set_defaults(func=_exec_aligned_lengths)
    add_al_parser_args(al_parser)

    cv_parser = subparsers.add_parser(
        "convert",
        help="Convert TRF BED file to other caller formats.")
    cv_parser.set_defaults(func=_exec_convert)
    add_cv_parser_args(cv_parser)

    vs_parser = subparsers.add_parser(
        "visualize",
        help="Start a web server to visualize aligned lengths of (repeat) a repeat region.")
    vs_parser.set_defaults(func=_exec_viz_server)
    add_vs_parser_args(vs_parser)

    args = args or sys.argv[1:]
    p_args = parser.parse_args(args)

    if not getattr(p_args, "func", None):
        p_args = parser.parse_args(("--help",))

    return p_args.func(p_args)


if __name__ == "__main__":
    exit(main())
