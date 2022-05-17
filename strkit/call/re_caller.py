import multiprocessing as mp
import sys

import strkit.constants as tc

from itertools import repeat
from typing import Optional

from .repeathmm import call_repeathmm
from .straglr import preprocess_lines_straglr, call_straglr
from .tandem_genotypes import call_tandem_genotypes

__all__ = [
    "re_call_all_alleles",
]


def _bound_param(param: int, min_val: int, max_val: int, flag_name: str):
    new_param = min(max(param, min_val), max_val)
    if new_param != param:
        sys.stderr.write(f"Warning: adjusting --{flag_name} to {new_param}")
    return new_param


def re_call_all_alleles(contig: Optional[str] = None,
                        sex_chr: Optional[str] = None,
                        bootstrap_iterations: int = 100,
                        min_reads: int = 4,
                        min_allele_reads: int = 2,
                        read_bias_corr_min: int = 4,
                        caller: str = tc.CALLER_TANDEM_GENOTYPES,
                        processes: int = 1) -> int:
    if caller not in tc.CALL_SUPPORTED_CALLERS:
        sys.stderr.write(f"Error: invalid caller '{caller}'")
        return 1

    n_proc = _bound_param(processes, 1, 512, "processes")
    min_reads = _bound_param(min_reads, 2, 512, "min-reads")
    min_allele_reads = _bound_param(min_allele_reads, 1, 512, "min-allele-reads")

    lines = [ls for ls in (line.strip() for line in sys.stdin) if ls]

    if caller == tc.CALLER_TANDEM_GENOTYPES:
        fn = call_tandem_genotypes

    elif caller == tc.CALLER_STRAGLR:
        fn = call_straglr

        # Need to group lines by locus since straglr does read-level stuff
        lines = preprocess_lines_straglr(lines)

    else:  # caller == tc.CALLER_REPEATHMM:
        fn = call_repeathmm

    args_iter = zip(
        repeat(contig),
        repeat(sex_chr),
        repeat(bootstrap_iterations),
        repeat(min_reads),
        repeat(min_allele_reads),
        repeat(read_bias_corr_min),
        lines,
    )

    sys.stderr.write(f"[DEBUG] Starting caller on {len(lines)} lines with {n_proc} processes\n")

    with mp.Pool(n_proc) as p:
        i = 0
        # noinspection PyTypeChecker
        for new_line in p.imap(fn, args_iter, chunksize=1):  # chunksize=1 seems fastest for some reason??
            sys.stdout.write(new_line)
            i += 1
            if i % 100 == 0:
                sys.stderr.write(f"[INFO] Processed {i} loci\n")
                sys.stderr.flush()
                sys.stdout.flush()

    return 0
