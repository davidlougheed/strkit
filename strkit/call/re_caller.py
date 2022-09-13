from __future__ import annotations

import multiprocessing as mp
import numpy as np
import sys

import strkit.constants as tc

from itertools import repeat
from typing import Optional

from strkit.exceptions import ParamError
from strkit.logger import logger
from .repeathmm import call_repeathmm
from .straglr import preprocess_lines_straglr, call_straglr
from .tandem_genotypes import call_tandem_genotypes

__all__ = [
    "re_call_all_alleles",
]


def _bound_param(param: int, min_val: int, max_val: int, flag_name: str):
    new_param = min(max(param, min_val), max_val)
    if new_param != param:
        logger.warning(f"Adjusting --{flag_name} to {new_param}")
    return new_param


def re_call_all_alleles(
    contig: Optional[str] = None,
    sex_chr: Optional[str] = None,
    bootstrap_iterations: int = 100,
    min_reads: int = 4,
    min_allele_reads: int = 2,
    read_bias_corr_min: int = 4,
    caller: str = tc.CALLER_TANDEM_GENOTYPES,
    processes: int = 1,
    seed: Optional[int] = None,
) -> None:
    if caller not in tc.CALL_SUPPORTED_CALLERS:
        raise ParamError(f"Invalid caller '{caller}'")

    n_proc = _bound_param(processes, 1, 512, "processes")
    min_reads = _bound_param(min_reads, 2, 512, "min-reads")
    min_allele_reads = _bound_param(min_allele_reads, 1, 512, "min-allele-reads")

    # Seed the random number generator if a seed is provided, for replicability
    rng = np.random.default_rng(seed=seed)

    def _rng_gen():
        while True:
            yield rng.integers(0, 4096).item()

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
        _rng_gen(),
        lines,
    )

    logger.debug(f"Starting caller on {len(lines)} lines with {n_proc} processes")

    with mp.Pool(n_proc) as p:
        i = 0
        # noinspection PyTypeChecker
        for new_line in p.imap(fn, args_iter, chunksize=1):  # chunksize=1 seems fastest for some reason??
            sys.stdout.write(new_line)
            i += 1
            if i % 100 == 0:
                logger.info(f"Processed {i} loci")
                sys.stderr.flush()
                sys.stdout.flush()
