from __future__ import annotations

# Disable OpenMP multithreading since it adds enormous overhead when multiprocessing
import os
os.environ["OMP_NUM_THREADS"] = "1"

# ----------------------------------------------------------------------------------------------------------------------

import numpy as np
import statistics

from sklearn.exceptions import ConvergenceWarning
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import normalize
from warnings import simplefilter

from typing import Iterable, Optional, Union

import strkit.constants as cc

__all__ = [
    "RepeatCounts",
    "get_n_alleles",
    "call_alleles",
]

RepeatCounts = Union[list[int], tuple[int, ...], list[float], tuple[float, ...]]


# K-means convergence errors - we expect convergence to some extent with homozygous alleles
simplefilter("ignore", category=ConvergenceWarning)

# TODO: parameterize
small_allele_min = 8
expansion_ratio = 5
N_GM_INIT = 3


def _calculate_cis(samples, force_int: bool = False, ci: str = "95") -> np.array:
    r = {
        "95": (2.5, 97.5),
        "99": (0.5, 99.5),
    }[ci]
    # TODO: enable for numpy >=1.22
    # percentiles = np.percentile(samples, r, axis=1, method="interpolated_inverted_cdf").transpose()
    percentiles = np.percentile(samples, r, axis=1, interpolation="nearest").transpose()
    return np.rint(percentiles).astype(np.int32) if force_int else percentiles


def get_n_alleles(default_n_alleles: int, sample_sex_chroms: Optional[str], contig: str) -> Optional[int]:
    if contig in cc.M_CHROMOSOME_NAMES:
        return 1

    if contig in cc.SEX_CHROMOSOMES:
        if sample_sex_chroms is None:
            return None
        if contig in cc.X_CHROMOSOME_NAMES:
            return sample_sex_chroms.count("X")
        if contig in cc.Y_CHROMOSOME_NAMES:
            return sample_sex_chroms.count("Y")

    return default_n_alleles


# noinspection PyUnresolvedReferences
def call_alleles(
    repeats_fwd: RepeatCounts,
    repeats_rev: RepeatCounts,
    read_weights_fwd: Optional[Iterable[float]],
    read_weights_rev: Optional[Iterable[float]],
    bootstrap_iterations: int,
    min_reads: int,
    min_allele_reads: int,
    n_alleles: int,
    separate_strands: bool,
    read_bias_corr_min: int,
    gm_filter_factor: int,
    hq: bool,
    force_int: bool,
    seed: Optional[int],
) -> Optional[dict]:
    fwd_strand_reads = np.array(repeats_fwd)
    rev_strand_reads = np.array(repeats_rev)

    fwd_len = fwd_strand_reads.shape[0]
    rev_len = rev_strand_reads.shape[0]

    fwd_strand_weights = np.array(
        read_weights_fwd if read_weights_fwd is not None else np.array(([1/fwd_len] * fwd_len) if fwd_len else []))
    rev_strand_weights = np.array(
        read_weights_rev if read_weights_rev is not None else np.array(([1/rev_len] * rev_len) if rev_len else []))

    assert fwd_strand_reads.shape == fwd_strand_weights.shape
    assert rev_strand_reads.shape == rev_strand_weights.shape

    combined_reads = np.concatenate((fwd_strand_reads, rev_strand_reads), axis=None)
    combined_weights = np.concatenate((fwd_strand_weights, rev_strand_weights), axis=None)
    combined_len = combined_reads.shape[-1]

    if combined_len < min_reads:
        return None

    def _na_length_list():
        return [list() for _ in range(n_alleles)]

    allele_samples = np.array(_na_length_list(), dtype=np.float32)
    allele_weight_samples = np.array(_na_length_list(), dtype=np.float32)
    allele_stdev_samples = np.array(_na_length_list(), dtype=np.float32)
    sample_peaks = np.array([], dtype=np.int32)

    rng: np.random.Generator = np.random.default_rng(seed=seed)

    # Perform a number of bootstrap iterations to get a 95% CI and more accurate estimate of repeat counts / differences

    if separate_strands and fwd_len >= read_bias_corr_min and rev_len >= read_bias_corr_min:
        target_length: int = max(fwd_len, rev_len)

        # Resample original sample, correcting for imbalances between
        # forward and reverse-strand reads along the way
        # (if we've passed the coverage threshold)

        fwd_strand_samples = rng.choice(
            fwd_strand_reads,
            size=(bootstrap_iterations, target_length),
            replace=True,
            p=fwd_strand_weights,
        )

        rev_strand_samples = rng.choice(
            rev_strand_reads,
            size=(bootstrap_iterations, target_length),
            replace=True,
            p=rev_strand_weights,
        )

        concat_samples = np.sort(
            np.concatenate((fwd_strand_samples, rev_strand_samples), axis=1),
            kind="stable")

    else:
        concat_samples = np.sort(
            rng.choice(
                combined_reads,
                size=(bootstrap_iterations, combined_len),
                replace=True,
                p=combined_weights,
            ) if bootstrap_iterations > 1 else np.array([combined_reads]),
            kind="stable")

    cache = {}

    for i in range(bootstrap_iterations):
        # Fit Gaussian mixture model to the resampled data

        sample = concat_samples[i, :]
        sample_t = tuple(sample)
        sample_rs = sample.reshape(-1, 1)

        g = None

        if sample_t in cache:
            g = cache[sample_t]
        else:
            n_components = n_alleles
            while n_components > 0:
                g = GaussianMixture(
                    n_components=n_components,
                    # init_params="kmeans",  TODO: parameterize
                    init_params="k-means++",
                    covariance_type="spherical",
                    n_init=N_GM_INIT,
                    random_state=rng.integers(0, 4096).item(),
                ).fit(sample_rs)

                means_and_weights = np.append(g.means_.transpose(), g.weights_.reshape(1, -1), axis=0)

                # Filter out peaks that aren't supported by ~min_allele_reads reads by probability, with some delta to
                # allow for peaks supported by "most of a read".
                mw_filter_1 = means_and_weights[1, :] > ((min_allele_reads - 0.1) / concat_samples.shape[0])

                # Filter out any peaks below some threshold using this magic constant filter factor
                # - Exception: Large expansions can have very few supporting reads due to quirks of sequencing beyond
                #   just chance/read length distribution; if we have 2 alleles and the large one is a lot bigger than
                #   the small one, don't apply this filter
                # - Discard anything below a specific weight threshold and resample means based on remaining weights
                #   to fill in the gap. E.g. below 1 / (5 * num alleles) - i.e. 5 times less than we expect with equal
                #   sharing in the worst case where it represents just one allele
                if n_components > 2 or (n_components == 2 and (not hq or (
                        means_and_weights[0, -1] < expansion_ratio * max(means_and_weights[0, 0],
                                                                         small_allele_min)))):
                    mw_filter_2 = means_and_weights[1, :] > (1 / (gm_filter_factor * n_alleles))
                else:
                    mw_filter_2 = means_and_weights[1, :] > np.finfo(np.float32).eps

                mw_filter = mw_filter_1 & mw_filter_2
                n_useless = np.size(mw_filter) - np.count_nonzero(mw_filter)
                if not n_useless:
                    # No useless components left to remove, so escape
                    break
                n_components -= n_useless

            cache[sample_t] = g

        if not g:
            # Could not fit any Gaussian mixture; skip this allele
            return None

        # Keep track of how many alleles were found for
        sample_peaks = np.append(sample_peaks, g.means_.shape[0])

        means_and_weights = np.append(g.means_.transpose(), g.weights_.reshape(1, -1), axis=0)

        means = means_and_weights[0, :]
        weights = means_and_weights[1, :]
        stdevs = np.sqrt(g.covariances_)
        n_to_resample = n_alleles - means.shape[0]

        if n_to_resample:
            # Re-sample means if any are removed, based on weights (re-normalized), to match total # of alleles
            resampled_indices = rng.choice(
                np.arange(len(means)),
                size=n_to_resample,
                p=normalize(weights.reshape(1, -1), norm="l1").flatten())
            resampled_means = np.append(means, means[resampled_indices])
            resampled_weights = np.append(weights, weights[resampled_indices])
            resampled_stdevs = np.append(stdevs, stdevs[resampled_indices])
        else:
            resampled_means = means
            resampled_weights = weights
            resampled_stdevs = stdevs

        argsorted_means = np.argsort(resampled_means, axis=0, kind="stable")
        sorted_allele_estimates = resampled_means[argsorted_means].reshape(-1, 1)
        sorted_allele_weight_estimates = resampled_weights[argsorted_means].reshape(-1, 1)
        sorted_allele_stdev_estimates = resampled_stdevs[argsorted_means].reshape(-1, 1)

        allele_samples = np.append(allele_samples, sorted_allele_estimates, axis=1)
        allele_weight_samples = np.append(allele_weight_samples, sorted_allele_weight_estimates, axis=1)
        allele_stdev_samples = np.append(allele_stdev_samples, sorted_allele_stdev_estimates, axis=1)

    # Calculate 95% and 99% confidence intervals for each allele from the bootstrap distributions.
    allele_samples.sort(axis=1, kind="stable")
    allele_cis_95 = _calculate_cis(allele_samples, force_int=force_int, ci="95")
    allele_cis_99 = _calculate_cis(allele_samples, force_int=force_int, ci="99")
    allele_weight_samples.sort(axis=1, kind="stable")
    allele_stdev_samples.sort(axis=1, kind="stable")
    sample_peaks.sort(kind="stable")  # To make mode consistent, given same set of peak #s

    # TODO: Calculate CIs based on Gaussians from allele samples instead? Ask someone...
    #  - Could take median of 2.5 percentiles and 97.5 percentiles from Gaussians instead, median of means

    # Report the median estimates (TODO: ???)
    # and the confidence intervals.

    medians_of_means = np.percentile(allele_samples, 50, axis=1, method="interpolated_inverted_cdf")
    medians_of_means_final = medians_of_means
    if force_int:
        medians_of_means_final = np.rint(medians_of_means).astype(np.int32)
    medians_of_weights = np.percentile(allele_weight_samples, 50, axis=1, interpolation="nearest")
    medians_of_stdevs = np.percentile(allele_stdev_samples, 50, axis=1, interpolation="nearest")
    modal_n_peaks = statistics.mode(sample_peaks).item()

    return {
        "call": medians_of_means_final.flatten(),
        "call_95_cis": allele_cis_95,
        "call_99_cis": allele_cis_99,

        "peaks": medians_of_means.flatten(),  # Don't round, so we can recover original Gaussian model
        # TODO: Do we want to report this, or the stdev of the final peaks as determined by the bootstrapped
        #  means and variances????????????????
        "peak_weights": medians_of_weights.flatten(),
        "peak_stdevs": medians_of_stdevs.flatten(),
        # TODO: should be ok to use this, because resample gets put at end, vertically (3rd allele in a 3-ploid case)
        #  so taking the first 2 alleles still works in terms of stdev/mean estimates? I think?
        #  Not quite, cause it's sorted...
        #  --> Only do the peak assignment with 1/2 peaks, which is the majority of human situations
        "modal_n_peaks": modal_n_peaks,
    }
