from __future__ import annotations

# Disable OpenMP/other multithreading since it adds enormous overhead when multiprocessing
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

# ----------------------------------------------------------------------------------------------------------------------

import logging  # For type hinting
import numpy as np
import statistics

from sklearn.exceptions import ConvergenceWarning
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import normalize
from warnings import simplefilter

from numpy.typing import NDArray
from typing import Iterable, Literal, Optional, TypedDict, Union

import strkit.constants as cc

from .params import CallParams
from .utils import get_new_seed

__all__ = [
    "RepeatCounts",
    "CallDict",
    "get_n_alleles",
    "call_alleles",
]

RepeatCounts = Union[list[int], tuple[int, ...], list[float], tuple[float, ...], NDArray[np.int_], NDArray[np.float_]]


# K-means convergence errors - we expect convergence to some extent with homozygous alleles
simplefilter("ignore", category=ConvergenceWarning)

# TODO: parameterize
small_allele_min = 8
expansion_ratio = 5
N_GM_INIT = 3

WEIGHT_1_0 = np.array([[1.0]])
FLOAT_32_EPSILON = np.finfo(np.float32).eps

CI_PERCENTILE_RANGES = {
    "95": (2.5, 97.5),
    "99": (0.5, 99.5),
}


def _array_as_int(n: Union[NDArray[np.int_], NDArray[np.float_]]) -> NDArray[np.int32]:
    return np.rint(n).astype(np.int32)


def _calculate_cis(samples, ci: str = Literal["95", "99"]) -> Union[NDArray[np.int32], NDArray[np.float_]]:
    percentiles = np.percentile(
        samples, CI_PERCENTILE_RANGES[ci], axis=1, method="interpolated_inverted_cdf"
    ).transpose()
    return _array_as_int(percentiles)


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


def na_length_list(n_alleles: int):
    return [list() for _ in range(n_alleles)]


def fit_gmm(
    rng: np.random.Generator,
    sample: np.array,
    n_alleles: int,
    allele_filter: float,
    hq: bool,
    gm_filter_factor: int,
    init_params: Literal["kmeans", "k-means++"] = "k-means++",  # TODO: parameterize outside
) -> Optional[object]:
    sample_rs = sample.reshape(-1, 1)
    g: Optional[object] = None

    n_components = n_alleles
    while n_components > 0:
        if n_components == 1:  # Don't need to do the full fit for a single peak, just calculate the parameters
            # I've confirmed this gives an ~identical result to fitting a GMM with one parameter.
            fake_g: object = type("", (), {})()
            fake_g.means_ = np.array([[np.mean(sample_rs)]])
            fake_g.weights_ = WEIGHT_1_0
            fake_g.covariances_ = np.array([[np.var(sample_rs)]])
            return fake_g

        g = GaussianMixture(
            n_components=n_components,
            init_params=init_params,
            covariance_type="spherical",
            n_init=N_GM_INIT,
            random_state=get_new_seed(rng),
        ).fit(sample_rs)

        # noinspection PyUnresolvedReferences
        means_and_weights = np.append(g.means_.transpose(), g.weights_.reshape(1, -1), axis=0)

        # Filter out peaks that aren't supported by ~min_allele_reads reads by probability, with some delta to
        # allow for peaks supported by "most of a read".
        mw_filter_1 = means_and_weights[1, :] > allele_filter

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
            mw_filter_2 = means_and_weights[1, :] > (1 / (gm_filter_factor * n_components))
        else:
            mw_filter_2 = means_and_weights[1, :] > FLOAT_32_EPSILON

        mw_filter = mw_filter_1 & mw_filter_2
        n_useless = np.size(mw_filter) - np.count_nonzero(mw_filter)
        if not n_useless:
            # No useless components left to remove, so return the GMM
            return g
        n_components -= n_useless

    return g


class CallDict(TypedDict):
    call: Union[NDArray[np.int32], NDArray[np.float_]]
    call_95_cis: Union[NDArray[np.int32], NDArray[np.float_]]  # 2D arrays
    call_99_cis: Union[NDArray[np.int32], NDArray[np.float_]]  # 2D arrays
    peaks: NDArray[np.float_]
    peak_weights: NDArray[np.float_]
    peak_stdevs: NDArray[np.float_]
    modal_n_peaks: int


def make_read_weights(read_weights: Optional[Iterable[float]], num_reads: int) -> NDArray[np.float_]:
    return np.array(
        read_weights if read_weights is not None else np.array(([1/num_reads] * num_reads) if num_reads else []))


def call_alleles(
    repeats_fwd: RepeatCounts,
    repeats_rev: RepeatCounts,
    read_weights_fwd: Optional[Iterable[float]],
    read_weights_rev: Optional[Iterable[float]],
    params: CallParams,
    min_reads: int,
    n_alleles: int,
    separate_strands: bool,
    read_bias_corr_min: int,
    gm_filter_factor: int,
    seed: Optional[int],
    logger_: logging.Logger,
    debug_str: str,
) -> Optional[CallDict]:
    fwd_strand_reads = np.array(repeats_fwd)
    rev_strand_reads = np.array(repeats_rev)

    fwd_len = fwd_strand_reads.shape[0]
    rev_len = rev_strand_reads.shape[0]

    fwd_strand_weights = make_read_weights(read_weights_fwd, fwd_len)
    rev_strand_weights = make_read_weights(read_weights_rev, rev_len)

    assert fwd_strand_reads.shape == fwd_strand_weights.shape
    assert rev_strand_reads.shape == rev_strand_weights.shape

    combined_reads = np.concatenate((fwd_strand_reads, rev_strand_reads), axis=None)
    combined_weights = np.concatenate((fwd_strand_weights, rev_strand_weights), axis=None)
    combined_len = combined_reads.shape[-1]

    if combined_len < min_reads:
        return None

    # If the locus/allele only has one value, don't bother bootstrapping
    if np.unique(combined_reads).shape[0] == 1:
        logger_.debug(f"{debug_str} - skipping bootstrap / GMM fitting for allele(s) (single value)")
        cn = combined_reads[0]

        call = _array_as_int(np.full(n_alleles, cn))
        call_cis = _array_as_int(np.full((n_alleles, 2), cn))

        peaks: NDArray[np.float_] = call.astype(np.float_)

        return {
            "call": call,
            "call_95_cis": call_cis,
            "call_99_cis": call_cis,
            "peaks": peaks,
            "peak_weights": np.full(n_alleles, 1.0 / n_alleles),
            "peak_stdevs": np.full(n_alleles, 0.0),
            "modal_n_peaks": 1,  # 1 peak, since we have 1 value
        }

    nal = na_length_list(n_alleles)
    allele_samples = np.array(nal, dtype=np.float32)
    allele_weight_samples = np.array(nal, dtype=np.float32)
    allele_stdev_samples = np.array(nal, dtype=np.float32)
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
            size=(params.num_bootstrap, target_length),
            replace=True,
            p=fwd_strand_weights,
        )

        rev_strand_samples = rng.choice(
            rev_strand_reads,
            size=(params.num_bootstrap, target_length),
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
                size=(params.num_bootstrap, combined_len),
                replace=True,
                p=combined_weights,
            ) if params.num_bootstrap > 1 else np.array([combined_reads]),
            kind="stable")

    gmm_cache = {}

    def _get_fitted_gmm(s: Union[NDArray[np.int_], NDArray[np.float_]]) -> Optional[object]:
        if (s_t := tuple(s)) not in gmm_cache:
            # Fit Gaussian mixture model to the resampled data
            gmm_cache[s_t] = fit_gmm(rng, s, n_alleles, allele_filter, params.hq, gm_filter_factor)

        return gmm_cache[s_t]

    # Filter out peaks that aren't supported by ~min_allele_reads reads by probability, with some delta to
    # allow for peaks supported by "most of a read".
    allele_filter = (params.min_allele_reads - 0.1) / concat_samples.shape[0]

    for i in range(params.num_bootstrap):
        sample = concat_samples[i, :]

        g: Optional[object] = _get_fitted_gmm(sample)
        if not g:
            # Could not fit any Gaussian mixture; skip this allele
            return None

        # Keep track of how many alleles were found for
        # noinspection PyUnresolvedReferences
        sample_peaks = np.append(sample_peaks, g.means_.shape[0])

        # noinspection PyUnresolvedReferences
        means_and_weights = np.append(g.means_.transpose(), g.weights_.reshape(1, -1), axis=0)

        means = means_and_weights[0, :]
        weights = means_and_weights[1, :]
        # noinspection PyUnresolvedReferences
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
    allele_samples_argsort = allele_samples.argsort(axis=1, kind="stable")
    allele_samples = np.take_along_axis(allele_samples, allele_samples_argsort, axis=1)
    allele_cis_95 = _calculate_cis(allele_samples, ci="95")
    allele_cis_99 = _calculate_cis(allele_samples, ci="99")
    allele_weight_samples = np.take_along_axis(allele_weight_samples, allele_samples_argsort, axis=1)
    allele_stdev_samples = np.take_along_axis(allele_stdev_samples, allele_samples_argsort, axis=1)

    sample_peaks.sort(kind="stable")  # To make mode consistent, given same set of peak #s

    # TODO: Calculate CIs based on Gaussians from allele samples instead? Ask someone...
    #  - Could take median of 2.5 percentiles and 97.5 percentiles from Gaussians instead, median of means

    # Report the median estimates and the confidence intervals.
    #  - we choose nearest for median rather than interpolating, so we can get real corresponding weights and stdevs.

    median_idx = allele_samples.shape[1] // 2  #
    medians_of_means = allele_samples[:, median_idx]
    medians_of_means_final = np.rint(medians_of_means).astype(np.int32)
    peak_weights = allele_weight_samples[:, median_idx]
    peak_stdevs = allele_stdev_samples[:, median_idx]
    modal_n_peaks: int = statistics.mode(sample_peaks).item()

    return {
        "call": medians_of_means_final.flatten(),
        "call_95_cis": allele_cis_95,
        "call_99_cis": allele_cis_99,

        "peaks": medians_of_means.flatten(),  # Don't round, so we can recover original Gaussian model
        "peak_weights": peak_weights.flatten(),
        "peak_stdevs": peak_stdevs.flatten(),
        # TODO: should be ok to use this, because resample gets put at end, vertically (3rd allele in a 3-ploid case)
        #  so taking the first 2 alleles still works in terms of stdev/mean estimates? I think?
        #  Not quite, cause it's sorted...
        #  --> Only do the peak assignment with 1/2 peaks, which is the majority of human situations
        "modal_n_peaks": modal_n_peaks,
    }
