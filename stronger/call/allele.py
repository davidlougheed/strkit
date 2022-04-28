import numpy as np

from sklearn.exceptions import ConvergenceWarning
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import normalize
from warnings import simplefilter

from typing import Iterable, List, Optional, Tuple, Union

import stronger.constants as cc

__all__ = [
    "RepeatCounts",
    "get_n_alleles",
    "call_alleles",
]

RepeatCounts = Union[List[int], Tuple[int, ...], List[float], Tuple[float, ...]]


# K-means convergence errors - we expect convergence to some extent with homozygous alleles
simplefilter("ignore", category=ConvergenceWarning)


def _calculate_cis(samples, force_int: bool = False, ci: str = "95") -> np.array:
    # TODO: enable for numpy >=1.22
    r = {
        "95": (2.5, 97.5),
        "99": (0.5, 99.5),
    }[ci]
    # percentiles = np.percentile(samples, (2.5, 97.5), axis=1, method="interpolated_inverted_cdf")
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
def call_alleles(repeats_fwd: RepeatCounts,
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
                 force_int: bool) -> Tuple[Optional[np.array], Optional[np.array], Optional[np.array]]:
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
        return None, None, None

    allele_samples = np.array(
        [list() for _ in range(n_alleles)],
        dtype=np.int32 if force_int else np.float32)

    # Perform a number of bootstrap iterations to get a 95% CI and more accurate estimate of repeat counts / differences

    if separate_strands and fwd_len >= read_bias_corr_min and rev_len >= read_bias_corr_min:
        target_length: int = max(fwd_len, rev_len)

        # Resample original sample, correcting for imbalances between
        # forward and reverse-strand reads along the way
        # (if we've passed the coverage threshold)

        fwd_strand_samples = np.random.choice(
            fwd_strand_reads,
            size=(bootstrap_iterations, target_length),
            replace=True,
            p=fwd_strand_weights,
        )

        rev_strand_samples = np.random.choice(
            rev_strand_reads,
            size=(bootstrap_iterations, target_length),
            replace=True,
            p=rev_strand_weights,
        )

        concat_samples = np.sort(np.concatenate((fwd_strand_samples, rev_strand_samples), axis=1))

    else:
        concat_samples = np.sort(
            np.random.choice(
                combined_reads,
                size=(bootstrap_iterations, combined_len),
                replace=True,
                p=combined_weights,
            ))

    cache = {}

    for i in range(bootstrap_iterations):
        # Fit Gaussian mixture model to the resampled data

        sample = concat_samples[i, :]
        sample_t = tuple(sample)

        if sample_t in cache:
            g = cache[sample_t]
        else:
            g = GaussianMixture(
                n_components=n_alleles,
                init_params="kmeans",
                covariance_type="spherical",
                max_iter=100,
            ).fit(sample.reshape(-1, 1))
            cache[sample_t] = g

        means_and_weights = np.append(g.means_.transpose(), g.weights_.reshape(1, -1), axis=0)

        # Filter out peaks that aren't supported by ~min_allele_reads reads by probability, with some delta to allow
        # for peaks supported by "most of a read".
        mw_filter_1 = means_and_weights[1, :] > ((min_allele_reads - 0.1) / concat_samples.shape[0])

        # Filter out any peaks below some threshold using this magic constant filter factor
        mw_filter_2 = means_and_weights[1, :] > (1 / (gm_filter_factor * n_alleles))

        filtered_means_and_weights = means_and_weights[:, mw_filter_1 & mw_filter_2]

        if filtered_means_and_weights.shape[1] < 1:
            return None, None, None

        filtered_means = filtered_means_and_weights[0, :]
        filtered_weights = filtered_means_and_weights[1, :]

        resampled_means = np.append(
            filtered_means,
            np.random.choice(
                filtered_means,
                size=n_alleles - filtered_means.shape[0],
                p=normalize(filtered_weights.reshape(1, -1), norm="l1").flatten(),  # TODO: is this normalization ok?
            ),
        )

        # TODO: Re-sample means if any are removed, based on weights (re-normalized), to match total # of alleles

        # TODO: discard anything below a specific weight threshold and resample means based on remaining weights
        #  to fill in the gap. Maybe below 1 / (5 * num alleles) - i.e. 5 times less than we expect with equal sharing
        #  in the worst case where it represents just one allele

        sorted_allele_estimates = np.sort(resampled_means, axis=0).reshape(-1, 1)

        if not np.issubdtype(combined_reads[0], np.floating):  # TODO: Add force_int
            sorted_allele_estimates = np.rint(sorted_allele_estimates).astype(np.int32)

        allele_samples = np.append(allele_samples, sorted_allele_estimates, axis=1)

    # Calculate 95% and 99% confidence intervals for each allele from the bootstrap distributions.
    allele_samples.sort(axis=1)
    allele_cis_95 = _calculate_cis(allele_samples, force_int=force_int, ci="95")
    allele_cis_99 = _calculate_cis(allele_samples, force_int=force_int, ci="99")

    # TODO: Calculate CIs based on Gaussians from allele samples instead? Ask someone...
    #  - Could take median of 2.5 percentiles and 97.5 percentiles from Gaussians instead, median of means

    # Report the median estimates (TODO: ???)
    # and the confidence intervals.

    # TODO: enable for numpy >=1.22
    # median = np.percentile(allele_samples, 50, axis=1, method="interpolated_inverted_cdf")
    median = np.percentile(allele_samples, 50, axis=1, interpolation="nearest")
    if force_int:
        median = np.rint(median).astype(np.int32)

    return median.flatten(), allele_cis_95, allele_cis_99
