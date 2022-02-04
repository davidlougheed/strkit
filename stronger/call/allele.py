import numpy as np

from sklearn.exceptions import ConvergenceWarning
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import normalize
from warnings import simplefilter

from typing import List, Optional, Tuple, Union

Alleles = Union[List[int], Tuple[int, ...], List[float], Tuple[float, ...]]


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


# noinspection PyUnresolvedReferences
def call_allele(allele_1: Alleles,
                allele_2: Alleles,
                bootstrap_iterations: int,
                min_reads: int,
                min_allele_reads: int,
                n_alleles: int,
                separate_strands: bool,
                read_bias_corr_min: int,
                gm_filter_factor: int,
                force_int: bool) -> Tuple[Optional[np.array], Optional[np.array], Optional[np.array]]:
    fwd_strand_reads = np.array(allele_1)
    rev_strand_reads = np.array(allele_2)

    combined = np.concatenate((fwd_strand_reads, rev_strand_reads), axis=None)
    combined_len = combined.shape[0]

    if combined_len < min_reads:
        return None, None, None

    allele_samples = np.array(
        [list() for _ in range(n_alleles)],
        dtype=np.int32 if force_int else np.float32)

    # Perform a number of bootstrap iterations to get a 95% CI and more accurate estimate of repeat counts / differences

    for _ in range(bootstrap_iterations):
        fwd_len = fwd_strand_reads.shape[0]
        rev_len = rev_strand_reads.shape[0]
        target_length: int = max(fwd_len, rev_len)

        # Resample original sample, correcting for imbalances between
        # forward and reverse-strand reads along the way
        # (if we've passed the coverage threshold)
        if separate_strands and fwd_len >= read_bias_corr_min and rev_len >= read_bias_corr_min:
            fwd_strand_sample = np.random.choice(fwd_strand_reads, size=target_length, replace=True)
            rev_strand_sample = np.random.choice(rev_strand_reads, size=target_length, replace=True)

            print(fwd_strand_sample, rev_strand_sample)

            concat_samples = np.concatenate((fwd_strand_sample, rev_strand_sample), axis=None)
        else:
            concat_samples = np.random.choice(combined, size=target_length*2, replace=True)

        # Fit Gaussian mixture model to the resampled data
        g = GaussianMixture(
            n_components=n_alleles,
            init_params="kmeans",
            # weights_init=[1/n_alleles]*n_alleles,
            #
            # TODO
            # We assume (under no mosaicism hypothesis) that the error model is the same between repeat sizes
            # There is probably a better way to do this, but otherwise it'll basically never call homozygous alleles.
            # covariance_type="tied",

            # n_init=10,
            # means_init=np.random.choice(concat_samples, size=n_alleles).reshape(-1, 1),
            max_iter=100,
        ).fit(concat_samples.reshape(-1, 1))

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

        if not np.issubdtype(combined[0], np.floating):  # TODO: Add force_int
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
