from __future__ import annotations

# Disable OpenMP/other multithreading since it adds enormous overhead when multiprocessing
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

# ----------------------------------------------------------------------------------------------------------------------

import numpy as np

from dataclasses import dataclass
from numpy.typing import NDArray
from numpy.random import Generator
from sklearn.exceptions import ConvergenceWarning
from sklearn.mixture import GaussianMixture
from typing import Literal
from warnings import simplefilter

from .utils import get_new_seed

__all__ = [
    "GMMInitParamsMethod",
    "GMMParams",
    "make_single_gaussian",
    "make_already_fitted_gmm",
]


# K-means convergence errors - we expect convergence to some extent with homozygous alleles
simplefilter("ignore", category=ConvergenceWarning)


GMMInitParamsMethod = Literal["kmeans", "k-means++"]

WEIGHT_1_0 = np.array([[1.0]])


@dataclass(frozen=True)
class GMMParams:
    # sklearn model parameters
    init_params_method: GMMInitParamsMethod
    n_init: int
    # other parameters for fit_gmm
    #  - for deciding pre-GMM-fit to only have one peak based on read copy number frequency patterns
    pre_filter_factor: int
    #  - for filtering out peaks which have a weight lower than (filter_factor * n_components) for a GMM of n_components
    filter_factor: int

    def make_fitted_gmm(self, n_components: int, sample_rs: NDArray, rng: Generator):
        return GaussianMixture(
            n_components=n_components,
            init_params=self.init_params_method,
            covariance_type="spherical",
            # For relatively simple cases (total of two different copy numbers in the sample set), force n_init to be 1
            # even if we have a different parameter value to save on needless computation.
            n_init=self.n_init,  # if sample_rs[-1][0] - sample_rs[0][0] > 1.0 else 1,  TODO: re-enable if ok
            random_state=get_new_seed(rng),
        ).fit(sample_rs)


def make_single_gaussian(sample_rs: NDArray) -> object:
    # Don't need to do the full fit for a single peak, just calculate the parameters.
    # I've confirmed this gives an ~identical result to fitting a GMM with one parameter.
    fake_g: object = type("", (), {})()
    fake_g.means_ = np.array([[np.mean(sample_rs)]])
    fake_g.weights_ = WEIGHT_1_0
    fake_g.covariances_ = np.array([[np.var(sample_rs)]])
    return fake_g


def make_already_fitted_gmm(
    n_components: int, peaks: NDArray[np.float_], stdevs: NDArray[np.float_], weights: NDArray[np.float_]
) -> GaussianMixture:
    g: GaussianMixture = GaussianMixture(n_components=n_components, covariance_type="spherical")
    g.means_ = peaks.reshape(-1, 1)
    g.covariances_ = stdevs ** 2
    g.weights_ = weights
    # https://github.com/scikit-learn/scikit-learn/blob/1.5.0/sklearn/mixture/_gaussian_mixture.py#L347
    g.precisions_cholesky_ = 1.0 / stdevs
    return g
