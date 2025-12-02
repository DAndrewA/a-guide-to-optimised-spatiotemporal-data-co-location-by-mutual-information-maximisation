"""Author: Andrew Martin,
Creation date: 2/12/25

Script handling common functionality for generating samples from two polynomial-based distributions, based on cumulative inverse sampling.
"""

import numpy as np
from scipy.special import erf as ERF

# define the polynomial functions for the normalised probability distributions on the interval [0,1]
F_a_of_x = np.polynomial.Polynomial([0,0,-6,8,-3])
F_b_of_x = np.polynomial.Polynomial([0,0,-18,52,-55,20])

def x_of_u_given_F_of_x(u: np.ndarray, F_of_x: np.polynomial.Polynomial, imag_threshold: float = 0.01) -> (np.ndarray, np.ndarray):
    """Given a polynomial F_of_x = f(x), and values of u=F(x), evaluate x=F^{-1}(x)."""
    possible_x_of_u = np.asarray([
        (np.polynomial.Polynomial([u_i]) + F_of_x).roots()
        for u_i in u
    ])

    in_correct_real_range = (possible_x_of_u.real >= 0) & (possible_x_of_u.real <= 1)
    in_correct_imag_range = (np.abs(possible_x_of_u.imag) < imag_threshold)

    indices = (in_correct_real_range & in_correct_imag_range)

    # sometimes, two roots of a polynomial saify the requirements simultaneously => they are likely degenerate roots.
    # Only take the first root in this case.
    first_true_idx = np.argmax(indices, axis=1)
    has_true = np.any(indices, axis=1)
    indices_filtered = np.zeros_like(indices, dtype=bool)
    rows_with_true = np.where(has_true)[0]
    indices_filtered[rows_with_true, first_true_idx[rows_with_true]] = True
        
    x_of_u = possible_x_of_u[indices_filtered].real
    return u[indices_filtered.any(axis=1)], x_of_u


def gen_uv_samples_mixture(n_samples: int, rhostar: float, mixture: float, rng: np.random.Generator | None = None):
    """Given a number of samples, a covariance, and a mixture ratio, generate samples proportionally mixing independent and gaussian dependent samples, and inverse transform them to obtain samples with uniform marginal distributions.
    """
    assert rhostar >= -1 and rhostar <= 1
    assert 0 <= mixture and mixture <= 1
    assert isinstance(n_samples, int)
    if rng is None:
        rng = np.random.default_rng()

    n_gaussian = int(np.ceil(n_samples * mixture))

    cov_gaussian = np.array([
        [1, rhostar],
        [rhostar, 1]
    ])
    xy_gaussian = rng.multivariate_normal(mean=[0,0], cov=cov_gaussian, size=n_gaussian)
    uv_gaussian = 0.5*( 1 + ERF(xy_gaussian / np.sqrt(2)) )

    n_independent = n_samples - n_gaussian
    uv_independent = rng.uniform(low=0, high=1, size=(n_independent, 2))

    uv = np.concat([uv_gaussian, uv_independent])
    return uv


def gen_xy_samples_mixture(n_samples: int, rhostar: float, mixture: float, rng: np.random.Generator | None = None):
    """Given a number of samples, covariance and mixture ratio, generate samples from inverse sampling F_[a,b]_of_x.
    """
    uv = gen_uv_samples_mixture(
        n_samples=n_samples,
        rhostar=rhostar,
        mixture=mixture,
        rng=rng,
    )
    u = uv[:,0]
    v = uv[:,1]

    u_a, x_a = x_of_u_given_F_of_x(u=u, F_of_x=F_a_of_x)
    v_b, y_b = x_of_u_given_F_of_x(u=v, F_of_x=F_b_of_x)

    return (u_a, v_b), (x_a, y_b)
