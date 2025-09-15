"""Author: Andrew Martin
Creation date: 3/7/25

Script handling conversion of SLURM array indices to pairs of (R, tau)
"""

from dataclasses import dataclass
import numpy as np
import datetime as dt

INDEX_FUNCTIONS = dict()


def register_index_function(f):
    INDEX_FUNCTIONS[f.__name__] = f
    return f


@dataclass(frozen=True, kw_only=True)
class Parametrisation:
    distance_km: float
    tau: dt.timedelta


class InvalidIndexError(ValueError):
    pass


@register_index_function
def original(index: int) -> Parametrisation:
    """Original indexing function used in scripts until 3/7/25
    VALID RANGE: [0,449]
    """
    if index < 0 or index > 449:
        raise InvalidIndexError(index)

    R_km = int(np_arange(5,151,5)[index % 30])
    tau_s = int(np_arange(300,4501,300)[index // 30])
    return Parametrisation(
        distance_km = float(R_km),
        tau = dt.timedelta(seconds=int(tau_s))
    )


@register_index_function
def extended_tau_meshgrid(index: int) -> Parametrisation:
    """Extends the range of tau values used with varying resolution along the tau axis
    VALID RANGE: [0,492]
    """
    if index < 0 or index >= 493:
        raise InvalidIndexError(index)
    R = np.concatenate([
        np.arange(5,50,5),
        np.arange(50,100,10),
        np.arange(100,150,25),
        [150]
    ])
    tau = np.concatenate([
        [300],
        np.arange(600,4800, 600),
        np.arange(4800, 21600, 1200),
        np.arange(21600, 43200, 3600),
        [43200]
    ])
    R_, tau_ = np.meshgrid(R, tau)
    [R_km, tau_s] = np.array([
        R_.flatten(),
        tau_.flatten()
    ])[:, index]
    return Parametrisation(
        distance_km = float(R_km),
        tau = dt.timedelta(seconds=int(tau_s))
    )


@register_index_function
def R_150km_tau_172800s(index: int) -> Parametrisation:
    """Extends the range of tau values used to up to two days, resolution varying across the allowed durations"""
    if index < 0 or index >= 26*17:
        raise InvalidIndexError(index)
    R = np.concatenate([
        np.arange(5,50,5),
        np.arange(50,100,10),
        np.arange(100,150,25),
        [150]
    ])

    tau = np.concatenate([
        np.arange(300,1800, 300),        # 5,10,...,30 minutes 
        np.arange(1800, 5400, 900),      # 30,45,...,90 minutes
        np.arange(5400, 10800, 1800),    # 90,120,150,180 minutes
        np.arange(10800, 21600, 3600),   # 3,4,5,6, hours
        np.arange(21600, 43200, 7200),   # 6,8,10,12 hours
        np.arange(43200, 64800, 10800),  # 12,15,18 hours
        np.arange(64800, 172800, 21600), # 18,24,30,...,48 hours
        [172800]
    ])
    R_, tau_ = np.meshgrid(R, tau)
    [R_km, tau_s] = np.array([
        R_.flatten(),
        tau_.flatten()
    ])[:, index]
    return Parametrisation(
        distance_km = float(R_km),
        tau = dt.timedelta(seconds=int(tau_s))
    )


@register_index_function
def R_500km_tau_172800s(index: int) -> Parametrisation:
    """Extends the range of R values checked to 500km and tau to two days"""
    R = np.concatenate([
        [5],
        np.arange(10, 100, 10),
        np.logspace(2, np.log10(500), 20)
    ])

    tau = np.concatenate([
        np.arange(300,1800, 300),        # 5,10,...,30 minutes 
        np.arange(1800, 5400, 900),      # 30,45,...,90 minutes
        np.arange(5400, 10800, 1800),    # 90,120,150,180 minutes
        np.arange(10800, 21600, 3600),   # 3,4,5,6, hours
        np.arange(21600, 43200, 7200),   # 6,8,10,12 hours
        np.arange(43200, 64800, 10800),  # 12,15,18 hours
        np.arange(64800, 172800, 21600), # 18,24,30,...,48 hours
        [172800]
    ])

    if index < 0 or index >= R.size * tau.size:
        raise InvalidIndexError(index)

    R_, tau_ = np.meshgrid(R, tau)
    [R_km, tau_s] = np.array([
        R_.flatten(),
        tau_.flatten()
    ])[:, index]
    return Parametrisation(
        distance_km = float(R_km),
        tau = dt.timedelta(seconds=int(tau_s))
    )
