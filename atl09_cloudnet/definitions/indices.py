"""Author: Andrew Martin
Creation date: 3/7/25

Script handling conversion of SLURM array indices to pairs of (R, tau)
"""

from dataclasses import dataclass
import numpy as np
import datetime as dt

INDEX_FUNCTIONS = dict()


@dataclass(frozen=True, kw_only=True)
class Parametrisation:
    distance_km: float
    tau: dt.timedelta


class InvalidIndexError(ValueError):
    pass


def register_index_function_product(name: str, R: np.ndarray, tau: np.ndarray):
    MAX_INDEX = R.size * tau.size
    def index_function(index: int) -> Parametrisation:
        if index < 0 or index >= MAX_INDEX:
            raise InvalidIndexError(index)
        
        Ri, taui = divmod(index, tau.size)
        
        R_km = float(R[Ri])
        tau_s = int(tau[taui])

        return Parametrisation(
            distance_km = R_km,
            tau = dt.timedelta(seconds=tau_s)
        )
    index_function.MAX_INDEX = MAX_INDEX
    index_function.__name__ = name
    INDEX_FUNCTIONS[name] = index_function


def register_index_function(name: str, R: np.ndarray, tau: np.ndarray):
    assert np.asarray(R).shape == np.asarray(tau).shape
    MAX_index = np.asarray(R).size
    def index_function(index: int) -> Parametrisation:
        if index < 0 or index >= MAX_index:
            raise InvalidIndexError(index)

        R_km = float(R[index])
        tau_s = int(tau[index])

        return Parametrisation(
            distance_km = R_km,
            tau = dt.timedelta(seconds=tau_s)
        )
    index_function.MAX_INDEX = MAX_INDEX
    index_function.__name__ = name
    INDEX_FUNCTIONS[name] = index_function

    

register_index_function_product(
    name = "original",
    R = np.arange(5,151,5),
    tau = np.arange(300,4501,300),
)


register_index_function_product(
    name = "extended_tau_meshgrid",
    R = np.concatenate([
        np.arange(5,50,5),
        np.arange(50,100,10),
        np.arange(100,150,25),
        [150]
    ]),
    tau = np.concatenate([
        [300],
        np.arange(600,4800, 600),
        np.arange(4800, 21600, 1200),
        np.arange(21600, 43200, 3600),
        [43200]
    ]),
)


register_index_function_product(
    name = "R_150km_tau_172800s",
    R = np.concatenate([
        np.arange(5,50,5),
        np.arange(50,100,10),
        np.arange(100,150,25),
        [150]
    ]),
    tau = np.concatenate([
        np.arange(300,1800, 300),        # 5,10,...,30 minutes 
        np.arange(1800, 5400, 900),      # 30,45,...,90 minutes
        np.arange(5400, 10800, 1800),    # 90,120,150,180 minutes
        np.arange(10800, 21600, 3600),   # 3,4,5,6, hours
        np.arange(21600, 43200, 7200),   # 6,8,10,12 hours
        np.arange(43200, 64800, 10800),  # 12,15,18 hours
        np.arange(64800, 172800, 21600), # 18,24,30,...,48 hours
        [172800]
    ]),
)


register_index_function_product(
    name="R_500km_tau_172800s",
    R = np.concatenate([
        [5],
        np.arange(10, 100, 10),
        np.logspace(2, np.log10(500), 20)
    ]),
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
)

R_lit_opt = [50,50,500,500,100]
tau_lit_opt = [1800,172800,1800,172800,3600]

register_index_function(
    "R_tau_extremal_lit_opt_ny-alesund",
    R = R_lit_opt + [60]
    tau = tau_lit_opt + [21600],
)

register_index_function(
    "R_tau_extremal_lit_opt_hyytiala",
    R = R_lit_opt + [180.9]
    tau = tau_lit_opt + [28800],
)

register_index_function(
    "R_tau_extremal_lit_opt_juelich",
    R = R_lit_opt + [214.3]
    tau = tau_lit_opt + [43200],
)

register_index_function(
    "R_tau_extremal_lit_opt_munich",
    R = R_lit_opt + [140.3]
    tau = tau_lit_opt + [14400],
)
