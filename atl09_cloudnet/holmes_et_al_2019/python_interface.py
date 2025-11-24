import ctypes
import numpy as np
import os
from dataclasses import dataclass
from typing import Self

fpath_so = os.path.join(
    os.path.dirname(__file__),
    "libMIxnyn.so"
)
LIB = ctypes.CDLL(fpath_so)

LIB.MIxnyn.argtypes = [
    ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),  # double **x
    ctypes.c_int,                                      # int dimx
    ctypes.c_int,                                      # int dimy
    ctypes.c_int,                                      # int K
    ctypes.c_int,                                      # int N
    ctypes.POINTER(ctypes.c_double)                   # double *MI
]
LIB.MIxnyn.restype = None



@dataclass(frozen=True, kw_only=True)
class MIEstimate:
    MI: float
    std: float
    N: int
    n_splits: int
    sigma_i: list[int]
    M: int
    K: int

    @property
    def ddof(self): 
        return self.M

    @classmethod
    def from_XYKMn_with_RNG(cls: Self, X: np.ndarray, Y: np.ndarray, K: int, M: int, n_splits: int, n_samples: int, RNG: np.random.Generator) -> Self:
        """Function that given input data X and Y, parameter K, degrees of freedom M, and number of splits per degree of freedom n_splits, will compute an MIEstimate instance.
        Also checks that the expected number of samples are provided (guards against transposed data input)"""
    # Start with assertions on the validity of the input data
    dimx, Nx = [int(v) for v in X.shape]
    dimy, Ny = [int(v) for v in Y.shape]

    assert Nx == Ny, f"Number of samples in X and Y are unqual: ({Nx=}) != ({Ny=})"
    assert Nx == n_samples, f"Number of samples given is not the expected number of samples, {Nx} != ({n_samples=})"

    assert isinstance(M, int), f"{M=} must be an integer"
    assert isinstance(n_splits, int), f"{n_splits=} must be an integer"
    assert n_splits > 3, f"{n_splits=} must be in Z/{1,2,3}, integer greater than 3"
    assert int( n_samples // (n_splits+1) ) > 3, f"Too few {n_samples=} for given ddof {M=} to compute std of MI estimator."  

    assert isinstance(K,int), f"{K=} must be integer"
    assert isinstance(N, int), f"{N=} must be integer"
    assert isinstance(dimx, int)
    assert isinstance(dimy, int)

    assert isinstance(RNG, np.random.Generator), f"{type(RNG)=} should be instance of {np.random.Generator}"

    XY = np.concatenate((X,Y),axis=0),

    MI =  _unsafe_call_MI_xnyn(
        X=X,
        Y=Y,
        K=K
    )

    # generate n_splits non-overlapping subsets of the input data, M times, and use to compute sigma_i
    sigma_i = list()
    n_i = list()
    for _ in range(M):
        MI_for_current_iteration = list()
        for XY_shuffled in np.array_split( 
                ary = RNG.permutation(XY, axis=1), 
                indices_or_sections = n_splits,
                axis = 1
            ):
            MI_for_current_iteration.append(
                _unsafe_call_MI_xnyn_from_XY(
                    XY = XY_shuffled,
                    dimx = dimx,
                    dimy = dimy,
                    K = K,
                )
            )
        sigma_i.append(np.std(MI_for_current_iteration))
        n_i.append(n_splits)

    sigma_i = np.asarray(sigma_i)
    n_i = np.asarray(n_i)
    # computation is Eq. (8)/N from (10.1103/PhysRevE.100.022404). Communicationas with Holmes shows that this is the correct formulation of the formula, as 1. the chi2 distribution has an additional factor of 1/2 in the exponentiation, and 2. x~sigma_i is poorly defined, so a probability density based on a value a_i sigma_i^2 / B needs to be used, making use of the Jacobian |d(a_i sigma_i^2/B)/d(sigma_i^2)|
    a_i_div_N = (n_i - 1) / n_i
    k_i = n_i - 1
    std_MI = np.sqrt(
        np.sum(
            a_i_div_N * np.power(sigma_i, 2)
        )
        / np.sum(k_i)
    )

    return MIEstimate(
        MI = MI,
        std = std_MI,
        N = n_samples,
        n_splits = n_splits,
        sigma_i = list(sigma_i),
        M = M,
        K = K,
    )
            


def call_MI_xnyn(X: np.ndarray, Y:np.ndarray, K:int, n_samples: int) -> float:
    """A safe function to call LIB.MI_xnyn, that ensures certain assertions on the data are met before the computation can proceed"""
    dimx, Nx = [int(v) for v in X.shape]
    dimy, Ny = [int(v) for v in Y.shape]

    assert Nx == Ny, f"Number of samples in X and Y are unqual: ({Nx=}) != ({Ny=})"
    assert Nx == n_samples, f"Number of samples given is not the expected number of samples, {Nx} != ({n_samples=})"

    assert isinstance(K,int), f"{K=} must be integer"
    assert isinstance(N, int), f"{N=} must be integer"
    assert isinstance(dimx, int)
    assert isinstance(dimy, int)

    return _unsafe_call_MI_xnyn(
        X=X,
        Y=Y,
        K=K
    )

def _unsafe_call_MI_xnyn(X: np.ndarray, Y: np.ndarray, K: int) -> float:
    """Function that calls LIB.MIxnyn without any assertions on the input data"""
    dimx, N = [int(v) for v in X.shape]
    dimy, _ = [int(v) for v in Y.shape]

    input_x = np.ascontiguousarray(
        np.concatenate((X,Y),axis=0),
        dtype=np.float64
    )

    MI = ctypes.c_double(0.0)

    x_ptrs = (ctypes.POINTER(ctypes.c_double) * input_x.shape[0])()
    for i in range(input_x.shape[0]):
        x_ptrs[i] = input_x[i].ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    LIB.MIxnyn(
        x_ptrs,
        dimx,
        dimy, 
        K,
        N,
        ctypes.byref(MI)
    )
    return MI.value

def _unsafe_call_MI_xnyn_from_XY(XY: np.ndarray, dimx: int, dimy: int, K: int) -> float:
    """Function that calls LIB.MIxnyn without assertions on the input data, where X and Y have already been correctly concatenated"""
    input_x = np.ascontiguousarray(
        XY,
        dtype=np.float64
    )
    ndims, N = input_x.shape 

    MI = ctypes.c_double(0.0)

    x_ptrs = (ctypes.POINTER(ctypes.c_double) * ndims)()
    for i in range(ndims):
        x_ptrs[i] = input_x[i].ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    LIB.MIxnyn(
        x_ptrs,
        dimx,
        dimy, 
        K,
        N,
        ctypes.byref(MI)
    )
    return MI.value
