import ctypes
import numpy as np

LIB = ctypes.CDLL("./libMIxnyn.so")

LIB.MIxnyn.argtypes = [
    ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),  # double **x
    ctypes.c_int,                                      # int dimx
    ctypes.c_int,                                      # int dimy
    ctypes.c_int,                                      # int K
    ctypes.c_int,                                      # int N
    ctypes.POINTER(ctypes.c_double)                   # double *MI
]
LIB.MIxnyn.restype = None


def call_MI_xnyn(X: np.ndarray,Y:np.ndarray,K:int) -> np.ndarray:
    dimx, Nx = [int(v) for v in X.shape]
    dimy, Ny = [int(v) for v in Y.shape]

    assert Nx == Ny
    N = Nx

    assert isinstance(K,int)
    assert isinstance(N, int)
    assert isinstance(dimx, int)
    assert isinstance(dimy, int)

    input_x = np.ascontiguousarray(
        np.concatenate((X,Y),axis=0),
        dtype=np.float64
    )

    MI = ctypes.c_double(0.0)

    # Convert x to array of pointers (double **)
    x_ptrs = (ctypes.POINTER(ctypes.c_double) * input_x.shape[0])()
    for i in range(input_x.shape[0]):
        x_ptrs[i] = input_x[i].ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    print(dimx, dimy, K, N)
    print(x_ptrs)
    print(input_x.shape)
    print(MI)
    # Call the function
    LIB.MIxnyn(
        x_ptrs,
        dimx,
        dimy, 
        K,
        N,
        ctypes.byref(MI)
    )
    print(f"Result: {MI=}")
    print(f"{MI.value=}")
    return MI.value

