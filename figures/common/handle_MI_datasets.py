"""Common functions for loading in the mutual information data"""

from .common import DIR_ROOT

import os
import xarray as xr

HEIGHT_SLICE = slice(1000,10_000)
R_slice = slice(25,None,None)

FPATH_MI = os.path.join(
    DIR_ROOT, "MI", "MI_merged.nc"
)


K = 10
R_slice = slice(15,None,None)

def get_MI_with_subsetting(**subset) -> xr.Dataset:
    """Returns the MI_merged dataset, subset by any additional specified parameters.
    NOTE: The dataset is pre-subset by K=10 and R=R_slice.
    """
    return (
        xr.open_dataset(FPATH_MI)
            .drop_dims("height")
            .sel(dict(K=K, R_km=R_slice))
            .sel(subset)
    )

