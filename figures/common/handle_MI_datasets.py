"""Common functions for loading in the mutual information data"""

import os
import xarray as xr

HEIGHT_SLICE = slice(1000,10_000)
R_slice = slice(25,None,None)

RESULTS_DIR = "/work/scratch-nopw2/eeasm/MI/old1_K20_no_std"


def load_dataset_from_site(site: str) -> xr.Dataset:
    return (
        xr.open_dataset(
            os.path.join(
                RESULTS_DIR,
                f"MI_{site}.nc"
            )
        )
            .sel(height=HEIGHT_SLICE).mean(dim="height")
            .sel(R_km=R_slice)
    )


def MI_from_ds(ds: xr.Dataset) -> xr.DataArray:
    return ds.MI_holmes_total

def N_events_from_ds(ds: xr.Dataset) -> xr.DataArray:
    return ds.N_events

def N_profiles_from_ds(ds: xr.Dataset) -> xr.DataArray:
    return ds.N_profiles

def ds_at_max_value(ds: xr.Dataset, value: xr.DataArray) -> xr.Dataset:
    """return the dataset ds, subset to the argmax of the dataarray value"""
    argmax = value.argmax(...)
    return ds.isel(argmax)