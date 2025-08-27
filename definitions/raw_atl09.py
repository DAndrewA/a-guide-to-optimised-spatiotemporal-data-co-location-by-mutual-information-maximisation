"""Author: Andrew Martin
Creation date: 22/08/2025

Class definition for RawATL09 to handle data from .h5 files
"""

from __future__ import annotations

#from collocation import RadiusDuration
#TODO: remove cyclic dependency between RadiusDuration and RawATL09
class RadiusDuration: pass # defined for type checking reasons
from sat_val_framework.implement import RawData

import xarray as xr 
import os
import datetime as dt


_ATL09_required_dims = (
    "time_index", "height", "ds_layers",
)
_ATL09_required_vars = (
    "time",
    "latitude", "longitude",
    "layer_top", "layer_bot", "layer_attr", "layer_con", "layer_conf_dens", "layer_flag",
    "cloud_flag_atm", 
    "surface_sig", "surface_thresh", "surface_height", 
    # TODO: identify others
)



def raw_groups_from_fpath(fpath: str) -> xr.Dataset:
    return [
        xr.open_dataset(
            fpath,
            engine = "h5netcdf", # TODO: include h5netcdf dependency in pyproject
            group = f"/profile_{p}/high_rate"
        )
        for p in (1,2,3,)
    ]

def safe_raw_groups_from_fpath(fpath: str) -> xr.Dataset | None:
    known_failure_types = ()
    try:
        data = raw_groups_from_fpath(fpath)
    except known_failure_types as ex:
        print("safe_raw_groups_from_fpath:", "Known error type:", type(ex), ex)
        data = None
    except Exception as e:
        raise e
    return data



def delta_time_as_variable(ds: xr.Dataset) -> xr.Dataset:
    return (
        ds
            .reset_index("delta_time")
            .reset_coords("delta_time")
            .rename_dims({"delta_time": "time_index"})
    )

def safe_delta_time_as_variable(ds: xr.Dataset) -> xr.Dataset | None:
    known_failure_types = (ValueError,)
    try:
        data = delta_time_as_variable(ds)
    except known_failure_types as ex:
        print("safe_delta_time_as_variable:", "Known error type:", type(ex), ex)
        data = None
    except Exception as e:
        raise e
    return data



def concat_profiles_with_padding(ds_objs: list[xr.Dataset]) -> xr.Dataset:
    longest_time_dim = max([
        ds.time_index.size
        for ds in ds_objs
    ])
    pad_sizes = [
        longest_time_dim - ds.time_index.size
        for ds in ds_objs
    ]
    return xr.concat(
        [
            ds.pad(
                pad_width = {"time_index": (0,pad_size)}
            )
            for ds, pad_size in zip(ds_objs, pad_sizes)
        ],
        dim = "profile"
    )

def safe_concat_profiles_with_padding(ds_objs: list[xr.Dataset]) -> xr.Dataset | None:
    known_failure_types = ()
    try:
        data = concat_profiles_with_padding(ds_objs)
    except known_failure_types as ex:
        print("safe_concat_profiles_with_padding:", "Known error type:", type(ex), ex)
        data = None
    except Exception as e:
        raise e
    return data



def set_extrapolated_time_as_coords(ds: xr.Dataset) -> xr.Dataset:
    return (
        ds
            .assign_coords(
                {"time_index": 
                    (ds.delta_time
                        .interpolate_na(
                            dim="time_index", 
                            method="linear", 
                            fill_value="extrapolate"
                        )
                        
                    )
                }
            )
    )

def rename_atl09_vars_dims(ds: xr.Dataset) -> xr.Dataset:
    return (
        ds
            .rename_dims({
                "ds_va_bin_h": "height",
            })
            .rename({
                "delta_time": "time",
                "ds_va_bin_h": "height"
            })
    )

def set_cloudmask(ds: xr.Dataset) -> xr.Dataset:
    data = ds.copy()
    data["cloudmask"] = (
        (
            (data["layer_bot"] <= data["height"]) &
            (data["layer_top"] >= data["height"])
        )
            .any(dim="ds_layers")
    )
    return data



def None_in(_iter):
    return any((v is None for v in _iter))

def safe_load_atl09_from_fpath(fpath: str) -> xr.Dataset | None:
    raw_groups = safe_raw_groups_from_fpath(fpath)
    if raw_groups is None:
        print("safe_load_atl09_from_fpath:", "None in raw_groups")
        return None
    groups = [
        safe_delta_time_as_variable(group)
        for group in raw_groups
    ]
    if None_in(groups):
        print("safe_load_atl09_from_fpath:", "None in groups")
        return None
    data = safe_concat_profiles_with_padding(groups)
    if data is None:
        print("safe_concat_profiles_with_padding:", "None returned after concatenation")
        return None
    # apply time interpolation and heighht renaming
    data = set_extrapolated_time_as_coords(data)
    data = rename_atl09_vars_dims(data)
    data = set_cloudmask(data)
    return data



class RawATL09(RawData):
    def assert_on_creation(self):
        # TODO: check dimensions in data
        # TODO: check required variables are in data
        for var in _ATL09_required_vars:
            assert var in self.data.variables, f"{var} not in ATL09 dataset"
        for dim in _ATL09_required_dims:
            assert dim in self.data.dims, f"{dim} not in ATL09 dataset dimensions"
        

    @classmethod
    def from_qualified_file(cls, fpath: str) -> Self | None:
        data = safe_load_atl09_from_fpath(fpath)
        if fpath is None:
            return None
        return cls(data=data, metadata=fpath)

    @classmethod
    def from_collocation_event_and_parameters(cls, event: CollocationEvent, parameters: CollocationParameters) -> Self | None:
        assert isinstance(parameters, RadiusDuration | None), f"collocation subsetting of {cls} not implemented for parameter type {type(parameters)}"
        # load the data safely from the 
        atle = event.event_atl09
        individual_data = [
            data
            for fpath in (atle.fpath1, atle.fpath2)
            if fpath is not None
            if (data := safe_load_atl09_from_fpath(fpath)) is not None
        ]
        if not individual_data:
            print(f"No data available for event={atle}")
            return None
        data = xr.combine_nested(
            individual_data,
            concat_dim = "time_index"
        )
        
        # apply subsetting based on parameters
        if parameters is None:
            return cls(data=data, metadata=event)

        if isinstance(parameters, RadiusDuration):
            raise NotImplementedError() # TODO: need to have a think about how collocation subsetting is handled between events and parameters



    def homogenise_to(self, H: Type[HomogenisedData]) -> H:
        raise NotImplementedError(f"Type {type(self)} does not implement .homngenise_to(self, H: Type[HomogenisedData])")
