"""Author: Andrew Martin
Creation date: 22/08/2025

Class definition for RawATL09 to handle data from .h5 files
"""

from sat_val_framework.implement import RawData

import xarray as xr 
import os
import datetime as dt


_Cloudnet_required_vars = (
    "longitude", "latitude",
    "classification_bits", "quality_bits",
)



def _load_cloudnet_from_file(fpath: str) -> xr.Dataset:
    """Function that, given a fully qualified file path to Cloudnet data, loads it and returns an xarray dataset"""
    ds = xr.open_dataset(fpath)
    # TODO: add quality and category bits to the dataset
    # TODO: ds = _include_cloudnet_cloudmask(ds)
    return ds



def _include_cloudnet_cloudmask(ds: xr.Dataset) -> xr.Dataset:
    data = ds.copy()
    data["cloudmask"] = (
        data.liq | (data.ice & data.frozen)
    ).rename("cloudmask")
    return data



class RawCloudnet(RawData):
    def assert_on_creation(self):
        # TODO: check dimensions in data
        # TODO: check required variables are in data
        for var in _Cloudnet_required_vars:
            assert var in self.data.variables, f"{var} not in Cloudnet dataset"

    @classmethod
    def from_qualified_file(cls, fpath: str) -> Self:
        data = _load_cloudnet_from_file(fpath=fpath)
        return cls(data = data, metadata = fpath)


    @classmethod
    def from_collocation_event_and_parameters(cls, event: CollocationEvent, parameters: CollocationParameters) -> Self:
        raise NotImplementedError(f"Type {type(self)} does not implement .from_collocation_event_and_parameters(cls, event: CollocationEvent, parameters: CollocationParameters)")

    def homogenise_to(self, H: Type[HomogenisedData]) -> H:
        raise NotImplementedError(f"Type {type(self)} does not implement .homngenise_to(self, H: Type[HomogenisedData])")
        
