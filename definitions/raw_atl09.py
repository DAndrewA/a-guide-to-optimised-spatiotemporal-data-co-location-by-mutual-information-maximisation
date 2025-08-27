"""Author: Andrew Martin
Creation date: 22/08/2025

Class definition for RawATL09 to handle data from .h5 files
"""

from __future__ import annotations

from .collocation import RadiusDuration
from sat_val_framework.implement import RawData

import xarray as xr 
import os
import datetime as dt


_ATL09_required_vars = (
    "latitude", "longitude",
    "layer_top", "layer_bot",
    "layer_attr", "layer_conf_dens",
    "surface_sig", "surface_thresh",
    "cloud_flag_atm", 
    # TODO: identify others
)


def _open_atl09_from_group(fpath: str) -> xr.Dataset:
    """Opens an ATL09 .h5 file and extracts the strong-beam profile groups as datasets"""
    # data = xr.merge(
    #     [
    #         xr.open_dataset(fpath, engine="h5cdf", group=f"profile_{p}/high_rate")
    #         for p in (1,2,3)
    #     ],
    #     dim = (("profile"), (1,2,3))
    # )
    return data
    
    
    
    
    raise NotImplementedError()






class RawATL09(RawData):
    def assert_on_creation(self):
        # TODO: check dimensions in data
        # TODO: check required variables are in data
        for var in _ATL09_required_vars:
            assert var in self.data.variables, f"{var} not in ATL09 dataset"
        

    @classmethod
    def from_qualified_file(cls, fpath: str) -> Self:
        data = _open_atl09_from_group(fpath)
        return cls(data=data, metadata=fpath)

    @classmethod
    def from_collocation_event_and_parameters(cls, event: CollocationEvent, parameters: CollocationParameters) -> Self:
        assert isinstance(parameters, RadiusDuration | None), f"collocation subsetting of {cls} not implemented for parameter type {type(parameters)}"

        atle = event.atl09_event
        fpaths = [
            f
            for f in [atle.fpath1, atle.fpath2]
            if f is not None
        ]
        #TODO: likely a combine="by_coords"-type jobby
        data = xr.merge(
            [
                _open_atl09_from_group(fpath)
                for fpath in fpaths
            ],
            dim = "time_index"
        )
        # apply subsetting based on parameters
        if parameters is None:
            pass
        if isinstance(parameters, RadiusDuration):
            pass # TODO: need to have a think about how collocation subsetting is handled between events and parameters



    def homogenise_to(self, H: Type[HomogenisedData]) -> H:
        raise NotImplementedError(f"Type {type(self)} does not implement .homngenise_to(self, H: Type[HomogenisedData])")
