"""Author: Andrew Martin
Creation date: 22/08/2025

Class definition for RawATL09 to handle data from .h5 files
"""

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

class RawATL09(RawData):
    def assert_on_creation(self):
        # TODO: check dimensions in data
        # TODO: check required variables are in data
        for var in _ATL09_required_vars:
            assert var in self.data.variables, f"{var} not in ATL09 dataset"
        
