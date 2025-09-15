"""Author: Andrew Martin
Creation date: 22/08/2025

Class definition for RawATL09 to handle data from .h5 files
"""

from __future__ import annotations

from sat_val_framework.implement import (
    RawDataSubsetter,
    RawDataEvent, 
    RawData,
    RawMetadata,
    HomogenisedData
)

from . import vcf

from typing import ClassVar
from dataclasses import dataclass
from pandas import Timestamp
import xarray as xr 
import numpy as np
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
CONF_DENS_THRESHOLD = 0.4


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
                {"time": 
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
                #"delta_time": "time",
                "ds_va_bin_h": "height"
            })
    )

def set_feature_mask(ds: xr.Dataset) -> xr.Dataset:
    data = ds.copy()
    clouds = (
        (data["layer_bot"] <= data["height"]) &
        (data["layer_top"] >= data["height"]) &
        ( data["layer_attr"].isin([1,]) ) # 1: cloud
    ).any(dim="ds_layers")
    aerosols = (
        (data["layer_bot"] <= data["height"]) &
        (data["layer_top"] >= data["height"]) &
        ( data["layer_attr"].isin([2,]) ) # 2: aerosol
    ).any(dim="ds_layers")
    attenuated = (
        (data["surface_sig"] < data["surface_thresh"]) &
        (data["height"] < data["layer_bot"].min(dim="ds_layers"))
    )
    
    feature_mask = xr.zeros_like(clouds, dtype=np.uint8)
    #NOTE: uint8 only supports 8 unique features
    features = (clouds, aerosols, attenuated)
    for i, feature in enumerate(features):
        feature_mask += (feature * 1<<i).astype(np.uint8)

    feature_mask = feature_mask.assign_attrs(
        description = """Feature mask from the ATL09 data. Bits are set to indicate the presence of a certain feature in the mask.
1<<0: Cloud
1<<1: Aerosol
1<<2: Attenuation"""
    )

    data["feature_mask"] = feature_mask
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
    data = set_feature_mask(data)
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
        if data is None:
            return None
        return cls(
            data=data, 
            metadata=RawMetadata(
                loader = fpath,
                subsetter = []
            )
        )

    @classmethod
    def from_collocation_event_and_parameters(cls, event: ATL09Event, parameters: RawDataSubsetter) -> Self | None:
        assert isinstance(event, ATL09Event), f"{type(event)=} must be ATL09 event"
        assert isinstance(parameters, ATL09_COLLOCATION_SUBSET_TYPES | None), f"collocation subsetting of {cls} not implemented for parameter type {type(parameters)}"
        
        metadata = RawMetadata(
            loader = event,
            subsetter = list()
        )

        # load the data safely from the 
        individual_data = [
            data
            for fpath in (event.fpath1, event.fpath2)
            if fpath is not None
            if (data := safe_load_atl09_from_fpath(fpath)) is not None
        ]
        if not individual_data:
            print(f"No data available for event={event}")
            return None
        data = xr.combine_nested(
            individual_data,
            concat_dim = "time_index"
        )
        
        # apply subsetting based on parameters
        raw_data_ob = cls(data=data, metadata=metadata)
        if parameters is None:
            return raw_data_ob
        return parameters.subset(raw_data_ob)

    @property
    def n_profiles(self) -> int:
        invalid_profiles_mask = self.data.feature_mask.isnull().all(dim="height")
        return int(
            invalid_profiles_mask.size - invalid_profiles_mask.sum()
        )


    def perform_qc(self) -> Self:
        # in this instance, qc is performed on the featuremask by ensuring that 
        #raise NotImplementedError()
        ds = self.data
        ds["qc_feature_mask"] = xr.where(
            (
                (ds["layer_bot"] <= ds["height"]) &
                (ds["layer_top"] >= ds["height"]) &
                (ds["layer_conf_dens"] <= CONF_DENS_THRESHOLD)
            ).any(dim="ds_layers"),
            0,
            ds.feature_mask,
        )
        return self


    def _homogenise_to_VCF(self, H: Type[vcf.VCF]) -> vcf.VCF:
        vertical_cloud_fraction = (
            (self.data.feature_mask == 1)
                .mean(dim=["time_index","profile"])
                .interp_like(H.CG.lin_interp_z)
                .rename("VCF")
        )
        return H(
            data = vertical_cloud_fraction,
            metadata = self.metadata
        )

    def homogenise_to(self, H: Type[HomogenisedData]) -> H:
        assert issubclass(H, HomogenisedData), f"{H} must be a subclass of {HomogenisedData}"
        if issubclass(H, vcf.VCF):
            return self._homogenise_to_VCF(H)
        raise TypeError(f"{type(self)} does not implement homogenise_to({H})")




@dataclass(kw_only=True, frozen=True)
class DistanceFromLocation(RawDataSubsetter):
    RDT: ClassVar = RawATL09
    distance_km: float
    latitude: float
    longitude: float

    # TODO: decide on a qc threshold
    MINIMUM_REQUIRED_PROFILES: int = 17 # floor (5 km of required data) / (0.28 km per profile)

    def _haversine_distance(self, lat_sat, lon_sat, a=6378):
        """Function to calculate the haversine distance between pairs of latitude and longitude coordinates and a fixed latitude and longitude location (lat0, lon0)
        NOTE: angles should be provided in decimal degrees
        NOTE: the value of R defines the units that the result is given in (in units per radian). The default is km along the surface of the Earth.
        """
        lat_satr, lon_satr = np.deg2rad(lat_sat), np.deg2rad(lon_sat)
        lat0r, lon0r = np.deg2rad(self.latitude), np.deg2rad(self.longitude)
        alphar = np.arccos(
            np.sin(lat0r) * np.sin(lat_satr) +
            np.cos(lat0r) * np.cos(lat_satr) * np.cos(lon0r - lon_satr)
        )
        s = a * alphar
        return s


    def get_distance_to_location(self, raw_data: RawATL09) -> xr.DataArray:
        stack = {"time_index_profile": ("time_index", "profile")}
        new_data = raw_data.data.copy().stack(stack)
        d2s = self._haversine_distance(
            lat_sat = new_data["latitude"],
            lon_sat = new_data["longitude"],
        )
        return d2s.rename("distance_to_site")
        
    def get_time_closest_approach(self, raw_data: RawATL09) -> dt.datetime:
        d2s = self.get_distance_to_location(raw_data)
        t0_dt64 = d2s.time.isel(
            time_index_profile = d2s.argmin()
        ).values
        return Timestamp(t0_dt64).to_pydatetime()

    def subset(self, raw_data: RawATL09) -> RawATL09:
        stack = {"time_index_profile": ("time_index", "profile")}
        new_data = raw_data.data.copy().stack(stack)
        d2s = self.get_distance_to_location(raw_data)
        #print("d2s range:", d2s.min().values, d2s.max().values)
        valid_subset = d2s <= self.distance_km
        if sum(valid_subset) < self.MINIMUM_REQUIRED_PROFILES:
            #TODO: implement a subsetting error
            print("Whelp, this is boring")
            #raise SubsetError(f"Insufficient profiles for {raw_data.metadata.loader} when subsetting with {self}")

        new_data = new_data.where(valid_subset, drop=True).unstack()
        # reset the time coordinates, even if data for a given column is null
        new_data["time"] = new_data.time.interpolate_na(
            dim="time_index",
            method="linear",
            fill_value="extrapolate",
        )
        # update metadata to indicate subsetting
        new_metadata = raw_data.metadata
        new_metadata.subsetter.append(self)

        return RawATL09(
            data = new_data,
            metadata = new_metadata
        )
        


@dataclass(frozen=True, kw_only=True)
class ATL09Event(RawDataEvent):
    """Class handling minimum required information to load a file(s) containing a collcoation event

    ATTRIBUTES:
        fpath1 (str): fully qualified path to a .h5 file containing ATL09 data to be loaded.
        min_separation_km (float | None) the minimum separation between the ATL09 data and Cloudnet site during the event. (optional)
    """
    RDT: ClassVar = RawATL09
    fpath: str
    min_separation_km: float | None



#NOTE: at bottom of file so that RawATL09 can be defined before DistanceFromLocation and ATL09Event require it in their definitions, but so that DistanceFromLocation can be defined before ATL09_COLLOCATION_SUBSET_TYPES
# include multiple implemented collocation types as a union type, rather than a tuple of types
ATL09_COLLOCATION_SUBSET_TYPES = (DistanceFromLocation | DistanceFromLocation)

