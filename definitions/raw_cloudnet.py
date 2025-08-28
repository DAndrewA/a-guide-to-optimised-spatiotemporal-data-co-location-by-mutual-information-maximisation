"""Author: Andrew Martin
Creation date: 22/08/2025

Class definition for RawCloudnet to handle raw data from .nc files.
"""

from __future__ import annotations

from sat_val_framework.implement import (
    RawData,
    RawMetadata,
    RawDataSubsetter,
    CollocationEvent,
)

import cloudnet_decode

from dataclasses import dataclass
from pandas import Timestamp
import xarray as xr 
import numpy as np
import os
import datetime as dt


_Cloudnet_required_vars = (
    "longitude", "latitude",
)


def _dates_from_time_slice(tslice: slice) -> list[dt.date]:
    current_date = tslice.start.date()
    end = tslice.stop
    dates = list()
    one_day = dt.timedelta(days=1)
    while current_date <= end.date():
        dates.append(current_date)
        current_date += one_day
    return dates



def _load_cloudnet_from_file(fpath: str) -> xr.Dataset:
    """Function that, given a fully qualified file path to Cloudnet data, loads it and returns an xarray dataset"""
    ds = xr.open_dataset(fpath)
    ds = cloudnet_decode.decode_cat_ds(ds)
    ds = _include_cloudnet_cloudmask(ds)
    return ds


#def _


def _include_cloudnet_cloudmask(ds: xr.Dataset) -> xr.Dataset:
    data = ds.copy()
    data["cloudmask"] = (
        data.liq | (data.falling & data.cold)
    ).rename("cloudmask")
    return data



@dataclass(kw_only=True, frozen=True)
class Duration(RawDataSubsetter):
    duration: dt.timedelta

    def subset(self, raw_data: RawCloudnet) -> RawCloudnet:
        new_data = raw_data.data.copy()
        new_metadata = raw_data.metadata
        # if no closest approach time is known, cannot subset by window centered on it
        if not isinstance(new_metadata.loader, CloudnetEvent):
            print(new_metadata, self, f"No closest approach time known, so no subsetting")
            return raw_data

        new_metadata.subsetter.append(self)

        center = new_metadata.loader.closest_approach_time
        tslice = slice(
            center - 0.5*self.duration,
            center + 0.5*self.duration
        )
        # see if more data needs loading or if current spans sufficient time
        min_time = Timestamp(new_data.time.min().values).to_pydatetime()
        max_time = Timestamp(new_data.time.max().values).to_pydatetime()
        tslice_in_time_range = (
            (tslice.start >= min_time) &
            (tslice.stop <= max_time)
        )
        if tslice_in_time_range:
            new_data = new_data.sel(time=tslice)
            return RawCloudnet(
                data = new_data,
                metadata = new_metadata
            )
        # if full time slice not already contained, load more data!
        dates_to_load = _dates_from_time_slice(tslice)
        fpaths = [
            os.path.join(
                new_metadata.loader.root_dir,
                RawCloudnet._fname_from_datetime(
                    datetime = date,
                    site = new_metadata.loader.site
                )
            )
            for date in dates_to_load
        ]
        new_data = xr.concat(
            [
                _load_cloudnet_from_file(fpath)
                for fpath in fpaths
            ],
            dim="time"
        ).sel(time=tslice)
        return RawCloudnet(
            data = new_data,
            metadata = new_metadata
        )



@dataclass(frozen=True, kw_only=True)
class CloudnetEvent(CollocationEvent):
    """Class handling minimum required information to load file containing a collocation event.

    ATTRIBUTES:
        closest_approach_time (dt.datetime): datetime object describing the time of closest approach of ICESat-2 to the Cloudnet site. Used to center the time window that is loaded.
        root_dir (str): Path of the directory containing the .nc Cloudnet files to be loaded.
        site (str): Name of the site for which Cloudnet data is loaded. This is used in generating filenames to be loaded.
    """
    closest_approach_time: dt.datetime
    root_dir: str
    site: str



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
    def from_collocation_event_and_parameters(cls, event: CloudnetEvent, parameters: Duration | None) -> Self:
        assert isinstance(event, CloudnetEvent), f"{type(event)=} must be CloudnetEvent"
        #assert isinstance(parameters, None | RadiusDuration), f"{type(parameters)=} subsetting is not implemented for {cls}"
        metadata = RawMetadata(
            loader = event,
            subsetter = list()
        )

        t_0 = event.closest_approach_time
        site = event.site
        root_dir = event.root_dir 

        # initially load single day of data, then apply subsetter if present
        fpath = os.path.join(
            root_dir,
            cls._fname_from_datetime(
                datetime=t_0, 
                site=site
            )
        )
        data = _load_cloudnet_from_file(fpath = fpath)
        raw_data = cls(
            data = data,
            metadata = metadata
        )
        if parameters is None:
            return raw_data
        else:
            return parameters.subset(raw_data)

    def perform_qc(self) -> Self:
        #TODO: implement
        return self

    def homogenise_to(self, H: Type[HomogenisedData]) -> H:
        raise NotImplementedError(f"Type {type(self)} does not implement .homngenise_to(self, H: Type[HomogenisedData])")
        

    @staticmethod
    def _fname_from_datetime(datetime: dt.datetime, site: str) -> str:
        """Create a valid Cloudnet classficiation filename from a datetime object and a site string.
        NOTE: this is not a fully qualified filepath. Prepend the directory for this to be the case.
        """
        # TODO: ensure that handling of "summit" is correct
        fname = f"{datetime:%Y%m%d}_{site}_categorize.nc"
        if site == "summit":
            fname = "{datetime:%Y%m%d}_categorize.nc"
        return fname
