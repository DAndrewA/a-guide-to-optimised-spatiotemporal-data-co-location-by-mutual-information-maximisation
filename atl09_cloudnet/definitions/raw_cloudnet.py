"""Author: Andrew Martin
Creation date: 22/08/2025

Class definition for RawCloudnet to handle raw data from .nc files.
"""

from __future__ import annotations

from sat_val_framework.implement import (
    RawData,
    RawMetadata,
    RawDataSubsetter,
    RawDataEvent,
    HomogenisedData,
)

from . import (
    cloudnet_decode, 
    vcf
)

from typing import ClassVar
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


def _safe_load_cloudnet_from_file(fpath: str) -> xr.Dataset | None:
    known_failure_types = (FileNotFoundError,)
    try:
        data = _load_cloudnet_from_file(fpath)
    except known_failure_types as ex:
        print("Known failure type:", type(ex), ex)
        data = None
    except Exception as e:
        raise e
    return data


#def _


def _include_cloudnet_cloudmask(ds: xr.Dataset) -> xr.Dataset:
    data = ds.copy()
    data["cloudmask"] = (
        data.liq | (data.falling & data.cold)
    ).rename("cloudmask")
    return data



class RawCloudnet(RawData):
    def assert_on_creation(self):
        # TODO: check dimensions in data
        # TODO: check required variables are in data
        for var in _Cloudnet_required_vars:
            assert var in self.data.variables, f"{var} not in Cloudnet dataset"

    @classmethod
    def from_qualified_file(cls, fpath: str) -> Self | None:
        data = _safe_load_cloudnet_from_file(fpath=fpath)
        if data is None:
            return None
        return cls(
            data = data, 
            metadata = RawMetadata(
                loader = fpath,
                subsetter = []
            )
        )


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
        data = _safe_load_cloudnet_from_file(fpath = fpath)
        if data is None:
            return None
        raw_data = cls(
            data = data,
            metadata = metadata
        )
        if parameters is None:
            return raw_data
        else:
            return parameters.subset(raw_data)

    def perform_qc(self) -> Self:
        #TODO: actually implement the qc procedure
        self.data["qc_cloudmask"] = (
            self.data.cloudmask
        )
        return self


    def _homogenise_to_VCF(self, H: Type[vcf.VCF]) -> vcf.VCF:
        vertical_cloud_fraction = (
            (self.data.qc_cloudmask)
                .mean(dim="time")
                .interp_like(H.CG.lin_interp_z)
                .fillna(0) # the interpolation can introduce NaNs at height < cloudnet.altitude, which mess up VCF assert_on_creation
                .rename("VCF")
        )
        return H(
            data = vertical_cloud_fraction,
            metadata = self.metadata
        )

    @property
    def n_profiles(self) -> int:
        return int(self.data.time.size)


    def homogenise_to(self, H: Type[HomogenisedData]) -> H:
        assert issubclass(H, HomogenisedData), f"{H} must be a subclass of {HomogenisedData}"
        if issubclass(H, vcf.VCF):
            return self._homogenise_to_VCF(H)
        raise TypeError(f"{type(self)} does not implement homogenise_to({H})")

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



@dataclass(kw_only=True, frozen=True)
class Duration(RawDataSubsetter):
    RDT: ClassVar = RawCloudnet
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
                data
                for fpath in fpaths
                if (data := _safe_load_cloudnet_from_file(fpath)) is not None
            ],
            dim="time"
        ).sel(time=tslice)
        return RawCloudnet(
            data = new_data,
            metadata = new_metadata
        )



@dataclass(frozen=True, kw_only=True)
class CloudnetEvent(RawDataEvent):
    """Class handling minimum required information to load file containing a collocation event.

    ATTRIBUTES:
        closest_approach_time (dt.datetime): datetime object describing the time of closest approach of ICESat-2 to the Cloudnet site. Used to center the time window that is loaded.
        root_dir (str): Path of the directory containing the .nc Cloudnet files to be loaded.
        site (str): Name of the site for which Cloudnet data is loaded. This is used in generating filenames to be loaded.
    """
    RDT: ClassVar = RawCloudnet
    closest_approach_time: dt.datetime
    root_dir: str
    site: str

