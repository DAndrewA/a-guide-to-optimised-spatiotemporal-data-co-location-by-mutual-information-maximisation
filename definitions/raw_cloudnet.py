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


_Cloudnet_required_vars = (
    "longitude", "latitude",
    "classification_bits", "quality_bits",
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
    # TODO: add quality and category bits to the dataset
    # TODO: ds = _include_cloudnet_cloudmask(ds)
    return ds


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
    def from_qualified_file(cls, fpath: str) -> Self:
        data = _load_cloudnet_from_file(fpath=fpath)
        return cls(data = data, metadata = fpath)


    @classmethod
    def from_collocation_event_and_parameters(cls, event: CollocationEvent, parameters: CollocationParameters | None) -> Self:
        assert isinstance(parameters, None | RadiusDuration), f"{type(parameters)=} subsetting is not implemented for {cls}"
        cne = event.cloudnet_event
        t_0 = cne.closest_approach_time
        site = cne.site
        root_dir = cne.root_dir 

        # initial implementation is for parameters is None
        if parameters is None:
            fpath = os.path.join(
                root_dir,
                cls._fname_from_datetime(
                    datetime=t_0, 
                    site=site
                )
            )
            data = _load_cloudnet_from_file(fpath = fpath)

        elif isinstance(parameters, RadiusDuration):
            tau = parameters.duration
            time_slice = slice(t_0 - 0.5*tau, t_0 + 0.5*tau)
            fpaths = [
                os.path.join(
                    root_dir,
                    fname
                )
                for fname in [
                    cls._fname_from_datetime(
                        datetime = date,
                        site = site
                    )
                    for date in _dates_from_time_slice(time_slice)
                ]
            ]
            data = xr.merge(
                [
                    cls.from_qualified_file(fpath)
                    for fpath in fpaths
                ],
                #TODO: assess if this needs special handling to avoid errors with model_time, etc
            )

        return cls(data=data, metadata=event)
         


    def homogenise_to(self, H: Type[HomogenisedData]) -> H:
        raise NotImplementedError(f"Type {type(self)} does not implement .homngenise_to(self, H: Type[HomogenisedData])")
        

    @staticmethod
    def _fname_from_datetime(datetime: dt.datetime, site: str) -> str:
        """Create a valid Cloudnet classficiation filename from a datetime object and a site string.
        NOTE: this is not a fully qualified filepath. Prepend the directory for this to be the case.
        """
        # TODO: ensure that handling of "summit" is correct
        fname = f"{datetime:%Y%m%d}_{site}_classfification.nc"
        if site == "summit":
            fname = "{datetime:%Y%m%d}_classfification.nc"
        return fname
