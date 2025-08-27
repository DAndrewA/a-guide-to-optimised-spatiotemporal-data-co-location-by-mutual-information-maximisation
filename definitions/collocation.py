"""Author: Andrew Martin
Creation date: 25/08/2025

Class definitions for RadiusDuration, ATL09Event, CloudnetEvent, CloudnetATL09Event and SchemeCloudnetATL09RadiusDuration.
"""


from sat_val_framework.implement import (
    CollocationParameters, 
    CollocationEvent, 
    CollocationScheme,
)
from sat_val_framework import CollocationEventList
from .raw_cloudnet import RawCloudnet
from .raw_atl09 import RawATL09

import datetime as dt
import os



def _one_of_each_raw_type(raw_data1, raw_data2) -> bool:
    """Returns True if both arguments span the required_types"""
    required_types = [RawATL09, RawCloudnet]
    first_type_is = [type(raw_data1) == t for t in required_types]
    second_type_is = [type(second) == t for t in required_types]
    both_types_included = all([ 
        bool1 & bool2  
        for bool1, bool2 in zip(first_type_is, second_type_is)
    ])
    return both_types_included

class RadiusDuration(CollocationParameters):
    radius_km: float
    tau: dt.timedelta
    longitude: float
    latitude: float

    def apply_collocation_subsetting(self, raw_data: RawData) -> RawData:
        if isinstance(raw_data, RawCloudnet):
            # TODO: subset based on temporal criteria
            pass
        elif isinstance(raw_data, RawATL09):
            # TODO: subset based on spatial criteria
            pass
        else:
            raise TypeError(f"{type(raw_data)} is not of type RawATL09 or RawCloudnet")
        return raw_data

    @staticmethod
    def calculate_collocation_criteria(raw_data1: RawData, raw_data2: RawData) -> tuple[RawData, RawData]:
        raise NotImplementedError(f"Type {type(self)} does not implement .calculate_collocation_criteria(self, raw_data1: RawData, raw_data2: RawData)")




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

class ATL09Event(CollocationEvent):
    """Class handling minimum required information to load a file(s) containing a collcoation event

    ATTRIBUTES:
        fpath1 (str): fully qualified path to a .h5 file containing ATL09 data to be loaded.
        fpath2 (str | None) fully qualified path to a .h5 file containing ATL09 data to be loaded. Use in the event that a collocation event spans a granule boundary and is therefore split across files.
        latitude (float): Latitude in decimal degrees of the associated Cloudnet site.
        longitude (float): Longitude in decimal degrees of the associated Cloudnet site.
    """
    fpath1: str
    fpath2: str | None
    latitude: float | None
    longitude: float | None

class CloudnetATL09Event(CollocationEvent):
    """Class handling a CollocationEvent between Cloudnet and ATL09 data.

    ATTRIBUTES:
        event_atl09 (ATL09Event): The information required to load ATL09 data associated with the given collocation event.
        event_cloudnet (CloudnetEvent): The information required to load Cloudnet data associated with the given collocation event.
    """
    event_atl09: ATL09Event
    event_cloudnet: CloudnetEvent



def _are_atl09_orbits_subsequent(fname1: str, fname2: str) -> bool:
    """Returns True if the filenames associated with ATL09 data are from subsequent orbits"""
    __, __, __, rgt_cycle1, __, __ = fname1.split("_")
    __, __, __, rgt_cycle2, __, __ = fname2.split("_")
    rgt1, cycle1 = int(rgt_cycle1[:4]), int(rgt_cycle1[4:6])
    rgt2, cycle2 = int(rgt_cycle2[:4]), int(rgt_cycle2[4:6])


    processed_ATL09_20190126121018_04440201_006_02.h5

    condition1 = (cycle1 == cycle2) and (rgt1 + 1 == rgt2)
    condition2 = (cycle1 + 1 == cycle2) and (rgt1 == 1387) and (rgt2 == 1)
    return condition1 or condition2



class SchemeCloudnetATL09RadiusDuration(CollocationScheme):
    """Collocation scheme to find collocation events between Cloudnet and ATL09 data, based on a radius-duration scheme.
    The radius-duration scheme is described by selecting ATL09 data that falls within a given radius R of the Cloudnet site, and selecting Cloudnet data within a temporal window of duration tau, centered on the time of closest approach.

    STATIC METHODS:
        get_matches_from_raw_directories(file_list_atl09: list[str], file_list_cloudnet: list[str]) -> CollocationEventList: from a provided list of ATL09 files and Cloudnet files, identifies all available matches and returns a list of CloudnetATL09Event instances.
    """

    @staticmethod
    def get_matches_from_raw_directories(file_list_atl09: list[str], dir_cloudnet: str, cloudnet_site: str) -> CollocationEventList:
        """From a provided list of ATL09 files and Cloudnet directory and sitename, identifies all available matches and returns a list of CloudnetATL09Event instances. 
        
        The premise for determining if an ATL09 collocation event is spread across multiple files is that, typically, for a randomly placed event, the granule will not be split and therefore all data will be obtained from a single file.
        However, in the event that a Cloudnet site is close to a granule-change latitude, then data will be split across files with subsequent orbit numbers.
        In the event that there are two collocation events on subsequent orbits, both consisting of data across two files, the harmony-subset data describes a time when the file is valid for.
        The two files belonging to the same event will have the smallest absolute timedelta between them, allowing us to identify events in long strings of subsequent orbit-number files.
        """
        skip_this_one = False
        event_list = list()
        for (fname_atl1, fname_atl2) in zip(file_list_atl09, file_list_atl09[1:] + [None]):
            if skip_this_one or fname_atl1 is None:
                # skip this one will be set if subsequent files are used for a collocation event
                skip_this_one = False
                continue
            
            atl09_event_args = {
                "fpath1": fname_atl1, 
                "fpath2": None,
                "latitude": None,
                "longitude": None,
            }
            if fname_atl2 is not None:
                if _are_atl09_orbits_subsequent(fname_atl1, fname_atl2):
                    skip_this_one = True
                    atl09_event_args["fpath2"] = fname_atl2

            atl09_event = ATL09Event(
                **atl09_event_args 
            )

            # use the given event to load the ATL09 data
            raw_atl09 = RawATL09.from_collocation_event_and_parameters(
                event = atl09_event, 
                parameters = None
            )
            
            # obtain a crude estimate of t_0 from the median time in the ATL09 data
            crude_t0 = raw_atl09.data.time.median()
            cloudnet_fname = RawCloudnet._fname_from_datetime(datetime = crude_t0, site=cloudnet_site)
            # check file exists during collocation event
            if not os.path.isfile( os.path.join(dir_cloudnet, cloudnet_fname) ):
                print(f"{cloudnet_fname} not found for ATL09 event")
                continue
            cloudnet_event_args = {
                "root_dir": dir_cloudnet,
                "site": cloudnet_site,
                "closest_approach_time": crude_t0
            }

            # load the Cloudnet data using the crude t_0 value.
            # NOTE: this only produces minor errors for collocation events +-5 minutes around midnight (i.e. 0.7% of cases assuming random uniform distribution of local overpass times.)
            crude_cloudnet_event = CloudnetEvent(
                **cloudnet_event_args
            )
            crude_raw_cloudnet = RawCloudnet.from_collocation_event_and_parameters(
                event = cloudnet_event,
                parameters = None 
            )

            # RadiusDuration implements calculating the collocation criteria from RawData instances
            # TODO: consolidate call signature for calculate_collocation_criteria, either change calling here, or signature definition
            raw_atl09, crude_raw_cloudnet = RadiusDuration.calculate_collocation_criteria(
                raw_atl09 = raw_atl09,
                raw_cloudnet = crude_raw_cloudnet,
            )
            
            # TODO: how am I going to consider RVPolarstern as a moving platform?
            atl09_event_args.update(
                latitude = crude_raw_cloudnet.data.latitude.mean(), 
                longitude = crude_raw_cloudnet.data.longitude.mean()
            )
            cloudnet_event_args.update(
                closest_approach_time = (
                    raw_atl09["time"]
                        .isel(
                            raw_atl09["distance_to_site"].argmin(dims="...", skipna=True)
                        )
                        .median()
                )
            )
            event_list.append(
                CloudnetATL09Event(
                    atl09_event = ATL09Event(**atl09_event_args),
                    cloudnet_event = CloudnetEvent(**cloudnet_event_args)
                )
            )

        # after all ATL09 filenames are checked:
        return CollocationEventList(
            data = event_list
        )


