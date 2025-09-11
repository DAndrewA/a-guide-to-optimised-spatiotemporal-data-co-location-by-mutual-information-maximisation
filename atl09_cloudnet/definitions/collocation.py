"""Author: Andrew Martin
Creation date: 25/08/2025

Class definitions for RadiusDuration, ATL09Event, CloudnetEvent, CloudnetATL09Event and SchemeCloudnetATL09RadiusDuration.
"""


from sat_val_framework.implement import (
    JointParameters, 
    CollocationEvent, 
    CollocationScheme,
    RawData, # importyed for type hinting
)

from sat_val_framework import CollocationEventList
from .raw_cloudnet import RawCloudnet, Duration, CloudnetEvent
from .raw_atl09 import RawATL09, DistanceFromLocation, ATL09Event

from dataclasses import dataclass
from pandas import Timestamp
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



class RadiusDuration(JointParameters):
    RAW_DATA_TYPES = (RawATL09, RawCloudnet)
    


class CollocationCloudnetATL09(CollocationEvent):
    """Class handling collocation event between ICESat-2 ATL09 and Cloudnet data."""



def _are_atl09_orbits_subsequent(fname1: str, fname2: str) -> bool:
    """Returns True if the filenames associated with ATL09 data are from subsequent orbits"""
    # processed_ATL09_20190126121018_04440201_006_02.h5
    #                           rgt -^^^^
    #                             cycle -^^
    __, __, rgt_cycle1, __, __ = fname1.split("_")
    __, __, rgt_cycle2, __, __ = fname2.split("_")
    rgt1, cycle1 = int(rgt_cycle1[:4]), int(rgt_cycle1[4:6])
    rgt2, cycle2 = int(rgt_cycle2[:4]), int(rgt_cycle2[4:6])


    condition1 = (cycle1 == cycle2) and (rgt1 + 1 == rgt2)
    condition2 = (cycle1 + 1 == cycle2) and (rgt1 == 1387) and (rgt2 == 1)
    return condition1 or condition2



@dataclass(kw_only=True, frozen=True)
class SchemeCloudnetATL09RadiusDuration(CollocationScheme):
    """Collocation scheme to find collocation events between Cloudnet and ATL09 data, based on a radius-duration scheme.
    The radius-duration scheme is described by selecting ATL09 data that falls within a given radius R of the Cloudnet site, and selecting Cloudnet data within a temporal window of duration tau, centered on the time of closest approach.

    STATIC METHODS:
        get_matches_from_raw_directories(file_list_atl09: list[str], file_list_cloudnet: list[str]) -> CollocationEventList: from a provided list of ATL09 files and Cloudnet files, identifies all available matches and returns a list of CloudnetATL09Event instances.
    """
    R_max_km: float
    min_required_atl09_profiles: int = 17 # floor (5 km of required data) / (0.28 km per profile)

    def get_matches_from_fpath_lists(self, file_list_atl09: list[str], dir_cloudnet: str, cloudnet_site: str) -> CollocationEventList:
        """From lists of file paths for ATL09 and Cloudnet data, identify collocation events and store the information using CollocationCloudnetATL09 instance
        
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
            print(f"Finding collocation for {fname_atl1}")
            
            atl09_event_args = {
                "fpath1": fname_atl1, 
                "fpath2": None,
                "min_separation_km": None,
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
            if raw_atl09 is None:
                print(f"No ATL09 data loaded, REJECTING")
                continue
            
            # obtain a crude estimate of t_0 from the median time in the ATL09 data
            crude_t0 = Timestamp(raw_atl09.data.time.mean().values).to_pydatetime()
            cloudnet_fname = RawCloudnet._fname_from_datetime(datetime = crude_t0, site=cloudnet_site)
            # check file exists during collocation event
            if not os.path.isfile( os.path.join(dir_cloudnet, cloudnet_fname) ):
                print(f"{cloudnet_fname} not found for ATL09 event, REJECTING")
                continue
            cloudnet_event_args = {
                "root_dir": dir_cloudnet,
                "site": cloudnet_site,
                "closest_approach_time": crude_t0
            }

            # load the Cloudnet data using the crude t_0 value.
            # NOTE: may produces minor errors for collocation events +-5 minutes around midnight (i.e. 0.7% of cases assuming random uniform distribution of local overpass times.)
            crude_cloudnet_event = CloudnetEvent(
                **cloudnet_event_args
            )
            crude_raw_cloudnet = RawCloudnet.from_collocation_event_and_parameters(
                event = crude_cloudnet_event,
                parameters = None 
            )
            
            d2s_finder = DistanceFromLocation(
                distance_km=self.R_max_km,
                longitude = crude_raw_cloudnet.data.longitude.mean().values,
                latitude = crude_raw_cloudnet.data.latitude.mean().values
            )
            t0 = d2s_finder.get_time_closest_approach(raw_atl09)
            d2s = d2s_finder.get_distance_to_location(raw_atl09)
            # if there are insufficient profiles
            if (d2s <= self.R_max_km).sum() < self.min_required_atl09_profiles:
                print(f"{atl09_event=} has fewer than {self.min_required_atl09_profiles} profiles within {self.R_max_km} km of {cloudnet_site}. REJECTING")
                continue

            d2s_min = float(d2s.min())
            atl09_event_args.update(
                min_separation_km = d2s_min
            )
             
            cloudnet_event_args.update(
                closest_approach_time = t0
            )

            collocation_event = CollocationEvent({
                RawATL09: ATL09Event(**atl09_event_args),
                RawCloudnet: CloudnetEvent(**cloudnet_event_args)
            })
            print("EVENT FOUND:",collocation_event)

            event_list.append(collocation_event)

        # after all ATL09 filenames are checked:
        return CollocationEventList(
            data = event_list
        )


