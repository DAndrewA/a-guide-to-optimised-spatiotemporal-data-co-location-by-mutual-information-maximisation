"""Author: Andrew Martin
Creation date: 29/8/2025

Script to handle finding collocation events between ATL09 and Cloudnet data.
"""


from ..definitions.collocation import (
    DistanceFromLocation, Duration, RadiusDuration,
    RawATL09, RawCloudnet
)
from ..definitions import (
    vcf,
    indices
)
from sat_val_framework import CollocationEventList

from typing import Callable, Any
from dataclasses import dataclass
import datetime as dt
import os
import xarray as xr

SITES = {
    "ny-alesund": dict(lat=78.923, lon=11.922), 
    "hyytiala": dict(lat=61.844, lon=24.287), 
    "juelich": dict(lat=50.908, lon=6.413), 
    "munich": dict(lat=48.148, lon=11.573),
}




@dataclass(kw_only=True, frozen=True)
class Args:
    fpath_pickle: str
    site: str
    dir_output: str
    joint_params: RadiusDuration


def parse_args() -> None | Args:
    import argparse
    parser = argparse.ArgumentParser(description="Collocate ATL09 and CloudNet data")

    # Required arguments (no default values in dataclass)
    parser.add_argument(
        '--pickle-dir',
        type=str,
        required=True,
        help='Directory path for the pickled CollocationEventList files'
    )

    parser.add_argument(
        '--site',
        type=str,
        choices = SITES.keys(),
        required=True,
        help='Site identifier'
    )

    parser.add_argument(
        '--output-dir',
        type=str,
        required=True,
        help='Directory path to location for saving output pickle files'
    )

    parser.add_argument(
        '--job-array-index',
        type=int,
        required=True,
        help='Job array index'
    )

    parser.add_argument(
        "--index-function",
        choices = indices.INDEX_FUNCTIONS.keys(),
        required=True,
        help="Name of the function mapping job array indices to (R,tau) pairs"
    )


    parsed_args = parser.parse_args()
    
    site = parsed_args.site

    # get R, tau from array index
    index_function = indices.INDEX_FUNCTIONS[parsed_args.index_function]
    R_tau = index_function(parsed_args.job_array_index)
    joint_params = RadiusDuration({
        RawATL09: DistanceFromLocation(
            distance_km = R_tau.distance_km,
            latitude = SITES[site]["lat"],
            longitude = SITES[site]["lon"]
        ),
        RawCloudnet: Duration(duration = R_tau.tau)
    })

    # Convert argparse Namespace to Args dataclass
    # Note: argparse uses dashes, but dataclass uses underscores
    fpath_pickle = os.path.join(
        parsed_args.pickle_dir, f"collocation_events_{site}.pkl" 
    )
    assert os.path.isfile(fpath_pickle)
    assert os.path.isdir(parsed_args.output_dir)

    return Args(
        fpath_pickle = fpath_pickle,
        site = site,
        dir_output = parsed_args.output_dir,
        joint_params = joint_params
    )



def pair_vcf_dataset_from_collocated_homogenised_data(coll_H):
    """From vcf.VCF instances, join the [RawATL09 | RawCloudnet] elements in a single xarray dataset.
    A deep copy is made before returning, so that no stray references to the collocated homogenised data are held onto.
    """
    vcf_atl09 = coll_H[RawATL09].data
    vcf_cloudnet = coll_H[RawCloudnet].data
    return xr.merge([
        vcf_atl09.rename("vcf_atl09"), vcf_cloudnet.rename("vcf_cloudnet")
    ]).copy(deep=True)



def make_safe(func: Callable[[Any],Any]) -> Callable[[Any],Any|None]:
    """Make a provided function safe, by wrapping it in a try except statement.
    If no exceptions are thrown, the return is [value, True].
    If an exception is thrown, the return is [message, False].
    """
    def safe_func(*args, **kwargs):
        try:
            value = func(*args, **kwargs)
            success = True
        except Exception as e:
            value = "\n".join((
                f"make_safe({func.__qualname__})",
                f"Thrown {type(e)}:",
                str(e)
            ))
            success = False
        return (value, success)
    return safe_func

# lambda to filter None from generator expressions
filter_none = lambda generator: (v for v in generator if v is not None)



def main(args: Args):
    all_collocation_events = CollocationEventList.from_file(args.fpath_pickle)
    print(f"Loaded {args.fpath_pickle}, with {len(all_collocation_events)} events in total")


    R_km = args.joint_params[RawATL09].distance_km
    tau_s = int( args.joint_params[RawCloudnet].duration.total_seconds() )
    fpath_out = os.path.join(
        args.dir_output,
        f"vcfs-per-event_{args.site}_{R_km}km_{tau_s:06}s.nc"
    )


    # forgive me for the attrocities that are about to be committed to code
    events_to_load = (
        (
            event, 
            print("\n",f"LOADING",event)
        )[0]
        for event in all_collocation_events
    )

    # reject events when the minimum ICESat-2 separation is gretaer than R_km:
    N_reject_minimum_R = 0
    events_to_load = filter_none((
        event

        if event[RawATL09].min_separation_km <= R_km
        else (
            None,
            print(f"ABORTING due to min_separation_km={event[RawATL09].min_separation_km} > {R_km}"),
            N_reject_minimum_R := N_reject_minimum_R + 1
        )[0]

        for event in events_to_load
    ))

    # load an event to a CollocatedRawData instance, and catch exceptions in loading
    N_fail_on_load = 0
    collocated_raw_data = filter_none((
        raw_data[0] # make_safe value

        if (
            raw_data := make_safe(collocation_event.load_with_joint_parameters)(args.joint_params)
        )[1] # make_safe success value

        else (
            None,
            print(raw_data[0]), # make_safe message
            print(f"ABORTING due to failed loading"),
            N_fail_on_load:=N_fail_on_load+1,
        )[0]

        for collocation_event in events_to_load
    ))

    # filter for any of the elements of CollocatedRawData instances being None
    N_fail_with_raw_data_None = 0
    collocated_raw_data = filter_none((
        raw_data

        if not any(
            value is None
            for value in (raw_data).values()
        )
        else (
            None,
            print(f"ABORTING: None value found in collocated raw data, { {k: (type(v) is None) for k,v in raw_data.items()}= } (likely no data after subsetting)"),
            N_fail_with_raw_data_None:=N_fail_with_raw_data_None+1
        )[0]

        for raw_data in collocated_raw_data
    ))

    # homogenise data, and catch exceptions
    N_fail_homogenisation = 0
    collocated_homogenised_data_and_n_profiles = filter_none((
        (
            collocated_h[0], # make_safe value
            b:={k: rd.n_profiles for k, rd in raw_data_pair.items()},
        )

        if (
            collocated_h := make_safe(raw_data_pair.homogenise_to)(vcf.VCF_240m)
        )[1] # make_safe success value
        else (
            None,
            print(collocated_h[0]), # make_safe message
            print(f"ABORTING: error thrown in homogenisation"),
            N_fail_homogenisation:=N_fail_homogenisation+1,
        )[0]

        for raw_data_pair in collocated_raw_data
    ))

    # create vcf data arrays from succesfully homogenised data
    N_success = 0
    n_profiles = {
        RawATL09: list(),
        RawCloudnet: list(),
    }
    paired_vcfs = (
        (
            pair_vcf_dataset_from_collocated_homogenised_data(coll_H),
            N_success := N_success + 1,
            a := n_profiles[RawATL09].append( num_profiles_dict[RawATL09] ),
            b := n_profiles[RawCloudnet].append( num_profiles_dict[RawCloudnet] ),
            print("Number of vertical profiles in event:", f"{num_profiles_dict}"),
            print(f"SUCCESS: homogenised data loaded to xarray DataArray"),
        )[0]

        for coll_H, num_profiles_dict in collocated_homogenised_data_and_n_profiles
    )

    vcfs_per_event = xr.concat(
        paired_vcfs,
        dim="collocation_event"
    )
    vcfs_per_event["n_profiles_atl09"] = ( ("collocation_event",) , n_profiles[RawATL09] )
    vcfs_per_event["n_profiles_cloudnet"] = ( ("collocation_event",) , n_profiles[RawCloudnet] )

    print(f"SAVING TO {fpath_out=}")
    vcfs_per_event.to_netcdf(fpath_out)
    print(f"SAVING SUCCESS")

    # print additional information like the number of successes and failures
    print("Additional information")
    for label, value in zip(
        ("loading", "none after laoding", "homogenissaton", "success", "total"),
        (N_fail_on_load, N_fail_with_raw_data_None, N_fail_homogenisation, N_success, len(all_collocation_events))
    ):
        print(f"{label:>25} | {value}")

    print(n_profiles)

    return



if __name__ == "__main__":
    args = parse_args()
    print(args)
    main(args=args)
    print("\n\n\n\n\n")
    print("SCRIPT SUCCESS")
