"""Author: Andrew Martin
Creation date: 30/8/2025

Script to handle loading netcdf vcfs_per_event files, and compute the mutual information between the multi-dimensional retrievals that are the VCFs as a function of height using the mutual information estimator from Holmes et al. (2019)

"""


from ..definitions import indices
from ..holmes_et_al_2019 import python_interface as holmes

from dataclasses import dataclass
import xarray as xr
import numpy as np
import os



SITES = ("ny-alesund", "hyytiala", "juelich", "munich",)



@dataclass(kw_only=True, frozen=True)
class Args:
    dir_vcfs: str
    site: str
    index_function_name: str
    dir_out: str
    K: int



def parse_args() -> Args:
    import argparse
    parser = argparse.ArgumentParser(description="Collocate ATL09 and CloudNet data")

    # Required arguments (no default values in dataclass)
    parser.add_argument(
        "--dir-vcfs",
        type=str,
        required=True,
        help="Directory path for the vcfs-per-event netcdf files"
    )

    parser.add_argument(
        "--site",
        type=str,
        choices = SITES,
        required=True,
        help="Site identifier"
    )

    parser.add_argument(
        "--dir-out",
        type=str,
        required=True,
        help="Directory path to location for saving output Mutual-Information-cube netcdf files"
    )

    parser.add_argument(
        "--index-function",
        choices = indices.INDEX_FUNCTIONS.keys(),
        required=True,
        help="Name of the function mapping job array indices to (R,tau) pairs"
    )

    parser.add_argument(
        "-K",
        type=int,
        required=True,
        help="The order of the Holmes estimator (finds the kth nearest neighbour)"
    )

    parsed_args = parser.parse_args()
    
    assert os.path.isdir( dir_vcfs := parsed_args.dir_vcfs)
    assert os.path.isdir( dir_out := parsed_args.dir_out )

    return Args(
        dir_vcfs = dir_vcfs,
        site = parsed_args.site,
        dir_out = dir_out,
        index_function_name = parsed_args.index_function,
        K = parsed_args.K
    )



def iterate_index_func(index_func):
    """Generator that stops when InvalidIndexError is raised"""
    i = 0
    while True:
        try:
            yield index_func(i)
            i += 1
        except indices.InvalidIndexError:
            return  # Stops the generator



def parameterisation_to_netcdf_fname(param: indices.Parametrisation, site: str) -> str:
    """takes a Parametrisation output of an index function and converts it to a filename according to the outputs from compute_vcfs_per_event.py
    """
    R_km = param.distance_km
    tau_s = int( param.tau.total_seconds() )
    return f"vcfs-per-event_{site}_{R_km}km_{tau_s:06}s.nc"



def main(args: Args):
    """For a given site, load all of the vcfs_per_event datasets, and calculate their MI using the Holmes estimator.
    Then, combine all datasets (including N_events and N_profiles) into a dataset with dimensions (R_km, tau_s).
    """

    index_function = indices.INDEX_FUNCTIONS[args.index_function_name]
    fpath_out = os.path.join(
        args.dir_out,
        f"MI-holmes_{args.site}.nc"
    )

    MI_das = list()
    for n, params in enumerate(iterate_index_func(index_function)):
        fpath = os.path.join(
            args.dir_vcfs,
            parameterisation_to_netcdf_fname(params, args.site)
        )
        print(f"LOADING {fpath=}", end=", ")
        ds = xr.load_dataset(fpath)
        print(f"loading-success")

        R_km = params.distance_km
        tau_s = int( params.tau.total_seconds() )

        X = ds.vcf_atl09.transpose("height","collocation_event").data
        Y = ds.vcf_cloudnet.transpose("height","collocation_event").data
        MI = holmes.call_MI_xnyn(X=X, Y=Y, K=args.K)

        MI_ds = xr.DataArray(
            data = MI,
            coords = {
                "R_km": [R_km],
                "tau_s": [tau_s]
            },
            dims=("R_km", "tau_s",)
        ).rename("mutual_information").to_dataset()


        MI_ds["n_collocation_events"] = (("R_km","tau_s",), [[ds.collocation_event.size]])
        MI_ds["n_profiles"] = (("R_km","tau_s",), [[ (ds.n_profiles_atl09*ds.n_profiles_cloudnet).sum() ]])
        print(f"dataset-generation-success")
        MI_das.append(MI_ds)

    # cmobine all individual mutual information datasets by their R and tau coordinates
    print(f"COMBINING alll datasets", end=", ")
    MI_cube = xr.combine_by_coords(MI_das)
    print(f"combining-success")

    print(f"SAVING to {fpath_out}", end=", ")
    MI_cube.to_netcdf(fpath_out)
    print("saving-success")



if __name__ == "__main__":
    args = parse_args()
    print(args)
    main(args)
    print("\n"*3)
    print("SCRIPT SUCCESS")
