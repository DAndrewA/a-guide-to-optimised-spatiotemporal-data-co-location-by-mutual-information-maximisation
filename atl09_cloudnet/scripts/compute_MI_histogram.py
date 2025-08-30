"""Author: Andrew Martin
Creation date: 30/8/25

Script to handle loading netcdf files containing vcfs_per_event, and compute mean vertically-averaged (mutual information as a function of height) using histograms of VCF at each height.
This method is inferior to the estimator outlined in Holmes et al. (2019), but significantly easier to implement.
"""

from ..definitions import indices

from dataclasses import dataclass
import xarray as xr
import numpy as np
import os
from xhistogram.xarray import histogram as xhist



SITES = ("ny-alesund", "hyytiala", "juelich", "munich",)
DEFAULT_BINS = 100



@dataclass(kw_only=True, frozen=True)
class Args:
    dir_vcfs: str
    site: str
    index_function_name: str
    dir_out: str



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

    parsed_args = parser.parse_args()
    
    assert os.path.isdir( dir_vcfs := parsed_args.dir_vcfs)
    assert os.path.isdir( dir_out := parsed_args.dir_out )

    return Args(
        dir_vcfs = dir_vcfs,
        site = parsed_args.site,
        dir_out = dir_out,
        index_function_name = parsed_args.index_function
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



# TODO: assess impact of choice of BINS_VCF 
def counts_from_raw_vcfs_per_profile(ds: xr.Dataset, BINS_VCF: int | np.ndarray = DEFAULT_BINS) -> xr.DataArray:
    """Function that maps a loaded vcfs_per_event dataset to a histogram of joint vcf values per height.
    """
    counts = xhist(
        ds.vcf_atl09,
        ds.vcf_cloudnet,
        bins=[BINS_VCF, BINS_VCF],
        dim=["collocation_event"],
    )
    return counts



def fill_non_finite(da: xr.DataArray, fill_val = 0) -> xr.DataArray:
    """Fills nan and inf values in xarray DataArrays with fill_val"""
    return xr.where(
        np.isfinite(da),
        da,
        fill_val
    ).copy()



def mutual_information_with_height(counts: xr.DataArray, dim: list[str] = ["vcf_atl09_bin", "vcf_cloudnet_bin"]) -> xr.DataArray:
    """Directly computes the mutual information as a function of height between the ATL09 and Cloudnet VCF distributions.
    """
    probability = counts / counts.sum(dim=dim)
    log_term = fill_non_finite(
        np.log2(probability)
        - np.log2(probability.sum(dim=dim[0]))
        - np.log2(probability.sum(dim=dim[1]))
    )
    return (
        (probability * log_term)
            .sum(dim=dim)
            .rename("mutual_information")
    )



def main(args: Args):
    """For a given site, load all of the vcfs_per_event datasets, and calculate the MI as a function of height.
    Then, combine all of the datasets (including N_events and N_profiles) into a dataset cube, with dimensions (height, R_km, tau_s).
    """
    index_function = indices.INDEX_FUNCTIONS[args.index_function_name]
    fpath_out = os.path.join(
        args.dir_out,
        f"MI-cube-by-histogram_{args.site}.nc"
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

        print("calculating MI", end=", ")
        MI_ds = (
            mutual_information_with_height(
                counts = counts_from_raw_vcfs_per_profile(ds)
            )
                .expand_dims(
                    dim={
                        "R_km": [R_km],
                        "tau_s": [tau_s]
                    }
                )
        ).to_dataset()
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
    print(f"\n"*3)
    print("SCRIPT SUCCESS")
