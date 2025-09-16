"""Author: Andrew Martin
Creation date: 16/9/2025

Script to handle computationof Mutual Information for the Vertical Cloud Fraction datasets using the following methods:
    + Holmes estimator of full (2x50)-dimensional dataset
    + Holmes estimator for 50 x 2-dimensional estimation per height level
    + Histogram with equal bin widths for 50 x 2-dimensional per height estimation
    + Histogram with equal bin counts for 50 x 2-dimensional per height estimation
"""

from ..definitions import indices
from ..holmes_et_al_2019 import python_interface as holmes
from xhistogram.xarray import histogram as xhist

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
    n_bins: int
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
        "--N-bins",
        type=int,
        required=True,
        help="The number of bins when computing histograms."
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
        n_bins = parsed_args.N_bins,
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
    return f"vcfs-per-event_{site}_{R_km:3.3f}km_{tau_s:06}s.nc"



# TODO: assess impact of choice of BINS_VCF 
def counts_from_raw_vcfs_per_profile(ds: xr.Dataset, bins_atl09: np.ndarray, bins_cloudnet) -> xr.DataArray:
    """Function that maps a loaded vcfs_per_event dataset to a histogram of joint vcf values per height.
    """
    counts = xhist(
        ds.vcf_atl09,
        ds.vcf_cloudnet,
        bins=[bins_atl09, bins_cloudnet],
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



def mutual_information_with_height_from_vcfs(ds: xr.Dataset, bins1: np.ndarray, bins2: np.ndarray, dim: list[str] = ["vcf_atl09", "vcf_cloudnet"]) -> xr.DataArray:
    counts = counts_from_raw_vcfs_per_profile(
        ds=ds,
        bins_atl09=bins1,
        bins_cloudnet=bins2
    )
    MI = mutual_information_with_height(
        counts = counts,
        dim = [
            dimname+"_bin"
            for dimname in dim
        ]
    )
    return MI



def main(args: Args):
    """For a given site, loads all of the vcfs per event datasets, and calculates the 4 MI values for each.
    Then, combines all datasets (including N_events and N_profiles) into a dataset with dimensions (R_km, tau_s, height).
    """

    index_function = indices.INDEX_FUNCTIONS[args.index_function_name]
    fpath_out = os.path.join(
        args.dir_out,
        f"MI_{args.site}.nc"
    )
    print(f"{fpath_out=}")

    MI_datasets = list()

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
        
        # data required for Holmes estimators
        X = ds.vcf_atl09.transpose("height","collocation_event").data
        Y = ds.vcf_cloudnet.transpose("height","collocation_event").data

        # 1) Holmes - total
        print("MI holmes total: ", end="")
        MI_holmes_total = holmes.call_MI_xnyn(X=X, Y=Y, K=args.K)
        print("success")

        # 2) Holmes - per height
        print("MI holmes per height: ", end="")
        MI_holmes_per_height = np.array([
            holmes.call_MI_xnyn(X=Xrow, Y=Yrow, K=args.K)
            for Xrow, Yrow in zip(X,Y)
        ])
        print("success")

        # bins required for histograms
        lift_degeneracy = lambda a: np.unique(a)
        BINS_equal_width = np.linspace(0,1, args.nbins+1)
        BINS_equal_counts_atl09 = lift_degeneracy(np.quantile(X, BINS_equal_width))
        BINS_equal_counts_cloudnet = lift_degeneracy(np.quantile(Y, BINS_equal_width))

        # 3) Histogram - equal width
        print("MI hist equal width: ", end="")
        MI_hist_equal_width = (
            mutual_information_with_height_from_vcfs(
                ds = ds,
                bins1 = BINS_equal_width,
                bins2 = BINS_equal_width,
            )
                .rename("MI_hist_equal_width")
        )
        print("sucess")

        # 4) Histogram - equal count
        print("MI hist equal count: ", end="")
        MI_hist_equal_count = (
            mutual_information_with_height_from_vcfs(
                ds = ds,
                bins1 = BINS_equal_counts_atl09,
                bins2 = BINS_equal_counts_cloudnet,
            )
                .rename("MI_hist_equal_count")
        )
        print("success")

        # create results dataset
        print(f"creating results dataset: ", end="")
        results_dataset = xr.Dataset(
            data_vars = {
                "MI_holmes_total": MI_holmes_total,  # float requires no extra dims
                "MI_holmes_per_height": (["height"], MI_holmes_per_height), # array requires dims specifying
                "MI_hist_equal_width": MI_hist_equal_width, # data array requires no further dims specifying
                "MI_hist_equal_count": MI_hist_equal_count, # (above)
                "N_events": int(ds.collocation_event.size), 
                "N_profiles": int(
                    (ds.n_profiles_atl09 * ds.n_profiles_cloudnet).sum()
                ),
            },
            coords = {
                "height": ds.coords["height"]
            }
        ).expand_dims(
            dim = {
                "R_km": [R_km],
                "tau_s": [tau_s],
            }
        )
        print("suceess")

        MI_datasets.append(results_dataset)

    print("\n"*3)
    # combine all individual results datasets by their R and tau coordinates
    print("COMBINING datasets: ", end="")
    MI_cube = xr.combine_by_coords(MI_datasets)
    print("success")

    print(f"SAVING dataset to {fpath_out=}")
    MI_cube.to_netcdf(fpath_out)
    print("SAVING SUCCESS")



if __name__ == "__main__":
    args = parse_args()
    print(args)
    main(args)
    print("SCRIPT SUCCESS")
