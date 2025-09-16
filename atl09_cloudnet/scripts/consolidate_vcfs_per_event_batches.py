"""Author: Andrew Martin
Creation date: 16/9/2025

Consolidates batches of datasets conatining vcfs per event, by concatting datasets along the collocation_event dimension.
"""

from ..definitions import indices

import os
import xarray as xr

from dataclasses import dataclass
from typing import Callable



SITES = ("ny-alesund", "hyytiala", "juelich", "munich")



@dataclass
class Args:
    site: str
    dir_vcfs: str
    index_function: Callable[[int], indices.Parametrisation]



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
        "--index-function",
        choices = indices.INDEX_FUNCTIONS.keys(),
        required=True,
        help="Name of the function mapping job array indices to (R,tau) pairs"
    )
    parsed_args = parser.parse_args()
    
    assert os.path.isdir( dir_vcfs := parsed_args.dir_vcfs)

    return Args(
        dir_vcfs = dir_vcfs,
        site = parsed_args.site,
        index_function = indices.INDEX_FUNCTIONS[parsed_args.index_function],
    )



def concat_batch_from_fpaths(fpaths: list[str]) -> xr.Dataset:
    return xr.concat(
        [
            xr.load_dataset(fpath)
            for fpath in sorted(fpaths)
        ],
        dim="collocation_event"
    )



def get_batch_from_params(params: indices.Parametrisation, fnames: list[str]) -> list[str]:
    """Given an indexed parametrisation, and a set of filenames, obtain the subset of the provided filenames that match the parametrisation for the given index.
    """
    R_km = params.distance_km
    tau_s = int( params.tau.total_seconds() )

    middle_of_fname = f"{R_km:3.3f}km_{tau_s:06}s"
    print(f"{middle_of_fname=}")
    get_middle_of_fname = lambda fname: "_".join(fname.split("_")[2:4])
    batch = [
        fname 
        for fname in fnames
        if get_middle_of_fname(fname) == middle_of_fname
    ]
    return batch



def get_filenames_for_site(site: str, dir_vcfs: str) -> list[str]:
    """Given a site and the directory containing the vcfs_per_event batch files, obtain a list of allfilenames associated with a given site.
    """
    fnames_site = [
        fname
        for fname in os.listdir(dir_vcfs)
        if (fname.split("_")[1] == site) and (fname.split("_")[4] == "batch")
    ]
    return fnames_site



def iterate_index_func(index_func):
    """Generator that stops when InvalidIndexError is raised"""
    i = 0
    while True:
        try:
            yield index_func(i)
            i += 1
        except indices.InvalidIndexError:
            return  # Stops the generator



def main(args: Args):
    print(f"Finding fnames for {args.site}")
    fnames_site = get_filenames_for_site(
        site = args.site,
        dir_vcfs = args.dir_vcfs
    )
    print(f"{len(fnames_site)=}")

    for i, params in enumerate(iterate_index_func(args.index_function)):
        print("\n")
        print(i, params)
        print("Finding fnames for batch: ", end="")
        batch_fnames = get_batch_from_params(
            params = params,
            fnames = fnames_site
        )
        print(f"{len(batch_fnames)} found")

        print(f"Concatting dataset: ", end="")
        concatted_dataset = concat_batch_from_fpaths(
            fpaths = [
                os.path.join(
                    args.dir_vcfs,
                    fname
                )
                for fname in batch_fnames
            ]
        )
        print("success")

        R_km = params.distance_km
        tau_s = int( params.tau.total_seconds() )
        fpath_out = os.path.join(
            args.dir_vcfs,
            f"vcfs-per-event_{args.site}_{R_km:3.3f}km_{tau_s:06}s.nc",
        )

        print(f"SAVING to {fpath_out}")
        concatted_dataset.to_netcdf(fpath_out)
        print("SAVING SUCCESS")
        print("")



if __name__ == "__main__":
    args = parse_args()
    print(args)
    main(args)
    print("SCRIPT SUCCESS")
