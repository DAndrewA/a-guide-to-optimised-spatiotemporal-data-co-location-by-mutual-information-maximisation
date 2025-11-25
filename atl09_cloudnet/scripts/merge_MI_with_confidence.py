"""Author: Andrew Martin
creation date: 14/10/25

Script to handle merging of MI datasets for each site, for each value of R, tau and K.

Can be ran independently of SLURM submissions.
"""

import xarray as xr
import os
from itertools import product
from typing import Callable
from dataclasses import dataclass
from tqdm import tqdm

from ..definitions import indices


SITES = ("ny-alesund", "hyytiala", "juelich", "munich")
K = (1,2,4,10,20)


@dataclass
class Args:
    dir_mi: str
    index_function_name: Callable[[int], indices.Parametrisation]


def parse_args() -> Args:
    import argparse
    parser = argparse.ArgumentParser(description="Collocate ATL09 and CloudNet data")

    # Required arguments (no default values in dataclass)
    parser.add_argument(
        "--dir-MI",
        type=str,
        required=True,
        help="Directory path containing the individual MI datasets"
    )

    parser.add_argument(
        "--index-function",
        choices = indices.INDEX_FUNCTIONS.keys(),
        required=True,
        help="Name of the function mapping job array indices to (R,tau) pairs"
    )

    parsed_args = parser.parse_args()
    
    assert os.path.isdir( dir_mi := parsed_args.dir_MI)

    return Args(
        dir_mi = dir_mi,
        index_function_name = parsed_args.index_function
    )



def main(args: Args):
    index_function = indices.INDEX_FUNCTIONS[args.index_function_name]
    datasets = list()
    for site, k, jobnum in tqdm(product(SITES, K, range(index_function.MAX_INDEX))):
        fpath_load = os.path.join(
            args.dir_mi,
            f"MI_{site}_K{k}_job{jobnum}.nc"
        )
        if not os.path.exists(fpath_load):
            print(f"Cannot find {fpath_load=}")
            continue
        #print(f"Loading {fpath_load=}", end="")
        datasets.append(
            xr.load_dataset(fpath_load).expand_dims(
                dim = {
                    "site": [site]
                }
            )
        )
        #print("  success")

    print("Combining datasets", end="")
    total_ds = xr.combine_by_coords(
        datasets
    )
    print("success")
    total_ds.to_netcdf(
        os.path.join(
            args.dir_mi,
            "MI_merged.nc"
        )
    )




if __name__ == "__main__":
    args = parse_args()
    print(args)
    main(args)
    print("SUCCESS")


