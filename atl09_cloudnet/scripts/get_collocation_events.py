"""Author: Andrew Martin
Creation date: 29/8/2025

Script to handle finding collocation events between ATL09 and Cloudnet data.
"""


from ..definitions.collocation import SchemeCloudnetATL09RadiusDuration

from dataclasses import dataclass
import os


SCHEME = SchemeCloudnetATL09RadiusDuration(
    R_max_km=150
)



@dataclass(kw_only=True, frozen=True)
class Args:
    dir_atl09: str
    dir_cloudnet: str
    site: str
    output_dir: str
    job_array_index: int
    slice_length: int = 120
    R_min_km: float = 150
    minimum_required_profiles_within_R: int = 17


def parse_args() -> None | Args:
    import argparse
    parser = argparse.ArgumentParser(description="Collocate ATL09 and CloudNet data")

    # Required arguments (no default values in dataclass)
    parser.add_argument(
        '--dir-atl09',
        type=str,
        required=True,
        help='Directory path for ATL09 data'
    )

    parser.add_argument(
        '--dir-cloudnet',
        type=str,
        required=True,
        help='Directory path for CloudNet data'
    )

    parser.add_argument(
        '--site',
        type=str,
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
        '--slice-length',
        type=int,
        default = 120,
        help='Number of ATL09 files to be processed'
    )

    # Arguments with default values
    parser.add_argument(
        '--R-min-km',
        type=float,
        default=150,
        help='Minimum radius in kilometers (default: 150)'
    )

    parser.add_argument(
        '--minimum-required-profiles-within-R',
        type=int,
        default=17,
        help='Minimum required profiles within radius (default: 17)'
    )

    parsed_args = parser.parse_args()

    # Convert argparse Namespace to Args dataclass
    # Note: argparse uses dashes, but dataclass uses underscores
    assert os.path.isdir(parsed_args.dir_atl09)
    assert os.path.isdir(parsed_args.dir_cloudnet)

    return Args(
        dir_atl09=parsed_args.dir_atl09,
        dir_cloudnet=parsed_args.dir_cloudnet,
        site=parsed_args.site,
        output_dir=parsed_args.output_dir,
        job_array_index=parsed_args.job_array_index,
        slice_length = parsed_args.slice_length,
        R_min_km=parsed_args.R_min_km,
        minimum_required_profiles_within_R=parsed_args.minimum_required_profiles_within_R
    )



def main(args: Args):
    """Function that, given the ATL09 directory, indices for which files should be selected, Cloudnet directory and site, and output locations, finds the collocation events between ATL09 and Cloudnet data, and saves the collocation event lists to dedicated pickle files.
    """
    SLICE_SIZE = args.slice_length
    atl09_indices = slice(
        SLICE_SIZE * args.job_array_index,
        SLICE_SIZE * (args.job_array_index + 1)
    )
    fpaths_atl09 = [
        os.path.join(args.dir_atl09, fname)
        for fname in sorted(os.listdir(args.dir_atl09))
        if fname[-3:] == ".h5"
    ][atl09_indices]

    if not fpaths_atl09:
        print(f"No ATL09 files found for indices {atl09_indices} in directory {args.dir_atl09}.")
        print("TERMINATING")
        return

    print(f"Finding collocation events for {fpaths_atl09}")

    collocation_event_list = SCHEME.get_matches_from_fpath_lists(
        file_list_atl09 = fpaths_atl09,
        dir_cloudnet = args.dir_cloudnet,
        cloudnet_site = args.site
    )

    print(f"{len(collocation_event_list)} collocation events found:")
    for event in collocation_event_list:
        print("===>      ", event, "\n")

    out_fpath = os.path.join(
        args.output_dir,
        f"collocation-events_{args.site}_{args.job_array_index:04}.pkl"
    )

    print(f"attempting to save to {out_fpath=}")
    collocation_event_list.to_file(out_fpath)


if __name__ == "__main__":
    args = parse_args()
    print(args)
    main(args=args)
    print("SCRIPT SUCCESS")
