"""Author: Andrew Martin
Creation date: 29/08/2025

Script that allows for collocation event list pickle files to be combined within a directory for a given site
"""


from sat_val_framework import CollocationEventList
from ..definitions.collocation import *

from dataclasses import dataclass
import os


@dataclass(frozen=True, kw_only=True)
class Args:
    pickle_dir: str
    site: str

def parse_args() -> Args:
    import argparse
    parser = argparse.ArgumentParser(description="Consolidate collocation event lists for a given site")

    # Required arguments (no default values in dataclass)
    parser.add_argument(
        '--site',
        type=str,
        required=True,
        help='Site identifier'
    )

    parser.add_argument(
        '--pickle-dir',
        type=str,
        required=True,
        help='Directory containing collocation event list pickle files to be combined.'
    )

    parsed_args = parser.parse_args()

    # Convert argparse Namespace to Args dataclass
    # Note: argparse uses dashes, but dataclass uses underscores
    assert os.path.isdir(parsed_args.pickle_dir)

    return Args(
        site = parsed_args.site,
        pickle_dir = parsed_args.pickle_dir
    )




def main(args: Args):
    pickle_files_for_site = sorted([
        fname
        for fname in os.listdir(args.pickle_dir)
        if fname[-4:] == ".pkl" and fname.split("_")[1] == args.site
    ])
    print(f"Loading from {pickle_files_for_site}")

    collocation_event_lists = list()

    for fname in pickle_files_for_site:
        fpath = os.path.join(args.pickle_dir, fname)

        with open(fpath, "rb") as f:
            collocation_event_list = CollocationEventList.from_file(fpath)
            collocation_event_lists.append(collocation_event_list)

    all_collocation_events = CollocationEventList([
        event
        for event_list in collocation_event_lists
        for event in event_list
    ])

    output_fname = f"collocation_events_{args.site}.pkl"
    output_fpath = os.path.join(args.pickle_dir, output_fname)

    print(f"{len(all_collocation_events)=}")
    print(f"Saving collocation events to {output_fpath}")
    all_collocation_events.to_file(output_fpath)


if __name__ == "__main__":
    args = parse_args()
    print(args)
    main(args)
    print("SCRIPT SUCCESS")




