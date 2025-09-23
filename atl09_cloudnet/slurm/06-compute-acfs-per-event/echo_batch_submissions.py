"""Author: Andrew Martin
Creation date: 23/9/2025

Script to handle creation of submission commands for batching compute_acfs_per_event
"""

import pickle
import os

SCRATCH = os.environ["SCRATCH"]
PICKLE_DIR = os.path.join(SCRATCH, "get_collocation_events_pickles")
OUTPUT_DIR = os.path.join(SCRATCH, "acfs_per_event")
SCRATCH_DIR = os.path.join(SCRATCH, "acfs_per_event_oe")


SITES = ("ny-alesund", "hyytiala", "juelich", "munich")
LENGTHS = (7977, 3156, 2480, 2304)

BATCH_SIZE = 1000
MAX_INDEX = 6

for site, N_events in zip(SITES, LENGTHS):
    n_batches = (N_events // BATCH_SIZE) + 1
    max_array_index = n_batches * MAX_INDEX

    job_name = f"compute_acfs_per_event_{site}"

    slurm_command=" ".join([
        "sbatch", 
        "--partition=standard",
        "--qos=standard",
        "--account=icecaps",
        "--time=08:00:00",
        "--mem=20G",
        f"--job-name={job_name}",
        f"-o {SCRATCH_DIR}/{job_name}_%a.o",
        f"-e {SCRATCH_DIR}/{job_name}_%a.e",
        f"--array=0-{max_array_index}",
        f"job_compute_acfs_per_event.sh {site}",
        f"--batch --batch-size {BATCH_SIZE}"
    ])

    print(r"#" + f"{site}")
    print(slurm_command)
    print("")


