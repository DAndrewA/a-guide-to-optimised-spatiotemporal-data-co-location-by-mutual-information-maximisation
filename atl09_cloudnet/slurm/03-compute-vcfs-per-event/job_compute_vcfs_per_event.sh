#!/bin/bash
script_name=atl09_cloudnet.scripts.compute_vcfs_per_event
run_from_dir=../../..

pickle_dir=$SCRATCH/get_collocation_events_pickles
output_dir=$SCRATCH/vcfs_per_event

site=$1

mamba activate overpass_analysis_again

echo $site $SLURM_ARRAY_TASK_ID

(cd $run_from_dir && time python -m $script_name --pickle-dir $pickle_dir --site $site --output-dir $output_dir --job-array-index=$SLURM_ARRAY_TASK_ID --index-function R_150km_tau_172800s)

