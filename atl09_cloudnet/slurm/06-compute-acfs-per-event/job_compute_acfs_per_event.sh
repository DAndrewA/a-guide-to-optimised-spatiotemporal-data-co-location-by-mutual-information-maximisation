#!/bin/bash
script_name=atl09_cloudnet.scripts.compute_acfs_per_event
run_from_dir=../../..

pickle_dir=$SCRATCH/get_collocation_events_pickles
output_dir=$SCRATCH/acfs_per_event

site=$1
shift 1
remaining_arguments=$@

mamba activate overpass_analysis_again

echo $site $SLURM_ARRAY_TASK_ID
echo "reminaing arguments" $remaining_arguments

(cd $run_from_dir && time python -m $script_name --pickle-dir $pickle_dir --site $site --output-dir $output_dir --job-array-index=$SLURM_ARRAY_TASK_ID --index-function R_tau_extremal_lit_opt_$site $remaining_arguments)

