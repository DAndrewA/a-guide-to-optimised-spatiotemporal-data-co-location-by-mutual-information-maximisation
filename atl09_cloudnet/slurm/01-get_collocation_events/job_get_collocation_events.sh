#!/bin/bash
script_name=atl09_cloudnet.scripts.get_collocation_events
run_from_dir=../../..

sites_dir=$MI_MAXIMISATION_RAW_DATA_DIRECTORY/sites
output_dir=$MI_MAXIMISATION_RESULTS_DIRECTORY/get_collocation_events_pickles

site=$1

mamba activate overpass_analysis_again

echo $site $SLURM_ARRAY_TASK_ID

(cd $run_from_dir && time python -m $script_name --dir-atl09 $sites_dir/$site/atl09 --dir-cloudnet $sites_dir/$site/cloudnet --site $site --output-dir $output_dir --job-array-index=$SLURM_ARRAY_TASK_ID --R-min-km 510)
