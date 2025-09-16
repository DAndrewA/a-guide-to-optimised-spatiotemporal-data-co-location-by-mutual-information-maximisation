#!/bin/bash
script_name=atl09_cloudnet.scripts.consolidate_vcfs_per_event_batches
run_from_dir=../../..

vcfs_dir=$SCRATCH/vcfs_per_event

site=$1

mamba activate overpass_analysis_again

echo $site

(cd $run_from_dir && time python -m $script_name --dir-vcfs $SCRATCH/vcfs_per_event --site $site --index-function R_500km_tau_172800s)

