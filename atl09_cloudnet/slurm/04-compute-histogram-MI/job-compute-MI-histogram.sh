#!/bin/bash

run_from_dir=../../..
script_name="atl09_cloudnet.scripts.compute_MI_histogram"

dir_vcfs=$SCRATCH/vcfs_per_event
out_dir=$SCRATCH/MI_histogram
site=$1

mkdir $out_dir

mamba activate overpass_analysis_again

echo "site = ${site}"
echo "dir_vcfs = ${dir_vcfs}"
echo "dir_out = ${dir_out}"

(cd $run_from_dir && \
    time python -m $script_name \
        --dir-vcfs $dir_vcfs \
        --dir-out $out_dir \
        --site $site \
        --index-function "R_150km_tau_172800s"
)
