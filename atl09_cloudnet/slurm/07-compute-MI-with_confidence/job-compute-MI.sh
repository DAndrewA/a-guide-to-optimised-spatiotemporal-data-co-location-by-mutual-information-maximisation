#!/bin/bash

run_from_dir=../../..
script_name="atl09_cloudnet.scripts.compute_MI_with_confidence"

dir_vcfs=$MI_MAXIMISATION_RESULTS_DIRECTORY/vcfs_per_event
out_dir=$MI_MAXIMISATION_RESULTS_DIRECTORY/MI
site=$1

mkdir $out_dir

mamba activate overpass_analysis_again

echo "site = ${site}"
echo "dir_vcfs = ${dir_vcfs}"
echo "dir_out = ${out_dir}"

(cd $run_from_dir && \
    time python -m $script_name \
        --dir-vcfs $dir_vcfs \
        --dir-out $out_dir \
        --site $site \
        --index-function "R_500km_tau_172800s" \
        -K 20 \
        --n-bootstraps 100 \
        --n-splits 10
)
