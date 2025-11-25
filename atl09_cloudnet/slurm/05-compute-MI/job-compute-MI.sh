#!/bin/bash

run_from_dir=../../..
script_name="atl09_cloudnet.scripts.compute_MI_with_confidence"

dir_vcfs=$SCRATCH/vcfs_per_event
out_dir=$SCRATCH/MI
site=$1
K=$2

mkdir $out_dir

mamba activate overpass_analysis_again

echo "site = ${site}"
echo "K = ${K}"
echo "dir_vcfs = ${dir_vcfs}"
echo "dir_out = ${dir_out}"
echo "array-index = ${SLURM_ARRAY_TASK_ID}"

(cd $run_from_dir && \
    time python -m $script_name \
        --dir-vcfs $dir_vcfs \
        --dir-out $out_dir \
        --site $site \
        --index-function "R_500km_tau_172800s" \
        --job-array-index $SLURM_ARRAY_TASK_ID \
        -K $K \
        --n-bootstraps 100 \
        --n-splits 10 \
        --n-B-repeats 20 \
)
