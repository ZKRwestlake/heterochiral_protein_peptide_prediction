#!/bin/sh
# Source the shell profile before enabling set -e: profiles can contain commands
# that exit non-zero, which set -e would otherwise treat as fatal.
. ~/.bashrc
set -e
conda activate CHAI_ENV
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

dirs=(`ls -d af_output/*/`)

for dir in ${dirs[@]};do
    echo $dir
    python collect_best_scored_model.py $dir
    python collect_most_confident_model.py $dir
done
