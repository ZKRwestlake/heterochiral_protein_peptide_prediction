#!/bin/sh
set -e
. ~/.bashrc
conda activate CHAI_ENV
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

dirs=(`ls -d prediction/*/`)

for dir in ${dirs[@]};do
    echo $dir
    python collect_best_scored_model.py $dir
    python collect_most_confident_model.py $dir
done
