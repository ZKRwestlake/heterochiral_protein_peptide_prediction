#!/bin/sh 
#SBATCH -p test,cpu,gpu_L40
#SBATCH -J dockq
#SBATCH -o dockq.log
#SBATCH -e dockq.err
#SBATCH -c 6
##SBATCH -a 1
##SBATCH --mem 1g
# set -e

. ~/.bashrc
conda activate CHAI_ENV
# Use the env's libstdc++ so scipy can be imported on older system glibc.
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

dirs=(`ls -d af_output/*/`)

for dir in ${dirs[@]};do
    echo $dir
    python dockq_iptm.py $dir
    python dockq_ddg.py $dir
done
