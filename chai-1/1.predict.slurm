#!/bin/bash							
##SBATCH -q huge
#SBATCH -J c1
#SBATCH -p gpu_L40,gpu
#SBATCH -o c1.log
#SBATCH -e c1.err
#SBATCH --gres=gpu:1
#SBATCH -a 1
#SBATCH -c 2
##SBATCH --mem 1g

. ~/.bashrc

dir=`sed -n ${SLURM_ARRAY_TASK_ID}p dirs.list`

python predict.py $dir
