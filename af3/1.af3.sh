#!/bin/bash
#SBATCH -J af3
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -p gpu_L40
#SBATCH --gres=gpu:1
#SBATCH -o af3.log
#SBATCH -e af3.err

pwd=$PWD
mkdir -p af_output

cd /public/software-share/lulab/alphafold3
singularity exec \
 --nv \
 --bind $pwd/af_input:/root/af_input \
 --bind $pwd/af_output:/root/af_output \
 --bind /data/alphafold_3.0.1:/root/models \
 --bind /data/alphafold_3.0.1/alphafold3_data/public_databases:/root/public_databases \
 /data/alphafold_3.0.1/alphafold3.sif \
 python run_alphafold.py \
 --model_dir=/root/models \
 --db_dir=/root/public_databases \
 --input_dir=/root/af_input \
 --output_dir=/root/af_output \
