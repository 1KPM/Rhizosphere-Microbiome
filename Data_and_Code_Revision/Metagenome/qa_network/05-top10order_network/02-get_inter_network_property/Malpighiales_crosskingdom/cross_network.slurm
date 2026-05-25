#!/bin/bash
#SBATCH --job-name=network_analysis
#SBATCH --partition=56c503g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=16G
#SBATCH --account=wangxiaolin_lab
#SBATCH --error=%x_%j.err
#SBATCH --output=%x_%j.out

source ~/.bashrc
conda activate r_env
Rscript Malpighiales_cross.R
Rscript Plot.R
