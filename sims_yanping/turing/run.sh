#!/bin/bash

#SBATCH --job-name="nextSimTest"
#SBATCH -D /home/ypei/obsReLOOP
#SBATCH --output=logs/nextSimTest_%A_%a.out
#SBATCH --error=logs/nextSimTest_%A_%a.err
#SBATCH --array=1-50
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --partition=short
#SBATCH --exclude=compute-3-01,compute-4-01

module purge
module load r

mkdir -p logs
mkdir -p results/next

Rscript /home/ypei/obsReLOOP/yhatRF.R $SLURM_ARRAY_TASK_ID 50
