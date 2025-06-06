#!/bin/bash
#SBATCH --job-name=rxn
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --time=0-1:00:00

conda activate /ihome/kjohnson/ska31/.conda/envs/deepmdGPU 

python run_neb.py
crc-job-stats
