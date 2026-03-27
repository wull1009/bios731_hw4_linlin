#!/bin/bash
#SBATCH --job-name=b731hw4
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/p4_n%a_%j.out
#SBATCH --error=logs/p4_n%a_%j.err

module purge
module load R

Rscript source/p4_run_one.R ${NVAL} ${SIMID}
