#!/bin/bash
#SBATCH --job-name=p4post
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH -o logs/p4_post_%j.out
#SBATCH -e logs/p4_post_%j.err

module load R/4.4.0
cd /projects/wrobel/users/linlin/bios731/bios731_hw4_linlin

Rscript source/p4_aggregate_results.R
Rscript source/p4_make_plots.R
