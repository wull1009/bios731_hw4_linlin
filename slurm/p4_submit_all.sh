#!/bin/bash

for N in 100 1000 10000
do
  for SIMID in $(seq 1 500)
  do
    sbatch --export=ALL,NVAL=${N},SIMID=${SIMID} slurm/p4_job.sh
  done
done
