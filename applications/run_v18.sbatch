#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --time=1-00:00:00
# usage: sbatch -c <cores> --mem=<mem> -J <name> run_v18.sbatch <app>   (round 12, applications/v18_split_half.R)
export DFDR_LIB=${DFDR_LIB:-/path/to/analysis_root/icb/v12_pkg/Rlib6}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p /path/to/analysis_root/icb/v18_split/logs
cd /path/to/analysis_root/icb/v18_split
Rscript /path/to/analysis_root/signature-decoy-fdr/applications/v18_split_half.R "$1" "$SLURM_CPUS_PER_TASK" > logs/v18_$1.log 2>&1
