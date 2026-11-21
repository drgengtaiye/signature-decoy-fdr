#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --time=1-00:00:00
# usage: sbatch -c <cores> --mem=<mem> -J <name> run_v19.sbatch <script> <TYPE> [stage]   (round 13, applications/v19_*.R)
export DFDR_LIB=${DFDR_LIB:-/path/to/analysis_root/icb/v12_pkg/Rlib6}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p /path/to/analysis_root/icb/v19_audit/logs
cd /path/to/analysis_root/icb/v19_audit
Rscript /path/to/analysis_root/signature-decoy-fdr/applications/$1.R "$2" "$SLURM_CPUS_PER_TASK" "${3:-all}" > logs/$1_$2_${3:-all}.log 2>&1
