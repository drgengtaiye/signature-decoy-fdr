#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --time=2-00:00:00
# usage: sbatch -c <cores> --mem=<mem> -J <name> run_v16.sbatch <app>   (round 11, applications/v16_external.R)
export DFDR_LIB=${DFDR_LIB:-/groups/adv2105_gp/yichen/Yi/com/gastric/icb/v12_pkg/Rlib6}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p /groups/adv2105_gp/yichen/Yi/com/gastric/icb/v16_external/logs
cd /groups/adv2105_gp/yichen/Yi/com/gastric/icb/v16_external
/users/yc4229/.conda/envs/r45/bin/Rscript /groups/adv2105_gp/yichen/Yi/com/gastric/signature-decoy-fdr/applications/v16_external.R "$1" "$SLURM_CPUS_PER_TASK" > logs/v16_$1.log 2>&1
