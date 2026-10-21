#!/bin/sh
# Target-decoy calibration v3: own-bank check on 4 decoys per target in the 12 primary cohorts.
cd /groups/adv2105_gp/yichen/Yi/com/gastric/icb || exit 1
mkdir -p calib_v3/logs
: > calib_v3/logs/exit_own.txt
export OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6
for c in ICB_Mariathasan ICB_Braun ICB_Liu ICB_Limagne1 ICB_Puch ICB_Van_Allen \
         ICB_Gide ICB_Riaz ICB_Kim ICB_Miao1 ICB_Jung ICB_Hugo ; do echo "$c" ; done |
xargs -P 12 -I{} sh -c 'Rscript calib_ownbank_v3.R {} 300 4 > calib_v3/logs/{}_own.log 2>&1; echo "$? {}" >> calib_v3/logs/exit_own.txt'
echo "ALL DONE $(date)" >> calib_v3/logs/exit_own.txt
