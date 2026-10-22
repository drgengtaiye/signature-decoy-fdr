#!/bin/bash
# the 12 response cohorts, one process each
cd /groups/adv2105_gp/yichen/Yi/com/gastric/icb
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
printf '%s\n' ICB_Braun ICB_Gide ICB_Hugo ICB_Jung ICB_Kim ICB_Limagne1 ICB_Liu \
  ICB_Mariathasan ICB_Miao1 ICB_Puch ICB_Riaz ICB_Van_Allen \
| xargs -P 12 -I{} sh -c 'Rscript v10_store.R {} 1000 > v10_store/logs/{}.log 2>&1; echo "$? {}" >> v10_store/logs/exit.txt'
echo "ALL DONE $(date)" >> v10_store/logs/exit.txt
