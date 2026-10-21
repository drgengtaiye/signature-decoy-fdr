#!/bin/bash
# 13 eligible cohorts (the 12 of the main analysis + VanDenEnde as a sensitivity analysis), 1,000 draws each
cd /groups/adv2105_gp/yichen/Yi/com/gastric/icb
export OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6
mkdir -p null_v2/logs
printf '%s\n' ICB_Braun ICB_Gide ICB_Hugo ICB_Jung ICB_Kim ICB_Limagne1 ICB_Liu \
  ICB_Mariathasan ICB_Miao1 ICB_Puch ICB_Riaz ICB_Van_Allen ICB_VanDenEnde \
| xargs -P 13 -I{} sh -c 'Rscript null_matched_v3_coexpr.R {} 1000 > null_v2/logs/{}.log 2>&1; echo "$? {}" >> null_v2/logs/exit.txt'
echo "ALL DONE $(date)" >> null_v2/logs/exit.txt
