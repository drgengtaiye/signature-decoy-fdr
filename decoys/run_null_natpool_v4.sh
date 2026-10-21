#!/bin/sh
# Natural-pool N1/N2 banks for the 52 signatures (300 draws), then the decoys read against them.
cd /groups/adv2105_gp/yichen/Yi/com/gastric/icb || exit 1
mkdir -p null_v4nat/logs
: > null_v4nat/logs/exit.txt
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
for c in ICB_Mariathasan ICB_Braun ICB_Liu ICB_Limagne1 ICB_Puch ICB_Van_Allen \
         ICB_Gide ICB_Riaz ICB_Kim ICB_Miao1 ICB_Jung ICB_Hugo ; do echo "$c" ; done |
xargs -P 12 -I{} sh -c 'Rscript null_natpool_v4.R {} 300 > null_v4nat/logs/{}.log 2>&1 && Rscript calib_read_natbank.R {} >> null_v4nat/logs/{}.log 2>&1; echo "$? {}" >> null_v4nat/logs/exit.txt'
echo "ALL DONE $(date)" >> null_v4nat/logs/exit.txt
