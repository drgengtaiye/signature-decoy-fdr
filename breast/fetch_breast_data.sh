#!/bin/bash
# Download public breast-cancer cohorts (GEO series matrices + GPL annotations, GSE96058, METABRIC)
# Each download goes into its own directory under $BASE; files already present are not re-fetched.
# Nothing downloaded is ever executed.  Usage: fetch_breast_data.sh [BASE]
set -u
BASE=${1:-/path/to/analysis_root/breast}
GEO=$BASE/download_geo
MB=$BASE/download_metabric
CURL="/usr/bin/curl -sS -L --fail --retry 4 --retry-delay 20 --max-time 3600"   # system curl: the conda curl's TLS handshake to ftp.ncbi times out
FTP=https://ftp.ncbi.nlm.nih.gov/geo

get () {  # $1 dest dir, $2 url
  local d=$1 u=$2 f
  f=$d/$(basename "${u%%\?*}")
  if [ -s "$f" ]; then echo "have $f"; return 0; fi
  for a in 1 2 3; do
    $CURL -o "$f.part" "$u" && [ -s "$f.part" ] && mv "$f.part" "$f" && echo "OK $f" && return 0
    sleep 30
  done
  rm -f "$f.part"; echo "FAILED $u"; return 1
}
stub () { echo "$1" | sed -E 's/[0-9]{1,3}$/nnn/'; }
gpl_annot () { local s; s=$(stub "$1"); echo "$FTP/platforms/$s/$1/annot/$1.annot.gz"; }

# GSE : platform : matrix file(s) : extra supplementary text files
declare -A PLAT=( [GSE2034]=GPL96 [GSE7390]=GPL96 [GSE11121]=GPL96 [GSE2990]=GPL96 [GSE1456]=GPL96
                  [GSE3494]=GPL96 [GSE4922]=GPL96 [GSE20685]=GPL570 [GSE21653]=GPL570 [GSE42568]=GPL570
                  [GSE17705]=GPL96 [GSE12093]=GPL96 [GSE6532]=GPL96 )
declare -A MAT=( [GSE1456]=GSE1456-GPL96_series_matrix.txt.gz [GSE3494]=GSE3494-GPL96_series_matrix.txt.gz
                 [GSE4922]=GSE4922-GPL96_series_matrix.txt.gz
                 [GSE6532]="GSE6532-GPL96_series_matrix.txt.gz GSE6532-GPL570_series_matrix.txt.gz" )
declare -A SUP=( [GSE7390]=GSE7390_transbig2006affy_demo.txt.gz [GSE2990]=GSE2990_suppl_info.txt
                 [GSE4922]=GSE4922_Clinical_file_for_both_Uppsala_Singapore_Samples.txt
                 [GSE6532]="GSE6532_LUMINAL_demo.txt.gz GSE6532_LUMINAL_README.txt.gz" )

for g in GSE2034 GSE7390 GSE11121 GSE2990 GSE1456 GSE3494 GSE4922 GSE20685 GSE21653 GSE42568 GSE17705 GSE12093 GSE6532; do
  d=$GEO/$g; mkdir -p "$d"; s=$(stub $g)
  for m in ${MAT[$g]:-${g}_series_matrix.txt.gz}; do get "$d" "$FTP/series/$s/$g/matrix/$m"; done
  for x in ${SUP[$g]:-}; do get "$d" "$FTP/series/$s/$g/suppl/$x"; done
  get "$d" "$(gpl_annot ${PLAT[$g]})"
done

# GSE96058 (SCAN-B): two platform series matrices (characteristics) + gene-level expression CSV
d=$GEO/GSE96058; mkdir -p "$d"
get "$d" "$FTP/series/GSE96nnn/GSE96058/matrix/GSE96058-GPL11154_series_matrix.txt.gz"
get "$d" "$FTP/series/GSE96nnn/GSE96058/matrix/GSE96058-GPL18573_series_matrix.txt.gz"
get "$d" "$FTP/series/GSE96nnn/GSE96058/suppl/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"

# METABRIC from the cBioPortal datahub (GitHub LFS media; the S3 tarball returns 403)
mkdir -p "$MB"
LFS=https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_metabric
RAW=https://raw.githubusercontent.com/cBioPortal/datahub/master/public/brca_metabric   # small non-LFS files
for f in data_clinical_patient.txt data_clinical_sample.txt data_mrna_illumina_microarray.txt; do
  get "$MB" "$LFS/$f"
done
for f in meta_mrna_illumina_microarray.txt meta_clinical_patient.txt meta_study.txt; do
  get "$MB" "$RAW/$f"
done
echo DONE
