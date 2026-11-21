#!/bin/bash
# Download each series' samples with their data tables, and each platform table, through the GEO query page (the GEO
# FTP server is not reachable from this cluster), gzipped, each series in its own directory.
R=/path/to/analysis_root/external/raw
get() {  # $1 accession, $2 dest
  [ -s "$2" ] && return 0
  for try in 1 2 3; do
    curl -s --max-time 7200 "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$1&targ=$3&form=text&view=data" | gzip > "$2.part" \
      && [ $(zcat "$2.part" | head -c 100000 | wc -c) -gt 1000 ] && mv "$2.part" "$2" && return 0
    sleep 30
  done
  echo "FAILED $1" >&2
}
while read gse gpls; do
  mkdir -p $R/$gse
  get $gse $R/$gse/${gse}_gsm_data.soft.gz gsm
  for g in ${gpls//,/ }; do get $g $R/$gse/${g}_table.soft.gz self; done
  echo "$gse done $(date -u +%H:%M:%S) $(du -sh $R/$gse | cut -f1)"
done
