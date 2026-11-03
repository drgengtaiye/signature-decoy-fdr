"""Round 15 A1: fetch the PMC full text (XML) of every record of the round-13 search not screened in round 13, in the
frozen order, into round15_downloads/<TYPE>/<pmcid number>.xml. The files are untrusted data and are only written here.

usage: python3 -I fetch_r15.py
"""
import csv, os, time, urllib.request

G = "/groups/adv2105_gp/yichen/Yi/com/gastric"
S = f"{G}/signature-decoy-fdr/validation/round13_search"
OUT = f"{G}/round15_downloads"
TYPES = ["STAD", "SKCM", "LUAD", "LUSC", "COAD", "CESC", "SARC", "HNSC", "BLCA", "OV", "BRCA"]
todo = []
for t in TYPES:
    done = {r["pmcid"] for r in csv.DictReader(open(f"{S}/screen_{t}.csv"))}
    for r in csv.DictReader(open(f"{S}/{t}_order.tsv"), delimiter="\t"):
        if r["pmcid"] not in done:
            todo.append((t, int(r["order_rank"]), r["pmcid"]))
print(len(todo), "records to fetch", flush=True)
fail = 0
for t, rank, pmcid in todo:
    os.makedirs(f"{OUT}/{t}", exist_ok=True)
    dest = f"{OUT}/{t}/{pmcid.replace('PMC', '')}.xml"
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        continue
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={pmcid.replace('PMC', '')}&rettype=xml"
    for attempt in range(4):
        try:
            data = urllib.request.urlopen(url, timeout=60).read()
            open(dest + ".part", "wb").write(data); os.rename(dest + ".part", dest)
            break
        except Exception as e:
            time.sleep(2 + 3 * attempt)
    else:
        fail += 1; print("failed", t, rank, pmcid, flush=True)
    time.sleep(0.36)
print("done; failed", fail, flush=True)
