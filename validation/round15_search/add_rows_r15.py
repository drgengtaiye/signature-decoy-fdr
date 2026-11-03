"""Append screening rows to a round-15 CSV. Rows come as JSON objects, one per line, on stdin; the columns are fixed,
a missing key is written empty, an unknown key is an error, and a PMCID already in the file is refused.

usage: python3 -I add_rows_r15.py <out.csv> [--backfill] < rows.jsonl
"""
import csv, json, os, sys

COLS = ["type", "order_rank", "pmcid", "pmid", "year", "journal", "title", "eligible", "reason", "theme", "endpoint",
        "training", "validation_cohorts", "n_genes", "genes", "coefs", "formula_text", "source_location", "notes",
        "trains_tcga_model", "validation_used"]
BACK = ["type", "order_rank", "pmcid", "trains_tcga_model", "validation_used", "notes"]
out = sys.argv[1]
cols = BACK if "--backfill" in sys.argv else COLS
have = set()
if os.path.exists(out):
    have = {r["pmcid"] for r in csv.DictReader(open(out))}
rows = [json.loads(l) for l in sys.stdin if l.strip()]
for r in rows:
    bad = set(r) - set(cols)
    if bad:
        sys.exit(f"unknown columns {sorted(bad)} in {r.get('pmcid')}")
    if r.get("pmcid") in have:
        sys.exit(f"{r.get('pmcid')} already recorded; nothing written")
    for k in ("eligible", "trains_tcga_model"):
        if k in r and str(r[k]).upper() not in ("TRUE", "FALSE"):
            sys.exit(f"{k} must be TRUE or FALSE in {r.get('pmcid')}")
new = not os.path.exists(out)
with open(out, "a", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    if new:
        w.writeheader()
    for r in rows:
        w.writerow({k: r.get(k, "") for k in cols}); have.add(r["pmcid"])
print(f"{len(rows)} rows written; {len(have)} in {out}")
