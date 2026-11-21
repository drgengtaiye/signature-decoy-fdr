"""Round 15 B1: per type, the validation cohorts named (validation_used) by the records that train a TCGA model, from
round 15's screen (screen_r15_final.csv) and the round-13 backfill (r13_backfill.csv); candidates are those named by at
least 5 records, or the 5 most named if fewer reach 5 (ties included), plus the round-11 cohort. TCGA is excluded.

usage: python3 -I candidates_r15.py
"""
import csv, re
from collections import Counter

S = "/path/to/analysis_root/signature-decoy-fdr/validation/round15_search"
TYPES = ["STAD", "SKCM", "LUAD", "LUSC", "COAD", "CESC", "SARC", "HNSC", "BLCA", "OV", "BRCA"]
R11 = dict(STAD="GSE84437", SKCM="GSE65904", LUAD="GSE72094", LUSC="GSE157011", COAD="GSE39582", CESC="GSE44001",
           SARC="GSE21050", HNSC="GSE65858", BLCA="GSE13507", OV="GSE32062", BRCA="GSE96058")


def tokens(s):
    out = set()
    for x in str(s).split(";"):
        x = x.strip()
        if not x or x.upper().startswith("TCGA"):
            continue
        m = re.match(r"^(GSE\d+|E-[A-Z]{4}-\d+)", x.upper().replace(" ", ""))
        out.add(m.group(1) if m else x)
    return out


rows = [r for r in csv.DictReader(open(f"{S}/screen_r15_final.csv"))] + \
       [r for r in csv.DictReader(open(f"{S}/r13_backfill.csv"))]
cnt = {t: Counter() for t in TYPES}
nrec = Counter()
for r in rows:
    if str(r.get("trains_tcga_model", "")).upper() != "TRUE":
        continue
    nrec[r["type"]] += 1
    for a in tokens(r.get("validation_used", "")):
        cnt[r["type"]][a] += 1
out = []
for t in TYPES:
    c = cnt[t]
    keep = {a for a, n in c.items() if n >= 5}
    if len(keep) < 5:
        ranked = c.most_common()
        if ranked:
            cut = ranked[min(4, len(ranked) - 1)][1]
            keep |= {a for a, n in ranked if n >= cut}
    keep.add(R11[t])
    for a in sorted(keep, key=lambda a: (-c[a], a)):
        out.append(dict(type=t, cohort=a, records_naming=c[a], records_training_tcga_model=nrec[t],
                        round11=a == R11[t]))
with open("candidates_r15.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
for t in TYPES:
    print(t, nrec[t], [(o["cohort"], o["records_naming"]) for o in out if o["type"] == t])
print(len(out), "candidates")
