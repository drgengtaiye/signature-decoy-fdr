"""Round 15 D2, secondary: ER status per sample from the GEO sample annotations (read as text). Writes er_status_r15.csv: cohort, sample, er."""
import re, sys
import pandas as pd
G = "/groups/adv2105_gp/yichen/Yi/com/gastric"
KEYS = {"BRCA_GSE20711": "er status", "BRCA_GSE21653": "er ihc", "BRCA_GSE42568": "er_status", "BRCA_GSE96058": "er status"}
DIRS = {"BRCA_GSE96058": f"{G}/external/cohorts", }
rows = []
for co, key in KEYS.items():
    gse = co.split("_")[1]
    d = DIRS.get(co, f"{G}/external/cohorts_r15")
    cl = pd.read_csv(f"{d}/{co}_clin.tsv", sep="\t", dtype=str)
    txt = open(f"{G}/external/raw_r15/{gse}/{gse}_samples_brief.txt", errors="replace").read()
    for block in re.split(r"(?m)^\^SAMPLE = ", txt)[1:]:
        gsm = block.split()[0]
        title = re.search(r"(?m)^!Sample_title = (.*)$", block)
        m = re.search(rf"(?mi)^!Sample_characteristics_ch1 = {re.escape(key)}: *(.*)$", block)
        val = m.group(1).strip() if m else ""
        rows.append(dict(cohort=co, gsm=gsm, title=title.group(1).strip() if title else "", raw=val))
    print(co, "clin samples", len(cl), "first ids", cl["sample"].head(2).tolist())
R = pd.DataFrame(rows)
print(R.groupby(["cohort", "raw"]).size().to_string())
out = []
for co, x in R.groupby("cohort"):
    d = DIRS.get(co, f"{G}/external/cohorts_r15")
    cl = pd.read_csv(f"{d}/{co}_clin.tsv", sep="\t", dtype=str)
    key = "gsm" if cl["sample"].isin(x.gsm).mean() > 0.5 else "title"
    y = x[x[key].isin(cl["sample"])]
    print(co, "matched on", key, len(y), "of", len(cl))
    val = y.raw.str.lower().map(lambda v: 1 if v in ("1", "pos", "positive", "er+") else (0 if v in ("0", "neg", "negative", "er-") else None))
    out.append(pd.DataFrame(dict(cohort=co, sample=y[key], er=val)).dropna())
O = pd.concat(out); O["er"] = O.er.astype(int); O.to_csv("er_status_r15.csv", index=False)
ev = []
for co, x in O.groupby("cohort"):
    d = DIRS.get(co, f"{G}/external/cohorts_r15"); cl = pd.read_csv(f"{d}/{co}_clin.tsv", sep="\t")
    m = cl.merge(x, on="sample"); ev.append(dict(cohort=co, annotated=len(m), er_pos=int((m.er == 1).sum()), events_er_pos=int(m[m.er == 1].event.sum())))
print(pd.DataFrame(ev).to_string(index=False))
