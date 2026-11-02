import json, glob, re, sys, csv
d = sys.argv[1]; rows = []
for f in sorted(glob.glob(d + "/chunk_*.json")):
    r = json.load(open(f))["result"]
    for k in r.get("uids", []):
        x = r[k]
        rows.append(dict(gse=x.get("accession"), n=int(x.get("n_samples") or 0), gpl=x.get("gpl"), date=x.get("pdat"),
                         taxon=x.get("taxon"), gdstype=x.get("gdstype"), title=x.get("title", ""), summary=x.get("summary", "")))
drug = re.compile(r"pembrolizumab|nivolumab|atezolizumab|durvalumab|avelumab|ipilimumab|tremelimumab|camrelizumab|sintilimab|"
                  r"tislelizumab|toripalimab|cemiplimab|spartalizumab|anti-?PD-?L?1|anti-?CTLA-?4|PD-?L?1 (blockade|inhibit)|"
                  r"checkpoint (blockade|inhibit)|\bICIs?\b|\bICB\b|\bICT\b", re.I)
resp = re.compile(r"respon|RECIST|benefit|resist|outcome|predict", re.I)
neg = re.compile(r"single[- ]cell|scRNA|single cell|mouse|murine|\bmice\b|cell line|organoid|PBMC|peripheral blood|plasma|serum|"
                 r"T cell clone|TCR|CAR[- ]T|in vitro|xenograft|methylation|ChIP|ATAC|spatial", re.I)
keep = []
for x in rows:
    txt = x["title"] + " " + x["summary"]
    if "Homo sapiens" not in (x["taxon"] or ""): continue
    if x["n"] < 25: continue
    if not drug.search(txt) or not resp.search(txt): continue
    x["neg"] = ";".join(sorted(set(m.group(0).lower() for m in neg.finditer(txt))))
    keep.append(x)
w = csv.DictWriter(open(d + "/shortlist.csv", "w"), fieldnames=list(keep[0].keys())); w.writeheader(); w.writerows(keep)
print(len(rows), "series;", len(keep), "shortlisted;", sum(1 for x in keep if not x["neg"]), "without exclusion words")
