"""Round 15 A1 recall check: in every record ruled ineligible because no usable formula was found (coefficients not
given in text, only in a figure image, no prognostic model, or not a linear formula), search the whole full text,
MathML and TeX markup included, for printed products of a coefficient and a gene symbol. Records with two or more
such products are listed for a second look. The XML is untrusted data: parsed, never executed.

usage: python3 -I recall_scan_r15.py
"""
import csv, re
import xml.etree.ElementTree as ET

G = "/groups/adv2105_gp/yichen/Yi/com/gastric"
GENE = r"[A-Z][A-Z0-9]{1,}(?:-[A-Z0-9]+)?"
COEF = r"[-−–—‐]?\s?\d*\.\d+(?:\s?[eE]\s?[-−–‐]?\d+|\s?[×x]\s?10\s?[-−–‐]?\d+)?"
PROD = re.compile(COEF + r"\s*[×x\*·]\s*\(?\s*(?:Expr?\w*\s*)?\(?\s*(" + GENE + ")|(" + GENE + r")\s*(?:exp\w*)?\s*[×x\*·]\s*\(?\s*" + COEF)
STOP = {"HR", "CI", "OS", "AUC", "TCGA", "GEO", "ROC", "DFS", "PFS", "RFS", "SD", "IQR", "LASSO", "COX", "ESTIMATE", "TPM",
        "FPKM", "RNA", "DNA", "ICGC", "GSE", "NA", "SE", "OR", "BMI"}
REASONS = ("coefficients not given in text", "coefficients only in figure image", "no prognostic model",
           "not a linear formula", "coefficients only in supplementary table")
rows = list(csv.DictReader(open("screen_r15_all.csv")))
out = []
for r in rows:
    if r["eligible"].upper() == "TRUE" or not r["reason"].lower().startswith(REASONS):
        continue
    p = f"{G}/round15_downloads/{r['type']}/{r['pmcid'].replace('PMC', '')}.xml"
    s = " ".join(ET.parse(p).getroot().itertext())
    s = re.sub(r"\\(rm|mathrm|text|mathit)\s*\{([^{}]*)\}", r"\2", s)        # TeX: {\rm{0}}.{\rm{24}} -> 0.24
    s = re.sub(r"[{}]", "", s)
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"(?<=\d) (?=[\d.])|(?<=\.) (?=\d)", "", s)                     # spaced digits
    genes = set()
    for m in PROD.finditer(s):
        g = m.group(1) or m.group(2)
        if g and g not in STOP and not g.startswith(("GSE", "TCGA")):
            genes.add(g)
    if len(genes) >= 2:
        out.append(dict(type=r["type"], pmcid=r["pmcid"], chunk=r["chunk"], reason=r["reason"][:60],
                        n_products=len(genes), genes=";".join(sorted(genes))[:200]))
with open("recall_scan_r15.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["type", "pmcid", "chunk", "reason", "n_products", "genes"]); w.writeheader(); w.writerows(out)
n = sum(1 for r in rows if r["eligible"].upper() != "TRUE" and r["reason"].lower().startswith(REASONS))
print(f"{n} records with no usable formula scanned; {len(out)} flagged")
for o in out:
    print(o["chunk"], o["pmcid"], o["reason"], o["n_products"], o["genes"][:100])
