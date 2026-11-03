"""Round 15 A1-A2: merge the readers' screening files and check them by script (round 13 step 4, done centrally).

For every chunk: each record of the chunk has exactly one row. For every eligible model: genes and coefficients are
lists of the same length as n_genes (2-50), every coefficient parses as a number after the round-13 normalization
(amendment 2 rule 3: minus-like characters read as a minus, brackets removed), and every gene and every coefficient
can be found in the verbatim formula_text. validation_cohorts is filled from the accessions the script found in the
full text (as round 13 recorded them). Writes screen_r15_all.csv and check_r15_report.txt.

usage: python3 -I check_r15.py
"""
import csv, glob, json, re

MINUS = "−–—˗‐‑‒－"
COHORT_ACC = re.compile(r"^(GSE\d+|E-(MTAB|GEOD|TABM)-\d+|METABRIC|CGGA|PRJ(NA|EB|DB)\d+|EGA[SD]\d+|phs\d+|ICGC|"
                        r"[A-Z]{2,4}-(JP|AU|CA|FR|CN|KR|UK|EU|DE|SG|IN|IT|ES|BR|MX|SA|RU))$")


def norm(s):
    s = str(s)
    for c in MINUS:
        s = s.replace(c, "-")
    return re.sub(r"\s+", "", s.replace("(", "").replace(")", "").replace("[", "").replace("]", ""))


def num(s):
    s = norm(s).replace("+", "")
    m = re.fullmatch(r"(-?\d*\.?\d+)(?:[eE](-?\d+)|[×x\*]10\^?(-?\d+))?", s)
    if not m:
        return None
    v = float(m.group(1))
    e = m.group(2) or m.group(3)
    return v * 10 ** int(e) if e else v


snip = {}
for f in glob.glob("snippets/c[01][0-9]_*.jsonl"):
    if "backfill" in f:
        continue
    for l in open(f):
        d = json.loads(l)
        snip[(d["type"], d["pmcid"])] = (f.split("/")[-1][:-6], d)
rows, report = [], []
for f in sorted(glob.glob("screen_r15_c*.csv")):
    chunk = f[len("screen_r15_"):-4]
    rr = list(csv.DictReader(open(f)))
    want = {k for k, (c, _) in snip.items() if c == chunk}
    got = [(r["type"], r["pmcid"]) for r in rr]
    if len(set(got)) != len(got):
        report.append(f"{chunk}: repeated rows {sorted({g for g in got if got.count(g) > 1})}")
    if set(got) != want:
        report.append(f"{chunk}: {len(want - set(got))} records without a row, {len(set(got) - want)} rows not in the chunk")
    for r in rr:
        r["chunk"] = chunk
        d = snip.get((r["type"], r["pmcid"]), (None, {}))[1]
        r["validation_cohorts"] = ";".join(a for a in d.get("accessions_any", []) if COHORT_ACC.match(a))
        rows.append(r)
elig = [r for r in rows if r["eligible"].upper() == "TRUE"]
problems = []
for r in elig:
    g = [x.strip() for x in r["genes"].split(";") if x.strip()]
    c = [x.strip() for x in r["coefs"].split(";") if x.strip()]
    ft = norm(r["formula_text"]).upper()
    msg = []
    if len(g) != len(c):
        msg.append(f"{len(g)} genes vs {len(c)} coefficients")
    if str(len(g)) != str(r["n_genes"]).strip():
        msg.append(f"n_genes {r['n_genes']} vs {len(g)} listed")
    if not 2 <= len(g) <= 50:
        msg.append(f"{len(g)} genes")
    bad = [x for x in c if num(x) is None]
    if bad:
        msg.append(f"unparsed coefficients {bad}")
    miss_g = [x for x in g if norm(x).upper() not in ft]
    if miss_g:
        msg.append(f"genes not in formula_text {miss_g}")
    miss_c = [x for x in c if num(x) is not None and norm(x).lstrip("+-").lstrip("0") not in ft.replace("+", "")]
    if miss_c:
        msg.append(f"coefficients not in formula_text {miss_c}")
    if len(set(g)) != len(g):
        msg.append("repeated gene")
    if msg:
        problems.append(f"{r['chunk']} {r['type']} {r['pmcid']}: " + "; ".join(msg))
pm = {}
for r in rows:
    pm.setdefault(r["pmcid"], []).append(r["type"])
dup = {k: v for k, v in pm.items() if len(v) > 1}
cols = list(rows[0].keys()) if rows else []
with open("screen_r15_all.csv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=cols); w.writeheader(); w.writerows(rows)
by = {}
for r in rows:
    b = by.setdefault(r["type"], [0, 0, 0]); b[0] += 1; b[1] += r["eligible"].upper() == "TRUE"
    b[2] += r["trains_tcga_model"].upper() == "TRUE"
out = [f"rows {len(rows)}; eligible {len(elig)}; PMCIDs under two or more types: {dup}"]
out += [f"  {t}: screened {b[0]}, eligible {b[1]}, trains a TCGA model {b[2]}" for t, b in sorted(by.items())]
out += ["coverage:"] + report + [f"eligible models with a problem: {len(problems)}"] + problems
open("check_r15_report.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
