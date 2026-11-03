"""Round 15 A2: agreement of the second, independent extraction with the first (round 13 A4's comparison), for the
10% sample of newly eligible models.
usage: python3 -I compare_second_pass_r15.py"""
import csv, re
DASH = "−–—˗‐‑‒﹣－"
def nums(s):
    out = []
    for x in str(s).split(";"):
        x = re.sub(r"\s+", "", x.translate({ord(c): "-" for c in DASH})).strip("[]()").lstrip("+")
        x = re.sub(r"(?<=\d),(?=\d)", "", x); x = re.sub(r"(×|x|\*)10\^?", "e", x)
        try: out.append(float(x))
        except ValueError: out.append(None)
    return out
genes = lambda s: [re.sub(r"\s+", "", g).upper() for g in str(s).split(";") if g.strip()]
first = {(r["type"], r["pmcid"]): r for r in csv.DictReader(open("screen_r15_all.csv"))}
rows = []
for r in csv.DictReader(open("second_pass_r15.csv")):
    a = first[(r["type"], r["pmcid"])]
    g1, g2 = genes(a["genes"]), genes(r["genes"]) if r["found"].upper() == "TRUE" else []
    c1, c2 = nums(a["coefs"]), nums(r["coefs"]) if g2 else []
    d1, d2 = dict(zip(g1, c1)), dict(zip(g2, c2))
    same_set = set(g1) == set(g2)
    same_coefs = same_set and all(d1[g] is not None and d2[g] is not None and abs(d1[g] - d2[g]) < 1e-9 for g in g1)
    rows.append(dict(type=r["type"], pmcid=r["pmcid"], found2=r["found"], n1=len(g1), n2=len(g2),
                     same_gene_order=g1 == g2, same_gene_set=same_set, same_coefs=same_coefs, genes1=";".join(g1),
                     genes2=";".join(g2), coefs1=a["coefs"], coefs2=r["coefs"] if g2 else "", notes2=r["notes"]))
with open("second_pass_agreement_r15.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(len(rows), "sampled;", sum(r["same_gene_set"] for r in rows), "same genes;", sum(r["same_coefs"] for r in rows),
      "same genes and coefficients")
for r in rows:
    if not r["same_coefs"]:
        print("--", r["type"], r["pmcid"]); print("   1:", r["genes1"], "|", r["coefs1"]); print("   2:", r["genes2"], "|", r["coefs2"])
