"""Plan A4: agreement of the second, independent extraction with the first, for the sampled models.
usage: python3 -I compare_second_pass.py <second_pass.csv> <round13_search dir> <out.csv>"""
import sys, re, pandas as pd
SP, DIR, OUT = sys.argv[1:4]
DASH = "−–—˗‐‑‒﹣－"
def nums(s):
    out = []
    for x in str(s).split(";"):
        x = re.sub(r"\s+", "", x.translate({ord(c): "-" for c in DASH})).strip("[]()")
        try: out.append(float(x))
        except ValueError: out.append(None)
    return out
genes = lambda s: [re.sub(r"\s+", "", g).upper() for g in str(s).split(";") if g.strip()]
sp = pd.read_csv(SP); rows = []
for r in sp.itertuples():
    f = pd.read_csv(f"{DIR}/screen_{r.type}.csv"); a = f[f.pmcid == r.pmcid].iloc[0]
    g1, g2 = genes(a.genes), genes(r.genes) if str(r.found).upper() == "TRUE" else []
    c1, c2 = nums(a.coefs), nums(r.coefs) if g2 else []
    same_genes = g1 == g2; same_set = set(g1) == set(g2)
    d1 = dict(zip(g1, c1)); d2 = dict(zip(g2, c2))
    same_coefs = same_set and all(d1[g] is not None and d2[g] is not None and abs(d1[g] - d2[g]) < 1e-9 for g in g1)
    rows.append(dict(type=r.type, pmcid=r.pmcid, found2=r.found, n1=len(g1), n2=len(g2), same_gene_order=same_genes,
                     same_gene_set=same_set, same_coefs=same_coefs, genes1=";".join(g1), genes2=";".join(g2),
                     coefs1=a.coefs, coefs2=r.coefs if g2 else "", notes2=r.notes))
R = pd.DataFrame(rows); R.to_csv(OUT, index=False)
print(len(R), "sampled;", int(R.same_gene_set.sum()), "same genes;", int(R.same_coefs.sum()), "same genes and coefficients")
print(R[~R.same_coefs][["type", "pmcid", "genes1", "genes2", "coefs1", "coefs2"]].to_string())
