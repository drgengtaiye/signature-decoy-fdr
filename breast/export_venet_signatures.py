"""Export the 48 sets of the Venet et al. 2011 bundle (cancer.signatures) as a long CSV, parsed with the pure-Python
rdata package (never with R).

usage: <python with rdata> -I export_venet_signatures.py <signatures.Rda> <venet_cancer_signatures.csv>
"""
import sys, csv, rdata
d = rdata.read_rda(sys.argv[1])
obj = d["cancer.signatures"]
with open(sys.argv[2], "w", newline="") as f:
    w = csv.writer(f); w.writerow(["signature", "pmid", "entrez", "symbol"])
    for k, x in obj.items():
        ids = [int(v) for v in x["sig"]]; sy = [str(s) for s in x["symb"]]
        pm = str(x["pmid"][0]) if len(x["pmid"]) else ""
        assert len(ids) == len(sy)
        for i, s in zip(ids, sy): w.writerow([str(k), pm, i, s])
print(len(obj))
