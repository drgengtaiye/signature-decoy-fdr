import sys, statistics as st
sys.path.insert(0, "/path/to/analysis_root/signature-decoy-fdr/validation/round15_cohorts/tmp_b2")
from soft_meta import load
from collections import Counter
m = load("GSE72970")
for c in ("0", "1"):
    t = [float(d["os"]) for d in m.values() if d["os censored"] == c]
    print("GSE72970 os censored", c, len(t), "median os", st.median(t))
    t = [float(d["pfs"]) for d in m.values() if d["pfs censored"] == c]
    print("GSE72970 pfs censored", c, len(t), "median pfs", st.median(t))
m = load("GSE14333")
sub = []
for g, d in m.items():
    v = "Location: " + d["Location"]
    kv = dict(x.split(":", 1) for x in [p.strip() for p in v.split(";")] if ":" in x)
    kv = {k.strip(): x.strip() for k, x in kv.items()}
    sub.append(kv)
print(Counter(tuple(sorted(k.keys())) for k in sub))
print(Counter(k["Location"] for k in sub))
print(Counter((k["DukesStage"], k["DFS_Cens"] if k["DFS_Cens"] in ("0","1") else "NA") for k in sub))
print(Counter(k["AdjCTX"] for k in sub))
