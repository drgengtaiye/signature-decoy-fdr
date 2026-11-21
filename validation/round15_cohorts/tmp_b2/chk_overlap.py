import sys, statistics as st, re
sys.path.insert(0, "/path/to/analysis_root/signature-decoy-fdr/validation/round15_cohorts/tmp_b2")
from soft_meta import load
m = load("GSE14333")
for c in ("0", "1"):
    t = []
    for d in m.values():
        kv = {p.split(":",1)[0].strip(): p.split(":",1)[1].strip() for p in ("Location: "+d["Location"]).split(";") if ":" in p}
        if kv["DFS_Cens"] == c: t.append(float(kv["DFS_Time"]))
    print("GSE14333 DFS_Cens", c, len(t), "median", st.median(t))
a = {d["(title)"] for d in load("GSE14333").values()}
b = {d["(title)"]: d.get("corresponding to gse17536") for d in load("GSE29621").values()}
print("GSE14333 x GSE29621 titles", len(a & set(b)))
s17536 = set(load("GSE17536").keys())
print("GSE29621 mapped GSMs in GSE17536:", sum(1 for v in b.values() if v in s17536), "of", len(b))
both = [t for t in a & set(b) if b[t] in s17536]
print("GSE14333 titles that map via GSE29621 to GSE17536:", len(both))
# GSE41271 vs GSE42127 titles
t41 = {re.sub(r"\s+", "", d["(title)"]) for d in load("GSE41271").values()}
t42 = {re.sub(r"^human lung cancer sample ", "", d["(title)"]).replace(" ", "") for d in load("GSE42127").values()}
print("GSE41271 x GSE42127 titles", len(t41 & t42), len(t41), len(t42))
# GSE29013 vs GSE42127/41271 (MDA)
# GSE32894 vs GSE32548
u1 = {d["(source)"].split()[-1] for d in load("GSE32894").values()}
u2 = {d["(source)"].split()[-1] for d in load("GSE32548").values()}
print("GSE32894 x GSE32548 UC ids", len(u1 & u2))
# GSE21050 vs GSE71118 titles
s1 = {d["(title)"].split("_")[0] for d in load("GSE21050").values()}
s2 = {d["(title)"] for d in load("GSE71118").values()}
print("GSE21050 x GSE71118 S-ids", len(s1 & s2), len(s1), len(s2))
# GSE17537 vs GSE38832 : no ids. GSE58812 vs GSE103091
g1 = set(load("GSE58812").keys()); g2 = set(load("GSE103091").keys())
print("GSE58812 GSMs in GSE103091", len(g1 & g2), len(g1))
g1 = set(load("GSE39582").keys()); g2 = set(load("GSE40967").keys())
print("GSE39582 GSMs in GSE40967", len(g1 & g2), len(g1))
g1 = set(load("GSE17536").keys()) | set(load("GSE17537").keys()); g2 = set(load("GSE17538").keys())
print("GSE17536+17537 GSMs in GSE17538", len(g1 & g2), len(g1))
g1 = set(load("GSE84433").keys()) | set(load("GSE84426").keys()); g2 = set(load("GSE84437").keys())
print("GSE84433+84426 GSMs in GSE84437", len(g1 & g2), len(g1))
g1 = set(load("GSE62254").keys()); g2 = set(load("GSE66229").keys())
print("GSE62254 in GSE66229", len(g1 & g2))
g1 = set(load("GSE15459").keys()); g2 = set(load("GSE15460").keys())
print("GSE15459 in GSE15460", len(g1 & g2))
