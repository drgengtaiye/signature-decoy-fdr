import sys, re
sys.path.insert(0, "/path/to/analysis_root/signature-decoy-fdr/validation/round15_cohorts/tmp_b2")
from soft_meta import load
from collections import Counter
m = load("GSE96058")
t = [d["(title)"] for d in m.values()]
print(Counter(re.sub(r"\d+", "#", x) for x in t).most_common(10))
print(Counter(d["(platform)"] for d in m.values() if "repl" in d["(title)"].lower()))
ids = Counter(d["scan-b external id"].split(".")[0] for d in m.values())
print("distinct Q ids", len(ids), Counter(ids.values()))
# GSE20711
m = load("GSE20711")
print(Counter((d["(source)"], d.get("e.os")) for d in m.values()))
# GSE65858 / GSE13507 / GSE32062 / GSE21050 / GSE65904 time ranges
