import sys
sys.path.insert(0, "/path/to/analysis_root/signature-decoy-fdr/validation/round15_cohorts/tmp_b2")
from soft_meta import load
a = load("GSE11969"); b = load("GSE13213")
ka = {(d.get("Age"), d.get("Sex"), d.get("Survival (days)")) for d in a.values() if d.get("Histology") == "AD"}
kb = {(d.get("Age"), d.get("Sex"), d.get("Survival (days)")) for d in b.values()}
print("GSE11969 AD x GSE13213 (age,sex,days) matches:", len(ka & kb), len(ka), len(kb))
# GSE17537 vs GSE38832 : match on stage + dss? GSE17537 has no dss; skip
# GSE20685 n with time max; GSE21050 time max
for g, f in (("GSE21050", "time"), ("GSE71118", "time"), ("GSE65858", "os"), ("GSE20685", "follow_up_duration (years)"), ("GSE50081", "survival time"), ("GSE103479", "overall survival time"), ("GSE72970", "os"), ("GSE140082", "final_ostm"), ("GSE38832", "dss_time (disease specific survival time, months)"), ("GSE31684", "survival.months"), ("GSE41613", "fu time"), ("GSE32894", "time_to_dod_(months)"), ("GSE26939", "survival_months"), ("GSE30219", "follow-up time (months)"), ("GSE42127", "overall survival months"), ("GSE17536", "overall survival follow-up time"), ("GSE21653", "dfs time (months)"), ("GSE58812", "os (days)"), ("GSE42568", "overall survival time_days")):
    v = []
    for d in load(g).values():
        try: v.append(float(d.get(f, "x").split()[0]))
        except Exception: pass
    v.sort()
    print(g, f, "n", len(v), "min", v[0], "median", v[len(v)//2], "max", v[-1])
