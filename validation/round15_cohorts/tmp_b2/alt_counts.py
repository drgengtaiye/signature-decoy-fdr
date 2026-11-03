"""Sensitivity counts for the uncertain filters (not used for eligibility)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import screen_b2 as sb  # re-runs the screen; harmless (rewrites the same outputs)
D, Y = sb.D, sb.Y
def c(**k):
    k.setdefault("filters", {}); return sb.count(k)
print("GSE26193 all histologies", c(gse="GSE26193", time="os time (years)", event="os event", event_map={"1": 1, "0": 0}, time_scale=Y))
print("GSE140082 all", c(gse="GSE140082", time="final_ostm", event="final_osid", event_map={"1": 1, "0": 0}, time_scale=D))
print("GSE63885 all", c(gse="GSE63885", time="os - overall survival [days]", event="clinical status at last follow-up (awd - alive with disease; dod - dead from disease; ned - no evidence of disease)", event_map={"dod": 1, "ned": 0, "awd": 0}, time_scale=D))
print("GSE72970 colon+rectosigmoid", c(gse="GSE72970", time="os", event="os censored", event_map={"1": 1, "0": 0}, filters={"tumor location": ["Left colon", "Right colon", "Caecum", "Transverse colon", "Rectum-sigmoid junction"]}))
print("GSE72970 all", c(gse="GSE72970", time="os", event="os censored", event_map={"1": 1, "0": 0}))
print("GSE103479 all", c(gse="GSE103479", time="overall survival time", event="status alive.dead", event_map={"dead": 1, "alive": 0}))
print("GSE14333 all sites", c(gse="GSE14333", composite="Location", time="DFS_Time", event="DFS_Cens", event_map={"0": 1, "1": 0}))
print("GSE31210 authors' exclusion applied", c(gse="GSE31210", time="days before death/censor", event="death", event_map={"dead": 1, "alive": 0}, filters={"tissue": ["primary lung tumor"], "exclude for prognosis analysis due to incomplete resection or adjuvant therapy": ["none"]}, time_scale=D))
print("GSE21050 incl time 0 (raw)", "see zero count")
