"""Planted-signal replicates summarized as ratios of totals: per scenario,
statistic and truth definition, the realized false discovery proportion sum(V)/sum(R) and each estimator's sum(E)/sum(R).
    python v10_calib_signal_totals.py   -> icb/v10_calib_signal_totals.csv"""
import pandas as pd
D = "/path/to/analysis_root/icb"
S = pd.read_csv(f"{D}/v10_calib_signal.csv")
rows = []
for k, g in S.groupby(["scenario", "beta", "cohorts", "stat", "truth"]):
    r = dict(zip(["scenario", "beta", "cohorts", "stat", "truth"], k), replicates=len(g), mean_calls=g.R.mean(),
             mean_false=g.V.mean(), realized_fdp=g.V.sum() / g.R.sum())
    for e in ("coex A", "random A", "coex B", "coex A, closest"):
        r[f"estimate {e}"] = g[f"E {e}"].sum() / g[f"R {e}"].sum()
    rows.append(r)
pd.DataFrame(rows).to_csv(f"{D}/v10_calib_signal_totals.csv", index=False)
print(pd.DataFrame(rows)[lambda x: x.truth == "10%"].round(3).to_string())
