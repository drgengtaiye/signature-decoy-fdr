"""Test B of validation/PLAN_three_tests.md: is the library one axis?

    python tb2_test.py   -> icb/tb_results.csv
"""
import functools, warnings; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd

D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
R = pd.read_csv(f"{D}/tb_axis_scores.csv")
dev = pd.read_csv(f"{D}/native_response_matched_v2.csv"); dev["eff"] = dev.AUC - dev.rand_med
best = dev.groupby("signature").eff.mean().idxmax()
print(f"comparator fixed in development: {best} (mean development effect "
      f"{dev.groupby('signature').eff.mean().max():+.4f})")

P = R.drop_duplicates("cohort")[["cohort", "n", "n_sig", "var_pc1", "auc_pc1"]].set_index("cohort")
B = R[R.signature == best].set_index("cohort").auc_sig
d = P.join(B.rename("auc_best")).dropna()
d["diff"] = d.auc_pc1 - d.auc_best
print(f"\ncohorts {len(d)} | PC1 explains {100*d.var_pc1.mean():.0f}% of the variance of the 52 scores on average")

rng = np.random.default_rng(20260923)
bs = np.array([d["diff"].to_numpy()[rng.integers(0, len(d), len(d))].mean() for _ in range(4000)])
lo = np.percentile(bs, 5)
print(f"\nPRIMARY  mean paired difference (axis - {best}) {d['diff'].mean():+.4f}, "
      f"one-sided 95% lower bound {lo:+.4f}")
print(f"  -> {'POSITIVE (non-inferior)' if lo > -0.02 else 'negative'} by the pre-specified margin of -0.02")

print("\nSECONDARY")
rank = R.groupby("cohort").apply(lambda g: (g.auc_sig > g.auc_pc1.iloc[0]).sum() + 1)
print(f"  the axis ranks {rank.median():.0f} of {int(P.n_sig.median())} signatures in the median cohort "
      f"(range {rank.min()}-{rank.max()})")
med = R.groupby("cohort").auc_sig.median()
print(f"  the axis beats the median signature in {int((P.auc_pc1 > med.reindex(P.index)).sum())} of {len(P)} cohorts "
      f"(all seven, where the comparator is missing in three)")
print(f"  the axis against the best available signature in each cohort: mean difference "
      f"{(P.auc_pc1 - R.groupby('cohort').auc_sig.max().reindex(P.index)).mean():+.4f}")
print("\nby cohort")
d["rank_pc1"] = rank.reindex(d.index)
d["median_sig"] = R.groupby("cohort").auc_sig.median().reindex(d.index)
print(d[["n", "auc_pc1", "auc_best", "diff", "median_sig", "rank_pc1"]].round(3).to_string())
d.to_csv(f"{D}/tb_results.csv")
