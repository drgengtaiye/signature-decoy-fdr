"""Analysis 2 of validation/PLAN_calibration_round5.md: why the response benchmark (Fig. 3) calls 20 signatures
and the N0 bank of the target-decoy calibration (Fig. 4) calls 21.

    python v10_null_decomp.py   -> icb/v10_null_decomp.csv, icb/v10_null_decomp_calls.csv

Four constructions of the folded-AUC null, the same observed AUCs throughout, 12 cohorts:
  (i)   ladder rung, unsigned weighted-mean null, 200 draws        native_response_matched_v2.csv (as published)
  (ii)  ladder rung, signed null for the four +1/-1 signatures     v10_signed_response.csv patched in
  (iii) exact size, signed, 1,000 draws                            null_v2/<cohort>_v2.rds, N0 (Fig. 4), via v10_export_n0.R
  (iv)  exact size, signed, 200 draws                              500 random subsamples of (iii)
A call is two-sided Wilcoxon P < 0.05 over cohorts (Fig. 3) and P < 0.05 with a positive median (Fig. 4).
"""
import glob, numpy as np, pandas as pd
from scipy.stats import wilcoxon

D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
EX = "ICB_VanDenEnde"
rng = np.random.default_rng(20261021)

R = pd.read_csv(f"{D}/native_response_matched_v2.csv"); R = R[R.cohort != EX]
S = pd.read_csv(f"{D}/v10_signed_response.csv"); S = S[S.cohort != EX].set_index(["cohort", "signature"])
obs = R.set_index(["cohort", "signature"]).AUC
med_i = R.set_index(["cohort", "signature"]).rand_med
med_ii = med_i.copy(); med_ii.loc[S.index] = S.rand_med

N0 = pd.read_csv(f"{D}/v10_n0_bank.csv.gz")                            # v10_export_n0.R
N0 = N0[N0.cohort != EX].set_index(["cohort", "signature"])
bank = {k: np.maximum(v, 1 - v) for k, v in zip(N0.index, N0.to_numpy(float))}   # folded, as Fig. 4 reads it
idx = obs.index
assert all(k in bank for k in idx), "bank missing a test"
B = np.stack([bank[k] for k in idx])                                   # tests x 1000
med_iii = pd.Series(np.median(B, 1), idx)


def calls(med):
    d = (obs - med).rename("d").reset_index()
    out = {}
    for s_, g in d.groupby("signature"):
        p = wilcoxon(g.d).pvalue; m = np.median(g.d)
        out[s_] = (p, m)
    t = pd.DataFrame(out, index=["p", "med"]).T
    return t, int((t.p < .05).sum()), int(((t.p < .05) & (t.med > 0)).sum())


rows, per = [], {}
for lab, med in [("(i) ladder, unsigned, 200", med_i), ("(ii) ladder, signed, 200", med_ii),
                 ("(iii) exact, signed, 1000", med_iii)]:
    t, c2, c1 = calls(med); per[lab] = t
    rows.append(dict(construction=lab, called_two_sided=c2, called_positive=c1))
    print(f"{lab:28s} two-sided {c2:2d}   positive median {c1:2d}")

NSUB = 500
sub2, sub1, freq = [], [], {}
for b in range(NSUB):
    j = np.stack([rng.choice(1000, 200, replace=False) for _ in range(len(idx))])
    med = pd.Series(np.median(np.take_along_axis(B, j, 1), 1), idx)
    t, c2, c1 = calls(med); sub2.append(c2); sub1.append(c1)
    for s_ in t.index[(t.p < .05)]:
        freq[s_] = freq.get(s_, 0) + 1
sub2, sub1 = np.array(sub2), np.array(sub1)
rows.append(dict(construction="(iv) exact, signed, 200 (500 subsamples)",
                 called_two_sided=float(np.median(sub2)), called_positive=float(np.median(sub1)),
                 two_sided_lo=float(np.quantile(sub2, .025)), two_sided_hi=float(np.quantile(sub2, .975)),
                 share_20=float(np.mean(sub2 == 20)), share_21=float(np.mean(sub2 == 21)),
                 share_ge21=float(np.mean(sub2 >= 21))))
print(f"(iv) 200-draw subsamples: two-sided median {np.median(sub2):.0f} (95% {np.quantile(sub2,.025):.0f}-"
      f"{np.quantile(sub2,.975):.0f}); 20 in {100*np.mean(sub2==20):.0f}%, 21 in {100*np.mean(sub2==21):.0f}%")

C = pd.concat({k: v for k, v in per.items()}, axis=1)
C.columns = [f"{a}|{b}" for a, b in C.columns]
C["call_freq_200draw_subsamples"] = pd.Series(freq).reindex(C.index).fillna(0) / NSUB
diff = C[[c for c in C.columns if c.endswith("|p")]].lt(.05).nunique(axis=1) > 1
print("\nsignatures whose call differs between constructions:")
print(C[diff | (C.call_freq_200draw_subsamples.between(.02, .98))].round(4).to_string())
pd.DataFrame(rows).to_csv(f"{D}/v10_null_decomp.csv", index=False)
C.reset_index(names="signature").to_csv(f"{D}/v10_null_decomp_calls.csv", index=False)
