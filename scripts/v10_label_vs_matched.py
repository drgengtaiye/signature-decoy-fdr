"""label_vs_matched.csv with the signed weighted-mean null (Additional file 1: Note S1, "The matched null against a
permutation of the labels"). Every 12-cohort benchmark row with the per-cohort 95th percentile of the folded AUC under
label permutation (label_p95, unchanged) and whether the signature clears each bar. Run on the v4 table it must
reproduce the v4 file.     python v10_label_vs_matched.py   -> icb/label_vs_matched.csv (v4 copy in _pre_v10/)"""
import os, shutil, numpy as np, pandas as pd
D = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
if not os.path.exists(f"{D}/_pre_v10/label_vs_matched.csv"):
    shutil.copy2(f"{D}/label_vs_matched.csv", f"{D}/_pre_v10/label_vs_matched.csv")
OLD = pd.read_csv(f"{D}/_pre_v10/label_vs_matched.csv")
LAB = OLD.groupby("cohort").label_p95.first()
assert (OLD.groupby("cohort").label_p95.nunique() == 1).all()
def build(tab):
    R = pd.read_csv(tab); R = R[R.cohort != "ICB_VanDenEnde"][list(OLD.columns[:11])].copy()
    R["label_p95"] = R.cohort.map(LAB); R["beats_matched"] = R.AUC > R.rand_p95; R["beats_label"] = R.AUC > R.label_p95
    return R.reset_index(drop=True)
chk = build(f"{D}/_pre_v10/native_response_matched_v2.csv")
assert chk[["cohort", "signature"]].equals(OLD[["cohort", "signature"]])
assert np.allclose(chk.select_dtypes("number"), OLD.select_dtypes("number")) and (chk.beats_matched == OLD.beats_matched).all()
NEW = build(f"{D}/native_response_matched_v2.csv"); NEW.to_csv(f"{D}/label_vs_matched.csv", index=False)
for lab, T in (("v4", OLD), ("v5", NEW)):
    m, l = T.beats_matched, T.beats_label
    only = T[m & ~l]
    print(f"{lab}: {len(T)} tests; matched {m.sum()}, label {l.sum()}, both {(m & l).sum()}; matched only {len(only)} "
          f"at rungs {sorted(only.null_size.unique())}")
