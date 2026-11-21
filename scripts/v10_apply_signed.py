"""Analysis 1 of validation/PLAN_calibration_round5.md: put the signed nulls of v10_signed_null.R into the three
benchmark tables.

    python v10_apply_signed.py   -> native_response_matched_v2.csv, native_prognosis_matched.csv,
                                    native_icbsurv_matched.csv (originals in icb/_pre_v10/)

For EMT_Thompson, IMPRES_Auslander, IRG_Yang and MPS_Guijarro, rand_med, rand_p95, rand_mean and rand_sd are
replaced by those of their signed null; every other row and every observed statistic is unchanged. A column
null_signed marks the replaced rows, and ladder_med, ladder_p95, ladder_mean and ladder_sd keep the unsigned
rung of every row, which is what the closed-form null (cf2_fit.py) models.
"""
import numpy as np, pandas as pd

D = "/path/to/analysis_root/icb"
COLS = ["rand_med", "rand_p95", "rand_mean", "rand_sd"]
for tab, sig, keys, stat in [("native_response_matched_v2.csv", "response", ["cohort", "signature"], "AUC"),
                             ("native_prognosis_matched.csv", "prognosis", ["cancer", "endpoint", "signature"], "C"),
                             ("native_icbsurv_matched.csv", "icbsurv", ["cohort", "endpoint", "signature"], "C")]:
    T = pd.read_csv(f"{D}/_pre_v10/{tab}")
    S = pd.read_csv(f"{D}/v10_signed_{sig}.csv")
    T = T.set_index(keys); S = S.set_index(keys)
    assert S.index.isin(T.index).all()
    assert np.allclose(T.loc[S.index, stat], S[stat], atol=1e-12), "observed statistic moved"
    assert np.allclose(T.loc[S.index, "rand_med"], S.old_med) and np.allclose(T.loc[S.index, "rand_sd"], S.old_sd)
    T["null_signed"] = False
    for c in COLS:
        T[c.replace("rand_", "ladder_")] = T[c]
    T.loc[S.index, COLS] = S[COLS].to_numpy()
    T.loc[S.index, "null_signed"] = True
    T = T.reset_index()[list(pd.read_csv(f"{D}/_pre_v10/{tab}", nrows=0).columns) + ["null_signed"] + [c.replace("rand_", "ladder_") for c in COLS]]
    T.to_csv(f"{D}/{tab}", index=False)
    print(f"{tab}: {len(S)} of {len(T)} rows given the signed null")
