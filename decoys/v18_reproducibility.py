"""Sections A to E of validation/PLAN_round12_reproducibility_by_type.md: split-half reliability of a cancer type's
family value of rho, its sampling SD, whether the round-11 TCGA-external discordance exceeds sampling noise, what the
round-11 test could detect, and single-signature claims in half cohorts.

    python v18_reproducibility.py   -> icb/v18_split/*.csv

Inputs: icb/v18_split/rho_split_{tcga,external}.csv.gz and *_sizes.csv (applications/v18_split_half.R), the round-10
TCGA values (icb/v15_content/atlas_by_type.csv) and the round-11 external values (icb/v17_external/per_cohort.csv).
"""
import numpy as np, pandas as pd
from scipy.stats import chi2, spearmanr, rankdata

ICB = "/path/to/analysis_root/icb"
OUT = f"{ICB}/v18_split"
FAM = pd.read_csv(f"{ICB}/v5_signature_clusters.csv").set_index("signature").family
INF, CON, OTH = "Inflamed/IFN", "Confounder", "Other"
LAB = {INF: "inflamed", CON: "confounder", OTH: "other"}
TYPE = lambda co: co.split("_")[0]
pd.set_option("display.width", 200)


def load(prefix):
    r = pd.read_csv(f"{OUT}/rho_split_{prefix}.csv.gz")
    r["family"] = r.target.map(FAM)
    sz = pd.read_csv(f"{OUT}/{prefix}_sizes.csv").set_index("cohort")
    return r, sz


T, TS = load("tcga")
E, ES = load("external")

# ------------------------------------------------ check: whole-cohort values reproduce rounds 10 and 11
atlas = pd.read_csv(f"{ICB}/v15_content/atlas_by_type.csv")
atlas = atlas[atlas.setting == "tcga_OS"].set_index("cohort")
pc = pd.read_csv(f"{ICB}/v17_external/per_cohort.csv").set_index("cohort")
chk = []
for name, R, ref in (("tcga", T, atlas), ("external", E, pc)):
    f = R[R.half == "full"].groupby(["cohort", "family"]).rho.mean().unstack()
    for fam, lab in LAB.items():
        d = (f[fam] - ref.loc[f.index, f"rho_{lab}"]).abs().max()
        chk.append({"setting": name, "family": lab, "max_abs_diff": d})
chk = pd.DataFrame(chk); print(chk)
chk.to_csv(f"{OUT}/check_rho_full.csv", index=False)
# AMENDMENT_2_round12: cell by cell; a mismatch is accepted only where the stored folded values of the target and one of
# its co-expressed decoys differ by less than 1e-13 without being equal (a tie that 15-digit storage split)
cellref = {"tcga": (pd.read_csv(f"{ICB}/v15_content/cells_tcga_OS.csv.gz"), f"{ICB}/v12_pkg/tcga_OS_units.csv.gz"),
           "external": (pd.read_csv(f"{ICB}/v17_external/cells_external_survival.csv.gz"),
                        f"{ICB}/v16_external/external_surv_units.csv.gz")}
mism = []
for name, R in (("tcga", T), ("external", E)):
    ref, ufile = cellref[name]
    x = R[R.half == "full"].merge(ref[["target", "cohort", "rho"]], on=["target", "cohort"], suffixes=("", "_ref"))
    assert len(x) == (R.half == "full").sum()
    bad = x[(x.rho - x.rho_ref).abs() > 1e-10]
    if len(bad):
        U = pd.read_csv(ufile)
        for _, b in bad.iterrows():
            u = U[(U.target == b.target) & (U.cohort == b.cohort) & ((U.kind == "target") | (U["class"] == "coex"))]
            tf = u[u.kind == "target"].fold_obs.iloc[0]; g = (u[u.kind == "decoy"].fold_obs - tf).abs()
            ok = bool(((g > 0) & (g < 1e-13)).any())
            mism.append({"setting": name, "target": b.target, "cohort": b.cohort, "rho_round12": b.rho,
                         "rho_stored": b.rho_ref, "split_tie": ok})
MM = pd.DataFrame(mism, columns=["setting", "target", "cohort", "rho_round12", "rho_stored", "split_tie"])
MM.to_csv(f"{OUT}/check_rho_cells_mismatch.csv", index=False); print(MM)
assert MM.split_tie.all(), "a whole-cohort rho differs from rounds 10/11 for a reason other than a split tie"


def halves(R, fam):
    """Family mean of rho per split, half and cohort: arrays A[b, t], B[b, t] and the whole-cohort value."""
    x = R[R.family == fam].groupby(["half", "split", "cohort"]).rho.mean()
    full = x.loc["full"].droplevel("split")
    A = x.loc["A"].unstack("cohort")[full.index]; B = x.loc["B"].unstack("cohort")[full.index]
    return A.values, B.values, full


def rowcorr(A, B):
    a = A - A.mean(1, keepdims=True); b = B - B.mean(1, keepdims=True)
    return (a * b).sum(1) / np.sqrt((a ** 2).sum(1) * (b ** 2).sum(1))


def rowspear(A, B):
    return rowcorr(np.apply_along_axis(rankdata, 1, A), np.apply_along_axis(rankdata, 1, B))


sb = lambda r: 2 * r / (1 + r)

# ------------------------------------------------ A, B: reliability over types and SD of each type's value
rel, tv, keep = [], [], {}
for setting, R, SZ, seed in (("tcga", T, TS, 916100), ("external", E, ES, 916100)):
    for fam, lab in LAB.items():
        A, B, full = halves(R, fam)
        r = np.median(rowcorr(A, B)); rs = np.median(rowspear(A, B))
        rng = np.random.default_rng(seed); nt = A.shape[1]; bt = []
        for _ in range(2000):
            i = rng.integers(0, nt, nt)
            if np.unique(i).size < 3:
                continue
            bt.append(sb(np.median(rowcorr(A[:, i], B[:, i]))))
        sd = (A - B).std(0, ddof=1) / 2
        rv = 1 - np.mean(sd ** 2) / full.var(ddof=1)
        rel.append({"setting": setting, "family": lab, "types": nt, "r_half_median": r, "R_spearman_brown": sb(r),
                    "R_lo": np.nanpercentile(bt, 2.5), "R_hi": np.nanpercentile(bt, 97.5),
                    "spearman_half_median": rs, "R_spearman_brown_from_spearman": sb(rs), "R_variance": rv,
                    "var_types": full.var(ddof=1), "mean_sd_full2": np.mean(sd ** 2)})
        for j, co in enumerate(full.index):
            tv.append({"setting": setting, "cohort": co, "type": TYPE(co), "family": lab, "value": full[co],
                       "sd_full": sd[j], "sd_half": np.concatenate([A[:, j], B[:, j]]).std(ddof=1),
                       "n": SZ.loc[co, "n"], "events": SZ.loc[co, "events"]})
        keep[(setting, lab)] = (A, B, full, sd)
REL = pd.DataFrame(rel); TV = pd.DataFrame(tv)
REL.to_csv(f"{OUT}/reliability.csv", index=False); TV.to_csv(f"{OUT}/type_values.csv", index=False)
print(REL.round(3).to_string())

# ------------------------------------------------ C: discordance against sampling noise (the 11 types of round 11)
tr = []
for lab in ("inflamed", "confounder"):
    t = TV[(TV.setting == "tcga") & (TV.family == lab)].set_index("type")
    e = TV[(TV.setting == "external") & (TV.family == lab)].set_index("type")
    ty = e.index
    D = e.value - t.loc[ty, "value"]; v = e.sd_full ** 2 + t.loc[ty, "sd_full"] ** 2
    w = 1 / v; mu = (w * D).sum() / w.sum(); Q = (w * (D - mu) ** 2).sum(); df = len(ty) - 1
    tr.append({"family": lab, "types": len(ty), "offset": mu, "Q": Q, "df": df, "P": chi2.sf(Q, df),
               "I2": max(0.0, (Q - df) / Q), "rms_D_centered": np.sqrt(np.mean((D - D.mean()) ** 2)),
               "rms_sampling_sd": np.sqrt(v.mean())})
    pd.DataFrame({"type": ty, "external": e.value.values, "tcga": t.loc[ty, "value"].values, "D": D.values,
                  "sd_ext": e.sd_full.values, "sd_tcga": t.loc[ty, "sd_full"].values, "z": ((D - mu) / np.sqrt(v)).values}
                 ).to_csv(f"{OUT}/transport_{lab}_by_type.csv", index=False)
TR = pd.DataFrame(tr); TR.to_csv(f"{OUT}/transport_test.csv", index=False); print(TR.round(4).to_string())

# ------------------------------------------------ D: what the round-11 primary could detect under transportability
t = TV[(TV.setting == "tcga") & (TV.family == "inflamed")].set_index("type")
e = TV[(TV.setting == "external") & (TV.family == "inflamed")].set_index("type")
ty = e.index; tt = t.loc[ty]
Rt = float(np.clip(REL.query("setting == 'tcga' and family == 'inflamed'").R_variance.iloc[0], 0, 1))
m = tt.value.mean(); theta = m + Rt * (tt.value.values - m)
rng = np.random.default_rng(916101); S = 10_000
Tsim = theta + rng.normal(size=(S, len(ty))) * tt.sd_full.values
Esim = theta + rng.normal(size=(S, len(ty))) * e.sd_full.values
rho_sim = rowspear(Tsim, Esim)
rng2 = np.random.default_rng(916102); base = np.arange(1, len(ty) + 1, dtype=float)
perm = np.sort(np.array([spearmanr(rng2.permutation(base), base).statistic for _ in range(100_000)]))
crit = np.percentile(perm, 95)
p_sim = (1 + (perm.size - np.searchsorted(perm, rho_sim - 1e-12, side="left"))) / (1 + perm.size)
SIM = pd.DataFrame([{"types": len(ty), "R_tcga_used": Rt, "spearman_median": np.median(rho_sim),
                     "spearman_lo": np.percentile(rho_sim, 2.5), "spearman_hi": np.percentile(rho_sim, 97.5),
                     "share_at_or_below_obs": np.mean(rho_sim <= -0.28 + 1e-12), "crit_95": crit,
                     "power": np.mean(p_sim < 0.05)}])
SIM.to_csv(f"{OUT}/simulation_round11_power.csv", index=False); print(SIM.round(4).to_string())

# ------------------------------------------------ E1: sampling SD against events
x = TV[TV.family == "inflamed"]
b, a = np.polyfit(np.log(x.events), np.log(x.sd_full), 1)
ev = lambda s: float(np.exp((np.log(s) - a) / b))
E1 = pd.DataFrame([{"family": "inflamed", "cohorts": len(x), "slope": b, "intercept": a,
                    "events_for_sd_0.10": ev(0.10), "events_for_sd_0.05": ev(0.05)}])
for lab in ("confounder", "other"):
    y = TV[TV.family == lab]; bb, aa = np.polyfit(np.log(y.events), np.log(y.sd_full), 1)
    E1.loc[len(E1)] = {"family": lab, "cohorts": len(y), "slope": bb, "intercept": aa,
                       "events_for_sd_0.10": float(np.exp((np.log(0.10) - aa) / bb)),
                       "events_for_sd_0.05": float(np.exp((np.log(0.05) - aa) / bb))}
E1.to_csv(f"{OUT}/sd_against_events.csv", index=False); print(E1.round(3).to_string())

# ------------------------------------------------ E2: single-signature claims (rho of 0.95 or more) across halves
sl = []
for setting, R in (("tcga", T), ("external", E)):
    h = R[R.half != "full"].pivot_table(index=["split", "cohort", "target"], columns="half", values="rho")
    hi = h >= 0.95 - 1e-12
    g = hi.groupby(level="split")
    both = g.apply(lambda z: (z.A & z.B).sum()); na = g.A.sum(); nb = g.B.sum()
    rep = (both / na.replace(0, np.nan) + both / nb.replace(0, np.nan)) / 2
    r = h.groupby(level="split").apply(lambda z: np.corrcoef(z.A, z.B)[0, 1])
    full = R[R.half == "full"]
    sl.append({"setting": setting, "cells": int(h.groupby(level="split").size().iloc[0]),
               "share_cells_ge95_half": float(np.mean(np.r_[h.A.values, h.B.values] >= 0.95 - 1e-12)),
               "claims_per_half_mean": float((na.mean() + nb.mean()) / 2),
               "repeat_rate": float(rep.mean()), "repeat_rate_lo": float(rep.quantile(0.025)),
               "repeat_rate_hi": float(rep.quantile(0.975)), "pearson_cells_median": float(r.median()),
               "share_cells_ge95_full": float((full.rho >= 0.95 - 1e-12).mean())})
SLv = pd.DataFrame(sl); SLv.to_csv(f"{OUT}/signature_level.csv", index=False); print(SLv.round(3).to_string())

print(TV[TV.family == "inflamed"].round(3).to_string())

# ------------------------------------------------ S (AMENDMENT_3, post hoc): bootstrap SDs in place of SD_full
import os
if os.path.exists(f"{OUT}/rho_boot_tcga.csv.gz") and os.path.exists(f"{OUT}/rho_boot_external.csv.gz"):
    sdb = []
    for setting in ("tcga", "external"):
        b = pd.read_csv(f"{OUT}/rho_boot_{setting}.csv.gz"); b["family"] = b.target.map(FAM).map(LAB)
        g = b.groupby(["cohort", "family", "boot"]).rho.mean().groupby(["cohort", "family"])
        sdb.append(pd.DataFrame({"sd_boot": g.std(ddof=1), "mean_boot": g.mean(), "nboot": g.size()})
                   .reset_index().assign(setting=setting))
    TVB = TV.merge(pd.concat(sdb), on=["setting", "cohort", "family"], how="left")
    assert TVB.sd_boot.notna().all() and (TVB.nboot == 200).all()
    TVB.to_csv(f"{OUT}/type_values_boot.csv", index=False)
    relb = []
    for (setting, lab), q in TVB.groupby(["setting", "family"]):
        relb.append({"setting": setting, "family": lab, "types": len(q), "var_types": q.value.var(ddof=1),
                     "mean_sd_boot2": np.mean(q.sd_boot ** 2),
                     "R_variance_boot": 1 - np.mean(q.sd_boot ** 2) / q.value.var(ddof=1),
                     "R_variance_prespecified": 1 - np.mean(q.sd_full ** 2) / q.value.var(ddof=1)})
    RELB = pd.DataFrame(relb); RELB.to_csv(f"{OUT}/reliability_boot.csv", index=False); print(RELB.round(3).to_string())
    trb = []
    for lab in ("inflamed", "confounder"):
        t = TVB[(TVB.setting == "tcga") & (TVB.family == lab)].set_index("type")
        e = TVB[(TVB.setting == "external") & (TVB.family == lab)].set_index("type")
        ty = e.index
        D = e.value - t.loc[ty, "value"]; v = e.sd_boot ** 2 + t.loc[ty, "sd_boot"] ** 2
        w = 1 / v; mu = (w * D).sum() / w.sum(); Q = (w * (D - mu) ** 2).sum(); df = len(ty) - 1
        trb.append({"family": lab, "types": len(ty), "offset": mu, "Q": Q, "df": df, "P": chi2.sf(Q, df),
                    "I2": max(0.0, (Q - df) / Q), "rms_sampling_sd": np.sqrt(v.mean())})
        pd.DataFrame({"type": ty, "external": e.value.values, "tcga": t.loc[ty, "value"].values, "D": D.values,
                      "sd_ext_boot": e.sd_boot.values, "sd_tcga_boot": t.loc[ty, "sd_boot"].values,
                      "z": ((D - mu) / np.sqrt(v)).values}).to_csv(f"{OUT}/transport_{lab}_by_type_boot.csv", index=False)
    TRB = pd.DataFrame(trb); TRB.to_csv(f"{OUT}/transport_test_boot.csv", index=False); print(TRB.round(4).to_string())
    t = TVB[(TVB.setting == "tcga") & (TVB.family == "inflamed")].set_index("type")
    e = TVB[(TVB.setting == "external") & (TVB.family == "inflamed")].set_index("type")
    ty = e.index; tt = t.loc[ty]
    Rtb = float(np.clip(RELB.query("setting == 'tcga' and family == 'inflamed'").R_variance_boot.iloc[0], 0, 1))
    m = tt.value.mean(); theta = m + Rtb * (tt.value.values - m)
    rng = np.random.default_rng(916101)
    Tsim = theta + rng.normal(size=(S, len(ty))) * tt.sd_boot.values
    Esim = theta + rng.normal(size=(S, len(ty))) * e.sd_boot.values
    rho_b = rowspear(Tsim, Esim)
    p_b = (1 + (perm.size - np.searchsorted(perm, rho_b - 1e-12, side="left"))) / (1 + perm.size)
    SIMB = pd.DataFrame([{"types": len(ty), "R_tcga_used": Rtb, "spearman_median": np.median(rho_b),
                          "spearman_lo": np.percentile(rho_b, 2.5), "spearman_hi": np.percentile(rho_b, 97.5),
                          "share_at_or_below_obs": np.mean(rho_b <= -0.28 + 1e-12), "crit_95": crit,
                          "power": np.mean(p_b < 0.05)}])
    SIMB.to_csv(f"{OUT}/simulation_round11_power_boot.csv", index=False); print(SIMB.round(4).to_string())
    e1 = []
    for lab in ("inflamed", "confounder", "other"):
        y = TVB[TVB.family == lab]; bb, aa = np.polyfit(np.log(y.events), np.log(y.sd_boot), 1)
        e1.append({"family": lab, "cohorts": len(y), "slope": bb, "intercept": aa,
                   "events_for_sd_0.10": float(np.exp((np.log(0.10) - aa) / bb)),
                   "events_for_sd_0.05": float(np.exp((np.log(0.05) - aa) / bb))})
    E1B = pd.DataFrame(e1); E1B.to_csv(f"{OUT}/sd_against_events_boot.csv", index=False); print(E1B.round(3).to_string())
    print(TVB[TVB.family == "inflamed"][["setting", "cohort", "value", "sd_full", "sd_half", "sd_boot", "mean_boot",
                                         "events"]].round(3).to_string())
