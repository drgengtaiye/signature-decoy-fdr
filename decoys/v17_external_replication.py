"""Section C of validation/PLAN_round11_external_replication.md: rho (a signature's place among its co-expressed decoys)
in the independent survival cohorts and in GSE176307, against the TCGA values of round 10.

    python v17_external_replication.py   -> icb/v17_external/*.csv

Inputs: icb/v16_external/external_surv_units.csv.gz, external_resp_units.csv.gz (applications/v16_external.R) and
icb/v15_content/atlas_by_type.csv (round 10).
"""
import os, numpy as np, pandas as pd
from scipy.stats import wilcoxon, spearmanr

ICB = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
V16, OUT = f"{ICB}/v16_external", f"{ICB}/v17_external"
os.makedirs(OUT, exist_ok=True)
FAM = pd.read_csv(f"{ICB}/v5_signature_clusters.csv").set_index("signature")
CL = FAM.cluster_fixed
INF, CON = "Inflamed/IFN", "Confounder"
FAMS = [INF, "Other", CON]


def cells(path):
    """Per target x cohort: Delta, structure S, content G, rho and clearance (as decoys/v15_content_structure.py)."""
    d = pd.read_csv(path)
    d["delta"] = d.fold_obs - d.fold_med
    d["clear"] = d.fold_obs > d.fold_p95
    t = d[d.kind == "target"].set_index(["target", "cohort"])[["n", "delta", "clear"]]
    out = t.copy()
    for cl in ["coex", "rand"]:
        m = d[(d.kind == "decoy") & (d["class"] == cl)].merge(t.delta.rename("tdelta"), left_on=["target", "cohort"],
                                                               right_index=True)
        m["below"] = (m.delta < m.tdelta) + 0.5 * (m.delta == m.tdelta)
        g = m.groupby(["target", "cohort"]).agg(S=("delta", "median"), rho=("below", "mean"), dclear=("clear", "mean"),
                                                nd=("delta", "size"))
        out = out.join(g.add_suffix("" if cl == "coex" else "_rand"))
    out["G"] = out.delta - out.S
    out = out.reset_index()
    out["family"] = out.target.map(FAM.family)
    return out


def boot_ci(fn, seed, B=10_000):
    """Percentile interval of fn(list of signatures) over resamples of the 13 correlation clusters."""
    r = np.random.default_rng(seed); cl = CL.unique()
    mem = {k: list(CL.index[CL == k]) for k in cl}
    v = []
    for _ in range(B):
        x = fn([s for k in r.choice(cl, len(cl), replace=True) for s in mem[k]])
        if np.isfinite(x):
            v.append(x)
    return np.percentile(v, [2.5, 97.5])


def fam_mean(c, cohort, fam):
    x = c[(c.cohort == cohort) & (c.family == fam)].set_index("target").rho
    def fn(sigs):
        y = x.reindex(sigs).dropna()
        return y.mean() if len(y) else np.nan
    return x.mean(), fn


S = cells(f"{V16}/external_surv_units.csv.gz")
S.to_csv(f"{OUT}/cells_external_survival.csv.gz", index=False)
TYPE = lambda co: co.split("_")[0]
TC = pd.read_csv(f"{ICB}/v15_content/atlas_by_type.csv")
TC = TC[TC.setting == "tcga_OS"].set_index("cohort")

# ------------------------------------------------ per cohort: family means with cluster intervals, single-cohort FDR
rows = []
for co in sorted(S.cohort.unique()):
    x = S[S.cohort == co]
    R = int(x.clear.sum()); ef = x.dclear.sum()
    r = {"cohort": co, "type": TYPE(co), "n": int(x.n.iloc[0]), "calls": R, "exp_false": ef,
         "fdr": 100 * ef / R if R else np.nan, "decoy_clear": x.dclear.mean(), "rho_all": x.rho.mean(),
         "tcga_rho_inflamed": TC.loc[TYPE(co), "rho_inflamed"], "tcga_rho_confounder": TC.loc[TYPE(co), "rho_confounder"]}
    for f, lab in ((INF, "inflamed"), (CON, "confounder"), ("Other", "other")):
        m, fn = fam_mean(S, co, f)
        lo, hi = boot_ci(fn, 916001)
        r.update({f"rho_{lab}": m, f"rho_{lab}_lo": lo, f"rho_{lab}_hi": hi,
                  f"calls_{lab}": int(x[x.family == f].clear.sum())})
    rows.append(r)
PC = pd.DataFrame(rows)
PC.to_csv(f"{OUT}/per_cohort.csv", index=False)

# ------------------------------------------------ 1. primary: Spearman over types, one-sided permutation P
ty = PC.groupby("type").agg(ext=("rho_inflamed", "mean"), tcga=("tcga_rho_inflamed", "first"),
                            cohorts=("cohort", "size")).reset_index()
obs = spearmanr(ty.ext, ty.tcga).statistic
rng = np.random.default_rng(916000); B = 100_000
perm = np.array([spearmanr(rng.permutation(ty.ext.values), ty.tcga).statistic for _ in range(B)])
p1 = (1 + np.sum(perm >= obs - 1e-12)) / (1 + B)
ty.to_csv(f"{OUT}/by_type.csv", index=False)
prim = {"types": len(ty), "spearman": obs, "p_one_sided": p1,
        "pearson": np.corrcoef(ty.ext, ty.tcga)[0, 1]}

# ------------------------------------------------ 2. stomach and melanoma
pred = []
for co, side in [(c, "below") for c in PC.cohort if c.startswith("STAD")] + \
                [(c, "above") for c in PC.cohort if c.startswith("SKCM")]:
    q = PC.set_index("cohort").loc[co]
    est, lo, hi = q.rho_inflamed, q.rho_inflamed_lo, q.rho_inflamed_hi
    holds = est < 0.5 if side == "below" else est > 0.5
    conf = hi < 0.5 if side == "below" else lo > 0.5
    pred.append({"test": f"{co} inflamed rho {side} 0.5", "estimate": est, "lo": lo, "hi": hi, "holds": holds,
                 "confirmed": conf})

# ------------------------------------------------ 3. response: GSE176307, and bladder response against survival
Rr = cells(f"{V16}/external_resp_units.csv.gz")
Rr.to_csv(f"{OUT}/cells_external_response.csv.gz", index=False)
m, fn = fam_mean(Rr, "GSE176307", INF); lo, hi = boot_ci(fn, 916001)
pred.append({"test": "GSE176307 inflamed rho above 0.5 (response)", "estimate": m, "lo": lo, "hi": hi,
             "holds": m > 0.5, "confirmed": lo > 0.5})
xr = Rr[(Rr.cohort == "GSE176307") & (Rr.family == INF)].set_index("target").rho
bl = [c for c in PC.cohort if c.startswith("BLCA")][0]
xs = S[(S.cohort == bl) & (S.family == INF)].set_index("target").rho
dd = (xr - xs.reindex(xr.index)).dropna()
def fnd(sigs):
    y = dd.reindex(sigs).dropna()
    return y.mean() if len(y) else np.nan
lo, hi = boot_ci(fnd, 916001)
pred.append({"test": f"bladder inflamed rho, response (GSE176307) minus survival ({bl})", "estimate": dd.mean(),
             "lo": lo, "hi": hi, "holds": dd.mean() > 0, "confirmed": lo > 0})
for f, lab in ((CON, "confounder"), ("Other", "other")):
    mm, fn2 = fam_mean(Rr, "GSE176307", f); l2, h2 = boot_ci(fn2, 916001)
    pred.append({"test": f"GSE176307 {lab} rho (response; descriptive)", "estimate": mm, "lo": l2, "hi": h2,
                 "holds": np.nan, "confirmed": np.nan})
PR = pd.DataFrame(pred)
PR.to_csv(f"{OUT}/predictions.csv", index=False)

# ------------------------------------------------ 4. library level over the external survival cohorts
def sig_level(c):
    def one(x):
        p = wilcoxon(x.G).pvalue if (x.G != 0).any() else 1.0
        return pd.Series({"ncoh": len(x), "rho_bar": x.rho.mean(), "med_G": x.G.median(), "p_content": p,
                          "content_call": bool((p < 0.05) and (x.G.median() > 0))})
    s = c.groupby("target")[["G", "rho"]].apply(one)
    s["family"] = s.index.map(FAM.family)
    return s
SL = sig_level(S)
calls = pd.read_csv(f"{V16}/external_surv_calls.csv.gz")
ct = calls[(calls.kind == "target") & (calls.stat == "fold")].set_index("target").call.astype(bool)
SL["called"] = ct.reindex(SL.index)
SL.to_csv(f"{OUT}/signature_external_survival.csv")
lib = []
for f in FAMS + ["all"]:
    x = SL if f == "all" else SL[SL.family == f]
    lib.append({"family": f, "n": len(x), "rho_bar": x.rho_bar.mean(), "content_calls": int(x.content_call.sum()),
                "calls": int(x.called.sum())})
LB = pd.DataFrame(lib)
LB.to_csv(f"{OUT}/library_external_survival.csv", index=False)
FD = pd.read_csv(f"{V16}/external_surv_fold_fdr.csv")
pd.DataFrame([prim]).to_csv(f"{OUT}/primary.csv", index=False)
pd.set_option("display.width", 200)
print(pd.DataFrame([prim])); print(PR); print(LB); print(FD[["called", "fdr", "fdr_lo", "fdr_hi", "pi0", "fdr_pi0"]])
print(PC[["cohort", "n", "calls", "fdr", "rho_inflamed", "rho_inflamed_lo", "rho_inflamed_hi", "tcga_rho_inflamed",
          "rho_confounder"]].round(3))
