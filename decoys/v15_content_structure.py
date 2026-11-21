"""Analyses 1-4 of validation/PLAN_round10_content_and_structure.md: what a signature achieves over its null, split into
what modules of non-signature genes with its structure achieve (structure) and what its own genes add (content), by
endpoint, signature family and cancer type.

    python v15_content_structure.py   -> icb/v15_content/*.csv

Inputs are the stored per-cohort outputs of the round-7 package runs; nothing is recomputed from expression.
"""
import os, numpy as np, pandas as pd
from scipy.stats import wilcoxon

ICB = "/path/to/analysis_root/icb"
PKG, OUT = f"{ICB}/v12_pkg", f"{ICB}/v15_content"
os.makedirs(OUT, exist_ok=True)

FAM = pd.read_csv(f"{ICB}/v5_signature_clusters.csv").set_index("signature")
ICB_TUM = {"ICB_Gide": "melanoma", "ICB_Hugo": "melanoma", "ICB_Liu": "melanoma", "ICB_Riaz": "melanoma",
           "ICB_Van_Allen": "melanoma", "ICB_Puch": "melanoma", "ICB_Braun": "kidney", "ICB_Miao1": "kidney",
           "ICB_Jung": "lung", "ICB_Limagne1": "lung", "ICB_Mariathasan": "bladder", "ICB_Kim": "gastric"}
TCGA_TUM = {"SKCM": "melanoma", "KIRC": "kidney", "LUAD": "lung", "LUSC": "lung", "BLCA": "bladder", "STAD": "gastric"}
FAMS = ["Inflamed/IFN", "Other", "Confounder"]


def cells(name, stat="fold", drop=None):
    """Per target x cohort: the target's Delta, structure S, content G, place among decoys rho, and clearance."""
    d = pd.read_csv(f"{PKG}/{name}_units.csv.gz")
    d["delta"] = d[f"{stat}_obs"] - d[f"{stat}_med"]
    d["clear"] = d[f"{stat}_obs"] > d[f"{stat}_p95"]
    if drop is not None:  # origin tests: removed for the target and its decoys
        d = d.merge(drop.assign(_x=1), on=["target", "cohort"], how="left")
        d = d[d._x.isna()].drop(columns="_x")
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
    return out.reset_index()


def sig_level(c):
    """Per target: mean rho over cohorts, median S and G, and the content call (signed-rank of G, P < 0.05, median > 0)."""
    def one(x):
        p = wilcoxon(x.G).pvalue if (x.G != 0).any() else 1.0
        return pd.Series({"ncoh": len(x), "rho_bar": x.rho.mean(), "rho_bar_rand": x.rho_rand.mean(),
                          "med_delta": x.delta.median(), "med_S": x.S.median(), "med_G": x.G.median(), "p_content": p,
                          "content_call": (p < 0.05) and (x.G.median() > 0)})
    s = c.groupby("target").apply(one)
    s["content_call"] = s.content_call.astype(bool)
    return s


def calls(name, stat="fold", origin=False):
    f = f"{PKG}/{name}_origin_calls.csv.gz" if origin else f"{PKG}/{name}_calls.csv.gz"
    c = pd.read_csv(f)
    c = c[(c.kind == "target") & (c.stat == stat)]
    return c.set_index("target").call.astype(bool)


SETTINGS = {  # setting: (units name, statistic, cohorts kept or None)
    "response": ("response", "fold", None),
    "tcga_OS": ("tcga_OS", "fold", None),
    "tcga_DSS": ("tcga_DSS", "fold", None),
    "tcga_PFI": ("tcga_PFI", "fold", None),
    "tcga_OS_matched": ("tcga_OS", "fold", list(TCGA_TUM)),
    "response_dir": ("response", "dir", None),
}
C, S = {}, {}
for k, (nm, st, keep) in SETTINGS.items():
    c = cells(nm, st)
    if keep is not None:
        c = c[c.cohort.isin(keep)]
    C[k] = c
    S[k] = sig_level(c)
    S[k]["family"] = S[k].index.map(FAM.family)
    S[k]["called"] = calls(nm, st).reindex(S[k].index)
    c.to_csv(f"{OUT}/cells_{k}.csv.gz", index=False)
    S[k].to_csv(f"{OUT}/signature_{k}.csv")

# ---------------------------------------------------------------- 1. content by endpoint and family (primary)
rng = np.random.default_rng(915000)
CL = FAM.cluster_fixed


def contrast(a, b):
    """Mean over inflamed minus mean over confounder of C(s) = rho_bar(s, a) - rho_bar(s, b)."""
    Cs = S[a].rho_bar - S[b].rho_bar
    fam = FAM.family.reindex(Cs.index)
    obs = Cs[fam == "Inflamed/IFN"].mean() - Cs[fam == "Confounder"].mean()
    # P: random reassignment of the two labels among the 28 signatures
    v = Cs[fam.isin(["Inflamed/IFN", "Confounder"])].values
    n_inf = int((fam == "Inflamed/IFN").sum())
    B = 100_000
    perm = np.empty(B)
    for i in range(B):
        o = rng.permutation(len(v))
        perm[i] = v[o[:n_inf]].mean() - v[o[n_inf:]].mean()
    p = (1 + np.sum(np.abs(perm) >= abs(obs) - 1e-12)) / (1 + B)
    return Cs, obs, p


def cluster_boot(stat_fn, seed, B=10_000):
    """Percentile interval of stat_fn(list of signatures) over resamples of the 13 correlation clusters."""
    r = np.random.default_rng(seed)
    cl = CL.unique()
    members = {k: list(CL.index[CL == k]) for k in cl}
    vals = []
    for _ in range(B):
        sigs = [s for k in r.choice(cl, len(cl), replace=True) for s in members[k]]
        v = stat_fn(sigs)
        if np.isfinite(v):
            vals.append(v)
    return np.percentile(vals, [2.5, 97.5]), len(vals)


rows = []
for a, b in [("response", "tcga_OS"), ("response", "tcga_OS_matched"), ("response", "tcga_DSS"),
             ("response", "tcga_PFI")]:
    Cs, obs, p = contrast(a, b)

    def fn(sigs, Cs=Cs):
        x = Cs.reindex(sigs); f = FAM.family.reindex(sigs).values
        i, c = x[f == "Inflamed/IFN"], x[f == "Confounder"]
        return i.mean() - c.mean() if len(i) and len(c) else np.nan
    (lo, hi), nb = cluster_boot(fn, 915001)
    fam = FAM.family.reindex(Cs.index)
    rows.append({"a": a, "b": b, "contrast_inflamed_minus_confounder": obs, "p_perm": p, "lo": lo, "hi": hi,
                 "boot_ok": nb, "C_inflamed": Cs[fam == "Inflamed/IFN"].mean(), "C_other": Cs[fam == "Other"].mean(),
                 "C_confounder": Cs[fam == "Confounder"].mean(), "C_all": Cs.mean(),
                 "paired_wilcoxon_p_all52": wilcoxon(S[a].rho_bar, S[b].rho_bar.reindex(S[a].index)).pvalue})
pd.DataFrame(rows).to_csv(f"{OUT}/primary_contrast.csv", index=False)

fr = []
for k in SETTINGS:
    s = S[k]
    for f in FAMS + ["all"]:
        sel = s.index if f == "all" else s.index[s.family == f]

        def fn(sigs, k=k, f=f):
            x = S[k].reindex(sigs)
            x = x if f == "all" else x[x.family == f]
            return x.rho_bar.mean() if len(x) else np.nan
        (lo, hi), _ = cluster_boot(fn, 915002)
        cc = C[k][C[k].target.isin(sel)]
        fr.append({"setting": k, "family": f, "n": len(sel), "rho_bar": s.loc[sel, "rho_bar"].mean(), "lo": lo,
                   "hi": hi, "rho_bar_rand": s.loc[sel, "rho_bar_rand"].mean(),
                   "content_calls": int(s.loc[sel, "content_call"].sum()), "calls": int(s.loc[sel, "called"].sum()),
                   "med_S": cc.S.median(), "med_G": cc.G.median(), "med_delta": cc.delta.median(),
                   "decoy_clear": cc.dclear.mean(), "rand_clear": cc.dclear_rand.mean(),
                   "target_clear": cc.clear.mean()})
pd.DataFrame(fr).to_csv(f"{OUT}/family_by_setting.csv", index=False)

# ---------------------------------------------------------------- 2. the same cancer types
tr = []
for tum in ["melanoma", "kidney", "lung", "bladder", "gastric"]:
    for k, cmap in [("response", ICB_TUM), ("tcga_OS", TCGA_TUM)]:
        cc = C[k][C[k].cohort.map(cmap) == tum].copy()
        cc["family"] = cc.target.map(FAM.family)
        for f in FAMS + ["all"]:
            x = cc if f == "all" else cc[cc.family == f]
            tr.append({"tumor": tum, "setting": k, "family": f, "cohorts": cc.cohort.nunique(),
                       "rho": x.rho.mean(), "med_G": x.G.median(), "decoy_clear": x.dclear.mean()})
pd.DataFrame(tr).to_csv(f"{OUT}/matched_tumor_types.csv", index=False)

# ---------------------------------------------------------------- 3. pan-cancer atlas
at = []
for k in ["tcga_OS", "response"]:
    cc = C[k].copy()
    cc["family"] = cc.target.map(FAM.family)
    cc.pivot_table(index="target", columns="cohort", values="rho").to_csv(f"{OUT}/atlas_rho_{k}.csv")
    for co, x in cc.groupby("cohort"):
        R = int(x.clear.sum())
        ef = x.dclear.sum()
        at.append({"setting": k, "cohort": co, "n": int(x.n.iloc[0]), "rho_all": x.rho.mean(),
                   "rho_inflamed": x[x.family == "Inflamed/IFN"].rho.mean(),
                   "rho_confounder": x[x.family == "Confounder"].rho.mean(),
                   "rho_other": x[x.family == "Other"].rho.mean(), "calls": R, "exp_false": ef,
                   "fdr": 100 * ef / R if R else np.nan, "decoy_clear": x.dclear.mean(),
                   "calls_inflamed": int(x[x.family == "Inflamed/IFN"].clear.sum()),
                   "calls_confounder": int(x[x.family == "Confounder"].clear.sum())})
pd.DataFrame(at).to_csv(f"{OUT}/atlas_by_type.csv", index=False)

# ---------------------------------------------------------------- 4. libraries
origin = pd.read_csv(f"{PKG}/breast_origin_excluded_tests.csv")
LIB = [("response", "response", None, False, "response_fdr.csv", None),
       ("tcga_OS", "tcga_OS", None, False, "tcga_OS_fdr.csv", None),
       ("tcga_DSS", "tcga_DSS", None, False, "tcga_DSS_fdr.csv", None),
       ("tcga_PFI", "tcga_PFI", None, False, "tcga_PFI_fdr.csv", None),
       ("breast_pc1", "breast_pc1", None, False, "breast_pc1_fdr.csv", None),
       ("breast_pc1_origin_out", "breast_pc1", origin, True, "breast_origin_fdr.csv", ("pc1", "D and C")),
       ("breast_mean", "breast_mean", None, False, "breast_mean_fdr.csv", None),
       ("hallmark_OS", "hallmark_OS", None, False, "hallmark_OS_fdr.csv", None)]
lr = []
for lab, nm, drop, org, ff, key in LIB:
    cc = cells(nm, "fold", drop)
    s = sig_level(cc)
    s["called"] = calls(nm, "fold", org).reindex(s.index)
    fd = pd.read_csv(f"{PKG}/{ff}")
    fd = fd[(fd.stat == "fold") & (fd["class"] == "coex") & (fd.boot == "targets")]
    if key is not None:
        fd = fd[(fd.score == key[0]) & (fd.rule == key[1])]
    fd = fd.iloc[0]
    lr.append({"library": lab, "m": len(s), "called": int(s.called.sum()), "rho_bar_all": s.rho_bar.mean(),
               "rho_bar_called": s[s.called].rho_bar.mean(), "rho_bar_rand_all": s.rho_bar_rand.mean(),
               "content_calls": int(s.content_call.sum()), "med_S": cc.S.median(), "med_G": cc.G.median(),
               "decoy_clear": cc.dclear.mean(), "rand_clear": cc.dclear_rand.mean(), "target_clear": cc.clear.mean(),
               "fdr_reported": fd.fdr, "pi0_reported": fd.pi0, "called_reported": int(fd.called)})
    s.to_csv(f"{OUT}/library_signature_{lab}.csv")
pd.DataFrame(lr).to_csv(f"{OUT}/libraries.csv", index=False)
print("done")
