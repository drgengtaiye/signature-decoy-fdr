"""Acceptance test of check_models() (validation/TOOL_check_models_acceptance.md): the outputs of
applications/v24_check_models_repro.R against round 15 (icb/v20_round15).

    python3 -I v24_check_models_repro.py [TYPE ...]     -> icb/v24_tool/repro_report.csv, printed summary
"""
import sys
import numpy as np, pandas as pd

G = "/groups/adv2105_gp/yichen/Yi/com/gastric"
A, T24 = f"{G}/icb/v20_round15", f"{G}/icb/v24_tool"
TYPES = ["STAD", "SKCM", "LUAD", "LUSC", "COAD", "CESC", "SARC", "HNSC", "BLCA", "OV", "BRCA"]
EXTID = dict(STAD="STAD_GSE84437", SKCM="SKCM_GSE65904", LUAD="LUAD_GSE72094", LUSC="LUSC_GSE157011",
             COAD="COAD_GSE39582", CESC="CESC_GSE44001", SARC="SARC_GSE21050", HNSC="HNSC_GSE65858",
             BLCA="BLCA_GSE13507", OV="OV_GSE32062", BRCA="BRCA_GSE96058")
KNOWN_B5 = {("BRCA_GSE96058", "BRCA_PMC11332736")}      # stated in the acceptance note before the run
rep = []


def check(tp, crit, what, ok, detail=""):
    rep.append(dict(type=tp, criterion=crit, check=what, passed=bool(ok), detail=detail))


def same_num(a, b, atol, rtol=0.0):
    a, b = np.asarray(a, float), np.asarray(b, float)
    both_na = np.isnan(a) & np.isnan(b)
    close = np.isclose(a, b, atol=atol, rtol=rtol)
    bad = ~(both_na | close)
    return int(bad.sum()), float(np.nanmax(np.abs(a - b))) if np.isfinite(a - b).any() else 0.0


def one(tp):
    R = pd.read_csv(f"{A}/models_{tp}.csv.gz")
    DI = pd.read_csv(f"{A}/models_decoy_info_{tp}.csv.gz")
    S = pd.read_csv(f"{T24}/repro_scores_{tp}.csv.gz")
    D = pd.read_csv(f"{T24}/repro_decoys_{tp}.csv.gz")
    Un = pd.read_csv(f"{T24}/repro_units_{tp}.csv.gz")
    ext = EXTID[tp]
    # 1. decoys
    check(tp, 1, "decoy names and order", list(D.decoy) == list(DI.decoy), f"{len(D)} vs {len(DI)}")
    check(tp, 1, "decoy genes", (D.genes.values == DI.genes.values).all() if len(D) == len(DI) else False)
    Rd = R[R.kind == "decoy"].set_index("decoy")
    check(tp, 1, "decoy coefficients (8 significant digits)", (D.set_index("decoy").coefs == Rd.weights.loc[D.decoy]).all())
    check(tp, 1, "decoy directions", (D.set_index("decoy").dirn == Rd.dirn.loc[D.decoy]).all())
    # 2. part A in TCGA and the external cohort
    for co, pre in ((tp, "tcga"), (ext, "ext")):
        x = S[S.cohort == co].copy()
        x["key"] = np.where(x.kind == "model", x.id, x.id)
        r = R.assign(key=np.where(R.kind == "target", R.model, R.decoy)).set_index("key")
        x = x.set_index("key")
        check(tp, 2, f"{pre}: same models and decoys", set(x.index) == set(r.index), f"{len(x)} vs {len(r)}")
        r = r.loc[x.index]
        nb, mx = same_num(x.cdir, r[f"{pre}_cdir"], 1e-10)
        check(tp, 2, f"{pre}: directed C-index within 1e-10", nb == 0, f"{nb} differ, max |diff| {mx:.2e}")
        nb, mx = same_num(x.p, r[f"{pre}_p"], 0, 1e-8)
        check(tp, 2, f"{pre}: log-rank P within 1e-8 relative", nb == 0, f"{nb} differ, max |diff| {mx:.2e}")
        nb, _ = same_num(x.worse, r[f"{pre}_worse"], 0)
        check(tp, 2, f"{pre}: high group worse", nb == 0, f"{nb} differ")
        if pre == "ext":
            nb, _ = same_num(x.validated, r.validated, 0)
            check(tp, 2, "ext: validated identical", nb == 0, f"{nb} differ")
    # 3. per model (audit_models_units)
    AU = pd.read_csv(f"{A}/audit_models_units.csv")
    AU = AU[AU.type == tp].set_index("model")
    u = Un[Un.cohort == ext].set_index("model")
    check(tp, 3, "same models", set(u.index) == set(AU.index), f"{len(u)} vs {len(AU)}")
    AU = AU.loc[u.index]
    has = AU.n_coex.notna()
    check(tp, 3, "models with decoys", (has == (u.n_coex > 0)).all(), f"{int(has.sum())} with decoys")
    for c in ("d_coex", "d_rand", "n_coex", "rho"):
        nb, mx = same_num(u.loc[has, c], AU.loc[has, c], 0)
        check(tp, 3, f"{c} identical", nb == 0, f"{nb} differ, max |diff| {mx:.2e}")
    tv = ((u.train_p < 0.05) & (u.train_worse == 1)).astype(int)
    check(tp, 3, "TCGA validated identical", (tv == AU.tcga_validated).all() and (u.train_validated == AU.tcga_validated).all())
    # 5. part B5 (published models and decoys)
    B5 = pd.read_csv(f"{A}/cohort_scores_{tp}.csv.gz")
    B5["kind"] = B5["kind"].fillna("null")
    B5 = B5[B5.kind.isin(["target", "decoy"])].copy()
    B5["kind"] = B5.kind.replace({"target": "model"})
    m = B5.merge(S[S.cohort != tp], on=["cohort", "kind", "id"], how="outer", suffixes=("_b", "_f"), indicator=True)   # part B did not score TCGA
    check(tp, 5, "same rows", (m._merge == "both").all(), m._merge.value_counts().to_dict().__repr__())
    m = m[m._merge == "both"]
    check(tp, 5, "scorable identical", (m.scorable_b.astype(bool) == m.scorable_f.astype(bool)).all())
    nb, _ = same_num(m.validated_b, m.validated_f, 0)
    check(tp, 5, "validated identical", nb == 0, f"{nb} differ")
    known = np.array([(c, i) in KNOWN_B5 for c, i in zip(m.cohort, m.id)])
    nb, mx = same_num(m.cdir_b[~known], m.cdir_f[~known], 1e-6)
    check(tp, 5, "directed C-index within 1e-6 (known exception aside)", nb == 0, f"{nb} differ, max |diff| {mx:.2e}")
    if known.any():
        k = m[known]
        check(tp, 5, "known exception: as stated", True,
              f"{k.id.iloc[0]} in {k.cohort.iloc[0]}: part B {k.measured_b.iloc[0]} genes, C {k.cdir_b.iloc[0]:.6f}; "
              f"function {k.measured_f.iloc[0]} genes, C {k.cdir_f.iloc[0]:.6f}; validated {k.validated_b.iloc[0]:.0f} / {k.validated_f.iloc[0]:.0f}")
    return u.assign(type=tp)


types = sys.argv[1:] or TYPES
U = pd.concat([one(tp) for tp in types])
if set(types) == set(TYPES):
    # 4. estimates
    PM = pd.read_csv(f"{G}/signature-decoy-fdr/validation/published_models_r15.csv").set_index("model_id")
    U = U[U.n_coex > 0].join(PM.uses_external)
    AM = pd.read_csv(f"{A}/audit_models.csv")
    for pop, X in (("primary: independent of the external cohort", U[~U.uses_external.astype(bool)]),
                   ("all models with decoys", U)):
        for dc in ("coex", "rand"):
            est = dict(M=len(X), validated=X.validated.sum(), validated_share=X.validated.mean(),
                       decoy_share=X[f"d_{dc}"].sum() / len(X), fdr=X[f"d_{dc}"].sum() / X.validated.sum())
            if dc == "coex":
                est["rho"] = X.rho.mean()
            ref = AM[(AM.population == pop) & (AM.decoys == dc)].set_index("quantity").estimate
            for q, v in est.items():
                d = abs(v - ref[q])
                check("all", 4, f"{pop}, {dc}: {q}", d < 1e-12, f"function {v:.10g}, round 15 {ref[q]:.10g}, |diff| {d:.1e}")
out = pd.DataFrame(rep)
out.to_csv(f"{T24}/repro_report{'' if set(types) == set(TYPES) else '_' + '_'.join(types)}.csv", index=False)
print(out.groupby("criterion").passed.agg(["sum", "size"]))
fail = out[~out.passed]
print("FAILED:" if len(fail) else "all checks passed")
print(fail.to_string() if len(fail) else "")
print(out[out.check.str.startswith("known")].detail.to_string())
