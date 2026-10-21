"""Within each arm of each randomized trial, against the same matched random-gene-set null.

The benchmark compares response under treatment with survival without it, in different
patients on different platforms. This script makes the comparison inside the trials, where
the patients, the disease, the platform and the endpoint are the same and only the treatment
differs: for every signature it fits progression-free survival within the checkpoint arm and
within the comparator arm, and refers each fit to random gene sets of the signature's own size
drawn from the same genes and scored the same way in the same patients.

    python v8_within_arm_null.py    -> pancan/v8_within_arm_null.csv (+ _rungs.csv)

Pre-specified in validation/PLAN_icb_survival.md's companion paragraph: the statistic is the
share of immune signature tests clearing the null's one-sided 5th percentile in the benefit
direction, in each arm, with the confounder controls as the internal comparison.
"""
import warnings, functools, sys
warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
import v7_lib as L

P = L.P
LADDER = [4, 7, 12, 20, 32, 50, 80, 130, 260, 440]
NRAND = 200
rng = np.random.default_rng(8)


def matrices(trial):
    if trial == "IMmotion150":
        C = pd.read_csv(f"{P}/v4_immotion150_clinical_api.csv", index_col=0)
        M = pd.read_csv(f"{P}/ici/rcc_iatlas_immotion150_2018.csv", index_col=0)
        arm = (C["ICI_RX"].astype(str) == "Atezolizumab").astype(int)
    else:
        C = pd.read_csv(f"{P}/v2_prince_clinical_api.csv", index_col=0)
        M = pd.read_csv(f"{P}/ici/paad_iatlas_prince_2022.csv", index_col=0)
        arm = C["ICI_RX"].astype(str).str.contains("Nivolumab", case=False, na=False).astype(int)
    keep = [s for s in M.columns if s in C.index]
    C, M, arm = C.loc[keep], M[keep], arm.loc[keep]
    t, e = L._pfs(C)
    ok = t.notna() & (t > 0) & arm.notna()
    Z = ((M.T - M.mean(1)) / M.std(1).replace(0, 1)).T
    return Z.loc[:, ok[ok].index], arm[ok].astype(int), t[ok].astype(float), e[ok].astype(int)


def coef(score, t, e):
    """Cox coefficient per standard deviation of the score, within one arm."""
    s = np.asarray(score, float)
    sd = s.std()
    if not np.isfinite(sd) or sd == 0:
        return np.nan
    b, V, ok = L.cox_fit(((s - s.mean()) / sd)[:, None], np.asarray(t, float), np.asarray(e, int))
    return float(b[0]) if ok else np.nan


SIGS = L.signatures()
rows, rungs = [], []
for trial in ("IMmotion150", "PRINCE"):
    Z, arm, t, e = matrices(trial)
    pool = list(Z.index)
    print(f"\n=== {trial}: n={len(arm)} (checkpoint {int(arm.sum())}, comparator {int((1-arm).sum())}), "
          f"events {int(e.sum())}, pool {len(pool)} genes")
    for a_, arm_name in ((1, "checkpoint"), (0, "comparator")):
        m = (arm == a_).values
        ta, ea = t.values[m], e.values[m]
        null = {}
        for k in LADDER:
            if k > len(pool):
                continue
            v = []
            for _ in range(NRAND):
                g = rng.choice(pool, k, replace=False)
                c = coef(Z.loc[g].values[:, m].mean(0), ta, ea)
                if np.isfinite(c):
                    v.append(c)
            v = np.array(v)
            null[k] = v
            rungs.append(dict(trial=trial, arm=arm_name, size=k, n_draws=len(v),
                              med=float(np.median(v)), p05=float(np.quantile(v, .05)),
                              sd=float(v.std(ddof=1))))
        for nm, gl in SIGS.items():
            g = [x for x in gl if x in Z.index]
            if len(g) < 3:
                continue
            fam = "Immune" if nm in L.IMM else "Control" if nm.startswith("CTRL_") else "Other"
            if fam == "Other":
                continue
            k = min(LADDER, key=lambda x: abs(x - len(g)))
            c = coef(Z.loc[g].values[:, m].mean(0), ta, ea)
            v = null[k]
            rows.append(dict(trial=trial, arm=arm_name, signature=nm, family=fam, n_genes=len(g),
                             null_size=k, coef=c, hr=float(np.exp(c)),
                             null_med=float(np.median(v)), null_p05=float(np.quantile(v, .05)),
                             null_sd=float(v.std(ddof=1)),
                             z=(c - float(np.median(v))) / float(v.std(ddof=1)),
                             clears=bool(c < np.quantile(v, .05))))
R = pd.DataFrame(rows); R.to_csv(f"{P}/v8_within_arm_null.csv", index=False)
pd.DataFrame(rungs).to_csv(f"{P}/v8_within_arm_null_rungs.csv", index=False)

print("\n=== share of tests clearing the null's 5th percentile in the benefit direction")
for fam in ("Immune", "Control"):
    g = R[R.family == fam]
    for arm_name in ("checkpoint", "comparator"):
        h = g[g.arm == arm_name]
        print(f"  {fam:8s} {arm_name:11s} {int(h.clears.sum()):2d} of {len(h):2d} "
              f"({100*h.clears.mean():4.1f}%)  median HR {np.exp(h.coef.median()):.2f}  median z {h.z.median():+.2f}")
print("\n=== per trial, immune signatures")
for trial in ("IMmotion150", "PRINCE"):
    for arm_name in ("checkpoint", "comparator"):
        h = R[(R.trial == trial) & (R.arm == arm_name) & (R.family == "Immune")]
        print(f"  {trial:12s} {arm_name:11s} {int(h.clears.sum()):2d} of {len(h):2d}  "
              f"median HR {np.exp(h.coef.median()):.2f}  median z {h.z.median():+.2f}")
