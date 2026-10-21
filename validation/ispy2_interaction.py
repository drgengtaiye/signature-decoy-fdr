"""Registration, secondary 5: randomized score x arm interaction on pCR in I-SPY2.

For each signature, pCR ~ score + arm + score:arm (logistic; score standardized within the
contrast). A hit is a checkpoint-favoring interaction (coefficient > 0) at Wald P < 0.05.
The count of hits is referred to the no-interaction simulation: pCR drawn from the fitted
pCR ~ PC1..PCK + arm model (K = 9, the registered value; the ladder 3, 4, 15, 20 is reported
as sensitivity), every signature refitted, 2,000 replicates.
Two libraries: the trial library of the survival trials (17 immune + 3 control gene sets, mean
gene z-score, pancan/v7_lib.py) and the 50 locked signatures (native scores)."""
import sys, numpy as np, pandas as pd
from scipy.stats import norm
from joblib import Parallel, delayed
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
from v7_lib import signatures, IMM, pc_design
V = "/groups/adv2105_gp/yichen/Yi/com/gastric/validation"; G = "/groups/adv2105_gp/yichen/Yi/com/gastric"
NREP = 2000; LADDER = [3, 4, 9, 15, 20]; KREG = 9
LOCK = pd.read_csv(f"{G}/icb/v3_LOCKED_predictions.csv").set_index("signature")


def logit_fit(X, y, iters=30):
    """Newton-Raphson logistic regression for many designs at once is unnecessary here; one
    design at a time, returning coefficients and Wald standard errors."""
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        eta = np.clip(X @ b, -30, 30); p = 1 / (1 + np.exp(-eta)); W = p * (1 - p)
        H = X.T @ (X * W[:, None]) + 1e-8 * np.eye(X.shape[1]); g = X.T @ (y - p)
        step = np.linalg.solve(H, g); b += step
        if np.max(np.abs(step)) < 1e-9: break
    eta = np.clip(X @ b, -30, 30); p = 1 / (1 + np.exp(-eta)); W = p * (1 - p)
    V_ = np.linalg.inv(X.T @ (X * W[:, None]) + 1e-8 * np.eye(X.shape[1]))
    return b, np.sqrt(np.diag(V_)), p


def interaction(s, arm, y):
    X = np.column_stack([np.ones_like(s), s, arm, s * arm])
    b, se, _ = logit_fit(X, y)
    return b[3], 2 * norm.sf(abs(b[3] / se[3]))


def hits(S, arm, y, mask):
    n = 0; out = []
    for i in range(S.shape[0]):
        b, p = interaction(S[i], arm, y); out.append((b, p))
        if mask[i] and np.isfinite(p) and p < .05 and b > 0: n += 1
    return n, out


def run(contrast, lib):
    X = pd.read_csv(f"{V}/cohorts/{contrast}.expr.tsv.gz", sep="\t", index_col=0)
    C = pd.read_csv(f"{V}/cohorts/{contrast}.clin.tsv", sep="\t"); C["sample"] = C["sample"].astype(str)
    X.columns = [str(c) for c in X.columns]; X = X[C["sample"]]
    arm, y = C.arm.to_numpy(float), C.y.to_numpy(float)
    if lib == "trial library":
        Z = ((X.T - X.mean(1)) / X.std(1).replace(0, 1)).T
        names, fams, rows = [], [], []
        for nm, gl in signatures().items():
            g = [x for x in gl if x in Z.index]
            if len(g) < 3: continue
            fam = "Immune" if nm in IMM else "Control" if nm.startswith("CTRL_") else "Other"
            if fam == "Other": continue
            names.append(nm); fams.append(fam); rows.append(Z.loc[g].mean(0).values)
        S = np.vstack(rows); mask = np.array([f == "Immune" for f in fams])
    else:
        N = pd.read_csv(f"{V}/cohorts/{contrast}.native_scores.csv", index_col=0); N.index = N.index.astype(str)
        N = N.loc[C["sample"], [s for s in LOCK.index if s in N.columns]]
        names = list(N.columns); fams = ["immune-attributed" if LOCK.immune_attributed[s] == 1 else "other" for s in names]
        S = N.to_numpy(float).T; mask = np.ones(len(names), bool)
    S = (S - S.mean(1, keepdims=True)) / (S.std(1, keepdims=True) + 1e-9)
    obs, per = hits(S, arm, y, mask)
    res = dict(contrast=contrast, library=lib, n=len(y), n_arm=int(arm.sum()), pcr_arm=int(y[arm == 1].sum()),
               pcr_ctrl=int(y[arm == 0].sum()), n_ctrl=int((arm == 0).sum()), tested=int(mask.sum()), hits=obs)
    for K in LADDER:
        pcs, cumvar = pc_design(S, K)
        Xr = np.column_stack([np.ones(len(y)), pcs, arm]); _, _, p0 = logit_fit(Xr, y)

        def one(r):
            rng = np.random.default_rng(10_000 * K + r)
            ys = (rng.random(len(y)) < p0).astype(float)
            return hits(S, arm, ys, mask)[0]
        null = np.array(Parallel(n_jobs=48)(delayed(one)(r) for r in range(NREP)))
        res[f"P_K{K}"] = float((np.sum(null >= obs) + 1) / (NREP + 1))
        res[f"null_med_K{K}"] = float(np.median(null)); res[f"null_p95_K{K}"] = float(np.quantile(null, .95))
        res[f"var_K{K}"] = float(cumvar[min(K, len(cumvar)) - 1])
    per_df = pd.DataFrame([dict(contrast=contrast, library=lib, signature=n_, family=f_, beta=b, p=p) for n_, f_, (b, p) in zip(names, fams, per)])
    return res, per_df


if __name__ == "__main__":
    R, P = [], []
    for contrast in ["ISPY2_pembro", "ISPY2_durva"]:
        for lib in ["trial library", "locked signatures"]:
            r, p = run(contrast, lib); R.append(r); P.append(p)
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
    pd.DataFrame(R).to_csv(f"{V}/ispy2_interaction_summary.csv", index=False)
    pd.concat(P).to_csv(f"{V}/ispy2_interaction_by_signature.csv", index=False)
