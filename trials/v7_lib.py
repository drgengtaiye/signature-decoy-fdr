"""Machinery shared by the v7 interaction-null scripts.

v6 ran four nulls for the score x treatment interaction and described three of them as
retaining the main effects. Two of those descriptions are wrong:

  * permuting arm labels within strata of the score still breaks the treatment main
    effect whenever beta_arm != 0, because the outcomes stay attached to the ORIGINAL
    assignment;
  * permuting the score within arm breaks the score's own main effect whenever
    beta_score != 0.

The hypothesis the paper needs is

    H0: beta(score x arm) = 0,  with beta_arm free and beta_score free.

This module provides the pieces for nulls that actually target it: a fast Cox fitter, a
Breslow baseline for model-based simulation, per-subject score residuals, and the
efficient (nuisance-projected) score for the interaction term, which is what lets a
multiplier bootstrap hold both main effects fixed while resampling only the interaction.

Everything here is self-consistent: the baseline hazard is computed from the same
uncentred linear predictor that is used to simulate from it, so no centring convention
has to be matched against an external library.
"""
import warnings, subprocess, functools
warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd

P = "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan"

IMM = {'T_cells', 'CD8_Tcells', 'Cytotoxic', 'NK_cells', 'B_cells', 'Plasma_cells',
       'Exhaustion', 'Checkpoints', 'GEP18_TIS (Ayers 2017)', 'IFNG_6gene (Ayers 2017)',
       'TLS (Cabrita 2020)', 'ESTIMATE_Immune', 'barriernet_Immune', 'APM', 'Treg', 'DC',
       'Chemokine12 (Messina 2012)', 'TLS_9gene (Prabhakaran)', 'CYT (Rooney 2015)'}


# --------------------------------------------------------------------------- data
def signatures():
    raw = subprocess.run(
        ["Rscript", "-e", f'source("{P}/signatures.R");SIGS<-SIGS[!sapply(SIGS,is.null)];'
         'for(n in names(SIGS)) cat(n,"\\t",paste(SIGS[[n]],collapse=","),"\\n",sep="")'],
        capture_output=True, text=True).stdout
    return {k.strip(): [x for x in v.split(",") if x]
            for k, v in (l.split("\t", 1) for l in raw.strip().split("\n") if "\t" in l)}


def _pfs(C):
    t = pd.to_numeric(C["PFS_MONTHS"], errors="coerce")
    s = C["PFS_STATUS"].astype(str)
    e = pd.to_numeric(s.str.extract(r"^\s*([01])\s*:")[0], errors="coerce")
    if e.isna().all():
        e = s.str.contains("1|PROGRESS", case=False, regex=True).astype(float)
    return t, e.fillna(0).astype(int)


def build(trial, SIGS=None, family_filter=True):
    """Scores, arm, time and event for one trial, exactly as v4/v5/v6 assembled them.

    family_filter=True reproduces the 20 signatures (17 immune + 3 control) that enter
    the interaction analysis; False keeps all 31 that clear the gene-presence filter.
    """
    SIGS = SIGS or signatures()
    if trial == "IMmotion150":
        C = pd.read_csv(f"{P}/v4_immotion150_clinical_api.csv", index_col=0)
        M = pd.read_csv(f"{P}/ici/rcc_iatlas_immotion150_2018.csv", index_col=0)
        arm = (C["ICI_RX"].astype(str) == "Atezolizumab").astype(int)
    else:
        C = pd.read_csv(f"{P}/v2_prince_clinical_api.csv", index_col=0)
        M = pd.read_csv(f"{P}/ici/paad_iatlas_prince_2022.csv", index_col=0)
        arm = C["ICI_RX"].astype(str).str.contains("Nivolumab", case=False, na=False).astype(int)
    keep = [s for s in M.columns if s in C.index]
    C, M = C.loc[keep], M[keep]
    arm = arm.loc[keep]
    t, e = _pfs(C)
    ok = t.notna() & (t > 0) & arm.notna()
    Z = ((M.T - M.mean(1)) / M.std(1).replace(0, 1)).T
    arm, t, e = arm[ok], t[ok], e[ok]
    idx = arm.index

    names, fams, scores = [], [], []
    for nm, gl in SIGS.items():
        g = [x for x in gl if x in Z.index]
        if len(g) < 3:
            continue
        fam = "Immune" if nm in IMM else "Control" if nm.startswith("CTRL_") else "Other"
        if family_filter and fam == "Other":
            continue
        names.append(nm); fams.append(fam)
        scores.append(Z.loc[g, idx].mean(0).values)
    S = np.vstack(scores)
    S = (S - S.mean(1, keepdims=True)) / (S.std(1, keepdims=True) + 1e-9)   # z per signature
    return dict(trial=trial, names=names, fams=np.array(fams), S=S,
                arm=arm.values.astype(int), t=t.values.astype(float),
                e=e.values.astype(int), n=len(idx))


# --------------------------------------------------------------------------- Cox
def _risk_sums(Xs, w, first_idx, need_S2=False):
    """Suffix sums over the sorted-by-time design: S0, S1 (and S2) on the risk sets."""
    S0 = np.cumsum(w[::-1])[::-1]
    S1 = np.cumsum((w[:, None] * Xs)[::-1], axis=0)[::-1]
    if not need_S2:
        return S0[first_idx], S1[first_idx], None
    XW = w[:, None] * Xs
    S2 = np.cumsum((XW[:, :, None] * Xs[:, None, :])[::-1], axis=0)[::-1]
    return S0[first_idx], S1[first_idx], S2[first_idx]


def _grad_hess(Xs, es, inv, first_idx, d, w, p, ties):
    """Score vector and observed information at the current beta.

    PFS is recorded in months and ties heavily (IMmotion150 has 164 events at 100
    distinct times, up to 6 at once), so the tie convention is not cosmetic: under
    Breslow one PRINCE signature crosses P = 0.05 and the published count moves from
    5/17 to 4/17. Efron is what the published observed statistics used, so it is the
    default here; Breslow stays available because the baseline-hazard simulator has to
    invert the same convention it was built from.
    """
    S0f = np.cumsum(w[::-1])[::-1]
    S1f = np.cumsum((w[:, None] * Xs)[::-1], axis=0)[::-1]
    XW = w[:, None] * Xs
    S2f = np.cumsum((XW[:, :, None] * Xs[:, None, :])[::-1], axis=0)[::-1]
    S0R, S1R, S2R = S0f[first_idx], S1f[first_idx], S2f[first_idx]
    xe = (es[:, None] * Xs).sum(0)

    if ties == "breslow":
        xbar = S1R / S0R[:, None]
        grad = xe - (d[:, None] * xbar).sum(0)
        I = (d[:, None, None] * (S2R / S0R[:, None, None]
                                 - xbar[:, :, None] * xbar[:, None, :])).sum(0)
        return grad, I

    # Efron: each of the d deaths at a tied time sees the risk set less an averaged
    # fraction of the death set, so the denominators walk down in d steps.
    nt = len(first_idx)
    ev = es > 0
    iv, wv, Xv = inv[ev], w[ev], Xs[ev]
    S0D = np.bincount(iv, weights=wv, minlength=nt)
    S1D = np.column_stack([np.bincount(iv, weights=wv * Xv[:, k], minlength=nt)
                           for k in range(p)])
    S2D = np.zeros((nt, p, p))
    for k in range(p):
        for l in range(k, p):
            v = np.bincount(iv, weights=wv * Xv[:, k] * Xv[:, l], minlength=nt)
            S2D[:, k, l] = v
            S2D[:, l, k] = v

    grad = xe.copy()
    I = np.zeros((p, p))
    for r in range(int(d.max()) if d.size else 0):
        m = d > r
        if not m.any():
            break
        frac = r / d[m]
        A0 = S0R[m] - frac * S0D[m]
        A1 = S1R[m] - frac[:, None] * S1D[m]
        A2 = S2R[m] - frac[:, None, None] * S2D[m]
        xb = A1 / A0[:, None]
        grad -= xb.sum(0)
        I += (A2 / A0[:, None, None] - xb[:, :, None] * xb[:, None, :]).sum(0)
    return grad, I


def cox_fit(X, t, e, ridge=0.0, tol=1e-9, maxit=60, ties="efron"):
    """Newton-Raphson Cox fit. Returns (beta, V, ok).

    V is the inverse observed information, so sqrt(diag(V)) are the model-based standard
    errors and the Wald test matches a standard Cox routine.
    """
    X = np.asarray(X, float)
    n, p = X.shape
    order = np.argsort(t, kind="mergesort")
    Xs, ts, es = X[order], t[order], e[order]
    uniq, first_idx, inv = np.unique(ts, return_index=True, return_inverse=True)
    d = np.bincount(inv, weights=es, minlength=len(uniq))

    beta = np.zeros(p)
    H = None
    for _ in range(maxit):
        eta = Xs @ beta
        w = np.exp(eta - eta.max())
        grad, I = _grad_hess(Xs, es, inv, first_idx, d, w, p, ties)
        H = -I
        if ridge:
            grad = grad - ridge * beta
            H = H - ridge * np.eye(p)
        try:
            step = np.linalg.solve(H, -grad)
        except np.linalg.LinAlgError:
            return None, None, False
        if not np.all(np.isfinite(step)):
            return None, None, False
        step = np.clip(step, -5.0, 5.0)          # keep the first iterations sane
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break
    try:
        V = np.linalg.inv(-H)
    except np.linalg.LinAlgError:
        return None, None, False
    if not np.all(np.isfinite(V)) or np.any(np.diag(V) <= 0):
        return None, None, False
    return beta, V, True


from scipy.stats import norm


def cox_interaction(score, arm, t, e, ties="efron"):
    """Wald coefficient and P for the interaction in  PFS ~ score + arm + score:arm."""
    X = np.column_stack([score, arm, score * arm])
    b, V, ok = cox_fit(X, t, e, ties=ties)
    if not ok:
        return np.nan, np.nan
    z = b[2] / np.sqrt(V[2, 2])
    return float(b[2]), float(2 * norm.sf(abs(z)))


def cox_arm_slope(score, arm, t, e, value, ties="efron"):
    """log HR per SD of score inside one arm."""
    m = arm == value
    sc = score[m]
    sc = (sc - sc.mean()) / (sc.std() + 1e-9)
    b, V, ok = cox_fit(sc[:, None], t[m], e[m], ties=ties)
    return float(b[0]) if ok else np.nan


# ------------------------------------------------------- residuals / efficient score
def score_residuals(X, t, e, beta):
    """Per-subject score residuals U (n x p) for a Cox model held at `beta`.

    U_ik = e_i (x_ik - xbar_k(t_i)) - exp(eta_i) * sum_{t_j <= t_i} (x_ik - xbar_k(t_j)) dL0_j

    At the MLE these sum to zero over subjects; away from it their column sums are the
    score, which is what the interaction test uses.
    """
    X = np.asarray(X, float)
    order = np.argsort(t, kind="mergesort")
    Xs, ts, es = X[order], t[order], e[order]
    uniq, first_idx, inv = np.unique(ts, return_index=True, return_inverse=True)
    d = np.bincount(inv, weights=es, minlength=len(uniq))

    eta = Xs @ beta
    w = np.exp(eta - eta.max())
    S0 = np.cumsum(w[::-1])[::-1]
    S1 = np.cumsum((w[:, None] * Xs)[::-1], axis=0)[::-1]
    S0u, S1u = S0[first_idx], S1[first_idx]
    xbar_u = S1u / S0u[:, None]
    dL0 = d / S0u
    cum_dL0 = np.cumsum(dL0)
    cum_xbar = np.cumsum(xbar_u * dL0[:, None], axis=0)

    k = inv
    U = (es[:, None] * (Xs - xbar_u[k])
         - w[:, None] * (Xs * cum_dL0[k][:, None] - cum_xbar[k]))
    out = np.empty_like(U)
    out[order] = U
    return out


def cox_information(X, t, e, beta):
    """Observed information at `beta` (Breslow)."""
    X = np.asarray(X, float)
    order = np.argsort(t, kind="mergesort")
    Xs, ts, es = X[order], t[order], e[order]
    uniq, first_idx, inv = np.unique(ts, return_index=True, return_inverse=True)
    d = np.bincount(inv, weights=es, minlength=len(uniq))
    eta = Xs @ beta
    w = np.exp(eta - eta.max())
    S0, S1, S2 = _risk_sums(Xs, w, first_idx, need_S2=True)
    xbar = S1 / S0[:, None]
    return (d[:, None, None] * (S2 / S0[:, None, None]
                                - xbar[:, :, None] * xbar[:, None, :])).sum(0)


def efficient_interaction_score(score, arm, t, e):
    """Everything a multiplier bootstrap needs for ONE signature's interaction test.

    Fits the REDUCED model  PFS ~ score + arm  (so both main effects sit at their
    observed size), then evaluates the three-column design at beta = (b_s, b_arm, 0) and
    projects the interaction's score residuals off the two nuisance ones.

    Returns dict with
      u      per-subject efficient score residual (n,)  - sums to the interaction score
      T      the interaction score statistic, sum(u)
      Ieff   efficient information (Schur complement)
      z      robust standardised statistic T / sqrt(sum u^2)
      b1     one-step coefficient estimate T / Ieff, comparable to the Cox MLE
    """
    s = np.asarray(score, float)
    Xr = np.column_stack([s, arm])
    # Breslow throughout this routine: the score residuals and the information below are
    # Breslow quantities, so the beta they are evaluated at has to be the Breslow fit.
    # The multiplier bootstrap only needs observed and null to share a convention, which
    # they do, and its reported z is compared against its own null rather than a table.
    br, Vr, ok = cox_fit(Xr, t, e, ties="breslow")
    if not ok:
        return None
    X = np.column_stack([s, arm, s * arm])
    beta0 = np.array([br[0], br[1], 0.0])
    U = score_residuals(X, t, e, beta0)
    I = cox_information(X, t, e, beta0)
    # empirical projection of the interaction residual on the nuisance residuals
    Un, Ui = U[:, :2], U[:, 2]
    try:
        B = np.linalg.lstsq(Un, Ui, rcond=None)[0]
    except np.linalg.LinAlgError:
        return None
    u = Ui - Un @ B
    T = float(u.sum())                      # == sum(Ui), since the nuisance scores vanish
    try:
        Ieff = float(I[2, 2] - I[2, :2] @ np.linalg.solve(I[:2, :2], I[:2, 2]))
    except np.linalg.LinAlgError:
        return None
    ss = float((u ** 2).sum())
    if not (np.isfinite(Ieff) and Ieff > 0 and ss > 0):
        return None
    return dict(u=u, T=T, Ieff=Ieff, z=T / np.sqrt(ss), b1=T / Ieff,
                se1=1.0 / np.sqrt(Ieff), b_red=br)


# --------------------------------------------------------------------- model-based null
def reduced_simulator(Xr, t, e, ridge=0.0):
    """Fit  PFS ~ Xr  and return a closure that simulates from it under H0.

    The baseline cumulative hazard is Breslow on the SAME uncentred linear predictor the
    simulator then uses, so the two cannot drift apart. Observed follow-up is reused as
    the censoring time, which keeps the censoring distribution and the number at risk
    close to the real trial.
    """
    # Breslow, to match the baseline hazard this routine builds and then inverts. The
    # simulated data are only a data-generating process satisfying H0; the tests applied
    # to them use whatever convention the caller's statistic uses.
    b, V, ok = cox_fit(Xr, t, e, ridge=ridge, ties="breslow")
    if not ok:
        return None
    order = np.argsort(t, kind="mergesort")
    ts, es = t[order], e[order]
    uniq, first_idx, inv = np.unique(ts, return_index=True, return_inverse=True)
    d = np.bincount(inv, weights=es, minlength=len(uniq))
    eta = np.asarray(Xr, float) @ b
    w = np.exp(eta)
    ws = w[order]
    S0 = np.cumsum(ws[::-1])[::-1][first_idx]
    dL0 = d / S0
    H0 = np.cumsum(dL0)
    keep = d > 0
    h_t, h_H = uniq[keep], H0[keep]

    def sim(rng):
        target = -np.log(rng.random(len(w))) / np.maximum(w, 1e-12)
        j = np.searchsorted(h_H, target, side="left")
        tstar = np.where(j < len(h_t), h_t[np.minimum(j, len(h_t) - 1)], np.inf)
        t_new = np.minimum(tstar, t)
        return t_new, (tstar <= t).astype(int)

    return dict(sim=sim, beta=b, n_par=Xr.shape[1])


def pc_design(S, npc):
    """Top-`npc` principal components of the signature matrix, standardised."""
    X = S.T - S.T.mean(0)
    U, sv, Vt = np.linalg.svd(X, full_matrices=False)
    npc = min(npc, len(sv))
    pcs = U[:, :npc] * sv[:npc]
    return pcs / (pcs.std(0) + 1e-9), (sv ** 2 / (sv ** 2).sum()).cumsum()


def perm_p(null, observed, tail):
    v = np.asarray(null, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan
    hit = (v >= observed).sum() if tail == "greater" else (v <= observed).sum()
    return float((hit + 1) / (v.size + 1))
