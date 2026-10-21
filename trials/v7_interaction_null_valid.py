"""Nulls that actually correspond to  H0: beta(score x arm) = 0, main effects free.

v6 ran four nulls and grouped three of them as "retaining the main effects". Two of those
do no such thing:

  restricted  permuting arm labels inside strata of the score still breaks the TREATMENT
              main effect whenever beta_arm != 0 - the outcomes stay attached to the
              original assignment, so the null population has no treatment effect;
  scoreperm   permuting the score inside arm breaks the SCORE main effect whenever
              beta_score != 0.

Each therefore tests a different, stronger hypothesis than the one the paper argues for,
and neither belongs in a sentence with the model-based bootstrap. This script replaces
them with nulls that hold both main effects and break only the interaction, and reports
what each one actually preserves.

Nulls
  persig      per-signature model-based bootstrap. For signature j, fit PFS ~ s_j + arm,
              simulate from it, refit with the interaction. This is the exact null for
              that signature's interaction coefficient with BOTH of its own main effects
              intact. It is per-signature only: one replicate is not simultaneously null
              for the other signatures, so it cannot calibrate the count statistic.
  pcK         joint model-based bootstrap simulating from PFS ~ PC1..PCK + arm. One
              replicate is null for every signature at once, so it does calibrate the
              count and the pooled contrast. K controls how much of the score space keeps
              its main effect: K=20 is the full rank of the signature matrix, i.e. no
              subspace restriction at all. This is the subspace-sensitivity ladder.
  wild        multiplier bootstrap of the efficient interaction score. Both main effects
              are pinned at their reduced-model estimates and nothing is permuted; the
              per-patient efficient score residuals are re-weighted by iid Rademacher
              multipliers. The SAME multipliers are applied to all signatures within a
              trial, so the between-signature correlation is preserved exactly. Free of
              any parametric assumption about the baseline hazard.
  freedman    permute the interaction regressor after projecting out score and arm, with
              the same permutation shared across signatures. Keeps the reduced design
              bolted to the outcomes; assumes the residualised regressor is exchangeable.

What each preserves is recorded in the output table rather than asserted in prose.

Outputs
  pancan/v7_interaction_null_valid.csv   observed, null quantiles and P per statistic
  pancan/v7_persignature_exact.csv       per-signature exact bootstrap P
  pancan/v7_null_properties.csv          what every null does and does not preserve
"""
import sys, time, functools
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
from v7_lib import (P, build, signatures, cox_fit, cox_interaction, cox_arm_slope,
                    efficient_interaction_score, reduced_simulator, pc_design, perm_p)

NREP  = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
NJOBS = 48
PC_LADDER = [3, 4, 9, 15, 20]
TRIALS = ("IMmotion150", "PRINCE")


# ------------------------------------------------------------------ statistics
def pass_stats(S, imm, arm, t, e):
    """Count of immune signatures with interaction P<0.05 and beta<0, plus the 17
    per-signature arm-slope differences that feed the pooled contrast."""
    hits = 0
    diffs = []
    for i in range(S.shape[0]):
        b, p = cox_interaction(S[i], arm, t, e)
        if imm[i]:
            if np.isfinite(p) and p < .05 and b < 0:
                hits += 1
            diffs.append(cox_arm_slope(S[i], arm, t, e, 1) - cox_arm_slope(S[i], arm, t, e, 0))
    return hits, np.array(diffs, float)


# ------------------------------------------------------------------ main
def main():
    SIGS = signatures()
    D = {tr: build(tr, SIGS) for tr in TRIALS}
    DALL = {tr: build(tr, SIGS, family_filter=False) for tr in TRIALS}
    rows, persig, props = [], [], []

    obs_hits, obs_diffs, E = {}, {}, {}
    for tr in TRIALS:
        d = D[tr]
        imm = d["fams"] == "Immune"
        h, df = pass_stats(d["S"], imm, d["arm"], d["t"], d["e"])
        obs_hits[tr], obs_diffs[tr] = h, df
        E[tr] = [efficient_interaction_score(d["S"][i], d["arm"], d["t"], d["e"])
                 for i in range(d["S"].shape[0])]
        print(f"{tr}: n={d['n']} events={d['e'].sum()} immune={imm.sum()} "
              f"observed hits {h}/{imm.sum()} median arm-slope diff {np.median(df):+.4f}")
    obs_pooled = float(np.median(np.concatenate([obs_diffs[tr] for tr in TRIALS])))
    print(f"pooled 34-test observed median log-HR difference {obs_pooled:+.4f}\n")

    # ---------------------------------------------------------- 1. per-signature exact
    print("=== per-signature model-based bootstrap (exact null for each signature) ===")
    for tr in TRIALS:
        d = DALL[tr]
        S, arm, t, e = d["S"], d["arm"], d["t"], d["e"]

        def one_sig(i):
            sc = S[i]
            b_obs, p_obs = cox_interaction(sc, arm, t, e)
            sim = reduced_simulator(np.column_stack([sc, arm]), t, e)
            if sim is None or not np.isfinite(b_obs):
                return None
            rng = np.random.default_rng(1000 + i)
            nb = []
            for _ in range(NREP):
                tn, en = sim["sim"](rng)
                nb.append(cox_interaction(sc, arm, tn, en)[0])
            nb = np.array(nb, float)
            return dict(trial=tr, signature=d["names"][i], family=d["fams"][i],
                        beta=b_obs, wald_p=p_obs,
                        boot_p_two_sided=perm_p(np.abs(nb[np.isfinite(nb)]), abs(b_obs), "greater"),
                        boot_p_checkpoint=perm_p(nb, b_obs, "less"),
                        null_sd=float(np.nanstd(nb)), n_rep=int(np.isfinite(nb).sum()))

        out = Parallel(n_jobs=NJOBS)(delayed(one_sig)(i) for i in range(S.shape[0]))
        out = [o for o in out if o]
        persig += out
        T = pd.DataFrame(out)
        immT = T[T.family == "Immune"]
        nh = int(((immT.boot_p_checkpoint < .05) & (immT.beta < 0)).sum())
        print(f"  {tr}: {len(T)} signatures tested; immune hits under each signature's own "
              f"exact null = {nh}/{len(immT)}  (Wald-based observed {obs_hits[tr]}/{len(immT)})")
        rows.append(dict(trial=tr, method="persig_exact_bootstrap", statistic="hit_immune",
                         observed=nh, null_median=np.nan, null_p05=np.nan, null_p95=np.nan,
                         n_rep=NREP, perm_p=np.nan,
                         note="count of signatures significant against their OWN exact null; "
                              "not a calibrated count statistic"))

    # ---------------------------------------------------------- 2. joint model-based, PC ladder
    print("\n=== joint model-based bootstrap across the PC ladder ===")
    pooled_null = {}
    for K in PC_LADDER:
        per_trial_diffs = {}
        for tr in TRIALS:
            d = D[tr]
            S, arm, t, e = d["S"], d["arm"], d["t"], d["e"]
            imm = d["fams"] == "Immune"
            pcs, cum = pc_design(S, K)
            sim = reduced_simulator(np.column_stack([pcs, arm]), t, e)
            if sim is None:
                print(f"  {tr} K={K}: reduced model failed"); continue

            def one(seed):
                rng = np.random.default_rng(seed)
                tn, en = sim["sim"](rng)
                return pass_stats(S, imm, arm, tn, en)

            t0 = time.time()
            seeds = np.random.SeedSequence(90000 + K).spawn(NREP)
            out = Parallel(n_jobs=NJOBS)(delayed(one)(s) for s in seeds)
            hits = np.array([o[0] for o in out], float)
            dif = np.vstack([o[1] for o in out])
            per_trial_diffs[tr] = dif
            p = perm_p(hits, obs_hits[tr], "greater")
            var = cum[min(K, len(cum)) - 1]
            print(f"  {tr} K={K:2d} (score variance retained {var:.3f}): "
                  f"hits obs {obs_hits[tr]} null med {np.median(hits):.1f} "
                  f"95th {np.quantile(hits,.95):.1f}  P={p:.4f}   [{(time.time()-t0)/60:.1f} min]")
            rows.append(dict(trial=tr, method=f"joint_parametric_pc{K}", statistic="hit_immune",
                             observed=obs_hits[tr], null_median=float(np.median(hits)),
                             null_p05=float(np.quantile(hits, .05)),
                             null_p95=float(np.quantile(hits, .95)), n_rep=NREP, perm_p=p,
                             note=f"score variance retained {var:.3f}"))
        if len(per_trial_diffs) == len(TRIALS):
            med = np.median(np.hstack([per_trial_diffs[tr] for tr in TRIALS]), axis=1)
            p = perm_p(med, obs_pooled, "less")
            pooled_null[f"joint_parametric_pc{K}"] = med
            print(f"  pooled 34-test median log-HR difference: obs {obs_pooled:+.4f} "
                  f"null med {np.median(med):+.4f} 5th {np.quantile(med,.05):+.4f}  P={p:.4f}")
            rows.append(dict(trial="pooled", method=f"joint_parametric_pc{K}",
                             statistic="median_log_diff_34tests", observed=obs_pooled,
                             null_median=float(np.median(med)),
                             null_p05=float(np.quantile(med, .05)),
                             null_p95=float(np.quantile(med, .95)), n_rep=NREP, perm_p=p,
                             note="arm-stratified refits on simulated outcomes"))

    # ---------------------------------------------------------- 3. multiplier bootstrap
    print("\n=== multiplier (wild) bootstrap of the efficient interaction score ===")
    rng = np.random.default_rng(4242)
    b1_null, z_obs_all, b1_obs_all = {}, {}, {}
    for tr in TRIALS:
        d = D[tr]
        imm = d["fams"] == "Immune"
        ok = [i for i in range(d["S"].shape[0]) if E[tr][i] is not None]
        Umat = np.vstack([E[tr][i]["u"] for i in ok])                 # nsig x n
        Ieff = np.array([E[tr][i]["Ieff"] for i in ok])
        ss = np.sqrt((Umat ** 2).sum(1))
        G = rng.choice([-1.0, 1.0], size=(NREP, d["n"]))              # SHARED across signatures
        Tn = G @ Umat.T                                               # NREP x nsig
        zn = Tn / ss
        b1n = Tn / Ieff
        z_obs = np.array([E[tr][i]["z"] for i in ok])
        b1_obs = np.array([E[tr][i]["b1"] for i in ok])
        immok = imm[ok]
        obs_h = int((z_obs[immok] < -1.96).sum())
        null_h = (zn[:, immok] < -1.96).sum(1).astype(float)
        p = perm_p(null_h, obs_h, "greater")
        print(f"  {tr}: robust-z hits obs {obs_h}/{immok.sum()} "
              f"(Wald-based {obs_hits[tr]}/{immok.sum()}) null med {np.median(null_h):.1f} "
              f"95th {np.quantile(null_h,.95):.1f}  P={p:.4f}")
        rows.append(dict(trial=tr, method="wild_multiplier", statistic="hit_immune",
                         observed=obs_h, null_median=float(np.median(null_h)),
                         null_p05=float(np.quantile(null_h, .05)),
                         null_p95=float(np.quantile(null_h, .95)), n_rep=NREP, perm_p=p,
                         note="hits defined by robust score z < -1.96, same rule observed and null"))
        b1_null[tr] = b1n[:, immok]
        b1_obs_all[tr] = b1_obs[immok]
    med_obs = float(np.median(np.concatenate([b1_obs_all[tr] for tr in TRIALS])))
    med_null = np.median(np.hstack([b1_null[tr] for tr in TRIALS]), axis=1)
    p = perm_p(med_null, med_obs, "less")
    pooled_null["wild_multiplier"] = med_null
    print(f"  pooled 34-test median one-step interaction coefficient: obs {med_obs:+.4f} "
          f"null med {np.median(med_null):+.4f} 5th {np.quantile(med_null,.05):+.4f}  P={p:.4f}")
    rows.append(dict(trial="pooled", method="wild_multiplier",
                     statistic="median_interaction_coef_34tests", observed=med_obs,
                     null_median=float(np.median(med_null)),
                     null_p05=float(np.quantile(med_null, .05)),
                     null_p95=float(np.quantile(med_null, .95)), n_rep=NREP, perm_p=p,
                     note="one-step interaction coefficient; tracks the arm-slope difference at r=0.99"))

    # ---------------------------------------------------------- 4. Freedman-Lane
    print("\n=== Freedman-Lane permutation of the residualised interaction regressor ===")
    fl_null = {}
    for tr in TRIALS:
        d = D[tr]
        S, arm, t, e, n = d["S"], d["arm"], d["t"], d["e"], d["n"]
        imm = d["fams"] == "Immune"
        R = []
        for i in range(S.shape[0]):
            sc = S[i]
            Dm = np.column_stack([sc, arm, np.ones(n)])
            v = sc * arm
            R.append(v - Dm @ np.linalg.lstsq(Dm, v, rcond=None)[0])
        R = np.vstack(R)

        def one(seed):
            r = np.random.default_rng(seed)
            pi = r.permutation(n)                       # SHARED across signatures
            hits, coefs = 0, []
            for i in range(S.shape[0]):
                b, V, ok = cox_fit(np.column_stack([S[i], arm, R[i][pi]]), t, e)
                if not ok:
                    if imm[i]:
                        coefs.append(np.nan)
                    continue
                z = b[2] / np.sqrt(V[2, 2])
                if imm[i]:
                    if abs(z) > 1.96 and b[2] < 0:
                        hits += 1
                    coefs.append(b[2])
            return hits, np.array(coefs, float)

        t0 = time.time()
        out = Parallel(n_jobs=NJOBS)(delayed(one)(s) for s in
                                     np.random.SeedSequence(555).spawn(NREP))
        hits = np.array([o[0] for o in out], float)
        coefs = np.vstack([o[1] for o in out])
        p = perm_p(hits, obs_hits[tr], "greater")
        print(f"  {tr}: hits obs {obs_hits[tr]} null med {np.median(hits):.1f} "
              f"95th {np.quantile(hits,.95):.1f}  P={p:.4f}  [{(time.time()-t0)/60:.1f} min]")
        rows.append(dict(trial=tr, method="freedman_lane", statistic="hit_immune",
                         observed=obs_hits[tr], null_median=float(np.median(hits)),
                         null_p05=float(np.quantile(hits, .05)),
                         null_p95=float(np.quantile(hits, .95)), n_rep=NREP, perm_p=p,
                         note="shared permutation across signatures"))
        fl_null[tr] = coefs
    obs_int = np.concatenate([
        [cox_interaction(D[tr]["S"][i], D[tr]["arm"], D[tr]["t"], D[tr]["e"])[0]
         for i in range(D[tr]["S"].shape[0]) if D[tr]["fams"][i] == "Immune"]
        for tr in TRIALS])
    med_obs_i = float(np.median(obs_int))
    med_fl = np.nanmedian(np.hstack([fl_null[tr] for tr in TRIALS]), axis=1)
    p = perm_p(med_fl, med_obs_i, "less")
    pooled_null["freedman_lane"] = med_fl
    print(f"  pooled 34-test median interaction coefficient: obs {med_obs_i:+.4f} "
          f"null med {np.median(med_fl):+.4f} 5th {np.quantile(med_fl,.05):+.4f}  P={p:.4f}")
    rows.append(dict(trial="pooled", method="freedman_lane",
                     statistic="median_interaction_coef_34tests", observed=med_obs_i,
                     null_median=float(np.median(med_fl)),
                     null_p05=float(np.quantile(med_fl, .05)),
                     null_p95=float(np.quantile(med_fl, .95)), n_rep=NREP, perm_p=p, note=""))

    # ---------------------------------------------------------- what each null preserves
    props = [
        dict(method="arm_label_permutation", targets_interaction_null=False,
             keeps_treatment_main_effect=False, keeps_score_main_effect=True,
             keeps_signature_correlation=True,
             tests="sharp null of no treatment effect of any kind",
             valid_for_H0_interaction=False),
        dict(method="restricted_arm_permutation", targets_interaction_null=False,
             keeps_treatment_main_effect=False, keeps_score_main_effect=True,
             keeps_signature_correlation=True,
             tests="no treatment effect within strata of the score",
             valid_for_H0_interaction=False),
        dict(method="score_within_arm_permutation", targets_interaction_null=False,
             keeps_treatment_main_effect=True, keeps_score_main_effect=False,
             keeps_signature_correlation=True,
             tests="no association between score and outcome of any kind",
             valid_for_H0_interaction=False),
        dict(method="joint_parametric_pcK", targets_interaction_null=True,
             keeps_treatment_main_effect=True, keeps_score_main_effect="within K-dim subspace",
             keeps_signature_correlation=True,
             tests="beta(score x arm) = 0", valid_for_H0_interaction=True),
        dict(method="persig_exact_bootstrap", targets_interaction_null=True,
             keeps_treatment_main_effect=True, keeps_score_main_effect=True,
             keeps_signature_correlation=False,
             tests="beta(score x arm) = 0 for one signature", valid_for_H0_interaction=True),
        dict(method="wild_multiplier", targets_interaction_null=True,
             keeps_treatment_main_effect=True, keeps_score_main_effect=True,
             keeps_signature_correlation=True,
             tests="beta(score x arm) = 0, no baseline-hazard model",
             valid_for_H0_interaction=True),
        dict(method="freedman_lane", targets_interaction_null=True,
             keeps_treatment_main_effect=True, keeps_score_main_effect=True,
             keeps_signature_correlation=True,
             tests="beta(score x arm) = 0, assumes residualised regressor exchangeable",
             valid_for_H0_interaction=True),
    ]

    pd.DataFrame(rows).to_csv(f"{P}/v7_interaction_null_valid.csv", index=False)
    pd.DataFrame(persig).to_csv(f"{P}/v7_persignature_exact.csv", index=False)
    pd.DataFrame(props).to_csv(f"{P}/v7_null_properties.csv", index=False)
    print("\n" + pd.DataFrame(rows).to_string(index=False))
    print(f"\nwrote v7_interaction_null_valid.csv, v7_persignature_exact.csv, "
          f"v7_null_properties.csv")


if __name__ == "__main__":
    main()
