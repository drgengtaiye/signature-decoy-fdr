"""NOT pre-registered (post hoc, labelled as such): number of the 17 immune gene sets whose
interaction points in the checkpoint-favoring direction (coefficient > 0), and their median
coefficient, referred to the same no-interaction simulation as the registered count."""
import sys, numpy as np, pandas as pd
from joblib import Parallel, delayed
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/validation"); sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan")
from ispy2_interaction import logit_fit, interaction, V
from v7_lib import signatures, IMM, pc_design
out = []
for contrast in ["ISPY2_pembro", "ISPY2_durva"]:
    X = pd.read_csv(f"{V}/cohorts/{contrast}.expr.tsv.gz", sep="\t", index_col=0); X.columns = [str(c) for c in X.columns]
    C = pd.read_csv(f"{V}/cohorts/{contrast}.clin.tsv", sep="\t"); C["sample"] = C["sample"].astype(str); X = X[C["sample"]]
    arm, y = C.arm.to_numpy(float), C.y.to_numpy(float)
    Z = ((X.T - X.mean(1)) / X.std(1).replace(0, 1)).T; rows, imm = [], []
    for nm, gl in signatures().items():
        g = [x for x in gl if x in Z.index]
        if len(g) < 3: continue
        fam = "Immune" if nm in IMM else "Control" if nm.startswith("CTRL_") else "Other"
        if fam == "Other": continue
        rows.append(Z.loc[g].mean(0).values); imm.append(fam == "Immune")
    S = np.vstack(rows); S = (S - S.mean(1, keepdims=True)) / (S.std(1, keepdims=True) + 1e-9); imm = np.array(imm)
    def stats(yv):
        b = np.array([interaction(S[i], arm, yv)[0] for i in np.where(imm)[0]])
        return int((b > 0).sum()), float(np.median(b))
    obs = stats(y)
    for K in [3, 4, 9, 15, 20]:
        pcs, _ = pc_design(S, K); _, _, p0 = logit_fit(np.column_stack([np.ones(len(y)), pcs, arm]), y)
        def one(r):
            rng = np.random.default_rng(20_000 * K + r); return stats((rng.random(len(y)) < p0).astype(float))
        null = np.array(Parallel(n_jobs=48)(delayed(one)(r) for r in range(2000)))
        out.append(dict(contrast=contrast, K=K, positive=obs[0], of=int(imm.sum()), median_beta=obs[1],
                        P_positive=float((np.sum(null[:, 0] >= obs[0]) + 1) / 2001),
                        P_median_beta=float((np.sum(null[:, 1] >= obs[1]) + 1) / 2001),
                        null_positive_median=float(np.median(null[:, 0]))))
        print(out[-1], flush=True)
pd.DataFrame(out).to_csv(f"{V}/ispy2_direction_posthoc.csv", index=False)
