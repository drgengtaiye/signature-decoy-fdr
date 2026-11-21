"""Tables read by Figure 3 and Supplementary Figure 1 (icb/v13_fig/).

    python scripts/figprep_tables.py [out dir]      (default: icb/v13_fig)

threshold_curve.csv   signatures called and the estimate (co-expressed decoys, ten sets) as the per-unit P threshold
                      varies, folded and directed (Fig. 3)
convergence_sets.csv  the estimate pooled over the first k of the ten decoy sets (Fig. 3s)
eigen_by_cohort.csv   eigenvalues of the 52 x 52 score correlation in each of the 12 checkpoint cohorts (Supplementary
                      Fig. 1)
cohort_similarity.csv correlation between cohorts of the pairwise signature correlations (Supplementary Fig. 1)
family_cluster_call.csv  each signature's family, correlation cluster (average linkage on 1 - |r| of the pooled
                      Fisher-z correlation, cut at 0.5), cluster size and folded response call (Supplementary Fig. 1)
most_shared_genes.csv the 12 genes in most signatures and the signatures holding them (Supplementary Fig. 1)

Inputs: icb/v12_pkg/response_calls.csv.gz, v4_standardised_effects.csv, v3_dissociation.csv, calib_v3/ICB_*_v3.csv,
scores/*_sig.csv.gz, SignatureSets-main/data-raw/IO/*.csv. Nothing is resampled.
"""
import glob, os, sys
from collections import Counter
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cns as K
import figstyle as fs

I = f"{K.D}/icb"
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{I}/v13_fig"
os.makedirs(OUT, exist_ok=True)
RC = pd.read_csv(f"{I}/v12_pkg/response_calls.csv.gz")

# ---------------------------------------------------------------- Figure 3: threshold curve and pooling of decoy sets
TH = np.exp(np.linspace(np.log(0.001), np.log(0.2), 40))
curves, conv = {}, {}
for st in ("fold", "dir"):
    c_ = RC[RC.stat == st]; t_ = c_[c_.kind == "target"]; d_ = c_[(c_.kind == "decoy") & (c_["class"] == "coex")]
    if st == "dir":
        d_ = d_[d_.dirn.to_numpy() == d_.target.map(t_.set_index("target").dirn).to_numpy()]
    Rk, Fk = [], []
    for th in TH:
        r_ = int(((t_.p < th) & (t_.med_delta > 0)).sum())
        e_ = d_.assign(k=(d_.p < th) & (d_.med_delta > 0)).groupby("target").k.mean().sum()
        Rk.append(r_); Fk.append(100 * e_ / r_ if r_ else np.nan)
    curves[st] = (np.array(Rk), np.array(Fk))
    r_ = int(t_.call.sum())                     # the estimate sums the decoy call rate over all 52 targets
    conv[st] = [100 * d_[d_.set <= k].groupby("target").call.mean().sum() / r_ for k in range(1, 11)]
pd.DataFrame({"threshold": TH, "called_fold": curves["fold"][0], "fdr_fold": curves["fold"][1], "called_dir": curves["dir"][0],
              "fdr_dir": curves["dir"][1]}).to_csv(f"{OUT}/threshold_curve.csv", index=False)
pd.DataFrame({"sets_pooled": range(1, 11), "fdr_fold": conv["fold"], "fdr_dir": conv["dir"]}).to_csv(
    f"{OUT}/convergence_sets.csv", index=False)

# ---------------------------------------------------------------- Supplementary Figure 1: the library's correlation
EFF = pd.read_csv(f"{I}/v4_standardised_effects.csv")
X = pd.read_csv(f"{I}/v3_dissociation.csv")[["signature", "family", "n_genes", "method"]].merge(
    EFF[["signature", "z_resp", "z_prog"]], on="signature")
DIRN = pd.concat([pd.read_csv(f, usecols=["unit", "kind", "dirn"]) for f in sorted(glob.glob(f"{I}/calib_v3/ICB_*_v3.csv"))])
DIRN = DIRN[DIRN.kind == "target"].drop_duplicates("unit").rename(columns={"unit": "signature"})[["signature", "dirn"]]
X = X.merge(DIRN, on="signature", how="left"); assert X.dirn.notna().all() and len(X) == 52
CALL = RC[(RC.kind == "target") & (RC.stat == "fold")].set_index("target").call
FAMORD = ["Inflamed/IFN", "Other", "Confounder"]
fam_of = X.set_index("signature").family

cohorts = sorted(os.path.basename(p)[:-11] for p in glob.glob(f"{I}/scores/*_sig.csv.gz"))
cohorts = [c for c in cohorts if c != "ICB_VanDenEnde"]
sigs = sorted(EFF.signature); S_IDX = {s: i for i, s in enumerate(sigs)}; nS = len(sigs)
zsum = np.zeros((nS, nS)); wsum = np.zeros((nS, nS))
for co in cohorts:
    S = pd.read_csv(f"{I}/scores/{co}_sig.csv.gz"); cols = [c for c in S.columns if c in S_IDX]
    r = np.clip(np.nan_to_num(np.corrcoef(S[cols].to_numpy(float), rowvar=False), nan=0.0), -.999999, .999999)
    idx = [S_IDX[c] for c in cols]; w = max(len(S) - 3, 1)
    zsum[np.ix_(idx, idx)] += np.arctanh(r) * w; wsum[np.ix_(idx, idx)] += w
Rm = np.tanh(np.divide(zsum, wsum, out=np.zeros_like(zsum), where=wsum > 0)); np.fill_diagonal(Rm, 1.0)
dist = 1 - np.abs(Rm); np.fill_diagonal(dist, 0); dist = (dist + dist.T) / 2
Z = linkage(squareform(dist, checks=False), method="average"); cl = fcluster(Z, t=.5, criterion="distance")
CLS = pd.Series(cl, index=sigs); csize = CLS.value_counts()
iu = np.triu_indices(nS, 1)

order = []                                      # family, then response effect, as the figure's circular overview
for f in FAMORD:
    order += X[X.family == f].sort_values("z_resp", ascending=False).signature.tolist()

EIG, RV = {}, {}
for co in cohorts:
    S = pd.read_csv(f"{I}/scores/{co}_sig.csv.gz"); cols = [c for c in sigs if c in S.columns]
    A = S[cols].to_numpy(float); A = A[:, np.nanstd(A, 0) > 0]
    r = np.corrcoef(A, rowvar=False)
    EIG[co] = np.sort(np.linalg.eigvalsh(r))[::-1]
    full = pd.DataFrame(np.nan, index=sigs, columns=sigs)
    keep = [c for c, ok in zip(cols, np.nanstd(S[cols].to_numpy(float), 0) > 0) if ok]
    full.loc[keep, keep] = r
    RV[co] = full.to_numpy()[iu]
corder = sorted(cohorts, key=lambda c: (fs.ICB_TUM[c], K.COHORT[c]))
SIM = np.array([[pd.Series(RV[a]).corr(pd.Series(RV[b])) for b in corder] for a in corder])
pd.DataFrame({co: pd.Series(ev) for co, ev in EIG.items()}).rename_axis("component").to_csv(f"{OUT}/eigen_by_cohort.csv")
pd.DataFrame(SIM, index=corder, columns=corder).to_csv(f"{OUT}/cohort_similarity.csv")

genes = {}
for s in X.signature:
    f = f"{I}/SignatureSets-main/data-raw/IO/{s}.csv"
    if os.path.exists(f):
        genes[s] = set(pd.read_csv(f, encoding="latin-1").get("gene_name", pd.Series(dtype=str)).dropna())
cnt_g = Counter(g for s_ in genes.values() for g in s_)
top = sorted(cnt_g.items(), key=lambda kv: (-kv[1], kv[0]))[:12]       # ties broken by name, reproducibly
pd.DataFrame({"signature": order, "family": [fam_of[s] for s in order], "cluster": [int(CLS[s]) for s in order],
              "cluster_size": [int(csize[CLS[s]]) for s in order], "called_response_folded": [bool(CALL[s]) for s in order]}
             ).to_csv(f"{OUT}/family_cluster_call.csv", index=False)
pd.DataFrame([{"gene": g, "n_signatures": n_, "signatures": ";".join(s for s in order if s in genes and g in genes[s])}
              for g, n_ in top]).to_csv(f"{OUT}/most_shared_genes.csv", index=False)
print("written to", OUT)
