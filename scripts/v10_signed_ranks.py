"""Analysis 1 of validation/PLAN_calibration_round5.md, for the label-permutation scripts: add the signed null
sets of the four +1/-1 weighted-mean signatures to the stored rank matrices.

    python v10_signed_ranks.py   -> icb/scores/<cohort>_ranks.npz (originals kept in icb/_pre_v10/scores/)

perm_inference.py and v5_perm_family.py read each signature's null from scores/<cohort>_ranks.npz through
its rung key. Here the 200 signed sets of each affected signature (v10_signed_null.R response, written to
scores/<cohort>_signed.csv.gz) are ranked, appended under the key 'wmean_signed|<signature>', and that
signature's rung is pointed at them. Nothing else in the file changes, and both scripts run unchanged.
"""
import glob, os, shutil, numpy as np, pandas as pd
from scipy.stats import rankdata

D = "/path/to/analysis_root/icb"
os.makedirs(f"{D}/_pre_v10/scores", exist_ok=True)
for f in sorted(glob.glob(f"{D}/scores/*_signed.csv.gz")):
    co = os.path.basename(f)[:-len("_signed.csv.gz")]
    npz, bak = f"{D}/scores/{co}_ranks.npz", f"{D}/_pre_v10/scores/{co}_ranks.npz"
    if not os.path.exists(bak):
        shutil.copy2(npz, bak)
    z = dict(np.load(bak, allow_pickle=True))
    order = pd.read_csv(f"{D}/scores/{co}_sig.csv.gz", usecols=["sample"])["sample"].astype(str).to_numpy()
    S = pd.read_csv(f).set_index("sample"); S.index = S.index.astype(str)
    S = S.loc[order]                                     # the patients in the order the ranks were built
    keys = np.array([c.rsplit("|", 1)[0] for c in S.columns])
    sig_names = list(z["sig_names"]); sig_rung = z["sig_rung"].astype(object).copy()
    for k in np.unique(keys):
        sig_rung[sig_names.index(k.split("|", 1)[1])] = k
    z["ranks_n"] = np.hstack([z["ranks_n"], rankdata(S.to_numpy(float), axis=0)])
    z["rung_key"] = np.concatenate([z["rung_key"].astype(object), keys.astype(object)])
    z["sig_rung"] = sig_rung
    np.savez_compressed(npz, **z)
    print(f"{co}: +{S.shape[1]} signed sets for {', '.join(sorted(set(k.split('|')[1] for k in keys)))}")
