"""Do immune signatures predict benefit from checkpoint blockade, or from treatment?

The interaction test says the signature effect differs between arms. It does not say what
the effect is inside each arm, and that is the question a clinician asks: if these markers
also rank patients in the control arm, they are prognostic markers of a treatable tumour
rather than predictors of checkpoint benefit. Both randomised trials have a
non-checkpoint arm, so the model can simply be fitted twice.

This replaces the cell-line drug analysis from the earlier gastric work, in which none of
11 compounds survived FDR correction and the strongest association ran opposite to the
one reported.
"""
import warnings, functools, subprocess; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import requests, numpy as np, pandas as pd
from lifelines import CoxPHFitter
from scipy.stats import wilcoxon, mannwhitneyu

API = "https://www.cbioportal.org/api"
P = "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan"
S = requests.Session()
TRIALS = [("rcc_iatlas_immotion150_2018", "NON_ICI_RX", "not_sunitinib", "PFS",
           "IMmotion150", "atezolizumab", "sunitinib"),
          ("paad_iatlas_prince_2022", "ICI_RX", "is_nivolumab", "PFS",
           "PRINCE", "nivolumab", "non-checkpoint")]

raw = subprocess.run(["Rscript","-e",
  f'source("{P}/signatures.R");'
  'SIGS<-SIGS[!sapply(SIGS,is.null)];'
  'for(n in names(SIGS)) cat(n,"\\t",paste(SIGS[[n]],collapse=","),"\\n",sep="")'],
  capture_output=True, text=True).stdout
SIGS = {k.strip(): [x for x in v.split(",") if x]
        for k, v in (l.split("\t",1) for l in raw.strip().split("\n") if "\t" in l)}
IMM = {'T_cells','CD8_Tcells','Cytotoxic','NK_cells','B_cells','Plasma_cells','Exhaustion',
       'Checkpoints','GEP18_TIS (Ayers 2017)','IFNG_6gene (Ayers 2017)','TLS (Cabrita 2020)',
       'ESTIMATE_Immune','barriernet_Immune','APM','Treg','DC','Chemokine12 (Messina 2012)',
       'TLS_9gene (Prabhakaran)','CYT (Rooney 2015)'}

def clinical(st):
    out = {}
    for kind in ("SAMPLE","PATIENT"):
        d = S.get(f"{API}/studies/{st}/clinical-data?clinicalDataType={kind}&projection=DETAILED",
                  timeout=180).json()
        key = "sampleId" if kind=="SAMPLE" else "patientId"
        df = pd.DataFrame([{"k":r[key],"a":r["clinicalAttributeId"],"v":r["value"]} for r in d])
        out[kind] = df.pivot_table(index="k", columns="a", values="v", aggfunc="first")
    smap = {x["sampleId"]:x["patientId"] for x in S.get(f"{API}/studies/{st}/samples",timeout=180).json()}
    C = out["SAMPLE"]
    for c in out["PATIENT"].columns:
        if c not in C.columns: C[c] = [out["PATIENT"][c].get(smap.get(s_)) for s_ in C.index]
    return C

rows = []
for st, armcol, rule, ep, name, ici_lab, ctl_lab in TRIALS:
    C = clinical(st); M = pd.read_csv(f"{P}/ici/{st}.csv", index_col=0)
    keep = [s for s in M.columns if s in C.index]; C = C.loc[keep]; M = M[keep]
    a = C[armcol].astype(str)
    ici = (~a.str.contains("Sunitinib", na=False)).astype(int) if rule=="not_sunitinib" \
          else a.str.contains("Nivolumab", case=False, na=False).astype(int)
    t = pd.to_numeric(C[f"{ep}_MONTHS"], errors="coerce")
    e = pd.to_numeric(C[f"{ep}_STATUS"].astype(str).str.extract(r"^\s*([01])\s*:")[0],
                      errors="coerce").fillna(0).astype(int)
    ok = t.notna() & (t>0) & ici.notna()
    Z = ((M.T - M.mean(1))/M.std(1).replace(0,1)).T
    print(f"\n=== {name} [{ep}]  {ici_lab} n={int(ici[ok].sum())}  "
          f"{ctl_lab} n={int((1-ici[ok]).sum())} ===")
    for k, gl in SIGS.items():
        g = [x for x in gl if x in Z.index]
        if len(g) < 3: continue
        v = Z.loc[g].mean(0)
        r = dict(trial=name, signature=k,
                 family="Immune" if k in IMM else ("Control" if k.startswith("CTRL_") else "Other"))
        for arm, lab in [(1,"ici"), (0,"ctl")]:
            m = ok & (ici==arm)
            df = pd.DataFrame({"T":t[m].values, "E":e[m].values,
                               "s":((v[m]-v[m].mean())/(v[m].std()+1e-9)).values})
            try:
                f = CoxPHFitter().fit(df,"T","E")
                r[f"hr_{lab}"] = float(np.exp(f.params_["s"])); r[f"p_{lab}"] = float(f.summary.loc["s","p"])
            except Exception: r[f"hr_{lab}"], r[f"p_{lab}"] = np.nan, np.nan
        rows.append(r)
R = pd.DataFrame(rows).dropna(subset=["hr_ici","hr_ctl"])
R.to_csv(f"{P}/arm_specific.csv", index=False)

print("\n" + "="*72)
for name in R.trial.unique():
    for fam in ["Immune","Control"]:
        g = R[(R.trial==name)&(R.family==fam)]
        if not len(g): continue
        print(f"{name:12s} {fam:8s} n={len(g):2d}   "
              f"checkpoint arm: median HR {g.hr_ici.median():.2f}, {int((g.p_ici<.05).sum())} at P<0.05   |   "
              f"control arm: median HR {g.hr_ctl.median():.2f}, {int((g.p_ctl<.05).sum())} at P<0.05")
imm = R[R.family=="Immune"]
print(f"\nacross both trials, immune signatures (n={len(imm)}):")
print(f"  benefit in the checkpoint arm (HR<1):  {int((imm.hr_ici<1).sum())}/{len(imm)}")
print(f"  benefit in the control arm (HR<1):     {int((imm.hr_ctl<1).sum())}/{len(imm)}")
w = wilcoxon(np.log(imm.hr_ici), np.log(imm.hr_ctl))
print(f"  paired difference in log HR:  P = {w.pvalue:.4f}")
