"""Treatment-interaction replication across every randomised ICI trial we can access."""
import warnings, functools, subprocess; warnings.filterwarnings('ignore')
print = functools.partial(print, flush=True)
import requests, numpy as np, pandas as pd
from lifelines import CoxPHFitter
from scipy.stats import binomtest

API = "https://www.cbioportal.org/api"
P = "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan"
rng = np.random.default_rng(0); NRAND = 300
S = requests.Session()

# (study, how to define the ICI arm, endpoint)
TRIALS = [
 ("rcc_iatlas_immotion150_2018", "NON_ICI_RX", "not_sunitinib", "PFS"),
 ("paad_iatlas_prince_2022",     "ICI_RX",     "is_nivolumab",  "PFS"),
 ("paad_iatlas_prince_2022",     "ICI_RX",     "is_nivolumab",  "OS"),
]

raw = subprocess.run(["Rscript","-e",
  f'source("{P}/signatures.R");suppressMessages(library(estimate));data(SI_geneset);'
  'SIGS[["ESTIMATE_Immune"]]<-{g<-as.character(unlist(SI_geneset["ImmuneSignature",-1]));g[g!=""]};'
  'SIGS[["ESTIMATE_Stromal"]]<-{g<-as.character(unlist(SI_geneset["StromalSignature",-1]));g[g!=""]};'
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
    for kind in ("SAMPLE", "PATIENT"):
        d = S.get(f"{API}/studies/{st}/clinical-data?clinicalDataType={kind}&projection=DETAILED",
                  timeout=180).json()
        key = "sampleId" if kind == "SAMPLE" else "patientId"
        df = pd.DataFrame([{"k": r[key], "a": r["clinicalAttributeId"], "v": r["value"]} for r in d])
        out[kind] = df.pivot_table(index="k", columns="a", values="v", aggfunc="first")
    smap = {x["sampleId"]: x["patientId"] for x in S.get(f"{API}/studies/{st}/samples", timeout=180).json()}
    C = out["SAMPLE"]
    for c in out["PATIENT"].columns:
        if c not in C.columns:
            C[c] = [out["PATIENT"][c].get(smap.get(s_)) for s_ in C.index]
    return C

allrows = []
for st, armcol, rule, ep in TRIALS:
    C = clinical(st)
    M = pd.read_csv(f"{P}/ici/{st}.csv", index_col=0)
    keep = [s for s in M.columns if s in C.index]
    C = C.loc[keep]; M = M[keep]
    if armcol not in C.columns: print(f"{st}: no {armcol}"); continue
    a = C[armcol].astype(str)
    ici = (~a.str.contains("Sunitinib", na=False)).astype(int) if rule == "not_sunitinib" \
          else a.str.contains("Nivolumab", case=False, na=False).astype(int)
    tcol, scol = (f"{ep}_MONTHS", f"{ep}_STATUS")
    if tcol not in C.columns: print(f"{st}: no {tcol}"); continue
    t = pd.to_numeric(C[tcol], errors="coerce")
    # cBioPortal encodes status as "0:Not Progressed" / "1:Progressed"; take the leading digit
    e = pd.to_numeric(C[scol].astype(str).str.extract(r"^\s*([01])\s*:")[0], errors="coerce")
    ok = t.notna() & (t > 0) & e.notna() & ici.notna()
    e = e.fillna(0).astype(int)
    if ok.sum() < 50 or ici[ok].sum() < 15 or (1-ici[ok]).sum() < 15:
        print(f"{st} {ep}: insufficient ({ok.sum()}, arms {int(ici[ok].sum())}/{int((1-ici[ok]).sum())})"); continue
    Z = ((M.T - M.mean(1)) / M.std(1).replace(0,1)).T
    idx = ok[ok].index

    def inter(v):
        df = pd.DataFrame({"T": t[ok].values, "E": e[ok].values,
                           "s": (v-np.mean(v))/(np.std(v)+1e-9), "arm": ici[ok].values})
        df["sx"] = df.s * df.arm
        try:
            f = CoxPHFitter().fit(df, "T", "E")
            return float(f.params_["sx"]), float(f.summary.loc["sx","p"])
        except Exception: return np.nan, np.nan

    rows = []
    for k, gl in SIGS.items():
        g = [x for x in gl if x in Z.index]
        if len(g) < 3: continue
        b, p = inter(Z.loc[g][idx].mean(0).values)
        rows.append(dict(trial=st, endpoint=ep, signature=k, coef=b, p=p,
                         family="Immune" if k in IMM else ("Control" if k.startswith("CTRL_") else "Other")))
    R = pd.DataFrame(rows).dropna()
    pool = list(Z.index)
    rp = [inter(Z.loc[rng.choice(pool,48,replace=False)][idx].mean(0).values)[1] for _ in range(NRAND)]
    rp = np.array([x for x in rp if np.isfinite(x)]); rr = (rp<.05).mean()
    imm = R[R.family=="Immune"]
    hit = int((imm.p<.05).sum())
    bt = binomtest(hit, len(imm), max(rr,1e-3), alternative="greater").pvalue
    print(f"\n=== {st}  [{ep}]  n={int(ok.sum())}  ICI={int(ici[ok].sum())} ctrl={int((1-ici[ok]).sum())} "
          f"events={int(e[ok].sum())} ===")
    print(f"  immune signatures p<0.05 : {hit}/{len(imm)} ({hit/len(imm)*100:.0f}%)")
    print(f"  random 48-gene sets      : {rr*100:.1f}%")
    print(f"  binomial vs random       : p = {bt:.5f}")
    print(f"  immune coef sign: {int((imm.coef<0).sum())}/{len(imm)} negative (= benefit in ICI arm)")
    ctl = R[R.family=="Control"]
    if len(ctl): print(f"  controls  coef sign: {int((ctl.coef<0).sum())}/{len(ctl)} negative")
    print(R.nsmallest(6,"p")[["signature","coef","p","family"]].round(4).to_string(index=False))
    allrows.append(R)

A = pd.concat(allrows, ignore_index=True)
A.to_csv(f"{P}/interaction_all_trials.csv", index=False)
print("\n=== combined: immune vs control interaction direction, all trials ===")
for f in ["Immune","Control"]:
    g = A[A.family==f]
    print(f"  {f:8s} n={len(g):3d}  median coef={g.coef.median():+.3f}  "
          f"{(g.coef<0).mean()*100:.0f}% negative  p<0.05 in {(g.p<.05).mean()*100:.0f}%")
