"""Build each registered validation cohort as a log-scale genes x samples matrix and a label file.

Labels: 1 = CR or PR (responder), 0 = SD or PD (non-responder); one pre-treatment sample per
patient. Chemo-immunotherapy cohorts use the registered endpoint (major pathological response
for GSE207422, CR/PR against SD for GSE235919). Nothing here looks at any signature.
"""
import gzip, io, re, glob, os, numpy as np, pandas as pd
V = "/groups/adv2105_gp/yichen/Yi/com/gastric/validation"; RAW = f"{V}/raw"; OUT = f"{V}/cohorts"
os.makedirs(OUT, exist_ok=True)


def soft(acc):
    out = []
    for r in open(f"{V}/{acc}_gsm.soft", encoding="utf-8", errors="ignore").read().split("^SAMPLE = ")[1:]:
        d = {"gsm": r.split("\n", 1)[0].strip()}
        for l in r.splitlines():
            if l.startswith("!Sample_title"): d["title"] = l.split(" = ", 1)[1]
            if l.startswith("!Sample_description"): d.setdefault("desc", l.split(" = ", 1)[1])
            if l.startswith("!Sample_characteristics_ch1 = ") and ": " in l:
                k, v = l.split(" = ", 1)[1].split(": ", 1); d[k] = v
        out.append(d)
    return pd.DataFrame(out)


def recist(v):
    v = str(v).strip().upper()
    if v in ("CR", "PR", "COMPLETE RESPONSE", "PARTIAL RESPONSE"): return 1
    if v in ("SD", "PD", "STABLE DISEASE", "PROGRESSIVE DISEASE"): return 0
    return np.nan


def save(name, M, y, group, note):
    M = M.loc[:, y.index]; M = M[~M.index.duplicated()]
    empty = M.columns[M.isna().all()]                 # deposited without any expression value
    if len(empty): note += f"; {len(empty)} sample(s) without expression dropped"
    M = M.drop(columns=empty); y = y.drop(index=empty)
    M.to_csv(f"{OUT}/{name}.expr.tsv.gz", sep="\t", compression="gzip")
    pd.DataFrame({"sample": y.index, "y": y.values.astype(int)}).to_csv(f"{OUT}/{name}.clin.tsv", sep="\t", index=False)
    print(f"{name:12s} {group:8s} genes {M.shape[0]:6d} samples {M.shape[1]:4d} responders {int(y.sum()):3d} non {int((1 - y).sum()):3d} | {note}")
    return dict(cohort=name, group=group, genes=M.shape[0], n=M.shape[1], n_resp=int(y.sum()), note=note)


def median_of_ratios(C):
    C = C.clip(lower=0); L = np.log(C.replace(0, np.nan)); ref = L.mean(axis=1)
    ok = ref.notna() & L.notna().all(axis=1)
    sf = np.exp((L[ok].sub(ref[ok], axis=0)).median(axis=0))
    return np.log2(C.div(sf, axis=1) + 1)


rows = []
# ---- primary: GSE176307, metastatic urothelial, Salmon gene TPM
S = soft("GSE176307"); S["y"] = S["io.response"].map(recist); S = S.dropna(subset=["y"])
M = pd.read_csv(f"{RAW}/GSE176307_salmon_tpm_gene.matrix.tsv.gz", sep="\t", index_col=0)
S = S[S.desc.isin(M.columns)]
y = pd.Series(S.y.values, index=S.desc.values)
rows.append(save("GSE176307", np.log2(M + 1), y, "primary", "Salmon TPM, log2(TPM+1)"))

# ---- panel: GSE161537 (NSCLC, HTG), decimal commas, patients as columns
S = soft("GSE161537"); S["y"] = S["best response on immunotherapy (recist)"].map(recist); S = S.dropna(subset=["y"])
M = pd.read_csv(f"{RAW}/GSE161537_nivobio_log2cpm.csv.gz", sep=";", decimal=",", index_col=0)
M.columns = [str(c).strip() for c in M.columns]
S = S[S["patient id"].isin(M.columns)]
rows.append(save("GSE161537", M, pd.Series(S.y.values, index=S["patient id"].values), "panel", "HTG EdgeSeq, log2 CPM"))

# ---- panel: GSE159067 (HNSCC, HTG)
S = soft("GSE159067"); S["y"] = S["best response on immunotherapy (recist)"].map(recist); S = S.dropna(subset=["y"])
M = pd.read_csv(f"{RAW}/GSE159067_IHN_log2cpm_data.txt.gz", sep="\t", index_col=0)
M = M[~M.index.str.startswith(("NEG_CTRL", "POS_CTRL", "GDNA", "ERCC"))]
S = S[S.desc.isin(M.columns)]
rows.append(save("GSE159067", M, pd.Series(S.y.values, index=S.desc.values), "panel", "HTG EdgeSeq, log2 CPM"))

# ---- panel: GSE303300 (cervical, HTG), baseline sample, response assessed before chemoradiation
S = soft("GSE303300"); S = S[S["sample timepoint"] == "baseline"].copy()
S["y"] = S["recist within 1week before rct"].map(recist); S = S.dropna(subset=["y"])
X = pd.read_excel(f"{RAW}/GSE303300_htg_log2cpm_ano.xlsx").set_index("Sample Name")
S["lib"] = S.desc.str.replace("Library name: ", "", regex=False).str.strip()
S = S[S.lib.isin(X.columns)]
rows.append(save("GSE303300", X, pd.Series(S.y.values, index=S.lib.values), "panel", "HTG EdgeSeq, log2 CPM, baseline"))

# ---- panel: GSE215868 (melanoma, NanoString IO 360 RCC): positive-control then housekeeping scaling
S = soft("GSE215868"); S["y"] = S["bor"].map(recist); S = S.dropna(subset=["y"])
cnt, hk, pos = {}, {}, {}
for f in sorted(glob.glob(f"{RAW}/GSE215868/*.txt.gz")):
    gsm = os.path.basename(f).split("_")[0]
    txt = gzip.open(f, "rt", errors="ignore").read()
    cs = txt.split("<Code_Summary>")[1].split("</Code_Summary>")[0].strip().splitlines()
    t = pd.read_csv(io.StringIO("\n".join(cs)))
    cnt[gsm] = t[t.CodeClass == "Endogenous"].set_index("Name").Count
    hk[gsm] = t[t.CodeClass == "Housekeeping"].set_index("Name").Count
    pos[gsm] = t[t.CodeClass == "Positive"].set_index("Name").Count
C = pd.DataFrame(cnt); H = pd.DataFrame(hk); P = pd.DataFrame(pos)
gm = lambda D: np.exp(np.log(D.clip(lower=1)).mean(axis=0))
pf = gm(P); pf = pf.mean() / pf; C = C * pf; H = H * pf
hf = gm(H); hf = hf.mean() / hf; C = C * hf
S = S[S.gsm.isin(C.columns)]
rows.append(save("GSE215868", np.log2(C + 1), pd.Series(S.y.values, index=S.gsm.values), "panel", "NanoString IO 360, positive + housekeeping normalized, log2"))

# ---- panel: GSE93157 (NSCLC/HNSCC/melanoma, NanoString PanCancer Immune), raw counts; SAMPLE k = Patientk
S = soft("GSE93157"); S["y"] = S["best.resp"].map(recist); S = S.dropna(subset=["y"])
raw = gzip.open(f"{RAW}/GSE93157_raw_data_values.txt.gz", "rt").read().splitlines()
raw = [l for l in raw if not l.startswith(("#", '"#')) and l.strip("\t ")]
C = pd.read_csv(io.StringIO("\n".join(raw)), sep="\t", index_col=0).dropna(axis=1, how="all")
C.columns = [re.sub(r"SAMPLE\s*", "Patient", c) for c in C.columns]
C = C.apply(pd.to_numeric, errors="coerce").dropna(how="all")
S = S[S.title.isin(C.columns)]
rows.append(save("GSE93157", median_of_ratios(C), pd.Series(S.y.values, index=S.title.values), "panel", "NanoString PanCancer Immune, median-of-ratios, log2"))

# ---- chemo-immunotherapy (exploratory): GSE207422 pre-treatment, major pathological response
m = pd.read_excel(f"{RAW}/GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx")
m = m[m.Resource.astype(str).str.contains("Pre", case=False)]
m["y"] = m["Pathologic Response"].map(lambda v: 1 if str(v).startswith("MPR") else (0 if str(v) == "NMPR" else np.nan)); m = m.dropna(subset=["y"])
M = pd.read_csv(f"{RAW}/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz", sep="\t", index_col=0)
m = m[m.Sample.isin(M.columns)]
rows.append(save("GSE207422", M, pd.Series(m.y.values, index=m.Sample.values), "chemoICB", "log2 TPM; MPR against none"))

# ---- chemo-immunotherapy (exploratory): GSE235919, CR/PR against SD
S = soft("GSE235919"); S["y"] = S["recist criterion"].map(recist); S = S.dropna(subset=["y"])
M = pd.read_csv(f"{RAW}/GSE235919_TPM.txt.gz", sep="\t", index_col=0)
S["col"] = S.desc.str.replace("-MEDI", "", regex=False).str.replace(r"^(\d\dBM)(\d+)$", r"\1-\2", regex=True)
S = S[S.col.isin(M.columns)]
rows.append(save("GSE235919", np.log2(M + 1), pd.Series(S.y.values, index=S.col.values), "chemoICB", "log2(TPM+1); CR/PR against SD"))

pd.DataFrame(rows).to_csv(f"{OUT}/cohort_manifest.csv", index=False)
