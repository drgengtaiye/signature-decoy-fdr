"""I-SPY2 randomized contrasts (registration, secondary 5): arm and pCR for
  ISPY2_pembro  GSE194040: paclitaxel + pembrolizumab against the HER2-negative paclitaxel controls
  ISPY2_durva   GSE173839: durvalumab + olaparib with paclitaxel against its controls (evaluable pCR)
Expression is the deposited gene-level log2 matrix, restricted to the patients of the contrast."""
import numpy as np, pandas as pd
V = "/groups/adv2105_gp/yichen/Yi/com/gastric/validation"; RAW = f"{V}/raw"; OUT = f"{V}/cohorts"


def soft(acc):
    out = []
    for r in open(f"{V}/{acc}_gsm.soft", encoding="utf-8", errors="ignore").read().split("^SAMPLE = ")[1:]:
        d = {"gsm": r.split("\n", 1)[0].strip()}
        for l in r.splitlines():
            if l.startswith("!Sample_characteristics_ch1 = ") and ": " in l:
                k, v = l.split(" = ", 1)[1].split(": ", 1); d[k] = v
        out.append(d)
    return pd.DataFrame(out)


S = soft("GSE194040")
S = S[S.arm.isin(["Paclitaxel", "Paclitaxel + Pembrolizumab"]) & (S.her2 == "0")].copy()
M = pd.read_csv(f"{RAW}/GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt.gz", sep="\t", index_col=0)
M.columns = [str(c) for c in M.columns]
S = S[S["patient id"].isin(M.columns)].drop_duplicates("patient id")
C = pd.DataFrame({"sample": S["patient id"].values, "arm": (S.arm == "Paclitaxel + Pembrolizumab").astype(int).values,
                  "y": S.pcr.astype(int).values, "hr": S.hr.astype(int).values})
M[C["sample"]].to_csv(f"{OUT}/ISPY2_pembro.expr.tsv.gz", sep="\t", compression="gzip")
C.to_csv(f"{OUT}/ISPY2_pembro.clin.tsv", sep="\t", index=False)
print("ISPY2_pembro", M.shape[0], "genes;", pd.crosstab(C.arm, C.y).to_dict())

B = pd.read_csv(f"{RAW}/GSE173839_ISPY2_DurvalumabOlaparibArm_biomarkers.csv.gz")
B = B[B["pCR.status"].isin([0, 1])]
M = pd.read_csv(f"{RAW}/GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_meanCol_geneLevel_n105.txt.gz", sep="\t", index_col=0)
M.columns = [str(c) for c in M.columns]; B["ResearchID"] = B.ResearchID.astype(str)
B = B[B.ResearchID.isin(M.columns)]
C = pd.DataFrame({"sample": B.ResearchID.values, "arm": (B.Arm == "durvalumab/olaparib").astype(int).values,
                  "y": B["pCR.status"].astype(int).values, "hr": B["HR.status..0.HR....1..HR.."].astype(int).values})
M[C["sample"]].to_csv(f"{OUT}/ISPY2_durva.expr.tsv.gz", sep="\t", compression="gzip")
C.to_csv(f"{OUT}/ISPY2_durva.clin.tsv", sep="\t", index=False)
print("ISPY2_durva", M.shape[0], "genes;", pd.crosstab(C.arm, C.y).to_dict())
