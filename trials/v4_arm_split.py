"""Is the interaction specific to checkpoint blockade, or to the regimen? (reviewer fix)

IMmotion150 randomised three arms: atezolizumab alone, atezolizumab + bevacizumab, and
sunitinib. The published analysis pools the two atezolizumab arms, so the interaction it
reports could in principle reflect the regimen (anti-VEGF included) rather than
checkpoint blockade itself. cBioPortal records the arms separately
(ICI_RX = Atezolizumab with NON_ICI_RX = Bevacizumab or missing), so the trial can be
split:

    atezolizumab monotherapy vs sunitinib     - checkpoint blockade, no anti-VEGF
    atezolizumab + bevacizumab vs sunitinib   - checkpoint blockade plus anti-VEGF
    pooled (as published)                     - both atezolizumab arms together

Only the first isolates checkpoint blockade, and it is the smallest, so its event count
is printed beside its P values: an absence of significance there is an absence of power,
not evidence against the claim.

PRINCE is also reported: what its two arms actually contain is printed from the data.

Output: pancan/v4_arm_split.csv (one row per signature x contrast, plus summary rows).
"""
import warnings, functools; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import subprocess
import numpy as np, pandas as pd
from lifelines import CoxPHFitter

P = "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan"

IMM = {'T_cells', 'CD8_Tcells', 'Cytotoxic', 'NK_cells', 'B_cells', 'Plasma_cells',
       'Exhaustion', 'Checkpoints', 'GEP18_TIS (Ayers 2017)', 'IFNG_6gene (Ayers 2017)',
       'TLS (Cabrita 2020)', 'ESTIMATE_Immune', 'barriernet_Immune', 'APM', 'Treg', 'DC',
       'Chemokine12 (Messina 2012)', 'TLS_9gene (Prabhakaran)', 'CYT (Rooney 2015)'}


def signatures():
    raw = subprocess.run(
        ["Rscript", "-e", f'source("{P}/signatures.R");SIGS<-SIGS[!sapply(SIGS,is.null)];'
         'for(n in names(SIGS)) cat(n,"\\t",paste(SIGS[[n]],collapse=","),"\\n",sep="")'],
        capture_output=True, text=True).stdout
    return {k.strip(): [x for x in v.split(",") if x]
            for k, v in (l.split("\t", 1) for l in raw.strip().split("\n") if "\t" in l)}


def pfs_fields(C):
    t = pd.to_numeric(C["PFS_MONTHS"], errors="coerce")
    s = C["PFS_STATUS"].astype(str)
    e = pd.to_numeric(s.str.extract(r"^\s*([01])\s*:")[0], errors="coerce")
    if e.isna().all():
        e = s.str.contains("1|PROGRESS", case=False, regex=True).astype(float)
    return t, e.fillna(0).astype(int)


def interaction(score, arm, t, e):
    df = pd.DataFrame({"T": t.values, "E": e.values,
                       "s": (score - np.mean(score)) / (np.std(score) + 1e-9),
                       "arm": arm.values})
    df["sx"] = df.s * df.arm
    try:
        f = CoxPHFitter().fit(df, "T", "E")
        return float(f.params_["sx"]), float(f.summary.loc["sx", "p"])
    except Exception:
        return np.nan, np.nan


SIGS = signatures()
rows = []

# ---------------------------------------------------------------- IMmotion150
C = pd.read_csv(f"{P}/v4_immotion150_clinical_api.csv", index_col=0)
M = pd.read_csv(f"{P}/ici/rcc_iatlas_immotion150_2018.csv", index_col=0)
keep = [s for s in M.columns if s in C.index]
C, M = C.loc[keep], M[keep]
Z = ((M.T - M.mean(1)) / M.std(1).replace(0, 1)).T
t, e = pfs_fields(C)

atezo = C["ICI_RX"].astype(str) == "Atezolizumab"
bev = C["NON_ICI_RX"].astype(str).str.contains("Bevacizumab", na=False)
suni = C["NON_ICI_RX"].astype(str).str.contains("Sunitinib", na=False)
groups = {"atezo_mono_vs_sunitinib": (atezo & ~bev, suni),
          "atezo_bev_vs_sunitinib": (atezo & bev, suni),
          "pooled_atezo_vs_sunitinib": (atezo, suni)}
print("IMmotion150 arms: atezolizumab mono n=%d, atezolizumab+bev n=%d, sunitinib n=%d"
      % (int((atezo & ~bev).sum()), int((atezo & bev).sum()), int(suni.sum())))

for label, (tr, ct) in groups.items():
    sel = (tr | ct) & t.notna() & (t > 0)
    arm = tr[sel].astype(int)
    tt, ee = t[sel], e[sel]
    n1, n0, ev = int(arm.sum()), int((1 - arm).sum()), int(ee.sum())
    print(f"\n--- {label}: treated n={n1}, comparator n={n0}, events={ev}")
    stats = []
    for name, gl in SIGS.items():
        g = [x for x in gl if x in Z.index]
        if len(g) < 3:
            continue
        b, p = interaction(Z.loc[g, sel[sel].index].mean(0).values, arm, tt, ee)
        if not np.isfinite(p):
            continue
        fam = ("Immune" if name in IMM else
               "Control" if name.startswith("CTRL_") else "Other")
        rows.append(dict(row_type="signature", trial="IMmotion150", contrast=label,
                         signature=name, family=fam, n_genes=len(g), n_treated=n1,
                         n_comparator=n0, events=ev, coef=b, p=p))
        stats.append((fam, b, p))
    S = pd.DataFrame(stats, columns=["fam", "b", "p"])
    for fam in ("Immune", "Control"):
        g = S[S.fam == fam]
        if not len(g):
            continue
        rows.append(dict(row_type="summary", trial="IMmotion150", contrast=label,
                         family=fam, n_signatures=len(g), n_treated=n1, n_comparator=n0,
                         events=ev, n_negative=int((g.b < 0).sum()),
                         n_p05=int((g.p < .05).sum()),
                         n_p05_checkpoint=int(((g.p < .05) & (g.b < 0)).sum()),
                         median_coef=float(g.b.median())))
        print(f"    {fam:8s} n={len(g):2d}: {int((g.b<0).sum())} negative, "
              f"{int(((g.p<.05)&(g.b<0)).sum())} at P<0.05 favouring checkpoint, "
              f"median coef {g.b.median():+.3f}")

# ---------------------------------------------------------------- PRINCE
C = pd.read_csv(f"{P}/v2_prince_clinical_api.csv", index_col=0)
M = pd.read_csv(f"{P}/ici/paad_iatlas_prince_2022.csv", index_col=0)
keep = [s for s in M.columns if s in C.index]
C, M = C.loc[keep], M[keep]
Z = ((M.T - M.mean(1)) / M.std(1).replace(0, 1)).T
t, e = pfs_fields(C)
nivo = C["ICI_RX"].astype(str).str.contains("Nivolumab", case=False, na=False)
print("\nPRINCE arm contents as recorded in the data:")
print("  ICI_RX      :", C["ICI_RX"].value_counts(dropna=False).to_dict())
print("  ICI_TARGET  :", C["ICI_TARGET"].value_counts(dropna=False).to_dict())
print("  the data record no chemotherapy backbone column; PRINCE gave gemcitabine/")
print("  nab-paclitaxel in every arm, the comparator adding the CD40 agonist sotigalimab")

sel = t.notna() & (t > 0)
arm = nivo[sel].astype(int)
tt, ee = t[sel], e[sel]
n1, n0, ev = int(arm.sum()), int((1 - arm).sum()), int(ee.sum())
print(f"--- nivolumab_vs_non_checkpoint: treated n={n1}, comparator n={n0}, events={ev}")
stats = []
for name, gl in SIGS.items():
    g = [x for x in gl if x in Z.index]
    if len(g) < 3:
        continue
    b, p = interaction(Z.loc[g, sel[sel].index].mean(0).values, arm, tt, ee)
    if not np.isfinite(p):
        continue
    fam = ("Immune" if name in IMM else
           "Control" if name.startswith("CTRL_") else "Other")
    rows.append(dict(row_type="signature", trial="PRINCE",
                     contrast="nivolumab_vs_non_checkpoint", signature=name, family=fam,
                     n_genes=len(g), n_treated=n1, n_comparator=n0, events=ev,
                     coef=b, p=p))
    stats.append((fam, b, p))
S = pd.DataFrame(stats, columns=["fam", "b", "p"])
for fam in ("Immune", "Control"):
    g = S[S.fam == fam]
    if not len(g):
        continue
    rows.append(dict(row_type="summary", trial="PRINCE",
                     contrast="nivolumab_vs_non_checkpoint", family=fam,
                     n_signatures=len(g), n_treated=n1, n_comparator=n0, events=ev,
                     n_negative=int((g.b < 0).sum()), n_p05=int((g.p < .05).sum()),
                     n_p05_checkpoint=int(((g.p < .05) & (g.b < 0)).sum()),
                     median_coef=float(g.b.median())))
    print(f"    {fam:8s} n={len(g):2d}: {int((g.b<0).sum())} negative, "
          f"{int(((g.p<.05)&(g.b<0)).sum())} at P<0.05 favouring checkpoint, "
          f"median coef {g.b.median():+.3f}")

R = pd.DataFrame(rows)
R.to_csv(f"{P}/v4_arm_split.csv", index=False)
print(f"\nwrote v4_arm_split.csv ({len(R)} rows)")
