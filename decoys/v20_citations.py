"""Round 15 A4, exploratory: citations per year (PubMed Central cited-by counts, retrieved once on one day) for models with rho >= 0.95
against rho < 0.95, among models with decoys. -> icb/v20_round15/citations.csv, citations_summary.csv"""
import json, time, urllib.parse, urllib.request
import numpy as np, pandas as pd
from scipy import stats
G = "/path/to/analysis_root"; O = f"{G}/icb/v20_round15"
U = pd.read_csv(f"{O}/audit_models_units.csv"); U = U[U.n_coex.notna()].copy()
ids = sorted({str(p).replace("PMC", "") for p in U.pmcid.dropna()})
print(len(U), "models with decoys,", len(ids), "papers")
base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
counts = {}
for i in range(0, len(ids), 100):
    chunk = ids[i:i + 100]
    q = urllib.parse.urlencode({"dbfrom": "pmc", "linkname": "pmc_pmc_citedby", "retmode": "json", "tool": "signature-audit"} ) + "".join(f"&id={c}" for c in chunk)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(base, data=q.encode()), timeout=120) as r:
                J = json.load(r)
            break
        except Exception as e:  # noqa: BLE001
            print("retry", e); time.sleep(3 * (attempt + 1))
    else:
        raise SystemExit("elink failed")
    for ls in J["linksets"]:
        pid = ls["ids"][0]
        n = 0
        for db in ls.get("linksetdbs", []):
            if db["linkname"] == "pmc_pmc_citedby":
                n = len(db["links"])
        counts[str(pid)] = n
    time.sleep(0.4)
U["pid"] = U.pmcid.astype(str).str.replace("PMC", "")
U["citedby"] = U.pid.map(counts)
RET = 2026 + 10.0 / 12    # retrieved 2026-11
U["cites_per_year"] = U.citedby / (RET - (U.year + 0.5))
U[["model", "type", "pmcid", "year", "rho", "uses_external", "validated", "citedby", "cites_per_year"]].to_csv(f"{O}/citations.csv", index=False)
rows = []
for lab, S in (("all models with decoys", U), ("independent", U[~U.uses_external.astype(bool)])):
    hi, lo = S[S.rho >= 0.95].cites_per_year.dropna(), S[S.rho < 0.95].cites_per_year.dropna()
    p = stats.mannwhitneyu(hi, lo, alternative="two-sided").pvalue
    rows.append(dict(population=lab, n_hi=len(hi), n_lo=len(lo), median_hi=hi.median(), median_lo=lo.median(), mean_hi=hi.mean(), mean_lo=lo.mean(), p=p))
R = pd.DataFrame(rows); R.to_csv(f"{O}/citations_summary.csv", index=False); print(R.round(3).to_string(index=False))
print("missing counts:", U.citedby.isna().sum())
