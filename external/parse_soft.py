"""Turn a series' downloaded samples (GEO query page, SOFT text with data tables) and its platform table into
<out>/<cohort>_expr.tsv.gz (gene symbol x sample; the probe with the largest mean kept per symbol; log2 when the
values are unlogged) and <out>/<cohort>_clin.tsv (sample, time, event), with the fields and filters of cohorts.json.
Read as text only; run with python -I.

usage: python -I parse_soft.py <cohorts.json> <cohort> <raw dir> <out dir>
"""
import sys, os, json, gzip, io, re
import numpy as np, pandas as pd

CFG, CID, RAW, OUT = sys.argv[1:5]
cfg = json.load(open(CFG))[CID]
gse = cfg["gse"]; d = os.path.join(RAW, gse)


def platform_symbols(gpl):
    with gzip.open(os.path.join(d, f"{gpl}_table.soft.gz"), "rt", encoding="utf-8", errors="ignore") as f:
        txt = f.read()
    a = txt.index("!platform_table_begin"); b = txt.index("!platform_table_end")
    t = pd.read_csv(io.StringIO(txt[a:b].split("\n", 1)[1]), sep="\t", dtype=str, low_memory=False)
    col = cfg.get("symbol_col")
    if col is None:
        for c in t.columns:
            if c.strip().lower() in ("gene symbol", "symbol", "gene_symbol", "ilmn_gene", "genesymbol", "gene_name"):
                col = c; break
    if col is None and "gene_assignment" in t.columns:
        col = "gene_assignment"
    s = t.set_index(t.columns[0])[col].astype(str)
    if col == "gene_assignment":
        s = s.str.split(" // ").str[1]
    s = s.str.split(r" /// |///|;|,").str[0].str.strip()
    s = s[s.notna() & ~s.isin(["", "nan", "---", "NA"])]
    return s, col


samples, meta = {}, {}
# sample annotations: the brief view saved at screening (external/meta/<GSE>_gsm.txt); data tables: the data view
MD = os.path.join(os.path.dirname(os.path.abspath(RAW)), "meta")
with open(os.path.join(MD, f"{gse}_gsm.txt"), encoding="utf-8", errors="ignore") as f:
    cur, chars = None, {}
    for line in f:
        if line.startswith("^SAMPLE"):
            if cur:
                meta[cur] = chars
            cur = line.split("=", 1)[1].strip(); chars = {}
        elif cur and line.startswith("!Sample_characteristics_ch1"):
            v = line.split("=", 1)[1].strip()
            if ":" in v:
                k, x = v.split(":", 1); k = k.strip()
                kk, i = k, 2
                while kk in chars:
                    kk = f"{k}#{i}"; i += 1
                chars[kk] = x.strip()
        elif cur and line.startswith("!Sample_source_name_ch1"):
            chars["(source)"] = line.split("=", 1)[1].strip()
        elif cur and line.startswith("!Sample_platform_id"):
            chars["(platform)"] = line.split("=", 1)[1].strip()
    if cur:
        meta[cur] = chars
with gzip.open(os.path.join(d, f"{gse}_gsm_data.soft.gz"), "rt", encoding="utf-8", errors="ignore") as f:
    cur, chars, plat, in_tab, ids, vals = None, {}, None, False, [], []
    for line in f:
        if line.startswith("^SAMPLE"):
            cur = line.split("=", 1)[1].strip(); chars = {}; plat = None
        elif cur and line.startswith("!Sample_characteristics_ch1"):
            v = line.split("=", 1)[1].strip()
            if ":" in v:
                k, x = v.split(":", 1); k = k.strip()
                kk, i = k, 2
                while kk in chars:
                    kk = f"{k}#{i}"; i += 1
                chars[kk] = x.strip()
        elif cur and line.startswith("!Sample_source_name_ch1"):
            chars["(source)"] = line.split("=", 1)[1].strip()
        elif cur and line.startswith("!Sample_platform_id"):
            plat = line.split("=", 1)[1].strip()
        elif line.startswith("!sample_table_begin"):
            in_tab = True; header = None; ids, vals = [], []
        elif line.startswith("!sample_table_end"):
            in_tab = False
            samples[cur] = pd.Series(pd.to_numeric(pd.Series(vals), errors="coerce").values, index=ids)
        elif in_tab:
            p = line.rstrip("\n").split("\t")
            if header is None:
                header = p; vi = header.index("VALUE") if "VALUE" in header else 1; continue
            ids.append(p[0]); vals.append(p[vi] if len(p) > vi else "")
meta = pd.DataFrame(meta).T
keep = pd.Series(True, index=meta.index)
for col, allowed in cfg.get("filters", {}).items():
    v = meta[col].astype(str).str.strip().str.lower()
    keep &= v.isin([a.lower() for a in allowed]) if isinstance(allowed, list) else v.str.startswith(allowed.lower())
t = pd.to_numeric(meta[cfg["time"]].astype(str).str.extract(r"([-\d.]+)")[0], errors="coerce") * cfg.get("time_scale", 1.0)
ev = meta[cfg["event"]].astype(str).str.strip().str.lower().map({str(k).lower(): v for k, v in cfg["event_map"].items()})
keep &= t.notna() & ev.notna() & (t > 0)
sel = list(meta.index[keep])
gpl = cfg["gpl"]
sym, symcol = platform_symbols(gpl)
X = pd.DataFrame({s: samples[s] for s in sel})
X.index = X.index.astype(str)
X = X.loc[X.index.intersection(sym.index)]
mx = float(np.nanmax(X.values))
if mx > 50:
    X = np.log2(X.clip(lower=1))
X = X[X.notna().all(axis=1)]
X["_s"] = sym.reindex(X.index).values
X["_m"] = X.drop(columns="_s").mean(axis=1)
X = X.sort_values("_m", ascending=False).drop_duplicates("_s").set_index("_s").drop(columns="_m")
X.index.name = "symbol"
os.makedirs(OUT, exist_ok=True)
X.to_csv(os.path.join(OUT, f"{CID}_expr.tsv.gz"), sep="\t", compression="gzip")
clin = pd.DataFrame({"sample": sel, "time": t[sel].values, "event": ev[sel].astype(int).values})
clin.to_csv(os.path.join(OUT, f"{CID}_clin.tsv"), sep="\t", index=False)
print(json.dumps({"cohort": CID, "gse": gse, "gpl": gpl, "symbol_col": symcol, "samples_in_series": len(meta),
                  "patients": len(sel), "events": int(clin.event.sum()), "genes": int(X.shape[0]), "max_raw": mx,
                  "endpoint": cfg["endpoint"]}))
