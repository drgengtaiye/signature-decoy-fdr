"""Round 15 (validation/AMENDMENT_3_round15...): parse_soft.py of round 11 with the round-15 layout (sample annotations
in <raw>/<GSE>/<GSE>_samples_brief.txt) and three config keys: "composite" (all annotations in one characteristics
line, split on ';' into key: value pairs), "time_start" (time = the date in "time" minus the date in "time_start", in
days, then time_scale), and samples restricted to the config's platform(s).

Turn a series' downloaded samples (GEO query page, SOFT text with data tables) and its platform table into
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
            if c.strip().lower() in ("gene symbol", "symbol", "gene_symbol", "ilmn_gene", "genesymbol", "gene_name", "gene.symbol", "orf"):
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
with open(os.path.join(d, f"{gse}_samples_brief.txt"), encoding="utf-8", errors="ignore") as f:
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
if cfg.get("composite"):                         # e.g. GSE14333: "Location: Left; DukesStage: B; DFS_Time: 22.38; ..."
    ck = cfg["composite"]
    for i, v in meta[ck].items():
        for j, part in enumerate(str(v).split(";")):
            if j == 0:
                meta.at[i, ck] = part.strip()
            elif ":" in part:
                a, b = part.split(":", 1); meta.at[i, a.strip()] = b.strip()
gpls = cfg["gpl"] if isinstance(cfg["gpl"], list) else [cfg["gpl"]]
keep = meta["(platform)"].isin(gpls) if "(platform)" in meta else pd.Series(True, index=meta.index)
for col, allowed in cfg.get("filters", {}).items():
    v = meta[col].astype(str).str.strip().str.lower()
    keep &= v.isin([a.lower() for a in allowed]) if isinstance(allowed, list) else v.str.startswith(allowed.lower())
if cfg.get("time_start"):                        # e.g. GSE41271: last follow-up date minus date of surgery
    t = (pd.to_datetime(meta[cfg["time"]], errors="coerce") - pd.to_datetime(meta[cfg["time_start"]], errors="coerce")).dt.days
    t = t.astype(float) * cfg.get("time_scale", 1.0)
else:
    t = pd.to_numeric(meta[cfg["time"]].astype(str).str.extract(r"([-\d.]+)")[0], errors="coerce") * cfg.get("time_scale", 1.0)
ev = meta[cfg["event"]].astype(str).str.strip().str.lower().map({str(k).lower(): v for k, v in cfg["event_map"].items()})
keep &= t.notna() & ev.notna() & (t > 0)
sel = list(meta.index[keep])
gpl = gpls[0]
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
