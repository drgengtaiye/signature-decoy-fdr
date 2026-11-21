#!/usr/bin/env python3
"""Prepare public breast-cancer expression cohorts with a survival endpoint.

DATA PREPARATION ONLY (plan PLAN_round6_applications.md, section 3).  Nothing in this
script relates expression to outcome: no survival model, concordance, Kaplan-Meier
estimate or gene-outcome correlation is computed anywhere.  The only expression-based
computations are sample-sample correlations used to find duplicated arrays/patients.

Order of the pipeline (all scripts in this directory; none is run from a download directory):
  1. fetch_breast_data.sh            GEO series matrices + GPL annot + suppl. clinical text, GSE96058, METABRIC
  2. Rscript export_orgHs_symbols.R  symbol/alias -> Entrez tables from the installed org.Hs.eg.db
  3. prep_breast_cohorts.py          this script
  4. python3 -I verify_breast_cohorts.py   format checks of the outputs

Run (the .Rda files of the Venet bundle are parsed with the pure-Python `rdata`
package, never with R):

    RDENV=<a Python virtual environment with rdata, numpy and pandas>/bin/python
    $RDENV -I prep_breast_cohorts.py [--base /path/to/analysis_root/breast]

Inputs (defaults under --base): download_venet/, download_geo/<GSE>/, download_metabric/,
annotation/orgHs_{symbol,alias}2entrez.tsv (written by export_orgHs_symbols.R from the
installed org.Hs.eg.db), and the local TCGA files given by --tcga-expr / --tcga-surv.
Outputs: <base>/cohorts/<ID>_expr.tsv.gz, <ID>_clin.tsv, manifest.csv, plus
overlap_diagnostic.csv, duplicates_within.csv and prep_log.txt.
"""
import argparse
import collections
import gzip
import io
import re
import sys
import warnings

import numpy as np
import pandas as pd

MIN_N, MIN_EVENTS, MIN_GENES = 100, 25, 10000
DAYS_PER_MONTH = 365.25 / 12.0
DUP_R = 0.95            # within-cohort sample-sample r above which two arrays are taken as one tumour
CROSS_R_REPORT = 0.70   # cross-cohort pairs reported in overlap_diagnostic.csv
ENDPOINT_ORDER = ["DMFS", "RFS", "DFS", "DFI", "DSS", "OS"]   # DMFS > RFS/DFS(/DFI) > DSS > OS

LOG = []


def log(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s, flush=True)


# --------------------------------------------------------------------------- generic readers
def read_series_matrix(path, want_table=True):
    """Parse a GEO series_matrix.txt.gz (text). Returns (meta DataFrame indexed by GSM, expr)."""
    rows = collections.OrderedDict()
    table = []
    started = False
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                started = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if started:
                if want_table:
                    table.append(line)
                continue
            if line.startswith("!Sample_"):
                parts = line.rstrip("\n").split("\t")
                key, k, i = parts[0], parts[0], 1
                while k in rows:
                    i += 1
                    k = f"{key}#{i}"
                rows[k] = [v.strip('"') for v in parts[1:]]
    raw = pd.DataFrame(rows)
    recs = []
    for _, r in raw.iterrows():
        d = {"gsm": r["!Sample_geo_accession"], "title": r["!Sample_title"],
             "source": r.get("!Sample_source_name_ch1", ""), "description": r.get("!Sample_description", "")}
        for c in raw.columns:
            if c.startswith("!Sample_characteristics"):
                v = r[c]
                if ": " in v:
                    kk, val = v.split(": ", 1)
                    d.setdefault(kk.strip(), val.strip())
        recs.append(d)
    meta = pd.DataFrame(recs).set_index("gsm", drop=False)
    expr = None
    if want_table and len(table) > 1:
        expr = pd.read_csv(io.StringIO("".join(table)), sep="\t", index_col=0,
                           na_values=["null", "NA", "", "NaN", "nan"], low_memory=False)
        expr.index = expr.index.astype(str).str.strip('"')
        expr.columns = [c.strip('"') for c in expr.columns]
        expr = expr.astype(np.float64)
    return meta, expr


def read_gpl_annot(path):
    """GEO GPLxxx.annot.gz -> Series probe -> Entrez (str, may contain '///' or be empty)."""
    with gzip.open(path, "rt", errors="replace") as f:
        lines = f.readlines()
    b = next(i for i, l in enumerate(lines) if l.startswith("!platform_table_begin"))
    e = next((i for i, l in enumerate(lines) if l.startswith("!platform_table_end")), len(lines))
    tab = pd.read_csv(io.StringIO("".join(lines[b + 1:e])), sep="\t", dtype=str, usecols=["ID", "Gene ID"])
    return tab.set_index("ID")["Gene ID"].fillna("")


def maybe_log2(expr):
    """log2-transform when the deposited values are clearly not logged (max > 50)."""
    mx = np.nanmax(expr.values)
    if mx > 50:
        pos = expr.values[expr.values > 0]
        floor = float(pos.min())
        n_floor = int((expr.values <= 0).sum())
        out = np.log2(expr.clip(lower=floor))
        return out, f"deposited values unlogged (max {mx:.4g}); log2-transformed" + (
            f" after flooring {n_floor} non-positive values at {floor:.3g}" if n_floor else "")
    return expr, f"deposited values already log scale (max {mx:.3g}); used as deposited"


def collapse_max_mean(expr, gene_ids):
    """expr rows -> genes: keep, for each Entrez ID, the row with the largest mean across samples."""
    df = expr.copy()
    df["__g"] = gene_ids.values
    df = df[df["__g"].notna()]
    df["__m"] = df.drop(columns="__g").mean(axis=1, skipna=True)
    df = df.sort_values(["__g", "__m"], ascending=[True, False], kind="mergesort")
    df = df.drop_duplicates("__g", keep="first")
    out = df.drop(columns=["__m"]).set_index("__g")
    out.index = out.index.astype(np.int64)
    out.index.name = "entrez"
    return out.sort_index()


def probes_to_genes(expr, probe2entrez):
    ann = probe2entrez.reindex(expr.index).fillna("").astype(str).str.strip()
    single = ann.str.fullmatch(r"\d+")
    n_multi = int(ann.str.contains("///").sum())
    n_none = int((ann == "").sum())
    keep = expr.loc[single.values]
    genes = collapse_max_mean(keep, ann[single].astype(np.int64))
    note = (f"{len(expr)} probes: {n_none} without Entrez, {n_multi} multi-Entrez dropped; "
            f"{int(single.sum())} single-Entrez probes -> {len(genes)} genes (max-mean probe)")
    return genes, note


class SymbolMapper:
    """Gene symbol -> Entrez via org.Hs.eg.db: official symbol first, else an alias that maps to one gene."""

    def __init__(self, annot_dir):
        sym = pd.read_csv(f"{annot_dir}/orgHs_symbol2entrez.tsv", sep="\t", dtype=str)
        ali = pd.read_csv(f"{annot_dir}/orgHs_alias2entrez.tsv", sep="\t", dtype=str)
        s = sym.groupby("symbol")["gene_id"].agg(lambda x: sorted(set(x)))
        a = ali.groupby("alias_symbol")["gene_id"].agg(lambda x: sorted(set(x)))
        self.sym, self.ali = s.to_dict(), a.to_dict()

    def map(self, symbols):
        out, how = [], collections.Counter()
        for x in symbols:
            x = str(x)
            ids = self.sym.get(x)
            if ids is None and x.upper() != x:
                ids = self.sym.get(x.upper())
            if ids is not None:
                if len(ids) == 1:
                    out.append(int(ids[0])); how["official"] += 1
                else:
                    out.append(np.nan); how["official_multi_dropped"] += 1
                continue
            ids = self.ali.get(x)
            if ids is not None and len(ids) == 1:
                out.append(int(ids[0])); how["alias_unique"] += 1
            elif ids is not None:
                out.append(np.nan); how["alias_multi_dropped"] += 1
            else:
                out.append(np.nan); how["unmapped"] += 1
        return pd.Series(out, index=range(len(symbols)), dtype="float64"), dict(how)


def make_clin(samples, time, event, unit, endpoint):
    """Build the clinical table: positive, non-missing time (months) and event 0/1."""
    t = pd.to_numeric(pd.Series(time, index=samples), errors="coerce").astype(float)
    e = pd.to_numeric(pd.Series(event, index=samples), errors="coerce").astype(float)
    factor = {"days": 1.0 / DAYS_PER_MONTH, "months": 1.0, "years": 12.0}[unit]
    df = pd.DataFrame({"sample": list(samples), "time_months": (t * factor).values, "event": e.values})
    ok = df.time_months.notna() & (df.time_months > 0) & df.event.isin([0, 1])
    dropped = {"missing_or_nonpositive_time_or_event": int((~ok).sum())}
    df = df[ok].copy()
    df["event"] = df["event"].astype(int)
    df["endpoint"] = endpoint
    return df.reset_index(drop=True), dropped


def new_cohort(cid, source, platform):
    return dict(id=cid, source=source, platform=platform, expr=None, clin=None, keys=None,
                notes=[], reason=None, endpoint=None, n_probe_note="", log_note="")


# --------------------------------------------------------------------------- Venet bundle
def read_venet(path):
    import rdata  # pure-Python parser; never R
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parsed = rdata.parser.parse_file(path)
        study = rdata.conversion.convert(parsed)["study"]
    d = study["data"]
    cols = [str(c) for c in d.coords["dim_1"].values]
    expr = pd.DataFrame(np.asarray(d.values, dtype=np.float64), index=np.asarray(study["genes"]).astype(np.int64),
                        columns=cols)
    # 'samples' is a character matrix (annotation fields x samples); rebuild it from the parsed object
    samp = parsed.object.value[0].value[1]
    a, attrs = samp.attributes, {}
    while a is not None and a.info.type.name == "LIST":
        tag = a.tag.referenced_object if a.tag.info.type.name == "REF" else a.tag
        attrs[tag.value.value.decode()] = a.value[0]
        a = a.value[1] if len(a.value) > 1 else None
    dn = attrs["dimnames"].value
    rown = [x.value.decode() if x.value is not None else None for x in dn[0].value]
    coln = [x.value.decode() if x.value is not None else None for x in dn[1].value]
    vals = [x.value.decode() if x.value is not None else None for x in samp.value]
    ann = pd.DataFrame(np.array(vals, dtype=object).reshape(len(coln), len(rown)), index=coln, columns=rown)
    surv = pd.DataFrame(np.asarray(study["survival"].values, dtype=float), index=cols, columns=["time", "status"])
    gd = study["geneDescr"]
    assert np.array_equal(np.asarray(study["genes"]).astype(np.int64), gd["geneId"].values.astype(np.int64))
    assert list(ann.index) == cols
    return expr, ann, surv


def load_nki(a):
    c = new_cohort("NKI", "Venet et al. 2011 Dataset S1 (expression-data-nki-*.Rda; van de Vijver 2002)",
                   "Agilent/Rosetta Hu25K two-colour (log10 ratio)")
    f_rfs = f"{a.venet}/expression-data-nki-rfs.Rda"
    f_os = f"{a.venet}/expression-data-nki-os.Rda"
    e1, ann, s_rfs = read_venet(f_rfs)
    e2, ann2, s_os = read_venet(f_os)
    assert list(e1.columns) == list(e2.columns) and np.allclose(e1.values, e2.values)
    # the bundle's 'survival' slot: rfs file = t.rfs/365 (years) with e.rfs; os file = t.os/365 with e.os
    t_dm = pd.to_numeric(ann["t.dmfs"]); t_rfs = pd.to_numeric(ann["t.rfs"])
    e_dm = pd.to_numeric(ann["e.dmfs"]); e_rfs = pd.to_numeric(ann["e.rfs"])
    assert np.allclose(s_rfs["status"].values, e_rfs.values) and np.allclose(s_rfs["time"].values, t_rfs.values / 365, atol=1e-6)
    assert np.allclose(s_os["status"].values, pd.to_numeric(ann["e.os"]).values)
    c["notes"].append("bundle files: *-rfs.Rda survival slot = t.rfs/365 (y), e.rfs; *-os.Rda = t.os/365, e.os; "
                      "sample annotation also carries t.dmfs/e.dmfs (days), used here (DMFS first in hierarchy)")
    if np.allclose(t_dm, t_rfs) and np.array_equal(e_dm.values, e_rfs.values):
        c["notes"].append("in this bundle t.dmfs/e.dmfs are identical to t.rfs/e.rfs for all 295 patients")
    assert (ann["series"] == "NKI2").all() and (ann["primary.anatomical.site"] == "1").all()
    c["expr"] = e1 * np.log2(10.0)
    c["log_note"] = ("deposited values are log10 ratios (range -2..2); multiplied by log2(10) to give log2 ratios "
                     "(linear rescaling)")
    assert not e1.index.duplicated().any()
    c["n_probe_note"] = f"gene-level in bundle ({len(e1)} Entrez genes, one row per gene as deposited)"
    c["clin"], dr = make_clin(e1.columns, t_dm.values, e_dm.values, "days", "DMFS")
    c["endpoint"] = "DMFS"
    c["keys"] = pd.Series({s: "NKI:" + s for s in e1.columns})
    return c


def loi_key(name):
    ser, code = name.split("_", 1)
    if ser in ("KIT", "KIU"):
        return "UPP:" + code          # KIT/KIU codes are Uppsala tumour codes (same as GSE3494/GSE4922 'X' + code)
    return f"{ser}:{code}"


def load_loi(a):
    c = new_cohort("Loi_GSE6532", "Venet et al. 2011 Dataset S1 (expression-data-loi-*.Rda; Loi 2007, GSE6532)",
                   "Affymetrix HG-U133A+B (GPL96/GPL97) and HG-U133 Plus 2 (GPL570, GUYT); gene level, RMA per "
                   "series, probes median-centred")
    e_rfs, ann, s_rfs = read_venet(f"{a.venet}/expression-data-loi-rfs.Rda")
    e_os, ann_os, s_os = read_venet(f"{a.venet}/expression-data-loi-os.Rda")
    shared = [s for s in e_os.columns if s in e_rfs.columns]
    assert len(shared) == e_os.shape[1] and np.allclose(e_rfs[shared].values, e_os[shared].values)
    demo = pd.read_csv(f"{a.geo}/GSE6532/GSE6532_LUMINAL_demo.txt.gz", sep="\t", na_values=["NA"]).set_index("samplename")
    dm = demo.reindex(e_rfs.columns)
    assert dm.index.isin(demo.index).all()
    # verify which bundle file holds what, against GEO's GSE6532_LUMINAL_demo.txt (t.* in days)
    ok_rfs = np.allclose(s_rfs["time"].values, dm["t.rfs"].values / 365.25, atol=2e-3) and \
        np.array_equal(s_rfs["status"].values, dm["e.rfs"].fillna(0).values)
    dmo = demo.reindex(e_os.columns)
    ok_dm = np.allclose(s_os["time"].values, dmo["t.dmfs"].values / 365.25, atol=2e-3) and \
        np.array_equal(s_os["status"].values, dmo["e.dmfs"].values)
    log("Loi bundle check: rfs file == GEO t.rfs/e.rfs:", ok_rfs, "; os file == GEO t.dmfs/e.dmfs:", ok_dm)
    assert ok_dm
    c["notes"].append("bundle '*-os.Rda' survival (fields time/event.disease.related.death, 380 pts) is identical to "
                      "GEO GSE6532_LUMINAL_demo t.dmfs/e.dmfs, i.e. DMFS (not overall survival); '*-rfs.Rda' "
                      f"(393 pts) = t.rfs/e.rfs (match {ok_rfs}). DMFS used (first in hierarchy); bundle expression "
                      "(identical in both files) used")
    c["expr"] = e_rfs
    c["log_note"] = ("log2 RMA, normalised per series and median-centred per probe by the depositors (README of "
                     "GSE6532 LUMINAL); used as deposited")
    assert not e_rfs.index.duplicated().any()
    c["n_probe_note"] = f"gene-level in bundle ({len(e_rfs)} Entrez genes, one row per gene as deposited)"
    c["clin"], dr = make_clin(e_rfs.columns, dm["t.dmfs"].values, dm["e.dmfs"].values, "days", "DMFS")
    c["endpoint"] = "DMFS"
    c["keys"] = pd.Series({s: loi_key(s) for s in e_rfs.columns})
    ser = pd.Series([s.split("_")[0] for s in c["clin"]["sample"]]).value_counts().to_dict()
    c["notes"].append(f"sites in final set: {ser}")
    return c


# --------------------------------------------------------------------------- GEO series
def geo_expr(c, a, gse, matrix_name, gpl):
    meta, raw = read_series_matrix(f"{a.geo}/{gse}/{matrix_name}")
    raw, c["log_note"] = maybe_log2(raw)
    p2e = read_gpl_annot(f"{a.geo}/{gse}/{gpl}.annot.gz")
    genes, c["n_probe_note"] = probes_to_genes(raw, p2e)
    return meta, genes


def load_geo(a, gse):
    plat = {"GSE20685": "GPL570", "GSE21653": "GPL570", "GSE42568": "GPL570"}.get(gse, "GPL96")
    mat = {"GSE1456": "GSE1456-GPL96_series_matrix.txt.gz", "GSE3494": "GSE3494-GPL96_series_matrix.txt.gz",
           "GSE4922": "GSE4922-GPL96_series_matrix.txt.gz"}.get(gse, f"{gse}_series_matrix.txt.gz")
    pname = {"GPL96": "Affymetrix HG-U133A (GPL96)", "GPL570": "Affymetrix HG-U133 Plus 2.0 (GPL570)"}[plat]
    c = new_cohort(gse, f"GEO {gse} ({mat})", pname + (" [GPL96 part of GPL96+GPL97 series]" if "-GPL96" in mat else ""))
    meta, genes = geo_expr(c, a, gse, mat, plat)
    m = meta
    title = m["title"]
    unit, endpoint, t, e = None, None, None, None
    if gse == "GSE2034":
        c["reason"] = ("no survival time/event in GEO: characteristics hold only 'bone relapses'; the clinical sheet "
                       "referred to in the series summary is not among the GEO files (suppl = RAW.tar only)")
        keys = "ROT2034:" + title
    elif gse == "GSE7390":
        endpoint, unit, t, e = "DMFS", "days", m["t.dmfs"], m["e.dmfs"]
        site = title.str.extract(r"^VDX([A-Z]+)_(.+)$")
        assert site[0].notna().all()
        keys = pd.Series([("UPP:" + cd) if s == "KIU" else f"{s}:{cd}" for s, cd in zip(site[0], site[1])], index=m.index)
        c["notes"].append("TRANSBIG; titles VDX<site>_<code>; KIU codes are Uppsala tumour codes, OXFU = Oxford "
                          "untreated (as in GSE6532/GSE2990)")
    elif gse == "GSE11121":
        endpoint, unit, t, e = "DMFS", "months", m["t.dmfs"], m["e.dmfs"]
        keys = "MAINZ:" + title
    elif gse == "GSE2990":
        sup = pd.read_csv(f"{a.geo}/GSE2990/GSE2990_suppl_info.txt", sep="\t", lineterminator="\r",
                          na_values=["NA"], dtype={"geo_accn": str})
        sup = sup[sup["geo_accn"].astype(str).str.startswith("GSM")].set_index("geo_accn").reindex(m.index)
        endpoint, unit, t, e = "DMFS", "years", sup["time.dmfs"], sup["event.dmfs"]
        keys = pd.Series([loi_key(x) for x in title], index=m.index)
        c["notes"].append("outcome from GEO supplementary GSE2990_suppl_info.txt (time.dmfs in years)")
    elif gse == "GSE1456":
        endpoint, unit, t, e = "RFS", "years", m["SURV_RELAPSE"], m["RELAPSE"]
        keys = "STO:" + title.str.replace(r" \(HG-U133A\)$", "", regex=True)
    elif gse == "GSE3494":
        c["reason"] = ("no survival in GEO (characteristics: 'Breast Tumor tissue' only; suppl = RAW.tar only); "
                       "its Uppsala patients are the Uppsala cohort of GSE4922, which carries DFS")
        keys = "UPP:" + title.str.extract(r"^X(\S+)")[0]
    elif gse == "GSE4922":
        endpoint, unit = "DFS", "years"
        t = m["DFS TIME (yrs)"]
        e = m["DFS EVENT (0=censored; 1=event defined as any type of recurrence (local, regional or distant) or death from breast cancer"]
        idx = title.str.replace(r"^A ", "", regex=True)          # titles 'A X100B08' (Uppsala), 'A NUH.BT022' (Singapore)
        keys = pd.Series([("UPP:" + x[1:]) if re.fullmatch(r"X\d+[A-Z]\d+", x) else "SGP:" + x for x in idx], index=m.index)
        c["notes"].append("DFS event = any recurrence (local, regional, distant) or breast-cancer death; the 40 "
                          "Singapore samples have no outcome")
    elif gse == "GSE20685":
        endpoint, unit = "DMFS", "years"
        ev = pd.to_numeric(m["event_metastasis"], errors="coerce")
        tm = pd.to_numeric(m["time_to_metastasis (years)"], errors="coerce")
        fu = pd.to_numeric(m["follow_up_duration (years)"], errors="coerce")
        assert ((ev == 1) == tm.notna()).all()
        t = np.where(ev == 1, tm, fu)
        e = ev
        keys = "TWN:" + title
        c["notes"].append("DMFS time = time_to_metastasis (years) for patients with event_metastasis=1, else "
                          "follow_up_duration (years); time 0 (metastatic at diagnosis) dropped as non-positive")
    elif gse == "GSE21653":
        endpoint, unit, t, e = "DFS", "months", m["dfs time (months)"], m["dfs evt"]
        keys = "MRS:" + title
    elif gse == "GSE42568":
        tumour = m["tissue"].eq("breast cancer")
        c["notes"].append(f"{int((~tumour).sum())} normal-breast samples dropped")
        m = m[tumour]
        title = m["title"]
        endpoint, unit, t, e = "RFS", "days", m["relapse free survival time_days"], m["relapse free survival event"]
        keys = "DUB:" + title.str.replace("Breast cancer, ", "", regex=False)
    elif gse == "GSE17705":
        endpoint, unit, t, e = "DMFS", "years", m["event time (years)"], m["distant relapse (1=dr, 0 censored)"]
        lab = m["profiling lab"]
        pid = m["patient id"].astype(str)
        keys = pd.Series([f"{l}:{p}" for l, p in zip(lab, pid)], index=m.index)   # JBI ids resolved later vs Loi
        c["notes"].append("endpoint 'distant relapse' with 'event time (years)' (taken as DMFS); profiling labs "
                          f"{lab.value_counts().to_dict()}")
    elif gse == "GSE12093":
        endpoint, unit, t, e = "DFS", "months", m["DFS.time"], m["DFS.status"]
        keys = "ROT12093:" + title
    else:
        raise ValueError(gse)
    if not isinstance(keys, pd.Series):
        keys = pd.Series(keys)
    keys.index = m.index
    c["keys"] = keys
    genes = genes[[g for g in m.index if g in genes.columns]]
    c["expr"] = genes
    if endpoint is not None:
        c["clin"], dr = make_clin(m.index, np.asarray(t), np.asarray(e), unit, endpoint)
        c["endpoint"] = endpoint
    return c


def load_gse96058(a):
    c = new_cohort("GSE96058", "GEO GSE96058 (SCAN-B; gene_expression_3273_samples_and_136_replicates_transformed.csv)",
                   "RNA-seq, Illumina HiSeq 2000 (GPL11154) + NextSeq 500 (GPL18573); log2(FPKM+0.1)")
    d = f"{a.geo}/GSE96058"
    metas = [read_series_matrix(f"{d}/GSE96058-{p}_series_matrix.txt.gz", want_table=False)[0]
             for p in ("GPL11154", "GPL18573")]
    meta = pd.concat(metas)
    meta["instrument"] = meta["instrument model"]
    expr = pd.read_csv(f"{d}/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz",
                       index_col=0, keep_default_na=False, engine="c")
    expr = expr.astype(np.float64)
    rep = expr.columns.str.endswith("repl")
    c["notes"].append(f"{int(rep.sum())} technical replicates ('...repl' columns) dropped")
    expr = expr.loc[:, ~rep]
    t2g = dict(zip(meta["title"], meta["gsm"]))
    assert all(x in t2g for x in expr.columns)
    expr.columns = [t2g[x] for x in expr.columns]
    meta = meta.set_index("gsm", drop=False)
    meta = meta[~meta["title"].str.endswith("repl")]
    assert meta["source"].eq("Primary breast tumor").all()
    c["log_note"] = "log2(FPKM + 0.1) as deposited (max %.3g); used as deposited" % np.nanmax(expr.values)
    mapper = SymbolMapper(a.annot)
    ent, how = mapper.map(list(expr.index))
    genes = collapse_max_mean(expr.reset_index(drop=True), ent)
    c["n_probe_note"] = f"{len(expr)} gene symbols -> Entrez via org.Hs.eg.db {how}; -> {len(genes)} genes (max-mean)"
    c["expr"] = genes[[g for g in meta.index if g in genes.columns]]
    c["clin"], dr = make_clin(meta.index, meta["overall survival days"].values, meta["overall survival event"].values,
                              "days", "OS")
    c["endpoint"] = "OS"
    c["keys"] = "SCANB:" + meta["title"]
    c["notes"].append("only overall survival is deposited (no DMFS/RFS/DSS)")
    return c


def load_metabric(a):
    c = new_cohort("METABRIC", "cBioPortal datahub brca_metabric (data_mrna_illumina_microarray.txt, data_clinical_*.txt)",
                   "Illumina HT-12 v3 (log2 intensity)")
    d = a.metabric
    pat = pd.read_csv(f"{d}/data_clinical_patient.txt", sep="\t", comment="#", dtype={"PATIENT_ID": str})
    smp = pd.read_csv(f"{d}/data_clinical_sample.txt", sep="\t", comment="#", dtype=str)
    m = pd.read_csv(f"{d}/data_mrna_illumina_microarray.txt", sep="\t")
    ent = pd.to_numeric(m["Entrez_Gene_Id"], errors="coerce")
    expr = m.drop(columns=["Hugo_Symbol", "Entrez_Gene_Id"]).astype(np.float64)
    n_rows = len(expr)
    genes = collapse_max_mean(expr, ent)
    c["n_probe_note"] = f"{n_rows} rows (gene level, cBioPortal Entrez_Gene_Id; {int(ent.isna().sum())} without Entrez) -> {len(genes)} genes (max-mean)"
    expr, c["log_note"] = maybe_log2(genes)
    smp = smp.set_index("SAMPLE_ID")
    cols = [s for s in expr.columns if s in smp.index and smp.loc[s, "SAMPLE_TYPE"] == "Primary"]
    assert len(set(smp.loc[cols, "PATIENT_ID"])) == len(cols)
    expr = expr[cols]
    pat = pat.set_index("PATIENT_ID")
    p = pat.reindex(smp.loc[cols, "PATIENT_ID"].values)
    ev = p["RFS_STATUS"].astype(str).str.extract(r"^(\d):")[0]
    c["clin"], dr = make_clin(cols, p["RFS_MONTHS"].values, ev.values, "months", "RFS")
    c["endpoint"] = "RFS"
    c["expr"] = expr
    c["keys"] = pd.Series({s: "MB:" + s for s in cols})
    c["notes"].append("RFS per cBioPortal: event = loco-regional relapse, distant relapse or disease-specific death "
                      "(RFS_STATUS '1:Recurred'), time RFS_MONTHS")
    return c


def load_tcga(a):
    c = new_cohort("TCGA_BRCA", "UCSC Xena TCGA BRCA HiSeqV2 (local BRCA.gz) + TCGA-CDR (Liu 2018) survival table",
                   "RNA-seq Illumina HiSeq (log2(RSEM normalized count + 1))")
    expr = pd.read_csv(a.tcga_expr, sep="\t", index_col=0, keep_default_na=False)
    expr.columns = [x.replace(".", "-") for x in expr.columns]
    prim = [x for x in expr.columns if x[13:15] == "01"]
    n_all = expr.shape[1]
    expr = expr[prim].astype(np.float64)
    pts = pd.Series([x[:12] for x in prim])
    assert not pts.duplicated().any()
    mapper = SymbolMapper(a.annot)
    ent, how = mapper.map(list(expr.index))
    genes = collapse_max_mean(expr.reset_index(drop=True), ent)
    c["n_probe_note"] = f"{len(expr)} gene symbols -> Entrez via org.Hs.eg.db {how}; -> {len(genes)} genes (max-mean)"
    genes, c["log_note"] = maybe_log2(genes)
    c["expr"] = genes
    s = pd.read_csv(a.tcga_surv, sep="\t", low_memory=False)
    s = s[s["cancer type abbreviation"] == "BRCA"].copy()
    s = s[s["sample"].isin(prim)].set_index("sample").reindex(prim)
    red = s["Redaction"].eq("Redacted")
    counts = {}
    for ep in ("DFI", "DSS", "OS", "PFI"):
        tt = pd.to_numeric(s[ep + ".time"], errors="coerce"); ee = pd.to_numeric(s[ep], errors="coerce")
        ok = tt.notna() & ee.notna() & (tt > 0) & ~red
        counts[ep] = (int(ok.sum()), int(ee[ok].sum()))
    c["notes"].append(f"{n_all} Xena columns, {len(prim)} primary tumours (sample type 01); {int(red.sum())} redacted "
                      f"patients dropped; usable n/events by endpoint {counts}")
    tt = pd.to_numeric(s["DFI.time"], errors="coerce").where(~red)
    ee = pd.to_numeric(s["DFI"], errors="coerce").where(~red)
    if counts["DFI"][0] >= MIN_N and counts["DFI"][1] >= MIN_EVENTS:
        c["clin"], dr = make_clin(prim, tt.values, ee.values, "days", "DFI")
        c["endpoint"] = "DFI"
        c["notes"].append("endpoint DFI (TCGA-CDR disease-free interval; defined only for patients disease-free "
                          "after initial treatment) as the RFS/DFS-class endpoint; DMFS not available")
    else:
        c["clin"], dr = make_clin(prim, s["DSS.time"].where(~red).values, s["DSS"].where(~red).values, "days", "DSS")
        c["endpoint"] = "DSS"
    c["keys"] = pd.Series({x: "TCGA:" + x[:12] for x in prim})
    return c


# --------------------------------------------------------------------------- duplicates / overlaps
def centred_top(expr, genes=None, ntop=3000):
    x = expr if genes is None else expr.reindex(genes)
    x = x[np.isfinite(x).all(axis=1)]
    x = x.sub(x.median(axis=1), axis=0)
    if genes is None:
        x = x.loc[x.var(axis=1).sort_values(ascending=False).index[:ntop]]
    return x


def zcols(x):
    v = x.values.astype(np.float64)
    v = (v - v.mean(0)) / v.std(0)
    return v


def within_duplicates(c):
    """Pairs of arrays within a cohort (all its expression samples, with or without outcome) with r > DUP_R
    on gene-centred top-variance genes."""
    x = centred_top(c["expr"])
    z = zcols(x)
    r = (z.T @ z) / z.shape[0]
    np.fill_diagonal(r, -1)
    iu = np.argwhere(np.triu(r > DUP_R, 1))
    cols = list(x.columns)
    nn = np.sort(r, axis=1)[:, -1]
    return [(cols[i], cols[j], float(r[i, j])) for i, j in iu], float(np.median(nn)), float(np.quantile(nn, 0.99))


def cross_diagnostic(cohorts, out):
    """Report cross-cohort sample pairs with high correlation (gene-centred within cohort, shared genes)."""
    arr = [c for c in cohorts if c["expr"] is not None and "Affymetrix" in c["platform"]]
    shared = sorted(set.intersection(*[set(c["expr"].index) for c in arr]))
    blocks, labs = [], []
    for c in arr:
        x = c["expr"].reindex(shared)
        x = x.sub(x.median(axis=1), axis=0)
        blocks.append(x)
    allx = pd.concat(blocks, axis=1)
    sd = pd.concat([b.std(axis=1).rank(pct=True) for b in blocks], axis=1).median(axis=1)
    sel = sd.sort_values(ascending=False).index[:3000]
    for c, b in zip(arr, blocks):
        labs += [(c["id"], s, c["keys"].get(s, "")) for s in b.columns]
    z = zcols(allx.loc[sel].fillna(0))
    r = (z.T @ z) / z.shape[0]
    coh = np.array([l[0] for l in labs])
    rows = []
    for i in range(len(labs)):
        mask = coh > coh[i]
        js = np.where(mask & (r[i] > CROSS_R_REPORT))[0]
        for j in js:
            rows.append(dict(cohort_a=labs[i][0], sample_a=labs[i][1], key_a=labs[i][2], cohort_b=labs[j][0],
                             sample_b=labs[j][1], key_b=labs[j][2], r=round(float(r[i, j]), 4),
                             same_key=labs[i][2] == labs[j][2]))
    # also r for every key-matched pair, whatever its value
    kidx = collections.defaultdict(list)
    for i, l in enumerate(labs):
        kidx[l[2]].append(i)
    for k, ii in kidx.items():
        for p in range(len(ii)):
            for q in range(p + 1, len(ii)):
                i, j = ii[p], ii[q]
                if coh[i] != coh[j] and r[i, j] <= CROSS_R_REPORT:
                    rows.append(dict(cohort_a=labs[i][0], sample_a=labs[i][1], key_a=k, cohort_b=labs[j][0],
                                     sample_b=labs[j][1], key_b=k, r=round(float(r[i, j]), 4), same_key=True))
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    np.fill_diagonal(r, np.nan)
    other_max = np.nanmax(np.where(coh[:, None] != coh[None, :], r, np.nan), axis=1)
    log(f"cross-cohort diagnostic over {len(arr)} Affymetrix cohorts, {len(shared)} shared genes, 3000 used; "
        f"median of each sample's best cross-cohort r = {np.nanmedian(other_max):.3f}; pairs reported: {len(df)}")
    return df


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default="/path/to/analysis_root/breast")
    ap.add_argument("--tcga-expr", default="/path/to/analysis_root/pancan/expr/BRCA.gz")
    ap.add_argument("--tcga-surv", default="/path/to/analysis_root/pancan/Survival_SupplementalTable_S1_20171025_xena_sp")
    ap.add_argument("--skip-write", action="store_true", help="build everything but do not write cohort files")
    a = ap.parse_args()
    a.venet = f"{a.base}/download_venet/ploscb-venet-dumont-detours/data"
    a.geo = f"{a.base}/download_geo"
    a.metabric = f"{a.base}/download_metabric"
    a.annot = f"{a.base}/annotation"
    a.out = f"{a.base}/cohorts"
    import os
    os.makedirs(a.out, exist_ok=True)

    geo_ids = ["GSE2034", "GSE7390", "GSE11121", "GSE2990", "GSE1456", "GSE3494", "GSE4922", "GSE20685",
               "GSE21653", "GSE42568", "GSE17705", "GSE12093"]
    cohorts = []
    for fn in (load_nki, load_loi, load_tcga):
        log("loading", fn.__name__)
        cohorts.append(fn(a))
    for g in geo_ids:
        log("loading", g)
        cohorts.append(load_geo(a, g))
    log("loading GSE96058")
    cohorts.append(load_gse96058(a))
    log("loading METABRIC")
    cohorts.append(load_metabric(a))
    by = {c["id"]: c for c in cohorts}

    # GSE17705: JBI-profiled samples carry institutional ids; resolve against Loi/GSE6532 sample names
    loi_names = set(by["Loi_GSE6532"]["expr"].columns) | set(
        pd.read_csv(f"{a.geo}/GSE6532/GSE6532_LUMINAL_demo.txt.gz", sep="\t")["samplename"])
    k = by["GSE17705"]["keys"].copy()
    nmatch = collections.Counter()
    for s, v in k.items():
        lab, pid = v.split(":", 1)
        if lab != "JBI":
            continue
        hits = [n for n in (f"OXFT_{pid}", f"KIT_{pid}", f"GUYT_{pid}") if n in loi_names]
        if len(hits) == 1:
            k[s] = loi_key(hits[0]); nmatch[hits[0].split("_")[0]] += 1
        elif len(hits) > 1:
            nmatch["ambiguous"] += 1
        else:
            nmatch["unmatched"] += 1
    by["GSE17705"]["keys"] = k
    log("GSE17705 JBI ids matched to GSE6532 sample names:", dict(nmatch))
    by["GSE17705"]["notes"].append(f"JBI-profiled patient ids match GSE6532 names {dict(nmatch)}")

    # within-cohort: duplicated patient keys, then near-identical arrays
    dup_rows = []
    for c in cohorts:
        if c["expr"] is None:
            continue
        keys = c["keys"]
        if c["clin"] is not None:
            cl = c["clin"]
            cl = cl[cl["sample"].isin(c["expr"].columns)]
            kk = cl["sample"].map(keys)
            dk = kk.duplicated(keep="first")
            if dk.any():
                c["notes"].append(f"{int(dk.sum())} samples with an already-seen patient key dropped")
                cl = cl[~dk.values]
            c["clin"] = cl.reset_index(drop=True)
        if c["clin"] is None or len(c["clin"]) < 3:
            continue
        pairs, med_nn, q99 = within_duplicates(c)
        log(f"{c['id']}: within-cohort nearest-neighbour r median {med_nn:.3f}, 99th pct {q99:.3f}; pairs r>{DUP_R}: {len(pairs)}")
        drop = set()
        for s1, s2, r in pairs:
            same = keys.get(s1) == keys.get(s2)
            dup_rows.append(dict(cohort=c["id"], sample_a=s1, sample_b=s2, r=round(r, 4), same_key=same,
                                 a_has_endpoint=s1 in set(c["clin"]["sample"]), b_has_endpoint=s2 in set(c["clin"]["sample"])))
            if not same:
                drop |= {s1, s2}
        drop &= set(c["clin"]["sample"])
        if drop:
            c["clin"] = c["clin"][~c["clin"]["sample"].isin(drop)].reset_index(drop=True)
            c["notes"].append(f"near-identical arrays (r>{DUP_R}) annotated as different patients: the member(s) with "
                              f"an endpoint dropped ({sorted(drop)}; see duplicates_within.csv)")
    pd.DataFrame(dup_rows, columns=["cohort", "sample_a", "sample_b", "r", "same_key", "a_has_endpoint",
                                    "b_has_endpoint"]).to_csv(
        f"{a.out}/duplicates_within.csv", index=False)

    cross = cross_diagnostic(cohorts, f"{a.out}/overlap_diagnostic.csv")

    # cross-cohort overlap: patients kept once, in the larger (eligible n) included series
    def eligible(c):
        return 0 if c["clin"] is None or c["reason"] else len(c["clin"])
    for c in cohorts:
        c["n_before_overlap"] = eligible(c)
    order = sorted(cohorts, key=lambda c: -eligible(c))
    claimed = {}
    kept_for = collections.defaultdict(collections.Counter)
    for c in order:
        if c["clin"] is None or c["reason"]:
            continue
        kk = c["clin"]["sample"].map(c["keys"])
        hit = kk.isin(claimed.keys())
        if hit.any():
            src = collections.Counter(claimed[x] for x in kk[hit])
            for k2, v2 in src.items():
                kept_for[k2][c["id"]] += v2
            c["notes"].append(f"{int(hit.sum())} patients also in larger series {dict(src)} removed here (kept there)")
            c["clin"] = c["clin"][~hit.values].reset_index(drop=True)
        n, ev = len(c["clin"]), int(c["clin"]["event"].sum())
        x = c["expr"][list(c["clin"]["sample"])]
        n_ok = int((np.isfinite(x).all(axis=1) & (x.std(axis=1) > 0)).sum())
        c["n_genes_ok"] = n_ok
        if n == 0 and hit.any():
            c["reason"] = "all patients with the endpoint are contained in a larger included series"
        elif n < MIN_N or ev < MIN_EVENTS or n_ok < MIN_GENES:
            c["reason"] = f"fails inclusion: n={n} (min {MIN_N}), events={ev} (min {MIN_EVENTS}), genes={n_ok} (min {MIN_GENES})"
        if not c["reason"]:
            for s, key in zip(c["clin"]["sample"], c["clin"]["sample"].map(c["keys"])):
                claimed[key] = c["id"]

    for cid, cnt in kept_for.items():
        by[cid]["notes"].append(f"patients also present in smaller series {dict(cnt)} are kept here and removed there")
    current = set(pd.read_csv(f"{a.annot}/orgHs_symbol2entrez.tsv", sep="\t", dtype=str)["gene_id"].astype(np.int64))
    for c in cohorts:
        if c["expr"] is not None:
            nc = int((~c["expr"].index.isin(list(current))).sum())
            c["notes"].append(f"{nc} of {len(c['expr'])} Entrez IDs are not current in org.Hs.eg.db (EGSOURCEDATE "
                              "2025-Sep24; discontinued/replaced IDs left as deposited)")
        if c["clin"] is not None:
            c["notes"].append(f"{c['n_before_overlap']} patients with the endpoint before cross-series overlap removal")

    # write
    man = []
    for c in cohorts:
        incl = c["reason"] is None and c["clin"] is not None
        n = len(c["clin"]) if c["clin"] is not None else 0
        ev = int(c["clin"]["event"].sum()) if c["clin"] is not None else 0
        ng = len(c["expr"]) if c["expr"] is not None else 0
        notes = "; ".join([c["log_note"], c["n_probe_note"]] + c["notes"] +
                          ([f"{c.get('n_genes_ok')} genes finite and non-constant in final samples"] if c.get("n_genes_ok") else []))
        man.append(dict(cohort=c["id"], source=c["source"], platform=c["platform"], n_samples=n, n_events=ev,
                        endpoint=c["endpoint"] or "", n_genes=ng, included=incl,
                        reason=c["reason"] or "meets inclusion (>=100 patients, >=25 events, >=10,000 genes)",
                        notes=notes))
        if incl and not a.skip_write:
            cl = c["clin"][["sample", "time_months", "event", "endpoint"]]
            cl.to_csv(f"{a.out}/{c['id']}_clin.tsv", sep="\t", index=False, float_format="%.4f")
            x = c["expr"][list(cl["sample"])]
            x.index.name = "entrez"
            x.to_csv(f"{a.out}/{c['id']}_expr.tsv.gz", sep="\t", float_format="%.6g", na_rep="NA",
                     compression={"method": "gzip", "compresslevel": 6})
        log(f"{c['id']}: included={incl} n={n} events={ev} endpoint={c['endpoint']} genes={ng} reason={c['reason']}")
    pd.DataFrame(man).to_csv(f"{a.out}/manifest.csv", index=False)
    with open(f"{a.out}/prep_log.txt", "w") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
