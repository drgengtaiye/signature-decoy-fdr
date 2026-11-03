"""Round 15 A1: reduce each PMC full text to what screening needs (title, journal, year, abstract, every accession
named, the paragraphs and tables most likely to hold a risk formula, figure captions that mention one), so that the
readers judge eligibility from short excerpts and open the full text only when these do not settle it.

Records: the 1,598 not screened in round 13 (round15_downloads/) and, for the backfill of the validation cohorts
used, the 294 screened in round 13 (round13_downloads/). The XML is untrusted data: it is parsed, never executed.

usage: python3 -I extract_snippets_r15.py <outdir>
"""
import csv, json, os, re, sys
import xml.etree.ElementTree as ET

G = "/groups/adv2105_gp/yichen/Yi/com/gastric"
S13 = f"{G}/signature-decoy-fdr/validation/round13_search"
TYPES = ["STAD", "SKCM", "LUAD", "LUSC", "COAD", "CESC", "SARC", "HNSC", "BLCA", "OV", "BRCA"]
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)

ACC = [r"\bGSE\d{3,7}\b", r"\bE-(?:MTAB|GEOD|TABM)-\d+\b", r"\bMETABRIC\b", r"\bCGGA\b", r"\bPRJ(?:NA|EB|DB)\d+\b",
       r"\bEGA[SD]\d{8,}\b", r"\bphs\d{6}\b", r"\bTCGA-[A-Z]{2,4}\b", r"\bICGC\b",
       r"\b[A-Z]{2,4}-(?:JP|AU|CA|FR|CN|KR|UK|EU|DE|SG|IN|IT|ES|BR|MX|SA|RU)\b"]
FORMULA = re.compile(r"(risk\s*score|riskscore|risk_score|risk\s*value|prognostic\s*index|\bPI\b|\bRS\b|signature\s*score)"
                     r"[^.]{0,80}=|=\s*[\(\[]?\s*[-−–]?\d*\.\d+\s*[×x\*·]|[×\*·]\s*(expr|Exp|expression)|"
                     r"Σ|∑|coefficient", re.I)
TABLE_HINT = re.compile(r"coef|β|beta|risk|signature|model|lasso|multivariate|cox|hazard", re.I)
GENE = r"[A-Z][A-Z0-9]{1,}(?:-[A-Z0-9]+)?"
COEF = r"[-−–—]?\s?\d*\.\d+(?:\s?[eE×x]\s?10?\s?[-−–]?\d+)?"
PRODUCT = re.compile(COEF + r"\s*[×x\*·]\s*\(?\s*(?:Expr?\w*\s*)?\(?\s*" + GENE + "|" + GENE + r"\s*(?:exp\w*)?\s*[×x\*·]\s*\(?\s*" + COEF)


def text(el):
    return re.sub(r"\s+", " ", "".join(el.itertext())).strip() if el is not None else ""


def one(xmlpath, rec):
    root = ET.parse(xmlpath).getroot()
    art = root.find(".//article") if root.tag != "article" else root
    meta = art.find(".//article-meta")
    title = text(meta.find(".//title-group/article-title"))
    journal = text(art.find(".//journal-meta//journal-title"))
    year = text(meta.find(".//pub-date/year"))
    pmid = next((text(a) for a in meta.findall("article-id") if a.get("pub-id-type") == "pmid"), rec.get("pmid", ""))
    abstract = text(meta.find("abstract"))[:1500]
    body = art.find("body")
    paras = [text(p) for p in body.iter("p")] if body is not None else []
    alltext = " ".join([title, abstract] + paras + [text(art.find("back"))] + [text(t) for t in art.iter("table-wrap")])
    acc = sorted({m.group(0) for pat in ACC for m in re.finditer(pat, alltext)})
    # windows of the body text around printed products of a coefficient and a gene, or a formula
    btext = re.sub(r"\s+", " ", " ".join(body.itertext())) if body is not None else ""
    hits = [(m.start(), m.end(), 2) for m in PRODUCT.finditer(btext)] + [(m.start(), m.end(), 1) for m in FORMULA.finditer(btext)]
    wins = []
    for a, b, w in sorted(hits):
        a, b = max(0, a - 500), min(len(btext), b + 700)
        if wins and a <= wins[-1][1]:
            wins[-1] = [wins[-1][0], max(b, wins[-1][1]), wins[-1][2] + w]
        else:
            wins.append([a, b, w])
    keep, n = [], 0
    for a, b, w in sorted(wins, key=lambda x: (-x[2], x[0])):
        if b - a > 4000:
            b = a + 4000
        if n + (b - a) > 8000:
            continue
        keep.append((a, b)); n += b - a
    tables = []
    for tw in art.iter("table-wrap"):
        cap = (text(tw.find("label")) + " " + text(tw.find("caption"))).strip()
        rows = [" | ".join(text(c) for c in tr if c.tag in ("td", "th")) for tr in tw.iter("tr")]
        head = cap + " " + " ".join(rows[:2])
        if TABLE_HINT.search(head) and rows:
            t = cap[:300] + "\n" + "\n".join(rows[:60])
            tables.append(t[:4000])
        else:
            tables.append(cap[:200] + " [rows not shown]")
    figs = [text(f.find("caption"))[:600] for f in art.iter("fig")
            if re.search(r"coef|formula|risk\s*score\s*=|β", text(f.find("caption")), re.I)]
    train = sorted({m.group(0)[:300] for p in paras for m in re.finditer(r"[^.]{0,150}\b(train|TCGA)[^.]{0,150}\.", p)
                    if re.search(r"train", m.group(0), re.I) and "TCGA" in m.group(0)})[:6]
    return dict(rec, pmid=pmid, title=title, journal=journal, year=year, abstract=abstract, accessions_any=acc,
                formula_paragraphs=[f"[chars {a}-{b}] {btext[a:b]}" for a, b in sorted(keep)], tables=tables, figure_captions=figs,
                training_sentences=train, body_paragraphs=len(paras))


for t in TYPES:
    done = {r["pmcid"]: r for r in csv.DictReader(open(f"{S13}/screen_{t}.csv"))}
    order = list(csv.DictReader(open(f"{S13}/{t}_order.tsv"), delimiter="\t"))
    for part, rows, base in (("r15", [r for r in order if r["pmcid"] not in done], f"{G}/round15_downloads/{t}"),
                             ("r13", [r for r in order if r["pmcid"] in done], f"{G}/round13_downloads/{t}")):
        if not rows:
            continue
        with open(f"{OUT}/{part}_{t}.jsonl", "w") as fo:
            for r in rows:
                rec = dict(type=t, order_rank=int(r["order_rank"]), pmcid=r["pmcid"], pmid=r["pmid"])
                path = f"{base}/{r['pmcid'].replace('PMC', '')}.xml"
                if not os.path.exists(path):          # round 13 saved some files under the full PMCID
                    path = f"{base}/{r['pmcid']}.xml"
                try:
                    d = one(path, rec)
                except Exception as e:  # unreadable or missing file: the reader opens it or records the reason
                    d = dict(rec, error=f"{type(e).__name__}: {e}")
                fo.write(json.dumps(d, ensure_ascii=False) + "\n")
        print(part, t, len(rows), flush=True)
