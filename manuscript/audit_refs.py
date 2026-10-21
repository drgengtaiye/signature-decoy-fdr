#!/usr/bin/env python3
"""Audit the manuscript bibliography against PubMed.

The author blocks are scrambled against the titles, but the in-text citations match the
TITLES (e.g. [73] is cited for Benjamini-Hochberg and ref 73's title is that paper), so the
titles anchor the numbering and only the author lists need replacing. Each reference is
identified by its journal/year/volume/first-page tail, which pins the article without
relying on the author/title boundary that is exactly what is corrupted.
"""
import json, re, time, urllib.parse, urllib.request, sys

DOCX = "/groups/adv2105_gp/yichen/Yi/self/shanghai/6gastirc/data/raw/submission_GM/MANUSCRIPT_GM1015.docx"
OUT  = "/groups/adv2105_gp/yichen/tmp/claude-3874/-groups-adv2105-gp-yichen-Yi-self-shanghai-6gastirc/30d20b73-d629-4090-b43c-c8ca139220a1/scratchpad/refs"
EU   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def get(url):
    for _ in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            time.sleep(1.5)
    return ""

from docx import Document
ps = [p.text.strip() for p in Document(DOCX).paragraphs]
refs = {}
for t in ps:
    m = re.match(r"^(\d{1,3})\.\s+(.*)$", t)
    if m and 1 <= int(m.group(1)) <= 96:
        refs[int(m.group(1))] = re.sub(r"\s+", " ", m.group(2)).strip()
print(f"parsed {len(refs)} references from the docx", file=sys.stderr)

TAIL = re.compile(r"\.\s*([^.]{2,60})\.\s*(\d{4})\s*;\s*([0-9]+)\s*(?:\(\d+\))?\s*:\s*(e?[0-9a-z]+)")
def parse_tail(t):
    ms = list(TAIL.finditer(t))
    if not ms: return None
    m = ms[-1]
    return dict(journal=m.group(1).strip(), year=m.group(2), vol=m.group(3), page=m.group(4),
                head=t[:m.start()].strip())

rows = []
for n in sorted(refs):
    raw = refs[n]
    tl = parse_tail(raw)
    rec = dict(ref=n, raw=raw, pmid="", pm_authors="", pm_title="", pm_journal="", status="")
    if not tl:
        rec["status"] = "no-journal-tail"; rows.append(rec); continue
    rec.update(journal=tl["journal"], year=tl["year"], vol=tl["vol"], page=tl["page"])
    q = f'{tl["journal"]}[ta] AND {tl["year"]}[dp] AND {tl["vol"]}[vi] AND {tl["page"]}[pg]'
    j = get(f'{EU}/esearch.fcgi?db=pubmed&retmode=json&retmax=3&term={urllib.parse.quote(q)}')
    time.sleep(0.4)
    ids = []
    try: ids = json.loads(j)["esearchresult"]["idlist"]
    except Exception: pass
    if not ids:
        MAN = {"J R Stat Soc B": "Benjamini Y, Hochberg Y",
               "J Am Stat Assoc": "Efron B",
               "J Bus Econ Stat": "Freedman D, Lane D"}
        for jk, au in MAN.items():
            if jk.lower() in tl["journal"].lower():
                rec["pm_authors"] = au; rec["pm_journal"] = tl["journal"]
                rec["status"] = "authors-OK" if raw.startswith(au.split()[0]) else "AUTHORS-WRONG"
                rec["pmid"] = "manual"
                break
        else:
            rec["status"] = "not-in-pubmed"
        rows.append(rec); continue
    rec["pmid"] = ids[0]
    x = get(f'{EU}/efetch.fcgi?db=pubmed&retmode=xml&id={ids[0]}')
    time.sleep(0.4)
    au = re.findall(r"<LastName>([^<]+)</LastName>\s*<ForeName>[^<]*</ForeName>\s*<Initials>([^<]*)</Initials>", x)
    if not au:
        au = [(a, b) for a, b in re.findall(r"<LastName>([^<]+)</LastName>.*?<Initials>([^<]*)</Initials>", x, re.S)]
    rec["pm_authors"] = ", ".join(f"{l} {i}" for l, i in au)
    mt = re.search(r"<ArticleTitle>(.*?)</ArticleTitle>", x, re.S)
    rec["pm_title"] = re.sub(r"<[^>]+>", "", mt.group(1)).strip() if mt else ""
    mj = re.search(r"<ISOAbbreviation>([^<]+)</ISOAbbreviation>", x)
    rec["pm_journal"] = mj.group(1) if mj else ""
    # does the manuscript's author block start with the authoritative first author?
    first = au[0][0] if au else ""
    rec["status"] = "authors-OK" if first and raw.lower().startswith(first.lower()) else "AUTHORS-WRONG"
    rows.append(rec)
    print(f"  {n:>3d} {rec['status']:<14s} pmid={rec['pmid']:<9s} {first}", file=sys.stderr)

json.dump(rows, open(f"{OUT}/audit.json", "w"), indent=1, ensure_ascii=False)
c = {}
for r in rows: c[r["status"]] = c.get(r["status"], 0) + 1
print("\n=== audit summary ===")
for k, v in sorted(c.items(), key=lambda kv: -kv[1]): print(f"  {k:<16s} {v}")
print(f"\nwrote {OUT}/audit.json")
