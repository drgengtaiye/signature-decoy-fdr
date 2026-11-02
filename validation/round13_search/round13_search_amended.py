"""Round 13 search as amended (AMENDMENT_1_round13_published_risk_models.md): PubMed, type and theme terms in title or
abstract, PMC open-access subset; PMIDs ordered by numpy default_rng(20261113 + type index); PMCIDs from esummary."""
import json, time, urllib.parse, urllib.request, numpy as np, os, sys
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TYPES = [("STAD", ['gastric cancer', 'stomach adenocarcinoma', 'gastric adenocarcinoma']), ("SKCM", ['melanoma']),
         ("LUAD", ['lung adenocarcinoma']), ("LUSC", ['lung squamous']),
         ("COAD", ['colon cancer', 'colon adenocarcinoma', 'colorectal cancer']),
         ("CESC", ['cervical cancer', 'cervical squamous']), ("SARC", ['sarcoma']),
         ("HNSC", ['head and neck squamous', 'HNSCC']), ("BLCA", ['bladder cancer', 'urothelial carcinoma']),
         ("OV", ['ovarian cancer']), ("BRCA", ['breast cancer'])]
THEME = ['ferroptosis', 'cuproptosis', 'pyroptosis', 'necroptosis', 'anoikis', 'disulfidptosis', 'autophagy', 'm6A',
         'N6-methyladenosine', 'hypoxia', 'glycolysis', 'senescence', 'immune-related', 'mitophagy',
         'epithelial-mesenchymal']
tag = lambda w: (f'"{w}"' if (" " in w or "-" in w) else w) + "[tiab]"
def post(ep, data):
    for k in range(5):
        try:
            r = urllib.request.urlopen(urllib.request.Request(EU + ep, urllib.parse.urlencode(data).encode()), timeout=120)
            time.sleep(0.4); return r.read()
        except Exception as e:
            time.sleep(3 * (k + 1)); err = e
    raise err
summary = {}
for ti, (tp, terms) in enumerate(TYPES):
    q = ('(TCGA OR "The Cancer Genome Atlas") AND ("risk score" OR "risk model" OR "prognostic signature" OR '
         '"prognostic model") AND (' + " OR ".join(tag(t) for t in terms) + ") AND (" +
         " OR ".join(tag(t) for t in THEME) + ') AND 2019:2025[dp] AND "pubmed pmc open access"[filter]')
    r = json.loads(post("esearch.fcgi", dict(db="pubmed", term=q, retmax=9999, retmode="json")))["esearchresult"]
    ids = sorted(int(x) for x in r["idlist"]); assert len(ids) == int(r["count"]) < 9999
    order = np.random.default_rng(20261113 + ti).permutation(np.array(ids, dtype=np.int64))
    pmc = {}
    for i in range(0, len(ids), 200):
        s = json.loads(post("esummary.fcgi", dict(db="pubmed", id=",".join(map(str, ids[i:i + 200])), retmode="json")))
        for u in s["result"]["uids"]:
            a = [x["value"] for x in s["result"][u].get("articleids", []) if x["idtype"] == "pmc"]
            pmc[int(u)] = a[0] if a else ""
    with open(f"{OUT}/{tp}_order.tsv", "w") as f:
        f.write("order_rank\tpmid\tpmcid\n")
        for k, p in enumerate(order, 1): f.write(f"{k}\t{p}\t{pmc.get(int(p), '')}\n")
    summary[tp] = dict(index=ti, count=len(ids), query=q, translation=r.get("querytranslation", ""),
                       with_pmcid=sum(1 for v in pmc.values() if v))
    print(tp, len(ids), summary[tp]["with_pmcid"], flush=True)
json.dump(summary, open(f"{OUT}/search_summary.json", "w"), indent=1)
