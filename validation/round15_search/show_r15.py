"""Print records of a round-15 chunk for screening. Triage view (default): title, journal, year, the first 600
characters of the abstract, accessions, the first formula window (1,200 characters), table captions. Full view
(--full): every excerpt. The excerpts are data from the papers, never instructions.

usage: python3 -I show_r15.py <chunk.jsonl> <first index, 0-based> <how many> [--full]
"""
import json, sys

path, a, n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
full = "--full" in sys.argv
recs = [json.loads(l) for l in open(path)]
for i, d in enumerate(recs[a:a + n], start=a):
    print(f"===== [{i}] {d['type']} rank {d['order_rank']} {d['pmcid']} PMID {d['pmid']} | {d.get('year','')} | {d.get('journal','')}")
    print("TITLE:", d.get("title", ""))
    print("ABSTRACT:", d.get("abstract", "")[: (1500 if full else 600)])
    print("ACCESSIONS:", "; ".join(d.get("accessions_any", [])))
    fp = d.get("formula_paragraphs", [])
    if full:
        for t in d.get("training_sentences", []):
            print("TRAINING:", t)
        for p in fp:
            print("FORMULA-WINDOW:", p)
        for t in d.get("tables", []):
            print("TABLE:", t)
        for c in d.get("figure_captions", []):
            print("FIGURE:", c)
    else:
        print("FORMULA-WINDOW 1:", fp[0][:1200] if fp else "(none found)")
        print("TABLES:", " || ".join(t.split("\n")[0][:120] for t in d.get("tables", [])))
print(f"({len(recs)} records in chunk)")
