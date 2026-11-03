"""Round 15 A1-A2: apply the rules fixed in AMENDMENT_1_round15 to the readers' screening rows, before any new model
is scored, and write screen_r15_final.csv (every record, with the decision that stands and why it changed).

usage: python3 -I apply_rules_r15.py
"""
import csv

# rule 4: the main model is the one the title or abstract presents; if they present none or several, the first printed
# survival model (re-read for every row whose notes record a choice between models)
# rule 5: coefficients fitted to a non-survival endpoint; rule 6: a coefficient with more than one reading;
# rule 2: a gene whose NCBI gene type is not protein-coding; rule 7: one paper, one entry
CHANGES = {
    ("LUAD", "PMC9302251"): ("FALSE", "not trained on TCGA: main model (title, abstract) is the 7-gene GEO-cohort signature"),
    ("HNSC", "PMC9948032"): ("FALSE", "not a linear formula with per-gene coefficients: main model (title) is a composite of "
                                      "two within-cohort z-scored scores"),
    ("LUAD", "PMC8182055"): ("FALSE", "no prognostic model: coefficients fitted to a non-survival endpoint (lymph-node "
                                      "metastasis, LASSO-logistic)"),
    ("LUAD", "PMC12686163"): ("FALSE", "no prognostic model: coefficients fitted to a non-survival endpoint (diagnostic "
                                       "LASSO-logistic)"),
    ("COAD", "PMC9513614"): ("FALSE", "coefficient not legible: DLAT printed '08,399', more than one reading"),
    ("LUAD", "PMC9252528"): ("FALSE", "not protein-coding: SFTA3 (NCBI gene type ncRNA)"),
    ("LUAD", "PMC8697852"): ("FALSE", "not protein-coding: SFTA3 (NCBI gene type ncRNA)"),
    ("LUSC", "PMC9715950"): ("FALSE", "not protein-coding: DIRC1, SPATA8 (NCBI gene type ncRNA)"),
    ("LUAD", "PMC6953332"): ("FALSE", "already in the table: the same model is round-13 LUSC PMC6953332 (rule 7)"),
    ("LUSC", "PMC12321723"): ("FALSE", "already in the table: the same model enters under LUAD (rule 7)"),
}
# rule 3: a comma between digits of a printed number is a digit-group separator
COEF_FIX = {("LUAD", "PMC9730711"): ("0.105,367", "0.105367", "0.134,121", "0.134121")}
# rule 8: weights that are the paper's printed hazard ratios: scored as printed, left out in a sensitivity analysis
HR_WEIGHTS = {("STAD", "PMC7982000"), ("CESC", "PMC11319947")}

rows = list(csv.DictReader(open("screen_r15_all.csv")))
seen = set()
for r in rows:
    k = (r["type"], r["pmcid"])
    r["eligible_reader"], r["reason_reader"] = r["eligible"], r["reason"]
    r["changed_by_amendment"] = ""
    if k in CHANGES:
        assert r["eligible"].upper() == "TRUE", k
        r["eligible"], r["reason"] = CHANGES[k]; r["changed_by_amendment"] = "AMENDMENT_1_round15"; seen.add(k)
    if k in COEF_FIX:
        f = COEF_FIX[k]
        for a, b in zip(f[0::2], f[1::2]):
            assert a in r["coefs"], (k, a)
            r["coefs"] = r["coefs"].replace(a, b)
        r["changed_by_amendment"] = "AMENDMENT_1_round15"; seen.add(k)
    r["hr_weights"] = "TRUE" if k in HR_WEIGHTS else "FALSE"
    seen |= {k} & HR_WEIGHTS
assert seen == set(CHANGES) | set(COEF_FIX) | HR_WEIGHTS, set(CHANGES) | set(COEF_FIX) | HR_WEIGHTS - seen
cols = list(rows[0].keys())
with open("screen_r15_final.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
print(len(rows), "records;", sum(r["eligible"].upper() == "TRUE" for r in rows), "eligible after the amendment;",
      sum(r["changed_by_amendment"] != "" for r in rows), "rows changed")
