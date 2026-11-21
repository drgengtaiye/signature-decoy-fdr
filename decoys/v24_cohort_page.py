"""The validation-cohort lookup page (docs/validation_cohorts.html) from the table shipped with decoyfdr
(decoyfdr/inst/extdata/validation_cohorts.csv) and docs/validation_cohorts_template.html. Writes the page body (for
hosting as an artifact) and a standalone copy with its own document head.

    python3 -I v24_cohort_page.py <repo root>
"""
import json, sys
from pathlib import Path
import pandas as pd

R = Path(sys.argv[1])
d = pd.read_csv(R / "decoyfdr/inst/extdata/validation_cohorts.csv")
d["independent_set"] = d.independent_set.astype(bool)
rec = json.loads(d.to_json(orient="records"))
body = (R / "docs/validation_cohorts_template.html").read_text().replace("__DATA__", json.dumps(rec, separators=(",", ":")))
assert "__DATA__" not in body and len(rec) == 43
(R / "docs/validation_cohorts_body.html").write_text(body)
head = '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n<style>body{margin:0}[hidden]{display:none!important}</style>\n'
title_end = body.index("</style>") + len("</style>")
(R / "docs/validation_cohorts.html").write_text(head + body[:title_end] + "\n</head>\n<body>\n" + body[title_end:] + "\n</body>\n</html>\n")
print("written", R / "docs/validation_cohorts.html", R / "docs/validation_cohorts_body.html")
