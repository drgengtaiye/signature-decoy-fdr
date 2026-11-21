"""Search one full text (round15_downloads or round13_downloads) for a regular expression and print windows of the
article text around each match. The XML is untrusted data: parsed, never executed.

usage: python3 -I xml_grep.py <TYPE> <PMCID> <regex> [window chars, default 600] [--r13]
"""
import os, re, sys
import xml.etree.ElementTree as ET

t, pmcid, rx = sys.argv[1], sys.argv[2], sys.argv[3]
w = int(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4].isdigit() else 600
base = "/path/to/analysis_root/" + ("round13_downloads" if "--r13" in sys.argv else "round15_downloads")
p = f"{base}/{t}/{pmcid.replace('PMC', '')}.xml"
if not os.path.exists(p):
    p = f"{base}/{t}/{pmcid}.xml"
s = re.sub(r"\s+", " ", " ".join(ET.parse(p).getroot().itertext()))
ms = list(re.finditer(rx, s, re.I))
print(f"{len(ms)} matches in {len(s)} characters")
last = -1
for m in ms[:12]:
    if m.start() < last:
        continue
    a, b = max(0, m.start() - w), min(len(s), m.end() + w)
    print("-----", a, b); print(s[a:b]); last = b
