#!/usr/bin/env python
"""Checks for one redrawn sheet. Usage: python -I pixl_nc_check.py Fig9 [old PDF]
Prints size, density, fonts, type sizes, overlapping or bleeding text, and the numbers the old sheet printed that this one lost.
Writes figs_pub/pixl/<sheet>_v.png at 200 dpi: open it with Read after every build."""
import collections
import re
import sys
import fitz

D = "/groups/adv2105_gp/yichen/Yi/com/gastric/figs_pub"
OLD = {"Fig1": "Figure1_GB10", "Fig2": "Figure2_GB10", "Fig3": "Figure3_GB10", "Fig4": "Figure4_GB10", "Fig5": "Figure5_GB11",
       "Fig6": "Figure6_GB10", "Fig7": "Figure7_GB10", "Fig8": "Figure8_GB14", "Fig9": "Figure9_NC16"}
n = sys.argv[1]
p, q = f"{D}/pixl/{n}.pdf", (sys.argv[2] if len(sys.argv) > 2 else f"{D}/{OLD.get(n, n)}.pdf")
d = fitz.open(p); pg = d[0]
W, H = pg.rect.width / 72 * 25.4, pg.rect.height / 72 * 25.4
spans = []
for b in pg.get_text("dict")["blocks"]:
    for l in b.get("lines", []):
        for s in l["spans"]:
            if s["text"].strip():
                spans.append((fitz.Rect(s["bbox"]), s["text"], s["size"], s["font"], l["bbox"][1]))
letters = sorted({s[1].strip() for s in spans if abs(s[2] - 8) < 0.05 and "Bold" in s[3] and re.fullmatch(r"[a-z]", s[1].strip())})
print(f"== {n}: {W:.1f} x {H:.1f} mm {'OK' if W <= 180.5 and H <= 240.5 else 'TOO BIG, limit 180 x 240'}; pages {len(d)}")
print(f"  panels {len(letters)} ({''.join(letters)}); {len(letters) / (W * H / 10000):.1f} per 100 cm2 (target 3.0-4.4; 13-19 panels)")
print("  fonts", sorted({f[3] for f in pg.get_fonts()}))
sizes = collections.Counter()
for s in spans:
    sizes[round(s[2], 1)] += len(s[1])
print("  sizes", dict(sorted(sizes.items())), "| under 4.95:", sum(v for k, v in sizes.items() if k < 4.95))
over, bleed = [], []
for i, a in enumerate(spans):
    if a[0].x0 < 1 or a[0].y0 < 1 or a[0].x1 > pg.rect.width - 1 or a[0].y1 > pg.rect.height - 1:
        bleed.append(a[1])
    for b in spans[i + 1:]:
        r = a[0] & b[0]
        if r.is_empty or abs(a[4] - b[4]) < 0.5:
            continue
        # glyph boxes are taller than the ink: shrink by a quarter of the height before calling it a collision
        h = min(a[0].height, b[0].height)
        if r.height > 0.55 * h and r.width > 0.5:
            over.append((a[1], b[1]))
print(f"  text overlaps {len(over)}", over[:12])
print(f"  off-page {len(bleed)}", bleed[:8])
try:
    num = lambda x: collections.Counter(w[4].strip("(),;:") for w in x.get_text("words") if re.search(r"\d", w[4]))  # noqa: E731
    A, B = num(fitz.open(q)[0]), num(pg)
    print("  numbers only on the old sheet (must be axis ticks or key ends only):", dict(A - B))
    print("  numbers only here:", dict(list((B - A).items())[:80]))
except Exception as e:  # noqa: BLE001
    print("  no old sheet to compare:", e)
pg.get_pixmap(dpi=200).save(f"{D}/pixl/{n}_v.png")
print("  view", f"{D}/pixl/{n}_v.png")
