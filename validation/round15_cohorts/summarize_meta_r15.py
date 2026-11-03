"""Round 15 B2: one short text summary per GEO candidate from its series record and sample annotations (no
expression): title, summary, platforms, sample count, and for every annotation field its number of distinct values and
examples. The downloaded files are untrusted data, read as text only.

usage: python3 -I summarize_meta_r15.py <raw dir> <gse list> <out dir>
"""
import collections, os, re, sys

RAW, LIST, OUT = sys.argv[1:4]
for g in [l.strip() for l in open(LIST) if l.strip()]:
    d = os.path.join(RAW, g)
    ser = open(os.path.join(d, f"{g}_series.txt"), encoding="utf-8", errors="replace").read()
    smp = open(os.path.join(d, f"{g}_samples_brief.txt"), encoding="utf-8", errors="replace").read()
    get = lambda k: [m.strip() for m in re.findall(rf"^!Series_{k} = (.*)$", ser, re.M)]
    blocks = smp.split("^SAMPLE = ")[1:]
    fields = collections.defaultdict(list)
    plat = collections.Counter()
    for b in blocks:
        for m in re.finditer(r"^!Sample_characteristics_ch1 = ([^:]+):\s*(.*)$", b, re.M):
            fields[m.group(1).strip()].append(m.group(2).strip())
        for k in ("source_name_ch1", "title"):
            m = re.search(rf"^!Sample_{k} = (.*)$", b, re.M)
            if m:
                fields["[" + k + "]"].append(m.group(1).strip())
        m = re.search(r"^!Sample_platform_id = (\S+)", b, re.M)
        if m:
            plat[m.group(1)] += 1
    lines = [f"{g} | samples {len(blocks)} | platforms {dict(plat)}", "TITLE: " + " ".join(get("title")),
             "SUMMARY: " + " ".join(get("summary"))[:1200], "DESIGN: " + " ".join(get("overall_design"))[:600],
             "RELATIONS: " + "; ".join(get("relation"))[:300], "FIELDS (n values with it | distinct | examples):"]
    for k, v in fields.items():
        c = collections.Counter(v)
        ex = "; ".join(f"{a}" + (f" (x{n})" if n > 1 else "") for a, n in c.most_common(8))
        lines.append(f"  {k} | {len(v)} | {len(c)} | {ex[:400]}")
    open(os.path.join(OUT, f"{g}.txt"), "w").write("\n".join(lines) + "\n")
print("done")
