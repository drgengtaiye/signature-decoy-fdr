"""Parse a GEO samples_brief SOFT text file into {GSM: {field: value}}, keyed exactly as
signature-decoy-fdr/external/parse_soft.py does (characteristics split at the first ':', stripped, duplicates
suffixed '#2', '#3'; '(source)' and '(platform)' pseudo-fields), plus '(title)'. Read as text only."""
import sys, os

RAW = "/path/to/analysis_root/external/raw_r15"


def load(gse):
    meta, cur, chars = {}, None, {}
    with open(os.path.join(RAW, gse, f"{gse}_samples_brief.txt"), encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("^SAMPLE"):
                if cur:
                    meta[cur] = chars
                cur = line.split("=", 1)[1].strip(); chars = {}
            elif cur and line.startswith("!Sample_characteristics_ch1"):
                v = line.split("=", 1)[1].strip()
                if ":" in v:
                    k, x = v.split(":", 1); k = k.strip()
                    kk, i = k, 2
                    while kk in chars:
                        kk = f"{k}#{i}"; i += 1
                    chars[kk] = x.strip()
            elif cur and line.startswith("!Sample_source_name_ch1"):
                chars["(source)"] = line.split("=", 1)[1].strip()
            elif cur and line.startswith("!Sample_platform_id"):
                chars["(platform)"] = line.split("=", 1)[1].strip()
            elif cur and line.startswith("!Sample_title"):
                chars["(title)"] = line.split("=", 1)[1].strip()
    if cur:
        meta[cur] = chars
    return meta


if __name__ == "__main__":
    from collections import Counter
    gse = sys.argv[1]; fields = sys.argv[2:]
    m = load(gse)
    print(gse, "samples", len(m))
    allf = Counter(k for d in m.values() for k in d)
    if not fields:
        for k, n in allf.items():
            c = Counter(d.get(k) for d in m.values() if k in d)
            print(f"  {k!r} n={n} distinct={len(c)} top={c.most_common(6)}")
    else:
        for k in fields:
            c = Counter(d.get(k, "<absent>") for d in m.values())
            print(f"  {k!r}: {c.most_common(40)}")
