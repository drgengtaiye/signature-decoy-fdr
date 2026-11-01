"""Sample annotations of GEO series through the GEO query page (brief view: no expression values), summarized per
series: platform(s), samples and the characteristic fields with their most common values.
usage: python -I screen_meta2.py <outdir> GSE..."""
import sys, os, re, time, collections, urllib.request
OUT = sys.argv[1]
def fetch(url, dest):
    if not os.path.exists(dest):
        data = urllib.request.urlopen(url, timeout=300).read()
        open(dest, "wb").write(data)
    return open(dest, encoding="utf-8", errors="ignore").read()
for acc in sys.argv[2:]:
    try:
        s = fetch(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=text&view=brief",
                  os.path.join(OUT, f"{acc}_self.txt"))
        g = fetch(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=gsm&form=text&view=brief",
                  os.path.join(OUT, f"{acc}_gsm.txt"))
    except Exception as e:
        print(f"## {acc} ERROR {e}"); continue
    title = re.findall(r"!Series_title = (.*)", s)
    plats = collections.Counter(re.findall(r"!Sample_platform_id = (\S+)", g))
    n = len(re.findall(r"^\^SAMPLE", g, re.M))
    fields = collections.defaultdict(list)
    for rec in g.split("^SAMPLE")[1:]:
        for k, v in re.findall(r"!Sample_characteristics_ch1 = ([^:\n]+):\s*(.*)", rec):
            fields[k.strip()].append(v.strip())
        src = re.findall(r"!Sample_source_name_ch1 = (.*)", rec)
        if src:
            fields["(source)"].append(src[0].strip())
    print(f"## {acc} n={n} {dict(plats)} | {title[0][:110] if title else ''}")
    for k, v in fields.items():
        c = collections.Counter(v)
        print(f"     {k[:38]:38s} n={len(v):4d} uniq={len(c):4d} {str(c.most_common(4))[:140]}")
    time.sleep(0.4)
