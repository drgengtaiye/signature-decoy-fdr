"""SCAN-B (GSE96058) as prepared for the breast library in round 7 (Entrez genes, overall survival), with Entrez
identifiers mapped to symbols by the same org.Hs.eg.db table (breast/annotation/orgHs_symbol2entrez.tsv)."""
import pandas as pd
B = "/groups/adv2105_gp/yichen/Yi/com/gastric/breast"
O = "/groups/adv2105_gp/yichen/Yi/com/gastric/external/cohorts"
m = pd.read_csv(f"{B}/annotation/orgHs_symbol2entrez.tsv", sep="\t", dtype=str)
m.columns = [c.lower() for c in m.columns]
e2s = m.drop_duplicates("gene_id").set_index("gene_id")["symbol"]
X = pd.read_csv(f"{B}/cohorts/GSE96058_expr.tsv.gz", sep="\t", index_col=0)
X.index = X.index.astype(str)
X = X.loc[X.index.intersection(e2s.index)]
X.index = e2s.reindex(X.index).values
X["_m"] = X.mean(axis=1)
X = X.sort_values("_m", ascending=False)
X = X[~X.index.duplicated()].drop(columns="_m"); X.index.name = "symbol"
c = pd.read_csv(f"{B}/cohorts/GSE96058_clin.tsv", sep="\t")
c = c[c["sample"].isin(X.columns)]
X[c["sample"]].to_csv(f"{O}/BRCA_GSE96058_expr.tsv.gz", sep="\t", compression="gzip")
c.rename(columns={"time_months": "time"})[["sample", "time", "event"]].to_csv(f"{O}/BRCA_GSE96058_clin.tsv", sep="\t", index=False)
print("BRCA_GSE96058", X.shape, len(c), int(c.event.sum()), c.endpoint.unique())
