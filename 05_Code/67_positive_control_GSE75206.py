# ============================================================
# 67_positive_control_GSE75206.py —— B-6：正对照分析
# GSE75206：小鼠海马 BaP 70 mg/kg/天 × 3 天 vs 橄榄油对照（4 vs 4）
# 流程与主分析一致：探针→基因（Ensembl Agilent-028005 映射）→
# 基因水平中位数折叠 → Welch t 检验（BH）→ 与 B-1 69 个可评估靶点取交集。
# 输出：results/16_Orthology/positive_control_GSE75206/
# ============================================================
import csv, gzip, os, statistics

ROOT = os.getcwd()
OUT = os.path.join(ROOT, "results", "16_Orthology", "positive_control_GSE75206")
os.makedirs(OUT, exist_ok=True)

MATRIX = os.path.join(ROOT, "data", "01_GEO_bulk", "GSE75206_series_matrix.txt.gz")
PROBE_MAP = os.path.join(ROOT, "data", "01_GEO_bulk", "GPL10787_ensembl_probe_map.csv")
EVAL = os.path.join(ROOT, "results", "16_Orthology", "evaluable_targets_CLEAN.csv")
TIER = os.path.join(ROOT, "results", "16_Orthology", "bap_targets_364_tier.csv")

# ---------- 1) 读取系列矩阵 ----------
import pandas as pd
from scipy import stats

with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as f:
    lines = f.readlines()
title_line = next(l for l in lines if l.startswith("!Sample_title"))
title_vals = [x.strip('"\n') for x in title_line.rstrip("\n").split("\t")][1:]
gse_line = next(l for l in lines if l.startswith("!Sample_geo_accession"))
gse_vals = [x.strip('"\n') for x in gse_line.rstrip("\n").split("\t")][1:]
gsm2title = dict(zip(gse_vals, title_vals))
start = next(i for i, l in enumerate(lines) if l.startswith('"ID_REF"'))
end = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
header = [x.strip('"\n') for x in lines[start].rstrip("\n").split("\t")]
rows = []
for l in lines[start+1:end]:
    parts = [x.strip('"') for x in l.rstrip("\n").split("\t")]
    rows.append(parts)
mx = pd.DataFrame(rows, columns=header)
mx["ID_REF"] = mx["ID_REF"].str.replace('"', "")
for c in header[1:]:
    mx[c] = pd.to_numeric(mx[c], errors="coerce")
print("matrix:", mx.shape)

# ---------- 2) 探针→基因 ----------
pm = pd.read_csv(PROBE_MAP)
pm.columns = ["gene_symbol", "probe_id"]
merged = mx.merge(pm, left_on="ID_REF", right_on="probe_id", how="inner")
print("probes mapped:", merged.shape[0])

sample_cols = header[1:]
grp = ["BaP" if "BaP" in gsm2title.get(c, "") else "Control" for c in sample_cols]
print("groups:", dict(zip(sample_cols, grp)))

# 基因水平（中位数折叠）
gene_expr = merged.groupby("gene_symbol")[sample_cols].median()
gene_expr.to_csv(os.path.join(OUT, "gene_expression_median.csv"))

# ---------- 3) Welch t 检验 + BH ----------
res = []
for g, row in gene_expr.iterrows():
    x_bap = row[[c for c, gr in zip(sample_cols, grp) if gr == "BaP"]].astype(float)
    x_ctrl = row[[c for c, gr in zip(sample_cols, grp) if gr == "Control"]].astype(float)
    diff = x_bap.mean() - x_ctrl.mean()
    p = stats.ttest_ind(x_bap, x_ctrl, equal_var=False).pvalue
    res.append({"gene": g, "log2FC_BaP_vs_Control": diff,
                "welch_p": p, "n_bap": int(x_bap.notna().sum()),
                "n_control": int(x_ctrl.notna().sum())})
deg = pd.DataFrame(res)
deg["BH_p"] = stats.false_discovery_control(deg["welch_p"].fillna(1), method="bh")
deg = deg.sort_values("welch_p")
deg.to_csv(os.path.join(OUT, "deg_results.csv"), index=False)

# ---------- 4) 与 69 可评估靶点交集 + 经典 BaP/AhR 基因状态 ----------
ev = pd.read_csv(EVAL)
ev_set = set(ev["matrix_symbol"].str.upper())
deg["is_evaluable_target"] = deg["gene"].str.upper().isin(ev_set)

cand_deg = deg[(deg["BH_p"] < 0.05) & (deg["is_evaluable_target"])]
canonical = ["Ahr", "Cyp1a1", "Cyp1a2", "Cyp1b1", "Ahrr", "Nqo1", "Nfe2l2",
             "Gstm1", "Maoa", "Hprt"]
canon = deg[deg["gene"].str.lower().isin([x.lower() for x in canonical])]

summary = []
summary.append("B-6 positive control: GSE75206")
summary.append(f"matrix probes: {mx.shape[0]}; mapped to genes: {merged.shape[0]}")
summary.append(f"gene-level features: {gene_expr.shape[0]}")
summary.append(f"genes with BH<0.05 (BaP vs control): {(deg['BH_p']<0.05).sum()}")
summary.append(f"evaluable-target DEGs (BH<0.05): {len(cand_deg)}")
summary.append("candidate DEGs among 69 evaluable targets:")
summary.append(cand_deg[["gene","log2FC_BaP_vs_Control","welch_p","BH_p"]].to_string(index=False))
summary.append("\ncanonical BaP/AhR battery status:")
summary.append(canon[["gene","log2FC_BaP_vs_Control","welch_p","BH_p"]].to_string(index=False))
summary.append("\ntier A targets present in this dataset:")
ta = deg[deg["gene"].str.upper().isin(["AHR","CYP1A2","CYP1B1","MAOA","HPRT"])]
summary.append(ta[["gene","log2FC_BaP_vs_Control","welch_p","BH_p"]].to_string(index=False))

with open(os.path.join(OUT, "positive_control_summary.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(summary))
cand_deg.to_csv(os.path.join(OUT, "candidate_deg_intersection.csv"), index=False)
canon.to_csv(os.path.join(OUT, "canonical_baP_genes.csv"), index=False)
print("\n".join(summary))
