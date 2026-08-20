# ============================================================
# 81_review_supplements.py —— BMC 投稿前评审 M2/M3/M4/M5 补充分析
# 输出（results/16_Orthology/review_supplements/）：
#   - S2_permutation_v7_meanvar.csv  主分析（5,043）均值-方差分层置换
#   - permutation_strata_counts.csv  分层变量分布
#   - power_GSE75206.csv             正对照功效（n=4 vs 4）
#   - S6_validation_ci_v7.csv        逐基因 Welch 95% CI（11 候选）
#   - S13_sample_manifest.csv        样本级清单（GSE58485+GSE41345）
#   - S14_full_universe_validation.csv 完整宇宙 6 候选方向一致性
# 全部为只读输入上的确定性计算（seed=2026）。
# ============================================================
import csv
import math
import os
import random

import numpy as np
import pandas as pd

ROOT = r"E:\sheng xin\tbi"
OUT = os.path.join(ROOT, "results", "16_Orthology", "review_supplements")
os.makedirs(OUT, exist_ok=True)
SEED = 2026
N_PERM = 10000


def read_csv_p(path):
    return pd.read_csv(path, encoding="utf-8-sig")


# ---------------- M2: primary (5,043) mean-variance stratified permutation ----
# 使用 80_wgcna_diagnostics_primary.R 导出的真实主分析矩阵（未校正，
# 与 70_single_cohort_main_v7.R prep() 完全一致），而非 ComBat 校正矩阵。
expr = read_csv_p(os.path.join(ROOT, "results", "16_Orthology", "v7_single_cohort",
                               "expression_GSE58485_only_primary.csv"))
expr = expr.set_index(expr.columns[0])
expr.columns = [str(c) for c in expr.columns]
gse_cols = [c for c in expr.columns if c.startswith("GSM1412")]
expr24 = expr[gse_cols]

key = pd.read_csv(os.path.join(ROOT, "results", "16_Orthology", "v7_single_cohort",
                               "key_genes_GSE58485_only.csv"))["gene"].tolist()
cand = pd.read_csv(os.path.join(ROOT, "results", "16_Orthology", "v7_single_cohort",
                                "candidates_GSE58485_only.csv"))["gene"].tolist()
universe = list(expr24.index)
key_set = set(key)

gene_mean = expr24.mean(axis=1)
gene_var = expr24.var(axis=1, ddof=1)


def decile(vals):
    q = vals.quantile(np.arange(0, 1.001, 0.1))
    q.iloc[0] -= 1e-12
    return pd.cut(vals, bins=q, labels=False, include_lowest=True)


mean_dec = decile(gene_mean)
var_dec = decile(gene_var)
strata_df = pd.DataFrame({"gene": universe, "mean_dec": mean_dec, "var_dec": var_dec})
strata_counts = strata_df.pivot_table(index="mean_dec", columns="var_dec",
                                      values="gene", aggfunc="count")
strata_counts.to_csv(os.path.join(OUT, "permutation_strata_counts.csv"))

rng = random.Random(SEED)
K = len(cand)


def perm_uniform():
    overlaps = []
    for _ in range(N_PERM):
        s = rng.sample(universe, K)
        overlaps.append(sum(1 for g in s if g in key_set))
    return np.array(overlaps)


def perm_meanvar():
    cand_strata = strata_df[strata_df["gene"].isin(cand)]
    pools = []
    for _, row in cand_strata.iterrows():
        pool = strata_df["gene"][
            (strata_df["mean_dec"] == row["mean_dec"]) &
            (strata_df["var_dec"] == row["var_dec"])
        ].tolist()
        pools.append(pool)
    overlaps = []
    for _ in range(N_PERM):
        s = [rng.choice(p) for p in pools]
        overlaps.append(sum(1 for g in s if g in key_set))
    return np.array(overlaps)


def report(path, label, overlaps, exp):
    obs = K
    emp_p = (1 + int((overlaps >= obs).sum())) / (N_PERM + 1)
    return {
        "path": label, "K": K, "universe": len(universe),
        "key_genes": len(key_set), "observed": obs,
        "expected": round(exp, 6), "enrichment": round(obs / exp, 3) if exp else None,
        "empirical_P": emp_p,
        "null_95_lo": int(np.quantile(overlaps, 0.025)),
        "null_95_hi": int(np.quantile(overlaps, 0.975)),
    }
    print(label, "obs", obs, "expected", round(exp, 4), "enrich", round(obs / exp, 2),
          "empP", emp_p)


ov_unif = perm_uniform()
ov_mv = perm_meanvar()
exp_unif = K * len(key_set) / len(universe)
cs = strata_df[strata_df["gene"].isin(cand)]
exp_mv = sum(
    len(strata_df["gene"][
        (strata_df["mean_dec"] == r["mean_dec"]) &
        (strata_df["var_dec"] == r["var_dec"]) &
        strata_df["gene"].isin(key_set)
    ]) / max(1, len(strata_df["gene"][
        (strata_df["mean_dec"] == r["mean_dec"]) &
        (strata_df["var_dec"] == r["var_dec"])
    ]))
    for _, r in cs.iterrows()
) if len(cand) else 0
rows_perm = [report(os.path.join(OUT, "S2_permutation_v7_meanvar.csv"), "uniform",
                    ov_unif, exp_unif),
             report(os.path.join(OUT, "S2_permutation_v7_meanvar.csv"),
                    "meanvar_stratified", ov_mv, exp_mv)]
pd.DataFrame(rows_perm).to_csv(os.path.join(OUT, "S2_permutation_v7_meanvar.csv"),
                               index=False)

# ---------------- M3: power for GSE75206 (n=4 vs 4) ----------------
n = 4
alpha = 0.05
power = 0.80
t_alpha = 1.943  # t_{0.975, df=6} for two-sample with n=4+4
t_beta = 0.889  # t_{0.80, df=6}
se_factor = math.sqrt(2.0 / n)
min_d = (t_alpha + t_beta) * se_factor  # in pooled-SD units

# empirical SD anchor from the primary analysis matrix (log2 units)
sd_median = float(gene_var.apply(math.sqrt).median())
sd_q75 = float(gene_var.apply(math.sqrt).quantile(0.75))
rows = []
for sigma in [0.5, 1.0, 1.5, sd_median, sd_q75]:
    rows.append({
        "n_per_group": n, "alpha": alpha, "power": power,
        "pooled_sd_log2": round(sigma, 4),
        "min_detectable_log2FC": round(min_d * sigma, 4),
        "anchor": "assumed" if sigma in (0.5, 1.0, 1.5) else
        ("primary_matrix_median_SD" if sigma == sd_median else "primary_matrix_Q75_SD"),
    })
pd.DataFrame(rows).to_csv(os.path.join(OUT, "power_GSE75206.csv"), index=False)
print("power rows written; median SD =", round(sd_median, 3), "Q75 =", round(sd_q75, 3))

# ---------------- M4: validation per-gene Welch 95% CIs (primary 11) ----------------
cand_stat = pd.read_csv(os.path.join(ROOT, "results", "16_Orthology", "v7_single_cohort",
                                     "candidate_stats_GSE58485_only.csv"))
genes11 = cand_stat["gene"].str.upper().tolist()
disc = dict(zip(cand_stat["gene"].str.upper(), cand_stat["logFC"].astype(float)))


def welch_ci(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    nx, ny = len(x), len(y)
    mx, my = x.mean(), y.mean()
    vx, vy = x.var(ddof=1), y.var(ddof=1)
    se = math.sqrt(vx / nx + vy / ny)
    nu = (vx / nx + vy / ny) ** 2 / (
        (vx / nx) ** 2 / (nx - 1) + (vy / ny) ** 2 / (ny - 1)
    )
    t = _t_ppf(0.975, nu)
    diff = mx - my
    tstat = diff / se if se > 0 else float("nan")
    p = t_pval(tstat, nu) if se > 0 else float("nan")
    return diff, diff - t * se, diff + t * se, p


def _t_ppf(q, df):
    # approximate Student t quantile via scipy if available, else normal approx
    try:
        from scipy.stats import t as tdist
        return float(tdist.ppf(q, df))
    except Exception:
        from scipy.stats import norm
        z = norm.ppf(q)
        return z + (z ** 3 + z) / (4 * df)


def t_pval(tval, df):
    try:
        from scipy.stats import t as tdist
        return float(tdist.sf(abs(tval), df) * 2)
    except Exception:
        from scipy.stats import norm
        return float(norm.sf(abs(tval)) * 2)


# GSE128543
g128 = pd.read_csv(os.path.join(ROOT, "data", "01_GEO_bulk", "GSE128543_expr.csv"))
g128 = g128.set_index(g128.columns[0])
ph128 = pd.read_csv(os.path.join(ROOT, "data", "01_GEO_bulk", "GSE128543_pheno.csv"),
                    index_col=0, encoding="utf-8-sig")
title128 = ph128["title"].astype(str).str.lower()
grp128 = np.where(title128.str.contains(r"(?<![a-z])tbi(?![a-z])|\bcci\b", regex=True) &
                  ~title128.str.contains("sham"), "TBI",
                  np.where(title128.str.contains("sham"), "Sham", "other"))
grp128 = pd.Series(grp128, index=ph128.index)
is_wt = ph128.get("genotype:ch1", pd.Series(index=ph128.index)).astype(str).str.lower().str.strip().isin(["wt", "wildtype"]) | \
        title128.str.contains(r"\bwt\b|wildtype", regex=True)
keep128 = is_wt & grp128.isin(["TBI", "Sham"])
cols128 = [c for c in g128.columns if c in ph128.index[keep128]]

# GSE205958
fpkm = pd.read_excel(os.path.join(ROOT, "data", "01_GEO_bulk",
                                  "GSE205958_genes_fpkm_expression.xlsx"))
fpkm_cols = [c for c in fpkm.columns if str(c).startswith("FPKM.")]
assert len(fpkm_cols) == 10
fpkm5 = fpkm.set_index("gene_name")
fpkm5 = fpkm5[[str(c) for c in fpkm_cols]]
grp205 = ["TBI"] * 5 + ["Sham"] * 5


def gene_row(df, gene):
    idx = [i for i in df.index if str(i).upper() == gene]
    return df.loc[idx[0]] if idx else None


rows_ci = []
for g in genes11:
    for cohort, df, groups, cols in [
        ("GSE128543_WT_only", g128, grp128[keep128], cols128),
        ("GSE205958", fpkm5, grp205, fpkm5.columns),
    ]:
        r = gene_row(df, g)
        if r is None:
            rows_ci.append(dict(validation=cohort, gene=g, diff=np.nan, ci_lo=np.nan,
                                ci_hi=np.nan, t_p=np.nan, direction_consistent=np.nan))
            continue
        x_tbi = [float(r[c]) for c, grp in zip(cols, groups) if grp == "TBI"]
        x_sham = [float(r[c]) for c, grp in zip(cols, groups) if grp == "Sham"]
        diff, lo, hi, p = welch_ci(x_tbi, x_sham)
        rows_ci.append(dict(validation=cohort, gene=g, diff=round(diff, 4),
                            ci_lo=round(lo, 4), ci_hi=round(hi, 4),
                            t_p=round(p, 6),
                            direction_consistent=(math.copysign(1, diff) == math.copysign(1, disc[g]))))
pd.DataFrame(rows_ci).to_csv(os.path.join(OUT, "S6_validation_ci_v7.csv"), index=False)
print("CI rows:", len(rows_ci))

# ---------------- M5: sample manifest ----------------
comp = read_csv_p(os.path.join(ROOT, "results", "00_Data_Prep", "sample_composition.csv"))
ph584 = pd.read_csv(os.path.join(ROOT, "data", "01_GEO_bulk", "GSE58485_pheno.csv"),
                    index_col=0, encoding="utf-8-sig")
ph413 = pd.read_csv(os.path.join(ROOT, "data", "01_GEO_bulk", "GSE41345_pheno.csv"),
                    index_col=0, encoding="utf-8-sig")


def exclusion_reason(row):
    if row["cohort"] == "GSE58485":
        if row["exclude_drug"]:
            return "cyclophosphamide i.p. treatment (drug-exposed)"
        return ""
    if row["exclude_drug"]:
        return "Ex-4 treatment (drug-exposed)"
    return ""


manifest = []
for _, r in comp.iterrows():
    ph = ph584 if r["cohort"] == "GSE58485" else ph413
    gt = ph.loc[r["sample"], "genotype:ch1"] if "genotype:ch1" in ph.columns else ""
    tx_col = ("treatment_protocol_ch1" if "treatment_protocol_ch1" in ph.columns
              else "treatment protocol:ch1" if "treatment protocol:ch1" in ph.columns else None)
    tx = ph.loc[r["sample"], tx_col] if tx_col else ""
    manifest.append({
        "gsm": r["sample"], "cohort": r["cohort"], "group": r["group"],
        "genotype": gt, "treatment_protocol": tx,
        "include_full": r["include_full"], "include_clean": r["include_clean"],
        "exclusion_reason": exclusion_reason(r),
    })
pd.DataFrame(manifest).to_csv(os.path.join(OUT, "S13_sample_manifest.csv"), index=False)
print("manifest rows:", len(manifest))

# ---------------- Full-universe 6 candidates: direction consistency ----------------
fu_cand = pd.read_csv(os.path.join(ROOT, "results", "16_Orthology", "v7_full_universe",
                                   "candidates_full_universe.csv"))
fu_genes = fu_cand["gene"].str.upper().tolist()
fu_deg = pd.read_csv(os.path.join(ROOT, "results", "16_Orthology", "v7_full_universe",
                                  "deg_all_full_universe.csv"))
fu_disc = dict(zip(fu_deg["gene"].str.upper(), fu_deg["logFC"].astype(float)))
fu_rows = []
for g in fu_genes:
    for cohort, df, groups, cols in [
        ("GSE128543_WT_only", g128, grp128[keep128], cols128),
        ("GSE205958", fpkm5, grp205, fpkm5.columns),
    ]:
        r = gene_row(df, g)
        if r is None:
            fu_rows.append(dict(validation=cohort, gene=g, diff=np.nan, t_p=np.nan,
                                direction_consistent=np.nan))
            continue
        x_tbi = [float(r[c]) for c, grp in zip(cols, groups) if grp == "TBI"]
        x_sham = [float(r[c]) for c, grp in zip(cols, groups) if grp == "Sham"]
        diff, lo, hi, p = welch_ci(x_tbi, x_sham)
        fu_rows.append(dict(validation=cohort, gene=g, diff=round(diff, 4),
                            t_p=round(p, 6),
                            direction_consistent=(math.copysign(1, diff) ==
                                                   math.copysign(1, fu_disc.get(g, 0)))))
pd.DataFrame(fu_rows).to_csv(os.path.join(OUT, "S14_full_universe_validation.csv"), index=False)
print("full-universe validation rows:", len(fu_rows))
print("done")
