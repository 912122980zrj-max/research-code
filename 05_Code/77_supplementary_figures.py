# ============================================================
# 77_supplementary_figures_v7.py
# v7 缺失补充图的可复现生成：
#   S1a_algorithm_benchmark_v7.pdf  —— 11 算法 CV AUC（D3_ml_benchmark_v7/model_comparison.csv）
#   S2_direction_replication_v7.pdf —— 方向复制散点（validation_direction_v7.csv）
# 输入均为 v7 已产出结果文件；本脚本只做可视化，不重新计算统计量。
# 输出：figures/Figure_S_Supplementary/
# ============================================================
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = r"E:\sheng xin\tbi"
OUT = os.path.join(ROOT, "figures", "Figure_S_Supplementary")
os.makedirs(OUT, exist_ok=True)


def read_csv_rows(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return list(csv.reader(fh))


def spearman(x, y):
    """Spearman rank correlation without scipy dependency."""
    def rank(v):
        order = np.argsort(np.argsort(v))
        # average ranks for ties
        s = np.sort(v)
        r = np.empty(len(v))
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and s[j + 1] == s[i]:
                j += 1
            r[order[i : j + 1]] = (i + j) / 2.0 + 1
            i = j + 1
        return r

    rx, ry = rank(np.asarray(x, float)), rank(np.asarray(y, float))
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = np.sqrt((rx**2).sum() * (ry**2).sum())
    return float((rx * ry).sum() / denom) if denom > 0 else float("nan")


# ---------------- S1a: 11-algorithm CV AUC benchmark (v7) ----------------
comp_path = os.path.join(
    ROOT, "results", "07_ML_Signature", "D3_ml_benchmark_v7", "model_comparison.csv"
)
rows = read_csv_rows(comp_path)
header, data = rows[0], rows[1:]
algos = [r[0] for r in data]
aucs = np.array([float(r[2]) for r in data])
sds = np.array([float(r[3]) if r[3] else 0.0 for r in data])

order = np.argsort(aucs)[::-1]
algos_s = [algos[i] for i in order]
aucs_s = aucs[order]
sds_s = sds[order]

fig, ax = plt.subplots(figsize=(9, 4.2))
colors = ["#C0392B" if a < 1.0 else "#2E86C1" for a in aucs_s]
bars = ax.bar(range(len(algos_s)), aucs_s, yerr=sds_s, capsize=4,
              color=colors, edgecolor="black", linewidth=0.6)
ax.set_xticks(range(len(algos_s)))
ax.set_xticklabels(algos_s, rotation=35, ha="right", fontsize=9)
ax.set_ylim(0.90, 1.02)
ax.set_ylabel("CV AUC (repeated stratified 5-fold)")
ax.set_title("Eleven-algorithm benchmark, 11 candidates, n=24 (descriptive only)")
for b, a in zip(bars, aucs_s):
    ax.text(b.get_x() + b.get_width() / 2, a + 0.004, f"{a:.3f}",
            ha="center", va="bottom", fontsize=8)
ax.axhline(0.973, color="grey", linestyle="--", linewidth=0.8)
ax.text(len(algos_s) - 0.5, 0.976, "StepGLM minimum 0.973", ha="right",
        fontsize=8, color="grey")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "S1a_algorithm_benchmark.pdf"))
plt.close(fig)

# ---------------- S2: per-gene direction replication (v7) ----------------
val_path = os.path.join(
    ROOT, "results", "08_External_Validation", "v7_single_cohort",
    "validation_direction_v7.csv",
)
rows = read_csv_rows(val_path)
header, data = rows[0], rows[1:]
cohorts = ["GSE128543_WT_only", "GSE205958"]
labels = {
    "GSE128543_WT_only": "GSE128543 (CCI 48 h, WT-only, n=14)",
    "GSE205958": "GSE205958 (CCI 7 d, RNA-seq, n=10)",
}
gene_order = sorted({r[1] for r in data})
gene_num = {g: i + 1 for i, g in enumerate(gene_order)}
fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.2), sharex=True, sharey=True)
for ax, cohort in zip(axes, cohorts):
    sub = [r for r in data if r[0] == cohort]
    x = np.array([float(r[5]) for r in sub])   # discovery log2FC
    y = np.array([float(r[2]) for r in sub])   # validation mean diff
    ok = np.array([r[4].strip().upper() == "TRUE" for r in sub])
    genes = [r[1] for r in sub]
    rho = spearman(x, y)
    ax.axhline(0, color="grey", linewidth=0.8)
    ax.axvline(0, color="grey", linewidth=0.8)
    ax.scatter(x[ok], y[ok], s=55, color="#27AE60", edgecolor="black",
               linewidth=0.5, zorder=3, label="Direction consistent")
    ax.scatter(x[~ok], y[~ok], s=70, color="#C0392B", edgecolor="black",
               linewidth=0.5, zorder=4, label="Opposite direction")
    for xi, yi, g in zip(x, y, genes):
        ax.annotate(str(gene_num[g]), (xi, yi), textcoords="offset points",
                    xytext=(0, -13), fontsize=7.5, ha="center",
                    fontweight="bold")
    lim = max(np.max(np.abs(x)), np.max(np.abs(y)), 1.0) + 0.4
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Discovery log2FC (GSE58485-only)")
    ax.set_ylabel("Validation mean difference (TBI - sham)")
    ax.set_title(f"{labels[cohort]}\nSpearman $\\rho$={rho:.3f}")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.25)
legend_handles = [
    plt.Line2D([0], [0], marker="o", linestyle="", markersize=0,
               label=f"{n}: {g}")
    for g, n in sorted(gene_num.items(), key=lambda kv: kv[1])
]
fig.legend(handles=legend_handles, loc="lower center", ncol=6,
           fontsize=7.5, frameon=False, columnspacing=1.2)
fig.suptitle("Per-gene direction replication of the 11 candidate genes",
             fontsize=12)
fig.tight_layout(rect=(0, 0.08, 1, 0.94))
fig.savefig(os.path.join(OUT, "S2_direction_replication.pdf"))
plt.close(fig)

print("wrote S1a_algorithm_benchmark.pdf and S2_direction_replication.pdf")
