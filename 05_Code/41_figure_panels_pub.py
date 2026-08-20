#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""补生成两个缺失面板（Figure 3D 模块-性状热图、Figure 5A 算法 AUC），
统一 Arial + 9/8/7 字号层级，导出可编辑 PDF（2026-08-09 部署规范）。"""

import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# 统一字体：强制 Arial
plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.labelsize"] = 8
plt.rcParams["xtick.labelsize"] = 7
plt.rcParams["ytick.labelsize"] = 7
plt.rcParams["legend.fontsize"] = 7
plt.rcParams["axes.titlesize"] = 9

ROOT = r"E:\sheng xin\tbi"


def read_csv(path):
    with open(path, encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


# ---------- Figure 3D：模块-性状相关热图 ----------
rows = read_csv(os.path.join(ROOT, "results", "04_WGCNA", "module_trait.csv"))
mods = [r["module"].replace("ME", "") for r in rows]
rs = [float(r["r"]) for r in rows]
ps = [float(r["P"]) for r in rows]
n_mods = len(mods)

fig, ax = plt.subplots(figsize=(3.2, 2.4 + 0.045 * n_mods))
im = ax.imshow([[r] for r in rs], cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
txt_size = max(5.0, min(8.0, 180.0 / n_mods))
for i, (r, p) in enumerate(zip(rs, ps)):
    star = "*" if p < 0.05 else ""
    ax.text(0, i, f"{r:.3f}{star}", ha="center", va="center",
            fontsize=txt_size, color="black")
ax.set_xticks([0])
ax.set_xticklabels(["TBI"], fontsize=8)
ax.set_yticks(range(len(mods)))
ax.set_yticklabels(mods, fontsize=6.5)
cbar = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
cbar.set_label("Pearson r", fontsize=8)
cbar.ax.tick_params(labelsize=7)
fig.tight_layout()
out3 = os.path.join(ROOT, "figures", "Figure_03_Discovery", "module_trait_heatmap.pdf")
fig.savefig(out3)
print("saved", out3)

# ---------- Figure 5A：11 算法 CV AUC（raw CLEAN，描述性） ----------
rows5 = read_csv(os.path.join(ROOT, "results", "07_ML_Signature",
                              "D3_ml_benchmark", "model_comparison.csv"))
raw = [r for r in rows5 if r["matrix"] == "raw_CLEAN"]
raw.sort(key=lambda r: float(r["cv_auc_mean"]), reverse=True)
algos = [r["algorithm"] for r in raw]
aucs = [float(r["cv_auc_mean"]) for r in raw]
sds = [float(r["cv_auc_sd"]) for r in raw]

fig, ax = plt.subplots(figsize=(5.2, 2.8))
ax.barh(algos, aucs, color="#4C72B0", edgecolor="black", linewidth=0.4,
        xerr=sds, error_kw=dict(elinewidth=0.5, capsize=1.5, ecolor="grey"))
for y, v in enumerate(aucs):
    ax.text(v + 0.008, y, f"{v:.3f}", va="center", fontsize=7)
ax.set_xlabel("CV AUC (raw CLEAN, descriptive)", fontsize=8)
ax.set_xlim(0.8, 1.08)
ax.tick_params(axis="both", labelsize=7)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
out5 = os.path.join(ROOT, "figures", "Figure_05_ML", "Figure5A_algorithm_auc.pdf")
fig.savefig(out5)
print("saved", out5)
