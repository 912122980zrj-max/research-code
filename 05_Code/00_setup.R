# ============================================================
# 00_setup.R —— 安装并加载本流水线所需全部 R 包
# 运行：Rscript scripts/00_setup.R   （R >= 4.3）
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# --- Bioconductor 包 ---
bioc_pkgs <- c(
  "GEOquery",        # GEO 数据下载
  "limma",           # 差异表达
  "sva",             # ComBat 批次校正
  "WGCNA",           # 加权共表达网络
  "clusterProfiler", # GO/KEGG 富集
  "org.Mm.eg.db",    # 小鼠注释
  "org.Hs.eg.db",    # 人类注释
  "Seurat",          # 单细胞
  "harmony",         # 单细胞批次整合
  "biomaRt",         # 小鼠->人同源映射（对接结构准备）
  "preprocessCore",  # CIBERSORT 分位数标准化（Windows 上如安装失败可跳过）
  "ComplexHeatmap"   # 热图（可选）
)
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing Bioc package: ", pkg)
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

# --- CRAN 包 ---
cran_pkgs <- c(
  "dplyr", "tidyr", "readr", "stringr", "data.table",
  "ggplot2", "ggrepel", "RColorBrewer", "pheatmap",
  "caret", "glmnet", "e1071", "randomForest", "xgboost",
  "kernlab", "mboost", "MASS", "pROC", "doParallel",
  "shapviz", "jsonlite", "rvest", "VennDiagram", "reshape2", "httr",
  "dynamicTreeCut", "cluster", "ggvenn", "remotes"
)
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing CRAN package: ", pkg)
    install.packages(pkg)
  }
}

# --- 单细胞通讯与虚拟敲除（GitHub / Bioconductor）---
if (!requireNamespace("CellChat", quietly = TRUE)) {
  message("Installing CellChat from GitHub (耗时较长)...")
  remotes::install_github("jinworks/CellChat", dependencies = TRUE, upgrade = "never")
}
if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
  message("Installing scTenifoldKnk from GitHub...")
  remotes::install_github("cailab-tamu/scTenifoldKnk", upgrade = "never")
}

# --- 输出已可用包清单 ---
pkgs <- c(bioc_pkgs, cran_pkgs, "CellChat", "scTenifoldKnk")
ip <- rownames(installed.packages())
message("=== 安装检查 ===")
for (pkg in pkgs) {
  cat(sprintf("%-20s %s\n", pkg, ifelse(pkg %in% ip, "OK", "MISSING")))
}
cat("setup finished.\n")
