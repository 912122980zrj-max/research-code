# 00_check_installed.R
# 用途：检测 TBI 项目所需的 R 包安装状态，输出“已装/缺失 + 版本”清单。
# 用法：Rscript 00_check_installed.R

options(timeout = 600)

cat("R:", R.version.string, "\n")
cat("libPaths:", paste(.libPaths(), collapse = " | "), "\n\n")

check_pkg <- function(p) {
  ok <- requireNamespace(p, quietly = TRUE)
  ver <- if (ok) tryCatch(as.character(packageVersion(p)), error = function(e) NA)
  cat(sprintf("%-20s %-8s %s\n", p, ifelse(ok, "OK", "MISSING"), ifelse(ok, ver, "")))
}

targets <- c(
  # 基础设施
  "BiocManager", "remotes",
  # 已安装核对（清单 §0）
  "Seurat", "WGCNA", "limma", "sva", "clusterProfiler", "org.Mm.eg.db",
  "caret", "glmnet", "pROC", "shapviz", "ggvenn", "harmony", "presto",
  "CellChat", "scTenifoldKnk", "scTenifoldNet", "BiocNeighbors",
  "assorthead", "enrichR", "pheatmap", "patchwork",
  # D3（注意：Mime 仓库已迁移，包名实际为 Mime1）
  "Mime1", "rms", "rmda",
  # D5
  "MuSiC", "BisqueRNA",
  # D7
  "DoubletFinder", "SoupX"
)

for (p in targets) check_pkg(p)
