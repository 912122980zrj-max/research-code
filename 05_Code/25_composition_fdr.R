# ============================================================
# 25_composition_fdr.R —— P0-1：单细胞组成检验补 BH-FDR
# 依据：FINAL_REVIEW_PUBLISHABILITY.md §4/P0-1（H5）
# 作用：为 06b 与 D7 主版本的样本级 Wilcoxon 表追加 BH 校正列；
#       不重算检验，仅由既有 p 值派生（8 类 × 1 次校正）。
# 输入：results/11_SingleCell/celltype_composition_stats.csv
#       results/11_SingleCell/D7_doublet/composition_stats_d7.csv
# 输出：同文件原位追加 wilcox_p_bh 列
# ============================================================

add_bh <- function(path) {
  d <- read.csv(path, check.names = FALSE, fileEncoding = "UTF-8")
  if (!"wilcox_p" %in% colnames(d)) stop(paste("no wilcox_p in", path))
  d$wilcox_p_bh <- round(p.adjust(d$wilcox_p, method = "BH"), 5)
  write.csv(d, path, row.names = FALSE, fileEncoding = "UTF-8")
  d
}

a <- add_bh("results/11_SingleCell/celltype_composition_stats.csv")
cat("06b stats:\n")
print(a[, c("celltype", "wilcox_p", "wilcox_p_bh")])

b <- add_bh("results/11_SingleCell/D7_doublet/composition_stats_d7.csv")
cat("D7 stats:\n")
print(b[, c("celltype", "wilcox_p", "wilcox_p_bh")])

cat("BH-FDR columns added.\n")
