# ============================================================
# 17_deconvolution_analysis.R —— 去卷积结果统计、验证与发表级绘图
# 输入：data/05_CIBERSORTx/deconvolution_*.csv（MuSiC/BisqueRNA）
#       data/05_CIBERSORTx/CIBERSORTx_Job7_Results.csv（网页版，S-mode 未生效）
#       data/05_CIBERSORTx/pseudobulk_loo_validation.csv
#       data/01_GEO_bulk/GSE104687_pheno.csv
# 输出：results/09_Immune_Deconvolution/
# ============================================================
suppressMessages({
  library(data.table); library(ggplot2)
})
options(warn = 1)

OUT <- "results/09_Immune_Deconvolution"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source("scripts/theme_pub.R")

ph <- fread("data/01_GEO_bulk/GSE104687_pheno.csv")
ph <- data.frame(sample = as.character(ph$sample),
                 group = as.character(ph$group), stringsAsFactors = FALSE)

read_est <- function(f) {
  m <- as.matrix(fread(f), rownames = 1)
  m
}

music8 <- read_est("data/05_CIBERSORTx/deconvolution_music_8ct.csv")
bisq8  <- read_est("data/05_CIBERSORTx/deconvolution_bisque_8ct.csv")
music6 <- read_est("data/05_CIBERSORTx/deconvolution_music_6ct.csv")
bisq6  <- read_est("data/05_CIBERSORTx/deconvolution_bisque_6ct.csv")
ct_names <- c("Neuron", "Oligodendrocyte", "OPC", "Microglia",
              "Astrocyte", "Endothelial", "T_NK", "Macrophage")
cib7_raw <- read_est("data/05_CIBERSORTx/CIBERSORTx_Job7_Results.csv")
cib7 <- cib7_raw[, intersect(ct_names, colnames(cib7_raw)), drop = FALSE]
cat("MuSiC8:", dim(music8), "| Bisque8:", dim(bisq8),
    "| MuSiC6:", dim(music6), "| Bisque6:", dim(bisq6),
    "| CIB7:", dim(cib7), "\n")

stopifnot(all(rownames(music8) %in% ph$sample))

methods <- list(
  MuSiC_8ct = music8, BisqueRNA_8ct = bisq8,
  MuSiC_6ct = music6, BisqueRNA_6ct = bisq6,
  CIBERSORTx_Job7 = cib7
)

# ---------- 1) 组间比较（Wilcoxon + BH） ----------
stat_rows <- list()
for (mname in names(methods)) {
  m <- methods[[mname]]
  for (ct in colnames(m)) {
    g <- ph$group[match(rownames(m), ph$sample)]
    a <- m[g == "Sham", ct]; b <- m[g == "TBI", ct]
    wt <- suppressWarnings(wilcox.test(a, b))
    stat_rows[[length(stat_rows) + 1]] <- data.frame(
      method = mname, celltype = ct,
      n_sham = sum(g == "Sham"), n_tbi = sum(g == "TBI"),
      mean_sham = mean(a), sd_sham = sd(a), median_sham = median(a),
      mean_tbi = mean(b), sd_tbi = sd(b), median_tbi = median(b),
      p_wilcox = wt$p.value, stringsAsFactors = FALSE)
  }
}
stats <- rbindlist(stat_rows)
stats[, fdr := p.adjust(p_wilcox, method = "BH"), by = method]
stats[, note := ""]
stats[method == "CIBERSORTx_Job7", note := "P值来自perm=0任务，不可用"]
fwrite(stats, file.path(OUT, "deconvolution_group_stats.csv"))
cat("\n== 组间比较（p<0.05 的条目）==\n")
print(stats[p_wilcox < 0.05, .(method, celltype, mean_sham, mean_tbi,
                               p_wilcox, fdr)])

# ---------- 2) 方法间一致性（Pearson） ----------
cor_rows <- list()
for (a_name in c("MuSiC_8ct", "BisqueRNA_8ct")) {
  for (b_name in c("CIBERSORTx_Job7", "MuSiC_8ct", "BisqueRNA_8ct")) {
    if (a_name == b_name) next
    A <- methods[[a_name]]; B <- methods[[b_name]]
    shared_ct <- intersect(colnames(A), colnames(B))
    for (ct in shared_ct) {
      x <- A[, ct]; y <- B[, ct]
      r <- cor(x, y, method = "pearson")
      cor_rows[[length(cor_rows) + 1]] <- data.frame(
        method_a = a_name, method_b = b_name, celltype = ct,
        pearson_r = r, n = length(x), stringsAsFactors = FALSE)
    }
  }
}
corr <- rbindlist(cor_rows)
corr <- corr[(method_a %in% c("MuSiC_8ct", "BisqueRNA_8ct") &
                method_b == "CIBERSORTx_Job7") |
               (method_a == "MuSiC_8ct" & method_b == "BisqueRNA_8ct")]
fwrite(corr, file.path(OUT, "deconvolution_method_corr.csv"))
cat("\n== 方法一致性 ==\n")
print(corr)

# ---------- 3) LOO 伪 bulk 验证 ----------
loo <- fread("data/05_CIBERSORTx/pseudobulk_loo_validation.csv")
loo_types <- intersect(sub("^true_", "", colnames(loo)[grepl("^true_", colnames(loo))]),
                       sub("^est_", "", colnames(loo)[grepl("^est_", colnames(loo))]))
val_rows <- list()
for (ct in loo_types) {
  tv <- loo[[paste0("true_", ct)]]
  ev <- loo[[paste0("est_", ct)]]
  ok <- !is.na(ev)
  r <- if (sum(ok) >= 3) cor(tv[ok], ev[ok]) else NA
  rmse <- sqrt(mean((tv[ok] - ev[ok])^2))
  val_rows[[length(val_rows) + 1]] <- data.frame(
    celltype = ct, n = sum(ok), pearson_r = r, rmse = rmse,
    stringsAsFactors = FALSE)
}
val <- rbindlist(val_rows)
fwrite(val, file.path(OUT, "pseudobulk_loo_stats.csv"))
cat("\n== LOO 验证 ==\n")
print(val)

# ---------- 4) 汇总表（组均值，供报告引用） ----------
sum_rows <- list()
for (mname in names(methods)) {
  m <- methods[[mname]]
  for (ct in colnames(m)) {
    g <- ph$group[match(rownames(m), ph$sample)]
    sum_rows[[length(sum_rows) + 1]] <- data.frame(
      method = mname, celltype = ct,
      mean_sham = mean(m[g == "Sham", ct]),
      mean_tbi = mean(m[g == "TBI", ct]),
      mean_all = mean(m[, ct]), stringsAsFactors = FALSE)
  }
}
fwrite(rbindlist(sum_rows), file.path(OUT, "deconvolution_summary_means.csv"))

# ---------- 5) 发表级图形 ----------
okabe_ito <- c("#0072B2", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#D55E00", "#CC79A7", "#999999")
names(okabe_ito) <- c("Neuron", "Oligodendrocyte", "OPC", "Microglia",
                       "Astrocyte", "Endothelial", "T_NK", "Macrophage")

make_long <- function(m, mname) {
  dt <- as.data.table(m, keep.rownames = "sample")
  melt(dt, id.vars = "sample", variable.name = "celltype",
       value.name = "proportion")
}

box_by_group <- function(m, mname, filebase, method_key, types_keep = NULL) {
  dt <- make_long(m, mname)
  dt$group <- ph$group[match(dt$sample, ph$sample)]
  if (!is.null(types_keep)) dt <- dt[celltype %in% types_keep]
  dt$celltype <- factor(dt$celltype, levels = names(okabe_ito)[names(okabe_ito) %in% unique(dt$celltype)])
  stats_dt <- read.csv(file.path(OUT, "deconvolution_group_stats.csv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
  stats_dt <- stats_dt[stats_dt$method == method_key, c("celltype", "fdr")]
  stats_dt$lab <- sprintf("FDR P = %.2f%s", stats_dt$fdr,
                          ifelse(stats_dt$fdr < 0.05, " *", ""))
  p <- ggplot(dt, aes(x = group, y = proportion, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.85) +
    geom_jitter(width = 0.18, size = 0.35, alpha = 0.25, color = "grey20") +
    facet_wrap(~ celltype, nrow = 2, scales = "free_y") +
    geom_text(data = stats_dt, aes(x = Inf, y = Inf, label = lab),
              hjust = 1.1, vjust = 1.4, size = 2.8, inherit.aes = FALSE) +
    scale_fill_manual(values = c(Sham = "#EFEFEF", TBI = "#6BAED6")) +
    labs(x = NULL, y = "Estimated proportion",
         title = mname) +
    theme_bw(base_size = 10) +
    theme(legend.position = "top",
          strip.background = element_rect(fill = "grey92", colour = NA),
          panel.grid.minor = element_blank())
  ggsave_pub(file.path(OUT, paste0(filebase, ".pdf")), p,
         width = 7.2, height = 4.8)
  p
}

cat("\n== 绘图 ==\n")
box_by_group(music8, "MuSiC (8-class input; Macrophage excluded by algorithm)",
             "fig1_music_8ct_by_group", method_key = "MuSiC_8ct")
box_by_group(bisq8, "BisqueRNA (8 classes)",
             "fig2_bisque_8ct_by_group", method_key = "BisqueRNA_8ct")
box_by_group(cib7, "CIBERSORTx web Job7 (uncorrected; S-mode not applied)",
             "fig3_cibersortx_job7_by_group", method_key = "CIBERSORTx_Job7")

# 方法对比散点
comp_dt <- make_long(music8, "MuSiC")
idx <- cbind(match(comp_dt$sample, rownames(cib7)),
             match(as.character(comp_dt$celltype), colnames(cib7)))
comp_dt$cibersortx <- cib7[idx]
comp_dt <- comp_dt[!is.na(cibersortx)]
comp_dt$celltype <- factor(comp_dt$celltype, levels = names(okabe_ito)[names(okabe_ito) %in% unique(comp_dt$celltype)])
p4 <- ggplot(comp_dt, aes(x = cibersortx, y = proportion, color = celltype)) +
  geom_point(size = 0.4, alpha = 0.35) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  facet_wrap(~ celltype, nrow = 2, scales = "free") +
  scale_color_manual(values = okabe_ito, guide = "none") +
  labs(x = "CIBERSORTx web Job7", y = "MuSiC",
       title = "MuSiC vs CIBERSORTx (7 shared types)") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        panel.grid.minor = element_blank())
ggsave_pub(file.path(OUT, "fig4_music_vs_cibersortx.pdf"), p4,
       width = 7.2, height = 4.8)

# LOO 验证散点
loo_dt <- melt(loo, id.vars = "sample",
               measure.vars = list(
                 true = grep("^true_", colnames(loo)),
                 est = grep("^est_", colnames(loo))),
               variable.name = "celltype")
loo_dt$celltype <- sub("^true_", "", loo_dt$celltype)
loo_dt <- loo_dt[!is.na(est)]
loo_dt$celltype <- factor(loo_dt$celltype, levels = names(okabe_ito)[names(okabe_ito) %in% unique(loo_dt$celltype)])
p5 <- ggplot(loo_dt, aes(x = true, y = est, color = celltype)) +
  geom_point(size = 1.6, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  facet_wrap(~ celltype, nrow = 2, scales = "free") +
  scale_color_manual(values = okabe_ito, guide = "none") +
  labs(x = "True proportion (pseudobulk)", y = "MuSiC estimate (LOO)",
       title = "Leave-one-sample-out pseudobulk validation") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        panel.grid.minor = element_blank())
ggsave_pub(file.path(OUT, "fig5_loo_validation.pdf"), p5,
       width = 7.2, height = 4.8)

cat("\n== 完成：", OUT, "==\n")
