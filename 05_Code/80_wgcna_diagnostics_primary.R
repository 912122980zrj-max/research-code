# ============================================================
# 80_wgcna_diagnostics_primary.R —— M8 应对方针：主分析 WGCNA 诊断
# 复现 70_single_cohort_main_v7.R 的 run_wgcna（top 8,000 MAD 上限、
# pickSoftThreshold、cutreeDynamic deepSplit=2），仅保存诊断输出：
#   results/16_Orthology/v7_single_cohort/wgcna_soft_threshold_GSE58485_only.csv
#   results/16_Orthology/v7_single_cohort/wgcna_module_sizes_GSE58485_only.csv
# 并核对 green 模块 r=0.978。
# ============================================================
suppressMessages({
  library(WGCNA)
  library(readr)
})

ROOT <- getwd()
OUT <- file.path(ROOT, "results", "16_Orthology", "v7_single_cohort")

comp <- read.csv(file.path(ROOT, "results", "00_Data_Prep",
                           "sample_composition.csv"), check.names = FALSE)
clean_samples <- comp$sample[comp$include_clean]          # n=33 (both cohorts)
clean24 <- comp$sample[comp$cohort == "GSE58485" & comp$include_clean]

# ---- 与 70_single_cohort_main_v7.R prep() 完全一致（未校正主矩阵）----
read_bulk <- function(gse) readRDS(file.path(ROOT, "data", "01_GEO_bulk",
                                             paste0(gse, ".rds")))
prep <- function() {
  gses <- c("GSE58485", "GSE41345")
  dat <- lapply(gses, read_bulk)
  names(dat) <- gses
  common <- Reduce(intersect, lapply(dat, function(d) rownames(d$expr)))
  expr_all <- do.call(cbind, lapply(gses, function(g) dat[[g]]$expr[common, , drop = FALSE]))
  expr_all <- expr_all[, intersect(clean_samples, colnames(expr_all)), drop = FALSE]
  cohort_vec <- comp$cohort[match(colnames(expr_all), comp$sample)]
  grp_map <- setNames(comp$group[match(colnames(expr_all), comp$sample)],
                      colnames(expr_all))
  keep_gene <- rowSums(is.na(expr_all)) < 0.2 * ncol(expr_all)
  expr_all <- expr_all[keep_gene, , drop = FALSE]
  cov_a <- rowMeans(expr_all[, cohort_vec == "GSE58485", drop = FALSE] > 0, na.rm = TRUE)
  cov_b <- rowMeans(expr_all[, cohort_vec == "GSE41345", drop = FALSE] > 0, na.rm = TRUE)
  keep_cov <- cov_a >= 0.10 & cov_b >= 0.10
  expr_all <- expr_all[keep_cov, , drop = FALSE]
  expr_all <- apply(expr_all, 2, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE); x
  })
  list(expr = expr_all, cohort = cohort_vec,
       group = ifelse(grp_map == "tbi", "TBI", "Sham"))
}
pp <- prep()
ex <- pp$expr[, pp$cohort == "GSE58485", drop = FALSE]
stopifnot(nrow(ex) == 5043, ncol(ex) == 24)

# 保存真实主分析矩阵（未校正；供源数据与后续分析使用）
write.csv(ex, file.path(OUT, "expression_GSE58485_only_primary.csv"))
cat("Hprt range:",
    paste(range(ex["Hprt", ]), collapse = " - "), "\n")

mad <- apply(ex, 1, mad)
top_genes <- names(sort(mad, decreasing = TRUE))[seq_len(min(8000, length(mad)))]
gene_var <- apply(ex[top_genes, ], 1, var, na.rm = TRUE)
top_genes <- top_genes[gene_var > 0]
datExpr <- t(ex[top_genes, ])
cat("WGCNA input genes:", ncol(datExpr), "\n")

sft <- pickSoftThreshold(datExpr, powerVector = 1:20, verbose = 0)
softPower <- sft$powerEstimate
if (is.na(softPower)) softPower <- 12
sft_tab <- data.frame(power = sft$fitIndices$Power,
                      SFT.R.sq = signif(sft$fitIndices$SFT.R.sq, 4),
                      slope = signif(sft$fitIndices$slope, 3),
                      mean.k = signif(sft$fitIndices$mean.k., 3),
                      median.k = signif(sft$fitIndices$median.k., 3))
write_csv(sft_tab, file.path(OUT, "wgcna_soft_threshold_GSE58485_only.csv"))
cat("softPower =", softPower, "\n")

adj <- adjacency(datExpr, power = softPower)
TOM <- TOMsimilarity(adj, verbose = 0)
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE)
dynamicColors <- labels2colors(dynamicMods)
names(dynamicColors) <- colnames(datExpr)

mod_size <- data.frame(module = names(table(dynamicColors)),
                       n_genes = as.integer(table(dynamicColors)))
mod_size <- mod_size[order(-mod_size$n_genes), ]
write_csv(mod_size, file.path(OUT, "wgcna_module_sizes_GSE58485_only.csv"))
print(head(mod_size, 40), row.names = FALSE)

MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
group_num <- as.numeric(ifelse(comp$group[match(clean24, comp$sample)] == "tbi",
                               "TBI", "Sham") == "TBI")
mod_trait <- cor(MEs, group_num, use = "p")
mod_p <- corPvalueStudent(mod_trait, nSamples = nrow(MEs))
sig_mods <- rownames(mod_p)[mod_p[, 1] < 0.05 & !grepl("^MEgrey$", rownames(mod_p))]
module_genes <- unique(colnames(datExpr)[dynamicColors %in% sub("^ME", "", sig_mods)])
mod_tab <- data.frame(module = rownames(mod_trait),
                      r = unname(mod_trait[, 1]),
                      p = unname(mod_p[, 1]),
                      significant = unname(mod_p[, 1] < 0.05 &
                                             !grepl("^MEgrey$", rownames(mod_trait))))
write_csv(mod_tab, file.path(OUT, "module_trait_GSE58485_only_diagnostics.csv"))
cat("significant modules:", paste(sub("^ME", "", sig_mods), collapse = ","), "\n")
cat("union module genes:", length(module_genes), "\n")
best <- rownames(mod_trait)[which.max(abs(mod_trait[, 1]))]
cat("best module:", sub("^ME", "", best),
    " r =", round(mod_trait[best, 1], 4),
    " P =", format(mod_p[best, 1], digits = 3), "\n")
