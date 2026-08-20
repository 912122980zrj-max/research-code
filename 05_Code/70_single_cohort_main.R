# ============================================================
# 70_single_cohort_main_v7.R —— v7 主分析：GSE58485-only
# 依据：2026-08-11 评估意见（意见一/意见二）：
#   “合并队列不可作为主分析；应以 GSE58485-only 为主发现队列，
#    合并分析降为异质性压力测试。”
# 数据口径与 63_sensitivity_matrix.R 的 GSE58485-only 路径完全一致：
#   CLEAN 合并矩阵（5,043 基因，公共基因合并 + 每队列覆盖过滤 + 中位数填补）
#   → 取 GSE58485 样本（n=24）→ limma ~group+genotype → WGCNA → 候选交集。
# 输出：results/16_Orthology/v7_single_cohort/（表 + 图）
# ============================================================
suppressMessages({
  library(limma)
  library(WGCNA)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggvenn)
})

ROOT <- getwd()
OUT <- file.path(ROOT, "results", "16_Orthology", "v7_single_cohort")
FIG <- file.path(ROOT, "figures", "Figure_03_Discovery")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts", "theme_pub.R"))

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(OUT, "v7_run.log"), append = TRUE)
}
if (file.exists(file.path(OUT, "v7_run.log"))) file.remove(file.path(OUT, "v7_run.log"))

GENO_CODE <- "wildtype"
COV_THRESHOLD <- 0.10

read_bulk <- function(gse) readRDS(file.path(ROOT, "data", "01_GEO_bulk",
                                             paste0(gse, ".rds")))

comp <- read.csv(file.path(ROOT, "results", "00_Data_Prep",
                           "sample_composition.csv"), check.names = FALSE)
clean_samples <- comp$sample[comp$include_clean]

# ---------- 与 03/63 完全一致的基础预处理 ----------
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
  keep_cov <- cov_a >= COV_THRESHOLD & cov_b >= COV_THRESHOLD
  expr_all <- expr_all[keep_cov, , drop = FALSE]
  expr_all <- apply(expr_all, 2, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE); x
  })
  list(expr = expr_all, cohort = cohort_vec,
       group = ifelse(grp_map == "tbi", "TBI", "Sham"))
}

# ---------- WGCNA（与 03/63 完全一致） ----------
run_wgcna <- function(expr_corrected, ph, label) {
  mad <- apply(expr_corrected, 1, mad)
  top_genes <- names(sort(mad, decreasing = TRUE))[seq_len(min(8000, length(mad)))]
  gene_var <- apply(expr_corrected[top_genes, ], 1, var, na.rm = TRUE)
  top_genes <- top_genes[gene_var > 0]
  datExpr <- t(expr_corrected[top_genes, ])
  sft <- pickSoftThreshold(datExpr, powerVector = 1:20, verbose = 0)
  softPower <- sft$powerEstimate
  if (is.na(softPower)) softPower <- 12
  adjacency <- adjacency(datExpr, power = softPower)
  TOM <- TOMsimilarity(adjacency, verbose = 0)
  dissTOM <- 1 - TOM
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                               deepSplit = 2, pamRespectsDendro = FALSE)
  dynamicColors <- labels2colors(dynamicMods)
  names(dynamicColors) <- colnames(datExpr)
  MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
  MEs <- MEList$eigengenes
  group_num <- as.numeric(ph$group == "TBI")
  mod_trait <- cor(MEs, group_num, use = "p")
  mod_p <- corPvalueStudent(mod_trait, nSamples = nrow(MEs))
  sig_mods <- rownames(mod_p)[mod_p[, 1] < 0.05 & !grepl("^MEgrey$", rownames(mod_p))]
  sig_colors <- sub("^ME", "", sig_mods)
  module_genes <- unique(colnames(datExpr)[dynamicColors %in% sig_colors])
  best <- rownames(mod_trait)[which.max(abs(mod_trait[, 1]))]
  best_genes <- unique(colnames(datExpr)[dynamicColors == sub("^ME", "", best)])
  list(module_genes = module_genes, best_module = sub("^ME", "", best),
       best_r = unname(mod_trait[best, 1]), best_p = unname(mod_p[best, 1]),
       softPower = softPower, colors = dynamicColors, MEs = MEs,
       mod_trait = mod_trait, mod_p = mod_p, best_genes = best_genes,
       geneTree = geneTree, datExpr = datExpr)
}

# ---------- DEG ----------
run_deg <- function(expr, ph, design = NULL) {
  if (is.null(design)) design <- model.matrix(~ group, data = ph)
  fit <- lmFit(expr, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = "groupTBI", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  deg <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.585, ]
  list(all = tt, deg = deg$gene)
}

# ---------- 靶点分级与可评估集 ----------
tier_tab <- read_csv(file.path(ROOT, "results", "16_Orthology",
                               "bap_targets_364_tier.csv"),
                     col_types = cols(.default = col_character()))
tier_map <- setNames(tier_tab$tier, toupper(tier_tab$gene_symbol))
eval_clean_orth <- read_csv(file.path(ROOT, "results", "16_Orthology",
                                      "evaluable_targets_CLEAN.csv"),
                            col_types = cols(.default = col_character()))
eval_clean_upper <- toupper(unique(eval_clean_orth$matrix_symbol))

candidates_tiered <- function(key_genes) {
  cand <- key_genes[toupper(key_genes) %in% eval_clean_upper]
  data.frame(gene = cand,
             tier = ifelse(is.na(tier_map[toupper(cand)]), "C",
                           unname(tier_map[toupper(cand)])),
             stringsAsFactors = FALSE)
}

# ============================================================
# 主分析：GSE58485-only
# ============================================================
log_msg("preparing GSE58485-only matrix (n=24)")
pp <- prep()
keep58 <- pp$cohort == "GSE58485"
ex58 <- pp$expr[, keep58, drop = FALSE]
ph58 <- data.frame(group = factor(pp$group[keep58], levels = c("Sham", "TBI")))
d58 <- read_bulk("GSE58485")
gt58 <- as.character(d58$pheno$`genotype:ch1`)[match(colnames(ex58),
                                                     rownames(d58$pheno))]
gt58[is.na(gt58) | gt58 == ""] <- GENO_CODE
ph58$gt <- factor(gt58)
design58 <- model.matrix(~ group + gt, data = ph58)

log_msg("running limma (design ~group+genotype)")
deg58 <- run_deg(ex58, ph58, design58)
n_deg <- length(deg58$deg)
log_msg("DEGs =", n_deg)

log_msg("running WGCNA (single cohort)")
wgcna58 <- run_wgcna(ex58, ph58, "GSE58485_only")
log_msg("best module =", wgcna58$best_module,
        " r =", round(wgcna58$best_r, 4),
        " p =", format(wgcna58$best_p, digits = 3))

key58 <- intersect(deg58$deg, wgcna58$module_genes)
cand58 <- candidates_tiered(key58)
n_key <- length(key58)
n_cand <- nrow(cand58)
log_msg("key genes =", n_key, " candidates =", n_cand)
log_msg(paste(cand58$gene, collapse = ","))

key58_best <- intersect(deg58$deg, wgcna58$best_genes)
cand58_best <- candidates_tiered(key58_best)
log_msg("best-module-only key genes =", length(key58_best),
        " candidates =", nrow(cand58_best))

# ---------- 保存表 ----------
write_csv(deg58$all, file.path(OUT, "deg_all_GSE58485_only.csv"))
write_csv(deg58$all[deg58$all$gene %in% deg58$deg, ],
          file.path(OUT, "deg_significant_GSE58485_only.csv"))
write_csv(data.frame(gene = wgcna58$module_genes),
          file.path(OUT, "module_genes_GSE58485_only.csv"))
write_csv(data.frame(gene = wgcna58$best_genes),
          file.path(OUT, "best_module_genes_GSE58485_only.csv"))
write_csv(data.frame(gene = key58), file.path(OUT, "key_genes_GSE58485_only.csv"))
write_csv(data.frame(gene = key58_best),
          file.path(OUT, "key_genes_best_module_GSE58485_only.csv"))
write_csv(cand58, file.path(OUT, "candidates_GSE58485_only.csv"))
write_csv(cand58_best, file.path(OUT, "candidates_best_module_GSE58485_only.csv"))

mod_tab <- data.frame(module = rownames(wgcna58$mod_trait),
                      r = unname(wgcna58$mod_trait[, 1]),
                      p = unname(wgcna58$mod_p[, 1]),
                      significant = unname(wgcna58$mod_p[, 1] < 0.05 &
                                             !grepl("^MEgrey$", rownames(wgcna58$mod_trait))))
write_csv(mod_tab, file.path(OUT, "module_trait_GSE58485_only.csv"))

# 候选统计表（logFC / P / adjP）
stat_rows <- deg58$all[match(cand58$gene, deg58$all$gene), ]
cand_stat <- data.frame(
  gene = cand58$gene,
  tier = cand58$tier,
  logFC = round(stat_rows$logFC, 3),
  P = format(stat_rows$P.Value, digits = 3, scientific = TRUE),
  adjP = format(stat_rows$adj.P.Val, digits = 3, scientific = TRUE),
  stringsAsFactors = FALSE)
write_csv(cand_stat, file.path(OUT, "candidate_stats_GSE58485_only.csv"))
print(cand_stat)

# ---------- 图：火山图 ----------
vol <- deg58$all
vol$sig <- vol$gene %in% deg58$deg
vol$cand <- vol$gene %in% cand58$gene
p_vol <- ggplot(vol, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = cand), size = 0.6, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#E6550D")) +
  theme_pub() + labs(x = "log2 fold change (TBI vs sham)",
                     y = expression(-log[10] ~ P),
                     title = "GSE58485-only DEGs (n=373; 11 candidates highlighted)") +
  theme(legend.position = "none")
ggsave(file.path(FIG, "volcano_deg_GSE58485_only.pdf"), p_vol,
       width = 6, height = 5)

# ---------- 图：模块-性状 ----------
mdf <- mod_tab[!grepl("^grey$", mod_tab$module), ]
mdf$label <- ifelse(mdf$significant, paste0("r=", round(mdf$r, 3),
                                            ", P=", format(mdf$p, digits = 2)), "")
p_mod <- ggplot(mdf, aes(reorder(module, r), r, fill = significant)) +
  geom_col() +
  geom_text(aes(label = label), hjust = -0.05, size = 2.8) +
  coord_flip() + scale_fill_manual(values = c("FALSE" = "grey75", "TRUE" = "#3182BD")) +
  theme_pub() + labs(x = NULL, y = "Module–trait correlation r (TBI)",
                     title = "GSE58485-only WGCNA modules") +
  theme(legend.position = "none")
ggsave(file.path(FIG, "module_trait_GSE58485_only.pdf"), p_mod,
       width = 7, height = 5)

# ---------- 图：候选 logFC ----------
cs <- cand_stat
cs$logFC <- as.numeric(cs$logFC)
p_cand <- ggplot(cs, aes(reorder(gene, logFC), logFC)) +
  geom_col(fill = "#E6550D") +
  coord_flip() +
  theme_pub() + labs(x = NULL, y = "log2 fold change (TBI vs sham)",
                     title = "11 candidate genes (GSE58485-only main path)")
ggsave(file.path(FIG, "figure3E_candidate_logFC_GSE58485_only.pdf"), p_cand,
       width = 6, height = 4.5)

# ---------- 图：Venn ----------
venn_list <- list(
  `Evaluable BaP targets` = eval_clean_upper,
  `DEGs` = toupper(deg58$deg),
  `Module key genes` = toupper(key58))
p_venn <- ggvenn(venn_list,
                 fill_color = c("#9ECAE1", "#A1D99B", "#FDBE85"),
                 stroke_size = 0.4, set_name_size = 4.2, text_size = 3.5)
ggsave(file.path(ROOT, "figures", "Figure_02_BaP_Targets",
                 "venn_candidates_v7.pdf"), p_venn, width = 8, height = 6)

# ---------- 汇总 ----------
writeLines(c(
  "v7 single-cohort main analysis (GSE58485-only)",
  paste0("generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("samples: ", ncol(ex58)),
  paste0("genes after coverage filter: ", nrow(ex58)),
  paste0("DEGs (|log2FC|>0.585, BH<0.05): ", n_deg),
  paste0("WGCNA softPower: ", wgcna58$softPower),
  paste0("best module: ", wgcna58$best_module,
         " r=", round(wgcna58$best_r, 4),
         " P=", format(wgcna58$best_p, digits = 3)),
  paste0("module genes (significant modules union): ", length(wgcna58$module_genes)),
  paste0("key genes (DEG∩module): ", n_key),
  paste0("candidates (key∩69 evaluable): ", n_cand,
         " -> ", paste(cand58$gene, collapse = ",")),
  paste0("best-module-only key genes: ", length(key58_best),
         " candidates: ", paste(cand58_best$gene, collapse = ","))
), file.path(OUT, "v7_summary.txt"))
log_msg("done")
