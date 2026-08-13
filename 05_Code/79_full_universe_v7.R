# ============================================================
# 79_full_universe_v7.R —— M1 应对方针：完整 GPL1261 宇宙重跑
# 依据：BMC 投稿前评审 M1（主分析基因宇宙 = 两平台公共基因∩覆盖过滤，
#       Ahr/Cyp1a1/Cyp1a2/Cyp1b1/Ahrr/Maoa 在 GPL1261 上可测但被双平台
#       覆盖过滤剔除）。
# 本脚本在“GSE58485 单队列完整基因矩阵（gene-level，27,250 基因）+
#       24 个 CLEAN 样本 + 单队列 ≥10% 覆盖过滤”上重跑主分析，
#       并输出 WGCNA 诊断（软阈值 R²、模块大小）与均值-方差分层置换。
# 输出：results/16_Orthology/v7_full_universe/
# ============================================================
suppressMessages({
  library(limma)
  library(WGCNA)
  library(dplyr)
  library(readr)
})

ROOT <- getwd()
OUT <- file.path(ROOT, "results", "16_Orthology", "v7_full_universe")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

COV_THRESHOLD <- 0.10
SEED <- 2026
N_PERM <- 10000

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(OUT, "run.log"), append = TRUE)
}
if (file.exists(file.path(OUT, "run.log"))) file.remove(file.path(OUT, "run.log"))

# ---------- 1. 完整 GSE58485 基因级矩阵 + 24 CLEAN 样本 ----------
comp <- read.csv(file.path(ROOT, "results", "00_Data_Prep",
                           "sample_composition.csv"), check.names = FALSE)
clean24 <- comp$sample[comp$cohort == "GSE58485" & comp$include_clean]
stopifnot(length(clean24) == 24)

ex_full <- read.csv(file.path(ROOT, "data", "01_GEO_bulk",
                             "GSE58485_expr.csv"),
                    row.names = 1, check.names = FALSE)
ex_full <- as.matrix(ex_full[, clean24, drop = FALSE])
ph <- read.csv(file.path(ROOT, "data", "01_GEO_bulk", "GSE58485_pheno.csv"),
               row.names = 1, check.names = FALSE)

# 覆盖过滤（单队列 ≥10% 非零）
keep <- rowMeans(ex_full > 0, na.rm = TRUE) >= COV_THRESHOLD
ex <- ex_full[keep, , drop = FALSE]
log_msg("full-universe genes before/after coverage filter:",
        nrow(ex_full), "->", nrow(ex))

# 表型：group + genotype（与 70 一致）
group <- ifelse(comp$group[match(clean24, comp$sample)] == "tbi", "TBI", "Sham")
gt <- as.character(ph$`genotype:ch1`)[match(clean24, rownames(ph))]
gt[is.na(gt) | gt == ""] <- "wildtype"
ph58 <- data.frame(group = factor(group, levels = c("Sham", "TBI")),
                   gt = factor(gt), row.names = clean24)
design58 <- model.matrix(~ group + gt, data = ph58)

# ---------- 2. limma DEG ----------
fit <- lmFit(ex, design58)
fit <- eBayes(fit)
tt <- topTable(fit, coef = "groupTBI", number = Inf, sort.by = "none")
tt$gene <- rownames(tt)
deg <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.585, ]
log_msg("DEGs =", nrow(deg))
write_csv(tt, file.path(OUT, "deg_all_full_universe.csv"))
write_csv(deg, file.path(OUT, "deg_significant_full_universe.csv"))

# ---------- 3. WGCNA（top 8,000 MAD；此处 8,000 上限生效）----------
mad <- apply(ex, 1, mad)
top_genes <- names(sort(mad, decreasing = TRUE))[seq_len(min(8000, length(mad)))]
gene_var <- apply(ex[top_genes, ], 1, var, na.rm = TRUE)
top_genes <- top_genes[gene_var > 0]
datExpr <- t(ex[top_genes, ])

sft <- pickSoftThreshold(datExpr, powerVector = 1:20, verbose = 0)
softPower <- sft$powerEstimate
if (is.na(softPower)) softPower <- 12
log_msg("softPower =", softPower)
sft_tab <- data.frame(power = sft$fitIndices$Power,
                      SFT.R.sq = signif(sft$fitIndices$SFT.R.sq, 4),
                      slope = signif(sft$fitIndices$slope, 3),
                      mean.k = signif(sft$fitIndices$mean.k., 3),
                      median.k = signif(sft$fitIndices$median.k., 3))
write_csv(sft_tab, file.path(OUT, "wgcna_soft_threshold_full_universe.csv"))

adj <- adjacency(datExpr, power = softPower)
TOM <- TOMsimilarity(adj, verbose = 0)
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE)
dynamicColors <- labels2colors(dynamicMods)
names(dynamicColors) <- colnames(datExpr)

MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
group_num <- as.numeric(ph58$group == "TBI")
mod_trait <- cor(MEs, group_num, use = "p")
mod_p <- corPvalueStudent(mod_trait, nSamples = nrow(MEs))
sig_mods <- rownames(mod_p)[mod_p[, 1] < 0.05 & !grepl("^MEgrey$", rownames(mod_p))]
sig_colors <- sub("^ME", "", sig_mods)
module_genes <- unique(colnames(datExpr)[dynamicColors %in% sig_colors])
best <- rownames(mod_trait)[which.max(abs(mod_trait[, 1]))]
best_genes <- unique(colnames(datExpr)[dynamicColors == sub("^ME", "", best)])

mod_size <- data.frame(module = names(table(dynamicColors)),
                       n_genes = as.integer(table(dynamicColors)))
mod_size <- mod_size[order(-mod_size$n_genes), ]
write_csv(mod_size, file.path(OUT, "wgcna_module_sizes_full_universe.csv"))

mod_tab <- data.frame(module = rownames(mod_trait),
                      r = unname(mod_trait[, 1]),
                      p = unname(mod_p[, 1]),
                      significant = unname(mod_p[, 1] < 0.05 &
                                             !grepl("^MEgrey$", rownames(mod_trait))))
write_csv(mod_tab, file.path(OUT, "module_trait_full_universe.csv"))
log_msg("modules:", paste(sub("^ME", "", sig_colors), collapse = ","),
        " best:", sub("^ME", "", best),
        " r =", round(mod_trait[best, 1], 4),
        " P =", format(mod_p[best, 1], digits = 3))

key <- intersect(deg$gene, module_genes)
key_best <- intersect(deg$gene, best_genes)

# ---------- 4. 可评估靶点（1:1 直系同源 ∩ 完整宇宙矩阵）----------
orth <- read_csv(file.path(ROOT, "results", "16_Orthology",
                           "orthology_mapping_full.csv"),
                 col_types = cols(.default = col_character()))
one2one <- orth[orth$orthology_class == "one2one" & orth$mapped == "TRUE", ]
one2one <- one2one[!duplicated(one2one$gene_symbol), ]
evaluable_syms <- unique(one2one$mouse_symbol)
evaluable <- evaluable_syms[evaluable_syms %in% rownames(ex)]
log_msg("1:1 ortholog targets:", nrow(one2one),
        " evaluable in full-universe matrix:", length(evaluable))

tier_map <- setNames(one2one$tier, toupper(one2one$gene_symbol))
mouse2human <- setNames(one2one$gene_symbol, one2one$mouse_symbol)

candidates <- function(kg) {
  cand <- kg[toupper(kg) %in% toupper(evaluable)]
  data.frame(gene = cand,
             human_symbol = unname(mouse2human[cand]),
             tier = ifelse(is.na(tier_map[toupper(cand)]), "C",
                           unname(tier_map[toupper(cand)])),
             stringsAsFactors = FALSE)
}
cand <- candidates(key)
cand_best <- candidates(key_best)
log_msg("key genes =", length(key), " candidates =", nrow(cand),
        " best-module candidates =", nrow(cand_best))
if (nrow(cand) > 0) log_msg("candidates:", paste(cand$gene, collapse = ","))

# ---------- 5. 置换：均匀 + 均值-方差分层 ----------
set.seed(SEED)
K <- nrow(cand)
universe <- rownames(ex)
key_set <- intersect(key, universe)

gene_mean <- rowMeans(ex, na.rm = TRUE)
gene_var <- apply(ex, 1, var, na.rm = TRUE)
mean_dec <- cut(gene_mean, breaks = quantile(gene_mean, probs = seq(0, 1, 0.1)),
                labels = FALSE, include.lowest = TRUE)
var_dec <- cut(gene_var, breaks = quantile(gene_var, probs = seq(0, 1, 0.1)),
               labels = FALSE, include.lowest = TRUE)
strata <- data.frame(gene = universe, mean_dec = mean_dec, var_dec = var_dec,
                     stringsAsFactors = FALSE)
strata_tab <- table(mean_dec, var_dec)
write.csv(as.data.frame.matrix(strata_tab), file.path(OUT, "strata_counts_full_universe.csv"))

perm_one <- function(cand_genes, mode = c("uniform", "meanvar")) {
  mode <- match.arg(mode)
  overlaps <- numeric(N_PERM)
  if (mode == "uniform") {
    for (i in seq_len(N_PERM)) {
      s <- sample(universe, K, replace = FALSE)
      overlaps[i] <- length(intersect(s, key_set))
    }
  } else {
    cand_strata <- strata[strata$gene %in% cand_genes, ]
    for (i in seq_len(N_PERM)) {
      s <- vapply(seq_len(nrow(cand_strata)), function(j) {
        pool <- strata$gene[strata$mean_dec == cand_strata$mean_dec[j] &
                            strata$var_dec == cand_strata$var_dec[j]]
        sample(pool, 1)
      }, character(1))
      overlaps[i] <- length(intersect(s, key_set))
    }
  }
  observed <- length(intersect(cand_genes, key_set))
  exp_val <- if (mode == "uniform") K * length(key_set) / length(universe) else
    mean(overlaps)
  data.frame(path = if (mode == "uniform") "uniform" else "meanvar_stratified",
             K = K, universe = length(universe), key_genes = length(key_set),
             observed = observed, expected = exp_val,
             enrichment = observed / exp_val,
             empirical_P = (1 + sum(overlaps >= observed)) / (N_PERM + 1),
             null_95_lo = quantile(overlaps, 0.025),
             null_95_hi = quantile(overlaps, 0.975))
}

if (K > 0) {
  perm_unif <- perm_one(cand$gene, "uniform")
  perm_mv <- perm_one(cand$gene, "meanvar")
  perm_res <- rbind(perm_unif, perm_mv)
} else {
  perm_res <- data.frame(path = c("uniform", "meanvar_stratified"), K = 0,
                         universe = length(universe), key_genes = length(key_set),
                         observed = 0, expected = 0, enrichment = NA,
                         empirical_P = NA, null_95_lo = 0, null_95_hi = 0)
}
write_csv(perm_res, file.path(OUT, "permutation_full_universe.csv"))
print(perm_res)

# ---------- 6. 保存汇总 ----------
write_csv(data.frame(gene = module_genes), file.path(OUT, "module_genes_full_universe.csv"))
write_csv(data.frame(gene = best_genes), file.path(OUT, "best_module_genes_full_universe.csv"))
write_csv(data.frame(gene = key), file.path(OUT, "key_genes_full_universe.csv"))
write_csv(data.frame(gene = key_best), file.path(OUT, "key_genes_best_module_full_universe.csv"))
write_csv(cand, file.path(OUT, "candidates_full_universe.csv"))
write_csv(cand_best, file.path(OUT, "candidates_best_module_full_universe.csv"))
write_csv(data.frame(mouse_symbol = evaluable), file.path(OUT, "evaluable_targets_full_universe.csv"))

# 关键基因的 AhR 电池状态（供手稿表述使用）
ahr_genes <- c("Ahr", "Cyp1a1", "Cyp1a2", "Cyp1b1", "Ahrr", "Nqo1", "Maoa")
ahr_tab <- data.frame(
  gene = ahr_genes,
  in_full_universe = ahr_genes %in% rownames(ex),
  in_deg = ahr_genes %in% deg$gene,
  in_module = ahr_genes %in% module_genes,
  in_key = ahr_genes %in% key,
  in_candidates = ahr_genes %in% cand$gene,
  stringsAsFactors = FALSE)
write_csv(ahr_tab, file.path(OUT, "ahr_battery_full_universe.csv"))
print(ahr_tab)

writeLines(c(
  "v7 M1 response: full GPL1261 universe rerun (GSE58485-only, gene-level)",
  paste0("generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("samples: ", ncol(ex)),
  paste0("genes full / after coverage: ", nrow(ex_full), " / ", nrow(ex)),
  paste0("DEGs: ", nrow(deg)),
  paste0("WGCNA softPower: ", softPower, "; SFT R^2 at chosen power: ",
         signif(sft$fitIndices$SFT.R.sq[sft$fitIndices$Power == softPower], 4)),
  paste0("best module: ", sub("^ME", "", best), " r=", round(mod_trait[best, 1], 4),
         " P=", format(mod_p[best, 1], digits = 3)),
  paste0("module genes: ", length(module_genes), "; key genes: ", length(key)),
  paste0("evaluable targets: ", length(evaluable)),
  paste0("candidates: ", nrow(cand), " -> ",
         if (nrow(cand)) paste(cand$gene, collapse = ",") else "none"),
  paste0("best-module candidates: ", nrow(cand_best))
), file.path(OUT, "summary_full_universe.txt"))
log_msg("done")
