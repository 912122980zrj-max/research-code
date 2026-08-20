# ============================================================
# 63_sensitivity_matrix.R —— B-3：敏感性汇总矩阵
# 依据：修改报告_BaP_TBI_最终版.md 阶段 B-3
#   每行：FULL / CLEAN / GSE58485-only / removeBatchEffect / SVA / LOCO；
#   列：候选基因 + Tier + 结论是否改变（含 GSE58485-only WGCNA 重跑、
#   removeBatchEffect/SVA 同套下游）。
# 输出：results/16_Orthology/sensitivity_matrix.csv 与 SENSITIVITY_REPORT.txt
# ============================================================

suppressMessages({
  library(limma)
  library(WGCNA)
  library(sva)
  library(dplyr)
  library(readr)
})

ROOT <- getwd()
OUT <- file.path(ROOT, "results", "16_Orthology")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(OUT, "sensitivity_run.log"), append = TRUE)
}
if (file.exists(file.path(OUT, "sensitivity_run.log"))) {
  file.remove(file.path(OUT, "sensitivity_run.log"))
}

GENO_CODE <- "wildtype"
COV_THRESHOLD <- 0.10

read_bulk <- function(gse) readRDS(file.path(ROOT, "data", "01_GEO_bulk",
                                             paste0(gse, ".rds")))

comp <- read.csv(file.path(ROOT, "results", "00_Data_Prep",
                           "sample_composition.csv"), check.names = FALSE)
clean_samples <- comp$sample[comp$include_clean]
full_samples  <- comp$sample[comp$include_full]

# ---------- 与 03 一致的基础预处理 ----------
prep <- function(use_clean = TRUE) {
  gses <- c("GSE58485", "GSE41345")
  dat <- lapply(gses, read_bulk)
  names(dat) <- gses
  common <- Reduce(intersect, lapply(dat, function(d) rownames(d$expr)))
  expr_all <- do.call(cbind, lapply(gses, function(g) {
    dat[[g]]$expr[common, , drop = FALSE]
  }))
  keep <- if (use_clean) clean_samples else full_samples
  expr_all <- expr_all[, intersect(keep, colnames(expr_all)), drop = FALSE]
  cohort_vec <- comp$cohort[match(colnames(expr_all), comp$sample)]
  grp_map <- setNames(comp$group[match(colnames(expr_all), comp$sample)],
                      colnames(expr_all))
  keep_gene <- rowSums(is.na(expr_all)) < 0.2 * ncol(expr_all)
  expr_all <- expr_all[keep_gene, , drop = FALSE]
  cov_a <- rowMeans(expr_all[, cohort_vec == "GSE58485", drop = FALSE] > 0,
                    na.rm = TRUE)
  cov_b <- rowMeans(expr_all[, cohort_vec == "GSE41345", drop = FALSE] > 0,
                    na.rm = TRUE)
  keep_cov <- cov_a >= COV_THRESHOLD & cov_b >= COV_THRESHOLD
  expr_all <- expr_all[keep_cov, , drop = FALSE]
  expr_all <- apply(expr_all, 2, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    x
  })
  list(expr = expr_all, cohort = cohort_vec,
       group = ifelse(grp_map == "tbi", "TBI", "Sham"))
}

# ---------- WGCNA（与 03 完全一致） ----------
run_wgcna <- function(expr_corrected, ph, label) {
  mad <- apply(expr_corrected, 1, mad)
  top_genes <- names(sort(mad, decreasing = TRUE))[
    seq_len(min(8000, length(mad)))]
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
  MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
  MEs <- MEList$eigengenes
  group_num <- as.numeric(ph$group == "TBI")
  mod_trait <- cor(MEs, group_num, use = "p")
  mod_p <- corPvalueStudent(mod_trait, nSamples = nrow(MEs))
  sig_mods <- rownames(mod_p)[mod_p[, 1] < 0.05 &
                                !grepl("^MEgrey$", rownames(mod_p))]
  sig_colors <- sub("^ME", "", sig_mods)
  module_genes <- unique(colnames(datExpr)[dynamicColors %in% sig_colors])
  best <- rownames(mod_trait)[which.max(abs(mod_trait[, 1]))]
  list(module_genes = module_genes,
       best_module = sub("^ME", "", best),
       best_r = unname(mod_trait[best, 1]),
       best_p = unname(mod_p[best, 1]),
       n_modules = length(dynamicColors),
       softPower = softPower)
}

# ---------- DEG ----------
run_deg <- function(expr, ph, design = NULL) {
  if (is.null(design)) {
    design <- model.matrix(~ group, data = ph)
  }
  fit <- lmFit(expr, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = "groupTBI", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  deg <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.585, ]
  list(all = tt, deg = deg$gene)
}

# ---------- 候选交集（沿用主分析的 union_targets.tsv 映射 + tier 标注） ----------
bap <- read.delim(file.path(ROOT, "data", "04_BaP_targets",
                            "union_targets.tsv"))
bap_upper <- toupper(unique(bap$gene_symbol))

tier_tab <- read_csv(file.path(OUT, "bap_targets_364_tier.csv"),
                     col_types = cols(.default = col_character()))
tier_map <- setNames(tier_tab$tier, toupper(tier_tab$gene_symbol))

# B-1 直系同源可评估集（matrix_symbol）
eval_clean_orth <- read_csv(file.path(OUT, "evaluable_targets_CLEAN.csv"),
                            col_types = cols(.default = col_character()))
eval_full_orth  <- read_csv(file.path(OUT, "evaluable_targets_FULL.csv"),
                            col_types = cols(.default = col_character()))
eval_clean_upper <- toupper(unique(eval_clean_orth$matrix_symbol))
eval_full_upper  <- toupper(unique(eval_full_orth$matrix_symbol))

candidates_tiered <- function(key_genes, eval_upper = eval_clean_upper) {
  cand <- key_genes[toupper(key_genes) %in% eval_upper]
  data.frame(gene = cand,
             tier = ifelse(is.na(tier_map[toupper(cand)]), "C",
                           unname(tier_map[toupper(cand)])),
             stringsAsFactors = FALSE)
}

# ============================================================
# 路径 1/2：CLEAN 与 FULL（读现有结果，不重跑）
# ============================================================
read_existing <- function(suffix = "") {
  key <- read.csv(file.path(ROOT, "results", "04_WGCNA",
                            paste0("tbi_key_genes", suffix, ".csv")))$gene
  mod <- read.csv(file.path(ROOT, "results", "04_WGCNA",
                            paste0("module_genes", suffix, ".csv")))$gene
  deg <- read.csv(file.path(ROOT, "results", "03_DEG_Analysis",
                            paste0("deg_significant", suffix, ".csv")))$gene
  cand <- candidates_tiered(key)
  list(key = key, mod = mod, deg = deg, cand = cand)
}

clean_res <- read_existing("")
full_res <- read_existing("_FULL")
# 用直系同源可评估集重算 CLEAN/FULL 候选（key genes 不变）
clean_res$cand <- candidates_tiered(clean_res$key, eval_clean_upper)
full_res$cand <- candidates_tiered(full_res$key, eval_full_upper)

# ============================================================
# 路径 3：GSE58485-only（含基因型协变量；单队列 WGCNA 重跑）
# ============================================================
log_msg("GSE58485-only WGCNA starting")
pp58 <- prep(TRUE)
keep58 <- pp58$cohort == "GSE58485"
ex58 <- pp58$expr[, keep58, drop = FALSE]
ph58 <- data.frame(group = factor(pp58$group[keep58], levels = c("Sham", "TBI")))
d58 <- read_bulk("GSE58485")
gt58 <- as.character(d58$pheno$`genotype:ch1`)[match(colnames(ex58),
                                                     rownames(d58$pheno))]
gt58[is.na(gt58) | gt58 == ""] <- GENO_CODE
ph58$gt <- factor(gt58)
design58 <- model.matrix(~ group + gt, data = ph58)
deg58 <- run_deg(ex58, ph58, design58)
wgcna58 <- run_wgcna(ex58, ph58, "GSE58485_only")
key58 <- intersect(deg58$deg, wgcna58$module_genes)
cand58 <- candidates_tiered(key58, eval_clean_upper)
write.csv(data.frame(gene = key58),
          file.path(OUT, "key_genes_GSE58485_only.csv"), row.names = FALSE)
log_msg("GSE58485-only done: DEGs=", length(deg58$deg),
        " modules=", length(wgcna58$module_genes),
        " key=", length(key58))

# ============================================================
# 路径 4：removeBatchEffect（CLEAN 样本，批次=cohort，协变量=group）
# ============================================================
log_msg("removeBatchEffect path starting")
pp_rbe <- prep(TRUE)
design_rbe <- model.matrix(~ group, data = data.frame(
  group = factor(pp_rbe$group, levels = c("Sham", "TBI"))))
ex_rbe <- removeBatchEffect(pp_rbe$expr, batch = pp_rbe$cohort,
                            design = design_rbe)
ph_rbe <- data.frame(group = factor(pp_rbe$group, levels = c("Sham", "TBI")))
deg_rbe <- run_deg(ex_rbe, ph_rbe)
wgcna_rbe <- run_wgcna(ex_rbe, ph_rbe, "removeBatchEffect")
key_rbe <- intersect(deg_rbe$deg, wgcna_rbe$module_genes)
cand_rbe <- candidates_tiered(key_rbe, eval_clean_upper)
write.csv(data.frame(gene = key_rbe),
          file.path(OUT, "key_genes_removeBatchEffect.csv"), row.names = FALSE)
log_msg("removeBatchEffect done: DEGs=", length(deg_rbe$deg),
        " key=", length(key_rbe))

# ============================================================
# 路径 5：SVA（CLEAN 样本；SV 作为协变量残差化后同套下游）
# ============================================================
log_msg("SVA path starting")
pp_sva <- prep(TRUE)
mod <- model.matrix(~ group, data = data.frame(
  group = factor(pp_sva$group, levels = c("Sham", "TBI"))))
mod0 <- model.matrix(~ 1, data = data.frame(group = pp_sva$group))
n_sv <- tryCatch(
  sva::num.sv(pp_sva$expr, mod, method = "leek"),
  error = function(e) 0
)
if (n_sv > 0) {
  sv_obj <- sva::sva(pp_sva$expr, mod, mod0, n.sv = n_sv)
  sv_mat <- sv_obj$sv
  ex_sva <- removeBatchEffect(pp_sva$expr, covariates = sv_mat,
                              design = mod)
} else {
  sv_mat <- NULL
  ex_sva <- pp_sva$expr
}
ph_sva <- data.frame(group = factor(pp_sva$group, levels = c("Sham", "TBI")))
deg_sva <- run_deg(ex_sva, ph_sva)
wgcna_sva <- run_wgcna(ex_sva, ph_sva, "SVA")
key_sva <- intersect(deg_sva$deg, wgcna_sva$module_genes)
cand_sva <- candidates_tiered(key_sva, eval_clean_upper)
write.csv(data.frame(gene = key_sva),
          file.path(OUT, "key_genes_SVA.csv"), row.names = FALSE)
log_msg("SVA done: n_sv=", n_sv, " DEGs=", length(deg_sva$deg),
        " key=", length(key_sva))

# ============================================================
# 汇总表
# ============================================================
row_of <- function(label, n_genes, degs, wgcna, key, cand, note = "") {
  data.frame(path = label,
             genes_after_filter = n_genes,
             n_deg = length(degs),
             best_module = ifelse(is.null(wgcna$best_module), "",
                                  wgcna$best_module),
             best_r = as.character(ifelse(is.null(wgcna$best_r), "",
                                          round(wgcna$best_r, 3))),
             best_p = as.character(ifelse(is.null(wgcna$best_p), "",
                                          format(wgcna$best_p, digits = 3))),
             n_key_genes = length(key),
             n_candidates = nrow(cand),
             candidate_genes = ifelse(nrow(cand) == 0, "-",
                                      paste(cand$gene, collapse = ",")),
             candidate_tiers = ifelse(nrow(cand) == 0, "-",
                                      paste(cand$tier, collapse = ",")),
             conclusion_changed = ifelse(
               setequal(cand$gene, clean_res$cand$gene),
               "No", "Yes"),
             note = note,
             stringsAsFactors = FALSE)
}

mat <- bind_rows(
  row_of("CLEAN (main)", nrow(pp58$expr), clean_res$deg, NULL,
         clean_res$key, clean_res$cand, "reference"),
  row_of("FULL", nrow(prep(FALSE)$expr), full_res$deg, NULL,
         full_res$key, full_res$cand, "includes drug-treated samples"),
  row_of("GSE58485-only", ncol(ex58), deg58$deg, wgcna58, key58, cand58,
         "single-cohort WGCNA rerun; genotype covariate included"),
  row_of("removeBatchEffect", nrow(ex_rbe), deg_rbe$deg, wgcna_rbe,
         key_rbe, cand_rbe, "limma removeBatchEffect (batch=cohort)"),
  row_of("SVA", nrow(ex_sva), deg_sva$deg, wgcna_sva, key_sva, cand_sva,
         paste0("sva n.sv=", n_sv))
)

# LOCO 行（ML 迁移，已有结果）
loco <- read.csv(file.path(ROOT, "results", "07_ML_Signature",
                           "cross_cohort_auc.csv"))
single <- read.csv(file.path(ROOT, "results", "07_ML_Signature",
                             "H3_cohort_sensitivity",
                             "single_cohort_model_auc.csv"))
mat$loco_clean <- ifelse(nrow(loco) > 0,
                         paste(loco$test_auc[grepl("CLEAN", loco$composition)],
                               collapse = "/"), "")
mat$loco_full <- ifelse(nrow(loco) > 0,
                        paste(loco$test_auc[grepl("FULL", loco$composition)],
                              collapse = "/"), "")
mat$single_cohort_cv_auc <- ifelse(nrow(single) > 0,
                                   paste(single$test_auc, collapse = "/"), "")

write_csv(mat, file.path(OUT, "sensitivity_matrix.csv"))

sink(file.path(OUT, "SENSITIVITY_REPORT.txt"))
cat("B-3 sensitivity matrix\n")
cat("generated:", format(Sys.time()), "\n")
cat("preprocessing: CLEAN samples, common genes, NA filter, coverage >=10%,",
    "median imputation; WGCNA top 8000 MAD (same as main).\n")
cat("candidate mapping: B-1 Ensembl Compara 1:1 orthology evaluable sets",
    "(CLEAN 69 / FULL 68); GSE58485-only/RBE/SVA use the CLEAN 69 set.\n\n")
print(as.data.frame(mat), row.names = FALSE)
sink()
log_msg("B-3 done")
print(as.data.frame(mat), row.names = FALSE)
