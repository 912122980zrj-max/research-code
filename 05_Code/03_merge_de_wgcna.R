# ============================================================
# 03_merge_de_wgcna.R —— 多队列合并、批次校正、DEG、WGCNA、
#                        BaP 靶点交集与 GO/KEGG 富集
# 运行：Rscript scripts/03_merge_de_wgcna.R
# 输入：data/01_GEO_bulk/*_expr.csv + *_pheno.csv（来自 01 脚本）
#       data/04_BaP_targets/union_targets.tsv（来自 02 脚本）
#       data/04_BaP_targets/union_targets_tiered.tsv（来自 02b 脚本，可选）
# 输出：results/00_Data_Prep、03_DEG_Analysis、04_WGCNA、05_Candidate_Intersection、06_GO_KEGG_Enrichment
#
# 2026-08-09 v4 定稿版：
#   - 在公共基因合并、NA 过滤之后、中位数填补与 ComBat 之前，增加
#     每队列表达覆盖过滤：两队列内均 >= COV_THRESHOLD（默认 0.10）样本
#     非零（log2 尺度下零值仍为 0，判定有效），消除平台覆盖伪影导致的
#     ComBat 后 PC1 残留队列结构；
#   - 报告过滤前后基因数与 ComBat 后 PC1 的 cohort eta²（验收 <=0.001）；
#   - GSE41345 无基因型列时默认赋 wildtype（C57BL/6 野生型假设）；可通过
#     环境变量 GENO_CODE=other 复现 other 编码敏感性口径；
#   - 环境变量 OUT_SUFFIX 可让敏感性口径写入带后缀文件而不覆盖主结果
#     （如 GENO_CODE=other OUT_SUFFIX=_GENO_OTHER Rscript scripts/03_merge_de_wgcna.R）；
#   - 环境变量 COV_THRESHOLD 可复现 0 / 0.5 等覆盖阈值敏感性。
#
# 历史说明：
#   - 2026-08-07 审查修复版（T1/T5/P2/P7/P8/P14）：
#   - FULL（42 样本）与 CLEAN（33 样本，剔除环磷酰胺 + Ex-4-mTBI）双口径同时运行；
#     CLEAN 为主结果（无后缀），FULL 输出带 _FULL 后缀归档对比。
#   - 输出样本构成表 sample_composition.csv 与 clean_vs_full_comparison.csv。
#   - GO/KEGG 富集增加 n_genes / is_1gene 列，并另存 Count>=3 的 *_keep.csv。
# ============================================================

suppressMessages({
  library(limma)
  library(sva)
  library(WGCNA)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(dplyr)
  library(readr)
})
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ---------- v4 敏感性环境变量 ----------
GENO_CODE <- Sys.getenv("GENO_CODE", unset = "wildtype")
COV_THRESHOLD <- as.numeric(Sys.getenv("COV_THRESHOLD", unset = "0.10"))
OUT_SUFFIX <- Sys.getenv("OUT_SUFFIX", unset = "")

# 与 01 脚本一致的分组推断（词边界，避免 uninjured 误判）
infer_group <- function(titles) {
  t <- tolower(titles)
  is_ctrl <- grepl("\\b(sham|control|uninjured|vehicle|naive)\\b", t)
  is_tbi <- grepl("\\b(tbi|mtbi|cci|injured|injury|injur|impact|fp)\\b", t) |
    grepl("fluid percussion", t)
  ifelse(is_tbi & !is_ctrl, "TBI",
         ifelse(is_ctrl & !is_tbi, "Sham", "Other"))
}

# ---------- 路径 ----------
if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT_PREP <- file.path(ROOT, "results", "00_Data_Prep")
OUT_DEG <- file.path(ROOT, "results", "03_DEG_Analysis")
OUT_WGCNA <- file.path(ROOT, "results", "04_WGCNA")
OUT_CAND <- file.path(ROOT, "results", "05_Candidate_Intersection")
OUT_GO <- file.path(ROOT, "results", "06_GO_KEGG_Enrichment")
for (od in c(OUT_PREP, OUT_DEG, OUT_WGCNA, OUT_CAND, OUT_GO)) dir.create(od, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts", "theme_pub.R"))

# ---------- 读入数据 ----------
read_bulk <- function(gse) {
  rds_file <- file.path(ROOT, "data", "01_GEO_bulk", paste0(gse, ".rds"))
  if (file.exists(rds_file)) {
    return(readRDS(rds_file))
  }
  expr_file <- file.path(ROOT, "data", "01_GEO_bulk", paste0(gse, "_expr.csv"))
  pheno_file <- file.path(ROOT, "data", "01_GEO_bulk", paste0(gse, "_pheno.csv"))
  if (!file.exists(expr_file)) return(NULL)
  ex <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE,
                           fileEncoding = "UTF-8"))
  ph <- read.csv(pheno_file, row.names = 1, check.names = FALSE,
                 fileEncoding = "UTF-8")
  list(expr = ex, pheno = ph)
}

# 发现队列：GSE58485 + GSE41345
DISCOVERY <- c("GSE58485", "GSE41345")
dat <- lapply(DISCOVERY, read_bulk)
names(dat) <- DISCOVERY
dat <- dat[!vapply(dat, is.null, logical(1))]
if (length(dat) < 1) stop("没有找到发现队列数据，请先运行 01_download_geo.R")

# ---------- 合并与分组 ----------
common_genes <- Reduce(intersect, lapply(dat, function(d) rownames(d$expr)))
expr_list <- lapply(names(dat), function(gse) {
  d <- dat[[gse]]
  d$expr[common_genes, , drop = FALSE]
})
expr_all <- do.call(cbind, expr_list)
pheno_all <- bind_rows(lapply(names(dat), function(gse) {
  p <- dat[[gse]]$pheno
  p$title <- as.character(p$title)
  p$group <- infer_group(p$title)
  p$cohort <- gse
  cols <- intersect(c("geo_accession", "title", "group", "cohort", "subseries"),
                    colnames(p))
  if ("genotype:ch1" %in% colnames(p)) {
    cols <- c(cols, "genotype:ch1")
  }
  p[, cols, drop = FALSE]
}))

pheno_all$group <- tolower(pheno_all$group)
if ("genotype:ch1" %in% colnames(pheno_all)) {
  pheno_all$genotype <- as.character(pheno_all$`genotype:ch1`)
} else {
  pheno_all$genotype <- GENO_CODE
}
pheno_all$genotype[is.na(pheno_all$genotype) | pheno_all$genotype == ""] <- GENO_CODE

# 药物处理亚组标记（审查 P2：环磷酰胺 / Ex-4-mTBI）
pheno_all$exclude_drug <- grepl("cyclophosphamide|ex-4-mtbi",
                                pheno_all$title, ignore.case = TRUE)

# 样本构成表（P14）
composition <- data.frame(
  sample = colnames(expr_all),
  cohort = pheno_all$cohort,
  group = pheno_all$group,
  title = pheno_all$title,
  genotype = pheno_all$genotype,
  exclude_drug = pheno_all$exclude_drug,
  include_full = pheno_all$group %in% c("sham", "tbi"),
  include_clean = pheno_all$group %in% c("sham", "tbi") & !pheno_all$exclude_drug,
  stringsAsFactors = FALSE
)
write.csv(composition, file.path(OUT_PREP, "sample_composition.csv"),
          row.names = FALSE)

# ---------- 单一口径分析主函数 ----------
run_mode <- function(suffix = "", exclude_drug = FALSE) {
  suffix <- paste0(suffix, OUT_SUFFIX)
  keep <- pheno_all$group %in% c("sham", "tbi")
  if (exclude_drug) keep <- keep & !pheno_all$exclude_drug
  ex <- expr_all[, keep, drop = FALSE]
  ph <- pheno_all[keep, , drop = FALSE]
  n_excluded_drug <- if (exclude_drug) {
    sum(pheno_all$group %in% c("sham", "tbi") & pheno_all$exclude_drug)
  } else {
    0
  }
  ph$group <- factor(ph$group, levels = c("sham", "tbi"))
  ph$cohort <- factor(ph$cohort)
  ph$genotype <- as.character(ph$genotype)

  report <- file.path(OUT_PREP, paste0("03_bulk_report", suffix, ".txt"))
  writeLines(paste0("03_merge_de_wgcna report", suffix), report)
  write_report <- function(...) {
    cat(..., "\n", sep = "", file = report, append = TRUE)
  }
  write_report("generated:", format(Sys.time()))
  write_report("cohorts:", paste(DISCOVERY, collapse = ","))
  write_report("samples:", ncol(ex), "| sham:", sum(ph$group == "sham"),
               "| tbi:", sum(ph$group == "tbi"),
               "| excluded_drug_samples:", n_excluded_drug)
  write_report("genotype covariate levels:", paste(unique(ph$genotype),
                                                    collapse = ","))
  write_report("GENO_CODE:", GENO_CODE)
  write_report("coverage threshold (per-cohort nonzero fraction):",
               COV_THRESHOLD)
  write_report("table of groups by cohort:")
  write_report(paste(capture.output(print(table(ph$cohort, ph$group))),
                     collapse = "\n"))

  # ---------- 预处理：过滤、缺失值、批次校正 ----------
  keep_gene <- rowSums(is.na(ex)) < 0.2 * ncol(ex)
  ex <- ex[keep_gene, ]
  write_report("genes after NA filter:", nrow(ex))

  # v4：每队列覆盖过滤（两队列内均 >= COV_THRESHOLD 样本非零）
  n_before_cov <- nrow(ex)
  cov_a <- rowMeans(ex[, ph$cohort == "GSE58485", drop = FALSE] > 0,
                    na.rm = TRUE)
  cov_b <- rowMeans(ex[, ph$cohort == "GSE41345", drop = FALSE] > 0,
                    na.rm = TRUE)
  keep_cov <- cov_a >= COV_THRESHOLD & cov_b >= COV_THRESHOLD
  ex <- ex[keep_cov, , drop = FALSE]
  write_report("genes before coverage filter:", n_before_cov)
  write_report("genes after coverage filter:", nrow(ex),
               "(removed:", n_before_cov - nrow(ex), ")")
  write_report("mean coverage of retained genes (GSE58485 / GSE41345):",
               paste0(round(mean(cov_a[keep_cov]) * 100, 1), "% / ",
                      round(mean(cov_b[keep_cov]) * 100, 1), "%"))

  ex <- apply(ex, 2, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    x
  })

  # ---------- PC1 上 cohort 的 eta²（批次残留指标） ----------
  eta2_cohort <- function(x, cohort) {
    aov_fit <- anova(lm(x ~ factor(cohort)))
    ss_between <- aov_fit$`Sum Sq`[1]
    ss_total <- sum(aov_fit$`Sum Sq`)
    ss_between / ss_total
  }
  pc <- prcomp(t(ex), scale. = TRUE)
  eta2_before <- eta2_cohort(pc$x[, 1], ph$cohort)
  pca_df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                       cohort = ph$cohort, group = ph$group)
  pct <- round(100 * summary(pc)$importance[2, 1:2], 1)
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cohort, shape = group)) +
    geom_point(size = 2.2, alpha = 0.9) +
    scale_color_manual(values = pal_cat[1:2], name = "Cohort") +
    scale_shape_manual(values = c(16, 17), name = "Group") +
    labs(title = "PCA before ComBat",
         x = paste0("PC1 (", pct[1], "%)"),
         y = paste0("PC2 (", pct[2], "%)")) +
    theme_pub()
  ggsave_pub(file.path(OUT_PREP, paste0("pca_before", suffix, ".pdf")), p_pca,
         width = 6, height = 5, dpi = 300)

  mod <- model.matrix(~ group, data = ph)
  expr_corrected <- ComBat(dat = ex, batch = ph$cohort,
                           mod = mod, par.prior = TRUE, prior.plots = FALSE)
  write.csv(expr_corrected, file.path(OUT_PREP, paste0("corrected_expr", suffix,
                                                  ".csv")),
            quote = FALSE)
  write.csv(ph, file.path(OUT_PREP, paste0("corrected_pheno", suffix, ".csv")),
            row.names = FALSE)

  pc2 <- prcomp(t(expr_corrected), scale. = TRUE)
  eta2_after <- eta2_cohort(pc2$x[, 1], ph$cohort)
  write_report("PC1 cohort eta2 before ComBat:", round(eta2_before, 6))
  write_report("PC1 cohort eta2 after ComBat:", round(eta2_after, 6),
               "| acceptance <=0.001 ->",
               ifelse(eta2_after <= 0.001, "PASS", "FAIL"))
  pca2_df <- data.frame(PC1 = pc2$x[, 1], PC2 = pc2$x[, 2],
                        cohort = ph$cohort, group = ph$group)
  pct2 <- round(100 * summary(pc2)$importance[2, 1:2], 1)
  p_pca2 <- ggplot(pca2_df, aes(x = PC1, y = PC2, color = cohort,
                                shape = group)) +
    geom_point(size = 2.2, alpha = 0.9) +
    scale_color_manual(values = pal_cat[1:2], name = "Cohort") +
    scale_shape_manual(values = c(16, 17), name = "Group") +
    labs(title = "PCA after ComBat",
         x = paste0("PC1 (", pct2[1], "%)"),
         y = paste0("PC2 (", pct2[2], "%)")) +
    theme_pub()
  ggsave_pub(file.path(OUT_PREP, paste0("pca_after", suffix, ".pdf")), p_pca2,
         width = 6, height = 5, dpi = 300)

  # ---------- 差异表达（limma） ----------
  design <- model.matrix(~ group + genotype, data = ph)
  fit <- lmFit(expr_corrected, design)
  fit <- eBayes(fit)
  deg_all <- topTable(fit, coef = "grouptbi", number = Inf, sort.by = "none")
  deg_all$gene <- rownames(deg_all)
  deg <- deg_all[deg_all$adj.P.Val < 0.05 & abs(deg_all$logFC) > 0.585, ]
  write_report("DEGs (|log2FC|>0.585, adjP<0.05):", nrow(deg))
  write_report("DEG up / down:",
               sum(deg$logFC > 0), "/", sum(deg$logFC < 0))
  write.csv(deg_all, file.path(OUT_DEG, paste0("deg_all", suffix, ".csv")),
            row.names = FALSE)
  write.csv(deg, file.path(OUT_DEG, paste0("deg_significant", suffix, ".csv")),
            row.names = FALSE)

  # ---------- WGCNA（top 8000 MAD 基因） ----------
  mad <- apply(expr_corrected, 1, mad)
  top_genes <- names(sort(mad, decreasing = TRUE))[seq_len(min(8000, length(mad)))]
  gene_var <- apply(expr_corrected[top_genes, ], 1, var, na.rm = TRUE)
  top_genes <- top_genes[gene_var > 0]
  datExpr <- t(expr_corrected[top_genes, ])

  powers <- c(1:20)
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
  softPower <- sft$powerEstimate
  if (is.na(softPower)) softPower <- 12
  write_report("WGCNA softPower:", softPower)

  adjacency <- adjacency(datExpr, power = softPower)
  TOM <- TOMsimilarity(adjacency, verbose = 0)
  dissTOM <- 1 - TOM
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                               deepSplit = 2, pamRespectsDendro = FALSE)
  dynamicColors <- labels2colors(dynamicMods)

  draw_wgcna_dendro(
    geneTree, dynamicColors,
    main = paste0("WGCNA gene dendrogram and module colors (",
                  ifelse(suffix == "", "CLEAN", sub("^_", "", suffix)), ")"),
    out_file = file.path(OUT_WGCNA, paste0("wgcna_dendro", suffix, ".pdf")))

  MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
  MEs <- MEList$eigengenes
  group_num <- as.numeric(ph$group == "tbi")
  mod_trait <- cor(MEs, group_num, use = "p")
  mod_p <- corPvalueStudent(mod_trait, nSamples = nrow(MEs))
  write_report("module-trait correlations (TBI):")
  write_report(paste(rownames(mod_trait), "r=", round(mod_trait[, 1], 3),
                     "P=", format(mod_p[, 1], digits = 3)))
  write.csv(data.frame(module = rownames(mod_trait),
                       r = round(mod_trait[, 1], 4),
                       P = mod_p[, 1],
                       significant = mod_p[, 1] < 0.05 &
                         !grepl("^MEgrey$", rownames(mod_trait))),
            file.path(OUT_WGCNA, paste0("module_trait", suffix, ".csv")),
            row.names = FALSE)

  best_mod <- rownames(mod_trait)[which.max(abs(mod_trait[, 1]))]
  best_mod_color <- sub("^ME", "", best_mod)
  write_report("best module:", best_mod, "r=", round(mod_trait[best_mod, 1], 3),
               "P=", format(mod_p[best_mod, 1], digits = 4))

  write.csv(data.frame(gene = colnames(datExpr), module = dynamicColors),
            file.path(OUT_WGCNA, paste0("module_colors", suffix, ".csv")),
            row.names = FALSE)

  sig_mods <- rownames(mod_p)[mod_p[, 1] < 0.05 & !grepl("^MEgrey$", rownames(mod_p))]
  sig_colors <- sub("^ME", "", sig_mods)
  write_report("TBI-associated modules (P<0.05):", paste(sig_colors, collapse = ","))
  module_genes <- unique(colnames(datExpr)[dynamicColors %in% sig_colors])
  write_report("module gene count (union of significant modules):",
               length(module_genes))
  write.csv(data.frame(gene = module_genes),
            file.path(OUT_WGCNA, paste0("module_genes", suffix, ".csv")),
            row.names = FALSE)

  # ---------- TBI 关键基因 = DEG ∩ 模块基因 ----------
  key_genes <- intersect(rownames(deg), module_genes)
  write_report("TBI key genes (DEG ∩ module):", length(key_genes))
  write.csv(data.frame(gene = key_genes),
            file.path(OUT_WGCNA, paste0("tbi_key_genes", suffix, ".csv")),
            row.names = FALSE)

  # ---------- 与 BaP 靶点取交集 ----------
  n_candidates <- NA_integer_
  n_candidates_ab <- NA_integer_
  bap_file <- file.path(ROOT, "data", "04_BaP_targets", "union_targets.tsv")
  if (file.exists(bap_file)) {
    bap <- read.delim(bap_file)
    bap_genes <- toupper(bap$gene_symbol)
    key_upper <- toupper(key_genes)
    candidate <- key_genes[key_upper %in% bap_genes]
    write_report("BaP targets:", length(unique(bap_genes)),
                 "| inferred BaP-TBI candidates:", length(candidate))
    write.csv(data.frame(gene = candidate),
              file.path(OUT_CAND, paste0("candidate_genes", suffix, ".csv")),
              row.names = FALSE)
    n_candidates <- length(candidate)

    # 分级靶点标注（T3，仅当 02b 已运行）
    tiered_file <- file.path(ROOT, "data", "04_BaP_targets",
                             "union_targets_tiered.tsv")
    if (file.exists(tiered_file) && length(candidate) > 0) {
      tier <- read.delim(tiered_file)
      tier$gene_symbol <- toupper(tier$gene_symbol)
      tier_order <- c(A = 0, B = 1, C = 2)
      tier$tier_num <- tier_order[tier$tier]
      best <- tier[order(tier$tier_num), ]
      best <- best[!duplicated(best$gene_symbol), ]
      cand_up <- toupper(candidate)
      m <- match(cand_up, best$gene_symbol)
      cand_tier <- data.frame(
        gene = candidate,
        tier = ifelse(is.na(m), NA_character_, best$tier[m]),
        source = ifelse(is.na(m), NA_character_,
                        vapply(m, function(i) {
                          paste(unique(tier$source[tier$gene_symbol ==
                                                     best$gene_symbol[i]]),
                                collapse = ";")
                        }, character(1))),
        score_detail = ifelse(is.na(m), NA_character_,
                              paste(best$score[m], best$detail[m], sep = "; ")),
        stringsAsFactors = FALSE
      )
      write.csv(cand_tier, file.path(OUT_CAND, paste0("candidate_genes_tiered",
                                                 suffix, ".csv")),
                row.names = FALSE)

      # 去别名后的候选分级表（第二轮审查通过条件 1：candidate_genes_tiered_v2.csv）
      dedup_file <- file.path(ROOT, "data", "04_BaP_targets",
                              "union_targets_tiered_dedup.tsv")
      if (file.exists(dedup_file)) {
        tier_d <- read.delim(dedup_file)
        tier_d$gene_symbol <- toupper(tier_d$gene_symbol)
        best_d <- tier_d[order(tier_order[tier_d$tier]), ]
        best_d <- best_d[!duplicated(best_d$gene_symbol), ]
        m_d <- match(cand_up, best_d$gene_symbol)
        cand_tier_d <- data.frame(
          gene = candidate,
          tier = ifelse(is.na(m_d), NA_character_, best_d$tier[m_d]),
          source = ifelse(is.na(m_d), NA_character_,
                          vapply(m_d, function(i) {
                            paste(unique(tier_d$source[tier_d$gene_symbol ==
                                                         best_d$gene_symbol[i]]),
                                  collapse = ";")
                          }, character(1))),
          score_detail = ifelse(is.na(m_d), NA_character_,
                                paste(best_d$score[m_d], best_d$detail[m_d],
                                      sep = "; ")),
          stringsAsFactors = FALSE
        )
        write.csv(cand_tier_d, file.path(OUT_CAND, paste0("candidate_genes_tiered_v2",
                                                     suffix, ".csv")),
                  row.names = FALSE)
      }

      cand_ab <- cand_tier[!is.na(cand_tier$tier) & cand_tier$tier %in% c("A", "B"),
                           c("gene", "tier", "source"), drop = FALSE]
      write.csv(cand_ab, file.path(OUT_CAND, paste0("candidate_genes_AB", suffix,
                                               ".csv")),
                row.names = FALSE)
      n_candidates_ab <- nrow(cand_ab)
      write_report("candidates with tier A/B:", n_candidates_ab)
    } else {
      write_report("union_targets_tiered.tsv 不存在，跳过分级交集标注")
    }

    # 三集合 Venn（BaP 靶点 / DEG / 模块基因）
    venn_list <- list(
      BaP = unique(bap_genes),
      DEG = unique(toupper(rownames(deg))),
      Module = unique(toupper(module_genes))
    )
    if (requireNamespace("ggvenn", quietly = TRUE)) {
      p_venn <- ggvenn::ggvenn(venn_list, show_elements = FALSE,
                               fill_color = pal_cat[c(4, 5, 6)],
                               stroke_size = 0.4,
                               set_name_size = 4.5)
      ggsave_pub(file.path(OUT_CAND, paste0("venn_candidates", suffix, ".pdf")),
             p_venn, width = 6, height = 5.5, dpi = 300)
    } else {
  pdf_pub(file.path(OUT_CAND, paste0("venn_candidates", suffix, ".pdf")))
      grid::grid.draw(VennDiagram::venn.diagram(venn_list, filename = NULL))
      dev.off()
    }

    # GO / KEGG 富集（T5：n_genes / is_1gene / Count>=3 过滤）
    if (length(candidate) > 0) {
      eg <- bitr(candidate, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)
      if (nrow(eg) > 0) {
        ego <- enrichGO(gene = eg$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP",
                        pAdjustMethod = "BH", pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2, readable = TRUE)
        if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
          godf <- as.data.frame(ego)
          godf$n_genes <- godf$Count
          godf$is_1gene <- godf$Count == 1
          write.csv(godf, file.path(OUT_GO, paste0("GO_BP_candidates", suffix,
                                                ".csv")),
                    row.names = FALSE)
          godf_keep <- godf[godf$Count >= 3, , drop = FALSE]
          write.csv(godf_keep, file.path(OUT_GO, paste0("GO_BP_keep", suffix,
                                                     ".csv")),
                    row.names = FALSE)
          write_report("GO BP terms (Count>=3):", nrow(godf_keep))
          p_go <- dotplot(ego, showCategory = 15) +
            scale_colour_gradient(low = "#084594", high = "#6BAED6",
                                  name = "p.adjust") +
            theme_pub()
          ggsave_pub(file.path(OUT_GO, paste0("GO_BP_bubble", suffix, ".pdf")),
                 p_go, width = 7, height = 6, dpi = 300)
        }
        ekegg <- NULL
        for (kegg_try in 1:3) {
          ekegg <- tryCatch(enrichKEGG(gene = eg$ENTREZID, organism = "mmu",
                                       pvalueCutoff = 0.05),
                            error = function(e) NULL)
          if (!is.null(ekegg)) break
          if (kegg_try < 3) Sys.sleep(5)
        }
        if (is.null(ekegg)) {
          write_report("KEGG enrichment failed after 3 attempts (network)")
        }
        if (!is.null(ekegg)) {
          kdf <- as.data.frame(ekegg)
          kdf$n_genes <- kdf$Count
          kdf$is_1gene <- kdf$Count == 1
          write.csv(kdf, file.path(OUT_GO, paste0("KEGG_candidates", suffix,
                                                ".csv")),
                    row.names = FALSE)
          kdf_keep <- kdf[kdf$Count >= 3, , drop = FALSE]
          write.csv(kdf_keep, file.path(OUT_GO, paste0("KEGG_keep", suffix,
                                                     ".csv")),
                    row.names = FALSE)
          write_report("KEGG candidate terms:", nrow(kdf))
          write_report("KEGG terms (Count>=3):", nrow(kdf_keep))
          if (nrow(kdf) > 0) {
            p_kegg <- dotplot(ekegg, showCategory = 15) +
              scale_colour_gradient(low = pal_blue[8], high = pal_blue[3],
                                    name = "p.adjust") +
              theme_pub()
            ggsave_pub(file.path(OUT_GO, paste0("KEGG_bubble", suffix, ".pdf")),
                   p_kegg, width = 7, height = 5, dpi = 300)
          }
        }
      }
    }
  } else {
    warning("未找到 data/04_BaP_targets/union_targets.tsv，跳过交集与富集")
  }

  list(n_samples = ncol(ex), n_deg = nrow(deg), softPower = softPower,
       best_module = best_mod_color,
       best_r = round(mod_trait[best_mod, 1], 4),
       best_p = format(mod_p[best_mod, 1], digits = 4),
       n_sig_mods = length(sig_colors), n_module_genes = length(module_genes),
       n_key_genes = length(key_genes), n_candidates = n_candidates,
       n_candidates_ab = n_candidates_ab)
}

# ---------- 运行 FULL 与 CLEAN ----------
metrics_full <- run_mode(suffix = "_FULL", exclude_drug = FALSE)
metrics_clean <- run_mode(suffix = "", exclude_drug = TRUE)

comparison <- data.frame(
  metric = c("n_samples", "n_deg", "softPower", "best_module", "best_r",
             "best_p", "n_sig_mods", "n_module_genes", "n_key_genes",
             "n_candidates", "n_candidates_ab"),
  full = unlist(metrics_full),
  clean = unlist(metrics_clean),
  stringsAsFactors = FALSE
)
if (nzchar(OUT_SUFFIX)) {
  write.csv(comparison,
            file.path(OUT_PREP, paste0("clean_vs_full_comparison",
                                       OUT_SUFFIX, ".csv")),
            row.names = FALSE)
} else {
  write.csv(comparison, file.path(OUT_PREP, "clean_vs_full_comparison.csv"),
            row.names = FALSE)
}

message("完成。CLEAN 主结果已按实验分目录输出（00_Data_Prep / 03_DEG_Analysis / 04_WGCNA / 05_Candidate_Intersection / 06_GO_KEGG_Enrichment）；FULL 归档带 _FULL 后缀。")
message("关键对比见 clean_vs_full_comparison", OUT_SUFFIX, ".csv。")
