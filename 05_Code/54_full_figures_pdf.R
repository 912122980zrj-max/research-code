# ============================================================
# 54_full_figures_pdf.R —— FULL 口径整体重跑 + 历史 PNG 补矢量 PDF
# 授权：用户 2026-08-09；背景：_archive_FULL 全线为 2026-08-07 尺度修复前的
# 旧结果，需在现行 log2/ComBat 矩阵上重跑 DEG→WGCNA→候选→GO→韦恩。
# 旧 FULL 文件先备份至各 *_FULL 目录的 _pre_log2/ 子目录。
# 输出（PDF 同步到 figures/Figure_S_Supplementary/）：
#   wgcna_softThreshold.pdf（CLEAN）、wgcna_softThreshold_FULL.pdf、
#   wgcna_dendro_FULL.pdf、GO_BP_bubble_FULL.pdf、venn_candidates_FULL.pdf
# ============================================================

suppressMessages({
  library(WGCNA)
  library(ggplot2)
  library(readr)
  library(limma)
  library(ggvenn)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})
allowWGCNAThreads()
source("scripts/theme_pub.R")

ROOT <- "E:/sheng xin/tbi"
OUT_PREP <- file.path(ROOT, "results", "00_Data_Prep")
OUT_DEG <- file.path(ROOT, "results", "03_DEG_Analysis", "_archive_FULL")
OUT_W <- file.path(ROOT, "results", "04_WGCNA")
OUT_WF <- file.path(OUT_W, "_archive_FULL")
OUT_GO <- file.path(ROOT, "results", "06_GO_KEGG_Enrichment", "_archive_FULL")
OUT_CF <- file.path(ROOT, "results", "05_Candidate_Intersection", "_archive_FULL")
FIG_S <- file.path(ROOT, "figures", "Figure_S_Supplementary")
dir.create(FIG_S, recursive = TRUE, showWarnings = FALSE)

# ---------- 0) 备份旧 FULL 文件 ----------
backup <- function(files) {
  for (f in files) {
    if (file.exists(f)) {
      dest_dir <- file.path(dirname(f), "_pre_log2")
      dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
      file.copy(f, file.path(dest_dir, basename(f)), overwrite = TRUE)
      cat("backed up:", f, "\n")
    }
  }
}
backup(list.files(OUT_DEG, full.names = TRUE))
backup(list.files(OUT_WF, full.names = TRUE))
backup(list.files(OUT_GO, full.names = TRUE))
backup(list.files(OUT_CF, full.names = TRUE))

# ---------- 1) FULL DEG（limma，现行 log2/ComBat 矩阵） ----------
exF <- as.matrix(read.csv(file.path(OUT_PREP, "corrected_expr_FULL.csv"),
                          row.names = 1, check.names = FALSE))
phF <- read.csv(file.path(OUT_PREP, "corrected_pheno_FULL.csv"),
                check.names = FALSE, fileEncoding = "UTF-8")
phF$group <- factor(phF$group, levels = c("sham", "tbi"))
designF <- model.matrix(~ group + genotype, data = phF)
fitF <- eBayes(lmFit(exF, designF))
degF_all <- topTable(fitF, coef = "grouptbi", number = Inf, sort.by = "none")
degF_all$gene <- rownames(degF_all)
degF_sig <- degF_all[degF_all$adj.P.Val < 0.05 & abs(degF_all$logFC) > 0.585, ]
write.csv(degF_all, file.path(OUT_DEG, "deg_all_FULL.csv"), row.names = FALSE)
write.csv(degF_sig, file.path(OUT_DEG, "deg_significant_FULL.csv"),
          row.names = FALSE)
cat("FULL DEG (log2 现行矩阵):", nrow(degF_sig),
    "（上调", sum(degF_sig$logFC > 0), "下调", sum(degF_sig$logFC < 0), "）\n")

# ---------- 2) FULL WGCNA（现行矩阵；重算 TOM） ----------
run_soft <- function(suffix, out_file) {
  ex <- as.matrix(read.csv(file.path(OUT_PREP,
                                     paste0("corrected_expr", suffix, ".csv")),
                           row.names = 1, check.names = FALSE))
  mad <- apply(ex, 1, mad)
  top <- names(sort(mad, decreasing = TRUE))[seq_len(min(8000, nrow(ex)))]
  top <- top[apply(ex[top, ], 1, var, na.rm = TRUE) > 0]
  datExpr <- t(ex[top, ])
  sft <- pickSoftThreshold(datExpr, powerVector = 1:20, verbose = 0)
  pdf_pub(out_file, width = 9, height = 4.5)
  fit <- sft$fitIndices
  par(mfrow = c(1, 2), family = "sans", ps = 8, mar = c(4.2, 4.2, 2, 1))
  plot(fit$Power, fit$SFT.R.sq, pch = 19, col = "#2166AC", cex = 1.2,
       xlab = "Soft threshold (power)",
       ylab = "Scale-free fit index (R^2)")
  text(fit$Power, fit$SFT.R.sq, labels = fit$Power, pos = 3, cex = 0.8)
  abline(h = 0.85, lty = 2, col = "grey50")
  plot(fit$Power, fit$mean.k., pch = 19, col = "#CC6677", cex = 1.2,
       log = "y", xlab = "Soft threshold (power)",
       ylab = "Mean connectivity")
  text(fit$Power, fit$mean.k., labels = fit$Power, pos = 3, cex = 0.8)
  dev.off()
  cat("soft threshold saved:", basename(out_file),
      "| powerEstimate =", sft$powerEstimate, "\n")
  invisible(sft)
}
run_soft("", file.path(OUT_W, "wgcna_softThreshold.pdf"))
sftF <- run_soft("_FULL", file.path(OUT_WF, "wgcna_softThreshold_FULL.pdf"))

madF <- apply(exF, 1, mad)
topF <- names(sort(madF, decreasing = TRUE))[seq_len(min(8000, nrow(exF)))]
topF <- topF[apply(exF[topF, ], 1, var, na.rm = TRUE) > 0]
datF <- t(exF[topF, ])
softPowerF <- if (is.na(sftF$powerEstimate)) 12 else sftF$powerEstimate
adjF <- adjacency(datF, power = softPowerF)
TOMF <- TOMsimilarity(adjF, verbose = 0)
geneTreeF <- hclust(as.dist(1 - TOMF), method = "average")
modsF <- cutreeDynamic(dendro = geneTreeF, distM = 1 - TOMF,
                       deepSplit = 2, pamRespectsDendro = FALSE)
colorsF <- labels2colors(modsF)
cat("FULL modules:", length(unique(colorsF)), "| sizes:\n")
print(table(colorsF))

MEsF <- moduleEigengenes(datF, colors = colorsF)$eigengenes
grpF <- as.numeric(phF$group == "tbi")
traitF <- cor(MEsF, grpF, use = "p")
traitP <- corPvalueStudent(traitF, nSamples = nrow(MEsF))
modTraitF <- data.frame(
  module = rownames(traitF),
  r = round(traitF[, 1], 4),
  P = traitP[, 1],
  significant = traitP[, 1] < 0.05 & !grepl("^MEgrey$", rownames(traitF))
)
write.csv(modTraitF, file.path(OUT_WF, "module_trait_FULL.csv"),
          row.names = FALSE)
sigModsF <- modTraitF$module[modTraitF$significant]
sigColsF <- sub("^ME", "", sigModsF)
moduleGenesF <- unique(colnames(datF)[colorsF %in% sigColsF])
write.csv(data.frame(gene = colnames(datF), module = colorsF),
          file.path(OUT_WF, "module_colors_FULL.csv"), row.names = FALSE)
write.csv(data.frame(gene = moduleGenesF),
          file.path(OUT_WF, "module_genes_FULL.csv"), row.names = FALSE)
keyF <- intersect(rownames(degF_sig), moduleGenesF)
write.csv(data.frame(gene = keyF),
          file.path(OUT_WF, "tbi_key_genes_FULL.csv"), row.names = FALSE)
cat("FULL significant modules:", paste(sigColsF, collapse = ","),
    "| module genes:", length(moduleGenesF),
    "| key genes:", length(keyF), "\n")

pdf_pub(file.path(OUT_WF, "wgcna_dendro_FULL.pdf"), width = 10, height = 6)
par(family = "sans")
plotDendroAndColors(geneTreeF, colorsF, "Dynamic Tree Cut",
                    dendroLabels = FALSE, cex.colorLabels = 0.9,
                    main = "FULL-pathway WGCNA gene dendrogram (log2)")
dev.off()
cat("FULL dendrogram saved\n")

# ---------- 3) FULL 候选交集与分级 ----------
bap <- read.delim(file.path(ROOT, "data", "04_BaP_targets",
                            "union_targets_tiered_dedup.tsv"),
                  stringsAsFactors = FALSE)
bap$gene_symbol <- toupper(bap$gene_symbol)
bap_t <- toupper(unique(bap$gene_symbol))
keyUp <- toupper(keyF)
candF <- keyF[keyUp %in% bap_t]
write.csv(data.frame(gene = candF),
          file.path(OUT_CF, "candidate_genes_FULL.csv"), row.names = FALSE)
tier_order <- c(A = 0, B = 1, C = 2)
bap_best <- bap[order(tier_order[bap$tier]), ]
bap_best <- bap_best[!duplicated(bap_best$gene_symbol), ]
mF <- match(toupper(candF), bap_best$gene_symbol)
candTierF <- data.frame(
  gene = candF,
  tier = ifelse(is.na(mF), NA, bap_best$tier[mF]),
  source = ifelse(is.na(mF), NA,
                  vapply(mF, function(i) paste(unique(bap_best$source[
                    bap_best$gene_symbol == bap_best$gene_symbol[i]]),
                    collapse = ";"), character(1))),
  score_detail = ifelse(is.na(mF), NA,
                        vapply(mF, function(i) paste(unique(bap_best$detail[
                          bap_best$gene_symbol == bap_best$gene_symbol[i]]),
                          collapse = ";"), character(1))),
  stringsAsFactors = FALSE
)
write.csv(candTierF, file.path(OUT_CF, "candidate_genes_tiered_FULL.csv"),
          row.names = FALSE)
write.csv(candTierF, file.path(OUT_CF, "candidate_genes_tiered_v2_FULL.csv"),
          row.names = FALSE)
write.csv(data.frame(gene = candF[!is.na(candTierF$tier) &
                                    candTierF$tier %in% c("A", "B")]),
          file.path(OUT_CF, "candidate_genes_AB_FULL.csv"), row.names = FALSE)
cat("FULL candidates (log2 现行矩阵):", length(candF), "\n")
print(candTierF)

# ---------- 4) FULL GO（14 候选口径改为新候选口径） ----------
egF <- AnnotationDbi::select(org.Mm.eg.db,
                             keys = candF,
                             columns = "ENTREZID", keytype = "SYMBOL")
egoF <- enrichGO(gene = unique(egF$ENTREZID), OrgDb = org.Mm.eg.db,
                 ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
keepF <- as.data.frame(egoF)
keepF <- keepF[keepF$Count >= 3, ]
write.csv(keepF, file.path(OUT_GO, "GO_BP_keep_FULL_log2.csv"),
          row.names = FALSE)
cat("FULL GO BP terms (Count>=3):", nrow(keepF), "\n")
if (nrow(keepF) > 0) {
  keepF$Description <- factor(keepF$Description,
                              levels = rev(keepF$Description[
                                order(keepF$p.adjust)]))
  p_go <- ggplot(keepF, aes(x = Count, y = Description,
                            size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#2166AC", high = "#F7F7F7", name = "adjP") +
    scale_size_continuous(range = c(2, 6)) +
    labs(x = "Gene count", y = NULL,
         title = "FULL-pathway GO BP (Count>=3, log2)") +
    theme_pub()
} else {
  p_go <- ggplot(data.frame(x = 1, y = 1,
                            lab = "no BP terms with Count>=3"),
                 aes(x, y, label = lab)) +
    geom_text(size = 7) + theme_void()
}
ggsave_pub(file.path(OUT_GO, "GO_BP_bubble_FULL.pdf"), p_go,
           width = 6.5, height = max(4, 0.45 * nrow(keepF) + 2))
cat("FULL GO bubble saved\n")

# ---------- 5) FULL 三集合韦恩 ----------
vennF <- list(DEG = toupper(rownames(degF_sig)),
              WGCNA_module = toupper(moduleGenesF),
              BaP_targets = bap_t)
p_venn <- ggvenn(vennF, fill_color = pal_cat[1:3], stroke_size = 0.4,
                 text_size = 3.2, set_name_size = 4)
ggsave_pub(file.path(OUT_CF, "venn_candidates_FULL.pdf"), p_venn,
           width = 6, height = 5.5)
cat("FULL venn intersection:",
    length(Reduce(intersect, vennF)), "\n")

# ---------- 6) 同步到 figures/Figure_S_Supplementary ----------
for (f in c(file.path(OUT_W, "wgcna_softThreshold.pdf"),
            file.path(OUT_WF, "wgcna_softThreshold_FULL.pdf"),
            file.path(OUT_WF, "wgcna_dendro_FULL.pdf"),
            file.path(OUT_GO, "GO_BP_bubble_FULL.pdf"),
            file.path(OUT_CF, "venn_candidates_FULL.pdf"))) {
  file.copy(f, file.path(FIG_S, basename(f)), overwrite = TRUE)
  cat("copied:", basename(f), "\n")
}
cat("== 54 FULL pipeline done ==\n")
