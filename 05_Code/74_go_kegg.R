# ============================================================
# 74_go_kegg_v7.R —— v7 主分析 GO（BP/CC/MF）+ KEGG（11 候选）
# 输入：results/16_Orthology/v7_single_cohort/candidates_GSE58485_only.csv
# 输出：results/06_GO_KEGG_Enrichment/v7_single_cohort/
# ============================================================
suppressMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

ROOT <- getwd()
OUT <- file.path(ROOT, "results", "06_GO_KEGG_Enrichment", "v7_single_cohort")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts", "theme_pub.R"))

cand <- read.csv(file.path(ROOT, "results", "16_Orthology", "v7_single_cohort",
                           "candidates_GSE58485_only.csv"))$gene
cat("candidates:", length(cand), "\n")

eg <- bitr(cand, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
cat("mapped:", nrow(eg), "\n")

ego <- enrichGO(gene = eg$ENTREZID, OrgDb = org.Mm.eg.db, ont = "ALL",
                pAdjustMethod = "BH", pvalueCutoff = 0.05,
                qvalueCutoff = 0.2, readable = TRUE)
if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  godf <- as.data.frame(ego)
  godf$n_genes <- godf$Count
  write_csv(godf, file.path(OUT, "GO_ALL_candidates.csv"))
  godf_keep <- godf[godf$Count >= 3, , drop = FALSE]
  write_csv(godf_keep, file.path(OUT, "GO_ALL_keep.csv"))
  cat("GO terms:", nrow(godf), " keep(Count>=3):", nrow(godf_keep), "\n")
  print(godf_keep[, c("ONTOLOGY", "Description", "Count", "p.adjust", "geneID")])
} else {
  cat("GO: no terms\n")
  writeLines("no GO terms", file.path(OUT, "GO_ALL_keep.csv"))
}

kk <- enrichKEGG(gene = eg$ENTREZID, organism = "mmu",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2)
if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
  kdf <- as.data.frame(kk)
  write_csv(kdf, file.path(OUT, "KEGG_candidates.csv"))
  cat("KEGG terms:", nrow(kdf), "\n")
  print(kdf[, c("Description", "Count", "p.adjust", "geneID")])
} else {
  cat("KEGG: no terms\n")
  writeLines("no KEGG terms", file.path(OUT, "KEGG_candidates.csv"))
}

writeLines(c(
  "v7 GO/KEGG (GSE58485-only 11 candidates)",
  paste0("generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("candidates: ", paste(cand, collapse = ",")),
  paste0("mapped ENTREZ: ", nrow(eg))
), file.path(OUT, "V7_GO_KEGG_REPORT.txt"))
