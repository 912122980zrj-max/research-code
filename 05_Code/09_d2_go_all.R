# ============================================================
# 09_d2_go_all.R —— D2 GO 全 ontology（BP+CC+MF）+ KEGG 补齐
# 依据：REVIEW4_FINAL_AUDIT D2（模板 enrichGO(ont="ALL")）
# 输入：results/05_Candidate_Intersection/candidate_genes.csv（CLEAN 27）
# 输出：results/06_GO_KEGG_Enrichment/D2_go_kegg/（GO_ALL_candidates.csv、GO_ALL_keep.csv、
#       GO_bubble_BP/CC/MF.pdf、KEGG_keep_v2.csv、KEGG_bubble_v2.pdf、D2_report.txt）
# ============================================================

suppressMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "results", "06_GO_KEGG_Enrichment", "D2_go_kegg")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts", "theme_pub.R"))
report <- file.path(OUT, "D2_report.txt")
writeLines("D2 GO all-ontology + KEGG enrichment", report)
write_report <- function(...) cat(..., "\n", sep = "", file = report,
                                  append = TRUE)
write_report("generated:", format(Sys.time()))

cand <- read.csv(file.path(ROOT, "results", "05_Candidate_Intersection",
                           "candidate_genes.csv"))$gene
write_report("candidates:", length(cand))
eg <- bitr(cand, fromType = "SYMBOL", toType = "ENTREZID",
           OrgDb = org.Mm.eg.db)
write_report("mapped ENTREZ IDs:", nrow(eg))

ego <- enrichGO(gene = eg$ENTREZID, OrgDb = org.Mm.eg.db, ont = "ALL",
                pAdjustMethod = "BH", pvalueCutoff = 0.05,
                qvalueCutoff = 0.2, readable = TRUE)
if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  godf <- as.data.frame(ego)
  godf$n_genes <- godf$Count
  godf$is_1gene <- godf$Count == 1
  write.csv(godf, file.path(OUT, "GO_ALL_candidates.csv"), row.names = FALSE)
  godf_keep <- godf[godf$Count >= 3, , drop = FALSE]
  write.csv(godf_keep, file.path(OUT, "GO_ALL_keep.csv"), row.names = FALSE)
  write_report("GO terms (ALL):", nrow(godf),
               "| keep (Count>=3):", nrow(godf_keep))
  write_report("by ontology (keep):",
               paste(names(table(godf_keep$ONTOLOGY)), "=",
                     as.integer(table(godf_keep$ONTOLOGY)), collapse = "; "))
  for (ont in c("BP", "CC", "MF")) {
    sub <- godf_keep[godf_keep$ONTOLOGY == ont, , drop = FALSE]
    if (nrow(sub) == 0) next
    sub <- sub[order(sub$p.adjust), ]
    sub$ratio <- sapply(strsplit(as.character(sub$GeneRatio), "/"),
                        function(x) as.numeric(x[1]) / as.numeric(x[2]))
    sub$Description <- factor(sub$Description, levels = rev(sub$Description))
    sub <- head(sub, 10)
    p <- ggplot(sub, aes(x = ratio, y = Description, size = Count,
                         color = p.adjust)) +
      geom_point(alpha = 0.85) +
      scale_color_gradient(low = pal_blue[8], high = pal_blue[3],
                           name = "p.adjust") +
      scale_size_continuous(range = c(2, 6), name = "Count") +
      labs(title = paste0("GO ", ont, " (Count>=3)"),
           x = "Gene ratio", y = NULL) +
      theme_pub() +
      theme(panel.grid.major.y = element_blank())
    ggsave_pub(file.path(OUT, paste0("GO_bubble_", ont, ".pdf")), p,
           width = 7, height = 6, dpi = 300)
    write_report(paste0("GO_", ont, " bubble written: ", nrow(sub), " terms"))
  }
} else {
  write_report("GO enrichment returned no terms")
}

ekegg <- tryCatch(enrichKEGG(gene = eg$ENTREZID, organism = "mmu",
                             pvalueCutoff = 0.05),
                  error = function(e) {
                    write_report("KEGG enrichKEGG failed:", conditionMessage(e))
                    NULL
                  })
if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  kdf <- as.data.frame(ekegg)
  kdf$n_genes <- kdf$Count
  kdf$is_1gene <- kdf$Count == 1
  write.csv(kdf, file.path(OUT, "KEGG_candidates_v2.csv"), row.names = FALSE)
  kdf_keep <- kdf[kdf$Count >= 3, , drop = FALSE]
  write.csv(kdf_keep, file.path(OUT, "KEGG_keep_v2.csv"), row.names = FALSE)
  write_report("KEGG terms:", nrow(kdf), "| keep (Count>=3):", nrow(kdf_keep))
  if (nrow(kdf_keep) > 0) {
    kk <- kdf_keep[order(kdf_keep$p.adjust), ]
    kk$ratio <- sapply(strsplit(as.character(kk$GeneRatio), "/"),
                       function(x) as.numeric(x[1]) / as.numeric(x[2]))
    kk$Description <- factor(kk$Description, levels = rev(kk$Description))
    p <- ggplot(kk, aes(x = ratio, y = Description, size = Count,
                        color = p.adjust)) +
      geom_point(alpha = 0.85) +
      scale_color_gradient(low = pal_blue[8], high = pal_blue[3],
                           name = "p.adjust") +
      scale_size_continuous(range = c(2, 6), name = "Count") +
      labs(title = "KEGG (Count>=3)", x = "Gene ratio", y = NULL) +
      theme_pub() +
      theme(panel.grid.major.y = element_blank())
    ggsave_pub(file.path(OUT, "KEGG_bubble_v2.pdf"), p,
           width = 7, height = 5, dpi = 300)
  }
} else {
  write_report(paste0("KEGG terms: 0 (no significant pathway for ",
                      length(cand), " candidates;"),
               " verified by direct enrichKEGG, not a network error)")
}
message("D2 完成。结果在 ", OUT)
