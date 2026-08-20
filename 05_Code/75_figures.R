# ============================================================
# 75_figures_v7.R —— v7 主图补充：置换背景 + 敏感性热图 + GO 气泡
# 输出：figures/Figure_03_Discovery/permutation_background_v7.pdf、
#       sensitivity_heatmap_v7.pdf；figures/Figure_04_Enrichment/GO_BP_bubble_v7.pdf
# ============================================================
suppressMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
})

ROOT <- getwd()
FIG3 <- file.path(ROOT, "figures", "Figure_03_Discovery")
FIG4 <- file.path(ROOT, "figures", "Figure_04_Enrichment")
source(file.path(ROOT, "scripts", "theme_pub.R"))

# ---------- 置换分布 ----------
perm <- readRDS(file.path(ROOT, "results", "16_Orthology", "v7_single_cohort",
                          "permutation_v7.rds"))
res <- perm$results
null_union <- attr(perm$results, "null")  # 可能为空；从 rds 重建
if (is.null(null_union)) {
  # 快速重建 union null（seed 与 71 一致）
  set.seed(2026)
  universe <- perm$universe
  key_upper <- toupper(read.csv(file.path(ROOT, "results", "16_Orthology",
                                          "v7_single_cohort",
                                          "key_genes_GSE58485_only.csv"))$gene)
  null_union <- sapply(1:10000, function(i)
    sum(toupper(sample(universe, 11)) %in% key_upper))
}
null_df <- data.frame(overlap = null_union)
p_perm <- ggplot(null_df, aes(overlap)) +
  geom_histogram(binwidth = 1, fill = "grey75", color = "grey30") +
  geom_vline(xintercept = 11, color = "#D73027", linewidth = 0.9) +
  annotate("text", x = 11, y = Inf, hjust = -0.1, vjust = 1.5,
           label = "observed = 11\n(expected 0.79; 13.9-fold; P<0.0001)",
           color = "#D73027", size = 3.2) +
  scale_x_continuous(breaks = 0:11, limits = c(-0.5, 12)) +
  theme_pub() +
  labs(x = "Overlap count under measurability-matched null (10,000 draws)",
       y = "Draws",
       title = "GSE58485-only: 11 candidates vs 364 key genes")
ggsave(file.path(FIG3, "permutation_background.pdf"), p_perm,
       width = 6.5, height = 4.5)

# ---------- 敏感性热图（v7 主行 = GSE58485-only） ----------
paths <- list(
  `GSE58485-only (main)` = c("Adam17","Ccr5","Chrm3","Hsd11b1","Lst1","Mmp2",
                             "Npy1r","Npy2r","Nqo2","P2ry6","Ripk1"),
  `Merged CLEAN (stress test)` = c("Ccr5","Hsd11b1","Lst1","Mmp2","Npy2r",
                                   "P2ry6","Ripk1"),
  `FULL (incl. drug-treated)` = c("Cyp1b1","Hsd11b1","Mmp2","P2ry6","Ripk1"),
  `removeBatchEffect` = c("Ccr5","Chrm3","Hsd11b1","Lst1","Mmp2","Npy2r",
                          "P2ry6","Ripk1"),
  `SVA` = c("Ccr5","Chrm3","Hsd11b1","Lst1","Mmp2","Npy2r","P2ry6","Ripk1"),
  `No ComBat` = character(0)
)
all_genes <- unique(unlist(paths))
heat <- do.call(rbind, lapply(names(paths), function(p) {
  data.frame(path = p, gene = all_genes,
             present = all_genes %in% paths[[p]],
             stringsAsFactors = FALSE)
}))
heat$path <- factor(heat$path, levels = names(paths))
heat$gene <- factor(heat$gene, levels = rev(all_genes))
p_heat <- ggplot(heat, aes(gene, path, fill = present)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#2166AC")) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "none") +
  labs(x = NULL, y = NULL,
       title = "Candidate sets across preprocessing/cohort paths")
ggsave(file.path(FIG3, "sensitivity_heatmap.pdf"), p_heat,
       width = 8, height = 3.5)

# ---------- GO BP 气泡（top 10） ----------
go <- read_csv(file.path(ROOT, "results", "06_GO_KEGG_Enrichment",
                         "v7_single_cohort", "GO_ALL_keep.csv"))
go_bp <- go %>% filter(ONTOLOGY == "BP") %>% arrange(p.adjust) %>% head(10)
p_go <- ggplot(go_bp, aes(reorder(Description, p.adjust), -log10(p.adjust),
                          size = Count)) +
  geom_point(color = "#2166AC", alpha = 0.8) +
  coord_flip() +
  theme_pub() +
  labs(x = NULL, y = expression(-log[10] ~ adjusted ~ P),
       size = "Gene count",
       title = "GO biological process enrichment (11 candidates, top 10)")
ggsave(file.path(FIG4, "GO_BP_bubble.pdf"), p_go,
       width = 7.5, height = 4.8)

writeLines("v7 figures generated", file.path(FIG3, "v7_figures.log"))
