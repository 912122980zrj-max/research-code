# ============================================================
# 65_v6_figures.R —— v6 主图新面板
# 依据：修改报告_BaP_TBI_最终版.md C-5（Fig3 背景检验、Fig4 敏感性）
# 输出：
#   figures/Figure_03_Discovery/permutation_background.pdf
#   figures/Figure_03_Discovery/sensitivity_heatmap.pdf
#   figures/Figure_04_Enrichment/figure4B_ahr_v6.pdf
#   figures/Figure_03_Discovery/venn_candidates_v6.pdf
# ============================================================

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggvenn)
  library(gridExtra)
})

ROOT <- getwd()
OUT <- file.path(ROOT, "results", "16_Orthology")

# ---------- Fig3F: permutation background ----------
dists <- readRDS(file.path(OUT, "permutation_distributions.rds"))
res <- read_csv(file.path(OUT, "permutation_results.csv"), show_col_types = FALSE)

p_union <- data.frame(overlap = dists$union_unstratified) %>%
  ggplot(aes(overlap)) +
  geom_histogram(binwidth = 1, fill = "#9ECAE1", color = "white") +
  geom_vline(xintercept = 7, color = "#D62728", linewidth = 0.9) +
  annotate("text", x = 6.2, y = Inf, label = "observed 7", hjust = 1,
           vjust = 1.5, color = "#D62728", size = 3.2) +
  labs(title = "Union-module path",
       subtitle = "expected 0.268 | 26.1-fold | P < 0.0001",
       x = "Overlap with 193 key genes", y = "Draws (10,000)") +
  theme_bw(base_size = 10)

p_black <- data.frame(overlap = dists$black_unstratified) %>%
  ggplot(aes(overlap)) +
  geom_histogram(binwidth = 1, fill = "#A1D99B", color = "white") +
  geom_vline(xintercept = 3, color = "#D62728", linewidth = 0.9) +
  annotate("text", x = 2.6, y = Inf, label = "observed 3", hjust = 1,
           vjust = 1.5, color = "#D62728", size = 3.2) +
  labs(title = "Black-module-only path",
       subtitle = "expected 0.0286 | 105.1-fold | P < 0.0001",
       x = "Overlap with 48 black-module key genes", y = "Draws (10,000)") +
  theme_bw(base_size = 10)

pdf(file.path(ROOT, "figures", "Figure_03_Discovery",
              "permutation_background.pdf"), width = 9, height = 4)
grid.arrange(p_union, p_black, ncol = 2)
dev.off()

# ---------- Fig4: sensitivity heatmap ----------
sm <- read_csv(file.path(OUT, "sensitivity_matrix.csv"), show_col_types = FALSE)
cand_sets <- list(
  CLEAN = c("Ccr5","Hsd11b1","Lst1","Mmp2","Npy2r","P2ry6","Ripk1"),
  FULL = c("Cyp1b1","Hsd11b1","Mmp2","P2ry6","Ripk1"),
  GSE58485_only = c("Adam17","Ccr5","Chrm3","Hsd11b1","Lst1","Mmp2",
                    "Npy1r","Npy2r","Nqo2","P2ry6","Ripk1"),
  removeBatchEffect = c("Ccr5","Chrm3","Hsd11b1","Lst1","Mmp2","Npy2r",
                        "P2ry6","Ripk1"),
  SVA = c("Ccr5","Chrm3","Hsd11b1","Lst1","Mmp2","Npy2r","P2ry6","Ripk1"),
  no_ComBat = character(0)
)
all_genes <- unique(unlist(cand_sets))
tier_tab <- read_csv(file.path(OUT, "bap_targets_364_tier.csv"),
                     col_types = cols(.default = col_character()))
tier_map <- setNames(tier_tab$tier, toupper(tier_tab$gene_symbol))
ht <- expand.grid(path = names(cand_sets), gene = all_genes,
                  stringsAsFactors = FALSE)
ht$present <- mapply(function(p, g) g %in% cand_sets[[p]],
                     ht$path, ht$gene)
ht$tier <- ifelse(ht$present, tier_map[toupper(ht$gene)], NA)
ht$tier[is.na(ht$tier)] <- ""
ht$path <- factor(ht$path, levels = rev(names(cand_sets)))
ht$gene <- factor(ht$gene, levels = all_genes)

p_sens <- ggplot(ht, aes(gene, path, fill = tier)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_manual(values = c(A = "#D62728", C = "#9ECAE1"),
                    na.value = "grey92", name = "Tier") +
  labs(x = NULL, y = NULL,
       title = "Candidate sets across sensitivity paths") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(ROOT, "figures", "Figure_03_Discovery",
              "sensitivity_heatmap.pdf"), width = 8, height = 4.5)
print(p_sens)
dev.off()

# ---------- Fig4B v6: AhR coverage status ----------
status <- data.frame(
  gene = c("AHR","CYP1A2","CYP1B1","MAOA","CYP1A1","AHRR","HPRT1 (Hprt)"),
  category = c("Tier A target","Tier A target","Tier A target",
               "Tier A target","AhR-responsive","AhR-responsive",
               "Tier A target (measured)"),
  coverage = c("not covered","not covered","not covered","not covered",
               "not covered","not covered","measured; not a DEG"),
  stringsAsFactors = FALSE
)
status$coverage[status$gene == "HPRT1 (Hprt)"] <- "measured; not a DEG"
tt <- tableGrob(status, rows = NULL, theme = gridExtra::ttheme_minimal(
  base_size = 10))
pdf(file.path(ROOT, "figures", "Figure_04_Enrichment",
              "figure4B_ahr_v6.pdf"), width = 8, height = 3.5)
gridExtra::grid.arrange(
  tt,
  top = grid::textGrob(
    "AhR axis coverage in CLEAN matrix\nFisher test not evaluable (array not covered)",
    gp = grid::gpar(fontsize = 11)))
dev.off()

# ---------- Fig2C v6: Venn with 66 evaluable targets ----------
eval66 <- read_csv(file.path(OUT, "evaluable_targets_v6_CLEAN.csv"),
                   show_col_types = FALSE)$gene
deg <- read.csv(file.path(ROOT, "results", "03_DEG_Analysis",
                          "deg_significant.csv"))$gene
mod <- read.csv(file.path(ROOT, "results", "04_WGCNA",
                          "module_genes.csv"))$gene
venn_list <- list(
  "Evaluable BaP targets (66)" = eval66,
  "CLEAN DEGs (198)" = deg,
  "Module genes (1,244)" = mod
)
pdf(file.path(ROOT, "figures", "Figure_03_Discovery",
              "venn_candidates_v6.pdf"), width = 8, height = 6)
print(ggvenn(venn_list, fill_color = c("#9ECAE1", "#A1D99B", "#FDBE85"),
             stroke_color = "white", text_size = 3.5,
             set_name_size = 4))
dev.off()

cat("v6 figures written\n")
