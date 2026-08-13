# ============================================================
# 30_figures_2_5d.R —— 补齐 Figure 2（BaP 靶点）与 Figure 5D（LOCO）
# 依据：FINAL_REVIEW_PUBLISHABILITY.md §3.14/P0-6
# 输出：figures/Figure_02_BaP_Targets/Figure2_BaP_targets.pdf
#       figures/Figure_05_ML/Figure5D_LOCO.pdf
# ============================================================

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

ROOT <- getwd()
source(file.path(ROOT, "scripts", "theme_pub.R"))

# ---------- 数据 ----------
targets <- read.delim(file.path(ROOT, "data", "04_BaP_targets",
                                "union_targets_tiered.tsv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
tier_col <- intersect(c("tier", "Tier", "evidence_tier"), colnames(targets))[1]
stopifnot(!is.na(tier_col))
targets$tier <- toupper(trimws(as.character(targets[[tier_col]])))
targets$tier <- ifelse(targets$tier %in% c("A", "B", "C"), targets$tier, "C")
gene_col_t <- intersect(c("gene_symbol", "gene", "Gene"), colnames(targets))[1]
targets$gene <- toupper(trimws(as.character(targets[[gene_col_t]])))
# 别名解析（与 02b target_summary_tiered.json 一致），按官方符号取最高 tier
alias_map <- c(BHLHE76 = "AHR", CYPIA2 = "CYP1A2")
targets$gene <- ifelse(targets$gene %in% names(alias_map),
                       unname(alias_map[targets$gene]), targets$gene)
tier_order <- c(A = 0, B = 1, C = 2)
targets <- targets[order(tier_order[targets$tier]), ]
targets <- targets[!duplicated(targets$gene), ]
bap_set <- toupper(unique(targets$gene))

cand <- read.csv(file.path(ROOT, "results", "05_Candidate_Intersection",
                           "candidate_genes_tiered_v2.csv"),
                 stringsAsFactors = FALSE)
cand$tier <- toupper(cand$tier)

deg <- read.csv(file.path(ROOT, "results", "03_DEG_Analysis",
                          "deg_significant.csv"), stringsAsFactors = FALSE)
mod_genes <- read.csv(file.path(ROOT, "results", "04_WGCNA",
                                "module_genes.csv"), stringsAsFactors = FALSE)
# module_genes.csv 结构可能是 module,gene
mod_col <- intersect(c("gene", "Gene", "genes"), colnames(mod_genes))[1]
mod_set <- toupper(unique(mod_genes[[mod_col]]))
deg_set <- toupper(unique(deg$gene))
# ---------- Figure 2A：靶点分级条形 ----------
tab_tier <- data.frame(tier = c("A", "B", "C"),
                       n = c(sum(targets$tier == "A"),
                             sum(targets$tier == "B"),
                             sum(targets$tier == "C")))
tab_tier$tier <- factor(tab_tier$tier, levels = c("A", "B", "C"))
tab_cand <- data.frame(tier = c("A", "B", "C"),
                       n = c(sum(cand$tier == "A"),
                             sum(cand$tier == "B"),
                             sum(cand$tier == "C")))
tab_cand$tier <- factor(tab_cand$tier, levels = c("A", "B", "C"))

p2a <- ggplot(tab_tier, aes(x = tier, y = n, fill = tier)) +
  geom_col(width = 0.65, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = n), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c(A = "#B2182B", B = "#F4A582", C = "#4393C3"),
                    guide = "none") +
  labs(x = "Evidence tier", y = "Predicted BaP targets (deduplicated)",
       title = "A  BaP target evidence tiers") +
  theme_pub(base_size = 11) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

p2b <- ggplot(tab_cand, aes(x = tier, y = n, fill = tier)) +
  geom_col(width = 0.65, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = n), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c(A = "#B2182B", B = "#F4A582", C = "#4393C3"),
                    guide = "none") +
  labs(x = "Evidence tier", y = "BaP × TBI candidate genes",
       title = paste0("B  Candidate tiers (n=", sum(cand$tier %in%
                         c("A", "B", "C")),
                      "; 3-gene signature all tier C)")) +
  theme_pub(base_size = 11) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)))

# ---------- Figure 2C：三集合交集（BaP / DEG / module） ----------
library(ggvenn)
ven <- list(
  `BaP targets` = bap_set,
  `DEG (CLEAN)` = deg_set,
  `Module genes` = mod_set
)
p2c <- ggvenn(ven, fill_color = c("#2166AC", "#B2182B", "#4393C3"),
              stroke_size = 0.4, text_size = 3.2,
              set_name_size = 4) +
  ggtitle("C  BaP targets ∩ DEG ∩ TBI-associated modules") +
  theme_pub(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, size = 11))

library(patchwork)
p2 <- (p2a | p2b) / p2c + plot_layout(heights = c(1, 1.2))
ggsave_pub(file.path(ROOT, "figures", "Figure_02_BaP_Targets",
                 "Figure2_BaP_targets.pdf"),
       p2, width = 10, height = 8, dpi = 300, bg = "white")

# ---------- Figure 5D：LOCO 汇总 ----------
loco <- read.csv(file.path(ROOT, "results", "07_ML_Signature",
                           "cross_cohort_auc.csv"), stringsAsFactors = FALSE)
loco$label <- paste0(loco$train_cohort, " → ", loco$test_cohort)
loco$comp_label <- paste0(loco$composition, "\n", loco$label)
loco$comp_label <- factor(loco$comp_label,
                          levels = loco$comp_label[order(loco$composition,
                                                         loco$test_auc)])
p5d <- ggplot(loco, aes(x = comp_label, y = test_auc)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
  geom_col(width = 0.6, fill = "#4393C3", colour = "grey30",
           linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", test_auc)), vjust = -0.4, size = 3.6) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Training → test cohort (fixed 5-gene core)",
       y = "Leave-one-cohort-out AUC",
       title = "D  Cross-cohort generalisation (LOCO, raw expression)") +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(size = 7, angle = 20, hjust = 1))
ggsave_pub(file.path(ROOT, "figures", "Figure_05_ML", "Figure5D_LOCO.pdf"),
       p5d, width = 7, height = 5, dpi = 300, bg = "white")

cat("Figure 2 and Figure 5D saved.\n")
cat("target tiers:", paste(names(table(targets$tier)),
                           table(targets$tier), sep = "=", collapse = ", "), "\n")
cat("candidate tiers:", paste(names(table(cand$tier)),
                              table(cand$tier), sep = "=", collapse = ", "), "\n")
cat("DEG:", length(deg_set), " module:", length(mod_set),
    " BaP:", length(bap_set),
    " triple intersection:", length(Reduce(intersect, ven)), "\n")
