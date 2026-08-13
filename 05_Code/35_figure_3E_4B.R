# ============================================================
# 35_figure_3E_4B.R —— Figure 3E（候选基因 log2FC 条形图）
#                      与 Figure 4B（AhR 反应基因状态图）
# 运行：Rscript scripts/35_figure_3E_4B.R
# 输出：figures/Figure_03_Discovery/figure3E_candidate_logFC.pdf
#       results/05_Candidate_Intersection/figure3E_candidate_logFC.pdf
#       figures/Figure_04_Enrichment/figure4B_ahr.pdf
#       results/06_GO_KEGG_Enrichment/figure4B_ahr.pdf
# 说明：2026-08-09 第三轮优化——加宽画布、y/x 轴留白，避免首末柱
#       数值标签贴边与旋转基因名被裁切（H 字母缺角）。
# ============================================================
suppressMessages({library(ggplot2); library(dplyr); library(readr)})
source("scripts/theme_pub.R")
ROOT <- getwd()

# ---------- Figure 3E: 候选基因 log2FC 条形图 ----------
deg <- read.csv(file.path(ROOT, "results/03_DEG_Analysis/deg_all.csv"),
                check.names = FALSE)
cand <- read.csv(file.path(ROOT,
                           "results/05_Candidate_Intersection/candidate_genes.csv"))$gene
cd <- deg[deg$gene %in% cand, c("gene", "logFC", "adj.P.Val")]
cd$gene <- factor(cd$gene, levels = cd$gene[order(cd$logFC)])
cd$dir <- ifelse(cd$logFC > 0, "Up", "Down")

p3e <- ggplot(cd, aes(x = gene, y = logFC, fill = dir)) +
  geom_col(width = 0.62, colour = "grey25", linewidth = 0.25) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_text(aes(label = sprintf("%.2f", logFC)),
            vjust = ifelse(cd$logFC > 0, -0.6, 1.6), size = 3.2) +
  scale_fill_manual(values = c(Up = pal_rdbu[6], Down = pal_rdbu[1]),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.20))) +
  scale_x_discrete(expand = expansion(mult = c(0.12, 0.12))) +
  labs(x = NULL, y = expression(log[2]~fold~change~"("*TBI~vs~Sham*")"),
       title = "E  Seven candidate genes (all tier C)") +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1,
                                   size = 9),
        plot.margin = margin(8, 8, 6, 8))
ggsave_pub(file.path(ROOT, "figures/Figure_03_Discovery",
                     "figure3E_candidate_logFC.pdf"),
           p3e, width = 6.2, height = 4.2, dpi = 300)
ggsave_pub(file.path(ROOT, "results/05_Candidate_Intersection",
                     "figure3E_candidate_logFC.pdf"),
           p3e, width = 6.2, height = 4.2, dpi = 300)
cat("Figure 3E saved\n")

# ---------- Figure 4B: AhR 反应基因可评估性 ----------
ahr <- data.frame(
  gene = c("CYP1A1", "CYP1B1", "AHRR"),
  in_candidates = c(FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)
bg <- unique(toupper(deg$gene))
ahr$in_background <- toupper(ahr$gene) %in% bg
ahr$label <- ifelse(ahr$in_candidates, "candidate",
                    ifelse(ahr$in_background, "background only",
                           "not measured\n(after coverage filter)"))
ahr$fill <- ifelse(ahr$in_candidates, "#B2182B",
                   ifelse(ahr$in_background, "#92C5DE", "#F4A582"))
ahr$gene <- factor(ahr$gene, levels = rev(c("CYP1A1", "CYP1B1", "AHRR")))
p4b <- ggplot(ahr, aes(x = 1, y = gene, fill = fill)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 3.2, colour = "white",
            fontface = "plain", lineheight = 0.9) +
  scale_fill_identity() +
  labs(title = "B  AhR-responsive genes: not measurable after coverage filter\n(0/3 hits; Fisher P = 1.00)",
       x = NULL, y = NULL) +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_text(size = 9, face = "plain", hjust = 0.5))
ggsave_pub(file.path(ROOT, "figures/Figure_04_Enrichment/figure4B_ahr.pdf"),
           p4b, width = 6.0, height = 2.8, dpi = 300)
ggsave_pub(file.path(ROOT, "results/06_GO_KEGG_Enrichment/figure4B_ahr.pdf"),
           p4b, width = 6.0, height = 2.8, dpi = 300)
cat("Figure 4B saved\n")
