# ============================================================
# 33_workflow_volcano.R —— Figure 1 流程图 + DEG 火山图（发表级）
# 运行：Rscript scripts/33_workflow_volcano.R
# 说明：与 theme_pub.R 色系一致（Paul Tol 分类色 / RdBu / 蓝系渐变）
# 输出：
#   figures/Figure_01_Workflow/Figure1_workflow.pdf (+pdf)
#   results/03_DEG_Analysis/volcano_deg.pdf (+pdf)
#   figures/Figure_03_Discovery/volcano_deg.pdf（副本，供图集）
# ============================================================
suppressMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(readr)
})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else {
  stop("请在项目根目录运行：E:\\sheng xin\\tbi")
}
source(file.path(ROOT, "scripts", "theme_pub.R"))

# ---------- 1) Figure 1 研究流程图 ----------
fig1_dir <- file.path(ROOT, "figures", "Figure_01_Workflow")
dir.create(fig1_dir, recursive = TRUE, showWarnings = FALSE)

# v4 动态数字（避免硬编码旧口径）
deg_v4 <- nrow(read.csv(file.path(ROOT, "results", "03_DEG_Analysis",
                                  "deg_significant.csv")))
cand_v4 <- nrow(read.csv(file.path(ROOT, "results",
                                   "05_Candidate_Intersection",
                                   "candidate_genes.csv")))
mt_v4 <- read.csv(file.path(ROOT, "results", "04_WGCNA",
                            "module_trait.csv"))
best_mod_v4 <- sub("^ME", "", mt_v4$module[which.max(abs(mt_v4$r))])

# 画布 0-100 x 0-100；box：x 中心、y 中心、宽、高
boxes <- data.frame(
  id = c("B1", "B2", "B3", "B4", "B5", "B6", "B7"),
  x  = c(50, 50, 50, 50, 50, 50, 82),
  y  = c(93, 78, 63, 48, 33, 18, 33),
  w  = c(64, 64, 64, 64, 64, 64, 28),
  h  = c( 9,  9,  9,  9,  9,  9, 10),
  fill = c(pal_cat[4], pal_cat[5], pal_cat[2], pal_cat[1], pal_cat[6], pal_cat[7], "white"),
  stringsAsFactors = FALSE
)

box_text <- c(
  "GEO datasets",
  "BaP target prediction & evidence tiering",
  "Bulk discovery  (CLEAN, n = 33)",
  "Machine-learning signature",
  "Dry-experiment modules",
  "Hypothesis-generating candidates",
  "Wet-lab validation (planned)"
)
box_sub <- c(
  "GSE58485 · GSE41345 · GSE128543 · GSE104687 · GSE209552",
  "ChEMBL / SEA / PharmMapper / SwissTargetPrediction  ->  364 genes (Tier A/B/C)",
  paste0("ComBat -> DEG (", deg_v4, ") -> WGCNA (", best_mod_v4,
         ") -> ", cand_v4, " candidates"),
  "P2RY6 / HSD11B1 / CCR5  ·  leave-one-cohort-out  ·  frozen external validation",
  "Deconvolution · CellChat · scTenifoldKnk · Astrocyte states · Docking",
  "BaP-target x TBI transcriptomic intersection  (exploratory)",
  "qPCR (CCI vs Sham)  ->  BaP x CCI"
)

arrows <- data.frame(
  x = c(50, 50, 50, 50, 50, 82, 50),
  y = c(88.5, 73.5, 58.5, 43.5, 28.5, 23, 33),
  xend = c(50, 50, 50, 50, 50, 82, 71),
  yend = c(82.5, 67.5, 52.5, 37.5, 22.5, 18, 33),
  lty = c(1, 1, 1, 1, 1, 2, 2)
)

lab_df <- data.frame(
  x = c(50, 50, 50, 50, 50, 50, 82),
  y = boxes$y + 2.2,
  label = box_text,
  sub = box_sub,
  size = c(4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.0),
  sub_size = c(2.9, 2.9, 2.9, 2.9, 2.9, 2.9, 2.7),
  stringsAsFactors = FALSE
)

p_flow <- ggplot() +
  # 连接箭头
  geom_segment(data = arrows,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.14, "inches"), type = "closed"),
               linewidth = 0.55, colour = "grey35", linetype = arrows$lty) +
  # 方框
  geom_rect(data = boxes,
            aes(xmin = x - w / 2, xmax = x + w / 2,
                ymin = y - h / 2, ymax = y + h / 2),
            fill = boxes$fill, colour = "grey25", linewidth = 0.5,
            alpha = 0.92) +
  # 主标题文字
  geom_text(data = lab_df, aes(x = x, y = y + 1.6, label = label),
            size = lab_df$size, fontface = "bold", colour = "grey15") +
  # 副标题文字
  geom_text(data = lab_df, aes(x = x, y = y - 2.4, label = sub),
            size = lab_df$sub_size, colour = "grey25") +
  coord_cartesian(xlim = c(0, 100), ylim = c(4, 100), clip = "off") +
  labs(title = "Study design and analytical workflow") +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14,
                              margin = margin(b = 6)),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave_pub(file.path(fig1_dir, "Figure1_workflow.pdf"), p_flow,
       width = 9, height = 10)
cat("Figure 1 workflow saved\n")

# ---------- 2) DEG 火山图 ----------
deg <- read.csv(file.path(ROOT, "results", "03_DEG_Analysis", "deg_all.csv"),
                check.names = FALSE)
stopifnot(all(c("gene", "logFC", "adj.P.Val") %in% colnames(deg)))

deg <- deg %>%
  mutate(
    neg_log10p = -log10(pmax(adj.P.Val, 1e-300)),
    status = case_when(
      adj.P.Val < 0.05 & logFC >  0.585 ~ "Up",
      adj.P.Val < 0.05 & logFC < -0.585 ~ "Down",
      TRUE ~ "NS"
    ),
    is_candidate = toupper(gene) %in% toupper(
      read.csv(file.path(ROOT, "results", "05_Candidate_Intersection",
                         "candidate_genes.csv"))$gene)
  )

up_n <- sum(deg$status == "Up"); down_n <- sum(deg$status == "Down")
cat("Up:", up_n, "| Down:", down_n, "| NS:", sum(deg$status == "NS"), "\n")

# 标注基因：显著性 top 12 + 强制包含签名/关注基因（若显著）
top_lab <- deg %>%
  filter(status != "NS") %>%
  arrange(desc(neg_log10p)) %>%
  slice_head(n = 12) %>%
  pull(gene)
focus <- c("P2RY6", "HSD11B1", "CCR5", "LST1", "NPY2R",
           "MMP2", "RIPK1", "CYP1B1")
focus_lab <- deg$gene[deg$gene %in% focus & deg$status != "NS"]
lab_genes <- unique(c(top_lab, focus_lab))
deg$label <- ifelse(deg$gene %in% lab_genes, deg$gene, NA_character_)

status_col <- c(Up = pal_rdbu[6], Down = pal_rdbu[1], NS = "grey75")

p_volcano <- ggplot(deg, aes(x = logFC, y = neg_log10p)) +
  geom_point(aes(colour = status), size = 0.9, alpha = 0.75) +
  geom_point(data = deg[deg$is_candidate, ],
             aes(x = logFC, y = neg_log10p),
             colour = "#B26A00", fill = "#B26A00", shape = 21,
             size = 2.0, alpha = 0.95, stroke = 0.25) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = 2, colour = "grey55",
             linewidth = 0.35) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey55",
             linewidth = 0.35) +
  scale_colour_manual(
    values = status_col,
    labels = c(Up = paste0("Up (", up_n, ")"),
               Down = paste0("Down (", down_n, ")"),
               NS = "Not significant"),
    name = NULL
  ) +
  geom_text_repel(
    data = deg[!is.na(deg$label), ],
    aes(label = label),
    size = 2.8, max.overlaps = 25, seed = 2026,
    segment.size = 0.25, segment.colour = "grey50",
    box.padding = 0.35, min.segment.length = 0.2
  ) +
  labs(
    x = expression(log[2]~fold~change~"("*TBI~vs~Sham*")"),
    y = expression(-log[10]~"("*adjusted~P*")"),
    title = "Differential expression (CLEAN, n = 33)"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = c(0.12, 0.88),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key.height = unit(0.45, "lines"))

out_deg <- file.path(ROOT, "results", "03_DEG_Analysis")
ggsave_pub(file.path(out_deg, "volcano_deg.pdf"), p_volcano,
       width = 7, height = 6)
fig3_dir <- file.path(ROOT, "figures", "Figure_03_Discovery")
dir.create(fig3_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(file.path(out_deg, "volcano_deg.pdf"),
          file.path(fig3_dir, "volcano_deg.pdf"), overwrite = TRUE)
cat("DEG volcano saved\n")

message("完成：Figure1_workflow 与 volcano_deg 已输出。")
