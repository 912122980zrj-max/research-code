# ============================================================
# 36_supplementary_figures.R —— 补充图 S5/S10/S11 生成
# 运行：Rscript scripts/36_supplementary_figures.R
# 输出：
#   Fig. S5  snRNA-seq QC 小提琴图 + 复核注释 dotplot
#   Fig. S10 星形胶质亚群 marker 热图 + 样本贡献堆叠条形
#   Fig. S11 版本对比（v1/v3/v4 候选漏斗与签名变化）
# ============================================================
suppressMessages({
  library(ggplot2); library(dplyr); library(readr)
  library(Seurat); library(pheatmap); library(patchwork)
})
source("scripts/theme_pub.R")
ROOT <- getwd()
FIG_S <- file.path(ROOT, "figures", "Figure_S_Supplementary")
dir.create(FIG_S, recursive = TRUE, showWarnings = FALSE)

# ---------- Fig. S5: snRNA-seq QC + dotplot ----------
cat("Loading seurat_object_reviewed.rds ...\n")
obj <- readRDS(file.path(ROOT, "results/11_SingleCell/seurat_object_reviewed.rds"))
md <- obj@meta.data
cat("meta cols:", paste(colnames(md), collapse = ", "), "\n")
sample_col <- intersect(c("sample", "orig.ident", "Sample", "sample_id"),
                        colnames(md))[1]
ct_col <- intersect(c("celltype_reviewed", "reviewed_celltype", "celltype", "cell_type",
                      "celltype_reviewed", "seurat_clusters", "RNA_snn_res.1.2"),
                    colnames(md))[1]
cat("sample col:", sample_col, "| celltype col:", ct_col, "\n")
md$sample <- as.character(md[[sample_col]])
md$celltype <- as.character(md[[ct_col]])
md$nFeature <- md$nFeature_RNA
md$nCount <- md$nCount_RNA
if (!"group" %in% colnames(md)) md$group <- "Unknown"
if (!"percent.mt" %in% colnames(md)) {
  mt_genes <- grep("^MT-", rownames(obj), value = TRUE)
  if (length(mt_genes) > 0) {
    md$percent.mt <- Matrix::colSums(GetAssayData(obj, "counts")[mt_genes, ]) /
      Matrix::colSums(GetAssayData(obj, "counts")) * 100
  } else {
    md$percent.mt <- 0
  }
}
obj@meta.data <- md

qc_long <- data.frame(
  sample = md$sample, group = md$group,
  nFeature = md$nFeature, nCount = md$nCount,
  percent_mt = md$percent.mt)
plot_vln <- function(var, ylab, title) {
  p <- ggplot(qc_long, aes(x = sample, y = .data[[var]], fill = group)) +
    geom_violin(scale = "width", width = 0.9) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    scale_fill_manual(values = c(Sham = pal_cat[4], TBI = pal_cat[6]),
                      name = NULL) +
    labs(x = NULL, y = ylab, title = title) +
    theme_pub(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     size = 6),
          legend.position = "top")
  p
}
p_qc1 <- plot_vln("nFeature", "nFeature_RNA", "snRNA-seq QC: nFeature by sample")
p_qc2 <- plot_vln("nCount", "nCount_RNA", "snRNA-seq QC: nCount by sample")
p_qc3 <- plot_vln("percent_mt", "Mitochondrial fraction (%)",
                  "snRNA-seq QC: mitochondrial fraction by sample")
qc_panel <- p_qc1 / p_qc2 / p_qc3 + plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "a")
ggsave_pub(file.path(ROOT, "results/11_SingleCell/figS5_qc_violin.pdf"),
           qc_panel,
           width = 9, height = 10, dpi = 300)
ggsave_pub(file.path(FIG_S, "figS5_qc_violin.pdf"),
           qc_panel,
           width = 9, height = 10, dpi = 300)
writeLines(sort(unique(as.character(md$sample))),
           file.path(ROOT, "results/11_SingleCell/sample_labels.txt"))
cat("Fig. S5 QC violin saved\n")

markers <- read.csv(file.path(ROOT, "results/11_SingleCell/cluster_markers.csv"),
                    check.names = FALSE)
top_feat <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 2) %>%
  pull(gene)
top_feat <- intersect(top_feat, rownames(obj))
if (length(top_feat) > 0) {
  p_dot <- DotPlot(obj, features = unique(top_feat),
                   group.by = ct_col) +
    theme_pub(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          legend.position = "right")
  ggsave_pub(file.path(FIG_S, "figS5_marker_dotplot.pdf"),
             p_dot + theme(axis.title.x = element_text(margin = margin(t = 6)),
                           plot.margin = margin(b = 12, unit = "pt")),
             width = 9, height = 6.0, dpi = 300)
  cat("Fig. S5 dotplot saved\n")
}

# ---------- Fig. S10: 星形胶质亚群诊断 ----------
cat("Loading astro_subset_v2.rds ...\n")
astro <- readRDS(file.path(ROOT,
                           "results/14_Astrocyte_Subclusters/D8_astro/astro_subset_v2.rds"))
am <- astro@meta.data
cluster_col <- intersect(c("seurat_clusters", "ident", "astro_clusters",
                           "RNA_snn_res.0.5"), colnames(am))[1]
if (is.null(cluster_col) || is.na(cluster_col)) cluster_col <- "seurat_clusters"
Idents(astro) <- cluster_col
cat("astro cluster col:", cluster_col, "| levels:",
    paste(levels(Idents(astro)), collapse = ","), "\n")

top3 <- read.csv(
  file.path(ROOT, "results/14_Astrocyte_Subclusters/D8_astro/astro_top10_markers.csv"),
  check.names = FALSE) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 3)
feat3 <- intersect(top3$gene, rownames(astro))
if (length(feat3) > 0) {
  expr_mat <- LayerData(astro, assay = "RNA", layer = "data")
  pdf_pub(file.path(FIG_S, "figS10_marker_heatmap.pdf"),
          width = 9, height = 6)
  pheatmap::pheatmap(
    as.matrix(expr_mat[feat3, ]),
    show_colnames = FALSE, cluster_cols = TRUE, cluster_rows = FALSE,
    color = colorRampPalette(pal_blue)(100),
    main = "Astrocyte subcluster markers (avg top3, log-normalized)",
    fontsize = 8)
  dev.off()
  cat("Fig. S10 marker heatmap saved\n")
}

per_samp <- read.csv(
  file.path(ROOT, "results/14_Astrocyte_Subclusters/D8_astro/astro_per_sample_counts.csv"),
  check.names = FALSE)
clusters <- setdiff(colnames(per_samp), c("sample", "group", "total_astro"))
long_s <- do.call(rbind, lapply(clusters, function(cl) {
  data.frame(sample = per_samp$sample, group = per_samp$group,
             cluster = cl, n = per_samp[[cl]], stringsAsFactors = FALSE)
}))
long_s$prop <- long_s$n / per_samp$total_astro[match(long_s$sample, per_samp$sample)]
p_stack <- ggplot(long_s, aes(x = sample, y = prop, fill = cluster)) +
  geom_col(position = "stack", width = 0.75, colour = "white", linewidth = 0.2) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = colorRampPalette(pal_cat)(length(clusters))) +
  labs(x = NULL, y = "Proportion of astrocyte nuclei",
       title = "Astrocyte subcluster composition by sample") +
  theme_pub(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 6, margin = margin(t = 1)),
        legend.position = "right")
ggsave_pub(file.path(FIG_S, "figS10_sample_composition.pdf"),
           p_stack, width = 8.5, height = 4.5, dpi = 300)
cat("Fig. S10 sample composition saved\n")

# ---------- Fig. S11: 版本对比 ----------
ver <- data.frame(
  version = factor(c("v1 (mixed scale)", "v3 (log2)", "v4 (coverage-filtered)"),
                   levels = c("v1 (mixed scale)", "v3 (log2)",
                              "v4 (coverage-filtered)")),
  DEG = c(1213, 570, 198),
  candidates = c(27, 12, 7),
  module = c("magenta", "yellow", "black"),
  signature = c("P2RY6/CCNA2/HSD11B1", "P2RY6/HSD11B1/NFE2L2",
                "P2RY6/HSD11B1/CCR5"),
  stringsAsFactors = FALSE)
p11 <- ggplot(ver, aes(x = version, y = candidates)) +
  geom_col(width = 0.55, fill = c("#D1E5F0", "#92C5DE", "#4393C3"),
           colour = "grey25", linewidth = 0.3) +
  geom_text(aes(label = paste0("n = ", candidates)),
            vjust = -0.5, size = 4) +
  geom_text(aes(y = 1, label = paste0("DEG ", DEG, "\nmodule ", module,
                                      "\nsignature ", signature)),
            hjust = 0.5, vjust = 1.0, size = 3.0, colour = "grey25",
            lineheight = 1.05) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(x = NULL, y = "BaP–TBI candidate genes",
       title = "Pipeline version comparison (candidate funnel and signature)") +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(size = 9))
ggsave_pub(file.path(FIG_S, "figS11_version_funnel.pdf"),
           p11, width = 7, height = 4.6, dpi = 300)
cat("Fig. S11 version funnel saved\n")
cat("ALL DONE\n")
