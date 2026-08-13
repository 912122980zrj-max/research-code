# ============================================================
# 07_publication_figures.R —— 发表级图形统一重绘
# 运行：Rscript scripts/07_publication_figures.R
# 说明（第二轮审查后）：
#   - 从现有结果文件重绘全部主图（不重跑分析、不依赖 KEGG 网络）；
#   - 统一使用 scripts/theme_pub.R 的主题与低饱和度配色；
#   - 覆盖主图：pca_before/after、wgcna_dendro、venn_candidates、
#     GO/KEGG_bubble、linear_attribution/beeswarm、umap_celltype(_reviewed)、
#     roe_heatmap、celltype_proportions、hub_featureplots、hub_module_score、
#     frozen_roc。
# ============================================================

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(WGCNA)
  library(Seurat)
  library(pheatmap)
  library(ggvenn)
  library(glmnet)
  library(caret)
  library(pROC)
})

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
OUT4 <- file.path(ROOT, "results", "07_ML_Signature")
OUT6 <- file.path(ROOT, "results", "11_SingleCell")
OUT8 <- file.path(ROOT, "results", "08_External_Validation")
source(file.path(ROOT, "scripts", "theme_pub.R"))
for (od in c(OUT_PREP, OUT_DEG, OUT_WGCNA, OUT_CAND, OUT_GO)) dir.create(od, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT4, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT6, recursive = TRUE, showWarnings = FALSE)

# 第三轮审查 §3.1 修复：PCA before = 原始合并矩阵（ComBat 前，与 03 预处理一致），
# PCA after = 直接读 03 保存的 corrected_expr*.csv（不重复校正）。
read_bulk_rds <- function(gse) {
  readRDS(file.path(ROOT, "data", "01_GEO_bulk", paste0(gse, ".rds")))
}
infer_group_titles <- function(titles) {
  t <- tolower(titles)
  is_ctrl <- grepl("\\b(sham|control|uninjured|vehicle|naive)\\b", t)
  is_tbi <- grepl("\\b(tbi|mtbi|cci|injured|injury|injur|impact|fp)\\b", t) |
    grepl("fluid percussion", t)
  ifelse(is_tbi & !is_ctrl, "TBI",
         ifelse(is_ctrl & !is_tbi, "Sham", "Other"))
}
build_raw_discovery <- function(clean = TRUE) {
  gses <- c("GSE58485", "GSE41345")
  dat <- lapply(gses, read_bulk_rds)
  names(dat) <- gses
  common <- Reduce(intersect, lapply(dat, function(d) rownames(d$expr)))
  expr_list <- lapply(gses, function(g) dat[[g]]$expr[common, , drop = FALSE])
  ex <- do.call(cbind, expr_list)
  ph <- do.call(rbind, lapply(gses, function(g) {
    p <- dat[[g]]$pheno
    data.frame(title = as.character(p$title),
               group = infer_group_titles(p$title),
               cohort = g, stringsAsFactors = FALSE)
  }))
  ph$group <- tolower(ph$group)
  keep <- ph$group %in% c("sham", "tbi")
  if (clean) {
    keep <- keep & !grepl("cyclophosphamide|ex-4-mtbi", ph$title,
                          ignore.case = TRUE)
  }
  ex <- ex[, keep, drop = FALSE]
  ph <- ph[keep, , drop = FALSE]
  keep_gene <- rowSums(is.na(ex)) < 0.2 * ncol(ex)
  ex <- ex[keep_gene, ]
  ex <- apply(ex, 2, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    x
  })
  list(expr = ex, pheno = ph)
}

cat("== 1/10 PCA（CLEAN + FULL）==\n")
draw_pca <- function(suffix = "") {
  raw <- build_raw_discovery(clean = (suffix == ""))
  cor_ex <- as.matrix(read.csv(file.path(OUT_PREP, paste0("corrected_expr",
                                                      suffix, ".csv")),
                               row.names = 1, check.names = FALSE))
  cor_ph <- read.csv(file.path(OUT_PREP, paste0("corrected_pheno", suffix,
                                            ".csv")),
                     check.names = FALSE, fileEncoding = "UTF-8")
  panels <- list(before = raw$expr, after = cor_ex)
  ph_use <- list(before = raw$pheno, after = cor_ph)
  for (tag in names(panels)) {
    mat <- panels[[tag]]
    ph <- ph_use[[tag]]
    pc <- prcomp(t(mat), scale. = TRUE)
    pct <- round(100 * summary(pc)$importance[2, 1:2], 1)
    df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                     cohort = factor(ph$cohort),
                     group = factor(ph$group, levels = c("sham", "tbi")))
    p <- ggplot(df, aes(x = PC1, y = PC2, color = cohort, shape = group)) +
      geom_point(size = 2.2, alpha = 0.9) +
      scale_color_manual(values = pal_cat[1:2], name = "Cohort") +
      scale_shape_manual(values = c(16, 17), name = "Group",
                         labels = c("Sham", "TBI")) +
      labs(title = paste0("PCA ", tag, " ComBat", suffix),
           x = paste0("PC1 (", pct[1], "%)"),
           y = paste0("PC2 (", pct[2], "%)")) +
      theme_pub()
    ggsave_pub(file.path(OUT_PREP, paste0("pca_", tag, suffix, ".pdf")), p,
           width = 6, height = 5, dpi = 300)
  }
}
draw_pca("")
  if (file.exists(file.path(OUT_PREP, "corrected_expr_FULL.csv"))) draw_pca("_FULL")

cat("== 2/10 WGCNA 树状图（CLEAN，重算 TOM）==\n")
ex <- as.matrix(read.csv(file.path(OUT_PREP, "corrected_expr.csv"),
                         row.names = 1, check.names = FALSE))
mad <- apply(ex, 1, mad)
top_genes <- names(sort(mad, decreasing = TRUE))[seq_len(min(8000, length(mad)))]
gv <- apply(ex[top_genes, ], 1, var, na.rm = TRUE)
top_genes <- top_genes[gv > 0]
datExpr <- t(ex[top_genes, ])
rep_lines <- readLines(file.path(OUT_PREP, "03_bulk_report.txt"),
                       warn = FALSE, encoding = "UTF-8")
softPower_v4 <- as.numeric(sub(".*WGCNA softPower:", "",
                               rep_lines[grep("WGCNA softPower:", rep_lines)]))
if (is.na(softPower_v4)) softPower_v4 <- 12
adj <- adjacency(datExpr, power = softPower_v4)
tom <- TOMsimilarity(adj, verbose = 0)
geneTree <- hclust(as.dist(1 - tom), method = "average")
mc <- read.csv(file.path(OUT_WGCNA, "module_colors.csv"))
colors_v <- mc$module[match(colnames(datExpr), mc$gene)]
if (any(is.na(colors_v))) colors_v[is.na(colors_v)] <- "grey"
draw_wgcna_dendro(
  geneTree, colors_v,
  main = "WGCNA gene dendrogram and module colors (CLEAN)",
  out_file = file.path(OUT_WGCNA, "wgcna_dendro.pdf"))

cat("== 3/10 Venn（CLEAN）==\n")
bap <- read.delim(file.path(ROOT, "data", "04_BaP_targets", "union_targets.tsv"))
deg <- read.csv(file.path(OUT_DEG, "deg_significant.csv"))
mod <- read.csv(file.path(OUT_WGCNA, "module_genes.csv"))
venn_list <- list(BaP = unique(toupper(bap$gene_symbol)),
                  DEG = unique(toupper(deg$gene)),
                  Module = unique(toupper(mod$gene)))
p_venn <- ggvenn(venn_list, show_elements = FALSE,
                 fill_color = c("#4393C3", "#B2182B", "#117733"),
                 stroke_size = 0.4, set_name_size = 4.5)
ggsave_pub(file.path(OUT_CAND, "venn_candidates.pdf"), p_venn,
       width = 6, height = 5.5, dpi = 300)

cat("== 4/10 GO/KEGG 气泡图 ==\n")
draw_dot <- function(csv_file, out_file, title) {
  if (!file.exists(csv_file)) {
    message("  [skip] ", csv_file)
    return(invisible(NULL))
  }
  d <- read.csv(csv_file)
  d$ratio <- sapply(strsplit(as.character(d$GeneRatio), "/"),
                    function(x) as.numeric(x[1]) / as.numeric(x[2]))
  d <- d[order(d$p.adjust), ]
  d$Description <- factor(d$Description,
                          levels = rev(d$Description))
  d <- head(d, 15)
  p <- ggplot(d, aes(x = ratio, y = Description, size = Count,
                     color = p.adjust)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "#084594", high = "#6BAED6",
                         name = "p.adjust") +
    scale_size_continuous(range = c(2.5, 6.5), name = "Count") +
    labs(title = title, x = "Gene ratio", y = NULL) +
    theme_pub() +
    theme(panel.grid.major.y = element_blank())
  ggsave_pub(out_file, p, width = 7, height = 6, dpi = 300)
}
draw_dot(file.path(OUT_GO, "GO_BP_keep.csv"),
         file.path(OUT_GO, "GO_BP_bubble.pdf"), "GO biological process (Count>=3)")
draw_dot(file.path(OUT_GO, "KEGG_keep.csv"),
         file.path(OUT_GO, "KEGG_bubble.pdf"), "KEGG pathway (Count>=3)")

cat("== 5/10 线性归因 ==\n")
cand <- toupper(read.csv(file.path(OUT_CAND, "candidate_genes.csv"))$gene)
ex <- as.matrix(read.csv(file.path(OUT_PREP, "corrected_expr.csv"),
                         row.names = 1, check.names = FALSE))
rn <- toupper(rownames(ex))
ex2 <- ex[rn %in% cand, , drop = FALSE]
rn2 <- rn[rn %in% cand]
if (anyDuplicated(rn2)) ex2 <- rowsum(ex2, group = rn2) else rownames(ex2) <- rn2
feat <- sort(rownames(ex2))
ex2 <- ex2[feat, , drop = FALSE]
coef_df <- read.csv(file.path(OUT4, "D3_ml_benchmark", "enet_coefficients.csv"))
beta <- rep(0, length(feat)); names(beta) <- feat
hit <- intersect(coef_df$gene, feat)
beta[hit] <- coef_df$coefficient[match(hit, coef_df$gene)]
Xs <- scale(t(ex2))
attr_mat <- sweep(Xs, 2, beta, "*")
mean_abs <- sort(colMeans(abs(attr_mat)), decreasing = TRUE)
imp_df <- data.frame(gene = factor(names(mean_abs),
                                   levels = rev(names(mean_abs))),
                     value = as.numeric(mean_abs))
p1 <- ggplot(imp_df, aes(x = gene, y = value)) +
  geom_col(fill = "#4A6FA5", width = 0.7) + coord_flip() +
  labs(title = "Linear attribution importance (mean |beta * scaled X|)",
       y = "mean |attribution|") +
  theme_pub(base_size = 11)
ggsave_pub(file.path(OUT4, "linear_attribution.pdf"), p1,
       width = 6, height = 5, dpi = 300)
long <- data.frame(gene = rep(colnames(attr_mat), each = nrow(attr_mat)),
                   attribution = as.vector(attr_mat),
                   expr = as.vector(Xs))
long$gene <- factor(long$gene, levels = names(mean_abs))
p2 <- ggplot(long, aes(x = attribution, y = gene, color = expr)) +
  geom_jitter(size = 0.6, alpha = 0.7, width = 0, height = 0.15) +
  scale_color_gradient2(low = pal_rdbu[1], mid = pal_rdbu[4],
                        high = pal_rdbu[6], name = "scaled expr") +
  labs(title = "Linear attribution (ENet approximation)",
       x = "attribution (beta * scaled X)") +
  theme_pub(base_size = 11)
ggsave_pub(file.path(OUT4, "linear_beeswarm.pdf"), p2,
       width = 7, height = 6, dpi = 300)

cat("== 6/10 单细胞 UMAP ==\n")
obj <- readRDS(file.path(OUT6, "seurat_object_reviewed.rds"))
if ("celltype_reviewed" %in% colnames(obj@meta.data)) {
  p <- DimPlot(obj, group.by = "celltype_reviewed", label = TRUE,
               repel = TRUE, cols = pal_cat) +
    theme_pub(base_size = 12) +
    labs(title = "UMAP (reviewed annotation)")
  ggsave_pub(file.path(OUT6, "umap_celltype_reviewed.pdf"), p,
         width = 8.5, height = 7, dpi = 300)
}
if ("celltype" %in% colnames(obj@meta.data)) {
  p0 <- DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE,
                cols = pal_cat) +
    theme_pub(base_size = 12) +
    labs(title = "UMAP (module-score annotation)")
  ggsave_pub(file.path(OUT6, "umap_celltype.pdf"), p0,
         width = 8.5, height = 7, dpi = 300)
}

cat("== 7/10 Ro/e 热图与按样本比例热图 ==\n")
roe_df <- read.csv(file.path(OUT6, "roe_celltype.csv"), check.names = FALSE)
roe_long <- do.call(rbind, lapply(names(roe_df)[-1], function(g) {
  data.frame(celltype = roe_df[[1]], group = g, roe = roe_df[[g]],
             stringsAsFactors = FALSE)
}))
roe_long$log2_roe <- log2(roe_long$roe + 0.01)
roe_long$celltype <- factor(roe_long$celltype, levels = rev(roe_df[[1]]))
p_roe <- ggplot(roe_long, aes(x = group, y = celltype, fill = log2_roe)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = sprintf("%.2f", log2_roe)), size = 3.0,
            colour = "grey20") +
  scale_fill_gradient2(low = pal_rdbu[1], mid = pal_rdbu[4],
                       high = pal_rdbu[6], midpoint = 0,
                       name = "log2(Ro/e)") +
  labs(x = NULL, y = NULL,
       title = "Ro/e (log2) cell composition (reviewed)") +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        legend.key.height = unit(0.35, "cm"))
ggsave_pub(file.path(OUT6, "roe_heatmap.pdf"), p_roe,
       width = 5.6, height = 4.4, dpi = 300)

counts <- read.csv(file.path(OUT6, "per_sample_celltype_counts.csv"),
                   check.names = FALSE)
cts <- setdiff(colnames(counts), c("sample", "group", "total_nuclei"))
prop <- as.data.frame(sapply(cts, function(ct) counts[[ct]] / counts$total_nuclei))
rownames(prop) <- counts$sample
prop_long2 <- reshape2::melt(as.matrix(prop))
colnames(prop_long2) <- c("sample", "celltype", "proportion")
prop_long2$group <- counts$group[match(prop_long2$sample, counts$sample)]
p_prop <- ggplot(prop_long2, aes(x = sample, y = celltype, fill = proportion)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  scale_fill_gradientn(colours = pal_blue, name = "Proportion") +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL, title = "Cell-type proportion per sample") +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.5),
        axis.text.y = element_text(size = 7),
        strip.text = element_text(size = 9),
        legend.key.height = unit(0.4, "cm"))
ggsave_pub(file.path(OUT6, "celltype_proportion_heatmap.pdf"), p_prop,
       width = 8.5, height = 4.6, dpi = 300)

cat("== 8/10 按样本比例箱线图 ==\n")
prop_long <- reshape2::melt(as.matrix(prop))
colnames(prop_long) <- c("sample", "celltype", "proportion")
prop_long$group <- counts$group[match(prop_long$sample, counts$sample)]
p_box <- ggplot(prop_long, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.8) +
  facet_wrap(~ celltype, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(Sham = pal_cat[4], TBI = pal_cat[6])) +
  labs(x = NULL, y = "Proportion of nuclei") +
  theme_pub(base_size = 10) +
  theme(legend.position = "none",
        strip.text = element_text(size = 9))
ggsave_pub(file.path(OUT6, "celltype_proportions.pdf"), p_box,
       width = 9, height = 7, dpi = 300)

cat("== 9/10 Hub 基因 FeaturePlot ==\n")
cand_v4 <- toupper(read.csv(file.path(OUT_CAND, "candidate_genes.csv"))$gene)
hub <- intersect(unique(c(cand_v4, "P2RY6", "HSD11B1", "CCR5")),
                 rownames(obj))
if (length(hub) > 0) {
  p_fp <- FeaturePlot(obj, features = hub, ncol = 3,
                      cols = c("#F7F7F7", "#2171B5")) &
    theme_pub(base_size = 10) &
    theme(legend.position = "right")
  ggsave_pub(file.path(OUT6, "hub_featureplots.pdf"), p_fp,
         width = 10, height = 4.2 * ceiling(length(hub) / 3), dpi = 300)
}
if ("HubScore1" %in% colnames(obj@meta.data)) {
  p_hs <- FeaturePlot(obj, features = "HubScore1",
                      cols = c("#F7F7F7", "#2171B5")) +
    theme_pub(base_size = 11)
  ggsave_pub(file.path(OUT6, "hub_module_score.pdf"), p_hs,
         width = 6.5, height = 5.5, dpi = 300)
}

cat("== 10/10 固定签名 ROC（GSE128543 WT-only）==\n")
read_bulk <- function(gse) {
  readRDS(file.path(ROOT, "data", "01_GEO_bulk", paste0(gse, ".rds")))
}
infer_group <- function(titles) {
  t <- tolower(titles)
  is_ctrl <- grepl("\\b(sham|control|uninjured|vehicle|naive)\\b|notbi|no_tbi|no-tbi", t)
  is_tbi <- grepl("(?<![a-z])tbi(?![a-z])|mtbi|\\bcci\\b|\\binjured\\b|\\binjury\\b|\\binjur\\b|\\bimpact\\b|\\bfp\\b",
                  t, perl = TRUE) | grepl("fluid percussion", t)
  ifelse(is_tbi & !is_ctrl, "TBI",
         ifelse(is_ctrl & !is_tbi, "Sham", "Other"))
}
build_raw <- function(gse, genes, clean = FALSE, wt_only = FALSE) {
  d <- read_bulk(gse)
  ex <- d$expr
  rn <- toupper(rownames(ex))
  ex2 <- ex[rn %in% genes, , drop = FALSE]
  rn2 <- rn[rn %in% genes]
  if (anyDuplicated(rn2)) ex2 <- rowsum(ex2, group = rn2) else rownames(ex2) <- rn2
  ex2 <- ex2[order(rownames(ex2)), , drop = FALSE]
  grp <- infer_group(as.character(d$pheno$title))
  if (clean) {
    excl <- grepl("cyclophosphamide|ex-4-mtbi", d$pheno$title,
                  ignore.case = TRUE)
    grp[excl] <- "Other"
  }
  if (wt_only) {
    gt <- if ("genotype:ch1" %in% colnames(d$pheno)) {
      as.character(d$pheno$`genotype:ch1`)
    } else rep("", nrow(d$pheno))
    is_wt <- grepl("^\\s*wt\\s*$|wildtype", gt, ignore.case = TRUE) |
      grepl("\\bwt\\b|wildtype", d$pheno$title, ignore.case = TRUE)
    grp[!is_wt] <- "Other"
  }
  df <- as.data.frame(t(ex2), check.names = FALSE)
  df$group <- factor(ifelse(grp == "TBI", "TBI", "Sham"),
                     levels = c("Sham", "TBI"))
  df[grp %in% c("TBI", "Sham"), ]
}
sig3 <- c("P2RY6", "HSD11B1", "CCR5")
core_full_file <- file.path(ROOT, "results", "05_Candidate_Intersection",
                            "candidate_genes_FULL.csv")
core6 <- if (file.exists(core_full_file)) {
  toupper(read.csv(core_full_file)$gene)
} else {
  c("CYP1B1", "HSD11B1", "MMP2", "P2RY6", "RIPK1")
}
test_genes <- unique(c(core6, sig3))
te <- build_raw("GSE128543", test_genes, wt_only = TRUE)
te3 <- te[, c(intersect(sig3, colnames(te)), "group")]
set.seed(2026)
# 冻结系数（log2 尺度，来自 04b）
fcoef <- read.csv(file.path(ROOT, "results", "08_External_Validation",
                            "frozen_coefficients.csv"))
cf3 <- fcoef[fcoef$model == "CLEAN_3gene_only", ]
genes3 <- cf3$term[cf3$term != "(Intercept)"]
lp <- cf3$coefficient[cf3$term == "(Intercept)"]
for (g in genes3) lp <- lp + cf3$coefficient[cf3$term == g] * te3[[g]]
prob <- plogis(lp)
roc_obj <- roc(te3$group, prob, levels = c("Sham", "TBI"), direction = "<",
               quiet = TRUE)
stopifnot(as.numeric(auc(roc_obj)) >= 0.5)  # 方向校验：TBI 应获得更高预测值
roc_df <- data.frame(fpr = 1 - roc_obj$specificities,
                     tpr = roc_obj$sensitivities)
roc_df <- roc_df[order(roc_df$fpr, roc_df$tpr), ]  # FPR 升序，保证单调
roc_df <- roc_df[!duplicated(roc_df), ]
p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_abline(intercept = 0, slope = 1, color = "grey70", linetype = 2) +
  geom_step(color = "#4A6FA5", linewidth = 1.0, direction = "hv") +
  annotate("text", x = 0.68, y = 0.25,
           label = paste0("AUC = ", round(as.numeric(auc(roc_obj)), 3)),
           size = 4.5) +
  labs(title = "Frozen 3-feature signature on GSE128543 WT-only (n=14)",
       x = "1 - specificity", y = "Sensitivity") +
  theme_pub()
ggsave_pub(file.path(OUT8, "frozen_roc.pdf"), p_roc,
       width = 5.5, height = 5, dpi = 300)

cat("全部主图已重绘（dpi=300，统一主题/配色）。\n")
