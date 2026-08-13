# ============================================================
# 18_d7_doublet_qc.R —— D7 双胞检测与 QC 重跑（GSE209552）
# 依据：scrna-seq-qc 技能（scDblFinder 按样本检测为必需方法）
# 流程：
#   1) 读取 seurat_object_reviewed.rds（raw counts + celltype_reviewed）
#   2) QC 指标（nFeature/nCount/percent.mt）复核
#   3) scDblFinder 按样本检测双胞（17 个样本分别运行）
#   4) 构建 passes_QC 标记（保留非双胞细胞），重跑
#      Normalize/HVG/Scale/PCA/Harmony/UMAP/FindClusters
#   5) 13 类 canonical marker 模块评分重注释新簇
#   6) 前后组成对比（按样本比例 + Wilcoxon）与 T_NK/Macrophage 结论对比
#   7) 新 UMAP 输出（交 M5 复核）
# 环境 RNA（SoupX）：仅当存在空滴/未过滤矩阵时才能可靠估计；
#   本项目 GEO 矩阵为过滤后计数（无空滴），故如实记录限制，不强行去污。
# 输出：results/11_SingleCell/D7_doublet_qc/
# ============================================================
suppressMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(SingleCellExperiment)
})
if (!requireNamespace("scDblFinder", quietly = TRUE)) {
  stop("scDblFinder 未安装，请先运行 BiocManager::install('scDblFinder')")
}

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "results", "11_SingleCell")
OUT_D7 <- file.path(OUT, "D7_doublet_qc")
dir.create(OUT_D7, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts", "theme_pub.R"))

LOG <- file.path(OUT_D7, "D7_run.log")
if (file.exists(LOG)) file.remove(LOG)
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(...), "\n")
  cat(msg)
  cat(msg, file = LOG, append = TRUE)
}

log_msg("== D7 开始 ==")
obj <- readRDS(file.path(OUT, "seurat_object_reviewed.rds"))
meta <- obj@meta.data
stopifnot(all(c("sample", "group", "celltype_reviewed") %in% colnames(meta)))
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
log_msg("输入细胞数: ", ncol(obj), " | 基因数: ", nrow(obj))

# ---------- 1) QC 指标 ----------
if (!"percent.mt" %in% colnames(meta)) {
  mt_genes <- grep("^MT-", rownames(obj), value = TRUE)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mt_genes)
  meta <- obj@meta.data
}
qc_df <- data.frame(
  sample = meta$sample, group = meta$group,
  nFeature_RNA = meta$nFeature_RNA,
  nCount_RNA = meta$nCount_RNA,
  percent_mt = meta$percent.mt,
  stringsAsFactors = FALSE)
write.csv(qc_df, file.path(OUT_D7, "qc_metrics_input.csv"), row.names = FALSE)
log_msg("QC 指标: nFeature 中位=", median(meta$nFeature_RNA),
        " nCount 中位=", median(meta$nCount_RNA),
        " percent.mt 中位=", round(median(meta$percent.mt), 2))

# ---------- 2) scDblFinder 按样本 ----------
log_msg("运行 scDblFinder（按样本，17 个样本）...")
samples <- sort(unique(meta$sample))
db_score <- rep(NA_real_, ncol(obj))
db_class <- rep(NA_character_, ncol(obj))
db_summary <- data.frame(sample = samples, n_cells = NA_integer_,
                         n_doublets = NA_integer_,
                         pct_doublets = NA_real_, stringsAsFactors = FALSE)
for (s in samples) {
  idx <- which(meta$sample == s)
  sce <- SingleCellExperiment(assays = list(counts = counts[, idx]))
  sce <- scDblFinder::scDblFinder(sce, BPPARAM = BiocParallel::SerialParam())
  db_score[idx] <- sce$scDblFinder.score
  db_class[idx] <- as.character(sce$scDblFinder.class)
  db_summary[db_summary$sample == s, c("n_cells", "n_doublets")] <-
    c(length(idx), sum(sce$scDblFinder.class == "doublet"))
  log_msg(s, ": n=", length(idx),
          " doublets=", sum(sce$scDblFinder.class == "doublet"))
}
db_summary$pct_doublets <- 100 * db_summary$n_doublets / db_summary$n_cells
meta$doublet_score <- db_score
meta$doublet_class <- db_class
meta$passes_QC <- !(meta$doublet_class == "doublet")
obj$doublet_score <- db_score
obj$doublet_class <- db_class
obj$passes_QC <- meta$passes_QC
write.csv(db_summary, file.path(OUT_D7, "doublet_summary_by_sample.csv"),
          row.names = FALSE)
log_msg("双胞总数: ", sum(!meta$passes_QC), "/", ncol(obj),
        " (", round(100 * mean(!meta$passes_QC), 2), "%)")

# 双胞在原始细胞类型中的分布（用于解释 T_NK/Macrophage 是否受双胞影响）
tab_dbl <- table(meta$celltype_reviewed, meta$doublet_class)
write.csv(as.data.frame.matrix(tab_dbl), file.path(OUT_D7, "doublet_by_celltype.csv"))
log_msg("双胞按细胞类型分布:")
print(tab_dbl)

# ---------- 3) 去双胞后重跑聚类 ----------
obj_f <- subset(obj, cells = colnames(obj)[meta$passes_QC])
log_msg("去双胞后细胞数: ", ncol(obj_f))
saveRDS(obj_f, file.path(OUT_D7, "seurat_object_doublet_filtered.rds"))

obj_f <- NormalizeData(obj_f)
obj_f <- FindVariableFeatures(obj_f, selection.method = "vst", nfeatures = 3000)
obj_f <- ScaleData(obj_f)
obj_f <- RunPCA(obj_f, npcs = 30, verbose = FALSE)
obj_f <- RunHarmony(obj_f, group.by.vars = "sample")
obj_f <- RunUMAP(obj_f, reduction = "harmony", dims = 1:20,
                 reduction.name = "umap_d7", reduction.key = "umapd7_")
obj_f <- FindNeighbors(obj_f, reduction = "harmony", dims = 1:20)
obj_f <- FindClusters(obj_f, resolution = 1.2)
log_msg("重聚类簇数: ", length(unique(as.character(Idents(obj_f)))))

# ---------- 4) 13 类 canonical marker 模块评分重注释 ----------
old_ms <- grep("^MS_", colnames(obj_f@meta.data), value = TRUE)
if (length(old_ms) > 0) {
  obj_f@meta.data <- obj_f@meta.data[, setdiff(colnames(obj_f@meta.data),
                                               old_ms), drop = FALSE]
}
human_canonical <- c(
  Microglia = "P2RY12|CX3CR1|TMEM119|CSF1R|TYROBP|AIF1|ITGAM|C1QA|C1QB|C1QC",
  Macrophage = "LYZ|CD68|CD163|MRC1|MSR1",
  Astrocyte = "GFAP|AQP4|SLC1A3|GJA1|S100B|ALDH1L1",
  Neuron = "SNAP25|RBFOX3|SYT1|MEG3|NRGN|STMN2",
  Oligodendrocyte = "MBP|MOG|PLP1|MOBP|CLDN11|MAG",
  OPC = "PDGFRA|CSPG4|OLIG1|OLIG2",
  Endothelial = "CLDN5|FLT1|PECAM1|VWF|EGFL7",
  Pericyte = "RGS5|PDGFRB|KCNJ8|NOTCH3",
  Fibroblast = "COL1A1|DCN|COL1A2|LUM",
  SMC = "ACTA2|MYH11|TAGLN",
  T_NK = "CD3D|CD3E|NKG7|CD8A|CD2",
  B = "CD79A|MS4A1|CD19",
  Monocyte = "CD14|FCGR3A|S100A8|S100A9"
)
for (ct in names(human_canonical)) {
  genes <- intersect(strsplit(human_canonical[[ct]], "\\|")[[1]], rownames(obj_f))
  if (length(genes) >= 2) {
    obj_f <- AddModuleScore(obj_f, features = list(genes), name = "D7MS",
                            assay = "RNA")
    colnames(obj_f@meta.data)[ncol(obj_f@meta.data)] <- paste0("D7MS_", ct)
  }
}
ms_cols <- grep("^D7MS_", colnames(obj_f@meta.data), value = TRUE)
cl <- as.character(Idents(obj_f))
cl_means <- aggregate(obj_f@meta.data[, ms_cols, drop = FALSE],
                      by = list(cluster = cl), FUN = mean)
assign_ct <- function(k) {
  row <- cl_means[cl_means$cluster == k, ms_cols, drop = FALSE]
  if (nrow(row) == 0) return(NA_character_)
  sub("^D7MS_", "", names(which.max(row))[1])
}
new_ct <- unname(vapply(sort(unique(cl)), assign_ct, character(1)))
map_new <- setNames(new_ct, sort(unique(cl)))
meta_f <- obj_f@meta.data
meta_f$celltype_d7 <- unname(map_new[cl])
obj_f$celltype_d7 <- meta_f$celltype_d7

# 簇证据表（模块评分均值 + 每簇细胞数/样本数）
ev <- do.call(rbind, lapply(sort(unique(cl)), function(k) {
  idx <- which(cl == k)
  data.frame(cluster = k, n_cells = length(idx),
             n_samples = length(unique(meta_f$sample[idx])),
             celltype_d7 = map_new[k],
             t(round(colMeans(meta_f[idx, ms_cols, drop = FALSE]), 3)),
             check.names = FALSE, stringsAsFactors = FALSE)
}))
write.csv(ev, file.path(OUT_D7, "cluster_evidence_d7.csv"), row.names = FALSE)
saveRDS(obj_f, file.path(OUT_D7, "seurat_object_d7_annotated.rds"))
log_msg("新注释类型: ", paste(sort(unique(meta_f$celltype_d7)), collapse = ", "))

# ---------- 5) 前后组成对比 ----------
prop_before <- prop.table(table(meta$sample, meta$celltype_reviewed), margin = 1)
prop_after <- prop.table(table(meta_f$sample, meta_f$celltype_d7), margin = 1)
shared_ct <- sort(union(colnames(prop_before), colnames(prop_after)))
comp_rows <- list()
for (ct in shared_ct) {
  b <- if (ct %in% colnames(prop_before)) prop_before[, ct] else
    setNames(rep(0, nrow(prop_before)), rownames(prop_before))
  a <- if (ct %in% colnames(prop_after)) prop_after[, ct] else
    setNames(rep(0, nrow(prop_after)), rownames(prop_after))
  smp <- intersect(names(b), names(a))
  g <- meta$group[match(smp, meta$sample)]
  wt_b <- suppressWarnings(wilcox.test(b[smp][g == "Sham"],
                                       b[smp][g == "TBI"]))
  wt_a <- suppressWarnings(wilcox.test(a[smp][g == "Sham"],
                                       a[smp][g == "TBI"]))
  comp_rows[[length(comp_rows) + 1]] <- data.frame(
    celltype = ct,
    mean_before_sham = mean(b[smp][g == "Sham"]),
    mean_before_tbi = mean(b[smp][g == "TBI"]),
    p_before = wt_b$p.value,
    mean_after_sham = mean(a[smp][g == "Sham"]),
    mean_after_tbi = mean(a[smp][g == "TBI"]),
    p_after = wt_a$p.value,
    stringsAsFactors = FALSE)
}
comp <- do.call(rbind, comp_rows)
comp$fdr_after <- p.adjust(comp$p_after, method = "BH")
write.csv(comp, file.path(OUT_D7, "composition_before_after.csv"),
          row.names = FALSE)
log_msg("组成对比表已写入 composition_before_after.csv")
print(comp)

# ---------- 6) UMAP 图 ----------
pal <- c("#0072B2", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
         "#D55E00", "#CC79A7", "#999999", "#A65628", "#F781BF",
         "#4DAF4A", "#FF7F00", "#377EB8")
names(pal) <- sort(unique(meta_f$celltype_d7))

pdf_pub(file.path(OUT_D7, "umap_d7_annotated.pdf"), width = 8, height = 6)
p <- DimPlot(obj_f, reduction = "umap_d7", group.by = "celltype_d7",
             cols = pal[intersect(names(pal), unique(meta_f$celltype_d7))],
             label = TRUE, repel = TRUE) +
  labs(title = paste0("After doublet QC (", ncol(obj_f), " nuclei)")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "right", panel.grid.minor = element_blank())
print(p)
dev.off()
ggsave_pub(file.path(OUT_D7, "umap_d7_annotated.pdf"), p,
       width = 8, height = 6, dpi = 300)
log_msg("UMAP 已保存（umap_d7_annotated.pdf/pdf）")

# ---------- 7) T_NK/Macrophage 结论对比 ----------
col_safe <- function(tab, ct) {
  if (ct %in% colnames(tab)) {
    tab[, ct, drop = FALSE]
  } else {
    setNames(rep(0, nrow(tab)), rownames(tab))
  }
}
tab_b <- table(meta$sample, meta$celltype_reviewed)
tab_a <- table(meta_f$sample, meta_f$celltype_d7)
tnk_before <- col_safe(tab_b, "T_NK")
tnk_after <- col_safe(tab_a, "T_NK")
macro_before <- col_safe(tab_b, "Macrophage")
macro_after <- col_safe(tab_a, "Macrophage")
g_before <- meta$group[match(rownames(tnk_before), meta$sample)]
g_after <- meta_f$group[match(rownames(tnk_after), meta_f$sample)]
summ <- data.frame(
  item = c("total_nuclei", "T_NK_n", "T_NK_sham_n", "T_NK_tbi_n",
           "Macrophage_n", "Macrophage_sham_n", "Macrophage_tbi_n"),
  before = c(ncol(obj),
             sum(tnk_before), sum(tnk_before[g_before == "Sham"]),
             sum(tnk_before[g_before == "TBI"]),
             sum(macro_before),
             sum(macro_before[g_before == "Sham"]),
             sum(macro_before[g_before == "TBI"])),
  after = c(ncol(obj_f),
            sum(tnk_after),
            sum(tnk_after[g_after == "Sham"]),
            sum(tnk_after[g_after == "TBI"]),
            sum(macro_after),
            sum(macro_after[g_after == "Sham"]),
            sum(macro_after[g_after == "TBI"]))
)
write.csv(summ, file.path(OUT_D7, "tnk_macrophage_before_after.csv"),
          row.names = FALSE)
log_msg("T_NK/Macrophage 前后对比:")
print(summ)

# 环境 RNA 限制说明
writeLines(c(
  "# D7 环境 RNA（SoupX）说明",
  "",
  "GEO 提供的 GSE209552 逐样本矩阵为过滤后计数矩阵（未见空滴/未过滤",
  "raw_feature_bc 矩阵）。SoupX 需要空滴数据估计 ambient profile；",
  "在无空滴数据时强行去污会引入未经检验的假设。因此 D7 只执行",
  "scDblFinder 双胞检测与 QC 重跑，环境 RNA 影响如实记为未评估限制。",
  "如需评估：需向 GSE209552 作者索取 raw_feature_bc_matrix 或重新下载",
  "原始 Cell Ranger 输出。"
), file.path(OUT_D7, "AMBIENT_RNA_LIMITATION.md"))

log_msg("== D7 完成 ==")
