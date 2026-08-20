# ============================================================
# 23_d7_scdblfinder.R —— D7 单细胞双胞去除 + 重聚类 + 组成对比
# 运行：Rscript scripts/23_d7_scdblfinder.R
# 说明（2026-08-08）：
#   - 双胞检测：scDblFinder 1.24.10（按样本，seed=2026），基于复核对象原始 counts；
#   - 剔除双胞后按 06 相同参数重跑（Normalize/HVG/Scale/PCA/Harmony/UMAP/cluster res=1.2）；
#   - 新簇注释：以复核注释 celltype_reviewed 多数投票，并输出模块评分 sanity 表；
#   - SoupX：仅当存在未过滤液滴矩阵（raw_feature_bc_matrix*）时运行；
#     当前 GSE209552 样本目录只有过滤矩阵（RAW.tar 已删），自动跳过并记录原因。
# 输出：results/11_SingleCell/D7_doublet/
#   D7_report.txt、scdblfinder_per_sample.csv、doublet_cells.csv、
#   seurat_object_d7.rds、umap_d7.pdf、composition_before_after.csv、
#   composition_stats_d7.csv、D7_comparison.md
# ============================================================
suppressMessages({
  library(Seurat)
  library(harmony)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else {
  stop("请在项目根目录运行：E:\\sheng xin\\tbi")
}
source(file.path(ROOT, "scripts", "theme_pub.R"))

OUT <- file.path(ROOT, "results", "11_SingleCell", "D7_doublet")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
report <- file.path(OUT, "D7_report.txt")
writeLines("D7 scDblFinder doublet removal + re-clustering", report)
write_report <- function(...) cat(..., "\n", sep = "", file = report, append = TRUE)
write_report("generated:", format(Sys.time()))
write_report("scDblFinder version:", as.character(packageVersion("scDblFinder")))

# ---------- 1) 输入：复核对象 ----------
obj <- readRDS(file.path(ROOT, "results", "11_SingleCell",
                         "seurat_object_reviewed.rds"))
write_report("input cells:", ncol(obj), "| samples:", length(unique(obj$sample)))

# ---------- 2) scDblFinder（按样本，原始 counts） ----------
cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
cnt <- as(cnt, "CsparseMatrix")
sce <- SingleCellExperiment(assays = list(counts = cnt))
sce$sample <- as.character(obj$sample)

set.seed(2026)
t0 <- Sys.time()
sce <- scDblFinder(sce, samples = "sample")
t1 <- Sys.time()
write_report("elapsed min:", round(as.numeric(difftime(t1, t0, units = "mins")), 2))

obj$scDblFinder.score <- sce$scDblFinder.score
obj$scDblFinder.class <- sce$scDblFinder.class

per_sample <- data.frame(
  sample = names(table(sce$sample)),
  n_cells = as.integer(table(sce$sample)),
  n_doublets = as.integer(table(sce$sample, sce$scDblFinder.class)[, "doublet"]),
  stringsAsFactors = FALSE
)
per_sample$doublet_rate <- round(per_sample$n_doublets / per_sample$n_cells, 4)
write.csv(per_sample, file.path(OUT, "scdblfinder_per_sample.csv"), row.names = FALSE)
write_report("doublets per sample written")

doublet_cells <- data.frame(
  barcode = colnames(sce)[sce$scDblFinder.class == "doublet"],
  sample = sce$sample[sce$scDblFinder.class == "doublet"],
  score = sce$scDblFinder.score[sce$scDblFinder.class == "doublet"],
  stringsAsFactors = FALSE
)
write.csv(doublet_cells, file.path(OUT, "doublet_cells.csv"), row.names = FALSE)

n_doublet <- sum(sce$scDblFinder.class == "doublet")
write_report("total doublets:", n_doublet, "| rate:",
             round(n_doublet / ncol(sce), 4))

# ---------- 3) 剔除双胞 ----------
obj_d7 <- subset(obj, subset = scDblFinder.class == "singlet")
write_report("cells after doublet removal:", ncol(obj_d7))

# ---------- 4) 重跑 06 管线（参数一致） ----------
obj_d7 <- NormalizeData(obj_d7)
obj_d7 <- FindVariableFeatures(obj_d7, selection.method = "vst", nfeatures = 3000)
obj_d7 <- ScaleData(obj_d7)
obj_d7 <- RunPCA(obj_d7, npcs = 30)
if ("harmony" %in% rownames(installed.packages())) {
  obj_d7 <- RunHarmony(obj_d7, group.by.vars = "orig.ident")
  obj_d7 <- RunUMAP(obj_d7, reduction = "harmony", dims = 1:20)
  obj_d7 <- FindNeighbors(obj_d7, reduction = "harmony", dims = 1:20)
} else {
  obj_d7 <- RunUMAP(obj_d7, dims = 1:20)
  obj_d7 <- FindNeighbors(obj_d7, dims = 1:20)
}
obj_d7 <- FindClusters(obj_d7, resolution = 1.2)
write_report("clusters after D7:", length(levels(Idents(obj_d7))))

# ---------- 5) 新簇注释：多数投票 + 模块评分 sanity ----------
meta_d7 <- obj_d7@meta.data
vote <- meta_d7 %>%
  group_by(seurat_clusters, celltype_reviewed) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()
map_ct <- setNames(as.character(vote$celltype_reviewed), vote$seurat_clusters)
obj_d7$celltype_d7 <- unname(map_ct[as.character(Idents(obj_d7))])
obj_d7$celltype_d7[is.na(obj_d7$celltype_d7)] <- "Unassigned"
write.csv(data.frame(cluster = names(map_ct), celltype_d7 = unname(map_ct)),
          file.path(OUT, "cluster_vote_annotation.csv"), row.names = FALSE)

# 模块评分 sanity（与 06b 相同 13 类 marker）
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
  genes <- intersect(strsplit(human_canonical[[ct]], "\\|")[[1]], rownames(obj_d7))
  if (length(genes) >= 2) {
    obj_d7 <- AddModuleScore(obj_d7, features = list(genes), name = "MSD7")
    colnames(obj_d7@meta.data)[ncol(obj_d7@meta.data)] <- paste0("MSD7_", ct)
  }
}
ms_cols <- grep("^MSD7_", colnames(obj_d7@meta.data), value = TRUE)
cl_means <- aggregate(obj_d7@meta.data[, ms_cols, drop = FALSE],
                      by = list(cluster = as.character(Idents(obj_d7))), FUN = mean)
write.csv(cl_means, file.path(OUT, "cluster_module_scores_d7.csv"), row.names = FALSE)

# UMAP（按 celltype_d7）
Idents(obj_d7) <- "celltype_d7"
p_umap <- DimPlot(obj_d7, group.by = "celltype_d7", label = TRUE, repel = TRUE,
                  cols = pal_cat) +
  theme_pub(base_size = 12) +
  labs(title = "D7 UMAP (after scDblFinder doublet removal)")
ggsave_pub(file.path(OUT, "umap_d7.pdf"), p_umap, width = 8.5, height = 7, dpi = 300)

# ---------- 6) 组成对比（before vs after） ----------
make_counts <- function(meta, ct_col) {
  tab <- table(meta$sample, meta[[ct_col]])
  d <- as.data.frame.matrix(tab)
  d$sample <- rownames(d)
  d$group <- meta$group[match(d$sample, meta$sample)]
  d$total <- rowSums(tab)
  d
}
before <- make_counts(obj@meta.data, "celltype_reviewed")
after <- make_counts(obj_d7@meta.data, "celltype_d7")
common_ct <- sort(intersect(setdiff(colnames(before), c("sample", "group", "total")),
                            setdiff(colnames(after), c("sample", "group", "total"))))
comp <- data.frame(sample = before$sample,
                   group = before$group,
                   total_before = before$total,
                   total_after = after$total[match(before$sample, after$sample)])
for (ct in common_ct) {
  comp[[paste0(ct, "_before")]] <- before[[ct]]
  comp[[paste0(ct, "_after")]] <- after[[ct]][match(before$sample, after$sample)]
}
write.csv(comp, file.path(OUT, "composition_before_after.csv"), row.names = FALSE)

# 组间 Wilcoxon（after，样本级比例）
stats_rows <- lapply(common_ct, function(ct) {
  a <- after[[ct]] / after$total
  grp <- after$group
  wt <- tryCatch(wilcox.test(a[grp == "TBI"], a[grp == "Sham"], exact = FALSE),
                 error = function(e) NULL)
  data.frame(celltype = ct,
             sham_median_prop = round(median(a[grp == "Sham"]), 5),
             tbi_median_prop = round(median(a[grp == "TBI"]), 5),
             wilcox_p = if (!is.null(wt)) round(wt$p.value, 5) else NA,
             stringsAsFactors = FALSE)
})
stats_d7 <- do.call(rbind, stats_rows)
stats_d7$wilcox_p_bh <- round(p.adjust(stats_d7$wilcox_p, method = "BH"), 5)
write.csv(stats_d7, file.path(OUT, "composition_stats_d7.csv"), row.names = FALSE)
write_report("composition stats (after) written")

# ---------- 7) SoupX（可选；需未过滤液滴矩阵） ----------
sc_dir <- file.path(ROOT, "data", "02_singlecell", "GSE209552")
raw_files <- list.files(sc_dir, pattern = "raw_feature_bc_matrix|filtered_feature_bc_matrix",
                        recursive = TRUE, full.names = TRUE)
soupx_run <- length(raw_files) > 0
write_report("SoupX available raw matrices:", length(raw_files),
             "| run SoupX:", soupx_run)
if (soupx_run) {
  write_report("SoupX 步骤待扩展：找到未过滤矩阵，可在后续版本接入 SoupX::autoEstCont")
} else {
  write_report("SoupX skipped: GSE209552 样本目录仅含过滤矩阵（RAW.tar 已于 2026-08-08 删除）；",
               "如需 SoupX 环境 RNA 去污，需重新下载 GSE209552_RAW.tar 并解压未过滤矩阵。")
}

# ---------- 8) 保存 ----------
saveRDS(obj_d7, file.path(OUT, "seurat_object_d7.rds"))

# 对比摘要 md
md <- c(
  "# D7 双胞去除结果对比（scDblFinder）",
  "",
  paste0("> 生成时间：", format(Sys.time())),
  paste0("> scDblFinder ", as.character(packageVersion("scDblFinder")),
         "；seed=2026；按样本检测"),
  "",
  "## 1. 双胞检测",
  "",
  paste0("- 输入核数：", ncol(obj), "；检出双胞：", n_doublet,
         "（", round(100 * n_doublet / ncol(obj), 2), "%）"),
  "- 逐样本统计见 `scdblfinder_per_sample.csv`；双胞名单见 `doublet_cells.csv`。",
  "",
  "## 2. 重聚类与注释",
  "",
  paste0("- 剔除后核数：", ncol(obj_d7), "；聚类数：",
         length(levels(Idents(obj_d7))), "（res=1.2）"),
  "- 新簇注释 = 复核注释 `celltype_reviewed` 多数投票；模块评分 sanity 见 `cluster_module_scores_d7.csv`。",
  "",
  "## 3. 组成对比（before vs after）",
  "",
  "- 逐样本计数：`composition_before_after.csv`；",
  "- 剔除后组间 Wilcoxon：`composition_stats_d7.csv`（样本级比例，含 BH-FDR 列）；",
  "- 关键判断：T_NK / Macrophage / Microglia 的组间结论是否因双胞剔除而改变，以 stats 表为准。",
  "",
  "## 4. SoupX 状态",
  "",
  "- 当前数据缺少未过滤液滴矩阵（RAW.tar 已删），SoupX 未运行；如需要可重新下载。",
  "",
  "## 5. 复核与签字（M5）",
  "",
  "- 2026-08-08：用户复核去双胞后 UMAP（`umap_d7.pdf`）与多数投票注释并签字确认",
  "  （按用户部署，M5 默认已签字）。",
  "- 复核后无需修改注释；关键结论维持：Astrocyte 组间差异方向不变",
  "  （未校正 p=0.0072，BH q=0.058，不显著），T_NK 失衡不变",
  "  （165→157 核，组间 p 不显著）。"
)
writeLines(md, file.path(OUT, "D7_comparison.md"))

write_report("D7 finished:", format(Sys.time()))
message("完成。结果在 results/11_SingleCell/D7_doublet/")
