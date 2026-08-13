# ============================================================
# 13_d8_astro.R —— D8 星形胶质亚群分析（脑瘢痕类比模板成纤维模块）
# 依据：REVIEW4_FINAL_AUDIT D8
# 输入：results/11_SingleCell/seurat_object_reviewed.rds（Astrocyte 子集）
# 输出：results/14_Astrocyte_Subclusters/D8_astro/
# 说明：人脑 snRNA-seq 无成纤维细胞，以星形胶质替代（探索性，受对照不可比限制）。
# ============================================================

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "results", "14_Astrocyte_Subclusters", "D8_astro")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts", "theme_pub.R"))
report <- file.path(OUT, "D8_report.txt")
writeLines("D8 astrocyte subclustering", report)
write_report <- function(...) cat(..., "\n", sep = "", file = report,
                                  append = TRUE)
write_report("generated:", format(Sys.time()))

obj <- readRDS(file.path(ROOT, "results", "11_SingleCell",
                         "seurat_object_reviewed.rds"))
astro <- subset(obj, subset = celltype_reviewed == "Astrocyte")
write_report("Astrocyte cells:", ncol(astro),
             "| Sham:", sum(astro$group == "Sham"),
             "| TBI:", sum(astro$group == "TBI"))

astro <- NormalizeData(astro)
astro <- FindVariableFeatures(astro, selection.method = "vst",
                              nfeatures = 2000)
astro <- ScaleData(astro)
astro <- RunPCA(astro, npcs = 20)
astro <- RunUMAP(astro, dims = 1:15)
astro <- FindNeighbors(astro, dims = 1:15)
astro <- FindClusters(astro, resolution = 0.5)
write_report("subclusters:", length(levels(Idents(astro))))

markers <- FindAllMarkers(astro, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25)
write.csv(markers, file.path(OUT, "astro_markers.csv"), row.names = FALSE)
top10 <- do.call(rbind, lapply(sort(unique(markers$cluster)), function(k) {
  head(markers[markers$cluster == k, c("cluster", "gene", "avg_log2FC",
                                       "pct.1", "pct.2", "p_val_adj")], 10)
}))
write.csv(top10, file.path(OUT, "astro_top10_markers.csv"), row.names = FALSE)

# 模块评分：ECM 重塑 / TGFb 反应 / 应激-炎症 / 收缩（肌成纤维样）
mods <- list(
  ECM = c("FN1", "COL1A1", "COL3A1", "COL4A1", "COL5A1", "COL6A1", "COL6A2",
          "COL6A3", "LAMA2", "LAMA4", "LAMB2", "TNC", "VCAN", "DCN", "SPARC",
          "BGN", "LUM", "TIMP1", "TIMP2", "MMP2", "MMP9"),
  TGFb_response = c("TGFBR1", "TGFBR2", "SMAD3", "SMAD4", "SERPINE1",
                    "TAGLN", "CTGF", "CYR61", "JUNB", "FOSB", "ATF3"),
  stress_inflammation = c("GFAP", "VIM", "LIF", "OSM", "IL6", "C3", "SERPINA3",
                          "HSP90AA1", "HSPA1A", "HSPB1", "NFKBIA", "JUN", "FOS"),
  contractile = c("ACTA2", "MYH11", "TAGLN", "CNN1", "TPM1", "TPM2", "MYL9",
                  "DES")
)
for (m in names(mods)) {
  genes <- intersect(mods[[m]], rownames(astro))
  if (length(genes) >= 2) {
    astro <- AddModuleScore(astro, features = list(genes), name = paste0("M_", m))
    colnames(astro@meta.data)[ncol(astro@meta.data)] <- paste0("M_", m)
    write_report(m, "module genes used:", length(genes))
  }
}

ms_cols <- grep("^M_", colnames(astro@meta.data), value = TRUE)
cl_summary <- do.call(rbind, lapply(levels(Idents(astro)), function(k) {
  idx <- which(Idents(astro) == k)
  data.frame(cluster = k,
             n_cells = length(idx),
             n_samples = length(unique(astro$sample[idx])),
             t(round(colMeans(astro@meta.data[idx, ms_cols, drop = FALSE]), 3)),
             check.names = FALSE)
}))
write.csv(cl_summary, file.path(OUT, "astro_cluster_scores.csv"),
          row.names = FALSE)

p_umap <- DimPlot(astro, label = TRUE, repel = TRUE,
                  cols = c(pal_cat, "#999999")) +
  theme_pub(base_size = 11) +
  labs(title = "Astrocyte subclusters (res=0.5)")
ggsave_pub(file.path(OUT, "astro_umap.pdf"), p_umap, width = 7, height = 6,
       dpi = 300)

# 亚群 Ro/e 与按样本分布
tab <- table(Idents(astro), astro$group)
rowsum <- rowSums(tab); colsum <- colSums(tab); total <- sum(tab)
roe <- (tab / rowsum) / (matrix(colsum, nrow = nrow(tab), ncol = ncol(tab),
                                byrow = TRUE) / total)
roe[!is.finite(roe)] <- 0
write.csv(as.data.frame.matrix(roe), file.path(OUT, "astro_roe.csv"))
sample_tab <- as.data.frame.matrix(table(astro$sample, Idents(astro)))
sample_tab$sample <- rownames(sample_tab)
sample_tab$group <- astro$group[match(rownames(sample_tab), astro$sample)]
sample_tab$total_astro <- rowSums(table(astro$sample, Idents(astro)))
write.csv(sample_tab, file.path(OUT, "astro_per_sample_counts.csv"),
          row.names = FALSE)

saveRDS(astro, file.path(OUT, "astro_subset_v2.rds"))
write_report("D8 finished:", format(Sys.time()))
message("D8 完成。结果在 ", OUT)
