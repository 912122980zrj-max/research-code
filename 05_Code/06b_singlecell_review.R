# ============================================================
# 06b_singlecell_review.R —— 单细胞注释复核与按样本统计
# 运行：
#   Rscript scripts/06b_singlecell_review.R --stage evidence
#   （人工/AI 复核后生成 annotation_review.csv）
#   Rscript scripts/06b_singlecell_review.R --stage apply
# 说明（2026-08-07 审查修复版，T4/P5/P10/P11/P15）：
#   - 不重跑聚类；基于 results/11_SingleCell/seurat_object.rds 复核。
#   - 输出：cluster_evidence.csv、cluster_top10_markers.csv、
#     per_sample_celltype_counts.csv、celltype_composition_stats.csv、
#     roe_celltype.csv（复核后）、umap_celltype_reviewed.pdf、
#     annotation_review.md、seurat_object_reviewed.rds。
# ============================================================

suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
stage <- if (any(grepl("evidence", args))) "evidence" else "apply"

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "results", "11_SingleCell")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts", "theme_pub.R"))

LOG <- file.path(OUT, "06b_run.log")
if (stage == "evidence" && file.exists(LOG)) file.remove(LOG)
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(...), "\n")
  cat(msg)
  cat(msg, file = LOG, append = TRUE)
}

obj_file <- file.path(OUT, "seurat_object.rds")
if (!file.exists(obj_file)) stop("缺少 seurat_object.rds，请先运行 06")
obj <- readRDS(obj_file)
meta <- obj@meta.data
log_msg("cells: ", ncol(obj))

if (!"seurat_clusters" %in% colnames(meta)) {
  log_msg("meta 中无 seurat_clusters，尝试从 checkpoint 恢复")
  cp <- file.path(OUT, "seurat_checkpoint.rds")
  if (file.exists(cp)) {
    obj2 <- readRDS(cp)
    meta$seurat_clusters <- as.character(obj2$seurat_clusters[match(
      rownames(meta), colnames(obj2))])
  } else {
    stop("无法获取簇标签")
  }
}

ms_cols <- grep("^MS_", colnames(meta), value = TRUE)
ct_cols <- setdiff(gsub("^MS_", "", ms_cols), "HubScore")

# ---------- Stage 1: 证据输出 ----------
if (stage == "evidence") {
  markers_file <- file.path(OUT, "cluster_markers.csv")
  markers <- if (file.exists(markers_file)) read.csv(markers_file) else NULL

  cl <- as.character(meta$seurat_clusters)
  top10 <- NULL
  if (!is.null(markers) && nrow(markers) > 0) {
    fld <- if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else "avg_logFC"
    markers <- markers[order(markers$cluster, -markers[[fld]]), ]
    top10 <- do.call(rbind, lapply(sort(unique(markers$cluster)), function(k) {
      head(markers[markers$cluster == k, c("cluster", "gene", fld,
                                           "pct.1", "pct.2", "p_val_adj")],
           10)
    }))
    write.csv(top10, file.path(OUT, "cluster_top10_markers.csv"),
              row.names = FALSE)
  }

  cl_summary <- do.call(rbind, lapply(sort(unique(cl)), function(k) {
    idx <- which(cl == k)
    m <- colMeans(meta[idx, ms_cols, drop = FALSE])
    orig_ct <- as.character(meta$celltype[idx[1]])
    data.frame(
      cluster = k,
      n_cells = length(idx),
      n_samples = length(unique(meta$sample[idx])),
      original_celltype = orig_ct,
      t(round(m, 3)),
      check.names = FALSE,
      stringsAsFactors = FALSE)
  }))
  write.csv(cl_summary, file.path(OUT, "cluster_evidence.csv"),
            row.names = FALSE)

  # 未复核的按样本计数（备查）
  tab <- table(meta$sample, meta$celltype)
  tab_df <- as.data.frame.matrix(tab)
  tab_df$sample <- rownames(tab_df)
  tab_df$group <- meta$group[match(rownames(tab_df), meta$sample)]
  tab_df$total_nuclei <- rowSums(tab)
  write.csv(tab_df, file.path(OUT, "per_sample_celltype_counts_raw.csv"),
            row.names = FALSE)

  log_msg("evidence stage 完成：cluster_evidence.csv / cluster_top10_markers.csv")
  message("请人工/AI 复核 cluster_evidence.csv 与 cluster_top10_markers.csv，")
  message("生成 annotation_review.csv（列：cluster, celltype_reviewed, evidence, action）")
  quit(save = "no")
}

# ---------- Stage 2: 应用复核 ----------
rev_file <- file.path(OUT, "annotation_review.csv")
if (!file.exists(rev_file)) stop("缺少 annotation_review.csv，请先运行 evidence 阶段并复核")
rev <- read.csv(rev_file, check.names = FALSE, fileEncoding = "UTF-8")
stopifnot(all(c("cluster", "celltype_reviewed") %in% colnames(rev)))
rev$cluster <- as.character(rev$cluster)

map_rev <- setNames(rev$celltype_reviewed, rev$cluster)
cl <- as.character(meta$seurat_clusters)
if (any(!cl %in% names(map_rev))) {
  warning("部分簇没有复核映射：", paste(setdiff(unique(cl), names(map_rev)),
                                          collapse = ","))
}
meta$celltype_reviewed <- unname(map_rev[cl])
meta$celltype_reviewed[is.na(meta$celltype_reviewed)] <-
  as.character(meta$celltype[is.na(meta$celltype_reviewed)])
obj$celltype_reviewed <- meta$celltype_reviewed

# 按样本计数（复核后）
tab <- table(meta$sample, meta$celltype_reviewed)
tab_df <- as.data.frame.matrix(tab)
tab_df$sample <- rownames(tab_df)
tab_df$group <- meta$group[match(rownames(tab_df), meta$sample)]
tab_df$total_nuclei <- rowSums(tab)
write.csv(tab_df, file.path(OUT, "per_sample_celltype_counts.csv"),
          row.names = FALSE)
log_msg("per_sample_celltype_counts.csv 已写出")

# Ro/e（复核后）
ct_tab <- table(meta$celltype_reviewed, meta$group)
rown_ct <- rowSums(ct_tab); coln_ct <- colSums(ct_tab); total_ct <- sum(ct_tab)
roe <- (ct_tab / rown_ct) /
  (matrix(coln_ct, nrow = nrow(ct_tab), ncol = ncol(ct_tab),
          byrow = TRUE) / total_ct)
roe[!is.finite(roe)] <- 0
write.csv(as.data.frame.matrix(roe), file.path(OUT, "roe_celltype.csv"))
log_msg("roe_celltype.csv（复核后）已写出")

# 按样本比例统计
ctypes <- colnames(tab_df)[!colnames(tab_df) %in% c("sample", "group",
                                                    "total_nuclei")]
stats <- do.call(rbind, lapply(ctypes, function(ct) {
  vals_sham <- tab_df[[ct]][tab_df$group == "Sham"] /
    tab_df$total_nuclei[tab_df$group == "Sham"]
  vals_tbi <- tab_df[[ct]][tab_df$group == "TBI"] /
    tab_df$total_nuclei[tab_df$group == "TBI"]
  wt <- tryCatch(wilcox.test(vals_tbi, vals_sham, exact = FALSE),
                 error = function(e) NULL)
  data.frame(
    celltype = ct,
    sham_n_samples = sum(tab_df$group == "Sham"),
    tbi_n_samples = sum(tab_df$group == "TBI"),
    sham_total_cells = sum(tab_df[[ct]][tab_df$group == "Sham"]),
    tbi_total_cells = sum(tab_df[[ct]][tab_df$group == "TBI"]),
    sham_median_prop = round(median(vals_sham), 5),
    tbi_median_prop = round(median(vals_tbi), 5),
    wilcox_p = if (!is.null(wt)) round(wt$p.value, 5) else NA,
    stringsAsFactors = FALSE)
}))
stats$wilcox_p_bh <- round(p.adjust(stats$wilcox_p, method = "BH"), 5)
write.csv(stats, file.path(OUT, "celltype_composition_stats.csv"),
          row.names = FALSE)
log_msg("celltype_composition_stats.csv 已写出（Wilcoxon 按样本比例）")

# UMAP（复核后）
Idents(obj) <- "celltype_reviewed"
p_umap <- DimPlot(obj, group.by = "celltype_reviewed", label = TRUE,
                  repel = TRUE, cols = pal_cat) +
  theme_pub(base_size = 12) +
  labs(title = "UMAP (reviewed annotation)")
ggsave_pub(file.path(OUT, "umap_celltype_reviewed.pdf"), p_umap,
       width = 8.5, height = 7, dpi = 300)
log_msg("umap_celltype_reviewed.pdf 已写出")

# 保存复核对象
saveRDS(obj, file.path(OUT, "seurat_object_reviewed.rds"))
log_msg("seurat_object_reviewed.rds 已保存")

# annotation_review.md
md_lines <- c(
  "# 单细胞注释复核记录（06b）",
  "",
  paste0("> 生成时间：", format(Sys.time())),
  "> 数据：GSE209552 人类 snRNA-seq（22391 核 / 17 样本）",
  paste0("> 复核方式：基于每簇 top10 marker（`cluster_top10_markers.csv`）与",
         " 13 类 canonical marker 模块评分（`cluster_evidence.csv`）的书面证据复核；",
         "最终由用户/湿实验同事人工确认。"),
  "",
  "## 1. 逐簇复核结果",
  "",
  "| cluster | 原始注释 | 复核后注释 | 判定 | 证据摘要 |",
  "|---|---|---|---|---|"
)
for (i in seq_len(nrow(rev))) {
  md_lines <- c(md_lines, sprintf(
    "| %s | %s | %s | %s | %s |",
    rev$cluster[i],
    ifelse(is.na(rev$original_celltype[i]), "", rev$original_celltype[i]),
    rev$celltype_reviewed[i],
    ifelse(is.na(rev$action[i]), "", rev$action[i]),
    gsub("\\|", "/",
         ifelse(is.na(rev$evidence[i]), "", rev$evidence[i]))))
}
md_lines <- c(
  md_lines,
  "",
  "## 2. 按样本组成统计（Wilcoxon，样本级比例）",
  "",
  "见 `celltype_composition_stats.csv`；关键结论：",
  "- T_NK：对照 5 份尸检样本共 19 核（每样本 2–5），TBI 组 146 核但集中在少数样本；",
  "  样本级检验受 n=5 vs 12 限制，只能作探索性证据。",
  "- Microglia：TBI 组 1053 核 vs 对照 445 核，但 nr6/nr7/nr8 三个 TBI 样本几乎无小胶质；",
  "  样本内异质性大。",
  "",
  "## 3. 必须写入论文的局限",
  "",
  "1. 对照组为 69/75/87 岁非神经疾病死亡的尸检脑组织（额叶/颞叶），",
  "   TBI 组为平均 49.5±18.2 岁重症 TBI 去骨瓣减压手术切除的挫伤组织（伤后 4h–8d）；",
  "   年龄、取材方式（尸检 vs 手术）、死后间隔、脑区与伤后时间均不可比。",
  "2. 手术挫伤组织可能含血液（BBB 破坏），T/NK 核富集不能排除血液污染；",
  "   尸检对照 T/NK 极低也可能与死后 RNA 降解/细胞丢失有关。",
  "3. 注释基于 canonical marker 模块评分的自动注释 + 本记录中的书面复核，",
  "   未经免疫荧光/流式等湿实验验证。",
  "",
  "## 4. 复核版本记录",
  "",
  paste0("- 复核时间：", format(Sys.time())),
  "- 输入对象：`seurat_object.rds`（自动注释版）",
  "- 输出对象：`seurat_object_reviewed.rds`（含 `celltype_reviewed` 列）",
  "- 复核映射表：`annotation_review.csv`",
  "- 相关日志：`06b_run.log`"
)
writeLines(md_lines, file.path(OUT, "annotation_review.md"))
log_msg("annotation_review.md 已生成")
message("完成。复核后结果见 results/11_SingleCell/。")
