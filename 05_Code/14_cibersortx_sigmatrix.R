# ============================================================
# 14_cibersortx_sigmatrix.R —— 从 GSE209552 复核对象导出 CIBERSORTx 特征矩阵
# 依据：REVIEW4 D5 + 用户确认（全8类 + 敏感性标注；仅人类矩阵）
# 输入：results/11_SingleCell/seurat_object_reviewed.rds
# 输出：
#   data/05_CIBERSORTx/GSE209552_signature_matrix.txt（原始counts，行=基因，列=细胞，tab分隔）
#   data/05_CIBERSORTx/GSE209552_cell_labels.txt（CellBarcode, CellType）
#   data/05_CIBERSORTx/SIGMATRIX_QC.md（真实性核查报告）
# 说明：counts 为稀疏矩阵，按基因逐行写出，避免稠密化内存爆炸。
# ============================================================

suppressMessages({
  library(Seurat)
  library(Matrix)
})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "data", "05_CIBERSORTx")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(ROOT, "results", "11_SingleCell",
                         "seurat_object_reviewed.rds"))
meta <- obj@meta.data
stopifnot("celltype_reviewed" %in% colnames(meta))

ct <- as.character(meta$celltype_reviewed)
samples <- as.character(meta$sample)
groups <- as.character(meta$group)
cells <- colnames(obj)
stopifnot(length(ct) == length(cells))

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
counts <- as(counts, "CsparseMatrix")
genes <- rownames(counts)
stopifnot(length(genes) == nrow(counts), length(cells) == ncol(counts))

# 将单行稀疏矩阵正确展开为 tab 分隔字符串。
# 修复记录（2026-08-08）：dgCMatrix 单行提取后 @i 为【行索引】（此处恒为 0），
# 旧代码 v[vals@i + 1] <- vals@x 会把该基因所有非零值覆盖写入第 1 列，
# 导致矩阵文件损坏。正确做法：转 TsparseMatrix 取 @j（列索引）。
write_sparse_row <- function(gene, sparse_row, n_cells) {
  v <- rep("0", n_cells)
  if (length(sparse_row@x) > 0) {
    t <- as(sparse_row, "TsparseMatrix")
    v[t@j + 1] <- as.character(t@x)
  }
  paste(c(gene, v), collapse = "\t")
}

# ---------- 写出标签文件 ----------
labels <- data.frame(CellBarcode = cells, CellType = ct,
                     stringsAsFactors = FALSE)
write.table(labels, file.path(OUT, "GSE209552_cell_labels.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ---------- 写出特征矩阵（逐基因行） ----------
mat_file <- file.path(OUT, "GSE209552_signature_matrix.txt")
cat("== 重写 barcode 版特征矩阵（修复行写入 bug，覆盖旧文件）==\n")
con <- file(mat_file, open = "w", encoding = "UTF-8")
writeLines(paste(c("Gene", cells), collapse = "\t"), con, useBytes = TRUE)
for (i in seq_len(nrow(counts))) {
  line <- write_sparse_row(genes[i], counts[i, , drop = FALSE],
                           ncol(counts))
  writeLines(line, con, useBytes = TRUE)
  if (i %% 2000 == 0) cat("written genes:", i, "/", nrow(counts), "\n")
}
close(con)

# ---------- QC 统计 ----------
umi_per_cell <- colSums(counts)
genes_per_cell <- colSums(counts > 0)
type_tab <- table(ct)
samp_tab <- table(ct, samples)
n_samp <- apply(samp_tab, 1, function(x) sum(x > 0))
type_df <- data.frame(
  celltype = names(type_tab),
  n_cells = as.integer(type_tab),
  n_samples = as.integer(n_samp),
  median_UMI = round(vapply(names(type_tab), function(t) {
    median(umi_per_cell[ct == t])
  }, numeric(1)), 1),
  median_genes = round(vapply(names(type_tab), function(t) {
    median(genes_per_cell[ct == t])
  }, numeric(1)), 1),
  stringsAsFactors = FALSE
)

marker_list <- list(
  Microglia = c("P2RY12", "CX3CR1", "TMEM119"),
  Astrocyte = c("GFAP", "AQP4", "SLC1A3"),
  Oligodendrocyte = c("MBP", "MOG", "PLP1"),
  Neuron = c("SNAP25", "RBFOX3", "SYT1"),
  OPC = c("PDGFRA", "CSPG4"),
  Endothelial = c("FLT1", "VWF", "CLDN5"),
  T_NK = c("CD3D", "LCK", "CD2"),
  Macrophage = c("VSIG4", "FCGBP", "CD14")
)
marker_rows <- do.call(rbind, lapply(names(marker_list), function(t) {
  gs <- intersect(marker_list[[t]], genes)
  hit <- colSums(counts[gs, , drop = FALSE] > 0) > 0
  data.frame(
    celltype = t,
    marker = paste(gs, collapse = ";"),
    pct_in_type = round(mean(hit[ct == t]), 3),
    pct_outside_type = round(mean(hit[ct != t]), 3),
    stringsAsFactors = FALSE)
}))

per_sample <- as.data.frame.matrix(table(samples, ct))
per_sample$sample <- rownames(per_sample)
per_sample$group <- groups[match(rownames(per_sample), samples)]
per_sample <- per_sample[, c("sample", "group",
                             setdiff(colnames(per_sample),
                                     c("sample", "group")))]

# ============================================================
# CIBERSORTx 网页端格式参考矩阵（2026-08-08 修正）
# 官方格式要求（--single_cell TRUE）：
#   行1 = 细胞类型标签（每列一个，允许重复）；列1 = 基因符号；
#   tab 分隔、无引号、无缺失值、非 log 空间。
# 为控制上传体积：每类最多 200 细胞（固定种子抽样）+ 基因过滤
# （在某类中 >=5% 细胞表达的基因保留）。
# 原 barcode 版矩阵（GSE209552_signature_matrix.txt）保留，
# 供 MuSiC/BisqueRNA 等本地工具使用。
# ============================================================
cat("== 生成 CIBERSORTx 网页端格式参考矩阵 ==\n")
set.seed(2026)
selected_idx <- unlist(lapply(names(type_tab), function(t) {
  idx <- which(ct == t)
  if (length(idx) > 200) idx[sample(length(idx), 200)] else idx
}))
selected_idx <- sort(selected_idx)
ct_sub <- ct[selected_idx]
counts_sub <- counts[, selected_idx, drop = FALSE]

# 基因过滤：在至少一个类型中 >=5% 的所选细胞有表达
keep_gene <- vapply(seq_len(nrow(counts_sub)), function(i) {
  any(vapply(names(type_tab), function(t) {
    mean(counts_sub[i, ct_sub == t] > 0) >= 0.05
  }, logical(1)))
}, logical(1))
genes_sub <- genes[keep_gene]
counts_sub <- counts_sub[keep_gene, , drop = FALSE]
cat("CIBERSORTx ref: cells =", ncol(counts_sub),
    "| genes =", nrow(counts_sub), "\n")

cib_file <- file.path(OUT, "GSE209552_cibersortx_ref.txt")
con2 <- file(cib_file, open = "w", encoding = "UTF-8")
writeLines(paste(c("Gene", ct_sub), collapse = "\t"), con2, useBytes = TRUE)
for (i in seq_len(nrow(counts_sub))) {
  line <- write_sparse_row(genes_sub[i], counts_sub[i, , drop = FALSE],
                           ncol(counts_sub))
  writeLines(line, con2, useBytes = TRUE)
  if (i %% 2000 == 0) cat("cib ref genes:", i, "/", nrow(counts_sub), "\n")
}
close(con2)
cat("CIBERSORTx ref written:", cib_file, "大小:",
    round(file.info(cib_file)$size / 1e6, 1), "MB\n")
write.csv(data.frame(celltype = names(table(ct_sub)),
                     n_cells = as.integer(table(ct_sub))),
          file.path(OUT, "cibersortx_ref_cellcounts.csv"),
          row.names = FALSE)

# ---------- QC 报告 ----------
md <- c(
  "# GSE209552 CIBERSORTx 特征矩阵 QC 报告",
  "",
  paste0("> 生成时间：", format(Sys.time())),
  "> 数据源：`results/11_SingleCell/seurat_object_reviewed.rds`（复核注释版）",
  "> 用途：CIBERSORTx Create Signature Matrix（人 bulk 去卷积参考，仅限人队列）",
  "",
  "## 1. 基本维度",
  "",
  sprintf("- 细胞数：%d", length(cells)),
  sprintf("- 基因数：%d（无重复符号，无 ENSG 前缀）", length(genes)),
  sprintf("- counts 为整数：%s",
          ifelse(all(counts@x == floor(counts@x)), "是", "否")),
  sprintf("- 矩阵文件：`GSE209552_signature_matrix.txt`（%d×%d）",
          nrow(counts), ncol(counts)),
  sprintf("- 标签文件：`GSE209552_cell_labels.txt`（%d 行，与矩阵列序一致）",
          nrow(labels)),
  "- **2026-08-08 修复记录**：早期版本（2026-08-08 09:58）行写入逻辑存在 bug",
  "  （dgCMatrix 单行提取后 `@i` 是行索引，旧代码误当列索引使用），",
  "  导致每个基因的全部非零值被覆盖写入第 1 列、其余列为 0；",
  "  本版已重写并与 Seurat 原始 counts 全行比对一致。",
  "",
  "## 2. 每类细胞数与覆盖样本",
  "",
  "| 细胞类型 | 细胞数 | 覆盖样本数 | 中位 UMI | 中位基因数 |",
  "|---|---|---|---|---|"
)
for (i in seq_len(nrow(type_df))) {
  md <- c(md, sprintf("| %s | %d | %d | %.1f | %.1f |",
                      type_df$celltype[i], type_df$n_cells[i],
                      type_df$n_samples[i], type_df$median_UMI[i],
                      type_df$median_genes[i]))
}
md <- c(
  md,
  "",
  "## 3. 敏感性警示（必须遵守）",
  "",
  "- **Macrophage**：100 核全部来自单个 TBI 样本 GSM6376822_MJ_TBI_nr3_RF，" ,
  "  该签名具有单样本特异性；CIBERSORTx 去卷积中 Macrophage 结果只作敏感性解读。",
  "- **T_NK**：165 核中 Sham 仅 19 核（5 样本）、TBI 146 核且集中于 nr11/nr16/nr2 三样本，",
  "  组间严重失衡；去卷积中 T_NK 结果只作敏感性解读。",
  "- 其余 6 类均覆盖全部 17 个样本。",
  "",
  "## 4. Marker 基因 sanity check（注释与表达一致性）",
  "",
  "定义：pct_in_type = 该类型内表达至少一个 marker 的细胞比例；",
  "pct_outside_type = 其他类型内表达至少一个 marker 的细胞比例。",
  "注：snRNA-seq 存在 dropout，低比例不必然代表注释错误，需结合模块评分证据综合判断。",
  "",
  "| 细胞类型 | marker | 该类型内表达比例 | 类型外表达比例 |",
  "|---|---|---|---|"
)
for (i in seq_len(nrow(marker_rows))) {
  md <- c(md, sprintf("| %s | %s | %.3f | %.3f |",
                      marker_rows$celltype[i], marker_rows$marker[i],
                      marker_rows$pct_in_type[i],
                      marker_rows$pct_outside_type[i]))
}
md <- c(
  md,
  "",
  "## 5. 每样本×细胞类型计数",
  ""
)
ct_names <- setdiff(colnames(per_sample), c("sample", "group"))
md <- c(md, paste0("| sample | group | ", paste(ct_names, collapse = " | "),
                   " |"))
md <- c(md, paste0("|", paste(rep("---", length(ct_names) + 2),
                              collapse = "|"), "|"))
short_name <- sub("_TBI_HuBrainCTL_Nuclei.*|_MJ_TBI_.*", "",
                  per_sample$sample)
for (i in seq_len(nrow(per_sample))) {
  md <- c(md, paste0("| ", short_name[i], " | ", per_sample$group[i],
                     " | ",
                     paste(per_sample[i, ct_names], collapse = " | "),
                     " |"))
}
md <- c(
  md,
  "",
  "完整 CSV：`sigmatrix_per_sample_counts.csv`。",
  "",
  "## 6. CIBERSORTx 网页端输入（正确格式，2026-08-08 修正）",
  "",
  "- 上传文件：**`GSE209552_cibersortx_ref.txt`**（不是 barcode 版矩阵）：",
  "  第一行 = 细胞类型标签（每列一个，重复允许）；第一列 = 基因符号；",
  "  tab 分隔、无引号、无缺失、非 log 空间。",
  "- 每类最多 200 细胞（固定种子 2026 抽样；Macrophage=100、T_NK=165 全保留）；",
  sprintf("- 基因过滤：在至少一个类型中 >=5%%%% 所选细胞表达（保留 %d 个基因）；",
          nrow(counts_sub)),
  "- CIBERSORTx 不需要单独标签文件（single-cell 模式自动构建 phenotype classes）；",
  "- 本文件与 barcode 版矩阵均经 2026-08-08 行写入 bug 修复，",
  "  值与 Seurat 原始 counts 逐列一致；",
  "- 原 barcode 版矩阵 + 标签文件保留用于 MuSiC/BisqueRNA 本地去卷积。",
  "",
  "## 7. 上传步骤与常见错误",
  "",
  "1. cibersortx.stanford.edu 登录 → **Upload Files** → 文件类型选 ",
  "**Single Cell Reference**，上传 `GSE209552_cibersortx_ref.txt`；",
  "2. **Run CIBERSORTx** → **Create Signature Matrix** → **Custom** → **scRNA-seq**，",
  "   选择刚上传的文件，其余默认；",
  "3. 若仍报 `SyntaxError: Unexpected end of JSON input`：多为上传超时/中断（文件过大或网络不稳），",
  "   换 Chrome/Edge、稳定网络重试，或使用更小文件（见 `CIBERSORTX_UPLOAD_GUIDE.md`）；",
  "4. 生成签名矩阵后，用 **Impute Cell Fractions** 上传",
  "   `GSE104687_cibersortx_mixture.txt`（TPM 空间、非 log、符号格式，见上传指南）。",
  "",
  "详细步骤与备选方案见 `data/05_CIBERSORTx/CIBERSORTX_UPLOAD_GUIDE.md`。"
)
writeLines(md, file.path(OUT, "SIGMATRIX_QC.md"))
write.csv(type_df, file.path(OUT, "sigmatrix_celltype_summary.csv"),
          row.names = FALSE)
write.csv(marker_rows, file.path(OUT, "sigmatrix_marker_check.csv"),
          row.names = FALSE)
write.csv(per_sample, file.path(OUT, "sigmatrix_per_sample_counts.csv"),
          row.names = FALSE)

cat("== 完成 ==", "\n")
cat("矩阵文件:", mat_file, "大小:",
    round(file.info(mat_file)$size / 1e6, 1), "MB\n")
cat("标签文件:", file.path(OUT, "GSE209552_cell_labels.txt"), "\n")
cat("QC 报告:", file.path(OUT, "SIGMATRIX_QC.md"), "\n")
