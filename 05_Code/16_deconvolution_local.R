# ============================================================
# 16_deconvolution_local.R —— 本地去卷积（MuSiC + BisqueRNA）
# 数据：
#   单细胞参考：results/11_SingleCell/seurat_object_reviewed.rds
#                （22391 核，celltype_reviewed 复核注释，原始 counts）
#   bulk 混合：data/05_CIBERSORTx/GSE104687_cibersortx_mixture.txt
#                （GSE104687 ACT 队列，19230 基因 x 376 样本，TPM）
#   分组：data/01_GEO_bulk/GSE104687_pheno.csv（TBI=194, Sham=182）
# 输出（data/05_CIBERSORTx/）：
#   deconvolution_music_8ct.csv / deconvolution_bisque_8ct.csv
#   deconvolution_music_6ct.csv / deconvolution_bisque_6ct.csv
#   pseudobulk_loo_validation.csv（留一队列外验证）
# 说明：CIBERSORTx 网页版 S-mode 无法在本地复现（需 Docker+token），
#   本脚本提供无 token、可复现的替代/敏感性方法。
# ============================================================
suppressMessages({
  library(Seurat); library(Matrix); library(data.table)
  library(SingleCellExperiment); library(Biobase)
  library(MuSiC); library(BisqueRNA)
})

SMOKE <- as.logical(Sys.getenv("DECONV_SMOKE", "FALSE"))

cat("== 1. 读取单细胞参考 ==\n")
obj <- readRDS("results/11_SingleCell/seurat_object_reviewed.rds")
meta <- obj@meta.data
ct <- as.character(meta$celltype_reviewed)
samp <- as.character(meta$sample)
grp <- as.character(meta$group)
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
counts <- as(counts, "CsparseMatrix")
cat("sc dims:", nrow(counts), "x", ncol(counts), "\n")

cat("== 2. 读取 bulk 混合矩阵（TPM）==\n")
bulk <- fread("data/05_CIBERSORTx/GSE104687_cibersortx_mixture.txt",
              header = TRUE)
bulk_genes <- as.character(bulk[[1]])
bulk_mtx <- as.matrix(bulk[, -1, with = FALSE])
rownames(bulk_mtx) <- bulk_genes
rm(bulk)
cat("bulk dims:", nrow(bulk_mtx), "x", ncol(bulk_mtx), "\n")

pheno <- fread("data/01_GEO_bulk/GSE104687_pheno.csv")
ph <- data.frame(sample = as.character(pheno$sample),
                 group = as.character(pheno$group), stringsAsFactors = FALSE)
stopifnot(all(colnames(bulk_mtx) %in% ph$sample))
cat("groups:", paste(names(table(ph$group)), table(ph$group),
                     sep = "=", collapse = " "), "\n")

if (SMOKE) {
  set.seed(1)
  keep_c <- sample(seq_len(ncol(counts)), 800)
  counts <- counts[, keep_c, drop = FALSE]
  ct <- ct[keep_c]; samp <- samp[keep_c]; grp <- grp[keep_c]
  keep_b <- sample(seq_len(ncol(bulk_mtx)), 12)
  bulk_mtx <- bulk_mtx[, keep_b, drop = FALSE]
  cat("SMOKE MODE: sc =", ncol(counts), "cells, bulk =", ncol(bulk_mtx), "samples\n")
}

common <- intersect(rownames(counts), rownames(bulk_mtx))
cat("common genes:", length(common), "\n")
sc_mtx <- counts[common, , drop = FALSE]
bulk_c <- bulk_mtx[common, , drop = FALSE]
stopifnot(length(common) > 10000)

run_music <- function(sc_mat, bulk_mat, cl, sb, max_cell = NULL) {
  if (!is.null(max_cell)) {
    set.seed(2026)
    idx <- unlist(lapply(sort(unique(cl)), function(t) {
      ii <- which(cl == t)
      if (length(ii) > max_cell) ii[sample(length(ii), max_cell)] else ii
    }))
    sc_mat <- sc_mat[, idx, drop = FALSE]
    cl <- cl[idx]; sb <- sb[idx]
  }
  sce <- SingleCellExperiment(assays = list(counts = sc_mat))
  colData(sce)$celltype <- cl
  colData(sce)$sample <- sb
  res <- music_prop(bulk.mtx = bulk_mat, sc.sce = sce,
                    clusters = "celltype", samples = "sample",
                    verbose = FALSE)
  est <- res$Est.prop.weighted
  if (is.null(est)) stop("music_prop returned no Est.prop.weighted")
  est
}

run_bisque <- function(sc_mat, bulk_mat, cl, sb, max_cell = 200) {
  if (!is.null(max_cell)) {
    set.seed(2026)
    idx <- unlist(lapply(sort(unique(cl)), function(t) {
      ii <- which(cl == t)
      if (length(ii) > max_cell) ii[sample(length(ii), max_cell)] else ii
    }))
    sc_mat <- sc_mat[, idx, drop = FALSE]
    cl <- cl[idx]; sb <- sb[idx]
  }
  sc_mat <- as.matrix(sc_mat)  # ExpressionSet 不支持稀疏矩阵；抽样后体积可控
  sc_eset <- ExpressionSet(assayData = sc_mat)
  pData(sc_eset)$cellType <- cl
  pData(sc_eset)$SubjectName <- sb
  bulk_eset <- ExpressionSet(assayData = bulk_mat)
  res <- ReferenceBasedDecomposition(bulk.eset = bulk_eset,
                                     sc.eset = sc_eset, markers = NULL,
                                     cell.types = "cellType",
                                     subject.names = "SubjectName",
                                     use.overlap = FALSE, verbose = FALSE)
  t(res$bulk.props)
}

cat("== 3. MuSiC 8 类 ==\n")
estM8 <- run_music(sc_mtx, bulk_c, ct, samp)
cat("MuSiC 8ct dims:", nrow(estM8), "x", ncol(estM8), "\n")
write.csv(estM8, "data/05_CIBERSORTx/deconvolution_music_8ct.csv")

cat("== 4. BisqueRNA 8 类 ==\n")
estB8 <- run_bisque(sc_mtx, bulk_c, ct, samp)
cat("Bisque 8ct dims:", nrow(estB8), "x", ncol(estB8), "\n")
write.csv(estB8, "data/05_CIBERSORTx/deconvolution_bisque_8ct.csv")

cat("== 5. 6 类敏感性（剔除 Macrophage/T_NK）==\n")
keep6 <- !(ct %in% c("Macrophage", "T_NK"))
estM6 <- run_music(sc_mtx[, keep6, drop = FALSE], bulk_c,
                   ct[keep6], samp[keep6])
write.csv(estM6, "data/05_CIBERSORTx/deconvolution_music_6ct.csv")
estB6 <- run_bisque(sc_mtx[, keep6, drop = FALSE], bulk_c,
                    ct[keep6], samp[keep6])
write.csv(estB6, "data/05_CIBERSORTx/deconvolution_bisque_6ct.csv")
cat("6ct done\n")

cat("== 6. 留一样本外伪 bulk 验证（MuSiC，每类最多 200 细胞）==\n")
if (!SMOKE) {
  samples <- unique(samp)
  tb <- prop.table(table(ct, samp), margin = 2)
  true_prop_m <- matrix(as.numeric(tb), nrow = nrow(tb), ncol = ncol(tb),
                        dimnames = list(rownames(tb), colnames(tb)))
  est_loo <- matrix(NA, nrow = length(samples), ncol = nrow(true_prop_m),
                    dimnames = list(samples, rownames(true_prop_m)))
  pseudo <- sapply(samples, function(s) {
    Matrix::rowSums(counts[, samp == s, drop = FALSE])
  })
  pseudo <- as.matrix(pseudo)
  pseudo_common <- pseudo[common, , drop = FALSE]
  for (s in samples) {
    ref_idx <- which(samp != s)
    e <- run_music(sc_mtx[, ref_idx, drop = FALSE], pseudo_common,
                   ct[ref_idx], samp[ref_idx], max_cell = 200)
    est_loo[s, colnames(e)] <- e[s, ]
  }
  tp_df <- as.data.frame(t(true_prop_m), check.names = FALSE)
  colnames(tp_df) <- paste0("true_", colnames(tp_df))
  es_df <- as.data.frame(est_loo, check.names = FALSE)
  colnames(es_df) <- paste0("est_", colnames(es_df))
  out <- data.frame(sample = samples, tp_df, es_df, check.names = FALSE)
  write.csv(out, "data/05_CIBERSORTx/pseudobulk_loo_validation.csv",
            row.names = FALSE)
  cat("LOO validation written\n")
}

cat("== 完成 ==\n")
