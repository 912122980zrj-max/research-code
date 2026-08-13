suppressMessages({library(data.table); library(Seurat); library(Matrix)})

# 用途：修复后 CIBERSORTx 文件的金标准一致性核验（对照 Seurat 原始 counts）
# 1) ref 文件全矩阵（12792 x 1465）逐值比对；
# 2) barcode 版大矩阵比对已停用（GSE209552_signature_matrix.txt 已于 2026-08-08
#    目录整理时删除，1.5GB；CIBERSORTx 网页端以 ref 版为输入）。
# 输出：R-side 检查通过后返回 0；任一不一致即报错退出。
obj <- readRDS("results/11_SingleCell/seurat_object_reviewed.rds")
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
counts <- as(counts, "CsparseMatrix")
sel <- fread("data/05_CIBERSORTx/cibersortx_ref_selection_map.csv")
sel_idx <- sel$full_index
sel_bc <- sel$barcode
stopifnot(ncol(counts) == 22391, length(sel_idx) == 1465)

# ---- FULL ref vs Seurat (12,792 x 1,465) ----
ref <- fread("data/05_CIBERSORTx/GSE209552_cibersortx_ref.txt", header = TRUE)
stopifnot(nrow(ref) == 12792, ncol(ref) == 1466)
ref_genes <- as.character(ref[[1]])
m <- match(ref_genes, rownames(counts))
stopifnot(!anyNA(m), !anyDuplicated(ref_genes))
ref_mat <- as.matrix(ref[, -1, with = FALSE])
true_mat <- as.matrix(counts[m, sel_idx])
d <- which(ref_mat != true_mat, arr.ind = TRUE)
cat("ref full compare | dims:", nrow(ref_mat), "x", ncol(ref_mat),
    "| mismatches:", nrow(d), "| all.equal:", isTRUE(all.equal(ref_mat, true_mat)), "\n")
if (nrow(d) > 0) { print(head(d, 5)); stop("REF MISMATCH") }

cat("ALL R-SIDE CHECKS PASSED\n")
