# ============================================================
# 15_cibersortx_mixture.R —— 准备 GSE104687 CIBERSORTx 混合矩阵
# 用途：CIBERSORTx Impute Cell Fractions（人脑 ACT 队列，仅人队列）
# 输入：data/01_GEO_bulk/GSE104687_RIN_corrected_FPKM.csv.gz（原始 FPKM）
# 输出：data/05_CIBERSORTx/GSE104687_cibersortx_mixture.txt
#       行=基因符号，列=样本，tab 分隔，无引号、无缺失、非 log 空间
# 处理依据：
#   1) 原始 FPKM 实测 19230 基因 × 376 样本，无缺失、无重复基因，
#      最大值 18.56（<50）；直接上传会被 CIBERSORTx 误判为 log 空间并按 2^x 反变换。
#   2) 转 TPM：TPM = FPKM / colSums * 1e6（实测最大值约 337，>50，不会被误判）。
#   3) 单细胞签名矩阵为 CPM（CIBERSORTx 自动归一化生成），
#      建议运行时开启 S-mode 批次校正（针对 10x/dropout 型 scRNA 签名）。
# 注意：本文件不改动 data/01_GEO_bulk/GSE104687_expr.csv（04b 外部验证仍用原文件）。
# ============================================================
suppressMessages({library(data.table)})
src <- "data/01_GEO_bulk/GSE104687_RIN_corrected_FPKM.csv.gz"
out <- file.path("data", "05_CIBERSORTx", "GSE104687_cibersortx_mixture.txt")
stopifnot(file.exists(src))

ex_raw <- read.csv(gzfile(src), row.names = 1, check.names = FALSE,
                    stringsAsFactors = FALSE)
ex <- as.matrix(ex_raw)
storage.mode(ex) <- "double"
cat("FPKM dims:", nrow(ex), "genes x", ncol(ex), "samples\n")
stopifnot(!anyNA(ex), !any(is.infinite(ex)))
cat("FPKM max:", max(ex), "| colSums range:", round(min(colSums(ex)), 1),
    "-", round(max(colSums(ex)), 1), "\n")

tpm <- sweep(ex, 2, colSums(ex), "/") * 1e6
cat("TPM max:", round(max(tpm), 2), "| colSums range:",
    round(min(colSums(tpm)), 1), "-", round(max(colSums(tpm)), 1), "\n")

dt <- data.frame(Gene = rownames(tpm), round(tpm, 4), check.names = FALSE)
fwrite(dt, out, sep = "\t", quote = FALSE)
cat("mixture written:", out, "size:",
    round(file.info(out)$size / 1e6, 1), "MB\n")

# 与 CIBERSORTx 参考文件的基因重叠（CIBERSORTx 会给出 <50% 重叠警告）
ref_genes <- fread("data/05_CIBERSORTx/GSE209552_cibersortx_ref.txt",
                   header = TRUE, select = "Gene")[[1]]
ov <- intersect(ref_genes, rownames(ex))
cat("ref genes:", length(ref_genes), "| mixture genes:", nrow(ex),
    "| overlap:", length(ov), sprintf("(%.1f%%)", 100 * length(ov) / length(ref_genes)), "\n")
