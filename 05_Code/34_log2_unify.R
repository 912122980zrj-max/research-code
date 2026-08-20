# ============================================================
# 34_log2_unify.R —— M1：bulk 表达尺度统一为 log2
# 依据：REVIEW_RESPONSE_FINAL_REVIEW.md M1（P0）
# 尺度核查（2026-08-09）：
#   GSE58485  线性（GEO 原始中位≈93）→ log2(x+1)
#   GSE41345  log2（中位≈0，含负值）→ 保持
#   GSE128543 log2（中位≈5.8）→ 保持
#   GSE104687 原始 RIN_corrected_FPKM 文件已是 log 尺度
#             （中位 3.13、最大 12.8，与存储一致）→ 保持，不二次变换
# 输出：GSE58485 expr（CSV+RDS）原位 log2(x+1)；线性版备份于 _pre_log2/
# ============================================================

OUT <- "data/01_GEO_bulk"
BAK <- file.path(OUT, "_pre_log2")
dir.create(BAK, showWarnings = FALSE)

gse <- "GSE58485"
expr_csv <- file.path(OUT, paste0(gse, "_expr.csv"))
rds_file <- file.path(OUT, paste0(gse, ".rds"))

if (!file.exists(file.path(BAK, paste0(gse, "_expr_linear_fixed.csv")))) {
  file.copy(expr_csv, file.path(BAK, paste0(gse, "_expr_linear_fixed.csv")))
  file.copy(rds_file, file.path(BAK, paste0(gse, "_rds_linear_fixed.rds")))
}

ex <- as.matrix(read.csv(expr_csv, row.names = 1, check.names = FALSE,
                         fileEncoding = "UTF-8"))
cat(gse, "before log2: median", median(ex), " max", max(ex), "\n")
ex_new <- log2(ex + 1)
write.csv(ex_new, expr_csv, quote = FALSE, fileEncoding = "UTF-8")
d <- readRDS(rds_file)
d$expr <- ex_new
saveRDS(d, rds_file)
cat(gse, "after log2: median", round(median(ex_new), 4),
    " max", round(max(ex_new), 4), "\n")

# GSE104687 校验：当前存储 == 原始 RIN_corrected_FPKM（log 尺度）
ex104 <- as.matrix(read.csv(file.path(OUT, "GSE104687_expr.csv"),
                            row.names = 1, check.names = FALSE,
                            fileEncoding = "UTF-8"))
cat("GSE104687 kept as-is: median", round(median(ex104), 4),
    " max", round(max(ex104), 4), "\n")

md <- c(
  "# log2 尺度统一记录（LOG2_UNIFY）",
  "",
  paste("> 时间：", Sys.time()),
  "> 依据：REVIEW_RESPONSE_FINAL_REVIEW.md M1",
  "- GSE58485：线性（GEO 原始中位≈93）→ log2(x+1)，CSV/RDS 原位更新；",
  "  线性版备份于 `_pre_log2/GSE58485_expr_linear_fixed.csv` 与 RDS；",
  "- GSE41345 / GSE128543：原始即为 log2 尺度，保持；",
  "- GSE104687：原始 RIN_corrected_FPKM 文件经核对已是 log 尺度",
  "  （中位 3.13、最大 12.8，与存储完全一致），保持，不二次变换；",
  "- GSE205958：验证时用 log2(FPKM+1) 与训练尺度对齐（21 脚本已实现）。",
  "- 影响：DEG 阈值 |log2FC|>0.585 在统一 log2 尺度上恢复生物学含义；",
  "  03→04/04b/19/20/21/26/08/09/10/31 全链重跑后的数字为最终口径。"
)
writeLines(md, file.path(OUT, "LOG2_UNIFY.md"))
message("log2 统一完成。备份在 ", BAK)
