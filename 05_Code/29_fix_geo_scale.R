# ============================================================
# 29_fix_geo_scale.R —— 修复 01/01b 下载脚本的表达尺度错误
# 问题（2026-08-08 终审 ML 核查发现）：
#   01/01b 在 rowsum 后除以 rowSums(!is.na())（=样本数列），
#   导致每个数据集的表达矩阵整体缩小了 1/n_samples；
#   已验证四套矩阵均无 NA，故 stored = 原始探针求和 / ncol，
#   乘以 ncol 可精确恢复原始求和尺度。
# 处理：CSV 与 RDS 原位修正；旧文件备份到 _pre_scale_fix/。
# 依据：数据真实性核查 + FINAL_REVIEW_PUBLISHABILITY.md（写作数字须真实）。
# ============================================================

OUT <- "data/01_GEO_bulk"
BAK <- file.path(OUT, "_pre_scale_fix")
dir.create(BAK, showWarnings = FALSE)

datasets <- c(GSE58485 = 29, GSE41345 = 18, GSE128543 = 28,
              GSE104687 = 376)

for (gse in names(datasets)) {
  n <- datasets[[gse]]
  expr_csv <- file.path(OUT, paste0(gse, "_expr.csv"))
  rds_file <- file.path(OUT, paste0(gse, ".rds"))

  # 备份旧文件（可追溯）
  if (!file.exists(file.path(BAK, paste0(gse, "_expr.csv")))) {
    file.copy(expr_csv, file.path(BAK, paste0(gse, "_expr.csv")))
    if (file.exists(rds_file)) {
      file.copy(rds_file, file.path(BAK, paste0(gse, ".rds")))
    }
  }

  ex <- as.matrix(read.csv(expr_csv, row.names = 1, check.names = FALSE,
                           fileEncoding = "UTF-8"))
  stopifnot(ncol(ex) == n)
  if (anyNA(ex)) stop(paste(gse, "存在 NA，无法精确反推"))
  ex_new <- ex * n
  write.csv(ex_new, expr_csv, quote = FALSE, fileEncoding = "UTF-8")

  if (file.exists(rds_file)) {
    d <- readRDS(rds_file)
    d$expr <- ex_new
    saveRDS(d, rds_file)
  }
  cat(gse, ":", ncol(ex), "samples; scale corrected x", n, "\n")
}

md <- c(
  "# 表达尺度修复记录（SCALE_REPAIR）",
  "",
  paste("> 修复时间：", Sys.time()),
  "> 问题：01/01b 下载脚本在 rowsum 后除以样本数列，使各数据集表达整体缩小 1/n；",
  "> 验证：四套矩阵无 NA，stored = 原始探针求和 / ncol；",
  "> 处理：全部乘回 ncol（GSE58485×29、GSE41345×18、GSE128543×28、GSE104687×376），",
  "  恢复原始探针求和尺度；CSV 与 RDS 原位更新；",
  "> 备份：`_pre_scale_fix/`（旧 CSV + RDS，可恢复）；",
  "> 脚本修复：01/01b 已改为 rowsum 求和、不再除样本数（未来重下不再复现）；",
  "> 影响：DEG/WGCNA/ML 等发现链与外部验证按新尺度重跑（见 RESPONSE_FINAL_REVIEW.md）。"
)
writeLines(md, file.path(OUT, "SCALE_REPAIR.md"))
message("尺度修复完成，备份在 ", BAK)
