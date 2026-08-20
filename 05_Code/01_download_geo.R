# ============================================================
# 01_download_geo.R —— 下载批量转录组数据（GEO）
# 运行：Rscript scripts/01_download_geo.R
# 说明：
#   - 默认下载 GSE58485 / GSE41345 / GSE128543（发现+验证）；
#     GSE2871（大鼠，可选）默认跳过，把 config 里那行前面的 # 去掉即可下载。
#   - 输出到 data/01_GEO_bulk/：每个数据集一个 .csv 表达矩阵 + 样本表 + RDS。
#   - 样本分组通过标题文本自动推断（Sham/Control vs TBI/CCI/mTBI），
#     推断不准时请手工修改 data/01_GEO_bulk/*_pheno_manual.csv 并重跑本脚本。
# ============================================================

suppressMessages({
  library(GEOquery)
  library(dplyr)
  library(readr)
  library(stringr)
})

# ---------- 路径 ----------
if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv，请把工作目录设为 TBI_pipeline 根目录")
}
OUT <- file.path(ROOT, "data", "01_GEO_bulk")
RES <- file.path(ROOT, "results", "01_GEO_Download")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

# ---------- 读取数据集清单 ----------
meta <- read.delim(file.path(ROOT, "config", "datasets.tsv"),
                   comment.char = "#", header = FALSE,
                   col.names = c("gse", "species", "tissue", "platform",
                                 "n", "role", "note"), stringsAsFactors = FALSE)
meta$gse <- trimws(meta$gse)
# 可选下载开关：只有 role 属于需要下载的才下载（singlecell 不在此下载）
download_roles <- c("discovery", "validation", "optional")
meta <- meta[meta$role %in% download_roles, ]
meta$skip <- grepl("GSE209552|GSE160763|GSE269748|GSE254162|GSE254163|GSE254164|GSE160651", meta$gse)

# ---------- 辅助函数 ----------
extract_symbols <- function(eset) {
  fd <- fData(eset)
  cand <- c("Gene Symbol", "GeneSymbol", "Symbol", "ILMN_Gene",
            "GENE_SYMBOL", "gene_symbol", "GENE", "Gene", "gene_assignment")
  hit <- NULL
  for (nm in cand) {
    if (nm %in% colnames(fd)) { hit <- nm; break }
  }
  if (is.null(hit)) {
    message("  [warn] 未在 fData 中找到基因符号列，尝试平台注释映射")
    return(rownames(exprs(eset)))
  }
  # Affymetrix ST 平台的 gene_assignment 形如 "NM_021297 // Tlr4 // 21898 ..."
  if (hit == "gene_assignment") {
    raw <- as.character(fd[[hit]])
    sym <- vapply(strsplit(raw, "///"), function(x) {
      parts <- trimws(strsplit(x[1], "//")[[1]])
      if (length(parts) >= 2 && !is.na(parts[2]) && parts[2] != "---") {
        parts[2]
      } else {
        NA_character_
      }
    }, character(1))
    sym[is.na(sym) | sym == ""] <- rownames(exprs(eset))[is.na(sym) | sym == ""]
    return(sym)
  }
  sym <- as.character(fd[[hit]])
  sym[is.na(sym) | sym == ""] <- rownames(exprs(eset))[is.na(sym) | sym == ""]
  sym
}

infer_group <- function(titles) {
  t <- tolower(titles)
  is_ctrl <- str_detect(t, "\\b(sham|control|uninjured|vehicle|naive)\\b")
  is_tbi <- str_detect(t, "\\b(tbi|mtbi|cci|injured|injury|injur|impact|fp)\\b") |
    str_detect(t, "fluid percussion")
  grp <- ifelse(is_tbi & !is_ctrl, "TBI",
         ifelse(is_ctrl & !is_tbi, "Sham", "Other"))
  grp
}

download_one <- function(gse, force = FALSE) {
  expr_file <- file.path(OUT, paste0(gse, "_expr.csv"))
  pheno_file <- file.path(OUT, paste0(gse, "_pheno.csv"))
  rds_file <- file.path(OUT, paste0(gse, ".rds"))
  if (file.exists(expr_file) && !force) {
    if (file.exists(rds_file)) {
      d <- readRDS(rds_file)
      ph <- d$pheno
      ph$title <- as.character(ph$title)
      ph$group <- infer_group(ph$title)
      write.csv(ph, pheno_file, fileEncoding = "UTF-8")
      message("  [repair] ", gse, " 已用 RDS 重写分组表（UTF-8）")
    } else {
      message("  [skip] ", gse, " 已存在（无 RDS 可修复）")
    }
    return(invisible(NULL))
  }
  message("  [download] ", gse)
  gse_list <- getGEO(gse, GSEMatrix = TRUE, getGPL = TRUE)
  if (!is.list(gse_list)) gse_list <- list(gse_list)

  all_expr <- list(); all_pheno <- list()
  for (i in seq_along(gse_list)) {
    eset <- gse_list[[i]]
    sym <- extract_symbols(eset)
    ex <- exprs(eset)
    if (!is.numeric(ex)) {
      ex <- matrix(as.numeric(ex), nrow = nrow(ex), ncol = ncol(ex),
                   dimnames = dimnames(ex))
    }
    rownames(ex) <- sym
    # 基因内重复探针按求和合并（2026-08-08 修复：原除以样本数列的写法
    # 会整体缩放表达值；本行只做 rowsum 求和，不除样本数）
    ex <- rowsum(ex, group = rownames(ex), na.rm = TRUE)
    ex[is.na(ex)] <- NA
    ph <- pData(eset)
    ph$geo_accession <- rownames(ph)
    ph$title <- as.character(ph$title)
    ph$group <- infer_group(ph$title)
    ph$cohort <- gse
    ph$subseries <- if (length(gse_list) > 1) paste0(gse, ".", i) else gse
    all_expr[[i]] <- ex
    all_pheno[[i]] <- ph
  }

  if (length(all_expr) > 1) {
    # 超级系列（如 GSE58485 的两个子系列）：取共同基因合并
    common <- Reduce(intersect, lapply(all_expr, rownames))
    ex_merged <- do.call(cbind, lapply(all_expr, function(m) m[common, , drop = FALSE]))
    ph_merged <- do.call(rbind, all_pheno)
    saveRDS(list(expr = ex_merged, pheno = ph_merged), file.path(OUT, paste0(gse, ".rds")))
    write.csv(ex_merged, expr_file, quote = FALSE, fileEncoding = "UTF-8")
    write.csv(ph_merged, pheno_file, fileEncoding = "UTF-8")
    message("  [ok] ", gse, " 合并子系列后样本数 = ", ncol(ex_merged))
  } else {
    ex <- all_expr[[1]]; ph <- all_pheno[[1]]
    saveRDS(list(expr = ex, pheno = ph), file.path(OUT, paste0(gse, ".rds")))
    write.csv(ex, expr_file, quote = FALSE, fileEncoding = "UTF-8")
    write.csv(ph, pheno_file, fileEncoding = "UTF-8")
    message("  [ok] ", gse, " 样本数 = ", ncol(ex), " 基因数 = ", nrow(ex))
  }
}

# ---------- 主流程 ----------
for (i in seq_len(nrow(meta))) {
  if (meta$skip[i]) {
    message("[skip] ", meta$gse[i], "（可选/单细胞队列，本脚本不下载）")
    next
  }
  tryCatch(download_one(meta$gse[i]),
           error = function(e) message("  [ERROR] ", meta$gse[i], ": ", conditionMessage(e)))
}

# ---------- 摘要报告 ----------
report <- file.path(RES, "01_download_report.txt")
con <- file(report, "w", encoding = "UTF-8")
writeLines(c("GEO bulk download report", paste("generated:", Sys.time()), ""), con)
for (f in list.files(OUT, pattern = "_expr.csv$")) {
  m <- read.csv(file.path(OUT, f), row.names = 1, check.names = FALSE,
                fileEncoding = "UTF-8")
  ph <- read.csv(file.path(OUT, sub("_expr.csv", "_pheno.csv", f)),
                 row.names = 1, check.names = FALSE, fileEncoding = "UTF-8")
  writeLines(sprintf("%s : %d genes x %d samples ; groups: %s",
                     sub("_expr.csv", "", f), nrow(m), ncol(m),
                     paste(names(table(ph$group)), table(ph$group), sep = "=", collapse = ",")),
             con)
}
close(con)
message("报告已写入: ", report)
message("完成。请打开 data/01_GEO_bulk/*_pheno.csv 核对分组列 group，"
        , "如果推断有误，请手工修改该文件后重跑本脚本（或直接在第 3 步用自定义分组）。")
