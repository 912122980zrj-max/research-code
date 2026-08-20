# ============================================================
# 01b_download_validation.R —— 补下载外部验证数据（网络恢复后运行）
# 运行：Rscript scripts/01b_download_validation.R [--tries 8] [--sleep 60]
# 功能：
#   1) GSE128543（小鼠 CCI 48h）：带重试下载，并用 Affymetrix gene_assignment
#      把探针映射为基因符号，写入 data/01_GEO_bulk/（供 03/04/04b 使用）
#   2) GSE104687（人脑 ACT 队列）：
#      a. 优先尝试 getGEO 获取样本元数据（TBI/对照分组）
#      b. 自动寻找手动下载的 FPKM 文件
#         （E:\sheng xin\tbi\ 或项目 data/ 下的 GSE104687*FPKM*.csv.gz）
#      c. 找不到分组信息时，生成待填写的 sample_group 模板
# 建议：早上网络好时双击 run_validation.bat 一键执行
# ============================================================

suppressMessages({
  library(GEOquery)
  library(dplyr)
  library(readr)
})

# ---------- 参数与路径 ----------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) as.numeric(args[i + 1]) else default
}
TRIES <- get_arg("--tries", 8)
SLEEP <- get_arg("--sleep", 60)

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv，请把工作目录设为 TBI_pipeline 根目录")
}
OUT <- file.path(ROOT, "data", "01_GEO_bulk")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
options(timeout = 600)

# ---------- 与 01 脚本一致的辅助函数 ----------
extract_symbols <- function(eset) {
  fd <- fData(eset)
  cand <- c("Gene Symbol", "GeneSymbol", "Symbol", "ILMN_Gene",
            "GENE_SYMBOL", "gene_symbol", "GENE", "Gene", "gene_assignment")
  hit <- NULL
  for (nm in cand) {
    if (nm %in% colnames(fd)) { hit <- nm; break }
  }
  if (is.null(hit)) {
    message("  [warn] 未在 fData 中找到基因符号列，返回探针 ID")
    return(rownames(exprs(eset)))
  }
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
  is_ctrl <- grepl("\\b(sham|control|uninjured|vehicle|naive)\\b|notbi|no_tbi|no-tbi", t)
  is_tbi <- grepl("(?<![a-z])tbi(?![a-z])|mtbi|\\bcci\\b|\\binjured\\b|\\binjury\\b|\\binjur\\b|\\bimpact\\b|\\bfp\\b",
                  t, perl = TRUE) |
    grepl("fluid percussion", t)
  ifelse(is_tbi & !is_ctrl, "TBI",
         ifelse(is_ctrl & !is_tbi, "Sham", "Other"))
}

process_eset <- function(eset, gse, subseries = gse) {
  sym <- extract_symbols(eset)
  ex <- exprs(eset)
  if (!is.numeric(ex)) {
    ex <- matrix(as.numeric(ex), nrow = nrow(ex), ncol = ncol(ex),
                 dimnames = dimnames(ex))
  }
  rownames(ex) <- sym
  # 基因内重复探针按求和合并（2026-08-08 修复：不再除样本数列）
  ex <- rowsum(ex, group = rownames(ex), na.rm = TRUE)
  ex[is.na(ex)] <- NA
  ph <- pData(eset)
  ph$geo_accession <- rownames(ph)
  ph$title <- as.character(ph$title)
  ph$group <- infer_group(ph$title)
  ph$cohort <- gse
  ph$subseries <- subseries
  list(expr = ex, pheno = ph)
}

save_bulk <- function(d, gse) {
  saveRDS(d, file.path(OUT, paste0(gse, ".rds")))
  write.csv(d$expr, file.path(OUT, paste0(gse, "_expr.csv")),
            quote = FALSE, fileEncoding = "UTF-8")
  write.csv(d$pheno, file.path(OUT, paste0(gse, "_pheno.csv")),
            fileEncoding = "UTF-8")
  message("  [ok] ", gse, " 已保存：", nrow(d$expr), " 基因 x ",
          ncol(d$expr), " 样本")
}

# ---------- 1) GSE128543（带重试） ----------
download_gse128543 <- function() {
  rds_file <- file.path(OUT, "GSE128543.rds")
  if (file.exists(rds_file)) {
    message("GSE128543 已存在，跳过（删除 data/01_GEO_bulk/GSE128543.* 可强制重下）")
    return(invisible(TRUE))
  }
  message("开始下载 GSE128543（最多重试 ", TRIES, " 次，间隔 ", SLEEP, " 秒）")
  eset_list <- NULL
  for (i in seq_len(TRIES)) {
    ok <- tryCatch({
      g <- getGEO("GSE128543", GSEMatrix = TRUE, getGPL = TRUE)
      eset_list <- g
      TRUE
    }, error = function(e) {
      message("  第 ", i, "/", TRIES, " 次失败：", conditionMessage(e))
      FALSE
    })
    if (ok) break
    if (i < TRIES) {
      message("  ", SLEEP, " 秒后重试...")
      Sys.sleep(SLEEP)
    }
  }
  if (is.null(eset_list)) stop("GSE128543 下载失败（网络仍不稳定），请稍后再试")
  if (!is.list(eset_list)) eset_list <- list(eset_list)
  d <- process_eset(eset_list[[1]], "GSE128543")
  save_bulk(d, "GSE128543")
  print(table(d$pheno$group))
}

# ---------- 2) GSE104687（优先元数据；FPKM 手动文件自动识别） ----------
find_fpkm_file <- function() {
  cand <- c(
    file.path(ROOT, "data", "01_GEO_bulk", "GSE104687_RIN_corrected_FPKM.csv.gz"),
    file.path("E:/sheng xin/tbi", "GSE104687_RIN_corrected_FPKM.csv"),
    file.path(ROOT, "data", "01_GEO_bulk", "GSE104687_RIN_corrected_FPKM.csv.gz"),
    file.path(ROOT, "data", "01_GEO_bulk", "GSE104687_RIN_corrected_FPKM.csv"),
    file.path(ROOT, "data", "01_GEO_bulk", "GSE104687_RIN_corrected_FPKM.csv.gz")
  )
  if (dir.exists("E:/sheng xin/tbi")) {
    extra <- list.files("E:/sheng xin/tbi", pattern = "GSE104687.*FPKM.*csv(\\.gz)?$",
                        recursive = TRUE, full.names = TRUE)
    cand <- c(cand, extra)
  }
  cand <- cand[file.exists(cand)]
  cand <- cand[!grepl("FPKMm", cand, ignore.case = TRUE)]
  # 优先选择带 .gz 的压缩版（若同时存在），其次无后缀版本
  cand[order(!grepl("\\.gz$", cand))]
  unique(cand)
}

fetch_gse104687_meta <- function() {
  # 只取元数据（pData），网络好时成功；失败返回 NULL 不阻塞
  tryCatch({
    g <- getGEO("GSE104687", GSEMatrix = TRUE, getGPL = FALSE)
    es <- if (is.list(g)) g[[1]] else g
    pData(es)
  }, error = function(e) {
    message("  [warn] GSE104687 元数据下载失败（不影响 FPKM 使用）：",
            conditionMessage(e))
    NULL
  })
}

normalize_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  gsub("[^A-Z0-9]", "", x)
}

process_gse104687 <- function() {
  rds_file <- file.path(OUT, "GSE104687.rds")
  if (file.exists(rds_file)) {
    message("GSE104687 已存在，跳过")
    return(invisible(TRUE))
  }
  message("处理 GSE104687（人脑 ACT 队列，人类外部验证）...")
  meta <- fetch_gse104687_meta()
  meta_ph <- NULL
  if (!is.null(meta)) {
    meta_ph <- data.frame(
      sample = rownames(meta),
      geo_accession = rownames(meta),
      title = as.character(meta$title),
      group = infer_group(as.character(meta$title)),
      description = as.character(meta$description),
      stringsAsFactors = FALSE
    )
    write.csv(meta_ph, file.path(OUT, "GSE104687_meta_pheno.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
    message("  [meta] 已获取 GSE104687 元数据并保存 GSE104687_meta_pheno.csv")
  }

  fpkm <- find_fpkm_file()
  if (length(fpkm) == 0) {
    message("  未找到手动下载的 FPKM 文件。")
    message("  请下载 GSE104687_RIN_corrected_FPKM.csv.gz 放到 E:\\sheng xin\\tbi\\ 或 data/ 下，")
    message("  或把 series matrix 放进 data/01_GEO_bulk/ 后重跑本脚本。")
    return(invisible(FALSE))
  }
  message("  找到 FPKM 文件：", fpkm[1])
  ex_raw <- read.csv(gzfile(fpkm[1]), row.names = 1, check.names = FALSE,
                     stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  ex <- as.matrix(ex_raw)
  storage.mode(ex) <- "double"
  ex <- rowsum(ex, group = rownames(ex), na.rm = TRUE)
  ex <- ex / matrix(rowSums(!is.na(ex)), nrow = nrow(ex), ncol = ncol(ex),
                    byrow = FALSE)
  ex[is.na(ex)] <- NA
  message("  FPKM 矩阵：", nrow(ex), " 基因 x ", ncol(ex), " 样本")

  # 构建样本表：优先用 getGEO 元数据匹配，否则生成待填模板
  sample_ids <- colnames(ex)
  ph <- NULL
  if (!is.null(meta_ph)) {
    # ACT 数字样本 ID 存在 description 列（如 496100278）
    desc <- meta_ph$description
    desc[is.na(desc) | desc == ""] <- meta_ph$geo_accession[is.na(desc) | desc == ""]
    meta_ph$key <- normalize_id(desc)
    keys <- normalize_id(sample_ids)
    m <- match(keys, meta_ph$key)
    hit <- sum(!is.na(m)) / length(keys)
    message("  样本与元数据（description 列）匹配率：", round(hit * 100, 1), "%")
    if (hit > 0.5) {
      ph <- data.frame(
        sample = sample_ids,
        geo_accession = sample_ids,
        title = ifelse(is.na(m), sample_ids, meta_ph$title[m]),
        group = ifelse(is.na(m), "Other", meta_ph$group[m]),
        stringsAsFactors = FALSE
      )
    }
  }
  if (is.null(ph)) {
    # 尝试用户提供的 sample_group.csv（两列：sample,group）
    sg_file <- file.path(OUT, "GSE104687_sample_group.csv")
    if (file.exists(sg_file)) {
      sg <- read.csv(sg_file, fileEncoding = "UTF-8")
      sg$sample <- as.character(sg$sample)
      ph <- data.frame(sample = sample_ids,
                       title = sample_ids,
                       group = sg$group[match(normalize_id(sample_ids),
                                              normalize_id(sg$sample))],
                       stringsAsFactors = FALSE)
      ph$group[is.na(ph$group)] <- "Other"
    } else {
      write.csv(data.frame(sample = sample_ids, group = ""),
                sg_file, row.names = FALSE, fileEncoding = "UTF-8")
      message("  [注意] 无法自动获取分组，已生成模板：")
      message("  ", sg_file)
      message("  请按 ACT 队列信息填写 group（TBI / Sham），保存后重跑本脚本")
      return(invisible(FALSE))
    }
  }
  ph$cohort <- "GSE104687"
  d <- list(expr = ex, pheno = ph)
  save_bulk(d, "GSE104687")
  print(table(ph$group))
}

# ---------- 主流程 ----------
cat("===== 01b 外部验证数据补下载 =====", "\n")
download_gse128543()
process_gse104687()
cat("完成。接下来运行 04b_external_validation.R 或用 run_validation.bat 一键执行。\n")
