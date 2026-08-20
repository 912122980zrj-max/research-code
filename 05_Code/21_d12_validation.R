# ============================================================
# 21_d12_validation.R —— D12 外部验证扩充：GSE205958（小鼠 CCI 皮层，7 d，RNA-seq）
# 运行：Rscript scripts/21_d12_validation.R
# 前置：04b 已生成 frozen_coefficients.csv；GSE205958 FPKM xlsx 与系列矩阵已下载
# 说明（2026-08-08，用户授权选定 D12 队列）：
#   - 主口径：固定 CLEAN_3gene_only 冻结系数（P2RY6/HSD11B1/CCR5，v4 D4 top3），
#     log2(FPKM+1)（与训练 log2 尺度一致），
#     在 GSE205958（TBI 5 vs Sham/naive 5，n=10）上评估 AUC + Delong 95% CI +
#     训练阈值灵敏度/特异度；不重训、不插补缺失基因。
#   - 敏感性：CLEAN_5gene / CLEAN_7gene / FULL_5gene 冻结模型；
#     原始线性 FPKM（跨尺度，仅排序参考）；GSE128543 WT-only 方法学互查。
#   - 如实标注：RNA-seq FPKM 与发现队列微阵列表达为跨平台线性尺度近似，
#     7 d 亚急性时点与发现队列（3 d）及 GSE128543（48 h）不同。
# 输出：results/08_External_Validation/D12_GSE205958/
#   external_validation_v3.csv、EXPERIMENT_REPORT.md、frozen_roc_GSE205958.pdf
# ============================================================

suppressMessages({
  library(pROC)
  library(dplyr)
  library(readr)
  library(readxl)
  library(ggplot2)
})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "results", "08_External_Validation")
OUT_D12 <- file.path(OUT, "D12_GSE205958")
dir.create(OUT_D12, recursive = TRUE, showWarnings = FALSE)
theme_pub_file <- file.path(ROOT, "scripts", "theme_pub.R")
if (file.exists(theme_pub_file)) source(theme_pub_file)

# ---------- 1. 冻结系数 ----------
coef_file <- file.path(OUT, "frozen_coefficients.csv")
coefs <- read.csv(coef_file, stringsAsFactors = FALSE)
models_all <- unique(coefs$model)
cat("frozen models:", paste(models_all, collapse = "; "), "\n")

# ---------- 2. 发现队列训练数据（用于训练阈值） ----------
read_bulk <- function(gse) {
  rds_file <- file.path(ROOT, "data", "01_GEO_bulk", paste0(gse, ".rds"))
  if (file.exists(rds_file)) return(readRDS(rds_file))
  ef <- file.path(ROOT, "data", "01_GEO_bulk", paste0(gse, "_expr.csv"))
  pf <- file.path(ROOT, "data", "01_GEO_bulk", paste0(gse, "_pheno.csv"))
  if (!file.exists(ef)) return(NULL)
  list(expr = as.matrix(read.csv(ef, row.names = 1, check.names = FALSE,
                                 fileEncoding = "UTF-8")),
       pheno = read.csv(pf, row.names = 1, check.names = FALSE,
                        fileEncoding = "UTF-8"))
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

build_df <- function(gse, genes, clean = FALSE) {
  d <- read_bulk(gse)
  if (is.null(d)) return(NULL)
  ex <- d$expr
  rn <- toupper(rownames(ex))
  ex2 <- ex[rn %in% genes, , drop = FALSE]
  rn2 <- rn[rn %in% genes]
  if (anyDuplicated(rn2)) {
    ex2 <- rowsum(ex2, group = rn2)
  } else {
    rownames(ex2) <- rn2
  }
  ex2 <- ex2[order(rownames(ex2)), , drop = FALSE]
  grp <- infer_group(as.character(d$pheno$title))
  if (clean) {
    excl <- grepl("cyclophosphamide|ex-4-mtbi", d$pheno$title,
                  ignore.case = TRUE)
    grp[excl] <- "Other"
  }
  df <- as.data.frame(t(ex2), check.names = FALSE)
  df$group <- factor(ifelse(grp == "TBI", "TBI", "Sham"),
                     levels = c("Sham", "TBI"))
  df <- df[grp %in% c("TBI", "Sham"), ]
  attr(df, "features") <- intersect(rownames(ex2), genes)
  df
}

core_full_file <- file.path(ROOT, "results", "05_Candidate_Intersection",
                            "candidate_genes_FULL.csv")
core13 <- if (file.exists(core_full_file)) {
  toupper(read.csv(core_full_file)$gene)
} else {
  c("CYP1B1", "HSD11B1", "MMP2", "P2RY6", "PGF", "RIPK1")
}
sig3 <- c("P2RY6", "HSD11B1", "CCR5")  # v4 重跑后 D4 SHAP top3
cand_file <- file.path(ROOT, "results", "05_Candidate_Intersection",
                       "candidate_genes.csv")
cand28 <- if (file.exists(cand_file)) toupper(read.csv(cand_file)$gene) else NULL

disc_full13 <- rbind(build_df("GSE58485", core13, clean = FALSE),
                     build_df("GSE41345", core13, clean = FALSE))
disc_clean13 <- rbind(build_df("GSE58485", core13, clean = TRUE),
                      build_df("GSE41345", core13, clean = TRUE))
disc_clean28 <- if (!is.null(cand28)) {
  rbind(build_df("GSE58485", cand28, clean = TRUE),
        build_df("GSE41345", cand28, clean = TRUE))
} else NULL
disc_sig3 <- rbind(build_df("GSE58485", sig3, clean = TRUE),
                   build_df("GSE41345", sig3, clean = TRUE))

tr_data <- list(
  FULL_5gene = disc_full13,
  CLEAN_5gene = disc_clean13,
  CLEAN_3gene_only = disc_sig3,
  CLEAN_7gene = disc_clean28
)

# ---------- 3. GSE205958 数据 ----------
xlsx_candidates <- c(
  file.path(ROOT, "data", "01_GEO_bulk", "GSE205958_genes_fpkm_expression.xlsx"),
  file.path(ROOT, "data", "01_GEO_bulk", "GSE205958_fpkm_dl3.xlsx"),
  file.path(ROOT, "data", "01_GEO_bulk", "GSE205958_fpkm_dl.xlsx")
)
xlsx_file <- xlsx_candidates[file.exists(xlsx_candidates)][1]
if (is.na(xlsx_file)) stop("GSE205958 xlsx 未找到")

sheets <- excel_sheets(xlsx_file)
cat("xlsx sheets:", paste(sheets, collapse = "; "), "\n")
fpkm_sheet <- sheets[grepl("fpkm", sheets, ignore.case = TRUE)]
sheet <- if (length(fpkm_sheet) > 0) fpkm_sheet[1] else sheets[1]

# 通用读取：先看列名是否为两行表头（col_names=FALSE 模式自行判断）
raw_head <- read_excel(xlsx_file, sheet = sheet, n_max = 5,
                       col_names = FALSE)
cat("raw head cols:", ncol(raw_head), " first row:", 
    paste(as.character(raw_head[1, ][1:min(8, ncol(raw_head))]), collapse = " | "),
    "\n")

# 采用第一行作列名读取（若 xlsx 为两行表头，读取后按基因列+样本列清洗）
dat <- as.data.frame(read_excel(xlsx_file, sheet = sheet),
                     stringsAsFactors = FALSE, check.names = FALSE)
cat("xlsx dim:", nrow(dat), "x", ncol(dat), "\n")

# 识别基因符号列（优先 gene_name/symbol，避免误选 gene_id）
cn <- colnames(dat)
gene_col <- which(tolower(cn) %in% c("gene_name", "symbol", "genesymbol",
                                     "gene symbol"))[1]
if (is.na(gene_col)) {
  gene_col <- which(grepl("gene_name|symbol", cn, ignore.case = TRUE))[1]
}
if (is.na(gene_col)) {
  gene_col <- which(grepl("^gene$", cn, ignore.case = TRUE))[1]
}
if (is.na(gene_col)) gene_col <- 1
cat("gene column:", gene_col, "->", cn[gene_col], "\n")

gene_vec <- as.character(dat[[gene_col]])
gene_vec <- toupper(trimws(gene_vec))
# 表达列：FPKM.*（count.* 列不用于本验证，因需文库大小校正）
sample_cols <- which(grepl("^FPKM\\.", cn, ignore.case = TRUE))
sample_names <- sub("^FPKM\\.", "", cn[sample_cols])
cat("sample columns:", length(sample_names),
    paste(head(sample_names, 12), collapse = "; "), "\n")

expr <- sapply(sample_cols, function(j) {
  v <- suppressWarnings(as.numeric(dat[[j]]))
  v
})
rownames(expr) <- gene_vec
expr <- expr[!is.na(gene_vec) & gene_vec != "", , drop = FALSE]
# 基因内重复行取均值（修复：按重复行数除，而非样本数）
if (anyDuplicated(rownames(expr))) {
  dup_n <- table(rownames(expr))
  expr <- rowsum(expr, group = rownames(expr), na.rm = TRUE)
  expr <- expr / matrix(as.numeric(dup_n[rownames(expr)]),
                        nrow = nrow(expr), ncol = ncol(expr),
                        byrow = FALSE)
}
expr[is.na(expr)] <- 0
cat("expr matrix:", nrow(expr), "x", ncol(expr), "\n")

# 样本分组：从系列矩阵解析
series_gz <- file.path(ROOT, "data", "01_GEO_bulk",
                       "GSE205958_series_matrix.txt.gz")
parse_series <- function(gzfile) {
  con <- gzfile(gzfile, "rt", encoding = "UTF-8")
  lines <- readLines(con, warn = FALSE)
  close(con)
  get_vals <- function(key) {
    ln <- lines[grepl(paste0("^!", key, "\t"), lines)]
    if (length(ln) == 0) return(NULL)
    vals <- sub(paste0("^!", key, "\t"), "", ln)
    vals <- strsplit(vals, "\t")[[1]]
    gsub('^"|"$', "", trimws(vals))
  }
  get_treatment <- function() {
    ln <- lines[grepl("^!Sample_characteristics_ch1\t", lines)]
    vals <- lapply(ln, function(l) {
      v <- strsplit(sub("^!Sample_characteristics_ch1\t", "", l), "\t")[[1]]
      gsub('^"|"$', "", trimws(v))
    })
    for (v in vals) {
      if (any(grepl("^treatment:", v, ignore.case = TRUE))) {
        return(sub("^treatment:\\s*", "",
                   v[grepl("^treatment:", v, ignore.case = TRUE)],
                   ignore.case = TRUE))
      }
    }
    NULL
  }
  list(title = get_vals("Sample_title"),
       geo = get_vals("Sample_geo_accession"),
       treatment = get_treatment())
}
sm <- parse_series(series_gz)
if (is.null(sm$title)) stop("系列矩阵解析失败")
grp_from_title <- if (!is.null(sm$treatment)) {
  g <- ifelse(tolower(sm$treatment) %in% c("tbi", "traumatic brain injury"),
              "TBI", "Sham")
  g
} else {
  infer_group(sm$title)
}
cat("series matrix groups:", paste(names(table(grp_from_title)),
                                   table(grp_from_title), sep = "=",
                                   collapse = ", "), "\n")

# 将 xlsx 列名映射到 series 样本（支持 GSM 或 TBI1/Sham1 标签）
norm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
xlsx_norm <- norm(sample_names)
title_norm <- norm(sm$title)
geo_norm <- norm(sm$geo)
map_idx <- vapply(xlsx_norm, function(x) {
  hit <- which(geo_norm == x)
  if (length(hit) == 0) hit <- which(title_norm == x |
                                       endsWith(title_norm, x))
  if (length(hit) == 0) NA_integer_ else hit[1]
}, integer(1))
if (anyNA(map_idx)) {
  cat("unmapped sample columns:", paste(sample_names[is.na(map_idx)],
                                        collapse = "; "), "\n")
  keep <- !is.na(map_idx)
} else {
  keep <- rep(TRUE, length(sample_names))
}
expr <- expr[, keep, drop = FALSE]
sample_names <- sample_names[keep]
map_idx <- map_idx[keep]
grp <- grp_from_title[map_idx]
df_val <- as.data.frame(t(expr), check.names = FALSE)
df_val$group <- factor(grp, levels = c("Sham", "TBI"))
df_val <- df_val[df_val$group %in% c("Sham", "TBI"), ]
cat("validation samples:", nrow(df_val), " TBI:", sum(df_val$group == "TBI"),
    " Sham:", sum(df_val$group == "Sham"), "\n")

# ---------- 4. 冻结模型评估 ----------
eval_frozen <- function(model_label, df, tr_df, scale_label = "fpkm_raw") {
  cf <- coefs[coefs$model == model_label, ]
  genes <- cf$term[cf$term != "(Intercept)"]
  missing <- setdiff(genes, colnames(df))
  if (length(missing) > 0) {
    return(data.frame(
      model = model_label, scale = scale_label, n = nrow(df),
      n_tbi = sum(df$group == "TBI"), n_sham = sum(df$group == "Sham"),
      auc = NA, ci_low = NA, ci_high = NA,
      concordant_pairs = NA, pairs_total = NA,
      exact_ci_low = NA, exact_ci_high = NA, threshold = NA,
      sensitivity = NA, specificity = NA,
      missing_genes = paste(missing, collapse = ";"),
      note = "not evaluable: signature gene(s) absent in cohort", 
      stringsAsFactors = FALSE))
  }
  lp <- cf$coefficient[cf$term == "(Intercept)"]
  for (g in genes) {
    lp <- lp + cf$coefficient[cf$term == g] * df[[g]]
  }
  prob <- plogis(lp)
  roc_v <- roc(as.numeric(df$group == "TBI"), prob, quiet = TRUE)
  a <- as.numeric(auc(roc_v))
  ci <- tryCatch(as.numeric(ci.auc(roc_v, method = "delong")),
                 error = function(e) rep(NA_real_, 3))
  # AUC=1 时 Delong/bootstrap CI 退化；改报配对一致性比例的精确二项
  # （Clopper–Pearson）CI（2026-08-09 v3.3：原 prop.test 为带连续性校正的
  # Wilson 区间，0.8342 与精确二项 0.8628 不一致，已改用 binom.test）
  n_pairs <- sum(df$group == "TBI") * sum(df$group == "Sham")
  concord <- round(a * n_pairs)
  ci_exact <- tryCatch(
    binom.test(concord, n_pairs, conf.level = 0.95)$conf.int,
    error = function(e) c(NA_real_, NA_real_))
  # 训练阈值（同一冻结系数在训练数据上的 best 阈值，与 04b 逻辑一致）
  tr_df2 <- tr_df[, c(genes[genes %in% colnames(tr_df)], "group"), drop = FALSE]
  tr_df2 <- tr_df2[complete.cases(tr_df2), ]
  lp_tr <- cf$coefficient[cf$term == "(Intercept)"]
  for (g in genes) lp_tr <- lp_tr + cf$coefficient[cf$term == g] * tr_df2[[g]]
  roc_tr <- roc(as.numeric(tr_df2$group == "TBI"), plogis(lp_tr), quiet = TRUE)
  co <- coords(roc_tr, "best", ret = c("threshold", "sensitivity",
                                       "specificity"))
  thr <- as.numeric(co$threshold[1])
  pred_class <- ifelse(prob >= thr, "TBI", "Sham")
  tab <- table(factor(pred_class, levels = c("Sham", "TBI")), df$group)
  sens <- if ("TBI" %in% colnames(tab)) tab["TBI", "TBI"] / sum(tab[, "TBI"]) else NA
  spec <- if ("Sham" %in% colnames(tab)) tab["Sham", "Sham"] / sum(tab[, "Sham"]) else NA
  data.frame(
    model = model_label, scale = scale_label, n = nrow(df),
    n_tbi = sum(df$group == "TBI"), n_sham = sum(df$group == "Sham"),
    auc = round(a, 4), ci_low = round(ci[1], 4), ci_high = round(ci[3], 4),
    concordant_pairs = concord, pairs_total = n_pairs,
    exact_ci_low = round(ci_exact[1], 4),
    exact_ci_high = round(ci_exact[2], 4),
    threshold = round(thr, 4), sensitivity = round(sens, 4),
    specificity = round(spec, 4),
    missing_genes = "", note = "", stringsAsFactors = FALSE)
}

df_fpkm <- df_val
df_log2 <- df_val
for (g in intersect(unique(c(core13, cand28, sig3)),
                    colnames(df_val))) {
  df_log2[[g]] <- log2(df_val[[g]] + 1)
}

res_raw <- list()
res_log2 <- list()
for (m in intersect(models_all, names(tr_data))) {
  tr_df <- tr_data[[m]]
  if (is.null(tr_df)) next
  res_raw[[length(res_raw) + 1]] <- eval_frozen(m, df_fpkm, tr_df,
                                                scale_label = "fpkm_raw")
  res_log2[[length(res_log2) + 1]] <- eval_frozen(m, df_log2, tr_df,
                                                  scale_label = "log2_fpkm")
}
res_gse205958 <- do.call(rbind, c(res_raw, res_log2))
res_gse205958$validation <- "GSE205958"
res_gse205958$species <- "mouse"
res_gse205958$model_type <- "CCI cortex 7d RNA-seq"

# ---------- 5. GSE128543 WT-only 方法学互查 ----------
gs128_wt <- build_df("GSE128543", unique(c(core13, cand28, sig3)),
                     clean = FALSE)
if (!is.null(gs128_wt)) {
  gt <- if ("genotype:ch1" %in% colnames(read_bulk("GSE128543")$pheno)) {
    as.character(read_bulk("GSE128543")$pheno$`genotype:ch1`)
  } else rep("", nrow(gs128_wt))
  is_wt <- grepl("^\\s*wt\\s*$|wildtype", gt, ignore.case = TRUE) |
    grepl("\\bwt\\b|wildtype",
          as.character(read_bulk("GSE128543")$pheno$title), ignore.case = TRUE)
  grp128 <- infer_group(as.character(read_bulk("GSE128543")$pheno$title))
  grp128[!is_wt] <- "Other"
  gs128_wt$group <- factor(ifelse(grp128 == "TBI", "TBI", "Sham"),
                           levels = c("Sham", "TBI"))
  gs128_wt <- gs128_wt[grp128 %in% c("TBI", "Sham"), ]
  res_ref <- eval_frozen("CLEAN_3gene_only", gs128_wt,
                         tr_data$CLEAN_3gene_only, scale_label = "as_stored")
  res_ref$validation <- "GSE128543_WT_only"
  res_ref$species <- "mouse"
  res_ref$model_type <- "CCI cortex 48h RNA-seq (reference check)"
  cat("GSE128543 WT-only check AUC:", res_ref$auc,
      " (expect 1.0 per frozen_signature_validation.csv)\n")
} else {
  res_ref <- NULL
}

res_all <- rbind(res_gse205958, res_ref)
write.csv(res_all, file.path(OUT_D12, "external_validation_v3.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("===== external_validation_v3 =====\n")
print(res_all)

# ---------- 6. 主 ROC 图（CLEAN_3gene_only, log2_fpkm） ----------
main_row <- res_gse205958[res_gse205958$model == "CLEAN_3gene_only" &
                            res_gse205958$scale == "log2_fpkm", ]
if (nrow(main_row) == 1 && !is.na(main_row$auc)) {
  cf3 <- coefs[coefs$model == "CLEAN_3gene_only", ]
  genes3 <- cf3$term[cf3$term != "(Intercept)"]
  lp3 <- cf3$coefficient[cf3$term == "(Intercept)"]
  for (g in genes3) lp3 <- lp3 + cf3$coefficient[cf3$term == g] * df_fpkm[[g]]
  prob3 <- plogis(lp3)
  roc3 <- roc(df_fpkm$group, prob3, levels = c("Sham", "TBI"),
              direction = "<", quiet = TRUE)
  roc_df <- data.frame(fpr = 1 - roc3$specificities,
                       tpr = roc3$sensitivities)
  roc_df <- roc_df[order(roc_df$fpr, roc_df$tpr), ]  # FPR 升序，保证单调
  roc_df <- roc_df[!duplicated(roc_df), ]
  auc_txt <- sprintf("AUC = %.3f (95%% CI %.3f–%.3f)",
                     main_row$auc, main_row$ci_low, main_row$ci_high)
  p <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey50") +
    geom_step(colour = "#B2182B", linewidth = 1.1, direction = "hv") +
    annotate("text", x = 0.68, y = 0.22, label = auc_txt, size = 4,
             hjust = 0) +
    labs(x = "1 – Specificity", y = "Sensitivity",
         title = "Frozen 3-gene signature in GSE205958 (CCI cortex, 7 d)",
         subtitle = "TBI 5 vs Sham 5; log2(FPKM+1) scale") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1))
  if (exists("theme_pub", inherits = FALSE)) p <- p + theme_pub()
  ggsave_pub(file.path(OUT_D12, "frozen_roc_GSE205958.pdf"), p,
         width = 5.5, height = 5.5, dpi = 300, bg = "white")
  cat("ROC saved.\n")
}

# ---------- 7. 3 基因表达摘要 ----------
present3 <- intersect(sig3, colnames(df_fpkm))
if (length(present3) == 3) {
  sum_df <- data.frame(gene = present3,
                       mean_fpkm_tbi = vapply(present3, function(g)
                         mean(df_fpkm[[g]][df_fpkm$group == "TBI"]), numeric(1)),
                       mean_fpkm_sham = vapply(present3, function(g)
                         mean(df_fpkm[[g]][df_fpkm$group == "Sham"]), numeric(1)),
                       median_fpkm_tbi = vapply(present3, function(g)
                         median(df_fpkm[[g]][df_fpkm$group == "TBI"]), numeric(1)),
                       median_fpkm_sham = vapply(present3, function(g)
                         median(df_fpkm[[g]][df_fpkm$group == "Sham"]), numeric(1)))
  cf3 <- coefs[coefs$model == "CLEAN_3gene_only", ]
  sign3 <- setNames(sign(cf3$coefficient[cf3$term != "(Intercept)"]),
                    cf3$term[cf3$term != "(Intercept)"])
  sum_df$observed_direction <- vapply(present3, function(g) {
    ifelse(sum_df$mean_fpkm_tbi[sum_df$gene == g] >
             sum_df$mean_fpkm_sham[sum_df$gene == g],
           "TBI>Sham", "TBI<Sham")
  }, character(1))
  sum_df$training_sign <- vapply(present3, function(g) {
    if (g %in% names(sign3)) {
      if (sign3[[g]] > 0) "TBI>Sham (training positive)" else
        "TBI<Sham (training negative)"
    } else "zero weight in frozen fit"
  }, character(1))
  sum_df$direction_consistent <- vapply(present3, function(g) {
    if (g %in% names(sign3)) {
      observed <- sum_df$mean_fpkm_tbi[sum_df$gene == g] >
        sum_df$mean_fpkm_sham[sum_df$gene == g]
      observed == (sign3[[g]] > 0)
    } else NA
  }, logical(1))
  write.csv(sum_df, file.path(OUT_D12, "signature_gene_expression_GSE205958.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  print(sum_df)
}

# ---------- 8. 报告 ----------
md <- c(
  "# D12 外部验证报告：GSE205958（小鼠 CCI 皮层，7 d）",
  "",
  paste("> 生成时间：", Sys.time()),
  "> 队列：Bao W et al. Brain Res 2023 (PMID 36335996)；GSE205958",
  "> 模型：CCI；组织：大脑皮层；取材：TBI 后 7 天；对照：naïve（未手术）小鼠",
  "> 平台：Illumina NovaSeq 6000 (GPL24247)；mm10；StringTie+ballgown FPKM",
  "",
  "## 1. 选择依据",
  "",
  "- 同物种（小鼠）、同模型（CCI）、同组织（皮层）、全转录组 RNA-seq；",
  "- 独立发表队列（Brain Res 2023），非本项目训练/验证来源；",
  "- 对照为 naïve 未手术小鼠（无开颅），与 GSE128543 Sham 定义一致；",
  "- 时点为 TBI 后 7 天（亚急性），与发现队列 GSE58485（3 d）和",
  "  GSE128543（48 h）不同，属跨时点敏感性验证。",
  "",
  "## 2. 方法",
  "",
"- 主口径：固定 CLEAN_3gene_only 冻结系数（P2RY6/HSD11B1/CCR5），",
  "  不重训、不插补；log2(FPKM+1) 与训练 log2 尺度一致；",
  "  原始线性 FPKM 作跨尺度敏感性（仅排序参考）；",
  "- 统计：AUC + Delong 95% CI；AUC=1.0 时 Delong CI 退化，改报",
  "  TBI–Sham 配对一致性比例的精确二项（Clopper–Pearson）95% CI",
  "  （AUC 即该比例；配对非独立，区间仅描述性）；",
  "  灵敏度/特异度采用训练数据 best 阈值；",
  "- 方法学互查：同一套冻结系数在 GSE128543 WT-only 复算，应与既有",
  "  frozen_signature_validation.csv 一致。",
  "",
  "## 3. 结果",
  "",
  paste0("- 样本：TBI 5 vs Sham 5（n=10）；"),
  paste0("- 主口径（CLEAN_3gene_only, log2_fpkm）：AUC=",
         if (nrow(main_row) == 1 && !is.na(main_row$auc)) {
           sprintf("%.3f (95%% CI %.3f–%.3f)", main_row$auc,
                   main_row$ci_low, main_row$ci_high)
         } else "not evaluable"),
  paste0("- 主口径配对一致性（精确二项（Clopper–Pearson）95% CI）：",
         if (nrow(main_row) == 1 && !is.na(main_row$exact_ci_low)) {
           sprintf("%d/%d = %.3f (%.3f–%.3f)",
                   main_row$concordant_pairs, main_row$pairs_total,
                   main_row$auc, main_row$exact_ci_low,
                   main_row$exact_ci_high)
         } else "not available"),
  "- 其余模型与 log2 尺度见 `external_validation_v3.csv`；",
  "- 3 基因表达均值/中位数见 `signature_gene_expression_GSE205958.csv`。",
  "",
  "## 4. 判读与局限",
  "",
  "- n=10（5 vs 5）样本量极小，即使 AUC=1.0 抽样不确定性仍大，",
  "  所有 CI 均须按此样本量解读；",
  "- AUC=1.0 为完全分离，Delong/bootstrap CI 均退化；精确二项区间",
  "  仅反映 n=10（25 个非独立 TBI–Sham 配对）下的抽样不确定性；",
  "- 训练阈值迁移到 FPKM 尺度后灵敏度/特异度可能失去判别意义（如",
  paste0("  CLEAN_3gene_only 阈值 ",
         if (nrow(main_row) == 1 && !is.na(main_row$threshold)) {
           main_row$threshold
         } else "NA",
         " 下 GSE205958 全部样本被判为 TBI："),
  "  sens=1/spec=0），故以 AUC 排序能力为主要判据，阈值仅报告不作结论；",
  "- 逐基因方向（vs v4 冻结系数符号）：GSE205958 中 P2RY6↑、HSD11B1↓",
  "  与训练方向一致，CCR5 观察方向为↑而冻结系数为负（方向相反），",
  "  与 GSE128543 的 per-gene 方向结果一致（见 26 脚本产物）；",
  "  该基因在论文中必须如实写为不稳定，不作为方向性复制证据；",
  "- 跨平台（微阵列 vs RNA-seq FPKM）与跨时点（3 d/48 h vs 7 d）使",
  "  系数尺度仅近似，方向性证据 > 数值证据；",
  "- 对照为 naïve 小鼠（无手术/无麻醉对照），与发现队列 Sham 定义",
  "  （vehicle/假手术）存在差异，可能放大组间效应；",
  "- 本验证为假设生成链上的复制性检查，不构成因果或临床预测证据。",
  "",
  "## 5. 数据来源",
  "",
  "- GEO：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205958",
  "- 补充文件：GSE205958_genes_fpkm_expression.xlsx（FPKM）",
  "- 论文：Bao W, Sun Y, Lin Y, Yang X, Chen Z. Brain Res. 2023;1799:148149."
)
writeLines(md, file.path(OUT_D12, "EXPERIMENT_REPORT.md"),
           useBytes = FALSE)
cat("report written.\n")
message("D12 完成。结果见 results/08_External_Validation/D12_GSE205958/")
