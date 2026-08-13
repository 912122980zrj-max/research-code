# ============================================================
# 04b_external_validation.R —— 外部验证（固定签名优先 + 95% CI）
# 运行：Rscript scripts/04b_external_validation.R
# 前置：03 / 04 已运行；01b 已补下载 GSE128543 / GSE104687
# 说明（2026-08-07 审查修复版，T2/P3/P4/P9/P16）：
#   - 主验证 = 固定 3 特征签名（P2RY6/HSD11B1/CCR5，log2 口径，冻结系数
#     不重训；v4 为 D4 SHAP top3）
#     在 GSE128543 WT-only（n=14）与 GSE205958（n=10）上评估，报 AUC + CI +
#     灵敏度/特异度；28 样本（含 RBM5 KO）作敏感性分析。
#   - GSE104687 明确标注为 transportability 阴性证据（P2RY6/CCR5/ESR1
#     缺失，非人类验证）。
# 输出：results/07_ML_Signature/frozen_signature_validation.csv
#       results/07_ML_Signature/external_validation_v2.csv
#       results/07_ML_Signature/frozen_coefficients.csv
# ============================================================

suppressMessages({
  library(caret)
  library(glmnet)
  library(pROC)
  library(dplyr)
  library(readr)
})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "results", "08_External_Validation")
OUT_CAND <- file.path(ROOT, "results", "05_Candidate_Intersection")
OUT_PREP <- file.path(ROOT, "results", "00_Data_Prep")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------- 数据准备 ----------
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

build_df_raw <- function(gse, genes, clean = FALSE, wt_only = FALSE) {
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
  if (wt_only) {
    gt <- if ("genotype:ch1" %in% colnames(d$pheno)) {
      as.character(d$pheno$`genotype:ch1`)
    } else {
      rep("", nrow(d$pheno))
    }
    is_wt <- grepl("^\\s*wt\\s*$|wildtype", gt, ignore.case = TRUE) |
      grepl("\\bwt\\b|wildtype", d$pheno$title, ignore.case = TRUE)
    grp[!is_wt] <- "Other"
  }
  df <- as.data.frame(t(ex2), check.names = FALSE)
  df$group <- factor(ifelse(grp == "TBI", "TBI", "Sham"),
                     levels = c("Sham", "TBI"))
  df <- df[grp %in% c("TBI", "Sham"), ]
  attr(df, "features") <- intersect(rownames(ex2), genes)
  df
}

fit_glmnet <- function(tr_df) {
  set.seed(2026)
  train(group ~ ., data = tr_df, method = "glmnet",
        trControl = trainControl(method = "repeatedcv", number = 5,
                                 repeats = 5, classProbs = TRUE,
                                 summaryFunction = twoClassSummary),
        metric = "ROC", tuneLength = 3)
}

best_threshold <- function(model, tr_df) {
  prob <- predict(model, newdata = tr_df, type = "prob")[, "TBI"]
  roc_tr <- roc(as.numeric(tr_df$group == "TBI"), prob, quiet = TRUE)
  co <- coords(roc_tr, "best", ret = c("threshold", "sensitivity",
                                       "specificity"))
  c(threshold = as.numeric(co$threshold[1]),
    sensitivity = as.numeric(co$sensitivity[1]),
    specificity = as.numeric(co$specificity[1]))
}

eval_frozen <- function(model, test_df, label_training, validation, subset,
                        note = "") {
  prob <- predict(model, newdata = test_df, type = "prob")[, "TBI"]
  roc_te <- roc(as.numeric(test_df$group == "TBI"), prob, quiet = TRUE)
  a <- as.numeric(auc(roc_te))
  ci <- tryCatch(as.numeric(ci.auc(roc_te, method = "delong")),
                 error = function(e) rep(NA_real_, 3))
  thr <- best_threshold(model, attr(model, "tr_df"))
  pred_class <- ifelse(prob >= thr["threshold"], "TBI", "Sham")
  tab <- table(factor(pred_class, levels = c("Sham", "TBI")),
               test_df$group)
  sens <- if ("TBI" %in% colnames(tab)) tab["TBI", "TBI"] / sum(tab[, "TBI"]) else NA
  spec <- if ("Sham" %in% colnames(tab)) tab["Sham", "Sham"] / sum(tab[, "Sham"]) else NA
  data.frame(
    signature_training = label_training,
    validation = validation,
    subset = subset,
    n = nrow(test_df),
    n_tbi = sum(test_df$group == "TBI"),
    n_sham = sum(test_df$group == "Sham"),
    auc = round(a, 4),
    ci_low = round(ci[1], 4),
    ci_high = round(ci[3], 4),
    threshold = round(thr["threshold"], 4),
    sensitivity = round(sens, 4),
    specificity = round(spec, 4),
    note = note,
    stringsAsFactors = FALSE)
}

extract_coefs <- function(model, label) {
  cf <- coef(model$finalModel, s = model$bestTune$lambda)
  g <- cf@Dimnames[[1]][cf@i + 1]
  v <- as.numeric(cf@x)
  df <- data.frame(model = label, term = g, coefficient = v,
                   stringsAsFactors = FALSE)
  df
}

# ---------- 候选基因 ----------
cand28_file <- file.path(OUT_CAND, "candidate_genes.csv")
cand28 <- if (file.exists(cand28_file)) toupper(read.csv(cand28_file)$gene) else NULL
core14_file <- file.path(OUT_CAND, "candidate_genes_FULL.csv")
core14 <- if (file.exists(core14_file)) {
  toupper(read.csv(core14_file)$gene)
} else {
  c("CYP1B1", "HSD11B1", "MMP2", "P2RY6", "PGF", "RIPK1")
}
sig3 <- c("P2RY6", "HSD11B1", "CCR5")  # v4 重跑后 D4 SHAP top3（7 候选 ENet）

# ---------- 训练（raw 口径，FULL 与 CLEAN） ----------
disc_full14 <- rbind(build_df_raw("GSE58485", core14, clean = FALSE),
                     build_df_raw("GSE41345", core14, clean = FALSE))
disc_clean14 <- rbind(build_df_raw("GSE58485", core14, clean = TRUE),
                      build_df_raw("GSE41345", core14, clean = TRUE))
disc_clean28 <- if (!is.null(cand28)) {
  rbind(build_df_raw("GSE58485", cand28, clean = TRUE),
        build_df_raw("GSE41345", cand28, clean = TRUE))
} else NULL

fit14_full <- fit_glmnet(disc_full14[, c(intersect(core14, colnames(disc_full14)), "group")])
fit14_clean <- fit_glmnet(disc_clean14[, c(intersect(core14, colnames(disc_clean14)), "group")])
disc_sig3 <- rbind(build_df_raw("GSE58485", sig3, clean = TRUE),
                   build_df_raw("GSE41345", sig3, clean = TRUE))
fit3_clean <- fit_glmnet(disc_sig3[, c(intersect(sig3, colnames(disc_sig3)),
                                       "group")])
disc_sig2 <- rbind(build_df_raw("GSE58485", sig3[2:3], clean = TRUE),
                   build_df_raw("GSE41345", sig3[2:3], clean = TRUE))
fit2_clean <- NULL
if (!is.null(cand28)) {
  fit28_clean <- fit_glmnet(disc_clean28[, c(intersect(cand28, colnames(disc_clean28)), "group")])
  fit2_clean <- fit_glmnet(disc_sig2[, c(intersect(sig3[2:3],
                                                    colnames(disc_sig2)), "group")])
} else {
  fit28_clean <- NULL
  fit2_clean <- fit_glmnet(disc_sig2[, c(intersect(sig3[2:3],
                                                    colnames(disc_sig2)), "group")])
}

# 记录训练集，供阈值计算
attr(fit14_full, "tr_df") <- disc_full14[, c(intersect(core14, colnames(disc_full14)), "group")]
attr(fit14_clean, "tr_df") <- disc_clean14[, c(intersect(core14, colnames(disc_clean14)), "group")]
attr(fit3_clean, "tr_df") <- disc_sig3[, c(intersect(sig3, colnames(disc_sig3)),
                                           "group")]
attr(fit2_clean, "tr_df") <- disc_sig2[, c(intersect(sig3[2:3],
                                                      colnames(disc_sig2)), "group")]
if (!is.null(fit28_clean)) {
  attr(fit28_clean, "tr_df") <- disc_clean28[, c(intersect(cand28, colnames(disc_clean28)), "group")]
}

coefs <- rbind(
  extract_coefs(fit14_full, "FULL_5gene"),
  extract_coefs(fit14_clean, "CLEAN_5gene"),
  extract_coefs(fit3_clean, "CLEAN_3gene_only"),
  extract_coefs(fit2_clean, "CLEAN_2gene_HSD11B1_CCR5")
)
if (!is.null(fit28_clean)) {
  coefs <- rbind(coefs, extract_coefs(fit28_clean, "CLEAN_7gene"))
}
write.csv(coefs, file.path(OUT, "frozen_coefficients.csv"), row.names = FALSE)

# ---------- 验证数据 ----------
test_genes <- unique(c(core14, sig3,
                       if (!is.null(cand28)) cand28 else character(0)))
gs128_all <- build_df_raw("GSE128543", test_genes, clean = FALSE)
gs128_wt <- build_df_raw("GSE128543", test_genes, clean = FALSE,
                         wt_only = TRUE)
gs104 <- build_df_raw("GSE104687", test_genes, clean = FALSE)

res <- list()
if (!is.null(gs128_all)) {
  res[[length(res) + 1]] <- eval_frozen(
    fit14_full, gs128_all, "FULL_5gene_frozen", "GSE128543", "all (n=28, incl RBM5 KO)",
    note = "与审查§2.3固定签名诊断互查")
  res[[length(res) + 1]] <- eval_frozen(
    fit14_clean, gs128_all, "CLEAN_5gene_frozen", "GSE128543", "all (n=28, incl RBM5 KO)",
    note = "敏感性分析")
  res[[length(res) + 1]] <- eval_frozen(
    fit3_clean, gs128_all, "CLEAN_3gene_only_frozen", "GSE128543", "all (n=28, incl RBM5 KO)",
    note = "敏感性分析：严格3基因冻结模型")
}
if (!is.null(gs128_wt)) {
  res[[length(res) + 1]] <- eval_frozen(
    fit3_clean, gs128_wt, "CLEAN_3gene_only_frozen", "GSE128543", "WT-only (n=14)",
    note = "主验证：严格3基因冻结系数，同型CCI皮层，WT 7 vs 7；Sham无开颅")
  res[[length(res) + 1]] <- eval_frozen(
    fit14_full, gs128_wt, "FULL_5gene_frozen", "GSE128543", "WT-only (n=14)",
    note = "敏感性分析")
  res[[length(res) + 1]] <- eval_frozen(
    fit14_clean, gs128_wt, "CLEAN_5gene_frozen", "GSE128543", "WT-only (n=14)",
    note = "敏感性分析：5基因冻结模型")
  if (!is.null(fit28_clean)) {
    gs128_wt28 <- build_df_raw("GSE128543", cand28, clean = FALSE, wt_only = TRUE)
    if (!is.null(gs128_wt28) && length(intersect(cand28, colnames(gs128_wt28))) >= 3) {
      res[[length(res) + 1]] <- eval_frozen(
        fit28_clean, gs128_wt28, "CLEAN_7gene_frozen", "GSE128543", "WT-only (n=14)",
        note = "7候选基因口径敏感性分析")
    }
  }
}
if (!is.null(gs104) && !is.null(fit2_clean)) {
  gs104_2 <- gs104[, c(intersect(sig3[2:3], colnames(gs104)), "group")]
  gs104_2 <- gs104_2[complete.cases(gs104_2), ]
  if (length(intersect(sig3[2:3], colnames(gs104))) >= 2 &&
      nrow(gs104_2) > 3 && length(unique(gs104_2$group)) == 2) {
    res[[length(res) + 1]] <- eval_frozen(
      fit2_clean, gs104_2, "CLEAN_2gene_frozen", "GSE104687", "all (n=376)",
      note = "transportability阴性：跨物种/跨时相/跨平台；HSD11B1/CCR5 口径，非人类验证")
  }
  miss3 <- setdiff(sig3, colnames(gs104))
  res[[length(res) + 1]] <- data.frame(
    signature_training = "CLEAN_3gene_frozen",
    validation = "GSE104687",
    subset = "not evaluable",
    n = NA, n_tbi = NA, n_sham = NA, auc = NA, ci_low = NA, ci_high = NA,
    threshold = NA, sensitivity = NA, specificity = NA,
    note = paste0("固定 3 基因签名在 GSE104687 不可直接评估（缺失：",
                  paste(miss3, collapse = "/"), "）"),
    stringsAsFactors = FALSE)
}

frozen <- do.call(rbind, res)
write.csv(frozen, file.path(OUT, "frozen_signature_validation.csv"),
          row.names = FALSE)
cat("===== 固定签名验证 =====\n")
print(frozen)

# ---------- 外部验证 v2（共同特征重训，标注口径；保留旧功能但明确标签） ----------
validate_cohort_v2 <- function(gse, label, use_sig = FALSE, wt_only = FALSE) {
  genes <- if (use_sig) sig3 else core14
  df <- build_df_raw(gse, genes, clean = FALSE, wt_only = wt_only)
  if (is.null(df)) return(NULL)
  common_feat <- intersect(genes, colnames(df))
  min_feat <- if (use_sig) 2 else 3
  if (length(common_feat) < min_feat) return(NULL)
  df <- df[, c(common_feat, "group")]
  df <- df[complete.cases(df), ]
  if (nrow(df) < 3 || length(unique(df$group)) < 2) return(NULL)
  tr_pool <- if (use_sig) disc_sig3 else disc_full14
  tr_df <- tr_pool[, c(common_feat, "group")]
  tr_df <- tr_df[complete.cases(tr_df), ]
  model <- tryCatch(fit_glmnet(tr_df), error = function(e) NULL)
  if (is.null(model)) return(NULL)
  prob <- tryCatch(predict(model, newdata = df, type = "prob")[, "TBI"],
                   error = function(e) NULL)
  if (is.null(prob)) return(NULL)
  roc_te <- roc(as.numeric(df$group == "TBI"), prob, quiet = TRUE)
  a <- as.numeric(auc(roc_te))
  ci <- tryCatch(as.numeric(ci.auc(roc_te, method = "delong")),
                 error = function(e) rep(NA_real_, 3))
  data.frame(
    validation = gse,
    model = label,
    subset = if (wt_only) "WT-only" else "all",
    auc = round(a, 4),
    ci_low = round(ci[1], 4),
    ci_high = round(ci[3], 4),
    n = nrow(df),
    n_tbi = sum(df$group == "TBI"),
    n_sham = sum(df$group == "Sham"),
    n_features = length(common_feat),
    features = paste(common_feat, collapse = ";"),
    note = if (gse == "GSE104687") {
      "transportability阴性：非人类验证"
    } else {
      "共同特征重训（跨平台标准做法），非固定签名"
    },
    stringsAsFactors = FALSE)
}

val_list <- list()
for (gse in c("GSE128543", "GSE104687")) {
  if (is.null(read_bulk(gse))) {
    message("[skip] ", gse, " 数据不在 data/01_GEO_bulk/，请先运行 01b")
    next
  }
  if (gse == "GSE128543") {
    val_list[[length(val_list) + 1]] <- validate_cohort_v2(
      gse, "ENet common features (FULL-trained)", use_sig = FALSE)
    val_list[[length(val_list) + 1]] <- validate_cohort_v2(
      gse, "ENet common features (FULL-trained)", use_sig = FALSE, wt_only = TRUE)
    val_list[[length(val_list) + 1]] <- validate_cohort_v2(
    gse, "3-gene retrained (CLEAN discovery)", use_sig = TRUE)
    val_list[[length(val_list) + 1]] <- validate_cohort_v2(
    gse, "3-gene retrained (CLEAN discovery)", use_sig = TRUE, wt_only = TRUE)
  } else {
    val_list[[length(val_list) + 1]] <- validate_cohort_v2(
      gse, "ENet common features (FULL-trained)", use_sig = FALSE)
    val_list[[length(val_list) + 1]] <- validate_cohort_v2(
    gse, "3-gene retrained (CLEAN discovery)", use_sig = TRUE)
  }
}
val2 <- do.call(rbind, val_list)
if (!is.null(val2)) {
  write.csv(val2, file.path(OUT, "external_validation_v2.csv"), row.names = FALSE)
  cat("===== 外部验证 v2（共同特征重训） =====\n")
  print(val2)
} else {
  cat("没有可用的验证队列。\n")
}

message("完成。固定签名验证见 frozen_signature_validation.csv；")
message("重训口径见 external_validation_v2.csv；系数见 frozen_coefficients.csv。")
