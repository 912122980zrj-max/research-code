# ============================================================
# 04_ml_shap.R —— 多算法机器学习 + 跨队列诊断 + 线性归因
# 运行：Rscript scripts/04_ml_shap.R
# 前置：03 已运行（CLEAN 主结果）
# 说明（2026-08-07 审查修复版，T2/P3/P12/P13/P16）：
#   - 主模型在 ComBat 校正矩阵（CLEAN n=33）上训练，
#     使用 cohort x group 分层的重复 5 折 CV。
#   - 新增 leave-one-cohort-out 跨队列 AUC（raw 口径，FULL 与 CLEAN）。
#   - "SHAP" 更名为 linear attribution（线性模型系数归因），
#     不再使用不准确的 SHAP 命名。
#   - 方法学声明：ComBat 批次参数使用全样本估计（含验证折），
#     分层重复 CV 仅为描述性指标；泛化性以 leave-one-cohort-out 与固定签名验证为准。
# 输出：results/07_ML_Signature/
# ============================================================

suppressMessages({
  library(caret)
  library(glmnet)
  library(pROC)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(e1071)
  library(randomForest)
  library(xgboost)
  library(kernlab)
  library(mboost)
  library(MASS)
})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "results", "07_ML_Signature")
OUT_CAND <- file.path(ROOT, "results", "05_Candidate_Intersection")
OUT_PREP <- file.path(ROOT, "results", "00_Data_Prep")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts", "theme_pub.R"))

# ---------- 读入 CLEAN 校正矩阵与候选基因 ----------
cand_file <- file.path(OUT_CAND, "candidate_genes.csv")
if (!file.exists(cand_file)) stop("缺少 candidate_genes.csv，请先运行 03")
cand <- toupper(read.csv(cand_file)$gene)

corrected_expr <- as.matrix(read.csv(file.path(OUT_PREP, "corrected_expr.csv"),
                                     row.names = 1, check.names = FALSE))
pheno <- read.csv(file.path(OUT_PREP, "corrected_pheno.csv"),
                  check.names = FALSE, fileEncoding = "UTF-8")

rn <- toupper(rownames(corrected_expr))
keep_r <- rn %in% cand
ex <- corrected_expr[keep_r, , drop = FALSE]
rn2 <- rn[keep_r]
if (anyDuplicated(rn2)) {
  ex <- rowsum(ex, group = rn2)
} else {
  rownames(ex) <- rn2
}
feat <- rownames(ex)
feat <- feat[order(feat)]
ex <- ex[feat, , drop = FALSE]

df_main <- as.data.frame(t(ex), check.names = FALSE)
df_main$group <- factor(ifelse(pheno$group == "tbi", "TBI", "Sham"),
                        levels = c("Sham", "TBI"))
df_main <- df_main[complete.cases(df_main), ]
cohort_vec <- as.character(pheno$cohort[match(rownames(df_main),
                                              colnames(ex))])
if (any(is.na(cohort_vec))) {
  cohort_vec <- as.character(pheno$cohort)
}

cat("features:", length(feat), "| samples:", nrow(df_main),
    "| TBI:", sum(df_main$group == "TBI"),
    "| Sham:", sum(df_main$group == "Sham"), "\n")

# ---------- 分层重复 CV（cohort x group 内分层） ----------
make_folds <- function(group, cohort, k = 5, repeats = 5, seed = 2026) {
  n <- length(group)
  strata <- interaction(cohort, group, drop = TRUE)
  train_idx <- list()
  test_idx <- list()
  r <- 0
  for (rep in seq_len(repeats)) {
    set.seed(seed + rep)
    fold <- integer(n)
    for (st in levels(strata)) {
      idx <- which(strata == st)
      fold[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
    }
    for (f in seq_len(k)) {
      r <- r + 1
      test_idx[[r]] <- which(fold == f)
      train_idx[[r]] <- which(fold != f)
    }
  }
  list(index = train_idx, indexOut = test_idx)
}

folds <- make_folds(df_main$group, cohort_vec, k = 5, repeats = 5, seed = 2026)
ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 5,
                     classProbs = TRUE, summaryFunction = twoClassSummary,
                     index = folds$index, indexOut = folds$indexOut,
                     savePredictions = "final", allowParallel = FALSE)

# 注：xgboost 2.x 与 caret 存在已知兼容问题（ROC 全为 NA），暂用 glmboost 代表 boosting 族
methods <- c("glmnet", "rf", "svmRadial", "glmboost", "lda")
model_list <- list()
for (m in methods) {
  cat("training:", m, "\n")
  set.seed(2026)
  model_list[[m]] <- tryCatch(
    train(group ~ ., data = df_main, method = m, trControl = ctrl,
          metric = "ROC", tuneLength = 3),
    error = function(e) {
      message("  [warn] ", m, " 失败: ", conditionMessage(e)); NULL
    }
  )
}

roc_summary <- sapply(model_list, function(x) {
  if (is.null(x)) return(NA_real_)
  max(x$results$ROC, na.rm = TRUE)
})
cat("CV AUC (stratified, ComBat-corrected):\n"); print(round(roc_summary, 3))
cat("注意：分层 CV 为描述性指标（ComBat 用全样本估计批次参数），",
    "泛化性以 leave-one-cohort-out 与固定签名验证为准。\n")
write.csv(data.frame(method = names(roc_summary), cv_auc = roc_summary),
          file.path(OUT, "model_auc.csv"), row.names = FALSE)

best_method <- names(which.max(roc_summary))
best_model <- model_list[[best_method]]
cat("best method:", best_method, "\n")

# 最终模型：弹性网络（与模板一致）
if (!is.null(model_list[["glmnet"]])) {
  final_model <- model_list[["glmnet"]]
  cat("final model: glmnet (elastic net), bestTune alpha =",
      final_model$bestTune$alpha, "\n")
  cf <- coef(final_model$finalModel, s = final_model$bestTune$lambda)
  coef_df <- data.frame(gene = cf@Dimnames[[1]][cf@i + 1],
                        coefficient = as.numeric(cf@x))
  coef_df <- coef_df[coef_df$gene != "(Intercept)" & coef_df$coefficient != 0, ]
  coef_df <- coef_df[order(-abs(coef_df$coefficient)), ]
  write.csv(coef_df, file.path(OUT, "enet_coefficients_legacy_04.csv"),
            row.names = FALSE)
} else {
  final_model <- best_model
  coef_df <- data.frame(gene = feat, coefficient = NA)
}

# ---------- 稳定性分析（5 个随机种子，glmnet 特征保留频率） ----------
retention <- lapply(1:5, function(seed) {
  fd <- make_folds(df_main$group, cohort_vec, k = 5, repeats = 1, seed = seed)
  tr <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     index = fd$index, indexOut = fd$indexOut)
  set.seed(seed)
  m <- tryCatch(train(group ~ ., data = df_main, method = "glmnet",
                      trControl = tr, metric = "ROC", tuneLength = 3),
                error = function(e) NULL)
  if (is.null(m)) return(character(0))
  cf <- coef(m$finalModel, s = m$bestTune$lambda)
  g <- cf@Dimnames[[1]][cf@i + 1]
  g[g != "(Intercept)"]
})
ret_tab <- sort(table(unlist(retention)), decreasing = TRUE)
write.csv(data.frame(gene = names(ret_tab), retain_seeds = as.integer(ret_tab)),
          file.path(OUT, "stability_glmnet.csv"), row.names = FALSE)

# ---------- Linear attribution（线性模型系数归因，非 SHAP） ----------
Xs <- scale(df_main[, feat, drop = FALSE])
if (exists("coef_df") && nrow(coef_df) > 0) {
  beta <- rep(0, length(feat)); names(beta) <- feat
  hit <- intersect(coef_df$gene, feat)
  beta[hit] <- coef_df$coefficient[match(hit, coef_df$gene)]
  attr_mat <- sweep(Xs, 2, beta, "*")
  mean_abs <- sort(colMeans(abs(attr_mat)), decreasing = TRUE)
  write.csv(data.frame(gene = names(mean_abs), mean_abs_attribution =
                         as.numeric(mean_abs)),
            file.path(OUT, "linear_attribution_legacy_04.csv"), row.names = FALSE)

  imp_df <- data.frame(gene = factor(names(mean_abs),
                                     levels = rev(names(mean_abs))),
                       value = as.numeric(mean_abs))
  p1 <- ggplot(imp_df, aes(x = gene, y = value)) +
    geom_col(fill = "#4A6FA5", width = 0.7) + coord_flip() +
    labs(title = "Linear attribution importance (mean |beta * scaled X|)",
         y = "mean |attribution|") +
    theme_pub(base_size = 11)
  ggsave_pub(file.path(OUT, "linear_attribution_legacy_04.pdf"), p1,
         width = 6, height = 5, dpi = 300)

  long <- data.frame(
    sample = rep(seq_len(nrow(attr_mat)), ncol(attr_mat)),
    gene = rep(colnames(attr_mat), each = nrow(attr_mat)),
    attribution = as.vector(attr_mat),
    expr = as.vector(Xs)
  )
  long$gene <- factor(long$gene, levels = names(mean_abs))
  p2 <- ggplot(long, aes(x = attribution, y = gene, color = expr)) +
    geom_jitter(size = 0.6, alpha = 0.7, width = 0, height = 0.15) +
    scale_color_gradient2(low = pal_rdbu[1], mid = pal_rdbu[4],
                          high = pal_rdbu[6],
                          name = "scaled expr") +
    labs(title = "Linear attribution (ENet approximation)",
         x = "attribution (beta * scaled X)") +
    theme_pub(base_size = 11)
  ggsave_pub(file.path(OUT, "linear_beeswarm_legacy_04.pdf"), p2,
         width = 7, height = 6, dpi = 300)
}

# ---------- 跨队列 leave-one-cohort-out（raw 口径，FULL + CLEAN） ----------
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

core14_file <- file.path(OUT_CAND, "candidate_genes_FULL.csv")
core14 <- if (file.exists(core14_file)) {
  toupper(read.csv(core14_file)$gene)
} else {
  c("CYP1B1", "HSD11B1", "MMP2", "P2RY6", "PGF", "RIPK1")
}

build_df_raw <- function(gse, genes, clean = FALSE) {
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

loco_rows <- list()
for (comp in c("FULL", "CLEAN")) {
  for (dir_ in list(c("GSE58485", "GSE41345"), c("GSE41345", "GSE58485"))) {
    tr_gse <- dir_[1]; te_gse <- dir_[2]
    tr_df <- build_df_raw(tr_gse, core14, clean = (comp == "CLEAN"))
    te_df <- build_df_raw(te_gse, core14, clean = (comp == "CLEAN"))
    if (is.null(tr_df) || is.null(te_df)) next
    common <- intersect(attr(tr_df, "features"), attr(te_df, "features"))
    if (length(common) < 3) next
    tr_df <- tr_df[, c(common, "group")]
    te_df <- te_df[, c(common, "group")]
    tr_df <- tr_df[complete.cases(tr_df), ]
    te_df <- te_df[complete.cases(te_df), ]
    if (length(unique(te_df$group)) < 2 || nrow(te_df) < 3) next
    set.seed(2026)
    m <- tryCatch(
      train(group ~ ., data = tr_df, method = "glmnet",
            trControl = trainControl(method = "repeatedcv", number = 5,
                                     repeats = 5, classProbs = TRUE,
                                     summaryFunction = twoClassSummary),
            metric = "ROC", tuneLength = 3),
      error = function(e) NULL)
    if (is.null(m)) next
    prob <- tryCatch(predict(m, newdata = te_df, type = "prob")[, "TBI"],
                     error = function(e) NULL)
    if (is.null(prob)) next
    test_auc <- tryCatch(auc(roc(as.numeric(te_df$group == "TBI"), prob)),
                         error = function(e) NA_real_)
    loco_rows[[length(loco_rows) + 1]] <- data.frame(
      composition = comp,
      train_cohort = tr_gse,
      test_cohort = te_gse,
      n_train = nrow(tr_df),
      n_test = nrow(te_df),
      train_cv_auc = max(m$results$ROC, na.rm = TRUE),
      test_auc = as.numeric(test_auc),
      n_features = length(common),
      stringsAsFactors = FALSE)
  }
}
if (length(loco_rows) > 0) {
  loco <- do.call(rbind, loco_rows)
  write.csv(loco, file.path(OUT, "cross_cohort_auc.csv"), row.names = FALSE)
  cat("===== leave-one-cohort-out AUC =====\n")
  print(loco)
} else {
  warning("未生成跨队列 AUC（数据缺失）")
}

message("完成。结果在 results/07_ML_Signature/（分层 CV AUC、系数、线性归因、跨队列 AUC、稳定性）。")
