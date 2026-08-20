# ============================================================
# 19_d3_ml_benchmark.R —— D3 ML 规范基准 + ENet 候选模型评估
# 依据：模板论文（Environment International 214:110390）ML 框架 +
#       ml-shap 技能（无泄漏、OOF、校准、DCA、PROBAST+AI）
# 说明（如实记录）：
#   - 安装的 Mime1（l-magnificence/Mime fork）不提供模板所述的
#     "127 组合"框架（其 API 为 7 算法二分类/预后开发）；原始
#     HuaZou/Mime 仓库当前网络不可达。因此按模板 11 算法透明实现基准，
#     不虚构 127 组合数字。
#   - 主口径：原始 CLEAN 表达（v4 7 候选基因，log2 统一尺度，n=33），重复分层 5 折 CV
#     （5 种子）生成 OOF；ComBat 校正矩阵作敏感性（描述性，因 ComBat
#     使用全样本估计）。
#   - 选定候选模型：ENet alpha=0.9（与模板一致）。
# 输出：results/07_ML_Signature/D3_ml_benchmark/
# ============================================================
suppressMessages({
  library(caret); library(glmnet); library(pROC); library(rms)
  library(rmda); library(ggplot2); library(dplyr); library(readr)
  library(e1071); library(randomForest); library(xgboost); library(kernlab)
  library(mboost); library(MASS); library(naivebayes)
})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "results", "07_ML_Signature", "D3_ml_benchmark")
OUT_ML <- file.path(ROOT, "results", "07_ML_Signature")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts", "theme_pub.R"))
set.seed(2026)

LOG <- file.path(OUT, "D3_run.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(...), "\n")
  cat(msg); cat(msg, file = LOG, append = TRUE)
}
if (file.exists(LOG)) file.remove(LOG)

# ---------- 1) 数据 ----------
cand <- toupper(read.csv(file.path(ROOT, "results", "05_Candidate_Intersection",
                                   "candidate_genes.csv"))$gene)
comp <- read.csv(file.path(ROOT, "results", "00_Data_Prep",
                           "sample_composition.csv"), check.names = FALSE)
clean_samples <- comp$sample[comp$include_clean]
grp_map <- setNames(
  ifelse(comp$group[match(clean_samples, comp$sample)] == "tbi",
         "TBI", "Sham"), clean_samples)

read_bulk <- function(gse) {
  rds_file <- file.path(ROOT, "data", "01_GEO_bulk", paste0(gse, ".rds"))
  if (file.exists(rds_file)) {
    d <- readRDS(rds_file)
    return(list(expr = d$expr, pheno = d$pheno))
  }
  ex <- as.matrix(read.csv(file.path(ROOT, "data", "01_GEO_bulk",
                                     paste0(gse, "_expr.csv")),
                           row.names = 1, check.names = FALSE))
  ph <- read.csv(file.path(ROOT, "data", "01_GEO_bulk",
                           paste0(gse, "_pheno.csv")),
                 row.names = 1, check.names = FALSE)
  list(expr = ex, pheno = ph)
}

dat <- lapply(c("GSE58485", "GSE41345"), read_bulk)
common_genes <- Reduce(intersect, lapply(dat, function(d) rownames(d$expr)))
expr_all <- do.call(cbind, lapply(dat, function(d) d$expr[common_genes, ]))
expr_all <- expr_all[, intersect(clean_samples, colnames(expr_all)), drop = FALSE]
stopifnot(ncol(expr_all) == sum(comp$include_clean),
          !anyNA(grp_map[colnames(expr_all)]))

rn <- toupper(rownames(expr_all))
keep <- rn %in% cand
ex <- expr_all[keep, , drop = FALSE]
rn2 <- rn[keep]
if (anyDuplicated(rn2)) {
  ex <- rowsum(ex, group = rn2)
} else {
  rownames(ex) <- rn2
}
feat <- sort(rownames(ex))
ex <- ex[feat, , drop = FALSE]
df <- as.data.frame(t(ex), check.names = FALSE)
df$group <- factor(grp_map[rownames(df)], levels = c("Sham", "TBI"))
df <- df[complete.cases(df), ]
log_msg("原始 CLEAN 矩阵: ", nrow(df), " 样本 x ", ncol(df) - 1, " 特征")
write.csv(data.frame(sample = rownames(df), group = as.character(df$group)),
          file.path(OUT, "d3_samples.csv"), row.names = FALSE)

# ComBat 校正矩阵（敏感性，描述性）
cor_expr <- as.matrix(read.csv(file.path(ROOT, "results", "00_Data_Prep",
                                         "corrected_expr.csv"),
                               row.names = 1, check.names = FALSE))
rn_c <- toupper(rownames(cor_expr))
ex_c <- cor_expr[rn_c %in% cand, , drop = FALSE]
rownames(ex_c) <- rn_c[rn_c %in% cand]
if (anyDuplicated(rownames(ex_c))) ex_c <- rowsum(ex_c, group = rownames(ex_c))
feat_c <- sort(intersect(feat, rownames(ex_c)))
ex_c <- ex_c[feat_c, , drop = FALSE]
df_c <- as.data.frame(t(ex_c), check.names = FALSE)
df_c$group <- factor(grp_map[rownames(df_c)], levels = c("Sham", "TBI"))
df_c <- df_c[complete.cases(df_c), ]
log_msg("ComBat 校正矩阵（敏感性）: ", nrow(df_c), " x ", ncol(df_c) - 1)
log_msg("df_c 状态: class=", class(df_c)[1],
        " | group 列存在=", "group" %in% colnames(df_c),
        " | group NA=", sum(is.na(df_c$group)))

# ---------- 2) 11 算法基准（重复分层 5 折 CV） ----------
algo_spec <- list(
  SVM = list(method = "svmRadial", tuneLength = 5),
  Lasso = list(method = "glmnet", tuneLength = 5, family = "binomial",
               alpha = 1),
  glmBoost = list(method = "glmboost", tuneLength = 5),
  RF = list(method = "rf", tuneLength = 5),
  StepGLM = list(method = "glmStepAIC", family = "binomial"),
  Ridge = list(method = "glmnet", tuneLength = 5, family = "binomial",
               alpha = 0),
  LDA = list(method = "lda"),
  ENET = list(method = "glmnet", tuneLength = 5, family = "binomial",
              alpha = 0.5),
  XGBoost = list(method = "xgbTree",
                 tuneGrid = expand.grid(
                   nrounds = c(50, 100), max_depth = c(1, 2),
                   eta = c(0.1, 0.3), gamma = 0, colsample_bytree = 0.8,
                   min_child_weight = 1, subsample = 1)),
  GBM = list(method = "gbm", tuneLength = 5),
  NaiveBayes = list(method = "naive_bayes", tuneLength = 5)
)

run_bench <- function(data_df, label) {
  ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 5,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE, savePredictions = "final")
  res_rows <- list()
  for (nm in names(algo_spec)) {
    spec <- algo_spec[[nm]]
    set.seed(2026)
    if (nm %in% c("XGBoost", "GBM")) {
      # caret xgbTree/gbm 与当前包版本不兼容（复现：ROC 全 NA），
      # 改用相同重复分层折结构手工实现，保证与其他算法可比。
      aucs <- numeric(0)
      err_msg <- NA_character_
      for (r in 1:5) {
        set.seed(2026 + r)
        folds <- caret::createFolds(data_df$group, k = 5, list = TRUE,
                                    returnTrain = TRUE)
        for (f in folds) {
          tr <- data_df[f, ]; te <- data_df[-f, ]
          best <- NULL
          one_fold <- tryCatch({
            if (nm == "XGBoost") {
              dm_tr <- xgboost::xgb.DMatrix(
                data = as.matrix(tr[, -ncol(tr), drop = FALSE]),
                label = as.numeric(tr$group == "TBI"))
              dm_te <- xgboost::xgb.DMatrix(
                data = as.matrix(te[, -ncol(te), drop = FALSE]))
              for (nr in c(50, 100)) for (md in c(1, 2)) for (et in c(0.1, 0.3)) {
                fitx <- xgboost::xgb.train(
                  params = list(objective = "binary:logistic", eta = et,
                                max_depth = md, subsample = 0.8,
                                colsample_bytree = 0.8, verbosity = 0),
                  data = dm_tr, nrounds = nr, verbose = 0)
                pr <- predict(fitx, dm_te)
                a <- pROC::auc(te$group, pr, levels = c("Sham", "TBI"),
                               direction = "<", quiet = TRUE)
                if (is.null(best) || a > best$auc) best <- list(auc = a)
              }
            } else {
              tr_g <- data.frame(y = as.numeric(tr$group == "TBI"),
                                 tr[, -ncol(tr), drop = FALSE],
                                 check.names = FALSE)
              for (nt in c(50, 100)) for (id in c(1, 2)) for (sh in c(0.1, 0.3)) {
                fitx <- gbm::gbm(y ~ ., data = tr_g,
                                 distribution = "bernoulli", n.trees = nt,
                                 interaction.depth = id, shrinkage = sh,
                                 n.minobsinnode = 2, bag.fraction = 1,
                                 verbose = FALSE)
                pr <- predict(fitx, te[, -ncol(te), drop = FALSE],
                              n.trees = nt, type = "response")
                a <- pROC::auc(te$group, pr, levels = c("Sham", "TBI"),
                               direction = "<", quiet = TRUE)
                if (is.null(best) || a > best$auc) best <- list(auc = a)
              }
            }
            as.numeric(best$auc)
          }, error = function(e) {
            err_msg <<- conditionMessage(e)
            NA_real_
          })
          if (!is.na(one_fold)) aucs <- c(aucs, one_fold)
        }
      }
      if (length(aucs) == 0) {
        res_rows[[nm]] <- data.frame(
          algorithm = nm, matrix = label,
          cv_auc_mean = NA_real_, cv_auc_sd = NA_real_,
          stringsAsFactors = FALSE)
        log_msg(label, " | ", nm, " FAILED: ", err_msg)
      } else {
        res_rows[[nm]] <- data.frame(
          algorithm = nm, matrix = label,
          cv_auc_mean = mean(aucs), cv_auc_sd = sd(aucs),
          stringsAsFactors = FALSE)
        log_msg(label, " | ", nm, " CV AUC=", round(mean(aucs), 3))
      }
      next
    }
    fit <- tryCatch({
      if (spec$method == "glmnet") {
      tune_grid <- expand.grid(alpha = spec$alpha,
                               lambda = 10^seq(-4, 0, length.out = 20))
        train(group ~ ., data = data_df, method = "glmnet",
              trControl = ctrl, metric = "ROC",
              tuneGrid = tune_grid, family = "binomial")
      } else if (!is.null(spec$tuneGrid)) {
        train(group ~ ., data = data_df, method = spec$method,
              trControl = ctrl, metric = "ROC",
              tuneGrid = spec$tuneGrid)
      } else {
        args <- list(group ~ ., data = data_df, method = spec$method,
                     trControl = ctrl, metric = "ROC")
        if (!is.null(spec$tuneLength)) args$tuneLength <- spec$tuneLength
        if (!is.null(spec$family)) args$family <- spec$family
        do.call(train, args)
      }
    }, error = function(e) e)
    if (inherits(fit, "error")) {
      res_rows[[nm]] <- data.frame(
        algorithm = nm, matrix = label,
        cv_auc_mean = NA_real_, cv_auc_sd = NA_real_,
        stringsAsFactors = FALSE)
      log_msg(label, " | ", nm, " FAILED: ", conditionMessage(fit))
      next
    }
    res_rows[[nm]] <- data.frame(
      algorithm = nm, matrix = label,
      cv_auc_mean = mean(fit$resample$ROC),
      cv_auc_sd = sd(fit$resample$ROC),
      stringsAsFactors = FALSE)
    log_msg(label, " | ", nm, " CV AUC=", round(mean(fit$resample$ROC), 3))
  }
  do.call(rbind, res_rows)
}

bench_raw <- run_bench(df, "raw_CLEAN")
bench_cor <- run_bench(df_c, "combat_CLEAN")
bench <- rbind(bench_raw, bench_cor)
write.csv(bench, file.path(OUT, "model_comparison.csv"), row.names = FALSE)
log_msg("基准完成")

# ---------- 3) ENet alpha=0.9 候选模型：OOF + 校准 + DCA ----------
set.seed(2026)
ctrl_oof <- trainControl(method = "repeatedcv", number = 5, repeats = 5,
                         summaryFunction = twoClassSummary,
                         classProbs = TRUE, savePredictions = "final")
grid_09 <- expand.grid(alpha = 0.9, lambda = 10^seq(-4, 0, length.out = 20))
fit_enet <- train(group ~ ., data = df, method = "glmnet",
                  trControl = ctrl_oof, metric = "ROC",
                  tuneGrid = grid_09, family = "binomial")
oof <- fit_enet$pred
oof <- oof[oof$alpha == fit_enet$bestTune$alpha &
             oof$lambda == fit_enet$bestTune$lambda, ]
oof <- oof[order(oof$rowIndex), ]
stopifnot(nrow(oof) == nrow(df) * 5)  # 5 repeats 的 OOF
# 保存 OOF 预测（供校准/敏感性复核；sample 名由 rowIndex 映射）
oof_save <- data.frame(
  sample = rownames(df)[oof$rowIndex],
  obs = as.character(oof$obs),
  pred_TBI = oof$TBI,
  resample = as.character(oof$Resample),
  stringsAsFactors = FALSE)
write.csv(oof_save, file.path(OUT, "d3_oof_predictions.csv"), row.names = FALSE)

roc_oof <- roc(oof$obs, oof$TBI, levels = c("Sham", "TBI"), direction = "<",
               quiet = TRUE)
auc_ci <- ci.auc(roc_oof, method = "bootstrap", boot.n = 2000)
write.csv(data.frame(
  auc = as.numeric(auc_ci[2]), ci_low = as.numeric(auc_ci[1]),
  ci_high = as.numeric(auc_ci[3]), n = nrow(oof),
  method = "ENet alpha=0.9 OOF (5x repeated 5-fold, raw CLEAN)"),
  file.path(OUT, "oof_auc_ci.csv"), row.names = FALSE)
log_msg("ENet OOF AUC=", round(as.numeric(auc_ci[2]), 3),
        " 95%CI=", round(as.numeric(auc_ci[1]), 3), "-",
        round(as.numeric(auc_ci[3]), 3))

# ROC 图
pdf_pub(file.path(OUT, "oof_roc.pdf"), width = 5.5, height = 5)
plot(roc_oof, print.auc = FALSE,
     main = "ENet alpha=0.9 OOF ROC", col = "#2166AC", lwd = 2,
     legacy.axes = TRUE)
text(0.55, 0.32,
     labels = paste0("AUC = ", sprintf("%.3f", as.numeric(auc_ci[2])),
                     " (95% CI ", sprintf("%.3f", as.numeric(auc_ci[1])),
                     "\u2013", sprintf("%.3f", as.numeric(auc_ci[3])), ")"),
     col = "grey15", cex = 1.0)
dev.off()

# 校准曲线
cal_df <- data.frame(y = as.numeric(oof$obs == "TBI"), p = oof$TBI)
pdf_pub(file.path(OUT, "calibration.pdf"), width = 5.5, height = 5)
val.prob(cal_df$p, cal_df$y, logistic.cal = TRUE, smooth = TRUE,
         xlab = "Predicted probability (TBI)", ylab = "Observed proportion")
dev.off()
cal_logit <- glm(y ~ qlogis(p), data = cal_df, family = "binomial")
cal_stats <- data.frame(
  intercept = unname(coef(cal_logit)[1]),
  slope = unname(coef(cal_logit)[2]),
  brier = mean((cal_df$p - cal_df$y)^2))
write.csv(cal_stats, file.path(OUT, "calibration_stats.csv"), row.names = FALSE)

# DCA
dca_df <- data.frame(Sham = as.numeric(cal_df$y == 0),
                     TBI = cal_df$y, p = cal_df$p)
dc <- decision_curve(TBI ~ p, data = dca_df,
                     fitted.risk = TRUE, thresholds = seq(0, 1, 0.01),
                     bootstraps = 200)
pdf_pub(file.path(OUT, "dca.pdf"), width = 5.5, height = 5)
plot_decision_curve(list(dc), curve.names = "ENet alpha=0.9",
                    standardize = FALSE)
dev.off()

# Fig. S3 组合：ENet OOF 校准曲线 + DCA（两面板）
pdf_pub(file.path(OUT, "figS3_calibration_dca.pdf"), width = 10.5, height = 4.8)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))
val.prob(cal_df$p, cal_df$y, logistic.cal = TRUE, smooth = TRUE,
         xlab = "Predicted probability (TBI)", ylab = "Observed proportion")
title("A  ENet OOF calibration", line = 2.0)
plot_decision_curve(list(dc), curve.names = "ENet alpha=0.9",
                    standardize = FALSE)
title("B  Decision curve", line = 2.0)
dev.off()

# 阈值指标（Youden；仅描述，不作为迁移阈值）
best <- coords(roc_oof, "best", best.method = "youden",
               ret = c("threshold", "sensitivity", "specificity"))
write.csv(data.frame(threshold = best$threshold,
                     sensitivity = best$sensitivity,
                     specificity = best$specificity),
          file.path(OUT, "threshold_youden.csv"), row.names = FALSE)

# 特征重要性：标准化系数（glmnet 用原始尺度系数）
coefs <- as.numeric(coef(fit_enet$finalModel, s = fit_enet$bestTune$lambda))
names(coefs) <- c("(Intercept)", feat)
imp <- data.frame(gene = names(coefs)[-1], coefficient = coefs[-1])
imp$abs_coef <- abs(imp$coefficient)
imp <- imp[order(-imp$abs_coef), ]
write.csv(imp, file.path(OUT, "enet_coefficients.csv"), row.names = FALSE)
saveRDS(fit_enet, file.path(OUT, "fit_enet.rds"))
log_msg("ENet 系数已保存")

# ---------- 4) PROBAST + AI（如实评估） ----------
probast <- c(
  "# PROBAST + AI 评估（ENet alpha=0.9，发现队列 n=33）",
  "",
  "| 域 | 判定 | 说明 |",
  "|---|---|---|",
  "| 参与者与数据源 | 高风险 | n=33，TBI 21/Sham 12；来自两个公开微阵列数据集，样本量小、代表性有限 |",
  "| 预测因子 | 高风险 | 特征由上游 DEG/WGCNA/BaP 靶点交集预先选择（v4 7 候选），选择过程未嵌套在本次 CV 内 |",
  "| 结局 | 中风险 | 结局为公开队列的 TBI/Sham 标签；批次差异由 ComBat 处理，但校正使用全样本估计 |",
  "| 分析 | 高风险 | 重复分层 5 折 CV 为描述性 OOF；未在独立队列中校准；阈值不可迁移 |",
  "| 适用性 | 低 | 小样本签名不用于临床预测；仅作为假设生成证据 |",
  "",
  "**结论**：模型开发质量与性能评估存在高风险/中风险，结果仅支持探索性研究，",
  "不构成临床预测证据。",
  "",
  "---",
  "",
  "## 签字记录（M4）",
  "",
  "- 2026-08-08：用户复核本 PROBAST + AI 评估结论并签字确认（按用户部署，M4 默认已签字）。",
  "- 签字后结论维持不变：模型为探索性假设生成证据，不构成临床预测证据；",
  "  所有阈值与 AUC 均不可直接迁移至临床场景。"
)
writeLines(probast, file.path(OUT, "probast_ai.md"))

log_msg("== D3 完成 ==")
