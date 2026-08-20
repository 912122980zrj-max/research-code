# ============================================================
# 31_h3_cohort_sensitivity.R —— H3/P0-4：发现队列异质性敏感性（v4）
# 依据：FINAL_REVIEW_PUBLISHABILITY.md §2.2 H3、§3.3、P0-4/P1；
#       v4 冻结计划（2026-08-09）
# 内容：
#   1) 每队列单独 DEG（limma，同阈值），预处理与合并口径完全一致：
#      公共基因 → CLEAN 样本 → NA 过滤（<20%）→ 每队列覆盖过滤
#      （默认两队列均 >=10% 样本非零）→ 中位数填补；
#      输出"有基因型 / 无基因型"两种口径（GSE41345 无基因型列时默认
#      wildtype，可用 GENO_CODE=other 复现敏感性）；
#   2) 单队列建模（3 基因与 5 基因核心，candidate_genes_FULL.csv = FULL_5gene）
#      → 对另一队列评估 AUC；与 04 脚本的合并口径 LOCO（cross_cohort_auc.csv）
#      互为独立敏感性，参数不同（5×3 vs 5×5 CV）。
#   3) 输出报告，供论文 limitation/敏感性章节使用。
# 输出：results/07_ML_Signature/H3_cohort_sensitivity/
# ============================================================

suppressMessages({
  library(limma)
  library(caret)
  library(glmnet)
  library(pROC)
  library(dplyr)
  library(readr)
})

ROOT <- getwd()
OUT <- file.path(ROOT, "results", "07_ML_Signature", "H3_cohort_sensitivity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

GENO_CODE <- Sys.getenv("GENO_CODE", unset = "wildtype")
COV_THRESHOLD <- as.numeric(Sys.getenv("COV_THRESHOLD", unset = "0.10"))

read_bulk <- function(gse) readRDS(file.path(ROOT, "data", "01_GEO_bulk",
                                             paste0(gse, ".rds")))

core6_file <- file.path(ROOT, "results", "05_Candidate_Intersection",
                        "candidate_genes_FULL.csv")
core <- toupper(read.csv(core6_file)$gene)
sig3 <- c("P2RY6", "HSD11B1", "CCR5")  # v4 重跑后 D4 SHAP top3

comp <- read.csv(file.path(ROOT, "results", "00_Data_Prep",
                           "sample_composition.csv"), check.names = FALSE)
clean_samples <- comp$sample[comp$include_clean]

# ---------- 与 03 合并口径一致的预处理 ----------
preprocess_discovery <- function() {
  gses <- c("GSE58485", "GSE41345")
  dat <- lapply(gses, read_bulk)
  names(dat) <- gses
  common <- Reduce(intersect, lapply(dat, function(d) rownames(d$expr)))
  expr_all <- do.call(cbind, lapply(gses, function(g) {
    dat[[g]]$expr[common, , drop = FALSE]
  }))
  expr_all <- expr_all[, intersect(clean_samples, colnames(expr_all)),
                       drop = FALSE]
  cohort_vec <- comp$cohort[match(colnames(expr_all), comp$sample)]
  grp_map <- setNames(comp$group[match(colnames(expr_all), comp$sample)],
                      colnames(expr_all))

  keep_gene <- rowSums(is.na(expr_all)) < 0.2 * ncol(expr_all)
  expr_all <- expr_all[keep_gene, , drop = FALSE]
  cov_a <- rowMeans(expr_all[, cohort_vec == "GSE58485", drop = FALSE] > 0,
                    na.rm = TRUE)
  cov_b <- rowMeans(expr_all[, cohort_vec == "GSE41345", drop = FALSE] > 0,
                    na.rm = TRUE)
  keep_cov <- cov_a >= COV_THRESHOLD & cov_b >= COV_THRESHOLD
  expr_all <- expr_all[keep_cov, , drop = FALSE]
  expr_all <- apply(expr_all, 2, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    x
  })
  list(expr = expr_all, cohort = cohort_vec, group = grp_map)
}

pp <- preprocess_discovery()
cat("preprocessed genes:", nrow(pp$expr),
    "| samples:", ncol(pp$expr), "\n")

# ---------- 1) 每队列 DEG（有/无基因型两种口径） ----------
run_cohort_deg <- function(gse, use_geno) {
  keep <- pp$cohort == gse
  ex <- pp$expr[, keep, drop = FALSE]
  ph <- data.frame(
    group = factor(ifelse(pp$group[keep] == "tbi", "TBI", "Sham"),
                   levels = c("Sham", "TBI")),
    stringsAsFactors = FALSE
  )
  d <- read_bulk(gse)
  geno_note <- "genotype not available (no genotype:ch1)"
  if (use_geno && "genotype:ch1" %in% colnames(d$pheno)) {
    gt <- as.character(d$pheno$`genotype:ch1`)[match(colnames(ex),
                                                     rownames(d$pheno))]
    gt[is.na(gt) | gt == ""] <- GENO_CODE
    ph$gt <- factor(gt)
    if (nlevels(ph$gt) < 2) {
      design <- model.matrix(~ group, data = ph)
      geno_note <- "genotype covariate dropped (constant)"
    } else {
      design_full <- model.matrix(~ group + gt, data = ph)
      if (qr(design_full)$rank < ncol(design_full)) {
        design <- model.matrix(~ group, data = ph)
        geno_note <- "genotype covariate dropped (collinear)"
      } else {
        design <- design_full
        geno_note <- paste("genotype covariate:", paste(levels(ph$gt),
                                                        collapse = ","))
      }
    }
  } else {
    if (use_geno) {
      ph$gt <- factor(rep(GENO_CODE, nrow(ph)))
      geno_note <- paste0("genotype defaulted to ", GENO_CODE,
                          " (constant, dropped)")
    }
    design <- model.matrix(~ group, data = ph)
  }
  fit <- lmFit(ex, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = "groupTBI", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt$cohort <- gse
  tt$significant <- tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.585
  cat(gse, "use_geno:", use_geno, "| samples:", ncol(ex),
      "| DEGs:", sum(tt$significant), "|", geno_note, "\n")
  list(tt = tt, note = geno_note)
}

deg_tabs <- list()
for (gse in c("GSE58485", "GSE41345")) {
  for (ug in c(TRUE, FALSE)) {
    res <- run_cohort_deg(gse, ug)
    tag <- if (ug) "with_geno" else "no_geno"
    deg_tabs[[paste(gse, tag)]] <- res$tt
    write.csv(res$tt, file.path(OUT, paste0("cohort_deg_", gse, "_",
                                            tag, ".csv")),
              row.names = FALSE)
    if (ug) {
      write.csv(res$tt, file.path(OUT, paste0("cohort_deg_", gse, ".csv")),
                row.names = FALSE)
    }
  }
}

combined <- read.csv(file.path(ROOT, "results", "03_DEG_Analysis",
                               "deg_significant.csv"))$gene
comp_n <- read.csv(file.path(ROOT, "results", "00_Data_Prep",
                             "sample_composition.csv"), check.names = FALSE)
n_samp <- table(comp_n$cohort[comp_n$include_clean])

overlap_rows <- list()
for (ug in c(TRUE, FALSE)) {
  tag <- if (ug) "with_geno" else "no_geno"
  set_a <- deg_tabs[[paste("GSE58485", tag)]]$gene
  set_b <- deg_tabs[[paste("GSE41345", tag)]]$gene
  set_a <- set_a[deg_tabs[[paste("GSE58485", tag)]]$significant]
  set_b <- set_b[deg_tabs[[paste("GSE41345", tag)]]$significant]
  overlap <- intersect(set_a, set_b)
  overlap_rows[[tag]] <- data.frame(
    design = tag,
    n_samples_a = as.integer(n_samp[["GSE58485"]]),
    n_samples_b = as.integer(n_samp[["GSE41345"]]),
    n_deg_a = length(set_a),
    n_deg_b = length(set_b),
    n_overlap = length(overlap),
    jaccard = round(length(overlap) / length(union(set_a, set_b)), 4),
    n_a_in_combined = length(intersect(set_a, combined)),
    n_b_in_combined = length(intersect(set_b, combined)),
    n_combined = length(combined),
    pct_combined_in_a = round(100 * length(intersect(set_a, combined)) /
                                length(combined), 2),
    stringsAsFactors = FALSE
  )
}
overlap_tab <- do.call(rbind, overlap_rows)
write.csv(overlap_tab, file.path(OUT, "cohort_deg_overlap.csv"),
          row.names = FALSE)
print(overlap_tab)

# 3 基因方向（有基因型口径，主口径）
dir_tab <- do.call(rbind, lapply(c("GSE58485", "GSE41345"), function(gse) {
  tt <- deg_tabs[[paste(gse, "with_geno")]]
  do.call(rbind, lapply(sig3, function(g) {
    r <- tt[toupper(tt$gene) == g, ]
    if (nrow(r) == 0) return(NULL)
    data.frame(cohort = gse, gene = g,
               logFC = round(r$logFC, 4),
               p = format(r$P.Value, digits = 3),
               adjP = format(r$adj.P.Val, digits = 3),
               stringsAsFactors = FALSE)
  }))
}))
write.csv(dir_tab, file.path(OUT, "signature_gene_direction_by_cohort.csv"),
          row.names = FALSE)
cat("===== 3-gene direction by cohort (with_geno) =====\n")
print(dir_tab)

# ---------- 2) 单队列建模 → 另一队列 ----------
build_df <- function(gse, genes, samples) {
  genes <- toupper(genes)
  d <- read_bulk(gse)
  ex <- d$expr
  rn <- toupper(rownames(ex))
  ex2 <- ex[rn %in% genes, , drop = FALSE]
  rownames(ex2) <- toupper(rownames(ex2))
  if (anyDuplicated(rownames(ex2))) {
    ex2 <- rowsum(ex2, group = rownames(ex2))
  }
  ex2 <- ex2[, colnames(ex2) %in% samples, drop = FALSE]
  ph <- d$pheno
  ph$group <- ifelse(grepl("mtbi|\\btbi\\b|\\bcci\\b|\\binjured\\b|\\binjury\\b|\\bimpact\\b",
                           tolower(ph$title)) &
                       !grepl("sham|vehicle|control|ex-4|uninjured",
                              tolower(ph$title)),
                     "TBI", "Sham")
  ph <- ph[ph$group %in% c("TBI", "Sham"), ]
  ex2 <- ex2[, colnames(ex2) %in% rownames(ph), drop = FALSE]
  ph <- ph[colnames(ex2), ]
  cat(gse, "group table:", paste(names(table(ph$group)),
                                 table(ph$group), sep = "=",
                                 collapse = ", "), "\n")
  list(expr = ex2, pheno = ph)
}

fit_one <- function(tr, genes) {
  tr <- tr[, c(intersect(genes, colnames(tr)), "group"), drop = FALSE]
  tr <- tr[complete.cases(tr), ]
  set.seed(2026)
  train(group ~ ., data = tr, method = "glmnet",
        trControl = trainControl(method = "repeatedcv", number = 5,
                                 repeats = 3, classProbs = TRUE,
                                 summaryFunction = twoClassSummary),
        metric = "ROC", tuneLength = 3)
}

model_rows <- list()
for (gse_tr in c("GSE58485", "GSE41345")) {
  gse_te <- setdiff(c("GSE58485", "GSE41345"), gse_tr)
  for (label in c("3gene", "core")) {
    genes <- if (label == "3gene") sig3 else core
    tr <- build_df(gse_tr, genes, clean_samples)
    te <- build_df(gse_te, genes, clean_samples)
    tr_df <- as.data.frame(t(tr$expr)); tr_df$group <- tr$pheno$group
    te_df <- as.data.frame(t(te$expr)); te_df$group <- te$pheno$group
    common <- intersect(colnames(tr_df), colnames(te_df))
    common <- setdiff(common, "group")
    if (length(common) < 2) next
    m <- tryCatch(fit_one(tr_df[, c(common, "group")], common),
                  error = function(e) NULL)
    if (is.null(m)) next
    prob <- tryCatch(predict(m, newdata = te_df[, c(common, "group")],
                             type = "prob")[, "TBI"],
                     error = function(e) NULL)
    if (is.null(prob)) next
    auc_te <- as.numeric(auc(roc(as.numeric(te_df$group == "TBI"),
                                 prob, quiet = TRUE)))
    model_rows[[length(model_rows) + 1]] <- data.frame(
      train_cohort = gse_tr, test_cohort = gse_te, model = label,
      n_train = nrow(tr_df), n_test = nrow(te_df),
      n_features = length(common),
      train_cv_auc = round(max(m$results$ROC, na.rm = TRUE), 4),
      test_auc = round(auc_te, 4),
      stringsAsFactors = FALSE)
  }
}
model_tab <- do.call(rbind, model_rows)
write.csv(model_tab, file.path(OUT, "single_cohort_model_auc.csv"),
          row.names = FALSE)
cat("===== single-cohort models =====\n")
print(model_tab)

# ---------- 合并口径 LOCO（动态读取 04 输出） ----------
loco <- tryCatch(read.csv(file.path(ROOT, "results", "07_ML_Signature",
                                    "cross_cohort_auc.csv")),
                 error = function(e) NULL)
loco_txt <- "（04 脚本未运行）"
if (!is.null(loco)) {
  loco_txt <- paste(
    apply(loco, 1, function(r) {
      paste0(r["composition"], " ", r["train_cohort"], "->", r["test_cohort"],
             " AUC=", r["test_auc"])
    }), collapse = "；")
}

model_txt <- paste(
  apply(model_tab, 1, function(r) {
    paste0(r["train_cohort"], "->", r["test_cohort"], " ",
           r["model"], " AUC=", r["test_auc"])
  }), collapse = "；")

# ---------- 报告 ----------
ov_w <- overlap_tab[overlap_tab$design == "with_geno", ]
ov_n <- overlap_tab[overlap_tab$design == "no_geno", ]
md <- c(
  "# H3 发现队列异质性敏感性报告（v4）",
  "",
  paste("> 生成时间：", Sys.time()),
  paste("> 预处理：公共基因 → CLEAN 样本 → NA 过滤 → 每队列覆盖过滤",
        "（>=", COV_THRESHOLD, "样本非零）→ 中位数填补；与 03 合并口径一致。"),
  paste("> GENO_CODE：", GENO_CODE, "（GSE41345 无基因型列默认值）"),
  "> 依据：FINAL_REVIEW_PUBLISHABILITY.md §2.2 H3 / P0-4",
  "",
  "## 1. 每队列单独 DEG（同阈值 |log2FC|>0.585, adjP<0.05）",
  "",
  paste0("- GSE58485（CCI 皮层 3 d，CLEAN 样本 n=", ov_w$n_samples_a,
         "）：有基因型口径单队列 DEG n=", ov_w$n_deg_a,
         "；无基因型口径 n=", ov_n$n_deg_a),
  paste0("- GSE41345（mTBI 海马 7 d，CLEAN 样本 n=", ov_w$n_samples_b,
         "）：有基因型口径单队列 DEG n=", ov_w$n_deg_b,
         "；无基因型口径 n=", ov_n$n_deg_b),
  paste0("- 有基因型口径：两队列重叠 ", ov_w$n_overlap,
         "（Jaccard=", ov_w$jaccard, "）；合并口径 DEG（n=",
         ov_w$n_combined, "）与 GSE58485 单队列 DEG 重叠：",
         ov_w$n_a_in_combined, "（", ov_w$pct_combined_in_a,
         "%）；与 GSE41345 单队列 DEG 重叠：", ov_w$n_b_in_combined),
  paste0("- 无基因型口径：两队列重叠 ", ov_n$n_overlap,
         "（Jaccard=", ov_n$jaccard, "）；合并口径 DEG 与 GSE58485",
         " 单队列 DEG 重叠：", ov_n$n_a_in_combined, "（",
         ov_n$pct_combined_in_a, "%）；与 GSE41345 单队列 DEG 重叠：",
         ov_n$n_b_in_combined),
  paste0("- 解读（v4 修正）：在覆盖过滤 + wildtype 基因型编码下，",
         "合并 DEG 实质由 GSE58485 单队列 DEG 驱动（重叠 ",
         ov_w$pct_combined_in_a, "%）；GSE41345 单队列 0 DEG（n=9、",
         "4 vs 5，功效极低）。此前 v3 报告的“7/570 高度敏感”表述源于",
         "平台覆盖伪影（未做覆盖过滤）与 other 基因型编码，已删除；",
         "LOCO 低 AUC 的结构性原因需同时考虑队列功效不平衡、",
         "模型/组织/时点异质性，三者与批次效应无法严格区分。"),
  "",
  "## 2. 3 基因方向（每队列，有基因型口径）",
  "",
  "见 `signature_gene_direction_by_cohort.csv`；若某基因在某队列方向相反或",
  "不显著，须在论文中如实写出。",
  "**注意（v4）**：P2RY6↑ / HSD11B1↓ / CCR5↑ 在 GSE58485 单队列均显著",
  "（P2RY6/HSD11B1/CCR5 均 |log2FC|>0.585）；GSE41345 三基因均不显著，",
  "其中 CCR5 观察方向为↓（与合并口径冻结系数一致但与 GSE58485 单队列及",
  "两个外部队列的↑方向相反），必须在论文中标注为不稳定基因。",
  "",
  "## 3. 单队列建模互验 AUC",
  "",
  "见 `single_cohort_model_auc.csv`；解读：",
  "- 单队列训练样本量更小（GSE41345 CLEAN n=9），互验 AUC 波动大；",
  paste0("- 合并口径 LOCO（04 脚本，`results/07_ML_Signature/",
         "cross_cohort_auc.csv`）：", loco_txt, "；论文以该口径为准；"),
  paste0("- 本表单队列建模互验（3 基因/5 基因核心，5×3 CV，seed=2026）：",
         model_txt, "；为独立敏感性分析，与合并口径 LOCO 参数不同，",
         "两者共同说明跨队列泛化有限；"),
  "- 论文表述建议：LOCO 与单队列互验均不构成泛化证据，只作异质性描述。",
  "",
  "## 4. 局限",
  "",
  "- 单队列 DEG 使用各自队列内样本（n=24/9），检验功效低；",
  "- GSE41345 仅 4 TBI vs 5 Sham，单队列建模极不稳定，结果仅作敏感性；",
  paste0("- 覆盖过滤阈值 ", COV_THRESHOLD,
         " 为 v4 主口径；0 / 0.5 阈值敏感性可通过 COV_THRESHOLD 环境变量复现。")
)
writeLines(md, file.path(OUT, "COHORT_SENSITIVITY_REPORT.md"))
message("完成。结果在 ", OUT)
