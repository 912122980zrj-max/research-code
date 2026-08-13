# ============================================================
# 26_external_validation_paper.R —— P0-2：论文口径外部验证摘要
# 依据：FINAL_REVIEW_PUBLISHABILITY.md §4/P0-2、§6.1（H2）
# 作用：
#   1) 逐基因方向一致性表（GSE128543 WT-only；GSE205958）；
#   2) 论文口径汇总表：只报 AUC + 方向；删除灵敏/特异度；
#      明确 CI 不可用（n=14）与阈值不可迁移；
#   3) 输出写作建议文本。
# 输出：results/08_External_Validation/
#   per_gene_direction_validation.csv、external_validation_paper_summary.csv、
#   EXTERNAL_VALIDATION_PAPER_NOTE.md
# ============================================================

suppressMessages({
  library(dplyr)
  library(readr)
  library(readxl)
})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else {
  stop("请在项目根目录运行")
}
OUT <- file.path(ROOT, "results", "08_External_Validation")

sig3 <- c("P2RY6", "HSD11B1", "CCR5")  # v4 重跑后 D4 SHAP top3

# 期望方向 = CLEAN_3gene_only 冻结系数的符号（先运行 04b）
coef3 <- read.csv(file.path(OUT, "frozen_coefficients.csv"),
                  stringsAsFactors = FALSE)
cf3 <- coef3[coef3$model == "CLEAN_3gene_only", ]
sign_map <- setNames(sign(cf3$coefficient[cf3$term != "(Intercept)"]),
                     cf3$term[cf3$term != "(Intercept)"])

read_bulk <- function(gse) {
  readRDS(file.path(ROOT, "data", "01_GEO_bulk", paste0(gse, ".rds")))
}

wt_test <- function(x_tbi, x_sham) {
  v <- tryCatch(t.test(x_tbi, x_sham)$p.value, error = function(e) NA_real_)
  v
}

# ---------- GSE128543 WT-only ----------
d128 <- read_bulk("GSE128543")
ex128 <- d128$expr
rn <- toupper(rownames(ex128))
ex128 <- ex128[rn %in% sig3, , drop = FALSE]
rownames(ex128) <- toupper(rownames(ex128))
gt <- if ("genotype:ch1" %in% colnames(d128$pheno)) {
  as.character(d128$pheno$`genotype:ch1`)
} else rep("", ncol(ex128))
is_wt <- grepl("^\\s*wt\\s*$|wildtype", gt, ignore.case = TRUE) |
  grepl("\\bwt\\b|wildtype", d128$pheno$title, ignore.case = TRUE)
title <- tolower(d128$pheno$title)
grp <- ifelse(grepl("(?<![a-z])tbi(?![a-z])|\\bcci\\b", title, perl = TRUE) &
                !grepl("sham", title), "TBI",
              ifelse(grepl("sham", title), "Sham", NA))
keep <- is_wt & grp %in% c("TBI", "Sham")
dir128 <- do.call(rbind, lapply(sig3, function(g) {
  if (!g %in% rownames(ex128)) {
    return(data.frame(validation = "GSE128543_WT_only", gene = g,
                      mean_tbi = NA, mean_sham = NA, diff = NA,
                      t_p = NA, direction_consistent = NA,
                      stringsAsFactors = FALSE))
  }
  x_tbi <- ex128[g, keep & grp == "TBI"]
  x_sham <- ex128[g, keep & grp == "Sham"]
  diff <- mean(x_tbi) - mean(x_sham)
  data.frame(validation = "GSE128543_WT_only", gene = g,
             mean_tbi = round(mean(x_tbi), 4),
             mean_sham = round(mean(x_sham), 4),
             diff = round(diff, 4),
             t_p = round(wt_test(x_tbi, x_sham), 6),
             direction_consistent =
               if (g %in% names(sign_map)) sign(diff) == sign_map[[g]] else NA,
             stringsAsFactors = FALSE)
}))

# ---------- GSE205958（直接读 FPKM xlsx，列序 TBI1-5, sham1-5） ----------
xlsx_file <- file.path(ROOT, "data", "01_GEO_bulk",
                       "GSE205958_genes_fpkm_expression.xlsx")
dat5 <- as.data.frame(read_excel(xlsx_file), stringsAsFactors = FALSE,
                      check.names = FALSE)
cn5 <- colnames(dat5)
fpkm_cols <- which(grepl("^FPKM\\.", cn5))
stopifnot(length(fpkm_cols) == 10)
rn5 <- toupper(as.character(dat5$gene_name))
ex205 <- as.matrix(dat5[, fpkm_cols])
rownames(ex205) <- rn5
ex205 <- ex205[!is.na(rn5) & rn5 != "", , drop = FALSE]
if (anyDuplicated(rownames(ex205))) {
  dup_n5 <- table(rownames(ex205))
  ex205 <- rowsum(ex205, group = rownames(ex205))
  ex205 <- ex205 / matrix(as.numeric(dup_n5[rownames(ex205)]),
                          nrow = nrow(ex205), ncol = ncol(ex205),
                          byrow = FALSE)
}
tbi_idx <- 1:5
sham_idx <- 6:10
dir205 <- do.call(rbind, lapply(sig3, function(g) {
  x_tbi <- ex205[g, tbi_idx]
  x_sham <- ex205[g, sham_idx]
  diff <- mean(x_tbi) - mean(x_sham)
  data.frame(validation = "GSE205958", gene = g,
             mean_tbi = round(mean(x_tbi), 4),
             mean_sham = round(mean(x_sham), 4),
             diff = round(diff, 4),
             t_p = round(wt_test(x_tbi, x_sham), 6),
             direction_consistent =
               if (g %in% names(sign_map)) sign(diff) == sign_map[[g]] else NA,
             stringsAsFactors = FALSE)
}))

dir_tab <- rbind(dir128, dir205)
write.csv(dir_tab, file.path(OUT, "per_gene_direction_validation.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("===== per-gene direction =====\n")
print(dir_tab)

# ---------- 论文口径汇总 ----------
fsv <- read.csv(file.path(OUT, "frozen_signature_validation.csv"),
                stringsAsFactors = FALSE)
row104 <- fsv[fsv$validation == "GSE104687" &
                fsv$signature_training == "CLEAN_2gene_frozen", ]
auc104 <- if (nrow(row104) == 1) row104$auc[1] else NA_real_
ci104 <- if (nrow(row104) == 1) {
  sprintf("%.4f–%.4f", row104$ci_low[1], row104$ci_high[1])
} else "NA"
summary_tab <- data.frame(
  validation = c("GSE128543 (WT-only)", "GSE205958",
                 "GSE104687 (transportability)"),
  model = c("Frozen 3-gene (CLEAN; P2RY6/HSD11B1/CCR5)",
            "Frozen 3-gene (CLEAN; P2RY6/HSD11B1/CCR5)",
            "not evaluable (P2RY6/CCR5 absent)"),
  n = c(14, 10, 376),
  n_tbi = c(7, 5, 194),
  n_sham = c(7, 5, 182),
  auc = c(1.0, 1.0, round(auc104, 4)),
  ci = c("not interpretable (n=14, perfect separation)",
         "exact binomial (Clopper–Pearson) 95% CI 0.863–1.000 (25/25 non-independent concordant pairs; descriptive)",
         ci104),
  per_gene_direction = c(
    "2/3 evaluable consistent (P2RY6 up, HSD11B1 down); CCR5 opposite in both external cohorts",
    "2/3 evaluable consistent (P2RY6 up, HSD11B1 down); CCR5 opposite in both external cohorts",
    "P2RY6/CCR5 absent; only HSD11B1 evaluable -> not evaluable under v4 signature"),
  interpretation = c(
    "Rank-order AUC + direction only; threshold NOT transferable; sensitivity/specificity NOT reported",
    "Exploratory replication; small n; cross-platform FPKM vs microarray; 7 d timepoint",
    "v4 signature not evaluable in this human aging chronic cohort; non-human validation excluded"),
  stringsAsFactors = FALSE
)
write.csv(summary_tab, file.path(OUT, "external_validation_paper_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("===== paper summary =====\n")
print(summary_tab)

# ---------- 写作建议 ----------
md <- c(
  "# 外部验证论文口径说明（P0-2）",
  "",
  paste("> 生成时间：", Sys.time()),
  "> 依据：FINAL_REVIEW_PUBLISHABILITY.md §4/P0-2、§6.1；本文件供写作线程直接引用。",
  "",
  "## 1. 必写表述",
  "",
  "- GSE128543（WT-only，n=14，7 vs 7）：固定 3 特征签名（P2RY6/HSD11B1/CCR5，",
  "  v4 D4 SHAP top3）的**排序 AUC=1.0**；",
  "  **不报灵敏度/特异度**（训练阈值跨平台迁移无意义：原 sens=1/spec=0）；",
  "  Delong/bootstrap CI 在完全分离时退化，**CI 不可用**；",
  "- 逐基因方向：P2RY6↑、HSD11B1↓ 在 GSE128543 与 GSE205958 均一致",
  "  （2/3 可评估）；CCR5 观察方向均为↑，而 v4 冻结系数为负（方向相反），",
  "  必须如实写为“不稳定，不作为方向性复制证据”；",
  "- GSE205958（n=10，5 vs 5）：AUC=1.0（25/25 配对，精确二项（Clopper–Pearson）",
  "  95% CI 0.863–1；配对非独立，区间描述性）；",
  "  探索性复制，跨平台（微阵列→FPKM）、跨时点（3 d/48 h→7 d）、",
  "  naïve 对照；",
  "- GSE104687：v4 固定 3 基因签名不可评估（P2RY6/CCR5 在该 FPKM 队列缺失）；",
  "  该队列只能作为跨物种/跨时相局限性的“不可评估”证据。",
  "",
  "## 2. 禁止表述",
  "",
  "- “外部验证显著”/“AUC=1.0 证明签名有效”；",
  "- 报告 sens/spec 或阈值迁移；",
  "- 将 GSE104687 称为验证失败/成功——只称“跨物种/跨队列迁移性阴性”。",
  "",
  "## 3. 数据文件",
  "",
  "- `per_gene_direction_validation.csv`（本脚本重算，可复现）；",
  "- `external_validation_paper_summary.csv`（论文表格底稿）；",
  "- 原始冻结验证（含阈值/sens/spec，仅归档不引用）：`frozen_signature_validation.csv`、",
  "  `external_validation_v2.csv`、`D12_GSE205958/external_validation_v3.csv`。"
)
writeLines(md, file.path(OUT, "EXTERNAL_VALIDATION_PAPER_NOTE.md"))
message("完成。结果在 ", OUT)
