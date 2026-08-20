# ============================================================
# 20_d4_shap.R —— D4 真 SHAP 解释（ENet alpha=0.9）
# 依据：模板论文 SHAP 四件套 + ml-shap 技能
# 方法：ENet 为线性模型，SHAP 值可精确计算：
#   phi_j = beta_j * (x_j - E[x_j])，baseline = b0 + sum(beta_j * E[x_j])
#   预测类别为 TBI（log-odds 空间）。
# 说明（2026-08-08 终审修复）：SHAP 解释 **D3 保存的同一 ENet 拟合**
#   （scripts/19_d3_ml_benchmark.R 输出 fit_enet.rds，caret repeated CV 选 λ，
#   seed=2026），不再另做 cv.glmnet 重拟合，避免系数口径不一致；
#   属样本内归因，不做因果解释。
# 输出：results/07_ML_Signature/D4_shap/
# ============================================================
suppressMessages({
  library(glmnet); library(shapviz); library(ggplot2)
  library(dplyr); library(readr)
})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "results", "07_ML_Signature", "D4_shap")
OUT_D3 <- file.path(ROOT, "results", "07_ML_Signature", "D3_ml_benchmark")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts", "theme_pub.R"))
set.seed(2026)

# ---------- 数据（与 D3 主口径一致：原始 CLEAN，v4 7 候选，log2 尺度） ----------
cand <- toupper(read.csv(file.path(ROOT, "results", "05_Candidate_Intersection",
                                   "candidate_genes.csv"))$gene)
comp <- read.csv(file.path(ROOT, "results", "00_Data_Prep",
                           "sample_composition.csv"), check.names = FALSE)
clean_samples <- comp$sample[comp$include_clean]
grp_map <- setNames(
  ifelse(comp$group[match(clean_samples, comp$sample)] == "tbi",
         "TBI", "Sham"), clean_samples)
read_bulk <- function(gse) {
  d <- readRDS(file.path(ROOT, "data", "01_GEO_bulk", paste0(gse, ".rds")))
  list(expr = d$expr, pheno = d$pheno)
}
dat <- lapply(c("GSE58485", "GSE41345"), read_bulk)
common_genes <- Reduce(intersect, lapply(dat, function(d) rownames(d$expr)))
expr_all <- do.call(cbind, lapply(dat, function(d) d$expr[common_genes, ]))
expr_all <- expr_all[, intersect(clean_samples, colnames(expr_all)), drop = FALSE]
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
x <- t(ex)
y <- as.numeric(grp_map[rownames(x)] == "TBI")
stopifnot(!anyNA(y))

# ---------- 加载 D3 最终 ENet（同一拟合） ----------
fit_file <- file.path(OUT_D3, "fit_enet.rds")
if (!file.exists(fit_file)) stop("缺少 fit_enet.rds，请先运行 19_d3_ml_benchmark.R")
fit_enet <- readRDS(fit_file)
stopifnot(fit_enet$bestTune$alpha == 0.9)
b <- as.numeric(coef(fit_enet$finalModel, s = fit_enet$bestTune$lambda))
names(b) <- c("(Intercept)", feat)
write.csv(data.frame(gene = names(b)[-1], coefficient = b[-1]),
          file.path(OUT, "enet_coefficients.csv"), row.names = FALSE)

# ---------- 精确线性 SHAP ----------
mu <- colMeans(x)
phi <- sweep(sweep(x, 2, mu, "-"), 2, b[-1], "*")
baseline <- b[1] + sum(b[-1] * mu)
link_pred <- as.numeric(predict(fit_enet$finalModel,
                                s = fit_enet$bestTune$lambda,
                                newx = x, type = "link"))
max_dev <- max(abs(rowSums(phi) + baseline - link_pred))
cat("SHAP 重构最大偏差:", max_dev, "\n")
stopifnot(max_dev < 1e-6)

shp <- shapviz(phi, X = as.data.frame(x, check.names = FALSE),
               baseline = baseline)
saveRDS(shp, file.path(OUT, "shap_object.rds"))
write.csv(cbind(sample = rownames(x), as.data.frame(phi, check.names = FALSE)),
          file.path(OUT, "shap_values.csv"), row.names = FALSE)

# ---------- 四件套 ----------
top3 <- names(sort(colMeans(abs(phi)), decreasing = TRUE))[1:3]
cat("top3 genes:", paste(top3, collapse = ", "), "\n")

make_plot <- function(gg, filebase, w = 7, h = 5.5) {
  ggsave_pub(file.path(OUT, paste0(filebase, ".pdf")), gg, width = w, height = h)
}

# 1) 全局重要性（柱状 + beeswarm）
p1 <- sv_importance(shp, kind = "both", max_display = 12) +
  labs(title = "SHAP importance (ENet alpha=0.9, class=TBI)")
make_plot(p1, "shap_importance_beeswarm", w = 8.5, h = 6)

# 2) beeswarm 单独
p2 <- sv_importance(shp, kind = "beeswarm", max_display = 15) +
  labs(title = "SHAP beeswarm (class=TBI)")
make_plot(p2, "shap_beeswarm", w = 7, h = 5.5)

# 3) 依赖图（top3）
for (v in top3) {
  p3 <- sv_dependence(shp, v = v) +
    labs(title = paste0("SHAP dependence: ", v))
  make_plot(p3, paste0("shap_dependence_", v))
}

# 4) waterfall：代表性样本（高置信 TBI / 高置信 Sham / 边界样本）
pr <- as.numeric(predict(fit_enet$finalModel,
                         s = fit_enet$bestTune$lambda,
                         newx = x, type = "response"))
names(pr) <- rownames(x)
rep_tbi <- rownames(x)[which(y == 1)[which.max(pr[y == 1])]]
rep_sham <- rownames(x)[which(y == 0)[which.min(pr[y == 0])]]
border_cands <- names(sort(abs(pr - 0.5)))
border <- border_cands[!border_cands %in% c(rep_tbi, rep_sham)][1]
write.csv(data.frame(role = c("rep_tbi", "rep_sham", "border"),
                     sample = c(rep_tbi, rep_sham, border)),
          file.path(OUT, "waterfall_samples.csv"), row.names = FALSE)
cat("waterfall samples:", rep_tbi, "/", rep_sham, "/", border, "\n")
for (sid in c(rep_tbi, rep_sham, border)) {
  p4 <- sv_waterfall(shp, row_id = sid) +
    labs(title = paste0("SHAP waterfall: ", sid,
                        " (p_TBI=", round(pr[sid], 3), ")"))
  make_plot(p4, paste0("shap_waterfall_", sid))
}

writeLines(c(
  "# D4 SHAP 说明",
  "",
  "- 模型：ENet alpha=0.9（D3 最终拟合 fit_enet.rds，caret repeated CV 选 lambda，",
  "  seed=2026），原始 CLEAN 矩阵（n=33，v4 7 特征）；",
  "- SHAP 为精确线性 SHAP（log-odds 空间），类别=TBI；",
  "- 背景集：发现队列全样本（样本内归因，不解释为因果）；",
  "- 重构检查：|rowSums(phi)+baseline - link_pred| < 1e-6。"
), file.path(OUT, "SHAP_NOTES.md"))

cat("== D4 完成 ==", "\n")
