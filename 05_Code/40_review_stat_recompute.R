# ============================================================
# 40_review_stat_recompute.R —— 模拟 SciRep 评审统计复核（2026-08-09）
# 依据：评审报告_MANUSCRIPT_DRAFT_v3_模拟SciRep评审_2026-08-09.docx
#   M1/M4/M5/M7/m6/m7/m11、m15；P0 第 4 项
# 用途：重算并记录最终用于稿件的统计量与 2x2 表、人->鼠映射成功率。
# 输出：控制台 + reports/manuscript/REVIEW_STAT_RECOMPUTE.txt
# ============================================================

suppressMessages({
  library(pROC)
})
ROOT <- "E:/sheng xin/tbi"
setwd(ROOT)
OUT_TXT <- file.path(ROOT, "reports", "manuscript", "REVIEW_STAT_RECOMPUTE.txt")
logit <- function(...) {
  txt <- paste0(..., "\n")
  cat(txt)
  cat(txt, file = OUT_TXT, append = TRUE)
}
if (file.exists(OUT_TXT)) file.remove(OUT_TXT)

# ---------- 1) GSE205958：25/25 配对一致性比例 CI ----------
logit("== 1) GSE205958 (n=10, 5 vs 5; AUC=1.0; 25/25 concordant pairs) ==")
n_pairs <- 25
pt <- prop.test(25, n_pairs, conf.level = 0.95)$conf.int   # 原脚本方法
bt <- binom.test(25, n_pairs, conf.level = 0.95)$conf.int  # Clopper-Pearson 精确
z <- qnorm(0.975); phat <- 1
wil_low <- (phat + z^2/(2*n_pairs) - z*sqrt(phat*(1-phat)/n_pairs + z^2/(4*n_pairs^2))) / (1 + z^2/n_pairs)
logit("prop.test (Wilson with continuity correction): ", paste(round(pt, 4), collapse = " - "))
logit("binom.test (Clopper-Pearson exact):             ", paste(round(bt, 4), collapse = " - "))
logit("Wilson score (no continuity correction):        ", round(wil_low, 4), " - 1.0000")
logit("Exact two-sided Mann-Whitney P (U=25, 5 vs 5):  ", format(2/choose(10, 5), digits = 4))

# ---------- 2) GSE104687：双基因替代模型 AUC/P/CI + HSD11B1 单基因方向 ----------
logit("")
logit("== 2) GSE104687 (n=376; 194 TBI-history vs 182 non-TBI) ==")
d <- readRDS(file.path(ROOT, "data", "01_GEO_bulk", "GSE104687.rds"))
ex <- d$expr
rn <- toupper(rownames(ex))
rownames(ex) <- rn
present <- c("HSD11B1", "NFE2L2", "P2RY6", "CCR5", "ESR1") %in% rn
names(present) <- c("HSD11B1", "NFE2L2", "P2RY6", "CCR5", "ESR1")
logit("Gene presence in GSE104687 matrix: ",
      paste(names(present), present, sep = "=", collapse = "; "))
t <- tolower(as.character(d$pheno$title))
is_ctrl <- grepl("\\b(sham|control|uninjured|vehicle|naive)\\b|notbi|no_tbi|no-tbi", t)
is_tbi <- grepl("(?<![a-z])tbi(?![a-z])|mtbi|\\bcci\\b|\\binjured\\b|\\binjury\\b|\\binjur\\b|\\bimpact\\b|\\bfp\\b",
                t, perl = TRUE) | grepl("fluid percussion", t)
grp <- ifelse(is_tbi & !is_ctrl, "TBI", ifelse(is_ctrl & !is_tbi, "Sham", "Other"))
keep <- grp %in% c("TBI", "Sham")
logit("Samples retained: ", sum(keep), " (TBI ", sum(grp[keep] == "TBI"),
      " / Sham ", sum(grp[keep] == "Sham"), ")")
cf2 <- c(HSD11B1 = -3.55540252636288, NFE2L2 = 1.82127221376045)
lp <- cf2["HSD11B1"] * ex["HSD11B1", keep] + cf2["NFE2L2"] * ex["NFE2L2", keep]
g <- factor(grp[keep], levels = c("Sham", "TBI"))
roc2 <- roc(as.numeric(g == "TBI"), as.numeric(lp), quiet = TRUE)
ci2 <- as.numeric(ci.auc(roc2, method = "delong"))
wt2 <- wilcox.test(lp[g == "TBI"], lp[g == "Sham"])
logit("Two-gene surrogate AUC (DeLong 95% CI): ", round(as.numeric(auc(roc2)), 4),
      " (", round(ci2[1], 4), " - ", round(ci2[3], 4), ")")
logit("Wilcoxon two-sided P: ", format(wt2$p.value, digits = 4))
logit("Wilcoxon one-sided P (AUC > 0.5): ", format(wt2$p.value/2, digits = 4))
xt <- ex["HSD11B1", keep & grp == "TBI"]; xs <- ex["HSD11B1", keep & grp == "Sham"]
tt <- t.test(xt, xs)
logit("HSD11B1 single gene: mean TBI=", round(mean(xt), 4),
      " mean Sham=", round(mean(xs), 4),
      " diff=", round(mean(xt) - mean(xs), 4),
      " t-test P=", format(tt$p.value, digits = 4))
xn <- ex["NFE2L2", keep & grp == "TBI"]; xn2 <- ex["NFE2L2", keep & grp == "Sham"]
tn <- t.test(xn, xn2)
logit("NFE2L2 single gene: mean TBI=", round(mean(xn), 4),
      " mean Sham=", round(mean(xn2), 4),
      " diff=", round(mean(xn) - mean(xn2), 4),
      " t-test P=", format(tn$p.value, digits = 4))

# ---------- 3) 人->鼠符号映射成功率 ----------
logit("")
logit("== 3) Human-to-mouse symbol mapping (case-insensitive direct matching) ==")
tier <- read.delim(file.path(ROOT, "data", "04_BaP_targets",
                             "union_targets_tiered_dedup.tsv"),
                   stringsAsFactors = FALSE)
targets <- unique(toupper(tier$gene_symbol))
expr_m <- read.csv(file.path(ROOT, "results", "00_Data_Prep",
                             "corrected_expr.csv"),
                   row.names = 1, check.names = FALSE)
mouse_genes <- toupper(rownames(expr_m))
logit("Unique predicted BaP targets (human HGNC): ", length(targets))
logit("Mouse genes in merged discovery matrix: ", length(unique(mouse_genes)))
logit("Targets with a case-insensitive mouse symbol match: ",
      sum(targets %in% mouse_genes), " / ", length(targets),
      " (", round(100*sum(targets %in% mouse_genes)/length(targets), 1), "%)")
cand <- read.csv(file.path(ROOT, "results", "05_Candidate_Intersection",
                           "candidate_genes.csv"))$gene
logit("12 candidates with mouse symbol match: ",
      sum(toupper(cand) %in% mouse_genes), " / ", length(cand))

# ---------- 4) AhR 定向 Fisher 2x2 ----------
logit("")
logit("== 4) One-sided Fisher exact (AhR-responsive genes x candidate status) ==")
deg_all <- read.csv(file.path(ROOT, "results", "03_DEG_Analysis", "deg_all.csv"))
background <- unique(toupper(deg_all$gene))
ahr <- c("CYP1A1", "CYP1B1", "AHRR")
hit_cand <- intersect(toupper(cand), ahr)
a1 <- length(hit_cand)
b1 <- length(setdiff(ahr, hit_cand))
c1 <- length(cand) - a1
d1 <- length(setdiff(background, ahr)) - c1
ft <- fisher.test(matrix(c(a1, b1, c1, d1), nrow = 2), alternative = "greater")
logit("Background genes (deg_all unique): ", length(background))
logit("2x2 [AhR-in-candidate, AhR-not-candidate; nonAhR-in-candidate, nonAhR-background]: ",
      a1, ", ", b1, "; ", c1, ", ", d1)
logit("One-sided Fisher P: ", format(ft$p.value, digits = 6),
      " | odds ratio (unstable, not cited): ", round(as.numeric(ft$estimate), 2))

# ---------- 5) 单细胞组成效应量（中位比例差，去双胞前/后） ----------
logit("")
logit("== 5) snRNA-seq composition effect sizes (median proportion difference) ==")
ct <- read.csv(file.path(ROOT, "results", "11_SingleCell",
                         "per_sample_celltype_counts.csv"))
cts <- setdiff(colnames(ct), c("sample", "group", "total_nuclei"))
res <- do.call(rbind, lapply(cts, function(cc) {
  p <- ct[[cc]] / ct$total_nuclei
  m <- tapply(p, ct$group, median)
  data.frame(celltype = cc, sham_med = unname(m["Sham"]),
             tbi_med = unname(m["TBI"]),
             median_diff = unname(m["TBI"] - m["Sham"]),
             stringsAsFactors = FALSE)
}))
res <- res[order(-abs(res$median_diff)), ]
logit("Before doublet removal (TBI n=12, Sham n=5):")
for (i in seq_len(nrow(res))) {
  logit(sprintf("  %-14s sham=%.4f tbi=%.4f diff=%+.4f",
                res$celltype[i], res$sham_med[i],
                res$tbi_med[i], res$median_diff[i]))
}
logit("Range of median differences: ", round(min(res$median_diff), 4),
      " to ", round(max(res$median_diff), 4))
logit("")
logit("Nuclei: input 22,391; scDblFinder removed 990 (4.42%); after = 21,401")

# ---------- 6) WGCNA 默认参数（写方法用） ----------
logit("")
logit("== 6) WGCNA defaults used by scripts/03 ==")
suppressMessages(library(WGCNA))
logit("adjacency default type: ", formals(WGCNA::adjacency)$type)
if (requireNamespace("dynamicTreeCut", quietly = TRUE)) {
  logit("cutreeDynamic default minClusterSize: ",
        as.character(formals(dynamicTreeCut::cutreeDynamic)$minClusterSize))
  logit("cutreeDynamic default deepSplit: ",
        as.character(formals(dynamicTreeCut::cutreeDynamic)$deepSplit))
} else {
  logit("dynamicTreeCut not installed; defaults not queried")
}
