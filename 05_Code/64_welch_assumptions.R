# ============================================================
# 64_welch_assumptions.R —— A-6：逐基因方向检验假设检查
# 依据：修改报告_BaP_TBI_最终版.md A-6
#   方法中写明 Welch's t-test；补方差齐性与正态性假设检查及稳健性说明。
# 输出：results/08_External_Validation/assumption_checks.csv
#       results/08_External_Validation/ASSUMPTION_CHECKS.md
# ============================================================

suppressMessages({
  library(dplyr)
  library(readr)
  library(readxl)
})

ROOT <- getwd()
OUT <- file.path(ROOT, "results", "08_External_Validation")

read_bulk <- function(gse) readRDS(file.path(ROOT, "data", "01_GEO_bulk",
                                             paste0(gse, ".rds")))

checks <- list()
add <- function(x) checks[[length(checks) + 1]] <<- x

# ---------- GSE128543 WT-only ----------
d128 <- read_bulk("GSE128543")
ex128 <- d128$expr
rn <- toupper(rownames(ex128))
ex128 <- ex128[rn %in% c("P2RY6", "HSD11B1", "CCR5"), , drop = FALSE]
rownames(ex128) <- toupper(rownames(ex128))
title <- tolower(d128$pheno$title)
grp <- ifelse(grepl("(?<![a-z])tbi(?![a-z])|\\bcci\\b", title, perl = TRUE) &
                !grepl("sham", title), "TBI",
              ifelse(grepl("sham", title), "Sham", NA))
gt <- if ("genotype:ch1" %in% colnames(d128$pheno)) {
  as.character(d128$pheno$`genotype:ch1`)
} else rep("", ncol(ex128))
is_wt <- grepl("^\\s*wt\\s*$|wildtype", gt, ignore.case = TRUE) |
  grepl("\\bwt\\b|wildtype", d128$pheno$title, ignore.case = TRUE)
keep <- is_wt & grp %in% c("TBI", "Sham")

# ---------- GSE205958（与 26 脚本一致：FPKM xlsx，列序 TBI1-5, sham1-5） ----------
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

for (g in c("P2RY6", "HSD11B1", "CCR5")) {
  # GSE128543
  x_tbi <- ex128[g, keep & grp == "TBI"]
  x_sham <- ex128[g, keep & grp == "Sham"]
  sw_tbi <- tryCatch(shapiro.test(as.numeric(x_tbi))$p.value,
                     error = function(e) NA_real_)
  sw_sham <- tryCatch(shapiro.test(as.numeric(x_sham))$p.value,
                      error = function(e) NA_real_)
  vf <- tryCatch(var.test(as.numeric(x_tbi), as.numeric(x_sham))$p.value,
                 error = function(e) NA_real_)
  wt <- tryCatch(t.test(as.numeric(x_tbi), as.numeric(x_sham))$p.value,
                 error = function(e) NA_real_)
  add(data.frame(validation = "GSE128543_WT_only", gene = g,
                 n_tbi = length(x_tbi), n_sham = length(x_sham),
                 shapiro_tbi_p = sw_tbi, shapiro_sham_p = sw_sham,
                 var_test_p = vf, welch_t_p = wt))
  # GSE205958
  x_tbi2 <- as.numeric(ex205[g, 1:5])
  x_sham2 <- as.numeric(ex205[g, 6:10])
  sw_tbi2 <- tryCatch(shapiro.test(x_tbi2)$p.value, error = function(e) NA_real_)
  sw_sham2 <- tryCatch(shapiro.test(x_sham2)$p.value, error = function(e) NA_real_)
  vf2 <- tryCatch(var.test(x_tbi2, x_sham2)$p.value, error = function(e) NA_real_)
  wt2 <- tryCatch(t.test(x_tbi2, x_sham2)$p.value, error = function(e) NA_real_)
  add(data.frame(validation = "GSE205958", gene = g,
                 n_tbi = length(x_tbi2), n_sham = length(x_sham2),
                 shapiro_tbi_p = sw_tbi2, shapiro_sham_p = sw_sham2,
                 var_test_p = vf2, welch_t_p = wt2))
}

res <- bind_rows(checks)
write_csv(res, file.path(OUT, "assumption_checks.csv"))

md <- c(
  "# Welch t-test assumption checks (A-6)",
  paste("> generated:", Sys.time()),
  "",
  "Per-gene direction testing used Welch's t-test (two-sided, unequal variance).",
  "Assumption checks: Shapiro–Wilk normality per group; F-test (var.test) for",
  "variance homogeneity. Welch's t-test is robust to variance inequality and",
  "mild non-normality at these small sample sizes; P values in the manuscript",
  "are the Welch t-test P values. For n=5 per group, normality tests have",
  "limited power and are reported descriptively.",
  "",
  "```",
  capture.output(print(as.data.frame(res), row.names = FALSE)),
  "```"
)
writeLines(md, file.path(OUT, "ASSUMPTION_CHECKS.md"))
cat("assumption checks written\n")
print(as.data.frame(res), row.names = FALSE)
