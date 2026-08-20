# ============================================================
# 71_permutation_v7.R —— v7 主分析置换背景检验（GSE58485-only）
# 依据：修改报告 B-2 设计；主零假设从 5,043 个可测基因中随机抽 K
#       重复 10,000 次（seed=2026）；表达广度分层变体；
#       union 路径（K=11, key=364）与 best-module（green）路径（K=3, key=109）。
# 输出：results/16_Orthology/v7_single_cohort/permutation_v7.csv / .rds / REPORT
# ============================================================
suppressMessages({
  library(dplyr)
  library(readr)
})

ROOT <- getwd()
OUT <- file.path(ROOT, "results", "16_Orthology", "v7_single_cohort")
set.seed(2026)
N_PERM <- 10000

# ---------- 5,043 基因宇宙 + GSE58485-only 表达广度 ----------
comp <- read.csv(file.path(ROOT, "results", "00_Data_Prep",
                           "sample_composition.csv"), check.names = FALSE)
clean_samples <- comp$sample[comp$include_clean]
read_bulk <- function(gse) readRDS(file.path(ROOT, "data", "01_GEO_bulk",
                                             paste0(gse, ".rds")))
gses <- c("GSE58485", "GSE41345")
dat <- lapply(gses, read_bulk)
names(dat) <- gses
common <- Reduce(intersect, lapply(dat, function(d) rownames(d$expr)))
expr_all <- do.call(cbind, lapply(gses, function(g) dat[[g]]$expr[common, , drop = FALSE]))
expr_all <- expr_all[, intersect(clean_samples, colnames(expr_all)), drop = FALSE]
cohort_vec <- comp$cohort[match(colnames(expr_all), comp$sample)]
keep_gene <- rowSums(is.na(expr_all)) < 0.2 * ncol(expr_all)
expr_all <- expr_all[keep_gene, , drop = FALSE]
cov_a <- rowMeans(expr_all[, cohort_vec == "GSE58485", drop = FALSE] > 0, na.rm = TRUE)
cov_b <- rowMeans(expr_all[, cohort_vec == "GSE41345", drop = FALSE] > 0, na.rm = TRUE)
keep_cov <- cov_a >= 0.10 & cov_b >= 0.10
expr58 <- expr_all[keep_cov, cohort_vec == "GSE58485", drop = FALSE]
universe <- rownames(expr58)
N <- length(universe)
stopifnot(N == 5043)

breadth <- rowMeans(expr58 > 0, na.rm = TRUE)
names(breadth) <- universe
decile <- cut(rank(breadth, ties.method = "first"), breaks = 10, labels = FALSE)
names(decile) <- universe

# ---------- 固定 key / candidate 集（DEG 与模块定义不进零假设） ----------
key_union <- read.csv(file.path(OUT, "key_genes_GSE58485_only.csv"))$gene
key_best  <- read.csv(file.path(OUT, "key_genes_best_module_GSE58485_only.csv"))$gene
cand_union <- read.csv(file.path(OUT, "candidates_GSE58485_only.csv"))$gene
cand_best  <- read.csv(file.path(OUT, "candidates_best_module_GSE58485_only.csv"))$gene
stopifnot(length(key_union) == 364, length(key_best) == 109,
          length(cand_union) == 11, length(cand_best) == 3)

perm_one <- function(k_genes, key_set, stratified = FALSE,
                     cand = NULL, n_perm = N_PERM) {
  key_upper <- toupper(key_set)
  obs <- sum(toupper(k_genes) %in% key_upper)
  expected <- length(k_genes) * length(key_set) / N
  null <- integer(n_perm)
  if (!stratified) {
    for (i in seq_len(n_perm)) {
      null[i] <- sum(toupper(sample(universe, length(k_genes))) %in% key_upper)
    }
  } else {
    cand_dec <- decile[toupper(names(decile)) %in% toupper(cand)]
    for (i in seq_len(n_perm)) {
      picked <- vapply(cand_dec, function(d) {
        pool <- universe[decile == d]
        sample(pool, 1)
      }, character(1))
      null[i] <- sum(toupper(picked) %in% key_upper)
    }
  }
  emp_p <- (1 + sum(null >= obs)) / (n_perm + 1)
  ci <- quantile(null, c(0.025, 0.975))
  hp <- phyper(obs - 1, length(key_set), N - length(key_set),
               length(k_genes), lower.tail = FALSE)
  data.frame(n_perm = n_perm, universe = N,
             key_genes = length(key_set), k_drawn = length(k_genes),
             observed = obs, expected = expected,
             enrichment_fold = obs / expected,
             empirical_p = emp_p,
             ci95_lo = unname(ci[1]), ci95_hi = unname(ci[2]),
             observed_percentile = mean(null <= obs),
             hypergeom_p = hp)
}

res <- bind_rows(
  perm_one(cand_union, key_union, FALSE) %>% mutate(path = "union_unstratified"),
  perm_one(cand_union, key_union, TRUE, cand_union) %>% mutate(path = "union_stratified"),
  perm_one(cand_best, key_best, FALSE) %>% mutate(path = "best_module_unstratified"),
  perm_one(cand_best, key_best, TRUE, cand_best) %>% mutate(path = "best_module_stratified")
) %>% select(path, everything())

write_csv(res, file.path(OUT, "permutation_v7.csv"))
saveRDS(list(results = res, universe = universe, breadth = breadth),
        file.path(OUT, "permutation_v7.rds"))

writeLines(c(
  "v7 permutation background test (GSE58485-only main path)",
  paste0("generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("seed=2026; draws=", N_PERM, "; universe=", N, " genes"),
  capture.output(print(as.data.frame(res), digits = 6))
), file.path(OUT, "PERMUTATION_V7_REPORT.txt"))
print(as.data.frame(res), digits = 6)
