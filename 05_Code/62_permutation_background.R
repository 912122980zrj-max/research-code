# ============================================================
# 62_permutation_background.R —— B-2：双层置换背景检验
# 依据：修改报告_BaP_TBI_最终版.md 阶段 B-2
#   set.seed(2026)；主零假设：从 5,043 个可测基因中随机抽 K 重复 10,000 次；
#   表达广度分层抽样；black-only 路径；报期望值、富集倍数、经验 P、
#   95% 区间、观测百分位；DEG/模块定义不进零假设。
# 输出：results/16_Orthology/permutation_results.csv 及
#       permutation_distributions.rds、PERMUTATION_REPORT.txt
# ============================================================

suppressMessages({
  library(dplyr)
  library(readr)
})

ROOT <- getwd()
OUT <- file.path(ROOT, "results", "16_Orthology")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(OUT, "permutation_run.log"), append = TRUE)
}
if (file.exists(file.path(OUT, "permutation_run.log"))) {
  file.remove(file.path(OUT, "permutation_run.log"))
}

set.seed(2026)
N_PERM <- 10000

ex <- read.csv(file.path(ROOT, "results", "00_Data_Prep",
                        "corrected_expr.csv"),
               row.names = 1, check.names = FALSE)
universe <- rownames(ex)
N <- length(universe)

key_union <- read.csv(file.path(ROOT, "results", "04_WGCNA",
                                "tbi_key_genes.csv"))$gene
key_black <- read.csv(file.path(ROOT, "results", "04_WGCNA",
                                "D1_wgcna_sensitivity",
                                "best_module_key_genes.csv"))$gene

cand_union <- c("Ccr5", "Hsd11b1", "Lst1", "Mmp2", "Npy2r", "P2ry6", "Ripk1")
cand_black <- c("Hsd11b1", "Mmp2", "P2ry6")

# 表达广度（非零样本比例）
breadth <- rowMeans(ex > 0, na.rm = TRUE)
names(breadth) <- universe

decile <- cut(rank(breadth, ties.method = "first"),
              breaks = 10, labels = FALSE)
names(decile) <- universe

perm_one <- function(k_genes, key_set, stratified = FALSE,
                     cand = NULL, n_perm = N_PERM) {
  key_upper <- toupper(key_set)
  obs <- sum(toupper(k_genes) %in% key_upper)
  expected <- length(k_genes) * length(key_set) / N
  null <- integer(n_perm)
  if (!stratified) {
    for (i in seq_len(n_perm)) {
      null[i] <- sum(toupper(sample(universe, length(k_genes))) %in%
                       key_upper)
    }
  } else {
    cand_dec <- decile[toupper(names(decile)) %in% toupper(cand)]
    for (i in seq_len(n_perm)) {
      picked <- vapply(cand_dec, function(d) {
        pool <- universe[decile == d]
        if (length(pool) < 2) pool <- universe
        sample(pool, 1)
      }, character(1))
      null[i] <- sum(toupper(picked) %in% key_upper)
    }
  }
  emp_p <- (1 + sum(null >= obs)) / (n_perm + 1)
  ci <- quantile(null, probs = c(0.025, 0.975), names = TRUE)
  percentile <- mean(null <= obs)
  list(obs = obs, expected = expected,
       enrichment = if (expected > 0) obs / expected else NA_real_,
       emp_p = emp_p, ci_lo = unname(ci[1]), ci_hi = unname(ci[2]),
       percentile = percentile, null = null)
}

run <- list(
  union_unstratified = perm_one(cand_union, key_union, FALSE, cand_union),
  union_stratified   = perm_one(cand_union, key_union, TRUE, cand_union),
  black_unstratified = perm_one(cand_black, key_black, FALSE, cand_black),
  black_stratified   = perm_one(cand_black, key_black, TRUE, cand_black)
)

res <- lapply(names(run), function(nm) {
  r <- run[[nm]]
  data.frame(path = nm, n_perm = N_PERM, universe = N,
             key_genes = if (grepl("black", nm)) length(key_black)
               else length(key_union),
             k_drawn = if (grepl("black", nm)) length(cand_black)
               else length(cand_union),
             observed = r$obs, expected = r$expected,
             enrichment_fold = r$enrichment,
             empirical_p = r$emp_p,
             ci95_lo = r$ci_lo, ci95_hi = r$ci_hi,
             observed_percentile = r$percentile,
             stringsAsFactors = FALSE)
}) %>% bind_rows()

write_csv(res, file.path(OUT, "permutation_results.csv"))
saveRDS(lapply(run, function(r) r$null),
        file.path(OUT, "permutation_distributions.rds"))

# 超几何参考（单侧，观测数 = 交集）
hg <- function(k, key_len, n) {
  phyper(k - 1, key_len, N - key_len, n, lower.tail = FALSE)
}
res$hypergeom_p <- c(
  hg(run$union_unstratified$obs, length(key_union), length(cand_union)),
  hg(run$union_stratified$obs, length(key_union), length(cand_union)),
  hg(run$black_unstratified$obs, length(key_black), length(cand_black)),
  hg(run$black_stratified$obs, length(key_black), length(cand_black))
)
write_csv(res, file.path(OUT, "permutation_results.csv"))

sink(file.path(OUT, "PERMUTATION_REPORT.txt"))
cat("B-2 double-layer permutation background test\n")
cat("generated:", format(Sys.time()), "\n")
cat("seed: 2026; n_perm:", N_PERM, "; universe: 5,043 measurable genes\n")
cat("key genes (union modules):", length(key_union),
    "; key genes (black-only):", length(key_black), "\n")
cat("observed candidates (union): 7; (black-only): 3\n\n")
print(as.data.frame(res), row.names = FALSE)
sink()

log_msg("B-2 done")
print(res)
