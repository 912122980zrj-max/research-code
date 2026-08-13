# ============================================================
# 08_d1_wgcna_sensitivity.R —— D1 WGCNA 模块口径敏感性分析
# 依据：REVIEW4_FINAL_AUDIT D1（与模板"单一最佳模块"口径对齐）
# 输入：results/00_Data_Prep/corrected_expr.csv、results/04_WGCNA/module_colors.csv、
#       deg_significant.csv、module_genes.csv、union_targets_tiered_dedup.tsv
# 输出：results/04_WGCNA/D1_*（新文件，不覆盖主结果）
# ============================================================

suppressMessages({library(dplyr); library(readr)})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "results", "04_WGCNA", "D1_wgcna_sensitivity")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
report <- file.path(OUT, "D1_report.txt")
write_report <- function(...) cat(..., "\n", sep = "", file = report,
                                  append = TRUE)
writeLines("D1 WGCNA module-criteria sensitivity", report)
write_report("generated:", format(Sys.time()))

deg <- read.csv(file.path(ROOT, "results", "03_DEG_Analysis",
                          "deg_significant.csv"))
deg_genes <- deg$gene
write_report("DEGs:", length(deg_genes))

mc <- read.csv(file.path(ROOT, "results", "04_WGCNA", "module_colors.csv"))
mu <- read.csv(file.path(ROOT, "results", "04_WGCNA", "module_genes.csv"))
mt <- read.csv(file.path(ROOT, "results", "04_WGCNA", "module_trait.csv"))
best_mod <- sub("^ME", "", mt$module[which.max(mt$r)])
write_report("best module (max |r|):", best_mod)
magenta_genes <- mc$gene[mc$module == "magenta"]
best_mod_genes <- mc$gene[mc$module == best_mod]
union_genes <- mu$gene
write_report("magenta module genes:", length(magenta_genes))
write_report(paste0(best_mod, " module genes:"), length(best_mod_genes))
write_report("13-module union genes:", length(union_genes))

bap <- read.delim(file.path(ROOT, "data", "04_BaP_targets",
                            "union_targets_tiered_dedup.tsv"))
bap_genes <- unique(toupper(bap$gene_symbol))

key_magenta <- intersect(deg_genes, magenta_genes)
key_best <- intersect(deg_genes, best_mod_genes)
key_union <- intersect(deg_genes, union_genes)
cand_magenta <- key_magenta[toupper(key_magenta) %in% bap_genes]
cand_best <- key_best[toupper(key_best) %in% bap_genes]
cand_union <- key_union[toupper(key_union) %in% bap_genes]
write_report("key genes (magenta-only):", length(key_magenta))
write_report(paste0("key genes (", best_mod, "-only):"), length(key_best))
write_report("key genes (13-module union):", length(key_union))
write_report("BaP-TBI candidates (magenta-only):", length(cand_magenta))
write_report(paste0("BaP-TBI candidates (", best_mod, "-only):"),
             length(cand_best))
write_report("BaP-TBI candidates (13-module union):", length(cand_union))

write.csv(data.frame(gene = key_magenta), file.path(OUT, "magenta_key_genes.csv"),
          row.names = FALSE)
write.csv(data.frame(gene = cand_magenta), file.path(OUT, "magenta_candidates.csv"),
          row.names = FALSE)
write.csv(data.frame(gene = key_best), file.path(OUT, "best_module_key_genes.csv"),
          row.names = FALSE)
write.csv(data.frame(gene = cand_best), file.path(OUT, "best_module_candidates.csv"),
          row.names = FALSE)
write.csv(data.frame(gene = cand_union), file.path(OUT, "union_candidates.csv"),
          row.names = FALSE)

# 3 特征签名在两种口径中的保留情况（v4 重跑后：P2RY6/HSD11B1/CCR5）
sig3 <- c("P2ry6", "Hsd11b1", "Ccr5")
sig3_up <- toupper(sig3)
hit_magenta <- sig3_up %in% toupper(cand_magenta)
hit_best <- sig3_up %in% toupper(cand_best)
hit_union <- sig3_up %in% toupper(cand_union)
write_report("3-gene signature in magenta-only:",
             paste(sig3[hit_magenta], collapse = ","),
             "| missing:", paste(sig3[!hit_magenta], collapse = ","))
write_report(paste0("3-gene signature in ", best_mod, "-only:"),
             paste(sig3[hit_best], collapse = ","),
             "| missing:", paste(sig3[!hit_best], collapse = ","))
write_report("3-gene signature in union:",
             paste(sig3[hit_union], collapse = ","),
             "| missing:", paste(sig3[!hit_union], collapse = ","))

# 两口径候选差集
write_report(paste0("candidates in union but not ", best_mod, ":"),
             paste(setdiff(cand_union, cand_best), collapse = ","))
write_report(paste0("candidates in ", best_mod, " but not union:"),
             paste(setdiff(cand_best, cand_union), collapse = ","))
write_report("note: magenta 在 log2 重跑后不再为 TBI 显著模块（0 基因），",
             "magenta-only 口径仅作历史对照。")

message("D1 完成。结果在 ", OUT)
