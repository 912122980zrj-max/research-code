# ============================================================
# 02c_tier_dedup.R —— 从 union_targets_tiered_dedup.tsv 生成去别名候选表
# 运行：Rscript scripts/02c_tier_dedup.R
# 输出：results/05_Candidate_Intersection/candidate_genes_tiered_v2.csv（CLEAN）
#       results/05_Candidate_Intersection/candidate_genes_tiered_v2_FULL.csv（FULL）
# 说明：第二轮审查通过条件 1；03 脚本已内置相同逻辑，本脚本供当前结果直接更新。
# ============================================================

tier_order <- c(A = 0, B = 1, C = 2)

make_v2 <- function(cand_file, out_file) {
  if (!file.exists(cand_file)) {
    message("跳过：", cand_file)
    return(invisible(NULL))
  }
  cand <- read.csv(cand_file)$gene
  tier_d <- read.delim("data/04_BaP_targets/union_targets_tiered_dedup.tsv")
  tier_d$gene_symbol <- toupper(tier_d$gene_symbol)
  best <- tier_d[order(tier_order[tier_d$tier]), ]
  best <- best[!duplicated(best$gene_symbol), ]
  cand_up <- toupper(cand)
  m <- match(cand_up, best$gene_symbol)
  out <- data.frame(
    gene = cand,
    tier = ifelse(is.na(m), NA_character_, best$tier[m]),
    source = ifelse(is.na(m), NA_character_,
                    vapply(m, function(i) {
                      paste(unique(tier_d$source[tier_d$gene_symbol ==
                                                   best$gene_symbol[i]]),
                            collapse = ";")
                    }, character(1))),
    score_detail = ifelse(is.na(m), NA_character_,
                          paste(best$score[m], best$detail[m], sep = "; ")),
    stringsAsFactors = FALSE
  )
  write.csv(out, out_file, row.names = FALSE)
  message("写出：", out_file, "（", nrow(out), " 行；Tier A/B = ",
          sum(out$tier %in% c("A", "B"), na.rm = TRUE), "）")
}

make_v2("results/05_Candidate_Intersection/candidate_genes.csv",
        "results/05_Candidate_Intersection/candidate_genes_tiered_v2.csv")
make_v2("results/05_Candidate_Intersection/candidate_genes_FULL.csv",
        "results/05_Candidate_Intersection/candidate_genes_tiered_v2_FULL.csv")
