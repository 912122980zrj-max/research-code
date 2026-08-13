# ============================================================
# 61_orthology_mapping.R —— B-1：Ensembl Compara 直系同源映射
# 依据：修改报告_BaP_TBI_最终版.md 阶段 B-1
#   输出全部 364 行映射表 + Ensembl release；1:1 进主分析、1:多单列报告；
#   199/165 账目、摘要百分比、超几何 n 随结果更新；Hprt 探针检测情况说明。
# 输出：
#   results/16_Orthology/
#     bap_targets_364_tier.csv         （364 行：human symbol + final tier）
#     orthology_mapping_full.csv       （364 行 × 映射明细，含 orthology type）
#     orthology_summary.txt
#     evaluable_targets_CLEAN.csv / evaluable_targets_FULL.csv
#     candidates_orthology_CLEAN.csv / candidates_orthology_FULL.csv
# ============================================================

suppressMessages({
  library(biomaRt)
  library(dplyr)
  library(readr)
})

ROOT <- getwd()
OUT <- file.path(ROOT, "results", "16_Orthology")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(OUT, "orthology_run.log"), append = TRUE)
}
if (file.exists(file.path(OUT, "orthology_run.log"))) {
  file.remove(file.path(OUT, "orthology_run.log"))
}

log_msg("B-1 orthology mapping started")

# ---------- 1) 重建 364 基因 final tier 表 ----------
rows <- read_tsv(file.path(ROOT, "data", "04_BaP_targets",
                          "union_targets_tiered_dedup.tsv"),
                 col_types = cols(.default = col_character()))
rank_tier <- c(A = 3, B = 2, C = 1)
best_tier <- tibble()
for (g in unique(rows$official_symbol)) {
  sub <- rows[rows$official_symbol == g, ]
  rt <- rank_tier[sub$tier]
  rt[is.na(rt)] <- 0
  ti <- names(rt)[which.max(rt)]
  if (max(rt) == 0) ti <- "C"
  best_tier <- bind_rows(best_tier, tibble(
    gene_symbol = g,
    tier = ti,
    sources = paste(unique(sub$source), collapse = ";")
  ))
}
stopifnot(nrow(best_tier) == 364)
stopifnot(all(sort(table(best_tier$tier)) ==
              sort(c(A = 5, B = 23, C = 336))))
write_csv(best_tier, file.path(OUT, "bap_targets_364_tier.csv"))
log_msg("364-gene tier table rebuilt: A=5 B=23 C=336")

# ---------- 2) Ensembl 版本记录 ----------
ens_version <- tryCatch({
  con <- url("https://rest.ensembl.org/info/software/version/?content-type=application/json",
             open = "rb")
  on.exit(close(con))
  rawToChar(readBin(con, "raw", n = 200))
}, error = function(e) NA_character_)
log_msg("Ensembl REST software version:", ens_version)

# ---------- 3) 直系同源查询（优先 REST 抓取表，回退 biomaRt） ----------
RAW_TSV <- file.path(OUT, "orthology_raw_ensembl.tsv")
if (file.exists(RAW_TSV)) {
  map_raw <- read_tsv(RAW_TSV, col_types = cols(.default = col_character()))
  map_raw <- map_raw %>%
    rename(human_symbol = human_symbol,
           human_ensembl = human_ensembl,
           mouse_ensembl = mouse_ensembl,
           mouse_symbol = mouse_symbol,
           orthology_type = orthology_type)
  log_msg("REST orthology table loaded:", nrow(map_raw), "rows")
} else {
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  log_msg("Human mart host:", human@host)
  attrs <- c("hgnc_symbol", "ensembl_gene_id",
             "mmusculus_homolog_ensembl_gene",
             "mmusculus_homolog_associated_gene_name",
             "mmusculus_homolog_orthology_type",
             "mmusculus_homolog_orthology_confidence")
  map_raw <- getBM(attributes = attrs,
                   filters = "hgnc_symbol",
                   values = best_tier$gene_symbol,
                   mart = human)
  log_msg("biomaRt rows returned:", nrow(map_raw))
  map_raw <- map_raw %>%
    rename(human_symbol = hgnc_symbol,
           human_ensembl = ensembl_gene_id,
           mouse_ensembl = mmusculus_homolog_ensembl_gene,
           mouse_symbol = mmusculus_homolog_associated_gene_name,
           orthology_type = mmusculus_homolog_orthology_type)
}

# ---------- 4) 1:1 / 1:多 分类 ----------
if ("hgnc_symbol" %in% colnames(map_raw)) {
  map <- map_raw %>%
    rename(human_symbol = hgnc_symbol,
           human_ensembl = ensembl_gene_id,
           mouse_ensembl = mmusculus_homolog_ensembl_gene,
           mouse_symbol = mmusculus_homolog_associated_gene_name,
           orthology_type = mmusculus_homolog_orthology_type) %>%
    filter(!is.na(mouse_ensembl) & mouse_ensembl != "" &
             !is.na(mouse_symbol) & mouse_symbol != "")
} else {
  map <- map_raw %>%
    filter(!is.na(mouse_ensembl) & mouse_ensembl != "" &
             !is.na(mouse_symbol) & mouse_symbol != "")
}

map <- map %>%
  group_by(human_symbol) %>%
  mutate(orthology_class = ifelse(
    all(orthology_type == "ortholog_one2one"),
    "one2one", "one2many")) %>%
  ungroup()

# 每个 human gene 是否恰有一个 mouse ortholog
n_mouse <- map %>% count(human_symbol, name = "n_mouse_orthologs")
map <- map %>% left_join(n_mouse, by = "human_symbol")
map$orthology_class[map$n_mouse_orthologs != 1] <- "one2many"

# ---------- 5) 合并 tier 并输出全表 ----------
full <- best_tier %>%
  left_join(map, by = c("gene_symbol" = "human_symbol")) %>%
  mutate(mapped = !is.na(mouse_symbol),
         orthology_class = ifelse(is.na(orthology_class), "unmapped",
                                  orthology_class))

# ---------- 6) 与小鼠矩阵比对（CLEAN / FULL） ----------
ex_clean <- read.csv(file.path(ROOT, "results", "00_Data_Prep",
                               "corrected_expr.csv"),
                     row.names = 1, check.names = FALSE)
ex_full <- read.csv(file.path(ROOT, "results", "00_Data_Prep",
                              "corrected_expr_FULL.csv"),
                    row.names = 1, check.names = FALSE)
rn_clean <- toupper(rownames(ex_clean))
rn_full <- toupper(rownames(ex_full))

# 经证实的 Ensembl 显示名 -> GEO 矩阵符号别名
# HPRT1 的 1:1 直系同源为 ENSMUSG00000025630；Ensembl 显示名 "Hprt1"，
# GEO 矩阵注释为 "Hprt"（CLEAN 33/33 非零），两者为同一基因。
alias_map <- c(Hprt1 = "Hprt")
full$matrix_symbol <- full$mouse_symbol
full$alias_note <- ""
for (i in seq_len(nrow(full))) {
  if (!is.na(full$mouse_symbol[i]) &&
      toupper(full$mouse_symbol[i]) %in% toupper(names(alias_map))) {
    hit <- names(alias_map)[toupper(names(alias_map)) ==
                              toupper(full$mouse_symbol[i])][1]
    full$matrix_symbol[i] <- unname(alias_map[[hit]])
    full$alias_note[i] <- paste0(
      "Ensembl display name ", full$mouse_symbol[i],
      " for ", full$mouse_ensembl[i],
      "; GEO matrix symbol ", full$matrix_symbol[i],
      " (verified in corrected_expr.csv)")
  }
}

full$in_matrix_CLEAN <- toupper(full$matrix_symbol) %in% rn_clean
full$in_matrix_FULL  <- toupper(full$matrix_symbol) %in% rn_full

eval_clean <- full %>% filter(mapped & orthology_class == "one2one" &
                                in_matrix_CLEAN)
eval_full  <- full %>% filter(mapped & orthology_class == "one2one" &
                                in_matrix_FULL)

write_csv(full, file.path(OUT, "orthology_mapping_full.csv"))
write_csv(eval_clean, file.path(OUT, "evaluable_targets_CLEAN.csv"))
write_csv(eval_full, file.path(OUT, "evaluable_targets_FULL.csv"))

# ---------- 7) 候选交集（沿用现有 key genes） ----------
key_clean <- read.csv(file.path(ROOT, "results", "04_WGCNA",
                                "tbi_key_genes.csv"))$gene
key_full  <- read.csv(file.path(ROOT, "results", "04_WGCNA",
                                "tbi_key_genes_FULL.csv"))$gene
key_black <- read.csv(file.path(ROOT, "results", "04_WGCNA",
                                "D1_wgcna_sensitivity",
                                "best_module_key_genes.csv"))$gene

cand_clean <- eval_clean %>%
  filter(toupper(mouse_symbol) %in% toupper(key_clean)) %>%
  arrange(tier, gene_symbol)
cand_full  <- eval_full %>%
  filter(toupper(mouse_symbol) %in% toupper(key_full)) %>%
  arrange(tier, gene_symbol)
cand_black <- eval_clean %>%
  filter(toupper(mouse_symbol) %in% toupper(key_black)) %>%
  arrange(tier, gene_symbol)

write_csv(cand_clean, file.path(OUT, "candidates_orthology_CLEAN.csv"))
write_csv(cand_full, file.path(OUT, "candidates_orthology_FULL.csv"))
write_csv(cand_black, file.path(OUT, "candidates_orthology_black_only.csv"))

# ---------- 8) 汇总（独特基因计数） ----------
uniq_mapped <- unique(full$gene_symbol[full$mapped])
uniq_one2one <- unique(full$gene_symbol[full$mapped &
                                          full$orthology_class == "one2one"])
uniq_one2many <- unique(full$gene_symbol[full$mapped &
                                           full$orthology_class == "one2many"])
summ <- list(
  ensembl_version = ens_version,
  human_genes = nrow(best_tier),
  mapped_genes = length(uniq_mapped),
  unmapped_genes = nrow(best_tier) - length(uniq_mapped),
  one2one_genes = length(uniq_one2one),
  one2many_genes = length(uniq_one2many),
  evaluable_CLEAN = nrow(eval_clean),
  evaluable_FULL = nrow(eval_full),
  tier_evaluable_CLEAN = paste(names(table(eval_clean$tier)),
                               as.integer(table(eval_clean$tier)),
                               collapse = "; "),
  tier_evaluable_FULL = paste(names(table(eval_full$tier)),
                              as.integer(table(eval_full$tier)),
                              collapse = "; "),
  candidates_CLEAN = paste(cand_clean$mouse_symbol, collapse = ","),
  candidates_FULL = paste(cand_full$mouse_symbol, collapse = ","),
  candidates_black_only = paste(cand_black$mouse_symbol, collapse = ",")
)

# 汇总内嵌 version（REST 抓取表自带 release）
if (file.exists(file.path(OUT, "orthology_raw_ensembl.tsv"))) {
  rel <- unique(read_tsv(file.path(OUT, "orthology_raw_ensembl.tsv"),
                         col_types = cols(.default = col_character()))$ensembl_release)
  summ$ensembl_version <- paste(rel, collapse = ";")
}

sink(file.path(OUT, "orthology_summary.txt"))
cat("B-1 orthology mapping summary\n")
cat("generated:", format(Sys.time()), "\n")
for (nm in names(summ)) cat(nm, "=", summ[[nm]], "\n")
cat("\nTier A detail:\n")
ta <- full %>% filter(tier == "A")
print(ta %>% select(gene_symbol, tier, sources, mouse_symbol,
                    orthology_type, orthology_class,
                    in_matrix_CLEAN, in_matrix_FULL),
      row.names = FALSE)
cat("\nHprt probe/expression note:\n")
hprt_row <- which(toupper(rownames(ex_clean)) == "HPRT")
if (length(hprt_row)) {
  x <- as.numeric(ex_clean[hprt_row[1], ])
  cat("Hprt in CLEAN matrix: yes; n=", length(x),
      "; non-zero=", sum(x > 0, na.rm = TRUE),
      "; range=", paste(round(range(x, na.rm = TRUE), 3), collapse = ".."),
      "\n", sep = "")
} else {
  cat("Hprt in CLEAN matrix: no\n")
}
sink()

log_msg("B-1 done: mapped=", summ$mapped_genes,
        " one2one=", summ$one2one_genes,
        " one2many=", summ$one2many_genes,
        " evaluable_CLEAN=", summ$evaluable_CLEAN,
        " evaluable_FULL=", summ$evaluable_FULL,
        " candidates_CLEAN=", summ$candidates_CLEAN)
