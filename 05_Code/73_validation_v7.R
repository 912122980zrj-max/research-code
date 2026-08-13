# ============================================================
# 73_validation_v7.R —— v7 主分析外部验证（11 候选）
# GSE128543（WT-only，n=14）与 GSE205958（n=10）：
#   逐基因 Welch t 检验方向一致性 + Spearman 秩相关；
# 人类 GSE104687：11 个人类同源符号覆盖检查。
# 输出：results/08_External_Validation/v7_single_cohort/
# ============================================================
suppressMessages({
  library(dplyr)
  library(readr)
  library(readxl)
})

ROOT <- getwd()
OUT <- file.path(ROOT, "results", "08_External_Validation", "v7_single_cohort")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cand_stat <- read_csv(file.path(ROOT, "results", "16_Orthology", "v7_single_cohort",
                                "candidate_stats_GSE58485_only.csv"),
                      col_types = cols(.default = col_character()))
genes <- toupper(cand_stat$gene)
disc_lfc <- setNames(as.numeric(cand_stat$logFC), genes)
disc_sign <- setNames(sign(disc_lfc), genes)

read_bulk <- function(gse) readRDS(file.path(ROOT, "data", "01_GEO_bulk",
                                             paste0(gse, ".rds")))
wt_test <- function(x_tbi, x_sham) {
  tryCatch(t.test(x_tbi, x_sham)$p.value, error = function(e) NA_real_)
}

# ---------- GSE128543 WT-only ----------
d128 <- read_bulk("GSE128543")
ex128 <- d128$expr
rn <- toupper(rownames(ex128))
ex128 <- ex128[rn %in% genes, , drop = FALSE]
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
cat("GSE128543 WT samples kept:", sum(keep), "\n")

dir128 <- do.call(rbind, lapply(genes, function(g) {
  if (!g %in% rownames(ex128)) {
    return(data.frame(validation = "GSE128543_WT_only", gene = g,
                      diff = NA, t_p = NA, direction_consistent = NA,
                      stringsAsFactors = FALSE))
  }
  x_tbi <- ex128[g, keep & grp == "TBI"]
  x_sham <- ex128[g, keep & grp == "Sham"]
  diff <- mean(x_tbi) - mean(x_sham)
  data.frame(validation = "GSE128543_WT_only", gene = g,
             diff = round(diff, 4),
             t_p = round(wt_test(x_tbi, x_sham), 6),
             direction_consistent = sign(diff) == disc_sign[[g]],
             stringsAsFactors = FALSE)
}))

# ---------- GSE205958 ----------
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
  ex205 <- rowsum(ex205, group = rownames(ex205))
}

dir205 <- do.call(rbind, lapply(genes, function(g) {
  if (!g %in% rownames(ex205)) {
    return(data.frame(validation = "GSE205958", gene = g,
                      diff = NA, t_p = NA, direction_consistent = NA,
                      stringsAsFactors = FALSE))
  }
  x_tbi <- ex205[g, 1:5]
  x_sham <- ex205[g, 6:10]
  diff <- mean(x_tbi) - mean(x_sham)
  data.frame(validation = "GSE205958", gene = g,
             diff = round(diff, 4),
             t_p = round(wt_test(x_tbi, x_sham), 6),
             direction_consistent = sign(diff) == disc_sign[[g]],
             stringsAsFactors = FALSE)
}))

val <- bind_rows(dir128, dir205) %>%
  left_join(cand_stat %>% select(gene, logFC, adjP) %>%
              mutate(gene = toupper(gene)), by = "gene")
write_csv(val, file.path(OUT, "validation_direction_v7.csv"))

cat("direction consistency:\n")
print(val %>% group_by(validation) %>%
        summarise(present = sum(!is.na(diff)),
                  consistent = sum(direction_consistent %in% TRUE,
                                   na.rm = TRUE),
                  .groups = "drop"))

# Spearman：发现 log2FC 与验证 diff（两队列合并，仅共有基因）
for (v in c("GSE128543_WT_only", "GSE205958")) {
  sub <- val %>% filter(validation == v, !is.na(diff))
  if (nrow(sub) >= 4) {
    sp <- cor.test(disc_lfc[sub$gene], sub$diff, method = "spearman")
    cat(v, "genes:", nrow(sub), "spearman rho =", round(sp$estimate, 3),
        "P =", format(sp$p.value, digits = 3), "\n")
  } else {
    cat(v, "genes:", nrow(sub), "too few for rank correlation\n")
  }
}

# ---------- 人类 GSE104687 覆盖 ----------
hsa <- c("ADAM17", "CCR5", "CHRM3", "HSD11B1", "LST1", "MMP2",
         "NPY1R", "NPY2R", "NQO2", "P2RY6", "RIPK1")
d104 <- readRDS(file.path(ROOT, "data", "01_GEO_bulk", "GSE104687.rds"))
rn_h <- toupper(rownames(d104$expr))
cov_h <- data.frame(gene = hsa,
                    in_human = hsa %in% rn_h)
write_csv(cov_h, file.path(OUT, "human_coverage_v7.csv"))
print(cov_h)
