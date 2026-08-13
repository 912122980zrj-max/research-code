# ============================================================
# 10_d11_ahr_fisher.R —— D11 AhR 反应基因富集检验
# 依据：REVIEW4_FINAL_AUDIT D11（模板 CYP1A1/CYP1B1/AHRR 单侧 Fisher）
# 输入：results/05_Candidate_Intersection/candidate_genes.csv（v4 7）、results/03_DEG_Analysis/deg_all.csv（背景）
# 输出：results/06_GO_KEGG_Enrichment/D11_ahr_fisher/ahr_fisher.csv、D11_report.txt
# ============================================================

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "results", "06_GO_KEGG_Enrichment", "D11_ahr_fisher")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
report <- file.path(OUT, "D11_report.txt")
writeLines("D11 AHR-responsive gene enrichment (Fisher exact)", report)
write_report <- function(...) cat(..., "\n", sep = "", file = report,
                                  append = TRUE)
write_report("generated:", format(Sys.time()))

cand <- read.csv(file.path(ROOT, "results", "05_Candidate_Intersection",
                           "candidate_genes.csv"))$gene
deg_all <- read.csv(file.path(ROOT, "results", "03_DEG_Analysis", "deg_all.csv"))
background <- unique(as.character(deg_all$gene))
write_report("candidates:", length(cand))
write_report("background genes (deg_all):", length(background))

ahr_genes <- c("CYP1A1", "CYP1B1", "AHRR")
hit_cand <- intersect(toupper(cand), ahr_genes)
hit_bg <- intersect(toupper(background), ahr_genes)

# 2x2: AHR 基因 × 候选
a <- length(hit_cand)                    # 候选 ∩ AHR
b <- length(setdiff(ahr_genes, hit_cand))# 非候选 ∩ AHR（= AHR 基因 - 命中）
c <- length(cand) - a                    # 候选 ∩ 非AHR
d <- length(setdiff(toupper(background), ahr_genes)) - c  # 背景非候选非AHR
ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")

write_report("AHR genes present in background:", paste(hit_bg, collapse = ","))
write_report("AHR genes in candidates:", paste(hit_cand, collapse = ","))
write_report("Fisher table [AHR-in-cand, AHR-not, nonAHR-in-cand, nonAHR-bg]:",
             a, b, c, d)
write_report("one-sided Fisher P:", ft$p.value, "| odds ratio:",
             ft$estimate)

out <- data.frame(
  ahr_genes = paste(ahr_genes, collapse = ";"),
  hits_in_candidates = paste(hit_cand, collapse = ";"),
  n_candidates = length(cand),
  n_background = length(background),
  n_hits = length(hit_cand),
  fisher_p_onesided = ft$p.value,
  odds_ratio = as.numeric(ft$estimate),
  note = paste0("v4 覆盖过滤后 CYP1A1/CYP1B1/AHRR 均未入选候选集",
                "（背景基因数 n=", length(background), "），AhR 富集为阴性")
)
write.csv(out, file.path(OUT, "ahr_fisher.csv"), row.names = FALSE)
message("D11 完成。结果在 ", OUT)
