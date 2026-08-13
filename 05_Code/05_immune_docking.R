# ============================================================
# 05_immune_docking.R —— 免疫浸润（mini-CIBERSORT）+ 相关分析
#                        + 分子对接结构准备（AlphaFold/PDB）
# 运行：Rscript scripts/05_immune_docking.R
# 输入：results/00_Data_Prep/corrected_expr.csv + corrected_pheno.csv
#       results/07_ML_Signature/enet_coefficients.csv（hub 基因）
#       data/03_LM22/LM22.txt（CIBERSORT 官网注册后下载，可选）
# 输出：results/05_immune_docking/
# ============================================================

suppressMessages({
  library(e1071)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(jsonlite)
})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "results", "05_immune_docking")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---------- mini-CIBERSORT（nu-SVR；与 CIBERSORT 原文同思路） ----------
run_cibersort <- function(expr, sig, perm = 0, QN = TRUE) {
  common <- intersect(rownames(expr), rownames(sig))
  if (length(common) < 100) {
    warning("签名矩阵与表达矩阵共有基因过少: ", length(common))
    return(NULL)
  }
  E <- expr[common, , drop = FALSE]
  S <- as.matrix(sig[common, , drop = FALSE])
  if (QN && requireNamespace("preprocessCore", quietly = TRUE)) {
    E <- preprocessCore::normalize.quantiles(E)
    S <- preprocessCore::normalize.quantiles(S)
  }
  frac <- matrix(0, nrow = ncol(S), ncol = ncol(E),
                 dimnames = list(colnames(S), colnames(E)))
  pval <- rep(1, ncol(E)); names(pval) <- colnames(E)
  for (j in seq_len(ncol(E))) {
    y <- E[, j]
    sv <- tryCatch(svm(x = S, y = y, type = "nu-regression", kernel = "linear",
                       nu = 0.5, cost = 10, scale = FALSE, cachesize = 500),
                   error = function(e) NULL)
    if (is.null(sv)) next
    w <- as.numeric(t(sv$coefs) %*% sv$SV)
    w <- pmax(w, 0)
    if (sum(w) > 0) frac[, j] <- w / sum(w)
    if (perm > 0) {
      fit <- S %*% w
      r0 <- suppressWarnings(cor(fit, y))
      rp <- sapply(seq_len(perm), function(k) {
        yp <- sample(y)
        svp <- tryCatch(svm(x = S, y = yp, type = "nu-regression",
                            kernel = "linear", nu = 0.5, cost = 10,
                            scale = FALSE, cachesize = 500), error = function(e) NULL)
        if (is.null(svp)) return(NA_real_)
        wp <- as.numeric(t(svp$coefs) %*% svp$SV)
        suppressWarnings(cor(S %*% pmax(wp, 0), yp))
      })
      pval[j] <- (1 + sum(rp >= r0, na.rm = TRUE)) / (1 + perm)
    }
  }
  list(fractions = frac, pvalue = pval)
}

# ---------- 1) 免疫浸润（LM22 存在时运行） ----------
lm22_file <- file.path(ROOT, "data", "03_LM22", "LM22.txt")
expr_file <- file.path(ROOT, "results", "00_Data_Prep", "corrected_expr.csv")
pheno_file <- file.path(ROOT, "results", "00_Data_Prep", "corrected_pheno.csv")

if (file.exists(lm22_file) && file.exists(expr_file)) {
  lm22 <- as.matrix(read.delim(lm22_file, row.names = 1, check.names = FALSE))
  expr <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE))
  rownames(expr) <- toupper(rownames(expr))
  rownames(lm22) <- toupper(rownames(lm22))
  res <- run_cibersort(expr, lm22, perm = 0)
  if (!is.null(res)) {
    write.csv(res$fractions, file.path(OUT, "cibersort_fractions.csv"))
    write.csv(data.frame(sample = names(res$pvalue), pvalue = res$pvalue),
              file.path(OUT, "cibersort_pvalue.csv"), row.names = FALSE)

    # 组间差异（Wilcoxon）
    ph <- read.csv(pheno_file, check.names = FALSE)
    ph$group <- factor(ifelse(tolower(ph$group) == "tbi", "TBI", "Sham"),
                       levels = c("Sham", "TBI"))
    f <- as.data.frame(t(res$fractions))
    f$sample <- rownames(f)
    f <- merge(f, ph[, c("sample", "group")], by = "sample")
    diff_tab <- lapply(setdiff(colnames(f), c("sample", "group")), function(ct) {
      p <- tryCatch(wilcox.test(reformulate("group", ct), data = f)$p.value,
                    error = function(e) NA_real_)
      data.frame(cell = ct, pvalue = p,
                 mean_sham = mean(f[[ct]][f$group == "Sham"]),
                 mean_tbi = mean(f[[ct]][f$group == "TBI"]))
    })
    diff_df <- do.call(rbind, diff_tab)
    diff_df <- diff_df[order(diff_df$pvalue), ]
    write.csv(diff_df, file.path(OUT, "immune_diff.csv"), row.names = FALSE)

    # hub 基因与免疫细胞相关
    coef_file <- file.path(ROOT, "results", "07_ML_Signature",
                           "enet_coefficients_legacy_04.csv")
    if (file.exists(coef_file)) {
      hub <- read.csv(coef_file)$gene
      hub <- intersect(hub, rownames(expr))
      if (length(hub) > 0) {
        e_sub <- expr[hub, rownames(f), drop = FALSE]
        cor_mat <- cor(t(e_sub), t(res$fractions[, rownames(f)]), method = "spearman")
        write.csv(cor_mat, file.path(OUT, "hub_immune_cor.csv"))
        if (requireNamespace("pheatmap", quietly = TRUE)) {
          pheatmap::pheatmap(cor_mat, cluster_cols = TRUE, cluster_rows = TRUE,
                             main = "hub gene vs immune cell (Spearman)",
                             filename = file.path(OUT, "hub_immune_heatmap.pdf"))
        }
      }
    }
  }
} else {
  message("未找到 LM22.txt（CIBERSORT 签名矩阵），跳过免疫浸润。")
  message("请在 cibersortx.stanford.edu 注册后下载 LM22.txt，放到 data/03_LM22/ 再重跑。")
}

# ---------- 2) 对接结构准备（hub 基因 -> 人源同源 -> UniProt -> AlphaFold/PDB） ----------
coef_file <- file.path(ROOT, "results", "07_ML_Signature",
                       "enet_coefficients_legacy_04.csv")
if (!file.exists(coef_file)) {
  message("未找到 enet_coefficients.csv，跳过对接结构准备。")
  quit(save = "no")
}
hub <- read.csv(coef_file)$gene
cat("hub genes:", paste(hub, collapse = ", "), "\n")

# 小鼠 -> 人同源（biomaRt；失败则原样输出并提示人工核对）
human_orth <- function(mouse_genes) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    message("  [warn] 未安装 biomaRt，使用原基因名（请人工核对是否为人类基因）")
    return(setNames(mouse_genes, mouse_genes))
  }
  tryCatch({
    mm <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    hs <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    out <- biomaRt::getLDS(attributes = "external_gene_name",
                           filters = "external_gene_name", values = mouse_genes,
                           mart = mm, attributesL = "external_gene_name",
                           martL = hs)
    setNames(toupper(out[, 2]), toupper(out[, 1]))
  }, error = function(e) {
    message("  [warn] biomaRt 查询失败: ", conditionMessage(e))
    setNames(mouse_genes, mouse_genes)
  })
}

uniprot_from_gene <- function(gene) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/search?query=gene_exact:",
                gene, "%20AND%20taxonomy_id:9606&fields=accession&format=json&size=1")
  tryCatch({
    d <- jsonlite::fromJSON(url)
    if (length(d$results) > 0) d$results$primaryAccession[1] else NA_character_
  }, error = function(e) NA_character_)
}

alphafold_pdb <- function(acc) {
  url <- paste0("https://alphafold.ebi.ac.uk/api/prediction/", acc)
  tryCatch({
    d <- jsonlite::fromJSON(url)
    d$pdbUrl[1]
  }, error = function(e) NA_character_)
}

orth <- human_orth(hub)
struct_rows <- lapply(seq_along(hub), function(i) {
  mouse_g <- hub[i]; human_g <- unname(orth[[mouse_g]])
  acc <- if (!is.na(human_g)) uniprot_from_gene(human_g) else NA_character_
  af <- if (!is.na(acc)) alphafold_pdb(acc) else NA_character_
  data.frame(mouse_gene = mouse_g, human_gene = human_g,
             uniprot = acc, alphafold_pdb_url = af)
})
struct_tab <- do.call(rbind, struct_rows)
write.csv(struct_tab, file.path(OUT, "hub_structure_table.csv"), row.names = FALSE)

# 下载 AlphaFold PDB（可选；AutoDock Vina 对接需在本地进行）
pdb_dir <- file.path(ROOT, "data", "docking")
dir.create(pdb_dir, showWarnings = FALSE)
for (i in seq_len(nrow(struct_tab))) {
  if (!is.na(struct_tab$alphafold_pdb_url[i])) {
    dst <- file.path(pdb_dir, paste0(struct_tab$uniprot[i], "_AF.pdb"))
    if (!file.exists(dst)) {
      tryCatch(download.file(struct_tab$alphafold_pdb_url[i], dst, mode = "wb", quiet = TRUE),
               error = function(e) message("  [warn] 下载失败: ", struct_tab$uniprot[i]))
    }
  }
}

message("完成。结果在 results/05_immune_docking/；PDB 结构在 data/docking/。")
message("对接请使用 AutoDock Vina 1.2（配体 BaP：SMILES C1=CC=C2C3=C4C(=CC2=C1)C=CC5=C4C(=CC=C5)C=C3）。")
