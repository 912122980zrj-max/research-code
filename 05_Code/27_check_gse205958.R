# QC：核对 GSE205958 FPKM 原始值（诊断 21 与 26 的 10x 差异）
suppressMessages(library(readxl))
p <- "data/01_GEO_bulk/GSE205958_genes_fpkm_expression.xlsx"
dat <- as.data.frame(read_excel(p), stringsAsFactors = FALSE,
                     check.names = FALSE)
cn <- colnames(dat)
fpkm_cols <- grep("^FPKM\\.", cn)
gene <- toupper(as.character(dat$gene_name))
for (g in c("P2RY6", "CCNA2", "HSD11B1")) {
  idx <- which(gene == g)
  cat(g, "rows:", length(idx), "\n")
  if (length(idx)) {
    m <- as.matrix(dat[idx, fpkm_cols])
    cat("  mean per column (raw):", round(colMeans(m), 6), "\n")
    cat("  sum per column (raw):", round(colSums(m), 6), "\n")
  }
}
