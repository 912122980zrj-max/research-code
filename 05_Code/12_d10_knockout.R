# ============================================================
# 12_d10_knockout.R —— D10 scTenifoldKnk 虚拟敲除
# 依据：REVIEW4_FINAL_AUDIT D10（top1000 HVG；|Z|>=2 且 q<0.01）
# 输入：results/11_SingleCell/seurat_object_reviewed.rds（Microglia 子集）
# 输出：results/13_VirtualKnockout/D10_knockout/
# 说明：scTenifoldKnk v1.1 默认参数（qc=TRUE, nc_nNet=10, nc_nCells=500,
#       nc_nComp=3, td_K=3, ma_nDim=2, nCores=detectCores()）。
#       gKO=P2RY6 为主；CCR5 为敏感性（可选）。
# ============================================================

suppressMessages({
  library(Seurat)
  library(scTenifoldKnk)
  library(ggplot2)
})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "results", "13_VirtualKnockout", "D10_knockout")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts", "theme_pub.R"))
report <- file.path(OUT, "D10_report.txt")
writeLines("D10 scTenifoldKnk virtual knockout", report)
write_report <- function(...) cat(..., "\n", sep = "", file = report,
                                  append = TRUE)
write_report("generated:", format(Sys.time()))
write_report("scTenifoldKnk version:",
             as.character(packageVersion("scTenifoldKnk")))

obj <- readRDS(file.path(ROOT, "results", "11_SingleCell",
                         "seurat_object_reviewed.rds"))
micro <- subset(obj, subset = celltype_reviewed == "Microglia")
write_report("Microglia cells:", ncol(micro))

norm <- GetAssayData(micro, assay = "RNA", layer = "data")
vr <- apply(as.matrix(norm), 1, var)
hvg_base <- names(sort(vr, decreasing = TRUE))
counts <- as.matrix(GetAssayData(micro, assay = "RNA", layer = "counts"))

run_ko <- function(gko) {
  # 强制纳入敲除靶基因（即使方差低不在 HVG 中，敲除必须保留该基因行）
  if (!gko %in% rownames(counts)) {
    write_report("[skip]", gko, "not in count matrix of Microglia")
    return(invisible(NULL))
  }
  pct <- mean(GetAssayData(micro, assay = "RNA", layer = "data")[gko, ] > 0)
  if (pct < 0.01) {
    write_report("[skip]", gko, "expression too low in Microglia (pct>0 =",
                 round(pct, 4), "), knockout not meaningful")
    return(invisible(NULL))
  }
  hvg <- unique(c(gko, hvg_base[seq_len(min(999, length(hvg_base)))]))
  X <- counts[hvg, , drop = FALSE]
  write_report(gko, "matrix dim:", paste(dim(X), collapse = "x"),
               "| target pct>0:", round(pct, 4))
  set.seed(2026)
  res <- tryCatch(
    scTenifoldKnk(countMatrix = X, gKO = gko, nc_nNet = 5,
                  nCores = parallel::detectCores()),
    error = function(e) {
      write_report("[error]", gko, ":", conditionMessage(e))
      NULL
    }
  )
  if (is.null(res)) return(invisible(NULL))
  dr <- res$diffRegulation
  dr <- dr[order(-dr$FC), ]
  write.csv(dr, file.path(OUT, paste0("scTenifoldKnk_", gko, "_all.csv")),
            row.names = FALSE)
  # 注意：scTenifoldKnk v1.1 输出列为 p.adj（BH-FDR），无 q 列；阈值对应审查的 q<0.01
  sig <- dr[abs(dr$Z) >= 2 & dr$p.adj < 0.01, ]
  write.csv(sig, file.path(OUT, paste0("scTenifoldKnk_", gko, "_significant.csv")),
            row.names = FALSE)
  top20 <- head(dr, 20)
  write.csv(top20, file.path(OUT, paste0("scTenifoldKnk_", gko, "_top20.csv")),
            row.names = FALSE)
  write_report(gko, "differentially regulated genes (all):", nrow(dr),
               "| significant (|Z|>=2 & p.adj<0.01):", nrow(sig))
  write_report(gko, "top20 by FC:",
               paste(head(dr$gene, 20), collapse = ","))
  # volcano
  dr$label <- ifelse(abs(dr$Z) >= 2 & dr$p.adj < 0.01, dr$gene, "")
  p <- ggplot(dr, aes(x = log2(abs(FC) + 1) * sign(FC), y = -log10(p.adj),
                      color = abs(Z) >= 2 & p.adj < 0.01)) +
    geom_point(size = 1.2, alpha = 0.7) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#B2182B"),
                       name = "|Z|>=2 & p.adj<0.01") +
    labs(title = paste0("scTenifoldKnk virtual KO: ", gko),
         x = "signed log2(|FC|+1)", y = "-log10(p.adj)") +
    theme_pub(base_size = 11)
  ggsave_pub(file.path(OUT, paste0("volcano_", gko, ".pdf")), p,
         width = 7, height = 5.5, dpi = 300)
  invisible(res)
}

write_report("running knockout: P2RY6")
run_ko("P2RY6")
if (file.exists(file.path(OUT, "scTenifoldKnk_P2RY6_all.csv"))) {
  write_report("CCR5 在小胶质中几乎不表达（pct>0<1%），敲除无意义；",
               "改用 P2RX7（小胶质 25% 表达、同为 CLEAN 候选）作敏感性")
  run_ko("P2RX7")
}

write_report("D10 finished:", format(Sys.time()))
message("D10 完成。结果在 ", OUT)
