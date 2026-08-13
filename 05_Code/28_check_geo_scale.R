# QC：核对 01 下载脚本的缩放是否引入偏差
# 对比 GEO 原始 exprs 与本地存储 expr（GSE58485/GSE41345）
suppressMessages({
  library(GEOquery)
})

compare <- function(gse, probe_gene) {
  cat("=====", gse, "=====\n")
  g <- getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE)
  if (!is.list(g)) g <- list(g)
  raw_list <- lapply(g, function(e) exprs(e))
  stored <- read.csv(file.path("data/01_GEO_bulk", paste0(gse, "_expr.csv")),
                     row.names = 1, check.names = FALSE, fileEncoding = "UTF-8")
  cat("stored dim:", dim(stored), " raw dims:",
      paste(sapply(raw_list, function(m) paste(dim(m), collapse = "x")),
            collapse = " ; "), "\n")
  # 找原始矩阵中目标基因（符号需从平台注释获取；这里用存储基因名直接匹配原始探针名困难，
  # 改为比较整体分布：raw 合并后按相同基因名聚合前的列均值比例）
  raw_all <- do.call(cbind, raw_list)
  # 对原始矩阵按行名匹配存储基因名（若行名就是符号）
  common <- intersect(rownames(stored), rownames(raw_all))
  cat("common rownames (symbols):", length(common), "\n")
  if (length(common) > 5) {
    g1 <- common[1:5]
    for (gn in g1) {
      rv <- raw_all[gn, ]
      sv <- stored[gn, ]
      # 样本顺序可能不同；直接比较中位数比值
      ratio <- median(as.numeric(rv)) / median(as.numeric(sv))
      cat(gn, "raw median:", round(median(as.numeric(rv)), 4),
          " stored median:", round(median(as.numeric(sv)), 4),
          " ratio:", round(ratio, 3), "\n")
    }
  }
}

compare("GSE58485", NULL)
compare("GSE41345", NULL)
