# 03_inspect_exports.R —— 调试：查看 BisqueRNA / DoubletFinder 能否加载及其导出函数

cat("== BisqueRNA ==", "\n")
r <- tryCatch({
  library(BisqueRNA)
  cat("loaded OK\n")
  print(ls("package:BisqueRNA"))
}, error = function(e) cat("LOAD ERROR:", conditionMessage(e), "\n"))

cat("\n== DoubletFinder ==", "\n")
r2 <- tryCatch({
  library(DoubletFinder)
  cat("loaded OK\n")
  print(ls("package:DoubletFinder"))
}, error = function(e) cat("LOAD ERROR:", conditionMessage(e), "\n"))
