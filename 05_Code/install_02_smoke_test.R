# 02_smoke_test.R
# 用途：对本次安装的关键包执行 library() 加载 + 关键函数存在性冒烟测试。
# 用法：Rscript 02_smoke_test.R

smoke <- function(pkg, funs = NULL) {
  ok <- tryCatch({
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    if (!is.null(funs)) {
      all(vapply(funs, function(f) exists(f, where = asNamespace(pkg), inherits = FALSE), logical(1)))
    } else TRUE
  }, error = function(e) {
    cat(sprintf("  %-16s LOAD FAIL: %s\n", pkg, conditionMessage(e)))
    FALSE
  })
  ver <- tryCatch(as.character(packageVersion(pkg)), error = function(e) "-")
  cat(sprintf("  %-16s %-5s %s\n", pkg, if (ok) "OK" else "FAIL", ver))
  invisible(ok)
}

cat("Smoke test start:", format(Sys.time()), "\n")
smoke("rms", "lrm")
smoke("rmda", "decision_curve")
smoke("Mime1", c("ML.Dev.Prog.Sig", "ML.Corefeature.Prog.Screen"))
smoke("CoxBoost")
smoke("fastAdaboost")
smoke("GSEABase")
smoke("GSVA", "gsva")
smoke("cancerclass")
smoke("mixOmics", "plsda")
smoke("sparrow")
smoke("ComplexHeatmap", "Heatmap")
smoke("TOAST")
smoke("MuSiC", "music_prop")
smoke("BisqueRNA", "ReferenceBasedDecomposition")
smoke("DoubletFinder", "doubletFinder")
smoke("SoupX", "load10X")
cat("Smoke test done:", format(Sys.time()), "\n")
