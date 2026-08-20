# 01_check_availability.R
# 用途：检查目标包及其依赖在当前 CRAN / Bioconductor 3.22 镜像中的可用性。
# 用法：Rscript 01_check_availability.R

options(timeout = 900)

cat("R:", R.version.string, "\n\n")

cran_repos <- c(
  CRAN = "https://mirrors.aliyun.com/CRAN/",
  CRAN2 = "https://mirrors.cloud.tencent.com/CRAN/",
  CRAN3 = "https://cloud.r-project.org"
)
bioc_repos <- c(
  BioCsoft = "https://bioconductor.org/packages/3.22/bioc",
  BioCann  = "https://bioconductor.org/packages/3.22/data/annotation",
  BioCexp  = "https://bioconductor.org/packages/3.22/data/experiment"
)

check_avail <- function(label, repo_urls) {
  cat("---", label, "---\n")
  db <- tryCatch(
    available.packages(repos = repo_urls, type = "source"),
    error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(db)) return(invisible(NULL))
  cat("总包数:", nrow(db), "\n")
  for (p in c(
    "rms", "rmda", "BisqueRNA", "limSolve", "BART", "ROCit", "compareC",
    "plsRcox", "randomForestSRC", "superpc", "survivalROC", "survivalsvm",
    "forestploter", "Ckmeans.1d.dp", "Boruta", "mixtools", "ggbreak",
    "miscTools", "compositions", "nnls", "MCMCpack", "fields", "ROCR",
    "KernSmooth", "TOAST", "GSEABase", "GSVA", "cancerclass", "mixOmics",
    "sparrow", "sva", "ComplexHeatmap", "SingleCellExperiment", "Biobase",
    "MuSiC", "CoxBoost", "fastAdaboost"
  )) {
    cat(sprintf("  %-20s %s\n", p, if (p %in% rownames(db)) "OK" else "MISSING"))
  }
  invisible(NULL)
}

check_avail("CRAN（阿里云/腾讯/官方，自动回退）", cran_repos)
check_avail("Bioc 3.22（官方源）", bioc_repos)
