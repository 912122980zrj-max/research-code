options(timeout = 6000)
options(BioC_mirror = "https://bioconductor.org")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
}
BiocManager::install("scDblFinder", update = FALSE, ask = FALSE)
cat("DONE scDblFinder:", requireNamespace("scDblFinder", quietly = TRUE), "\n")
