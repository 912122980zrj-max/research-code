# ============================================================================
# 一键安装 TBI 项目缺失的 R 包（Windows / R 4.5.3 / Bioc 3.22）
# ----------------------------------------------------------------------------
# 安装对象（2026-08-08 实测缺失）：
#   D3  : rms, rmda(CRAN 归档), Mime(GitHub, 包名 Mime1)
#   D5  : MuSiC(GitHub), BisqueRNA(GitHub)
#   D7  : DoubletFinder(GitHub), SoupX(GitHub)
# 策略：
#   1. CRAN 包：阿里云 -> 腾讯 -> USTC -> 官方 cloud 自动回退
#   2. Bioc 包：官方 Bioconductor 3.22（实测 USTC/TUNA 镜像当前不可用）
#   3. GitHub 包：remotes::install_github -> codeload 下载源码本地安装
#   4. rmda：已从 CRAN 移除（2026-05-04），从 CRAN 归档装 rmda_1.6.tar.gz
#   5. 已安装的包自动跳过；全部完成后验证 library() 并生成 INSTALL_REPORT.txt
#
# 用法：
#   Rscript install_missing_packages.R
#   或双击同目录 install_missing_packages.bat
# ============================================================================

options(timeout = 900)
options(warn = 1)

# ---------- 路径与日志 ----------
args0 <- commandArgs(FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_dir <- if (length(file_arg)) {
  normalizePath(dirname(sub("^--file=", "", file_arg[1])))
} else {
  normalizePath(getwd())
}

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_path <- file.path(script_dir, paste0("INSTALL_LOG_", stamp, ".txt"))
report_path <- file.path(script_dir, "INSTALL_REPORT.txt")
log_con <- file(log_path, open = "wt", encoding = "UTF-8")

log <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = " "))
  cat(msg, "\n")
  writeLines(msg, log_con)
  flush(log_con)
}

on.exit({
  if (isOpen(log_con)) close(log_con)
})

is_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

pkg_version <- function(pkg) {
  tryCatch(as.character(packageVersion(pkg)), error = function(e) NA_character_)
}

write_report <- function(lines) {
  writeLines(c(
    "# TBI 缺失 R 包一键安装报告",
    paste0("# 生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("# R 版本：", R.version.string),
    "",
    lines,
    ""
  ), report_path, useBytes = TRUE)
  log("报告已写入:", report_path)
}

# ---------- 镜像配置 ----------
CRAN_MIRRORS <- c(
  "https://mirrors.aliyun.com/CRAN/",
  "https://mirrors.cloud.tencent.com/CRAN/",
  "https://mirrors.ustc.edu.cn/CRAN/",
  "https://cloud.r-project.org"
)

BIOC_SOFT <- "https://bioconductor.org/packages/3.22/bioc"
BIOC_REPOS <- c(
  BioCsoft = BIOC_SOFT,
  BioCann  = "https://bioconductor.org/packages/3.22/data/annotation",
  BioCexp  = "https://bioconductor.org/packages/3.22/data/experiment"
)

# ---------- 安装函数 ----------
install_cran <- function(pkgs) {
  pkgs <- setdiff(pkgs, installed.packages()[, "Package"])
  if (!length(pkgs)) return(invisible(TRUE))
  log(">> CRAN 安装:", paste(pkgs, collapse = ", "))
  for (p in pkgs) {
    if (is_installed(p)) next
    for (m in CRAN_MIRRORS) {
      log("   尝试镜像:", m)
      ok <- tryCatch({
        install.packages(p, repos = c(CRAN = m), quiet = TRUE)
        is_installed(p)
      }, error = function(e) {
        log("   失败:", conditionMessage(e))
        FALSE
      })
      if (ok) {
        log("   成功:", p, pkg_version(p))
        break
      }
    }
    if (!is_installed(p)) log("   !! 未能安装:", p)
  }
  invisible(TRUE)
}

install_bioc <- function(pkgs) {
  pkgs <- setdiff(pkgs, installed.packages()[, "Package"])
  if (!length(pkgs)) return(invisible(TRUE))
  log(">> BiocManager 安装:", paste(pkgs, collapse = ", "))
  for (p in pkgs) {
    if (is_installed(p)) next
    ok <- tryCatch({
      BiocManager::install(p, update = FALSE, ask = FALSE, quiet = TRUE)
      is_installed(p)
    }, error = function(e) {
      log("   BiocManager 失败:", conditionMessage(e))
      FALSE
    })
    if (ok) {
      log("   成功:", p, pkg_version(p))
      next
    }
    # 回退：官方源下载 tar.gz 本地源码安装
    db <- tryCatch(
      available.packages(repos = BIOC_REPOS, type = "source"),
      error = function(e) NULL
    )
    if (!is.null(db) && p %in% rownames(db)) {
      ver <- db[p, "Version"]
      url <- sprintf("%s/src/contrib/%s_%s.tar.gz", BIOC_SOFT, p, ver)
      dest <- file.path(tempdir(), basename(url))
      dl <- tryCatch({
        download.file(url, dest, mode = "wb", quiet = TRUE)
        TRUE
      }, error = function(e) {
        log("   下载失败:", conditionMessage(e))
        FALSE
      })
      if (dl) {
        tryCatch({
          install.packages(dest, repos = NULL, type = "source", quiet = TRUE)
        }, error = function(e) log("   源码安装失败:", conditionMessage(e)))
      }
    }
    if (!is_installed(p)) log("   !! 未能安装:", p)
  }
  invisible(TRUE)
}

install_github_pkg <- function(repo, pkg, ref = "master") {
  if (is_installed(pkg)) {
    log("已安装（跳过）:", pkg, pkg_version(pkg))
    return(invisible(TRUE))
  }
  log(">> GitHub 安装:", repo, "->", pkg, "(ref:", ref, ")")
  ok <- tryCatch({
    remotes::install_github(
      repo = repo,
      ref = ref,
      upgrade = "never",
      build_vignettes = FALSE,
      quiet = FALSE
    )
    is_installed(pkg)
  }, error = function(e) {
    log("   remotes 安装失败:", conditionMessage(e))
    FALSE
  })
  if (ok) {
    log("   成功:", pkg, pkg_version(pkg))
    return(invisible(TRUE))
  }
  # 回退：codeload 下载源码 tar.gz 后本地安装
  tar_name <- sprintf("%s_%s.tar.gz", sub("/", "-", repo), ref)
  url <- sprintf("https://codeload.github.com/%s/tar.gz/refs/heads/%s", repo, ref)
  dest <- file.path(tempdir(), tar_name)
  log("   回退下载:", url)
  dl <- tryCatch({
    download.file(url, dest, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    log("   下载失败:", conditionMessage(e))
    FALSE
  })
  if (dl) {
    tryCatch({
      install.packages(dest, repos = NULL, type = "source", quiet = TRUE)
    }, error = function(e) log("   源码安装失败:", conditionMessage(e)))
  }
  if (!is_installed(pkg)) log("   !! 未能安装:", pkg)
  invisible(is_installed(pkg))
}

install_rmda_archive <- function() {
  if (is_installed("rmda")) {
    log("已安装（跳过）: rmda", pkg_version("rmda"))
    return(invisible(TRUE))
  }
  log(">> rmda：从 CRAN 归档安装 1.6")
  install_cran(c("caret", "pander", "MASS"))
  # rmda 依赖旧包 reshape（非 reshape2），已从 CRAN 移除，需从归档安装
  if (!is_installed("reshape")) {
    log("   reshape：从 CRAN 归档安装 0.8.9")
    url_r <- "https://cran.r-project.org/src/contrib/Archive/reshape/reshape_0.8.9.tar.gz"
    dest_r <- file.path(tempdir(), "reshape_0.8.9.tar.gz")
    dl_r <- tryCatch({
      download.file(url_r, dest_r, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) {
      log("   下载失败:", conditionMessage(e))
      FALSE
    })
    if (dl_r) {
      tryCatch({
        install.packages(dest_r, repos = NULL, type = "source", quiet = TRUE)
      }, error = function(e) log("   源码安装失败:", conditionMessage(e)))
    }
  }
  if (!is_installed("reshape")) log("   !! 未能安装 reshape（rmda 将无法安装）")
  url <- "https://cran.r-project.org/src/contrib/Archive/rmda/rmda_1.6.tar.gz"
  dest <- file.path(tempdir(), "rmda_1.6.tar.gz")
  dl <- tryCatch({
    download.file(url, dest, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    log("   下载失败:", conditionMessage(e))
    FALSE
  })
  if (dl) {
    tryCatch({
      install.packages(dest, repos = NULL, type = "source", quiet = TRUE)
    }, error = function(e) log("   源码安装失败:", conditionMessage(e)))
  }
  if (!is_installed("rmda")) log("   !! 未能安装: rmda")
  invisible(is_installed("rmda"))
}

# ---------- 开始安装 ----------
log("================================================================")
log("TBI 缺失 R 包一键安装开始")
log("R:", R.version.string)
log("库目录:", paste(.libPaths(), collapse = " | "))
log("================================================================")

# 基础设施（通常已装，确保存在）
install_cran(c("BiocManager", "remotes", "CoxBoost"))
if (!is_installed("CoxBoost")) install_github_pkg("binderh/CoxBoost", "CoxBoost", "master")

# ---- D3：机器学习规范化 ----
log("---- D3 组 ----")
install_cran("rms")
install_rmda_archive()

mime_bioc_deps <- c(
  "GSEABase", "GSVA", "cancerclass", "mixOmics", "sparrow", "sva",
  "ComplexHeatmap"
)
install_bioc(mime_bioc_deps)
install_github_pkg("souravc83/fastAdaboost", "fastAdaboost", "master")
install_github_pkg("l-magnificence/Mime", "Mime1", "main")

# ---- D5：免疫去卷积 ----
log("---- D5 组 ----")
install_bioc("TOAST")
install_cran(c("nnls", "MCMCpack", "limSolve"))
install_github_pkg("xuranw/MuSiC", "MuSiC", "master")
install_github_pkg("cozygene/bisque", "BisqueRNA", "master")

# ---- D7：双胞去除与环境 RNA 去污 ----
log("---- D7 组 ----")
install_github_pkg("chris-mcginnis-ucsf/DoubletFinder", "DoubletFinder", "master")
install_github_pkg("constantAmateur/SoupX", "SoupX", "master")

# ---------- 验证 ----------
log("================================================================")
log("安装结束，开始验证 ...")

verify <- c(
  "rms", "rmda", "Mime1", "CoxBoost", "fastAdaboost",
  "GSEABase", "GSVA", "cancerclass", "mixOmics", "sparrow", "ComplexHeatmap",
  "TOAST", "MuSiC", "BisqueRNA", "limSolve",
  "DoubletFinder", "SoupX"
)

report_lines <- character()
all_ok <- TRUE
for (p in verify) {
  ok <- is_installed(p)
  ver <- if (ok) pkg_version(p) else "-"
  status <- if (ok) "OK" else "FAIL"
  log(sprintf("%-18s %-5s %s", p, status, ver))
  report_lines <- c(
    report_lines,
    sprintf("| %-18s | %-5s | %s |", p, status, ver)
  )
  if (!ok) all_ok <- FALSE
}

report_lines <- c(
  "## 验证结果",
  "",
  "| 包 | 状态 | 版本 |",
  "|----|------|------|",
  report_lines,
  "",
  sprintf("**总体：%s**", if (all_ok) "全部安装成功" else "存在失败项，请查看上方日志"),
  "",
  "## 重要说明",
  "1. Mime 的 GitHub 仓库已迁移为 l-magnificence/Mime，包名是 **Mime1**，",
  "   加载用 `library(Mime1)`（不是 library(Mime)）；D3 脚本需相应调整。",
  "2. rmda 已于 2026-05-04 从 CRAN 移除，本次由 CRAN 归档安装 1.6。",
  "3. MuSiC / BisqueRNA 已不在 Bioconductor 发布版中，本次从官方 GitHub 安装。",
  "4. USTC/TUNA 的 Bioc 镜像当前实测不可用（404/403），本次 Bioc 包使用官方源。",
  "5. 完整安装日志：", log_path
)
write_report(report_lines)

log("完成。报告：", report_path)
