# 生成 reports/METHODS_REPRODUCIBILITY.md（P0-8/9）
si <- capture.output(sessionInfo())
pkgs <- c("Seurat", "harmony", "scDblFinder", "CellChat", "limma", "sva",
          "WGCNA", "clusterProfiler", "caret", "glmnet", "pROC", "shapviz",
          "MuSiC", "BisqueRNA", "scTenifoldKnk", "readxl", "dplyr", "ggplot2")
vers <- vapply(pkgs, function(p) {
  if (requireNamespace(p, quietly = TRUE)) {
    as.character(packageVersion(p))
  } else "NOT INSTALLED"
}, character(1))

order <- c(
  "01_download_geo.R", "01b_download_validation.R", "02_bap_targets.py",
  "02b_tier_targets.py", "02c_tier_dedup.R", "03_merge_de_wgcna.R",
  "04_ml_shap.R", "04b_external_validation.R", "05_immune_docking.R",
  "06_singlecell.R", "06b_singlecell_review.R", "07_publication_figures.R",
  "08_d1_wgcna_sensitivity.R", "09_d2_go_all.R", "10_d11_ahr_fisher.R",
  "11_d9_cellchat.R", "12_d10_knockout.R", "13_d8_astro.R",
  "14_cibersortx_sigmatrix.R", "15_cibersortx_mixture.R",
  "16_deconvolution_local.R", "17_deconvolution_analysis.R",
  "18_d7_doublet_qc.R", "19_d3_ml_benchmark.R", "20_d4_shap.R",
  "21_d12_validation.R",
  "22_d8_annotation_review.R", "23_d7_scdblfinder.R",
  "25_composition_fdr.R", "26_external_validation_paper.R",
  "29_fix_geo_scale.R", "30_figures_2_5d.R", "31_h3_cohort_sensitivity.R",
  "32_sessioninfo.R", "33_workflow_volcano.R", "34_log2_unify.R",
  "58_figure1_workflow_v5.R", "59_build_manuscripts_docx.py",
  "60_build_submission_aux.py", "61_orthology_mapping.R",
  "62_permutation_background.R", "63_sensitivity_matrix.R")

md <- c(
  "# 方法学与可复现性说明（P0-8/9）",
  "",
  paste("> 生成时间：", Sys.time()),
  "> 依据：FINAL_REVIEW_PUBLISHABILITY.md §7 P0-8/9",
  "",
  "## 1. 运行环境与关键包版本",
  "",
  paste0("- R：", R.version.string),
  paste0("- 平台：", R.version$platform),
  paste0("- 关键包版本：", paste(pkgs, vers, sep = " = ", collapse = "；")),
  "",
  "## 2. 完整 sessionInfo",
  "",
  "```r",
  si,
  "```",
  "",
  "## 3. 脚本运行顺序（数据准备 → 分析 → 复核）",
  "",
  paste0(seq_along(order), ". `", order, "`", collapse = "\n"),
  "",
  "## 4. 关键参数速查",
  "",
  "- DEG：limma（线性尺度表达；设计 ~group+genotype，队列内奇异时降级 ~group）；",
  "  阈值 adj.P<0.05 & |logFC|>0.585（与主脚本一致）；",
  "- ComBat：par.prior=TRUE，mod=~group（全样本估计，CV 为描述性）；",
  "- WGCNA：top 8000 MAD 基因，softPower 12（CLEAN），模块-性状 Pearson；",
  "- ML：7 候选，raw CLEAN n=33；D3 11 算法基准 + ENet alpha=0.9 OOF；",
  "  D4 SHAP 解释 D3 同一拟合（fit_enet.rds）；seed=2026；",
  "- scDblFinder 1.24.10：按 17 样本分别检测，seed=2026；",
  "- CellChat v2.2.0.9001：Secreted Signaling，min.cells=3，raw.use=TRUE；",
  "- scTenifoldKnk：nc_nNet=5（参数已如实记录）；",
  "- 外部验证：固定冻结系数、不重训、不插补；GSE205958 主口径 FPKM 原值。",
  "",
  "## 5. 数据口径",
  "",
  "- CLEAN 主结果：GSE58485(24)+GSE41345(9)，n=33（剔除环磷酰胺与 Ex-4-mTBI）；",
  "- 表达尺度修复：见 `data/01_GEO_bulk/SCALE_REPAIR.md`（01/01b 曾除以样本数，",
  "  已乘回 ncol 恢复原始探针求和尺度并重跑全部 bulk 链）；",
  "- config/datasets.tsv 为数据集清单唯一权威。"
)
writeLines(md, "reports/METHODS_REPRODUCIBILITY.md")
cat("METHODS_REPRODUCIBILITY.md written\n")
