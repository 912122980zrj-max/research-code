# ============================================================
# 22_d8_annotation_review.R —— M9：D8 星形胶质亚群命名复核
# 运行：Rscript scripts/22_d8_annotation_review.R
# 前置：13_d8_astro.R 已生成 astro_subset_v2.rds / markers / scores
# 说明（2026-08-08）：
#   - D8 原输出只有簇编号（0–8）与模块评分，无正式命名；
#   - 本脚本基于 top10 markers + 4 个模块评分 + 按样本比例统计，
#     给出证据驱动的命名、置信度与警示，并将命名写入 RDS 元数据；
#   - 明确标注单样本/双样本驱动的簇（2/4/5）与低置信簇（8），
#     不做超出证据的生物学命名。
# 输出（原位，不新建打包）：
#   astro_cluster_annotation.csv / astro_subcluster_group_stats.csv /
#   astro_purity_check.csv / astro_umap_named.pdf /
#   ASTROCYTE_ANNOTATION_REVIEW.md / 更新 astro_subset_v2.rds
# ============================================================

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "results", "14_Astrocyte_Subclusters", "D8_astro")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts", "theme_pub.R"))

astro <- readRDS(file.path(OUT, "astro_subset_v2.rds"))
astro$subcluster <- as.character(astro$RNA_snn_res.0.5)
cat("clusters:", paste(sort(unique(astro$subcluster)), collapse = ", "), "\n")
cat("cells:", ncol(astro), " Sham:", sum(astro$group == "Sham"),
    " TBI:", sum(astro$group == "TBI"), "\n")

# ---------- 1. 命名表（证据驱动，含置信度与警示） ----------
ann <- data.frame(
  cluster = as.character(0:8),
  subcluster_name = c(
    "GFAP-high / mitochondrial-rich astrocyte",
    "Homeostatic astrocyte (SLC1A2-high)",
    "OSMR/ABCC9-high reactive-like astrocyte",
    "TBI-associated NAMPT-high reactive-like astrocyte",
    "APOE/CST3-high, blood/ambient-enriched (HBA2+)",
    "SLC1A2-high homeostatic-like, sample-enriched",
    "Reactive astrocyte (C3+/ECM+)",
    "Proliferative/ECM-remodeling astrocyte (GINS3/ADAMTS10-high)",
    "Low-confidence neuronal-like/ambient (EDIL3+)"
  ),
  short_label = c(
    "C0 GFAP-high/mt",
    "C1 Homeostatic",
    "C2 OSMR/ABCC9+",
    "C3 NAMPT+ TBI",
    "C4 APOE+ blood/ambient",
    "C5 SLC1A2+ sample-enriched",
    "C6 Reactive C3+/ECM+",
    "C7 GINS3+ proliferative/ECM",
    "C8 low-conf neuronal-like"
  ),
  evidence_markers = c(
    "GFAP, BCL6, SORBS1, FAM107A, MT-ND/CO genes",
    "SLC1A2, NRXN1, NKAIN3, NTRK2, GPM6A, CTNNA2",
    "OSMR, ABCC9, LYVE1, HSPA4L, CNTNAP3, PTGES3",
    "NAMPT, SCG2, LINGO1, RASGEF1B, ELL2, TESC",
    "MTRNR2L1, HBA2, APOE, CST3, CLU, GAPDH, GADD45B",
    "ARL17B, TRPM3, SLC1A2, ABLIM3, FAM198B-AS1",
    "C3, TNC, WWTR1, NAMPT, GAP43, CADPS; ECM/TGFb/stress/contractile scores top",
    "GINS3, ADAMTS10, AIF1L, DPYSL5, ETV5, TRAF4, MALAT1",
    "SLC24A2, GRIN2B, IL1RAPL1, KIRREL3, GALNT13, EDIL3, ANK3"
  ),
  confidence = c(
    "medium",
    "medium-high",
    "medium-low",
    "medium",
    "low (biological subtype not supported)",
    "medium-low",
    "high",
    "low",
    "low"
  ),
  caveat = c(
    "Sham-skewed (RO/E 1.40); mt 基因可能反映代谢状态或 ambient，5 个样本",
    "组间相对均衡；含部分神经元样基因（NRXN1/GPM6A），不排除 ambient",
    "400 个细胞但仅来自 2 个 Sham 样本（donor 驱动）；LYVE1 提示血管周/ambient",
    "TBI-skewed（RO/E 3.42，9 个样本但 nr20/nr21 贡献大）",
    "148 个细胞中 146 个来自单一样本 nr11_LFP；HBA2/GAPDH 提示血液/ambient，不应作为亚群结论",
    "141 个细胞中 138 个来自单个 Sham 样本 502T",
    "TBI-skewed（RO/E 3.40，12 个样本）；C3/TNC 为反应性星形经典证据",
    "84 个细胞中 80 个来自单一样本 nr3_RF；样本级 Wilcoxon p=0.72（BH 0.72），不可作群体性结论；snRNA 增殖基因可能含 ambient/伪影",
    "12 个样本但每样本细胞数少；神经元样基因可能为注释/ambient 不确定性"
  ),
  stringsAsFactors = FALSE
)
astro$subcluster_name <- ann$subcluster_name[match(astro$subcluster,
                                                   ann$cluster)]
astro$subcluster_short <- ann$short_label[match(astro$subcluster,
                                                ann$cluster)]
write.csv(ann, file.path(OUT, "astro_cluster_annotation.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ---------- 2. 按样本组间比例统计（Wilcoxon + BH） ----------
sample_tab <- as.data.frame.matrix(table(astro$sample, astro$subcluster))
sample_tab$sample <- rownames(sample_tab)
sample_tab$group <- astro$group[match(rownames(sample_tab), astro$sample)]
sample_tab$total_astro <- rowSums(sample_tab[, as.character(0:8)])
stats_list <- lapply(as.character(0:8), function(k) {
  prop <- sample_tab[[k]] / sample_tab$total_astro
  sham <- prop[sample_tab$group == "Sham"]
  tbi <- prop[sample_tab$group == "TBI"]
  p <- tryCatch(wilcox.test(tbi, sham)$p.value, error = function(e) NA_real_)
  data.frame(cluster = k,
             n_sham_samples = sum(sample_tab$group == "Sham"),
             n_tbi_samples = sum(sample_tab$group == "TBI"),
             n_cells_sham = sum(sample_tab[[k]][sample_tab$group == "Sham"]),
             n_cells_tbi = sum(sample_tab[[k]][sample_tab$group == "TBI"]),
             median_prop_sham = round(median(sham), 4),
             median_prop_tbi = round(median(tbi), 4),
             mean_prop_sham = round(mean(sham), 4),
             mean_prop_tbi = round(mean(tbi), 4),
             wilcox_p = round(p, 5),
             stringsAsFactors = FALSE)
})
stats <- do.call(rbind, stats_list)
stats$wilcox_p_bh <- round(p.adjust(stats$wilcox_p, method = "BH"), 5)
stats <- merge(stats, ann[, c("cluster", "subcluster_name")],
               by = "cluster", all.x = TRUE)
write.csv(stats, file.path(OUT, "astro_subcluster_group_stats.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("===== group stats =====\n")
print(stats)

# ---------- 3. 纯度核查（每簇典型标记的 log-normalized 均值） ----------
markers_check <- c(
  Astrocyte = "GFAP", Astrocyte2 = "AQP4", Astrocyte3 = "SLC1A2",
  Astrocyte4 = "ALDH1L1", Astrocyte5 = "GJA1",
  Microglia = "P2RY12", Microglia2 = "CX3CR1",
  Oligodendrocyte = "MBP", Oligo2 = "PLP1",
  Neuron = "SNAP25", Neuron2 = "RBFOX3",
  T_NK = "CD3D", Endothelial = "FLT1", Endo2 = "CLDN5"
)
present <- intersect(markers_check, rownames(astro))
ex <- as.matrix(GetAssayData(astro, assay = "RNA", layer = "data"))
purity <- do.call(rbind, lapply(as.character(0:8), function(k) {
  idx <- which(astro$subcluster == k)
  m <- rowMeans(ex[present, idx, drop = FALSE])
  data.frame(cluster = k, t(m), check.names = FALSE)
}))
write.csv(purity, file.path(OUT, "astro_purity_check.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("===== purity check =====\n")
print(purity)

# ---------- 4. 命名 UMAP ----------
p_umap <- DimPlot(astro, group.by = "subcluster_short",
                  label = TRUE, repel = TRUE,
                  cols = c(pal_cat, "#999999")) +
  theme_pub(base_size = 11) +
  labs(title = "Astrocyte subclusters (reviewed names, res=0.5)") +
  theme(legend.position = "right")
ggsave_pub(file.path(OUT, "astro_umap_named.pdf"), p_umap,
       width = 9, height = 6, dpi = 300, bg = "white")

# ---------- 5. 复核报告 ----------
md <- c(
  "# D8 星形胶质亚群命名复核（M9）",
  "",
  paste("> 生成时间：", Sys.time()),
  paste0("> 输入：`astro_subset_v2.rds`（", ncol(astro),
         " 个 Astrocyte 核；Sham ", sum(astro$group == "Sham"),
         "，TBI ", sum(astro$group == "TBI"), "）"),
  "> 依据：FindAllMarkers top10 + 4 个模块评分 + 按样本比例 Wilcoxon/BH",
  "",
  "## 1. 结论摘要",
  "",
  "- 簇 6 为证据最强的**反应性星形胶质（C3+/ECM+，TBI 富集，12 个样本）**；",
  "- 簇 1 为**稳态/原浆样（SLC1A2-high）**，组间相对均衡；",
  "- 簇 3 为 TBI 富集的反应性样（NAMPT-high，9 个样本，置信度中）；",
  "- 簇 0 为 GFAP-high/mt-rich（Sham 偏倚，置信度中）；",
  "- **簇 2/4/5/7 由 1–2 个样本驱动，不可作为群体性亚群结论**：",
  "  簇 2（400 细胞，2 个 Sham 样本）、簇 4（148 细胞中 146 来自 nr11_LFP，",
  "  含 HBA2/GAPDH 血液/ambient 信号）、簇 5（141 细胞中 138 来自 502T）、",
  "  簇 7（84 细胞中 80 来自 nr3_RF，样本级 p=0.72）；",
  "- 簇 8 为低置信（12 个样本但每样本细胞少，神经元样基因，疑似 ambient）。",
  "",
  "## 2. 逐簇命名与证据",
  "",
  "见 `astro_cluster_annotation.csv`（cluster / name / evidence / confidence / caveat）。",
  "",
  "## 3. 按样本组间比例（Wilcoxon，BH 校正）",
  "",
  "见 `astro_subcluster_group_stats.csv`。要点：",
  paste0("- 仅当 n_cells_sham 与 n_cells_tbi 均 > 30 且样本数 ≥ 5 时，",
         "比例差异才具初步可解释性；否则按样本驱动处理。"),
  "- RO/E（`astro_roe.csv`）为细胞层面比值，受零计数影响大，",
  "  本复核一律以样本级比例统计为准。",
  "",
  "## 4. 纯度核查",
  "",
  "见 `astro_purity_check.csv`。所有簇均保留星形标记（GFAP/AQP4/SLC1A2），",
  "但簇 4（HBA2+）与簇 8（神经元样基因）存在 ambient/污染信号，",
  "论文中不应将其作为独立星形亚群功能结论。",
  "",
  "## 5. 对论文/湿实验的影响",
  "",
  "- 建议正文只写：星形胶质内存在**反应性亚群（C3+/TNC+，TBI 富集）**",
  "  与**稳态亚群（SLC1A2-high）**的区分；",
  "- 簇 2/4/5/7 以“样本驱动簇”表述或在敏感性分析中剔除；",
  "- 若用于验证/去卷积参考，只采用 6/1（必要时 0/3）作为稳健代表；",
  "  样本级统计中仅簇 6 通过 BH<0.05（p=0.0047，BH=0.042）。",
  "",
  "## 6. 局限",
  "",
  "- 原始 13_d8_astro.R 未显式 set.seed，簇编号以 RDS 冻结结果为准；",
  "- snRNA-seq 中增殖/血红蛋白/神经元样基因可能来自 ambient，未做 SoupX（无空滴矩阵）；",
  "- 对照（尸检）与 TBI（手术）组织的可比性限制延续自单细胞主分析。"
)
writeLines(md, file.path(OUT, "ASTROCYTE_ANNOTATION_REVIEW.md"))

# ---------- 6. 更新 RDS（追加命名列，默认身份仍为簇编号）与报告 ----------
Idents(astro) <- factor(astro$subcluster, levels = as.character(0:8))
saveRDS(astro, file.path(OUT, "astro_subset_v2.rds"))
write("M9 naming review done: 2026-08-08",
      file.path(OUT, "D8_report.txt"), append = TRUE)
message("M9 完成。结果在 ", OUT)
