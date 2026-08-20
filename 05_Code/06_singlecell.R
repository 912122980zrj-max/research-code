# ============================================================
# 06_singlecell.R —— 单细胞分析（Seurat + Harmony + CellChat
#                    + scTenifoldKnk 虚拟敲除）
# 运行：Rscript scripts/06_singlecell.R
# 输入：
#   A. 优先：用户提供的已整合 Seurat 对象 data/02_singlecell/integrated.rds
#      （例如人源 GSE209552 的整合对象，或 CEREBRI 图谱子集）
#   B. 自动：小鼠 GSE160763（10x 原始数据，GEO 补充文件，约数百 MB）
# 输出：results/11_SingleCell/
# 注意：原始 10x 数据很大，建议先用 GEO 下载或使用论文补充数据；
#       本脚本只做"分析骨架"，参数请按实际数据调整。
# ============================================================

suppressMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "results", "11_SingleCell")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts", "theme_pub.R"))

# 运行日志：手动跑的时候每一步都会记录到 results/11_SingleCell/06_run.log
LOG <- file.path(OUT, "06_run.log")
if (file.exists(LOG)) file.remove(LOG)
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(...), "\n")
  cat(msg)
  cat(msg, file = LOG, append = TRUE)
}
options(warn = 1)

# 物种：小鼠用 CellChatDB.mouse，人类用 CellChatDB.human
# 用法：Rscript scripts/06_singlecell.R human   （跑 GSE209552 人类数据时）
args <- commandArgs(trailingOnly = TRUE)
arg_sp <- tolower(args)
if (any(arg_sp %in% c("human", "h"))) {
  SPECIES <- "human"
} else if (any(arg_sp %in% c("mouse", "m"))) {
  SPECIES <- "mouse"
} else if (identical(Sys.getenv("SPECIES"), "human")) {
  SPECIES <- "human"
} else if (identical(Sys.getenv("SPECIES"), "mouse")) {
  SPECIES <- "mouse"
} else {
  human_dir0 <- file.path(ROOT, "data", "02_singlecell", "GSE209552")
  mouse_dir0 <- file.path(ROOT, "data", "02_singlecell", "GSE160763")
  if (dir.exists(human_dir0) && !dir.exists(mouse_dir0)) {
    SPECIES <- "human"
  } else if (dir.exists(mouse_dir0) && !dir.exists(human_dir0)) {
    SPECIES <- "mouse"
  } else {
    stop("请指定分析物种：Rscript scripts/06_singlecell.R human（GSE209552）或 mouse（GSE160763）")
  }
}
cat("SPECIES:", SPECIES, "\n")

# ---------- 1) 数据加载（本地优先，不再强制联网下载） ----------
find_10x_dirs <- function(base) {
  if (!dir.exists(base)) return(character(0))
  bc <- list.files(base, pattern = "barcodes\\.tsv(\\.gz)?$",
                   recursive = TRUE, full.names = TRUE)
  dirs <- unique(dirname(bc))
  dirs[vapply(dirs, function(d) {
    any(grepl("matrix\\.mtx(\\.gz)?$", list.files(d, full.names = TRUE)))
  }, logical(1))]
}

load_10x_to_seurat <- function(sc_dir) {
  dirs <- find_10x_dirs(sc_dir)
  if (length(dirs) == 0) {
    stop("在 ", sc_dir, " 下未找到 10x 矩阵（barcodes/matrix），请检查目录结构")
  }
  cat("找到 10x 样本目录:", length(dirs), "\n")
  objs <- lapply(dirs, function(d) {
    # Seurat 5 的 Read10X 只识别 features.tsv.gz；10x v2 的 genes.tsv.gz 需改名
    feat <- file.path(d, "features.tsv.gz")
    genes <- file.path(d, "genes.tsv.gz")
    if (!file.exists(feat) && file.exists(genes)) {
      file.rename(genes, feat)
      cat("  [fix] ", basename(d), ": genes.tsv.gz -> features.tsv.gz\n", sep = "")
    }
    counts <- Read10X(d)
    name <- basename(d)
    # min.cells=0：保留全部基因，保证各样本特征一致，合并后 JoinLayers 才能生效
    o <- CreateSeuratObject(counts = counts, project = name,
                            min.cells = 0, min.features = 200)
    o$sample <- name
    o
  })
  if (length(objs) > 1) merge(objs[[1]], objs[-1]) else objs[[1]]
}

rds_file <- file.path(ROOT, "data", "02_singlecell", "integrated.rds")
sc_dir <- file.path(ROOT, "data", "02_singlecell")
if (!dir.exists(sc_dir) && dir.exists(file.path(ROOT, "data", "singlescell"))) {
  sc_dir <- file.path(ROOT, "data", "singlescell")
  rds_file <- file.path(ROOT, "data", "singlescell", "integrated.rds")
}
# select dataset dir by species to avoid mixing human and mouse samples
dataset_dir <- if (SPECIES == "human") {
  file.path(sc_dir, "GSE209552")
} else {
  file.path(sc_dir, "GSE160763")
}
if (!dir.exists(dataset_dir)) {
  cat("[warn] dataset dir not found: ", dataset_dir,
      " | fall back to scanning whole data/02_singlecell\n", sep = "")
  dataset_dir <- sc_dir
}
if (file.exists(rds_file)) {
  cat("加载用户提供的整合对象:", rds_file, "\n")
  obj <- readRDS(rds_file)
} else if (length(find_10x_dirs(dataset_dir)) > 0) {
  cat("使用本地 10x 数据（", sc_dir, "）\n")
  obj <- load_10x_to_seurat(dataset_dir)
} else {
  cat("未找到 integrated.rds 或本地 10x 目录。\n")
  cat("请用以下任一方式准备单细胞数据（原始 10x 压缩包通常 >1 GB，建议用浏览器下载）：\n")
  cat("  A) 下载 GSE160763_RAW.tar：\n")
  cat("     https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160763/suppl/\n")
  cat("     解压到 data/02_singlecell/GSE160763/（脚本自动递归识别各样本 10x 目录）\n")
  cat("  B) 使用论文处理数据（如 GSE209552 人源 / CEREBRI 小鼠图谱），\n")
  cat("     保存为 data/02_singlecell/integrated.rds 后重跑\n")
  stop("单细胞数据缺失：请按上述说明准备 data/02_singlecell/ 后重跑")
}

# 样本-分组映射（TBI/Sham），用于组成分析与 CellChat
# Seurat v5: after merge(), join layers so downstream functions work reliably
if (exists("JoinLayers", mode = "function")) {
  obj[["RNA"]] <- tryCatch(JoinLayers(obj[["RNA"]]),
                           error = function(e) {
                             log_msg("JoinLayers 失败：", conditionMessage(e))
                             obj[["RNA"]]
                           })
}
log_msg("RNA layers: ", paste(Layers(obj, assay = "RNA"), collapse = ","))

map_file <- file.path(ROOT, "data", "02_singlecell", "sample_group.csv")
if (!file.exists(map_file) && dir.exists(file.path(ROOT, "data", "singlescell"))) {
  map_file <- file.path(ROOT, "data", "singlescell", "sample_group.csv")
}
if (!"group" %in% colnames(obj@meta.data)) {
  if (file.exists(map_file)) {
    mp <- read.csv(map_file, fileEncoding = "UTF-8")
    id_vec <- if ("sample" %in% colnames(obj@meta.data)) {
      as.character(obj$sample)
    } else {
      rownames(obj@meta.data)
    }
    obj$group <- mp$group[match(id_vec, as.character(mp$sample))]
    cat("已从 sample_group.csv 映射分组\n")
    # 仅保留 TBI / Sham 组（排除 PLX5622 等处理组）
    n_before <- ncol(obj)
    obj <- subset(obj, subset = group %in% c("TBI", "Sham"))
    cat("过滤后保留细胞数：", ncol(obj), "（原 ", n_before, "）\n")
  } else {
    obj$group <- "Sham"
    cat("[warn] 未提供 sample_group.csv（sample,group 两列），group 暂设 Sham；",
        "请补充分组后重跑组成分析/CellChat 部分\n")
  }
}

# ---------- 2) QC、标准化、聚类 ----------
mt_pat <- ifelse(SPECIES == "human", "^MT-", "^mt-")
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pat)
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                percent.mt < 20)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)
log_msg("PCA done, cells: ", ncol(obj))
if ("harmony" %in% rownames(installed.packages())) {
  obj <- RunHarmony(obj, group.by.vars = "orig.ident")
  log_msg("Harmony done")
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:20)
} else {
  obj <- RunUMAP(obj, dims = 1:20)
  obj <- FindNeighbors(obj, dims = 1:20)
}
obj <- FindClusters(obj, resolution = 1.2)
log_msg("clustering done, clusters: ", length(levels(Idents(obj))))

# FindAllMarkers 前再次确保 layers 已合并（subset 可能重新分层）
if (exists("JoinLayers", mode = "function")) {
  obj[["RNA"]] <- tryCatch(JoinLayers(obj[["RNA"]]),
                           error = function(e) obj[["RNA"]])
}
log_msg("RNA layers before markers: ", paste(Layers(obj, assay = "RNA"), collapse = ","))

if (!requireNamespace("presto", quietly = TRUE)) {
  log_msg("提示：未安装 presto，FindAllMarkers 会很慢；可运行 BiocManager::install('presto')")
}
markers <- tryCatch(
  FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5),
  error = function(e) {
    log_msg("FindAllMarkers 出错（已跳过，不影响注释）：", conditionMessage(e))
    data.frame()
  }
)
log_msg("cluster markers rows: ", nrow(markers))
write.csv(markers, file.path(OUT, "cluster_markers.csv"), row.names = FALSE)
# 检查点：聚类+marker 完成后保存，后续出错可从这里继续
saveRDS(obj, file.path(OUT, "seurat_checkpoint.rds"))
log_msg("checkpoint saved")

# ---------- 3) 细胞类型注释（小鼠默认；人类请换 marker） ----------
mouse_canonical <- c(
  Microglia = "P2ry12|Cx3cr1|Tmem119|Csf1r|Tyrobp|Aif1|Itgam|C1qa|C1qb|C1qc",
  Macrophage = "Lyz2|Cd68|Cd163|Mrc1|Msr1",
  Astrocyte = "Gfap|Aqp4|Slc1a3|Gja1|S100b|Aldh1l1",
  Neuron = "Snap25|Rbfox3|Syt1|Meg3|Nrgn|Stmn2",
  Oligodendrocyte = "Mbp|Mog|Plp1|Mobp|Cldn11|Mag",
  OPC = "Pdgfra|Cspg4|Olig1|Olig2",
  Endothelial = "Cldn5|Flt1|Pecam1|Vwf|Egfl7",
  Pericyte = "Rgs5|Pdgfrb|Kcnj8|Notch3",
  Fibroblast = "Col1a1|Dcn|Col1a2|Lum",
  SMC = "Acta2|Myh11|Tagln",
  T_NK = "Cd3d|Cd3e|Nkg7|Cd8a|Cd2",
  B = "Cd79a|Ms4a1|Cd19",
  Monocyte = "Cd14|Fcgr3a|Ly6c2|S100a8"
)
human_canonical <- c(
  Microglia = "P2RY12|CX3CR1|TMEM119|CSF1R|TYROBP|AIF1|ITGAM|C1QA|C1QB|C1QC",
  Macrophage = "LYZ|CD68|CD163|MRC1|MSR1",
  Astrocyte = "GFAP|AQP4|SLC1A3|GJA1|S100B|ALDH1L1",
  Neuron = "SNAP25|RBFOX3|SYT1|MEG3|NRGN|STMN2",
  Oligodendrocyte = "MBP|MOG|PLP1|MOBP|CLDN11|MAG",
  OPC = "PDGFRA|CSPG4|OLIG1|OLIG2",
  Endothelial = "CLDN5|FLT1|PECAM1|VWF|EGFL7",
  Pericyte = "RGS5|PDGFRB|KCNJ8|NOTCH3",
  Fibroblast = "COL1A1|DCN|COL1A2|LUM",
  SMC = "ACTA2|MYH11|TAGLN",
  T_NK = "CD3D|CD3E|NKG7|CD8A|CD2",
  B = "CD79A|MS4A1|CD19",
  Monocyte = "CD14|FCGR3A|S100A8|S100A9"
)
canonical <- if (SPECIES == "human") human_canonical else mouse_canonical
# 细胞类型注释：每个类型计算模块评分（AddModuleScore），簇取均值后按最高分注释
for (ct in names(canonical)) {
  genes <- intersect(strsplit(canonical[[ct]], "\\|")[[1]], rownames(obj))
  if (length(genes) >= 2) {
    obj <- AddModuleScore(obj, features = list(genes), name = "MS", assay = "RNA")
    # AddModuleScore 会给列名加数字后缀（如 MS1），这里重命名为 MS_<细胞类型>
    colnames(obj@meta.data)[ncol(obj@meta.data)] <- paste0("MS_", ct)
  }
}
ms_cols <- grep("^MS_", colnames(obj@meta.data), value = TRUE)
cl_means <- aggregate(obj@meta.data[, ms_cols, drop = FALSE],
                      by = list(cluster = as.character(Idents(obj))), FUN = mean)
assign_ct <- function(cl) {
  row <- cl_means[cl_means$cluster == cl, ms_cols, drop = FALSE]
  if (nrow(row) == 0) return(NA_character_)
  best <- names(which.max(row))[1]
  sub("^MS_", "", best)
}
celltype_vec <- vapply(levels(Idents(obj)), assign_ct, character(1))[as.character(Idents(obj))]
obj$celltype <- unname(celltype_vec)
obj$celltype <- factor(obj$celltype, levels = names(canonical))
obj$celltype[is.na(obj$celltype)] <- "Neuron"
Idents(obj) <- "celltype"
log_msg("cell type annotation done (module score)")

p_umap0 <- DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE,
                   cols = pal_cat) +
  theme_pub(base_size = 12) +
  labs(title = "UMAP (module-score annotation)")
ggsave_pub(file.path(OUT, "umap_celltype.pdf"), p_umap0,
       width = 8.5, height = 7, dpi = 300)

# ---------- 4) 组间组成（Ro/e） ----------
if (!"group" %in% colnames(obj@meta.data)) {
  message("meta.data 中没有 group 列（TBI/Sham），请从原始注释补充后重跑")
} else {
  tab <- table(obj$celltype, obj$group)
  rown <- rowSums(tab); coln <- colSums(tab); total <- sum(tab)
  roe <- (tab / rown) / (matrix(coln, nrow = nrow(tab), ncol = ncol(tab),
                                byrow = TRUE) / total)
  roe[!is.finite(roe)] <- 0
  write.csv(as.data.frame.matrix(roe), file.path(OUT, "roe_celltype.csv"))
  roe_mat <- log2(roe + 0.01)
  if (length(unique(as.vector(roe_mat))) > 1) {
    png_pub(file.path(OUT, "roe_heatmap.pdf"), width = 1000, height = 900, res = 150)
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      tryCatch(
        pheatmap::pheatmap(
          roe_mat, cluster_cols = FALSE,
          color = colorRampPalette(pal_rdbu)(100),
          border_color = NA,
          fontsize = 7, fontfamily = "Arial",
          main = "Ro/e (log2) cell composition"),
        error = function(e) log_msg("Ro/e 热图失败（跳过绘图）：", conditionMessage(e))
      )
    }
    dev.off()
  } else {
    log_msg("Ro/e 热图跳过：数值无差异")
  }
}

# ---------- 5) hub 基因定位 ----------
coef_file <- file.path(ROOT, "results", "07_ML_Signature",
                       "enet_coefficients_legacy_04.csv")
if (file.exists(coef_file)) {
  hub <- intersect(read.csv(coef_file)$gene, rownames(obj))
  if (length(hub) > 0) {
    png_pub(file.path(OUT, "hub_featureplots.pdf"), width = 1500,
        height = 300 * ceiling(length(hub) / 4), res = 150)
    print(FeaturePlot(obj, features = hub, ncol = 4,
                      cols = c("#F7F7F7", "#2171B5")))
    dev.off()
    obj <- AddModuleScore(obj, features = list(hub), name = "HubScore")
    png_pub(file.path(OUT, "hub_module_score.pdf"), width = 1000, height = 900, res = 150)
    print(FeaturePlot(obj, features = "HubScore1",
                      cols = c("#F7F7F7", "#2171B5")))
    dev.off()
  }
}

# ---------- 6) CellChat（免疫-基质串扰，TGF-β 通路） ----------
if (requireNamespace("CellChat", quietly = TRUE)) {
  suppressMessages(library(CellChat))
  if (SPECIES == "mouse") db <- CellChatDB.mouse else db <- CellChatDB.human
  run_chat <- function(sub) {
    cc <- createCellChat(object = sub, group.by = "celltype", assay = "RNA")
    cc@DB <- db
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc, type = "triMean")
    cc <- filterCommunication(cc, min.cells = 10)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    cc
  }
  if ("group" %in% colnames(obj@meta.data)) {
    objs <- SplitObject(obj, split.by = "group")
    chats <- lapply(objs, function(o) tryCatch(run_chat(o), error = function(e) NULL))
    chats <- chats[!vapply(chats, is.null, logical(1))]
    if (length(chats) >= 2) {
      merged <- mergeCellChat(chats, add.names = names(chats))
      pdf_pub(file.path(OUT, "cellchat_tgfb.pdf"), width = 10, height = 8)
      print(netVisual_aggregate(chats[[1]], signaling = "TGFb",
                                layout = "circle", vertex.label.cex = 0.7))
      print(netVisual_aggregate(chats[[2]], signaling = "TGFb",
                                layout = "circle", vertex.label.cex = 0.7))
      dev.off()
      saveRDS(merged, file.path(OUT, "cellchat_merged.rds"))
    } else {
      message("CellChat：分组对象不足，跳过对比分析")
    }
  }
}

# ---------- 7) scTenifoldKnk 虚拟敲除（目标细胞类群） ----------
if (requireNamespace("scTenifoldKnk", quietly = TRUE)) {
  target_ct <- ifelse("Microglia" %in% levels(Idents(obj)), "Microglia",
                      as.character(levels(Idents(obj))[1]))
  cat("虚拟敲除目标细胞类群:", target_ct, "\n")
  sub <- subset(obj, idents = target_ct)
  X <- as.matrix(GetAssayData(sub, assay = "RNA", slot = "data"))
  # scTenifoldKnk 的输入/输出格式以包内帮助为准（不同版本有差异）
  knk <- tryCatch(
    scTenifoldKnk::scTenifoldKnk(X = X, k = 3, q = 0.05, nCores = 1),
    error = function(e) {
      message("  [warn] scTenifoldKnk 运行失败: ", conditionMessage(e),
              "（请查看 help(scTenifoldKnk) 调整参数）")
      NULL
    }
  )
  if (!is.null(knk)) {
    dr <- knk$diffRegulation
    dr <- dr[order(dr$FC, decreasing = TRUE), ]
    write.csv(dr, file.path(OUT, paste0("scTenifoldKnk_", target_ct, ".csv")),
              row.names = FALSE)
  }
}

log_msg("analysis complete")
saveRDS(obj, file.path(OUT, "seurat_object.rds"))
message("完成。结果在 results/11_SingleCell/（UMAP、Ro/e、FeaturePlot、CellChat、虚拟敲除）。")
