# ============================================================
# 11_d9_cellchat.R —— D9 CellChat 细胞通讯（模板参数对齐）
# 依据：REVIEW4_FINAL_AUDIT D9
# 输入：results/11_SingleCell/seurat_object_reviewed.rds
# 输出：results/12_CellCommunication/D9_cellchat/
# 说明：CellChat v2.2.0.9001 已移除 projectData()（模板 v1 显式调用），
#       v2 将 PPI 投影作为内部处理，本脚本如实记录该差异。
# ============================================================

suppressMessages({
  library(Seurat)
  library(CellChat)
})

if (file.exists("config/datasets.tsv")) {
  ROOT <- getwd()
} else if (file.exists("../config/datasets.tsv")) {
  ROOT <- ".."
} else {
  stop("找不到 config/datasets.tsv")
}
OUT <- file.path(ROOT, "results", "12_CellCommunication", "D9_cellchat")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts", "theme_pub.R"))
report <- file.path(OUT, "D9_report.txt")
writeLines("D9 CellChat communication analysis", report)
write_report <- function(...) cat(..., "\n", sep = "", file = report,
                                  append = TRUE)
write_report("generated:", format(Sys.time()))
write_report("CellChat version:", as.character(packageVersion("CellChat")))
write_report("note: projectData() 在 CellChat v2 中已移除，PPI 投影为 v2 内部处理")

obj <- readRDS(file.path(ROOT, "results", "11_SingleCell",
                         "seurat_object_reviewed.rds"))
write_report("cells:", ncol(obj))
write_report("groups:", paste(names(table(obj$group)), collapse = ","))

db <- CellChatDB.human
db_use <- subsetDB(db, search = "Secreted Signaling")
write_report("CellChatDB.human secreted signaling interactions:",
             nrow(db_use$interaction))

run_chat <- function(grp) {
  sub <- subset(obj, subset = group == grp)
  cc <- createCellChat(object = sub, group.by = "celltype_reviewed",
                       assay = "RNA")
  cc@DB <- db_use
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, raw.use = TRUE, population.size = TRUE)
  cc <- filterCommunication(cc, min.cells = 3)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc)
  cc
}

chats <- list()
for (grp in c("Sham", "TBI")) {
  write_report("running CellChat for group:", grp, "...")
  chats[[grp]] <- tryCatch(run_chat(grp), error = function(e) {
    write_report("  [error]", grp, ":", conditionMessage(e))
    NULL
  })
}
chats <- chats[!vapply(chats, is.null, logical(1))]
if (length(chats) < 1) stop("CellChat 两个分组均失败")

for (grp in names(chats)) {
  cc <- chats[[grp]]
  df_net <- subsetCommunication(cc)
  write_report(grp, "interactions (pairwise, min.cells=3):", nrow(df_net))
  pw <- as.data.frame(cc@net$count)
  write.csv(pw, file.path(OUT, paste0("net_count_", grp, ".csv")))
  path_tab <- table(df_net$pathway_name)
  write_report(grp, "top pathways:",
               paste(names(sort(path_tab, decreasing = TRUE))[1:10],
                     collapse = ","))
  write.csv(df_net, file.path(OUT, paste0("communication_", grp, ".csv")),
            row.names = FALSE)
  # TGFb 网络
  if ("TGFb" %in% unique(df_net$pathway_name)) {
    png_pub(file.path(OUT, paste0("TGFb_network_", grp, ".pdf")),
        width = 1200, height = 1200, res = 300)
    tryCatch(netVisual_aggregate(cc, signaling = "TGFb",
                                 layout = "circle", vertex.label.cex = 0.7),
             error = function(e) write_report("  TGFb plot error:", 
                                               conditionMessage(e)))
    dev.off()
    lr <- tryCatch(subsetCommunication(cc, signaling = "TGFb"),
                   error = function(e) {
                     write_report("  TGFb LR error:", conditionMessage(e))
                     NULL
                   })
    write.csv(lr, file.path(OUT, paste0("TGFb_LR_", grp, ".csv")),
              row.names = FALSE)
    write_report(grp, "TGFb L-R pairs:", nrow(lr))
  } else {
    write_report(grp, "TGFb pathway not detected")
  }
}

if (length(chats) == 2) {
  merged <- mergeCellChat(chats, add.names = names(chats))
  saveRDS(merged, file.path(OUT, "cellchat_merged.rds"))
  ok_heat <- FALSE
  for (measure in c("weight", "count")) {
    png_pub(file.path(OUT, paste0("pathway_heatmap_diff_", measure, ".pdf")),
        width = 1200, height = 1000, res = 300)
    ok_heat <- tryCatch({
      netVisual_heatmap(merged, measure = measure)
      TRUE
    }, error = function(e) {
      write_report("  merged heatmap error (", measure, "):",
                   conditionMessage(e))
      FALSE
    })
    dev.off()
    if (ok_heat) break
  }
  if (!ok_heat) {
    # 回退：逐组 pathway heatmap
    for (grp in names(chats)) {
      png_pub(file.path(OUT, paste0("pathway_heatmap_", grp, ".pdf")),
          width = 1200, height = 1000, res = 300)
      tryCatch(netVisual_heatmap(chats[[grp]], measure = "weight"),
               error = function(e) write_report("  per-group heatmap error (",
                                                 grp, "):", conditionMessage(e)))
      dev.off()
    }
  }
  ok_role <- FALSE
  png_pub(file.path(OUT, "signaling_role_heatmap.pdf"), width = 1200,
      height = 900, res = 300)
  ok_role <- tryCatch({
    netAnalysis_signalingRole_heatmap(merged, pattern = "outgoing")
    TRUE
  }, error = function(e) {
    write_report("  merged role heatmap error:", conditionMessage(e))
    FALSE
  })
  dev.off()
  if (!ok_role) {
    for (grp in names(chats)) {
      png_pub(file.path(OUT, paste0("signaling_role_heatmap_", grp, ".pdf")),
          width = 1200, height = 900, res = 300)
      tryCatch(netAnalysis_signalingRole_heatmap(chats[[grp]],
                                                 pattern = "outgoing"),
               error = function(e) write_report("  per-group role heatmap error (",
                                                 grp, "):", conditionMessage(e)))
      dev.off()
    }
  }
  write_report("merged cellchat saved")
}

write_report("D9 finished:", format(Sys.time()))
message("D9 完成。结果在 ", OUT)
