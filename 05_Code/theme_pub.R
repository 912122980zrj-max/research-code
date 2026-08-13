# ============================================================
# theme_pub.R —— 发表级绘图主题与低饱和度配色（Scientific Reports / Nature 风格）
# 供全部 R 绘图脚本 source() 使用。
# 部署规范（2026-08-09，组图统一）：
#   - 字体：强制无衬线 Arial 族。注意本机 R/cairo 对字面 "Arial" 存在间歇性
#     "invalid font type" 缺陷（已复现），而 family="sans" 在 Windows cairo 下
#     稳定解析并内嵌 ArialMT（已用 pypdf 核验 /BaseFont=ArialMT、文字可提取）。
#     因此代码统一用 "sans"，产出 PDF 的嵌入字体即 Arial。
#   - 字号层级：图内主标签/图例标题 9 pt；坐标轴标签 8 pt；刻度数字 7 pt；
#     子图字母（a/b/c）一律不在 R 中绘制，交给 AI 排版；
#   - 导出：全部用 cairo_pdf（矢量、文字可编辑），不导 png/jpg；
#     包装函数 ggsave_pub / png_pub / pdf_pub 统一处理。
# ============================================================

# 干净、克制的期刊主题：细边框、浅灰网格、无花哨背景
theme_pub <- function(base_size = 8, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
      panel.border = element_rect(colour = "grey40", fill = NA,
                                  linewidth = 0.6),
      axis.ticks = element_line(colour = "grey40", linewidth = 0.4),
      axis.title = element_text(size = 8, colour = "black"),
      axis.text = element_text(size = 7, colour = "black"),
      plot.title = element_text(size = 9, face = "plain", hjust = 0.5),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = 9, colour = "black"),
      legend.text = element_text(size = 7, colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text = element_text(size = 9, colour = "black"),
      complete = TRUE
    )
}

# 统一导出包装：ggsave 强制 cairo_pdf（矢量、Arial 可编辑文字）
ggsave_pub <- function(filename, plot = ggplot2::last_plot(),
                       width = 7, height = 7, dpi = NULL, bg = "white", ...) {
  if (!grepl("\\.pdf$", filename)) {
    filename <- sub("\\.[^.]+$", ".pdf", filename)
  }
  grDevices::cairo_pdf(filename, width = width, height = height,
                       bg = bg, family = "sans")
  print(plot)
  grDevices::dev.off()
  invisible(filename)
}

# 统一导出包装：直接 png() 设备改为 cairo_pdf（px 尺寸自动换算为英寸）
png_pub <- function(filename, width = 7, height = 7, res = 300,
                    units = "in", type = "cairo", bg = "white",
                    pointsize = 12, ...) {
  if (is.null(units) || tolower(units) == "px") {
    width <- width / res
    height <- height / res
  }
  grDevices::cairo_pdf(filename = filename, width = width, height = height,
                       bg = bg, pointsize = pointsize, family = "sans", ...)
}

# 统一导出包装：直接 pdf() 设备改为 cairo_pdf（默认 Arial）
pdf_pub <- function(filename, ..., family = "sans") {
  grDevices::cairo_pdf(filename = filename, ..., family = family)
}

# 全局强制主题：任何脚本内 theme_bw(...) 也统一为 Arial + 9/8/7 字号层级
theme_bw <- function(base_size = 8, base_family = "sans", ...) {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family, ...) +
    theme(
      axis.title = element_text(size = 8, colour = "black"),
      axis.text = element_text(size = 7, colour = "black"),
      plot.title = element_text(size = 9, face = "plain", hjust = 0.5),
      legend.title = element_text(size = 9, colour = "black"),
      legend.text = element_text(size = 7, colour = "black"),
      strip.text = element_text(size = 9, colour = "black")
    )
}

# Paul Tol 柔和分类色（8 色，色盲友好，低饱和）
pal_cat <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77",
             "#CC6677", "#AA4499", "#882255")

# RdBu 柔和渐变（组间对比/热图）
pal_rdbu <- c("#2166AC", "#4393C3", "#92C5DE", "#F7F7F7", "#D6604D", "#B2182B")

# 蓝系渐变（连续值：富集 p 值、表达水平）
pal_blue <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6",
              "#4292C6", "#2171B5", "#084594")

# WGCNA 树状图 + 模块颜色条（修复：补齐标题、Height 纵轴、Genes 横轴与图例）
# 用法：draw_wgcna_dendro(geneTree, colors_vec, main = "...", out_file = "...pdf")
# geneTree：hclust 对象；colors_vec：与 geneTree$labels 等长的模块颜色向量
draw_wgcna_dendro <- function(geneTree, colors_vec, main = "",
                              out_file = "wgcna_dendro.pdf",
                              width = 10, height = 6.2) {
  n <- length(geneTree$order)  # 叶节点数（hclust labels 可能为空，order 始终有效）
  stopifnot(length(colors_vec) == n)
  ord <- geneTree$order
  col_ord <- colors_vec[ord]
  uniq_cols <- unique(col_ord)
  uniq_names <- unique(names(colors_vec))
  if (is.null(uniq_names) || length(uniq_names) != length(uniq_cols)) {
    # 无名字时用颜色本身做图例
    uniq_names <- uniq_cols
  }
  # 按出现顺序取图例（与颜色条一致）
  first_pos <- match(uniq_cols, col_ord)
  legend_df <- data.frame(col = uniq_cols, name = uniq_names,
                          pos = first_pos, stringsAsFactors = FALSE)
  legend_df <- legend_df[order(legend_df$pos), ]

  pdf_pub(out_file, width = width, height = height)
  layout(matrix(c(1, 2, 3), nrow = 3), heights = c(4.6, 1.2, 1.2))

  # 上：树状图（含标题与 Height 纵轴）
  par(mar = c(0.5, 4.5, 3.2, 1.5))
  plot(geneTree, main = main, ylab = "Height", xlab = NA,
       sub = NA, labels = FALSE,
       cex.main = 1.15, cex.lab = 1.0, cex.axis = 0.85,
       las = 1, axes = TRUE)
  box(col = "grey55", lwd = 0.6)

  # 下：模块颜色条 + Genes 横轴
  par(mar = c(3.2, 4.5, 0.3, 1.5))
  plot(0, type = "n", xlim = c(1, n), ylim = c(0, 1),
       axes = FALSE, xlab = NA, ylab = NA)
  for (i in seq_len(n)) {
    rect(i - 0.5, 0, i + 0.5, 1, col = col_ord[i], border = NA)
  }
  box(col = "grey55", lwd = 0.6)
  axis(1, at = c(1, n), labels = c("1", as.character(n)),
       cex.axis = 0.8, col.axis = "grey25")
  mtext("Genes (dendrogram order)", side = 1, line = 2.0, cex = 0.85,
        col = "grey15")
  mtext("Module colors", side = 4, line = 0.2, cex = 0.85, col = "grey15")

  # 图例（模块名 → 颜色），放在颜色条下方独立区域
  par(mar = c(0.5, 4.5, 0.2, 1.5))
  plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
       xlab = NA, ylab = NA)
  legend("center", legend = legend_df$name, fill = legend_df$col,
         border = "grey40", ncol = 6, cex = 0.62, bty = "n",
         x.intersp = 0.6, y.intersp = 0.9, text.width = 0.16)
  dev.off()
  invisible(out_file)
}
