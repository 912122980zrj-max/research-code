# 52_ko_volcano_pdf.R —— 仅重绘 scTenifoldKnk 火山图（不重跑分析），输出 Arial 矢量 PDF
suppressMessages({ library(ggplot2); library(dplyr) })
source("scripts/theme_pub.R")
ROOT <- "E:/sheng xin/tbi"
KO <- file.path(ROOT, "results", "13_VirtualKnockout", "D10_knockout")
FIG <- file.path(ROOT, "figures", "Figure_09_Comm_KO_Astro")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

for (g in c("P2RY6", "P2RX7")) {
  d <- read.csv(file.path(KO, paste0("scTenifoldKnk_", g, "_all.csv")),
                stringsAsFactors = FALSE)
  d$sig <- d$p.adj < 0.01 & abs(d$Z) >= 2
  p <- ggplot(d, aes(x = Z, y = -log10(p.adj + 1e-300))) +
    geom_point(aes(colour = sig), size = 0.8, alpha = 0.7) +
    scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "#CC6677"),
                        labels = c("FALSE" = "not significant",
                                   "TRUE" = "|Z|>=2 & adjP<0.01")) +
    labs(x = "Z score", y = "-log10(adjusted P)", colour = NULL) +
    theme_pub()
  ggsave_pub(file.path(KO, paste0("volcano_", g, ".pdf")), p,
             width = 5.5, height = 4)
  file.copy(file.path(KO, paste0("volcano_", g, ".pdf")),
            file.path(FIG, paste0("volcano_", g, ".pdf")), overwrite = TRUE)
  cat("saved volcano", g, "\n")
}
