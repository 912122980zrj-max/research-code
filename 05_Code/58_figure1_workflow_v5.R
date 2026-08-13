# ============================================================
# 58_figure1_workflow_v5.R
# Figure 1 (study design) for the v4 fixed numbers.
# v5.1: single-line header, larger/uniform fonts, straight
# fan-out arrows into the validation row (no bottom elbows).
# ============================================================

suppressMessages({
  library(ggplot2)
})

out_png <- file.path("figures", "Figure_01_Workflow", "Figure1_workflow_v5.png")
out_pdf <- file.path("figures", "Figure_01_Workflow", "Figure1_workflow_v5.pdf")

# ---------------- helper: rounded box ----------------
box_df <- function(id, x0, x1, y0, y1, title, lines, fill, edge,
                   title_col = "#17324D", body_col = "#222222",
                   dashed = FALSE, bold = FALSE) {
  data.frame(
    id = id, x0 = x0, x1 = x1, y0 = y0, y1 = y1,
    title = title, lines = I(list(lines)), fill = fill, edge = edge,
    title_col = title_col, body_col = body_col,
    dashed = dashed, bold = bold,
    stringsAsFactors = FALSE
  )
}

boxes <- rbind(
  # ---- Phase 1: inputs ----
  box_df("in1", 2, 30, 68, 79,
         "Discovery bulk cohorts",
         c("GSE58485  cortex CCI 3 d  (n = 24)",
           "GSE41345  hippocampus mTBI 7 d  (n = 9)",
           "CLEAN n = 33  (21 TBI / 12 sham)"),
         "#EAF1FB", "#3B6EA5"),
  box_df("in2", 35, 63, 68, 79,
         "BaP target prediction",
         c("ChEMBL · SEA · PharmMapper ·",
           "SwissTargetPrediction",
           "364 unique targets (A = 5, B = 23, C = 336)"),
         "#EAF1FB", "#3B6EA5"),
  box_df("in3", 68, 98, 68, 79,
         "External / contextual datasets",
         c("GSE128543  CCI 48 h, WT n = 14",
           "GSE205958  CCI 7 d RNA-seq, n = 10",
           "GSE104687  human aging, n = 376",
           "GSE209552  snRNA-seq, n = 17"),
         "#EAF1FB", "#3B6EA5"),

  # ---- Phase 2: preprocessing ----
  box_df("pp1", 2, 30, 55, 66,
         "Merged common-gene matrix",
         c("10,717 genes",
           "NA filter + per-cohort coverage filter",
           "(>=10% non-zero in both cohorts)"),
         "#FFFFFF", "#6B7A90"),
  box_df("pp2", 35, 63, 55, 66,
         "Coverage filter + ComBat",
         c("5,043 genes retained",
           "ComBat (batch = cohort, mod = ~group)",
           "PC1 cohort eta^2 0.9999 -> 0.0001"),
         "#FFFFFF", "#6B7A90"),
  box_df("pp3", 68, 98, 55, 66,
         "Target evaluability",
         c("human -> mouse symbol mapping",
           "65 / 364 targets evaluable (17.9%)",
           "tier A = 0, B = 6, C = 59"),
         "#FFFFFF", "#6B7A90"),

  # ---- Phase 3: discovery ----
  box_df("d1", 2, 28, 42, 52,
         "Differential expression (limma)",
         c("198 DEGs (177 up / 21 down)",
           "|log2FC| > 0.585, BH adjP < 0.05"),
         "#E5F0FF", "#2F5DA8"),
  box_df("d2", 34, 60, 42, 52,
         "WGCNA co-expression",
         c("top 8,000 MAD genes, softPower 12",
           "black module r = 0.827 (P = 2.9e-9)",
           "10 TBI-associated modules"),
         "#E5F0FF", "#2F5DA8"),
  box_df("d3", 2, 98, 31, 40.5,
         "Candidate intersection",
         c("198 DEG ∩ 1,244 module genes = 193 key genes",
           "∩ 65 evaluable BaP targets = 7 candidates (all tier C)"),
         "#D6E4FF", "#1F4E9C", bold = TRUE),

  # ---- Phase 4: enrichment + ML ----
  box_df("e1", 2, 28, 20, 30,
         "Enrichment & AhR test",
         c("GO BP 5 terms; KEGG 0",
           "AhR Fisher 0 / 3 (P = 1.00)",
           "not measurable after coverage filter"),
         "#FFF6E5", "#C07F00"),
  box_df("m1", 34, 60, 20, 30,
         "Machine learning & explanation",
         c("11 algorithms; ENet OOF",
           "AUC = 0.952 (0.916–0.974)",
           "SHAP top3: P2RY6, HSD11B1, CCR5"),
         "#FFF6E5", "#C07F00"),
  box_df("m2", 66, 98, 20, 30,
         "Frozen 3-gene signature",
         c("P2RY6 up, HSD11B1 down, CCR5 down",
           "non-nested selection: descriptive only",
           "(PROBAST+AI high risk)"),
         "#FFE9C2", "#B26A00", bold = TRUE),

  # ---- Phase 5: validation ----
  box_df("v1", 2, 30, 10.5, 19,
         "Same-model CCI cortex",
         c("GSE128543 (n = 14) & GSE205958 (n = 10)",
           "rank AUC = 1.0",
           "per-gene direction 2/3 (CCR5 unstable)"),
         "#EDF6E7", "#4E8C2E"),
  box_df("v2", 35, 61, 10.5, 19,
         "Leave-one-cohort-out",
         c("CLEAN 0.60 / 0.50",
           "FULL 0.70 / 0.50",
           "no robust transfer across heterogeneity"),
         "#FBEFE3", "#C07F00"),
  box_df("v3", 66, 98, 10.5, 19,
         "Human aging cohort",
         c("GSE104687 (n = 376)",
           "not evaluable",
           "(P2RY6 / CCR5 absent)"),
         "#FBEFE3", "#C07F00")
)

# ---------------- arrows (straight fan-out into validation) ----------------
arrows <- data.frame(
  x  = c(16, 49, 83, 30, 49, 49, 83, 15, 47, 15, 47, 82, 28, 60,
         66, 70, 86, 82),
  y  = c(68, 68, 68, 60.5, 55, 55, 55, 42, 42, 31, 31, 31, 25, 25,
         20, 20, 20, 10.5),
  xend = c(16, 49, 83, 35, 15, 47, 83, 15, 47, 15, 47, 82, 34, 66,
           30, 61, 86, 82),
  yend = c(66, 66, 66, 60.5, 52, 52, 40.5, 40.5, 40.5, 30, 30, 30, 25, 25,
           19, 19, 19, 8.5),
  color = c("#6B7A90", "#6B7A90", "#6B7A90", "#2F5DA8", "#2F5DA8", "#2F5DA8", "#2F5DA8",
            "#2F5DA8", "#2F5DA8", "#C07F00", "#C07F00", "#C07F00", "#C07F00", "#2F5DA8",
            "#4E8C2E", "#C07F00", "#C07F00", "#8A94A0"),
  dashed = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
             FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
             FALSE, FALSE, FALSE, TRUE),
  head = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
           TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
           TRUE, TRUE, TRUE, TRUE)
)

# ---------------- phase bands ----------------
phases <- data.frame(
  y0 = c(67, 54, 30, 19, 9.5),
  y1 = c(80, 67, 53, 31, 20),
  label = c("DATA & TARGET PREDICTION", "PREPROCESSING",
            "DISCOVERY & CANDIDATES", "ENRICHMENT & ML", "VALIDATION"),
  lab_y = c(73.5, 60.5, 41.5, 25, 14.75)
)

p <- ggplot() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 85), expand = FALSE) +
  theme_void() +

  # single concise header
  annotate("text", x = 50, y = 82.5, hjust = 0.5, size = 6.0,
           fontface = "bold", family = "sans", colour = "#14263F",
           label = "Predicted BaP targets and acute TBI transcriptomes: study design") +

  # phase bands
  geom_rect(data = phases,
            aes(xmin = 0.4, xmax = 99.6, ymin = y0, ymax = y1),
            fill = "#7F8CA3", alpha = 0.055, inherit.aes = FALSE) +
  geom_text(data = phases, aes(x = 1.25, y = lab_y, label = label),
            angle = 90, hjust = 0.5, size = 3.2, fontface = "bold",
            family = "sans", colour = "#5A6B82", inherit.aes = FALSE) +

  # contextual strip
  annotate("rect", xmin = 2, xmax = 98, ymin = 2.5, ymax = 8.5,
           fill = "#F5F6F7", colour = "#8A94A0", linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = 50, y = 7.4, hjust = 0.5, size = 3.4, fontface = "bold",
           family = "sans", colour = "#333333",
           label = "Contextual analyses (exploratory)") +
  annotate("text", x = 50, y = 5.2, hjust = 0.5, size = 2.9,
           family = "sans", colour = "#444444",
           label = "deconvolution: all negative   ·   CellChat: single-sample artifact   ·   virtual knockout: target-only   ·   docking: structural plausibility only")

# ---------------- boxes ----------------
for (i in seq_len(nrow(boxes))) {
  b <- boxes[i, ]
  lt <- if (b$dashed) "dashed" else "solid"
  p <- p +
    annotate("rect", xmin = b$x0, xmax = b$x1, ymin = b$y0, ymax = b$y1,
             fill = b$fill, colour = b$edge, linewidth = if (b$bold) 1.1 else 0.7,
             linetype = lt)
}
for (i in seq_len(nrow(boxes))) {
  b <- boxes[i, ]
  nline <- length(b$lines[[1]])
  title_h <- 4.8
  body_top <- b$y1 - title_h
  line_h <- min(2.6, (body_top - b$y0) / nline)
  title_size <- if (b$bold) 3.9 else 3.7
  p <- p +
    annotate("text", x = (b$x0 + b$x1) / 2, y = b$y1 - 1.9, hjust = 0.5,
             size = title_size, fontface = "bold", family = "sans",
             colour = b$title_col, label = b$title)
  for (j in seq_len(nline)) {
    y <- body_top - (j - 0.5) * line_h
    p <- p + annotate("text", x = (b$x0 + b$x1) / 2, y = y, hjust = 0.5,
                      size = 3.0, family = "sans", colour = b$body_col,
                      label = b$lines[[1]][j])
  }
}

# ---------------- arrows ----------------
for (i in seq_len(nrow(arrows))) {
  a <- arrows[i, ]
  ar <- if (a$head) arrow(length = unit(0.14, "cm"), type = "closed") else NULL
  p <- p + annotate("segment",
                    x = a$x, y = a$y, xend = a$xend, yend = a$yend,
                    colour = a$color, linewidth = 0.75,
                    linetype = if (a$dashed) "dashed" else "solid",
                    arrow = ar)
}

ggsave(out_png, p, width = 16, height = 13.2, dpi = 240, bg = "white")
ggsave(out_pdf, p, width = 16, height = 13.2, dpi = 300, bg = "white")
cat("Saved:", out_png, "\n", out_pdf, "\n")
