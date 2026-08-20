# ============================================================
# 82_figure2_v7.R —— 按 v7 数据重绘 Figure 2（单张合成 PDF）
# 面板 A：364 预测靶点证据分级（A=5, B=23, C=336）
# 面板 B：GSE58485-only 主路径 11 个候选（全 tier C）
# 面板 C：69 可评估靶点 / 373 DEGs / 364 关键基因 三集合韦恩图（交集 11）
# 输出：figures/Figure_02_BaP_Targets/Figure2_v7.pdf（单页）
# ============================================================
suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggvenn)
  library(patchwork)
})

ROOT <- getwd()
source(file.path(ROOT, "scripts", "theme_pub.R"))

# ---------- 数据 ----------
tier_tab <- read_csv(file.path(ROOT, "results", "16_Orthology",
                               "bap_targets_364_tier.csv"),
                     col_types = cols(.default = col_character()))
tab_tier <- data.frame(
  tier = c("A", "B", "C"),
  n = c(sum(tier_tab$tier == "A", na.rm = TRUE),
        sum(tier_tab$tier == "B", na.rm = TRUE),
        sum(tier_tab$tier == "C", na.rm = TRUE)))
tab_tier$tier <- factor(tab_tier$tier, levels = c("A", "B", "C"))
stopifnot(sum(tab_tier$n) == 364)

cand <- read_csv(file.path(ROOT, "results", "16_Orthology", "v7_single_cohort",
                           "candidates_GSE58485_only.csv"),
                 col_types = cols(.default = col_character()))
cand$tier <- toupper(cand$tier)
tab_cand <- data.frame(
  tier = c("A", "B", "C"),
  n = c(sum(cand$tier == "A"), sum(cand$tier == "B"), sum(cand$tier == "C")))
tab_cand$tier <- factor(tab_cand$tier, levels = c("A", "B", "C"))
stopifnot(nrow(cand) == 11, all(cand$tier == "C"))

eval_upper <- toupper(unique(read_csv(file.path(ROOT, "results", "16_Orthology",
                                                "evaluable_targets_v6_CLEAN.csv"),
                                      col_types = cols(.default = col_character()))$gene))
deg_upper <- toupper(unique(read_csv(file.path(ROOT, "results", "16_Orthology",
                                               "v7_single_cohort",
                                               "deg_significant_GSE58485_only.csv"))$gene))
key_upper <- toupper(unique(read_csv(file.path(ROOT, "results", "16_Orthology",
                                               "v7_single_cohort",
                                               "key_genes_GSE58485_only.csv"))$gene))
stopifnot(length(eval_upper) == 69, length(deg_upper) == 373, length(key_upper) == 364)

# ---------- 面板 A / B ----------
p2a <- ggplot(tab_tier, aes(x = tier, y = n, fill = tier)) +
  geom_col(width = 0.65, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = n), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c(A = "#B2182B", B = "#F4A582", C = "#4393C3"),
                    guide = "none") +
  labs(x = "Evidence tier", y = "Predicted BaP targets (deduplicated)",
       title = "BaP target evidence tiers") +
  theme_pub(base_size = 11) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

p2b <- ggplot(tab_cand, aes(x = tier, y = n, fill = tier)) +
  geom_col(width = 0.65, colour = "grey30", linewidth = 0.3) +
  geom_text(aes(label = n), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c(A = "#B2182B", B = "#F4A582", C = "#4393C3"),
                    guide = "none") +
  labs(x = "Evidence tier", y = "BaP \u00d7 TBI candidate genes",
       title = "Candidate tiers (n = 11; all tier C)") +
  theme_pub(base_size = 11) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)))

# ---------- 面板 C：三集合韦恩图 ----------
venn_list <- list(`Evaluable BaP targets` = eval_upper,
                  DEGs = deg_upper,
                  `Module key genes` = key_upper)
p2c <- ggvenn(venn_list,
              fill_color = c("#9ECAE1", "#A1D99B", "#FDBE85"),
              stroke_size = 0.4, text_size = 3.5, set_name_size = 4.2) +
  ggtitle("BaP targets \u2229 DEG \u2229 TBI-associated modules") +
  theme_pub(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, size = 11))

# ---------- 合成 ----------
p2 <- (p2a | p2b) / p2c + plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 11, face = "bold"))
ggsave_pub(file.path(ROOT, "figures", "Figure_02_BaP_Targets", "Figure2.pdf"),
           p2, width = 10, height = 8, bg = "white")
writeLines(c("Figure2.pdf regenerated",
             paste0("targets: ", paste(tab_tier$n, collapse = "/")),
             paste0("candidates: ", paste(tab_cand$n, collapse = "/")),
             paste0("venn sets: 69/373/364; intersection = ",
                    length(intersect(intersect(eval_upper, deg_upper), key_upper)))),
           file.path(ROOT, "figures", "Figure_02_BaP_Targets", "figure2_v7.log"))
