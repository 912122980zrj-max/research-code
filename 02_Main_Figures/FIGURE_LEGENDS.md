# 图题与图件清单（ 定稿口径）

> 更新日期：2026-08-13。数字与 `reports/manuscript/MANUSCRIPT_FINAL__*.md` 一致。
>  依据 2026-08-11 评估意见重构：主分析 = GSE58485-only（单均质队列），
> 合并分析 = 压力测试；叙事定位 = 网络毒理学靶点–疾病交集方法基准检验。
> 主图 4 张（每图一个合成 PDF），补充图 S1–S5，补充表 S1–S14。
> 单独上传用 TIFF（300 DPI、Deflate 无损压缩）：Figure1.tiff–Figure4.tiff。

## 主图

**Figure 1. Study design / 研究设计。** 预测 BaP 靶点并按证据分级（A/B/C），Ensembl Compara 1:1 直系同源映射；主分析使用单一均质发现队列（GSE58485，皮层 CCI 3 d），经覆盖过滤、差异表达、WGCNA 与候选交集后，将观测交集与可测量性匹配置换零假设比较；合并队列与替代校正路径作为压力测试；候选方向在两个同型队列中检验；单细胞、情境与正对照分析为探索性。文件：`Fig1_StudyDesign.png`（单图全宽）+ `Fig1_StudyDesign.pdf`（单页矢量版，由 PNG 生成）；可编辑源文件 `Fig1_StudyDesign.drawio`。

**Figure 2. Predicted BaP targets and candidate funnel / BaP 预测靶点与候选漏斗。**（A）364 个独特靶点的证据分级（A=5、B=23、C=336；Tier C 占 92.3%）。（B）11 个候选基因的分级分布（全部 Tier C；A/B=0）。（C）三集合韦恩图：可评估（直系同源）BaP 靶点（69）、GSE58485-only DEG（373）、WGCNA 模块关键基因（364），交集 n=11。文件：`Fig2.pdf`（单张合成 PDF，A/B/C 三面板；由 `82_figure2.R` 生成）。

**Figure 3. Primary single-cohort analysis / 主分析（单队列）。**（A）373 个 DEG 火山图（上调 295、下调 78；11 个候选以橙色高亮；|log2FC|>0.585，BH adjP<0.05）。（B）WGCNA 模块–性状相关（green r=0.978，P=1.7×10⁻¹⁶）。（C）11 个候选基因的发现队列 log2FC（全部 Tier C）。（D）置换背景检验：10,000 次抽样的重叠计数零分布（seed=2026；从 5,043 个可测基因抽 11 个 vs 364 个关键基因），观测值 11 以红线标注（期望 0.79；13.9 倍；经验 P<0.0001）。文件：`Fig3.pdf`（单张合成 PDF，A–D 四面板）。

**Figure 4. Enrichment and sensitivity summary / 富集与敏感性汇总。**（A）11 个候选基因 GO 生物过程富集前 10 条（Count≥3，BH<0.05）。（B）AhR 反应基因与 Tier A 靶点未被芯片覆盖（Fisher 检验不可评估，不报告 P 值）。（C）各路径候选集合（行：GSE58485-only、green 单模块、合并 CLEAN、FULL、removeBatchEffect、SVA、未校正），按存在与否着色；主分析行高亮。文件：`Fig4.pdf`（单张合成 PDF，A/B/C 三面板）。

## 补充图 S1–S5（提交包文件位于 `03_Supplementary/Figures/`）

| 编号 | 内容 | 文件 |
|---|---|---|
| S1 | 11 种算法 CV AUC 基准；ENet OOF ROC/校准/DCA（，11 基因） | `S1a_algorithm_benchmark.pdf`；`S1b_oof_roc.pdf`；`S1c_calibration_DCA.pdf` |
| S2 | 方向复制：逐基因验证 log2FC vs 发现（GSE128543、GSE205958） | `S2_direction_replication.pdf` |
| S3 | snRNA-seq QC、注释、组成与标记热图 | `S3a_QC_violin.pdf`；`S3b_marker_dotplot.pdf`；`S3c_marker_heatmap.pdf`；`S3d_sample_composition.pdf` |
| S4 | 去卷积完整结果（MuSiC/BisqueRNA/CIBERSORTx + LOO） | `S4a_MuSiC_proportions.pdf`；`S4b_BisqueRNA_proportions.pdf`；`S4c_CIBERSORTx_proportions.pdf`；`S4d_LOO_validation.pdf` |
| S5 | CellChat、虚拟敲除与分子对接 | `S5a_TGFb_network_Sham.pdf`；`S5b_TGFb_network_TBI.pdf`；`S5c_KO_volcano_P2RY6.pdf`；`S5d_KO_volcano_P2RX7.pdf`；`S5e_docking_energy_summary.pdf`；`S5f–S5o_pose_*.png`（10 个受体） |

## 补充表 S1–S14（提交包文件位于 `03_Supplementary/Tables/`）

> 工作簿：`Supplementary_Tables_S1-S12.xlsx`；原始 CSV 位于 `Raw_CSV/`，
> 与工作簿各页一一对应（详见 `03_Supplementary/Text/SUPPLEMENTARY_LEGENDS_EN.md`）。

| 编号 | 内容 | 文件 |
|---|---|---|
| S1 | 直系同源映射（364 行；Ensembl release；1:1 vs 1:多；Hprt 探针说明） | `results/16_Orthology/orthology_mapping_full.csv` |
| S2 | 置换检验结果（ 主分析 + v6 合并参考） | `results/16_Orthology/single_cohort/permutation_.csv`；`permutation_results.csv` |
| S3 | 敏感性矩阵 + 阈值扫描 + 模块口径与基因型/未校正敏感性 | `results/16_Orthology/single_cohort/`；`sensitivity_matrix.csv` |
| S4 | ML 基准（，11 基因）；ENet OOF/校准；单队列 DEG 重叠 | `results/07_ML_Signature/D3_ml_benchmark/` |
| S5 | GO/KEGG 富集（11 候选）；AhR 覆盖检查 | `results/06_GO_KEGG_Enrichment/single_cohort/` |
| S6 | Welch 假设检查；逐基因方向复制；既往研究 [11] hub 基因 DEG 统计 | `results/08_External_Validation/` |
| S7 | snRNA-seq 组成统计 | `results/11_SingleCell/` |
| S8 | 对接得分与对照 | `results/10_Molecular_Docking/` |
| S9 | 星形亚群统计 | `results/14_Astrocyte_Subclusters/D8_astro/` |
| S10 | 正对照检索（GSE75206 结果 + 候选系列） | `results/16_Orthology/positive_control_GSE75206/`；`geo_positive_control_search.csv` |
| S11 | 去卷积方法一致性矩阵 | `results/09_Immune_Deconvolution/deconvolution_method_corr.csv` |
| S12 | 非 BaP 化合物对照（ChEMBL 单源；芘/菲；0 个可评估靶点） | `results/16_Orthology/compound_control/compound_control_results.csv` |
| S13 | 样本级清单（全部 47 个样本：队列、基因型、处理、分组、纳入标记与排除原因） | `Raw_CSV/S13_sample_manifest.csv` |
| S14 | 完整平台（GPL1261，27,250 基因）敏感性分析：汇总、候选、AhR 电池、置换与方向一致性 | `Raw_CSV/S14_full_universe_*.csv` |
