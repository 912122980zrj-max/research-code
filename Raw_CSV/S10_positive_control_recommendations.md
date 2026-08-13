# B-6 正对照检索与执行结果

- 更新：2026-08-11 07:30（网络窗口内完成检索与初步分析；未超止损）
- 检索：NCBI GEO E-utilities（`benzo[a]pyrene` + `Mus musculus` + `gse[ETYP]`，三组检索式）
- 检索结果：52 个唯一系列；脑相关 1 个（`results/16_Orthology/geo_positive_control_search.csv` 全表）

## 已验证正对照（主推荐）

| 项目 | 内容 |
|---|---|
| 系列 | **GSE75206**（GDS 200075206；GPL10787 Agilent-028005 SurePrint G3 Mouse GE 8x60K） |
| 设计 | 成年雄性 Muta Mouse，海马；BaP 70 mg/kg/天 × 3 天口服 vs 橄榄油对照；24 h 后取材 |
| 分组 | Control n=4（GSM1945106/8/11/13）、BaP n=4（GSM1945107/9/10/12）——每组 ≥3 例 ✔ |
| 处理数据 | `GSE75206_series_matrix.txt.gz`（探针水平，LOWESS 归一化 log2 比值） |
| 探针注释 | Ensembl BioMart `agilent_sureprint_g3_8x60k`（42,107 探针 → 30,456 基因） |
| 文献 | PubMed 27195522（Chepelev et al.，2016；原文结论：海马未见经典 BaP 靶点（异生物代谢）效应） |

## GSE75206 正对照分析结果（67_positive_control_GSE75206.py）

- 与主分析同一映射/检验思路（基因水平中位数折叠 → Welch t 检验 → BH）。
- 全部基因 BH<0.05：**0 个**（n=4 vs 4，功效低）；69 个可评估靶点中 68 个可检测（Hprt 以 Ensembl 符号 Hprt1 存在），**BaP–对照 DEG ∩ 69 靶点 = 0**。
- 经典 BaP/AhR 电池（本数据集内）：Cyp1a2（FC −0.15，原始 P=0.040，BH=0.769）、Ahr（FC −0.06，P=0.087）、Nqo1（FC +0.34，P=0.098）、Cyp1a1（FC −0.08，P=0.450）——均未达到 BH 显著，方向也不一致；与原文“海马中经典 BaP 靶点无效应”一致。
- 已知阳性对照（原文 RT-PCR 验证）：Grina（FC +0.91）、Grin2a（FC +0.40）方向与原文一致（本芯片检验未达显著）。
- 解读：直接脑暴露的正对照也未显示 BaP 靶点–脑转录组强桥接；该结果支持主研究“交集微弱、AhR 轴在脑中不可见”的结论，同时如实反映正对照功效有限（n=4 vs 4）与组织特异性（海马）。

## 次级候选（未逐一核实分组，使用前需确认）

见 `geo_positive_control_search.csv`（52 行）。示例（n_samples 来自 GEO 摘要，分组需再核实）：

| 系列线索 | 主题 | n_samples |
|---|---|---|
| 200337058 | 小鼠 MOVAS 细胞 BaP 暴露 circRNA | 6 |
| 200230356 / 200213947 / 200213946 / 200213945 | 产前 BaP 暴露卵巢生殖细胞 | 6–12 |
| 200317094–200317185 | 野火烟尘颗粒（含 BaP 途径）卵巢雄激素化 | 6–11 |

## 止损与兜底

- 检索+分析耗时：约 1.5 小时（< 4 小时检索止损 + 半天跑数止损）✔
- CTD 弱替代未启用；如后续需要更大功效，可再评估肺/肝 BaP 暴露系列（本检索 CSV 内）。
