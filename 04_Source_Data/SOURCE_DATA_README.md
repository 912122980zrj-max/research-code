# Source Data (v7)

This folder contains the processed data underlying the v7 manuscript
(`01_Manuscript/MANUSCRIPT_v7_BMC_EN.md`). The primary analysis is
GSE58485-only (single homogeneous cohort); merged-cohort files are the stress
test and are labeled CLEAN/merged.

## Primary analysis (GSE58485-only, n=24)

- `expression_GSE58485_only_primary.csv` — primary analysis expression matrix,
  5,043 genes × 24 samples (uncorrected; identical to the matrix produced by
  `70_single_cohort_main_v7.R` prep(); subset to the GSE58485 CLEAN samples).
- `corrected_phenotype_GSE58485_only.csv` — sample metadata for the 24
  primary samples (group, genotype coding).
- `deg_significant_GSE58485_only.csv` — limma results (373 DEGs; 295 up,
  78 down; |log2FC|>0.585, BH adjP<0.05; design ~group+genotype).
- `key_genes_GSE58485_only.csv` — 364 key genes (DEG ∩ significant WGCNA
  modules).
- `candidates_GSE58485_only.csv` — 11 tier-C candidates (key genes ∩ 69
  orthology-evaluable BaP targets).
- `candidate_stats_GSE58485_only.csv` — log2FC/P/adjP for the 11 candidates.
- `validation_direction_v7.csv` — per-gene direction replication
  (GSE128543 11/11; GSE205958 10/11).

## Merged-cohort stress test (CLEAN n=33 / FULL n=42)

- `corrected_expression_CLEAN.csv` — ComBat-corrected merged matrix used by the
  merged stress-test path, 5,043 genes × 33 samples (24 GSE58485 + 9 GSE41345);
  note this is NOT the primary matrix (the primary path uses the uncorrected
  matrix above).
- `corrected_phenotype_CLEAN.csv` — sample metadata (cohort, group,
  genotype coding).
- `sample_composition.csv` — discovery cohort composition (CLEAN/FULL flags).
- `deg_all.csv` / `deg_significant_CLEAN.csv` — merged-path limma results
  (198 DEGs under ComBat).
- `module_genes.csv`, `module_trait.csv`, `tbi_key_genes.csv` — merged-path
  WGCNA outputs (black r=0.827; 193 key genes).
- `candidate_genes_tiered_v2.csv` — merged-path candidates (7, all tier C).
- `evaluable_targets_v6_CLEAN.csv` — orthology-evaluable target list used
  by both paths.

## Single-cell and contextual analyses (exploratory)

- `singlecell_per_sample_counts.csv`, `singlecell_roe.csv`,
  `singlecell_cluster_markers.csv` — GSE209552 composition, Ro/e, and markers.
- `deconvolution_group_stats.csv`, `pseudobulk_loo_stats.csv` — deconvolution
  group statistics.
- `docking_best_poses.csv`, `docking_nearest_residues.csv` — AutoDock Vina
  best poses and interacting residues (10-receptor panel).

All files are also available in the project `results/` tree. Raw GEO data are
public; accession list in Table 1 of the manuscript. v6-era files (SHAP,
frozen-signature validation, per-gene direction validation) were moved to
`reports/_archive_docs/manuscript_history/submission_legacy_20260811/`.
