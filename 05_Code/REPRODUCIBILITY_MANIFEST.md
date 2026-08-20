# Reproducibility Manifest 

## Script order

01→02→02b→02c→03→04→04b→06→06b→07→08→09→10→11→12→13→14→15→16→17→18→19→20→21→22→23→25→26→27→28→29→30→31→32→33→34→35→36→40→41→52→53→54→58→59→60→61→61b→62→63→64→65→66→67→68→70→71→72→73→74→75→76→77→78

Numbered scripts are in this folder (`scripts/*.R`, `scripts/*.py`); helper
`theme_pub.R` and install scripts are prefixed `install_*`. Scripts 01–68 are
the pre- and stress-test pipeline (merged-cohort analyses are retained as
stress tests); scripts 70–78 are the the primary and submission pipeline:

- `70_single_cohort_main.R` — GSE58485-only primary analysis (limma,
  WGCNA, candidate intersection; seed-independent except WGCNA soft-threshold
  selection).
- `71_permutation.R` — measurability-matched permutation null
  (seed=2026, 10,000 draws).
- `72_ml_benchmark.R` — descriptive 11-algorithm benchmark (seed=2026,
  5 repeats × 5-fold CV; PROBAST+AI high risk; supplementary only).
- `73_validation.R` — direction replication (GSE128543/GSE205958) and
  human ortholog coverage.
- `74_go_kegg.R` — GO/KEGG enrichment for the 11 candidates.
- `75_figures.R` — main-text figures (volcano, module–trait, permutation,
  sensitivity heatmap, GO bubble).
- `76_compound_control.py` — ChEMBL single-source pyrene/phenanthrene control
  (S12).
- `77_supplementary_figures.py` — S1a algorithm benchmark and S2 direction
  replication figures.
- `78_renumber_references.py` — BMC first-appearance reference renumbering
  (mapping in `reports/reference_renumber_map.csv`).

## Environment

- R 4.5.3 (Windows), Bioconductor packages: sva, org.Mm.eg.db, scDblFinder,
  harmony, Seurat, CellChat, scTenifoldKnk, clusterProfiler, pROC, etc.
- Python 3.12: matplotlib/PIL (figure rendering), python-docx (document
  assembly), openpyxl (supplementary workbook).
- Seeds: permutation 2026 (10,000 draws); ML 2026; scDblFinder 2026;
  docking seeds 2026/42/8888.
- Full `sessionInfo()` outputs are in the source tree
  `scripts/install/INSTALL_REPORT.txt` and project reports.

## Data

All raw data are public GEO series: GSE58485, GSE41345, GSE128543,
GSE205958, GSE104687, GSE209552, GSE75206 (see `datasets.tsv`).
GSE2871 and GSE160763 were downloaded but excluded (see manuscript Methods
4.1 and supplementary methods S6).

## Figure generation

- R/Python figures are exported as vector PDFs (`theme_pub.R` /
  `ggsave_pub`).
- Figure 1 workflow: `58_figure1_workflow.R` renders
  `Figure1_workflow_v5.png/.pdf`; editable source `Figure1_workflow_v5.drawio`.
-  main figures: `75_figures.R`; supplementary figures S1a/S2:
  `77_supplementary_figures.py`.
- Manuscript DOCX: `59_build_manuscripts_docx.py` (EN/ZH masters); the
  BMC-formatted manuscript is generated from
  `submission/01_Manuscript/MANUSCRIPT_SR_MAIN.md` with pandoc.
- Submission assembly: `scripts/_assemble_submission.ps1`;
  supplementary workbook: `68_build_supplementary.py` ( sheets included).
