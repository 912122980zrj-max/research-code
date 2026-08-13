# research-code — BaP–TBI network toxicology benchmarking (v1.1.0)

Analysis code for the manuscript:

> Quantifying the target–disease intersection in network toxicology: a
> transcriptomic stress test of predicted benzo[a]pyrene targets in
> traumatic brain injury (submitted to *BMC Bioinformatics*).

## Requirements

- R ≥ 4.5 (developed on R 4.5.3, Windows) with Bioconductor packages:
  `sva`, `org.Mm.eg.db`, `scDblFinder`, `harmony`, `Seurat`, `CellChat`,
  `scTenifoldKnk`, `clusterProfiler`, `pROC`; CRAN packages: `WGCNA`,
  `limma`, `caret`, `glmnet`, `xgboost`, `dplyr`, `readr`, `ggplot2`,
  `ggvenn`, `patchwork`, `openxlsx`.
- Python ≥ 3.9 with `requests`, `pandas`, `numpy`, `matplotlib`,
  `Pillow`, `python-docx`, `openpyxl`, `pymupdf` (see `requirements.txt`).

## Installation

### R

Run the provided installer (installs missing CRAN/Bioconductor packages):

```bat
install_install_missing_packages.bat
```

or from R:

```r
source("install_install_missing_packages.R")
```

Check the installation:

```r
source("install_00_check_installed.R")
source("install_01_check_availability.R")
source("install_02_smoke_test.R")
```

### Python

```bash
pip install -r requirements.txt
```

## Data

All raw data are public GEO series (see `datasets.tsv`): GSE58485 (primary),
GSE41345 (merged stress test), GSE128543 and GSE205958 (direction
consistency), GSE104687 (human transportability), GSE209552 (single-cell),
GSE75206 (positive control). Download with:

```r
source("01_download_geo.R")        # primary/merged cohorts
source("01b_download_validation.R") # validation cohorts
```

Raw GEO files are not included in this archive; processed matrices are in
the submission package (`04_Source_Data/`).

## Run order

v7 primary/submission pipeline (recommended):

```
70_single_cohort_main_v7.R → 71_permutation_v7.R → 72_ml_benchmark_v7.R →
73_validation_v7.R → 74_go_kegg_v7.R → 75_figures_v7.R →
76_compound_control.py → 77_supplementary_figures_v7.py →
82_figure2_v7.R → 83_composite_figures.py → 78_renumber_references.py
```

The full script order (01→…→78) and version notes are in
`REPRODUCIBILITY_MANIFEST.md`.

## Reproducibility

- Random sampling uses fixed seeds (permutation and cross-validation:
  seed = 2026; scDblFinder and docking seeds recorded in the manifest).
- DEG and module definitions were fixed before all permutation draws.
- Software versions, output paths, and timestamps are recorded in
  `REPRODUCIBILITY_MANIFEST.md`.

## License

Creative Commons Attribution 4.0 International (CC BY 4.0).
