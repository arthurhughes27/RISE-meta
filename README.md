# RISE-Meta
Repository for the development and implementation of the **Rank-Based Identification of High-Dimensional Surrogate Markers Meta-Analysis (RISE-Meta)** method.
The methodology is described in the paper *Meta-Analysis of High-Dimensional Surrogate Markers*, available soon on arXiv.

---

## Quick start

1. Clone this repository.
2. Download the required raw data (described below) and place it in `data-raw/`.
3. Create the required project structure:
```r
source("analysis/create_project_structure.R")
```
4. Install dependencies and recreate software environment
```r
renv::restore()
```
5. Run the full pipeline:
```r
source("analysis/analysis_master.R")
```

---

### Reproducing the software environment

This repository uses [`renv`](https://rstudio.github.io/renv/) for package management. To restore the exact package versions used in this analysis:

```r
renv::restore()
```

The complete package environment is recorded in `renv.lock`.

---

## Data access

No raw or processed data is stored in this repository.

To run the project:
1. Download required raw data (see below)
2. Place all files in `data-raw/`
3. Run scripts in `analysis/`

---

### Low-dimensional application

Data is accessed via the **Surrogate** R package:
```r
install.packages("Surrogate")
```

---

### High-dimensional application

#### ImmuneSpace data

1. Go to https://immunespace.org
2. Click **Resources → Immune Signature Projects**
3. Open **Project 2** (Immune Signatures 2)
4. Click **Proceed to study data** (account required)
5. Open **Clinical and Assay Data**

Download:
- `all_noNorm_eset.rds`

From the **Assay** section download:
- Neutralizing antibody titer
- Hemagglutination inhibition (HAI)

Export each as:
1. Click dataset
2. Click **Export**
3. Select **Excel Workbook (.xlsx)**
4. Download file

Expected filenames:
- `neut_ab_titer_2025-01-10_01-13-22.xlsx`
- `hai_2025-01-10_01-13-41.xlsx`

---

### Gene set data

- **Blood Transcriptional Modules (BTMs)**
  https://github.com/shuzhao-li/BTM/tree/master/BTM/datasets
  File: `BTM_for_GSEA_20131008.gmt`

- **BloodGen3 Modules (BG3M)**
  DOI: 10.1093/bioinformatics/btab121
  File: `Suppl_File_1_BIOINF.xls`

---

## Project structure

- `data-raw/`
  Raw data files (not stored in this repository)

- `data/`
  Processed data (not stored in this repository)

- `R/`
  Helper functions

- `analysis/`
  Analysis scripts:
  - preprocessing
  - descriptive
  - application
  - simulation

- `output/`
  Outputs (not stored in this repository):
  - `results/`
  - `figures/`

- `manuscript/`
  Manuscript-related files:
  - `diagrams/` manually created figures

---

## Analysis

### Full pipeline

```r
source("analysis/analysis_master.R")
```

Runs:
- package installation
- project folder structure
- preprocessing
- descriptive analysis
- real-data application
- simulation study (may take several days)

---

### Preprocessing

```r
source("analysis/preprocessing/preprocessing_master.R")
```

Steps:
- clinical data preprocessing
- gene expression preprocessing
- immune response preprocessing
- gene set preprocessing
- dataset merging

Outputs of preprocessing go in `data/`

---

### Descriptive analysis

```r
source("analysis/descriptive/descriptive_master.R")
```

Outputs:
- visual description of dataset
- data for flowchart of samples after each preprocessing step

Saved in:
- `output/results/descriptive/`
- `output/figures/descriptive/`

---

### Application

```r
source("analysis/application/application_master.R")
```

Saved in:
- `output/figures/application/`

---

### Simulation

```r
source("analysis/simulation/simulation_master.R")
```

Saved in:
- `output/results/simulation/`
- `output/figures/simulation/`

---

## Citation

The methodology is described in:

*Meta-Analysis of High-Dimensional Surrogate Markers*

and implemented in the R package SurrogateRank.

(arXiv preprint coming soon)

--- 

## Software and version information

Below is a summary of the output of 

```{r}
sessionInfo()
```

after loading all the libraries necessary for this repository. 

## Software and version information

- **OS:** Ubuntu 24.04.4 LTS (x86_64-pc-linux-gnu)
- **R version:** 4.5.3 (2026-03-11)
- **BLAS/LAPACK:** OpenBLAS (LAPACK version 3.12.0)
- **Time zone:** Europe/Paris

### Attached packages

```
knitr_1.51          kableExtra_1.4.0    digest_0.6.39       viridis_0.6.5
viridisLite_0.4.3   scales_1.4.0        janitor_2.2.1       Biobase_2.70.0
BiocGenerics_0.56.0 generics_0.1.4      readxl_1.4.5        GSA_1.03.3
renv_1.2.2          Surrogate_3.4.1     patchwork_1.3.2     boot_1.3-32
lubridate_1.9.5     forcats_1.0.1       stringr_1.6.0       purrr_1.2.2
readr_2.2.0         tidyr_1.3.2         tibble_3.3.1        ggplot2_4.0.3
tidyverse_2.0.0     SurrogateRank_3.0   dplyr_1.2.1         fs_2.1.0
```

### Packages loaded via namespace

```
gtable_0.3.6        xfun_0.57           lattice_0.22-9      tzdb_0.5.0
vctrs_0.7.3         tools_4.5.3         Rdpack_2.6.6        parallel_4.5.3
sandwich_3.1-1      pbmcapply_1.5.1     pkgconfig_2.0.3     Matrix_1.7-5
RColorBrewer_1.1-3  S7_0.2.2            lifecycle_1.0.5     compiler_4.5.3
farver_2.1.2        textshaping_1.0.5   maxLik_1.5-2.2      snakecase_0.11.1
htmltools_0.5.9     pillar_1.11.1       nloptr_2.2.1        MASS_7.3-65
reformulas_0.4.4    nlme_3.1-169        tidyselect_1.2.1    stringi_1.8.7
splines_4.5.3       miscTools_0.6-30    fastmap_1.2.0       grid_4.5.3
cli_3.6.6           magrittr_2.0.5      survival_3.8-6      withr_3.0.2
timechange_0.4.0    rmarkdown_2.31      lme4_2.0-1          gridExtra_2.3
cellranger_1.1.0    zoo_1.8-15          hms_1.1.4           evaluate_1.0.5
rbibutils_2.4.1     rlang_1.2.0         Rcpp_1.1.1          glue_1.8.1
xml2_1.5.2          BiocManager_1.30.27 svglite_2.2.2       rstudioapi_0.18.0
minqa_1.2.8         R6_2.6.1            systemfonts_1.3.2
```

> Full version details are recorded in `renv.lock`. Run `renv::restore()` to reproduce the exact package environment.