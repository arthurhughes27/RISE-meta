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

4. Install required packages:

```r
source("analysis/install_dependencies.R")
```

5. Run the full pipeline:

```r
source("analysis/analysis_master.R")
```

---

## Requirements

The following R packages are required:

- Biobase
- digest
- dplyr
- forcats
- fs
- ggplot2
- GSA
- janitor
- kableExtra
- knitr
- patchwork
- purrr
- readr
- readxl
- scales
- stringr
- Surrogate
- SurrogateRank
- tibble
- tidyr
- tidyverse
- viridis

Install them with:

```r
source("analysis/install_dependencies.R")
```

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