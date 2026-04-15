# RISE-Meta

Repository for the development and implementation of the **Rank-Based Identification of High-Dimensional Surrogate Markers Meta-Analysis (RISE-Meta)** method.

The methodology is described in the paper "Meta-Analysis of High-Dimensional Surrogate Markers", available soon on arxiv. 

## Project structure

- `data-raw/`  
  Raw data files downloaded by the user (not stored in this repository)

- `data/`  
  Processed data generated during preprocessing

- `R/`  
  Helper functions used throughout the project

- `analysis/`  
  Scripts used to generate results and figures for the manuscript  
  Subfolders:
  - `preprocessing/`
  - `descriptive/`
  - `application/`
  - `simulation/`

- `output/`  
  Project outputs:
  - `results/` (raw results)
  - `figures/` (manuscript figures)

## Data access

No raw or processed data is stored in this repository.

To run the project:
1. Download the required raw data (see below)
2. Place all files in `data-raw/`
3. Run the scripts in `analysis/`

---

## Low-dimensional application

Data is accessed via the **Surrogate** R package:

```r
install.packages("Surrogate")
```

---

## High-dimensional application

### ImmuneSpace data

1. Go to: https://immunespace.org  
2. Click **Resources → Immune Signature Projects**  
3. Scroll to **Project 2** and open the Immune Signatures 2 page  
4. Click **Proceed to study data** (account required)  
5. Open **Clinical and Assay Data**

Download:
- `all_noNorm_eset.rds`

Also download assay data from the **Assay** section:
- Neutralizing antibody titer 
- Hemagglutination inhibition (HAI)

For each:
1. Click the dataset  
2. Click **Export**  
3. Select **Excel Workbook (.xlsx)**  
4. Download the file

This should download files named:
- `neut_ab_titer_2025-01-10_01-13-22.xlsx`
- `hai_2025-01-10_01-13-41.xlsx`

---

### Gene set data

- **Blood Transcriptional Modules (BTMs)**  
  Available from the author's github
  https://github.com/shuzhao-li/BTM/tree/master/BTM/datasets  
  File: `BTM_for_GSEA_20131008.gmt`

- **BloodGen3 Modules (BG3M)**  
  Available from the supplementary files in the original publication
  DOI: 10.1093/bioinformatics/btab121  
  File: `Suppl_File_1_BIOINF.xls`

---

## Analysis

Run the full pipeline:

```r
source("analysis/analysis_master.R")
```

This executes all components:
- installation of required packages
- preprocessing
- descriptive analysis
- real-data application
- simulation study (note: this can run for up to a few days on a standard pc)

You can also run each component separately.

---

## Preprocessing

```r
source("analysis/preprocessing/preprocessing_master.R")
```

Steps:
- preprocess clinical data (study IDs, vaccine info, variable names)
- preprocess gene expression data
- preprocess immune response data
- preprocess gene set data
- merge datasets

Outputs are written to `data/`.

---

## Descriptive analysis

```r
source("analysis/descriptive/descriptive_master.R")
```

Outputs:
- study descriptions
- sample flow summaries

Saved in:
- `output/results/descriptive/`
- `output/figures/descriptive/`

---

## Application

```r
source("analysis/application/application_master.R")
```

Runs all real-data analyses.

Saved in:
- `output/figures/application/`

---

## Simulation

```r
source("analysis/simulation/simulation_master.R")
```

Runs simulation study analyses.

Saved in:
- `output/results/simulation/`
- `output/figures/simulation/`

---

## Note

Simulation may take several hours to multiple days on a standard computer.

# Dependencies 

The following R packages are required to run the files in this repository, and will be installed upon running

```r
source("analysis/install_dependencies.R")
```

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