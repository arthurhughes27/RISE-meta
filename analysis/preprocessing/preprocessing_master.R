# Master preprocessing script: runs all preprocessing steps in order
# Each script reads from and writes to the data/ directory

source("analysis/preprocessing/preprocessing_clinical.R")
source("analysis/preprocessing/preprocessing_expression.R")
source("analysis/preprocessing/preprocessing_immuneresponse.R")
source("analysis/preprocessing/preprocessing_GSA.R")
source("analysis/preprocessing/preprocessing_merging.R")
