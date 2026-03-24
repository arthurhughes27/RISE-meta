# Master preprocessing script: runs all preprocessing steps in order
# Each script reads from and writes to the data/ directory

source("analysis/preprocessing/preprocessing_clinical.R")
source("analysis/preprocessing/preprocessing_expression.R")
source("analysis/preprocessing/preprocessing_immuneresponse.R")
source("analysis/preprocessing/preprocessing_BTM.R")
source("analysis/preprocessing/preprocessing_BG3M.R")
source("analysis/preprocessing/preprocessing_merging.R")
