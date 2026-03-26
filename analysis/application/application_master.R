# Master script to run all application scripts
# Main application script
source("analysis/application/application_main.R")

# Classic applications to low dimensional data
source("analysis/application/application_ARMD.R")
source("analysis/application/application_ovarian.R")

# Supplementary scripts for sensitivity analysis
source("analysis/application/application_supplementary_MetaAnalysisSpec.R")
source("analysis/application/application_supplementary_epsilon.R")
source("analysis/application/application_supplementary_geneset.R")

