# Master script to run all preprocessing, descriptive, data analysis and simulation scripts

# Script to install any missing dependencies
source(fs::path("analysis", "install_dependencies.R"))

# Script to ensure project structure is in place
source(fs::path("analysis", "create_project_structure.R"))

# Master script for preprocessing files
source(fs::path("analysis", "preprocessing", "preprocessing_master.R"))

# Master script for data description files
source(fs::path("analysis", "descriptive", "descriptive_master.R"))

# Master script for data analysis files
source(fs::path("analysis", "application", "application_master.R"))

# Master script for simulation files
source(fs::path("analysis", "simulation", "simulation_master.R"))