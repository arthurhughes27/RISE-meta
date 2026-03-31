# Master simulation script: runs all main parametric simulation scripts
main_folder_path = fs::path("analysis", "simulation", "main")

source(fs::path(main_folder_path, "simulation_parametric_calibration_main.R"))
source(fs::path(main_folder_path, "simulation_parametric_power_main.R"))
source(fs::path(main_folder_path, "simulation_parametric_calibration_invalidStrength.R"))
source(fs::path(main_folder_path, "simulation_parametric_power_validStrength.R"))