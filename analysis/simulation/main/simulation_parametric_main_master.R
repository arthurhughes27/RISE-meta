# Master simulation script: runs all main parametric simulation scripts
source(fs::path("analysis", "simulation", "main", "simulation_parametric_calibration_main.R"))
source(fs::path("analysis", "simulation", "main", "simulation_parametric_power_main.R"))
source(fs::path("analysis", "simulation", "main", "simulation_parametric_calibration_invalidStrength.R"))
source(fs::path("analysis", "simulation", "main", "simulation_parametric_power_validStrength.R"))