# Master simulation script: runs all simulation scripts
source(fs::path("analysis", "simulation", "simulation_1_calibration.R"))
source(fs::path("analysis", "simulation", "simulation_2_power.R"))
source(fs::path("analysis", "simulation", "simulation_3_calibration_varEstim.R"))
source(fs::path("analysis", "simulation", "simulation_4_calibration_sampleSize.R"))
source(fs::path("analysis", "simulation", "simulation_5_power_sampleSize.R"))
source(fs::path("analysis", "simulation", "simulation_6_calibration_FE_RE.R"))
source(fs::path("analysis", "simulation", "simulation_nonparametric_1_nstudies.R"))
source(fs::path("analysis", "simulation", "simulation_nonparametric_2_nsamples.R"))
