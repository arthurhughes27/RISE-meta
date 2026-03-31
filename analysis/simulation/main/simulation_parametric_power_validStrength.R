# Script to derive results for parametric simulation of trial-level surrogate effects
# In this script, we generate only valid surrogates, in order to examine the false positive rate 

# Libraries
library(SurrogateRank)
library(tidyverse)
library(ggplot2)
library(viridis)
library(dplyr)
library(tibble)
library(scales)

# Folder to store results and figures
simulation_figures_folder = fs::path("output", "figures", "simulation", "main")
simulation_results_folder = fs::path("output", "results", "simulation", "main")

# Number of independent markers to generate trial-level effects for
J = 10000

# Value of epsilon defining the validity region
epsilon <- 0.1

# Grid of number of studies to generate
M_grid <- seq(5, 25, 5)

# Grid of number of samples within each study to generate
nm_grid = c(50)

# Grid of maximum between-study variability values
u_tau_max_vals <- c(epsilon/100, epsilon/10, epsilon, epsilon*10, epsilon*100)

# Grid of maximum within-study variability values
u_nu_max_vals <- c(epsilon/100, epsilon/10, epsilon, epsilon*10, epsilon*100)

# Grid of alpha values for the calibration plot
alpha_grid <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5)

# Grid of values for valid surrogate mean generation
mu_valid_vals <- seq(-epsilon, epsilon, epsilon/10)

# Initialise data frame WITHOUT alpha
results <- expand.grid(
  M = M_grid,
  nm = nm_grid,
  u_tau_max = u_tau_max_vals,
  u_nu_max = u_nu_max_vals,
  mu_valid = mu_valid_vals,
  stringsAsFactors = FALSE
)

# Initialise list for storage
results_list <- vector("list", nrow(results))

for (i in seq_len(nrow(results))) {
  
  M <- results$M[i]
  sample_sizes <- rep(results$nm[i], M)
  u_tau_max <- results$u_tau_max[i]
  u_nu_max <- results$u_nu_max[i]
  mu_valid = results$mu_valid[i]
  
  # --- DATA GENERATION ---
  data <- simulate.multi.study.surrogates(
    epsilon = epsilon,
    M = M,
    sample_sizes = sample_sizes,
    J = J,
    prop_valid = 1, # 0% valid surrogates
    u_tau_min = 0,
    u_tau_max = u_tau_max,
    u_nu_min = 0,
    u_nu_max = u_nu_max,
    prop_invalid_under = 0.5,
    valid_mean_discrete = mu_valid
  )
  
  # Storage for p-values
  p_vals <- numeric(J)
  
  for (j in seq_len(J)) { # For each marker
    resj <- delta.reml.meta(
      delta = data$delta[, j],
      sd.delta = data$sd.delta[, j],
      epsilon = epsilon,
      alternative = "two.sided",
      test = "knha",
      meta.analysis.method = "RE"
    )
    p_vals[j] <- resj$results$p
  }
  
  # Derive FPR for each value of alpha
  tpr_vec <- colMeans(outer(p_vals, alpha_grid, "<"), na.rm = TRUE)
  
  # Store results 
  results_list[[i]] <- tibble(
    M = M,
    nm = results$nm[i],
    u_tau_max = u_tau_max,
    u_nu_max = u_nu_max,
    alpha = alpha_grid,
    tpr = tpr_vec,
    mu_valid = mu_valid
  )
}

# Extract results into a dataframe
results <- bind_rows(results_list)

# Save results
saveRDS(results, file = fs::path(simulation_results_folder, "simulation_parametric_calibration_validStrength.rds"))

results = readRDS(fs::path(simulation_results_folder, "simulation_parametric_calibration_validStrength.rds"))

# Define breaks and labels in terms of epsilon
x_breaks <- c(-epsilon, -epsilon/2, 0, epsilon/2, epsilon)
x_labels <- c(
  expression(-epsilon),
  expression(-epsilon/2),
  "0",
  expression(epsilon/2),
  expression(epsilon)
)

u_tau_max_fixed <- epsilon /10
u_nu_max_fixed <- epsilon

results_plot = results %>% 
  filter(alpha == 0.05,
         u_tau_max == u_tau_max_fixed,
         u_nu_max == u_nu_max_fixed)

p1 <- ggplot(results_plot, aes(
  x = mu_valid,
  y = tpr,
  color = factor(M),
  group = M
)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.5) +
  labs(title = "Empirical power as a function of valid surrogate mean by number of studies",
       x = expression(mu[valid]),
       y = "Empirical Power",
       color = "Number of studies"
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels
  ) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(size = 22, hjust = 0.5),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

p1

ggsave(
  filename = "power_plot_validStrength.pdf",
  path = simulation_figures_folder,
  plot = p1,
  width = 35,
  height = 18,
  units = "cm"
)

rm(list = ls())
