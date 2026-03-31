# Script to derive results for parametric simulation of trial-level surrogate effects
# In this script, we generate only invalid surrogates, in order to examine the false positive rate 

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
M_grid <- c(10)

# Grid of number of samples within each study to generate
nm_grid = c(50)

# Grid of maximum between-study variability values
u_tau_max_vals <- c(epsilon/100, epsilon/10, epsilon, epsilon*10, epsilon*100)

# Grid of maximum within-study variability values
u_nu_max_vals <- c(epsilon/100, epsilon/10, epsilon, epsilon*10, epsilon*100)

# Grid of alpha values for the calibration plot
alpha_grid <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5)

# Grid of values for invalid surrogate mean generation
abs_mu_invalid_vals <- c(seq(epsilon, 0.3, 0.05))

# Initialise data frame WITHOUT alpha
results <- expand.grid(
  M = M_grid,
  nm = nm_grid,
  u_tau_max = u_tau_max_vals,
  u_nu_max = u_nu_max_vals,
  abs_mu_invalid = abs_mu_invalid_vals,
  stringsAsFactors = FALSE
)

# Initialise list for storage
results_list <- vector("list", nrow(results))

for (i in seq_len(nrow(results))) {
  
  M <- results$M[i]
  sample_sizes <- rep(results$nm[i], M)
  u_tau_max <- results$u_tau_max[i]
  u_nu_max <- results$u_nu_max[i]
  abs_mu_invalid = results$abs_mu_invalid[i]
  
  mu_invalid_vec <- c(-abs_mu_invalid, abs_mu_invalid)
  
  # --- DATA GENERATION ---
  data <- simulate.multi.study.surrogates(
    epsilon = epsilon,
    M = M,
    sample_sizes = sample_sizes,
    J = J,
    prop_valid = 0, # 0% valid surrogates
    u_tau_min = 0,
    u_tau_max = u_tau_max,
    u_nu_min = 0,
    u_nu_max = u_nu_max,
    prop_invalid_under = 0.5,
    invalid_mean_discrete = mu_invalid_vec
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
  fpr_vec <- colMeans(outer(p_vals, alpha_grid, "<"), na.rm = TRUE)
  
  # Store results 
  results_list[[i]] <- tibble(
    M = M,
    nm = results$nm[i],
    u_tau_max = u_tau_max,
    u_nu_max = u_nu_max,
    alpha = alpha_grid,
    fpr = fpr_vec,
    abs_mu_invalid = abs_mu_invalid
  )
}

# Extract results into a dataframe
results <- bind_rows(results_list)

# Save results
saveRDS(results, file = fs::path(simulation_results_folder, "simulation_parametric_calibration_invalidStrength.rds"))

results = readRDS(fs::path(simulation_results_folder, "simulation_parametric_calibration_invalidStrength.rds"))

# Calibration plot faceted by u_tau_max

epsilon <- 0.1

u_tau_max_vals <- c(epsilon / 10, epsilon, epsilon * 10)
u_nu_max_fixed <- epsilon / 10

# Helper labels for parsed facet strips
tau_levels <- c(
  "u[tau*','*max] == epsilon/100",
  "u[tau*','*max] == epsilon/10",
  "u[tau*','*max] == epsilon",
  "u[tau*','*max] == 10*epsilon",
  "u[tau*','*max] == 100*epsilon"
)

results_plot <- results %>% 
  filter(
    u_tau_max %in% u_tau_max_vals,
    u_nu_max == u_nu_max_fixed
  ) %>%
  mutate(
    tau_lab = case_when(
      dplyr::near(u_tau_max, epsilon / 100) ~ "u[tau*','*max] == epsilon/100",
      dplyr::near(u_tau_max, epsilon / 10)  ~ "u[tau*','*max] == epsilon/10",
      dplyr::near(u_tau_max, epsilon)       ~ "u[tau*','*max] == epsilon",
      dplyr::near(u_tau_max, epsilon * 10)  ~ "u[tau*','*max] == 10*epsilon",
      dplyr::near(u_tau_max, epsilon * 100) ~ "u[tau*','*max] == 100*epsilon"
    ),
    tau_lab = factor(tau_lab, levels = tau_levels)
  )

# Get unique abs_mu_invalid values for color scale
abs_mu_vals <- sort(unique(results_plot$abs_mu_invalid))
n_colors <- length(abs_mu_vals)

min_val <- 1e-4
max_val <- 1

# Custom percent labels: 0.0001 -> 0.01%, 0.2 -> 20%, 1 -> 100%
fmt_percent <- function(x) {
  lab <- scales::label_percent(accuracy = 0.01)(x)
  sub("\\.?0+%$", "%", lab)
}

p1 <- ggplot(results_plot, aes(x = alpha, y = fpr, color = factor(abs_mu_invalid))) +
  geom_ribbon(
    data = tibble(x = c(min_val, max_val)),
    aes(x = x, ymin = x, ymax = max_val),
    inherit.aes = FALSE,
    fill = "grey80",
    alpha = 0.5
  ) +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.6) +
  scale_x_log10(
    limits = c(min_val, max_val),
    breaks = c(alpha_grid, max_val),
    labels = fmt_percent
  ) +
  scale_y_log10(
    limits = c(min_val, max_val),
    breaks = c(alpha_grid, max_val),
    labels = fmt_percent
  ) +
  scale_color_viridis(
    discrete = TRUE,
    option = "D",
    begin = 0, end = 0.9
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  facet_wrap(~ tau_lab, labeller = label_parsed, nrow = 1) +
  labs(
    title = "Calibration plots across invalid surrogate means and between-study heterogeneity",
    x = "Nominal FPR",
    y = "Observed FPR",
    color = expression(abs(mu[invalid]))
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey80", color = NA),
    strip.text = element_text(face = "bold", size = 30),
    axis.title = element_text(size = 35),
    plot.title = element_text(hjust = 0.5, size = 45)
  )

p1

ggsave(
  filename = "calibration_plot_invalidStrength.pdf",
  path = simulation_figures_folder,
  plot = p1,
  width = 75,
  height = 18,
  units = "cm"
)

rm(list = ls())
