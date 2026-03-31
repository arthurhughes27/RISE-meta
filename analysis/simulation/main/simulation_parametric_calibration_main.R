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
M_grid <- c(3, 10, 25)

# Grid of number of samples within each study to generate
nm_grid = c(10, 50, 250)

# Grid of maximum between-study variability values
u_tau_max_vals <- c(epsilon/100, epsilon/10, epsilon, epsilon*10, epsilon*100)

# Grid of maximum within-study variability values
u_nu_max_vals <- c(epsilon/100, epsilon/10, epsilon, epsilon*10, epsilon*100)

# Grid of alpha values for the calibration plot
alpha_grid <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5)

# Initialise data frame WITHOUT alpha
results <- expand.grid(
  M = M_grid,
  nm = nm_grid,
  u_tau_max = u_tau_max_vals,
  u_nu_max = u_nu_max_vals,
  stringsAsFactors = FALSE
)

# Initialise list for storage
results_list <- vector("list", nrow(results))

for (i in seq_len(nrow(results))) {
  
  M <- results$M[i]
  sample_sizes <- rep(results$nm[i], M)
  u_tau_max <- results$u_tau_max[i]
  u_nu_max <- results$u_nu_max[i]
  
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
    invalid_at_boundary = TRUE # Worst case-scenario for invalid surrogates
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
    fpr = fpr_vec
  )
}

# Extract results into a dataframe
results <- bind_rows(results_list)

# Save results
saveRDS(results, file = fs::path(simulation_results_folder, "simulation_parametric_calibration_main.rds"))

results = readRDS(fs::path(simulation_results_folder, "simulation_parametric_calibration_main.rds"))

# Plot 1: Grid of calibration plots for number of studies and sample sizes
fixed_u_tau_max = epsilon/10
fixed_u_nu_max = epsilon*100

# Keep only the fixed heterogeneity setting
results_plot <- results %>%
  filter(
   u_tau_max == fixed_u_tau_max, 
   u_nu_max == fixed_u_nu_max
  ) %>%
  mutate(
    M = factor(M, levels = sort(unique(M))),
    nm = factor(nm, levels = sort(unique(nm)))
  )

# Define plotting limits
min_val <- 1e-5
max_val <- 1

# Custom percent labels
fmt_percent <- function(x) {
  lab <- scales::label_percent(accuracy = 0.01)(x)
  sub("\\.?0+%$", "%", lab)
}

p1 <- ggplot(results_plot, aes(
  x = alpha,
  y = fpr,
  group = 1
)) +
  # Shaded grey region above y=x
  geom_ribbon(
    data = tibble(x = c(min_val, max_val)),
    aes(x = x, ymin = x, ymax = max_val),
    inherit.aes = FALSE,
    fill = "grey80",
    alpha = 0.5
  ) +
  geom_line(size = 2, alpha = 0.6) +
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
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(
    title = "Calibration plots across number of studies and sample sizes",
    x = "Nominal FPR",
    y = "Observed FPR"
  ) +
  facet_grid(
    nm ~ M,
    labeller = labeller(
      M = function(x) paste0("M = ", x),
      nm = function(x) paste0("n = ", x)
    )
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey80", color = NA),
    strip.text = element_text(face = "bold", size = 24),
    plot.title = element_text(size = 50, hjust = 0.5),
    axis.title = element_text(size = 45)
  )

p1

ggsave(
  filename = "calibration_grid_samples.pdf",
  path = simulation_figures_folder,
  plot = p1,
  width = 75,
  height = 35,
  units = "cm"
)


# Plot 2: Grid of calibration plots for between- and within-study variability
fixed_M  = 10
fixed_nm = 10

tau_vals <- c(epsilon / 10, epsilon, epsilon * 10)
nu_vals  <- c(epsilon, epsilon * 10, epsilon * 100)

# Prepare data + create parsed labels
results_plot2 <- results %>%
  filter(
    M == fixed_M,
    nm == fixed_nm,
    u_tau_max %in% tau_vals,
    u_nu_max %in% nu_vals
  ) %>%
  mutate(
    u_tau_lab = factor(
      paste0("u[tau*','*max] == ", format(u_tau_max, scientific = TRUE)),
      levels = paste0("u[tau*','*max] == ", format(tau_vals, scientific = TRUE))
    ),
    u_nu_lab = factor(
      paste0(
        "u[nu*','*max]/n == ",
        format(u_nu_max / fixed_nm, scientific = TRUE)
      ),
      levels = paste0(
        "u[nu*','*max]/n == ",
        format(nu_vals / fixed_nm, scientific = TRUE)
      )
    )
  )

# Define plotting limits
min_val <- 1e-5
max_val <- 1

# Custom percent labels
fmt_percent <- function(x) {
  lab <- scales::label_percent(accuracy = 0.01)(x)
  sub("\\.?0+%$", "%", lab)
}

# Plot
p2 <- ggplot(results_plot2, aes(
  x = alpha,
  y = fpr,
  group = 1
)) +
  geom_ribbon(
    data = tibble(x = c(min_val, max_val)),
    aes(x = x, ymin = x, ymax = max_val),
    inherit.aes = FALSE,
    fill = "grey80",
    alpha = 0.5
  ) +
  geom_line(size = 2, alpha = 0.6) +
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
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(
    title = "Calibration plots across between- and within-study variability",
    x = "Nominal FPR",
    y = "Observed FPR"
  ) +
  facet_grid(
    u_nu_lab ~ u_tau_lab,
    labeller = label_parsed
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey80", color = NA),
    strip.text = element_text(face = "bold", size = 24),
    plot.title = element_text(size = 50, hjust = 0.5),
    axis.title = element_text(size = 45)
  )

p2

# Save
ggsave(
  filename = "calibration_grid_heterogeneity.pdf",
  path = simulation_figures_folder,
  plot = p2,
  width = 75,
  height = 35,
  units = "cm"
)

rm(list = ls())
