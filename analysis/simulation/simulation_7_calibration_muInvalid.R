library(SurrogateRank)
library(tidyverse)
library(fs)
library(ggplot2)
library(scales)
library(viridis)

simulation_figures_folder <- fs::path("output", "figures", "simulation")
simulation_results_folder <- fs::path("output", "results", "simulation")

J <- 1000
epsilon <- 0.1

# Define a grid of nominal significance levels
alpha_grid <- c(0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2)

Ms <- c(25)
tau_max_vals <- c(epsilon/10, epsilon, epsilon*10)
nu_max_vals <- c(epsilon/10, epsilon, epsilon*10)
fixed_sample_size <- 100
abs_mu_invalid_grid <- c(seq(epsilon, 0.2, 0.02), 0.3)

# Expand results for all alpha values
results <- expand.grid(
  M = Ms,
  abs_mu_invalid = abs_mu_invalid_grid,
  alpha = alpha_grid,
  stringsAsFactors = FALSE
) %>%
  mutate(fpr_obs = NA_real_)

for (i in seq_len(nrow(results))) {
  M <- results$M[i]
  u_tau_max <- tau_max_vals
  u_nu_max <- nu_max_vals
  
  mu_invalid_vec <- c(-results$abs_mu_invalid[i], results$abs_mu_invalid[i])
  sample_sizes <- rep(fixed_sample_size, M)
  
  # Simulate data once per configuration of M and abs_mu_invalid
  data <- simulate.multi.study.surrogates(
    epsilon = epsilon,
    M = M,
    sample_sizes = sample_sizes,
    J = J,
    prop_valid = 0,
    u_tau_min = 0,
    u_tau_max = u_tau_max,
    u_nu_min = 0,
    u_nu_max = u_nu_max,
    prop_invalid_under = 0.5,
    invalid_at_boundary = FALSE, 
    invalid_mean_discrete = mu_invalid_vec
  )
  
  # Compute p-values for all replicates
  p_vals <- numeric(J)
  for (j in seq_len(J)) {
    resj <- delta.reml.meta(
      delta = data$delta[, j],
      sd.delta = data$sd.delta[, j],
      epsilon = epsilon,
      alpha = results$alpha[i],  # current nominal alpha
      alternative = "two.sided",
      test = "knha",
      meta.analysis.method = "RE"
    )
    p_vals[j] <- resj$results$p
  }
  
  # Observed FPR at this nominal alpha
  results$fpr_obs[i] <- mean(p_vals < results$alpha[i], na.rm = TRUE)
}

# Optional: reshape for plotting
results_long <- results %>%
  rename(nominal_alpha = alpha, observed_fpr = fpr_obs)


saveRDS(results, file = fs::path(simulation_results_folder, "simulation_7_calibration.rds"))

# Get unique abs_mu_invalid values for color scale
abs_mu_vals <- sort(unique(results_long$abs_mu_invalid))
n_colors <- length(abs_mu_vals)

min_val <- 1e-4
max_val <- 1

# custom percent labels: 0.0001 -> 0.01%, 0.2 -> 20%, 1 -> 100%
fmt_percent <- function(x) {
  lab <- scales::label_percent(accuracy = 0.01)(x)
  sub("\\.?0+%$", "%", lab)
}

p1 <- ggplot(results_long, aes(x = nominal_alpha, y = observed_fpr, color = factor(abs_mu_invalid))) +
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
  labs(
    x = "Nominal FPR",
    y = "Observed FPR",
    color = expression(abs(mu))
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

p1

ggsave(
  filename = "nominal_observed_mu.pdf",
  path = simulation_figures_folder,
  plot = p1,
  width = 28,
  height = 18,
  units = "cm"
)

# rm(list = ls())
