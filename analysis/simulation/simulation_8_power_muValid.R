library(SurrogateRank)
library(tidyverse)
library(fs)
library(ggplot2)
library(scales)
library(viridis)

simulation_figures_folder <- fs::path("output", "figures", "simulation")
simulation_results_folder <- fs::path("output", "results", "simulation")

J <- 10000
epsilon <- 0.1

# Define a grid of nominal significance levels
alpha_grid <- c(0.05)

Ms <- c(5, 10, 15, 20, 25)
tau_max_vals <- c(epsilon/5)
nu_max_vals <- c(epsilon/10)
fixed_sample_size <- 100
mu_valid_grid <- seq(-epsilon, epsilon, epsilon/10)

# Expand results for all alpha values
results <- expand.grid(
  M = Ms,
  mu_valid = mu_valid_grid,
  alpha = alpha_grid,
  stringsAsFactors = FALSE
) %>%
  mutate(tpr_obs = NA_real_)

for (i in seq_len(nrow(results))) {
  M <- results$M[i]
  u_tau_max <- tau_max_vals
  u_nu_max <- nu_max_vals
  
  mu_valid_vec <- results$mu_valid[i]
  sample_sizes <- rep(fixed_sample_size, M)
  
  # Simulate data once per configuration of M and mu_valid
  data <- simulate.multi.study.surrogates(
    epsilon = epsilon,
    M = M,
    sample_sizes = sample_sizes,
    J = J,
    prop_valid = 1,
    u_tau_min = 0,
    u_tau_max = u_tau_max,
    u_nu_min = 0,
    u_nu_max = u_nu_max,
    prop_invalid_under = 0.5,
    invalid_at_boundary = FALSE, 
    valid_mean_discrete = mu_valid_vec
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
  
  # Observed tpr at this nominal alpha
  results$tpr_obs[i] <- mean(p_vals < results$alpha[i], na.rm = TRUE)
}

# Optional: reshape for plotting
results_long <- results %>%
  rename(nominal_alpha = alpha, observed_tpr = tpr_obs)


saveRDS(results, file = fs::path(simulation_results_folder, "simulation_8_power_muValid.rds"))

# Define breaks and labels in terms of epsilon
x_breaks <- c(-epsilon, -epsilon/2, 0, epsilon/2, epsilon)
x_labels <- c(
  expression(-epsilon),
  expression(-epsilon/2),
  "0",
  expression(epsilon/2),
  expression(epsilon)
)

p1 <- ggplot(results_long, aes(
  x = mu_valid,
  y = observed_tpr,
  color = factor(M),
  group = M
)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.5) +
  labs(title = "Empirical power as a function of distribution mean by number of studies",
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
  filename = "powerPlot.pdf",
  path = simulation_figures_folder,
  plot = p1,
  width = 30,
  height = 18,
  units = "cm"
)

# rm(list = ls())
