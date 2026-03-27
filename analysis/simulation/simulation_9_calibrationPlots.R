library(SurrogateRank)
library(tidyverse)
library(ggplot2)
library(viridis)
library(dplyr)
library(tibble)
library(scales)

simulation_figures_folder = fs::path("output", "figures", "simulation")
simulation_results_folder = fs::path("output", "results", "simulation")

J <- 10000
epsilon <- 0.1
alpha <- 0.05

Ms <- c(3, 10, 25)
tau_max_vals <- c(0.001, 0.01, 0.1, 1)
nu_max_vals <- c(epsilon)
fixed_sample_size = 50
alpha_grid <- c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2)

results <- expand.grid(
  M = Ms,
  u_tau_max = tau_max_vals,
  u_nu_max = nu_max_vals,
  stringsAsFactors = FALSE,
  alpha = alpha_grid
) %>%
  mutate(fpr = NA_real_)

for (i in seq_len(nrow(results))) {
  M <- results$M[i]
  u_tau_max <- results$u_tau_max[i]
  u_nu_max <- results$u_nu_max[i]
  
  sample_sizes <- rep(fixed_sample_size, M)
  
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
    invalid_at_boundary = TRUE
  )
  
  p_vals <- numeric(J)
  for (j in seq_len(J)) {
    resj <- delta.reml.meta(
      delta = data$delta[, j],
      sd.delta = data$sd.delta[, j],
      epsilon = epsilon,
      alpha = alpha,
      alternative = "two.sided",
      test = "knha",
      meta.analysis.method = "RE"
    )
    p_vals[j] <- resj$results$p
  }
  
  results$fpr[i] <- mean(p_vals < results$alpha[i], na.rm = TRUE)
}

# Define plotting limits
min_val <- 1e-5
max_val <- 1

# Custom percent labels
fmt_percent <- function(x) {
  lab <- scales::label_percent(accuracy = 0.01)(x)
  sub("\\.?0+%$", "%", lab)
}

# tol = min(alpha_grid)
# 
# results = results %>% 
#   mutate(fpr = ifelse(fpr == 0, fpr + tol, fpr)) 

# Ensure u_tau_max is treated as factor for coloring
results$u_tau_max <- factor(results$u_tau_max)

p1 <- ggplot(results, aes(
  x = alpha,
  y = fpr,
  color = u_tau_max,
  group = u_tau_max
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
  scale_color_viridis(
    discrete = TRUE,
    option = "D",
    begin = 0, end = 0.9
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(
    title = "Calibration across number of studies and between-study variance",
    x = "Nominal FPR",
    y = "Observed FPR",
    color = expression(u[tau[max]])
  ) +
  facet_grid(. ~ M, labeller = labeller(M = function(x) paste0("M = ", x))) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey80", color = NA), 
    strip.text = element_text(face = "bold", size = 30),
    plot.title = element_text(size = 35, hjust = 0.5),
    axis.title = element_text(size = 30)
  )

p1

# Add a label column for the facets
results <- results %>%
  mutate(u_tau_max_label = paste0("max Tau = ", u_tau_max))  # \u03C4 is the tau symbol

# Then plot
p2 <- ggplot(results, aes(
  x = alpha,
  y = fpr,
  color = factor(M),
  group = factor(M)
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
  scale_color_viridis(
    discrete = TRUE,
    option = "D",
    begin = 0, end = 0.9
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(
    title = "Calibration across number of studies and between-study variance",
    x = "Nominal FPR",
    y = "Observed FPR",
    color = "Number of studies (M)"
  ) +
  facet_grid(
    . ~ u_tau_max_label
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey80", color = NA),
    strip.text = element_text(face = "bold", size = 30),
    plot.title = element_text(size = 35, hjust = 0.5),
    axis.title = element_text(size = 30)
  )

p2

ggsave(
  filename = "calibration_grid_byM.pdf",
  path = simulation_figures_folder,
  plot = p1,
  width = 75,
  height = 18,
  units = "cm"
)

ggsave(
  filename = "calibration_grid_byTau.pdf",
  path = simulation_figures_folder,
  plot = p2,
  width = 75,
  height = 18,
  units = "cm"
)

# rm(list = ls())
