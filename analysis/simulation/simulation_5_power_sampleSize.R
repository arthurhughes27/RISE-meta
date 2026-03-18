library(SurrogateRank)
library(tidyverse)

simulation_figures_folder = fs::path("output", "figures", "simulation")
simulation_results_folder = fs::path("output", "results", "simulation")

J <- 10000
M <- 5
epsilon <- 0.1
alpha <- 0.05

n_vals       <- c(10, 25, 50, 75, 100)
tau_max_vals <- c(0.001, 0.01, 0.05, 0.1, 1, 5, 10)
nu_max_vals  <- c(0.001, 0.01, 0.05, 0.1, 1, 5, 10)

results <- expand.grid(
  n         = n_vals,
  u_tau_max = tau_max_vals,
  u_nu_max  = nu_max_vals,
  stringsAsFactors = FALSE
) %>%
  mutate(tpr = NA_real_)

for (i in seq_len(nrow(results))) {
  n         <- results$n[i]
  u_tau_max <- results$u_tau_max[i]
  u_nu_max  <- results$u_nu_max[i]

  sample_sizes <- rep(n, M)

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
    prop_invalid_under = 0.5
  )

  p_vals <- numeric(J)
  for (j in seq_len(J)) {
    resj <- delta.reml.meta(
      delta = data$delta[, j],
      sd.delta = data$sd.delta[, j],
      epsilon = epsilon,
      alpha = alpha,
      alternative = "two.sided",
      test = "knha"
    )
    p_vals[j] <- resj$results$p
  }

  results$tpr[i] <- mean(p_vals < alpha, na.rm = TRUE)
}


saveRDS(results, file = fs::path(simulation_results_folder, "simulation_5_power_sampleSize.rds"))

p1 <- ggplot(results, aes(
  x = factor(u_tau_max, levels = sort(unique(u_tau_max))),
  y = tpr,
  color = factor(u_nu_max),
  group = factor(u_nu_max)
)) +
  geom_line(size = 1.2, alpha = 0.65) +
  geom_point(size = 3, alpha = 0.65) +
  facet_wrap(
    ~ n,
    nrow = 1,
    labeller = labeller(n = function(x) paste("n =", x))
  ) +
  scale_color_manual(
    values = scales::viridis_pal(option = "D")(length(unique(results$u_nu_max))),
    name = "Max within-study variance"
  ) +
  ylim(0, 1) +
  labs(
    x = "Maximum between-trial variance",
    y = "Power",
    title = "Power - true positive rate across per-study sample sizes (N trials = 5)"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(size = 26, hjust = 0.5),
    panel.spacing = unit(5, "lines"),
    strip.background = element_rect(fill = "#f0f0f0", color = "black"),
    strip.text = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p1

ggsave(
  filename = "power_lineplot_sampleSize.pdf",
  path = simulation_figures_folder,
  plot = p1,
  width = 60,
  height = 18,
  units = "cm"
)

rm(list = ls())