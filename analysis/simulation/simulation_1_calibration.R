library(SurrogateRank)
library(tidyverse)

simulation_figures_folder = fs::path("output", "figures", "simulation")

J <- 10000
epsilon <- 0.1
alpha <- 0.05

Ms <- c(3, 10, 25)
tau_max_vals <- c(0.001, 0.01, 0.05, 0.1, 1, 5, 10)
nu_max_vals <- c(0.001, 0.01, 0.05, 0.1, 1, 5, 10)

results <- expand.grid(M = Ms, u_tau_max = tau_max_vals, u_nu_max = nu_max_vals,
                       stringsAsFactors = FALSE) %>%
  mutate(fpr = NA_real_)

for (i in seq_len(nrow(results))) {
  M <- results$M[i]
  u_tau_max <- results$u_tau_max[i]
  u_nu_max <- results$u_nu_max[i]
  
  sample_sizes <- seq(25, 25 * M, 25)
  
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
      sample.sizes = sample_sizes
    )
    p_vals[j] <- resj$results$p
  }
  
  results$fpr[i] <- mean(p_vals < alpha, na.rm = TRUE)
}

p1 = ggplot(results, aes(y = factor(u_tau_max), x = factor(u_nu_max), fill = fpr)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", fpr)), size = 4) +
  facet_wrap(~ M, labeller = labeller(M = function(x) paste("N trials =", x))) +
  scale_fill_gradientn(
    colours = c("#2166AC", "white", "#B2182B"),
    values = scales::rescale(c(0, alpha, 0.15)),
    limits = c(0, 0.15)
  ) +
  labs(
    x = "Maximum between-trial variance",
    y = "Maximum within-study variance",
    fill = "FPR",
    title = "Calibration - false positive rate across different settings"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 30, hjust = 0.5),
    panel.spacing = unit(5, "lines"),            # increase spacing between facets
    strip.background = element_rect(fill = "#f0f0f0", color = "black"),  # gray background for facet labels
    strip.text = element_text(face = "bold", size = 16)                   # bold facet text
  )

p1

ggsave(
  filename = "calibration_tileplot.pdf",
  path = simulation_figures_folder,
  plot = p1,
  width = 40,
  height = 18,
  units = "cm"
)


p2 <- ggplot(results, aes(x = factor(u_tau_max, levels = sort(unique(u_tau_max))),
                          y = fpr,
                          color = factor(u_nu_max),
                          group = factor(u_nu_max))) +
  geom_line(size = 1.2, alpha = 0.65) +
  geom_point(size = 3, alpha = 0.65) +
  facet_wrap(~ M, labeller = labeller(M = function(x) paste("N trials =", x))) +
  scale_color_manual(
    values = scales::viridis_pal(option = "D")(length(unique(results$u_nu_max))),
    name = "Max within-study variance"
  ) +
  ylim(0, 0.2) +
  geom_hline(yintercept = alpha, linetype = "dashed", color = "black", alpha = 0.7) +
  labs(
    x = "Maximum between-trial variance",
    y = "False Positive Rate (FPR)",
    title = "Calibration - false positive rate across different settings"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(size = 30, hjust = 0.5),
    panel.spacing = unit(5, "lines"),                            # wider separation between facets
    strip.background = element_rect(fill = "#f0f0f0", color = "black"),  # gray background for facet labels
    strip.text = element_text(face = "bold", size = 16),         # bold facet text
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"                                     # move legend to the right
  )

p2

ggsave(
  filename = "calibration_lineplot.pdf",
  path = simulation_figures_folder,
  plot = p2,
  width = 40,
  height = 18,
  units = "cm"
)
