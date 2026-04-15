# Script to derive results for parametric simulation of trial-level surrogate effects
# In this script, we generate only invalid surrogates, in order to examine the false positive rate 

seed = 01042026
set.seed(seed)

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
J = 100000

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
  data <- generate.example.data.highdim.multistudy(
    epsilon = epsilon,
    M = M,
    sample_sizes = sample_sizes,
    J = J,
    prop_valid = 1, # 100% valid surrogates
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
  tpr_vec <- colMeans(outer(p_vals, alpha_grid, "<"), na.rm = TRUE)
  
  # Store results 
  results_list[[i]] <- tibble(
    M = M,
    nm = results$nm[i],
    u_tau_max = u_tau_max,
    u_nu_max = u_nu_max,
    alpha = alpha_grid,
    tpr = tpr_vec
  )
}

# Extract results into a dataframe
results <- bind_rows(results_list)

# Save results
saveRDS(results, file = fs::path(simulation_results_folder, "simulation_parametric_power_main.rds"))

results = readRDS(fs::path(simulation_results_folder, "simulation_parametric_power_main.rds"))

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
  y = tpr,
  group = 1
)) +
  geom_line(size = 2, alpha = 0.6) +
  geom_point(size = 2, alpha = 0.6) +
  scale_x_log10(
    limits = c(min_val, max_val),
    breaks = c(alpha_grid, max_val),
    labels = fmt_percent
  ) +
  labs(
    title = "Empirical powervs significance level across number of studies and sample sizes",
    x = "Nominal FPR",
    y = "Empirical Power"
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
  ) + 
  ylim(0, 1) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red")

p1

ggsave(
  filename = "power_grid_samples.pdf",
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
  y = tpr,
  group = 1
)) +
  geom_line(size = 2, alpha = 0.6) +
  geom_point(size = 2, alpha = 0.6) +
  scale_x_log10(
    limits = c(min_val, max_val),
    breaks = c(alpha_grid, max_val),
    labels = fmt_percent
  ) +
  labs(
    title = "Empirical power vs significance level across between- and within-study variability",
    x = "Nominal FPR",
    y = "Empirical power"
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
  )  + 
  ylim(0, 1) + 
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red")

p2

# Save
ggsave(
  filename = "power_grid_heterogeneity.pdf",
  path = simulation_figures_folder,
  plot = p2,
  width = 75,
  height = 35,
  units = "cm"
)


# Plot 3 - line plots of power against between-study heterogeneity,
# faceted by number of studies, coloured by within-study heterogeneity

epsilon <- 0.1

u_tau_max_vals_plot3 <- c(epsilon / 100, epsilon / 10, epsilon, epsilon * 10, epsilon * 100)
u_nu_max_vals_plot3  <- c(epsilon / 100, epsilon / 10, epsilon, epsilon * 10, epsilon * 100)

fixed_nm <- 50

# Break labels for x axis
tau_levels <- c(
  "epsilon/100",
  "epsilon/10",
  "epsilon",
  "10*epsilon",
  "100*epsilon"
)

tau_labels <- c(
  expression(epsilon/100),
  expression(epsilon/10),
  expression(epsilon),
  expression(10*epsilon),
  expression(100*epsilon)
)

# Break labels for legend, without division by n
nu_levels <- c(
  "epsilon/100",
  "epsilon/10",
  "epsilon",
  "10*epsilon",
  "100*epsilon"
)

nu_labels <- c(
  expression(epsilon/100),
  expression(epsilon/10),
  expression(epsilon),
  expression(10*epsilon),
  expression(100*epsilon)
)

results_plot3 <- results %>% 
  filter(
    u_tau_max %in% u_tau_max_vals_plot3,
    u_nu_max %in% u_nu_max_vals_plot3,
    nm == fixed_nm,
    alpha == 0.05
  ) %>%
  mutate(
    tau_lab = case_when(
      dplyr::near(u_tau_max, epsilon / 100) ~ "epsilon/100",
      dplyr::near(u_tau_max, epsilon / 10)  ~ "epsilon/10",
      dplyr::near(u_tau_max, epsilon)       ~ "epsilon",
      dplyr::near(u_tau_max, epsilon * 10)  ~ "10*epsilon",
      dplyr::near(u_tau_max, epsilon * 100) ~ "100*epsilon"
    ),
    nu_lab = case_when(
      dplyr::near(u_nu_max, epsilon / 100) ~ "epsilon/100",
      dplyr::near(u_nu_max, epsilon / 10)  ~ "epsilon/10",
      dplyr::near(u_nu_max, epsilon)       ~ "epsilon",
      dplyr::near(u_nu_max, epsilon * 10)  ~ "10*epsilon",
      dplyr::near(u_nu_max, epsilon * 100) ~ "100*epsilon"
    ),
    tau_lab = factor(tau_lab, levels = tau_levels),
    nu_lab  = factor(nu_lab, levels = nu_levels)
  )

p3 <- ggplot(results_plot3, aes(
  x = tau_lab,
  y = tpr,
  color = nu_lab,
  group = nu_lab
)) +
  geom_line(size = 1.2, alpha = 0.65) +
  geom_point(size = 3, alpha = 0.65) +
  facet_wrap(~ M, labeller = labeller(
    M = function(x) paste("M =", x)
  )) +
  scale_x_discrete(
    breaks = tau_levels,
    labels = tau_labels
  ) +
  scale_color_manual(
    values = scales::viridis_pal(option = "D")(length(nu_levels)),
    breaks = nu_levels,
    labels = nu_labels,
    name = expression(u[nu*","*max])
  ) +
  labs(
    x = expression(u[tau*","*max]),
    y = "Empirical Power",
    title = "Empirical power across number of studies and heterogeneity settings"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(size = 30, hjust = 0.5),
    panel.spacing = unit(5, "lines"),
    strip.background = element_rect(fill = "#f0f0f0", color = "black"),
    strip.text = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p3

# Save
ggsave(
  filename = "power_lines.pdf",
  path = simulation_figures_folder,
  plot = p3,
  width = 40,
  height = 18,
  units = "cm"
)

rm(list = ls())
