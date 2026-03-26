library(Surrogate)
library(SurrogateRank)
library(ggplot2)
library(dplyr)

application_figures_folder <- fs::path("output", "figures", "application", "classic")

# Load data
data(ARMD)

# Keep only centers with at least 1 participant per arm and at least 5 overall
center_counts <- ARMD %>%
  group_by(Center, Treat) %>%
  summarise(n_patients = n(), .groups = "drop")

ok_centers <- center_counts %>%
  group_by(Center) %>%
  filter(all(n_patients >= 2) & sum(n_patients) >= 5) %>%
  pull(Center)

dat <- ARMD %>%
  filter(Center %in% ok_centers)

# Full joint model (Molenberghs/Buyse)
fit_full <- BifixedContCont(
  Dataset  = dat,
  Surr     = Diff24,
  True     = Diff52,
  Treat    = Treat,
  Trial.ID = Center,
  Pat.ID   = Id,
  Model    = "Full",
  Weighted = TRUE
)

# Trial-level surrogacy metrics
trial_R2 <- fit_full$Trial.R2$`R2 Trial`
trial_R  <- fit_full$Trial.R$`R Trial`

cat("Full-model trial-level R^2:", trial_R2, "\n")
cat("Full-model trial-level R  :", trial_R, "\n")

trial_effects <- fit_full$Results.Stage.1 %>%
  transmute(
    Center = as.factor(Trial),
    n      = Obs.per.trial,
    u.s    = Treatment.S,
    u.y    = Treatment.T
  )

# Adaptive legend
n_vals <- trial_effects$n
if (length(unique(n_vals)) == 1) {
  legend_breaks <- unique(n_vals)
  legend_labels <- as.character(unique(n_vals))
} else {
  min_val <- min(n_vals)
  max_val <- max(n_vals)
  median_val <- median(n_vals)
  mid_val <- round(median_val)
  if (mid_val <= min_val) mid_val <- min_val + 1
  if (mid_val >= max_val) mid_val <- max_val - 1
  legend_breaks <- c(min_val, mid_val, max_val)
  legend_labels <- as.character(legend_breaks)
}

# Full-model plot
plot.min.global <- min(trial_effects$u.s, trial_effects$u.y) - 2
jointModel_plot <- trial_effects %>%
  ggplot(aes(x = u.s, y = u.y)) +
  geom_point(aes(size = n), shape = 21, alpha = 0.5, stroke = 1, fill = "#6FB1EF") +
  geom_abline(slope = 1, intercept = 0, color = "#FF2128", linetype = "dashed", linewidth = 0.8, alpha = 0.5) +
  annotate("text", x = plot.min.global, y = max(trial_effects$u.y, na.rm = TRUE),
           label = paste0("Trial-level R2 = ", round(trial_R2, 2)),
           hjust = 0, vjust = 1, color = "red", size = 8) +
  scale_size_continuous(range = c(5, 20), breaks = legend_breaks, labels = legend_labels) +
  scale_x_continuous(limits = c(plot.min.global, max(trial_effects$u.s, na.rm = TRUE) + 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(plot.min.global, max(trial_effects$u.y, na.rm = TRUE) + 2), expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  labs(title = "Trial-level surrogacy: bivariate joint modelling approach",
       x = "Treatment effect on change in visual acuity at 6 months",
       y = "Treatment effect on change in visual acuity at 12 months",
       size = "Center N") +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 18),
        legend.position = "right")

jointModel_plot

# RISE analysis
yone    <- dat$Diff52[dat$Treat ==  1]
yzero   <- dat$Diff52[dat$Treat == -1]
sone    <- dat$Diff24[dat$Treat ==  1]
szero   <- dat$Diff24[dat$Treat == -1]
Centerone <- as.character(dat$Center[dat$Treat ==  1])
Centerzero <- as.character(dat$Center[dat$Treat == -1])

rise_fit <- rise.screen.meta(
  yone = yone, yzero = yzero, sone = sone, szero = szero,
  studyone = Centerone, studyzero = Centerzero,
  alpha = 0.05, 
  epsilon.study = 0.2, epsilon.meta = 0.2,
  p.correction = "none", return.study.similarity.plot = FALSE,
  paired.all = FALSE, test = "knha", meta.analysis.method = "RE"
)

gamma_df <- rise_fit[["screening.metrics.study"]]
n_vals <- gamma_df$n
if (length(unique(n_vals)) == 1) {
  legend_breaks <- unique(n_vals)
  legend_labels <- as.character(unique(n_vals))
} else {
  min_val <- min(n_vals)
  max_val <- max(n_vals)
  median_val <- median(n_vals)
  mid_val <- round(median_val)
  if (mid_val <= min_val) mid_val <- min_val + 1
  if (mid_val >= max_val) mid_val <- max_val - 1
  legend_breaks <- c(min_val, mid_val, max_val)
  legend_labels <- as.character(legend_breaks)
}

ccc <- rise_fit[["evaluation.metrics.meta"]][["ccc"]]

riseMeta_plot <- gamma_df %>%
  ggplot(aes(x = u.s, y = u.y)) +
  geom_point(aes(size = n), shape = 21, alpha = 0.5, stroke = 1, fill = "#6FB1EF") +
  geom_abline(slope = 1, intercept = 0, color = "#FF2128", linetype = "dashed", linewidth = 0.8, alpha = 0.5) +
  annotate("text", x = -0.1, y = 1.05, label = paste0("CCC = ", round(ccc, 2)),
           hjust = 0, vjust = 1, color = "red", size = 8) +
  scale_size_continuous(range = c(5, 20), breaks = legend_breaks, labels = legend_labels) +
  scale_x_continuous(limits = c(-0.1, 1.1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.1, 1.1), expand = c(0, 0)) +
  coord_fixed(ratio = 1) +
  labs(title = "Trial-level surrogacy: RISE-Meta",
       x = "Treatment effect on change in visual acuity at 6 months",
       y = "Treatment effect on change in visual acuity at 12 months",
       size = "Center N") +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 18),
        legend.position = "right")

riseMeta_plot

# Save these graphics with an informative name
ggsave(
  filename = paste0("jointModel_plot_ARMD.pdf"),
  path     = application_figures_folder,
  plot     = jointModel_plot,
  width    = 37,
  height   = 20,
  units    = "cm"
)

ggsave(
  filename = paste0("riseMeta_plot_ARMD.pdf"),
  path     = application_figures_folder,
  plot     = riseMeta_plot,
  width    = 37,
  height   = 20,
  units    = "cm"
)
