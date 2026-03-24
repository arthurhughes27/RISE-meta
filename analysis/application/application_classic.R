library(Surrogate)
library(SurrogateRank)
library(dplyr)
library(tidyverse)
library(knitr)
library(kableExtra)

data("ARMD")

# -------------------------------------------------------------------
# 1) Standardize names and coding once
# -------------------------------------------------------------------
armd <- ARMD %>%
  transmute(
    Id      = as.factor(Id),
    Center  = as.factor(Center),
    Treat   = as.integer(Treat),
    # -1 = control, 1 = treatment
    Treat01 = if_else(Treat == 1L, 1L, 0L),
    # 0/1 coding for Molenberghs
    Diff24  = Diff24,
    Diff52  = Diff52
  )

# -------------------------------------------------------------------
# 2) Keep only centers with at least 3 individuals in each arm
# -------------------------------------------------------------------
min_per_arm <- 3L

center_keep <- armd %>%
  group_by(Center, Treat) %>%
  summarise(n_patients = n_distinct(Id), .groups = "drop") %>%
  group_by(Center) %>%
  filter(all(n_patients >= min_per_arm)) %>%
  pull(Center)

armd_filtered <- armd %>%
  filter(Center %in% center_keep) %>%
  arrange(Center, Id)

# Optional: summary of sample sizes by center and treatment
armd_summary <- armd_filtered %>%
  group_by(Center, Treat) %>%
  summarise(n_patients = n_distinct(Id), .groups = "drop") %>%
  arrange(Center, Treat)

n_centers <- n_distinct(armd_filtered$Center)

# -------------------------------------------------------------------
# 3) RISE analysis
# -------------------------------------------------------------------
rise_input <- armd_filtered %>%
  mutate(study = as.character(Center))

yone    <- rise_input %>% filter(Treat ==  1) %>% pull(Diff52)
yzero   <- rise_input %>% filter(Treat == -1) %>% pull(Diff52)
sone    <- rise_input %>% filter(Treat ==  1) %>% pull(Diff24)
szero   <- rise_input %>% filter(Treat == -1) %>% pull(Diff24)
studyone <- rise_input %>% filter(Treat ==  1) %>% pull(study)
studyzero <- rise_input %>% filter(Treat == -1) %>% pull(study)

rise_fit <- rise.screen.meta(
  yone = yone,
  yzero = yzero,
  sone = sone,
  szero = szero,
  studyone = studyone,
  studyzero = studyzero,
  alpha = 0.05,
  epsilon.study = 0.2,
  epsilon.meta = 0.2,
  p.correction = "none",
  return.study.similarity.plot = FALSE,
  test = "knha",
  meta.analysis.method = "RE"
)

ccc = rise_fit[["evaluation.metrics.meta"]][["ccc"]]
plot_rise = rise_fit$gamma.s.plot$fit.plot
# plot_rise_forest = rise_fit$gamma.s.plot$forest.plot
# plot_rise_forest

# -------------------------------------------------------------------
# 4) Molenberghs / Buyse trial-level model
# -------------------------------------------------------------------
fit_molen_fixed <- BifixedContCont(
  Dataset  = armd_filtered,
  Surr     = Diff24,
  True     = Diff52,
  Treat    = Treat01,
  Trial.ID = Center,
  Pat.ID   = Id
)

trial_effects <- fit_molen_fixed$Results.Stage.1 %>%
  transmute(
    Center = as.factor(Trial),
    N      = Obs.per.trial,
    Alpha  = Treatment.S,
    # trial-level surrogate effect
    Beta   = Treatment.T    # trial-level true endpoint effect
  )

trial_fit <- TrialLevelMA(
  Alpha.Vector = trial_effects$Alpha,
  Beta.Vector  = trial_effects$Beta,
  N.Vector     = trial_effects$N,
  Weighted     = TRUE,
  Alpha        = 0.05
)

r2_trial <- trial_fit[["Trial.R2"]][["R2 Trial"]]


# -------------------------------------------------------------------
# 5) Quick outputs
# -------------------------------------------------------------------
print(plot_rise)
plot.TrialLevelMA(trial_fit)

ccc
r2_trial


# rm(list = ls())
