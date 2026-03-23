# Meta-analysis with RISE: application to HIPC influenza vaccine data
library(tidyverse)
library(SurrogateRank)
library(knitr)
library(kableExtra)

set.seed(08012025)

# Paths to processed data and output figures
processed_data_folder <- "data"
application_figures_folder <- fs::path("output", "figures", "application")

# Load merged gene expression and BTM gene set objects
hipc_merged_all_noNorm <- readRDS(fs::path(processed_data_folder, "hipc_merged_all_noNorm.rds"))
BTM <- readRDS(fs::path(processed_data_folder, "BTM_processed.rds"))

# Gene columns present in the data with no missing values
gene_names <- hipc_merged_all_noNorm %>%
  select(a1cf:zzz3) %>%
  select(where( ~ !any(is.na(.)))) %>%
  colnames()

# Day 1 post-vaccination; timepoints_to_keep includes baseline (day 0)
tp <- 1
timepoints_to_keep <- c(0, tp)

# Meta-analysis method
meta.analysis.method = "RE"

# nAb studies use a different response variable than HAI studies
nab_studies <- c("SDY80", "SDY180", "SDY1276", "SDY67")

# Build a unified pre/post response column, then filter to influenza participants
# with exactly one baseline and one post-vaccination measurement per person,
# and drop studies with fewer than 3 participants
hipc_merged_all_noNorm_filtered <- hipc_merged_all_noNorm %>%
  mutate(
    response_pre = ifelse(
      study_accession %in% nab_studies,
      immResp_mean_nAb_pre_value,
      immResp_mean_hai_pre_value
    ),
    response_post = ifelse(
      study_accession %in% nab_studies,
      immResp_mean_nAb_post_value,
      immResp_mean_hai_post_value
    )
  ) %>%
  filter(
    !is.na(immResp_MFC_anyAssay_log2_MFC),
    vaccine_name == "Influenza (IN)",
    study_time_collected %in% timepoints_to_keep
  ) %>%
  group_by(participant_id) %>%
  filter(sum(study_time_collected == 0) == 1,
         sum(study_time_collected == tp) == 1) %>%
  ungroup() %>%
  group_by(study_accession) %>%
  filter(length(unique(participant_id)) > 2) %>%
  ungroup() %>%
  select(
    participant_id,
    age_imputed,
    gender,
    race,
    study_accession,
    study_time_collected,
    response_pre,
    response_post,
    all_of(gene_names)
  ) %>%
  arrange(participant_id)

# Sample 66% of participants per study for training; remainder becomes test set
train_indices <- hipc_merged_all_noNorm_filtered %>%
  distinct(study_accession, participant_id) %>%
  group_by(study_accession) %>%
  slice_sample(prop = 0.66) %>%
  ungroup()

df_train <- hipc_merged_all_noNorm_filtered %>%
  semi_join(train_indices, by = c("study_accession", "participant_id"))

df_test <- hipc_merged_all_noNorm_filtered %>%
  anti_join(train_indices, by = c("study_accession", "participant_id"))

# BTM gene sets and their names (excluding top-level aggregates)
btm_filter <- which(BTM[["geneset.aggregates"]] != "NA")
BTM_genes <- BTM[["genesets"]][btm_filter] %>% unlist()
BTM_genes_names <- BTM[["geneset.names.descriptions"]][btm_filter]

# Summarise gene expression within each BTM gene set by applying FUN row-wise
aggregate_to_geneset <- function(df, genesets, geneset_names, FUN = mean) {
  df <- as.data.frame(df)
  out <- imap_dfc(genesets, function(genes, i) {
    present <- intersect(genes, colnames(df))
    if (length(present) == 0)
      return(NULL)
    mat <- df[, present, drop = FALSE]
    vec <- apply(mat, 1, function(r)
      FUN(r, na.rm = TRUE))
    # Replace NaN (all-NA rows) with proper NA
    vec[is.nan(vec)] <- NA_real_
    tibble::tibble(!!geneset_names[i] := vec)
  })
  as.data.frame(out)
}

# Extract RISE inputs (response vectors, gene-set matrices, study labels)
# from a data frame containing both baseline and post-vaccination rows
extract_rise_inputs <- function(df) {
  sone_raw  <- df %>% filter(study_time_collected > 0)  %>% select(any_of(BTM_genes))
  szero_raw <- df %>% filter(study_time_collected == 0) %>% select(any_of(BTM_genes))
  
  list(
    yone      = df %>% filter(study_time_collected > 0)  %>% pull(response_post),
    yzero     = df %>% filter(study_time_collected == 0) %>% pull(response_pre),
    sone      = aggregate_to_geneset(sone_raw, BTM[["genesets"]][btm_filter], BTM_genes_names),
    szero     = aggregate_to_geneset(szero_raw, BTM[["genesets"]][btm_filter], BTM_genes_names),
    studyone  = df %>% filter(study_time_collected > 0)  %>% pull(study_accession),
    studyzero = df %>% filter(study_time_collected == 0) %>% pull(study_accession)
  )
}

# ----- Screening on training data -----

train_inputs <- extract_rise_inputs(df_train)

# Screen for surrogate markers across studies using BH-corrected meta-analysis
rise_screen_result <- rise.screen.meta(
  yone                         = train_inputs$yone,
  yzero                        = train_inputs$yzero,
  sone                         = train_inputs$sone,
  szero                        = train_inputs$szero,
  studyone                     = train_inputs$studyone,
  studyzero                    = train_inputs$studyzero,
  alpha                        = 0.05,
  epsilon.meta.mode            = "user",
  power.want.s.study           = 0.8,
  epsilon.meta                 = 0.2,
  alternative                  = "two.sided",
  paired.all                   = TRUE,
  return.all.screen            = TRUE,
  epsilon.study                = 0.2,
  p.correction                 = "BH",
  show.pooled.effect           = TRUE,
  return.study.similarity.plot = FALSE,
  test                         = "knha",
  meta.analysis.method         = meta.analysis.method,
  n.cores                      = 5
)

# Extract per-marker screening metrics from result
screening_metrics <- rise_screen_result[["screening.metrics.meta"]]

# Format significant markers into a publication-ready table
screening_metrics_select <- screening_metrics %>%
  arrange(p.unadjusted) %>%
  mutate(mu.delta.ci = paste0(
    round(mu.delta, 3),
    " (",
    round(ci.delta.lower, 3),
    ", ",
    round(ci.delta.upper, 3),
    ")"
  )) %>%
  select(marker, mu.delta.ci, p.unadjusted, p.adjusted) %>%
  filter(p.adjusted < 0.05)

# Render LaTeX table of significant screening markers
kable(
  screening_metrics_select,
  format   = "latex",
  booktabs = TRUE,
  caption  = "Screening Metrics Table"
) %>%
  kable_styling(latex_options = "hold_position") %>%
  row_spec(0, bold = TRUE) %>%
  column_spec(1:ncol(screening_metrics_select))

# Forest plot and rank-correlation fit plot from screening
p1 <- rise_screen_result[["gamma.s.plot"]]$forest.plot
p1

p2 <- rise_screen_result[["gamma.s.plot"]]$fit.plot
p2

p3 = rise_screen_result[["gamma.s.plot"]]$screen.plot
p3

ggsave(
  filename = paste0("TIV_d1_screening_", meta.analysis.method, ".pdf"),
  path     = application_figures_folder,
  plot     = p3,
  width    = 40,
  height   = 18,
  units    = "cm"
)

# ----- Evaluation on test data -----

test_inputs <- extract_rise_inputs(df_test)

# Evaluate significant markers from screening on held-out test data
rise_evaluate_result <- rise.evaluate.meta(
  yone               = test_inputs$yone,
  yzero              = test_inputs$yzero,
  sone               = test_inputs$sone,
  szero              = test_inputs$szero,
  studyone           = test_inputs$studyone,
  studyzero          = test_inputs$studyzero,
  alpha              = 0.05,
  epsilon.meta       = 0.2,
  alternative        = "two.sided",
  paired.all         = TRUE,
  epsilon.study      = 0.2,
  p.correction       = "none",
  show.pooled.effect = TRUE,
  screening.weights  = rise_screen_result[["screening.weights"]],
  markers            = rise_screen_result[["significant.markers"]],
  test               = "knha",
  epsilon.meta.mode  = "user",
  power.want.s.study = 0.8,
  meta.analysis.method = meta.analysis.method
)

# Save evaluation forest plot
p4 <- rise_evaluate_result[["gamma.s.plot"]]$forest.plot

p4

ggsave(
  filename = paste0("TIV_evaluation_forest_", meta.analysis.method, ".pdf"),
  path     = application_figures_folder,
  plot     = p4,
  width    = 32,
  height   = 15,
  units    = "cm"
)

# Save evaluation rank-correlation fit plot
p5 <- rise_evaluate_result[["gamma.s.plot"]]$fit.plot

p5

ggsave(
  filename = paste0("TIV_evaluation_fitplot_", meta.analysis.method, ".pdf"),
  path     = application_figures_folder,
  plot     = p5,
  width    = 37,
  height   = 20,
  units    = "cm"
)
