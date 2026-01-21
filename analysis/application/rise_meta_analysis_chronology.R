# R script to illustrate a meta-analysis with RISE
library(tidyverse)
library(SurrogateRank)
library(grid)    # for unit()
library(scales)  # pretty formatting
library(cowplot)
library(metafor) # required for REML random-effects fit

set.seed(08012025)

# Directory containing engineered / processed data files
processed_data_folder <- "data"

# Paths to processed gene-level data and gene-set objects
p_load_expr_all_noNorm <- fs::path(processed_data_folder, "hipc_merged_all_noNorm.rds")

# Load data objects
hipc_merged_all_noNorm <- readRDS(p_load_expr_all_noNorm)

gene_names <- hipc_merged_all_noNorm %>%
  select(a1cf:zzz3) %>%
  select(where(~ !any(is.na(.)))) %>%
  colnames()

tp <- 2
timepoints_to_keep <- c(0, tp)

hipc_merged_all_noNorm_filtered <- hipc_merged_all_noNorm %>%
  mutate(response_pre = ifelse(study_accession == "SDY80", immResp_MFC_nAb_pre_value, immResp_MFC_hai_pre_value),
         response_post = ifelse(study_accession == "SDY80", immResp_MFC_nAb_post_value, immResp_MFC_hai_post_value)) %>% 
  filter(
    !is.na(immResp_MFC_anyAssay_log2_MFC),
    vaccine_name == "Influenza (IN)",
    study_time_collected %in% timepoints_to_keep,
    !(study_accession %in% c("SDY61", "SDY224"))
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
    study_accession,
    study_time_collected,
    response_pre,
    response_post,
    all_of(gene_names)
  ) %>%
  arrange(participant_id)

# test_studies = c("SDY224")

test_studies = hipc_merged_all_noNorm_filtered %>%
  pull(study_accession) %>%
  unique() %>%
  sample(size = 3)

df_train = hipc_merged_all_noNorm_filtered %>%
  # filter(
  #   # !(study_accession %in% test_studies),
  #        gender == "Male") %>% 
  group_by(study_accession) %>% 
  filter(length(unique(participant_id)) > 8) %>% 
  ungroup()

yone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(response_post)

yzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(response_pre)

sone = df_train %>%
  filter(study_time_collected > 0) %>%
  select(all_of(gene_names))

szero = df_train %>%
  filter(study_time_collected == 0) %>%
  select(all_of(gene_names))

studyone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession)

studyzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(study_accession)

rise.screen.meta.result = rise.screen.meta(
  yone,
  yzero,
  sone,
  szero,
  studyone,
  studyzero,
  alpha = 0.05,
  epsilon.study = 0.2,
  epsilon.meta = 0.2,
  p.correction = "none",
  n.cores = 12,
  alternative = "two.sided",
  paired.all = T,
  return.all.screen = F,
  return.all.weights = F,
  weight.mode = "diff.epsilon"
)

markers = rise.screen.meta.result[["significant.markers"]]

screening.weights = rise.screen.meta.result[["screening.weights"]]

df_test = hipc_merged_all_noNorm_filtered %>%
  filter(study_accession %in% test_studies)

yone = df_test %>%
  filter(study_time_collected > 0) %>%
  pull(response_post)

yzero = df_test %>%
  filter(study_time_collected == 0) %>%
  pull(response_pre)

sone = df_test %>%
  filter(study_time_collected > 0) %>%
  select(all_of(gene_names))

szero = df_test %>%
  filter(study_time_collected == 0) %>%
  select(all_of(gene_names))

studyone = df_test %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession)

studyzero = df_test %>%
  filter(study_time_collected == 0) %>%
  pull(study_accession)

rise.evaluate.meta.result = rise.evaluate.meta(
  yone,
  yzero,
  sone,
  szero,
  studyone,
  studyzero,
  alpha = 0.05,
  epsilon.meta = 0.2,
  alternative = "two.sided",
  paired.all = T,
  evaluate.weights = T,
  return.all.evaluate = T,
  epsilon.study = 0.2,
  p.correction = "BH",
  screening.weights = screening.weights,
  markers = markers
)

rise.screen.meta.result$gamma.s.plot
rise.evaluate.meta.result$gamma.s.plot

