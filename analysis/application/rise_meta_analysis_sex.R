# R script to illustrate a meta-analysis with RISE
library(tidyverse)
library(SurrogateRank)
library(grid)    # for unit()
library(scales)  # pretty formatting
library(cowplot)
library(metafor) # required for REML random-effects fit

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

tp <- 1
timepoints_to_keep <- c(0, tp)

hipc_merged_all_noNorm_filtered <- hipc_merged_all_noNorm %>%
  filter(
    !is.na(immResp_MFC_anyAssay_log2_MFC),
    vaccine_name == "Influenza (IN)",
    study_time_collected %in% timepoints_to_keep
  ) %>%
  group_by(participant_id) %>%
  filter(sum(study_time_collected == 0) == 1,
         sum(study_time_collected == tp) == 1) %>%
  ungroup() %>%
  select(
    participant_id,
    gender,
    study_accession,
    study_time_collected,
    immResp_MFC_anyAssay_pre_value,
    immResp_MFC_anyAssay_post_value,
    all_of(gene_names)
  ) %>%
  arrange(participant_id)

df_train = hipc_merged_all_noNorm_filtered %>%
  filter(gender == "Male")

yone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(immResp_MFC_anyAssay_post_value)

yzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(immResp_MFC_anyAssay_pre_value)

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
  filter(gender == "Female")

yone = df_test %>%
  filter(study_time_collected > 0) %>%
  pull(immResp_MFC_anyAssay_post_value)

yzero = df_test %>%
  filter(study_time_collected == 0) %>%
  pull(immResp_MFC_anyAssay_pre_value)

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
  alpha = 0.1,
  epsilon.meta = 0.2,
  alternative = "two.sided",
  paired.all = T,
  evaluate.weights = T,
  screening.weights = screening.weights,
  markers = markers
)
