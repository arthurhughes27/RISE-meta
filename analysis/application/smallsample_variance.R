# R script to illustrate a meta-analysis with RISE
library(tidyverse)
library(SurrogateRank)
library(knitr)
library(kableExtra)

set.seed(08012025)

# Directory containing engineered / processed data files
processed_data_folder <- "data"
# Folder to store images
application_figures_folder = fs::path("output", "figures", "application")

# Paths to processed gene-level data and gene-set objects
p_load_expr_all_noNorm <- fs::path(processed_data_folder, "hipc_merged_all_noNorm.rds")
p_load_BTM <- fs::path(processed_data_folder, "BTM_processed.rds")

# Load data objects
hipc_merged_all_noNorm <- readRDS(p_load_expr_all_noNorm)
BTM <- readRDS(p_load_BTM)

gene_names <- hipc_merged_all_noNorm %>%
  select(a1cf:zzz3) %>%
  select(where(~ !any(is.na(.)))) %>%
  colnames()

tp <- c(1)
timepoints_to_keep <- c(0, tp)

hipc_merged_all_noNorm_filtered <- hipc_merged_all_noNorm %>%
  filter(study_accession == "SDY180") %>%
  mutate(
    response_pre = ifelse(
      study_accession %in% c("SDY80", "SDY180", "SDY1276", "SDY67"),
      immResp_mean_nAb_pre_value,
      immResp_mean_hai_pre_value
    ),
    response_post = ifelse(
      study_accession %in% c("SDY80", "SDY180", "SDY1276", "SDY67"),
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

df_train = hipc_merged_all_noNorm_filtered

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

res = rise.screen(yone, yzero, sone, szero, epsilon = 0.1, paired = TRUE, alpha = 0.05, alternative = "less")
check1 = res[["screening.metrics"]]

res = rise.screen(yone, yzero, sone, szero, epsilon = 0.1, paired = TRUE, alpha = 0.5, alternative = "less",)
check2 = res[["screening.metrics"]]

epsilon = 0.1
paired = TRUE
alpha = 0.5
alternative = "less"
power.want.s = NULL
u.y.hyp = NULL
p.correction = "none"
n.cores = 1
return.all.screen = TRUE
return.all.weights = FALSE
weight.mode = "diff.epsilon"
normalise.weights = TRUE
verbose = TRUE

library(pbmcapply)
library(dplyr)
library(ggplot2)

