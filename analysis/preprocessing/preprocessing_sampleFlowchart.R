# R Script to extract sample and study sizes for visual representation of sample filtering

library(dplyr)

# Specify folder within folder root where the processed data lives
processed_data_folder = "data"

# Use fs::path() to specify the data path robustly
p_load <- fs::path(processed_data_folder, "hipc_merged_all_noNorm.rds")

# Load the merged data
hipc_merged_all_noNorm = readRDS(p_load)

# Base participants, studies
n1 = hipc_merged_all_noNorm %>%
  pull(participant_id) %>%
  unique() %>%
  length()

M1 = hipc_merged_all_noNorm %>%
  pull(study_accession) %>%
  unique() %>%
  length()

# Inactivated influenza only
hipc_merged_all_noNorm = hipc_merged_all_noNorm %>%
  filter(vaccine_name == "Influenza (IN)")

n2 = hipc_merged_all_noNorm %>%
  pull(participant_id) %>%
  unique() %>%
  length()

M2 = hipc_merged_all_noNorm %>%
  pull(study_accession) %>%
  unique() %>%
  length()

# Timepoints 0 and 1 only
hipc_merged_all_noNorm <- hipc_merged_all_noNorm %>%
  filter(study_time_collected %in% c(0, 1)) %>%
  group_by(participant_id) %>%
  filter(n_distinct(study_time_collected) == 2) %>%
  ungroup()

n3 = hipc_merged_all_noNorm %>%
  pull(participant_id) %>%
  unique() %>%
  length()

M3 = hipc_merged_all_noNorm %>%
  pull(study_accession) %>%
  unique() %>%
  length()

# Valid immune response at both pre and post for at least one of HAI or nAb
hipc_merged_all_noNorm <- hipc_merged_all_noNorm %>%
  filter((
    !is.na(immResp_mean_hai_pre_value) &
      !is.na(immResp_mean_hai_post_value)
  ) |
    (
      !is.na(immResp_mean_nAb_pre_value) &
        !is.na(immResp_mean_nAb_post_value)
    ))

n4 = hipc_merged_all_noNorm %>%
  pull(participant_id) %>%
  unique() %>%
  length()

M4 = hipc_merged_all_noNorm %>%
  pull(study_accession) %>%
  unique() %>%
  length()

# Remove studies with fewer than 6 total participants
hipc_merged_all_noNorm <- hipc_merged_all_noNorm %>%
  group_by(study_accession) %>%
  filter(length(unique(participant_id)) > 5) %>%
  ungroup()

n5 = hipc_merged_all_noNorm %>%
  pull(participant_id) %>%
  unique() %>%
  length()

M5 = hipc_merged_all_noNorm %>%
  pull(study_accession) %>%
  unique() %>%
  length()

n1
n2
n3
n4
n5

M1
M2
M3
M4
M5

n1-n2
n2-n3
n3-n4
n4-n5
