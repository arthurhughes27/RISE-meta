# File for preprocessing the raw immune response data files
## In particular, we are interested in deriving "maximum fold-change in immune response" (MFC) values in the first 28 (+-7) days after vaccination
## We use the fold-change to take into account baseline immunity
## We use the maximum of fold-change values to choose one measurement across multiple analytes, timepoints, or assays
## In this script, the final result will be an immune response dataframe with one row per unique participant identifier
## We will have one column for the MFC values for each of three assay types (neutralising antibodies, ELISA, and HAI)
## We will also create a column giving the maximum of any of the assay types

# Load packages
library(fs)
library(dplyr)
library(stringr)
library(janitor)
library(readxl)
library(purrr)


# Specify folder within folder root where the raw data lives
raw_data_folder = "data-raw"
processed_data_folder = "data"

# Use fs::path() to specify the data paths robustly
p_load_nAb <- fs::dir_ls(raw_data_folder, regexp = "neut_ab_titer.*\\.xlsx$")
p_load_hai <- fs::dir_ls(raw_data_folder, regexp = "hai.*\\.xlsx$")
p_load_clinical <- fs::path(processed_data_folder, "hipc_clinical.rds")

# HAI response data
response_hai = read_excel(p_load_hai) %>%
  clean_names()

# nAb response data
response_nAb = read_excel(p_load_nAb) %>%
  clean_names()

# Clinical data for filtration of participants
hipc_clinical_all_noNorm = readRDS(p_load_clinical)

# First, filter each dataframe to only contain information on participants for which we have gene expression data
# Find identifiers of participants with gene expression measurements
participants = hipc_clinical_all_noNorm %>%
  filter(vaccine_name %in% c("Influenza (IN)", "Influenza (LV)")) %>% 
  pull(participant_id) %>%
  unique()

# Filter immune response data by these participants
response_hai = response_hai %>%
  filter(participant_id %in% participants)

response_nAb = response_nAb %>%
  filter(participant_id %in% participants)

# Depending on the assay, there may be multiple viral strains or analytes measured.
# We rename some columns and add a column for the assay name, so that we can merge these data together.

response_hai = response_hai %>%
  rename(response_strain_analyte = virus) %>%
  mutate(assay = "hai")

response_nAb = response_nAb %>%
  rename(response_strain_analyte = virus) %>%
  mutate(assay = "nAb")

# Now merge all the raw response data
raw_response_influenzain = bind_rows(response_nAb, response_hai) %>%
  select(-gender) %>% 
  arrange(participant_id) %>%
  distinct()

# First get the studies and other clinical data corresponding to each participant id
hipc_studies = hipc_clinical_all_noNorm %>%
  select(participant_id,
         age_imputed,
         gender,
         study_accession,
         vaccine,
         vaccine_type,
         pathogen) %>%
  distinct()

# Merge the study names into the immune response data (it is not directly given)
raw_response_influenzain = merge(x = raw_response_influenzain,
                                 y = hipc_studies,
                                 by = "participant_id",
                                 all = F) %>%
  select(
    participant_id,
    study_accession,
    age_imputed, 
    gender,
    response_strain_analyte,
    assay,
    study_time_collected,
    value_preferred
  )

# Use fs::path() to specify the data path robustly
p_save <- fs::path(processed_data_folder, "raw_response_influenzain_all_noNorm.rds")

# Save dataframe
saveRDS(raw_response_influenzain, file = p_save)

# Now derive maximum fold-change values for each individual

# We intend to predict on immune response data at day 28 (plus or minus 7 days), so we filter to get the immune data at these timepoints or pre-vaccination.

raw_response_influenzain = raw_response_influenzain %>%
  filter((study_time_collected < 36 &
            study_time_collected > 20) | study_time_collected <= 0)

# 2) For each participant × response_strain_analyte × assay, keep one pre measurement:
#    - choose the measurement closest to day 0
#    - because all pre values are <= 0, the closest one is the largest study_time_collected
response_pre <- raw_response_influenzain %>%
  filter(study_time_collected <= 0) %>%
  group_by(participant_id, response_strain_analyte, assay) %>%
  slice_max(order_by = study_time_collected, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(study_time_collected = 0)

# 3) For each participant × response_strain_analyte × assay, keep one post measurement:
#    - choose the measurement closest to day 28
#    - if there is a tie, this keeps the earlier timepoint because of the arrange() order
response_post <- raw_response_influenzain %>%
  filter(between(study_time_collected, 21, 35)) %>%
  group_by(participant_id, response_strain_analyte, assay) %>%
  arrange(abs(study_time_collected - 28), study_time_collected, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(study_time_collected = 28)

# 4) Combine the selected pre and post measurements
raw_response_influenzain <- bind_rows(response_pre, response_post) %>%
  arrange(participant_id, response_strain_analyte, assay, study_time_collected)

# Now derive MFC values, taking most recent measurement as baseline in the case where there are two pre-vaccination measurements

pre_vax <- raw_response_influenzain %>%
  filter(study_time_collected <= 0) %>%
  group_by(participant_id, assay, response_strain_analyte) %>%
  slice_max(study_time_collected, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(response_MFC_pre_value = value_preferred,
         response_MFC_pre_time = study_time_collected) %>% 
  select(participant_id, study_accession, assay, response_strain_analyte, response_MFC_pre_time, response_MFC_pre_value)

post_vax <- raw_response_influenzain %>%
  filter(study_time_collected > 0) %>%
  rename(response_MFC_post_value = value_preferred,
         response_MFC_post_time = study_time_collected) %>% 
  select(participant_id, study_accession, assay, response_strain_analyte, response_MFC_post_time, response_MFC_post_value)

# Now we merge pre- and post- values to get a dataframe in wide format with one line per unique combination of participant/assay/strain-analyte
merged_vax <- full_join(
  pre_vax,
  post_vax,
  by = c("participant_id", "study_accession","response_strain_analyte", "assay"),
  relationship = "many-to-many"
) %>%
  select(
    participant_id,
    study_accession,
    assay,
    response_strain_analyte,
    response_MFC_pre_time,
    response_MFC_post_time,
    response_MFC_pre_value,
    response_MFC_post_value
  )

# Now we can exclude any rows with NA values in either pre- or post-.
merged_vax = merged_vax %>%
  filter(!is.na(response_MFC_pre_value),
         !is.na(response_MFC_post_value))

# Derive the log2 fold changes.
response_MFC_df <- merged_vax %>%
  mutate(
    # Define fold change
    response_MFC = response_MFC_post_value / response_MFC_pre_value,
    # If pre-vaccination value is 0, set fold change to post-vaccination value
    response_MFC = ifelse(
      is.infinite(response_MFC),
      response_MFC_post_value,
      response_MFC
    ),
    # Set fold change to 1 if both response_MFC_pre_value and response_MFC_post_value are 0
    response_MFC = ifelse(
      response_MFC_pre_value == 0 &
        response_MFC_post_value == 0,
      1,
      response_MFC
    )
  ) %>%
  # Calculate log2 fold change
  mutate(
    response_log2_MFC = ifelse(is.infinite(log2(response_MFC)), NA, log2(response_MFC)),
    response_log2_MFC_pre_value = log2(response_MFC_pre_value),
    response_log2_MFC_post_value = log2(response_MFC_post_value)
  ) %>%
  # Select relevant columns
  select(
    participant_id,
    study_accession,
    assay,
    response_strain_analyte,
    response_MFC_pre_time,
    response_MFC_post_time,
    response_MFC_pre_value,
    response_MFC_post_value,
    response_MFC, 
    response_log2_MFC_pre_value,
    response_log2_MFC_post_value,
    response_log2_MFC
  )



# Now, for each participant, we can derive the maximum of these fold changes.
## Across all assays
max_response_MFC_df_any <- response_MFC_df %>%
  group_by(participant_id) %>%
  slice_max(response_MFC, n = 1, with_ties = F) %>%
  ungroup() %>%
  rename(
    immResp_MFC_anyAssay_response_strain_analyte = response_strain_analyte,
    immResp_MFC_anyAssay_assay = assay,
    immResp_MFC_anyAssay_pre_time = response_MFC_pre_time,
    immResp_MFC_anyAssay_post_time = response_MFC_post_time,
    immResp_MFC_anyAssay_pre_value = response_MFC_pre_value,
    immResp_MFC_anyAssay_post_value = response_MFC_post_value,
    immResp_MFC_anyAssay_MFC = response_MFC,
    immResp_MFC_anyAssay_log2_pre_value = response_log2_MFC_pre_value,
    immResp_MFC_anyAssay_log2_post_value = response_log2_MFC_post_value,
    immResp_MFC_anyAssay_log2_MFC = response_log2_MFC
  )

## And within each assay
max_response_MFC_df_each <- response_MFC_df %>%
  group_by(participant_id, assay) %>%
  slice_max(response_MFC, n = 1, with_ties = F) %>%
  ungroup()

max_response_MFC_df_nAb = max_response_MFC_df_each %>%
  filter(assay == "nAb") %>%
  rename(
    immResp_MFC_nAb_response_strain_analyte = response_strain_analyte,
    immResp_MFC_nAb_pre_time = response_MFC_pre_time,
    immResp_MFC_nAb_post_time = response_MFC_post_time,
    immResp_MFC_nAb_pre_value = response_MFC_pre_value,
    immResp_MFC_nAb_post_value = response_MFC_post_value,
    immResp_MFC_nAb_MFC = response_MFC,
    immResp_MFC_nAb_log2_pre_value = response_log2_MFC_pre_value,
    immResp_MFC_nAb_log2_post_value = response_log2_MFC_post_value,
    immResp_MFC_nAb_log2_MFC = response_log2_MFC
  ) %>%
  select(-assay)

max_response_MFC_df_hai = max_response_MFC_df_each %>%
  filter(assay == "hai") %>%
  rename(
    immResp_MFC_hai_response_strain_analyte = response_strain_analyte,
    immResp_MFC_hai_pre_time = response_MFC_pre_time,
    immResp_MFC_hai_post_time = response_MFC_post_time,
    immResp_MFC_hai_pre_value = response_MFC_pre_value,
    immResp_MFC_hai_post_value = response_MFC_post_value,
    immResp_MFC_hai_MFC = response_MFC,
    immResp_MFC_hai_log2_pre_value = response_log2_MFC_pre_value,
    immResp_MFC_hai_log2_post_value = response_log2_MFC_post_value,
    immResp_MFC_hai_log2_MFC = response_log2_MFC
  ) %>%
  select(-assay)

mean_response_nAb <- response_MFC_df %>%
  filter(assay == "nAb") %>%
  group_by(participant_id, study_accession) %>%
  summarise(immResp_mean_nAb_pre_value  = mean(response_MFC_pre_value,  na.rm = TRUE),
            immResp_mean_nAb_post_value = mean(response_MFC_post_value, na.rm = TRUE),
            .groups = "drop") %>% 
  mutate(immResp_mean_nAb_log2_pre_value = log2(immResp_mean_nAb_pre_value),
         immResp_mean_nAb_log2_post_value = log2(immResp_mean_nAb_post_value),
         immResp_mean_nAb_FC = immResp_mean_nAb_post_value/immResp_mean_nAb_pre_value,
         immResp_mean_nAb_log2_FC = log2(immResp_mean_nAb_FC))

mean_response_hai <- response_MFC_df %>%
  filter(assay == "hai") %>%
  group_by(participant_id, study_accession) %>%
  summarise(immResp_mean_hai_pre_value  = mean(response_MFC_pre_value,  na.rm = TRUE),
            immResp_mean_hai_post_value = mean(response_MFC_post_value, na.rm = TRUE),
            .groups = "drop") %>% 
  mutate(immResp_mean_hai_log2_pre_value = log2(immResp_mean_hai_pre_value),
         immResp_mean_hai_log2_post_value = log2(immResp_mean_hai_post_value),
         immResp_mean_hai_FC = immResp_mean_hai_post_value/immResp_mean_hai_pre_value,
         immResp_mean_hai_log2_FC = log2(immResp_mean_hai_FC))

hipc_immResp <- list(
  max_response_MFC_df_any,
  max_response_MFC_df_nAb,
  max_response_MFC_df_hai,
  mean_response_nAb,
  mean_response_hai
) %>%
  reduce(full_join, by = c("participant_id", "study_accession")) %>%
  arrange(participant_id)

# ---------------------------
# Compute strain-wise z-scores, then participant-level mean/max of z-scores
# Produces columns:
# immResp_mean_nAb_std_pre_value, immResp_mean_nAb_std_post_value,
# immResp_max_nAb_std_pre_value,  immResp_max_nAb_std_post_value,
# immResp_mean_hai_std_pre_value, immResp_mean_hai_std_post_value,
# immResp_max_hai_std_pre_value,  immResp_max_hai_std_post_value
# ---------------------------
library(tidyr)

# 1) pivot pre/post into long time column and z-score within (study x assay x strain x time)
std_long <- response_MFC_df %>%
  select(participant_id, study_accession, assay, response_strain_analyte,
         response_MFC_pre_value, response_MFC_post_value) %>%
  pivot_longer(
    cols = c(response_MFC_pre_value, response_MFC_post_value),
    names_to = "time",
    values_to = "value"
  ) %>%
  mutate(time = ifelse(time == "response_MFC_pre_value", "pre", "post")) %>%
  group_by(study_accession, assay, response_strain_analyte, time) %>%
  mutate(
    grp_mean = mean(value, na.rm = TRUE),
    grp_sd   = sd(value, na.rm = TRUE),
    value_z  = ifelse(is.na(value), NA_real_,
                      ifelse(is.na(grp_sd) | grp_sd == 0, 0, (value - grp_mean) / grp_sd))
  ) %>%
  ungroup() %>%
  select(participant_id, study_accession, assay, response_strain_analyte, time, value_z)

# 2) summarise per participant/study/assay/time: mean and max of z-scores across strains
std_summary <- std_long %>%
  group_by(participant_id, study_accession, assay, time) %>%
  summarise(
    mean_std = if(all(is.na(value_z))) NA_real_ else mean(value_z, na.rm = TRUE),
    max_std  = if(all(is.na(value_z))) NA_real_ else max(value_z, na.rm = TRUE),
    .groups = "drop"
  )

# 3) split into separate tibbles for each (transformation x assay x time) with requested names
immResp_mean_nAb_std_pre  <- std_summary %>%
  filter(assay == "nAb", time == "pre") %>%
  rename(immResp_mean_nAb_std_pre_value  = mean_std)

immResp_mean_nAb_std_post <- std_summary %>%
  filter(assay == "nAb", time == "post") %>%
  rename(immResp_mean_nAb_std_post_value = mean_std)

immResp_max_nAb_std_pre   <- std_summary %>%
  filter(assay == "nAb", time == "pre") %>%
  rename(immResp_max_nAb_std_pre_value   = max_std)

immResp_max_nAb_std_post  <- std_summary %>%
  filter(assay == "nAb", time == "post") %>%
  rename(immResp_max_nAb_std_post_value  = max_std)

immResp_mean_hai_std_pre  <- std_summary %>%
  filter(assay == "hai", time == "pre") %>%
  rename(immResp_mean_hai_std_pre_value  = mean_std)

immResp_mean_hai_std_post <- std_summary %>%
  filter(assay == "hai", time == "post") %>%
  rename(immResp_mean_hai_std_post_value = mean_std)

immResp_max_hai_std_pre   <- std_summary %>%
  filter(assay == "hai", time == "pre") %>%
  rename(immResp_max_hai_std_pre_value   = max_std)

immResp_max_hai_std_post  <- std_summary %>%
  filter(assay == "hai", time == "post") %>%
  rename(immResp_max_hai_std_post_value  = max_std)

# 4) join the eight tibbles into one table keyed by participant_id + study_accession,
#    keeping only the final eight value columns plus keys
std_cols <- list(
  immResp_mean_nAb_std_pre,
  immResp_mean_nAb_std_post,
  immResp_max_nAb_std_pre,
  immResp_max_nAb_std_post,
  immResp_mean_hai_std_pre,
  immResp_mean_hai_std_post,
  immResp_max_hai_std_pre,
  immResp_max_hai_std_post
) %>%
  purrr::reduce(function(x, y) full_join(x, y, by = c("participant_id", "study_accession"))) %>%
  select(
    participant_id,
    study_accession,
    immResp_mean_nAb_std_pre_value,
    immResp_mean_nAb_std_post_value,
    immResp_max_nAb_std_pre_value,
    immResp_max_nAb_std_post_value,
    immResp_mean_hai_std_pre_value,
    immResp_mean_hai_std_post_value,
    immResp_max_hai_std_pre_value,
    immResp_max_hai_std_post_value
  )

# 5) merge into hipc_immResp (preserve existing cols; append only the 8 standardized columns)
hipc_immResp <- hipc_immResp %>%
  left_join(std_cols, by = c("participant_id", "study_accession")) %>%
  arrange(participant_id)


# Use fs::path() to specify the data path robustly
p_save <- fs::path(processed_data_folder, "hipc_immResp.rds")

# Save dataframe
saveRDS(hipc_immResp, file = p_save)

rm(list = ls())
