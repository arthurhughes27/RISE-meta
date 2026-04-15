library(dplyr)
library(fs)

processed_data_folder <- "data"
results_folder = fs::path("output", "results", "descriptive")
p_load <- fs::path(processed_data_folder, "hipc_merged_all_noNorm.rds")

hipc_merged_all_noNorm <- readRDS(p_load)

# Helper to count unique individuals and studies
get_counts <- function(df) {
  tibble(
    individuals = n_distinct(df$participant_id),
    studies = n_distinct(df$study_accession)
  )
}

# Initialize tracking table
summary_df <- tibble(
  step = character(),
  individuals = integer(),
  difference_in_individuals = integer(),
  studies = integer()
)

# Step 1: starting dataset
counts1 <- get_counts(hipc_merged_all_noNorm)
summary_df <- bind_rows(
  summary_df,
  tibble(
    step = "Starting dataset",
    individuals = counts1$individuals,
    difference_in_individuals = NA_integer_,
    studies = counts1$studies
  )
)

# Step 2: inactivated influenza only
prev_n <- counts1$individuals

hipc_merged_all_noNorm <- hipc_merged_all_noNorm %>%
  filter(vaccine_name == "Influenza (IN)")

counts2 <- get_counts(hipc_merged_all_noNorm)
summary_df <- bind_rows(
  summary_df,
  tibble(
    step = "Filtered for inactivated influenza vaccine only",
    individuals = counts2$individuals,
    difference_in_individuals = prev_n - counts2$individuals,
    studies = counts2$studies
  )
)

# Step 3: timepoints 0 and 1 only, keeping participants with both
prev_n <- counts2$individuals

hipc_merged_all_noNorm <- hipc_merged_all_noNorm %>%
  filter(study_time_collected %in% c(0, 1)) %>%
  group_by(participant_id) %>%
  filter(n_distinct(study_time_collected) == 2) %>%
  ungroup()

counts3 <- get_counts(hipc_merged_all_noNorm)
summary_df <- bind_rows(
  summary_df,
  tibble(
    step = "Filtered for paired observations at timepoints 0 and 1",
    individuals = counts3$individuals,
    difference_in_individuals = prev_n - counts3$individuals,
    studies = counts3$studies
  )
)

# Step 4: valid immune response at both pre and post for at least one of HAI or nAb
prev_n <- counts3$individuals

hipc_merged_all_noNorm <- hipc_merged_all_noNorm %>%
  filter(
    (!is.na(immResp_mean_hai_pre_value) & !is.na(immResp_mean_hai_post_value)) |
      (!is.na(immResp_mean_nAb_pre_value) & !is.na(immResp_mean_nAb_post_value))
  )

counts4 <- get_counts(hipc_merged_all_noNorm)
summary_df <- bind_rows(
  summary_df,
  tibble(
    step = "Filtered for valid immune responses",
    individuals = counts4$individuals,
    difference_in_individuals = prev_n - counts4$individuals,
    studies = counts4$studies
  )
)

# Step 5: remove studies with fewer than 6 total participants
prev_n <- counts4$individuals

hipc_merged_all_noNorm <- hipc_merged_all_noNorm %>%
  group_by(study_accession) %>%
  filter(n_distinct(participant_id) > 5) %>%
  ungroup()

counts5 <- get_counts(hipc_merged_all_noNorm)
summary_df <- bind_rows(
  summary_df,
  tibble(
    step = "Filtered for studies with at least 6 participants",
    individuals = counts5$individuals,
    difference_in_individuals = prev_n - counts5$individuals,
    studies = counts5$studies
  )
)

# Write CSV
out_file <- fs::path(results_folder, "sampleFlowchart.csv")
write.csv(summary_df, out_file, row.names = FALSE)

summary_df

rm(list = ls())