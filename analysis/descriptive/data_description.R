# Descriptive summaries and figures for the HIPC influenza vaccine cohorts
library(tidyverse)
library(SurrogateRank)
library(grid)
library(scales)
library(cowplot)
library(metafor)
library(forcats)
library(stringr)

# Paths to processed data and output figures
processed_data_folder <- "data"
descriptive_figures_folder <- fs::path("output", "figures", "descriptive")

# Load merged gene expression dataset
hipc_merged_all_noNorm <- readRDS(fs::path(processed_data_folder, "hipc_merged_all_noNorm.rds"))
raw_response_influenzain = readRDS(fs::path(processed_data_folder, "raw_response_influenzain_all_noNorm.rds"))

# Gene columns present in the data with no missing values
gene_names <- hipc_merged_all_noNorm %>%
  select(a1cf:zzz3) %>%
  select(where(~ !any(is.na(.)))) %>%
  colnames()

# Inclusion criteria: inactivated influenza vaccine, MFC available, observations
# at days 0-7 only, exactly one baseline row and at least one post-vaccination row
hipc_merged_all_noNorm_filtered <- hipc_merged_all_noNorm %>%
  filter(
    !is.na(immResp_MFC_anyAssay_log2_MFC),
    vaccine_name == "Influenza (IN)",
    study_time_collected %in% c(0, 1, 2, 3, 7)
  ) %>%
  group_by(participant_id) %>%
  filter(
    sum(study_time_collected == 0) == 1,
    sum(study_time_collected > 0) > 0
  ) %>%
  ungroup() %>%
  select(
    participant_id,
    study_accession,
    arm_accession,
    age_imputed,
    study_time_collected,
    immResp_MFC_nAb_pre_value,
    immResp_MFC_nAb_post_value,
    immResp_MFC_hai_pre_value,
    immResp_MFC_hai_post_value,
    immResp_MFC_anyAssay_pre_value,
    immResp_MFC_anyAssay_post_value,
    gender,
    all_of(gene_names)
  ) %>%
  mutate(gender = droplevels(gender)) %>%
  arrange(participant_id)

# Canonical study order used consistently across all plots
study_order <- c(
  "SDY61",
  "SDY269",
  "SDY270",
  "SDY180",
  "SDY56",
  "SDY67",
  "SDY1119",
  "SDY224",
  "SDY404",
  "SDY63",
  "SDY400",
  "SDY1276",
  "SDY520",
  "SDY80",
  "SDY640"
)

# Apply canonical study ordering
hipc_merged_all_noNorm_filtered <- hipc_merged_all_noNorm_filtered %>%
  mutate(study_accession = factor(study_accession, levels = study_order))

# ===== Bubble plot: transcriptomic sample counts by study and timepoint =====

# Count observations per study and timepoint
counts <- hipc_merged_all_noNorm_filtered %>%
  filter(!is.na(study_time_collected)) %>%
  group_by(study_accession, study_time_collected) %>%
  summarise(n = n(), .groups = "drop")

# Convert timepoint to an ordered factor for correct x-axis ordering
time_levels <- counts %>%
  distinct(study_time_collected) %>%
  arrange(as.numeric(study_time_collected)) %>%
  pull(study_time_collected) %>%
  as.character()

counts <- counts %>%
  mutate(
    study_time_collected = factor(as.character(study_time_collected), levels = time_levels, ordered = TRUE),
    # Reverse study order so the first study appears at the top of the y-axis
    study_accession = factor(study_accession, levels = rev(levels(study_accession)))
  )

# Bubble size is proportional to raw count; labels show the count inside each bubble
p1 <- ggplot(counts, aes(x = study_time_collected, y = study_accession)) +
  geom_point(
    aes(size = n),
    fill = "#5062FF",
    shape = 21,
    colour = "black",
    alpha = 0.8,
    show.legend = FALSE
  ) +
  geom_text(aes(label = n), colour = "white", size = 3.5, vjust = 0.5, show.legend = TRUE) +
  scale_fill_identity(guide = "none") +
  scale_size_area(max_size = 28, guide = "none") +
  labs(
    x        = "Days post-vaccination",
    y        = "Study",
    title    = "Number of transcriptomic observations across time per study",
    subtitle = "Selected timepoints, healthy adults with baseline and post-vaccination data"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor   = element_blank(),
    axis.title         = element_text(size = 17),
    axis.text          = element_text(size = 15),
    plot.title         = element_text(size = 20, hjust = 0.5, face = "bold"),
    plot.subtitle      = element_text(size = 17, hjust = 0.5)
  )

p1

ggsave(
  filename = "transcriptomic_samples_influenzain.pdf",
  path     = descriptive_figures_folder,
  plot     = p1,
  width    = 30,
  height   = 26,
  units    = "cm"
)