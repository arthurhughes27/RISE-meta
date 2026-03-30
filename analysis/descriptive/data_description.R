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

# ===== Tile plot: antibody assay availability per study =====

# Reduce to one row per participant (baseline row) to avoid double-counting
first_tp <- hipc_merged_all_noNorm_filtered %>%
  group_by(participant_id) %>%
  slice_min(study_time_collected, with_ties = FALSE) %>%
  ungroup()

# Flag whether each participant has both pre and post values for each assay
participant_assays <- first_tp %>%
  transmute(
    participant_id,
    study_accession,
    nAb_both = !is.na(immResp_MFC_nAb_pre_value) & !is.na(immResp_MFC_nAb_post_value),
    HAI_both = !is.na(immResp_MFC_hai_pre_value)  & !is.na(immResp_MFC_hai_post_value)
  ) %>%
  pivot_longer(
    cols      = c(nAb_both, HAI_both),
    names_to  = "assay_raw",
    values_to = "has_both"
  ) %>%
  mutate(assay = recode(assay_raw, nAb_both = "nAb", HAI_both = "HAI")) %>%
  select(participant_id, study_accession, assay, has_both)

# Count participants with complete assay data per study
summary_df <- participant_assays %>%
  group_by(study_accession, assay) %>%
  summarise(n_participants = sum(has_both, na.rm = TRUE), .groups = "drop")

# Ensure all study x assay combinations are present so tiles render even when n = 0
summary_df <- expand_grid(
  study_accession = unique(first_tp$study_accession),
  assay           = c("nAb", "HAI")
) %>%
  left_join(summary_df, by = c("study_accession", "assay")) %>%
  mutate(
    n_participants  = replace_na(n_participants, 0),
    # Boolean fill: green if any participants have the assay, white otherwise
    has             = n_participants > 0,
    study_accession = factor(study_accession, levels = study_order),
    assay           = factor(assay, levels = c("nAb", "HAI"))
  )

p2 <- ggplot(summary_df, aes(x = study_accession, y = assay)) +
  geom_tile(aes(fill = has), colour = "grey70") +
  geom_text(aes(label = n_participants), size = 5) +
  scale_fill_manual(values = c("TRUE" = "#2ca25f", "FALSE" = "white"), guide = FALSE) +
  labs(
    x        = "Study accession",
    y        = "Assay",
    title    = "Availability of antibody assays per study",
    subtitle = "Participants with both baseline and day 28 (+- 7 days) post-vaccination measurements"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y   = element_text(size = 13),
    axis.title    = element_text(size = 14),
    panel.grid    = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA)
  )

p2

ggsave(
  filename = "antibody_samples_influenzain.pdf",
  path     = descriptive_figures_folder,
  plot     = p2,
  width    = 30,
  height   = 13,
  units    = "cm"
)

# ===== Stacked bar plot: gender distribution per study =====

# One row per participant at baseline
hipc_merged_all_noNorm_filtered_unique <- hipc_merged_all_noNorm_filtered %>%
  filter(study_time_collected == 0) %>%
  distinct(participant_id, .keep_all = TRUE)

# Proportion and raw count of each gender within each study
summary_df <- hipc_merged_all_noNorm_filtered %>%
  filter(study_time_collected == 0) %>%
  count(study_accession, gender, name = "n") %>%
  group_by(study_accession) %>%
  mutate(total = sum(n), prop = n / total) %>%
  ungroup()

# 100% stacked bars with raw counts annotated inside each segment
p3 <- ggplot(summary_df, aes(x = study_accession, y = n, fill = gender)) +
  geom_col(position = "fill", width = 0.7) +
  geom_text(
    aes(label = n),
    position  = position_fill(vjust = 0.5),
    size      = 5,
    color     = "white",
    fontface  = "bold"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x     = "Study",
    y     = "Proportion",
    title = "Gender distribution by study",
    fill  = "Gender"
  ) +
  theme_bw() +
  theme(
    text            = element_text(size = 16),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y     = element_text(size = 14),
    axis.title      = element_text(size = 16),
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 14)
  )

p3

ggsave(
  filename = "sex_distribution_influenzain.pdf",
  path     = descriptive_figures_folder,
  plot     = p3,
  width    = 30,
  height   = 13,
  units    = "cm"
)

# ===== Violin plot: age distribution per study =====

# Drop participants with missing imputed age
age_df <- hipc_merged_all_noNorm_filtered_unique %>%
  filter(!is.na(age_imputed)) %>%
  select(participant_id, study_accession, age_imputed) %>%
  mutate(study_accession = factor(study_accession, levels = study_order))

# Violin shows the full age distribution; jittered points show individual participants
p4 <- ggplot(age_df, aes(x = study_accession, y = age_imputed)) +
  geom_violin(trim = FALSE, width = 0.9, colour = "black", fill = "#5062FF") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.3, size = 1.8) +
  labs(
    x     = "Study",
    y     = "Age",
    title = "Age distribution by study"
  ) +
  theme_bw() +
  theme(
    text               = element_text(size = 16),
    plot.title         = element_text(hjust = 0.5, face = "bold", size = 20),
    plot.subtitle      = element_text(hjust = 0.5, size = 14),
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y        = element_text(size = 14),
    axis.title         = element_text(size = 16),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

p4

ggsave(
  filename = "age_distribution_influenzain.pdf",
  path     = descriptive_figures_folder,
  plot     = p4,
  width    = 30,
  height   = 13,
  units    = "cm"
)






# Bubble plot of antibody measurement availability

plot_antibody_availability <- function(dat, assay_name) {
  counts <- dat %>%
    filter(assay == assay_name) %>%
    filter(!is.na(study_time_collected)) %>%
    filter(!is.na(value_preferred)) %>%
    distinct(participant_id, study_accession, study_time_collected, response_strain_analyte) %>%
    group_by(study_accession, study_time_collected) %>%
    summarise(n = n_distinct(participant_id), .groups = "drop")
  
  time_levels <- counts %>%
    distinct(study_time_collected) %>%
    mutate(time_num = suppressWarnings(as.numeric(as.character(study_time_collected)))) %>%
    arrange(
      ifelse(is.na(time_num), Inf, time_num),
      as.character(study_time_collected)
    ) %>%
    pull(study_time_collected) %>%
    as.character()
  
  counts <- counts %>%
    mutate(
      study_time_collected = factor(as.character(study_time_collected),
                                    levels = time_levels, ordered = TRUE),
      study_accession = factor(study_accession,
                               levels = rev(levels(factor(study_accession))))
    )
  
  ggplot(counts, aes(x = study_time_collected, y = study_accession)) +
    geom_point(
      aes(size = n),
      fill = "#5062FF",
      shape = 21,
      colour = "black",
      alpha = 0.8,
      show.legend = FALSE
    ) +
    geom_text(aes(label = n), colour = "white", size = 3.5, vjust = 0.5) +
    scale_fill_identity(guide = "none") +
    scale_size_area(max_size = 28, guide = "none") +
    labs(
      x = "Days post-vaccination",
      y = "Study",
      title = paste0("Availability of ", assay_name, " antibody measurements across time"),
      subtitle = "Number of unique individuals with at least one non-missing value_preferred"
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
}

p5 <- plot_antibody_availability(raw_response_influenzain, "hai")
p6 <- plot_antibody_availability(raw_response_influenzain, "nAb")

p5
p6

ggsave(
  filename = "hai_bubble_plot.pdf",
  path     = descriptive_figures_folder,
  plot     = p5,
  width    = 30,
  height   = 20,
  units    = "cm"
)

ggsave(
  filename = "nAb_bubble_plot.pdf",
  path     = descriptive_figures_folder,
  plot     = p6,
  width    = 30,
  height   = 20,
  units    = "cm"
)



counts_combined <- raw_response_influenzain %>%
  filter(assay %in% c("hai", "nAb")) %>%
  filter(!is.na(study_time_collected)) %>%
  filter(!is.na(value_preferred)) %>%
  distinct(participant_id, study_accession, study_time_collected) %>%
  group_by(study_accession, study_time_collected) %>%
  summarise(n = n_distinct(participant_id), .groups = "drop")

time_levels <- counts_combined %>%
  distinct(study_time_collected) %>%
  mutate(time_num = suppressWarnings(as.numeric(as.character(study_time_collected)))) %>%
  arrange(ifelse(is.na(time_num), Inf, time_num), as.character(study_time_collected)) %>%
  pull(study_time_collected) %>%
  as.character()

study_levels <- counts_combined %>%
  distinct(study_accession) %>%
  pull(study_accession) %>%
  as.character()

counts_combined <- counts_combined %>%
  mutate(
    study_time_collected = factor(as.character(study_time_collected),
                                  levels = time_levels, ordered = TRUE),
    study_accession = factor(as.character(study_accession),
                             levels = rev(study_levels))
  )

p7 <- ggplot(counts_combined, aes(x = study_time_collected, y = study_accession)) +
  geom_point(
    aes(size = n),
    fill = "#5062FF",
    shape = 21,
    colour = "black",
    alpha = 0.8,
    show.legend = FALSE
  ) +
  geom_text(aes(label = n), colour = "white", size = 3.5, vjust = 0.5) +
  scale_fill_identity(guide = "none") +
  scale_size_area(max_size = 28, guide = "none") +
  labs(
    x = "Days post-vaccination",
    y = "Study",
    title = "Availability of antibody measurements across time",
    subtitle = "Unique participants with at least one non-missing measurement in either hai or nAb"
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

p7

ggsave(
  filename = "antibody_bubble_plot.pdf",
  path     = descriptive_figures_folder,
  plot     = p7,
  width    = 30,
  height   = 25,
  units    = "cm"
)
