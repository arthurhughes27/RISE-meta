# R script to illustrate a meta-analysis with RISE
library(tidyverse)
library(SurrogateRank)
library(grid)    # for unit()
library(scales)  # pretty formatting
library(cowplot)
library(metafor) # required for REML random-effects fit
library(forcats)
library(stringr)

# Directory containing engineered / processed data files
processed_data_folder <- "data"

# Folder to store images
descriptive_figures_folder = fs::path("output", "figures", "descriptive")

# Paths to processed gene-level data and gene-set objects
p_load_expr_all_noNorm <- fs::path(processed_data_folder, "hipc_merged_all_noNorm.rds")

# Load data objects
hipc_merged_all_noNorm <- readRDS(p_load_expr_all_noNorm)

# Extract gene names
gene_names <- hipc_merged_all_noNorm %>%
  select(a1cf:zzz3) %>%
  select(where(~ !any(is.na(.)))) %>%
  colnames()

# Filter dataframe
# Inclusion criteria : participants must
# be healthy adults (18+ adults, no type 2 diabetes),
# receive inactivated influenza vaccine,
# have at least one immune response assay available at baseline and day 28 (+-7 days)
# have a baseline measurement and at least one non-baseline measurement
# observations occur on or before day 14 post-vaccination

hipc_merged_all_noNorm_filtered <- hipc_merged_all_noNorm %>%
  filter(
    !is.na(immResp_MFC_anyAssay_log2_MFC),
    # baseline and post-vaccine immune response for at least one assay
    vaccine_name == "Influenza (IN)",
    # receipt of inactivated influenza vaccine!(cohort %in% c("young T2D", "old T2D")),
    # remove T2 diabetes cohorts
    study_time_collected %in% c(0, 1, 2, 3, 7)
  ) %>%
  group_by(participant_id) %>%
  filter(sum(study_time_collected == 0) == 1,
         sum(study_time_collected > 0) > 0) %>% # baseline and at least one post-vaccination observation
  ungroup() %>%
  select(
    participant_id,
    study_accession,
    arm_accession,
    age_imputed,
    study_time_collected,
    immResp_MFC_nAb_pre_value,
    immResp_MFC_nAb_post_value,
    immResp_MFC_elisa_pre_value,
    immResp_MFC_elisa_post_value,
    immResp_MFC_hai_pre_value,
    immResp_MFC_hai_post_value,
    immResp_MFC_anyAssay_pre_value,
    immResp_MFC_anyAssay_post_value,
    gender,
    all_of(gene_names)
  ) %>%
  mutate(gender = droplevels(gender)) %>%
  arrange(participant_id)

study_names = hipc_merged_all_noNorm_filtered %>%
  pull(study_accession) %>%
  unique() %>%
  sort()

study_descriptions = c(
  "SDY1119" = "Systems Biology of 2011 trivalent Influenza vaccine (TIV) in young and elderly individuals, healthy or with T2D",
  "SDY1276" = "Time series of global gene expression after trivalent influenza vaccination in humans ",
  "SDY180" = "Systems scale interactive exploration reveals quantitative and qualitative differences in response to 2009-2010 Fluzone influenza vaccine and pneumococcal vaccine",
  "SDY224" = "Immune Responses to Seasonal TIV 2010-2011 Influenza Vaccination in Humans",
  "SDY269" = "Systems Biology of 2008 Influenza Vaccination in Humans",
  "SDY270" = "Systems Biology of 2009 Influenza Vaccination in Humans",
  "SDY400" = "Immunologic and genomic signatures of influenza vaccine response - 2012",
  "SDY404" = "Immunologic and genomic signatures of influenza vaccine response - 2011",
  "SDY520" = "Immunologic and genomic signatures of influenza vaccine response - 2013",
  "SDY56" = "Systems Biology of 2010 trivalent Influenza vaccine (TIV) in young and elderly",
  "SDY61" = "Systems Biology of 2007 Influenza Vaccination in Humans",
  "SDY63" = "Immunologic and genomic signatures of influenza vaccine response - 2010",
  "SDY640" = "Immunologic and genomic signatures of influenza vaccine response - 2014",
  "SDY67" = "Bioinformatics Approach to 2010-2011 TIV Influenza A/H1N1 Vaccine Immune Profiling",
  "SDY80" = "Cellular and molecular characterization of the immune response in healthy NIH employees at baseline, and after immunization with the H1N1 or seasonal influenza vaccines"
)

study_dates = c(
  "SDY1119" = "2011-01-01",
  "SDY1276" = "2013-06-17",
  "SDY180" = "2010-01-01",
  "SDY224" = "2011-01-01",
  "SDY269" = "2009-01-01",
  "SDY270" = "2009-01-01",
  "SDY400" = "2012-01-01",
  "SDY404" = "2011-01-01",
  "SDY520" = "2013-10-10",
  "SDY56" = "2010-01-01",
  "SDY61" = "2008-01-01",
  "SDY63" = "2011-01-01",
  "SDY640" = "2014-10-01",
  "SDY67" = "2010-01-01",
  "SDY80" = "2014-01-01"
)

study_order = c(
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

hipc_merged_all_noNorm_filtered = hipc_merged_all_noNorm_filtered %>% 
  mutate(study_accession = factor(study_accession, levels = study_order))

# Bubble plot - number of transcriptomic samples across time

# Counts per vaccine x time
counts <- hipc_merged_all_noNorm_filtered %>%
  filter(!is.na(study_time_collected)) %>%
  group_by(study_accession, study_time_collected) %>%
  summarise(n = n(), .groups = "drop")

# Order the time points and make study_time_collected an ordered factor
time_levels <- counts %>%
  distinct(study_time_collected) %>%
  arrange(as.numeric(study_time_collected)) %>%
  pull(study_time_collected) %>%
  as.character()   # factor levels must be character

# Order the counts by the study time
counts <- counts %>%
  mutate(study_time_collected = factor(
    as.character(study_time_collected),
    levels = time_levels,
    ordered = TRUE
  ))

# Since the range of the number of samples is large, we make the bubble size proportional to the count/square root of count
counts <- counts %>%
  mutate(size_var = n,
         study_accession =  factor(study_accession, levels = rev(levels(study_accession)))) 

# Bubble plot
p1 <- ggplot(counts, aes(x = study_time_collected, y = study_accession)) +
  # points coloured by vaccine hex codes; shape 21 allows fill + black border
  geom_point(
    aes(size = size_var),
    fill = "#5062FF" ,
    shape = 21,
    colour = "black",
    alpha = 0.8,
    show.legend = FALSE
  ) +
  # numeric labels inside bubbles; label shows raw counts n and is white
  geom_text(
    aes(label = n),
    colour = "white",
    size = 3.5,
    vjust = 0.5,
    show.legend = T
  ) +
  scale_fill_identity(guide = "none") +     # use hex codes directly, no fill legend
  scale_size_area(max_size = 28, guide = "none") +
  labs(
    x = "Days post-vaccination",
    y = "Study",
    title = "Number of transcriptomic observations across time per study",
    subtitle = "Selected timepoints, healthy adults with baseline and post-vaccination data"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 15),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 17, hjust = 0.5)
  )

print(p1)

ggsave(
  filename = "transcriptomic_samples_influenzain.pdf",
  path = descriptive_figures_folder,
  plot = p1,
  width = 30,
  height = 26,
  units = "cm"
)


# Immune Assays 

# --- 1. Reduce to one row per participant: the first study_time_collected observed
first_tp <- hipc_merged_all_noNorm_filtered %>%
  group_by(participant_id) %>%
  slice_min(study_time_collected, with_ties = FALSE) %>%
  ungroup()

# --- 2. Compute per-participant logical flags: TRUE only if both pre AND post are present
participant_assays <- first_tp %>%
  transmute(
    participant_id,
    study_accession,
    nAb_both  = !is.na(immResp_MFC_nAb_pre_value)   & !is.na(immResp_MFC_nAb_post_value),
    HAI_both  = !is.na(immResp_MFC_hai_pre_value)    & !is.na(immResp_MFC_hai_post_value)
  ) %>%
  pivot_longer(
    cols = c(nAb_both, HAI_both),
    names_to = "assay_raw",
    values_to = "has_both"
  ) %>%
  mutate(
    assay = recode(assay_raw,
                   nAb_both  = "nAb",
                   HAI_both  = "HAI")
  ) %>%
  select(participant_id, study_accession, assay, has_both)

# --- 3. Summarise: number of distinct participants per study x assay with both measurements
summary_df <- participant_assays %>%
  group_by(study_accession, assay) %>%
  summarise(
    n_participants = sum(has_both, na.rm = TRUE), # logical -> integer sum counts TRUEs
    .groups = "drop"
  )

# ensure all study x assay combinations are present (so tiles render even when 0)
all_combos <- expand_grid(
  study_accession = unique(first_tp$study_accession),
  assay = c("nAb", "HAI")
)

summary_df <- all_combos %>%
  left_join(summary_df, by = c("study_accession", "assay")) %>%
  mutate(
    n_participants = replace_na(n_participants, 0),
    has = n_participants > 0
  )

# --- 4. Reorder studies by total measured (descending) for nicer plotting (optional)
summary_df <- summary_df %>%
  mutate(
    study_accession = factor(study_accession, levels = study_order),
    assay = factor(assay, levels = c("nAb", "HAI"))
  )

# --- 5. Plot: geom_tile with counts annotated
# Increase base_size for global text enlargement; fine-tune specific element sizes below.
p2 <- ggplot(summary_df, aes(x = study_accession, y = assay)) +
  geom_tile(aes(fill = has), colour = "grey70") +
  geom_text(aes(label = n_participants), size = 5) + # larger, bold labels
  scale_fill_manual(values = c("TRUE" = "#2ca25f", "FALSE" = "white"), guide = FALSE) +
  labs(
    x = "Study accession",
    y = "Assay",
    title = "Availability of antibody assays per study",
    subtitle = "Participants with both baseline and day 28 (+- 7 days) post-vaccination measurements"
  ) +
  theme_minimal(base_size = 16) + # global text size increase
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),   # centered + bold title
    plot.subtitle = element_text(hjust = 0.5, size = 14),              # centered subtitle
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 14),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA)
  )

print(p2)

ggsave(
  filename = "antibody_samples_influenzain.pdf",
  path = descriptive_figures_folder,
  plot = p2,
  width = 30,
  height = 13,
  units = "cm"
)


# Proportion of each gender per study

hipc_merged_all_noNorm_filtered_unique = hipc_merged_all_noNorm_filtered %>%
  filter(study_time_collected == 0) %>%
  distinct(participant_id, .keep_all = TRUE)

# Summarise counts by study and gender
summary_df <- hipc_merged_all_noNorm_filtered %>%
  filter(study_time_collected == 0) %>%
  count(study_accession, gender, name = "n") %>%
  group_by(study_accession) %>%
  mutate(total = sum(n), prop = n / total) %>%
  ungroup()

# Plot: one 100% stacked bar per study, annotated with raw counts
p3 <- ggplot(summary_df, aes(x = study_accession, y = n, fill = gender)) +
  geom_col(position = "fill", width = 0.7) +
  # place raw counts centered inside each stacked segment
  geom_text(
    aes(label = n),
    position = position_fill(vjust = 0.5),
    size = 5,
    # label size (large)
    color = "white",
    # choose white for contrast; adjust if needed
    fontface = "bold"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Study",
       y = "Proportion",
       title = "Gender distribution by study",
       fill = "Gender") +
  theme_bw() +
  theme(
    text = element_text(size = 16),
    # make all text large
    plot.title = element_text(
      hjust = 0.5,
      # center the title
      face = "bold",
      size = 20
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 14
    ),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

# Print the plot
print(p3)

ggsave(
  filename = "sex_distribution_influenzain.pdf",
  path = descriptive_figures_folder,
  plot = p3,
  width = 30,
  height = 13,
  units = "cm"
)


# 2) Optionally remove participants with missing age_imputed
age_df <- hipc_merged_all_noNorm_filtered_unique %>%
  filter(!is.na(age_imputed)) %>%
  select(participant_id, study_accession, age_imputed)

age_df <- age_df %>%
  mutate(study_accession = factor(study_accession, levels = study_order))

# 4) Create the violin + jitter plot
p4 <- ggplot(age_df, aes(x = study_accession, y = age_imputed)) +
  geom_violin(trim = FALSE, width = 0.9, colour = "black", fill = "#5062FF") +
  geom_jitter(
    aes(),
    width = 0.15,
    height = 0,
    alpha = 0.3,
    size = 1.8
  ) +
  labs(
    x = "Study",
    y = "Age",
    title = "Age distribution by study"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# 5) Print and save
print(p4)

# adjust filename/path as needed; descriptive_figures_folder assumed defined earlier
ggsave(
  filename = "age_distribution_influenzain.pdf",
  plot = p4,
  path = descriptive_figures_folder,
  width = 30,
  height = 13,
  units = "cm"
)


# Create age groups using 64 as cutoff.
age_grouped <- hipc_merged_all_noNorm_filtered_unique %>%
  mutate(
    age_group = case_when(
      !is.na(age_imputed) & age_imputed <= 64 ~ "<=64",
      !is.na(age_imputed) & age_imputed > 64  ~ ">64",
      TRUE ~ NA_character_
    )
  )

# Summarise counts by study and age_group (exclude missing age)
summary_df_age <- age_grouped %>%
  filter(!is.na(age_group)) %>%
  count(study_accession, age_group, name = "n") %>%
  group_by(study_accession) %>%
  mutate(total = sum(n), prop = n / total) %>%
  ungroup()

summary_df_age <- summary_df_age %>%
  mutate(study_accession = factor(study_accession, levels = study_order),
         age_group = factor(age_group, levels = c("<=64", ">64")))

# Plot: one 100% stacked bar per study, annotated with raw counts
p5 <- ggplot(summary_df_age, aes(x = study_accession, y = n, fill = age_group)) +
  geom_col(position = "fill", width = 0.7) +
  # place raw counts centered inside each stacked segment
  geom_text(
    aes(label = n),
    position = position_fill(vjust = 0.5),
    size = 5,
    color = "white",
    fontface = "bold"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Study",
    y = "Proportion",
    title = "Age category distribution by study — cutoff at 64 years",
    fill = "Age group"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

# Print and save
print(p5)

ggsave(
  filename = "agecat_distribution_influenzain.pdf",
  path = descriptive_figures_folder,
  plot = p5,
  width = 30,
  height = 13,
  units = "cm"
)

