# Simulation study with permuted dataset to examine the calibration of the approach
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

source("./R/generate_permuted_dataset.R")

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

alpha = 0.05

tp <- c(7)
timepoints_to_keep <- c(0, tp)

# Within-study sample size grid: smallest first, full dataset last (NA = no resampling)
sample_sizes <- c(10, 25, 35, NA_real_)
sample_sizes <- c(NA_real_)
size_labels  <- ifelse(is.na(sample_sizes), "Full", as.character(sample_sizes))


df <- hipc_merged_all_noNorm %>%
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
  # filter(length(unique(participant_id)) > max(sample_sizes, na.rm = TRUE)) %>%
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

n_studies_max = df$study_accession %>%
  unique() %>%
  length()

# Number-of-studies grid: 2 to 5, plus all studies (NA = no subsampling)
n_studies_grid   <- c(3, 4, 5, 8, n_studies_max)
n_studies_labels <- ifelse(is.na(n_studies_grid), "All", as.character(n_studies_grid))

# Number of permutation replicates
B <- 100

# Pre-generate seeds for each replicate from the master RNG so results are
# fully reproducible regardless of iteration order or future code changes.
perm_seeds <- sample.int(n = .Machine$integer.max, size = B)

# Collect results across the full grid into a flat list of data frames
results <- vector("list", length(n_studies_grid) * length(sample_sizes))
idx <- 1L

for (j in seq_along(n_studies_grid)) {
  n_st     <- n_studies_grid[j]
  n_st_lbl <- n_studies_labels[j]
  
  for (i in seq_along(sample_sizes)) {
    n_study <- sample_sizes[i]
    lbl     <- size_labels[i]
    
    cat(sprintf(
      "\n==== n_studies: %s | n_within_study: %s ====\n",
      n_st_lbl,
      lbl
    ))
    
    fpr_vec <- numeric(B)
    
    for (b in seq_len(B)) {
      df_perm_list <- generate_permuted_dataset(
        df             = df,
        gene_names     = gene_names,
        seed           = perm_seeds[b],
        n_within_study = if (is.na(n_study))
          NULL
        else
          n_study,
        n_studies      = if (is.na(n_st))
          NULL
        else
          n_st
      )
      
      yone      <- df_perm_list$yone
      yzero     <- df_perm_list$yzero
      sone      <- df_perm_list$sone
      szero     <- df_perm_list$szero
      studyone  <- df_perm_list$studyone
      studyzero <- df_perm_list$studyzero
      
      rise.screen.meta.result <- rise.screen.meta(
        yone,
        yzero,
        sone = sone,
        szero = szero,
        studyone,
        studyzero,
        alpha = alpha,
        epsilon.meta = 0.2,
        alternative = "two.sided",
        paired.all = TRUE,
        return.all.screen = TRUE,
        epsilon.study = 0.2,
        p.correction = "BH",
        show.pooled.effect = TRUE,
        return.study.similarity.plot = FALSE,
        n.cores = 6,
        return.fit.plot = FALSE,
        return.forest.plot = FALSE,
        normalise.weights = FALSE,
        weight.mode = "diff.epsilon",
        return.all.weights = FALSE,
        return.evaluate.results = FALSE
      )
      
      p_vals <- rise.screen.meta.result[["screening.metrics.meta"]]$p.unadjusted
      fpr_vec[b] <- mean(p_vals < alpha)
      
      cat(sprintf("Permutation %d / %d  |  FPR = %.4f\n", b, B, fpr_vec[b]))
    }
    
    cat(
      sprintf(
        "Mean FPR (n_studies=%s, n_within_study=%s): %.4f\n",
        n_st_lbl,
        lbl,
        mean(fpr_vec)
      )
    )
    
    results[[idx]] <- data.frame(fpr = fpr_vec,
                                 n_within_study = lbl,
                                 n_studies = n_st_lbl)
    idx <- idx + 1L
  }
}

# Convert to long data frame with ordered factors
fpr_df <- bind_rows(results)

results_folder = fs::path("output", "simulation")

# Use fs::path() to specify the data path robustly
p_save <- fs::path(results_folder, "simulation_nonparametric_nstudies.rds")

# Save dataframe
saveRDS(fpr_df, file = p_save)

# Faceted boxplot: x = within-study sample size, facets = number of studies
ggplot(fpr_df, aes(x = n_within_study, y = fpr)) +
  geom_boxplot() +
  geom_hline(
    yintercept = alpha,
    colour = "black",
    linetype = "dashed",
    alpha = 0.8
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap( ~ n_studies, labeller = label_both) +
  labs(
    x = "Within-study sample size",
    y = "False positive rate",
    title = sprintf("FPR distribution over %d permutations", B)
  )
