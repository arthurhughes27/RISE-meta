# Supplementary analyses: RISE meta-analysis across timepoints and marker types
library(tidyverse)
library(SurrogateRank)
library(knitr)
library(kableExtra)

set.seed(08012025)

# Paths to processed data and output figures
processed_data_folder <- "data"
application_figures_folder <- fs::path("output", "figures", "application")

# Load merged gene expression and BTM gene set objects
hipc_merged_all_noNorm <- readRDS(fs::path(processed_data_folder, "hipc_merged_all_noNorm.rds"))
BTM <- readRDS(fs::path(processed_data_folder, "BTM_processed.rds"))

# Gene columns present in the data with no missing values
gene_names <- hipc_merged_all_noNorm %>%
  select(a1cf:zzz3) %>%
  select(where(~ !any(is.na(.)))) %>%
  colnames()

# nAb studies use a different response variable than HAI studies
nab_studies <- c("SDY80", "SDY180", "SDY1276", "SDY67")

# BTM gene sets and their names (excluding top-level aggregates)
btm_filter <- which(BTM[["geneset.aggregates"]] != "NA")
BTM_genes <- BTM[["genesets"]][btm_filter] %>% unlist()
BTM_genes_names <- BTM[["geneset.names.descriptions"]][btm_filter]

# Summarise gene expression within each BTM gene set by applying FUN row-wise
aggregate_to_geneset <- function(df, genesets, geneset_names, FUN = mean) {
  df <- as.data.frame(df)
  out <- imap_dfc(genesets, function(genes, i) {
    present <- intersect(genes, colnames(df))
    if (length(present) == 0)
      return(NULL)
    mat <- df[, present, drop = FALSE]
    vec <- apply(mat, 1, function(r) FUN(r, na.rm = TRUE))
    # Replace NaN (all-NA rows) with proper NA
    vec[is.nan(vec)] <- NA_real_
    tibble::tibble(!!geneset_names[i] := vec)
  })
  as.data.frame(out)
}

# Extract RISE inputs (response vectors, gene-set matrices, study labels)
# from a data frame containing both baseline and post-vaccination rows;
# marker_cols controls whether BTM gene sets or individual genes are used
extract_rise_inputs <- function(df, marker_cols) {
  sone_raw  <- df %>% filter(study_time_collected > 0)  %>% select(any_of(marker_cols))
  szero_raw <- df %>% filter(study_time_collected == 0) %>% select(any_of(marker_cols))

  list(
    yone      = df %>% filter(study_time_collected > 0)  %>% pull(response_post),
    yzero     = df %>% filter(study_time_collected == 0) %>% pull(response_pre),
    sone      = aggregate_to_geneset(sone_raw,  BTM[["genesets"]][btm_filter], BTM_genes_names),
    szero     = aggregate_to_geneset(szero_raw, BTM[["genesets"]][btm_filter], BTM_genes_names),
    studyone  = df %>% filter(study_time_collected > 0)  %>% pull(study_accession),
    studyzero = df %>% filter(study_time_collected == 0) %>% pull(study_accession)
  )
}

# Filter the merged data to influenza participants at the specified timepoints,
# requiring exactly one baseline and one post-vaccination row per participant,
# and dropping studies with fewer than 3 participants
filter_hipc <- function(tp) {
  timepoints_to_keep <- c(0, tp)
  hipc_merged_all_noNorm %>%
    mutate(
      response_pre = ifelse(
        study_accession %in% nab_studies,
        immResp_mean_nAb_pre_value,
        immResp_mean_hai_pre_value
      ),
      response_post = ifelse(
        study_accession %in% nab_studies,
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
    filter(
      sum(study_time_collected == 0) == 1,
      sum(study_time_collected == tp) == 1
    ) %>%
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
}

# Build a forest plot of pooled effect estimates for the top N markers,
# coloured by significance category, with dashed epsilon equivalence bounds
make_screening_forest_plot <- function(screening_metrics, top_N = 15, epsilon_val = 0.2) {
  df_plot <- screening_metrics %>%
    arrange(p.unadjusted) %>%
    slice_head(n = top_N) %>%
    mutate(
      # Reverse factor order so top marker appears at the top of the y-axis
      marker = factor(marker, levels = rev(unique(marker))),
      sig_cat = case_when(
        p.adjusted < 0.05                         ~ "adj < 0.05",
        p.adjusted >= 0.05 & p.unadjusted < 0.05 ~ "unadj < 0.05 only",
        TRUE                                      ~ "ns"
      ),
      sig_cat = factor(sig_cat, levels = c("adj < 0.05", "unadj < 0.05 only", "ns"))
    )

  ggplot(df_plot, aes(x = mu.delta, y = marker)) +
    geom_segment(
      aes(x = ci.delta.lower, xend = ci.delta.upper, y = marker, yend = marker, color = sig_cat),
      size = 1.1, lineend = "round"
    ) +
    geom_point(aes(color = sig_cat), size = 4) +
    # Dashed lines mark the equivalence bounds (+/- epsilon)
    geom_vline(xintercept = c(-epsilon_val, epsilon_val), linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.25)) +
    scale_color_manual(
      values = c(
        "adj < 0.05"        = "#F20A0A",
        "unadj < 0.05 only" = "#0A63F2",
        "ns"                = "black"
      ),
      labels = c(
        "adj < 0.05"        = "BH-adjusted p < 0.05",
        "unadj < 0.05 only" = "Unadjusted p < 0.05 only",
        "ns"                = "Not significant"
      ),
      name = "Significance"
    ) +
    labs(
      x     = expression("Pooled effect " ~ mu[delta]),
      y     = NULL,
      title = glue::glue("Screening results: Top {top_N} markers by p-value")
    ) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title         = element_text(size = 25, face = "bold", hjust = 0.5),
      axis.text.y        = element_text(size = 13),
      axis.text.x        = element_text(size = 15),
      axis.title.x       = element_text(size = 30),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

# Format significant screening markers as a LaTeX table
make_screening_table <- function(screening_metrics) {
  screening_metrics_select <- screening_metrics %>%
    arrange(p.unadjusted) %>%
    mutate(mu.delta.ci = paste0(
      round(mu.delta, 3), " (",
      round(ci.delta.lower, 3), ", ",
      round(ci.delta.upper, 3), ")"
    )) %>%
    select(marker, mu.delta.ci, p.unadjusted, p.adjusted) %>%
    filter(p.adjusted < 0.05)

  kable(
    screening_metrics_select,
    format   = "latex",
    booktabs = TRUE,
    caption  = "Screening Metrics Table"
  ) %>%
    kable_styling(latex_options = "hold_position") %>%
    row_spec(0, bold = TRUE) %>%
    column_spec(1:ncol(screening_metrics_select))
}

# Run rise.screen.meta with standard settings and BH correction
run_screen <- function(inputs, p_correction = "BH") {
  rise.screen.meta(
    inputs$yone,
    inputs$yzero,
    sone                         = inputs$sone,
    szero                        = inputs$szero,
    inputs$studyone,
    inputs$studyzero,
    alpha                        = 0.05,
    epsilon.meta                 = 0.2,
    alternative                  = "two.sided",
    paired.all                   = TRUE,
    return.all.screen            = TRUE,
    epsilon.study                = 0.2,
    p.correction                 = p_correction,
    show.pooled.effect           = TRUE,
    return.study.similarity.plot = FALSE
  )
}

# ===== BTM analysis at day 2 (screening only, full data) =====

hipc_d2 <- filter_hipc(tp = 2)

# No train/test split for this supplementary screening: use all participants
inputs_d2 <- extract_rise_inputs(hipc_d2, BTM_genes)

screen_d2 <- run_screen(inputs_d2)

screening_metrics_d2 <- screen_d2[["screening.metrics.meta"]]

# Render LaTeX table of significant markers
make_screening_table(screening_metrics_d2)

# Forest and fit plots from the RISE screening output
p1 <- screen_d2[["gamma.s.plot"]]$forest.plot
p1

p2 <- screen_d2[["gamma.s.plot"]]$fit.plot
p2

# Custom screening forest plot for the top 15 markers
p3 <- make_screening_forest_plot(screening_metrics_d2, top_N = 15)
p3

ggsave(
  filename = "TIV_d2_screening.pdf",
  path     = application_figures_folder,
  plot     = p3,
  width    = 40,
  height   = 18,
  units    = "cm"
)

# ===== BTM analysis at day 1+2 (with train/test split and evaluation) =====

# Use day 1 as primary post timepoint; fall back to day 2 if day 1 is absent
tp         <- 1L
fallback_tp <- 2L
timepoints_to_keep <- c(0L, tp, fallback_tp)

# Build unified response columns and select the best available post timepoint
# per participant (day 1 preferred over day 2)
hipc_d1d2 <- hipc_merged_all_noNorm %>%
  mutate(
    response_pre = ifelse(
      study_accession %in% nab_studies,
      immResp_mean_nAb_pre_value,
      immResp_mean_hai_pre_value
    ),
    response_post = ifelse(
      study_accession %in% nab_studies,
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
  mutate(
    has_pre = sum(study_time_collected == 0L),
    has_tp1 = sum(study_time_collected == tp),
    has_tp2 = sum(study_time_collected == fallback_tp),
    # Prefer day 1; fall back to day 2 if day 1 is unavailable
    desired_post = case_when(
      has_tp1 > 0 ~ tp,
      has_tp2 > 0 ~ fallback_tp,
      TRUE        ~ NA_integer_
    )
  ) %>%
  # Keep only participants with exactly one baseline and one chosen post row
  filter(
    has_pre == 1,
    !is.na(desired_post),
    study_time_collected %in% c(0L, desired_post)
  ) %>%
  filter(
    sum(study_time_collected == 0L) == 1,
    sum(study_time_collected == desired_post) == 1
  ) %>%
  ungroup() %>%
  group_by(study_accession) %>%
  filter(n_distinct(participant_id) > 2) %>%
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

# Sample 66% of participants per study for training; remainder becomes test set
train_indices_d1d2 <- hipc_d1d2 %>%
  distinct(study_accession, participant_id) %>%
  group_by(study_accession) %>%
  slice_sample(prop = 0.66) %>%
  ungroup()

df_train_d1d2 <- hipc_d1d2 %>%
  semi_join(train_indices_d1d2, by = c("study_accession", "participant_id"))

df_test_d1d2 <- hipc_d1d2 %>%
  anti_join(train_indices_d1d2, by = c("study_accession", "participant_id"))

# Screen on training data
train_inputs_d1d2 <- extract_rise_inputs(df_train_d1d2, BTM_genes)
screen_d1d2 <- run_screen(train_inputs_d1d2)

screening_metrics_d1d2 <- screen_d1d2[["screening.metrics.meta"]]

# Render LaTeX table of significant markers
make_screening_table(screening_metrics_d1d2)

# Forest and fit plots from the RISE screening output
p1 <- screen_d1d2[["gamma.s.plot"]]$forest.plot
p1

p2 <- screen_d1d2[["gamma.s.plot"]]$fit.plot
p2

# Custom screening forest plot for the top 15 markers
p4 <- make_screening_forest_plot(screening_metrics_d1d2, top_N = 15)
p4

ggsave(
  filename = "TIV_d1d2_screening.pdf",
  path     = application_figures_folder,
  plot     = p4,
  width    = 40,
  height   = 18,
  units    = "cm"
)

# Evaluate significant markers from screening on held-out test data
test_inputs_d1d2 <- extract_rise_inputs(df_test_d1d2, BTM_genes)

evaluate_d1d2 <- rise.evaluate.meta(
  test_inputs_d1d2$yone,
  test_inputs_d1d2$yzero,
  sone               = test_inputs_d1d2$sone,
  szero              = test_inputs_d1d2$szero,
  test_inputs_d1d2$studyone,
  test_inputs_d1d2$studyzero,
  alpha              = 0.05,
  epsilon.meta       = 0.2,
  alternative        = "two.sided",
  paired.all         = TRUE,
  epsilon.study      = 0.2,
  p.correction       = "none",
  show.pooled.effect = TRUE,
  # Pass screening weights and selected markers from the training stage
  screening.weights  = screen_d1d2[["screening.weights"]],
  markers            = screen_d1d2[["significant.markers"]]
)

# Save evaluation forest plot
p4 <- evaluate_d1d2[["gamma.s.plot"]]$forest.plot
p4

ggsave(
  filename = "TIV_evaluation_forest_d1d2.pdf",
  path     = application_figures_folder,
  plot     = p4,
  width    = 32,
  height   = 15,
  units    = "cm"
)

# Save evaluation rank-correlation fit plot
p5 <- evaluate_d1d2[["gamma.s.plot"]]$fit.plot
p5

ggsave(
  filename = "TIV_evaluation_fitplot_d1d2.pdf",
  path     = application_figures_folder,
  plot     = p5,
  width    = 37,
  height   = 20,
  units    = "cm"
)

# ===== BTM analysis at day 7 (screening only, full data) =====

hipc_d7 <- filter_hipc(tp = 7)

# No train/test split: use all participants for this supplementary screening
inputs_d7 <- extract_rise_inputs(hipc_d7, BTM_genes)

screen_d7 <- run_screen(inputs_d7)

screening_metrics_d7 <- screen_d7[["screening.metrics.meta"]]

# Render LaTeX table of significant markers
make_screening_table(screening_metrics_d7)

# Forest and fit plots from the RISE screening output
p1 <- screen_d7[["gamma.s.plot"]]$forest.plot
p1

p2 <- screen_d7[["gamma.s.plot"]]$fit.plot
p2

# Custom screening forest plot for the top 15 markers
p3 <- make_screening_forest_plot(screening_metrics_d7, top_N = 15)
p3

ggsave(
  filename = "TIV_d7_screening.pdf",
  path     = application_figures_folder,
  plot     = p3,
  width    = 40,
  height   = 18,
  units    = "cm"
)

# Re-run screening without p-value correction for the forest plot figure
screen_d7_uncorrected <- run_screen(inputs_d7, p_correction = "none")

# Save the uncorrected forest plot
p1 <- screen_d7_uncorrected[["gamma.s.plot"]]$forest.plot
p1

ggsave(
  filename = "TIV_evaluation_forest_d7.pdf",
  path     = application_figures_folder,
  plot     = p1,
  width    = 32,
  height   = 15,
  units    = "cm"
)

# ===== Gene-wise screening (individual genes as markers, full data) =====
# The following four blocks repeat the same analysis at days 1, 2, 3, and 7,
# using individual genes rather than BTM gene set summaries as markers.

# --- Day 1 ---

hipc_d1_genes <- filter_hipc(tp = 1)

# Use individual gene columns directly as surrogate markers
inputs_d1_genes <- list(
  yone      = hipc_d1_genes %>% filter(study_time_collected > 0)  %>% pull(response_post),
  yzero     = hipc_d1_genes %>% filter(study_time_collected == 0) %>% pull(response_pre),
  sone      = hipc_d1_genes %>% filter(study_time_collected > 0)  %>% select(any_of(gene_names)),
  szero     = hipc_d1_genes %>% filter(study_time_collected == 0) %>% select(any_of(gene_names)),
  studyone  = hipc_d1_genes %>% filter(study_time_collected > 0)  %>% pull(study_accession),
  studyzero = hipc_d1_genes %>% filter(study_time_collected == 0) %>% pull(study_accession)
)

screen_d1_genes <- run_screen(inputs_d1_genes)

screening_metrics_d1_genes <- screen_d1_genes[["screening.metrics.meta"]]

# Render LaTeX table of significant markers
make_screening_table(screening_metrics_d1_genes)

# Forest and fit plots from the RISE screening output
p1 <- screen_d1_genes[["gamma.s.plot"]]$forest.plot
p1

p2 <- screen_d1_genes[["gamma.s.plot"]]$fit.plot
p2

# Custom screening forest plot for the top 20 genes
p3 <- make_screening_forest_plot(screening_metrics_d1_genes, top_N = 20)
p3

ggsave(
  filename = "TIV_d1_screening_genewise.pdf",
  path     = application_figures_folder,
  plot     = p3,
  width    = 40,
  height   = 18,
  units    = "cm"
)

# --- Day 2 ---

hipc_d2_genes <- filter_hipc(tp = 2)

inputs_d2_genes <- list(
  yone      = hipc_d2_genes %>% filter(study_time_collected > 0)  %>% pull(response_post),
  yzero     = hipc_d2_genes %>% filter(study_time_collected == 0) %>% pull(response_pre),
  sone      = hipc_d2_genes %>% filter(study_time_collected > 0)  %>% select(any_of(gene_names)),
  szero     = hipc_d2_genes %>% filter(study_time_collected == 0) %>% select(any_of(gene_names)),
  studyone  = hipc_d2_genes %>% filter(study_time_collected > 0)  %>% pull(study_accession),
  studyzero = hipc_d2_genes %>% filter(study_time_collected == 0) %>% pull(study_accession)
)

screen_d2_genes <- run_screen(inputs_d2_genes)

screening_metrics_d2_genes <- screen_d2_genes[["screening.metrics.meta"]]

# Render LaTeX table of significant markers
make_screening_table(screening_metrics_d2_genes)

# Forest and fit plots from the RISE screening output
p1 <- screen_d2_genes[["gamma.s.plot"]]$forest.plot
p1

p2 <- screen_d2_genes[["gamma.s.plot"]]$fit.plot
p2

# Custom screening forest plot for the top 20 genes
p3 <- make_screening_forest_plot(screening_metrics_d2_genes, top_N = 20)
p3

ggsave(
  filename = "TIV_d2_screening_genewise.pdf",
  path     = application_figures_folder,
  plot     = p3,
  width    = 40,
  height   = 18,
  units    = "cm"
)

# --- Day 3 ---

hipc_d3_genes <- filter_hipc(tp = 3)

inputs_d3_genes <- list(
  yone      = hipc_d3_genes %>% filter(study_time_collected > 0)  %>% pull(response_post),
  yzero     = hipc_d3_genes %>% filter(study_time_collected == 0) %>% pull(response_pre),
  sone      = hipc_d3_genes %>% filter(study_time_collected > 0)  %>% select(any_of(gene_names)),
  szero     = hipc_d3_genes %>% filter(study_time_collected == 0) %>% select(any_of(gene_names)),
  studyone  = hipc_d3_genes %>% filter(study_time_collected > 0)  %>% pull(study_accession),
  studyzero = hipc_d3_genes %>% filter(study_time_collected == 0) %>% pull(study_accession)
)

screen_d3_genes <- run_screen(inputs_d3_genes)

screening_metrics_d3_genes <- screen_d3_genes[["screening.metrics.meta"]]

# Render LaTeX table of significant markers
make_screening_table(screening_metrics_d3_genes)

# Forest and fit plots from the RISE screening output
p1 <- screen_d3_genes[["gamma.s.plot"]]$forest.plot
p1

p2 <- screen_d3_genes[["gamma.s.plot"]]$fit.plot
p2

# Custom screening forest plot for the top 20 genes
p3 <- make_screening_forest_plot(screening_metrics_d3_genes, top_N = 20)
p3

ggsave(
  filename = "TIV_d3_screening_genewise.pdf",
  path     = application_figures_folder,
  plot     = p3,
  width    = 40,
  height   = 18,
  units    = "cm"
)

# --- Day 7 ---

hipc_d7_genes <- filter_hipc(tp = 7)

inputs_d7_genes <- list(
  yone      = hipc_d7_genes %>% filter(study_time_collected > 0)  %>% pull(response_post),
  yzero     = hipc_d7_genes %>% filter(study_time_collected == 0) %>% pull(response_pre),
  sone      = hipc_d7_genes %>% filter(study_time_collected > 0)  %>% select(any_of(gene_names)),
  szero     = hipc_d7_genes %>% filter(study_time_collected == 0) %>% select(any_of(gene_names)),
  studyone  = hipc_d7_genes %>% filter(study_time_collected > 0)  %>% pull(study_accession),
  studyzero = hipc_d7_genes %>% filter(study_time_collected == 0) %>% pull(study_accession)
)

screen_d7_genes <- run_screen(inputs_d7_genes)

screening_metrics_d7_genes <- screen_d7_genes[["screening.metrics.meta"]]

# Render LaTeX table of significant markers
make_screening_table(screening_metrics_d7_genes)

# Forest and fit plots from the RISE screening output
p1 <- screen_d7_genes[["gamma.s.plot"]]$forest.plot
p1

p2 <- screen_d7_genes[["gamma.s.plot"]]$fit.plot
p2

# Custom screening forest plot for the top 20 genes
p3 <- make_screening_forest_plot(screening_metrics_d7_genes, top_N = 20)
p3

ggsave(
  filename = "TIV_d7_screening_genewise.pdf",
  path     = application_figures_folder,
  plot     = p3,
  width    = 40,
  height   = 18,
  units    = "cm"
)
