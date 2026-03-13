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
  # mutate(study_time_collected = ifelse(study_time_collected %in% c(2), 1, study_time_collected)) %>%
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

# For each study, randomly assign 66% of participants to training
train_indices <- hipc_merged_all_noNorm_filtered %>%
  distinct(study_accession, participant_id) %>%  # unique participants per study
  group_by(study_accession) %>%
  slice_sample(prop = 0.66) %>%                # sample 66% of participants per study
  ungroup()

# df_train: all rows corresponding to selected participants
df_train <- hipc_merged_all_noNorm_filtered %>%
  semi_join(train_indices, by = c("study_accession", "participant_id"))

# df_test: the remaining participants
df_test <- hipc_merged_all_noNorm_filtered %>%
  anti_join(train_indices, by = c("study_accession", "participant_id"))


innate_pathways = c(
  "Interferon/Antiviral Sensing",
  "Monocytes",
  "Antigen Presentation",
  "NK Cells",
  "Neutrophils",
  "Inflammatory/TLR/Chemokines",
  "DC Activation",
  "T Cells",
  "B Cells",
  "Plasma Cells",
  "ECM And Migration",
  "Platelets",
  "Cell Cycle", 
  "Signal Transduction",
  "Energy Metabolism"
)

BTM_genes = BTM[["genesets"]][which(BTM[["geneset.aggregates"]] %in% innate_pathways)] %>%
  unlist()

BTM_genes_names = BTM[["geneset.names.descriptions"]][which(BTM[["geneset.aggregates"]] %in% innate_pathways)]

yone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(response_post)

yzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(response_pre)

sone = df_train %>%
  filter(study_time_collected > 0) %>%
  select(any_of(BTM_genes))

szero = df_train %>%
  filter(study_time_collected == 0) %>%
  select(any_of(BTM_genes))

aggregate_to_geneset <- function(df, genesets, geneset_names, FUN = mean) {
  df <- as.data.frame(df)
  out <- imap_dfc(genesets, function(genes, i) {
    present <- intersect(genes, colnames(df))
    if (length(present) == 0)
      return(NULL)
    colname <- geneset_names[i]
    mat <- df[, present, drop = FALSE]
    vec <- apply(mat, 1, function(r)
      FUN(r, na.rm = TRUE))
    vec[is.nan(vec)] <- NA_real_
    tibble::tibble(!!colname := vec)
  })
  as.data.frame(out)
}

# aggregate sone and szero
sone <- aggregate_to_geneset(
  df = sone,
  genesets = BTM[["genesets"]][which(BTM[["geneset.aggregates"]] %in% innate_pathways)],
  geneset_names = BTM_genes_names,
  FUN = mean
)

szero <- aggregate_to_geneset(
  df = szero,
  genesets = BTM[["genesets"]][which(BTM[["geneset.aggregates"]] %in% innate_pathways)],
  geneset_names = BTM_genes_names,
  FUN = mean
)

studyone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession)

studyzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(study_accession)

rise.screen.meta.result = rise.screen.meta(
  yone,
  yzero,
  sone = sone,
  szero = szero,
  studyone,
  studyzero,
  alpha = 0.05,
  epsilon.meta = 0.2,
  alternative = "two.sided",
  paired.all = T,
  return.all.screen = T,
  epsilon.study = 0.2,
  p.correction = "BH",
  show.pooled.effect = T,
  return.study.similarity.plot = F
)

screening.metrics.meta = rise.screen.meta.result[["screening.metrics.meta"]]

screening.metrics.meta.select = screening.metrics.meta %>%
  arrange(p.unadjusted) %>%
  mutate(mu.delta.ci = paste0(
    round(mu.delta, 3),
    " (",
    round(ci.delta.lower, 3),
    ", ",
    round(ci.delta.upper, 3) ,
    ")"
  )) %>%
  select(marker, mu.delta.ci, p.unadjusted, p.adjusted)  %>%
  filter(p.adjusted < 0.05)

# Suppose your dataframe is screening.metrics.meta.select
kable(
  screening.metrics.meta.select,
  format = "latex",
  booktabs = TRUE,
  caption = "Screening Metrics Table"
) %>%
  kable_styling(latex_options = "hold_position") %>%
  row_spec(0, bold = TRUE) %>%              # make header bold
  column_spec(1:ncol(screening.metrics.meta.select))

p1 = rise.screen.meta.result[["gamma.s.plot"]]$forest.plot

p1

p2 = rise.screen.meta.result[["gamma.s.plot"]]$fit.plot

p2

yone = df_test %>%
  filter(study_time_collected > 0) %>%
  pull(response_post)

yzero = df_test %>%
  filter(study_time_collected == 0) %>%
  pull(response_pre)

sone = df_test %>%
  filter(study_time_collected > 0) %>%
  select(any_of(BTM_genes))

szero = df_test %>%
  filter(study_time_collected == 0) %>%
  select(any_of(BTM_genes))

# aggregate sone and szero
sone <- aggregate_to_geneset(
  df = sone,
  genesets = BTM[["genesets"]][which(BTM[["geneset.aggregates"]] %in% innate_pathways)],
  geneset_names = BTM_genes_names,
  FUN = mean
)

szero <- aggregate_to_geneset(
  df = szero,
  genesets = BTM[["genesets"]][which(BTM[["geneset.aggregates"]] %in% innate_pathways)],
  geneset_names = BTM_genes_names,
  FUN = mean
)

studyone = df_test %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession)

studyzero = df_test %>%
  filter(study_time_collected == 0) %>%
  pull(study_accession)

rise.evaluate.meta.result = rise.evaluate.meta(
  yone,
  yzero,
  sone = sone,
  szero = szero,
  studyone,
  studyzero,
  alpha = 0.05,
  epsilon.meta = 0.2,
  alternative = "two.sided",
  paired.all = T,
  epsilon.study = 0.2,
  p.correction = "none",
  show.pooled.effect = T,
  screening.weights = rise.screen.meta.result[["screening.weights"]],
  markers = rise.screen.meta.result[["significant.markers"]]
)

p3 = rise.evaluate.meta.result[["gamma.s.plot"]]$forest.plot

p3

p4 = rise.evaluate.meta.result[["gamma.s.plot"]]$fit.plot

p4
