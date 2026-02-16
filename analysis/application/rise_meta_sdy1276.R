# R script to illustrate a meta-analysis with RISE
library(tidyverse)
library(SurrogateRank)
library(grid)    # for unit()
library(scales)  # pretty formatting
library(cowplot)
library(metafor) # required for REML random-effects fit
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

BTM_interferon = which(BTM[["geneset.aggregates"]] == "Interferon/Antiviral Sensing")
BTM_filtered = BTM[["genesets"]][BTM_interferon] %>%
  unlist()

gene_names <- hipc_merged_all_noNorm %>%
  select(a1cf:zzz3) %>%
  select(where(~ !any(is.na(.)))) %>%
  colnames()

tp <- 1
timepoints_to_keep <- c(0, tp)

hipc_merged_all_noNorm_filtered <- hipc_merged_all_noNorm %>%
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

# Here, we do various analyses of the 222-gene signature identified in the RISE paper

# Load the gene list

path_weights = fs::path(
  "/home",
  "ah3",
  "Desktop",
  "Work",
  "PhD",
  "RISE-Project",
  "output",
  "application",
  "weights_influenzain_sdy1276_female.rds"
)

screening.weights = readRDS(path_weights) %>%
  mutate(marker = tolower(marker))

markers = screening.weights$marker

# Analysis 1 : SDY1276 males
hipc_merged_all_noNorm_filtered_copy = hipc_merged_all_noNorm_filtered %>%
  filter(study_accession == "SDY1276") %>%
  mutate(study_accession = ifelse(gender == "Female", "SDY1276_Female", "SDY1276_Male"))

df_train = hipc_merged_all_noNorm_filtered_copy %>%
  group_by(study_accession) %>%
  # filter(length(unique(participant_id)) > 8) %>%
  ungroup()

yone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(response_post)

yzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(response_pre)

sone = df_train %>%
  filter(study_time_collected > 0) %>%
  select(all_of(markers))

szero = df_train %>%
  filter(study_time_collected == 0) %>%
  select(all_of(markers))

studyone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession)

studyzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(study_accession)

rise.evaluate.meta.result = rise.evaluate.meta(
  yone,
  yzero,
  sone,
  szero,
  studyone,
  studyzero,
  alpha = 0.05,
  epsilon.meta = 0.2,
  alternative = "two.sided",
  paired.all = T,
  evaluate.weights = T,
  return.all.evaluate = T,
  epsilon.study = 0.2,
  p.correction = "BH",
  screening.weights = screening.weights,
  markers = markers,
  # weight.mode = "inverse.delta",
  show.pooled.effect = F
)

p1 = rise.evaluate.meta.result[["gamma.s.plot"]]$forest.plot

p1

ggsave(
  filename = "rise_influenzain_sdy1276_gender.pdf",
  path = application_figures_folder,
  plot = p1,
  width = 27,
  height = 13,
  units = "cm"
)

# Analysis two : across all studies, all participants
hipc_merged_all_noNorm_filtered_copy = hipc_merged_all_noNorm_filtered

df_train = hipc_merged_all_noNorm_filtered_copy %>%
  group_by(study_accession) %>%
  ungroup()

yone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(response_post)

yzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(response_pre)

sone = df_train %>%
  filter(study_time_collected > 0) %>%
  select(all_of(markers))

szero = df_train %>%
  filter(study_time_collected == 0) %>%
  select(all_of(markers))

studyone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession)

studyzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(study_accession)

rise.evaluate.meta.result = rise.evaluate.meta(
  yone,
  yzero,
  sone,
  szero,
  studyone,
  studyzero,
  alpha = 0.05,
  epsilon.meta = 0.2,
  alternative = "two.sided",
  paired.all = T,
  evaluate.weights = T,
  return.all.evaluate = T,
  epsilon.study = 0.2,
  p.correction = "BH",
  screening.weights = screening.weights,
  markers = markers,
  show.pooled.effect = F
)

p2 = rise.evaluate.meta.result[["gamma.s.plot"]]$forest.plot

p2

ggsave(
  filename = "rise_influenzain_crossstudy_all.pdf",
  path = application_figures_folder,
  plot = p2,
  width = 27,
  height = 18,
  units = "cm"
)

# Analysis three : across all studies, stratified by gender
hipc_merged_all_noNorm_filtered_copy = hipc_merged_all_noNorm_filtered %>%
  mutate(study_accession = ifelse(
    gender == "Female",
    paste0(study_accession, "_Female"),
    paste0(study_accession, "_Male")
  ))

df_train = hipc_merged_all_noNorm_filtered_copy %>%
  group_by(study_accession) %>%
  ungroup()

yone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(response_post)

yzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(response_pre)

sone = df_train %>%
  filter(study_time_collected > 0) %>%
  select(all_of(markers))

szero = df_train %>%
  filter(study_time_collected == 0) %>%
  select(all_of(markers))

studyone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession)

studyzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(study_accession)

rise.evaluate.meta.result = rise.evaluate.meta(
  yone,
  yzero,
  sone,
  szero,
  studyone,
  studyzero,
  alpha = 0.05,
  epsilon.meta = 0.2,
  alternative = "two.sided",
  paired.all = T,
  evaluate.weights = T,
  return.all.evaluate = T,
  epsilon.study = 0.2,
  p.correction = "BH",
  screening.weights = screening.weights,
  markers = markers,
  show.pooled.effect = F
)

p3 = rise.evaluate.meta.result[["gamma.s.plot"]]$forest.plot

p3

ggsave(
  filename = "rise_influenzain_crossstudy_gender.pdf",
  path = application_figures_folder,
  plot = p3,
  width = 27,
  height = 17,
  units = "cm"
)

# Analysis four : across all studies, stratified by age
hipc_merged_all_noNorm_filtered_copy = hipc_merged_all_noNorm_filtered %>%
  mutate(study_accession = ifelse(
    age_imputed > 64,
    paste0(study_accession, "_Older"),
    paste0(study_accession, "_Younger")
  ))

df_train = hipc_merged_all_noNorm_filtered_copy %>%
  group_by(study_accession) %>%
  ungroup()

yone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(response_post)

yzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(response_pre)

sone = df_train %>%
  filter(study_time_collected > 0) %>%
  select(all_of(markers))

szero = df_train %>%
  filter(study_time_collected == 0) %>%
  select(all_of(markers))

studyone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession)

studyzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(study_accession)

rise.evaluate.meta.result = rise.evaluate.meta(
  yone,
  yzero,
  sone,
  szero,
  studyone,
  studyzero,
  alpha = 0.05,
  epsilon.meta = 0.2,
  alternative = "two.sided",
  paired.all = T,
  evaluate.weights = T,
  return.all.evaluate = T,
  epsilon.study = 0.2,
  p.correction = "BH",
  screening.weights = screening.weights,
  markers = markers,
  show.pooled.effect = F
)

p4 = rise.evaluate.meta.result[["gamma.s.plot"]]$forest.plot

p4

ggsave(
  filename = "rise_influenzain_crossstudy_age.pdf",
  path = application_figures_folder,
  plot = p4,
  width = 27,
  height = 14,
  units = "cm"
)


# Analysis five : Overlap between study signatures at day 1
hipc_merged_all_noNorm_filtered_copy = hipc_merged_all_noNorm_filtered

df_train = hipc_merged_all_noNorm_filtered_copy %>%
  group_by(study_accession) %>%
  ungroup()

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

# sone = df_train %>%
#   filter(study_time_collected > 0) %>%
#   select(all_of(markers))
#
# szero = df_train %>%
#   filter(study_time_collected == 0) %>%
#   select(all_of(markers))

studyone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession)

studyzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(study_accession)

rise.screen.meta.result = rise.screen.meta(
  yone,
  yzero,
  sone,
  szero,
  studyone,
  studyzero,
  alpha = 0.05,
  epsilon.meta = 0.2,
  alternative = "two.sided",
  paired.all = T,
  return.all.screen = T,
  epsilon.study = 0.2,
  p.correction = "BH",
  show.pooled.effect = F,
  return.study.similarity.plot = TRUE
)

p5 = rise.screen.meta.result[["gamma.s.plot"]]$similarity.plots$upset.plot

p5

ggsave(
  filename = "rise_influenzain_crossstudy_similarities.pdf",
  path = application_figures_folder,
  plot = p5,
  width = 27,
  height = 17,
  units = "cm"
)

markers_study = rise.screen.meta.result[["screening.metrics.study"]]

# number of studies
n_studies <- n_distinct(markers_study$study)

# markers significant in all studies
markers_all_studies <- markers_study %>%
  filter(p_adjusted < 0.05) %>%
  group_by(marker) %>%
  summarise(n_studies_sig = n_distinct(study), .groups = "drop") %>%
  filter(n_studies_sig == n_studies) %>%
  pull(marker)

# Analysis six: RISE-meta applied to all genes and studies

hipc_merged_all_noNorm_filtered_copy = hipc_merged_all_noNorm_filtered

df_train = hipc_merged_all_noNorm_filtered_copy %>%
  group_by(study_accession) %>%
  ungroup()

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

studyone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession)

studyzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(study_accession)

rise.screen.meta.result = rise.screen.meta(
  yone,
  yzero,
  sone,
  szero,
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
  return.study.similarity.plot = TRUE
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
  select(marker, mu.delta.ci, p.unadjusted, p.adjusted) %>%
  head(10)

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


# Analysis seven: RISE-meta applied to 222 genes and studies

hipc_merged_all_noNorm_filtered_copy = hipc_merged_all_noNorm_filtered

df_train = hipc_merged_all_noNorm_filtered_copy %>%
  group_by(study_accession) %>%
  ungroup()

yone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(response_post)

yzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(response_pre)

sone = df_train %>%
  filter(study_time_collected > 0) %>%
  select(all_of(markers))

szero = df_train %>%
  filter(study_time_collected == 0) %>%
  select(all_of(markers))

studyone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession)

studyzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(study_accession)

rise.screen.meta.result = rise.screen.meta(
  yone,
  yzero,
  sone,
  szero,
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
  return.study.similarity.plot = TRUE
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


p6 = rise.screen.meta.result[["gamma.s.plot"]]$forest.plot

p6

ggsave(
  filename = "rise_influenzain_222_meta.pdf",
  path = application_figures_folder,
  plot = p6,
  width = 27,
  height = 17,
  units = "cm"
)

p7 = rise.screen.meta.result[["gamma.s.plot"]]$similarity.plots$upset.plot

p7

ggsave(
  filename = "rise_influenzain_222_crossstudy_similarities.pdf",
  path = application_figures_folder,
  plot = p7,
  width = 27,
  height = 17,
  units = "cm"
)

p13 = rise.screen.meta.result[["gamma.s.plot"]]$fit.plot

p13

ggsave(
  filename = "rise_influenzain_222_crossstudy_meta_fitplot.pdf",
  path = application_figures_folder,
  plot = p13,
  width = 37,
  height = 20,
  units = "cm"
)

# Analysis eight: pre-selecting BTM interferon/antiviral genes

hipc_merged_all_noNorm_filtered_copy = hipc_merged_all_noNorm_filtered

df_train = hipc_merged_all_noNorm_filtered_copy %>%
  group_by(study_accession) %>%
  ungroup()

yone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(response_post)

yzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(response_pre)

sone = df_train %>%
  filter(study_time_collected > 0) %>%
  select(any_of(BTM_filtered))

szero = df_train %>%
  filter(study_time_collected == 0) %>%
  select(any_of(BTM_filtered))

studyone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession)

studyzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(study_accession)

rise.screen.meta.result = rise.screen.meta(
  yone,
  yzero,
  sone,
  szero,
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
  return.study.similarity.plot = TRUE
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


p8 = rise.screen.meta.result[["gamma.s.plot"]]$forest.plot

p8

ggsave(
  filename = "rise_influenzain_inteferon_meta.pdf",
  path = application_figures_folder,
  plot = p8,
  width = 27,
  height = 17,
  units = "cm"
)

p9 = rise.screen.meta.result[["gamma.s.plot"]]$similarity.plots$upset.plot

p9

ggsave(
  filename = "rise_influenzain_inteferon_crossstudy_similarities.pdf",
  path = application_figures_folder,
  plot = p9,
  width = 27,
  height = 17,
  units = "cm"
)

# Analysis nine: BTM genes from innate immune pathways
innate_pathways = c(
  "Interferon/Antiviral Sensing",
  "Monocytes",
  "Antigen Presentation",
  "NK Cells",
  "Neutrophils",
  "Inflammatory/TLR/Chemokines",
  "DC Activation"
)

BTM_innate = which(BTM[["geneset.aggregates"]] %in% innate_pathways)
BTM_filtered = BTM[["genesets"]][BTM_innate] %>%
  unlist()

length(BTM_filtered)


hipc_merged_all_noNorm_filtered_copy = hipc_merged_all_noNorm_filtered

df_train = hipc_merged_all_noNorm_filtered_copy %>%
  group_by(study_accession) %>%
  ungroup()

yone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(response_post)

yzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(response_pre)

sone = df_train %>%
  filter(study_time_collected > 0) %>%
  select(any_of(BTM_filtered))

szero = df_train %>%
  filter(study_time_collected == 0) %>%
  select(any_of(BTM_filtered))

studyone = df_train %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession)

studyzero = df_train %>%
  filter(study_time_collected == 0) %>%
  pull(study_accession)

rise.screen.meta.result = rise.screen.meta(
  yone,
  yzero,
  sone,
  szero,
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
  return.study.similarity.plot = TRUE
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


# Analysis ten: geneset-level features

hipc_merged_all_noNorm_filtered_copy = hipc_merged_all_noNorm_filtered

innate_pathways = c(
  "Interferon/Antiviral Sensing",
  "Monocytes",
  "Antigen Presentation",
  "NK Cells",
  "Neutrophils",
  "Inflammatory/TLR/Chemokines",
  "DC Activation"
)

BTM_genes = BTM[["genesets"]][which(BTM[["geneset.aggregates"]] != "NA")] %>%
  unlist()

BTM_genes_names = BTM[["geneset.names.descriptions"]][which(BTM[["geneset.aggregates"]] != "NA")]

df_train = hipc_merged_all_noNorm_filtered_copy %>%
  group_by(study_accession) %>%
  ungroup()

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
sone_geneset <- aggregate_to_geneset(
  df = sone,
  genesets = BTM[["genesets"]][which(BTM[["geneset.aggregates"]] != "NA")],
  geneset_names = BTM_genes_names,
  FUN = mean
)

szero_geneset <- aggregate_to_geneset(
  df = szero,
  genesets = BTM[["genesets"]][which(BTM[["geneset.aggregates"]] != "NA")],
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
  sone = sone_geneset,
  szero = szero_geneset,
  studyone,
  studyzero,
  alpha = 0.05,
  epsilon.meta = 0.2,
  alternative = "two.sided",
  paired.all = T,
  return.all.screen = T,
  epsilon.study = 0.25,
  p.correction = "BH",
  show.pooled.effect = T,
  return.study.similarity.plot = TRUE
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

p10 = rise.screen.meta.result[["gamma.s.plot"]]$forest.plot

p10

ggsave(
  filename = "rise_influenzain_BTMmean_meta.pdf",
  path = application_figures_folder,
  plot = p10,
  width = 27,
  height = 17,
  units = "cm"
)

p11 = rise.screen.meta.result[["gamma.s.plot"]]$similarity.plots$upset.plot

p11

ggsave(
  filename = "rise_influenzain_BTMmean_crossstudy_similarities.pdf",
  path = application_figures_folder,
  plot = p11,
  width = 27,
  height = 17,
  units = "cm"
)

p12 = rise.screen.meta.result[["gamma.s.plot"]]$fit.plot

p12

ggsave(
  filename = "rise_influenzain_BTMmean_crossstudy_meta_fitplot.pdf",
  path = application_figures_folder,
  plot = p12,
  width = 37,
  height = 20,
  units = "cm"
)
