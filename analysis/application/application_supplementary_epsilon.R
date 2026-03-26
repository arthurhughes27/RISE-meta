# Script to run supplementary application exploring the impact of the choice of epsilon


# Define global hyperparameters for analysis
hyperparameter_list = list(
  # Hyperparameters for data pre-processing
  tp = 1,
  # Timepoint for gene expression
  screen.fraction = 1,
  # Fraction of data for screening
  seed = 08012025,
  # seed for random data splitting
  
  # Hyperparameters to define methodology
  meta.analysis.method = "RE",
  # meta analysis method (random or fixed effects)
  test = "knha",
  # method for variance estimation of pooled effect
  alternative = "two.sided",
  # form of alternative hypothesis
  epsilon.meta.mode = "mean.power",
  # choice of how to define epsilon
  paired.all = TRUE,
  # paired mode
  paired.studies = NULL,
  # which studies are paired
  evaluate.weights = TRUE, 
  # Whether to use weighting for evaluation stage
  
  # Numeric hyperparameters for testing procedure
  alpha = 0.05,
  # significance level
  power.want.s.study = 0.8,
  # within-study power for epsilon
  epsilon.meta = NULL,
  # fixed value for epsilon
  epsilon.study = NULL,
  # epsilon for within-study testing
  p.correction = "BH",
  # multiplicity correction for p-values
  u.y.hyp = NULL,
  # hypothesised effect size on y
  weight.mode = "diff.epsilon",
  # How to weight surrogates in combination
  normalise.weights = TRUE,
  # normalise weights for the combination
  
  # Hyperparameters to define which objects to return
  return.all.screen = TRUE,
  # returns all screening results
  show.pooled.effect = TRUE,
  # show pooled effect in forest plot?
  return.study.similarity.plot = FALSE,
  # return similarity of within-study analyses plot
  return.forest.plot = TRUE,
  # return forest plot for combined marker
  return.fit.plot = TRUE,
  # return fit plot for combined marker
  return.evaluate.results = TRUE,
  # return evaluation results for screening data
  return.screen.plot = TRUE, 
  # return screening plot
  return.all.weights = FALSE,
  # return weights for all predictors
  
  # Predictor transformation parameters
  aggregation_function = mean,
  # Function defining aggregation from gene to geneset level 
  geneset_definition = "BTM",
  # Argument stating the definition of the genesets (options are BTM or BG3M)
  
  
  # Other hyperparameters
  n.cores = 5,
  # number of cores for parallel computing
  screen.plot.topN = 15, # how many predictors to plot
  
  # Graphical parameters
  screen.plot.width = 40,
  screen.plot.height = 18,
  forest.plot.width = 32,
  forest.plot.height = 15,
  fit.plot.width = 37,
  fit.plot.height = 20
)

file_name_tag = paste0("_timepoint",
                       hyperparameter_list$tp,
                       "_method",
                       hyperparameter_list$meta.analysis.method,
                       "_test",
                       hyperparameter_list$test,
                       "_epsMode",
                       hyperparameter_list$epsilon.meta.mode,
                       ifelse(hyperparameter_list$epsilon.meta.mode == "user", 
                              paste0("_eps", hyperparameter_list$epsilon.meta), 
                              paste0("_power", hyperparameter_list$power.want.s.study)))


# Libraries
library(tidyverse)
library(SurrogateRank)

# Load internal functions
sapply(list.files("R/", pattern = "\\.R$", full.names = TRUE), source)

# Paths to processed data and output figures
processed_data_folder <- "data"
application_figures_folder <- fs::path("output", "figures", "application", "supplementary", "epsilon")
application_results_folder <- fs::path("output", "results", "application", "supplementary", "epsilon")


# Load merged gene expression and GS_list gene set objects
df <- readRDS(fs::path(processed_data_folder, "hipc_merged_all_noNorm.rds"))
GS_list <- readRDS(fs::path(processed_data_folder, paste0(hyperparameter_list$geneset_definition, "_processed.rds")))

preprocessed_data_list = preprocess_data(
  df = df,
  tp = hyperparameter_list$tp,
  screen.fraction = hyperparameter_list$screen.fraction,
  seed = hyperparameter_list$seed
)

df_train = preprocessed_data_list[["df.screen"]]
df_test = preprocessed_data_list[["df.evaluate"]]

predictor_names = df_train %>% 
  dplyr::select(a1cf:zzz3) %>% 
  colnames()

# ----- Screening on training data -----

train_inputs <- extract_rise_inputs(df_train, 
                                    predictor_names = predictor_names, 
                                    genesets = GS_list[["genesets"]], 
                                    geneset_names = GS_list[["geneset.names.descriptions"]],
                                    aggregation_function = hyperparameter_list$aggregation_function)

# ----- Run screening for the four requested specifications -----
spec_grid <- tibble::tribble(
  ~epsilon.meta.mode, ~epsilon.meta, ~power.want.s.study,
  "user",             0.1,          NA_real_,
  "user",             0.2,          NA_real_,
  "mean.power",       NA_real_,     0.8,
  "mean.power",       NA_real_,     0.9
)

run_one_spec <- function(epsilon.meta.mode, epsilon.meta, power.want.s.study) {
  
  if (length(power.want.s.study) == 1 && is.na(power.want.s.study)) {
    power.want.s.study <- NULL
  }
  
  rise_screen_result <- rise.screen.meta(
    yone                         = train_inputs$yone,
    yzero                        = train_inputs$yzero,
    sone                         = train_inputs$sone,
    szero                        = train_inputs$szero,
    studyone                     = train_inputs$studyone,
    studyzero                    = train_inputs$studyzero,
    alpha                        = hyperparameter_list$alpha,
    epsilon.meta.mode            = epsilon.meta.mode,
    power.want.s.study           = power.want.s.study,
    epsilon.meta                 = epsilon.meta,
    alternative                  = hyperparameter_list$alternative,
    paired.all                   = hyperparameter_list$paired.all,
    return.all.screen            = hyperparameter_list$return.all.screen,
    epsilon.study                = hyperparameter_list$epsilon.study,
    p.correction                 = hyperparameter_list$p.correction,
    show.pooled.effect           = hyperparameter_list$show.pooled.effect,
    return.study.similarity.plot = hyperparameter_list$return.study.similarity.plot,
    test                         = hyperparameter_list$test,
    meta.analysis.method         = hyperparameter_list$meta.analysis.method,
    n.cores                      = hyperparameter_list$n.cores,
    screen.plot.topN             = hyperparameter_list$screen.plot.topN,
    return.evaluate.results      = hyperparameter_list$return.evaluate.results,
    return.fit.plot              = hyperparameter_list$return.fit.plot,
    return.forest.plot           = hyperparameter_list$return.forest.plot, 
    normalise.weights            = hyperparameter_list$normalise.weights, 
    return.screen.plot           = hyperparameter_list$return.screen.plot, 
    weight.mode                  = hyperparameter_list$weight.mode, 
    return.all.weights           = hyperparameter_list$return.all.weights, 
    paired.studies               = hyperparameter_list$paired.studies, 
    u.y.hyp                      = hyperparameter_list$u.y.hyp
  )
  
  screen_output = extract_rise_outputs(screen_result = rise_screen_result)
  
  file_name_tag_spec = paste0("_timepoint",
                              hyperparameter_list$tp,
                              "_method",
                              hyperparameter_list$meta.analysis.method,
                              "_test",
                              hyperparameter_list$test,
                              "_epsMode",
                              epsilon.meta.mode,
                              ifelse(epsilon.meta.mode == "user", 
                                     paste0("_eps", epsilon.meta), 
                                     paste0("_power", power.want.s.study)))
  
  screen_plot_1 = screen_output$screen_plot
  screen_forest_1 = screen_output$screen_forest
  screen_fit_1 = screen_output$screen_fit
  
  if (!is.null(screen_plot_1)) {
    ggsave(
      filename = paste0("screen_plot", file_name_tag_spec, ".pdf"),
      path     = application_figures_folder,
      plot     = screen_plot_1,
      width    = hyperparameter_list$screen.plot.width,
      height   = hyperparameter_list$screen.plot.height,
      units    = "cm"
    )
  }
  
  if (!is.null(screen_forest_1)) {
    ggsave(
      filename = paste0("screen_forest", file_name_tag_spec, ".pdf"),
      path     = application_figures_folder,
      plot     = screen_forest_1,
      width    = hyperparameter_list$forest.plot.width,
      height   = hyperparameter_list$forest.plot.height,
      units    = "cm"
    )
  }
  
  if (!is.null(screen_fit_1)) {
    ggsave(
      filename = paste0("screen_fit", file_name_tag_spec, ".pdf"),
      path     = application_figures_folder,
      plot     = screen_fit_1,
      width    = hyperparameter_list$fit.plot.width,
      height   = hyperparameter_list$fit.plot.height,
      units    = "cm"
    )
  }
  
  tibble::tibble(
    epsilon.meta.mode = epsilon.meta.mode,
    epsilon.meta = epsilon.meta,
    epsilon.used = screen_output$epsilon,
    power.want.s.study = power.want.s.study,
    n_significant = screen_output$n_significant
  )
}

results_df <- purrr::pmap_dfr(spec_grid, run_one_spec)

results_df

rm(list = ls())
