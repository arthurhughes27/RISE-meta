# R script to perform meta-analysis with rise
# Modified to produce two separate forest plots:
#   - one for the Fixed-Effect (FE) inverse-variance meta-analysis
#   - one for the Random-Effects (RE) REML + Knapp-Hartung meta-analysis
# Added display of tau^2 and I^2 on the RE plot.
library(tidyverse)
library(SurrogateRank)
library(grid)    # for unit()
library(scales)  # pretty formatting
library(cowplot)
library(metafor) # required for REML random-effects fit

# Directory containing engineered / processed data files
processed_data_folder <- "data"

# Paths to processed gene-level data and gene-set objects
p_load_expr_all_noNorm <- fs::path(processed_data_folder, "hipc_merged_all_noNorm.rds")

# Load data objects
hipc_merged_all_noNorm <- readRDS(p_load_expr_all_noNorm)

# Gene and timepoint to illustrate
gene_name = "jchain"
tp = 7

# Timepoints of interest (numeric)
timepoints_to_keep <- c(0,tp)

# Filter to samples with non-missing immune response, Influenza vaccine,
# and collected at one of the specified timepoints.
hipc_merged_all_noNorm_filtered <- hipc_merged_all_noNorm %>%
  filter(
    !is.na(immResp_MFC_anyAssay_log2_MFC),
    vaccine_name == "Influenza (IN)",
    study_time_collected %in% timepoints_to_keep
  ) %>% 
  select(participant_id, 
         study_accession,
         study_time_collected,
         immResp_MFC_anyAssay_pre_value,
         immResp_MFC_anyAssay_post_value,
         all_of(gene_name))

# Extract the study names
study_names = hipc_merged_all_noNorm_filtered %>%
  filter(study_time_collected > 0) %>%
  pull(study_accession) %>%
  unique()

timepoints_of_interest = c(tp)

epsilon = 0.2

for (sdy in study_names) {
  # sdy = study_names[1]
  
  hipc_merged_all_noNorm_filtered_sdy = hipc_merged_all_noNorm_filtered %>%
    filter(study_time_collected %in% c(0, tp),
           study_accession == sdy) %>%
    group_by(participant_id) %>%
    filter(all(c(0, tp) %in% study_time_collected)) %>%
    ungroup()
  
  # Pre-vaccination response
  yzero_sdy = hipc_merged_all_noNorm_filtered_sdy %>%
    filter(study_time_collected == 0) %>%
    pull(immResp_MFC_anyAssay_pre_value)
  
  # Post-vaccination response
  yone_sdy = hipc_merged_all_noNorm_filtered_sdy %>%
    filter(study_time_collected == tp) %>%
    pull(immResp_MFC_anyAssay_post_value)
  
  # Pre-vaccination gene expression
  szero_sdy = hipc_merged_all_noNorm_filtered_sdy %>%
    filter(study_time_collected == 0) %>%
    pull(gene_name)
  
  # Post-vaccination gene expression
  sone_sdy = hipc_merged_all_noNorm_filtered_sdy %>%
    filter(study_time_collected == tp) %>%
    pull(gene_name)
  
  test_result_sdy = SurrogateRank::test.surrogate.extension(
    yone = yone_sdy,
    yzero = yzero_sdy,
    sone = sone_sdy,
    szero = szero_sdy,
    epsilon = epsilon,
    alternative = "two.sided",
    paired = T,
    alpha = 0.05
  )
  
  test_result_sdy$ci.delta.lower <- test_result_sdy$ci.delta[1]
  test_result_sdy$ci.delta.upper <- test_result_sdy$ci.delta[2]
  test_result_sdy$ci.delta <- NULL
  
  test_result_df = data.frame(test_result_sdy)
  test_result_df$study_accession = sdy
  test_result_df$n = 2*length(yone_sdy)
  
  test_result_df = test_result_df %>%
    select(
      study_accession,
      n,
      epsilon.used,
      u.y,
      u.s,
      delta.estimate,
      ci.delta.lower,
      ci.delta.upper,
      sd.delta,
      p.delta
    )
  
  if (sdy == study_names[1]) {
    all_results = test_result_df
  } else {
    all_results = rbind(all_results, test_result_df)
  }
}


# ----------------- 0. Data prep and pooled estimates ------------------------
# Use your all_results dataframe
df <- all_results %>%
  filter(sd.delta != 0)

meta_reml_tost <- function(df = NULL,
                           delta = NULL,
                           sd_delta = NULL,
                           vi = NULL,
                           epsilon = NULL,
                           alpha = 0.05,
                           tol = 1e-10,
                           verbose = FALSE) {
  # INPUTS:
  #  - df: optional data.frame with columns delta.estimate and sd.delta (or vi)
  #  - delta: numeric vector of study estimates (overrides df if provided)
  #  - sd_delta: numeric vector of standard errors (overrides df if provided)
  #  - vi: optional vector of variances (if provided, used directly)
  #  - epsilon: equivalence margin (required)
  #  - alpha: significance level for CIs and TOST
  #
  # OUTPUT: list with elements: tau2, mu, se_HK, ci_HK, se_conv, ci_conv,
  #         p_L, p_U, p_TOST, Q, I2, m, weights, details
  
  if (!is.null(df)) {
    if (!is.null(delta) || !is.null(sd_delta) || !is.null(vi)) {
      warning("df provided and also delta/sd_delta/vi provided. Using explicit vectors (delta/sd_delta/vi) if present.")
    }
    # prefer explicit args when present
    if (is.null(delta)) {
      if (!is.null(df$delta.estimate)) delta <- df$delta.estimate
      else stop("df provided but no column delta.estimate found and delta not supplied")
    }
    if (is.null(vi) && is.null(sd_delta)) {
      if (!is.null(df$sd.delta)) sd_delta <- df$sd.delta
      else if (!is.null(df$vi)) vi <- df$vi
      else stop("df provided but no sd.delta or vi column found and sd_delta/vi not supplied")
    }
  }
  
  if (is.null(delta)) stop("delta (trial-level estimates) must be supplied (either via df or delta argument)")
  delta <- as.numeric(delta)
  m <- length(delta)
  
  if (!is.null(vi)) {
    vi <- as.numeric(vi)
  } else {
    if (is.null(sd_delta)) stop("Either sd_delta or vi must be supplied (or df with sd.delta)")
    sd_delta <- as.numeric(sd_delta)
    vi <- sd_delta^2
  }
  
  if (length(vi) != m) stop("length(vi) must equal length(delta)")
  if (any(is.na(delta)) || any(is.na(vi))) stop("NA values in delta or vi are not allowed")
  if (is.null(epsilon)) stop("epsilon (equivalence margin) must be supplied")
  
  # -------------------------
  # REML estimation of tau^2
  # -------------------------
  # weighted mean function for a given tau2
  mu_of <- function(tau2) {
    w <- 1 / (vi + tau2)
    sum(w * delta) / sum(w)
  }
  
  # restricted score function: root gives REML tau^2
  score_fn <- function(tau2) {
    v <- vi + tau2
    w <- 1 / v
    mu <- sum(w * delta) / sum(w)
    sum((delta - mu)^2 / (v^2)) - sum(1 / v)
  }
  
  # negative restricted log-likelihood (minimize this if no root found)
  neg_restricted_loglik <- function(tau2) {
    v <- vi + tau2
    wsum <- sum(1 / v)
    mu  <- sum((1 / v) * delta) / wsum
    0.5 * ( sum(log(v)) + log(wsum) + sum((delta - mu)^2 / v) )
  }
  
  # choose an upper bound for tau2 search
  s2 <- if (m > 1) var(delta) else 0
  upper <- max(abs(s2 - mean(vi)), s2, max(vi), 1e-6)
  upper <- abs(upper) * 10 + 1e-6
  if (upper <= 0) upper <- 1e-6
  
  # try uniroot on [0, upper] if a sign change exists
  tau2_hat <- NA
  method_used <- NA
  converged <- FALSE
  details <- list()
  
  s0 <- score_fn(0)
  su <- score_fn(upper)
  if (verbose) cat(sprintf("REML score at endpoints: score(0)=%.6g, score(upper=%.6g)=%.6g\n", s0, upper, su))
  
  if (is.finite(s0) && is.finite(su) && s0 * su < 0) {
    ur <- tryCatch(
      uniroot(score_fn, lower = 0, upper = upper, tol = tol),
      error = function(e) e
    )
    if (!inherits(ur, "error")) {
      tau2_hat <- max(0, ur$root)
      method_used <- "score_uniroot"
      converged <- TRUE
      details$uniroot <- ur
    } else {
      if (verbose) message("uniroot failed; falling back to optimize on restricted log-likelihood")
    }
  }
  
  if (!converged) {
    opt <- optimize(neg_restricted_loglik, lower = 0, upper = upper, tol = tol)
    tau2_hat <- max(0, opt$minimum)
    method_used <- "optimize_neg_restricted_loglik"
    converged <- TRUE
    details$optimize <- opt
  }
  
  # -------------------------
  # pooled estimates & HK SE
  # -------------------------
  w_tau <- 1 / (vi + tau2_hat)
  mu_hat <- sum(w_tau * delta) / sum(w_tau)
  var_conv <- 1 / sum(w_tau)
  se_conv <- sqrt(var_conv)
  
  # Hartung-Knapp scaling factor q (as in your formula)
  if (m > 1) {
    q <- (1 / (m - 1)) * sum(w_tau * (delta - mu_hat)^2)
    # ensure q non-negative (numerical)
    q <- max(0, q)
    se_HK <- sqrt(q) * sqrt(var_conv)
  } else {
    q <- NA
    se_HK <- NA
  }
  
  # HK CI (uses t_{m-1})
  if (m > 1) {
    tcrit <- qt(1 - alpha / 2, df = m - 1)
    ci_HK <- c(mu_hat - tcrit * se_HK, mu_hat + tcrit * se_HK)
  } else {
    ci_HK <- c(NA, NA)
  }
  
  # Conventional normal-approx CI (for reference)
  zcrit <- qnorm(1 - alpha / 2)
  ci_conv <- c(mu_hat - zcrit * se_conv, mu_hat + zcrit * se_conv)
  
  # -------------------------
  # Two-one-sided test (TOST) using HK SE
  # -------------------------
  if (m > 1) {
    T_L <- (mu_hat + epsilon) / se_HK
    T_U <- (mu_hat - epsilon) / se_HK
    p_L <- 1 - pt(T_L, df = m - 1)   # p-value for H0L: mu <= -epsilon
    p_U <- pt(T_U, df = m - 1)       # p-value for H0U: mu >= +epsilon
    p_TOST <- max(p_L, p_U)
  } else {
    T_L <- T_U <- p_L <- p_U <- p_TOST <- NA
  }
  
  # -------------------------
  # Cochran's Q and I^2
  # -------------------------
  # Fixed-effect weights use vi only
  w0 <- 1 / vi
  w0_sum <- sum(w0)
  delta_FE <- sum(w0 * delta) / w0_sum
  Q <- sum(w0 * (delta - delta_FE)^2)
  # handle degenerate Q
  if (m > 1 && is.finite(Q) && Q > (m - 1)) {
    I2 <- max(0, (Q - (m - 1)) / Q) * 100
  } else {
    I2 <- max(0, (Q - (m - 1)) / max(Q, 1e-12)) * 100
    I2 <- max(0, I2)
  }
  
  # -------------------------
  # prepare output
  # -------------------------
  results <- list(
    m = m,
    tau2 = tau2_hat,
    mu = mu_hat,
    var_conv = var_conv,
    se_conv = se_conv,
    ci_conv = ci_conv,
    se_HK = se_HK,
    ci_HK = ci_HK,
    q_HK = q,
    T_L = T_L,
    T_U = T_U,
    p_L = p_L,
    p_U = p_U,
    p_TOST = p_TOST,
    Q = Q,
    I2 = I2,
    weights_tau = w_tau,
    weights_fe = w0,
    delta_FE = delta_FE,
    method = method_used,
    converged = converged,
    details = details
  )
  
  # small summary table per-study (optional)
  study_tbl <- data.frame(
    delta = delta,
    vi = vi,
    w_tau = w_tau,
    stringsAsFactors = FALSE
  )
  
  return(list(results = results, per_study = study_tbl))
}

out <- meta_reml_tost(df = df, epsilon = epsilon, alpha = 0.05, verbose = TRUE)

# ------------------ Validate inputs ------------------------------------------------
if (!exists("df")) stop("data.frame `df` not found. Provide df with study-level data.")
if (!exists("out")) stop("meta results `out` not found. Run meta_reml_tost(...) first.")
if (is.null(epsilon)) stop("equivalence margin `epsilon` is not defined in the environment.")

# ------------------ Prepare per-study data ----------------------------------------
# Ensure numeric columns present
df <- df %>% mutate(delta = as.numeric(delta.estimate))

# derive se and vi
if (!is.null(df$vi)) {
  df <- df %>% mutate(vi = as.numeric(vi),
                      sd = sqrt(vi))
} else if (!is.null(df$sd.delta)) {
  df <- df %>% mutate(sd = as.numeric(sd.delta),
                      vi = (sd)^2)
} else {
  stop("df must contain either 'vi' or 'sd.delta' (standard error) for each study.")
}

# per-study CI fallback if missing: normal approx 95% CI
if (is.null(df$ci.delta.lower) || any(is.na(df$ci.delta.lower))) {
  df <- df %>%
    mutate(ci.delta.lower = ifelse(is.na(ci.delta.lower),
                                   delta - qnorm(0.975) * sd,
                                   ci.delta.lower),
           ci.delta.upper = ifelse(is.na(ci.delta.upper),
                                   delta + qnorm(0.975) * sd,
                                   ci.delta.upper))
}

# per-study p-value: prefer provided p.delta, otherwise compute per-study TOST p (normal approx if n missing)
df <- df %>%
  mutate(
    p_L_study = (delta + epsilon) / sd,                    # T_L statistic (z or t)
    p_U_study = (delta - epsilon) / sd,
    pL = ifelse(!is.null(p.delta) & !is.na(p.delta), NA_real_, 1 - pnorm(p_L_study)),
    pU = ifelse(!is.null(p.delta) & !is.na(p.delta), NA_real_, pnorm(p_U_study)),
    p_TOST_study = ifelse(!is.null(p.delta) & !is.na(p.delta), p.delta, pmax(pL, pU)),
    label_pval = ifelse(is.na(p_TOST_study), "", formatC(p_TOST_study, format = "f", digits = 3)),
    label_n = ifelse(!is.null(n) & !is.na(n), as.character(n), "")
  )

# assign vertical positions: studies from top (k) down to 1; summary will be at y = 0
k <- nrow(df)
df <- df %>% mutate(y = seq(from = k, to = 1))

# ------------------ Extract RE summary from out -----------------------------------
res <- out$results
# res <- rma.uni(yi = df$delta, vi = df$vi, method = "REML", test = "knha")

mu_re       <- res$mu
ci_re_lower <- res$ci_HK[1]
ci_re_upper <- res$ci_HK[2]
se_hk       <- res$se_HK
var_hk      <- if (!is.na(se_hk)) se_hk^2 else NA_real_
p_tost_re   <- res$p_TOST
tau2_reml   <- res$tau2
I2_pct      <- res$I2
k_used      <- res$m

# ------------------ Random-effects weights & percent weights ----------------------
# compute RE weights using estimated tau2
w_re <- 1 / (df$vi + tau2_reml)
pct_w <- 100 * w_re / sum(w_re, na.rm = TRUE)

df <- df %>%
  mutate(weight_RE = w_re,
         pct_weight = pct_w,
         label_wgt = ifelse(!is.na(pct_weight), formatC(pct_weight, format = "f", digits = 1), "")
  )

# ------------------ Prepare summary (RE) row -------------------------------------
summary_row_re <- tibble::tibble(
  study_accession = "Summary",
  delta.estimate  = mu_re,
  ci.delta.lower  = ci_re_lower,
  ci.delta.upper  = ci_re_upper,
  sd.delta        = se_hk,
  p.delta         = p_tost_re,
  epsilon         = epsilon,
  var             = var_hk,
  weight          = NA_real_,
  pct_weight      = NA_real_,
  label_pval      = ifelse(is.na(p_tost_re), "", formatC(p_tost_re, format = "f", digits = 3)),
  label_wgt       = "",
  label_n         = "",
  y = 0
)

# Combine the study data with the summary row without renaming
plot_df_re <- bind_rows(
  df,             # df already has delta.estimate, ci.delta.lower, ci.delta.upper
  summary_row_re  # summary row
) %>%
  mutate(
    study_label = study_accession,
    is_summary  = (study_label == "Summary")
  )

# ------------------ Plot layout parameters --------------------------------------
base_text_size <- 16   # large, readable text
y_min <- -1
y_max <- k + 1
rel_w_left  <- 0.45
rel_w_mid   <- 1.10
rel_w_right <- 0.55

# x axis range
x_min <- -1
x_max <- 1

# ------------------ Left panel: study labels -------------------------------------
left_labels_re <- ggplot(plot_df_re, aes(y = y)) +
  geom_text(aes(x = 0, label = study_label,
                fontface = ifelse(is_summary, "bold", "plain")), 
            hjust = 0, size = 5) +
  scale_x_continuous(limits = c(-0.1, 0.95), expand = c(0, 0)) +
  scale_y_continuous(breaks = plot_df_re$y, limits = c(y_min, y_max), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = unit(c(1.2, 0.6, 1, 1.2), "lines"),
    plot.title = element_text(hjust = 0.5, size = base_text_size + 4, face = "bold")
  )

# ------------------ Middle panel: forest plot -----------------------------------
forest_mid_re <- ggplot(plot_df_re, aes(x = delta.estimate, y = y)) +
  # CIs
  geom_errorbarh(aes(xmin = ci.delta.lower, xmax = ci.delta.upper), height = 0.15, size = 0.8) +
  # study points
  geom_point(data = filter(plot_df_re, !is_summary), shape = 16, size = 3.5) +
  # summary point (distinct)
  geom_point(data = filter(plot_df_re, is_summary), shape = 5, size = 8) +
  # axis and limits
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0), breaks = seq(x_min, x_max, by = 0.5)) +
  scale_y_continuous(breaks = plot_df_re$y, limits = c(y_min, y_max), expand = c(0, 0)) +
  coord_cartesian(ylim = c(y_min, y_max), xlim = c(x_min, x_max), clip = "off") +
  labs(x = expression(delta), y = NULL,
       title = paste0("Random-effects meta-analysis of ", ifelse(exists("gene_name"), gene_name, ""), 
                      ifelse(exists("tp"), paste0(" at Day ", tp), ""))) +
  theme_minimal(base_size = base_text_size) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = base_text_size + 6),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = base_text_size + 4),
    axis.text.x = element_text(size = base_text_size),
    plot.margin = unit(c(1.2, 0.6, 1, 0.6), "lines")
  ) +
  # non-inferiority margins and center line
  geom_vline(xintercept = c(-epsilon, epsilon), linetype = "dashed", color = "red", size = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey60")

# ------------------ Right panel: p-value / weight / N ---------------------------
# prepare a display dataframe for the right panel aligned by y
right_df <- plot_df_re %>%
  mutate(
    display_pval = ifelse(is_summary, label_pval, ifelse(!is.na(label_pval) & label_pval != "", label_pval, "")),
    display_wgt = ifelse(is_summary, "", ifelse(!is.na(label_wgt) & label_wgt != "", paste0(label_wgt, "%"), "")),
    display_n   = ifelse(is_summary, "", label_n)
  )

right_table_re <- ggplot(right_df, aes(y = y)) +
  # column headers
  annotate("text", x = 1, y = k + 0.9, label = "p-value", fontface = "bold", hjust = 0, size = 5) +
  annotate("text", x = 2.1, y = k + 0.9, label = "Weight", fontface = "bold", hjust = 0, size = 5) +
  annotate("text", x = 3.3, y = k + 0.9, label = "N", fontface = "bold", hjust = 0, size = 5) +
  geom_text(aes(x = 1, label = display_pval), hjust = 0, size = 4.5) +
  geom_text(aes(x = 2.1, label = display_wgt), hjust = 0, size = 4.5) +
  geom_text(aes(x = 3.3, label = display_n), hjust = 0, size = 4.5) +
  scale_y_continuous(breaks = right_df$y, limits = c(y_min, y_max), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0.9, 4.0), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = unit(c(1.2, 1.2, 1, 0.6), "lines"))

# ------------------ Combine panels ------------------------------------------------
combined_re <- plot_grid(left_labels_re, forest_mid_re, right_table_re, nrow = 1,
                         rel_widths = c(rel_w_left, rel_w_mid, rel_w_right), align = "h")

# ------------------ Bottom info row (tau^2 and I^2) --------------------------------
tau2_txt <- ifelse(is.na(tau2_reml), "NA", formatC(tau2_reml, digits = 4, format = "f"))
I2_txt   <- ifelse(is.na(I2_pct), "NA", formatC(I2_pct, digits = 1, format = "f"))
k_txt    <- as.character(k_used)

info_text <- paste0("RE (REML):  τ² = ", tau2_txt, "   |   I² = ", I2_txt)

info_grob <- ggdraw() +
  draw_label(info_text, x = 0.5, y = 0.5, hjust = 0.5, size = base_text_size + 2)

final_re_plot <- plot_grid(combined_re, info_grob, ncol = 1, rel_heights = c(0.95, 0.05))

# ------------------ Print the final plot ---------------------------------------
print(final_re_plot)
