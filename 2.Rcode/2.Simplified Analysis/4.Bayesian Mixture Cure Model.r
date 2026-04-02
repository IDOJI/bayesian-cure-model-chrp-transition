# Configure Paths ---------------------------------------------------------
find_repo_root <- function(start_dir) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)

  repeat {
    has_repo_markers <- dir.exists(file.path(current_dir, ".git")) ||
      (
        dir.exists(file.path(current_dir, "0.Data")) &&
          dir.exists(file.path(current_dir, "2.Rcode")) &&
          dir.exists(file.path(current_dir, "3.Results files"))
      )

    if (has_repo_markers) {
      return(current_dir)
    }

    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      break
    }
    current_dir <- parent_dir
  }

  stop("Could not locate the repository root.", call. = FALSE)
}

command_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", command_args, value = TRUE)
search_start_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

repo_root <- find_repo_root(search_start_dir)

source_data_file <- Sys.getenv(
  "STAGE8A_SIMPLE_DATA_FILE",
  unset = Sys.getenv(
    "STAGE8A_SIMPLE_MERGED_DATA_PATH",
    unset = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/0.Data/2.Preprocessed data/Preprocessed_Merged_PNUH_SNUH_Data.csv"
  )
)
export_path <- Sys.getenv(
  "STAGE8A_SIMPLE_EXPORT_PATH",
  unset = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/1.Modeling/4.Bayesian Mixture Cure"
)
pnu_site_label <- Sys.getenv("STAGE8A_SIMPLE_PNU_SITE_LABEL", unset = "PNU")
snu_site_label <- Sys.getenv("STAGE8A_SIMPLE_SNU_SITE_LABEL", unset = "SNU")

horizon_years <- 1:10
curve_horizon_max_year <- 10
curve_step_year <- 0.05
detail_tick_step_year <- 0.5
time_origin_epsilon_year <- 1e-08
main_risk_scale <- "transition_only_main"
output_prefix <- "bayesian_mixture_cure"

stan_chains <- as.integer(Sys.getenv("STAGE8A_SIMPLE_STAN_CHAINS", unset = "4"))
stan_iter <- as.integer(Sys.getenv("STAGE8A_SIMPLE_STAN_ITER", unset = "2000"))
stan_warmup <- as.integer(Sys.getenv("STAGE8A_SIMPLE_STAN_WARMUP", unset = "1000"))
stan_thin <- as.integer(Sys.getenv("STAGE8A_SIMPLE_STAN_THIN", unset = "1"))
stan_seed <- as.integer(Sys.getenv("STAGE8A_SIMPLE_STAN_SEED", unset = "20260329"))
stan_adapt_delta <- as.numeric(Sys.getenv("STAGE8A_SIMPLE_STAN_ADAPT_DELTA", unset = "0.95"))
stan_max_treedepth <- as.integer(Sys.getenv("STAGE8A_SIMPLE_STAN_MAX_TREEDEPTH", unset = "12"))
stan_refresh <- as.integer(Sys.getenv("STAGE8A_SIMPLE_STAN_REFRESH", unset = "0"))
posterior_prediction_draws <- as.integer(Sys.getenv("STAGE8A_SIMPLE_POSTERIOR_PRED_DRAWS", unset = "400"))
reuse_existing_fit_rds <- identical(
  toupper(trimws(Sys.getenv("STAGE8A_SIMPLE_REUSE_EXISTING_FIT_RDS", unset = "TRUE"))),
  "TRUE"
)

plot_width_in <- 10
plot_height_in <- 6
plot_dpi <- 320
detail_plot_width_in <- 14
detail_plot_height_in <- 8

horizon_summary_file <- file.path(export_path, paste0(output_prefix, "_horizon_summary.csv"))
plot_source_file <- file.path(export_path, paste0(output_prefix, "_plot_source.csv"))
detail_annotation_file <- file.path(export_path, paste0(output_prefix, "_detail_annotation_table.csv"))
model_registry_file <- file.path(export_path, paste0(output_prefix, "_model_registry.csv"))
parameter_summary_file <- file.path(export_path, paste0(output_prefix, "_parameter_summary.csv"))
plot_rds_file <- file.path(export_path, paste0(output_prefix, "_plot_objects.rds"))
trace_plot_pdf_file <- file.path(export_path, paste0(output_prefix, "_trace_plots.pdf"))
fit_rds_path <- function(model_id, export_dir = export_path) {
  file.path(export_dir, paste0(output_prefix, "__", model_id, "__fit.rds"))
}

dataset_registry <- tibble::tibble(
  dataset = c("PNU", "SNU", "merged_no_site", "merged_site_adjusted"),
  dataset_label = c("PNU", "SNU", "Merged", "Merged (site-adjusted)")
)
dataset_order <- dataset_registry$dataset

dataset_palette <- c(
  "PNU" = "#1B4332",
  "SNU" = "#2A6F97",
  "merged_no_site" = "#C1666B",
  "merged_site_adjusted" = "#B8860B"
)

at_risk_bar_fill <- "#8FA3BF"
transition_bar_fill <- "#D98C6C"

if (anyNA(c(
  stan_chains,
  stan_iter,
  stan_warmup,
  stan_thin,
  stan_seed,
  stan_adapt_delta,
  stan_max_treedepth,
  posterior_prediction_draws
))) {
  stop("Stan and posterior-prediction control values must all be finite.", call. = FALSE)
}
if (stan_iter <= stan_warmup) {
  stop("`stan_iter` must be larger than `stan_warmup`.", call. = FALSE)
}
if (posterior_prediction_draws < 1L) {
  stop("`posterior_prediction_draws` must be at least 1.", call. = FALSE)
}

# Load Packages -----------------------------------------------------------
required_packages <- c("readr", "dplyr", "tibble", "ggplot2", "scales", "rstan")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install required packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(rstan)
})

options(stringsAsFactors = FALSE, scipen = 999)
rstan_options(auto_write = TRUE)
detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}
options(mc.cores = max(1L, min(stan_chains, detected_cores)))
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# Define Helpers ----------------------------------------------------------
assert_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

clip_prob <- function(x, eps = 1e-12) {
  pmin(pmax(x, eps), 1 - eps)
}

read_delimited_or_rds <- function(path) {
  assert_file_exists(path, "Input file")
  ext <- tolower(tools::file_ext(path))

  if (grepl("\\.csv\\.gz$", path, ignore.case = TRUE)) {
    return(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  }
  if (ext == "csv") {
    return(readr::read_csv(
      path,
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE,
      progress = FALSE
    ))
  }
  if (ext == "rds") {
    return(readRDS(path))
  }

  stop("Unsupported extension for file: ", path, call. = FALSE)
}

ensure_draw_matrix <- function(x, ncol_expected) {
  if (is.null(dim(x))) {
    return(matrix(as.numeric(x), ncol = ncol_expected))
  }

  if (length(dim(x)) != 2L) {
    stop("Posterior draws must be a vector or a 2D matrix.", call. = FALSE)
  }

  matrix(
    as.numeric(x),
    nrow = dim(x)[1L],
    ncol = ncol_expected,
    dimnames = dimnames(x)
  )
}

format_half_year_tick_labels <- function(x) {
  ifelse(
    abs(x - round(x)) < 1e-9,
    formatC(x, format = "f", digits = 0),
    formatC(x, format = "f", digits = 1)
  )
}

# Define Data Preparation -------------------------------------------------
standardize_known_site_labels <- function(df, pnu_label, snu_label) {
  if (!("site" %in% names(df))) {
    stop("Input data must contain a `site` column.", call. = FALSE)
  }

  df %>%
    mutate(
      site = trimws(as.character(site)),
      site = case_when(
        toupper(site) == toupper(pnu_label) ~ pnu_label,
        toupper(site) == toupper(snu_label) ~ snu_label,
        TRUE ~ site
      )
    )
}

prepare_analysis_dataset <- function(df, dataset_name, pnu_label, snu_label, site_mode = c("single", "merged")) {
  site_mode <- match.arg(site_mode)

  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("[%s] Missing required columns: %s", dataset_name, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  out <- tibble::as_tibble(df) %>%
    standardize_known_site_labels(pnu_label = pnu_label, snu_label = snu_label) %>%
    mutate(
      id = trimws(as.character(id)),
      sex_num = as.integer(safe_numeric(sex_num)),
      age_exact_entry = safe_numeric(age_exact_entry),
      days_followup = safe_numeric(days_followup),
      status_num = as.integer(safe_numeric(status_num))
    )

  if (nrow(out) == 0L) {
    stop(sprintf("[%s] Dataset has zero rows.", dataset_name), call. = FALSE)
  }
  if (anyNA(out[, required_cols])) {
    stop(sprintf("[%s] Missing values detected in required columns.", dataset_name), call. = FALSE)
  }
  if (any(out$id == "", na.rm = TRUE)) {
    stop(sprintf("[%s] Blank `id` values detected.", dataset_name), call. = FALSE)
  }
  if (any(out$site == "", na.rm = TRUE)) {
    stop(sprintf("[%s] Blank `site` values detected.", dataset_name), call. = FALSE)
  }
  if (any(!out$sex_num %in% c(0L, 1L), na.rm = TRUE)) {
    stop(sprintf("[%s] `sex_num` must be coded 0/1.", dataset_name), call. = FALSE)
  }
  if (any(!out$status_num %in% c(0L, 1L, 2L), na.rm = TRUE)) {
    stop(sprintf("[%s] `status_num` must be coded 0/1/2.", dataset_name), call. = FALSE)
  }
  if (any(out$days_followup < 0, na.rm = TRUE)) {
    stop(sprintf("[%s] Negative `days_followup` values detected.", dataset_name), call. = FALSE)
  }

  n_site_levels <- dplyr::n_distinct(out$site)
  if (site_mode == "single" && n_site_levels != 1L) {
    stop(sprintf("[%s] Single-cohort input must contain exactly one site.", dataset_name), call. = FALSE)
  }
  if (site_mode == "merged" && n_site_levels < 2L) {
    stop(sprintf("[%s] Merged input must contain at least two site levels.", dataset_name), call. = FALSE)
  }

  age_mean <- mean(out$age_exact_entry)
  age_sd <- stats::sd(out$age_exact_entry)
  if (is.na(age_sd) || age_sd <= 0) {
    stop(sprintf("[%s] `age_exact_entry` must have positive SD.", dataset_name), call. = FALSE)
  }

  out <- out %>%
    mutate(
      unique_person_id = paste(site, id, sep = "_"),
      time_year = pmax(days_followup / 365.25, time_origin_epsilon_year),
      age_s = (age_exact_entry - age_mean) / (2 * age_sd),
      event_main = as.integer(status_num == 1L),
      censor_main = as.integer(status_num %in% c(0L, 2L))
    )

  if (dplyr::n_distinct(out$unique_person_id) != nrow(out)) {
    stop(sprintf("[%s] `site + id` is not unique.", dataset_name), call. = FALSE)
  }

  out
}

build_analysis_datasets_from_source <- function(source_file, pnu_label, snu_label) {
  source_df <- read_delimited_or_rds(source_file)
  if (!inherits(source_df, "data.frame")) {
    stop("Merged preprocessed dataset could not be loaded as a data frame.", call. = FALSE)
  }

  source_df <- tibble::as_tibble(source_df) %>%
    standardize_known_site_labels(pnu_label = pnu_label, snu_label = snu_label)

  pnu_df <- source_df %>% filter(site == pnu_label)
  snu_df <- source_df %>% filter(site == snu_label)

  if (nrow(pnu_df) == 0L) {
    stop("No PNU rows were found after standardizing `site` labels.", call. = FALSE)
  }
  if (nrow(snu_df) == 0L) {
    stop("No SNU rows were found after standardizing `site` labels.", call. = FALSE)
  }

  list(
    PNU = prepare_analysis_dataset(
      df = pnu_df,
      dataset_name = "PNU",
      pnu_label = pnu_label,
      snu_label = snu_label,
      site_mode = "single"
    ),
    SNU = prepare_analysis_dataset(
      df = snu_df,
      dataset_name = "SNU",
      pnu_label = pnu_label,
      snu_label = snu_label,
      site_mode = "single"
    ),
    merged = prepare_analysis_dataset(
      df = source_df,
      dataset_name = "merged",
      pnu_label = pnu_label,
      snu_label = snu_label,
      site_mode = "merged"
    )
  )
}

summarize_analysis_datasets <- function(analysis_datasets) {
  bind_rows(lapply(names(analysis_datasets), function(dataset_name) {
    df <- analysis_datasets[[dataset_name]]
    tibble(
      source_dataset = dataset_name,
      n_subjects = nrow(df),
      n_transition = sum(df$event_main == 1L, na.rm = TRUE),
      n_censored = sum(df$censor_main == 1L, na.rm = TRUE),
      max_observed_followup_year = max(as.numeric(df$time_year), na.rm = TRUE)
    )
  }))
}

build_detail_annotation_table <- function(analysis_df, times) {
  tibble(
    time_horizon_year = as.numeric(times),
    n_at_risk = vapply(
      times,
      function(tt) sum(as.numeric(analysis_df$time_year) >= tt, na.rm = TRUE),
      integer(1)
    ),
    n_transition_cumulative = vapply(
      times,
      function(tt) sum(as.integer(analysis_df$event_main) == 1L & as.numeric(analysis_df$time_year) <= tt, na.rm = TRUE),
      integer(1)
    )
  )
}

# Define Priors and Model Specifications ---------------------------------
build_prior_artifacts <- function() {
  list(
    alpha_gp_vdw = -9.581369553169,
    mu_beta_inc_anchor = c(
      0.419871845822,
      0.907608052926,
      0.586202561451,
      0.466865123863,
      0.037997248763
    ),
    sd_beta_inc_anchor = c(
      0.132789397422,
      0.173731076538,
      0.191221553945,
      0.270393197518,
      0.302838606651
    ),
    alpha_inc_sd = 3.5,
    sd_gamma0 = 2.5,
    sd_gamma_lat = 1.0,
    sd_log_sigma = 0.50
  )
}

build_model_spec_registry <- function() {
  tibble(
    dataset = c("PNU", "SNU", "merged_no_site", "merged_site_adjusted"),
    dataset_label = c("PNU", "SNU", "Merged", "Merged (site-adjusted)"),
    source_dataset_key = c("PNU", "SNU", "merged", "merged"),
    model_id = c(
      "PNU-L0-LN-ANCH-NOSITE",
      "SNU-L0-LN-ANCH-NOSITE",
      "MERGED-I0-L0S0-LN-ANCH-NOSITE",
      "MERGED-I1-L0S1-LN-ANCH-N01"
    ),
    structural_model_id = c("PNU-L0", "SNU-L0", "MERGED-I0-L0S0", "MERGED-I1-L0S1"),
    site_adjustment_flag = c(FALSE, FALSE, FALSE, TRUE),
    incidence_site_indicator = c(FALSE, FALSE, FALSE, TRUE),
    latency_site_indicator = c(FALSE, FALSE, FALSE, TRUE),
    formula_anchor = c("base", "base", "base", "site_added")
  )
}

make_design_bundle <- function(df, model_row, prior_artifacts, snu_label) {
  z_i <- as.integer(df$sex_num)
  x20_i <- as.integer(df$age_exact_entry >= 20 & df$age_exact_entry < 30)
  x30_i <- as.integer(df$age_exact_entry >= 30)
  s_i <- as.integer(df$site == snu_label)

  X_inc <- cbind(
    sex_num = z_i,
    age20_29 = x20_i,
    age30plus = x30_i,
    sex_x_age20_29 = z_i * x20_i,
    sex_x_age30plus = z_i * x30_i
  )
  inc_is_site <- rep(0L, ncol(X_inc))

  mu_beta_inc <- prior_artifacts$mu_beta_inc_anchor
  sd_beta_inc <- prior_artifacts$sd_beta_inc_anchor

  if (isTRUE(model_row$incidence_site_indicator[[1L]])) {
    X_inc <- cbind(X_inc, site_SNU = s_i)
    mu_beta_inc <- c(mu_beta_inc, 0)
    sd_beta_inc <- c(sd_beta_inc, 1)
    inc_is_site <- c(inc_is_site, 1L)
  }

  a_i <- as.numeric(df$age_s)

  X_lat <- switch(
    model_row$structural_model_id[[1L]],
    "PNU-L0" = cbind(age_s = a_i, sex_num = z_i),
    "SNU-L0" = cbind(age_s = a_i, sex_num = z_i),
    "MERGED-I0-L0S0" = cbind(age_s = a_i, sex_num = z_i),
    "MERGED-I1-L0S1" = cbind(age_s = a_i, sex_num = z_i, site_SNU = s_i),
    stop("Unknown structural model id: ", model_row$structural_model_id[[1L]], call. = FALSE)
  )

  list(
    X_inc = unclass(as.matrix(X_inc)),
    inc_is_site = as.integer(inc_is_site),
    mu_beta_inc = as.numeric(mu_beta_inc),
    sd_beta_inc = as.numeric(sd_beta_inc),
    alpha_inc_prior_mean = as.numeric(prior_artifacts$alpha_gp_vdw),
    alpha_inc_prior_sd = as.numeric(prior_artifacts$alpha_inc_sd),
    X_lat = unclass(as.matrix(X_lat)),
    lat_is_site = as.integer(colnames(X_lat) == "site_SNU"),
    sd_gamma_lat = rep(prior_artifacts$sd_gamma_lat, ncol(X_lat)),
    sd_gamma0 = as.numeric(prior_artifacts$sd_gamma0),
    sd_log_sigma = as.numeric(prior_artifacts$sd_log_sigma),
    time = as.numeric(df$time_year),
    event = as.integer(df$event_main)
  )
}

# Define Stan Model -------------------------------------------------------
compile_bayesian_mixture_cure_stan_model <- function() {
  stan_code <- r"(
data {
  int<lower=1> N;
  vector<lower=1e-8>[N] time;
  int<lower=0, upper=1> event[N];
  int<lower=1> K_inc;
  matrix[N, K_inc] X_inc;
  array[K_inc] int<lower=0, upper=1> inc_is_site;
  int<lower=1> K_lat;
  matrix[N, K_lat] X_lat;
  array[K_lat] int<lower=0, upper=1> lat_is_site;
  real alpha_inc_prior_mean;
  real<lower=0> alpha_inc_prior_sd;
  vector[K_inc] mu_beta_inc;
  vector<lower=0>[K_inc] sd_beta_inc;
  real<lower=0> sd_gamma0;
  vector<lower=0>[K_lat] sd_gamma_lat;
  real<lower=0> sd_log_sigma;
}
parameters {
  real alpha_inc;
  vector[K_inc] beta_inc;
  real gamma0;
  vector[K_lat] gamma_lat;
  real log_sigma;
}
model {
  real sigma;
  sigma = exp(log_sigma);

  alpha_inc ~ normal(alpha_inc_prior_mean, alpha_inc_prior_sd);

  for (j in 1:K_inc) {
    if (inc_is_site[j] == 1) {
      beta_inc[j] ~ normal(0, 1);
    } else {
      beta_inc[j] ~ normal(mu_beta_inc[j], sd_beta_inc[j]);
    }
  }

  gamma0 ~ normal(0, sd_gamma0);
  for (j in 1:K_lat) {
    if (lat_is_site[j] == 1) {
      gamma_lat[j] ~ normal(0, 1);
    } else {
      gamma_lat[j] ~ normal(0, sd_gamma_lat[j]);
    }
  }

  log_sigma ~ normal(0, sd_log_sigma);

  for (i in 1:N) {
    real eta_inc;
    real pi_i;
    real mu_lat;
    real logS;
    real logf;

    eta_inc = alpha_inc + dot_product(row(X_inc, i), beta_inc);
    pi_i = inv_logit(eta_inc);
    mu_lat = gamma0 + dot_product(row(X_lat, i), gamma_lat);
    logS = normal_lccdf(log(time[i]) | mu_lat, sigma);
    logf = lognormal_lpdf(time[i] | mu_lat, sigma);

    if (event[i] == 1) {
      target += log(pi_i) + logf;
    } else {
      target += log_sum_exp(log1m(pi_i), log(pi_i) + logS);
    }
  }
}
)"

  rstan::stan_model(
    model_code = stan_code,
    model_name = "bayesian_lognormal_mixture_cure"
  )
}

fit_bayesian_lognormal_mixture_cure <- function(stan_model_compiled, design_bundle, seed_value) {
  stan_data <- list(
    N = nrow(design_bundle$X_inc),
    time = as.numeric(design_bundle$time),
    event = as.integer(design_bundle$event),
    K_inc = ncol(design_bundle$X_inc),
    X_inc = design_bundle$X_inc,
    inc_is_site = as.integer(design_bundle$inc_is_site),
    K_lat = ncol(design_bundle$X_lat),
    X_lat = design_bundle$X_lat,
    lat_is_site = as.integer(design_bundle$lat_is_site),
    alpha_inc_prior_mean = as.numeric(design_bundle$alpha_inc_prior_mean),
    alpha_inc_prior_sd = as.numeric(design_bundle$alpha_inc_prior_sd),
    mu_beta_inc = as.numeric(design_bundle$mu_beta_inc),
    sd_beta_inc = as.numeric(design_bundle$sd_beta_inc),
    sd_gamma0 = as.numeric(design_bundle$sd_gamma0),
    sd_gamma_lat = as.numeric(design_bundle$sd_gamma_lat),
    sd_log_sigma = as.numeric(design_bundle$sd_log_sigma)
  )

  rstan::sampling(
    object = stan_model_compiled,
    data = stan_data,
    chains = stan_chains,
    iter = stan_iter,
    warmup = stan_warmup,
    thin = stan_thin,
    seed = seed_value,
    refresh = stan_refresh,
    control = list(
      adapt_delta = stan_adapt_delta,
      max_treedepth = stan_max_treedepth
    )
  )
}

extract_prediction_draws <- function(fit, K_inc, K_lat, max_draws, seed_value) {
  ext <- rstan::extract(
    fit,
    pars = c("alpha_inc", "beta_inc", "gamma0", "gamma_lat", "log_sigma"),
    permuted = TRUE,
    inc_warmup = FALSE
  )

  total_draws <- length(ext$alpha_inc)
  set.seed(seed_value)
  keep_draw_idx <- if (total_draws <= max_draws) {
    seq_len(total_draws)
  } else {
    sort(sample(seq_len(total_draws), size = max_draws, replace = FALSE))
  }

  list(
    alpha_inc = as.numeric(ext$alpha_inc[keep_draw_idx]),
    beta_inc = ensure_draw_matrix(ext$beta_inc, ncol_expected = K_inc)[keep_draw_idx, , drop = FALSE],
    gamma0 = as.numeric(ext$gamma0[keep_draw_idx]),
    gamma_lat = ensure_draw_matrix(ext$gamma_lat, ncol_expected = K_lat)[keep_draw_idx, , drop = FALSE],
    log_sigma = as.numeric(ext$log_sigma[keep_draw_idx]),
    n_draws = length(keep_draw_idx)
  )
}

# Define Parameter Diagnostics and Trace Plots ----------------------------
build_parameter_names <- function(K_inc, K_lat) {
  c(
    "alpha_inc",
    paste0("beta_inc[", seq_len(K_inc), "]"),
    "gamma0",
    paste0("gamma_lat[", seq_len(K_lat), "]"),
    "log_sigma"
  )
}

extract_parameter_summary <- function(fit, model_row, K_inc, K_lat) {
  parameter_names <- build_parameter_names(K_inc = K_inc, K_lat = K_lat)
  fit_summary <- rstan::summary(fit, pars = parameter_names)$summary

  tibble(
    dataset = model_row$dataset[[1L]],
    dataset_label = model_row$dataset_label[[1L]],
    model_id = model_row$model_id[[1L]],
    site_adjustment_flag = as.logical(model_row$site_adjustment_flag[[1L]]),
    parameter = rownames(fit_summary),
    mean = fit_summary[, "mean"],
    sd = fit_summary[, "sd"],
    q025 = fit_summary[, "2.5%"],
    q50 = fit_summary[, "50%"],
    q975 = fit_summary[, "97.5%"],
    n_eff = fit_summary[, "n_eff"],
    rhat = fit_summary[, "Rhat"]
  )
}

make_model_registry_row <- function(
  model_row,
  analysis_df,
  design_bundle,
  elapsed_seconds,
  divergent_transitions,
  treedepth_hits,
  prediction_draws,
  parameter_summary_df
) {
  tibble(
    dataset = model_row$dataset[[1L]],
    dataset_label = model_row$dataset_label[[1L]],
    source_dataset_key = model_row$source_dataset_key[[1L]],
    model_id = model_row$model_id[[1L]],
    structural_model_id = model_row$structural_model_id[[1L]],
    formula_anchor = model_row$formula_anchor[[1L]],
    site_adjustment_flag = as.logical(model_row$site_adjustment_flag[[1L]]),
    incidence_site_indicator = as.logical(model_row$incidence_site_indicator[[1L]]),
    latency_site_indicator = as.logical(model_row$latency_site_indicator[[1L]]),
    risk_scale = main_risk_scale,
    fit_status = "ok",
    n_subjects = nrow(analysis_df),
    n_transition = sum(analysis_df$event_main == 1L, na.rm = TRUE),
    n_censored = sum(analysis_df$censor_main == 1L, na.rm = TRUE),
    K_inc = ncol(design_bundle$X_inc),
    K_lat = ncol(design_bundle$X_lat),
    stan_chains = stan_chains,
    stan_iter = stan_iter,
    stan_warmup = stan_warmup,
    stan_seed = NA_integer_,
    posterior_prediction_draws = as.integer(prediction_draws),
    elapsed_seconds = as.numeric(elapsed_seconds),
    divergent_transitions = as.integer(divergent_transitions),
    max_treedepth_hits = as.integer(treedepth_hits),
    max_rhat = max(parameter_summary_df$rhat, na.rm = TRUE),
    min_n_eff = min(parameter_summary_df$n_eff, na.rm = TRUE)
  )
}

write_trace_plot_pdf <- function(fit_records, output_file) {
  if (length(fit_records) == 0L) {
    return(invisible(FALSE))
  }

  grDevices::pdf(output_file, width = 11, height = 8.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  chain_cols <- c("#1B4332", "#2A6F97", "#C1666B", "#B8860B", "#6A4C93", "#7F5539")

  for (record in fit_records) {
    fit <- record$fit
    parameter_names <- build_parameter_names(K_inc = record$K_inc, K_lat = record$K_lat)
    arr <- as.array(fit, pars = parameter_names)
    n_chains <- dim(arr)[2L]

    chunk_starts <- seq.int(1L, length(parameter_names), by = 4L)
    for (chunk_start in chunk_starts) {
      chunk_idx <- chunk_start:min(chunk_start + 3L, length(parameter_names))
      par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))

      for (jj in chunk_idx) {
        param_name <- parameter_names[[jj]]
        trace_mat <- arr[, , param_name, drop = TRUE]
        if (is.null(dim(trace_mat))) {
          trace_mat <- matrix(trace_mat, ncol = 1L)
        }

        y_limits <- range(trace_mat, finite = TRUE)
        x_index <- seq_len(nrow(trace_mat))

        plot(
          x_index,
          trace_mat[, 1L],
          type = "l",
          col = chain_cols[[1L]],
          lwd = 1,
          xlab = "Post-warmup iteration",
          ylab = "Value",
          main = param_name,
          ylim = y_limits
        )
        if (n_chains > 1L) {
          for (cc in 2:n_chains) {
            lines(x_index, trace_mat[, cc], col = chain_cols[[cc]], lwd = 1)
          }
        }
        legend(
          "topright",
          legend = paste("Chain", seq_len(n_chains)),
          col = chain_cols[seq_len(n_chains)],
          lty = 1,
          lwd = 1,
          bty = "n",
          cex = 0.75
        )
      }

      if (length(chunk_idx) < 4L) {
        for (kk in seq_len(4L - length(chunk_idx))) {
          plot.new()
        }
      }

      mtext(
        paste0("Trace plots: ", record$model_id, " (", record$dataset, ")"),
        outer = TRUE,
        cex = 1.1,
        font = 2
      )
    }
  }

  invisible(TRUE)
}

# Define Prediction Helpers ----------------------------------------------
build_curve_prediction_df_from_draws <- function(model_row, design_bundle, draw_bundle, times) {
  X_inc <- design_bundle$X_inc
  X_lat <- design_bundle$X_lat
  n_subjects <- nrow(X_inc)

  eta_inc_mat <- X_inc %*% t(draw_bundle$beta_inc)
  eta_inc_mat <- sweep(eta_inc_mat, 2, draw_bundle$alpha_inc, "+")
  pi_mat <- clip_prob(plogis(eta_inc_mat))

  mu_lat_mat <- X_lat %*% t(draw_bundle$gamma_lat)
  mu_lat_mat <- sweep(mu_lat_mat, 2, draw_bundle$gamma0, "+")
  sigma_vec <- exp(draw_bundle$log_sigma)

  subject_susceptible_fraction <- rowMeans(pi_mat)

  bind_rows(lapply(times, function(tt) {
    if (tt <= 0) {
      subject_susceptible_only_survival <- rep(1, n_subjects)
      subject_overall_survival <- rep(1, n_subjects)
    } else {
      z_mat <- sweep(log(tt) - mu_lat_mat, 2, sigma_vec, "/")
      susceptible_only_survival_draws <- 1 - stats::pnorm(z_mat)
      overall_survival_draws <- (1 - pi_mat) + pi_mat * susceptible_only_survival_draws

      subject_susceptible_only_survival <- rowMeans(susceptible_only_survival_draws)
      subject_overall_survival <- rowMeans(overall_survival_draws)
    }

    tibble(
      dataset = model_row$dataset[[1L]],
      dataset_label = model_row$dataset_label[[1L]],
      plot_dataset_label = model_row$dataset_label[[1L]],
      model_id = model_row$model_id[[1L]],
      site_adjustment_flag = as.logical(model_row$site_adjustment_flag[[1L]]),
      risk_scale = main_risk_scale,
      time_horizon_year = as.numeric(tt),
      susceptible_fraction = mean(subject_susceptible_fraction),
      cure_fraction = mean(1 - subject_susceptible_fraction),
      overall_survival_prob = mean(subject_overall_survival),
      overall_risk_prob = mean(1 - subject_overall_survival),
      susceptible_only_survival_prob = mean(subject_susceptible_only_survival),
      susceptible_only_risk_prob = mean(1 - subject_susceptible_only_survival),
      n_subjects = as.integer(n_subjects)
    )
  }))
}

# Build Plot Helpers ------------------------------------------------------
build_plot_group_registry <- function() {
  tibble(
    plot_group = c("all_cohorts", "pnu_only", "snu_only", "merged_only"),
    group_title = c(
      "PNU, SNU, and merged cohort models",
      "PNU cohort",
      "SNU cohort",
      "Merged cohort models"
    ),
    dataset_keys = list(
      c("PNU", "SNU", "merged_no_site", "merged_site_adjusted"),
      "PNU",
      "SNU",
      c("merged_no_site", "merged_site_adjusted")
    ),
    include_pnu_reference = c(TRUE, TRUE, FALSE, FALSE)
  )
}

build_dataset_detail_plot_registry <- function() {
  tibble(
    plot_group = c("pnu_detail", "snu_detail", "merged_detail"),
    group_title = c(
      "PNU cohort",
      "SNU cohort",
      "Merged cohort models"
    ),
    dataset_keys = list(
      "PNU",
      "SNU",
      c("merged_no_site", "merged_site_adjusted")
    ),
    source_dataset = c("PNU", "SNU", "merged"),
    include_pnu_reference = c(TRUE, FALSE, FALSE)
  )
}

make_plot_output_file <- function(export_dir, metric_name, plot_group) {
  file.path(export_dir, paste0(output_prefix, "_", metric_name, "__", plot_group, ".png"))
}

remove_existing_detail_plot_outputs <- function(export_dir, annotation_file) {
  if (!dir.exists(export_dir)) {
    return(invisible(NULL))
  }

  existing_detail_pngs <- list.files(
    export_dir,
    pattern = paste0("^", output_prefix, "_.*with_counts__.*\\.png$"),
    full.names = TRUE
  )
  existing_targets <- unique(c(existing_detail_pngs, annotation_file[file.exists(annotation_file)]))

  if (length(existing_targets) > 0L) {
    unlink(existing_targets)
  }

  invisible(existing_targets)
}

make_curve_plot <- function(
  curve_df,
  value_col,
  y_label,
  title_text,
  pnu_reference_year = NA_real_,
  x_breaks = horizon_years,
  x_labels = scales::label_number(accuracy = 1),
  y_lower_limit = 0,
  caption_text = NULL
) {
  dataset_labels <- curve_df %>%
    distinct(dataset, plot_dataset_label) %>%
    arrange(match(dataset, dataset_order))

  yearly_points <- curve_df %>%
    filter(time_horizon_year %in% horizon_years)

  show_legend <- nrow(dataset_labels) > 1L
  subtitle_text <- "Cohort-level posterior means from direct Bayesian log-normal mixture cure model fits"

  if (is.finite(pnu_reference_year) && any(curve_df$dataset == "PNU")) {
    subtitle_text <- paste0(
      subtitle_text,
      sprintf(" | Dashed line = PNU max observed follow-up (%.2f years)", pnu_reference_year)
    )
  }

  plot_object <- ggplot(curve_df, aes(x = time_horizon_year, y = .data[[value_col]], color = dataset)) +
    geom_line(linewidth = 1.05) +
    geom_point(
      data = yearly_points,
      aes(x = time_horizon_year, y = .data[[value_col]], color = dataset),
      size = 1.8
    ) +
    scale_color_manual(
      values = dataset_palette,
      breaks = dataset_labels$dataset,
      labels = dataset_labels$plot_dataset_label
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      limits = c(0, curve_horizon_max_year),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      breaks = seq(0, 1, by = 0.2),
      limits = c(y_lower_limit, 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Years after cohort entry (k)",
      y = y_label,
      color = "Dataset",
      caption = caption_text
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = if (show_legend) "top" else "none",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 9, hjust = 0),
      plot.caption.position = "plot",
      plot.margin = margin(
        t = 5.5,
        r = 5.5,
        b = if (y_lower_limit < 0) 28 else 5.5,
        l = 5.5
      ),
      axis.text.x = element_text(size = if (length(x_breaks) > length(horizon_years)) 7 else 9)
    ) +
    coord_cartesian(clip = "off")

  if (is.finite(pnu_reference_year) && any(curve_df$dataset == "PNU")) {
    plot_object <- plot_object +
      geom_vline(
        xintercept = pnu_reference_year,
        linetype = "dashed",
        linewidth = 0.7,
        color = dataset_palette[["PNU"]]
      )
  }

  plot_object
}

make_detail_curve_plot <- function(
  curve_df,
  value_col,
  y_label,
  title_text,
  analysis_df,
  pnu_reference_year = NA_real_
) {
  detail_times <- seq(0, curve_horizon_max_year, by = detail_tick_step_year)
  detail_bar_x_limits <- c(
    -detail_tick_step_year * 0.4,
    curve_horizon_max_year + detail_tick_step_year * 0.4
  )

  curve_plot <- make_curve_plot(
    curve_df = curve_df,
    value_col = value_col,
    y_label = y_label,
    title_text = title_text,
    pnu_reference_year = pnu_reference_year,
    x_breaks = detail_times,
    x_labels = format_half_year_tick_labels,
    y_lower_limit = 0,
    caption_text = NULL
  ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5)
    )

  annotation_tbl <- build_detail_annotation_table(analysis_df, detail_times)

  at_risk_plot <- ggplot(annotation_tbl, aes(x = time_horizon_year, y = n_at_risk)) +
    geom_col(width = detail_tick_step_year * 0.72, fill = at_risk_bar_fill) +
    geom_text(aes(label = n_at_risk), vjust = -0.25, size = 2.8) +
    scale_x_continuous(
      breaks = detail_times,
      labels = format_half_year_tick_labels,
      limits = detail_bar_x_limits,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(x = NULL, y = "n at risk") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 5.5, b = 0, l = 5.5)
    )

  transition_plot <- ggplot(annotation_tbl, aes(x = time_horizon_year, y = n_transition_cumulative)) +
    geom_col(width = detail_tick_step_year * 0.72, fill = transition_bar_fill) +
    geom_text(aes(label = n_transition_cumulative), vjust = -0.25, size = 2.8) +
    scale_x_continuous(
      breaks = detail_times,
      labels = format_half_year_tick_labels,
      limits = detail_bar_x_limits,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(
      x = "Years after cohort entry (k)",
      y = "cumulative\ntransitions"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(size = 7),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5)
    )

  list(
    curve_plot = curve_plot,
    at_risk_plot = at_risk_plot,
    transition_plot = transition_plot,
    annotation_tbl = annotation_tbl
  )
}

save_plot_png <- function(plot_object, output_file, width = plot_width_in, height = plot_height_in) {
  ggplot2::ggsave(
    filename = output_file,
    plot = plot_object,
    width = width,
    height = height,
    dpi = plot_dpi,
    units = "in"
  )
}

save_stacked_plot_png <- function(plot_bundle, output_file, width = detail_plot_width_in, height = detail_plot_height_in) {
  grDevices::png(filename = output_file, width = width, height = height, units = "in", res = plot_dpi)
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()
  grid::pushViewport(
    grid::viewport(
      layout = grid::grid.layout(
        nrow = 3,
        ncol = 1,
        heights = grid::unit(c(3.8, 1.4, 1.4), "null")
      )
    )
  )

  print(plot_bundle$curve_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(plot_bundle$at_risk_plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(plot_bundle$transition_plot, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
}

# Load Data and Compile Model --------------------------------------------
analysis_datasets <- build_analysis_datasets_from_source(
  source_file = source_data_file,
  pnu_label = pnu_site_label,
  snu_label = snu_site_label
)
dataset_summary <- summarize_analysis_datasets(analysis_datasets)
pnu_max_observed_followup_year <- dataset_summary %>%
  filter(source_dataset == "PNU") %>%
  pull(max_observed_followup_year) %>%
  .[[1L]]

if (!is.finite(pnu_max_observed_followup_year)) {
  stop("Could not derive the PNU maximum observed follow-up in years.", call. = FALSE)
}

prior_artifacts <- build_prior_artifacts()
model_spec_registry <- build_model_spec_registry() %>%
  mutate(fit_rds_path = vapply(model_id, fit_rds_path, character(1)))

curve_time_grid <- sort(unique(round(c(
  0,
  seq(0, curve_horizon_max_year, by = curve_step_year),
  horizon_years
), 10)))

n_models_requiring_fit <- sum(!(reuse_existing_fit_rds & file.exists(model_spec_registry$fit_rds_path)))
if (n_models_requiring_fit > 0L) {
  message(
    "Compiling Bayesian log-normal mixture cure model for ",
    n_models_requiring_fit,
    " model(s) that still need fitting."
  )
  stan_model_compiled <- compile_bayesian_mixture_cure_stan_model()
} else {
  message("All fitted-model RDS files already exist. Stan compilation skipped.")
  stan_model_compiled <- NULL
}

# Fit Models and Build Cohort-Level Prediction Data -----------------------
curve_data_list <- vector("list", nrow(model_spec_registry))
parameter_summary_list <- vector("list", nrow(model_spec_registry))
model_registry_list <- vector("list", nrow(model_spec_registry))
fit_records <- vector("list", nrow(model_spec_registry))

for (ii in seq_len(nrow(model_spec_registry))) {
  model_row <- model_spec_registry[ii, , drop = FALSE]
  source_dataset_key <- model_row$source_dataset_key[[1L]]
  analysis_df <- analysis_datasets[[source_dataset_key]]
  seed_value <- stan_seed + ii

  message(
    "[", ii, "/", nrow(model_spec_registry), "] Preparing ",
    model_row$model_id[[1L]],
    " for dataset=",
    model_row$dataset[[1L]]
  )

  design_bundle <- make_design_bundle(
    df = analysis_df,
    model_row = model_row,
    prior_artifacts = prior_artifacts,
    snu_label = snu_site_label
  )

  one_fit_rds_path <- model_row$fit_rds_path[[1L]]
  reused_existing_fit <- FALSE
  fit <- NULL
  elapsed_seconds <- NA_real_

  if (isTRUE(reuse_existing_fit_rds) && file.exists(one_fit_rds_path)) {
    fit_try <- tryCatch(readRDS(one_fit_rds_path), error = function(e) e)
    if (!inherits(fit_try, "error") && inherits(fit_try, "stanfit")) {
      fit <- fit_try
      reused_existing_fit <- TRUE
      message("  reused existing fit RDS: ", basename(one_fit_rds_path))
    } else {
      message("  existing fit RDS could not be reused; refitting.")
    }
  }

  if (is.null(fit)) {
    if (is.null(stan_model_compiled)) {
      stop("Internal error: no compiled Stan model available for a required refit.", call. = FALSE)
    }

    started_at <- Sys.time()
    fit <- fit_bayesian_lognormal_mixture_cure(
      stan_model_compiled = stan_model_compiled,
      design_bundle = design_bundle,
      seed_value = seed_value
    )
    elapsed_seconds <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
    saveRDS(fit, one_fit_rds_path, compress = TRUE)
  }

  draw_bundle <- extract_prediction_draws(
    fit = fit,
    K_inc = ncol(design_bundle$X_inc),
    K_lat = ncol(design_bundle$X_lat),
    max_draws = posterior_prediction_draws,
    seed_value = stan_seed + 100000L + ii
  )

  curve_data_list[[ii]] <- build_curve_prediction_df_from_draws(
    model_row = model_row,
    design_bundle = design_bundle,
    draw_bundle = draw_bundle,
    times = curve_time_grid
  )

  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergent_transitions <- sum(vapply(sampler_params, function(x) sum(x[, "divergent__"] > 0), numeric(1)))
  treedepth_hits <- sum(vapply(sampler_params, function(x) sum(x[, "treedepth__"] >= stan_max_treedepth), numeric(1)))

  parameter_summary_list[[ii]] <- extract_parameter_summary(
    fit = fit,
    model_row = model_row,
    K_inc = ncol(design_bundle$X_inc),
    K_lat = ncol(design_bundle$X_lat)
  )
  model_registry_list[[ii]] <- make_model_registry_row(
    model_row = model_row,
    analysis_df = analysis_df,
    design_bundle = design_bundle,
    elapsed_seconds = elapsed_seconds,
    divergent_transitions = divergent_transitions,
    treedepth_hits = treedepth_hits,
    prediction_draws = draw_bundle$n_draws,
    parameter_summary_df = parameter_summary_list[[ii]]
  ) %>%
    mutate(
      stan_seed = as.integer(seed_value),
      reused_existing_fit_rds = reused_existing_fit,
      fit_rds_path = one_fit_rds_path
    )
  fit_records[[ii]] <- list(
    dataset = model_row$dataset[[1L]],
    model_id = model_row$model_id[[1L]],
    fit = fit,
    K_inc = ncol(design_bundle$X_inc),
    K_lat = ncol(design_bundle$X_lat)
  )

  message(
    "  ",
    if (isTRUE(reused_existing_fit)) {
      "reuse completed"
    } else {
      paste0("completed in ", formatC(elapsed_seconds, digits = 1L, format = "f"), " sec")
    },
    "; divergent=",
    divergent_transitions,
    ", treedepth_hits=",
    treedepth_hits,
    ", prediction_draws=",
    draw_bundle$n_draws
  )
}

curve_data_df <- bind_rows(curve_data_list) %>%
  mutate(time_horizon_year = round(as.numeric(time_horizon_year), 10)) %>%
  arrange(match(dataset, dataset_order), time_horizon_year)

parameter_summary_df <- bind_rows(parameter_summary_list) %>%
  arrange(match(dataset, dataset_order), model_id, parameter)

model_registry_df <- bind_rows(model_registry_list) %>%
  arrange(match(dataset, dataset_order), model_id)

plot_source_df <- curve_data_df %>%
  transmute(
    dataset,
    dataset_label,
    plot_dataset_label,
    model_id,
    site_adjustment_flag,
    risk_scale,
    time_horizon_year,
    susceptible_fraction,
    cure_fraction,
    overall_survival_prob,
    overall_risk_prob,
    susceptible_only_survival_prob,
    susceptible_only_risk_prob,
    n_subjects
  ) %>%
  arrange(match(dataset, dataset_order), time_horizon_year)

horizon_summary_df <- plot_source_df %>%
  filter(time_horizon_year %in% horizon_years) %>%
  mutate(horizon_year = as.integer(time_horizon_year)) %>%
  arrange(match(dataset, dataset_order), horizon_year)

# Build Plot Objects ------------------------------------------------------
plot_group_registry <- build_plot_group_registry()
detail_plot_group_registry <- build_dataset_detail_plot_registry()
plot_registry <- list()
detail_annotation_list <- list()

for (ii in seq_len(nrow(plot_group_registry))) {
  plot_group <- plot_group_registry$plot_group[[ii]]
  group_title <- plot_group_registry$group_title[[ii]]
  group_dataset_keys <- plot_group_registry$dataset_keys[[ii]]
  include_pnu_reference <- isTRUE(plot_group_registry$include_pnu_reference[[ii]])
  group_curve_df <- plot_source_df %>%
    filter(dataset %in% group_dataset_keys)

  if (nrow(group_curve_df) == 0L) {
    stop(sprintf("No plotting rows found for plot group `%s`.", plot_group), call. = FALSE)
  }

  pnu_reference_year <- if (include_pnu_reference) pnu_max_observed_followup_year else NA_real_

  survival_plot <- make_curve_plot(
    curve_df = group_curve_df,
    value_col = "overall_survival_prob",
    y_label = "Estimated survival probability",
    title_text = paste0("Estimated survival probability: ", group_title),
    pnu_reference_year = pnu_reference_year
  )

  risk_plot <- make_curve_plot(
    curve_df = group_curve_df,
    value_col = "overall_risk_prob",
    y_label = "Estimated risk probability (1 - survival)",
    title_text = paste0("Estimated risk probability: ", group_title),
    pnu_reference_year = pnu_reference_year
  )

  plot_registry[[paste0(plot_group, "__survival_curve")]] <- survival_plot
  plot_registry[[paste0(plot_group, "__risk_curve")]] <- risk_plot

  save_plot_png(
    survival_plot,
    make_plot_output_file(export_path, "estimated_survival_curve", plot_group)
  )
  save_plot_png(
    risk_plot,
    make_plot_output_file(export_path, "estimated_risk_curve", plot_group)
  )
}

remove_existing_detail_plot_outputs(
  export_dir = export_path,
  annotation_file = detail_annotation_file
)

for (ii in seq_len(nrow(detail_plot_group_registry))) {
  plot_group <- detail_plot_group_registry$plot_group[[ii]]
  group_title <- detail_plot_group_registry$group_title[[ii]]
  group_dataset_keys <- detail_plot_group_registry$dataset_keys[[ii]]
  source_dataset_name <- detail_plot_group_registry$source_dataset[[ii]]
  include_pnu_reference <- isTRUE(detail_plot_group_registry$include_pnu_reference[[ii]])
  analysis_df <- analysis_datasets[[source_dataset_name]]
  group_curve_df <- plot_source_df %>%
    filter(dataset %in% group_dataset_keys)

  if (nrow(group_curve_df) == 0L) {
    stop(sprintf("No plotting rows found for detail plot group `%s`.", plot_group), call. = FALSE)
  }

  pnu_reference_year <- if (include_pnu_reference) pnu_max_observed_followup_year else NA_real_

  detail_annotation_list[[ii]] <- build_detail_annotation_table(
    analysis_df = analysis_df,
    times = seq(0, curve_horizon_max_year, by = detail_tick_step_year)
  ) %>%
    mutate(
      plot_group = plot_group,
      source_dataset = source_dataset_name
    )

  survival_detail_plot_bundle <- make_detail_curve_plot(
    curve_df = group_curve_df,
    value_col = "overall_survival_prob",
    y_label = "Estimated survival probability",
    title_text = paste0("Estimated survival probability with risk counts: ", group_title),
    analysis_df = analysis_df,
    pnu_reference_year = pnu_reference_year
  )

  risk_detail_plot_bundle <- make_detail_curve_plot(
    curve_df = group_curve_df,
    value_col = "overall_risk_prob",
    y_label = "Estimated risk probability (1 - survival)",
    title_text = paste0("Estimated risk probability with risk counts: ", group_title),
    analysis_df = analysis_df,
    pnu_reference_year = pnu_reference_year
  )

  plot_registry[[paste0(plot_group, "__survival_curve_with_counts")]] <- survival_detail_plot_bundle
  plot_registry[[paste0(plot_group, "__risk_curve_with_counts")]] <- risk_detail_plot_bundle

  save_stacked_plot_png(
    survival_detail_plot_bundle,
    make_plot_output_file(export_path, "estimated_survival_curve_with_counts", plot_group),
    width = detail_plot_width_in,
    height = detail_plot_height_in
  )
  save_stacked_plot_png(
    risk_detail_plot_bundle,
    make_plot_output_file(export_path, "estimated_risk_curve_with_counts", plot_group),
    width = detail_plot_width_in,
    height = detail_plot_height_in
  )
}

detail_annotation_df <- bind_rows(detail_annotation_list)

# Write Outputs -----------------------------------------------------------
readr::write_csv(horizon_summary_df, horizon_summary_file)
readr::write_csv(plot_source_df, plot_source_file)
readr::write_csv(detail_annotation_df, detail_annotation_file)
readr::write_csv(model_registry_df, model_registry_file)
readr::write_csv(parameter_summary_df, parameter_summary_file)
saveRDS(plot_registry, plot_rds_file)
write_trace_plot_pdf(fit_records = fit_records, output_file = trace_plot_pdf_file)

message("Bayesian mixture cure export completed.")
message("Horizon summary CSV: ", normalizePath(horizon_summary_file, winslash = "/", mustWork = FALSE))
message("Plot source CSV: ", normalizePath(plot_source_file, winslash = "/", mustWork = FALSE))
message("Detail annotation CSV: ", normalizePath(detail_annotation_file, winslash = "/", mustWork = FALSE))
message("Model registry CSV: ", normalizePath(model_registry_file, winslash = "/", mustWork = FALSE))
message("Parameter summary CSV: ", normalizePath(parameter_summary_file, winslash = "/", mustWork = FALSE))
message("Plot RDS: ", normalizePath(plot_rds_file, winslash = "/", mustWork = FALSE))
message("Trace plot PDF: ", normalizePath(trace_plot_pdf_file, winslash = "/", mustWork = FALSE))
message("Per-model fit RDS files are stored in: ", normalizePath(export_path, winslash = "/", mustWork = FALSE))
