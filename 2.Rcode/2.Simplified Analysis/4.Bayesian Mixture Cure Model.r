# Configure Paths ---------------------------------------------------------
find_repo_root <- function(start_dir) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)

  repeat {
    stage8_reference <- file.path(current_dir, "2.Rcode", "stage8A_Bayesian transition-only cure.r")
    if (file.exists(stage8_reference)) {
      return(current_dir)
    }

    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      break
    }
    current_dir <- parent_dir
  }

  stop(
    "Could not locate the repository root containing `2.Rcode/stage8A_Bayesian transition-only cure.r`.",
    call. = FALSE
  )
}

command_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", command_args, value = TRUE)
search_start_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

repo_root <- find_repo_root(search_start_dir)
results_root <- file.path(repo_root, "3.Results files")

sys_name <- Sys.info()[["sysname"]]
dropbox_results_root <- switch(
  sys_name,
  "Darwin" = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  "Windows" = "C:/Users/clair/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  stop("Unsupported OS: ", sys_name)
)

merged_data_path <- Sys.getenv(
  "STAGE8A_SIMPLE_MERGED_DATA_PATH",
  unset = file.path(dropbox_results_root, "data", "MERGED_dataset3_pnu_snu.csv")
)
export_path <- Sys.getenv(
  "STAGE8A_SIMPLE_EXPORT_PATH",
  unset = results_root
)

pnu_site_label <- "PNU"
snu_site_label <- "SNU"

main_risk_scale <- "transition_only_main"
horizon_years <- 1:10
curve_horizon_max_year <- 10
curve_step_year <- 0.05
time_origin_epsilon_year <- 1e-08

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

horizon_summary_file <- file.path(
  export_path,
  "stage8A_simple_bayesian_transition_only_horizon_summary.csv"
)
curve_data_file <- file.path(
  export_path,
  "stage8A_simple_bayesian_transition_only_curve_data.csv"
)
fit_rds_path <- function(model_id, export_dir = export_path) {
  file.path(
    export_dir,
    paste0("stage8A_simple_bayesian_transition_only__", model_id, "__fit.rds")
  )
}
model_registry_file <- file.path(
  export_path,
  "stage8A_simple_bayesian_transition_only_model_registry.csv"
)
parameter_summary_file <- file.path(
  export_path,
  "stage8A_simple_bayesian_transition_only_parameter_summary.csv"
)
plot_rds_file <- file.path(
  export_path,
  "stage8A_simple_bayesian_transition_only_plot_objects.rds"
)
trace_plot_pdf_file <- file.path(
  export_path,
  "stage8A_simple_bayesian_transition_only_trace_plots.pdf"
)
overall_survival_png_file <- file.path(
  export_path,
  "stage8A_simple_bayesian_transition_only_overall_survival.png"
)
overall_risk_png_file <- file.path(
  export_path,
  "stage8A_simple_bayesian_transition_only_overall_risk.png"
)
susceptible_survival_png_file <- file.path(
  export_path,
  "stage8A_simple_bayesian_transition_only_susceptible_only_survival.png"
)
susceptible_risk_png_file <- file.path(
  export_path,
  "stage8A_simple_bayesian_transition_only_susceptible_only_risk.png"
)

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

dataset_palette <- c(
  "PNU" = "#1B4332",
  "SNU" = "#2A6F97",
  "merged_no_site" = "#C1666B",
  "merged_site_adjusted" = "#B8860B"
)

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
rstan_options(auto_write = FALSE)
detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}
options(mc.cores = max(1L, min(stan_chains, detected_cores)))
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# Define Helpers ----------------------------------------------------------
safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

clip_prob <- function(x, eps = 1e-12) {
  pmin(pmax(x, eps), 1 - eps)
}

read_delimited_or_rds <- function(path) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))
  if (grepl("\\.csv\\.gz$", path, ignore.case = TRUE)) {
    return(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  }
  if (ext == "csv") {
    return(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
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

# Define Data Preparation -------------------------------------------------
prepare_analysis_dataset <- function(df, dataset_name, pnu_label, snu_label) {
  if (is.null(df) || !inherits(df, "data.frame")) {
    stop("Input dataset is NULL or not a data frame for ", dataset_name, call. = FALSE)
  }

  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop("[", dataset_name, "] Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  out <- tibble::as_tibble(df) %>%
    mutate(
      id = trimws(as.character(id)),
      site = trimws(as.character(site)),
      site = case_when(
        toupper(site) == toupper(pnu_label) ~ pnu_label,
        toupper(site) == toupper(snu_label) ~ snu_label,
        TRUE ~ site
      ),
      sex_num = as.integer(safe_numeric(sex_num)),
      age_exact_entry = safe_numeric(age_exact_entry),
      days_followup = safe_numeric(days_followup),
      status_num = as.integer(safe_numeric(status_num))
    )

  if (nrow(out) == 0L) {
    stop("[", dataset_name, "] Dataset has zero rows.", call. = FALSE)
  }
  if (anyNA(out[required_cols])) {
    stop("[", dataset_name, "] Missing values detected in required backbone columns.", call. = FALSE)
  }
  if (any(out$id == "")) {
    stop("[", dataset_name, "] Blank `id` values detected.", call. = FALSE)
  }
  if (any(!out$sex_num %in% c(0L, 1L))) {
    stop("[", dataset_name, "] `sex_num` must be coded 0/1.", call. = FALSE)
  }
  if (any(!out$status_num %in% c(0L, 1L, 2L))) {
    stop("[", dataset_name, "] `status_num` must be coded 0/1/2.", call. = FALSE)
  }
  if (any(out$days_followup < 0)) {
    stop("[", dataset_name, "] Negative `days_followup` values detected.", call. = FALSE)
  }

  if (dataset_name != "merged") {
    out <- out %>% mutate(site = dataset_name)
  }

  age_mean <- mean(out$age_exact_entry)
  age_sd <- stats::sd(out$age_exact_entry)
  if (is.na(age_sd) || age_sd <= 0) {
    stop("[", dataset_name, "] `age_exact_entry` must have positive SD.", call. = FALSE)
  }

  out <- out %>%
    mutate(
      unique_person_id = paste(site, id, sep = "_"),
      time_year = pmax(days_followup / 365.25, time_origin_epsilon_year),
      age_s = (age_exact_entry - age_mean) / (2 * age_sd),
      # Transition-only event rule:
      # - event_main = 1 only when status_num == 1
      # - status_num == 0 or 2 are treated as censoring in this simplified script
      event_main = as.integer(status_num == 1L),
      censor_main = as.integer(status_num %in% c(0L, 2L))
    )

  if (dplyr::n_distinct(out$unique_person_id) != nrow(out)) {
    stop("[", dataset_name, "] `site + id` is not unique.", call. = FALSE)
  }
  if (dataset_name == "merged" && dplyr::n_distinct(out$site) < 2L) {
    stop("[merged] Merged dataset must contain at least two site levels.", call. = FALSE)
  }

  list(
    data = out,
    summary = tibble(
      source_dataset_key = dataset_name,
      n_subjects = nrow(out),
      n_transition = sum(out$event_main == 1L),
      n_censored = sum(out$censor_main == 1L)
    )
  )
}

load_transition_only_datasets <- function(merged_path, pnu_label, snu_label) {
  merged_raw <- read_delimited_or_rds(merged_path)
  if (!inherits(merged_raw, "data.frame")) {
    stop("Merged dataset could not be loaded.", call. = FALSE)
  }

  merged_raw <- tibble::as_tibble(merged_raw)
  if (!("site" %in% names(merged_raw))) {
    stop("Merged dataset must contain a `site` column.", call. = FALSE)
  }

  merged_raw <- merged_raw %>%
    mutate(
      site = trimws(as.character(site)),
      site = case_when(
        toupper(site) == toupper(pnu_label) ~ pnu_label,
        toupper(site) == toupper(snu_label) ~ snu_label,
        TRUE ~ site
      )
    )

  pnu_raw <- merged_raw %>% filter(site == pnu_label)
  snu_raw <- merged_raw %>% filter(site == snu_label)

  if (nrow(pnu_raw) == 0L) {
    stop("No rows found for PNU site label `", pnu_label, "`.", call. = FALSE)
  }
  if (nrow(snu_raw) == 0L) {
    stop("No rows found for SNU site label `", snu_label, "`.", call. = FALSE)
  }

  pnu_obj <- prepare_analysis_dataset(pnu_raw, "PNU", pnu_label, snu_label)
  snu_obj <- prepare_analysis_dataset(snu_raw, "SNU", pnu_label, snu_label)
  merged_obj <- prepare_analysis_dataset(merged_raw, "merged", pnu_label, snu_label)

  list(
    datasets = list(
      PNU = pnu_obj$data,
      SNU = snu_obj$data,
      merged = merged_obj$data
    ),
    dataset_summary = bind_rows(pnu_obj$summary, snu_obj$summary, merged_obj$summary)
  )
}

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

# Define Model Specifications ---------------------------------------------
build_model_spec_registry <- function() {
  # Retained main-model rule for the merged site-adjusted branch:
  # follow the Stage 8A reference selection that places `site` in both the
  # incidence part and the latency part (`site_in_both`).
  tibble(
    dataset = c("PNU", "SNU", "merged_no_site", "merged_site_adjusted"),
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
compile_stage8a_simple_stan_model <- function() {
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
    model_name = "stage8A_simple_transition_only_bayesian_lognormal_cure"
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

make_model_registry_row <- function(model_row, analysis_df, design_bundle, elapsed_seconds, divergent_transitions, treedepth_hits, prediction_draws, parameter_summary_df) {
  tibble(
    dataset = model_row$dataset[[1L]],
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

# Define Prediction Helpers -----------------------------------------------
build_curve_prediction_df_from_draws <- function(model_row, n_subjects, design_bundle, draw_bundle, times) {
  # The export tables are built by first computing subject-level predictions,
  # then averaging those subject-level posterior means to the cohort level.
  X_inc <- design_bundle$X_inc
  X_lat <- design_bundle$X_lat

  eta_inc_mat <- X_inc %*% t(draw_bundle$beta_inc)
  eta_inc_mat <- sweep(eta_inc_mat, 2, draw_bundle$alpha_inc, "+")
  pi_mat <- clip_prob(plogis(eta_inc_mat))

  mu_lat_mat <- X_lat %*% t(draw_bundle$gamma_lat)
  mu_lat_mat <- sweep(mu_lat_mat, 2, draw_bundle$gamma0, "+")
  sigma_vec <- exp(draw_bundle$log_sigma)

  subject_susceptible_fraction <- rowMeans(pi_mat)

  curve_rows <- lapply(times, function(tt) {
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
      model_id = model_row$model_id[[1L]],
      site_adjustment_flag = as.logical(model_row$site_adjustment_flag[[1L]]),
      risk_scale = main_risk_scale,
      time_year = as.numeric(tt),
      susceptible_fraction = mean(subject_susceptible_fraction),
      overall_survival_prob = mean(subject_overall_survival),
      overall_risk_prob = mean(1 - subject_overall_survival),
      susceptible_only_survival_prob = mean(subject_susceptible_only_survival),
      susceptible_only_risk_prob = mean(1 - subject_susceptible_only_survival),
      n_subjects = as.integer(n_subjects)
    )
  })

  bind_rows(curve_rows)
}

# Build Plots -------------------------------------------------------------
make_curve_plot <- function(curve_df, value_col, y_label, title_text) {
  yearly_points <- curve_df %>%
    filter(time_year %in% as.numeric(horizon_years))

  ggplot(curve_df, aes(x = time_year, y = .data[[value_col]], color = dataset)) +
    geom_line(linewidth = 1.05) +
    geom_point(
      data = yearly_points,
      aes(x = time_year, y = .data[[value_col]], color = dataset),
      size = 1.8
    ) +
    scale_color_manual(values = dataset_palette) +
    scale_x_continuous(
      breaks = 0:curve_horizon_max_year,
      limits = c(0, curve_horizon_max_year),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = title_text,
      subtitle = "Cohort-level means of subject-level posterior predictions from simplified Stage 8A Bayesian log-normal mixture cure models",
      x = "Years after cohort entry",
      y = y_label,
      color = "Dataset"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
}

save_plot_png <- function(plot_object, output_file) {
  ggplot2::ggsave(
    filename = output_file,
    plot = plot_object,
    width = plot_width_in,
    height = plot_height_in,
    dpi = plot_dpi,
    units = "in"
  )
}

# Load Data and Compile Model ---------------------------------------------
loaded_objects <- load_transition_only_datasets(
  merged_path = merged_data_path,
  pnu_label = pnu_site_label,
  snu_label = snu_site_label
)

analysis_datasets <- loaded_objects$datasets
dataset_summary <- loaded_objects$dataset_summary
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
    "Compiling simplified Stage 8A Bayesian log-normal mixture cure model for ",
    n_models_requiring_fit,
    " model(s) that still need fitting."
  )
  stan_model_compiled <- compile_stage8a_simple_stan_model()
} else {
  message("All model fit RDS files already exist. Stan compilation skipped.")
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
    n_subjects = nrow(analysis_df),
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
    if (isTRUE(reused_existing_fit)) "reuse completed" else paste0("completed in ", formatC(elapsed_seconds, digits = 1L, format = "f"), " sec"),
    "; divergent=",
    divergent_transitions,
    ", treedepth_hits=",
    treedepth_hits,
    ", prediction_draws=",
    draw_bundle$n_draws
  )
}

curve_data_df <- bind_rows(curve_data_list) %>%
  mutate(time_year = round(as.numeric(time_year), 10)) %>%
  arrange(match(dataset, model_spec_registry$dataset), time_year)
parameter_summary_df <- bind_rows(parameter_summary_list) %>%
  arrange(match(dataset, model_spec_registry$dataset), model_id, parameter)
model_registry_df <- bind_rows(model_registry_list) %>%
  arrange(match(dataset, model_spec_registry$dataset), model_id)

# Derive Horizon Summary --------------------------------------------------
horizon_lookup <- tibble(
  time_year = as.numeric(horizon_years),
  horizon_year = as.integer(horizon_years)
)

horizon_summary_df <- curve_data_df %>%
  inner_join(horizon_lookup, by = "time_year") %>%
  transmute(
    dataset = dataset,
    model_id = model_id,
    site_adjustment_flag = site_adjustment_flag,
    risk_scale = risk_scale,
    horizon_year = horizon_year,
    susceptible_fraction = susceptible_fraction,
    overall_survival_prob = overall_survival_prob,
    overall_risk_prob = overall_risk_prob,
    susceptible_only_survival_prob = susceptible_only_survival_prob,
    susceptible_only_risk_prob = susceptible_only_risk_prob,
    n_subjects = n_subjects
  ) %>%
  arrange(match(dataset, model_spec_registry$dataset), horizon_year)

# Build Plot Objects ------------------------------------------------------
overall_survival_plot <- make_curve_plot(
  curve_df = curve_data_df,
  value_col = "overall_survival_prob",
  y_label = "Overall survival probability",
  title_text = "Stage 8A simplified Bayesian mixture cure: overall survival"
)

overall_risk_plot <- make_curve_plot(
  curve_df = curve_data_df,
  value_col = "overall_risk_prob",
  y_label = "Overall cumulative risk",
  title_text = "Stage 8A simplified Bayesian mixture cure: overall cumulative risk"
)

susceptible_survival_plot <- make_curve_plot(
  curve_df = curve_data_df,
  value_col = "susceptible_only_survival_prob",
  y_label = "Susceptible-only survival probability",
  title_text = "Stage 8A simplified Bayesian mixture cure: susceptible-only survival"
)

susceptible_risk_plot <- make_curve_plot(
  curve_df = curve_data_df,
  value_col = "susceptible_only_risk_prob",
  y_label = "Susceptible-only cumulative risk",
  title_text = "Stage 8A simplified Bayesian mixture cure: susceptible-only cumulative risk"
)

plot_registry <- list(
  overall_survival_curve = overall_survival_plot,
  overall_risk_curve = overall_risk_plot,
  susceptible_only_survival_curve = susceptible_survival_plot,
  susceptible_only_risk_curve = susceptible_risk_plot
)

# Write Outputs -----------------------------------------------------------
readr::write_csv(horizon_summary_df, horizon_summary_file)
readr::write_csv(curve_data_df, curve_data_file)
readr::write_csv(model_registry_df, model_registry_file)
readr::write_csv(parameter_summary_df, parameter_summary_file)
saveRDS(plot_registry, plot_rds_file)
write_trace_plot_pdf(fit_records = fit_records, output_file = trace_plot_pdf_file)

save_plot_png(overall_survival_plot, overall_survival_png_file)
save_plot_png(overall_risk_plot, overall_risk_png_file)
save_plot_png(susceptible_survival_plot, susceptible_survival_png_file)
save_plot_png(susceptible_risk_plot, susceptible_risk_png_file)

message("Stage 8A simplified Bayesian transition-only cure export completed.")
message("Horizon summary CSV: ", normalizePath(horizon_summary_file, winslash = "/", mustWork = FALSE))
message("Curve data CSV: ", normalizePath(curve_data_file, winslash = "/", mustWork = FALSE))
message("Model registry CSV: ", normalizePath(model_registry_file, winslash = "/", mustWork = FALSE))
message("Parameter summary CSV: ", normalizePath(parameter_summary_file, winslash = "/", mustWork = FALSE))
message("Plot RDS: ", normalizePath(plot_rds_file, winslash = "/", mustWork = FALSE))
message("Trace plot PDF: ", normalizePath(trace_plot_pdf_file, winslash = "/", mustWork = FALSE))
