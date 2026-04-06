# Configure Paths ---------------------------------------------------------
find_repo_root <- function(start_dir) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)

  repeat {
    has_repo_markers <- dir.exists(file.path(current_dir, ".git")) ||
      (
        dir.exists(file.path(current_dir, "0.Data")) &&
          dir.exists(file.path(current_dir, "2.Rcode"))
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
block1_root_default <- "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/1.Block1"
dropbox_root <- Sys.getenv(
  "CHRPMIXTURE_DROPBOX_ROOT",
  unset = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model"
)
source_data_file <- Sys.getenv(
  "SIMPLE_MODEL_COMPARISON_SOURCE_DATA_FILE",
  unset = file.path(
    dropbox_root,
    "0.Data",
    "2.Preprocessed data",
    "Preprocessed_Merged_PNUH_SNUH_Data.csv"
  )
)
modeling_root <- Sys.getenv(
  "SIMPLE_MODELING_RESULTS_ROOT",
  unset = block1_root_default
)
results_integration_root <- Sys.getenv(
  "BLOCK1_RESULTS_INTEGRATION_ROOT",
  unset = file.path(block1_root_default, "5.Results Integration")
)
export_path <- Sys.getenv(
  "SIMPLE_MODEL_COMPARISON_TABLE_EXPORT_PATH",
  unset = file.path(results_integration_root, "simple_model_probability_comparison_tables")
)

stage3_input_file <- file.path(modeling_root, "1.KM", "simple_km_yearly_estimates.csv")
stage5_input_file <- file.path(
  modeling_root,
  "2.MLE No-Cure",
  "mle_nocure_lognormal_horizon_summary.csv"
)
stage5_fit_file <- file.path(
  modeling_root,
  "2.MLE No-Cure",
  "mle_nocure_lognormal_fitted_models.rds"
)
stage7_input_file <- file.path(
  modeling_root,
  "3.MLE Mixture Cure",
  "mle_mixture_cure_lognormal_horizon_summary.csv"
)
stage7_fit_file <- file.path(
  modeling_root,
  "3.MLE Mixture Cure",
  "mle_mixture_cure_lognormal_fitted_objects.rds"
)
stage8_input_candidates <- c(
  file.path(modeling_root, "4.Bayesian Mixture Cure", "bayesian_mixture_cure_horizon_summary.csv"),
  file.path(modeling_root, "4.Bayesian Mixture Cure", "bayesian_mixture_cure_prior_branch_horizon_summary.csv"),
  file.path(repo_root, "3.Results files", "stage8A_simple_bayesian_transition_only_horizon_summary.csv")
)

required_horizons <- 1:10
round_digits <- as.integer(Sys.getenv("SIMPLE_MODEL_COMPARISON_ROUND_DIGITS", unset = "6"))
stage5_ci_draws <- as.integer(Sys.getenv("SIMPLE_MODEL_COMPARISON_STAGE5_CI_DRAWS", unset = "4000"))
stage7_ci_draws <- as.integer(Sys.getenv("SIMPLE_MODEL_COMPARISON_STAGE7_CI_DRAWS", unset = "4000"))
delta_interval_draws <- as.integer(Sys.getenv("SIMPLE_MODEL_COMPARISON_DELTA_DRAWS", unset = "4000"))
pnu_site_label <- Sys.getenv("SIMPLE_MODEL_COMPARISON_PNU_SITE_LABEL", unset = "PNU")
snu_site_label <- Sys.getenv("SIMPLE_MODEL_COMPARISON_SNU_SITE_LABEL", unset = "SNU")
time_origin_epsilon_year <- 1e-08

# Load Packages -----------------------------------------------------------
required_packages <- c("readr", "dplyr", "tidyr", "tibble", "flexsurv")
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
  library(tidyr)
  library(tibble)
  library(flexsurv)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# Define Helpers ----------------------------------------------------------
assert_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

resolve_existing_path <- function(paths, label) {
  existing_paths <- paths[file.exists(paths)]

  if (length(existing_paths) == 0L) {
    stop(
      sprintf(
        "%s not found. Checked paths:\n- %s",
        label,
        paste(paths, collapse = "\n- ")
      ),
      call. = FALSE
    )
  }

  normalizePath(existing_paths[[1L]], winslash = "/", mustWork = TRUE)
}

read_csv_checked <- function(path, label) {
  assert_file_exists(path, label)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

read_rds_checked <- function(path, label) {
  assert_file_exists(path, label)
  readRDS(path)
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

normalize_dataset_variant <- function(x) {
  x <- trimws(as.character(x))

  dplyr::case_when(
    toupper(x) == "PNU" ~ "PNU",
    toupper(x) == "SNU" ~ "SNU",
    tolower(x) %in% c("merged", "merged_no_site") ~ "merged",
    tolower(x) == "merged_site_adjusted" ~ "merged_site_adjusted",
    TRUE ~ x
  )
}

dataset_group_from_variant <- function(dataset_variant) {
  dataset_variant <- normalize_dataset_variant(dataset_variant)

  dplyr::case_when(
    dataset_variant == "PNU" ~ "PNU",
    dataset_variant == "SNU" ~ "SNU",
    dataset_variant %in% c("merged", "merged_site_adjusted") ~ "Merged",
    TRUE ~ dataset_variant
  )
}

clip_probability <- function(x) {
  pmin(pmax(as.numeric(x), 0), 1)
}

clip_probability_open <- function(x, eps = 1e-12) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

duplicate_dataset_rows <- function(df, source_dataset_variant, target_dataset_variant) {
  duplicated_rows <- df %>%
    filter(dataset_variant == source_dataset_variant) %>%
    mutate(
      dataset_variant = target_dataset_variant,
      dataset_group = dataset_group_from_variant(target_dataset_variant)
    )

  bind_rows(df, duplicated_rows)
}

round_probability <- function(x) {
  round(as.numeric(x), digits = round_digits)
}

standardize_known_site_labels <- function(df, pnu_label, snu_label) {
  if (!("site" %in% names(df))) {
    stop("Input data must contain a `site` column.", call. = FALSE)
  }

  df %>%
    mutate(
      site = trimws(as.character(site)),
      site = dplyr::case_when(
        toupper(site) == toupper(pnu_label) ~ pnu_label,
        toupper(site) == toupper(snu_label) ~ snu_label,
        TRUE ~ site
      )
    )
}

validate_source_data <- function(df, pnu_label, snu_label) {
  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "Merged preprocessed data is missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  observed_sites <- sort(unique(trimws(as.character(df$site))))
  expected_sites <- c(pnu_label, snu_label)
  unexpected_sites <- setdiff(observed_sites, expected_sites)

  if (length(unexpected_sites) > 0L) {
    stop(
      sprintf(
        "Unexpected site labels found in merged input: %s",
        paste(unexpected_sites, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

split_site_dataset <- function(merged_df, site_label, dataset_name) {
  out <- merged_df %>%
    filter(toupper(trimws(as.character(site))) == toupper(site_label)) %>%
    mutate(site = site_label)

  if (nrow(out) == 0L) {
    stop(
      sprintf("[%s] No rows found in merged input for site label `%s`.", dataset_name, site_label),
      call. = FALSE
    )
  }

  out
}

prepare_analysis_dataset <- function(df, dataset_name, site_mode = c("single", "merged")) {
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
    mutate(
      id = trimws(as.character(id)),
      site = trimws(as.character(site)),
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
  if (any(!out$sex_num %in% c(0L, 1L), na.rm = TRUE)) {
    stop(sprintf("[%s] `sex_num` must be coded 0/1.", dataset_name), call. = FALSE)
  }
  if (any(!out$status_num %in% c(0L, 1L, 2L), na.rm = TRUE)) {
    stop(sprintf("[%s] `status_num` must be coded 0/1/2.", dataset_name), call. = FALSE)
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

  out %>%
    mutate(
      unique_person_id = paste(site, id, sep = "_"),
      time_year = pmax(days_followup / 365.25, time_origin_epsilon_year),
      time_year_model = pmax(days_followup / 365.25, time_origin_epsilon_year),
      age_s = (age_exact_entry - age_mean) / (2 * age_sd),
      event_main = as.integer(status_num == 1L),
      site = factor(as.character(site))
    )
}

build_analysis_datasets_from_source <- function(source_file, pnu_label, snu_label) {
  assert_file_exists(source_file, "Merged preprocessed data CSV")
  source_df <- readr::read_csv(
    source_file,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
  source_df <- standardize_known_site_labels(source_df, pnu_label, snu_label)
  validate_source_data(source_df, pnu_label, snu_label)

  list(
    PNU = prepare_analysis_dataset(
      split_site_dataset(source_df, pnu_label, "PNU"),
      dataset_name = "PNU",
      site_mode = "single"
    ),
    SNU = prepare_analysis_dataset(
      split_site_dataset(source_df, snu_label, "SNU"),
      dataset_name = "SNU",
      site_mode = "single"
    ),
    merged = prepare_analysis_dataset(
      source_df,
      dataset_name = "merged",
      site_mode = "merged"
    )
  )
}

trim_rhs <- function(rhs) {
  rhs <- gsub("\\s+", " ", trimws(as.character(rhs)))
  rhs[is.na(rhs)] <- NA_character_
  rhs
}

make_rhs_formula <- function(rhs) {
  rhs <- trim_rhs(rhs)
  if (is.na(rhs) || rhs == "") {
    stats::as.formula("~ 1")
  } else {
    stats::as.formula(paste("~", rhs))
  }
}

model_matrix_from_rhs <- function(df, rhs) {
  stats::model.matrix(make_rhs_formula(rhs), data = df)
}

lognormal_surv_density <- function(time, lp, log_sigma) {
  time <- pmax(as.numeric(time), time_origin_epsilon_year)
  lp <- as.numeric(lp)
  sigma <- exp(as.numeric(log_sigma))
  z <- (log(time) - lp) / sigma
  surv <- 1 - stats::pnorm(z)
  dens <- stats::dnorm(z) / (sigma * time)

  list(
    surv = pmin(pmax(surv, 0), 1),
    dens = pmax(dens, 1e-300)
  )
}

make_psd_matrix <- function(x, eigen_floor = 1e-10) {
  eig <- eigen(x, symmetric = TRUE)
  eig$values[eig$values < eigen_floor] <- eigen_floor
  eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
}

invert_hessian <- function(hessian_mat, eigen_floor = 1e-8) {
  sym_hessian <- 0.5 * (hessian_mat + t(hessian_mat))
  eig <- eigen(sym_hessian, symmetric = TRUE)
  eig$values[eig$values < eigen_floor] <- eigen_floor
  eig$vectors %*% diag(1 / eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
}

draw_mvn <- function(n_draws, mu, sigma, seed_value) {
  sigma_psd <- make_psd_matrix(sigma)
  eig <- eigen(sigma_psd, symmetric = TRUE)
  transform_mat <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow = length(eig$values))

  set.seed(seed_value)
  z <- matrix(stats::rnorm(length(mu) * n_draws), nrow = length(mu), ncol = n_draws)
  draws <- matrix(mu, nrow = length(mu), ncol = n_draws) + transform_mat %*% z
  t(draws)
}

summarize_draws <- function(x) {
  tibble(
    uncertainty_sd = stats::sd(x, na.rm = TRUE),
    interval_lower_95 = unname(stats::quantile(x, probs = 0.025, na.rm = TRUE, names = FALSE)),
    interval_median_50 = unname(stats::quantile(x, probs = 0.5, na.rm = TRUE, names = FALSE)),
    interval_upper_95 = unname(stats::quantile(x, probs = 0.975, na.rm = TRUE, names = FALSE))
  )
}

simulate_probability_draws_from_summary <- function(
    estimate_probability,
    interval_lower_95,
    interval_upper_95,
    interval_median_50 = NA_real_,
    n_draws,
    seed_value
) {
  center_probability <- if (is.finite(interval_median_50)) {
    interval_median_50
  } else {
    estimate_probability
  }

  center_probability <- clip_probability_open(center_probability)
  lower_probability <- suppressWarnings(as.numeric(interval_lower_95))
  upper_probability <- suppressWarnings(as.numeric(interval_upper_95))

  if (!is.finite(lower_probability) || !is.finite(upper_probability)) {
    return(rep(clip_probability(estimate_probability), n_draws))
  }

  lower_probability <- clip_probability_open(lower_probability)
  upper_probability <- clip_probability_open(upper_probability)

  if ((upper_probability - lower_probability) <= 1e-12) {
    return(rep(clip_probability(estimate_probability), n_draws))
  }

  logit_center <- stats::qlogis(center_probability)
  logit_lower <- stats::qlogis(lower_probability)
  logit_upper <- stats::qlogis(upper_probability)
  logit_sd <- (logit_upper - logit_lower) / (2 * stats::qnorm(0.975))

  if (!is.finite(logit_sd) || logit_sd <= 0) {
    return(rep(clip_probability(estimate_probability), n_draws))
  }

  set.seed(seed_value)
  clip_probability(stats::plogis(stats::rnorm(n_draws, mean = logit_center, sd = logit_sd)))
}

build_probability_draw_lookup <- function(long_df, n_draws) {
  key_df <- long_df %>%
    distinct(
      dataset_variant,
      dataset_group,
      model_family,
      model_scope,
      probability_type,
      horizon_year,
      estimate_probability,
      interval_type,
      interval_lower_95,
      interval_median_50,
      interval_upper_95
    ) %>%
    arrange(
      factor(dataset_group, levels = c("PNU", "SNU", "Merged")),
      factor(dataset_variant, levels = c("PNU", "SNU", "merged", "merged_site_adjusted")),
      model_family,
      model_scope,
      probability_type,
      horizon_year
    )

  draw_list <- lapply(seq_len(nrow(key_df)), function(ii) {
    simulate_probability_draws_from_summary(
      estimate_probability = key_df$estimate_probability[[ii]],
      interval_lower_95 = key_df$interval_lower_95[[ii]],
      interval_upper_95 = key_df$interval_upper_95[[ii]],
      interval_median_50 = key_df$interval_median_50[[ii]],
      n_draws = n_draws,
      seed_value = 980000L + ii
    )
  })

  key_df$probability_draws <- draw_list
  key_df
}

delta_interval_type_from_components <- function(interval_type_a, interval_type_b) {
  dplyr::case_when(
    interval_type_a == "95% CrI" & interval_type_b == "95% CrI" ~ "95% CrI",
    interval_type_a == "95% CI" & interval_type_b == "95% CI" ~ "95% CI",
    TRUE ~ "95% uncertainty interval"
  )
}

cleanup_existing_comparison_outputs <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    return(invisible(character()))
  }

  stale_files <- list.files(
    output_dir,
    pattern = ".*\\.csv$",
    full.names = TRUE
  )

  if (length(stale_files) == 0L) {
    return(invisible(character()))
  }

  removed_flag <- file.remove(stale_files)
  if (any(!removed_flag)) {
    stop(
      sprintf(
        "Failed to remove stale comparison-table outputs: %s",
        paste(basename(stale_files[!removed_flag]), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(stale_files[removed_flag])
}

validate_long_table <- function(df) {
  required_columns <- c(
    "dataset_variant",
    "dataset_group",
    "model_family",
    "model_scope",
    "probability_type",
    "horizon_year",
    "estimate_probability",
    "interval_type",
    "interval_lower_95",
    "interval_upper_95"
  )
  missing_columns <- setdiff(required_columns, names(df))

  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "Integrated long table is missing required columns: %s",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  duplicate_rows <- df %>%
    count(
      dataset_variant,
      model_family,
      model_scope,
      probability_type,
      horizon_year,
      name = "row_count"
    ) %>%
    filter(row_count > 1L)

  if (nrow(duplicate_rows) > 0L) {
    stop("Integrated long table contains duplicated rows.", call. = FALSE)
  }

  missing_horizons <- df %>%
    count(
      dataset_variant,
      dataset_group,
      model_family,
      model_scope,
      probability_type,
      name = "n_horizons"
    ) %>%
    filter(n_horizons != length(required_horizons))

  if (nrow(missing_horizons) > 0L) {
    stop("At least one dataset-model-probability series is missing yearly horizons 1-10.", call. = FALSE)
  }

  unexpected_horizons <- setdiff(sort(unique(df$horizon_year)), required_horizons)
  if (length(unexpected_horizons) > 0L) {
    stop("Unexpected yearly horizons detected in integrated long table.", call. = FALSE)
  }

  if (anyNA(df$estimate_probability)) {
    stop("Integrated long table contains missing probabilities.", call. = FALSE)
  }

  if (anyNA(df$interval_type) || anyNA(df$interval_lower_95) || anyNA(df$interval_upper_95)) {
    stop("Integrated long table contains missing interval metadata.", call. = FALSE)
  }

  invisible(TRUE)
}

model_row_label <- function(dataset_group, dataset_variant, model_family, model_scope) {
  if (identical(model_family, "KM benchmark")) {
    base_label <- "KM benchmark"
  } else if (identical(model_family, "Frequentist no-cure lognormal")) {
    base_label <- "Frequentist no-cure lognormal"
  } else if (identical(model_family, "Frequentist mixture cure")) {
    base_label <- if (identical(model_scope, "susceptible")) {
      "Frequentist mixture cure (susceptible only)"
    } else {
      "Frequentist mixture cure (overall)"
    }
  } else if (identical(model_family, "Bayesian transition-only cure")) {
    base_label <- if (identical(model_scope, "susceptible")) {
      "Bayesian transition-only cure (susceptible only)"
    } else {
      "Bayesian transition-only cure (overall)"
    }
  } else {
    base_label <- model_family
  }

  if (!identical(dataset_group, "Merged")) {
    return(base_label)
  }

  variant_suffix <- if (identical(dataset_variant, "merged_site_adjusted")) {
    "[merged(site-adj)]"
  } else {
    "[merged]"
  }

  paste(base_label, variant_suffix)
}

model_row_order <- function(dataset_group, dataset_variant, model_family, model_scope) {
  if (!identical(dataset_group, "Merged")) {
    return(dplyr::case_when(
      model_family == "KM benchmark" ~ 1,
      model_family == "Frequentist no-cure lognormal" ~ 2,
      model_family == "Frequentist mixture cure" & model_scope == "overall" ~ 3,
      model_family == "Frequentist mixture cure" & model_scope == "susceptible" ~ 4,
      model_family == "Bayesian transition-only cure" & model_scope == "overall" ~ 5,
      model_family == "Bayesian transition-only cure" & model_scope == "susceptible" ~ 6,
      TRUE ~ 999
    ))
  }

  dplyr::case_when(
    model_family == "KM benchmark" & dataset_variant == "merged" ~ 1,
    model_family == "KM benchmark" & dataset_variant == "merged_site_adjusted" ~ 2,
    model_family == "Frequentist no-cure lognormal" & dataset_variant == "merged" ~ 3,
    model_family == "Frequentist no-cure lognormal" & dataset_variant == "merged_site_adjusted" ~ 4,
    model_family == "Frequentist mixture cure" & model_scope == "overall" & dataset_variant == "merged" ~ 5,
    model_family == "Frequentist mixture cure" & model_scope == "susceptible" & dataset_variant == "merged" ~ 6,
    model_family == "Frequentist mixture cure" & model_scope == "overall" & dataset_variant == "merged_site_adjusted" ~ 7,
    model_family == "Frequentist mixture cure" & model_scope == "susceptible" & dataset_variant == "merged_site_adjusted" ~ 8,
    model_family == "Bayesian transition-only cure" & model_scope == "overall" & dataset_variant == "merged" ~ 9,
    model_family == "Bayesian transition-only cure" & model_scope == "susceptible" & dataset_variant == "merged" ~ 10,
    model_family == "Bayesian transition-only cure" & model_scope == "overall" & dataset_variant == "merged_site_adjusted" ~ 11,
    model_family == "Bayesian transition-only cure" & model_scope == "susceptible" & dataset_variant == "merged_site_adjusted" ~ 12,
    TRUE ~ 999
  )
}

compute_stage5_interval_lookup <- function(fit_path, analysis_datasets, target_horizons, n_draws) {
  fitted_models <- read_rds_checked(fit_path, "Stage 5 fitted models RDS")

  stage5_registry <- tibble(
    dataset_variant = c("PNU", "SNU", "merged", "merged_site_adjusted"),
    source_dataset = c("PNU", "SNU", "merged", "merged"),
    fit_key = c("PNU", "SNU", "merged", "merged_site_adjusted")
  )

  bind_rows(lapply(seq_len(nrow(stage5_registry)), function(ii) {
    row <- stage5_registry[ii, ]
    fit_key <- as.character(row$fit_key[[1L]])
    fit_obj <- fitted_models[[fit_key]]
    if (is.null(fit_obj)) {
      stop(sprintf("Stage 5 fitted model missing for `%s`.", fit_key), call. = FALSE)
    }

    analysis_df <- analysis_datasets[[as.character(row$source_dataset[[1L]])]]
    model_formula <- fit_obj$concat.formula
    if (is.null(model_formula) && !is.null(fit_obj$all.formulae$meanlog)) {
      model_formula <- fit_obj$all.formulae$meanlog
    }
    if (is.null(model_formula)) {
      stop(sprintf("Stage 5 model formula is unavailable for `%s`.", fit_key), call. = FALSE)
    }
    model_terms <- stats::delete.response(stats::terms(model_formula))
    pred_data <- analysis_df
    design_matrix <- stats::model.matrix(model_terms, data = pred_data)

    parameter_draws <- flexsurv::normboot.flexsurvreg(
      fit_obj,
      B = n_draws,
      newdata = pred_data[1, , drop = FALSE],
      raw = TRUE
    )
    parameter_draws <- as.matrix(parameter_draws)

    if (is.null(colnames(parameter_draws)) || !("sdlog" %in% colnames(parameter_draws))) {
      stop(sprintf("Stage 5 parameter draws are missing `sdlog` for `%s`.", fit_key), call. = FALSE)
    }

    location_draws <- parameter_draws[, setdiff(colnames(parameter_draws), "sdlog"), drop = FALSE]
    colnames(location_draws) <- ifelse(colnames(location_draws) == "meanlog", "(Intercept)", colnames(location_draws))

    missing_design_cols <- setdiff(colnames(design_matrix), colnames(location_draws))
    if (length(missing_design_cols) > 0L) {
      stop(
        sprintf(
          "Stage 5 location draws for `%s` are missing design-matrix columns: %s",
          fit_key,
          paste(missing_design_cols, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    location_draws <- location_draws[, colnames(design_matrix), drop = FALSE]
    sigma_draws <- pmax(as.numeric(parameter_draws[, "sdlog"]), 1e-8)
    sigma_matrix <- matrix(sigma_draws, nrow = nrow(design_matrix), ncol = nrow(parameter_draws), byrow = TRUE)
    lp_matrix <- design_matrix %*% t(location_draws)

    bind_rows(lapply(target_horizons, function(horizon_year) {
      z_matrix <- (log(pmax(as.numeric(horizon_year), time_origin_epsilon_year)) - lp_matrix) / sigma_matrix
      survival_draws <- clip_probability(colMeans(1 - stats::pnorm(z_matrix), na.rm = TRUE))
      risk_draws <- clip_probability(1 - survival_draws)

      bind_rows(
        summarize_draws(survival_draws) %>%
          mutate(
            dataset_variant = as.character(row$dataset_variant[[1L]]),
            dataset_group = dataset_group_from_variant(row$dataset_variant[[1L]]),
            model_family = "Frequentist no-cure lognormal",
            model_scope = "overall",
            probability_type = "survival",
            horizon_year = as.integer(horizon_year),
            interval_type = "95% CI"
          ),
        summarize_draws(risk_draws) %>%
          mutate(
            dataset_variant = as.character(row$dataset_variant[[1L]]),
            dataset_group = dataset_group_from_variant(row$dataset_variant[[1L]]),
            model_family = "Frequentist no-cure lognormal",
            model_scope = "overall",
            probability_type = "risk",
            horizon_year = as.integer(horizon_year),
            interval_type = "95% CI"
          )
      )
    }))
  })) %>%
    select(
      dataset_variant,
      dataset_group,
      model_family,
      model_scope,
      probability_type,
      horizon_year,
      interval_type,
      uncertainty_sd,
      interval_lower_95,
      interval_median_50,
      interval_upper_95
    )
}

compute_stage7_interval_lookup <- function(analysis_datasets, fit_path, target_horizons, n_draws) {
  fitted_objects <- read_rds_checked(fit_path, "Stage 7 fitted objects RDS")

  stage7_registry <- tibble(
    dataset_variant = c("PNU", "SNU", "merged", "merged_site_adjusted"),
    source_dataset = c("PNU", "SNU", "merged", "merged"),
    fit_key = c("PNU", "SNU", "merged_no_site", "merged_site_adjusted")
  )

  bind_rows(lapply(seq_len(nrow(stage7_registry)), function(ii) {
    row <- stage7_registry[ii, ]
    fit_key <- as.character(row$fit_key[[1L]])
    fit_obj <- fitted_objects[[fit_key]]
    if (is.null(fit_obj)) {
      stop(sprintf("Stage 7 fitted object missing for `%s`.", fit_key), call. = FALSE)
    }

    analysis_df <- analysis_datasets[[as.character(row$source_dataset[[1L]])]] %>%
      mutate(site = factor(as.character(site), levels = fit_obj$site_levels))

    X_inc <- model_matrix_from_rhs(analysis_df, fit_obj$incidence_rhs)
    X_lat <- model_matrix_from_rhs(analysis_df, fit_obj$latency_rhs)
    y_time <- pmax(as.numeric(analysis_df$time_year_model), time_origin_epsilon_year)
    y_event <- as.integer(analysis_df$event_main)

    neg_loglik <- function(par) {
      gamma <- par[seq_len(ncol(X_inc))]
      beta <- par[ncol(X_inc) + seq_len(ncol(X_lat))]
      log_sigma <- par[[length(par)]]

      uncured_prob_inner <- clip_probability_open(plogis(drop(X_inc %*% gamma)))
      lp_lat <- drop(X_lat %*% beta)
      dens_surv <- lognormal_surv_density(y_time, lp_lat, log_sigma)
      surv_uncured <- pmax(dens_surv$surv, 1e-300)
      dens_uncured <- pmax(dens_surv$dens, 1e-300)

      loglik <- sum(
        y_event * (log(uncured_prob_inner) + log(dens_uncured)) +
          (1 - y_event) * log((1 - uncured_prob_inner) + uncured_prob_inner * surv_uncured)
      )

      if (!is.finite(loglik)) {
        return(1e12)
      }

      -loglik
    }

    par_hat <- c(fit_obj$gamma, fit_obj$beta, fit_obj$log_sigma)
    hessian_mat <- tryCatch(
      stats::optimHess(par = par_hat, fn = neg_loglik),
      error = function(e) NULL
    )

    if (is.null(hessian_mat) || any(!is.finite(hessian_mat))) {
      stop(sprintf("Stage 7 Hessian computation failed for `%s`.", fit_key), call. = FALSE)
    }

    covariance_mat <- tryCatch(
      invert_hessian(hessian_mat),
      error = function(e) NULL
    )

    if (is.null(covariance_mat) || any(!is.finite(covariance_mat))) {
      stop(sprintf("Stage 7 covariance approximation failed for `%s`.", fit_key), call. = FALSE)
    }

    parameter_draws <- draw_mvn(
      n_draws = n_draws,
      mu = par_hat,
      sigma = covariance_mat,
      seed_value = 920000L + ii
    )

    gamma_draws <- parameter_draws[, seq_len(ncol(X_inc)), drop = FALSE]
    beta_draws <- parameter_draws[, ncol(X_inc) + seq_len(ncol(X_lat)), drop = FALSE]
    log_sigma_draws <- parameter_draws[, length(par_hat)]

    pi_matrix <- clip_probability_open(plogis(X_inc %*% t(gamma_draws)))
    latency_lp_matrix <- X_lat %*% t(beta_draws)
    sigma_matrix <- matrix(exp(log_sigma_draws), nrow = nrow(X_lat), ncol = nrow(parameter_draws), byrow = TRUE)

    bind_rows(lapply(target_horizons, function(horizon_year) {
      z_matrix <- (log(pmax(as.numeric(horizon_year), time_origin_epsilon_year)) - latency_lp_matrix) / sigma_matrix
      susceptible_survival_matrix <- 1 - stats::pnorm(z_matrix)
      overall_survival_draws <- clip_probability(colMeans((1 - pi_matrix) + pi_matrix * susceptible_survival_matrix, na.rm = TRUE))
      overall_risk_draws <- clip_probability(1 - overall_survival_draws)
      susceptible_survival_draws <- clip_probability(colMeans(susceptible_survival_matrix, na.rm = TRUE))
      susceptible_risk_draws <- clip_probability(1 - susceptible_survival_draws)

      bind_rows(
        summarize_draws(overall_survival_draws) %>%
          mutate(
            dataset_variant = as.character(row$dataset_variant[[1L]]),
            dataset_group = dataset_group_from_variant(row$dataset_variant[[1L]]),
            model_family = "Frequentist mixture cure",
            model_scope = "overall",
            probability_type = "survival",
            horizon_year = as.integer(horizon_year),
            interval_type = "95% CI"
          ),
        summarize_draws(overall_risk_draws) %>%
          mutate(
            dataset_variant = as.character(row$dataset_variant[[1L]]),
            dataset_group = dataset_group_from_variant(row$dataset_variant[[1L]]),
            model_family = "Frequentist mixture cure",
            model_scope = "overall",
            probability_type = "risk",
            horizon_year = as.integer(horizon_year),
            interval_type = "95% CI"
          ),
        summarize_draws(susceptible_survival_draws) %>%
          mutate(
            dataset_variant = as.character(row$dataset_variant[[1L]]),
            dataset_group = dataset_group_from_variant(row$dataset_variant[[1L]]),
            model_family = "Frequentist mixture cure",
            model_scope = "susceptible",
            probability_type = "survival",
            horizon_year = as.integer(horizon_year),
            interval_type = "95% CI"
          ),
        summarize_draws(susceptible_risk_draws) %>%
          mutate(
            dataset_variant = as.character(row$dataset_variant[[1L]]),
            dataset_group = dataset_group_from_variant(row$dataset_variant[[1L]]),
            model_family = "Frequentist mixture cure",
            model_scope = "susceptible",
            probability_type = "risk",
            horizon_year = as.integer(horizon_year),
            interval_type = "95% CI"
          )
      )
    }))
  })) %>%
    select(
      dataset_variant,
      dataset_group,
      model_family,
      model_scope,
      probability_type,
      horizon_year,
      interval_type,
      uncertainty_sd,
      interval_lower_95,
      interval_median_50,
      interval_upper_95
    )
}

build_comparison_table <- function(long_df, dataset_group, probability_type) {
  comparison_df <- long_df %>%
    filter(
      dataset_group == !!dataset_group,
      probability_type == !!probability_type
    ) %>%
    mutate(
      row_label = mapply(
        model_row_label,
        dataset_group = dataset_group,
        dataset_variant = dataset_variant,
        model_family = model_family,
        model_scope = model_scope,
        USE.NAMES = FALSE
      ),
      row_order = mapply(
        model_row_order,
        dataset_group = dataset_group,
        dataset_variant = dataset_variant,
        model_family = model_family,
        model_scope = model_scope,
        USE.NAMES = FALSE
      )
    )

  duplicate_rows <- comparison_df %>%
    count(row_label, horizon_year, name = "row_count") %>%
    filter(row_count > 1L)

  if (nrow(duplicate_rows) > 0L) {
    stop(
      sprintf(
        "Duplicated rows detected while building %s / %s comparison table.",
        dataset_group,
        probability_type
      ),
      call. = FALSE
    )
  }

  comparison_df %>%
    transmute(
      model = row_label,
      row_order = row_order,
      interval_type = interval_type,
      horizon_year = horizon_year,
      estimate_probability = round_probability(estimate_probability),
      interval_lower_95 = round_probability(interval_lower_95),
      interval_upper_95 = round_probability(interval_upper_95)
    ) %>%
    tidyr::pivot_longer(
      cols = c(estimate_probability, interval_lower_95, interval_upper_95),
      names_to = "value_type",
      values_to = "value"
    ) %>%
    mutate(
      value_type = dplyr::case_when(
        value_type == "estimate_probability" ~ "estimate",
        value_type == "interval_lower_95" ~ "lower_95",
        value_type == "interval_upper_95" ~ "upper_95",
        TRUE ~ value_type
      ),
      output_column = paste0(horizon_year, " year ", value_type)
    ) %>%
    select(model, row_order, interval_type, output_column, value) %>%
    tidyr::pivot_wider(
      names_from = output_column,
      values_from = value
    ) %>%
    arrange(row_order, model) %>%
    select(-row_order)
}

write_comparison_table <- function(table_df, export_dir, dataset_group, probability_type) {
  dataset_stub <- dplyr::case_when(
    dataset_group == "PNU" ~ "pnu",
    dataset_group == "SNU" ~ "snu",
    dataset_group == "Merged" ~ "merged",
    TRUE ~ tolower(dataset_group)
  )

  probability_stub <- dplyr::case_when(
    probability_type == "survival" ~ "survival_probability",
    probability_type == "risk" ~ "cumulative_risk_probability",
    TRUE ~ probability_type
  )

  output_file <- file.path(
    export_dir,
    sprintf("%s_%s_comparison_table.csv", dataset_stub, probability_stub)
  )

  readr::write_csv(table_df, output_file)
  invisible(output_file)
}

build_probability_difference_summary <- function(long_df, draw_lookup) {
  model_rank_lookup <- c(
    "KM benchmark" = 1,
    "Frequentist no-cure lognormal" = 2,
    "Frequentist mixture cure" = 3,
    "Bayesian transition-only cure" = 4
  )

  overall_df <- long_df %>%
    filter(model_scope == "overall") %>%
    mutate(model_rank = unname(model_rank_lookup[model_family])) %>%
    distinct(
      dataset_variant,
      dataset_group,
      probability_type,
      horizon_year,
      model_family,
      .keep_all = TRUE
    ) %>%
    left_join(
      draw_lookup %>%
        select(
          dataset_variant,
          dataset_group,
          model_family,
          model_scope,
          probability_type,
          horizon_year,
          probability_draws
        ),
      by = c(
        "dataset_variant",
        "dataset_group",
        "model_family",
        "model_scope",
        "probability_type",
        "horizon_year"
      )
    )

  pair_df <- overall_df %>%
    transmute(
      dataset_variant,
      dataset_group,
      probability_type,
      horizon_year,
      model_a = model_family,
      model_rank_a = model_rank,
      estimate_a = estimate_probability,
      interval_type_a = interval_type,
      uncertainty_sd_a = uncertainty_sd,
      interval_lower_95_a = interval_lower_95,
      interval_upper_95_a = interval_upper_95,
      interval_median_50_a = interval_median_50,
      probability_draws_a = probability_draws
    ) %>%
    inner_join(
      overall_df %>%
        transmute(
          dataset_variant,
          dataset_group,
          probability_type,
          horizon_year,
          model_b = model_family,
          model_rank_b = model_rank,
          estimate_b = estimate_probability,
          interval_type_b = interval_type,
          uncertainty_sd_b = uncertainty_sd,
          interval_lower_95_b = interval_lower_95,
          interval_upper_95_b = interval_upper_95,
          interval_median_50_b = interval_median_50,
          probability_draws_b = probability_draws
        ),
      by = c("dataset_variant", "dataset_group", "probability_type", "horizon_year"),
      relationship = "many-to-many"
    ) %>%
    filter(model_rank_a < model_rank_b)

  delta_summary_df <- bind_rows(lapply(seq_len(nrow(pair_df)), function(ii) {
    row <- pair_df[ii, ]
    draws_a <- row$probability_draws_a[[1L]]
    draws_b <- row$probability_draws_b[[1L]]
    delta_draws <- as.numeric(draws_b) - as.numeric(draws_a)

    summarize_draws(delta_draws) %>%
      transmute(
        delta_uncertainty_sd = uncertainty_sd,
        delta_interval_lower_95 = interval_lower_95,
        delta_interval_median_50 = interval_median_50,
        delta_interval_upper_95 = interval_upper_95
      )
  }))

  bind_cols(pair_df, delta_summary_df) %>%
    transmute(
      dataset_variant,
      dataset_group,
      probability_type,
      horizon_year,
      comparison = paste(model_b, "minus", model_a),
      model_a,
      model_b,
      estimate_a = round_probability(estimate_a),
      interval_type_a,
      uncertainty_sd_a = round_probability(uncertainty_sd_a),
      interval_lower_95_a = round_probability(interval_lower_95_a),
      interval_median_50_a = round_probability(interval_median_50_a),
      interval_upper_95_a = round_probability(interval_upper_95_a),
      estimate_b = round_probability(estimate_b),
      interval_type_b,
      uncertainty_sd_b = round_probability(uncertainty_sd_b),
      interval_lower_95_b = round_probability(interval_lower_95_b),
      interval_median_50_b = round_probability(interval_median_50_b),
      interval_upper_95_b = round_probability(interval_upper_95_b),
      delta_interval_type = delta_interval_type_from_components(interval_type_a, interval_type_b),
      delta_probability = round_probability(estimate_b - estimate_a),
      delta_uncertainty_sd = round_probability(delta_uncertainty_sd),
      delta_interval_lower_95 = round_probability(delta_interval_lower_95),
      delta_interval_median_50 = round_probability(delta_interval_median_50),
      delta_interval_upper_95 = round_probability(delta_interval_upper_95),
      abs_delta_probability = round_probability(abs(estimate_b - estimate_a))
    ) %>%
    arrange(dataset_group, dataset_variant, probability_type, horizon_year, comparison)
}

standardize_stage3 <- function(path) {
  stage3_df <- read_csv_checked(path, "Stage 3 yearly estimate file")

  bind_rows(
    stage3_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "KM benchmark",
        model_scope = "overall",
        probability_type = "survival",
        horizon_year = as.integer(horizon_year),
        estimate_probability = clip_probability(estimated_survival_probability),
        interval_type = "95% CI",
        uncertainty_sd = NA_real_,
        interval_lower_95 = clip_probability(estimated_survival_lower_95),
        interval_median_50 = NA_real_,
        interval_upper_95 = clip_probability(estimated_survival_upper_95)
      ),
    stage3_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "KM benchmark",
        model_scope = "overall",
        probability_type = "risk",
        horizon_year = as.integer(horizon_year),
        estimate_probability = clip_probability(estimated_risk_probability),
        interval_type = "95% CI",
        uncertainty_sd = NA_real_,
        interval_lower_95 = clip_probability(estimated_risk_lower_95),
        interval_median_50 = NA_real_,
        interval_upper_95 = clip_probability(estimated_risk_upper_95)
      )
  ) %>%
    duplicate_dataset_rows(
      source_dataset_variant = "merged",
      target_dataset_variant = "merged_site_adjusted"
    )
}

standardize_stage5 <- function(path, interval_lookup) {
  stage5_df <- read_csv_checked(path, "Stage 5 yearly estimate file")

  bind_rows(
    stage5_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset_version_key),
        dataset_group = dataset_group_from_variant(dataset_version_key),
        model_family = "Frequentist no-cure lognormal",
        model_scope = "overall",
        probability_type = "survival",
        horizon_year = as.integer(horizon_year),
        estimate_probability = clip_probability(estimated_survival_probability)
      ),
    stage5_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset_version_key),
        dataset_group = dataset_group_from_variant(dataset_version_key),
        model_family = "Frequentist no-cure lognormal",
        model_scope = "overall",
        probability_type = "risk",
        horizon_year = as.integer(horizon_year),
        estimate_probability = clip_probability(estimated_risk_probability)
      )
  ) %>%
    left_join(
      interval_lookup,
      by = c(
        "dataset_variant",
        "dataset_group",
        "model_family",
        "model_scope",
        "probability_type",
        "horizon_year"
      )
    )
}

standardize_stage7 <- function(path, interval_lookup) {
  stage7_df <- read_csv_checked(path, "Stage 7 horizon summary file")

  bind_rows(
    stage7_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "Frequentist mixture cure",
        model_scope = "overall",
        probability_type = "survival",
        horizon_year = as.integer(time_horizon_year),
        estimate_probability = clip_probability(overall_survival_prob)
      ),
    stage7_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "Frequentist mixture cure",
        model_scope = "overall",
        probability_type = "risk",
        horizon_year = as.integer(time_horizon_year),
        estimate_probability = clip_probability(overall_risk_prob)
      ),
    stage7_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "Frequentist mixture cure",
        model_scope = "susceptible",
        probability_type = "survival",
        horizon_year = as.integer(time_horizon_year),
        estimate_probability = clip_probability(susceptible_only_survival_prob)
      ),
    stage7_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "Frequentist mixture cure",
        model_scope = "susceptible",
        probability_type = "risk",
        horizon_year = as.integer(time_horizon_year),
        estimate_probability = clip_probability(susceptible_only_risk_prob)
      )
  ) %>%
    left_join(
      interval_lookup,
      by = c(
        "dataset_variant",
        "dataset_group",
        "model_family",
        "model_scope",
        "probability_type",
        "horizon_year"
      )
    )
}

standardize_stage8 <- function(path) {
  stage8_df <- read_csv_checked(path, "Stage 8 horizon summary file")
  if ("prior_branch" %in% names(stage8_df)) {
    stage8_df <- stage8_df %>%
      filter(prior_branch == "anchor_informed")
  }

  bind_rows(
    stage8_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "Bayesian transition-only cure",
        model_scope = "overall",
        probability_type = "survival",
        horizon_year = as.integer(horizon_year),
        estimate_probability = clip_probability(overall_survival_prob),
        interval_type = "95% CrI",
        uncertainty_sd = as.numeric(overall_survival_prob_sd),
        interval_lower_95 = clip_probability(overall_survival_prob_q025),
        interval_median_50 = clip_probability(overall_survival_prob_q50),
        interval_upper_95 = clip_probability(overall_survival_prob_q975)
      ),
    stage8_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "Bayesian transition-only cure",
        model_scope = "overall",
        probability_type = "risk",
        horizon_year = as.integer(horizon_year),
        estimate_probability = clip_probability(overall_risk_prob),
        interval_type = "95% CrI",
        uncertainty_sd = as.numeric(overall_risk_prob_sd),
        interval_lower_95 = clip_probability(overall_risk_prob_q025),
        interval_median_50 = clip_probability(overall_risk_prob_q50),
        interval_upper_95 = clip_probability(overall_risk_prob_q975)
      ),
    stage8_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "Bayesian transition-only cure",
        model_scope = "susceptible",
        probability_type = "survival",
        horizon_year = as.integer(horizon_year),
        estimate_probability = clip_probability(susceptible_only_survival_prob),
        interval_type = "95% CrI",
        uncertainty_sd = as.numeric(susceptible_only_survival_prob_sd),
        interval_lower_95 = clip_probability(susceptible_only_survival_prob_q025),
        interval_median_50 = clip_probability(susceptible_only_survival_prob_q50),
        interval_upper_95 = clip_probability(susceptible_only_survival_prob_q975)
      ),
    stage8_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "Bayesian transition-only cure",
        model_scope = "susceptible",
        probability_type = "risk",
        horizon_year = as.integer(horizon_year),
        estimate_probability = clip_probability(susceptible_only_risk_prob),
        interval_type = "95% CrI",
        uncertainty_sd = as.numeric(susceptible_only_risk_prob_sd),
        interval_lower_95 = clip_probability(susceptible_only_risk_prob_q025),
        interval_median_50 = clip_probability(susceptible_only_risk_prob_q50),
        interval_upper_95 = clip_probability(susceptible_only_risk_prob_q975)
      )
  )
}

# Execute Integration -----------------------------------------------------
resolved_stage8_input_file <- resolve_existing_path(
  stage8_input_candidates,
  "Stage 8 Bayesian horizon summary file"
)

analysis_datasets <- build_analysis_datasets_from_source(
  source_file = source_data_file,
  pnu_label = pnu_site_label,
  snu_label = snu_site_label
)
stage5_interval_lookup <- compute_stage5_interval_lookup(
  fit_path = stage5_fit_file,
  analysis_datasets = analysis_datasets,
  target_horizons = required_horizons,
  n_draws = stage5_ci_draws
)
stage7_interval_lookup <- compute_stage7_interval_lookup(
  analysis_datasets = analysis_datasets,
  fit_path = stage7_fit_file,
  target_horizons = required_horizons,
  n_draws = stage7_ci_draws
)

metadata_df <- tibble::tibble(
  input_label = c(
    "source_data",
    "km_benchmark",
    "mle_nocure_summary",
    "mle_nocure_fit",
    "mle_mixture_cure_summary",
    "mle_mixture_cure_fit",
    "bayesian_mixture_cure"
  ),
  input_path = c(
    normalizePath(source_data_file, winslash = "/", mustWork = TRUE),
    normalizePath(stage3_input_file, winslash = "/", mustWork = TRUE),
    normalizePath(stage5_input_file, winslash = "/", mustWork = TRUE),
    normalizePath(stage5_fit_file, winslash = "/", mustWork = TRUE),
    normalizePath(stage7_input_file, winslash = "/", mustWork = TRUE),
    normalizePath(stage7_fit_file, winslash = "/", mustWork = TRUE),
    resolved_stage8_input_file
  )
)

integrated_long_df <- bind_rows(
  standardize_stage3(stage3_input_file),
  standardize_stage5(stage5_input_file, stage5_interval_lookup),
  standardize_stage7(stage7_input_file, stage7_interval_lookup),
  standardize_stage8(resolved_stage8_input_file)
) %>%
  filter(horizon_year %in% required_horizons) %>%
  arrange(
    factor(dataset_group, levels = c("PNU", "SNU", "Merged")),
    factor(dataset_variant, levels = c("PNU", "SNU", "merged", "merged_site_adjusted")),
    model_family,
    model_scope,
    probability_type,
    horizon_year
  )

validate_long_table(integrated_long_df)

probability_draw_lookup <- build_probability_draw_lookup(
  long_df = integrated_long_df,
  n_draws = delta_interval_draws
)

cleanup_existing_comparison_outputs(export_path)

readr::write_csv(
  metadata_df,
  file.path(export_path, "input_file_registry.csv")
)

readr::write_csv(
  integrated_long_df %>%
    mutate(
      estimate_probability = round_probability(estimate_probability),
      uncertainty_sd = round_probability(uncertainty_sd),
      interval_lower_95 = round_probability(interval_lower_95),
      interval_median_50 = round_probability(interval_median_50),
      interval_upper_95 = round_probability(interval_upper_95)
    ),
  file.path(export_path, "integrated_model_probabilities_long.csv")
)

comparison_specs <- tibble::tibble(
  dataset_group = c("PNU", "PNU", "SNU", "SNU", "Merged", "Merged"),
  probability_type = c("survival", "risk", "survival", "risk", "survival", "risk")
)

output_registry <- comparison_specs %>%
  rowwise() %>%
  mutate(
    output_file = write_comparison_table(
      table_df = build_comparison_table(integrated_long_df, dataset_group, probability_type),
      export_dir = export_path,
      dataset_group = dataset_group,
      probability_type = probability_type
    )
  ) %>%
  ungroup()

readr::write_csv(
  output_registry,
  file.path(export_path, "output_table_registry.csv")
)

readr::write_csv(
  build_probability_difference_summary(integrated_long_df, probability_draw_lookup),
  file.path(export_path, "model_probability_difference_summary_long.csv")
)

message("Comparison tables written to: ", export_path)
