# Configure: Block 3 paths and constants ===============================
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

dropbox_project_root <- Sys.getenv(
  "BLOCK3_DROPBOX_PROJECT_ROOT",
  unset = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model"
)
source_data_file <- Sys.getenv(
  "BLOCK3_SOURCE_DATA_FILE",
  unset = file.path(
    dropbox_project_root,
    "0.Data",
    "2.Preprocessed data",
    "Preprocessed_Merged_PNUH_SNUH_Data.csv"
  )
)
export_path <- Sys.getenv(
  "BLOCK3_EXPORT_PATH",
  unset = file.path(dropbox_project_root, "3.Block3")
)
system_rscript <- Sys.getenv(
  "BLOCK3_SYSTEM_RSCRIPT",
  unset = "/Library/Frameworks/R.framework/Resources/bin/Rscript"
)
frequentist_cure_script <- file.path(
  repo_root,
  "2.Rcode",
  "3.Block3",
  "2.Frequentist Cure Engine.r"
)
bayesian_cure_script <- file.path(
  repo_root,
  "2.Rcode",
  "3.Block3",
  "3.Bayesian Cure Engine.r"
)

pnu_site_label <- Sys.getenv("BLOCK3_PNU_SITE_LABEL", unset = "PNU")
snu_site_label <- Sys.getenv("BLOCK3_SNU_SITE_LABEL", unset = "SNU")

common_horizons_year <- 1:10
curve_time_grid_year <- seq(0, 10, by = 0.05)
time_origin_epsilon_year <- 1e-08
plot_width_in <- 14
plot_height_in <- 9
plot_dpi <- 320

bayesian_stan_chains <- as.integer(Sys.getenv("BLOCK3_BAYESIAN_STAN_CHAINS", unset = "2"))
bayesian_stan_iter <- as.integer(Sys.getenv("BLOCK3_BAYESIAN_STAN_ITER", unset = "1000"))
bayesian_stan_warmup <- as.integer(Sys.getenv("BLOCK3_BAYESIAN_STAN_WARMUP", unset = "500"))
bayesian_stan_thin <- as.integer(Sys.getenv("BLOCK3_BAYESIAN_STAN_THIN", unset = "1"))
bayesian_stan_adapt_delta <- as.numeric(Sys.getenv("BLOCK3_BAYESIAN_STAN_ADAPT_DELTA", unset = "0.95"))
bayesian_stan_max_treedepth <- as.integer(Sys.getenv("BLOCK3_BAYESIAN_STAN_MAX_TREEDEPTH", unset = "12"))
bayesian_posterior_draws <- as.integer(Sys.getenv("BLOCK3_BAYESIAN_POSTERIOR_PRED_DRAWS", unset = "200"))
no_cure_ci_draws <- as.integer(Sys.getenv("BLOCK3_NO_CURE_CI_DRAWS", unset = "1000"))

scenario_registry <- tibble::tibble(
  scenario_id = c("full", "truncated_2y", "truncated_3y"),
  scenario_label = c("Full follow-up", "Truncated at 2 years", "Truncated at 3 years"),
  truncation_year = c(NA_real_, 2, 3)
)

dataset_registry <- tibble::tibble(
  dataset_scope = c("SNU", "merged"),
  dataset_label = c("SNU", "Merged"),
  site_mode = c("single", "merged"),
  formula_rhs = c("age_s + sex_num", "age_s + sex_num + site"),
  frequentist_dataset_key = c("SNU", "merged_site_adjusted"),
  bayesian_dataset_key = c("SNU", "merged_site_adjusted")
)

no_cure_family_registry <- tibble::tibble(
  model_family = c("exponential", "weibull", "lognormal", "loglogistic"),
  flexsurv_dist = c("exp", "weibull", "lnorm", "llogis")
)

frequentist_required_files <- c(
  "block3_frequentist_cure_horizon_summary.csv",
  "block3_frequentist_cure_fit_summary.csv",
  "block3_frequentist_cure_plot_source.csv"
)
bayesian_required_files <- c(
  "block3_bayesian_cure_horizon_summary.csv",
  "block3_bayesian_cure_posterior_draw_summary.csv",
  "block3_bayesian_cure_plot_source.csv"
)
frequentist_interpretation_files <- c(
  "block3_frequentist_cure_horizon_summary.csv",
  "block3_frequentist_cure_fit_summary.csv",
  "block3_frequentist_cure_fitted_objects.rds",
  "block3_frequentist_cure_plot_source.csv",
  "block3_frequentist_cure_detail_annotation_table.csv",
  "block3_frequentist_cure_latency_plot_source.csv"
)
bayesian_interpretation_files <- c(
  "block3_bayesian_cure_horizon_summary.csv",
  "block3_bayesian_cure_prior_branch_horizon_summary.csv",
  "block3_bayesian_cure_anchor_vs_neutral_delta_summary.csv",
  "block3_bayesian_cure_model_registry.csv",
  "block3_bayesian_cure_posterior_draw_summary.csv",
  "block3_bayesian_cure_plot_source.csv",
  "block3_bayesian_cure_detail_annotation_table.csv",
  "block3_bayesian_cure_latency_plot_source.csv",
  "block3_bayesian_cure_parameter_summary.csv",
  "block3_bayesian_cure_uncured_latency_summary.csv",
  "block3_bayesian_cure_latency_aft_effect_summary.csv",
  file.path("sub", "block3_bayesian_cure_trace_plots.pdf")
)
frequentist_interpretation_priority_files <- c(
  "block3_frequentist_cure_fit_summary.csv",
  "block3_frequentist_cure_detail_annotation_table.csv"
)
bayesian_interpretation_priority_files <- c(
  "block3_bayesian_cure_parameter_summary.csv",
  "block3_bayesian_cure_posterior_draw_summary.csv",
  "block3_bayesian_cure_model_registry.csv",
  "block3_bayesian_cure_anchor_vs_neutral_delta_summary.csv"
)
interpretation_priority_files <- c(
  frequentist_interpretation_priority_files,
  bayesian_interpretation_priority_files
)
frequentist_interpretation_subdir_files <- setdiff(
  frequentist_interpretation_files,
  frequentist_interpretation_priority_files
)
bayesian_interpretation_subdir_files <- setdiff(
  bayesian_interpretation_files,
  bayesian_interpretation_priority_files
)

interpretation_dir <- file.path(export_path, "1.Interpretation Diagnostics")
interpretation_frequentist_dir <- file.path(interpretation_dir, "Frequentist Cure")
interpretation_bayesian_dir <- file.path(interpretation_dir, "Bayesian Cure")
interpretation_readme_file <- file.path(interpretation_dir, "README.md")
supporting_table_dir <- file.path(export_path, "2.Supporting Tables")
technical_dir <- file.path(export_path, "3.Technical")
scenario_source_dir <- file.path(technical_dir, "scenario_sources")
supporting_model_dir <- file.path(technical_dir, "supporting_model_exports")
log_dir <- file.path(technical_dir, "logs")

km_horizon_file <- file.path(supporting_table_dir, "block3_km_horizon_summary.csv")
km_curve_file <- file.path(supporting_table_dir, "block3_km_curve_data.csv")
nocure_registry_file <- file.path(supporting_table_dir, "block3_nocure_model_registry.csv")
nocure_horizon_file <- file.path(supporting_table_dir, "block3_nocure_horizon_summary.csv")
nocure_curve_file <- file.path(supporting_table_dir, "block3_nocure_curve_data.csv")
freq_cure_horizon_file <- file.path(supporting_table_dir, "block3_frequentist_cure_horizon_summary.csv")
freq_cure_curve_file <- file.path(supporting_table_dir, "block3_frequentist_cure_curve_data.csv")
bayes_horizon_file <- file.path(supporting_table_dir, "block3_bayesian_cure_horizon_summary.csv")
bayes_curve_file <- file.path(supporting_table_dir, "block3_bayesian_cure_curve_data.csv")
probability_table_file <- file.path(export_path, "block3_probability_comparison_long.csv")
curve_table_file <- file.path(supporting_table_dir, "block3_model_curve_comparison_long.csv")
divergence_long_file <- file.path(supporting_table_dir, "block3_divergence_long.csv")
divergence_summary_file <- file.path(export_path, "block3_divergence_summary.csv")
cure_signal_file <- file.path(export_path, "block3_cure_signal_summary.csv")
standard_error_registry_file <- file.path(supporting_table_dir, "block3_standard_error_table_registry.csv")
standard_error_long_file <- file.path(supporting_table_dir, "block3_standard_error_long.csv")
scenario_manifest_file <- file.path(technical_dir, "block3_scenario_manifest.csv")
dataset_summary_file <- file.path(export_path, "block3_dataset_summary.csv")
metadata_registry_file <- file.path(technical_dir, "block3_metadata_registry.csv")
run_log_file <- file.path(technical_dir, "block3_run_log.csv")
summary_markdown_file <- file.path(export_path, "block3_summary.md")
manifest_file <- file.path(technical_dir, "block3_export_manifest.csv")
km_plot_file <- file.path(export_path, "block3_plot_km_risk_curves.png")
model_plot_file <- file.path(export_path, "block3_plot_model_risk_curves.png")

# Attach: packages and runtime state ===============================
required_packages <- c("readr", "dplyr", "tibble", "purrr", "tidyr", "ggplot2", "survival", "flexsurv", "scales")
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
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(survival)
  library(flexsurv)
  library(scales)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)
dir.create(interpretation_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(interpretation_frequentist_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(interpretation_bayesian_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supporting_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(technical_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(scenario_source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supporting_model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# Define: general helpers ===============================
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

assert_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

normalize_existing_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

remove_legacy_block3_paths <- function() {
  legacy_interpretation_dirs <- c(
    paste0("3", ".MLE Mixture Cure"),
    paste0("4", ".Bayesian Mixture Cure")
  )
  legacy_paths <- c(
    file.path(technical_dir, "supporting_stage_exports"),
    file.path(interpretation_dir, legacy_interpretation_dirs)
  )
  legacy_log_files <- Sys.glob(file.path(log_dir, "block3__*__stage*.log"))
  existing_paths <- legacy_paths[file.exists(legacy_paths) | dir.exists(legacy_paths)]

  if (length(existing_paths) > 0L) {
    unlink(existing_paths, recursive = TRUE, force = TRUE)
  }
  if (length(legacy_log_files) > 0L) {
    unlink(legacy_log_files, force = TRUE)
  }

  invisible(c(existing_paths, legacy_log_files))
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

clip_prob <- function(x, lower = 0, upper = 1) {
  pmin(pmax(as.numeric(x), lower), upper)
}

read_csv_text <- function(path, label) {
  assert_file_exists(path, label)
  readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
}

new_run_log <- function() {
  tibble(
    component = character(),
    status = character(),
    detail = character(),
    file_path = character()
  )
}

append_run_log <- function(run_log, component, status, detail = NA_character_, file_path = NA_character_) {
  bind_rows(
    run_log,
    tibble(
      component = component,
      status = status,
      detail = detail,
      file_path = file_path
    )
  )
}

make_scenario_source_file <- function(scenario_id) {
  file.path(scenario_source_dir, paste0("block3_source__", scenario_id, ".csv"))
}

make_model_target_dir <- function(scenario_id, model_folder_name) {
  file.path(supporting_model_dir, scenario_id, model_folder_name)
}

summarize_dataset_counts <- function(df, dataset_scope, scenario_id) {
  tibble(
    scenario_id = scenario_id,
    dataset_scope = dataset_scope,
    n_subjects = nrow(df),
    n_transition = sum(df$status_num == 1L, na.rm = TRUE),
    n_censored = sum(df$status_num == 0L, na.rm = TRUE),
    n_remission = sum(df$status_num == 2L, na.rm = TRUE),
    max_followup_year = max(df$days_followup, na.rm = TRUE) / 365.25,
    median_followup_year = stats::median(df$days_followup, na.rm = TRUE) / 365.25
  )
}

make_manifest_tbl <- function(file_paths, object_names, descriptions) {
  tibble(
    file_name = basename(file_paths),
    object_name = object_names,
    description = descriptions,
    file_path = normalize_existing_path(file_paths),
    portable_file_reference = basename(file_paths)
  )
}

make_env_assignment <- function(key, value) {
  paste0(key, "=", shQuote(as.character(value)))
}

standard_error_column_patterns <- c(
  "(^|_)sd$",
  "(^|_)se$",
  "std_error$",
  "std\\.error$",
  "stderr$",
  "posterior_sd$",
  "uncertainty_sd$",
  "robust_se$"
)

identify_standard_error_columns <- function(df) {
  col_names <- names(df)
  if (is.null(col_names) || length(col_names) == 0L) {
    return(character())
  }

  col_names[vapply(
    col_names,
    function(col_name) any(grepl(standard_error_column_patterns, col_name, ignore.case = TRUE)),
    logical(1)
  )]
}

identify_standard_error_id_columns <- function(df) {
  preferred_cols <- c(
    "scenario_id",
    "scenario_label",
    "dataset",
    "dataset_label",
    "source_dataset",
    "dataset_scope",
    "model_family",
    "model_name",
    "model_type",
    "horizon_year",
    "time_year",
    "curve_source",
    "comparison_type",
    "contrast_label",
    "site_counterfactual",
    "parameter",
    "term"
  )

  intersect(preferred_cols, names(df))
}

empty_standard_error_registry <- function() {
  tibble(
    source_object = character(),
    n_rows = integer(),
    n_columns = integer(),
    n_standard_error_columns = integer(),
    standard_error_columns = character()
  )
}

empty_standard_error_long <- function() {
  tibble(
    source_object = character(),
    row_id = integer(),
    standard_error_column = character(),
    standard_error_value = numeric(),
    standard_error_value_raw = character()
  )
}

build_standard_error_registry_entry <- function(df, source_name) {
  se_cols <- identify_standard_error_columns(df)

  tibble(
    source_object = source_name,
    n_rows = nrow(df),
    n_columns = ncol(df),
    n_standard_error_columns = length(se_cols),
    standard_error_columns = if (length(se_cols) > 0L) paste(se_cols, collapse = "|") else NA_character_
  )
}

build_standard_error_long_table <- function(df, source_name) {
  se_cols <- identify_standard_error_columns(df)
  if (length(se_cols) == 0L || nrow(df) == 0L) {
    return(empty_standard_error_long())
  }

  base_df <- tibble::as_tibble(df) %>%
    mutate(row_id = dplyr::row_number())
  id_cols <- identify_standard_error_id_columns(base_df)

  bind_rows(lapply(se_cols, function(se_col) {
    out <- tibble(
      source_object = source_name,
      row_id = base_df$row_id,
      standard_error_column = se_col,
      standard_error_value = suppressWarnings(as.numeric(base_df[[se_col]])),
      standard_error_value_raw = as.character(base_df[[se_col]])
    )

    if (length(id_cols) > 0L) {
      out <- bind_cols(out, base_df[, id_cols, drop = FALSE])
    }

    out
  }))
}

build_standard_error_export_bundle <- function(named_tables) {
  named_tables <- named_tables[vapply(named_tables, function(x) inherits(x, "data.frame"), logical(1))]

  registry_tbl <- bind_rows(lapply(names(named_tables), function(table_name) {
    build_standard_error_registry_entry(named_tables[[table_name]], table_name)
  }))
  long_tbl <- bind_rows(lapply(names(named_tables), function(table_name) {
    build_standard_error_long_table(named_tables[[table_name]], table_name)
  }))

  if (nrow(registry_tbl) == 0L) {
    registry_tbl <- empty_standard_error_registry()
  }
  if (nrow(long_tbl) == 0L) {
    long_tbl <- empty_standard_error_long()
  }

  list(
    registry = registry_tbl,
    long = long_tbl
  )
}

# Define: source-data preparation ===============================
standardize_known_site_labels <- function(df, pnu_label, snu_label) {
  if (!"site" %in% names(df)) {
    stop("Merged input must contain column `site`.", call. = FALSE)
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

read_source_data <- function(path, pnu_label, snu_label) {
  df <- read_csv_text(path, "Block 3 merged preprocessed data CSV")
  df <- standardize_known_site_labels(df, pnu_label, snu_label)
  validate_source_data(df, pnu_label, snu_label)
  df
}

apply_followup_truncation <- function(source_df, truncation_year = NA_real_) {
  out <- source_df %>%
    mutate(
      days_followup = safe_numeric(days_followup),
      status_num = as.integer(safe_numeric(status_num))
    )

  if (!is.finite(truncation_year)) {
    return(out)
  }

  cutoff_days <- truncation_year * 365.25

  out %>%
    mutate(
      truncated_flag = as.integer(days_followup > cutoff_days),
      days_followup = if_else(days_followup > cutoff_days, cutoff_days, days_followup),
      status_num = if_else(truncated_flag == 1L, 0L, status_num),
      status = if ("status" %in% names(.)) {
        if_else(truncated_flag == 1L, sprintf("censored_at_%sy", truncation_year), as.character(status))
      } else {
        NULL
      }
    ) %>%
    select(-tidyselect::any_of("truncated_flag"))
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

prepare_backbone_dataset <- function(df, dataset_name, site_mode = c("single", "merged")) {
  site_mode <- match.arg(site_mode)

  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("[%s] Missing required columns: %s", dataset_name, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  out <- df %>%
    mutate(
      id = trimws(as.character(id)),
      site = trimws(as.character(site)),
      sex_num = as.integer(safe_numeric(sex_num)),
      age_exact_entry = safe_numeric(age_exact_entry),
      days_followup = safe_numeric(days_followup),
      status_num = as.integer(safe_numeric(status_num))
    )

  if (nrow(out) == 0L) {
    stop(sprintf("[%s] Dataset has zero rows after loading.", dataset_name), call. = FALSE)
  }

  if (anyNA(out[, required_cols])) {
    stop(sprintf("[%s] Missing values detected in required columns.", dataset_name), call. = FALSE)
  }

  if (any(out$id == "", na.rm = TRUE)) {
    stop(sprintf("[%s] Blank `id` values detected.", dataset_name), call. = FALSE)
  }

  if (any(!out$sex_num %in% c(0L, 1L), na.rm = TRUE)) {
    stop(sprintf("[%s] `sex_num` must be coded as 0/1 only.", dataset_name), call. = FALSE)
  }

  if (any(!out$status_num %in% c(0L, 1L, 2L), na.rm = TRUE)) {
    stop(sprintf("[%s] `status_num` must be coded as 0/1/2 only.", dataset_name), call. = FALSE)
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

  age_sd <- stats::sd(out$age_exact_entry)
  if (is.na(age_sd) || age_sd <= 0) {
    stop(sprintf("[%s] `age_exact_entry` must have positive standard deviation.", dataset_name), call. = FALSE)
  }

  age_mean <- mean(out$age_exact_entry)

  out %>%
    mutate(
      unique_person_id = paste(site, id, sep = "_"),
      time_year = days_followup / 365.25,
      event_main = as.integer(status_num == 1L),
      censor_main = as.integer(status_num %in% c(0L, 2L)),
      age_s = (age_exact_entry - age_mean) / (2 * age_sd)
    )
}

build_analysis_datasets_from_source_df <- function(source_df, pnu_label, snu_label) {
  list(
    PNU = prepare_backbone_dataset(
      split_site_dataset(source_df, pnu_label, "PNU"),
      dataset_name = "PNU",
      site_mode = "single"
    ),
    SNU = prepare_backbone_dataset(
      split_site_dataset(source_df, snu_label, "SNU"),
      dataset_name = "SNU",
      site_mode = "single"
    ),
    merged = prepare_backbone_dataset(
      source_df,
      dataset_name = "merged",
      site_mode = "merged"
    )
  )
}

# Define: KM and no-cure helpers ===============================
extract_flexsurv_matrix <- function(summary_list, target_times) {
  if (!is.list(summary_list)) {
    summary_list <- list(summary_list)
  }

  t(
    vapply(
      summary_list,
      FUN.VALUE = numeric(length(target_times)),
      FUN = function(x) {
        xx <- as.data.frame(x)
        if ("time" %in% names(xx) && nrow(xx) != length(target_times)) {
          matched_idx <- vapply(
            target_times,
            FUN.VALUE = integer(1),
            FUN = function(tt) which.min(abs(xx$time - tt))
          )
          xx <- xx[matched_idx, , drop = FALSE]
        }
        if ("est" %in% names(xx)) {
          return(as.numeric(xx$est))
        }
        as.numeric(xx[[ncol(xx)]])
      }
    )
  )
}

predict_parametric_survival <- function(fit, newdata, target_times) {
  surv_list <- summary(
    fit,
    newdata = newdata,
    t = target_times,
    type = "survival",
    ci = FALSE,
    se = FALSE
  )

  extract_flexsurv_matrix(surv_list, target_times)
}

summarize_interval_draws <- function(x) {
  tibble(
    uncertainty_sd = stats::sd(x, na.rm = TRUE),
    lower_95 = unname(stats::quantile(x, probs = 0.025, na.rm = TRUE, names = FALSE)),
    median_50 = unname(stats::quantile(x, probs = 0.5, na.rm = TRUE, names = FALSE)),
    upper_95 = unname(stats::quantile(x, probs = 0.975, na.rm = TRUE, names = FALSE))
  )
}

extract_flexsurv_model_formula <- function(fit) {
  model_formula <- tryCatch(
    stats::formula(fit),
    error = function(e) NULL
  )

  if (is.null(model_formula) && !is.null(fit$concat.formula) && inherits(fit$concat.formula, "formula")) {
    model_formula <- fit$concat.formula
  }

  if (is.null(model_formula) && !is.null(fit$call$formula) && inherits(fit$call$formula, "formula")) {
    model_formula <- fit$call$formula
  }

  if (is.null(model_formula) && !is.null(fit$all.formulae)) {
    candidate_names <- c("meanlog", "scale", "rate", "location")
    available_candidates <- intersect(candidate_names, names(fit$all.formulae))
    if (length(available_candidates) > 0L) {
      model_formula <- fit$all.formulae[[available_candidates[[1L]]]]
    }
  }

  if (is.null(model_formula)) {
    stop("No-cure fit object is missing the model formula needed for interval calculation.", call. = FALSE)
  }

  model_formula
}

build_nocure_interval_summaries <- function(fit, df, target_times, n_draws) {
  model_formula <- extract_flexsurv_model_formula(fit)
  design_matrix <- stats::model.matrix(
    stats::delete.response(stats::terms(model_formula)),
    data = df
  )
  model_pars <- fit$dlist$pars
  summary_fn <- get("summary_fns", asNamespace("flexsurv"))(fit, "survival")

  sim_tbl <- flexsurv::normboot.flexsurvreg(
    fit,
    B = n_draws,
    X = design_matrix,
    tidy = TRUE
  ) %>%
    tibble::as_tibble()

  missing_pars <- setdiff(model_pars, names(sim_tbl))
  if (length(missing_pars) > 0L) {
    stop(
      sprintf(
        "No-cure interval draws are missing model parameters: %s",
        paste(missing_pars, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  bind_rows(lapply(target_times, function(tt) {
    time_value <- as.numeric(tt)

    if (time_value <= 0) {
      draw_tbl <- sim_tbl %>%
        distinct(repno) %>%
        transmute(
          repno = as.integer(repno),
          estimated_survival_probability = 1,
          estimated_risk_probability = 0
        )
    } else {
      fn_args <- c(
        list(
          t = rep(time_value, nrow(sim_tbl)),
          start = rep(0, nrow(sim_tbl))
        ),
        as.list(sim_tbl[model_pars]),
        fit$aux
      )

      survival_by_subject_draw <- do.call(summary_fn, fn_args)

      draw_tbl <- tibble(
        repno = as.integer(sim_tbl$repno),
        estimated_survival_probability = clip_prob(survival_by_subject_draw)
      ) %>%
        group_by(repno) %>%
        summarise(
          estimated_survival_probability = mean(estimated_survival_probability, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          estimated_risk_probability = clip_prob(1 - estimated_survival_probability)
        )
    }

    summarize_interval_draws(draw_tbl$estimated_survival_probability) %>%
      rename_with(~ paste0("estimated_survival_", .x), everything()) %>%
      bind_cols(
        summarize_interval_draws(draw_tbl$estimated_risk_probability) %>%
          rename_with(~ paste0("estimated_risk_", .x), everything())
      ) %>%
      mutate(time_year = time_value)
  }))
}

ensure_numeric_columns <- function(df, cols) {
  for (one_col in cols) {
    if (!(one_col %in% names(df))) {
      df[[one_col]] <- NA_real_
    }
    df[[one_col]] <- as.numeric(df[[one_col]])
  }

  df
}

build_km_exports <- function(df, scenario_id, scenario_label, dataset_scope, dataset_label, time_grid, horizon_years) {
  fit <- survival::survfit(
    survival::Surv(time_year, event_main) ~ 1,
    data = df,
    conf.type = "log-log"
  )

  horizon_summary <- summary(fit, times = horizon_years, extend = TRUE)
  curve_summary <- summary(fit, times = time_grid, extend = TRUE)

  km_horizon_tbl <- tibble(
    scenario_id = scenario_id,
    scenario_label = scenario_label,
    dataset_scope = dataset_scope,
    dataset_label = dataset_label,
    model_group = "KM",
    model_name = "km",
    model_label = "KM",
    interval_type = "95% CI",
    horizon_year = as.integer(horizon_years),
    estimated_survival_probability = clip_prob(horizon_summary$surv),
    estimated_risk_probability = clip_prob(1 - horizon_summary$surv),
    estimated_survival_uncertainty_sd = NA_real_,
    estimated_survival_lower_95 = clip_prob(horizon_summary$lower),
    estimated_survival_median_50 = clip_prob(horizon_summary$surv),
    estimated_survival_upper_95 = clip_prob(horizon_summary$upper),
    estimated_risk_uncertainty_sd = NA_real_,
    estimated_risk_lower_95 = clip_prob(1 - horizon_summary$upper),
    estimated_risk_median_50 = clip_prob(1 - horizon_summary$surv),
    estimated_risk_upper_95 = clip_prob(1 - horizon_summary$lower),
    lower_survival_probability = clip_prob(horizon_summary$lower),
    upper_survival_probability = clip_prob(horizon_summary$upper),
    n_at_risk = vapply(
      horizon_years,
      FUN.VALUE = integer(1),
      FUN = function(tt) sum(df$time_year >= tt, na.rm = TRUE)
    ),
    n_subjects = nrow(df),
    n_transition_events_total = sum(df$event_main == 1L, na.rm = TRUE),
    cure_fraction = NA_real_,
    susceptible_fraction = NA_real_
  )

  km_curve_tbl <- tibble(
    scenario_id = scenario_id,
    scenario_label = scenario_label,
    dataset_scope = dataset_scope,
    dataset_label = dataset_label,
    model_group = "KM",
    model_name = "km",
    model_label = "KM",
    interval_type = "95% CI",
    time_year = as.numeric(time_grid),
    estimated_survival_probability = clip_prob(curve_summary$surv),
    estimated_risk_probability = clip_prob(1 - curve_summary$surv),
    estimated_survival_uncertainty_sd = NA_real_,
    estimated_survival_lower_95 = clip_prob(curve_summary$lower),
    estimated_survival_median_50 = clip_prob(curve_summary$surv),
    estimated_survival_upper_95 = clip_prob(curve_summary$upper),
    estimated_risk_uncertainty_sd = NA_real_,
    estimated_risk_lower_95 = clip_prob(1 - curve_summary$upper),
    estimated_risk_median_50 = clip_prob(1 - curve_summary$surv),
    estimated_risk_upper_95 = clip_prob(1 - curve_summary$lower),
    lower_survival_probability = clip_prob(curve_summary$lower),
    upper_survival_probability = clip_prob(curve_summary$upper)
  )

  list(horizon = km_horizon_tbl, curve = km_curve_tbl)
}

fit_single_nocure_family <- function(df, formula_rhs, flexsurv_dist) {
  survival_formula <- stats::as.formula(
    paste0("survival::Surv(time_year_model, event_main) ~ ", formula_rhs)
  )

  tryCatch(
    flexsurv::flexsurvreg(
      formula = survival_formula,
      data = df,
      dist = flexsurv_dist
    ),
    error = function(e) e
  )
}

extract_nocure_fit_stats <- function(fit_obj, df, scenario_id, dataset_scope, model_family, formula_rhs) {
  if (inherits(fit_obj, "error")) {
    return(
      tibble(
        scenario_id = scenario_id,
        dataset_scope = dataset_scope,
        model_family = model_family,
        formula_rhs = formula_rhs,
        converged = FALSE,
        error_message = conditionMessage(fit_obj),
        n_subjects = nrow(df),
        n_transition_events_total = sum(df$event_main == 1L, na.rm = TRUE),
        n_parameters = NA_integer_,
        logLik = NA_real_,
        AIC = NA_real_,
        BIC = NA_real_
      )
    )
  }

  n_parameters <- length(stats::coef(fit_obj))
  loglik_value <- tryCatch(as.numeric(stats::logLik(fit_obj))[1], error = function(e) NA_real_)
  aic_value <- tryCatch(stats::AIC(fit_obj), error = function(e) NA_real_)
  bic_value <- if (!is.na(loglik_value)) {
    -2 * loglik_value + log(nrow(df)) * n_parameters
  } else {
    NA_real_
  }

  tibble(
    scenario_id = scenario_id,
    dataset_scope = dataset_scope,
    model_family = model_family,
    formula_rhs = formula_rhs,
    converged = TRUE,
    error_message = NA_character_,
    n_subjects = nrow(df),
    n_transition_events_total = sum(df$event_main == 1L, na.rm = TRUE),
    n_parameters = as.integer(n_parameters),
    logLik = loglik_value,
    AIC = aic_value,
    BIC = bic_value
  )
}

build_nocure_exports <- function(df, scenario_id, scenario_label, dataset_scope, dataset_label, formula_rhs, family_row, time_grid, horizon_years) {
  model_family <- family_row$model_family[[1]]
  flexsurv_dist <- family_row$flexsurv_dist[[1]]

  model_input <- df %>%
    mutate(
      site = factor(as.character(site)),
      time_year_model = pmax(as.numeric(time_year), time_origin_epsilon_year),
      event_main = as.integer(event_main)
    )

  fit_obj <- fit_single_nocure_family(model_input, formula_rhs, flexsurv_dist)
  fit_registry_tbl <- extract_nocure_fit_stats(fit_obj, model_input, scenario_id, dataset_scope, model_family, formula_rhs)

  if (inherits(fit_obj, "error")) {
    return(
      list(
        registry = fit_registry_tbl,
        horizon = tibble(),
        curve = tibble()
      )
    )
  }

  yearly_survival_mat <- predict_parametric_survival(
    fit = fit_obj,
    newdata = model_input,
    target_times = horizon_years
  )
  curve_survival_mat <- predict_parametric_survival(
    fit = fit_obj,
    newdata = model_input,
    target_times = time_grid
  )

  yearly_survival <- clip_prob(colMeans(yearly_survival_mat, na.rm = TRUE))
  curve_survival <- clip_prob(colMeans(curve_survival_mat, na.rm = TRUE))
  yearly_interval_tbl <- build_nocure_interval_summaries(
    fit = fit_obj,
    df = model_input,
    target_times = horizon_years,
    n_draws = no_cure_ci_draws
  ) %>%
    transmute(
      horizon_year = as.integer(round(time_year)),
      estimated_survival_uncertainty_sd,
      estimated_survival_lower_95,
      estimated_survival_median_50,
      estimated_survival_upper_95,
      estimated_risk_uncertainty_sd,
      estimated_risk_lower_95,
      estimated_risk_median_50,
      estimated_risk_upper_95
    )
  curve_interval_tbl <- build_nocure_interval_summaries(
    fit = fit_obj,
    df = model_input,
    target_times = time_grid,
    n_draws = no_cure_ci_draws
  ) %>%
    transmute(
      time_year = as.numeric(time_year),
      estimated_survival_uncertainty_sd,
      estimated_survival_lower_95,
      estimated_survival_median_50,
      estimated_survival_upper_95,
      estimated_risk_uncertainty_sd,
      estimated_risk_lower_95,
      estimated_risk_median_50,
      estimated_risk_upper_95
    )

  model_name <- paste0("no_cure_", model_family)
  model_label <- paste("No-cure", model_family)

  horizon_tbl <- tibble(
    scenario_id = scenario_id,
    scenario_label = scenario_label,
    dataset_scope = dataset_scope,
    dataset_label = dataset_label,
    model_group = "No-cure",
    model_name = model_name,
    model_label = model_label,
    interval_type = "95% CI",
    model_family = model_family,
    formula_rhs = formula_rhs,
    horizon_year = as.integer(horizon_years),
    estimated_survival_probability = yearly_survival,
    estimated_risk_probability = clip_prob(1 - yearly_survival),
    n_subjects = nrow(model_input),
    n_transition_events_total = sum(model_input$event_main == 1L, na.rm = TRUE),
    cure_fraction = NA_real_,
    susceptible_fraction = NA_real_
  ) %>%
    left_join(yearly_interval_tbl, by = "horizon_year")

  curve_tbl <- tibble(
    scenario_id = scenario_id,
    scenario_label = scenario_label,
    dataset_scope = dataset_scope,
    dataset_label = dataset_label,
    model_group = "No-cure",
    model_name = model_name,
    model_label = model_label,
    interval_type = "95% CI",
    model_family = model_family,
    formula_rhs = formula_rhs,
    time_year = as.numeric(time_grid),
    estimated_survival_probability = curve_survival,
    estimated_risk_probability = clip_prob(1 - curve_survival)
  ) %>%
    left_join(curve_interval_tbl, by = "time_year")

  list(
    registry = fit_registry_tbl,
    horizon = horizon_tbl,
    curve = curve_tbl
  )
}

# Define: model-script helpers ===============================
ensure_required_paths <- function(dir_path, relative_paths, label) {
  missing_paths <- relative_paths[!file.exists(file.path(dir_path, relative_paths))]
  if (length(missing_paths) > 0L) {
    stop(
      sprintf(
        "%s is missing required files in `%s`: %s",
        label,
        dir_path,
        paste(missing_paths, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

ensure_required_files <- function(dir_path, file_names, label) {
  ensure_required_paths(dir_path, file_names, label)
}

copy_selected_files <- function(source_dir, target_dir, relative_paths) {
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  missing_source_paths <- relative_paths[!file.exists(file.path(source_dir, relative_paths))]
  if (length(missing_source_paths) > 0L) {
    stop(
      sprintf(
        "Cannot copy stage outputs because source files are missing in `%s`: %s",
        source_dir,
        paste(missing_source_paths, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  purrr::walk(
    relative_paths,
    function(relative_path) {
      target_path <- file.path(target_dir, relative_path)
      dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
      file.copy(
        from = file.path(source_dir, relative_path),
        to = target_path,
        overwrite = TRUE
      )
    }
  )
}

run_model_script <- function(script_path, export_dir, env_pairs, log_file) {
  dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

  command_string <- paste(
    paste(env_pairs, collapse = " "),
    shQuote(system_rscript),
    shQuote(script_path),
    ">",
    shQuote(log_file),
    "2>&1"
  )

  status_code <- system(
    command = command_string,
    intern = FALSE,
    wait = TRUE
  )

  if (!identical(status_code, 0L)) {
    stop(
      sprintf(
        "Model script failed with status %s. Inspect log: %s",
        as.character(status_code),
        log_file
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

write_interpretation_readme <- function(path) {
  lines <- c(
    "# Block 3 Interpretation Diagnostics",
    "",
    "This folder bundles the full follow-up Block 3 frequentist and Bayesian cure files that are most useful for interpreting why cure-fraction estimates diverge, especially for the SNU full-follow-up comparison.",
    "",
    "Suggested reading order",
    "1. Start with the six key interpretation CSVs placed directly in this folder.",
    "2. Use the Frequentist Cure/ and Bayesian Cure/ subfolders only for supplemental full-follow-up files.",
    "3. Open the fitted-object RDS or trace plots only if deeper inspection is needed.",
    "3. Bayesian prior-branch and anchor-vs-neutral summaries for prior anchoring sensitivity.",
    "4. Bayesian trace plots PDF can be used for chain mixing and sampler behavior checks.",
    "",
    "Folder contents",
    "- Root-level CSVs: key interpretation files needed for direct review and upload.",
    "- Frequentist Cure/: supplemental full-follow-up frequentist mixture-cure files from Block 3.",
    "- Bayesian Cure/: supplemental full-follow-up Bayesian mixture-cure files and diagnostics from Block 3."
  )

  writeLines(lines, con = path)
}

ensure_interpretation_bundle <- function(full_frequentist_dir, full_bayesian_dir) {
  dir.create(interpretation_dir, recursive = TRUE, showWarnings = FALSE)
  copy_selected_files(full_frequentist_dir, interpretation_dir, frequentist_interpretation_priority_files)
  copy_selected_files(full_bayesian_dir, interpretation_dir, bayesian_interpretation_priority_files)
  copy_selected_files(full_frequentist_dir, interpretation_frequentist_dir, frequentist_interpretation_subdir_files)
  copy_selected_files(full_bayesian_dir, interpretation_bayesian_dir, bayesian_interpretation_subdir_files)

  unlink(file.path(interpretation_frequentist_dir, frequentist_interpretation_priority_files))
  unlink(file.path(interpretation_bayesian_dir, bayesian_interpretation_priority_files))

  ensure_required_paths(
    interpretation_dir,
    interpretation_priority_files,
    "Block 3 interpretation priority bundle"
  )
  ensure_required_paths(
    interpretation_frequentist_dir,
    frequentist_interpretation_subdir_files,
    "Block 3 interpretation frequentist supplemental bundle"
  )
  ensure_required_paths(
    interpretation_bayesian_dir,
    bayesian_interpretation_subdir_files,
    "Block 3 interpretation Bayesian supplemental bundle"
  )

  write_interpretation_readme(interpretation_readme_file)

  list(
    readme = interpretation_readme_file,
    interpretation_dir = interpretation_dir,
    frequentist_dir = interpretation_frequentist_dir,
    bayesian_dir = interpretation_bayesian_dir
  )
}

ensure_scenario_model_exports <- function(scenario_id, scenario_source_file) {
  frequentist_target_dir <- make_model_target_dir(scenario_id, "frequentist_cure")
  bayesian_target_dir <- make_model_target_dir(scenario_id, "bayesian_cure")

  frequentist_log_file <- file.path(log_dir, paste0("block3__", scenario_id, "__frequentist_cure.log"))
  bayesian_log_file <- file.path(log_dir, paste0("block3__", scenario_id, "__bayesian_cure.log"))

  if (!all(file.exists(file.path(frequentist_target_dir, frequentist_required_files)))) {
    run_model_script(
      script_path = frequentist_cure_script,
      export_dir = frequentist_target_dir,
      env_pairs = c(
        make_env_assignment("BLOCK3_FREQUENTIST_DATA_FILE", scenario_source_file),
        make_env_assignment("BLOCK3_FREQUENTIST_EXPORT_PATH", frequentist_target_dir),
        make_env_assignment("BLOCK3_FREQUENTIST_REUSE_EXISTING_FIT_RDS", "TRUE"),
        make_env_assignment("BLOCK3_PNU_SITE_LABEL", pnu_site_label),
        make_env_assignment("BLOCK3_SNU_SITE_LABEL", snu_site_label)
      ),
      log_file = frequentist_log_file
    )
  }

  if (!all(file.exists(file.path(bayesian_target_dir, bayesian_required_files)))) {
    run_model_script(
      script_path = bayesian_cure_script,
      export_dir = bayesian_target_dir,
      env_pairs = c(
        make_env_assignment("BLOCK3_BAYESIAN_DATA_FILE", scenario_source_file),
        make_env_assignment("BLOCK3_BAYESIAN_EXPORT_PATH", bayesian_target_dir),
        make_env_assignment("BLOCK3_BAYESIAN_REUSE_EXISTING_FIT_RDS", "TRUE"),
        make_env_assignment("BLOCK3_BAYESIAN_STAN_CHAINS", bayesian_stan_chains),
        make_env_assignment("BLOCK3_BAYESIAN_STAN_ITER", bayesian_stan_iter),
        make_env_assignment("BLOCK3_BAYESIAN_STAN_WARMUP", bayesian_stan_warmup),
        make_env_assignment("BLOCK3_BAYESIAN_STAN_THIN", bayesian_stan_thin),
        make_env_assignment("BLOCK3_BAYESIAN_STAN_ADAPT_DELTA", bayesian_stan_adapt_delta),
        make_env_assignment("BLOCK3_BAYESIAN_STAN_MAX_TREEDEPTH", bayesian_stan_max_treedepth),
        make_env_assignment("BLOCK3_BAYESIAN_POSTERIOR_PRED_DRAWS", bayesian_posterior_draws),
        make_env_assignment("BLOCK3_PNU_SITE_LABEL", pnu_site_label),
        make_env_assignment("BLOCK3_SNU_SITE_LABEL", snu_site_label)
      ),
      log_file = bayesian_log_file
    )
  }

  ensure_required_files(frequentist_target_dir, frequentist_required_files, sprintf("%s frequentist cure export", scenario_id))
  ensure_required_files(bayesian_target_dir, bayesian_required_files, sprintf("%s Bayesian cure export", scenario_id))

  list(
    frequentist_dir = frequentist_target_dir,
    bayesian_dir = bayesian_target_dir,
    frequentist_log = frequentist_log_file,
    bayesian_log = bayesian_log_file
  )
}

standardize_frequentist_outputs <- function(dir_path, scenario_id, scenario_label) {
  horizon_raw_tbl <- readr::read_csv(
    file.path(dir_path, "block3_frequentist_cure_horizon_summary.csv"),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    filter(dataset %in% dataset_registry$frequentist_dataset_key)

  horizon_raw_tbl <- ensure_numeric_columns(
    horizon_raw_tbl,
    c(
      "overall_survival_prob_sd",
      "overall_survival_prob_q025",
      "overall_survival_prob_q50",
      "overall_survival_prob_q975",
      "overall_risk_prob_sd",
      "overall_risk_prob_q025",
      "overall_risk_prob_q50",
      "overall_risk_prob_q975",
      "cure_fraction",
      "susceptible_fraction"
    )
  )

  horizon_tbl <- horizon_raw_tbl %>%
    mutate(
      dataset_scope = if_else(dataset == "merged_site_adjusted", "merged", dataset),
      dataset_label = if_else(dataset_scope == "merged", "Merged", dataset_scope),
      scenario_id = scenario_id,
      scenario_label = scenario_label,
      model_group = "Frequentist cure",
      model_name = "frequentist_mixture_cure_lognormal",
      model_label = "Frequentist mixture cure",
      interval_type = "95% CI",
      horizon_year = as.integer(time_horizon_year),
      estimated_survival_probability = as.numeric(overall_survival_prob),
      estimated_survival_uncertainty_sd = as.numeric(overall_survival_prob_sd),
      estimated_survival_lower_95 = as.numeric(overall_survival_prob_q025),
      estimated_survival_median_50 = as.numeric(overall_survival_prob_q50),
      estimated_survival_upper_95 = as.numeric(overall_survival_prob_q975),
      estimated_risk_probability = as.numeric(overall_risk_prob),
      estimated_risk_uncertainty_sd = as.numeric(overall_risk_prob_sd),
      estimated_risk_lower_95 = as.numeric(overall_risk_prob_q025),
      estimated_risk_median_50 = as.numeric(overall_risk_prob_q50),
      estimated_risk_upper_95 = as.numeric(overall_risk_prob_q975),
      n_subjects = NA_integer_
    ) %>%
    select(
      scenario_id,
      scenario_label,
      dataset_scope,
      dataset_label,
      model_group,
      model_name,
      model_label,
      interval_type,
      horizon_year,
      estimated_survival_probability,
      estimated_survival_uncertainty_sd,
      estimated_survival_lower_95,
      estimated_survival_median_50,
      estimated_survival_upper_95,
      estimated_risk_probability,
      estimated_risk_uncertainty_sd,
      estimated_risk_lower_95,
      estimated_risk_median_50,
      estimated_risk_upper_95,
      cure_fraction,
      susceptible_fraction,
      n_subjects
    )

  curve_raw_tbl <- readr::read_csv(
    file.path(dir_path, "block3_frequentist_cure_plot_source.csv"),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    filter(dataset %in% dataset_registry$frequentist_dataset_key)

  curve_raw_tbl <- ensure_numeric_columns(
    curve_raw_tbl,
    c(
      "overall_survival_prob_sd",
      "overall_survival_prob_q025",
      "overall_survival_prob_q50",
      "overall_survival_prob_q975",
      "overall_risk_prob_sd",
      "overall_risk_prob_q025",
      "overall_risk_prob_q50",
      "overall_risk_prob_q975"
    )
  )

  curve_tbl <- curve_raw_tbl %>%
    mutate(
      dataset_scope = if_else(dataset == "merged_site_adjusted", "merged", dataset),
      dataset_label = if_else(dataset_scope == "merged", "Merged", dataset_scope),
      scenario_id = scenario_id,
      scenario_label = scenario_label,
      model_group = "Frequentist cure",
      model_name = "frequentist_mixture_cure_lognormal",
      model_label = "Frequentist mixture cure",
      interval_type = "95% CI",
      time_year = as.numeric(time_horizon_year),
      estimated_survival_probability = as.numeric(overall_survival_prob),
      estimated_survival_uncertainty_sd = as.numeric(overall_survival_prob_sd),
      estimated_survival_lower_95 = as.numeric(overall_survival_prob_q025),
      estimated_survival_median_50 = as.numeric(overall_survival_prob_q50),
      estimated_survival_upper_95 = as.numeric(overall_survival_prob_q975),
      estimated_risk_probability = as.numeric(overall_risk_prob),
      estimated_risk_uncertainty_sd = as.numeric(overall_risk_prob_sd),
      estimated_risk_lower_95 = as.numeric(overall_risk_prob_q025),
      estimated_risk_median_50 = as.numeric(overall_risk_prob_q50),
      estimated_risk_upper_95 = as.numeric(overall_risk_prob_q975)
    ) %>%
    select(
      scenario_id,
      scenario_label,
      dataset_scope,
      dataset_label,
      model_group,
      model_name,
      model_label,
      interval_type,
      time_year,
      estimated_survival_probability,
      estimated_survival_uncertainty_sd,
      estimated_survival_lower_95,
      estimated_survival_median_50,
      estimated_survival_upper_95,
      estimated_risk_probability,
      estimated_risk_uncertainty_sd,
      estimated_risk_lower_95,
      estimated_risk_median_50,
      estimated_risk_upper_95
    )

  list(horizon = horizon_tbl, curve = curve_tbl)
}

standardize_bayesian_outputs <- function(dir_path, scenario_id, scenario_label) {
  horizon_raw_tbl <- readr::read_csv(
    file.path(dir_path, "block3_bayesian_cure_horizon_summary.csv"),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    filter(
      dataset %in% dataset_registry$bayesian_dataset_key,
      is.na(risk_scale) | risk_scale == "transition_only_main"
    )

  horizon_raw_tbl <- ensure_numeric_columns(
    horizon_raw_tbl,
    c(
      "overall_survival_prob_sd",
      "overall_survival_prob_q025",
      "overall_survival_prob_q50",
      "overall_survival_prob_q975",
      "overall_risk_prob_sd",
      "overall_risk_prob_q025",
      "overall_risk_prob_q50",
      "overall_risk_prob_q975",
      "cure_fraction",
      "susceptible_fraction",
      "n_subjects"
    )
  )

  horizon_tbl <- horizon_raw_tbl %>%
    mutate(
      dataset_scope = if_else(dataset == "merged_site_adjusted", "merged", dataset),
      dataset_label = if_else(dataset_scope == "merged", "Merged", dataset_scope),
      scenario_id = scenario_id,
      scenario_label = scenario_label,
      model_group = "Bayesian cure",
      model_name = "bayesian_mixture_cure",
      model_label = "Bayesian cure",
      interval_type = "95% CrI",
      horizon_year = as.integer(dplyr::coalesce(horizon_year, time_horizon_year)),
      estimated_survival_probability = as.numeric(overall_survival_prob),
      estimated_survival_uncertainty_sd = as.numeric(overall_survival_prob_sd),
      estimated_survival_lower_95 = as.numeric(overall_survival_prob_q025),
      estimated_survival_median_50 = as.numeric(overall_survival_prob_q50),
      estimated_survival_upper_95 = as.numeric(overall_survival_prob_q975),
      estimated_risk_probability = as.numeric(overall_risk_prob),
      estimated_risk_uncertainty_sd = as.numeric(overall_risk_prob_sd),
      estimated_risk_lower_95 = as.numeric(overall_risk_prob_q025),
      estimated_risk_median_50 = as.numeric(overall_risk_prob_q50),
      estimated_risk_upper_95 = as.numeric(overall_risk_prob_q975),
      n_subjects = as.integer(n_subjects)
    ) %>%
    select(
      scenario_id,
      scenario_label,
      dataset_scope,
      dataset_label,
      model_group,
      model_name,
      model_label,
      interval_type,
      horizon_year,
      estimated_survival_probability,
      estimated_survival_uncertainty_sd,
      estimated_survival_lower_95,
      estimated_survival_median_50,
      estimated_survival_upper_95,
      estimated_risk_probability,
      estimated_risk_uncertainty_sd,
      estimated_risk_lower_95,
      estimated_risk_median_50,
      estimated_risk_upper_95,
      cure_fraction,
      susceptible_fraction,
      n_subjects
    )

  curve_raw_tbl <- readr::read_csv(
    file.path(dir_path, "block3_bayesian_cure_plot_source.csv"),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    filter(
      dataset %in% dataset_registry$bayesian_dataset_key,
      is.na(risk_scale) | risk_scale == "transition_only_main"
    )

  curve_raw_tbl <- ensure_numeric_columns(
    curve_raw_tbl,
    c(
      "overall_survival_prob_sd",
      "overall_survival_prob_q025",
      "overall_survival_prob_q50",
      "overall_survival_prob_q975",
      "overall_risk_prob_sd",
      "overall_risk_prob_q025",
      "overall_risk_prob_q50",
      "overall_risk_prob_q975"
    )
  )

  curve_tbl <- curve_raw_tbl %>%
    mutate(
      dataset_scope = if_else(dataset == "merged_site_adjusted", "merged", dataset),
      dataset_label = if_else(dataset_scope == "merged", "Merged", dataset_scope),
      scenario_id = scenario_id,
      scenario_label = scenario_label,
      model_group = "Bayesian cure",
      model_name = "bayesian_mixture_cure",
      model_label = "Bayesian cure",
      interval_type = "95% CrI",
      time_year = as.numeric(time_horizon_year),
      estimated_survival_probability = as.numeric(overall_survival_prob),
      estimated_survival_uncertainty_sd = as.numeric(overall_survival_prob_sd),
      estimated_survival_lower_95 = as.numeric(overall_survival_prob_q025),
      estimated_survival_median_50 = as.numeric(overall_survival_prob_q50),
      estimated_survival_upper_95 = as.numeric(overall_survival_prob_q975),
      estimated_risk_probability = as.numeric(overall_risk_prob),
      estimated_risk_uncertainty_sd = as.numeric(overall_risk_prob_sd),
      estimated_risk_lower_95 = as.numeric(overall_risk_prob_q025),
      estimated_risk_median_50 = as.numeric(overall_risk_prob_q50),
      estimated_risk_upper_95 = as.numeric(overall_risk_prob_q975)
    ) %>%
    select(
      scenario_id,
      scenario_label,
      dataset_scope,
      dataset_label,
      model_group,
      model_name,
      model_label,
      interval_type,
      time_year,
      estimated_survival_probability,
      estimated_survival_uncertainty_sd,
      estimated_survival_lower_95,
      estimated_survival_median_50,
      estimated_survival_upper_95,
      estimated_risk_probability,
      estimated_risk_uncertainty_sd,
      estimated_risk_lower_95,
      estimated_risk_median_50,
      estimated_risk_upper_95
    )

  list(horizon = horizon_tbl, curve = curve_tbl)
}

# Define: comparison and reporting helpers ===============================
extract_full_vs_truncated_differences <- function(probability_tbl) {
  full_tbl <- probability_tbl %>%
    filter(scenario_id == "full") %>%
    select(
      dataset_scope,
      model_name,
      horizon_year,
      full_survival_probability = estimated_survival_probability,
      full_risk_probability = estimated_risk_probability,
      full_cure_fraction = cure_fraction,
      full_susceptible_fraction = susceptible_fraction
    )

  probability_tbl %>%
    filter(scenario_id != "full") %>%
    left_join(full_tbl, by = c("dataset_scope", "model_name", "horizon_year")) %>%
    mutate(
      survival_diff_vs_full = estimated_survival_probability - full_survival_probability,
      risk_diff_vs_full = estimated_risk_probability - full_risk_probability,
      abs_risk_diff_vs_full = abs(risk_diff_vs_full),
      cure_fraction_diff_vs_full = cure_fraction - full_cure_fraction,
      susceptible_fraction_diff_vs_full = susceptible_fraction - full_susceptible_fraction
    )
}

summarize_divergence <- function(divergence_long_tbl) {
  extract_point <- function(df, horizon_value, column_name) {
    out <- df %>%
      filter(horizon_year == horizon_value) %>%
      pull({{ column_name }})
    if (length(out) == 0L) NA_real_ else as.numeric(out[[1L]])
  }

  divergence_long_tbl %>%
    group_by(dataset_scope, scenario_id, scenario_label, model_group, model_name, model_label) %>%
    group_modify(
      ~ tibble(
        max_abs_risk_diff_vs_full = max(.x$abs_risk_diff_vs_full, na.rm = TRUE),
        risk_diff_vs_full_at_2y = extract_point(.x, 2L, risk_diff_vs_full),
        risk_diff_vs_full_at_3y = extract_point(.x, 3L, risk_diff_vs_full),
        risk_diff_vs_full_at_5y = extract_point(.x, 5L, risk_diff_vs_full),
        risk_diff_vs_full_at_10y = extract_point(.x, 10L, risk_diff_vs_full),
        cure_fraction_diff_vs_full = {
          tmp <- unique(.x$cure_fraction_diff_vs_full[!is.na(.x$cure_fraction_diff_vs_full)])
          if (length(tmp) == 0L) NA_real_ else as.numeric(tmp[[1L]])
        }
      )
    ) %>%
    ungroup()
}

classify_cure_signal <- function(freq_delta, bayes_delta, freq_trunc, bayes_trunc) {
  deltas <- c(freq_delta, bayes_delta)
  trunc_values <- c(freq_trunc, bayes_trunc)
  valid_delta <- deltas[!is.na(deltas)]
  valid_trunc <- trunc_values[!is.na(trunc_values)]

  if (length(valid_trunc) > 0L && all(valid_trunc <= 0.05)) {
    return("minimal_under_truncation")
  }
  if (length(valid_delta) > 0L && all(valid_delta <= -0.05)) {
    return("weakened_under_truncation")
  }
  if (length(valid_delta) > 0L && all(abs(valid_delta) < 0.05)) {
    return("largely_preserved")
  }
  "mixed_change"
}

build_cure_signal_summary <- function(probability_tbl) {
  cure_fraction_tbl <- probability_tbl %>%
    filter(model_name %in% c("frequentist_mixture_cure_lognormal", "bayesian_mixture_cure")) %>%
    select(dataset_scope, scenario_id, scenario_label, model_name, cure_fraction) %>%
    distinct()

  risk_10_tbl <- probability_tbl %>%
    filter(horizon_year == 10L) %>%
    select(dataset_scope, scenario_id, model_name, estimated_risk_probability)

  risk_10_wide <- risk_10_tbl %>%
    filter(model_name %in% c("no_cure_lognormal", "frequentist_mixture_cure_lognormal", "bayesian_mixture_cure")) %>%
    select(dataset_scope, scenario_id, model_name, estimated_risk_probability) %>%
    tidyr::pivot_wider(
      names_from = model_name,
      values_from = estimated_risk_probability,
      names_prefix = "risk10_"
    )

  full_cure_tbl <- cure_fraction_tbl %>%
    filter(scenario_id == "full") %>%
    select(dataset_scope, model_name, full_cure_fraction = cure_fraction)

  cure_fraction_tbl %>%
    filter(scenario_id != "full") %>%
    left_join(full_cure_tbl, by = c("dataset_scope", "model_name")) %>%
    mutate(cure_fraction_diff_vs_full = cure_fraction - full_cure_fraction) %>%
    select(-full_cure_fraction) %>%
    tidyr::pivot_wider(
      names_from = model_name,
      values_from = c(cure_fraction, cure_fraction_diff_vs_full),
      names_sep = "__"
    ) %>%
    left_join(risk_10_wide, by = c("dataset_scope", "scenario_id")) %>%
    mutate(
      cure_signal_interpretation = purrr::pmap_chr(
        list(
          cure_fraction_diff_vs_full__frequentist_mixture_cure_lognormal,
          cure_fraction_diff_vs_full__bayesian_mixture_cure,
          cure_fraction__frequentist_mixture_cure_lognormal,
          cure_fraction__bayesian_mixture_cure
        ),
        classify_cure_signal
      )
    )
}

make_km_plot <- function(km_curve_tbl) {
  scenario_levels <- scenario_registry$scenario_label

  plot_df <- km_curve_tbl %>%
    mutate(
      scenario_label = factor(scenario_label, levels = scenario_levels),
      dataset_label = factor(dataset_label, levels = dataset_registry$dataset_label)
    )

  plot_object <- ggplot(
    plot_df,
    aes(x = time_year, y = estimated_risk_probability, color = scenario_label)
  )

  if (all(c("estimated_risk_lower_95", "estimated_risk_upper_95") %in% names(plot_df))) {
    plot_object <- plot_object +
      geom_ribbon(
        aes(
          ymin = estimated_risk_lower_95,
          ymax = estimated_risk_upper_95,
          fill = scenario_label,
          group = scenario_label
        ),
        inherit.aes = TRUE,
        alpha = 0.12,
        color = NA
      )
  }

  plot_object +
    geom_step(linewidth = 0.9, alpha = 0.95) +
    facet_wrap(~ dataset_label, ncol = 1) +
    scale_fill_manual(
      values = c(
        "Full follow-up" = "#1B4332",
        "Truncated at 2 years" = "#D97706",
        "Truncated at 3 years" = "#2A6F97"
      ),
      guide = "none"
    ) +
    scale_color_manual(
      values = c(
        "Full follow-up" = "#1B4332",
        "Truncated at 2 years" = "#D97706",
        "Truncated at 3 years" = "#2A6F97"
      )
    ) +
    scale_x_continuous(breaks = 0:10, limits = c(0, 10)) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = "Block 3 KM risk curves: full vs truncated follow-up",
      x = "Years since cohort entry",
      y = "Transition risk",
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

make_model_plot <- function(model_curve_tbl) {
  scenario_levels <- scenario_registry$scenario_label

  model_order <- c(
    "No-cure exponential",
    "No-cure weibull",
    "No-cure lognormal",
    "No-cure loglogistic",
    "Frequentist mixture cure",
    "Bayesian cure"
  )

  plot_df <- model_curve_tbl %>%
    filter(model_group != "KM") %>%
    mutate(
      scenario_label = factor(scenario_label, levels = scenario_levels),
      dataset_label = factor(dataset_label, levels = dataset_registry$dataset_label),
      model_label = factor(model_label, levels = model_order)
    )

  plot_object <- ggplot(
    plot_df,
    aes(x = time_year, y = estimated_risk_probability, color = scenario_label)
  )

  if (all(c("estimated_risk_lower_95", "estimated_risk_upper_95") %in% names(plot_df))) {
    plot_object <- plot_object +
      geom_ribbon(
        aes(
          ymin = estimated_risk_lower_95,
          ymax = estimated_risk_upper_95,
          fill = scenario_label,
          group = scenario_label
        ),
        inherit.aes = TRUE,
        alpha = 0.12,
        color = NA
      )
  }

  plot_object +
    geom_line(linewidth = 0.85, alpha = 0.95) +
    facet_grid(model_label ~ dataset_label) +
    scale_fill_manual(
      values = c(
        "Full follow-up" = "#1B4332",
        "Truncated at 2 years" = "#D97706",
        "Truncated at 3 years" = "#2A6F97"
      ),
      guide = "none"
    ) +
    scale_color_manual(
      values = c(
        "Full follow-up" = "#1B4332",
        "Truncated at 2 years" = "#D97706",
        "Truncated at 3 years" = "#2A6F97"
      )
    ) +
    scale_x_continuous(breaks = 0:10, limits = c(0, 10)) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = "Block 3 model-based risk curves: full vs truncated follow-up",
      x = "Years since cohort entry",
      y = "Transition risk",
      color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

write_summary_markdown <- function(path, divergence_summary_tbl, cure_signal_tbl) {
  top_divergence_tbl <- divergence_summary_tbl %>%
    arrange(desc(max_abs_risk_diff_vs_full)) %>%
    slice_head(n = 8)

  lines <- c(
    "# Block 3 summary",
    "",
    "## Largest divergence from full follow-up",
    ""
  )

  if (nrow(top_divergence_tbl) == 0L) {
    lines <- c(lines, "- No divergence rows were generated.", "")
  } else {
    divergence_lines <- sprintf(
      "- %s / %s / %s: max |risk diff| = %.3f; 5y diff = %.3f; 10y diff = %.3f",
      top_divergence_tbl$dataset_scope,
      top_divergence_tbl$scenario_label,
      top_divergence_tbl$model_label,
      top_divergence_tbl$max_abs_risk_diff_vs_full,
      top_divergence_tbl$risk_diff_vs_full_at_5y,
      top_divergence_tbl$risk_diff_vs_full_at_10y
    )
    lines <- c(lines, divergence_lines, "")
  }

  lines <- c(lines, "## Cure-like signal summary", "")

  if (nrow(cure_signal_tbl) == 0L) {
    lines <- c(lines, "- No cure-like signal summary rows were generated.", "")
  } else {
    cure_lines <- sprintf(
      "- %s / %s: frequentist cure fraction = %.3f, Bayesian cure fraction = %.3f, interpretation = %s",
      cure_signal_tbl$dataset_scope,
      cure_signal_tbl$scenario_label,
      cure_signal_tbl$cure_fraction__frequentist_mixture_cure_lognormal,
      cure_signal_tbl$cure_fraction__bayesian_mixture_cure,
      cure_signal_tbl$cure_signal_interpretation
    )
    lines <- c(lines, cure_lines, "")
  }

  writeLines(lines, con = path)
}

# Build: source scenarios ===============================
run_log <- new_run_log()

source_df <- read_source_data(source_data_file, pnu_site_label, snu_site_label)
run_log <- append_run_log(
  run_log,
  component = "source_data",
  status = "loaded",
  detail = sprintf("Loaded source rows: %s", scales::comma(nrow(source_df))),
  file_path = source_data_file
)

scenario_source_registry <- purrr::pmap_dfr(
  scenario_registry,
  function(scenario_id, scenario_label, truncation_year) {
    scenario_df <- apply_followup_truncation(source_df, truncation_year)
    scenario_source_file <- make_scenario_source_file(scenario_id)
    readr::write_csv(scenario_df, scenario_source_file)

    tibble(
      scenario_id = scenario_id,
      scenario_label = scenario_label,
      truncation_year = truncation_year,
      source_file = scenario_source_file,
      n_rows = nrow(scenario_df),
      n_transition = sum(as.integer(safe_numeric(scenario_df$status_num)) == 1L, na.rm = TRUE)
    )
  }
)

run_log <- append_run_log(
  run_log,
  component = "scenario_sources",
  status = "written",
  detail = sprintf("Scenario source files created: %s", nrow(scenario_source_registry)),
  file_path = scenario_source_dir
)

# Build: KM and no-cure comparisons ===============================
analysis_results <- purrr::pmap(
  scenario_registry,
  function(scenario_id, scenario_label, truncation_year) {
    scenario_source_file <- scenario_source_registry %>%
      filter(scenario_id == !!scenario_id) %>%
      pull(source_file)

    scenario_df <- read_csv_text(scenario_source_file, paste("Scenario source", scenario_id))
    scenario_df <- standardize_known_site_labels(scenario_df, pnu_site_label, snu_site_label)
    analysis_datasets <- build_analysis_datasets_from_source_df(scenario_df, pnu_site_label, snu_site_label)

    dataset_summaries <- bind_rows(
      summarize_dataset_counts(analysis_datasets$SNU, "SNU", scenario_id),
      summarize_dataset_counts(analysis_datasets$merged, "merged", scenario_id)
    ) %>%
      left_join(
        scenario_registry %>% select(scenario_id, scenario_label),
        by = "scenario_id"
      ) %>%
      relocate(scenario_label, .after = scenario_id)

    km_exports <- purrr::pmap(
      dataset_registry,
      function(dataset_scope, dataset_label, site_mode, formula_rhs, frequentist_dataset_key, bayesian_dataset_key) {
        df <- analysis_datasets[[dataset_scope]]
        build_km_exports(
          df = df,
          scenario_id = scenario_id,
          scenario_label = scenario_label,
          dataset_scope = dataset_scope,
          dataset_label = dataset_label,
          time_grid = curve_time_grid_year,
          horizon_years = common_horizons_year
        )
      }
    )

    nocure_exports <- purrr::pmap(
      dataset_registry,
      function(dataset_scope, dataset_label, site_mode, formula_rhs, frequentist_dataset_key, bayesian_dataset_key) {
        df <- analysis_datasets[[dataset_scope]]
        purrr::map(
          seq_len(nrow(no_cure_family_registry)),
          function(ii) {
            build_nocure_exports(
              df = df,
              scenario_id = scenario_id,
              scenario_label = scenario_label,
              dataset_scope = dataset_scope,
              dataset_label = dataset_label,
              formula_rhs = formula_rhs,
              family_row = no_cure_family_registry[ii, , drop = FALSE],
              time_grid = curve_time_grid_year,
              horizon_years = common_horizons_year
            )
          }
        )
      }
    )

    list(
      dataset_summary = dataset_summaries,
      km_horizon = bind_rows(purrr::map(km_exports, "horizon")),
      km_curve = bind_rows(purrr::map(km_exports, "curve")),
      nocure_registry = bind_rows(unlist(purrr::map(nocure_exports, ~ purrr::map(.x, "registry")), recursive = FALSE)),
      nocure_horizon = bind_rows(unlist(purrr::map(nocure_exports, ~ purrr::map(.x, "horizon")), recursive = FALSE)),
      nocure_curve = bind_rows(unlist(purrr::map(nocure_exports, ~ purrr::map(.x, "curve")), recursive = FALSE))
    )
  }
)

dataset_summary_tbl <- bind_rows(purrr::map(analysis_results, "dataset_summary"))
km_horizon_tbl <- bind_rows(purrr::map(analysis_results, "km_horizon"))
km_curve_tbl <- bind_rows(purrr::map(analysis_results, "km_curve"))
nocure_registry_tbl <- bind_rows(purrr::map(analysis_results, "nocure_registry"))
nocure_horizon_tbl <- bind_rows(purrr::map(analysis_results, "nocure_horizon"))
nocure_curve_tbl <- bind_rows(purrr::map(analysis_results, "nocure_curve"))

run_log <- append_run_log(
  run_log,
  component = "km_and_nocure",
  status = "assembled",
  detail = sprintf(
    "KM rows = %s; no-cure horizon rows = %s",
    scales::comma(nrow(km_horizon_tbl)),
    scales::comma(nrow(nocure_horizon_tbl))
  ),
  file_path = export_path
)

# Build: frequentist and Bayesian cure exports ===============================
remove_legacy_block3_paths()

model_export_registry <- purrr::pmap_dfr(
  scenario_registry,
  function(scenario_id, scenario_label, truncation_year) {
    scenario_source_file <- scenario_source_registry %>%
      filter(scenario_id == !!scenario_id) %>%
      pull(source_file)

    model_dirs <- ensure_scenario_model_exports(scenario_id, scenario_source_file)

    tibble(
      scenario_id = scenario_id,
      scenario_label = scenario_label,
      frequentist_dir = model_dirs$frequentist_dir,
      bayesian_dir = model_dirs$bayesian_dir,
      frequentist_log = model_dirs$frequentist_log,
      bayesian_log = model_dirs$bayesian_log
    )
  }
)

if (nrow(model_export_registry) > 0L) {
  run_log <- append_run_log(
    run_log,
    component = "model_exports",
    status = "ready",
    detail = sprintf("Prepared Block 3 frequentist and Bayesian cure exports for %s scenarios.", nrow(model_export_registry)),
    file_path = supporting_model_dir
  )
}

full_model_dirs <- model_export_registry %>%
  filter(scenario_id == "full")

if (nrow(full_model_dirs) != 1L) {
  stop("Expected exactly one full-follow-up model export row.", call. = FALSE)
}

interpretation_bundle_dirs <- ensure_interpretation_bundle(
  full_frequentist_dir = full_model_dirs$frequentist_dir[[1L]],
  full_bayesian_dir = full_model_dirs$bayesian_dir[[1L]]
)
run_log <- append_run_log(
  run_log,
  component = "interpretation_diagnostics",
  status = "ready",
  detail = "Copied full-follow-up Block 3 cure diagnostics into the interpretation bundle.",
  file_path = interpretation_dir
)

frequentist_model_tbl <- purrr::pmap(
  model_export_registry,
  function(scenario_id, scenario_label, frequentist_dir, bayesian_dir, frequentist_log, bayesian_log) {
    standardize_frequentist_outputs(frequentist_dir, scenario_id, scenario_label)
  }
)
bayesian_model_tbl <- purrr::pmap(
  model_export_registry,
  function(scenario_id, scenario_label, frequentist_dir, bayesian_dir, frequentist_log, bayesian_log) {
    standardize_bayesian_outputs(bayesian_dir, scenario_id, scenario_label)
  }
)

freq_cure_horizon_tbl <- bind_rows(
  purrr::map(frequentist_model_tbl, "horizon")
)
freq_cure_curve_tbl <- bind_rows(
  purrr::map(frequentist_model_tbl, "curve")
)
bayes_horizon_tbl <- bind_rows(
  purrr::map(bayesian_model_tbl, "horizon")
)
bayes_curve_tbl <- bind_rows(
  purrr::map(bayesian_model_tbl, "curve")
)

# Assemble: combined tables, plots, and summaries ===============================
probability_tbl <- bind_rows(
  km_horizon_tbl %>%
    select(
      scenario_id,
      scenario_label,
      dataset_scope,
      dataset_label,
      model_group,
      model_name,
      model_label,
      interval_type,
      horizon_year,
      estimated_survival_probability,
      estimated_survival_uncertainty_sd,
      estimated_survival_lower_95,
      estimated_survival_median_50,
      estimated_survival_upper_95,
      estimated_risk_probability,
      estimated_risk_uncertainty_sd,
      estimated_risk_lower_95,
      estimated_risk_median_50,
      estimated_risk_upper_95,
      cure_fraction,
      susceptible_fraction,
      n_subjects
    ),
  nocure_horizon_tbl %>%
    select(
      scenario_id,
      scenario_label,
      dataset_scope,
      dataset_label,
      model_group,
      model_name,
      model_label,
      interval_type,
      horizon_year,
      estimated_survival_probability,
      estimated_survival_uncertainty_sd,
      estimated_survival_lower_95,
      estimated_survival_median_50,
      estimated_survival_upper_95,
      estimated_risk_probability,
      estimated_risk_uncertainty_sd,
      estimated_risk_lower_95,
      estimated_risk_median_50,
      estimated_risk_upper_95,
      cure_fraction,
      susceptible_fraction,
      n_subjects
    ),
  freq_cure_horizon_tbl,
  bayes_horizon_tbl
) %>%
  arrange(
    match(dataset_scope, dataset_registry$dataset_scope),
    match(scenario_id, scenario_registry$scenario_id),
    match(
      model_name,
      c(
        "km",
        "no_cure_exponential",
        "no_cure_weibull",
        "no_cure_lognormal",
        "no_cure_loglogistic",
        "frequentist_mixture_cure_lognormal",
        "bayesian_mixture_cure"
      )
    ),
    horizon_year
  )

curve_comparison_tbl <- bind_rows(
  km_curve_tbl,
  nocure_curve_tbl %>%
    select(
      scenario_id,
      scenario_label,
      dataset_scope,
      dataset_label,
      model_group,
      model_name,
      model_label,
      interval_type,
      time_year,
      estimated_survival_probability,
      estimated_survival_uncertainty_sd,
      estimated_survival_lower_95,
      estimated_survival_median_50,
      estimated_survival_upper_95,
      estimated_risk_probability,
      estimated_risk_uncertainty_sd,
      estimated_risk_lower_95,
      estimated_risk_median_50,
      estimated_risk_upper_95
    ),
  freq_cure_curve_tbl,
  bayes_curve_tbl
) %>%
  arrange(
    match(dataset_scope, dataset_registry$dataset_scope),
    match(scenario_id, scenario_registry$scenario_id),
    model_name,
    time_year
  )

divergence_long_tbl <- extract_full_vs_truncated_differences(probability_tbl)
divergence_summary_tbl <- summarize_divergence(divergence_long_tbl)
cure_signal_tbl <- build_cure_signal_summary(probability_tbl)

km_plot <- make_km_plot(km_curve_tbl)
model_plot <- make_model_plot(curve_comparison_tbl)

ggsave(km_plot_file, plot = km_plot, width = plot_width_in, height = 7, dpi = plot_dpi)
ggsave(model_plot_file, plot = model_plot, width = plot_width_in, height = plot_height_in, dpi = plot_dpi)

metadata_tbl <- tibble::tribble(
  ~metadata_group, ~metadata_name, ~metadata_value,
  "block", "block_name", "Block 3 short vs long follow-up direct comparison",
  "block", "spec_reference", "1.Model specifciation/spec.md",
  "inputs", "source_data_file", normalize_existing_path(source_data_file),
  "inputs", "export_path", normalize_existing_path(export_path),
  "outputs", "interpretation_diagnostics_path", normalize_existing_path(interpretation_dir),
  "outputs", "supporting_model_exports_path", normalize_existing_path(supporting_model_dir),
  "inputs", "system_rscript", normalize_existing_path(system_rscript),
  "scenario", "scenario_ids", paste(scenario_registry$scenario_id, collapse = ","),
  "scenario", "scenario_labels", paste(scenario_registry$scenario_label, collapse = " | "),
  "scenario", "truncation_years", paste(ifelse(is.na(scenario_registry$truncation_year), "full", scenario_registry$truncation_year), collapse = ","),
  "modeling", "no_cure_families", paste(no_cure_family_registry$model_family, collapse = ","),
  "modeling", "frequentist_cure_engine_script", normalize_existing_path(frequentist_cure_script),
  "modeling", "bayesian_cure_engine_script", normalize_existing_path(bayesian_cure_script),
  "modeling", "bayesian_stan_chains", as.character(bayesian_stan_chains),
  "modeling", "bayesian_stan_iter", as.character(bayesian_stan_iter),
  "modeling", "bayesian_stan_warmup", as.character(bayesian_stan_warmup),
  "modeling", "no_cure_ci_draws", as.character(no_cure_ci_draws),
  "time", "common_horizons_year", paste(common_horizons_year, collapse = ","),
  "time", "curve_time_step_year", format(diff(curve_time_grid_year)[[1L]], trim = TRUE),
  "event", "event_definition", "status_num == 1",
  "event", "main_censoring_definition", "status_num %in% c(0, 2)",
  "rule", "truncation_rule", "If follow-up exceeds the scenario cutoff, follow-up time is truncated to the cutoff and status_num is recoded to 0.",
  "rule", "merged_model_formula", "age_s + sex_num + site",
  "rule", "snu_model_formula", "age_s + sex_num"
)

block3_standard_error_bundle <- build_standard_error_export_bundle(list(
  dataset_summary_tbl = dataset_summary_tbl,
  km_horizon_tbl = km_horizon_tbl,
  km_curve_tbl = km_curve_tbl,
  nocure_registry_tbl = nocure_registry_tbl,
  nocure_horizon_tbl = nocure_horizon_tbl,
  nocure_curve_tbl = nocure_curve_tbl,
  freq_cure_horizon_tbl = freq_cure_horizon_tbl,
  freq_cure_curve_tbl = freq_cure_curve_tbl,
  bayes_horizon_tbl = bayes_horizon_tbl,
  bayes_curve_tbl = bayes_curve_tbl,
  probability_tbl = probability_tbl,
  curve_comparison_tbl = curve_comparison_tbl,
  divergence_long_tbl = divergence_long_tbl,
  divergence_summary_tbl = divergence_summary_tbl,
  cure_signal_tbl = cure_signal_tbl
))

write_summary_markdown(summary_markdown_file, divergence_summary_tbl, cure_signal_tbl)

readr::write_csv(scenario_source_registry, scenario_manifest_file)
readr::write_csv(dataset_summary_tbl, dataset_summary_file)
readr::write_csv(km_horizon_tbl, km_horizon_file)
readr::write_csv(km_curve_tbl, km_curve_file)
readr::write_csv(nocure_registry_tbl, nocure_registry_file)
readr::write_csv(nocure_horizon_tbl, nocure_horizon_file)
readr::write_csv(nocure_curve_tbl, nocure_curve_file)
readr::write_csv(freq_cure_horizon_tbl, freq_cure_horizon_file)
readr::write_csv(freq_cure_curve_tbl, freq_cure_curve_file)
readr::write_csv(bayes_horizon_tbl, bayes_horizon_file)
readr::write_csv(bayes_curve_tbl, bayes_curve_file)
readr::write_csv(probability_tbl, probability_table_file)
readr::write_csv(curve_comparison_tbl, curve_table_file)
readr::write_csv(divergence_long_tbl, divergence_long_file)
readr::write_csv(divergence_summary_tbl, divergence_summary_file)
readr::write_csv(cure_signal_tbl, cure_signal_file)
readr::write_csv(block3_standard_error_bundle$registry, standard_error_registry_file)
readr::write_csv(block3_standard_error_bundle$long, standard_error_long_file)
readr::write_csv(metadata_tbl, metadata_registry_file)

run_log <- append_run_log(
  run_log,
  component = "exports",
  status = "written",
  detail = "Block 3 summary tables and plots were written.",
  file_path = export_path
)
readr::write_csv(run_log, run_log_file)

interpretation_file_paths <- c(
  interpretation_bundle_dirs$readme,
  file.path(interpretation_bundle_dirs$interpretation_dir, interpretation_priority_files),
  file.path(interpretation_bundle_dirs$frequentist_dir, frequentist_interpretation_subdir_files),
  file.path(interpretation_bundle_dirs$bayesian_dir, bayesian_interpretation_subdir_files)
)
interpretation_object_names <- c(
  "interpretation_readme",
  "frequentist_interpretation_fit_summary",
  "frequentist_interpretation_detail_annotation",
  "bayesian_interpretation_parameter_summary",
  "bayesian_interpretation_posterior_draw_summary",
  "bayesian_interpretation_model_registry",
  "bayesian_interpretation_anchor_vs_neutral_delta",
  "frequentist_interpretation_horizon_summary",
  "frequentist_interpretation_fitted_objects",
  "frequentist_interpretation_plot_source",
  "frequentist_interpretation_latency_plot_source",
  "bayesian_interpretation_horizon_summary",
  "bayesian_interpretation_prior_branch_horizon_summary",
  "bayesian_interpretation_plot_source",
  "bayesian_interpretation_detail_annotation",
  "bayesian_interpretation_latency_plot_source",
  "bayesian_interpretation_uncured_latency_summary",
  "bayesian_interpretation_latency_aft_effect_summary",
  "bayesian_interpretation_trace_plots"
)
interpretation_descriptions <- c(
  "Guide to the Block 3 full-follow-up cure-model files that support interpretation.",
  "Full-follow-up frequentist mixture-cure fit summary copied to the top-level Block 3 interpretation folder.",
  "Full-follow-up frequentist mixture-cure detail annotation table copied to the top-level Block 3 interpretation folder.",
  "Full-follow-up Bayesian parameter summary copied to the top-level Block 3 interpretation folder.",
  "Full-follow-up Bayesian posterior-draw summary copied to the top-level Block 3 interpretation folder.",
  "Full-follow-up Bayesian model registry copied to the top-level Block 3 interpretation folder.",
  "Full-follow-up Bayesian anchor-vs-neutral delta summary copied to the top-level Block 3 interpretation folder.",
  "Full-follow-up frequentist mixture-cure horizon summary copied to the frequentist supplemental interpretation subfolder.",
  "Full-follow-up frequentist mixture-cure fitted-object RDS copied to the frequentist supplemental interpretation subfolder.",
  "Full-follow-up frequentist mixture-cure plot source copied to the frequentist supplemental interpretation subfolder.",
  "Full-follow-up frequentist mixture-cure latency plot source copied to the frequentist supplemental interpretation subfolder.",
  "Full-follow-up Bayesian mixture-cure horizon summary copied to the Bayesian supplemental interpretation subfolder.",
  "Full-follow-up Bayesian prior-branch horizon summary copied to the Bayesian supplemental interpretation subfolder.",
  "Full-follow-up Bayesian plot source copied to the Bayesian supplemental interpretation subfolder.",
  "Full-follow-up Bayesian detail annotation table copied to the Bayesian supplemental interpretation subfolder.",
  "Full-follow-up Bayesian latency plot source copied to the Bayesian supplemental interpretation subfolder.",
  "Full-follow-up Bayesian uncured-latency summary copied for Block 3 interpretation.",
  "Full-follow-up Bayesian latency AFT effect summary copied for Block 3 interpretation.",
  "Full-follow-up Bayesian trace plots copied for Block 3 interpretation."
)

manifest_tbl <- make_manifest_tbl(
  file_paths = c(
    interpretation_file_paths,
    scenario_manifest_file,
    dataset_summary_file,
    km_horizon_file,
    km_curve_file,
    nocure_registry_file,
    nocure_horizon_file,
    nocure_curve_file,
    freq_cure_horizon_file,
    freq_cure_curve_file,
    bayes_horizon_file,
    bayes_curve_file,
    probability_table_file,
    curve_table_file,
    divergence_long_file,
    divergence_summary_file,
    cure_signal_file,
    standard_error_registry_file,
    standard_error_long_file,
    metadata_registry_file,
    run_log_file,
    summary_markdown_file,
    km_plot_file,
    model_plot_file,
    manifest_file
  ),
  object_names = c(
    interpretation_object_names,
    "scenario_source_registry",
    "dataset_summary_tbl",
    "km_horizon_tbl",
    "km_curve_tbl",
    "nocure_registry_tbl",
    "nocure_horizon_tbl",
    "nocure_curve_tbl",
    "freq_cure_horizon_tbl",
    "freq_cure_curve_tbl",
    "bayes_horizon_tbl",
    "bayes_curve_tbl",
    "probability_tbl",
    "curve_comparison_tbl",
    "divergence_long_tbl",
    "divergence_summary_tbl",
    "cure_signal_tbl",
    "block3_standard_error_registry_tbl",
    "block3_standard_error_long_tbl",
    "metadata_tbl",
    "run_log",
    "summary_markdown",
    "km_plot",
    "model_plot",
    "manifest_tbl"
  ),
  descriptions = c(
    interpretation_descriptions,
    "Scenario source-file registry for full, truncated_2y, and truncated_3y inputs.",
    "Dataset-level counts and follow-up summaries for SNU and merged Block 3 scenarios.",
    "KM annual horizon summary table.",
    "KM curve data on the 0-10 year grid.",
    "No-cure model fit registry with convergence and information criteria.",
    "No-cure annual horizon summary table for four parametric families.",
    "No-cure curve data on the 0-10 year grid.",
    "Frequentist mixture-cure annual horizon summary table.",
    "Frequentist mixture-cure curve data on the 0-10 year grid.",
    "Bayesian cure annual horizon summary table.",
    "Bayesian cure curve data on the 0-10 year grid.",
    "Combined annual probability comparison table across KM, no-cure, frequentist cure, and Bayesian cure.",
    "Combined curve data for model-based Block 3 plots.",
    "Long-form truncated-vs-full divergence table at each common horizon.",
    "Scenario-by-model divergence summary from truncated analyses versus full follow-up.",
    "Descriptive cure-like signal summary across truncated scenarios.",
    "Registry of Block 3 tables that contain standard-error or uncertainty columns.",
    "Long-form Block 3 standard-error table assembled from the exported comparison sources.",
    "Block 3 metadata registry.",
    "Component-level run log for the Block 3 pipeline.",
    "Markdown summary of the largest divergences and cure-like signal changes.",
    "KM risk-curve comparison plot.",
    "Model-based risk-curve comparison plot.",
    "Manifest of Block 3 exported files."
  )
)

readr::write_csv(manifest_tbl, manifest_file)

message("Block 3 short-vs-long follow-up direct comparison export completed.")
