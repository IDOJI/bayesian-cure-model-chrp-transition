# Configure: direct input and output paths ===============================
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
results_root_default <- file.path(repo_root, "3.Results files")
block1_root_default <- "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/1.Block1"
dropbox_export_path_default <- file.path(block1_root_default, "2.MLE No-Cure")

source_data_file <- Sys.getenv(
  "STAGE5_SIMPLE_DATA_FILE",
  unset = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/0.Data/2.Preprocessed data/Preprocessed_Merged_PNUH_SNUH_Data.csv"
)
export_path <- Sys.getenv(
  "STAGE5_SIMPLE_EXPORT_PATH",
  unset = dropbox_export_path_default
)

pnu_site_label <- Sys.getenv("STAGE5_SIMPLE_PNU_SITE_LABEL", unset = "PNU")
snu_site_label <- Sys.getenv("STAGE5_SIMPLE_SNU_SITE_LABEL", unset = "SNU")
output_prefix <- "mle_nocure_lognormal"

horizon_summary_file <- file.path(export_path, paste0(output_prefix, "_horizon_summary.csv"))
plot_source_file <- file.path(export_path, paste0(output_prefix, "_plot_source.csv"))
fitted_models_rds_file <- file.path(export_path, paste0(output_prefix, "_fitted_models.rds"))
survival_plot_rds_file <- file.path(export_path, paste0(output_prefix, "_survival_plot.rds"))
risk_plot_rds_file <- file.path(export_path, paste0(output_prefix, "_risk_plot.rds"))
survival_plot_png_file <- file.path(export_path, paste0(output_prefix, "_survival_plot.png"))
risk_plot_png_file <- file.path(export_path, paste0(output_prefix, "_risk_plot.png"))
plot_bundle_rds_file <- file.path(export_path, paste0(output_prefix, "_plot_bundle.rds"))
detail_plot_bundle_rds_file <- file.path(export_path, paste0(output_prefix, "_detail_plot_bundle.rds"))

required_datasets <- c("PNU", "SNU", "merged")
common_horizons_year <- 1:10
time_origin_epsilon_year <- 1e-08
curve_step_year <- 0.05
curve_horizon_max_year <- 10
risk_table_step_year <- 0.5
detail_plot_width_in <- 14
detail_plot_height_in <- 9
stage5_ci_draws <- as.integer(Sys.getenv("STAGE5_SIMPLE_CI_DRAWS", unset = "4000"))

dataset_version_registry <- tibble::tibble(
  dataset_version_key = c("PNU", "SNU", "merged", "merged_site_adjusted"),
  dataset_version_label = c("PNU", "SNU", "Merged", "Merged (site-adjusted)"),
  source_dataset = c("PNU", "SNU", "merged", "merged"),
  formula_name = c("base", "base", "base", "site_added"),
  expected_site_branch = c("site_free", "site_free", "site_free", "site_adjusted"),
  expected_interaction_branch = rep("no_age_sex_interaction", 4L)
)

dataset_palette <- c(
  "PNU" = "#1B4332",
  "SNU" = "#2A6F97",
  "Merged" = "#C1666B",
  "Merged (site-adjusted)" = "#B8860B"
)

plot_group_registry <- list(
  all = list(dataset_version_keys = c("PNU", "SNU", "merged"), label = "PNU + SNU + Merged"),
  pnu = list(dataset_version_keys = "PNU", label = "PNU"),
  snu = list(dataset_version_keys = "SNU", label = "SNU"),
  merged = list(dataset_version_keys = "merged", label = "Merged")
)

required_packages <- c("readr", "dplyr", "tibble", "ggplot2", "scales", "survival", "flexsurv", "patchwork")
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
  library(survival)
  library(flexsurv)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# Define: helpers ===============================
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

assert_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

read_csv_checked <- function(path, label) {
  assert_file_exists(path, label)
  readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
}

cleanup_existing_stage5_outputs <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    return(invisible(character()))
  }

  sub_dir <- file.path(output_dir, "sub")
  if (dir.exists(sub_dir)) {
    unlink(sub_dir, recursive = TRUE, force = TRUE)
  }

  stale_files <- list.files(
    output_dir,
    pattern = "^(stage5_simple_lognormal_|mle_nocure_lognormal_).*\\.(csv|rds|png)$",
    full.names = TRUE
  )

  if (length(stale_files) == 0L) {
    return(invisible(character()))
  }

  removed_flag <- file.remove(stale_files)
  if (any(!removed_flag)) {
    stop(
      sprintf(
        "Failed to remove stale no-cure outputs: %s",
        paste(basename(stale_files[!removed_flag]), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(stale_files[removed_flag])
}

organize_stage5_outputs <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    return(invisible(character()))
  }

  sub_dir <- file.path(output_dir, "sub")
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)

  keep_files <- c(
    paste0(output_prefix, "_fitted_models.rds"),
    paste0(output_prefix, "_horizon_summary.csv"),
    paste0(output_prefix, "_plot_source.csv"),
    paste0(output_prefix, "_survival_plot.png"),
    paste0(output_prefix, "_survival_plot.rds"),
    paste0(output_prefix, "_risk_plot.png"),
    paste0(output_prefix, "_risk_plot.rds")
  )

  top_level_files <- list.files(
    output_dir,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE
  )
  top_level_files <- top_level_files[file.info(top_level_files)$isdir %in% FALSE]

  move_files <- top_level_files[!(basename(top_level_files) %in% keep_files)]

  if (length(move_files) == 0L) {
    return(invisible(character()))
  }

  moved_flag <- file.rename(move_files, file.path(sub_dir, basename(move_files)))
  if (any(!moved_flag)) {
    stop(
      sprintf(
        "Failed to move no-cure support files into `sub`: %s",
        paste(basename(move_files[!moved_flag]), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(file.path(sub_dir, basename(move_files[moved_flag])))
}

coerce_numeric_text <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

clip_probability_open <- function(x, eps = 1e-12) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

summarize_draws <- function(x) {
  tibble(
    uncertainty_sd = stats::sd(x, na.rm = TRUE),
    interval_lower_95 = unname(stats::quantile(x, probs = 0.025, na.rm = TRUE, names = FALSE)),
    interval_median_50 = unname(stats::quantile(x, probs = 0.5, na.rm = TRUE, names = FALSE)),
    interval_upper_95 = unname(stats::quantile(x, probs = 0.975, na.rm = TRUE, names = FALSE))
  )
}

build_parametric_interval_summaries <- function(fit, df, target_times, n_draws) {
  model_formula <- fit$concat.formula
  if (is.null(model_formula) && !is.null(fit$all.formulae$meanlog)) {
    model_formula <- fit$all.formulae$meanlog
  }
  if (is.null(model_formula)) {
    stop("No-cure fit object is missing the model formula needed for CI calculation.", call. = FALSE)
  }

  design_matrix <- stats::model.matrix(
    stats::delete.response(stats::terms(model_formula)),
    data = df
  )
  parameter_draws <- flexsurv::normboot.flexsurvreg(
    fit,
    B = n_draws,
    newdata = df[1, , drop = FALSE],
    raw = TRUE
  )
  parameter_draws <- as.matrix(parameter_draws)

  if (is.null(colnames(parameter_draws)) || !("sdlog" %in% colnames(parameter_draws))) {
    stop("No-cure parameter draws are missing `sdlog`.", call. = FALSE)
  }

  location_draws <- parameter_draws[, setdiff(colnames(parameter_draws), "sdlog"), drop = FALSE]
  colnames(location_draws) <- ifelse(colnames(location_draws) == "meanlog", "(Intercept)", colnames(location_draws))

  missing_design_cols <- setdiff(colnames(design_matrix), colnames(location_draws))
  if (length(missing_design_cols) > 0L) {
    stop(
      sprintf(
        "No-cure parameter draws are missing design-matrix columns: %s",
        paste(missing_design_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  location_draws <- location_draws[, colnames(design_matrix), drop = FALSE]
  sigma_draws <- pmax(as.numeric(parameter_draws[, "sdlog"]), 1e-8)
  sigma_matrix <- matrix(sigma_draws, nrow = nrow(design_matrix), ncol = nrow(parameter_draws), byrow = TRUE)
  lp_matrix <- design_matrix %*% t(location_draws)

  bind_rows(lapply(target_times, function(tt) {
    time_value <- as.numeric(tt)
    if (time_value <= 0) {
      survival_draws <- rep(1, nrow(parameter_draws))
    } else {
      z_matrix <- (log(pmax(time_value, time_origin_epsilon_year)) - lp_matrix) / sigma_matrix
      survival_draws <- clip_prob(colMeans(1 - stats::pnorm(z_matrix), na.rm = TRUE))
    }
    risk_draws <- clip_prob(1 - survival_draws)

    summarize_draws(survival_draws) %>%
      rename_with(~ paste0("estimated_survival_", .x), everything()) %>%
      bind_cols(
        summarize_draws(risk_draws) %>%
          rename_with(~ paste0("estimated_risk_", .x), everything())
      ) %>%
      mutate(time_year = time_value)
  }))
}

clip_prob <- function(x) {
  pmin(pmax(as.numeric(x), 0), 1)
}

parse_flag <- function(x, default = FALSE) {
  if (length(x) == 0L || is.na(x) || trimws(x) == "") {
    return(default)
  }

  normalized_x <- tolower(trimws(as.character(x)))
  if (normalized_x %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (normalized_x %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }

  default
}

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
      sex_num = as.integer(coerce_numeric_text(sex_num)),
      age_exact_entry = coerce_numeric_text(age_exact_entry),
      days_followup = coerce_numeric_text(days_followup),
      status_num = as.integer(coerce_numeric_text(status_num))
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

  if (any(out$site == "", na.rm = TRUE)) {
    stop(sprintf("[%s] Blank `site` values detected.", dataset_name), call. = FALSE)
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

  out <- out %>%
    mutate(
      unique_person_id = paste(site, id, sep = "_"),
      time_year = days_followup / 365.25,
      event_main = as.integer(status_num == 1L),
      right_censor_flag = as.integer(status_num == 0L),
      remission_flag = as.integer(status_num == 2L),
      censor_main = as.integer(status_num %in% c(0L, 2L)),
      sex_label = factor(if_else(sex_num == 0L, "Female", "Male"), levels = c("Female", "Male")),
      age_s = (age_exact_entry - age_mean) / (2 * age_sd)
    )

  if (nrow(out) != dplyr::n_distinct(out$unique_person_id)) {
    stop(sprintf("[%s] `site + id` is not unique within dataset.", dataset_name), call. = FALSE)
  }

  out
}

validate_analysis_datasets <- function(analysis_datasets) {
  if (!is.list(analysis_datasets) || !all(required_datasets %in% names(analysis_datasets))) {
    stop("Analysis dataset bundle must contain `PNU`, `SNU`, and `merged`.", call. = FALSE)
  }
  invisible(TRUE)
}

build_analysis_datasets_from_source <- function(source_file, pnu_label, snu_label) {
  source_df <- read_csv_checked(source_file, "Merged preprocessed data CSV")
  source_df <- standardize_known_site_labels(source_df, pnu_label, snu_label)
  validate_source_data(source_df, pnu_label, snu_label)

  analysis_datasets <- list(
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

  validate_analysis_datasets(analysis_datasets)
  analysis_datasets
}

normalize_dataset_label <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    toupper(x) == "PNU" ~ "PNU",
    toupper(x) == "SNU" ~ "SNU",
    tolower(x) == "merged" ~ "merged",
    TRUE ~ x
  )
}

build_formula_registry <- function() {
  bind_rows(
    tibble::tibble(
      dataset = c("PNU", "PNU", "SNU", "SNU"),
      formula_name = c("base", "interaction", "base", "interaction"),
      formula_label = c("Base", "Interaction", "Base", "Interaction"),
      formula_rhs = c(
        "age_s + sex_num",
        "age_s + sex_num + age_s:sex_num",
        "age_s + sex_num",
        "age_s + sex_num + age_s:sex_num"
      )
    ),
    tibble::tibble(
      dataset = c("merged", "merged", "merged", "merged"),
      formula_name = c("base", "interaction", "site_added", "site_interaction"),
      formula_label = c("Base", "Interaction", "Site-added", "Site + interaction"),
      formula_rhs = c(
        "age_s + sex_num",
        "age_s + sex_num + age_s:sex_num",
        "age_s + sex_num + site",
        "age_s + sex_num + age_s:sex_num + site"
      )
    )
  ) %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      formula_id = paste(dataset, formula_name, sep = "__"),
      formula_full = paste("~", formula_rhs),
      uses_site = grepl("\\bsite\\b", formula_rhs),
      uses_age_sex_interaction = grepl("age_s:sex_num", formula_rhs, fixed = TRUE),
      site_branch = if_else(uses_site, "site_adjusted", "site_free"),
      interaction_branch = if_else(uses_age_sex_interaction, "age_sex_interaction", "no_age_sex_interaction")
    ) %>%
    select(
      dataset,
      formula_id,
      formula_name,
      formula_label,
      formula_rhs,
      formula_full,
      uses_site,
      uses_age_sex_interaction,
      site_branch,
      interaction_branch
    )
}

derive_support_tier <- function(dataset_name, horizon_year) {
  h0 <- as.integer(horizon_year)

  if (dataset_name == "PNU") {
    if (h0 == 1L) return("primary_supported")
    if (h0 == 2L) return("sensitivity")
    return("projection")
  }

  if (dataset_name %in% c("SNU", "merged")) {
    if (h0 %in% c(1L, 2L)) return("primary_supported")
    if (h0 <= 5L) return("secondary")
    return("projection")
  }

  stop(sprintf("Unknown dataset `%s` in support-tier derivation.", dataset_name), call. = FALSE)
}

derive_horizon_evidence_class <- function(dataset_name, horizon_year) {
  h0 <- as.integer(horizon_year)

  if (dataset_name == "PNU") {
    if (h0 == 1L) return("directly_observed_data_supported")
    if (h0 == 2L) return("partly_model_dependent")
    return("mostly_extrapolated")
  }

  if (dataset_name %in% c("SNU", "merged")) {
    if (h0 %in% c(1L, 2L)) return("directly_observed_data_supported")
    if (h0 <= 5L) return("partly_model_dependent")
    return("mostly_extrapolated")
  }

  stop(sprintf("Unknown dataset `%s` in horizon-evidence derivation.", dataset_name), call. = FALSE)
}

derive_claim_restriction_flag <- function(horizon_evidence_class) {
  if (is.na(horizon_evidence_class)) {
    return(NA_character_)
  }

  if (horizon_evidence_class == "directly_observed_data_supported") {
    return("primary_claim_allowed")
  }
  if (horizon_evidence_class == "partly_model_dependent") {
    return("secondary_or_sensitivity_only")
  }
  "projection_only"
}

derive_interpretation_note <- function(dataset_name, horizon_year, support_tier, claim_restriction_flag) {
  h0 <- as.integer(horizon_year)

  if (support_tier == "primary_supported") {
    return("Primary supported horizon with comparatively direct follow-up support.")
  }
  if (dataset_name == "PNU" && h0 == 2L) {
    return("Sensitivity horizon for PNU; partly model-dependent and not for primary claims.")
  }
  if (support_tier == "secondary") {
    return("Secondary horizon with growing tail uncertainty; interpret with model-dependence explicitly acknowledged.")
  }
  if (claim_restriction_flag == "projection_only") {
    return("Projection-dominant horizon; mostly extrapolated and not eligible for primary claims.")
  }
  "Common comparison horizon retained for cross-model comparability."
}

build_horizon_registry <- function() {
  tibble::as_tibble(expand.grid(
    dataset = required_datasets,
    horizon_year = common_horizons_year,
    stringsAsFactors = FALSE
  )) %>%
    mutate(
      dataset_key = dataset,
      horizon_id = paste0("year_", horizon_year),
      horizon_days = horizon_year * 365.25,
      support_tier = mapply(derive_support_tier, dataset, horizon_year, USE.NAMES = FALSE),
      horizon_evidence_class = mapply(derive_horizon_evidence_class, dataset, horizon_year, USE.NAMES = FALSE),
      claim_restriction_flag = mapply(derive_claim_restriction_flag, horizon_evidence_class, USE.NAMES = FALSE),
      interpretation_note = mapply(
        derive_interpretation_note,
        dataset,
        horizon_year,
        support_tier,
        claim_restriction_flag,
        USE.NAMES = FALSE
      )
    ) %>%
    arrange(match(dataset, required_datasets), horizon_year)
}

validate_inputs <- function(analysis_datasets, formula_registry, horizon_registry) {
  validate_analysis_datasets(analysis_datasets)

  required_formula_cols <- c(
    "dataset",
    "formula_id",
    "formula_name",
    "formula_label",
    "formula_rhs",
    "site_branch",
    "interaction_branch"
  )
  missing_formula_cols <- setdiff(required_formula_cols, names(formula_registry))
  if (length(missing_formula_cols) > 0L) {
    stop(
      sprintf(
        "Formula registry is missing required columns: %s",
        paste(missing_formula_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  required_horizon_cols <- c(
    "dataset",
    "horizon_year",
    "horizon_days",
    "support_tier",
    "horizon_evidence_class",
    "claim_restriction_flag",
    "interpretation_note"
  )
  missing_horizon_cols <- setdiff(required_horizon_cols, names(horizon_registry))
  if (length(missing_horizon_cols) > 0L) {
    stop(
      sprintf(
        "Horizon registry is missing required columns: %s",
        paste(missing_horizon_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

prepare_dataset <- function(df, dataset_name) {
  required_cols <- c("unique_person_id", "site", "id", "sex_num", "age_exact_entry", "age_s")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("[%s] Dataset is missing required columns: %s", dataset_name, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!("time_year" %in% names(df))) {
    if (!("days_followup" %in% names(df))) {
      stop(sprintf("[%s] Dataset must contain `time_year` or `days_followup`.", dataset_name), call. = FALSE)
    }
    df$time_year <- as.numeric(df$days_followup) / 365.25
  }

  if (!("event_main" %in% names(df))) {
    if (!("status_num" %in% names(df))) {
      stop(sprintf("[%s] Dataset must contain `event_main` or `status_num`.", dataset_name), call. = FALSE)
    }
    df$event_main <- as.integer(as.numeric(df$status_num) == 1)
  }

  out <- df %>%
    mutate(
      unique_person_id = as.character(unique_person_id),
      site = factor(as.character(site)),
      id = as.character(id),
      sex_num = as.integer(sex_num),
      age_exact_entry = as.numeric(age_exact_entry),
      age_s = as.numeric(age_s),
      time_year = as.numeric(time_year),
      time_year_model = pmax(as.numeric(time_year), time_origin_epsilon_year),
      event_main = as.integer(event_main)
    )

  if (nrow(out) != dplyr::n_distinct(out$unique_person_id)) {
    stop(sprintf("[%s] `unique_person_id` must be unique.", dataset_name), call. = FALSE)
  }

  if (anyNA(out$time_year) || any(out$time_year < 0, na.rm = TRUE)) {
    stop(sprintf("[%s] Invalid `time_year` values detected.", dataset_name), call. = FALSE)
  }

  if (anyNA(out$event_main) || any(!out$event_main %in% c(0L, 1L))) {
    stop(sprintf("[%s] `event_main` must be coded as 0/1.", dataset_name), call. = FALSE)
  }

  if (anyNA(out$sex_num) || any(!out$sex_num %in% c(0L, 1L))) {
    stop(sprintf("[%s] `sex_num` must be coded as 0/1.", dataset_name), call. = FALSE)
  }

  out
}

select_model_specs <- function(formula_registry) {
  formula_tbl <- formula_registry %>%
    mutate(dataset = normalize_dataset_label(dataset))

  selected_specs <- dataset_version_registry %>%
    left_join(
      formula_tbl,
      by = c(
        "source_dataset" = "dataset",
        "formula_name" = "formula_name",
        "expected_site_branch" = "site_branch",
        "expected_interaction_branch" = "interaction_branch"
      )
    )

  if (anyNA(selected_specs$formula_id) || anyNA(selected_specs$formula_rhs)) {
    missing_keys <- selected_specs$dataset_version_key[is.na(selected_specs$formula_id) | is.na(selected_specs$formula_rhs)]
    stop(
      sprintf(
        "Could not resolve simplified Stage 5 formulas for dataset versions: %s",
        paste(unique(missing_keys), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  selected_specs %>%
    select(
      dataset_version_key,
      dataset_version_label,
      source_dataset,
      formula_id,
      formula_name,
      formula_label,
      formula_rhs,
      site_branch = expected_site_branch,
      interaction_branch = expected_interaction_branch
    )
}

get_horizon_table <- function(horizon_registry, dataset_name) {
  horizon_registry %>%
    filter(dataset == dataset_name) %>%
    transmute(
      source_dataset = as.character(dataset),
      horizon_year = as.integer(horizon_year),
      horizon_days = as.numeric(horizon_days),
      support_tier = as.character(support_tier),
      horizon_evidence_class = as.character(horizon_evidence_class),
      claim_restriction_flag = as.character(claim_restriction_flag),
      interpretation_note = as.character(interpretation_note)
    ) %>%
    distinct() %>%
    arrange(horizon_year)
}

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

fit_single_lognormal_model <- function(spec_row, analysis_datasets) {
  dataset_name <- spec_row$source_dataset[[1]]
  formula_rhs <- spec_row$formula_rhs[[1]]

  df <- prepare_dataset(analysis_datasets[[dataset_name]], dataset_name)

  surv_formula <- stats::as.formula(
    paste0("survival::Surv(time_year_model, event_main) ~ ", formula_rhs)
  )

  flexsurv::flexsurvreg(
    formula = surv_formula,
    data = df,
    dist = "lnorm"
  )
}

build_model_exports_from_fit <- function(spec_row, fit, analysis_datasets, horizon_registry, curve_times) {
  dataset_name <- spec_row$source_dataset[[1]]
  dataset_version_key <- spec_row$dataset_version_key[[1]]
  dataset_version_label <- spec_row$dataset_version_label[[1]]
  formula_rhs <- spec_row$formula_rhs[[1]]

  df <- prepare_dataset(analysis_datasets[[dataset_name]], dataset_name)
  horizon_tbl <- get_horizon_table(horizon_registry, dataset_name)

  yearly_survival_mat <- predict_parametric_survival(
    fit = fit,
    newdata = df,
    target_times = horizon_tbl$horizon_year
  )
  yearly_survival <- clip_prob(colMeans(yearly_survival_mat, na.rm = TRUE))

  curve_survival_mat <- predict_parametric_survival(
    fit = fit,
    newdata = df,
    target_times = curve_times
  )
  curve_survival <- clip_prob(colMeans(curve_survival_mat, na.rm = TRUE))

  yearly_interval_df <- build_parametric_interval_summaries(
    fit = fit,
    df = df,
    target_times = horizon_tbl$horizon_year,
    n_draws = stage5_ci_draws
  ) %>%
    transmute(
      horizon_year = as.integer(round(time_year)),
      estimated_survival_uncertainty_sd = estimated_survival_uncertainty_sd,
      estimated_survival_lower_95 = estimated_survival_interval_lower_95,
      estimated_survival_median_50 = estimated_survival_interval_median_50,
      estimated_survival_upper_95 = estimated_survival_interval_upper_95,
      estimated_risk_uncertainty_sd = estimated_risk_uncertainty_sd,
      estimated_risk_lower_95 = estimated_risk_interval_lower_95,
      estimated_risk_median_50 = estimated_risk_interval_median_50,
      estimated_risk_upper_95 = estimated_risk_interval_upper_95
    )

  curve_interval_df <- build_parametric_interval_summaries(
    fit = fit,
    df = df,
    target_times = c(0, curve_times),
    n_draws = stage5_ci_draws
  ) %>%
    transmute(
      time_year = as.numeric(time_year),
      estimated_survival_uncertainty_sd = estimated_survival_uncertainty_sd,
      estimated_survival_lower_95 = estimated_survival_interval_lower_95,
      estimated_survival_median_50 = estimated_survival_interval_median_50,
      estimated_survival_upper_95 = estimated_survival_interval_upper_95,
      estimated_risk_uncertainty_sd = estimated_risk_uncertainty_sd,
      estimated_risk_lower_95 = estimated_risk_interval_lower_95,
      estimated_risk_median_50 = estimated_risk_interval_median_50,
      estimated_risk_upper_95 = estimated_risk_interval_upper_95
    )

  yearly_summary <- horizon_tbl %>%
    mutate(
      dataset_version_key = dataset_version_key,
      dataset_version_label = dataset_version_label,
      source_dataset = dataset_name,
      model_family = "lognormal",
      formula_id = spec_row$formula_id[[1]],
      formula_name = spec_row$formula_name[[1]],
      formula_label = spec_row$formula_label[[1]],
      formula_rhs = formula_rhs,
      site_branch = spec_row$site_branch[[1]],
      interaction_branch = spec_row$interaction_branch[[1]],
      interval_type = "95% CI",
      estimated_survival_probability = yearly_survival,
      estimated_risk_probability = 1 - yearly_survival,
      n_subjects_averaged = nrow(df),
      n_transition_events_total = sum(df$event_main == 1L, na.rm = TRUE),
      max_followup_years = max(df$time_year, na.rm = TRUE)
    ) %>%
    left_join(yearly_interval_df, by = "horizon_year") %>%
    select(
      dataset_version_key,
      dataset_version_label,
      source_dataset,
      model_family,
      formula_id,
      formula_name,
      formula_label,
      formula_rhs,
      site_branch,
      interaction_branch,
      interval_type,
      horizon_year,
      horizon_days,
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
      n_subjects_averaged,
      n_transition_events_total,
      max_followup_years,
      support_tier,
      horizon_evidence_class,
      claim_restriction_flag,
      interpretation_note
    )

  curve_df <- tibble(
    dataset_version_key = dataset_version_key,
    dataset_version_label = dataset_version_label,
    source_dataset = dataset_name,
    model_family = "lognormal",
    formula_id = spec_row$formula_id[[1]],
    formula_name = spec_row$formula_name[[1]],
    formula_label = spec_row$formula_label[[1]],
    formula_rhs = formula_rhs,
    site_branch = spec_row$site_branch[[1]],
    interaction_branch = spec_row$interaction_branch[[1]],
    interval_type = "95% CI",
    time_year = c(0, curve_times),
    estimated_survival_probability = c(1, curve_survival),
    estimated_risk_probability = c(0, 1 - curve_survival),
    n_subjects_averaged = nrow(df),
    n_transition_events_total = sum(df$event_main == 1L, na.rm = TRUE)
  ) %>%
    left_join(curve_interval_df, by = "time_year")

  list(
    yearly_summary = yearly_summary,
    curve_df = curve_df
  )
}

validate_fitted_model_cache <- function(fitted_models, model_specs) {
  if (!is.list(fitted_models) || length(fitted_models) != nrow(model_specs)) {
    return(FALSE)
  }

  required_names <- model_specs$dataset_version_key
  if (!identical(sort(names(fitted_models)), sort(required_names))) {
    return(FALSE)
  }

  ordered_models <- fitted_models[required_names]
  all(vapply(ordered_models, inherits, logical(1), what = "flexsurvreg"))
}

get_plot_subset_levels <- function(dataset_version_keys) {
  dataset_version_registry$dataset_version_label[
    match(dataset_version_keys, dataset_version_registry$dataset_version_key)
  ]
}

filter_plot_data <- function(curve_df, yearly_df, dataset_version_keys) {
  label_levels <- get_plot_subset_levels(dataset_version_keys)

  curve_subset <- curve_df %>%
    filter(dataset_version_key %in% dataset_version_keys) %>%
    mutate(dataset_version_label = factor(dataset_version_label, levels = label_levels))

  yearly_subset <- yearly_df %>%
    filter(dataset_version_key %in% dataset_version_keys) %>%
    mutate(dataset_version_label = factor(dataset_version_label, levels = label_levels))

  if (nrow(curve_subset) == 0L || nrow(yearly_subset) == 0L) {
    stop("Requested plot subset produced no rows.", call. = FALSE)
  }

  list(
    curve_df = curve_subset,
    yearly_df = yearly_subset,
    label_levels = label_levels
  )
}

add_pnu_followup_line <- function(plot_object, dataset_version_keys, pnu_followup_end_year) {
  if (!("PNU" %in% dataset_version_keys) || is.na(pnu_followup_end_year) || !is.finite(pnu_followup_end_year)) {
    return(plot_object)
  }

  plot_object +
    geom_vline(
      xintercept = pnu_followup_end_year,
      linewidth = 0.7,
      linetype = "22",
      color = "#6B7280"
    )
}

format_halfyear_labels <- function(x) {
  vapply(
    x,
    FUN.VALUE = character(1),
    FUN = function(tt) {
      if (abs(tt - round(tt)) < 1e-8) {
        as.character(as.integer(round(tt)))
      } else {
        format(tt, nsmall = 1, trim = TRUE)
      }
    }
  )
}

build_dataset_followup_summary <- function(df) {
  time_grid <- seq(0, curve_horizon_max_year, by = risk_table_step_year)

  tibble(
    time_year = time_grid,
    n_at_risk = vapply(
      time_grid,
      FUN.VALUE = integer(1),
      FUN = function(tt) sum(df$time_year >= tt, na.rm = TRUE)
    ),
    n_transition_cum = vapply(
      time_grid,
      FUN.VALUE = integer(1),
      FUN = function(tt) sum(df$event_main == 1L & df$time_year <= tt, na.rm = TRUE)
    )
  )
}

prepare_followup_bar_data <- function(followup_summary_df) {
  bind_rows(
    followup_summary_df %>%
      transmute(
        time_year,
        metric_label = "At risk",
        value = as.numeric(n_at_risk)
      ),
    followup_summary_df %>%
      transmute(
        time_year,
        metric_label = "Cumulative transitions",
        value = as.numeric(n_transition_cum)
      )
  ) %>%
    mutate(
      metric_label = factor(metric_label, levels = c("At risk", "Cumulative transitions"))
    )
}

make_followup_annotation_plot <- function(followup_summary_df, dataset_version_keys, pnu_followup_end_year) {
  bar_df <- prepare_followup_bar_data(followup_summary_df)
  metric_palette <- c(
    "At risk" = "#94A3B8",
    "Cumulative transitions" = "#F59E0B"
  )

  plot_object <- ggplot(bar_df, aes(x = time_year, y = value, fill = metric_label)) +
    geom_col(width = risk_table_step_year * 0.72, show.legend = FALSE) +
    geom_text(
      aes(label = scales::comma(value)),
      vjust = -0.25,
      size = 2.8,
      color = "#111827"
    ) +
    facet_grid(metric_label ~ ., scales = "free_y", switch = "y") +
    scale_fill_manual(values = metric_palette, drop = FALSE) +
    scale_x_continuous(
      breaks = followup_summary_df$time_year,
      labels = format_halfyear_labels(followup_summary_df$time_year),
      limits = c(-risk_table_step_year / 2, curve_horizon_max_year + risk_table_step_year / 2),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      labels = scales::comma,
      expand = expansion(mult = c(0, 0.22))
    ) +
    labs(
      x = "Years after cohort entry (0.5-year intervals)",
      y = "Count"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "#E5E7EB", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 9),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 10),
      legend.position = "none",
      plot.margin = margin(0, 5.5, 5.5, 5.5)
    )

  add_pnu_followup_line(plot_object, dataset_version_keys, pnu_followup_end_year)
}

combine_main_plot_with_followup_table <- function(main_plot, followup_plot) {
  main_plot <- main_plot +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(5.5, 5.5, 0, 5.5)
    )

  patchwork::wrap_plots(
    main_plot,
    followup_plot,
    ncol = 1,
    heights = c(4.4, 2.6)
  )
}

make_probability_plot <- function(
  curve_df,
  yearly_df,
  dataset_version_keys,
  value_col,
  title_text,
  y_label,
  pnu_followup_end_year
) {
  plot_data <- filter_plot_data(curve_df, yearly_df, dataset_version_keys)
  legend_position <- if (length(dataset_version_keys) == 1L) "none" else "top"
  plot_caption <- if ("PNU" %in% dataset_version_keys) {
    sprintf(
      "Dashed vertical line marks the maximum observed follow-up in PNU (%.2f years).",
      pnu_followup_end_year
    )
  } else {
    NULL
  }

  plot_object <- ggplot(
    plot_data$curve_df,
    aes(x = time_year, y = .data[[value_col]], color = dataset_version_label)
  ) +
    geom_line(linewidth = 1.1) +
    geom_point(
      data = plot_data$yearly_df,
      aes(x = horizon_year, y = .data[[value_col]], color = dataset_version_label),
      size = 2
    ) +
    scale_color_manual(values = dataset_palette[plot_data$label_levels], drop = FALSE) +
    scale_x_continuous(
      breaks = common_horizons_year,
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
      subtitle = "Population-averaged predictions from subject-level log-normal no-cure models",
      x = "Years after cohort entry",
      y = y_label,
      color = "Dataset",
      caption = plot_caption
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = legend_position,
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 9)
    )

  add_pnu_followup_line(plot_object, dataset_version_keys, pnu_followup_end_year)
}

make_survival_plot <- function(curve_df, yearly_survival_df, dataset_version_keys, plot_group_label, pnu_followup_end_year) {
  make_probability_plot(
    curve_df = curve_df,
    yearly_df = yearly_survival_df,
    dataset_version_keys = dataset_version_keys,
    value_col = "estimated_survival_probability",
    title_text = paste0(plot_group_label, ": Estimated survival probability after cohort entry"),
    y_label = "Estimated survival probability",
    pnu_followup_end_year = pnu_followup_end_year
  )
}

make_risk_plot <- function(curve_df, yearly_risk_df, dataset_version_keys, plot_group_label, pnu_followup_end_year) {
  make_probability_plot(
    curve_df = curve_df,
    yearly_df = yearly_risk_df,
    dataset_version_keys = dataset_version_keys,
    value_col = "estimated_risk_probability",
    title_text = paste0(plot_group_label, ": Estimated 1 - survival probability after cohort entry"),
    y_label = "Estimated risk probability",
    pnu_followup_end_year = pnu_followup_end_year
  )
}

save_plot_png <- function(plot_object, output_file, width = 10, height = 6, dpi = 300) {
  ggplot2::ggsave(
    filename = output_file,
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    units = "in"
  )
}

save_plot_rds <- function(plot_object, output_file) {
  saveRDS(plot_object, output_file)
}

make_plot_file_base <- function(plot_type, plot_group_key) {
  if (plot_group_key == "all") {
    return(file.path(export_path, paste0(output_prefix, "_", plot_type, "_plot")))
  }

  file.path(
    export_path,
    paste0(output_prefix, "_", plot_type, "_plot_", plot_group_key)
  )
}

make_detail_plot_file_base <- function(plot_type, plot_group_key) {
  file.path(
    export_path,
    paste0(output_prefix, "_", plot_type, "_plot_", plot_group_key, "_with_risk_table")
  )
}

# Read: direct source data and build analysis inputs ===============================
message("Reading source data directly from: ", source_data_file)
analysis_datasets <- build_analysis_datasets_from_source(
  source_file = source_data_file,
  pnu_label = pnu_site_label,
  snu_label = snu_site_label
)
pnu_followup_end_year <- max(analysis_datasets$PNU$time_year, na.rm = TRUE)
formula_registry <- build_formula_registry()
horizon_registry <- build_horizon_registry()

validate_inputs(analysis_datasets, formula_registry, horizon_registry)
model_specs <- select_model_specs(formula_registry)

curve_times <- seq(curve_step_year, curve_horizon_max_year, by = curve_step_year)
reuse_existing_fitted_models <- parse_flag(
  Sys.getenv("STAGE5_SIMPLE_REUSE_FITTED_MODELS", unset = "true"),
  default = TRUE
)
force_refit_fitted_models <- parse_flag(
  Sys.getenv("STAGE5_SIMPLE_FORCE_REFIT", unset = "false"),
  default = FALSE
)

# Fit: simplified Stage 5 log-normal models ===============================
fitted_models <- NULL

if (isTRUE(reuse_existing_fitted_models) && !isTRUE(force_refit_fitted_models) && file.exists(fitted_models_rds_file)) {
  fitted_models <- tryCatch(readRDS(fitted_models_rds_file), error = function(e) NULL)
  if (validate_fitted_model_cache(fitted_models, model_specs)) {
    fitted_models <- fitted_models[model_specs$dataset_version_key]
    message("Loaded cached fitted models from: ", fitted_models_rds_file)
  } else {
    fitted_models <- NULL
    message("Cached fitted model RDS was unavailable or invalid; refitting models.")
  }
}

if (is.null(fitted_models)) {
  fitted_models <- lapply(
    seq_len(nrow(model_specs)),
    function(i) {
      fit_single_lognormal_model(
        spec_row = model_specs[i, , drop = FALSE],
        analysis_datasets = analysis_datasets
      )
    }
  )
  names(fitted_models) <- model_specs$dataset_version_key
  saveRDS(fitted_models, fitted_models_rds_file)
  message("Saved fitted models to: ", fitted_models_rds_file)
}

fit_results <- lapply(
  seq_len(nrow(model_specs)),
  function(i) {
    build_model_exports_from_fit(
      spec_row = model_specs[i, , drop = FALSE],
      fit = fitted_models[[model_specs$dataset_version_key[[i]]]],
      analysis_datasets = analysis_datasets,
      horizon_registry = horizon_registry,
      curve_times = curve_times
    )
  }
)
names(fit_results) <- model_specs$dataset_version_key

yearly_summary <- bind_rows(lapply(fit_results, `[[`, "yearly_summary")) %>%
  mutate(
    dataset_version_label = factor(
      dataset_version_label,
      levels = dataset_version_registry$dataset_version_label
    )
  ) %>%
  arrange(match(dataset_version_key, dataset_version_registry$dataset_version_key), horizon_year) %>%
  mutate(dataset_version_label = as.character(dataset_version_label))

curve_df <- bind_rows(lapply(fit_results, `[[`, "curve_df")) %>%
  mutate(
    dataset_version_label = factor(
      dataset_version_label,
      levels = dataset_version_registry$dataset_version_label
    )
  ) %>%
  arrange(match(dataset_version_key, dataset_version_registry$dataset_version_key), time_year) %>%
  mutate(dataset_version_label = as.character(dataset_version_label))

horizon_summary_df <- yearly_summary

survival_plot_bundle <- lapply(
  names(plot_group_registry),
  function(plot_group_key) {
    plot_group <- plot_group_registry[[plot_group_key]]
    make_survival_plot(
      curve_df = curve_df,
      yearly_survival_df = yearly_summary,
      dataset_version_keys = plot_group$dataset_version_keys,
      plot_group_label = plot_group$label,
      pnu_followup_end_year = pnu_followup_end_year
    )
  }
)
names(survival_plot_bundle) <- names(plot_group_registry)

risk_plot_bundle <- lapply(
  names(plot_group_registry),
  function(plot_group_key) {
    plot_group <- plot_group_registry[[plot_group_key]]
    make_risk_plot(
      curve_df = curve_df,
      yearly_risk_df = yearly_summary,
      dataset_version_keys = plot_group$dataset_version_keys,
      plot_group_label = plot_group$label,
      pnu_followup_end_year = pnu_followup_end_year
    )
  }
)
names(risk_plot_bundle) <- names(plot_group_registry)

detail_plot_group_keys <- setdiff(names(plot_group_registry), "all")
analysis_dataset_lookup <- c(pnu = "PNU", snu = "SNU", merged = "merged")

followup_summary_bundle <- lapply(
  detail_plot_group_keys,
  function(plot_group_key) {
    build_dataset_followup_summary(analysis_datasets[[analysis_dataset_lookup[[plot_group_key]]]])
  }
)
names(followup_summary_bundle) <- detail_plot_group_keys

plot_metadata_df <- yearly_summary %>%
  distinct(
    dataset_version_key,
    dataset_version_label,
    source_dataset,
    model_family,
    formula_id,
    formula_name,
    formula_label,
    formula_rhs,
    site_branch,
    interaction_branch,
    .keep_all = FALSE
  )

followup_plot_source_df <- bind_rows(lapply(
  detail_plot_group_keys,
  function(plot_group_key) {
    dataset_version_key <- plot_group_registry[[plot_group_key]]$dataset_version_keys[[1L]]

    followup_summary_bundle[[plot_group_key]] %>%
      mutate(
        plot_group_key = plot_group_key,
        dataset_version_key = dataset_version_key
      )
  }
)) %>%
  left_join(plot_metadata_df, by = "dataset_version_key") %>%
  transmute(
    component_type = "count_panel",
    plot_group_key = plot_group_key,
    dataset_version_key = dataset_version_key,
    dataset_version_label = dataset_version_label,
    source_dataset = source_dataset,
    model_family = model_family,
    formula_id = formula_id,
    formula_name = formula_name,
    formula_label = formula_label,
    formula_rhs = formula_rhs,
    site_branch = site_branch,
    interaction_branch = interaction_branch,
    interval_type = NA_character_,
    time_horizon_year = as.numeric(time_year),
    horizon_year = as.integer(ifelse(abs(time_year - round(time_year)) < 1e-8, round(time_year), NA)),
    estimated_survival_probability = NA_real_,
    estimated_survival_uncertainty_sd = NA_real_,
    estimated_survival_lower_95 = NA_real_,
    estimated_survival_median_50 = NA_real_,
    estimated_survival_upper_95 = NA_real_,
    estimated_risk_probability = NA_real_,
    estimated_risk_uncertainty_sd = NA_real_,
    estimated_risk_lower_95 = NA_real_,
    estimated_risk_median_50 = NA_real_,
    estimated_risk_upper_95 = NA_real_,
    n_subjects_averaged = NA_integer_,
    n_transition_events_total = NA_integer_,
    n_at_risk = as.integer(n_at_risk),
    n_transition_cumulative = as.integer(n_transition_cum)
  )

plot_source_df <- bind_rows(
  curve_df %>%
    transmute(
      component_type = "curve",
      plot_group_key = NA_character_,
      dataset_version_key = dataset_version_key,
      dataset_version_label = dataset_version_label,
      source_dataset = source_dataset,
      model_family = model_family,
      formula_id = formula_id,
      formula_name = formula_name,
      formula_label = formula_label,
      formula_rhs = formula_rhs,
      site_branch = site_branch,
      interaction_branch = interaction_branch,
      interval_type = as.character(interval_type),
      time_horizon_year = as.numeric(time_year),
      horizon_year = NA_integer_,
      estimated_survival_probability = as.numeric(estimated_survival_probability),
      estimated_survival_uncertainty_sd = as.numeric(estimated_survival_uncertainty_sd),
      estimated_survival_lower_95 = as.numeric(estimated_survival_lower_95),
      estimated_survival_median_50 = as.numeric(estimated_survival_median_50),
      estimated_survival_upper_95 = as.numeric(estimated_survival_upper_95),
      estimated_risk_probability = as.numeric(estimated_risk_probability),
      estimated_risk_uncertainty_sd = as.numeric(estimated_risk_uncertainty_sd),
      estimated_risk_lower_95 = as.numeric(estimated_risk_lower_95),
      estimated_risk_median_50 = as.numeric(estimated_risk_median_50),
      estimated_risk_upper_95 = as.numeric(estimated_risk_upper_95),
      n_subjects_averaged = as.integer(n_subjects_averaged),
      n_transition_events_total = as.integer(n_transition_events_total),
      n_at_risk = NA_integer_,
      n_transition_cumulative = NA_integer_
    ),
  horizon_summary_df %>%
    transmute(
      component_type = "yearly_point",
      plot_group_key = NA_character_,
      dataset_version_key = dataset_version_key,
      dataset_version_label = dataset_version_label,
      source_dataset = source_dataset,
      model_family = model_family,
      formula_id = formula_id,
      formula_name = formula_name,
      formula_label = formula_label,
      formula_rhs = formula_rhs,
      site_branch = site_branch,
      interaction_branch = interaction_branch,
      interval_type = as.character(interval_type),
      time_horizon_year = as.numeric(horizon_year),
      horizon_year = as.integer(horizon_year),
      estimated_survival_probability = as.numeric(estimated_survival_probability),
      estimated_survival_uncertainty_sd = as.numeric(estimated_survival_uncertainty_sd),
      estimated_survival_lower_95 = as.numeric(estimated_survival_lower_95),
      estimated_survival_median_50 = as.numeric(estimated_survival_median_50),
      estimated_survival_upper_95 = as.numeric(estimated_survival_upper_95),
      estimated_risk_probability = as.numeric(estimated_risk_probability),
      estimated_risk_uncertainty_sd = as.numeric(estimated_risk_uncertainty_sd),
      estimated_risk_lower_95 = as.numeric(estimated_risk_lower_95),
      estimated_risk_median_50 = as.numeric(estimated_risk_median_50),
      estimated_risk_upper_95 = as.numeric(estimated_risk_upper_95),
      n_subjects_averaged = as.integer(n_subjects_averaged),
      n_transition_events_total = as.integer(n_transition_events_total),
      n_at_risk = NA_integer_,
      n_transition_cumulative = NA_integer_
    ),
  followup_plot_source_df
) %>%
  mutate(
    component_type = factor(
      component_type,
      levels = c("curve", "yearly_point", "count_panel")
    )
  ) %>%
  arrange(
    match(dataset_version_key, dataset_version_registry$dataset_version_key),
    component_type,
    time_horizon_year
  ) %>%
  mutate(component_type = as.character(component_type))

detail_survival_plot_bundle <- lapply(
  detail_plot_group_keys,
  function(plot_group_key) {
    combine_main_plot_with_followup_table(
      main_plot = survival_plot_bundle[[plot_group_key]],
      followup_plot = make_followup_annotation_plot(
        followup_summary_df = followup_summary_bundle[[plot_group_key]],
        dataset_version_keys = plot_group_registry[[plot_group_key]]$dataset_version_keys,
        pnu_followup_end_year = pnu_followup_end_year
      )
    )
  }
)
names(detail_survival_plot_bundle) <- detail_plot_group_keys

detail_risk_plot_bundle <- lapply(
  detail_plot_group_keys,
  function(plot_group_key) {
    combine_main_plot_with_followup_table(
      main_plot = risk_plot_bundle[[plot_group_key]],
      followup_plot = make_followup_annotation_plot(
        followup_summary_df = followup_summary_bundle[[plot_group_key]],
        dataset_version_keys = plot_group_registry[[plot_group_key]]$dataset_version_keys,
        pnu_followup_end_year = pnu_followup_end_year
      )
    )
  }
)
names(detail_risk_plot_bundle) <- detail_plot_group_keys

survival_plot_obj <- survival_plot_bundle[["all"]]
risk_plot_obj <- risk_plot_bundle[["all"]]

# Export: simplified outputs ===============================
cleanup_existing_stage5_outputs(export_path)
saveRDS(fitted_models, fitted_models_rds_file)
readr::write_csv(horizon_summary_df, horizon_summary_file)
readr::write_csv(plot_source_df, plot_source_file)
save_plot_rds(survival_plot_obj, survival_plot_rds_file)
save_plot_rds(risk_plot_obj, risk_plot_rds_file)
save_plot_png(survival_plot_obj, survival_plot_png_file)
save_plot_png(risk_plot_obj, risk_plot_png_file)
save_plot_rds(
  list(
    survival = survival_plot_bundle,
    risk = risk_plot_bundle,
    pnu_followup_end_year = pnu_followup_end_year
  ),
  plot_bundle_rds_file
)
save_plot_rds(
  list(
    survival = detail_survival_plot_bundle,
    risk = detail_risk_plot_bundle,
    followup_summary = followup_summary_bundle
  ),
  detail_plot_bundle_rds_file
)

for (plot_group_key in setdiff(names(plot_group_registry), "all")) {
  survival_file_base <- make_plot_file_base("survival", plot_group_key)
  risk_file_base <- make_plot_file_base("risk", plot_group_key)

  save_plot_rds(survival_plot_bundle[[plot_group_key]], paste0(survival_file_base, ".rds"))
  save_plot_rds(risk_plot_bundle[[plot_group_key]], paste0(risk_file_base, ".rds"))
  save_plot_png(survival_plot_bundle[[plot_group_key]], paste0(survival_file_base, ".png"))
  save_plot_png(risk_plot_bundle[[plot_group_key]], paste0(risk_file_base, ".png"))

  detail_survival_file_base <- make_detail_plot_file_base("survival", plot_group_key)
  detail_risk_file_base <- make_detail_plot_file_base("risk", plot_group_key)

  save_plot_rds(detail_survival_plot_bundle[[plot_group_key]], paste0(detail_survival_file_base, ".rds"))
  save_plot_rds(detail_risk_plot_bundle[[plot_group_key]], paste0(detail_risk_file_base, ".rds"))
  save_plot_png(
    detail_survival_plot_bundle[[plot_group_key]],
    paste0(detail_survival_file_base, ".png"),
    width = detail_plot_width_in,
    height = detail_plot_height_in
  )
  save_plot_png(
    detail_risk_plot_bundle[[plot_group_key]],
    paste0(detail_risk_file_base, ".png"),
    width = detail_plot_width_in,
    height = detail_plot_height_in
  )
}

organize_stage5_outputs(export_path)

message("Saved simplified no-cure log-normal outputs to: ", export_path)
