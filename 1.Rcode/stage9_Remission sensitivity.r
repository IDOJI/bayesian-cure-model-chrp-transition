# 🔴 Configure: paths and stage-9 controls ===============================
stage1_path <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/stage1_Backbone lock'
export_path <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/Stage9_Remission sensitivity'

stage1_bundle_file <- file.path(stage1_path, "stage1_backbone_bundle.rds")
stage1_datasets_file <- file.path(stage1_path, "stage1_analysis_datasets.rds")
stage1_formula_registry_file <- file.path(stage1_path, "stage1_formula_registry.csv")
stage1_horizon_registry_file <- file.path(stage1_path, "stage1_horizon_registry.csv")
stage1_threshold_registry_file <- file.path(stage1_path, "stage1_threshold_registry.csv")
stage1_modeling_registry_file <- file.path(stage1_path, "stage1_modeling_registry.csv")
stage1_metadata_registry_file <- file.path(stage1_path, "stage1_metadata_registry.csv")
stage1_export_manifest_file <- file.path(stage1_path, "stage1_export_manifest.csv")

run_merged_site_adjusted_benchmark <- TRUE
run_cause_specific_models <- TRUE
run_subdistribution_models <- TRUE
save_fit_objects <- TRUE
make_plot_pdf <- TRUE
make_plot_png <- TRUE
reuse_existing_valid_outputs <- TRUE
overwrite_existing_plots <- TRUE

minimum_transition_events <- 5L
minimum_remission_events <- 2L
cox_ties <- "efron"
subdistribution_failcode <- 1L
subdistribution_cencode <- 0L

main_risk_scale <- "transition_only_main"
supplementary_risk_scale <- "transition_cif_competing"

plot_png_width <- 11
plot_png_height <- 8.5
plot_png_dpi <- 300
projection_plot_caption <- "Grey shading marks projection-dominant horizons; prioritize rows with claim_restriction_flag != 'projection_only' for primary or supported interpretation."

# 🔴 Initialize: packages, options, and runtime logging ===============================
required_packages <- c("survival", "dplyr", "readr", "tibble", "tidyr", "purrr", "ggplot2")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop(
    sprintf("Install required packages before running Stage 9: %s", paste(missing_packages, collapse = ", ")),
    call. = FALSE
  )
}

has_cmprsk <- requireNamespace("cmprsk", quietly = TRUE)

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE, scipen = 999)

if (!dir.exists(export_path) && !dir.create(export_path, recursive = TRUE, showWarnings = FALSE)) {
  stop(sprintf("Could not create export_path: %s", export_path), call. = FALSE)
}

progress_log_file <- file.path(export_path, "stage9_progress_log.txt")
stage_start_time <- Sys.time()

format_elapsed <- function(start_time) {
  elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (!is.finite(elapsed_sec)) {
    return("NA sec")
  }
  if (elapsed_sec < 60) {
    return(sprintf("%.1f sec", elapsed_sec))
  }
  if (elapsed_sec < 3600) {
    return(sprintf("%.1f min", elapsed_sec / 60))
  }
  sprintf("%.2f hr", elapsed_sec / 3600)
}

log_step <- function(..., level = "INFO") {
  line <- sprintf(
    "[%s] [%s] %s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    level,
    paste0(..., collapse = "")
  )
  message(line)
  cat(line, "\n", file = progress_log_file, append = TRUE)
  flush.console()
  invisible(line)
}

writeLines(character(0), progress_log_file)
log_step("Stage 9 started.")
if (isTRUE(run_subdistribution_models) && !isTRUE(has_cmprsk)) {
  log_step("cmprsk package not available; Fine-Gray subdistribution sensitivity will be skipped.", level = "WARN")
}


# 🔴 Define: scalar helpers, readers, and validators ===============================
## 🟠 Define: scalar helpers ===============================
normalize_existing_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

safe_divide <- function(num, den) {
  ifelse(is.na(den) | den <= 0, NA_real_, num / den)
}

collapse_unique_chr <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA_character_)
  }
  paste(sort(x), collapse = "|")
}

summarise_single_logical <- function(x, label = "logical value") {
  x <- unique(as.logical(x))
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA)
  }
  if (length(x) > 1L) {
    stop(sprintf("Expected a single `%s`, observed: %s", label, paste(x, collapse = ", ")), call. = FALSE)
  }
  x[[1]]
}

summarise_single_text <- function(x, label = "text value") {
  x <- unique(as.character(x))
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA_character_)
  }
  if (length(x) > 1L) {
    stop(sprintf("Expected a single `%s`, observed: %s", label, paste(x, collapse = ", ")), call. = FALSE)
  }
  x[[1]]
}

equal_or_both_na <- function(x, y) {
  (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
}

step_eval <- function(times, values, at, left_value = 0) {
  if (length(times) == 0L || length(values) == 0L) {
    return(rep(left_value, length(at)))
  }
  idx <- findInterval(at, times)
  out <- rep(left_value, length(at))
  keep <- idx > 0L
  out[keep] <- values[idx[keep]]
  out
}

compact_bind_rows <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) == 0L) {
    tibble::tibble()
  } else {
    dplyr::bind_rows(x)
  }
}

safe_write_csv <- function(df, path) {
  readr::write_csv(df, path)
  invisible(path)
}

safe_save_rds <- function(object, path) {
  saveRDS(object, path)
  invisible(path)
}

build_file_signature <- function(paths) {
  tibble(path = normalize_existing_path(paths)) %>%
    mutate(
      exists = file.exists(path),
      size_bytes = ifelse(exists, file.info(path)$size, NA_real_),
      modified_time = ifelse(exists, as.character(file.info(path)$mtime), NA_character_)
    )
}

## 🟠 Define: registry harmonizers ===============================
derive_support_tier <- function(dataset, horizon_year) {
  h <- as.integer(horizon_year)

  dplyr::case_when(
    dataset == "PNU" & h == 1L ~ "primary_supported",
    dataset == "PNU" & h == 2L ~ "sensitivity",
    dataset == "PNU" & h >= 3L ~ "projection",
    dataset %in% c("SNU", "merged") & h %in% c(1L, 2L) ~ "primary_supported",
    dataset %in% c("SNU", "merged") & h %in% 3L:5L ~ "secondary",
    dataset %in% c("SNU", "merged") & h >= 6L ~ "projection",
    TRUE ~ "projection"
  )
}

derive_horizon_evidence_class <- function(dataset, horizon_year) {
  h <- as.integer(horizon_year)

  dplyr::case_when(
    dataset == "PNU" & h == 1L ~ "directly_observed_data_supported",
    dataset == "PNU" & h == 2L ~ "partly_model_dependent",
    dataset == "PNU" & h >= 3L ~ "mostly_extrapolated",
    dataset %in% c("SNU", "merged") & h %in% c(1L, 2L) ~ "directly_observed_data_supported",
    dataset %in% c("SNU", "merged") & h %in% 3L:5L ~ "partly_model_dependent",
    dataset %in% c("SNU", "merged") & h >= 6L ~ "mostly_extrapolated",
    TRUE ~ "mostly_extrapolated"
  )
}

derive_claim_restriction_flag <- function(horizon_evidence_class, support_tier) {
  dplyr::case_when(
    horizon_evidence_class == "directly_observed_data_supported" ~ "primary_claim_allowed",
    horizon_evidence_class == "partly_model_dependent" ~ "secondary_or_sensitivity_only",
    horizon_evidence_class == "mostly_extrapolated" ~ "projection_only",
    support_tier == "projection" ~ "projection_only",
    TRUE ~ "secondary_or_sensitivity_only"
  )
}

derive_legacy_interpretation_tier <- function(support_tier) {
  dplyr::case_when(
    support_tier == "primary_supported" ~ "primary_supported",
    support_tier == "secondary" ~ "secondary",
    support_tier == "sensitivity" ~ "sensitivity",
    support_tier == "projection" ~ "projection",
    TRUE ~ "secondary"
  )
}

derive_interpretation_note <- function(dataset, horizon_year, support_tier, horizon_evidence_class, claim_restriction_flag) {
  h <- as.integer(horizon_year)

  if (support_tier == "primary_supported" && horizon_evidence_class == "directly_observed_data_supported") {
    return("Primary supported horizon with comparatively direct follow-up support.")
  }

  if (dataset == "PNU" && h == 2L) {
    return("Sensitivity horizon for PNU; partly model-dependent and not for primary claims.")
  }

  if (support_tier == "secondary" && horizon_evidence_class == "partly_model_dependent") {
    return("Secondary horizon with growing tail uncertainty; interpret with model dependence acknowledged.")
  }

  if (claim_restriction_flag == "projection_only") {
    return("Projection-dominant horizon; mostly extrapolated and not eligible for primary claims.")
  }

  "Common comparison horizon retained for cross-model comparability."
}

normalize_horizon_registry <- function(horizon_registry) {
  assert_required_columns(horizon_registry, c("dataset", "horizon_year"), "normalize_horizon_registry")

  out <- horizon_registry %>%
    mutate(
      dataset = as.character(dataset),
      horizon_year = as.integer(horizon_year),
      horizon_id = if ("horizon_id" %in% names(.)) as.character(horizon_id) else paste0("year_", horizon_year),
      horizon_days = if ("horizon_days" %in% names(.)) as.numeric(horizon_days) else horizon_year * 365.25,
      support_tier = if ("support_tier" %in% names(.)) as.character(support_tier) else derive_support_tier(dataset, horizon_year),
      support_tier_standard = if ("support_tier_standard" %in% names(.)) as.character(support_tier_standard) else support_tier,
      interpretation_tier = if ("interpretation_tier" %in% names(.)) as.character(interpretation_tier) else derive_legacy_interpretation_tier(support_tier),
      primary_supported_flag = if ("primary_supported_flag" %in% names(.)) as.logical(primary_supported_flag) else support_tier == "primary_supported",
      horizon_evidence_class = if ("horizon_evidence_class" %in% names(.)) as.character(horizon_evidence_class) else derive_horizon_evidence_class(dataset, horizon_year),
      claim_restriction_flag = if ("claim_restriction_flag" %in% names(.)) as.character(claim_restriction_flag) else derive_claim_restriction_flag(horizon_evidence_class, support_tier),
      interpretation_note = if ("interpretation_note" %in% names(.)) as.character(interpretation_note) else mapply(
        derive_interpretation_note,
        dataset,
        horizon_year,
        support_tier,
        horizon_evidence_class,
        claim_restriction_flag,
        USE.NAMES = FALSE
      )
    ) %>%
    mutate(
      horizon_support_label = support_tier,
      projection_flag = horizon_evidence_class == "mostly_extrapolated" | claim_restriction_flag %in% c("projection_only", "projection_plus_prior_sensitive"),
      supported_for_reporting = !projection_flag,
      reporting_status = ifelse(projection_flag, "projection_only", "supported")
    ) %>%
    arrange(match(dataset, c("PNU", "SNU", "merged")), horizon_year)

  out
}

## 🟠 Define: file readers ===============================
read_stage1_csv <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Required Stage 1 CSV file does not exist: %s", path), call. = FALSE)
  }
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

read_stage1_rds <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Required Stage 1 RDS file does not exist: %s", path), call. = FALSE)
  }
  readRDS(path)
}

safe_read_existing_csv <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
}

safe_read_existing_rds <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(readRDS(path), error = function(e) NULL)
}

## 🟠 Define: integrity helpers ===============================
assert_required_columns <- function(df, required_cols, context_label) {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(sprintf("[%s] Missing required columns: %s", context_label, paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

assert_prob_in_unit_interval <- function(x, label, tol = 1e-8) {
  bad <- which(!is.na(x) & (x < -tol | x > 1 + tol))
  if (length(bad) > 0L) {
    stop(sprintf("%s has %d values outside [0,1] tolerance.", label, length(bad)), call. = FALSE)
  }
  invisible(TRUE)
}

assert_prob_sum_close_to_one <- function(x, label, tol = 1e-6) {
  bad <- which(!is.na(x) & abs(x - 1) > tol)
  if (length(bad) > 0L) {
    stop(sprintf("%s has %d values not summing to 1 within tolerance.", label, length(bad)), call. = FALSE)
  }
  invisible(TRUE)
}

prepare_stage1_dataset <- function(df, dataset_name) {
  required_cols <- c(
    "unique_person_id", "id", "site", "sex_num", "age_exact_entry", "age_s",
    "days_followup", "time_year", "status_num", "event_main", "remission_flag", "censor_main"
  )
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(sprintf("[%s] Missing required Stage 1 columns: %s", dataset_name, paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  sex_label_source <- if ("sex_label" %in% names(df)) as.character(df$sex_label) else if ("sex_fact" %in% names(df)) as.character(df$sex_fact) else rep(NA_character_, nrow(df))
  status_label_source <- if ("status_label" %in% names(df)) as.character(df$status_label) else if ("status" %in% names(df)) as.character(df$status) else rep(NA_character_, nrow(df))

  out <- df %>%
    mutate(
      unique_person_id = as.character(unique_person_id),
      id = as.character(id),
      site = as.character(site),
      sex_num = as.integer(sex_num),
      sex_label = dplyr::case_when(
        !is.na(sex_label_source) ~ sex_label_source,
        sex_num == 0L ~ "Female",
        sex_num == 1L ~ "Male",
        TRUE ~ NA_character_
      ),
      age_exact_entry = as.numeric(age_exact_entry),
      age_s = as.numeric(age_s),
      days_followup = as.numeric(days_followup),
      time_year = as.numeric(time_year),
      status_num = as.integer(status_num),
      status_label = dplyr::case_when(
        !is.na(status_label_source) ~ status_label_source,
        status_num == 0L ~ "right_censoring",
        status_num == 1L ~ "transition",
        status_num == 2L ~ "remission",
        TRUE ~ NA_character_
      ),
      event_main = as.integer(event_main),
      remission_flag = as.integer(remission_flag),
      censor_main = as.integer(censor_main)
    )

  if (nrow(out) != dplyr::n_distinct(out$unique_person_id)) {
    stop(sprintf("[%s] `unique_person_id` must be unique.", dataset_name), call. = FALSE)
  }
  if (anyNA(out[, required_cols])) {
    stop(sprintf("[%s] Missing values detected in required Stage 1 columns.", dataset_name), call. = FALSE)
  }
  if (any(out$days_followup < 0, na.rm = TRUE)) {
    stop(sprintf("[%s] Negative `days_followup` detected.", dataset_name), call. = FALSE)
  }
  if (any(!out$sex_num %in% c(0L, 1L), na.rm = TRUE)) {
    stop(sprintf("[%s] `sex_num` must remain coded as 0/1.", dataset_name), call. = FALSE)
  }
  if (any(!out$status_num %in% c(0L, 1L, 2L), na.rm = TRUE)) {
    stop(sprintf("[%s] `status_num` must remain coded as 0/1/2.", dataset_name), call. = FALSE)
  }

  out
}

validate_stage1_backbone <- function(stage1_bundle, analysis_datasets, formula_registry, horizon_registry, threshold_registry, modeling_registry, metadata_registry) {
  required_dataset_names <- c("PNU", "SNU", "merged")
  observed_dataset_names <- names(analysis_datasets)

  if (!all(required_dataset_names %in% observed_dataset_names)) {
    stop(
      sprintf("Stage 1 analysis datasets must include PNU, SNU, merged. Observed: %s", paste(observed_dataset_names, collapse = ", ")),
      call. = FALSE
    )
  }

  invisible(lapply(required_dataset_names, function(dataset_name) prepare_stage1_dataset(analysis_datasets[[dataset_name]], dataset_name)))

  observed_horizons <- sort(unique(as.integer(horizon_registry$horizon_year)))
  if (!identical(observed_horizons, 1:10)) {
    stop("Stage 1 horizon registry must contain the common 1:10 year horizon grid.", call. = FALSE)
  }

  observed_thresholds <- sort(unique(as.numeric(threshold_registry$threshold)))
  if (length(observed_thresholds) == 0L || anyNA(observed_thresholds) || any(observed_thresholds <= 0 | observed_thresholds >= 1)) {
    stop("Stage 1 threshold registry must contain unique non-missing probabilities strictly between 0 and 1.", call. = FALSE)
  }

  required_formula_cols <- c(
    "dataset", "formula_id", "formula_name", "formula_label", "formula_rhs",
    "site_branch", "interaction_branch"
  )
  if (length(setdiff(required_formula_cols, names(formula_registry))) > 0L) {
    stop("Stage 1 formula registry is missing required columns.", call. = FALSE)
  }

  if (length(setdiff(c("dataset", "formula_id"), names(modeling_registry))) > 0L) {
    stop("Stage 1 modeling registry is missing required columns.", call. = FALSE)
  }

  if (length(setdiff(c("metadata_group", "metadata_name", "metadata_value"), names(metadata_registry))) > 0L) {
    stop("Stage 1 metadata registry is missing required columns.", call. = FALSE)
  }

  bundle_horizons <- sort(unique(as.integer(stage1_bundle$config$common_horizons_year)))
  if (!identical(bundle_horizons, observed_horizons)) {
    stop("Stage 1 bundle horizon vector does not match Stage 1 horizon registry.", call. = FALSE)
  }

  bundle_thresholds <- sort(unique(as.numeric(stage1_bundle$config$risk_thresholds)))
  if (!identical(bundle_thresholds, observed_thresholds)) {
    stop("Stage 1 bundle threshold vector does not match Stage 1 threshold registry.", call. = FALSE)
  }

  invisible(TRUE)
}


# 🔴 Define: component annotations and output validators ===============================
## 🟠 Define: component annotation rules ===============================
derive_component_annotations <- function(comparison_scope, comparison_family, analysis_variant, sensitivity_method = NA_character_) {
  comparison_scope <- as.character(comparison_scope)
  comparison_family <- as.character(comparison_family)
  analysis_variant <- as.character(analysis_variant)
  sensitivity_method <- as.character(sensitivity_method)

  tibble(
    remission_component_available = dplyr::case_when(
      analysis_variant == "main_censoring_mirror" ~ FALSE,
      comparison_family == "subdistribution" & analysis_variant == "remission_sensitive" ~ FALSE,
      comparison_family == "cause_specific" &
        analysis_variant == "remission_sensitive" &
        sensitivity_method == "cause_specific_cox_transition_only" ~ FALSE,
      analysis_variant == "remission_sensitive" ~ TRUE,
      TRUE ~ FALSE
    ),
    eventfree_component_available = dplyr::case_when(
      analysis_variant == "main_censoring_mirror" ~ TRUE,
      comparison_family == "subdistribution" & analysis_variant == "remission_sensitive" ~ FALSE,
      analysis_variant == "remission_sensitive" ~ TRUE,
      TRUE ~ FALSE
    ),
    component_mode = dplyr::case_when(
      analysis_variant == "main_censoring_mirror" ~ "transition_plus_eventfree_main_mirror",
      comparison_scope == "benchmark" & analysis_variant == "remission_sensitive" ~ "nonparametric_transition_remission_eventfree",
      comparison_family == "cause_specific" &
        analysis_variant == "remission_sensitive" &
        sensitivity_method == "cause_specific_cox_transition_only" ~ "transition_only_cause_specific_model",
      comparison_family == "cause_specific" & analysis_variant == "remission_sensitive" ~ "full_cause_specific_decomposition",
      comparison_family == "subdistribution" & analysis_variant == "remission_sensitive" ~ "transition_only_subdistribution_model",
      analysis_variant == "remission_sensitive" ~ "transition_remission_eventfree_decomposition_available",
      TRUE ~ "transition_plus_eventfree_main_mirror"
    ),
    component_note = dplyr::case_when(
      analysis_variant == "main_censoring_mirror" ~ "transition_plus_eventfree_main_mirror",
      comparison_family == "cause_specific" &
        analysis_variant == "remission_sensitive" &
        sensitivity_method == "cause_specific_cox_transition_only" ~ "transition_only_cause_specific_model",
      comparison_family == "subdistribution" & analysis_variant == "remission_sensitive" ~ "transition_only_subdistribution_model",
      analysis_variant == "remission_sensitive" ~ "transition_remission_eventfree_decomposition_available",
      TRUE ~ "transition_plus_eventfree_main_mirror"
    )
  )
}

bind_component_annotations <- function(df) {
  dplyr::bind_cols(
    df,
    derive_component_annotations(
      comparison_scope = df$comparison_scope,
      comparison_family = df$comparison_family,
      analysis_variant = df$analysis_variant,
      sensitivity_method = df$sensitivity_method
    )
  )
}

validate_component_annotations <- function(tbl, table_label) {
  if (nrow(tbl) == 0L) {
    return(invisible(TRUE))
  }

  expected <- derive_component_annotations(
    comparison_scope = tbl$comparison_scope,
    comparison_family = tbl$comparison_family,
    analysis_variant = tbl$analysis_variant,
    sensitivity_method = tbl$sensitivity_method
  )

  if (any(!equal_or_both_na(tbl$remission_component_available, expected$remission_component_available), na.rm = TRUE)) {
    stop(sprintf("[%s] remission_component_available does not match model structure.", table_label), call. = FALSE)
  }
  if (any(!equal_or_both_na(tbl$eventfree_component_available, expected$eventfree_component_available), na.rm = TRUE)) {
    stop(sprintf("[%s] eventfree_component_available does not match model structure.", table_label), call. = FALSE)
  }
  if (any(!equal_or_both_na(tbl$component_mode, expected$component_mode), na.rm = TRUE)) {
    stop(sprintf("[%s] component_mode does not match model structure.", table_label), call. = FALSE)
  }
  if (any(!equal_or_both_na(tbl$component_note, expected$component_note), na.rm = TRUE)) {
    stop(sprintf("[%s] component_note does not match model structure.", table_label), call. = FALSE)
  }

  invisible(TRUE)
}

## 🟠 Define: horizon reporting checks ===============================
validate_horizon_reporting_fields <- function(tbl, table_label) {
  if (nrow(tbl) == 0L) {
    return(invisible(TRUE))
  }

  required_cols <- c(
    "support_tier", "horizon_evidence_class", "claim_restriction_flag",
    "horizon_support_label", "supported_for_reporting", "projection_flag", "reporting_status"
  )
  assert_required_columns(tbl, required_cols, paste0("validate_horizon_reporting_fields:", table_label))

  expected_projection <- tbl$horizon_evidence_class == "mostly_extrapolated" |
    tbl$claim_restriction_flag %in% c("projection_only", "projection_plus_prior_sensitive")

  if (any(!equal_or_both_na(tbl$projection_flag, expected_projection), na.rm = TRUE)) {
    stop(sprintf("[%s] projection_flag inconsistent with horizon_evidence_class / claim_restriction_flag.", table_label), call. = FALSE)
  }

  invisible(TRUE)
}

## 🟠 Define: output-specific validators ===============================
validate_status_summary <- function(status_summary) {
  if (nrow(status_summary) == 0L) {
    return(invisible(TRUE))
  }
  required_cols <- c(
    "summary_level", "dataset", "group_variable", "group_value",
    "n_subjects", "n_status_0_right_censor", "n_status_1_transition", "n_status_2_remission"
  )
  assert_required_columns(status_summary, required_cols, "status_summary")
  invisible(TRUE)
}

validate_benchmark_state_probabilities <- function(tbl) {
  if (nrow(tbl) == 0L) {
    return(invisible(TRUE))
  }
  required_cols <- c(
    "dataset", "comparison_scope", "comparison_family", "comparison_id", "risk_scale",
    "horizon_year", "main_km_risk", "sensitivity_transition_risk", "sensitivity_remission_risk",
    "support_tier", "horizon_evidence_class", "claim_restriction_flag"
  )
  assert_required_columns(tbl, required_cols, "benchmark_state_probabilities")
  validate_horizon_reporting_fields(tbl, "benchmark_state_probabilities")
  if (!all(tbl$risk_scale %in% supplementary_risk_scale)) {
    stop("benchmark_state_probabilities must keep risk_scale == supplementary_risk_scale.", call. = FALSE)
  }
  assert_prob_in_unit_interval(tbl$main_km_risk, "benchmark_state_probabilities main_km_risk")
  assert_prob_in_unit_interval(tbl$sensitivity_transition_risk, "benchmark_state_probabilities sensitivity_transition_risk")
  assert_prob_in_unit_interval(tbl$sensitivity_remission_risk, "benchmark_state_probabilities sensitivity_remission_risk")
  invisible(TRUE)
}

validate_subject_predictions <- function(predictions_long, block_name) {
  if (nrow(predictions_long) == 0L) {
    return(invisible(TRUE))
  }

  required_cols <- c(
    "dataset", "comparison_scope", "comparison_family", "comparison_id",
    "prediction_level", "model_class", "model_family", "risk_scale",
    "analysis_variant", "evaluation_mode", "risk_estimator",
    "predicted_transition_risk", "predicted_remission_risk", "predicted_eventfree_prob",
    "risk_difference_from_main_signed", "risk_difference_from_main_absolute",
    "support_tier", "horizon_evidence_class", "claim_restriction_flag",
    "horizon_support_label", "supported_for_reporting", "projection_flag", "reporting_status",
    "remission_component_available", "eventfree_component_available", "component_mode", "component_note"
  )
  assert_required_columns(predictions_long, required_cols, paste0("validate_subject_predictions:", block_name))
  validate_horizon_reporting_fields(predictions_long, paste0("subject_predictions:", block_name))
  validate_component_annotations(predictions_long, paste0("subject_predictions:", block_name))

  if (!all(predictions_long$risk_scale %in% supplementary_risk_scale)) {
    stop(sprintf("[%s] risk_scale must remain supplementary remission-sensitive scale in Stage 9 outputs.", block_name), call. = FALSE)
  }

  assert_prob_in_unit_interval(predictions_long$predicted_transition_risk, paste0(block_name, " predicted_transition_risk"))
  if (!all(is.na(predictions_long$predicted_remission_risk))) {
    assert_prob_in_unit_interval(predictions_long$predicted_remission_risk, paste0(block_name, " predicted_remission_risk"))
  }
  if (!all(is.na(predictions_long$predicted_eventfree_prob))) {
    assert_prob_in_unit_interval(predictions_long$predicted_eventfree_prob, paste0(block_name, " predicted_eventfree_prob"))
  }

  main_rows <- predictions_long %>% filter(analysis_variant == "main_censoring_mirror")
  if (nrow(main_rows) > 0L) {
    if (anyNA(main_rows$predicted_eventfree_prob)) {
      stop(sprintf("[%s] Main mirror rows must have non-missing event-free probability.", block_name), call. = FALSE)
    }
    assert_prob_sum_close_to_one(
      main_rows$predicted_transition_risk + main_rows$predicted_eventfree_prob,
      paste0(block_name, " main transition + eventfree")
    )
    if (any(main_rows$predicted_remission_risk != 0, na.rm = TRUE)) {
      stop(sprintf("[%s] Main mirror rows must keep predicted_remission_risk == 0.", block_name), call. = FALSE)
    }
  }

  invisible(TRUE)
}

validate_risk_summary_table <- function(summary_tbl, table_label) {
  if (nrow(summary_tbl) == 0L) {
    return(invisible(TRUE))
  }

  required_cols <- c(
    "dataset", "comparison_scope", "comparison_family", "comparison_id",
    "model_class", "model_family", "risk_scale", "horizon_year",
    "support_tier", "horizon_evidence_class", "claim_restriction_flag",
    "mean_main_transition_risk", "mean_sensitivity_transition_risk",
    "mean_sensitivity_remission_risk", "mean_sensitivity_eventfree_prob"
  )
  assert_required_columns(summary_tbl, required_cols, table_label)
  validate_horizon_reporting_fields(summary_tbl, table_label)
  if (!all(summary_tbl$risk_scale %in% supplementary_risk_scale)) {
    stop(sprintf("[%s] risk_scale must remain supplementary remission-sensitive scale.", table_label), call. = FALSE)
  }
  assert_prob_in_unit_interval(summary_tbl$mean_main_transition_risk, paste0(table_label, " mean_main_transition_risk"))
  assert_prob_in_unit_interval(summary_tbl$mean_sensitivity_transition_risk, paste0(table_label, " mean_sensitivity_transition_risk"))
  if (!all(is.na(summary_tbl$mean_sensitivity_remission_risk))) {
    assert_prob_in_unit_interval(summary_tbl$mean_sensitivity_remission_risk, paste0(table_label, " mean_sensitivity_remission_risk"))
  }
  if (!all(is.na(summary_tbl$mean_sensitivity_eventfree_prob))) {
    assert_prob_in_unit_interval(summary_tbl$mean_sensitivity_eventfree_prob, paste0(table_label, " mean_sensitivity_eventfree_prob"))
  }
  invisible(TRUE)
}

validate_classification_table <- function(classification_tbl, table_label) {
  if (nrow(classification_tbl) == 0L) {
    return(invisible(TRUE))
  }

  required_cols <- c(
    "n_subjects", "n_event_by_horizon", "n_remission_by_horizon",
    "n_known_nonevent_by_horizon", "n_evaluable_by_horizon", "n_unknown_due_to_early_followup_loss",
    "n_predicted_positive", "n_predicted_negative",
    "n_predicted_positive_evaluable", "n_predicted_negative_evaluable",
    "tp_count", "fp_count", "tn_count", "fn_count",
    "positive_classification_rate", "positive_classification_rate_evaluable",
    "false_positive_rate", "false_positive_burden", "specificity", "ppv", "tpr",
    "event_prevalence_all", "event_prevalence_evaluable", "remission_prevalence_all",
    "net_benefit", "net_benefit_evaluable", "standardized_net_benefit",
    "support_tier", "horizon_evidence_class", "claim_restriction_flag",
    "risk_scale", "remission_component_available", "eventfree_component_available", "component_mode", "component_note"
  )
  assert_required_columns(classification_tbl, required_cols, paste0("validate_classification_table:", table_label))
  validate_horizon_reporting_fields(classification_tbl, paste0("classification:", table_label))
  validate_component_annotations(classification_tbl, paste0("classification:", table_label))

  if (!all(classification_tbl$risk_scale %in% supplementary_risk_scale)) {
    stop(sprintf("[%s] risk_scale must remain supplementary remission-sensitive scale.", table_label), call. = FALSE)
  }

  if (any(classification_tbl$n_subjects != classification_tbl$n_predicted_positive + classification_tbl$n_predicted_negative)) {
    stop(sprintf("[%s] n_subjects != predicted positive + negative.", table_label), call. = FALSE)
  }
  if (any(classification_tbl$n_subjects != classification_tbl$n_evaluable_by_horizon + classification_tbl$n_unknown_due_to_early_followup_loss)) {
    stop(sprintf("[%s] n_subjects != evaluable + unknown.", table_label), call. = FALSE)
  }
  if (any(classification_tbl$n_event_by_horizon != classification_tbl$tp_count + classification_tbl$fn_count)) {
    stop(sprintf("[%s] event count mismatch versus tp/fn.", table_label), call. = FALSE)
  }
  if (any(classification_tbl$n_known_nonevent_by_horizon != classification_tbl$fp_count + classification_tbl$tn_count)) {
    stop(sprintf("[%s] known nonevent mismatch versus fp/tn.", table_label), call. = FALSE)
  }
  if (any(classification_tbl$n_predicted_positive_evaluable != classification_tbl$tp_count + classification_tbl$fp_count)) {
    stop(sprintf("[%s] evaluable predicted positive mismatch versus tp/fp.", table_label), call. = FALSE)
  }
  if (any(classification_tbl$n_predicted_negative_evaluable != classification_tbl$tn_count + classification_tbl$fn_count)) {
    stop(sprintf("[%s] evaluable predicted negative mismatch versus tn/fn.", table_label), call. = FALSE)
  }

  for (nm in c(
    "positive_classification_rate", "positive_classification_rate_evaluable", "false_positive_rate",
    "false_positive_burden", "specificity", "ppv", "tpr", "event_prevalence_all",
    "event_prevalence_evaluable", "remission_prevalence_all"
  )) {
    assert_prob_in_unit_interval(classification_tbl[[nm]], paste0(table_label, " ", nm))
  }

  invisible(TRUE)
}

validate_delta_table <- function(delta_tbl, table_label) {
  if (nrow(delta_tbl) == 0L) {
    return(invisible(TRUE))
  }

  required_cols <- c(
    "comparison_scope", "comparison_family", "sensitivity_method",
    "risk_change_signed", "risk_change_absolute", "remission_risk_added",
    "false_positive_burden_change", "false_positive_per_100_change",
    "positive_classification_rate_change", "positive_classification_rate_evaluable_change",
    "false_positive_rate_change", "net_benefit_change", "net_benefit_evaluable_change",
    "standardized_net_benefit_change",
    "main_mean_predicted_transition_risk", "sensitivity_mean_predicted_transition_risk",
    "main_mean_predicted_remission_risk", "sensitivity_mean_predicted_remission_risk",
    "main_false_positive_burden", "sensitivity_false_positive_burden",
    "main_false_positive_per_100", "sensitivity_false_positive_per_100",
    "support_tier", "horizon_evidence_class", "claim_restriction_flag", "risk_scale"
  )
  assert_required_columns(delta_tbl, required_cols, paste0("validate_delta_table:", table_label))
  validate_horizon_reporting_fields(delta_tbl, paste0("delta:", table_label))
  if (!all(delta_tbl$risk_scale %in% supplementary_risk_scale)) {
    stop(sprintf("[%s] risk_scale must remain supplementary remission-sensitive scale.", table_label), call. = FALSE)
  }

  risk_abs_diff <- abs(delta_tbl$risk_change_absolute - abs(delta_tbl$risk_change_signed))
  risk_abs_diff[is.na(risk_abs_diff)] <- 0
  if (any(risk_abs_diff > 1e-10)) {
    stop(sprintf("[%s] risk_change_absolute != abs(risk_change_signed).", table_label), call. = FALSE)
  }

  invisible(TRUE)
}


validate_horizon_alignment <- function(tbl, horizon_registry, table_label) {
  if (nrow(tbl) == 0L) {
    return(invisible(TRUE))
  }
  expected <- horizon_registry %>%
    distinct(dataset, horizon_id, horizon_year)
  observed <- tbl %>%
    distinct(dataset, horizon_id, horizon_year)
  missing_pairs <- dplyr::anti_join(expected, observed, by = c("dataset", "horizon_id", "horizon_year"))
  unexpected_pairs <- dplyr::anti_join(observed, expected, by = c("dataset", "horizon_id", "horizon_year"))
  if (nrow(missing_pairs) > 0L) {
    stop(sprintf("[%s] Existing outputs are missing current Stage 1 horizon rows.", table_label), call. = FALSE)
  }
  if (nrow(unexpected_pairs) > 0L) {
    stop(sprintf("[%s] Existing outputs contain horizon rows not present in current Stage 1 registry.", table_label), call. = FALSE)
  }
  invisible(TRUE)
}

validate_threshold_alignment <- function(tbl, threshold_registry, table_label) {
  if (nrow(tbl) == 0L || !all(c("threshold_id", "threshold", "threshold_label") %in% names(tbl))) {
    return(invisible(TRUE))
  }
  expected <- threshold_registry %>%
    distinct(threshold_id, threshold, threshold_label)
  observed <- tbl %>%
    distinct(threshold_id, threshold, threshold_label)
  missing_rows <- dplyr::anti_join(expected, observed, by = c("threshold_id", "threshold", "threshold_label"))
  unexpected_rows <- dplyr::anti_join(observed, expected, by = c("threshold_id", "threshold", "threshold_label"))
  if (nrow(missing_rows) > 0L) {
    stop(sprintf("[%s] Existing outputs are missing current Stage 1 threshold rows.", table_label), call. = FALSE)
  }
  if (nrow(unexpected_rows) > 0L) {
    stop(sprintf("[%s] Existing outputs contain threshold rows not present in current Stage 1 registry.", table_label), call. = FALSE)
  }
  invisible(TRUE)
}

validate_formula_alignment <- function(tbl, formula_registry, table_label) {
  if (nrow(tbl) == 0L || !all(c("dataset", "formula_id") %in% names(tbl))) {
    return(invisible(TRUE))
  }
  expected <- formula_registry %>% distinct(dataset, formula_id)
  observed <- tbl %>% filter(!is.na(formula_id)) %>% distinct(dataset, formula_id)
  missing_rows <- dplyr::anti_join(expected, observed, by = c("dataset", "formula_id"))
  unexpected_rows <- dplyr::anti_join(observed, expected, by = c("dataset", "formula_id"))
  if (nrow(missing_rows) > 0L) {
    stop(sprintf("[%s] Existing outputs are missing formula IDs from the current Stage 1 registry.", table_label), call. = FALSE)
  }
  if (nrow(unexpected_rows) > 0L) {
    stop(sprintf("[%s] Existing outputs contain formula IDs not present in current Stage 1 registry.", table_label), call. = FALSE)
  }
  invisible(TRUE)
}

validate_export_manifest <- function(export_manifest, export_path) {
  if (nrow(export_manifest) == 0L) {
    return(invisible(TRUE))
  }

  missing_files <- export_manifest$file_name[!file.exists(file.path(export_path, export_manifest$file_name))]
  if (length(missing_files) > 0L) {
    stop(
      sprintf("Export manifest references files that were not created: %s", paste(missing_files, collapse = ", ")),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# 🔴 Define: factor helpers, status summaries, and benchmark engines ===============================
## 🟠 Define: panel helpers ===============================
make_dataset_factor <- function(x, dataset_order) {
  factor(as.character(x), levels = dataset_order)
}

make_threshold_factor <- function(x, threshold_registry) {
  factor(as.character(x), levels = unique(as.character(threshold_registry$threshold_label)))
}

make_panel_factor <- function(x, levels = NULL) {
  x_chr <- as.character(x)
  if (is.null(levels)) {
    levels <- unique(x_chr)
  }
  factor(x_chr, levels = levels)
}

make_benchmark_risk_panel_label <- function(dataset, comparison_id, assignment_value) {
  dplyr::case_when(
    as.character(comparison_id) == "overall" ~ paste(as.character(dataset), "Overall group-assigned benchmark", sep = " | "),
    TRUE ~ paste(as.character(dataset), paste("Site-adjusted benchmark:", as.character(assignment_value)), sep = " | ")
  )
}

make_benchmark_delta_panel_label <- function(dataset, comparison_id) {
  dplyr::case_when(
    as.character(comparison_id) == "overall" ~ paste(as.character(dataset), "Overall group-assigned benchmark", sep = " | "),
    TRUE ~ paste(as.character(dataset), "Site-adjusted group-assigned benchmark", sep = " | ")
  )
}

make_formula_panel_label <- function(dataset, formula_label) {
  paste(as.character(dataset), as.character(formula_label), sep = " | ")
}

make_benchmark_risk_panel_levels <- function(state_probabilities_tbl, dataset_order) {
  state_probabilities_tbl %>%
    distinct(dataset, comparison_id, assignment_value) %>%
    mutate(
      dataset = factor(as.character(dataset), levels = dataset_order),
      comparison_id = factor(as.character(comparison_id), levels = c("overall", "site_adjusted")),
      panel_label = make_benchmark_risk_panel_label(dataset, comparison_id, assignment_value)
    ) %>%
    arrange(dataset, comparison_id, assignment_value) %>%
    pull(panel_label) %>%
    unique()
}

make_benchmark_delta_panel_levels <- function(delta_tbl, dataset_order) {
  delta_tbl %>%
    distinct(dataset, comparison_id) %>%
    mutate(
      dataset = factor(as.character(dataset), levels = dataset_order),
      comparison_id = factor(as.character(comparison_id), levels = c("overall", "site_adjusted")),
      panel_label = make_benchmark_delta_panel_label(dataset, comparison_id)
    ) %>%
    arrange(dataset, comparison_id) %>%
    pull(panel_label) %>%
    unique()
}

make_formula_panel_levels <- function(summary_tbl, formula_registry, dataset_order) {
  candidate_levels <- formula_registry %>%
    filter(dataset %in% unique(as.character(summary_tbl$dataset))) %>%
    mutate(dataset = factor(as.character(dataset), levels = dataset_order)) %>%
    arrange(dataset, formula_id) %>%
    transmute(panel_label = make_formula_panel_label(dataset, formula_label)) %>%
    pull(panel_label) %>%
    unique()

  observed_levels <- summary_tbl %>%
    distinct(dataset, formula_label) %>%
    mutate(
      dataset = factor(as.character(dataset), levels = dataset_order),
      panel_label = make_formula_panel_label(dataset, formula_label)
    ) %>%
    arrange(dataset, formula_label) %>%
    pull(panel_label) %>%
    unique()

  candidate_levels[candidate_levels %in% observed_levels]
}

build_projection_rectangles <- function(df, panel_col = "panel_label") {
  if (nrow(df) == 0L || !all(c(panel_col, "horizon_year", "claim_restriction_flag") %in% names(df))) {
    return(tibble())
  }

  df %>%
    mutate(projection_flag = claim_restriction_flag %in% c("projection_only", "projection_plus_prior_sensitive")) %>%
    group_by(across(all_of(panel_col))) %>%
    summarise(
      has_projection = any(projection_flag, na.rm = TRUE),
      xmin = if (any(projection_flag, na.rm = TRUE)) min(horizon_year[projection_flag], na.rm = TRUE) - 0.5 else NA_real_,
      xmax = max(horizon_year, na.rm = TRUE) + 0.5,
      ymin = -Inf,
      ymax = Inf,
      .groups = "drop"
    ) %>%
    filter(has_projection) %>%
    select(-has_projection)
}

## 🟠 Define: status summary ===============================
build_status_summary <- function(analysis_datasets) {
  dataset_rows <- bind_rows(lapply(names(analysis_datasets), function(dataset_name) {
    df <- analysis_datasets[[dataset_name]]

    tibble(
      summary_level = "dataset",
      dataset = dataset_name,
      group_variable = "dataset",
      group_value = dataset_name,
      n_subjects = nrow(df),
      n_status_0_right_censor = sum(df$status_num == 0L),
      n_status_1_transition = sum(df$status_num == 1L),
      n_status_2_remission = sum(df$status_num == 2L),
      prop_status_0_right_censor = mean(df$status_num == 0L),
      prop_status_1_transition = mean(df$status_num == 1L),
      prop_status_2_remission = mean(df$status_num == 2L),
      median_followup_years = stats::median(df$time_year),
      max_followup_years = max(df$time_year),
      note = dplyr::case_when(
        dataset_name == "PNU" ~ "Right censoring and remission separated explicitly.",
        dataset_name == "SNU" ~ "Remission is structurally absent or near-absent if encoded as 2.",
        TRUE ~ "Merged cohort retains original site-specific status patterns."
      )
    )
  }))

  site_rows <- bind_rows(lapply(names(analysis_datasets), function(dataset_name) {
    df <- analysis_datasets[[dataset_name]]
    df %>%
      group_by(site) %>%
      summarise(
        n_subjects = n(),
        n_status_0_right_censor = sum(status_num == 0L),
        n_status_1_transition = sum(status_num == 1L),
        n_status_2_remission = sum(status_num == 2L),
        prop_status_0_right_censor = mean(status_num == 0L),
        prop_status_1_transition = mean(status_num == 1L),
        prop_status_2_remission = mean(status_num == 2L),
        median_followup_years = stats::median(time_year),
        max_followup_years = max(time_year),
        .groups = "drop"
      ) %>%
      mutate(
        summary_level = "dataset_by_site",
        dataset = dataset_name,
        group_variable = "site",
        group_value = as.character(site),
        note = "Site-specific status composition retained for remission sensitivity."
      ) %>%
      select(
        summary_level, dataset, group_variable, group_value,
        n_subjects,
        n_status_0_right_censor, n_status_1_transition, n_status_2_remission,
        prop_status_0_right_censor, prop_status_1_transition, prop_status_2_remission,
        median_followup_years, max_followup_years, note
      )
  }))

  bind_rows(dataset_rows, site_rows) %>%
    arrange(factor(dataset, levels = c("PNU", "SNU", "merged")), summary_level, group_value)
}

## 🟠 Define: benchmark helper metadata ===============================
make_benchmark_meta <- function(dataset_name, df, allow_site_adjusted = TRUE) {
  overall_row <- tibble(
    dataset = dataset_name,
    benchmark_core_id = "overall",
    benchmark_label = "Overall remission-sensitive group-assigned benchmark",
    assignment_variable = "overall",
    main_variant_id = "km_overall",
    main_variant_label = "Overall KM main mirror (group-assigned)",
    sensitivity_variant_id = "aj_overall",
    sensitivity_variant_label = "Overall Aalen-Johansen benchmark (group-assigned)",
    sensitivity_method = "nonparametric_competingrisk",
    group_value = "overall"
  )

  if (dataset_name != "merged" || !isTRUE(allow_site_adjusted)) {
    return(overall_row)
  }

  site_rows <- tibble(
    dataset = dataset_name,
    benchmark_core_id = "site_adjusted",
    benchmark_label = "Site-adjusted remission-sensitive group-assigned benchmark",
    assignment_variable = "site",
    main_variant_id = "km_site_adjusted",
    main_variant_label = "Site-adjusted KM main mirror (group-assigned)",
    sensitivity_variant_id = "aj_site_adjusted",
    sensitivity_variant_label = "Site-adjusted Aalen-Johansen benchmark (group-assigned)",
    sensitivity_method = "nonparametric_competingrisk",
    group_value = sort(unique(as.character(df$site)))
  )

  bind_rows(overall_row, site_rows)
}

## 🟠 Define: nonparametric engines ===============================
compute_km_curve <- function(time, event) {
  time <- as.numeric(time)
  event <- as.integer(event)

  if (length(time) != length(event)) {
    stop("`time` and `event` must have equal length.", call. = FALSE)
  }

  event_times <- sort(unique(time[event == 1L]))

  if (length(event_times) == 0L) {
    return(
      tibble(
        time = numeric(0),
        n_risk = numeric(0),
        d_event = numeric(0),
        km_survival = numeric(0),
        km_risk = numeric(0)
      )
    )
  }

  n_risk <- vapply(event_times, function(t) sum(time >= t), numeric(1))
  d_event <- vapply(event_times, function(t) sum(time == t & event == 1L), numeric(1))
  km_survival <- cumprod(pmax(0, 1 - safe_divide(d_event, n_risk)))

  tibble(
    time = event_times,
    n_risk = n_risk,
    d_event = d_event,
    km_survival = km_survival,
    km_risk = 1 - km_survival
  )
}

evaluate_km_at_horizons <- function(curve_df, horizons) {
  tibble(
    horizon_year = horizons,
    main_km_survival = step_eval(curve_df$time, curve_df$km_survival, horizons, left_value = 1),
    main_km_risk = step_eval(curve_df$time, curve_df$km_risk, horizons, left_value = 0)
  )
}

compute_aj_curve <- function(time, status) {
  time <- as.numeric(time)
  status <- as.integer(status)

  if (length(time) != length(status)) {
    stop("`time` and `status` must have equal length.", call. = FALSE)
  }

  event_times <- sort(unique(time[status %in% c(1L, 2L)]))

  if (length(event_times) == 0L) {
    return(
      tibble(
        time = numeric(0),
        n_risk = numeric(0),
        d_transition = numeric(0),
        d_remission = numeric(0),
        p_entry = numeric(0),
        transition_cif = numeric(0),
        remission_cif = numeric(0)
      )
    )
  }

  n_risk <- vapply(event_times, function(t) sum(time >= t), numeric(1))
  d_transition <- vapply(event_times, function(t) sum(time == t & status == 1L), numeric(1))
  d_remission <- vapply(event_times, function(t) sum(time == t & status == 2L), numeric(1))

  p_entry <- numeric(length(event_times))
  transition_cif <- numeric(length(event_times))
  remission_cif <- numeric(length(event_times))

  s_prev <- 1
  f_transition_prev <- 0
  f_remission_prev <- 0

  for (i in seq_along(event_times)) {
    y_i <- n_risk[i]
    if (!is.finite(y_i) || y_i <= 0) {
      p_entry[i] <- s_prev
      transition_cif[i] <- f_transition_prev
      remission_cif[i] <- f_remission_prev
      next
    }

    h_transition <- d_transition[i] / y_i
    h_remission <- d_remission[i] / y_i

    f_transition_prev <- f_transition_prev + s_prev * h_transition
    f_remission_prev <- f_remission_prev + s_prev * h_remission
    s_prev <- s_prev * (1 - h_transition - h_remission)

    p_entry[i] <- s_prev
    transition_cif[i] <- f_transition_prev
    remission_cif[i] <- f_remission_prev
  }

  tibble(
    time = event_times,
    n_risk = n_risk,
    d_transition = d_transition,
    d_remission = d_remission,
    p_entry = p_entry,
    transition_cif = transition_cif,
    remission_cif = remission_cif
  )
}

evaluate_aj_at_horizons <- function(curve_df, horizons) {
  tibble(
    horizon_year = horizons,
    sensitivity_eventfree_prob = step_eval(curve_df$time, curve_df$p_entry, horizons, left_value = 1),
    sensitivity_transition_risk = step_eval(curve_df$time, curve_df$transition_cif, horizons, left_value = 0),
    sensitivity_remission_risk = step_eval(curve_df$time, curve_df$remission_cif, horizons, left_value = 0)
  )
}

build_group_horizon_table <- function(df_group, dataset_name, meta_row, horizons_tbl) {
  km_curve <- compute_km_curve(df_group$time_year, as.integer(df_group$status_num == 1L))
  aj_curve <- compute_aj_curve(df_group$time_year, df_group$status_num)

  km_horizon <- evaluate_km_at_horizons(km_curve, horizons_tbl$horizon_year)
  aj_horizon <- evaluate_aj_at_horizons(aj_curve, horizons_tbl$horizon_year)

  horizons_tbl %>%
    select(
      horizon_id, horizon_year, horizon_days,
      interpretation_tier, primary_supported_flag, interpretation_note,
      support_tier, support_tier_standard, horizon_evidence_class, claim_restriction_flag,
      horizon_support_label, supported_for_reporting, projection_flag, reporting_status
    ) %>%
    left_join(km_horizon, by = "horizon_year") %>%
    left_join(aj_horizon, by = "horizon_year") %>%
    mutate(
      dataset = dataset_name,
      comparison_scope = "benchmark",
      comparison_family = "nonparametric",
      comparison_id = meta_row$benchmark_core_id[[1]],
      prediction_level = "group_assigned",
      model_class = "Stage9_benchmark",
      model_family = "nonparametric",
      risk_scale = supplementary_risk_scale,
      main_model_family = "KaplanMeier",
      sensitivity_model_family = "AalenJohansen",
      formula_id = NA_character_,
      formula_name = NA_character_,
      formula_label = NA_character_,
      site_branch = ifelse(meta_row$assignment_variable[[1]] == "site", "site_adjusted", "site_free"),
      interaction_branch = NA_character_,
      benchmark_core_id = meta_row$benchmark_core_id[[1]],
      benchmark_label = meta_row$benchmark_label[[1]],
      assignment_variable = meta_row$assignment_variable[[1]],
      assignment_value = meta_row$group_value[[1]],
      main_variant_id = meta_row$main_variant_id[[1]],
      main_variant_label = meta_row$main_variant_label[[1]],
      sensitivity_variant_id = meta_row$sensitivity_variant_id[[1]],
      sensitivity_variant_label = meta_row$sensitivity_variant_label[[1]],
      sensitivity_method = meta_row$sensitivity_method[[1]],
      n_subjects_in_group = nrow(df_group),
      n_transition_total = sum(df_group$status_num == 1L),
      n_remission_total = sum(df_group$status_num == 2L),
      n_right_censor_total = sum(df_group$status_num == 0L),
      n_observed_through_horizon = vapply(horizon_year, function(h) sum(df_group$time_year >= h), numeric(1)),
      n_transition_by_horizon = vapply(horizon_year, function(h) sum(df_group$status_num == 1L & df_group$time_year <= h), numeric(1)),
      n_remission_by_horizon = vapply(horizon_year, function(h) sum(df_group$status_num == 2L & df_group$time_year <= h), numeric(1)),
      n_right_censor_by_horizon = vapply(horizon_year, function(h) sum(df_group$status_num == 0L & df_group$time_year <= h), numeric(1)),
      risk_difference_signed = sensitivity_transition_risk - main_km_risk,
      risk_difference_absolute = abs(risk_difference_signed)
    ) %>%
    arrange(horizon_year)
}

build_benchmark_outputs <- function(df, dataset_name, horizons_tbl, allow_site_adjusted = TRUE) {
  meta_tbl <- make_benchmark_meta(dataset_name = dataset_name, df = df, allow_site_adjusted = allow_site_adjusted)

  comparison_ids <- unique(meta_tbl$benchmark_core_id)
  comparison_group_tables <- vector("list", length(comparison_ids))
  comparison_subject_tables <- vector("list", length(comparison_ids))
  fit_object_list <- list()

  for (i in seq_along(comparison_ids)) {
    comparison_id <- comparison_ids[i]
    meta_rows <- meta_tbl %>% filter(benchmark_core_id == comparison_id)

    group_tables_for_comparison <- vector("list", nrow(meta_rows))
    fit_objects_for_comparison <- vector("list", nrow(meta_rows))

    for (j in seq_len(nrow(meta_rows))) {
      meta_row <- meta_rows[j, , drop = FALSE]

      if (identical(meta_row$assignment_variable[[1]], "overall")) {
        df_group <- df
      } else {
        group_value <- as.character(meta_row$group_value[[1]])
        df_group <- df[df[[meta_row$assignment_variable[[1]]]] == group_value, , drop = FALSE]
      }

      if (nrow(df_group) == 0L) {
        next
      }

      group_tables_for_comparison[[j]] <- build_group_horizon_table(
        df_group = df_group,
        dataset_name = dataset_name,
        meta_row = meta_row,
        horizons_tbl = horizons_tbl
      )

      fit_objects_for_comparison[[j]] <- list(
        dataset = dataset_name,
        comparison_id = comparison_id,
        benchmark_label = meta_row$benchmark_label[[1]],
        assignment_variable = meta_row$assignment_variable[[1]],
        group_value = meta_row$group_value[[1]],
        km_curve = compute_km_curve(df_group$time_year, as.integer(df_group$status_num == 1L)),
        aj_curve = compute_aj_curve(df_group$time_year, df_group$status_num)
      )
    }

    group_horizon_tbl <- compact_bind_rows(group_tables_for_comparison)

    subject_base <- df %>%
      mutate(
        assignment_value = if (identical(meta_rows$assignment_variable[[1]], "overall")) {
          "overall"
        } else {
          as.character(.data[[meta_rows$assignment_variable[[1]]]])
        }
      ) %>%
      select(
        unique_person_id, id, site, sex_num, sex_label, age_exact_entry, age_s,
        days_followup, time_year, status_num, status_label, assignment_value
      )

    subject_joined <- subject_base %>%
      left_join(
        group_horizon_tbl %>%
          select(
            dataset, comparison_scope, comparison_family, comparison_id, prediction_level,
            model_class, risk_scale, main_model_family, sensitivity_model_family,
            formula_id, formula_name, formula_label,
            site_branch, interaction_branch, benchmark_core_id, benchmark_label,
            assignment_variable, assignment_value,
            main_variant_id, main_variant_label,
            sensitivity_variant_id, sensitivity_variant_label,
            sensitivity_method,
            horizon_id, horizon_year, horizon_days,
            interpretation_tier, primary_supported_flag,
            interpretation_note, support_tier, support_tier_standard,
            horizon_evidence_class, claim_restriction_flag,
            horizon_support_label, supported_for_reporting, projection_flag, reporting_status,
            main_km_survival, main_km_risk,
            sensitivity_eventfree_prob,
            sensitivity_transition_risk,
            sensitivity_remission_risk,
            risk_difference_signed,
            risk_difference_absolute
          ),
        by = c("assignment_value")
      )

    prediction_main <- subject_joined %>%
      transmute(
        dataset,
        comparison_scope,
        comparison_family,
        comparison_id,
        prediction_level,
        model_class,
        model_family = main_model_family,
        risk_scale,
        formula_id,
        formula_name,
        formula_label,
        site_branch,
        interaction_branch,
        benchmark_core_id,
        benchmark_label,
        assignment_variable,
        assignment_value,
        variant_id = main_variant_id,
        variant_label = main_variant_label,
        sensitivity_method,
        analysis_variant = "main_censoring_mirror",
        evaluation_mode = "main_censoring",
        risk_estimator = "KaplanMeier",
        unique_person_id,
        id,
        site,
        sex_num,
        sex_label,
        age_exact_entry,
        age_s,
        days_followup,
        time_year,
        status_num,
        status_label,
        horizon_id,
        horizon_year,
        horizon_days,
        interpretation_tier,
        primary_supported_flag,
        interpretation_note,
        support_tier,
        support_tier_standard,
        horizon_evidence_class,
        claim_restriction_flag,
        horizon_support_label,
        supported_for_reporting,
        projection_flag,
        reporting_status,
        predicted_transition_risk = main_km_risk,
        predicted_remission_risk = 0,
        predicted_eventfree_prob = main_km_survival,
        risk_difference_from_main_signed = 0,
        risk_difference_from_main_absolute = 0
      ) %>%
      bind_component_annotations()

    prediction_sensitive <- subject_joined %>%
      transmute(
        dataset,
        comparison_scope,
        comparison_family,
        comparison_id,
        prediction_level,
        model_class,
        model_family = sensitivity_model_family,
        risk_scale,
        formula_id,
        formula_name,
        formula_label,
        site_branch,
        interaction_branch,
        benchmark_core_id,
        benchmark_label,
        assignment_variable,
        assignment_value,
        variant_id = sensitivity_variant_id,
        variant_label = sensitivity_variant_label,
        sensitivity_method,
        analysis_variant = "remission_sensitive",
        evaluation_mode = "remission_sensitive",
        risk_estimator = "AalenJohansen",
        unique_person_id,
        id,
        site,
        sex_num,
        sex_label,
        age_exact_entry,
        age_s,
        days_followup,
        time_year,
        status_num,
        status_label,
        horizon_id,
        horizon_year,
        horizon_days,
        interpretation_tier,
        primary_supported_flag,
        interpretation_note,
        support_tier,
        support_tier_standard,
        horizon_evidence_class,
        claim_restriction_flag,
        horizon_support_label,
        supported_for_reporting,
        projection_flag,
        reporting_status,
        predicted_transition_risk = sensitivity_transition_risk,
        predicted_remission_risk = sensitivity_remission_risk,
        predicted_eventfree_prob = sensitivity_eventfree_prob,
        risk_difference_from_main_signed = risk_difference_signed,
        risk_difference_from_main_absolute = risk_difference_absolute
      ) %>%
      bind_component_annotations()

    comparison_group_tables[[i]] <- group_horizon_tbl
    comparison_subject_tables[[i]] <- bind_rows(prediction_main, prediction_sensitive)
    fit_object_list[[comparison_id]] <- fit_objects_for_comparison
  }

  output <- list(
    group_horizon_table = compact_bind_rows(comparison_group_tables),
    subject_predictions = compact_bind_rows(comparison_subject_tables),
    fit_objects = fit_object_list
  )

  validate_benchmark_state_probabilities(output$group_horizon_table)
  validate_subject_predictions(output$subject_predictions, "benchmark")
  output
}


# 🔴 Define: threshold evaluation and delta builders ===============================
## 🟠 Define: threshold-based classification ===============================
make_threshold_classification <- function(predictions_long, threshold_registry) {
  if (nrow(predictions_long) == 0L) {
    return(tibble())
  }

  group_cols <- c(
    "dataset", "comparison_scope", "comparison_family", "comparison_id",
    "prediction_level",
    "model_class", "model_family", "risk_scale",
    "formula_id", "formula_name", "formula_label",
    "site_branch", "interaction_branch",
    "benchmark_core_id", "benchmark_label",
    "assignment_variable",
    "variant_id", "variant_label",
    "sensitivity_method",
    "analysis_variant", "evaluation_mode", "risk_estimator",
    "horizon_id", "horizon_year", "horizon_days",
    "interpretation_tier", "primary_supported_flag",
    "interpretation_note", "support_tier", "support_tier_standard",
    "horizon_evidence_class", "claim_restriction_flag",
    "horizon_support_label", "supported_for_reporting", "projection_flag", "reporting_status",
    "remission_component_available", "eventfree_component_available", "component_mode", "component_note"
  )

  out <- predictions_long %>%
    tidyr::crossing(
      threshold_registry %>%
        select(threshold_id, threshold, threshold_label, positive_rule)
    ) %>%
    mutate(
      predicted_positive = predicted_transition_risk >= threshold,
      event_by_horizon = status_num == 1L & time_year <= horizon_year,
      remission_by_horizon = status_num == 2L & time_year <= horizon_year,
      known_nonevent_by_horizon = dplyr::case_when(
        evaluation_mode == "main_censoring" ~ (!event_by_horizon & time_year >= horizon_year),
        evaluation_mode == "remission_sensitive" ~ (remission_by_horizon | (!event_by_horizon & time_year >= horizon_year)),
        TRUE ~ FALSE
      ),
      evaluable_by_horizon = event_by_horizon | known_nonevent_by_horizon,
      unknown_due_to_early_followup_loss = !evaluable_by_horizon
    ) %>%
    group_by(across(all_of(group_cols)), threshold_id, threshold, threshold_label, positive_rule) %>%
    summarise(
      n_subjects = n(),
      n_event_by_horizon = sum(event_by_horizon),
      n_remission_by_horizon = sum(remission_by_horizon),
      n_known_nonevent_by_horizon = sum(known_nonevent_by_horizon),
      n_evaluable_by_horizon = sum(evaluable_by_horizon),
      n_unknown_due_to_early_followup_loss = sum(unknown_due_to_early_followup_loss),
      n_predicted_positive = sum(predicted_positive),
      n_predicted_negative = sum(!predicted_positive),
      n_predicted_positive_evaluable = sum(predicted_positive & evaluable_by_horizon),
      n_predicted_negative_evaluable = sum((!predicted_positive) & evaluable_by_horizon),
      tp_count = sum(predicted_positive & event_by_horizon),
      fp_count = sum(predicted_positive & known_nonevent_by_horizon),
      tn_count = sum(!predicted_positive & known_nonevent_by_horizon),
      fn_count = sum(!predicted_positive & event_by_horizon),
      mean_predicted_transition_risk = mean(predicted_transition_risk),
      min_predicted_transition_risk = min(predicted_transition_risk),
      max_predicted_transition_risk = max(predicted_transition_risk),
      mean_predicted_remission_risk = if (all(is.na(predicted_remission_risk))) NA_real_ else mean(predicted_remission_risk, na.rm = TRUE),
      mean_predicted_eventfree_prob = if (all(is.na(predicted_eventfree_prob))) NA_real_ else mean(predicted_eventfree_prob, na.rm = TRUE),
      mean_risk_difference_from_main_signed = mean(risk_difference_from_main_signed),
      mean_risk_difference_from_main_absolute = mean(risk_difference_from_main_absolute),
      .groups = "drop"
    ) %>%
    mutate(
      positive_classification_rate = safe_divide(n_predicted_positive, n_subjects),
      positive_classification_rate_evaluable = safe_divide(n_predicted_positive_evaluable, n_evaluable_by_horizon),
      false_positive_rate = safe_divide(fp_count, n_known_nonevent_by_horizon),
      false_positive_burden = safe_divide(fp_count, n_subjects),
      false_positive_burden_known_nonevent = safe_divide(fp_count, n_known_nonevent_by_horizon),
      false_positive_per_100 = 100 * false_positive_burden,
      specificity = safe_divide(tn_count, tn_count + fp_count),
      ppv = safe_divide(tp_count, tp_count + fp_count),
      tpr = safe_divide(tp_count, tp_count + fn_count),
      event_prevalence_all = safe_divide(n_event_by_horizon, n_subjects),
      event_prevalence_evaluable = safe_divide(n_event_by_horizon, n_evaluable_by_horizon),
      remission_prevalence_all = safe_divide(n_remission_by_horizon, n_subjects),
      net_benefit = safe_divide(tp_count, n_subjects) - safe_divide(fp_count, n_subjects) * (threshold / (1 - threshold)),
      net_benefit_evaluable = safe_divide(tp_count, n_evaluable_by_horizon) - safe_divide(fp_count, n_evaluable_by_horizon) * (threshold / (1 - threshold)),
      standardized_net_benefit = safe_divide(net_benefit, event_prevalence_all),
      evaluation_rule = dplyr::case_when(
        evaluation_mode == "main_censoring" ~ "event by horizon OR observed event-free through horizon; remission before horizon remains censored/unknown",
        evaluation_mode == "remission_sensitive" ~ "event by horizon OR remission by horizon OR observed transition-free through horizon; early right censoring before horizon remains unknown",
        TRUE ~ "unknown"
      )
    ) %>%
    arrange(
      factor(dataset, levels = c("PNU", "SNU", "merged")),
      comparison_scope, comparison_id,
      horizon_year, threshold,
      analysis_variant
    )

  validate_classification_table(out, "threshold_classification")
  out
}

## 🟠 Define: risk summaries ===============================
make_risk_summary_from_predictions <- function(predictions_long) {
  if (nrow(predictions_long) == 0L) {
    return(tibble())
  }

  key_cols <- c(
    "dataset", "comparison_scope", "comparison_family", "comparison_id",
    "prediction_level",
    "model_class",
    "formula_id", "formula_name", "formula_label",
    "site_branch", "interaction_branch",
    "benchmark_core_id", "benchmark_label",
    "assignment_variable",
    "sensitivity_method",
    "risk_scale",
    "horizon_id", "horizon_year", "horizon_days",
    "interpretation_tier", "primary_supported_flag",
    "interpretation_note", "support_tier", "support_tier_standard",
    "horizon_evidence_class", "claim_restriction_flag",
    "horizon_support_label", "supported_for_reporting", "projection_flag", "reporting_status",
    "unique_person_id"
  )

  main_tbl <- predictions_long %>%
    filter(analysis_variant == "main_censoring_mirror") %>%
    select(
      all_of(key_cols),
      main_model_family = model_family,
      main_variant_id = variant_id,
      main_variant_label = variant_label,
      main_risk_estimator = risk_estimator,
      main_transition_risk = predicted_transition_risk,
      main_remission_component_available = remission_component_available,
      main_eventfree_component_available = eventfree_component_available,
      main_component_mode = component_mode,
      main_component_note = component_note
    )

  sensitivity_tbl <- predictions_long %>%
    filter(analysis_variant == "remission_sensitive") %>%
    select(
      all_of(key_cols),
      sensitivity_model_family = model_family,
      sensitivity_variant_id = variant_id,
      sensitivity_variant_label = variant_label,
      sensitivity_risk_estimator = risk_estimator,
      sensitivity_transition_risk = predicted_transition_risk,
      sensitivity_remission_risk = predicted_remission_risk,
      sensitivity_eventfree_prob = predicted_eventfree_prob,
      risk_difference_from_main_signed,
      risk_difference_from_main_absolute,
      sensitivity_remission_component_available = remission_component_available,
      sensitivity_eventfree_component_available = eventfree_component_available,
      sensitivity_component_mode = component_mode,
      sensitivity_component_note = component_note
    )

  joined <- main_tbl %>%
    inner_join(sensitivity_tbl, by = key_cols)

  if (nrow(joined) != nrow(main_tbl) || nrow(joined) != nrow(sensitivity_tbl)) {
    stop("Main/sensitivity prediction join mismatch in `make_risk_summary_from_predictions()`.", call. = FALSE)
  }

  summary_cols <- setdiff(key_cols, "unique_person_id")

  out <- joined %>%
    group_by(across(all_of(summary_cols))) %>%
    summarise(
      main_model_family = collapse_unique_chr(main_model_family),
      main_variant_id = collapse_unique_chr(main_variant_id),
      main_variant_label = collapse_unique_chr(main_variant_label),
      main_risk_estimator = collapse_unique_chr(main_risk_estimator),
      sensitivity_model_family = collapse_unique_chr(sensitivity_model_family),
      sensitivity_variant_id = collapse_unique_chr(sensitivity_variant_id),
      sensitivity_variant_label = collapse_unique_chr(sensitivity_variant_label),
      sensitivity_risk_estimator = collapse_unique_chr(sensitivity_risk_estimator),
      n_subjects = n(),
      mean_main_transition_risk = mean(main_transition_risk),
      mean_sensitivity_transition_risk = mean(sensitivity_transition_risk),
      mean_sensitivity_remission_risk = if (all(is.na(sensitivity_remission_risk))) NA_real_ else mean(sensitivity_remission_risk, na.rm = TRUE),
      mean_sensitivity_eventfree_prob = if (all(is.na(sensitivity_eventfree_prob))) NA_real_ else mean(sensitivity_eventfree_prob, na.rm = TRUE),
      mean_risk_difference_signed = mean(risk_difference_from_main_signed),
      mean_risk_difference_absolute = mean(risk_difference_from_main_absolute),
      median_risk_difference_absolute = stats::median(risk_difference_from_main_absolute),
      max_risk_difference_absolute = max(risk_difference_from_main_absolute),
      main_remission_component_available = summarise_single_logical(main_remission_component_available, "main_remission_component_available"),
      main_eventfree_component_available = summarise_single_logical(main_eventfree_component_available, "main_eventfree_component_available"),
      main_component_mode = summarise_single_text(main_component_mode, "main_component_mode"),
      main_component_note = summarise_single_text(main_component_note, "main_component_note"),
      sensitivity_remission_component_available = summarise_single_logical(sensitivity_remission_component_available, "sensitivity_remission_component_available"),
      sensitivity_eventfree_component_available = summarise_single_logical(sensitivity_eventfree_component_available, "sensitivity_eventfree_component_available"),
      sensitivity_component_mode = summarise_single_text(sensitivity_component_mode, "sensitivity_component_mode"),
      sensitivity_component_note = summarise_single_text(sensitivity_component_note, "sensitivity_component_note"),
      .groups = "drop"
    ) %>%
    mutate(
      model_family = comparison_family
    ) %>%
    relocate(model_family, .after = model_class) %>%
    arrange(
      factor(dataset, levels = c("PNU", "SNU", "merged")),
      comparison_scope, comparison_id,
      horizon_year
    )

  validate_risk_summary_table(out, "risk_summary")
  out
}

## 🟠 Define: delta tables ===============================
make_delta_table <- function(classification_tbl) {
  if (nrow(classification_tbl) == 0L) {
    return(tibble())
  }

  join_keys <- c(
    "dataset", "comparison_scope", "comparison_family", "comparison_id",
    "prediction_level",
    "model_class", "risk_scale",
    "formula_id", "formula_name", "formula_label",
    "site_branch", "interaction_branch",
    "benchmark_core_id", "benchmark_label",
    "assignment_variable",
    "sensitivity_method",
    "horizon_id", "horizon_year", "horizon_days",
    "interpretation_tier", "primary_supported_flag",
    "interpretation_note", "support_tier", "support_tier_standard",
    "horizon_evidence_class", "claim_restriction_flag",
    "horizon_support_label",
    "supported_for_reporting", "projection_flag", "reporting_status",
    "threshold_id", "threshold", "threshold_label", "positive_rule"
  )

  main_tbl <- classification_tbl %>%
    filter(analysis_variant == "main_censoring_mirror") %>%
    select(
      all_of(join_keys),
      main_model_family = model_family,
      main_variant_id = variant_id,
      main_variant_label = variant_label,
      main_risk_estimator = risk_estimator,
      main_evaluation_mode = evaluation_mode,
      main_mean_predicted_transition_risk = mean_predicted_transition_risk,
      main_mean_predicted_remission_risk = mean_predicted_remission_risk,
      main_mean_predicted_eventfree_prob = mean_predicted_eventfree_prob,
      main_n_subjects = n_subjects,
      main_n_event_by_horizon = n_event_by_horizon,
      main_n_remission_by_horizon = n_remission_by_horizon,
      main_n_known_nonevent_by_horizon = n_known_nonevent_by_horizon,
      main_n_evaluable_by_horizon = n_evaluable_by_horizon,
      main_n_unknown_due_to_early_followup_loss = n_unknown_due_to_early_followup_loss,
      main_n_predicted_positive = n_predicted_positive,
      main_n_predicted_negative = n_predicted_negative,
      main_n_predicted_positive_evaluable = n_predicted_positive_evaluable,
      main_n_predicted_negative_evaluable = n_predicted_negative_evaluable,
      main_tp_count = tp_count,
      main_fp_count = fp_count,
      main_tn_count = tn_count,
      main_fn_count = fn_count,
      main_positive_classification_rate = positive_classification_rate,
      main_positive_classification_rate_evaluable = positive_classification_rate_evaluable,
      main_false_positive_rate = false_positive_rate,
      main_false_positive_burden = false_positive_burden,
      main_false_positive_per_100 = false_positive_per_100,
      main_specificity = specificity,
      main_ppv = ppv,
      main_tpr = tpr,
      main_net_benefit = net_benefit,
      main_net_benefit_evaluable = net_benefit_evaluable,
      main_standardized_net_benefit = standardized_net_benefit,
      main_remission_component_available = remission_component_available,
      main_eventfree_component_available = eventfree_component_available,
      main_component_mode = component_mode,
      main_component_note = component_note
    )

  sensitivity_tbl <- classification_tbl %>%
    filter(analysis_variant == "remission_sensitive") %>%
    select(
      all_of(join_keys),
      sensitivity_model_family = model_family,
      sensitivity_variant_id = variant_id,
      sensitivity_variant_label = variant_label,
      sensitivity_risk_estimator = risk_estimator,
      sensitivity_evaluation_mode = evaluation_mode,
      sensitivity_mean_predicted_transition_risk = mean_predicted_transition_risk,
      sensitivity_mean_predicted_remission_risk = mean_predicted_remission_risk,
      sensitivity_mean_predicted_eventfree_prob = mean_predicted_eventfree_prob,
      sensitivity_n_subjects = n_subjects,
      sensitivity_n_event_by_horizon = n_event_by_horizon,
      sensitivity_n_remission_by_horizon = n_remission_by_horizon,
      sensitivity_n_known_nonevent_by_horizon = n_known_nonevent_by_horizon,
      sensitivity_n_evaluable_by_horizon = n_evaluable_by_horizon,
      sensitivity_n_unknown_due_to_early_followup_loss = n_unknown_due_to_early_followup_loss,
      sensitivity_n_predicted_positive = n_predicted_positive,
      sensitivity_n_predicted_negative = n_predicted_negative,
      sensitivity_n_predicted_positive_evaluable = n_predicted_positive_evaluable,
      sensitivity_n_predicted_negative_evaluable = n_predicted_negative_evaluable,
      sensitivity_tp_count = tp_count,
      sensitivity_fp_count = fp_count,
      sensitivity_tn_count = tn_count,
      sensitivity_fn_count = fn_count,
      sensitivity_positive_classification_rate = positive_classification_rate,
      sensitivity_positive_classification_rate_evaluable = positive_classification_rate_evaluable,
      sensitivity_false_positive_rate = false_positive_rate,
      sensitivity_false_positive_burden = false_positive_burden,
      sensitivity_false_positive_per_100 = false_positive_per_100,
      sensitivity_specificity = specificity,
      sensitivity_ppv = ppv,
      sensitivity_tpr = tpr,
      sensitivity_net_benefit = net_benefit,
      sensitivity_net_benefit_evaluable = net_benefit_evaluable,
      sensitivity_standardized_net_benefit = standardized_net_benefit,
      sensitivity_remission_component_available = remission_component_available,
      sensitivity_eventfree_component_available = eventfree_component_available,
      sensitivity_component_mode = component_mode,
      sensitivity_component_note = component_note
    )

  out <- sensitivity_tbl %>%
    inner_join(main_tbl, by = join_keys) %>%
    mutate(
      risk_change_signed = sensitivity_mean_predicted_transition_risk - main_mean_predicted_transition_risk,
      risk_change_absolute = abs(risk_change_signed),
      remission_risk_added = sensitivity_mean_predicted_remission_risk - main_mean_predicted_remission_risk,
      false_positive_burden_change = sensitivity_false_positive_burden - main_false_positive_burden,
      false_positive_per_100_change = sensitivity_false_positive_per_100 - main_false_positive_per_100,
      positive_classification_rate_change = sensitivity_positive_classification_rate - main_positive_classification_rate,
      positive_classification_rate_evaluable_change = sensitivity_positive_classification_rate_evaluable - main_positive_classification_rate_evaluable,
      false_positive_rate_change = sensitivity_false_positive_rate - main_false_positive_rate,
      net_benefit_change = sensitivity_net_benefit - main_net_benefit,
      net_benefit_evaluable_change = sensitivity_net_benefit_evaluable - main_net_benefit_evaluable,
      standardized_net_benefit_change = sensitivity_standardized_net_benefit - main_standardized_net_benefit,
      fp_count_change = sensitivity_fp_count - main_fp_count,
      tp_count_change = sensitivity_tp_count - main_tp_count,
      n_evaluable_change = sensitivity_n_evaluable_by_horizon - main_n_evaluable_by_horizon,
      n_unknown_change = sensitivity_n_unknown_due_to_early_followup_loss - main_n_unknown_due_to_early_followup_loss,
      n_predicted_positive_evaluable_change = sensitivity_n_predicted_positive_evaluable - main_n_predicted_positive_evaluable,
      model_family = comparison_family
    ) %>%
    relocate(model_family, .after = model_class) %>%
    arrange(
      factor(dataset, levels = c("PNU", "SNU", "merged")),
      comparison_scope, comparison_id,
      horizon_year, threshold
    )

  if (nrow(out) != nrow(main_tbl) || nrow(out) != nrow(sensitivity_tbl)) {
    stop("Main/sensitivity classification join mismatch in `make_delta_table()`.", call. = FALSE)
  }

  validate_delta_table(out, "delta_table")
  out
}


# 🔴 Define: regression engines and prediction assemblers ===============================
## 🟠 Define: formula utilities ===============================
build_formula_object <- function(formula_rhs, event_expression) {
  stats::as.formula(paste0("survival::Surv(time_year, ", event_expression, ") ~ ", formula_rhs))
}

build_rhs_formula <- function(formula_rhs) {
  stats::as.formula(paste0("~", formula_rhs))
}

make_design_matrix <- function(df, formula_rhs, training_columns = NULL) {
  mm <- stats::model.matrix(build_rhs_formula(formula_rhs), data = df)
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }

  if (!is.null(training_columns)) {
    if (length(training_columns) == 0L) {
      return(matrix(numeric(0), nrow = nrow(df), ncol = 0L))
    }

    missing_cols <- setdiff(training_columns, colnames(mm))
    extra_cols <- setdiff(colnames(mm), training_columns)

    if (length(missing_cols) > 0L) {
      zero_block <- matrix(
        0,
        nrow = nrow(mm),
        ncol = length(missing_cols),
        dimnames = list(NULL, missing_cols)
      )
      mm <- cbind(mm, zero_block)
    }

    if (length(extra_cols) > 0L) {
      mm <- mm[, colnames(mm) %in% training_columns, drop = FALSE]
    }

    mm <- mm[, training_columns, drop = FALSE]
  }

  if (is.null(training_columns) && ncol(mm) > 0L) {
    zero_variance <- vapply(
      seq_len(ncol(mm)),
      function(j) {
        vj <- stats::var(mm[, j])
        is.na(vj) || vj == 0
      },
      logical(1)
    )
    if (any(zero_variance)) {
      mm <- mm[, !zero_variance, drop = FALSE]
    }
  }

  mm
}

## 🟠 Define: Cox fitters ===============================
fit_main_transition_cox <- function(df, formula_row) {
  transition_events <- sum(df$status_num == 1L)
  if (transition_events < minimum_transition_events) {
    return(list(ok = FALSE, reason = "too_few_transition_events", fit = NULL))
  }

  fit <- tryCatch(
    survival::coxph(
      formula = build_formula_object(formula_row$formula_rhs[[1]], "status_num == 1L"),
      data = df,
      ties = cox_ties,
      x = TRUE,
      model = TRUE,
      singular.ok = TRUE
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(ok = FALSE, reason = conditionMessage(fit), fit = NULL))
  }

  list(ok = TRUE, reason = NA_character_, fit = fit)
}

fit_remission_cox <- function(df, formula_row) {
  remission_events <- sum(df$status_num == 2L)
  if (remission_events < minimum_remission_events) {
    return(list(ok = FALSE, reason = "too_few_remission_events", fit = NULL))
  }

  fit <- tryCatch(
    survival::coxph(
      formula = build_formula_object(formula_row$formula_rhs[[1]], "status_num == 2L"),
      data = df,
      ties = cox_ties,
      x = TRUE,
      model = TRUE,
      singular.ok = TRUE
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(ok = FALSE, reason = conditionMessage(fit), fit = NULL))
  }

  list(ok = TRUE, reason = NA_character_, fit = fit)
}

baseline_hazard_table <- function(cox_fit) {
  if (is.null(cox_fit)) {
    return(tibble(time = numeric(0), cumhaz = numeric(0), jump = numeric(0)))
  }

  bh <- survival::basehaz(cox_fit, centered = FALSE)
  bh <- tibble(
    time = as.numeric(bh$time),
    cumhaz = as.numeric(bh$hazard)
  ) %>%
    group_by(time) %>%
    summarise(cumhaz = max(cumhaz), .groups = "drop") %>%
    arrange(time) %>%
    mutate(jump = cumhaz - lag(cumhaz, default = 0))

  bh
}

predict_main_cox_risk <- function(main_fit, newdata, horizons) {
  bh1 <- baseline_hazard_table(main_fit)
  lp1 <- as.numeric(stats::predict(main_fit, newdata = newdata, type = "lp"))
  exp_lp1 <- exp(lp1)
  cumhaz_at_horizons <- step_eval(bh1$time, bh1$cumhaz, horizons, left_value = 0)

  risk_matrix <- outer(exp_lp1, cumhaz_at_horizons, FUN = "*")
  risk_matrix <- 1 - exp(-risk_matrix)

  list(
    risk_matrix = risk_matrix,
    lp = lp1,
    baseline = bh1
  )
}

predict_cause_specific_cif <- function(main_fit, remission_fit, newdata, horizons) {
  bh1 <- baseline_hazard_table(main_fit)
  bh2 <- baseline_hazard_table(remission_fit)

  lp1 <- as.numeric(stats::predict(main_fit, newdata = newdata, type = "lp"))
  exp_lp1 <- exp(lp1)

  if (is.null(remission_fit)) {
    lp2 <- rep(NA_real_, nrow(newdata))
    exp_lp2 <- rep(0, nrow(newdata))
  } else {
    lp2 <- as.numeric(stats::predict(remission_fit, newdata = newdata, type = "lp"))
    exp_lp2 <- exp(lp2)
  }

  main_cumhaz_at_horizons <- step_eval(bh1$time, bh1$cumhaz, horizons, left_value = 0)
  main_risk_matrix <- 1 - exp(-outer(exp_lp1, main_cumhaz_at_horizons, FUN = "*"))

  event_times <- sort(unique(c(bh1$time, bh2$time)))
  n_subjects <- nrow(newdata)
  n_horizons <- length(horizons)

  transition_matrix <- matrix(0, nrow = n_subjects, ncol = n_horizons)
  remission_matrix <- matrix(0, nrow = n_subjects, ncol = n_horizons)
  eventfree_matrix <- matrix(1, nrow = n_subjects, ncol = n_horizons)

  if (length(event_times) == 0L) {
    return(list(
      main_risk_matrix = main_risk_matrix,
      transition_matrix = transition_matrix,
      remission_matrix = remission_matrix,
      eventfree_matrix = eventfree_matrix,
      baseline_transition = bh1,
      baseline_remission = bh2,
      lp_transition = lp1,
      lp_remission = lp2
    ))
  }

  jump_table <- tibble(time = event_times) %>%
    left_join(bh1 %>% select(time, jump_transition = jump), by = "time") %>%
    left_join(bh2 %>% select(time, jump_remission = jump), by = "time") %>%
    mutate(
      jump_transition = dplyr::coalesce(jump_transition, 0),
      jump_remission = dplyr::coalesce(jump_remission, 0)
    )

  for (i in seq_len(n_subjects)) {
    s_prev <- 1
    f_transition_prev <- 0
    f_remission_prev <- 0
    h_index <- 1L

    for (j in seq_len(nrow(jump_table))) {
      while (h_index <= n_horizons && horizons[h_index] < jump_table$time[j]) {
        transition_matrix[i, h_index] <- f_transition_prev
        remission_matrix[i, h_index] <- f_remission_prev
        eventfree_matrix[i, h_index] <- s_prev
        h_index <- h_index + 1L
      }

      dH1 <- jump_table$jump_transition[j] * exp_lp1[i]
      dH2 <- jump_table$jump_remission[j] * exp_lp2[i]
      dH_total <- dH1 + dH2

      if (is.finite(dH_total) && dH_total > 0) {
        event_mass <- s_prev * (1 - exp(-dH_total))
        f_transition_prev <- f_transition_prev + event_mass * (dH1 / dH_total)
        f_remission_prev <- f_remission_prev + event_mass * (dH2 / dH_total)
        s_prev <- s_prev * exp(-dH_total)
      }

      while (h_index <= n_horizons && horizons[h_index] <= jump_table$time[j]) {
        transition_matrix[i, h_index] <- f_transition_prev
        remission_matrix[i, h_index] <- f_remission_prev
        eventfree_matrix[i, h_index] <- s_prev
        h_index <- h_index + 1L
      }
    }

    while (h_index <= n_horizons) {
      transition_matrix[i, h_index] <- f_transition_prev
      remission_matrix[i, h_index] <- f_remission_prev
      eventfree_matrix[i, h_index] <- s_prev
      h_index <- h_index + 1L
    }
  }

  list(
    main_risk_matrix = main_risk_matrix,
    transition_matrix = transition_matrix,
    remission_matrix = remission_matrix,
    eventfree_matrix = eventfree_matrix,
    baseline_transition = bh1,
    baseline_remission = bh2,
    lp_transition = lp1,
    lp_remission = lp2
  )
}

## 🟠 Define: Fine-Gray fitters ===============================
fit_fine_gray_transition <- function(df, formula_row) {
  if (!isTRUE(has_cmprsk) || !isTRUE(run_subdistribution_models)) {
    return(list(ok = FALSE, reason = "cmprsk_not_available_or_disabled", fit = NULL, design_columns = character(0)))
  }

  transition_events <- sum(df$status_num == 1L)
  if (transition_events < minimum_transition_events) {
    return(list(ok = FALSE, reason = "too_few_transition_events", fit = NULL, design_columns = character(0)))
  }

  mm <- make_design_matrix(df, formula_row$formula_rhs[[1]])

  fit <- tryCatch(
    cmprsk::crr(
      ftime = df$time_year,
      fstatus = df$status_num,
      cov1 = mm,
      failcode = subdistribution_failcode,
      cencode = subdistribution_cencode
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(ok = FALSE, reason = conditionMessage(fit), fit = NULL, design_columns = colnames(mm)))
  }

  list(ok = TRUE, reason = NA_character_, fit = fit, design_columns = colnames(mm))
}

predict_fine_gray_transition <- function(fg_fit, formula_rhs, newdata, horizons, training_columns) {
  if (is.null(fg_fit)) {
    return(list(
      transition_matrix = matrix(0, nrow = nrow(newdata), ncol = length(horizons)),
      remission_matrix = matrix(NA_real_, nrow = nrow(newdata), ncol = length(horizons)),
      eventfree_matrix = matrix(NA_real_, nrow = nrow(newdata), ncol = length(horizons))
    ))
  }

  mm_new <- make_design_matrix(newdata, formula_rhs, training_columns = training_columns)
  pred_raw <- stats::predict(fg_fit, cov1 = mm_new)
  if (is.null(dim(pred_raw))) {
    pred_raw <- matrix(pred_raw, ncol = 2L)
  }

  pred_times <- as.numeric(pred_raw[, 1])
  pred_values <- pred_raw[, -1, drop = FALSE]

  transition_matrix <- matrix(0, nrow = nrow(newdata), ncol = length(horizons))
  for (i in seq_len(nrow(newdata))) {
    transition_matrix[i, ] <- step_eval(pred_times, pred_values[, i], horizons, left_value = 0)
  }

  list(
    transition_matrix = transition_matrix,
    remission_matrix = matrix(NA_real_, nrow = nrow(newdata), ncol = length(horizons)),
    eventfree_matrix = matrix(NA_real_, nrow = nrow(newdata), ncol = length(horizons))
  )
}

## 🟠 Define: subject-level prediction assembly ===============================
make_prediction_long <- function(
    df,
    dataset_name,
    formula_row,
    horizons_tbl,
    main_risk_matrix,
    sensitivity_transition_matrix,
    sensitivity_remission_matrix,
    sensitivity_eventfree_matrix,
    comparison_scope,
    comparison_family,
    model_class,
    main_model_family,
    sensitivity_model_family,
    main_variant_id,
    main_variant_label,
    sensitivity_variant_id,
    sensitivity_variant_label,
    sensitivity_method,
    main_risk_estimator,
    sensitivity_risk_estimator
) {
  index_tbl <- tidyr::expand_grid(
    row_id = seq_len(nrow(df)),
    horizon_index = seq_len(nrow(horizons_tbl))
  ) %>%
    mutate(
      unique_person_id = df$unique_person_id[row_id],
      id = df$id[row_id],
      site = df$site[row_id],
      sex_num = df$sex_num[row_id],
      sex_label = df$sex_label[row_id],
      age_exact_entry = df$age_exact_entry[row_id],
      age_s = df$age_s[row_id],
      days_followup = df$days_followup[row_id],
      time_year = df$time_year[row_id],
      status_num = df$status_num[row_id],
      status_label = df$status_label[row_id],
      horizon_id = horizons_tbl$horizon_id[horizon_index],
      horizon_year = horizons_tbl$horizon_year[horizon_index],
      horizon_days = horizons_tbl$horizon_days[horizon_index],
      interpretation_tier = horizons_tbl$interpretation_tier[horizon_index],
      primary_supported_flag = horizons_tbl$primary_supported_flag[horizon_index],
      interpretation_note = horizons_tbl$interpretation_note[horizon_index],
      support_tier = horizons_tbl$support_tier[horizon_index],
      support_tier_standard = horizons_tbl$support_tier_standard[horizon_index],
      horizon_evidence_class = horizons_tbl$horizon_evidence_class[horizon_index],
      claim_restriction_flag = horizons_tbl$claim_restriction_flag[horizon_index],
      horizon_support_label = horizons_tbl$horizon_support_label[horizon_index],
      supported_for_reporting = horizons_tbl$supported_for_reporting[horizon_index],
      projection_flag = horizons_tbl$projection_flag[horizon_index],
      reporting_status = horizons_tbl$reporting_status[horizon_index]
    )

  main_tbl <- index_tbl %>%
    transmute(
      dataset = dataset_name,
      comparison_scope = comparison_scope,
      comparison_family = comparison_family,
      comparison_id = formula_row$formula_id[[1]],
      prediction_level = "individualized",
      model_class = model_class,
      model_family = main_model_family,
      risk_scale = supplementary_risk_scale,
      formula_id = formula_row$formula_id[[1]],
      formula_name = formula_row$formula_name[[1]],
      formula_label = formula_row$formula_label[[1]],
      site_branch = formula_row$site_branch[[1]],
      interaction_branch = formula_row$interaction_branch[[1]],
      benchmark_core_id = NA_character_,
      benchmark_label = NA_character_,
      assignment_variable = "individualized_prediction",
      variant_id = main_variant_id,
      variant_label = main_variant_label,
      sensitivity_method = sensitivity_method,
      analysis_variant = "main_censoring_mirror",
      evaluation_mode = "main_censoring",
      risk_estimator = main_risk_estimator,
      unique_person_id,
      id,
      site,
      sex_num,
      sex_label,
      age_exact_entry,
      age_s,
      days_followup,
      time_year,
      status_num,
      status_label,
      horizon_id,
      horizon_year,
      horizon_days,
      interpretation_tier,
      primary_supported_flag,
      interpretation_note,
      support_tier,
      support_tier_standard,
      horizon_evidence_class,
      claim_restriction_flag,
      horizon_support_label,
      supported_for_reporting,
      projection_flag,
      reporting_status,
      predicted_transition_risk = main_risk_matrix[cbind(row_id, horizon_index)],
      predicted_remission_risk = 0,
      predicted_eventfree_prob = 1 - predicted_transition_risk,
      risk_difference_from_main_signed = 0,
      risk_difference_from_main_absolute = 0
    ) %>%
    bind_component_annotations()

  sensitivity_tbl <- index_tbl %>%
    transmute(
      dataset = dataset_name,
      comparison_scope = comparison_scope,
      comparison_family = comparison_family,
      comparison_id = formula_row$formula_id[[1]],
      prediction_level = "individualized",
      model_class = model_class,
      model_family = sensitivity_model_family,
      risk_scale = supplementary_risk_scale,
      formula_id = formula_row$formula_id[[1]],
      formula_name = formula_row$formula_name[[1]],
      formula_label = formula_row$formula_label[[1]],
      site_branch = formula_row$site_branch[[1]],
      interaction_branch = formula_row$interaction_branch[[1]],
      benchmark_core_id = NA_character_,
      benchmark_label = NA_character_,
      assignment_variable = "individualized_prediction",
      variant_id = sensitivity_variant_id,
      variant_label = sensitivity_variant_label,
      sensitivity_method = sensitivity_method,
      analysis_variant = "remission_sensitive",
      evaluation_mode = "remission_sensitive",
      risk_estimator = sensitivity_risk_estimator,
      unique_person_id,
      id,
      site,
      sex_num,
      sex_label,
      age_exact_entry,
      age_s,
      days_followup,
      time_year,
      status_num,
      status_label,
      horizon_id,
      horizon_year,
      horizon_days,
      interpretation_tier,
      primary_supported_flag,
      interpretation_note,
      support_tier,
      support_tier_standard,
      horizon_evidence_class,
      claim_restriction_flag,
      horizon_support_label,
      supported_for_reporting,
      projection_flag,
      reporting_status,
      predicted_transition_risk = sensitivity_transition_matrix[cbind(row_id, horizon_index)],
      predicted_remission_risk = sensitivity_remission_matrix[cbind(row_id, horizon_index)],
      predicted_eventfree_prob = sensitivity_eventfree_matrix[cbind(row_id, horizon_index)],
      risk_difference_from_main_signed = predicted_transition_risk - main_risk_matrix[cbind(row_id, horizon_index)],
      risk_difference_from_main_absolute = abs(risk_difference_from_main_signed)
    ) %>%
    bind_component_annotations()

  out <- bind_rows(main_tbl, sensitivity_tbl) %>%
    arrange(unique_person_id, analysis_variant, horizon_year)

  expected_rows <- 2L * nrow(df) * nrow(horizons_tbl)
  if (nrow(out) != expected_rows) {
    stop(
      sprintf("[%s:%s] Subject prediction row count mismatch. Expected %d, observed %d.", dataset_name, formula_row$formula_id[[1]], expected_rows, nrow(out)),
      call. = FALSE
    )
  }

  out
}

## 🟠 Define: modeling wrappers ===============================
run_cause_specific_block_for_dataset <- function(df, dataset_name, formula_registry, horizons_tbl, threshold_registry) {
  if (!isTRUE(run_cause_specific_models)) {
    return(list(
      risk_summary = tibble(),
      subject_predictions = tibble(),
      classification = tibble(),
      deltas = tibble(),
      fit_objects = list()
    ))
  }

  formula_rows <- formula_registry %>% filter(dataset == dataset_name)

  subject_prediction_list <- list()
  fit_object_list <- list()

  if (nrow(formula_rows) == 0L) {
    return(list(
      risk_summary = tibble(),
      subject_predictions = tibble(),
      classification = tibble(),
      deltas = tibble(),
      fit_objects = list()
    ))
  }

  for (i in seq_len(nrow(formula_rows))) {
    formula_row <- formula_rows[i, , drop = FALSE]

    log_step("Cause-specific mirror fitting: dataset=", dataset_name, ", formula_id=", formula_row$formula_id[[1]])

    main_fit_info <- fit_main_transition_cox(df, formula_row)
    if (!isTRUE(main_fit_info$ok)) {
      log_step(
        "Skipping cause-specific block for dataset=", dataset_name,
        ", formula_id=", formula_row$formula_id[[1]],
        " :: ", main_fit_info$reason,
        level = "WARN"
      )
      next
    }

    remission_fit_info <- fit_remission_cox(df, formula_row)
    if (!isTRUE(remission_fit_info$ok)) {
      log_step(
        "Remission Cox fit unavailable for dataset=", dataset_name,
        ", formula_id=", formula_row$formula_id[[1]],
        " :: ", remission_fit_info$reason,
        ". Transition CIF will be generated with remission hazard fixed at zero.",
        level = "WARN"
      )
    }

    cs_pred <- predict_cause_specific_cif(
      main_fit = main_fit_info$fit,
      remission_fit = remission_fit_info$fit,
      newdata = df,
      horizons = horizons_tbl$horizon_year
    )

    sensitivity_method_label <- if (isTRUE(remission_fit_info$ok)) {
      "cause_specific_cox"
    } else {
      "cause_specific_cox_transition_only"
    }

    subject_predictions <- make_prediction_long(
      df = df,
      dataset_name = dataset_name,
      formula_row = formula_row,
      horizons_tbl = horizons_tbl,
      main_risk_matrix = cs_pred$main_risk_matrix,
      sensitivity_transition_matrix = cs_pred$transition_matrix,
      sensitivity_remission_matrix = cs_pred$remission_matrix,
      sensitivity_eventfree_matrix = cs_pred$eventfree_matrix,
      comparison_scope = "cause_specific_formula",
      comparison_family = "cause_specific",
      model_class = "Stage9_cause_specific",
      main_model_family = "CoxPH",
      sensitivity_model_family = "CauseSpecificCox",
      main_variant_id = paste0(formula_row$formula_id[[1]], "__cox_main_mirror"),
      main_variant_label = paste("Main censoring mirror:", formula_row$formula_label[[1]]),
      sensitivity_variant_id = paste0(formula_row$formula_id[[1]], "__cause_specific"),
      sensitivity_variant_label = paste("Cause-specific sensitivity:", formula_row$formula_label[[1]]),
      sensitivity_method = sensitivity_method_label,
      main_risk_estimator = "CoxPH_main_mirror",
      sensitivity_risk_estimator = "CauseSpecificCox"
    )

    subject_prediction_list[[length(subject_prediction_list) + 1L]] <- subject_predictions
    fit_object_list[[formula_row$formula_id[[1]]]] <- list(
      dataset = dataset_name,
      formula = formula_row,
      main_fit = main_fit_info$fit,
      remission_fit = remission_fit_info$fit,
      baseline_transition = cs_pred$baseline_transition,
      baseline_remission = cs_pred$baseline_remission,
      lp_transition = cs_pred$lp_transition,
      lp_remission = cs_pred$lp_remission,
      remission_fit_status = if (isTRUE(remission_fit_info$ok)) "fitted" else remission_fit_info$reason
    )
  }

  all_predictions <- compact_bind_rows(subject_prediction_list)
  validate_subject_predictions(all_predictions, "cause_specific")

  all_classification <- make_threshold_classification(all_predictions, threshold_registry)
  all_risk_summary <- make_risk_summary_from_predictions(all_predictions)
  all_deltas <- make_delta_table(all_classification)

  list(
    risk_summary = all_risk_summary,
    subject_predictions = all_predictions,
    classification = all_classification,
    deltas = all_deltas,
    fit_objects = fit_object_list
  )
}

run_subdistribution_block_for_dataset <- function(df, dataset_name, formula_registry, horizons_tbl, threshold_registry) {
  formula_rows <- formula_registry %>% filter(dataset == dataset_name)

  if (!isTRUE(run_subdistribution_models) || !isTRUE(has_cmprsk) || nrow(formula_rows) == 0L) {
    return(list(
      risk_summary = tibble(),
      subject_predictions = tibble(),
      classification = tibble(),
      deltas = tibble(),
      fit_objects = list()
    ))
  }

  subject_prediction_list <- list()
  fit_object_list <- list()

  for (i in seq_len(nrow(formula_rows))) {
    formula_row <- formula_rows[i, , drop = FALSE]

    log_step("Fine-Gray fitting: dataset=", dataset_name, ", formula_id=", formula_row$formula_id[[1]])

    main_fit_info <- fit_main_transition_cox(df, formula_row)
    if (!isTRUE(main_fit_info$ok)) {
      log_step(
        "Skipping Fine-Gray block for dataset=", dataset_name,
        ", formula_id=", formula_row$formula_id[[1]],
        " :: ", main_fit_info$reason,
        level = "WARN"
      )
      next
    }

    fg_fit_info <- fit_fine_gray_transition(df, formula_row)
    if (!isTRUE(fg_fit_info$ok)) {
      log_step(
        "Skipping Fine-Gray sensitivity for dataset=", dataset_name,
        ", formula_id=", formula_row$formula_id[[1]],
        " :: ", fg_fit_info$reason,
        level = "WARN"
      )
      next
    }

    main_pred <- predict_main_cox_risk(
      main_fit = main_fit_info$fit,
      newdata = df,
      horizons = horizons_tbl$horizon_year
    )

    fg_pred <- predict_fine_gray_transition(
      fg_fit = fg_fit_info$fit,
      formula_rhs = formula_row$formula_rhs[[1]],
      newdata = df,
      horizons = horizons_tbl$horizon_year,
      training_columns = fg_fit_info$design_columns
    )

    subject_predictions <- make_prediction_long(
      df = df,
      dataset_name = dataset_name,
      formula_row = formula_row,
      horizons_tbl = horizons_tbl,
      main_risk_matrix = main_pred$risk_matrix,
      sensitivity_transition_matrix = fg_pred$transition_matrix,
      sensitivity_remission_matrix = fg_pred$remission_matrix,
      sensitivity_eventfree_matrix = fg_pred$eventfree_matrix,
      comparison_scope = "subdistribution_formula",
      comparison_family = "subdistribution",
      model_class = "Stage9_subdistribution",
      main_model_family = "CoxPH",
      sensitivity_model_family = "FineGray",
      main_variant_id = paste0(formula_row$formula_id[[1]], "__cox_main_mirror"),
      main_variant_label = paste("Main censoring mirror:", formula_row$formula_label[[1]]),
      sensitivity_variant_id = paste0(formula_row$formula_id[[1]], "__fine_gray"),
      sensitivity_variant_label = paste("Fine-Gray sensitivity:", formula_row$formula_label[[1]]),
      sensitivity_method = "fine_gray",
      main_risk_estimator = "CoxPH_main_mirror",
      sensitivity_risk_estimator = "FineGray"
    )

    subject_prediction_list[[length(subject_prediction_list) + 1L]] <- subject_predictions
    fit_object_list[[formula_row$formula_id[[1]]]] <- list(
      dataset = dataset_name,
      formula = formula_row,
      main_fit = main_fit_info$fit,
      fine_gray_fit = fg_fit_info$fit,
      design_columns = fg_fit_info$design_columns,
      baseline_main = main_pred$baseline
    )
  }

  all_predictions <- compact_bind_rows(subject_prediction_list)
  validate_subject_predictions(all_predictions, "subdistribution")

  all_classification <- make_threshold_classification(all_predictions, threshold_registry)
  all_risk_summary <- make_risk_summary_from_predictions(all_predictions)
  all_deltas <- make_delta_table(all_classification)

  list(
    risk_summary = all_risk_summary,
    subject_predictions = all_predictions,
    classification = all_classification,
    deltas = all_deltas,
    fit_objects = fit_object_list
  )
}


# 🔴 Define: cache-aware loaders and plot builders ===============================
## 🟠 Define: cache validators ===============================
try_load_existing_table <- function(path, validator = NULL, required_cols = NULL, label = basename(path)) {
  tbl <- safe_read_existing_csv(path)
  if (is.null(tbl)) {
    return(NULL)
  }

  ok <- tryCatch({
    if (!is.null(required_cols)) {
      assert_required_columns(tbl, required_cols, paste0("existing:", label))
    }
    if (!is.null(validator)) {
      validator(tbl)
    }
    TRUE
  }, error = function(e) {
    log_step("Existing file failed validation and will be recomputed: ", label, " :: ", conditionMessage(e), level = "WARN")
    FALSE
  })

  if (!isTRUE(ok)) {
    return(NULL)
  }

  tbl
}

load_existing_stage9_fit_objects <- function(path) {
  obj <- safe_read_existing_rds(path)
  if (is.null(obj) || !is.list(obj)) {
    return(NULL)
  }
  obj
}

load_existing_benchmark_block <- function(export_path, horizon_registry, threshold_registry) {
  block <- list(
    benchmark_state_probabilities = try_load_existing_table(
      file.path(export_path, "stage9_benchmark_state_probabilities.csv"),
      validator = validate_benchmark_state_probabilities,
      label = "stage9_benchmark_state_probabilities.csv"
    ),
    benchmark_subject_predictions = try_load_existing_table(
      file.path(export_path, "stage9_benchmark_subject_horizon_predictions.csv"),
      validator = function(x) validate_subject_predictions(x, "benchmark"),
      label = "stage9_benchmark_subject_horizon_predictions.csv"
    ),
    benchmark_classification = try_load_existing_table(
      file.path(export_path, "stage9_benchmark_classification.csv"),
      validator = function(x) validate_classification_table(x, "benchmark_classification"),
      label = "stage9_benchmark_classification.csv"
    ),
    benchmark_delta_vs_main = try_load_existing_table(
      file.path(export_path, "stage9_benchmark_delta_vs_main.csv"),
      validator = function(x) validate_delta_table(x, "benchmark_delta_vs_main"),
      label = "stage9_benchmark_delta_vs_main.csv"
    )
  )

  block_ok <- all(vapply(block, function(x) !is.null(x), logical(1)))
  if (!isTRUE(block_ok)) {
    return(NULL)
  }

  tryCatch({
    validate_horizon_alignment(block$benchmark_state_probabilities, horizon_registry, "benchmark_state_probabilities")
    validate_horizon_alignment(block$benchmark_subject_predictions, horizon_registry, "benchmark_subject_predictions")
    validate_horizon_alignment(block$benchmark_classification, horizon_registry, "benchmark_classification")
    validate_horizon_alignment(block$benchmark_delta_vs_main, horizon_registry, "benchmark_delta_vs_main")
    validate_threshold_alignment(block$benchmark_classification, threshold_registry, "benchmark_classification")
    validate_threshold_alignment(block$benchmark_delta_vs_main, threshold_registry, "benchmark_delta_vs_main")
    TRUE
  }, error = function(e) {
    log_step("Existing benchmark block does not match current Stage 1 registries and will be recomputed. :: ", conditionMessage(e), level = "WARN")
    FALSE
  }) -> aligned_ok

  if (!isTRUE(aligned_ok)) {
    return(NULL)
  }

  block
}

load_existing_formula_block <- function(export_path, prefix, horizon_registry, threshold_registry, formula_registry) {
  risk_summary_name <- paste0("stage9_", prefix, "_risk_summary.csv")
  subject_name <- paste0("stage9_", prefix, "_subject_horizon_predictions.csv")
  class_name <- paste0("stage9_", prefix, "_classification.csv")
  delta_name <- paste0("stage9_", prefix, "_delta_vs_main.csv")

  block <- list(
    risk_summary = try_load_existing_table(
      file.path(export_path, risk_summary_name),
      validator = function(x) validate_risk_summary_table(x, paste0(prefix, "_risk_summary")),
      label = risk_summary_name
    ),
    subject_predictions = try_load_existing_table(
      file.path(export_path, subject_name),
      validator = function(x) validate_subject_predictions(x, prefix),
      label = subject_name
    ),
    classification = try_load_existing_table(
      file.path(export_path, class_name),
      validator = function(x) validate_classification_table(x, paste0(prefix, "_classification")),
      label = class_name
    ),
    delta_vs_main = try_load_existing_table(
      file.path(export_path, delta_name),
      validator = function(x) validate_delta_table(x, paste0(prefix, "_delta_vs_main")),
      label = delta_name
    )
  )

  block_ok <- all(vapply(block, function(x) !is.null(x), logical(1)))
  if (!isTRUE(block_ok)) {
    return(NULL)
  }

  tryCatch({
    validate_horizon_alignment(block$risk_summary, horizon_registry, paste0(prefix, "_risk_summary"))
    validate_horizon_alignment(block$subject_predictions, horizon_registry, paste0(prefix, "_subject_predictions"))
    validate_horizon_alignment(block$classification, horizon_registry, paste0(prefix, "_classification"))
    validate_horizon_alignment(block$delta_vs_main, horizon_registry, paste0(prefix, "_delta_vs_main"))
    validate_threshold_alignment(block$classification, threshold_registry, paste0(prefix, "_classification"))
    validate_threshold_alignment(block$delta_vs_main, threshold_registry, paste0(prefix, "_delta_vs_main"))
    validate_formula_alignment(block$risk_summary, formula_registry, paste0(prefix, "_risk_summary"))
    validate_formula_alignment(block$subject_predictions, formula_registry, paste0(prefix, "_subject_predictions"))
    TRUE
  }, error = function(e) {
    log_step("Existing ", prefix, " block does not match current Stage 1 registries and will be recomputed. :: ", conditionMessage(e), level = "WARN")
    FALSE
  }) -> aligned_ok

  if (!isTRUE(aligned_ok)) {
    return(NULL)
  }

  block
}

## 🟠 Define: plot constructors ===============================
build_stage9_plot_objects <- function(
    benchmark_state_probabilities,
    benchmark_delta_vs_main,
    cause_specific_risk_summary,
    cause_specific_delta_vs_main,
    subdistribution_risk_summary,
    subdistribution_delta_vs_main,
    formula_registry,
    threshold_registry,
    dataset_order
) {
  plot_list <- list()

  if (nrow(benchmark_state_probabilities) > 0L) {
    benchmark_panel_levels <- make_benchmark_risk_panel_levels(
      state_probabilities_tbl = benchmark_state_probabilities,
      dataset_order = dataset_order
    )

    benchmark_plot_data <- benchmark_state_probabilities %>%
      select(
        dataset, comparison_id, benchmark_label, assignment_value,
        horizon_year, claim_restriction_flag,
        main_km_risk, sensitivity_transition_risk, sensitivity_remission_risk
      ) %>%
      pivot_longer(
        cols = c(main_km_risk, sensitivity_transition_risk, sensitivity_remission_risk),
        names_to = "curve_type",
        values_to = "value"
      ) %>%
      mutate(
        dataset = make_dataset_factor(dataset, dataset_order),
        panel_label = make_benchmark_risk_panel_label(dataset, comparison_id, assignment_value),
        panel_label = make_panel_factor(panel_label, levels = benchmark_panel_levels),
        curve_type = factor(
          curve_type,
          levels = c("main_km_risk", "sensitivity_transition_risk", "sensitivity_remission_risk"),
          labels = c("Main KM risk", "Transition CIF", "Remission CIF")
        )
      )

    benchmark_projection_rect <- build_projection_rectangles(benchmark_plot_data)

    plot_list$benchmark_risk_comparison <- ggplot(
      benchmark_plot_data,
      aes(x = horizon_year, y = value, group = curve_type, linetype = curve_type)
    ) +
      {
        if (nrow(benchmark_projection_rect) > 0L) {
          geom_rect(
            data = benchmark_projection_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey95"
          )
        } else {
          NULL
        }
      } +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      facet_wrap(~ panel_label, ncol = 2) +
      labs(
        title = "Stage 9 group-assigned benchmark risk comparison",
        x = "Horizon (years)",
        y = "Risk / cumulative incidence",
        linetype = "Curve",
        caption = projection_plot_caption
      ) +
      theme_bw()
  }

  if (nrow(benchmark_delta_vs_main) > 0L) {
    benchmark_delta_panel_levels <- make_benchmark_delta_panel_levels(
      delta_tbl = benchmark_delta_vs_main,
      dataset_order = dataset_order
    )

    benchmark_delta_plot_data <- benchmark_delta_vs_main %>%
      mutate(
        dataset = make_dataset_factor(dataset, dataset_order),
        comparison_id = factor(as.character(comparison_id), levels = c("overall", "site_adjusted")),
        threshold_label = make_threshold_factor(threshold_label, threshold_registry),
        panel_label = make_benchmark_delta_panel_label(dataset, comparison_id),
        panel_label = make_panel_factor(panel_label, levels = benchmark_delta_panel_levels)
      )

    benchmark_delta_projection_rect <- build_projection_rectangles(benchmark_delta_plot_data)

    plot_list$benchmark_false_positive_burden_change <- ggplot(
      benchmark_delta_plot_data,
      aes(x = horizon_year, y = false_positive_burden_change, group = 1)
    ) +
      {
        if (nrow(benchmark_delta_projection_rect) > 0L) {
          geom_rect(
            data = benchmark_delta_projection_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey95"
          )
        } else {
          NULL
        }
      } +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.5) +
      facet_grid(panel_label ~ threshold_label) +
      labs(
        title = "Group-assigned benchmark false-positive burden change versus main mirror",
        x = "Horizon (years)",
        y = "False-positive burden change",
        caption = projection_plot_caption
      ) +
      theme_bw()

    plot_list$benchmark_net_benefit_change <- ggplot(
      benchmark_delta_plot_data,
      aes(x = horizon_year, y = net_benefit_change, group = 1)
    ) +
      {
        if (nrow(benchmark_delta_projection_rect) > 0L) {
          geom_rect(
            data = benchmark_delta_projection_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey95"
          )
        } else {
          NULL
        }
      } +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.5) +
      facet_grid(panel_label ~ threshold_label) +
      labs(
        title = "Group-assigned benchmark net benefit change versus main mirror",
        x = "Horizon (years)",
        y = "Net benefit change",
        caption = projection_plot_caption
      ) +
      theme_bw()
  }

  if (nrow(cause_specific_risk_summary) > 0L) {
    cause_panel_levels <- make_formula_panel_levels(
      summary_tbl = cause_specific_risk_summary,
      formula_registry = formula_registry,
      dataset_order = dataset_order
    )

    cause_plot_data <- cause_specific_risk_summary %>%
      select(
        dataset, formula_label, horizon_year, claim_restriction_flag,
        sensitivity_remission_component_available,
        mean_main_transition_risk,
        mean_sensitivity_transition_risk,
        mean_sensitivity_remission_risk
      ) %>%
      pivot_longer(
        cols = c(mean_main_transition_risk, mean_sensitivity_transition_risk, mean_sensitivity_remission_risk),
        names_to = "curve_type",
        values_to = "value"
      ) %>%
      mutate(
        value = dplyr::case_when(
          curve_type == "mean_sensitivity_remission_risk" & !sensitivity_remission_component_available ~ NA_real_,
          TRUE ~ value
        ),
        dataset = make_dataset_factor(dataset, dataset_order),
        panel_label = make_formula_panel_label(dataset, formula_label),
        panel_label = make_panel_factor(panel_label, levels = cause_panel_levels),
        curve_type = factor(
          curve_type,
          levels = c("mean_main_transition_risk", "mean_sensitivity_transition_risk", "mean_sensitivity_remission_risk"),
          labels = c("Main Cox risk", "Cause-specific transition CIF", "Cause-specific remission CIF")
        )
      )

    cause_projection_rect <- build_projection_rectangles(cause_plot_data)

    plot_list$cause_specific_risk_comparison <- ggplot(
      cause_plot_data,
      aes(x = horizon_year, y = value, group = curve_type, linetype = curve_type)
    ) +
      {
        if (nrow(cause_projection_rect) > 0L) {
          geom_rect(
            data = cause_projection_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey95"
          )
        } else {
          NULL
        }
      } +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      facet_wrap(~ panel_label, ncol = 2) +
      labs(
        title = "Cause-specific individualized mirror risk comparison",
        x = "Horizon (years)",
        y = "Mean predicted risk / cumulative incidence",
        linetype = "Curve",
        caption = projection_plot_caption
      ) +
      theme_bw()
  }

  if (nrow(cause_specific_delta_vs_main) > 0L) {
    cause_delta_panel_levels <- make_formula_panel_levels(
      summary_tbl = cause_specific_delta_vs_main,
      formula_registry = formula_registry,
      dataset_order = dataset_order
    )

    cause_delta_plot_data <- cause_specific_delta_vs_main %>%
      mutate(
        dataset = make_dataset_factor(dataset, dataset_order),
        threshold_label = make_threshold_factor(threshold_label, threshold_registry),
        panel_label = make_formula_panel_label(dataset, formula_label),
        panel_label = make_panel_factor(panel_label, levels = cause_delta_panel_levels)
      )

    cause_delta_projection_rect <- build_projection_rectangles(cause_delta_plot_data)

    plot_list$cause_specific_false_positive_burden_change <- ggplot(
      cause_delta_plot_data,
      aes(x = horizon_year, y = false_positive_burden_change, group = 1)
    ) +
      {
        if (nrow(cause_delta_projection_rect) > 0L) {
          geom_rect(
            data = cause_delta_projection_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey95"
          )
        } else {
          NULL
        }
      } +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.5) +
      facet_grid(panel_label ~ threshold_label) +
      labs(
        title = "Cause-specific false-positive burden change versus main mirror",
        x = "Horizon (years)",
        y = "False-positive burden change",
        caption = projection_plot_caption
      ) +
      theme_bw()

    plot_list$cause_specific_net_benefit_change <- ggplot(
      cause_delta_plot_data,
      aes(x = horizon_year, y = net_benefit_change, group = 1)
    ) +
      {
        if (nrow(cause_delta_projection_rect) > 0L) {
          geom_rect(
            data = cause_delta_projection_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey95"
          )
        } else {
          NULL
        }
      } +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.5) +
      facet_grid(panel_label ~ threshold_label) +
      labs(
        title = "Cause-specific net benefit change versus main mirror",
        x = "Horizon (years)",
        y = "Net benefit change",
        caption = projection_plot_caption
      ) +
      theme_bw()
  }

  if (nrow(subdistribution_risk_summary) > 0L) {
    fg_panel_levels <- make_formula_panel_levels(
      summary_tbl = subdistribution_risk_summary,
      formula_registry = formula_registry,
      dataset_order = dataset_order
    )

    fg_plot_data <- subdistribution_risk_summary %>%
      select(
        dataset, formula_label, horizon_year, claim_restriction_flag,
        mean_main_transition_risk,
        mean_sensitivity_transition_risk
      ) %>%
      pivot_longer(
        cols = c(mean_main_transition_risk, mean_sensitivity_transition_risk),
        names_to = "curve_type",
        values_to = "value"
      ) %>%
      mutate(
        dataset = make_dataset_factor(dataset, dataset_order),
        panel_label = make_formula_panel_label(dataset, formula_label),
        panel_label = make_panel_factor(panel_label, levels = fg_panel_levels),
        curve_type = factor(
          curve_type,
          levels = c("mean_main_transition_risk", "mean_sensitivity_transition_risk"),
          labels = c("Main Cox risk", "Fine-Gray transition CIF")
        )
      )

    fg_projection_rect <- build_projection_rectangles(fg_plot_data)

    plot_list$fine_gray_risk_comparison <- ggplot(
      fg_plot_data,
      aes(x = horizon_year, y = value, group = curve_type, linetype = curve_type)
    ) +
      {
        if (nrow(fg_projection_rect) > 0L) {
          geom_rect(
            data = fg_projection_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey95"
          )
        } else {
          NULL
        }
      } +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      facet_wrap(~ panel_label, ncol = 2) +
      labs(
        title = "Fine-Gray individualized mirror risk comparison",
        x = "Horizon (years)",
        y = "Mean predicted transition risk",
        linetype = "Curve",
        caption = projection_plot_caption
      ) +
      theme_bw()
  }

  if (nrow(subdistribution_delta_vs_main) > 0L) {
    subdistribution_delta_panel_levels <- make_formula_panel_levels(
      summary_tbl = subdistribution_delta_vs_main,
      formula_registry = formula_registry,
      dataset_order = dataset_order
    )

    subdistribution_delta_plot_data <- subdistribution_delta_vs_main %>%
      mutate(
        dataset = make_dataset_factor(dataset, dataset_order),
        threshold_label = make_threshold_factor(threshold_label, threshold_registry),
        panel_label = make_formula_panel_label(dataset, formula_label),
        panel_label = make_panel_factor(panel_label, levels = subdistribution_delta_panel_levels)
      )

    subdistribution_delta_projection_rect <- build_projection_rectangles(subdistribution_delta_plot_data)

    plot_list$fine_gray_false_positive_burden_change <- ggplot(
      subdistribution_delta_plot_data,
      aes(x = horizon_year, y = false_positive_burden_change, group = 1)
    ) +
      {
        if (nrow(subdistribution_delta_projection_rect) > 0L) {
          geom_rect(
            data = subdistribution_delta_projection_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey95"
          )
        } else {
          NULL
        }
      } +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.5) +
      facet_grid(panel_label ~ threshold_label) +
      labs(
        title = "Fine-Gray false-positive burden change versus main mirror",
        x = "Horizon (years)",
        y = "False-positive burden change",
        caption = projection_plot_caption
      ) +
      theme_bw()

    plot_list$fine_gray_net_benefit_change <- ggplot(
      subdistribution_delta_plot_data,
      aes(x = horizon_year, y = net_benefit_change, group = 1)
    ) +
      {
        if (nrow(subdistribution_delta_projection_rect) > 0L) {
          geom_rect(
            data = subdistribution_delta_projection_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey95"
          )
        } else {
          NULL
        }
      } +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.5) +
      facet_grid(panel_label ~ threshold_label) +
      labs(
        title = "Fine-Gray net benefit change versus main mirror",
        x = "Horizon (years)",
        y = "Net benefit change",
        caption = projection_plot_caption
      ) +
      theme_bw()
  }

  plot_list
}

save_stage9_plots <- function(plot_list, export_path, make_plot_pdf, make_plot_png, overwrite_existing_plots) {
  png_files <- character(0)
  pdf_file <- NULL

  if (length(plot_list) == 0L) {
    return(list(pdf_file = pdf_file, png_files = png_files))
  }

  if (isTRUE(make_plot_pdf)) {
    pdf_file <- file.path(export_path, "stage9_plots.pdf")
    grDevices::pdf(pdf_file, width = 11, height = 8.5)
    for (nm in names(plot_list)) {
      print(plot_list[[nm]])
    }
    grDevices::dev.off()
  }

  if (isTRUE(make_plot_png)) {
    for (nm in names(plot_list)) {
      png_path <- file.path(export_path, paste0("stage9_plot__", nm, ".png"))
      if (!isTRUE(overwrite_existing_plots) && file.exists(png_path)) {
        png_files <- c(png_files, basename(png_path))
        next
      }
      ggplot2::ggsave(
        filename = png_path,
        plot = plot_list[[nm]],
        width = plot_png_width,
        height = plot_png_height,
        dpi = plot_png_dpi
      )
      png_files <- c(png_files, basename(png_path))
    }
  }

  list(
    pdf_file = if (!is.null(pdf_file)) basename(pdf_file) else NULL,
    png_files = png_files
  )
}


# 🔴 Define: metadata and manifest builders ===============================
## 🟠 Define: manifest builders ===============================
build_stage9_export_manifest <- function(plot_pdf_file, plot_png_files, save_fit_objects) {
  core_manifest <- tibble(
    file_name = c(
      "stage9_status_summary.csv",
      "stage9_benchmark_state_probabilities.csv",
      "stage9_benchmark_subject_horizon_predictions.csv",
      "stage9_benchmark_classification.csv",
      "stage9_benchmark_delta_vs_main.csv",
      "stage9_cause_specific_risk_summary.csv",
      "stage9_cause_specific_subject_horizon_predictions.csv",
      "stage9_cause_specific_classification.csv",
      "stage9_cause_specific_delta_vs_main.csv",
      "stage9_subdistribution_risk_summary.csv",
      "stage9_subdistribution_subject_horizon_predictions.csv",
      "stage9_subdistribution_classification.csv",
      "stage9_subdistribution_delta_vs_main.csv",
      "stage9_stage_metadata.csv",
      "stage9_progress_log.txt",
      "stage9_output_bundle.rds",
      "stage9_export_manifest.csv"
    ),
    object_name = c(
      "status_summary",
      "benchmark_state_probabilities",
      "benchmark_subject_predictions",
      "benchmark_classification",
      "benchmark_delta_vs_main",
      "cause_specific_risk_summary",
      "cause_specific_subject_predictions",
      "cause_specific_classification",
      "cause_specific_delta_vs_main",
      "subdistribution_risk_summary",
      "subdistribution_subject_predictions",
      "subdistribution_classification",
      "subdistribution_delta_vs_main",
      "stage9_metadata",
      "stage9_progress_log",
      "stage9_output_bundle",
      "stage9_export_manifest"
    ),
    role = c(
      "status count source of truth",
      "group-level group-assigned benchmark source of truth",
      "subject-level group-assigned benchmark predictions",
      "threshold-based benchmark evaluation",
      "benchmark change-versus-main table",
      "cause-specific horizon risk source of truth",
      "subject-level cause-specific predictions",
      "threshold-based cause-specific evaluation",
      "cause-specific change-versus-main table",
      "subdistribution horizon risk source of truth",
      "subject-level subdistribution predictions",
      "threshold-based subdistribution evaluation",
      "subdistribution change-versus-main table",
      "stage metadata",
      "runtime log",
      "full output bundle",
      "manifest"
    ),
    file_type = c(
      rep("csv", 14),
      "txt",
      "rds",
      "csv"
    )
  )

  out <- core_manifest

  if (!is.null(plot_pdf_file)) {
    out <- bind_rows(
      out,
      tibble(
        file_name = plot_pdf_file,
        object_name = "stage9_plot_pdf",
        role = "combined stage-9 plots PDF",
        file_type = "pdf"
      )
    )
  }

  if (isTRUE(save_fit_objects)) {
    out <- bind_rows(
      out,
      tibble(
        file_name = "stage9_fit_objects.rds",
        object_name = "stage9_fit_objects",
        role = "fit objects bundle",
        file_type = "rds"
      )
    )
  }

  if (length(plot_png_files) > 0L) {
    out <- bind_rows(
      out,
      tibble(
        file_name = plot_png_files,
        object_name = paste0("plot_png__", tools::file_path_sans_ext(plot_png_files)),
        role = "individual plot PNG",
        file_type = "png"
      )
    )
  }

  out
}

## 🟠 Define: metadata builders ===============================
build_stage9_metadata <- function(
    elapsed_value,
    stage1_path,
    stage1_bundle_file,
    stage1_datasets_file,
    stage1_formula_registry_file,
    stage1_horizon_registry_file,
    stage1_threshold_registry_file,
    stage1_modeling_registry_file,
    stage1_metadata_registry_file,
    stage1_export_manifest_file,
    stage1_file_signature,
    horizon_registry,
    threshold_registry,
    has_cmprsk,
    run_subdistribution_models,
    run_merged_site_adjusted_benchmark,
    make_plot_pdf,
    make_plot_png,
    save_fit_objects,
    block_reuse_registry,
    plot_pdf_file,
    plot_png_files
) {
  core_rows <- bind_rows(
    tibble(
      metadata_group = "stage",
      metadata_name = c("stage_name", "stage_role", "model_fitting_allowed"),
      metadata_value = c(
        "Stage 9",
        "Reanalyze remission as a dedicated sensitivity block without hiding remission inside censoring.",
        "Yes"
      )
    ),
    tibble(
      metadata_group = "inputs",
      metadata_name = c(
        "stage1_path", "stage1_bundle_file", "stage1_datasets_file", "stage1_formula_registry_file",
        "stage1_horizon_registry_file", "stage1_threshold_registry_file", "stage1_modeling_registry_file",
        "stage1_metadata_registry_file", "stage1_export_manifest_file"
      ),
      metadata_value = c(
        normalize_existing_path(stage1_path),
        normalize_existing_path(stage1_bundle_file),
        normalize_existing_path(stage1_datasets_file),
        normalize_existing_path(stage1_formula_registry_file),
        normalize_existing_path(stage1_horizon_registry_file),
        normalize_existing_path(stage1_threshold_registry_file),
        normalize_existing_path(stage1_modeling_registry_file),
        normalize_existing_path(stage1_metadata_registry_file),
        normalize_existing_path(stage1_export_manifest_file)
      )
    ),
    tibble(
      metadata_group = "event",
      metadata_name = c("main_event_definition", "main_censoring_definition", "remission_sensitive_rule", "stage9_risk_scale"),
      metadata_value = c(
        "status_num == 1",
        "status_num %in% c(0, 2)",
        "status_num == 2 treated as remission competing event in dedicated sensitivity block",
        supplementary_risk_scale
      )
    ),
    tibble(
      metadata_group = "horizons",
      metadata_name = c("horizon_vector", "projection_reporting_rule"),
      metadata_value = c(
        paste(sort(unique(horizon_registry$horizon_year)), collapse = ","),
        "claim_restriction_flag controls supported versus projection interpretation; projection horizons are shaded in Stage 9 plots"
      )
    ),
    tibble(
      metadata_group = "thresholds",
      metadata_name = "threshold_vector",
      metadata_value = paste(format(sort(unique(threshold_registry$threshold)), trim = TRUE, scientific = FALSE), collapse = ",")
    ),
    tibble(
      metadata_group = "methods",
      metadata_name = c(
        "benchmark_method",
        "benchmark_prediction_level",
        "cause_specific_method",
        "subdistribution_method",
        "merged_site_adjusted_benchmark",
        "component_annotation_rule",
        "reuse_existing_valid_outputs"
      ),
      metadata_value = c(
        "Nonparametric Kaplan-Meier main mirror versus Aalen-Johansen transition/remission cumulative incidence benchmark",
        "group_assigned benchmark classifier; not an individualized regression model",
        "Main censoring-based Cox mirror versus paired cause-specific Cox cumulative incidence",
        if (isTRUE(run_subdistribution_models) && has_cmprsk) {
          "Main censoring-based Cox mirror versus Fine-Gray transition subdistribution sensitivity"
        } else {
          "Skipped because Fine-Gray is disabled or cmprsk is unavailable"
        },
        if (isTRUE(run_merged_site_adjusted_benchmark)) "Enabled" else "Disabled",
        "Component availability fields reflect model structure rather than whether remission/event-free numeric columns are missing",
        if (isTRUE(reuse_existing_valid_outputs)) "Enabled" else "Disabled"
      )
    ),
    tibble(
      metadata_group = "runtime",
      metadata_name = c("cmprsk_available", "elapsed"),
      metadata_value = c(as.character(has_cmprsk), elapsed_value)
    ),
    tibble(
      metadata_group = "outputs",
      metadata_name = c(
        "status_summary_file",
        "benchmark_state_probabilities_file",
        "benchmark_subject_predictions_file",
        "benchmark_classification_file",
        "benchmark_delta_vs_main_file",
        "cause_specific_risk_summary_file",
        "cause_specific_subject_predictions_file",
        "cause_specific_classification_file",
        "cause_specific_delta_vs_main_file",
        "subdistribution_risk_summary_file",
        "subdistribution_subject_predictions_file",
        "subdistribution_classification_file",
        "subdistribution_delta_vs_main_file",
        "stage_metadata_file",
        "export_manifest_file",
        "output_bundle_file",
        "progress_log_file"
      ),
      metadata_value = c(
        "stage9_status_summary.csv",
        "stage9_benchmark_state_probabilities.csv",
        "stage9_benchmark_subject_horizon_predictions.csv",
        "stage9_benchmark_classification.csv",
        "stage9_benchmark_delta_vs_main.csv",
        "stage9_cause_specific_risk_summary.csv",
        "stage9_cause_specific_subject_horizon_predictions.csv",
        "stage9_cause_specific_classification.csv",
        "stage9_cause_specific_delta_vs_main.csv",
        "stage9_subdistribution_risk_summary.csv",
        "stage9_subdistribution_subject_horizon_predictions.csv",
        "stage9_subdistribution_classification.csv",
        "stage9_subdistribution_delta_vs_main.csv",
        "stage9_stage_metadata.csv",
        "stage9_export_manifest.csv",
        "stage9_output_bundle.rds",
        "stage9_progress_log.txt"
      )
    )
  )

  reuse_rows <- if (nrow(block_reuse_registry) > 0L) {
    block_reuse_registry %>%
      transmute(
        metadata_group = "reuse",
        metadata_name = paste0("reuse_", block_name),
        metadata_value = ifelse(reused, "TRUE", "FALSE")
      )
  } else {
    tibble()
  }

  signature_rows <- if (nrow(stage1_file_signature) > 0L) {
    stage1_file_signature %>%
      mutate(sig_id = seq_len(n())) %>%
      transmute(
        metadata_group = "stage1_file_signature",
        metadata_name = paste0("file_", sig_id),
        metadata_value = paste(
          path,
          paste0("exists=", exists),
          paste0("size=", size_bytes),
          paste0("mtime=", modified_time),
          sep = " | "
        )
      )
  } else {
    tibble()
  }

  plot_rows <- bind_rows(
    if (!is.null(plot_pdf_file)) {
      tibble(
        metadata_group = "outputs",
        metadata_name = "plot_pdf_file",
        metadata_value = plot_pdf_file
      )
    } else {
      NULL
    },
    if (length(plot_png_files) > 0L) {
      tibble(
        metadata_group = "outputs",
        metadata_name = paste0("plot_png_", seq_along(plot_png_files)),
        metadata_value = plot_png_files
      )
    } else {
      NULL
    },
    if (isTRUE(save_fit_objects)) {
      tibble(
        metadata_group = "outputs",
        metadata_name = "fit_objects_file",
        metadata_value = "stage9_fit_objects.rds"
      )
    } else {
      NULL
    }
  )

  bind_rows(core_rows, reuse_rows, signature_rows, plot_rows)
}


# 🔴 Read: Stage-1 backbone inputs and harmonize registries ===============================
## 🟠 Read: Stage 1 reusable artifacts ===============================
stage1_bundle <- read_stage1_rds(stage1_bundle_file)
analysis_datasets_raw <- read_stage1_rds(stage1_datasets_file)
formula_registry <- read_stage1_csv(stage1_formula_registry_file)
horizon_registry <- read_stage1_csv(stage1_horizon_registry_file)
threshold_registry <- read_stage1_csv(stage1_threshold_registry_file)
modeling_registry <- read_stage1_csv(stage1_modeling_registry_file)
metadata_registry <- read_stage1_csv(stage1_metadata_registry_file)
stage1_export_manifest <- if (file.exists(stage1_export_manifest_file)) read_stage1_csv(stage1_export_manifest_file) else tibble()

if (!is.null(stage1_bundle$config$main_risk_scale)) {
  main_risk_scale <- as.character(stage1_bundle$config$main_risk_scale)
}
if (!is.null(stage1_bundle$config$supplementary_risk_scale)) {
  supplementary_risk_scale <- as.character(stage1_bundle$config$supplementary_risk_scale)
}

analysis_datasets <- lapply(names(analysis_datasets_raw), function(dataset_name) {
  prepare_stage1_dataset(analysis_datasets_raw[[dataset_name]], dataset_name)
})
names(analysis_datasets) <- names(analysis_datasets_raw)

horizon_registry <- normalize_horizon_registry(horizon_registry)

threshold_registry <- threshold_registry %>%
  mutate(
    threshold = as.numeric(threshold),
    threshold_id = if ("threshold_id" %in% names(.)) as.character(threshold_id) else paste0("threshold_", row_number()),
    threshold_label = if ("threshold_label" %in% names(.)) as.character(threshold_label) else paste0(format(100 * threshold, trim = TRUE, scientific = FALSE), "%"),
    positive_rule = if ("positive_rule" %in% names(.)) as.character(positive_rule) else "Classify as high risk when predicted risk >= threshold"
  ) %>%
  arrange(threshold)

formula_registry <- formula_registry %>%
  mutate(
    dataset = as.character(dataset),
    formula_id = as.character(formula_id),
    formula_name = as.character(formula_name),
    formula_label = as.character(formula_label),
    formula_rhs = as.character(formula_rhs),
    site_branch = as.character(site_branch),
    interaction_branch = as.character(interaction_branch)
  ) %>%
  arrange(match(dataset, c("PNU", "SNU", "merged")), formula_id)

validate_stage1_backbone(
  stage1_bundle = stage1_bundle,
  analysis_datasets = analysis_datasets,
  formula_registry = formula_registry,
  horizon_registry = horizon_registry,
  threshold_registry = threshold_registry,
  modeling_registry = modeling_registry,
  metadata_registry = metadata_registry
)

dataset_order <- c("PNU", "SNU", "merged")
horizons_by_dataset <- split(horizon_registry, horizon_registry$dataset)
stage1_file_signature <- build_file_signature(
  c(
    stage1_bundle_file,
    stage1_datasets_file,
    stage1_formula_registry_file,
    stage1_horizon_registry_file,
    stage1_threshold_registry_file,
    stage1_modeling_registry_file,
    stage1_metadata_registry_file,
    stage1_export_manifest_file
  )
)

log_step(
  "Stage 1 backbone loaded. Datasets=",
  paste(names(analysis_datasets), collapse = ", "),
  "; thresholds=",
  paste(format(threshold_registry$threshold, trim = TRUE, scientific = FALSE), collapse = ", ")
)

# 🔴 Summarize: remission and censoring composition ===============================
status_summary <- build_status_summary(analysis_datasets)
validate_status_summary(status_summary)

# 🔴 Reuse or fit: benchmark, cause-specific, and Fine-Gray blocks ===============================
## 🟠 Read: existing fit object bundle when available ===============================
existing_fit_objects <- if (isTRUE(reuse_existing_valid_outputs) && isTRUE(save_fit_objects)) {
  load_existing_stage9_fit_objects(file.path(export_path, "stage9_fit_objects.rds"))
} else {
  NULL
}

block_reuse_registry <- tibble(block_name = character(0), reused = logical(0))

## 🟠 Resolve: benchmark block ===============================
benchmark_outputs <- NULL
benchmark_state_probabilities <- tibble()
benchmark_subject_predictions <- tibble()
benchmark_classification <- tibble()
benchmark_delta_vs_main <- tibble()
benchmark_fit_objects <- list(PNU = list(), SNU = list(), merged = list())

if (isTRUE(reuse_existing_valid_outputs)) {
  benchmark_existing <- load_existing_benchmark_block(export_path, horizon_registry = horizon_registry, threshold_registry = threshold_registry)
} else {
  benchmark_existing <- NULL
}

if (!is.null(benchmark_existing)) {
  log_step("Reusing existing valid benchmark outputs.")
  benchmark_state_probabilities <- benchmark_existing$benchmark_state_probabilities
  benchmark_subject_predictions <- benchmark_existing$benchmark_subject_predictions
  benchmark_classification <- benchmark_existing$benchmark_classification
  benchmark_delta_vs_main <- benchmark_existing$benchmark_delta_vs_main
  if (!is.null(existing_fit_objects$benchmark)) {
    benchmark_fit_objects <- existing_fit_objects$benchmark
  }
  block_reuse_registry <- bind_rows(block_reuse_registry, tibble(block_name = "benchmark", reused = TRUE))
} else {
  log_step("Building benchmark sensitivity outputs.")
  benchmark_outputs <- lapply(dataset_order, function(dataset_name) {
    log_step("Building benchmark sensitivity outputs for dataset=", dataset_name)
    build_benchmark_outputs(
      df = analysis_datasets[[dataset_name]],
      dataset_name = dataset_name,
      horizons_tbl = horizons_by_dataset[[dataset_name]],
      allow_site_adjusted = run_merged_site_adjusted_benchmark
    )
  })
  names(benchmark_outputs) <- dataset_order

  benchmark_state_probabilities <- compact_bind_rows(lapply(benchmark_outputs, `[[`, "group_horizon_table"))
  benchmark_subject_predictions <- compact_bind_rows(lapply(benchmark_outputs, `[[`, "subject_predictions"))
  benchmark_classification <- make_threshold_classification(benchmark_subject_predictions, threshold_registry)
  benchmark_delta_vs_main <- make_delta_table(benchmark_classification)

  benchmark_fit_objects <- list(
    PNU = benchmark_outputs[["PNU"]]$fit_objects,
    SNU = benchmark_outputs[["SNU"]]$fit_objects,
    merged = benchmark_outputs[["merged"]]$fit_objects
  )
  block_reuse_registry <- bind_rows(block_reuse_registry, tibble(block_name = "benchmark", reused = FALSE))
}

## 🟠 Resolve: cause-specific block ===============================
cause_specific_risk_summary <- tibble()
cause_specific_subject_predictions <- tibble()
cause_specific_classification <- tibble()
cause_specific_delta_vs_main <- tibble()
cause_specific_fit_objects <- list(PNU = list(), SNU = list(), merged = list())

if (isTRUE(reuse_existing_valid_outputs)) {
  cause_specific_existing <- load_existing_formula_block(export_path, "cause_specific", horizon_registry = horizon_registry, threshold_registry = threshold_registry, formula_registry = formula_registry)
} else {
  cause_specific_existing <- NULL
}

if (!is.null(cause_specific_existing)) {
  log_step("Reusing existing valid cause-specific outputs.")
  cause_specific_risk_summary <- cause_specific_existing$risk_summary
  cause_specific_subject_predictions <- cause_specific_existing$subject_predictions
  cause_specific_classification <- cause_specific_existing$classification
  cause_specific_delta_vs_main <- cause_specific_existing$delta_vs_main
  if (!is.null(existing_fit_objects$cause_specific)) {
    cause_specific_fit_objects <- existing_fit_objects$cause_specific
  }
  block_reuse_registry <- bind_rows(block_reuse_registry, tibble(block_name = "cause_specific", reused = TRUE))
} else {
  cause_specific_outputs <- lapply(dataset_order, function(dataset_name) {
    run_cause_specific_block_for_dataset(
      df = analysis_datasets[[dataset_name]],
      dataset_name = dataset_name,
      formula_registry = formula_registry,
      horizons_tbl = horizons_by_dataset[[dataset_name]],
      threshold_registry = threshold_registry
    )
  })
  names(cause_specific_outputs) <- dataset_order

  cause_specific_risk_summary <- compact_bind_rows(lapply(cause_specific_outputs, `[[`, "risk_summary"))
  cause_specific_subject_predictions <- compact_bind_rows(lapply(cause_specific_outputs, `[[`, "subject_predictions"))
  cause_specific_classification <- compact_bind_rows(lapply(cause_specific_outputs, `[[`, "classification"))
  cause_specific_delta_vs_main <- compact_bind_rows(lapply(cause_specific_outputs, `[[`, "deltas"))

  cause_specific_fit_objects <- list(
    PNU = cause_specific_outputs[["PNU"]]$fit_objects,
    SNU = cause_specific_outputs[["SNU"]]$fit_objects,
    merged = cause_specific_outputs[["merged"]]$fit_objects
  )
  block_reuse_registry <- bind_rows(block_reuse_registry, tibble(block_name = "cause_specific", reused = FALSE))
}

## 🟠 Resolve: subdistribution block ===============================
subdistribution_risk_summary <- tibble()
subdistribution_subject_predictions <- tibble()
subdistribution_classification <- tibble()
subdistribution_delta_vs_main <- tibble()
subdistribution_fit_objects <- list(PNU = list(), SNU = list(), merged = list())

if (isTRUE(reuse_existing_valid_outputs)) {
  subdistribution_existing <- load_existing_formula_block(export_path, "subdistribution", horizon_registry = horizon_registry, threshold_registry = threshold_registry, formula_registry = formula_registry)
} else {
  subdistribution_existing <- NULL
}

if (!is.null(subdistribution_existing)) {
  log_step("Reusing existing valid subdistribution outputs.")
  subdistribution_risk_summary <- subdistribution_existing$risk_summary
  subdistribution_subject_predictions <- subdistribution_existing$subject_predictions
  subdistribution_classification <- subdistribution_existing$classification
  subdistribution_delta_vs_main <- subdistribution_existing$delta_vs_main
  if (!is.null(existing_fit_objects$subdistribution)) {
    subdistribution_fit_objects <- existing_fit_objects$subdistribution
  }
  block_reuse_registry <- bind_rows(block_reuse_registry, tibble(block_name = "subdistribution", reused = TRUE))
} else {
  subdistribution_outputs <- lapply(dataset_order, function(dataset_name) {
    run_subdistribution_block_for_dataset(
      df = analysis_datasets[[dataset_name]],
      dataset_name = dataset_name,
      formula_registry = formula_registry,
      horizons_tbl = horizons_by_dataset[[dataset_name]],
      threshold_registry = threshold_registry
    )
  })
  names(subdistribution_outputs) <- dataset_order

  subdistribution_risk_summary <- compact_bind_rows(lapply(subdistribution_outputs, `[[`, "risk_summary"))
  subdistribution_subject_predictions <- compact_bind_rows(lapply(subdistribution_outputs, `[[`, "subject_predictions"))
  subdistribution_classification <- compact_bind_rows(lapply(subdistribution_outputs, `[[`, "classification"))
  subdistribution_delta_vs_main <- compact_bind_rows(lapply(subdistribution_outputs, `[[`, "deltas"))

  subdistribution_fit_objects <- list(
    PNU = subdistribution_outputs[["PNU"]]$fit_objects,
    SNU = subdistribution_outputs[["SNU"]]$fit_objects,
    merged = subdistribution_outputs[["merged"]]$fit_objects
  )
  block_reuse_registry <- bind_rows(block_reuse_registry, tibble(block_name = "subdistribution", reused = FALSE))
}

# 🔴 Check: final Stage-9 output integrity ===============================
validate_benchmark_state_probabilities(benchmark_state_probabilities)
validate_subject_predictions(benchmark_subject_predictions, "benchmark")
validate_classification_table(benchmark_classification, "benchmark_classification")
validate_delta_table(benchmark_delta_vs_main, "benchmark_delta_vs_main")

validate_subject_predictions(cause_specific_subject_predictions, "cause_specific")
validate_risk_summary_table(cause_specific_risk_summary, "cause_specific_risk_summary")
validate_classification_table(cause_specific_classification, "cause_specific_classification")
validate_delta_table(cause_specific_delta_vs_main, "cause_specific_delta_vs_main")

validate_subject_predictions(subdistribution_subject_predictions, "subdistribution")
validate_risk_summary_table(subdistribution_risk_summary, "subdistribution_risk_summary")
validate_classification_table(subdistribution_classification, "subdistribution_classification")
validate_delta_table(subdistribution_delta_vs_main, "subdistribution_delta_vs_main")

# 🔴 Export: core Stage-9 CSV tables ===============================
safe_write_csv(status_summary, file.path(export_path, "stage9_status_summary.csv"))
safe_write_csv(benchmark_state_probabilities, file.path(export_path, "stage9_benchmark_state_probabilities.csv"))
safe_write_csv(benchmark_subject_predictions, file.path(export_path, "stage9_benchmark_subject_horizon_predictions.csv"))
safe_write_csv(benchmark_classification, file.path(export_path, "stage9_benchmark_classification.csv"))
safe_write_csv(benchmark_delta_vs_main, file.path(export_path, "stage9_benchmark_delta_vs_main.csv"))

safe_write_csv(cause_specific_risk_summary, file.path(export_path, "stage9_cause_specific_risk_summary.csv"))
safe_write_csv(cause_specific_subject_predictions, file.path(export_path, "stage9_cause_specific_subject_horizon_predictions.csv"))
safe_write_csv(cause_specific_classification, file.path(export_path, "stage9_cause_specific_classification.csv"))
safe_write_csv(cause_specific_delta_vs_main, file.path(export_path, "stage9_cause_specific_delta_vs_main.csv"))

safe_write_csv(subdistribution_risk_summary, file.path(export_path, "stage9_subdistribution_risk_summary.csv"))
safe_write_csv(subdistribution_subject_predictions, file.path(export_path, "stage9_subdistribution_subject_horizon_predictions.csv"))
safe_write_csv(subdistribution_classification, file.path(export_path, "stage9_subdistribution_classification.csv"))
safe_write_csv(subdistribution_delta_vs_main, file.path(export_path, "stage9_subdistribution_delta_vs_main.csv"))

# 🔴 Plot: PDF and plot-specific PNG exports ===============================
stage9_plot_objects <- build_stage9_plot_objects(
  benchmark_state_probabilities = benchmark_state_probabilities,
  benchmark_delta_vs_main = benchmark_delta_vs_main,
  cause_specific_risk_summary = cause_specific_risk_summary,
  cause_specific_delta_vs_main = cause_specific_delta_vs_main,
  subdistribution_risk_summary = subdistribution_risk_summary,
  subdistribution_delta_vs_main = subdistribution_delta_vs_main,
  formula_registry = formula_registry,
  threshold_registry = threshold_registry,
  dataset_order = dataset_order
)

plot_export_info <- save_stage9_plots(
  plot_list = stage9_plot_objects,
  export_path = export_path,
  make_plot_pdf = make_plot_pdf,
  make_plot_png = make_plot_png,
  overwrite_existing_plots = overwrite_existing_plots
)

# 🔴 Build: fit-object bundle, metadata, and output bundle ===============================
stage9_fit_objects <- list(
  benchmark = benchmark_fit_objects,
  cause_specific = cause_specific_fit_objects,
  subdistribution = subdistribution_fit_objects
)

if (isTRUE(save_fit_objects)) {
  safe_save_rds(stage9_fit_objects, file.path(export_path, "stage9_fit_objects.rds"))
}

stage9_elapsed <- format_elapsed(stage_start_time)

stage9_metadata <- build_stage9_metadata(
  elapsed_value = stage9_elapsed,
  stage1_path = stage1_path,
  stage1_bundle_file = stage1_bundle_file,
  stage1_datasets_file = stage1_datasets_file,
  stage1_formula_registry_file = stage1_formula_registry_file,
  stage1_horizon_registry_file = stage1_horizon_registry_file,
  stage1_threshold_registry_file = stage1_threshold_registry_file,
  stage1_modeling_registry_file = stage1_modeling_registry_file,
  stage1_metadata_registry_file = stage1_metadata_registry_file,
  stage1_export_manifest_file = stage1_export_manifest_file,
  stage1_file_signature = stage1_file_signature,
  horizon_registry = horizon_registry,
  threshold_registry = threshold_registry,
  has_cmprsk = has_cmprsk,
  run_subdistribution_models = run_subdistribution_models,
  run_merged_site_adjusted_benchmark = run_merged_site_adjusted_benchmark,
  make_plot_pdf = make_plot_pdf,
  make_plot_png = make_plot_png,
  save_fit_objects = save_fit_objects,
  block_reuse_registry = block_reuse_registry,
  plot_pdf_file = plot_export_info$pdf_file,
  plot_png_files = plot_export_info$png_files
)

stage9_export_manifest <- build_stage9_export_manifest(
  plot_pdf_file = plot_export_info$pdf_file,
  plot_png_files = plot_export_info$png_files,
  save_fit_objects = save_fit_objects
)

stage9_output_bundle <- list(
  stage = "Stage 9",
  created_at = as.character(Sys.time()),
  session_info = utils::sessionInfo(),
  config = list(
    stage1_path = normalize_existing_path(stage1_path),
    export_path = normalize_existing_path(export_path),
    run_merged_site_adjusted_benchmark = run_merged_site_adjusted_benchmark,
    run_cause_specific_models = run_cause_specific_models,
    run_subdistribution_models = run_subdistribution_models,
    reuse_existing_valid_outputs = reuse_existing_valid_outputs,
    overwrite_existing_plots = overwrite_existing_plots,
    has_cmprsk = has_cmprsk,
    minimum_transition_events = minimum_transition_events,
    minimum_remission_events = minimum_remission_events,
    cox_ties = cox_ties,
    main_risk_scale = main_risk_scale,
    supplementary_risk_scale = supplementary_risk_scale
  ),
  stage1_inputs = list(
    bundle = stage1_bundle,
    formula_registry = formula_registry,
    horizon_registry = horizon_registry,
    threshold_registry = threshold_registry,
    modeling_registry = modeling_registry,
    metadata_registry = metadata_registry,
    stage1_export_manifest = stage1_export_manifest
  ),
  reuse_registry = block_reuse_registry,
  outputs = list(
    status_summary = status_summary,
    benchmark_state_probabilities = benchmark_state_probabilities,
    benchmark_subject_predictions = benchmark_subject_predictions,
    benchmark_classification = benchmark_classification,
    benchmark_delta_vs_main = benchmark_delta_vs_main,
    cause_specific_risk_summary = cause_specific_risk_summary,
    cause_specific_subject_predictions = cause_specific_subject_predictions,
    cause_specific_classification = cause_specific_classification,
    cause_specific_delta_vs_main = cause_specific_delta_vs_main,
    subdistribution_risk_summary = subdistribution_risk_summary,
    subdistribution_subject_predictions = subdistribution_subject_predictions,
    subdistribution_classification = subdistribution_classification,
    subdistribution_delta_vs_main = subdistribution_delta_vs_main,
    stage9_metadata = stage9_metadata,
    stage9_export_manifest = stage9_export_manifest,
    stage9_plot_pdf = plot_export_info$pdf_file,
    stage9_plot_png_files = plot_export_info$png_files,
    stage9_fit_objects_file = if (isTRUE(save_fit_objects)) "stage9_fit_objects.rds" else NULL
  )
)

safe_write_csv(stage9_metadata, file.path(export_path, "stage9_stage_metadata.csv"))
safe_write_csv(stage9_export_manifest, file.path(export_path, "stage9_export_manifest.csv"))
safe_save_rds(stage9_output_bundle, file.path(export_path, "stage9_output_bundle.rds"))

validate_export_manifest(stage9_export_manifest, export_path)

log_step("Stage 9 completed in ", stage9_elapsed, ".")
