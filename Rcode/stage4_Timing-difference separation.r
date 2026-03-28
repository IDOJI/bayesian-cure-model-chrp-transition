rm(list = ls())

# 🔴 Configure: paths and Stage 4 constants ===============================
data_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.데이터 분석 New/🟦Stage1/attachments'
export_path <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.데이터 분석 New/🟦Stage4_PNU–SNU timing-difference/attachments"

stage1_analysis_datasets_file <- file.path(data_path, "stage1_analysis_datasets.rds")
stage1_backbone_bundle_file <- file.path(data_path, "stage1_backbone_bundle.rds")
stage1_dataset_registry_file <- file.path(data_path, "stage1_dataset_registry.csv")
stage1_scaling_registry_file <- file.path(data_path, "stage1_scaling_registry.csv")
stage1_metadata_registry_file <- file.path(data_path, "stage1_metadata_registry.csv")
stage1_formula_registry_file <- file.path(data_path, "stage1_formula_registry.csv")
stage1_horizon_registry_file <- file.path(data_path, "stage1_horizon_registry.csv")
stage1_threshold_registry_file <- file.path(data_path, "stage1_threshold_registry.csv")

bootstrap_iterations <- 300L
bootstrap_seed <- 20260322L
bootstrap_min_success_rate <- 0.90

time_zero_epsilon <- 1e-08
piecewise_cuts_year <- c(1, 2)
hazard_band_breaks_year <- c(0, 1, 2, 5, 10)
supported_short_horizons_year <- c(1, 2)
late_tail_horizons_year <- c(5, 10)

piecewise_min_events_per_site <- 1L
piecewise_min_subjects_per_site <- 5L

pnu_site_label <- "PNU"
snu_site_label <- "SNU"

# 🔴 Attach: packages and session options ===============================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(survival)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# 🔴 Define: reusable Stage 4 helpers ===============================
## 🟠 Define: path readers and scalar utilities ===============================
normalize_existing_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Required CSV file does not exist: %s", path), call. = FALSE)
  }
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

safe_divide <- function(num, den) {
  ifelse(is.na(den) | den <= 0, NA_real_, num / den)
}

coalesce_character <- function(x, y) {
  ifelse(is.na(x) | x == "", y, x)
}

combine_reason_tokens <- function(tokens) {
  tokens <- tokens[!is.na(tokens) & nzchar(tokens)]
  if (length(tokens) == 0L) {
    return(NA_character_)
  }
  paste(unique(tokens), collapse = ";")
}

equal_numeric_or_both_na <- function(x, y) {
  (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
}

format_check_locations <- function(df, cols) {
  if (nrow(df) == 0L) {
    return("ok")
  }
  apply(
    X = df[, cols, drop = FALSE],
    MARGIN = 1,
    FUN = function(row) paste(row, collapse = "/")
  ) %>%
    paste(collapse = "; ")
}

make_portable_file_reference <- function(path) {
  basename(path)
}

make_reporting_status <- function(estimable_flag, reported_estimate_suppressed_flag, bootstrap_instability_flag = NA) {
  dplyr::case_when(
    dplyr::coalesce(as.logical(reported_estimate_suppressed_flag), FALSE) ~ "suppressed",
    !dplyr::coalesce(as.logical(estimable_flag), FALSE) ~ "suppressed",
    dplyr::coalesce(as.logical(bootstrap_instability_flag), FALSE) ~ "reported_with_bootstrap_instability_flag",
    TRUE ~ "reported"
  )
}

make_plotting_status <- function(estimable_flag, reported_estimate_suppressed_flag, bootstrap_instability_flag = NA) {
  dplyr::case_when(
    dplyr::coalesce(as.logical(reported_estimate_suppressed_flag), FALSE) ~ "suppressed",
    !dplyr::coalesce(as.logical(estimable_flag), FALSE) ~ "suppressed",
    dplyr::coalesce(as.logical(bootstrap_instability_flag), FALSE) ~ "bootstrap_unstable",
    TRUE ~ "stable"
  )
}

make_reporting_note <- function(reporting_status, availability_note, bootstrap_note = NA_character_, support_issue_reason = NA_character_) {
  dplyr::case_when(
    reporting_status == "reported_with_bootstrap_instability_flag" ~ coalesce_character(bootstrap_note, "bootstrap_instability_flagged"),
    reporting_status == "suppressed" ~ coalesce_character(support_issue_reason, availability_note),
    TRUE ~ availability_note
  )
}

read_stage1_inputs <- function() {
  backbone_bundle <- if (file.exists(stage1_backbone_bundle_file)) readRDS(stage1_backbone_bundle_file) else NULL
  analysis_datasets <- if (file.exists(stage1_analysis_datasets_file)) {
    readRDS(stage1_analysis_datasets_file)
  } else if (!is.null(backbone_bundle) && !is.null(backbone_bundle$datasets)) {
    backbone_bundle$datasets
  } else {
    stop("Could not find `stage1_analysis_datasets.rds` or datasets inside `stage1_backbone_bundle.rds`.", call. = FALSE)
  }
  
  list(
    backbone_bundle = backbone_bundle,
    analysis_datasets = analysis_datasets,
    dataset_registry = safe_read_csv(stage1_dataset_registry_file),
    scaling_registry = safe_read_csv(stage1_scaling_registry_file),
    metadata_registry = safe_read_csv(stage1_metadata_registry_file),
    formula_registry = safe_read_csv(stage1_formula_registry_file),
    horizon_registry = safe_read_csv(stage1_horizon_registry_file),
    threshold_registry = safe_read_csv(stage1_threshold_registry_file)
  )
}

## 🟠 Define: inherited backbone checks ===============================
validate_stage1_inputs <- function(stage1_inputs) {
  required_datasets <- c("PNU", "SNU", "merged")
  
  if (!is.list(stage1_inputs$analysis_datasets) || !all(required_datasets %in% names(stage1_inputs$analysis_datasets))) {
    stop("Stage 1 analysis datasets must be an R list containing PNU, SNU, and merged.", call. = FALSE)
  }
  
  required_cols <- c(
    "id", "site", "unique_person_id", "sex_num", "age_exact_entry", "age_s",
    "days_followup", "time_year", "status_num", "event_main", "censor_main",
    "right_censor_flag", "remission_flag"
  )
  
  for (dataset_name in required_datasets) {
    df <- stage1_inputs$analysis_datasets[[dataset_name]]
    
    if (!all(required_cols %in% names(df))) {
      missing_cols <- setdiff(required_cols, names(df))
      stop(sprintf("[%s] Missing required Stage 1 columns: %s", dataset_name, paste(missing_cols, collapse = ", ")), call. = FALSE)
    }
    
    if (nrow(df) == 0) {
      stop(sprintf("[%s] Stage 1 dataset contains zero rows.", dataset_name), call. = FALSE)
    }
    
    if (nrow(df) != dplyr::n_distinct(df$unique_person_id)) {
      stop(sprintf("[%s] `unique_person_id` must be unique.", dataset_name), call. = FALSE)
    }
    
    if (anyNA(df[, required_cols])) {
      stop(sprintf("[%s] Missing values detected in required Stage 1 columns.", dataset_name), call. = FALSE)
    }
  }
  
  merged_sites <- sort(unique(as.character(stage1_inputs$analysis_datasets$merged$site)))
  if (!identical(merged_sites, sort(c(pnu_site_label, snu_site_label)))) {
    stop(sprintf("Merged Stage 1 dataset must contain exactly sites `%s` and `%s`.", pnu_site_label, snu_site_label), call. = FALSE)
  }
  
  if (!all(c("dataset", "horizon_year", "interpretation_tier", "primary_supported_flag", "interpretation_note") %in% names(stage1_inputs$horizon_registry))) {
    stop("Stage 1 horizon registry is missing required columns.", call. = FALSE)
  }
  
  horizon_values <- sort(unique(as.integer(stage1_inputs$horizon_registry$horizon_year)))
  if (!identical(horizon_values, 1:10)) {
    stop("Stage 1 horizon registry must be locked to 1:10 years.", call. = FALSE)
  }
  
  threshold_values <- sort(unique(as.numeric(stage1_inputs$threshold_registry$threshold)))
  if (length(threshold_values) == 0 || anyNA(threshold_values)) {
    stop("Stage 1 threshold registry must contain at least one non-missing threshold.", call. = FALSE)
  }
  
  if (!all(c("dataset", "scaled_variable", "scaling_rule") %in% names(stage1_inputs$scaling_registry))) {
    stop("Stage 1 scaling registry is missing required columns.", call. = FALSE)
  }
  
  age_scaling_rows <- stage1_inputs$scaling_registry %>%
    filter(scaled_variable == "age_s", dataset %in% required_datasets)
  
  if (nrow(age_scaling_rows) < length(required_datasets)) {
    stop("Stage 1 scaling registry must contain `age_s` rows for PNU, SNU, and merged.", call. = FALSE)
  }
  
  if (!all(c("metadata_name", "metadata_value") %in% names(stage1_inputs$metadata_registry))) {
    stop("Stage 1 metadata registry is missing required columns.", call. = FALSE)
  }
  
  metadata_lookup <- setNames(stage1_inputs$metadata_registry$metadata_value, stage1_inputs$metadata_registry$metadata_name)
  
  if (!identical(unname(metadata_lookup[["analysis_time_variable"]]), "days_followup")) {
    stop("Stage 1 metadata registry must lock `analysis_time_variable` to `days_followup`.", call. = FALSE)
  }
  
  if (!identical(unname(metadata_lookup[["event_definition"]]), "status_num == 1")) {
    stop("Stage 1 metadata registry must lock the main event definition to `status_num == 1`.", call. = FALSE)
  }
  
  if (!identical(unname(metadata_lookup[["main_censoring_definition"]]), "status_num %in% c(0, 2)")) {
    stop("Stage 1 metadata registry must lock the main censoring definition to `status_num %in% c(0, 2)`.", call. = FALSE)
  }
  
  invisible(TRUE)
}

## 🟠 Define: support hierarchy helpers ===============================
normalize_interpretation_tier <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x == "primary-supported" ~ "primary_supported",
    TRUE ~ as.character(x)
  )
}

make_support_priority <- function(support_label) {
  dplyr::case_when(
    is.na(support_label) ~ NA_character_,
    support_label == "primary_supported" ~ "primary",
    support_label == "sensitivity" ~ "sensitivity",
    support_label == "secondary" ~ "secondary",
    support_label == "projection" ~ "projection",
    TRUE ~ "context"
  )
}

make_stage1_support_lookup <- function(horizon_registry) {
  horizon_registry %>%
    transmute(
      dataset = as.character(.data$dataset),
      horizon_year = as.integer(.data$horizon_year),
      interpretation_tier = as.character(.data$interpretation_tier),
      interpretation_note = as.character(.data$interpretation_note),
      primary_supported_flag = as.logical(.data$primary_supported_flag)
    ) %>%
    mutate(
      support_label = normalize_interpretation_tier(.data$interpretation_tier),
      support_priority = make_support_priority(.data$support_label)
    ) %>%
    arrange(match(.data$dataset, c("PNU", "SNU", "merged")), .data$horizon_year)
}

make_single_horizon_support <- function(support_lookup, dataset_value, horizon_value) {
  row <- support_lookup %>%
    filter(
      .data$dataset == .env$dataset_value,
      .data$horizon_year == as.integer(.env$horizon_value)
    )
  
  if (nrow(row) != 1L) {
    stop(
      sprintf(
        "Support lookup failed for dataset `%s` at horizon `%s`. Matched rows: %s",
        dataset_value,
        horizon_value,
        nrow(row)
      ),
      call. = FALSE
    )
  }
  
  row %>%
    transmute(
      support_label = .data$support_label,
      support_priority = .data$support_priority,
      primary_supported_flag = .data$primary_supported_flag
    )
}

make_conservative_observed_contrast_support <- function(horizon_year) {
  h0 <- as.integer(horizon_year)
  
  support_label <- dplyr::case_when(
    h0 == 1L ~ "primary_supported",
    h0 == 2L ~ "sensitivity",
    h0 >= 3L ~ "projection",
    TRUE ~ "projection"
  )
  
  tibble(
    support_label = support_label,
    support_priority = make_support_priority(support_label),
    primary_supported_flag = support_label == "primary_supported"
  )
}

make_common_restricted_window_flag <- function(horizon_year) {
  horizon_year %in% supported_short_horizons_year
}

make_support_basis <- function(source_key) {
  dplyr::case_when(
    source_key == "site_specific_observed_followup" ~ "site_specific_observed_followup",
    source_key == "separate_cohort_observed_support" ~ "separate_cohort_observed_support",
    source_key == "merged_pooled_observed_support" ~ "merged_pooled_observed_support",
    source_key == "merged_restricted_window_model_support" ~ "merged_restricted_window_model_support",
    source_key == "merged_piecewise_interval_support" ~ "merged_piecewise_interval_support",
    source_key == "observed_person_time_band_support" ~ "observed_person_time_band_support",
    TRUE ~ source_key
  )
}

## 🟠 Define: dataset preparation helpers ===============================
add_model_time_year <- function(df) {
  df %>%
    mutate(
      site = factor(as.character(site), levels = c(pnu_site_label, snu_site_label)),
      sex_num = as.integer(sex_num),
      event_main = as.integer(event_main),
      censor_main = as.integer(censor_main),
      right_censor_flag = as.integer(right_censor_flag),
      remission_flag = as.integer(remission_flag),
      model_time_year = pmax(as.numeric(time_year), time_zero_epsilon)
    )
}

summarize_dataset_analysis_totals <- function(df) {
  tibble(
    analysis_subject_n_total = nrow(df),
    analysis_transition_events_total = sum(df$event_main == 1L)
  )
}

summarize_horizon_window_counts <- function(df, horizon_year, require_full_observation = TRUE) {
  h0 <- as.numeric(horizon_year)
  max_followup_year <- max(df$model_time_year)
  full_observation_flag <- h0 <= max_followup_year
  site_chr <- as.character(df$site)
  
  if (require_full_observation && !full_observation_flag) {
    return(
      bind_cols(
        summarize_dataset_analysis_totals(df),
        tibble(
          at_risk_subject_n_total = NA_real_,
          transition_events_total = NA_real_,
          pnu_event_n = NA_real_,
          snu_event_n = NA_real_
        )
      )
    )
  }
  
  bind_cols(
    summarize_dataset_analysis_totals(df),
    tibble(
      at_risk_subject_n_total = sum(df$model_time_year >= h0),
      transition_events_total = sum(df$event_main == 1L & df$model_time_year <= h0),
      pnu_event_n = if (any(site_chr == pnu_site_label)) {
        sum(df$event_main == 1L & df$model_time_year <= h0 & site_chr == pnu_site_label)
      } else {
        NA_real_
      },
      snu_event_n = if (any(site_chr == snu_site_label)) {
        sum(df$event_main == 1L & df$model_time_year <= h0 & site_chr == snu_site_label)
      } else {
        NA_real_
      }
    )
  )
}

make_analysis_totals_lookup <- function(analysis_datasets) {
  bind_rows(lapply(names(analysis_datasets), function(dataset_name) {
    bind_cols(
      tibble(dataset = dataset_name),
      summarize_dataset_analysis_totals(analysis_datasets[[dataset_name]])
    )
  })) %>%
    arrange(match(dataset, c("PNU", "SNU", "merged")))
}

make_stage4_dataset_summary <- function(analysis_datasets, stage1_dataset_registry, support_lookup) {
  support_counts <- support_lookup %>%
    group_by(dataset) %>%
    summarise(
      n_primary_supported_horizons = sum(primary_supported_flag, na.rm = TRUE),
      .groups = "drop"
    )
  
  bind_rows(lapply(names(analysis_datasets), function(dataset_name) {
    df <- analysis_datasets[[dataset_name]]
    tibble(
      dataset = dataset_name,
      n_rows = nrow(df),
      n_unique_person_id = dplyr::n_distinct(df$unique_person_id),
      n_transition = sum(df$event_main),
      n_remission = sum(df$remission_flag),
      n_right_censoring = sum(df$right_censor_flag),
      n_main_censoring = sum(df$censor_main),
      person_time_years = sum(df$model_time_year),
      median_followup_years = stats::median(df$model_time_year),
      max_followup_years = max(df$model_time_year),
      stage4_role = dplyr::case_when(
        dataset_name == "merged" ~ "site_free_descriptive_plus_site_adjusted_contrast",
        TRUE ~ "separate_cohort_descriptive_and_contrast"
      )
    )
  })) %>%
    left_join(support_counts, by = "dataset") %>%
    left_join(
      stage1_dataset_registry %>%
        transmute(
          dataset = as.character(dataset),
          source_description = as.character(source_description),
          site_values = as.character(site_values)
        ),
      by = "dataset"
    ) %>%
    arrange(match(dataset, c("PNU", "SNU", "merged")))
}

## 🟠 Define: Kaplan-Meier risk helpers ===============================
extract_survfit_point <- function(fit, time_value, max_followup_year) {
  if (is.na(time_value) || time_value > max_followup_year) {
    return(tibble(
      survival = NA_real_,
      survival_conf_low = NA_real_,
      survival_conf_high = NA_real_,
      n_risk = NA_real_
    ))
  }
  
  ss <- summary(fit, times = time_value, extend = TRUE)
  
  survival_value <- if (length(ss$surv) == 0) 1 else as.numeric(ss$surv[1])
  lower_value <- if (length(ss$lower) == 0) survival_value else as.numeric(ss$lower[1])
  upper_value <- if (length(ss$upper) == 0) survival_value else as.numeric(ss$upper[1])
  n_risk_value <- if (length(ss$n.risk) == 0) sum(fit$n) else as.numeric(ss$n.risk[1])
  
  tibble(
    survival = survival_value,
    survival_conf_low = lower_value,
    survival_conf_high = upper_value,
    n_risk = n_risk_value
  )
}

compute_km_risk_value <- function(df, horizon_year) {
  fit <- survival::survfit(Surv(model_time_year, event_main) ~ 1, data = df)
  point <- extract_survfit_point(fit, horizon_year, max(df$model_time_year))
  if (is.na(point$survival)) {
    return(NA_real_)
  }
  1 - point$survival
}

compute_km_risk_trajectory <- function(df, dataset_name, support_lookup) {
  fit <- survival::survfit(Surv(model_time_year, event_main) ~ 1, data = df)
  max_followup_year <- max(df$model_time_year)
  
  purrr::map_dfr(1:10, function(h) {
    point <- extract_survfit_point(fit, h, max_followup_year)
    
    tibble(
      dataset = dataset_name,
      analysis_view = "site_free_observed",
      model_name = "KM",
      horizon_year = h,
      survival = point$survival,
      risk = ifelse(is.na(point$survival), NA_real_, 1 - point$survival),
      risk_conf_low = ifelse(is.na(point$survival_conf_high), NA_real_, 1 - point$survival_conf_high),
      risk_conf_high = ifelse(is.na(point$survival_conf_low), NA_real_, 1 - point$survival_conf_low),
      n_risk = point$n_risk,
      cumulative_transition_n = sum(df$event_main == 1L & df$model_time_year <= h),
      cumulative_main_censor_n = sum(df$censor_main == 1L & df$model_time_year <= h),
      cumulative_remission_n = sum(df$remission_flag == 1L & df$model_time_year <= h),
      max_followup_years = max_followup_year,
      available_within_observed_followup = h <= max_followup_year,
      estimable_flag = h <= max_followup_year,
      availability_note = dplyr::case_when(
        h <= max_followup_year ~ "within_observed_followup",
        TRUE ~ "beyond_max_observed_followup"
      )
    )
  }) %>%
    left_join(
      support_lookup %>%
        filter(dataset == dataset_name) %>%
        select(
          horizon_year,
          interpretation_tier,
          interpretation_note,
          support_label,
          support_priority,
          primary_supported_flag
        ),
      by = "horizon_year"
    ) %>%
    mutate(
      common_restricted_window_flag = make_common_restricted_window_flag(horizon_year),
      support_basis = dplyr::case_when(
        dataset == "merged" ~ make_support_basis("merged_pooled_observed_support"),
        TRUE ~ make_support_basis("site_specific_observed_followup")
      )
    )
}

build_restricted_dataset <- function(df, horizon_year) {
  df %>%
    mutate(
      window_year = horizon_year,
      time_restricted_year = pmin(model_time_year, horizon_year),
      event_restricted = as.integer(event_main == 1L & model_time_year <= horizon_year),
      site = factor(as.character(site), levels = c(pnu_site_label, snu_site_label))
    )
}

resample_rows <- function(df) {
  df[sample.int(nrow(df), size = nrow(df), replace = TRUE), , drop = FALSE]
}

resample_within_site <- function(df) {
  bind_rows(lapply(split(df, df$site), resample_rows))
}

compute_percentile_ci <- function(x, probs = c(0.025, 0.975)) {
  if (all(is.na(x))) {
    return(c(NA_real_, NA_real_))
  }
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, type = 6, names = FALSE))
}

## 🟠 Define: restricted Cox standardization helpers ===============================
fit_restricted_cox <- function(df, horizon_year, formula_rhs, model_name) {
  restricted_df <- tryCatch(
    build_restricted_dataset(df, horizon_year),
    error = function(e) NULL
  )
  
  if (is.null(restricted_df) || !is.data.frame(restricted_df) || nrow(restricted_df) == 0) {
    return(NULL)
  }
  
  fit_formula <- stats::as.formula(
    paste0("Surv(time_restricted_year, event_restricted) ~ ", formula_rhs)
  )
  
  fit_object <- suppressWarnings(
    tryCatch(
      survival::coxph(
        fit_formula,
        data = restricted_df,
        ties = "efron",
        x = TRUE,
        model = TRUE
      ),
      error = function(e) NULL
    )
  )
  
  if (is.null(fit_object)) {
    return(NULL)
  }
  
  coef_vec <- tryCatch(
    stats::coef(fit_object),
    error = function(e) NULL
  )
  
  if (is.null(coef_vec) || anyNA(coef_vec) || any(!is.finite(coef_vec))) {
    return(NULL)
  }
  
  list(
    fit = fit_object,
    data = restricted_df,
    horizon_year = horizon_year,
    model_name = model_name,
    formula_rhs = formula_rhs
  )
}

linear_predictor_from_cox <- function(fit, newdata) {
  design_terms <- stats::delete.response(stats::terms(fit))
  mm <- stats::model.matrix(design_terms, data = newdata)
  if ("(Intercept)" %in% colnames(mm)) {
    mm <- mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
  }
  coef_vec <- stats::coef(fit)
  mm <- mm[, names(coef_vec), drop = FALSE]
  as.numeric(mm %*% coef_vec)
}

extract_cumhaz_at_time <- function(fit, time_value) {
  bh <- suppressWarnings(survival::basehaz(fit, centered = FALSE))
  idx <- which(bh$time <= time_value)
  if (length(idx) == 0) {
    return(0)
  }
  as.numeric(bh$hazard[max(idx)])
}

make_empty_adjusted_supported_rows <- function(horizon_year, model_name, formula_rhs, n_boot, support_lookup) {
  h0 <- as.numeric(horizon_year)
  merged_support <- make_single_horizon_support(support_lookup, "merged", h0)
  
  tibble(
    dataset = c("merged", "merged", "PNU_vs_SNU"),
    analysis_view = "site_adjusted_standardized",
    source_block = "merged_restricted_cox_standardized",
    model_name = model_name,
    formula_rhs = formula_rhs,
    horizon_year = c(h0, h0, h0),
    window_year = c(h0, h0, h0),
    result_type = c("site_standardized_risk", "site_standardized_risk", "absolute_risk_difference"),
    site_counterfactual = c(pnu_site_label, snu_site_label, NA_character_),
    contrast_direction = c(NA_character_, NA_character_, "PNU_minus_SNU"),
    interval_label = NA_character_,
    estimate = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    n_risk = NA_real_,
    max_followup_years = NA_real_,
    available_within_observed_followup = NA,
    bootstrap_success = c(0L, 0L, 0L),
    bootstrap_iterations = c(n_boot, n_boot, n_boot),
    bootstrap_success_rate = rep(safe_divide(0L, n_boot), 3L),
    bootstrap_instability_flag = rep(ifelse(n_boot > 0, TRUE, NA), 3L),
    bootstrap_note = rep(ifelse(n_boot > 0, "bootstrap_fit_failed_all_replicates", NA_character_), 3L),
    estimable_flag = rep(FALSE, 3L),
    availability_note = rep("point_or_bootstrap_fit_failed", 3L),
    support_issue_reason = rep("point_or_bootstrap_fit_failed", 3L),
    analysis_subject_n_total = NA_real_,
    analysis_transition_events_total = NA_real_,
    at_risk_subject_n_total = NA_real_,
    transition_events_total = NA_real_,
    pnu_event_n = NA_real_,
    snu_event_n = NA_real_,
    support_label = rep(merged_support$support_label, 3L),
    support_priority = rep(merged_support$support_priority, 3L),
    primary_supported_flag = rep(merged_support$primary_supported_flag, 3L),
    common_restricted_window_flag = rep(make_common_restricted_window_flag(h0), 3L),
    support_basis = rep(make_support_basis("merged_restricted_window_model_support"), 3L),
    reported_estimate_suppressed_flag = rep(TRUE, 3L),
    output_label = "supported_short_horizon_contrast"
  )
}

compute_standardized_site_risks <- function(fit, reference_df, horizon_year, model_name, formula_rhs, support_lookup) {
  h0 <- as.numeric(horizon_year)
  merged_support <- make_single_horizon_support(support_lookup, "merged", h0)
  max_followup_years <- max(reference_df$model_time_year)
  n_risk_h0 <- sum(reference_df$model_time_year >= h0)
  support_counts <- summarize_horizon_window_counts(reference_df, h0, require_full_observation = TRUE)
  
  reference_df <- reference_df %>%
    mutate(site = factor(as.character(site), levels = c(pnu_site_label, snu_site_label)))
  
  H0_t <- tryCatch(
    extract_cumhaz_at_time(fit, h0),
    error = function(e) NA_real_
  )
  
  site_rows <- bind_rows(lapply(c(pnu_site_label, snu_site_label), function(site_value) {
    risk_estimate <- tryCatch({
      newdata_cf <- reference_df %>%
        mutate(site = factor(site_value, levels = c(pnu_site_label, snu_site_label)))
      
      lp <- linear_predictor_from_cox(fit, newdata_cf)
      
      if (is.na(H0_t) || !is.finite(H0_t) || anyNA(lp) || any(!is.finite(lp))) {
        NA_real_
      } else {
        mean(1 - exp(-H0_t * exp(lp)))
      }
    }, error = function(e) NA_real_)
    
    tibble(
      dataset = "merged",
      analysis_view = "site_adjusted_standardized",
      source_block = "merged_restricted_cox_standardized",
      model_name = model_name,
      formula_rhs = formula_rhs,
      horizon_year = h0,
      window_year = h0,
      result_type = "site_standardized_risk",
      site_counterfactual = site_value,
      contrast_direction = NA_character_,
      interval_label = NA_character_,
      estimate = risk_estimate,
      conf_low = NA_real_,
      conf_high = NA_real_,
      n_risk = n_risk_h0,
      max_followup_years = max_followup_years,
      available_within_observed_followup = h0 <= max_followup_years,
      bootstrap_success = NA_integer_,
      bootstrap_iterations = NA_integer_,
      bootstrap_success_rate = NA_real_,
      bootstrap_instability_flag = NA,
      bootstrap_note = NA_character_,
      estimable_flag = !is.na(risk_estimate),
      availability_note = dplyr::case_when(
        !is.na(risk_estimate) ~ "model_estimated",
        TRUE ~ "standardized_risk_not_estimable"
      ),
      support_issue_reason = dplyr::case_when(
        !is.na(risk_estimate) ~ NA_character_,
        TRUE ~ "standardized_risk_not_estimable"
      ),
      analysis_subject_n_total = support_counts$analysis_subject_n_total,
      analysis_transition_events_total = support_counts$analysis_transition_events_total,
      at_risk_subject_n_total = n_risk_h0,
      transition_events_total = support_counts$transition_events_total,
      pnu_event_n = support_counts$pnu_event_n,
      snu_event_n = support_counts$snu_event_n,
      support_label = merged_support$support_label,
      support_priority = merged_support$support_priority,
      primary_supported_flag = merged_support$primary_supported_flag,
      common_restricted_window_flag = make_common_restricted_window_flag(h0),
      support_basis = make_support_basis("merged_restricted_window_model_support"),
      reported_estimate_suppressed_flag = is.na(risk_estimate),
      output_label = "supported_short_horizon_contrast"
    )
  }))
  
  pnu_est <- site_rows %>%
    filter(site_counterfactual == pnu_site_label) %>%
    pull(estimate)
  pnu_est <- if (length(pnu_est) >= 1L) pnu_est[1] else NA_real_
  
  snu_est <- site_rows %>%
    filter(site_counterfactual == snu_site_label) %>%
    pull(estimate)
  snu_est <- if (length(snu_est) >= 1L) snu_est[1] else NA_real_
  
  contrast_row <- tibble(
    dataset = "PNU_vs_SNU",
    analysis_view = "site_adjusted_standardized",
    source_block = "merged_restricted_cox_standardized",
    model_name = model_name,
    formula_rhs = formula_rhs,
    horizon_year = h0,
    window_year = h0,
    result_type = "absolute_risk_difference",
    site_counterfactual = NA_character_,
    contrast_direction = "PNU_minus_SNU",
    interval_label = NA_character_,
    estimate = pnu_est - snu_est,
    conf_low = NA_real_,
    conf_high = NA_real_,
    n_risk = n_risk_h0,
    max_followup_years = max_followup_years,
    available_within_observed_followup = h0 <= max_followup_years,
    bootstrap_success = NA_integer_,
    bootstrap_iterations = NA_integer_,
    bootstrap_success_rate = NA_real_,
    bootstrap_instability_flag = NA,
    bootstrap_note = NA_character_,
    estimable_flag = !is.na(pnu_est) & !is.na(snu_est),
    availability_note = dplyr::case_when(
      !is.na(pnu_est) & !is.na(snu_est) ~ "model_estimated",
      TRUE ~ "contrast_not_estimable"
    ),
    support_issue_reason = dplyr::case_when(
      !is.na(pnu_est) & !is.na(snu_est) ~ NA_character_,
      TRUE ~ "contrast_not_estimable"
    ),
    analysis_subject_n_total = support_counts$analysis_subject_n_total,
    analysis_transition_events_total = support_counts$analysis_transition_events_total,
    at_risk_subject_n_total = n_risk_h0,
    transition_events_total = support_counts$transition_events_total,
    pnu_event_n = support_counts$pnu_event_n,
    snu_event_n = support_counts$snu_event_n,
    support_label = merged_support$support_label,
    support_priority = merged_support$support_priority,
    primary_supported_flag = merged_support$primary_supported_flag,
    common_restricted_window_flag = make_common_restricted_window_flag(h0),
    support_basis = make_support_basis("merged_restricted_window_model_support"),
    reported_estimate_suppressed_flag = is.na(pnu_est) | is.na(snu_est),
    output_label = "supported_short_horizon_contrast"
  )
  
  bind_rows(site_rows, contrast_row)
}

extract_supported_triplet <- function(rows) {
  out <- c(PNU = NA_real_, SNU = NA_real_, diff = NA_real_)
  
  if (is.null(rows) || !is.data.frame(rows) || nrow(rows) == 0) {
    return(out)
  }
  
  pnu_val <- rows %>%
    filter(
      result_type == "site_standardized_risk",
      as.character(site_counterfactual) == pnu_site_label
    ) %>%
    pull(estimate)
  
  snu_val <- rows %>%
    filter(
      result_type == "site_standardized_risk",
      as.character(site_counterfactual) == snu_site_label
    ) %>%
    pull(estimate)
  
  diff_val <- rows %>%
    filter(result_type == "absolute_risk_difference") %>%
    pull(estimate)
  
  if (length(pnu_val) >= 1L) out["PNU"] <- pnu_val[1]
  if (length(snu_val) >= 1L) out["SNU"] <- snu_val[1]
  if (length(diff_val) >= 1L) out["diff"] <- diff_val[1]
  
  out
}

bootstrap_adjusted_supported <- function(df, horizon_year, formula_rhs, model_name, n_boot, seed_value, support_lookup) {
  h0 <- as.numeric(horizon_year)
  empty_rows <- make_empty_adjusted_supported_rows(h0, model_name, formula_rhs, n_boot, support_lookup)
  
  reference_df <- tryCatch(
    build_restricted_dataset(df, h0),
    error = function(e) NULL
  )
  
  if (is.null(reference_df) || !is.data.frame(reference_df) || nrow(reference_df) == 0) {
    return(empty_rows)
  }
  
  reference_counts <- summarize_horizon_window_counts(reference_df, h0, require_full_observation = TRUE)
  reference_at_risk_n <- sum(reference_df$model_time_year >= h0)
  reference_max_followup <- max(reference_df$model_time_year)
  reference_available <- h0 <= reference_max_followup
  
  point_fit <- fit_restricted_cox(df, h0, formula_rhs, model_name)
  
  if (is.null(point_fit)) {
    return(
      empty_rows %>%
        mutate(
          n_risk = reference_at_risk_n,
          max_followup_years = reference_max_followup,
          available_within_observed_followup = reference_available,
          analysis_subject_n_total = reference_counts$analysis_subject_n_total,
          analysis_transition_events_total = reference_counts$analysis_transition_events_total,
          at_risk_subject_n_total = reference_at_risk_n,
          transition_events_total = reference_counts$transition_events_total,
          pnu_event_n = reference_counts$pnu_event_n,
          snu_event_n = reference_counts$snu_event_n
        )
    )
  }
  
  point_rows <- tryCatch(
    compute_standardized_site_risks(point_fit$fit, reference_df, h0, model_name, formula_rhs, support_lookup),
    error = function(e) NULL
  )
  
  if (is.null(point_rows) || nrow(point_rows) != 3L || any(!point_rows$estimable_flag)) {
    return(
      empty_rows %>%
        mutate(
          n_risk = reference_at_risk_n,
          max_followup_years = reference_max_followup,
          available_within_observed_followup = reference_available,
          analysis_subject_n_total = reference_counts$analysis_subject_n_total,
          analysis_transition_events_total = reference_counts$analysis_transition_events_total,
          at_risk_subject_n_total = reference_at_risk_n,
          transition_events_total = reference_counts$transition_events_total,
          pnu_event_n = reference_counts$pnu_event_n,
          snu_event_n = reference_counts$snu_event_n
        )
    )
  }
  
  if (n_boot <= 0) {
    return(
      point_rows %>%
        mutate(
          bootstrap_success = NA_integer_,
          bootstrap_iterations = n_boot,
          bootstrap_success_rate = NA_real_,
          bootstrap_instability_flag = NA,
          bootstrap_note = NA_character_
        )
    )
  }
  
  set.seed(seed_value)
  
  boot_matrix <- vapply(
    X = seq_len(n_boot),
    FUN = function(b) {
      boot_df <- resample_within_site(df)
      
      boot_fit <- tryCatch(
        fit_restricted_cox(boot_df, h0, formula_rhs, model_name),
        error = function(e) NULL
      )
      
      if (is.null(boot_fit)) {
        return(c(PNU = NA_real_, SNU = NA_real_, diff = NA_real_))
      }
      
      boot_rows <- tryCatch(
        compute_standardized_site_risks(boot_fit$fit, reference_df, h0, model_name, formula_rhs, support_lookup),
        error = function(e) NULL
      )
      
      extract_supported_triplet(boot_rows)
    },
    FUN.VALUE = c(PNU = NA_real_, SNU = NA_real_, diff = NA_real_)
  )
  
  ci_pnu <- compute_percentile_ci(boot_matrix["PNU", ])
  ci_snu <- compute_percentile_ci(boot_matrix["SNU", ])
  ci_diff <- compute_percentile_ci(boot_matrix["diff", ])
  
  success_pnu <- sum(!is.na(boot_matrix["PNU", ]))
  success_snu <- sum(!is.na(boot_matrix["SNU", ]))
  success_diff <- sum(!is.na(boot_matrix["diff", ]))
  
  point_rows %>%
    mutate(
      conf_low = case_when(
        result_type == "site_standardized_risk" & site_counterfactual == pnu_site_label ~ ci_pnu[1],
        result_type == "site_standardized_risk" & site_counterfactual == snu_site_label ~ ci_snu[1],
        result_type == "absolute_risk_difference" ~ ci_diff[1],
        TRUE ~ conf_low
      ),
      conf_high = case_when(
        result_type == "site_standardized_risk" & site_counterfactual == pnu_site_label ~ ci_pnu[2],
        result_type == "site_standardized_risk" & site_counterfactual == snu_site_label ~ ci_snu[2],
        result_type == "absolute_risk_difference" ~ ci_diff[2],
        TRUE ~ conf_high
      ),
      bootstrap_success = case_when(
        result_type == "site_standardized_risk" & site_counterfactual == pnu_site_label ~ success_pnu,
        result_type == "site_standardized_risk" & site_counterfactual == snu_site_label ~ success_snu,
        result_type == "absolute_risk_difference" ~ success_diff,
        TRUE ~ NA_integer_
      ),
      bootstrap_iterations = n_boot,
      bootstrap_success_rate = case_when(
        result_type == "site_standardized_risk" & site_counterfactual == pnu_site_label ~ safe_divide(success_pnu, n_boot),
        result_type == "site_standardized_risk" & site_counterfactual == snu_site_label ~ safe_divide(success_snu, n_boot),
        result_type == "absolute_risk_difference" ~ safe_divide(success_diff, n_boot),
        TRUE ~ NA_real_
      ),
      bootstrap_instability_flag = case_when(
        !is.na(bootstrap_success_rate) ~ bootstrap_success_rate < bootstrap_min_success_rate,
        TRUE ~ NA
      ),
      bootstrap_note = case_when(
        is.na(bootstrap_success_rate) ~ NA_character_,
        bootstrap_success_rate < bootstrap_min_success_rate ~ sprintf("bootstrap_success_below_%.2f", bootstrap_min_success_rate),
        TRUE ~ "bootstrap_success_acceptable"
      )
    )
}

bootstrap_observed_supported <- function(pnu_df, snu_df, horizon_year, n_boot, seed_value, support_lookup) {
  h0 <- as.numeric(horizon_year)
  
  fit_pnu <- survival::survfit(Surv(model_time_year, event_main) ~ 1, data = pnu_df)
  fit_snu <- survival::survfit(Surv(model_time_year, event_main) ~ 1, data = snu_df)
  
  max_followup_pnu <- max(pnu_df$model_time_year)
  max_followup_snu <- max(snu_df$model_time_year)
  
  point_pnu_detail <- extract_survfit_point(fit_pnu, h0, max_followup_pnu)
  point_snu_detail <- extract_survfit_point(fit_snu, h0, max_followup_snu)
  
  point_pnu <- ifelse(is.na(point_pnu_detail$survival), NA_real_, 1 - point_pnu_detail$survival)
  point_snu <- ifelse(is.na(point_snu_detail$survival), NA_real_, 1 - point_snu_detail$survival)
  point_diff <- point_pnu - point_snu
  
  pnu_support_counts <- summarize_horizon_window_counts(pnu_df, h0, require_full_observation = TRUE)
  snu_support_counts <- summarize_horizon_window_counts(snu_df, h0, require_full_observation = TRUE)
  
  pnu_support <- make_single_horizon_support(support_lookup, "PNU", h0)
  snu_support <- make_single_horizon_support(support_lookup, "SNU", h0)
  contrast_support <- make_conservative_observed_contrast_support(h0)
  
  support_label_vec <- c(
    pnu_support$support_label,
    snu_support$support_label,
    contrast_support$support_label
  )
  support_priority_vec <- c(
    pnu_support$support_priority,
    snu_support$support_priority,
    contrast_support$support_priority
  )
  primary_flag_vec <- c(
    pnu_support$primary_supported_flag,
    snu_support$primary_supported_flag,
    contrast_support$primary_supported_flag
  )
  common_window_vec <- rep(make_common_restricted_window_flag(h0), 3L)
  max_followup_vec <- c(max_followup_pnu, max_followup_snu, min(max_followup_pnu, max_followup_snu))
  available_vec <- c(h0 <= max_followup_pnu, h0 <= max_followup_snu, h0 <= min(max_followup_pnu, max_followup_snu))
  n_risk_vec <- c(point_pnu_detail$n_risk, point_snu_detail$n_risk, ifelse(anyNA(c(point_pnu_detail$n_risk, point_snu_detail$n_risk)), NA_real_, point_pnu_detail$n_risk + point_snu_detail$n_risk))
  support_basis_vec <- c(
    make_support_basis("site_specific_observed_followup"),
    make_support_basis("site_specific_observed_followup"),
    make_support_basis("separate_cohort_observed_support")
  )
  
  point_rows <- tibble(
    dataset = c("PNU", "SNU", "PNU_vs_SNU"),
    analysis_view = "site_free_observed",
    source_block = "observed_km_restricted",
    model_name = "KM",
    formula_rhs = NA_character_,
    horizon_year = c(h0, h0, h0),
    window_year = c(h0, h0, h0),
    result_type = c("observed_risk", "observed_risk", "absolute_risk_difference"),
    site_counterfactual = c(pnu_site_label, snu_site_label, NA_character_),
    contrast_direction = c(NA_character_, NA_character_, "PNU_minus_SNU"),
    interval_label = NA_character_,
    estimate = c(point_pnu, point_snu, point_diff),
    conf_low = NA_real_,
    conf_high = NA_real_,
    n_risk = n_risk_vec,
    max_followup_years = max_followup_vec,
    available_within_observed_followup = available_vec,
    bootstrap_success = NA_integer_,
    bootstrap_iterations = NA_integer_,
    bootstrap_success_rate = NA_real_,
    bootstrap_instability_flag = NA,
    bootstrap_note = NA_character_,
    estimable_flag = !is.na(c(point_pnu, point_snu, point_diff)),
    availability_note = case_when(
      !is.na(c(point_pnu, point_snu, point_diff)) ~ "within_observed_followup",
      TRUE ~ "not_estimable"
    ),
    support_issue_reason = case_when(
      !is.na(c(point_pnu, point_snu, point_diff)) ~ NA_character_,
      TRUE ~ "not_estimable"
    ),
    analysis_subject_n_total = c(
      pnu_support_counts$analysis_subject_n_total,
      snu_support_counts$analysis_subject_n_total,
      pnu_support_counts$analysis_subject_n_total + snu_support_counts$analysis_subject_n_total
    ),
    analysis_transition_events_total = c(
      pnu_support_counts$analysis_transition_events_total,
      snu_support_counts$analysis_transition_events_total,
      pnu_support_counts$analysis_transition_events_total + snu_support_counts$analysis_transition_events_total
    ),
    at_risk_subject_n_total = n_risk_vec,
    transition_events_total = c(
      pnu_support_counts$transition_events_total,
      snu_support_counts$transition_events_total,
      pnu_support_counts$transition_events_total + snu_support_counts$transition_events_total
    ),
    pnu_event_n = c(
      pnu_support_counts$transition_events_total,
      NA_real_,
      pnu_support_counts$transition_events_total
    ),
    snu_event_n = c(
      NA_real_,
      snu_support_counts$transition_events_total,
      snu_support_counts$transition_events_total
    ),
    support_label = support_label_vec,
    support_priority = support_priority_vec,
    primary_supported_flag = primary_flag_vec,
    common_restricted_window_flag = common_window_vec,
    support_basis = support_basis_vec,
    reported_estimate_suppressed_flag = is.na(c(point_pnu, point_snu, point_diff)),
    output_label = "supported_short_horizon_contrast"
  )
  
  if (n_boot <= 0) {
    return(
      point_rows %>%
        mutate(
          bootstrap_success = NA_integer_,
          bootstrap_iterations = n_boot,
          bootstrap_success_rate = NA_real_,
          bootstrap_instability_flag = NA,
          bootstrap_note = NA_character_
        )
    )
  }
  
  set.seed(seed_value)
  
  boot_matrix <- vapply(
    X = seq_len(n_boot),
    FUN = function(b) {
      boot_pnu <- resample_rows(pnu_df)
      boot_snu <- resample_rows(snu_df)
      
      pnu_risk <- tryCatch(
        compute_km_risk_value(boot_pnu, h0),
        error = function(e) NA_real_
      )
      snu_risk <- tryCatch(
        compute_km_risk_value(boot_snu, h0),
        error = function(e) NA_real_
      )
      
      c(PNU = pnu_risk, SNU = snu_risk, diff = pnu_risk - snu_risk)
    },
    FUN.VALUE = c(PNU = NA_real_, SNU = NA_real_, diff = NA_real_)
  )
  
  ci_pnu <- compute_percentile_ci(boot_matrix["PNU", ])
  ci_snu <- compute_percentile_ci(boot_matrix["SNU", ])
  ci_diff <- compute_percentile_ci(boot_matrix["diff", ])
  
  success_pnu <- sum(!is.na(boot_matrix["PNU", ]))
  success_snu <- sum(!is.na(boot_matrix["SNU", ]))
  success_diff <- sum(!is.na(boot_matrix["diff", ]))
  
  point_rows %>%
    mutate(
      conf_low = case_when(
        dataset == "PNU" ~ ci_pnu[1],
        dataset == "SNU" ~ ci_snu[1],
        dataset == "PNU_vs_SNU" ~ ci_diff[1],
        TRUE ~ conf_low
      ),
      conf_high = case_when(
        dataset == "PNU" ~ ci_pnu[2],
        dataset == "SNU" ~ ci_snu[2],
        dataset == "PNU_vs_SNU" ~ ci_diff[2],
        TRUE ~ conf_high
      ),
      bootstrap_success = case_when(
        dataset == "PNU" ~ success_pnu,
        dataset == "SNU" ~ success_snu,
        dataset == "PNU_vs_SNU" ~ success_diff,
        TRUE ~ NA_integer_
      ),
      bootstrap_iterations = n_boot,
      bootstrap_success_rate = case_when(
        dataset == "PNU" ~ safe_divide(success_pnu, n_boot),
        dataset == "SNU" ~ safe_divide(success_snu, n_boot),
        dataset == "PNU_vs_SNU" ~ safe_divide(success_diff, n_boot),
        TRUE ~ NA_real_
      ),
      bootstrap_instability_flag = case_when(
        !is.na(bootstrap_success_rate) ~ bootstrap_success_rate < bootstrap_min_success_rate,
        TRUE ~ NA
      ),
      bootstrap_note = case_when(
        is.na(bootstrap_success_rate) ~ NA_character_,
        bootstrap_success_rate < bootstrap_min_success_rate ~ sprintf("bootstrap_success_below_%.2f", bootstrap_min_success_rate),
        TRUE ~ "bootstrap_success_acceptable"
      )
    )
}

## 🟠 Define: piecewise site-effect helpers ===============================
make_piecewise_dataset <- function(df) {
  base_df <- df %>%
    transmute(
      unique_person_id,
      site = factor(as.character(site), levels = c(pnu_site_label, snu_site_label)),
      age_s,
      sex_num = as.integer(sex_num),
      model_time_year,
      event_main = as.integer(event_main)
    ) %>%
    as.data.frame()
  
  split_df <- survival::survSplit(
    Surv(model_time_year, event_main) ~ .,
    data = base_df,
    cut = piecewise_cuts_year,
    start = "tstart_year",
    end = "tstop_year",
    event = "event_piece",
    episode = "interval_id"
  )
  
  split_df %>%
    mutate(
      interval_label = factor(
        dplyr::case_when(
          interval_id == 1L ~ "0-1",
          interval_id == 2L ~ "1-2",
          TRUE ~ ">2"
        ),
        levels = c("0-1", "1-2", ">2")
      ),
      site_snu = as.integer(as.character(site) == snu_site_label),
      site_snu_int_0_1 = as.integer(interval_label == "0-1") * site_snu,
      site_snu_int_1_2 = as.integer(interval_label == "1-2") * site_snu,
      site_snu_int_gt_2 = as.integer(interval_label == ">2") * site_snu
    )
}

summarize_piecewise_interval_support <- function(split_df, min_events_per_site, min_subjects_per_site) {
  site_support <- split_df %>%
    mutate(
      interval_label = as.character(interval_label),
      site = as.character(site)
    ) %>%
    group_by(interval_label, site) %>%
    summarise(
      site_subject_n = dplyr::n_distinct(unique_person_id),
      site_event_n = sum(event_piece == 1L),
      site_person_years = sum(pmax(tstop_year - tstart_year, 0)),
      .groups = "drop"
    )
  
  site_template <- tidyr::expand_grid(
    interval_label = c("0-1", "1-2", ">2"),
    site = c(pnu_site_label, snu_site_label)
  )
  
  site_support_full <- site_template %>%
    left_join(site_support, by = c("interval_label", "site")) %>%
    mutate(
      site_subject_n = dplyr::coalesce(site_subject_n, 0L),
      site_event_n = dplyr::coalesce(site_event_n, 0L),
      site_person_years = dplyr::coalesce(site_person_years, 0)
    )
  
  interval_support <- site_support_full %>%
    pivot_wider(
      names_from = site,
      values_from = c(site_subject_n, site_event_n, site_person_years),
      names_sep = "__"
    ) %>%
    mutate(
      interval_label = factor(interval_label, levels = c("0-1", "1-2", ">2")),
      interval_start_year = case_when(
        interval_label == "0-1" ~ 0,
        interval_label == "1-2" ~ 1,
        TRUE ~ 2
      ),
      interval_end_year = case_when(
        interval_label == "0-1" ~ 1,
        interval_label == "1-2" ~ 2,
        TRUE ~ NA_real_
      ),
      pnu_subject_n = dplyr::coalesce(.data[[paste0("site_subject_n__", pnu_site_label)]], 0L),
      snu_subject_n = dplyr::coalesce(.data[[paste0("site_subject_n__", snu_site_label)]], 0L),
      pnu_event_n = dplyr::coalesce(.data[[paste0("site_event_n__", pnu_site_label)]], 0L),
      snu_event_n = dplyr::coalesce(.data[[paste0("site_event_n__", snu_site_label)]], 0L),
      pnu_person_years = dplyr::coalesce(.data[[paste0("site_person_years__", pnu_site_label)]], 0),
      snu_person_years = dplyr::coalesce(.data[[paste0("site_person_years__", snu_site_label)]], 0),
      at_risk_subject_n_total = pnu_subject_n + snu_subject_n,
      transition_events_total = pnu_event_n + snu_event_n,
      person_years_total = pnu_person_years + snu_person_years,
      support_ok = (
        pnu_event_n >= min_events_per_site &
          snu_event_n >= min_events_per_site &
          pnu_subject_n >= min_subjects_per_site &
          snu_subject_n >= min_subjects_per_site &
          pnu_person_years > 0 &
          snu_person_years > 0
      ),
      support_issue_reason = purrr::pmap_chr(
        list(pnu_event_n, snu_event_n, pnu_subject_n, snu_subject_n, pnu_person_years, snu_person_years),
        function(pnu_e, snu_e, pnu_n, snu_n, pnu_py, snu_py) {
          tokens <- character(0)
          if (pnu_e < min_events_per_site || snu_e < min_events_per_site) {
            tokens <- c(tokens, "zero_or_sparse_site_events")
          }
          if (pnu_n < min_subjects_per_site || snu_n < min_subjects_per_site) {
            tokens <- c(tokens, sprintf("site_subjects_below_%s", min_subjects_per_site))
          }
          if (pnu_py <= 0 || snu_py <= 0) {
            tokens <- c(tokens, "no_person_time_in_one_site")
          }
          combine_reason_tokens(tokens)
        }
      ),
      support_label = case_when(
        interval_label %in% c("0-1", "1-2") ~ "primary_supported",
        TRUE ~ "projection"
      ),
      support_priority = case_when(
        interval_label %in% c("0-1", "1-2") ~ "primary",
        TRUE ~ "projection"
      ),
      primary_supported_flag = interval_label %in% c("0-1", "1-2"),
      common_restricted_window_flag = interval_label %in% c("0-1", "1-2"),
      output_label = case_when(
        interval_label %in% c("0-1", "1-2") ~ "supported_short_horizon_contrast",
        TRUE ~ "late_horizon_tail_contrast"
      ),
      dataset = "merged",
      analysis_view = "site_adjusted",
      source_block = "merged_piecewise_interval_support",
      support_basis = make_support_basis("merged_piecewise_interval_support")
    ) %>%
    transmute(
      dataset, analysis_view, source_block,
      interval_label, interval_start_year, interval_end_year,
      pnu_subject_n, snu_subject_n, pnu_event_n, snu_event_n,
      pnu_person_years, snu_person_years,
      at_risk_subject_n_total, transition_events_total, person_years_total,
      estimable_flag = support_ok,
      support_issue_reason,
      support_label, support_priority, primary_supported_flag,
      common_restricted_window_flag, support_basis, output_label
    ) %>%
    arrange(match(interval_label, c("0-1", "1-2", ">2")))
  
  interval_support
}

fit_piecewise_site_model_bundle <- function(split_df, model_name, covariate_rhs, interval_support_df) {
  rhs_terms <- character(0)
  if (!is.null(covariate_rhs) && nzchar(trimws(covariate_rhs))) {
    rhs_terms <- c(rhs_terms, covariate_rhs)
  }
  rhs_terms <- c(
    rhs_terms,
    "site_snu_int_0_1",
    "site_snu_int_1_2",
    "site_snu_int_gt_2",
    "strata(interval_label)",
    "cluster(unique_person_id)"
  )
  
  fit_formula <- stats::as.formula(
    paste("Surv(tstart_year, tstop_year, event_piece) ~", paste(rhs_terms, collapse = " + "))
  )
  
  fit_object <- suppressWarnings(
    tryCatch(
      survival::coxph(fit_formula, data = split_df, ties = "efron"),
      error = function(e) NULL
    )
  )
  
  if (is.null(fit_object)) {
    return(NULL)
  }
  
  coef_table <- as.data.frame(summary(fit_object)$coefficients)
  coef_table$term <- rownames(coef_table)
  robust_col <- if ("robust se" %in% colnames(coef_table)) "robust se" else "se(coef)"
  p_col <- if ("Pr(>|z|)" %in% colnames(coef_table)) "Pr(>|z|)" else NA_character_
  
  term_template <- tibble(
    term = c("site_snu_int_0_1", "site_snu_int_1_2", "site_snu_int_gt_2"),
    interval_label = c("0-1", "1-2", ">2"),
    interval_start_year = c(0, 1, 2),
    interval_end_year = c(1, 2, NA_real_)
  )
  
  interval_support_join <- interval_support_df %>%
    select(
      interval_label,
      at_risk_subject_n_total, transition_events_total, person_years_total,
      pnu_subject_n, snu_subject_n, pnu_event_n, snu_event_n,
      pnu_person_years, snu_person_years,
      estimable_flag, support_issue_reason,
      support_label, support_priority, primary_supported_flag,
      common_restricted_window_flag, support_basis, output_label
    )
  
  summary_table <- term_template %>%
    left_join(
      coef_table %>%
        transmute(
          term = term,
          raw_estimate_log_hr = coef,
          raw_robust_se = .data[[robust_col]],
          raw_p_value = if (is.character(p_col)) .data[[p_col]] else NA_real_
        ),
      by = "term"
    ) %>%
    left_join(interval_support_join, by = "interval_label") %>%
    mutate(
      raw_estimate_hr = exp(raw_estimate_log_hr),
      raw_conf_low_hr = exp(raw_estimate_log_hr - 1.96 * raw_robust_se),
      raw_conf_high_hr = exp(raw_estimate_log_hr + 1.96 * raw_robust_se),
      fit_issue_reason = case_when(
        is.na(raw_estimate_log_hr) | !is.finite(raw_estimate_log_hr) | is.na(raw_robust_se) | !is.finite(raw_robust_se) ~ "coefficient_not_estimable_from_fit",
        !is.finite(raw_estimate_hr) | !is.finite(raw_conf_low_hr) | !is.finite(raw_conf_high_hr) ~ "coefficient_scale_unstable",
        TRUE ~ NA_character_
      ),
      support_issue_reason = case_when(
        !estimable_flag & !is.na(support_issue_reason) & !is.na(fit_issue_reason) ~ paste(support_issue_reason, fit_issue_reason, sep = ";"),
        !estimable_flag ~ support_issue_reason,
        estimable_flag & !is.na(fit_issue_reason) ~ fit_issue_reason,
        TRUE ~ NA_character_
      ),
      estimable_flag = dplyr::coalesce(estimable_flag, FALSE) &
        is.finite(raw_estimate_log_hr) &
        is.finite(raw_robust_se) &
        is.finite(raw_estimate_hr) &
        is.finite(raw_conf_low_hr) &
        is.finite(raw_conf_high_hr),
      estimate_log_hr = ifelse(estimable_flag, raw_estimate_log_hr, NA_real_),
      robust_se = ifelse(estimable_flag, raw_robust_se, NA_real_),
      estimate_hr = ifelse(estimable_flag, raw_estimate_hr, NA_real_),
      conf_low_hr = ifelse(estimable_flag, raw_conf_low_hr, NA_real_),
      conf_high_hr = ifelse(estimable_flag, raw_conf_high_hr, NA_real_),
      p_value = ifelse(estimable_flag, raw_p_value, NA_real_),
      reported_estimate_suppressed_flag = !estimable_flag,
      reporting_rule = case_when(
        estimable_flag ~ "reported",
        TRUE ~ "suppressed_for_inadequate_interval_support_or_fit_instability"
      ),
      raw_estimate_log_hr = ifelse(estimable_flag, raw_estimate_log_hr, NA_real_),
      raw_robust_se = ifelse(estimable_flag, raw_robust_se, NA_real_),
      raw_estimate_hr = ifelse(estimable_flag, raw_estimate_hr, NA_real_),
      raw_conf_low_hr = ifelse(estimable_flag, raw_conf_low_hr, NA_real_),
      raw_conf_high_hr = ifelse(estimable_flag, raw_conf_high_hr, NA_real_),
      raw_p_value = ifelse(estimable_flag, raw_p_value, NA_real_),
      dataset = "merged",
      analysis_view = dplyr::case_when(
        model_name == "piecewise_unadjusted" ~ "site_comparison_unadjusted",
        TRUE ~ "site_adjusted"
      ),
      source_block = "merged_piecewise_cox",
      model_name = model_name,
      formula_rhs = ifelse(nzchar(trimws(covariate_rhs)), covariate_rhs, NA_character_),
      contrast_direction = "SNU_vs_PNU",
      availability_note = case_when(
        estimable_flag ~ "interval_support_adequate",
        TRUE ~ coalesce_character(support_issue_reason, "interval_support_inadequate")
      ),
      n_subjects = dplyr::n_distinct(split_df$unique_person_id),
      n_events_total = sum(split_df$event_piece),
      max_followup_years = max(split_df$tstop_year),
      available_within_observed_followup = TRUE
    ) %>%
    select(
      dataset, analysis_view, source_block, model_name, formula_rhs,
      interval_label, interval_start_year, interval_end_year,
      contrast_direction,
      at_risk_subject_n_total, transition_events_total, person_years_total,
      pnu_subject_n, snu_subject_n, pnu_event_n, snu_event_n,
      pnu_person_years, snu_person_years,
      raw_estimate_log_hr, raw_robust_se, raw_estimate_hr, raw_conf_low_hr, raw_conf_high_hr, raw_p_value,
      estimate_log_hr, robust_se, estimate_hr, conf_low_hr, conf_high_hr, p_value,
      estimable_flag, reported_estimate_suppressed_flag, reporting_rule, support_issue_reason, availability_note,
      n_subjects, n_events_total, max_followup_years, available_within_observed_followup,
      support_label, support_priority, primary_supported_flag, common_restricted_window_flag, support_basis, output_label
    )
  
  list(
    fit = fit_object,
    summary_table = summary_table,
    model_name = model_name,
    formula_rhs = ifelse(nzchar(trimws(covariate_rhs)), covariate_rhs, NA_character_)
  )
}

## 🟠 Define: hazard-pattern helpers ===============================
compute_survival_at_time <- function(fit, time_value, max_followup_year) {
  if (time_value <= 0) {
    return(1)
  }
  point <- extract_survfit_point(fit, time_value, max_followup_year)
  point$survival
}

summarize_hazard_bands <- function(df, dataset_name, band_breaks, support_lookup) {
  fit <- survival::survfit(Surv(model_time_year, event_main) ~ 1, data = df)
  max_followup_year <- max(df$model_time_year)
  
  bands <- tibble(
    dataset = dataset_name,
    band_start_year = band_breaks[-length(band_breaks)],
    band_end_year = band_breaks[-1]
  ) %>%
    mutate(
      band_label = paste0(band_start_year, "-", band_end_year),
      support_ref_horizon_year = as.integer(band_end_year),
      output_label = if_else(band_end_year <= 2, "supported_short_horizon_contrast", "late_horizon_tail_contrast")
    )
  
  bands %>%
    rowwise() %>%
    mutate(
      observed_band_end_year = pmin(band_end_year, max_followup_year),
      band_full_width_year = band_end_year - band_start_year,
      observed_band_width_year = pmax(observed_band_end_year - band_start_year, 0),
      observed_band_fraction = safe_divide(observed_band_width_year, band_full_width_year),
      full_band_observed_flag = band_end_year <= max_followup_year,
      partial_band_flag = band_start_year < max_followup_year & band_end_year > max_followup_year,
      no_observed_time_flag = band_start_year >= max_followup_year,
      at_risk_start_n = sum(df$model_time_year >= band_start_year),
      transition_events = sum(df$event_main == 1L & df$model_time_year > band_start_year & df$model_time_year <= observed_band_end_year),
      remission_censor_n = sum(df$remission_flag == 1L & df$model_time_year > band_start_year & df$model_time_year <= observed_band_end_year),
      right_censor_n = sum(df$right_censor_flag == 1L & df$model_time_year > band_start_year & df$model_time_year <= observed_band_end_year),
      main_censor_n = sum(df$censor_main == 1L & df$model_time_year > band_start_year & df$model_time_year <= observed_band_end_year),
      person_years = sum(pmax(pmin(df$model_time_year, observed_band_end_year) - band_start_year, 0)),
      crude_hazard_per_100py = ifelse(person_years > 0, 100 * transition_events / person_years, NA_real_),
      survival_start = compute_survival_at_time(fit, band_start_year, max_followup_year),
      survival_observed_end = compute_survival_at_time(fit, observed_band_end_year, max_followup_year),
      survival_end = compute_survival_at_time(fit, band_end_year, max_followup_year),
      conditional_interval_risk = dplyr::if_else(
        full_band_observed_flag & !is.na(survival_start) & !is.na(survival_end) & survival_start > 0,
        1 - (survival_end / survival_start),
        NA_real_
      ),
      conditional_interval_risk_observed_portion = dplyr::if_else(
        observed_band_width_year > 0 & !is.na(survival_start) & !is.na(survival_observed_end) & survival_start > 0,
        1 - (survival_observed_end / survival_start),
        NA_real_
      ),
      max_followup_years = max_followup_year,
      available_within_observed_followup = band_end_year <= max_followup_year,
      estimable_flag = person_years > 0,
      availability_note = dplyr::case_when(
        no_observed_time_flag ~ "no_observed_time_in_interval",
        full_band_observed_flag ~ "full_band_observed",
        partial_band_flag ~ "partial_band_observed",
        person_years > 0 ~ "interval_person_time_available",
        TRUE ~ "no_person_time_in_interval"
      ),
      analysis_view = "site_free_observed",
      source_block = "interval_hazard_summary",
      support_basis = make_support_basis("observed_person_time_band_support")
    ) %>%
    ungroup() %>%
    left_join(
      support_lookup %>%
        filter(dataset == dataset_name) %>%
        select(horizon_year, support_label, support_priority, primary_supported_flag),
      by = c("support_ref_horizon_year" = "horizon_year")
    ) %>%
    mutate(
      common_restricted_window_flag = band_end_year <= 2
    ) %>%
    select(
      dataset, analysis_view, source_block,
      band_label, band_start_year, band_end_year,
      observed_band_end_year, band_full_width_year, observed_band_width_year, observed_band_fraction,
      full_band_observed_flag, partial_band_flag, no_observed_time_flag,
      at_risk_start_n, transition_events, remission_censor_n,
      right_censor_n, main_censor_n, person_years,
      crude_hazard_per_100py, conditional_interval_risk, conditional_interval_risk_observed_portion,
      max_followup_years, available_within_observed_followup,
      estimable_flag, availability_note,
      support_label, support_priority, primary_supported_flag,
      common_restricted_window_flag, support_basis, output_label
    )
}

## 🟠 Define: output-alignment checks ===============================
validate_supported_short_counts <- function(supported_df, risk_trajectory, analysis_totals_lookup, piecewise_interval_support) {
  analysis_lookup <- analysis_totals_lookup %>%
    rename(
      expected_analysis_subject_n_total = analysis_subject_n_total,
      expected_analysis_transition_events_total = analysis_transition_events_total
    )
  
  observed_noncontrast_check <- supported_df %>%
    filter(
      result_type == "observed_risk",
      source_block %in% c("observed_km_restricted", "observed_km_pooled_merged")
    ) %>%
    select(
      dataset, horizon_year,
      analysis_subject_n_total, analysis_transition_events_total,
      at_risk_subject_n_total, transition_events_total
    ) %>%
    left_join(
      risk_trajectory %>%
        select(dataset, horizon_year, n_risk, cumulative_transition_n),
      by = c("dataset", "horizon_year")
    ) %>%
    left_join(analysis_lookup, by = "dataset") %>%
    mutate(
      passed = equal_numeric_or_both_na(analysis_subject_n_total, expected_analysis_subject_n_total) &
        equal_numeric_or_both_na(analysis_transition_events_total, expected_analysis_transition_events_total) &
        equal_numeric_or_both_na(at_risk_subject_n_total, n_risk) &
        equal_numeric_or_both_na(transition_events_total, cumulative_transition_n)
    )
  
  contrast_check <- supported_df %>%
    filter(
      dataset == "PNU_vs_SNU",
      source_block == "observed_km_restricted",
      result_type == "absolute_risk_difference"
    ) %>%
    select(
      horizon_year,
      analysis_subject_n_total, analysis_transition_events_total,
      at_risk_subject_n_total, transition_events_total,
      pnu_event_n, snu_event_n
    ) %>%
    left_join(
      risk_trajectory %>%
        filter(dataset == "PNU") %>%
        transmute(
          horizon_year,
          expected_pnu_n_risk = n_risk,
          expected_pnu_transition_events = cumulative_transition_n
        ),
      by = "horizon_year"
    ) %>%
    left_join(
      risk_trajectory %>%
        filter(dataset == "SNU") %>%
        transmute(
          horizon_year,
          expected_snu_n_risk = n_risk,
          expected_snu_transition_events = cumulative_transition_n
        ),
      by = "horizon_year"
    ) %>%
    mutate(
      expected_analysis_subject_n_total = sum(analysis_lookup$expected_analysis_subject_n_total[match(c("PNU", "SNU"), analysis_lookup$dataset)]),
      expected_analysis_transition_events_total = sum(analysis_lookup$expected_analysis_transition_events_total[match(c("PNU", "SNU"), analysis_lookup$dataset)]),
      expected_at_risk_subject_n_total = expected_pnu_n_risk + expected_snu_n_risk,
      expected_transition_events_total = expected_pnu_transition_events + expected_snu_transition_events,
      passed = equal_numeric_or_both_na(analysis_subject_n_total, expected_analysis_subject_n_total) &
        equal_numeric_or_both_na(analysis_transition_events_total, expected_analysis_transition_events_total) &
        equal_numeric_or_both_na(at_risk_subject_n_total, expected_at_risk_subject_n_total) &
        equal_numeric_or_both_na(transition_events_total, expected_transition_events_total) &
        equal_numeric_or_both_na(pnu_event_n, expected_pnu_transition_events) &
        equal_numeric_or_both_na(snu_event_n, expected_snu_transition_events)
    )
  
  standardized_check <- supported_df %>%
    filter(
      source_block == "merged_restricted_cox_standardized",
      result_type %in% c("site_standardized_risk", "absolute_risk_difference")
    ) %>%
    mutate(
      expected_analysis_subject_n_total = analysis_lookup$expected_analysis_subject_n_total[match("merged", analysis_lookup$dataset)],
      expected_analysis_transition_events_total = analysis_lookup$expected_analysis_transition_events_total[match("merged", analysis_lookup$dataset)],
      passed = analysis_subject_n_total == expected_analysis_subject_n_total &
        analysis_transition_events_total == expected_analysis_transition_events_total &
        at_risk_subject_n_total == n_risk &
        transition_events_total == pnu_event_n + snu_event_n
    )
  
  piecewise_short_check <- supported_df %>%
    filter(result_type == "hazard_ratio", interval_label %in% c("0-1", "1-2")) %>%
    select(
      model_name, interval_label,
      analysis_subject_n_total, analysis_transition_events_total,
      at_risk_subject_n_total, transition_events_total,
      pnu_event_n, snu_event_n
    ) %>%
    left_join(
      piecewise_interval_support %>%
        select(
          interval_label,
          expected_at_risk_subject_n_total = at_risk_subject_n_total,
          expected_transition_events_total = transition_events_total,
          expected_pnu_event_n = pnu_event_n,
          expected_snu_event_n = snu_event_n
        ),
      by = "interval_label"
    ) %>%
    mutate(
      expected_analysis_subject_n_total = analysis_lookup$expected_analysis_subject_n_total[match("merged", analysis_lookup$dataset)],
      expected_analysis_transition_events_total = analysis_lookup$expected_analysis_transition_events_total[match("merged", analysis_lookup$dataset)],
      passed = equal_numeric_or_both_na(analysis_subject_n_total, expected_analysis_subject_n_total) &
        equal_numeric_or_both_na(analysis_transition_events_total, expected_analysis_transition_events_total) &
        equal_numeric_or_both_na(at_risk_subject_n_total, expected_at_risk_subject_n_total) &
        equal_numeric_or_both_na(transition_events_total, expected_transition_events_total) &
        equal_numeric_or_both_na(pnu_event_n, expected_pnu_event_n) &
        equal_numeric_or_both_na(snu_event_n, expected_snu_event_n)
    )
  
  failed_ids <- c(
    if (any(!observed_noncontrast_check$passed)) {
      sprintf("supported observed rows: %s", format_check_locations(observed_noncontrast_check %>% filter(!passed), c("dataset", "horizon_year")))
    } else {
      NULL
    },
    if (any(!contrast_check$passed)) {
      sprintf("supported observed contrast rows: %s", format_check_locations(contrast_check %>% filter(!passed), c("horizon_year")))
    } else {
      NULL
    },
    if (any(!standardized_check$passed)) {
      sprintf("supported standardized rows: %s", format_check_locations(standardized_check %>% filter(!passed), c("model_name", "horizon_year", "result_type")))
    } else {
      NULL
    },
    if (any(!piecewise_short_check$passed)) {
      sprintf("supported piecewise rows: %s", format_check_locations(piecewise_short_check %>% filter(!passed), c("model_name", "interval_label")))
    } else {
      NULL
    }
  )
  
  if (length(failed_ids) > 0L) {
    stop(sprintf("Supported short-horizon count alignment failed: %s", paste(failed_ids, collapse = " | ")), call. = FALSE)
  }
  
  invisible(TRUE)
}

validate_late_tail_counts <- function(late_df, risk_trajectory, analysis_totals_lookup, piecewise_interval_support) {
  analysis_lookup <- analysis_totals_lookup %>%
    rename(
      expected_analysis_subject_n_total = analysis_subject_n_total,
      expected_analysis_transition_events_total = analysis_transition_events_total
    )
  
  estimable_observed_check <- late_df %>%
    filter(result_type == "observed_risk", estimable_flag) %>%
    select(
      dataset, horizon_year,
      analysis_subject_n_total, analysis_transition_events_total,
      at_risk_subject_n_total, transition_events_total
    ) %>%
    left_join(
      risk_trajectory %>%
        select(dataset, horizon_year, n_risk, cumulative_transition_n),
      by = c("dataset", "horizon_year")
    ) %>%
    left_join(analysis_lookup, by = "dataset") %>%
    mutate(
      passed = equal_numeric_or_both_na(analysis_subject_n_total, expected_analysis_subject_n_total) &
        equal_numeric_or_both_na(analysis_transition_events_total, expected_analysis_transition_events_total) &
        equal_numeric_or_both_na(at_risk_subject_n_total, n_risk) &
        equal_numeric_or_both_na(transition_events_total, cumulative_transition_n)
    )
  
  nonestimable_observed_check <- late_df %>%
    filter(result_type == "observed_risk", !estimable_flag) %>%
    mutate(
      passed = is.na(at_risk_subject_n_total) &
        is.na(transition_events_total) &
        is.na(pnu_event_n) &
        is.na(snu_event_n)
    )
  
  piecewise_tail_check <- late_df %>%
    filter(result_type == "hazard_ratio", interval_label == ">2") %>%
    select(
      model_name, interval_label,
      analysis_subject_n_total, analysis_transition_events_total,
      at_risk_subject_n_total, transition_events_total,
      pnu_event_n, snu_event_n
    ) %>%
    left_join(
      piecewise_interval_support %>%
        select(
          interval_label,
          expected_at_risk_subject_n_total = at_risk_subject_n_total,
          expected_transition_events_total = transition_events_total,
          expected_pnu_event_n = pnu_event_n,
          expected_snu_event_n = snu_event_n
        ),
      by = "interval_label"
    ) %>%
    mutate(
      expected_analysis_subject_n_total = analysis_lookup$expected_analysis_subject_n_total[match("merged", analysis_lookup$dataset)],
      expected_analysis_transition_events_total = analysis_lookup$expected_analysis_transition_events_total[match("merged", analysis_lookup$dataset)],
      passed = equal_numeric_or_both_na(analysis_subject_n_total, expected_analysis_subject_n_total) &
        equal_numeric_or_both_na(analysis_transition_events_total, expected_analysis_transition_events_total) &
        equal_numeric_or_both_na(at_risk_subject_n_total, expected_at_risk_subject_n_total) &
        equal_numeric_or_both_na(transition_events_total, expected_transition_events_total) &
        equal_numeric_or_both_na(pnu_event_n, expected_pnu_event_n) &
        equal_numeric_or_both_na(snu_event_n, expected_snu_event_n)
    )
  
  failed_ids <- c(
    if (any(!estimable_observed_check$passed)) {
      sprintf("late estimable observed rows: %s", format_check_locations(estimable_observed_check %>% filter(!passed), c("dataset", "horizon_year")))
    } else {
      NULL
    },
    if (any(!nonestimable_observed_check$passed)) {
      sprintf("late non-estimable observed rows: %s", format_check_locations(nonestimable_observed_check %>% filter(!passed), c("dataset", "horizon_year")))
    } else {
      NULL
    },
    if (any(!piecewise_tail_check$passed)) {
      sprintf("late piecewise rows: %s", format_check_locations(piecewise_tail_check %>% filter(!passed), c("model_name", "interval_label")))
    } else {
      NULL
    }
  )
  
  if (length(failed_ids) > 0L) {
    stop(sprintf("Late-tail count alignment failed: %s", paste(failed_ids, collapse = " | ")), call. = FALSE)
  }
  
  invisible(TRUE)
}

validate_reporting_flags <- function(supported_df, late_df, piecewise_df) {
  supported_bootstrap_check <- supported_df %>%
    filter(!is.na(bootstrap_instability_flag), bootstrap_instability_flag) %>%
    mutate(
      passed = reporting_status == "reported_with_bootstrap_instability_flag" &
        plotting_status == "bootstrap_unstable" &
        plot_inclusion_flag
    )
  
  supported_suppressed_check <- supported_df %>%
    filter(reported_estimate_suppressed_flag) %>%
    mutate(
      passed = reporting_status == "suppressed" &
        plotting_status == "suppressed" &
        !plot_inclusion_flag
    )
  
  late_suppressed_check <- late_df %>%
    filter(reported_estimate_suppressed_flag) %>%
    mutate(
      passed = reporting_status == "suppressed" &
        plotting_status == "suppressed" &
        !plot_inclusion_flag
    )
  
  piecewise_suppressed_check <- piecewise_df %>%
    filter(reported_estimate_suppressed_flag) %>%
    mutate(
      passed = reporting_status == "suppressed" &
        plotting_status == "suppressed" &
        !plot_inclusion_flag
    )
  
  failed_ids <- c(
    if (nrow(supported_bootstrap_check) > 0L && any(!supported_bootstrap_check$passed)) {
      sprintf("supported bootstrap-flag rows: %s", format_check_locations(supported_bootstrap_check %>% filter(!passed), c("dataset", "horizon_year", "model_name", "result_type")))
    } else {
      NULL
    },
    if (nrow(supported_suppressed_check) > 0L && any(!supported_suppressed_check$passed)) {
      sprintf("supported suppressed rows: %s", format_check_locations(supported_suppressed_check %>% filter(!passed), c("dataset", "horizon_year", "model_name", "result_type")))
    } else {
      NULL
    },
    if (nrow(late_suppressed_check) > 0L && any(!late_suppressed_check$passed)) {
      sprintf("late suppressed rows: %s", format_check_locations(late_suppressed_check %>% filter(!passed), c("dataset", "horizon_year", "model_name", "result_type")))
    } else {
      NULL
    },
    if (nrow(piecewise_suppressed_check) > 0L && any(!piecewise_suppressed_check$passed)) {
      sprintf("piecewise suppressed rows: %s", format_check_locations(piecewise_suppressed_check %>% filter(!passed), c("model_name", "interval_label")))
    } else {
      NULL
    }
  )
  
  if (length(failed_ids) > 0L) {
    stop(sprintf("Reporting-flag alignment failed: %s", paste(failed_ids, collapse = " | ")), call. = FALSE)
  }
  
  invisible(TRUE)
}

## 🟠 Define: plot constructors ===============================
plot_risk_trajectory <- function(risk_trajectory_df) {
  ggplot(risk_trajectory_df, aes(x = horizon_year, y = risk, color = dataset, group = dataset)) +
    geom_line(na.rm = TRUE, linewidth = 0.8) +
    geom_point(na.rm = TRUE, size = 2) +
    scale_x_continuous(breaks = 1:10) +
    labs(
      title = "Stage 4 risk trajectory: observed KM risk by dataset",
      subtitle = "Values beyond each dataset's observed follow-up remain exported with availability flags and are omitted from plotted points.",
      x = "Horizon (years)",
      y = "Risk = 1 - KM survival",
      color = "Dataset"
    ) +
    theme_bw()
}

plot_supported_short <- function(supported_df) {
  dodge_width <- 0.35
  plot_df <- supported_df %>%
    filter(result_type == "absolute_risk_difference") %>%
    mutate(
      model_name = factor(model_name, levels = unique(model_name)),
      plotting_status = factor(plotting_status, levels = c("stable", "bootstrap_unstable", "suppressed"))
    )
  
  ggplot(plot_df, aes(x = factor(horizon_year), y = estimate, color = model_name, group = model_name)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_linerange(
      aes(ymin = conf_low, ymax = conf_high),
      position = position_dodge(width = dodge_width),
      na.rm = TRUE
    ) +
    geom_point(
      aes(shape = plotting_status),
      position = position_dodge(width = dodge_width),
      size = 2.5,
      na.rm = TRUE
    ) +
    geom_text(
      data = plot_df %>% filter(plotting_status == "bootstrap_unstable"),
      aes(label = "unstable"),
      position = position_dodge(width = dodge_width),
      vjust = -0.8,
      size = 3,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    scale_shape_manual(
      values = c(stable = 16, bootstrap_unstable = 1, suppressed = 4),
      drop = FALSE
    ) +
    labs(
      title = "Supported short-horizon contrast: PNU minus SNU absolute risk difference",
      subtitle = "Open points and `unstable` labels indicate rows kept for reporting but flagged for bootstrap instability; suppressed rows remain in the CSV only.",
      x = "Restricted horizon (years)",
      y = "Absolute risk difference",
      color = "Source",
      shape = "Reporting status"
    ) +
    theme_bw()
}

plot_piecewise_hr <- function(piecewise_df) {
  plot_df <- piecewise_df %>%
    filter(plot_inclusion_flag) %>%
    mutate(interval_label = factor(interval_label, levels = c("0-1", "1-2", ">2")))
  
  ggplot(plot_df, aes(x = interval_label, y = estimate_hr, color = model_name, group = model_name)) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
    geom_pointrange(aes(ymin = conf_low_hr, ymax = conf_high_hr), position = position_dodge(width = 0.3), na.rm = TRUE) +
    scale_y_log10() +
    labs(
      title = "Piecewise site effect in merged data: SNU versus PNU hazard ratio",
      subtitle = "Intervals flagged as non-estimable because of sparse or zero site-specific support are omitted from plotted point ranges and remain explicitly suppressed in the CSV.",
      x = "Time interval (years)",
      y = "Hazard ratio (log scale)",
      color = "Model"
    ) +
    theme_bw()
}

plot_piecewise_support <- function(piecewise_support_df) {
  plot_df <- piecewise_support_df %>%
    select(interval_label, pnu_event_n, snu_event_n) %>%
    pivot_longer(
      cols = c(pnu_event_n, snu_event_n),
      names_to = "site_name",
      values_to = "event_n"
    ) %>%
    mutate(
      interval_label = factor(interval_label, levels = c("0-1", "1-2", ">2")),
      site_name = factor(
        case_when(
          site_name == "pnu_event_n" ~ pnu_site_label,
          TRUE ~ snu_site_label
        ),
        levels = c(pnu_site_label, snu_site_label)
      )
    )
  
  ggplot(plot_df, aes(x = interval_label, y = event_n, fill = site_name)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, na.rm = TRUE) +
    labs(
      title = "Piecewise interval support: site-specific transition events",
      subtitle = "Large or undefined site-effect estimates after 2 years should be suppressed when one site has sparse event support.",
      x = "Time interval (years)",
      y = "Transition events",
      fill = "Site"
    ) +
    theme_bw()
}

plot_hazard_pattern <- function(hazard_df) {
  ggplot(hazard_df, aes(x = band_label, y = crude_hazard_per_100py, color = dataset, group = dataset)) +
    geom_line(na.rm = TRUE, linewidth = 0.8) +
    geom_point(na.rm = TRUE, size = 2) +
    labs(
      title = "Hazard-pattern summary: crude transition hazard per 100 person-years",
      subtitle = "Bands extending beyond observed follow-up are flagged in the CSV; crude hazards use observed person-time only.",
      x = "Time band (years)",
      y = "Crude hazard per 100 person-years",
      color = "Dataset"
    ) +
    theme_bw()
}

plot_late_tail <- function(late_df) {
  plot_df <- late_df %>%
    filter(result_type == "observed_risk", plot_inclusion_flag) %>%
    mutate(dataset = factor(dataset, levels = unique(dataset)))
  
  ggplot(plot_df, aes(x = factor(horizon_year), y = estimate, color = dataset, group = dataset)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_point(size = 2, na.rm = TRUE) +
    geom_line(na.rm = TRUE, linewidth = 0.8) +
    labs(
      title = "Late-horizon tail contrast: observed risk summaries",
      subtitle = "Rows beyond observed follow-up remain exported with availability flags; only reported observed-risk points are drawn.",
      x = "Horizon (years)",
      y = "Observed KM risk",
      color = "Dataset"
    ) +
    theme_bw()
}

# 🔴 Read: Stage 1 backbone inputs ===============================
stage1_inputs <- read_stage1_inputs()

# 🔴 Validate: inherited registries and prepared datasets ===============================
validate_stage1_inputs(stage1_inputs)

common_horizons_year <- sort(unique(as.integer(stage1_inputs$horizon_registry$horizon_year)))
risk_thresholds <- sort(unique(as.numeric(stage1_inputs$threshold_registry$threshold)))
support_lookup <- make_stage1_support_lookup(stage1_inputs$horizon_registry)

analysis_datasets <- stage1_inputs$analysis_datasets
analysis_datasets <- lapply(analysis_datasets, add_model_time_year)
analysis_totals_lookup <- make_analysis_totals_lookup(analysis_datasets)

stage4_dataset_summary <- make_stage4_dataset_summary(
  analysis_datasets = analysis_datasets,
  stage1_dataset_registry = stage1_inputs$dataset_registry,
  support_lookup = support_lookup
)

formula_registry <- stage1_inputs$formula_registry %>%
  mutate(
    dataset = as.character(dataset),
    formula_name = as.character(formula_name),
    formula_rhs = as.character(formula_rhs)
  )

merged_site_formula_registry <- formula_registry %>%
  filter(dataset == "merged", formula_name %in% c("site_added", "site_interaction")) %>%
  arrange(match(formula_name, c("site_added", "site_interaction")))

if (nrow(merged_site_formula_registry) != 2L) {
  stop("Stage 1 formula registry must contain `merged` site_added and site_interaction formulas.", call. = FALSE)
}

# 🔴 Compute: observed risk trajectories ===============================
risk_trajectory <- bind_rows(lapply(names(analysis_datasets), function(dataset_name) {
  compute_km_risk_trajectory(
    df = analysis_datasets[[dataset_name]],
    dataset_name = dataset_name,
    support_lookup = support_lookup
  )
})) %>%
  mutate(
    output_label = if_else(
      horizon_year %in% supported_short_horizons_year,
      "supported_short_horizon_contrast",
      if_else(horizon_year %in% late_tail_horizons_year, "late_horizon_tail_contrast", "stage4_context_only")
    )
  ) %>%
  arrange(match(dataset, c("PNU", "SNU", "merged")), horizon_year)

# 🔴 Estimate: short-horizon contrasts ===============================
pnu_data <- analysis_datasets[["PNU"]]
snu_data <- analysis_datasets[["SNU"]]
merged_data <- analysis_datasets[["merged"]]

observed_supported_rows <- bind_rows(lapply(supported_short_horizons_year, function(h) {
  bootstrap_observed_supported(
    pnu_df = pnu_data,
    snu_df = snu_data,
    horizon_year = h,
    n_boot = bootstrap_iterations,
    seed_value = bootstrap_seed + h * 1000L,
    support_lookup = support_lookup
  )
}))

merged_site_window_counts <- bind_rows(lapply(supported_short_horizons_year, function(h) {
  bind_cols(
    tibble(horizon_year = as.numeric(h)),
    summarize_horizon_window_counts(merged_data, h, require_full_observation = TRUE)
  )
}))

merged_site_free_context_rows <- risk_trajectory %>%
  filter(dataset == "merged", horizon_year %in% supported_short_horizons_year) %>%
  left_join(merged_site_window_counts, by = "horizon_year") %>%
  transmute(
    dataset = dataset,
    analysis_view = analysis_view,
    source_block = "observed_km_pooled_merged",
    model_name = model_name,
    formula_rhs = NA_character_,
    horizon_year = horizon_year,
    window_year = horizon_year,
    result_type = "observed_risk",
    site_counterfactual = NA_character_,
    contrast_direction = NA_character_,
    interval_label = NA_character_,
    estimate = risk,
    conf_low = risk_conf_low,
    conf_high = risk_conf_high,
    n_risk = n_risk,
    max_followup_years = max_followup_years,
    available_within_observed_followup = available_within_observed_followup,
    bootstrap_success = NA_integer_,
    bootstrap_iterations = NA_integer_,
    bootstrap_success_rate = NA_real_,
    bootstrap_instability_flag = NA,
    bootstrap_note = NA_character_,
    estimable_flag = estimable_flag,
    availability_note = availability_note,
    support_issue_reason = NA_character_,
    analysis_subject_n_total = analysis_subject_n_total,
    analysis_transition_events_total = analysis_transition_events_total,
    at_risk_subject_n_total = n_risk,
    transition_events_total = transition_events_total,
    pnu_event_n = NA_real_,
    snu_event_n = NA_real_,
    support_label = support_label,
    support_priority = support_priority,
    primary_supported_flag = primary_supported_flag,
    common_restricted_window_flag = common_restricted_window_flag,
    support_basis = make_support_basis("merged_pooled_observed_support"),
    reported_estimate_suppressed_flag = !estimable_flag,
    output_label = "supported_short_horizon_contrast"
  )

restricted_model_outputs <- list()
restricted_supported_rows_list <- list()
restricted_index <- 1L

for (i in seq_len(nrow(merged_site_formula_registry))) {
  model_name <- paste0("restricted_", merged_site_formula_registry$formula_name[i])
  formula_rhs <- merged_site_formula_registry$formula_rhs[i]
  
  for (h in supported_short_horizons_year) {
    restricted_model_outputs[[paste(model_name, h, sep = "__")]] <- fit_restricted_cox(
      df = merged_data,
      horizon_year = h,
      formula_rhs = formula_rhs,
      model_name = model_name
    )
    
    restricted_supported_rows_list[[restricted_index]] <- bootstrap_adjusted_supported(
      df = merged_data,
      horizon_year = h,
      formula_rhs = formula_rhs,
      model_name = model_name,
      n_boot = bootstrap_iterations,
      seed_value = bootstrap_seed + i * 10000L + h * 100L,
      support_lookup = support_lookup
    )
    
    restricted_index <- restricted_index + 1L
  }
}

restricted_supported_rows <- bind_rows(restricted_supported_rows_list)

# 🔴 Fit: merged piecewise site-effect models ===============================
piecewise_merged_data <- make_piecewise_dataset(merged_data)

piecewise_interval_support <- summarize_piecewise_interval_support(
  split_df = piecewise_merged_data,
  min_events_per_site = piecewise_min_events_per_site,
  min_subjects_per_site = piecewise_min_subjects_per_site
)

piecewise_model_specs <- tibble(
  model_name = c("piecewise_unadjusted", "piecewise_site_added", "piecewise_site_interaction"),
  covariate_rhs = c("", "age_s + sex_num", "age_s + sex_num + age_s:sex_num")
)

piecewise_model_objects <- list()
piecewise_rows_list <- list()

for (i in seq_len(nrow(piecewise_model_specs))) {
  model_bundle <- fit_piecewise_site_model_bundle(
    split_df = piecewise_merged_data,
    model_name = piecewise_model_specs$model_name[i],
    covariate_rhs = piecewise_model_specs$covariate_rhs[i],
    interval_support_df = piecewise_interval_support
  )
  
  piecewise_model_objects[[piecewise_model_specs$model_name[i]]] <- model_bundle
  piecewise_rows_list[[i]] <- if (is.null(model_bundle)) tibble() else model_bundle$summary_table
}

piecewise_site_effect <- bind_rows(piecewise_rows_list) %>%
  arrange(match(model_name, piecewise_model_specs$model_name), match(interval_label, c("0-1", "1-2", ">2"))) %>%
  mutate(
    reporting_status = dplyr::case_when(
      reported_estimate_suppressed_flag ~ "suppressed",
      TRUE ~ "reported"
    ),
    plotting_status = dplyr::case_when(
      reported_estimate_suppressed_flag ~ "suppressed",
      TRUE ~ "stable"
    ),
    plot_inclusion_flag = estimable_flag & !reported_estimate_suppressed_flag,
    reporting_note = dplyr::case_when(
      reported_estimate_suppressed_flag ~ coalesce_character(support_issue_reason, availability_note),
      TRUE ~ availability_note
    )
  )

# 🔴 Summarize: interval hazards and tail contrasts ===============================
hazard_pattern_summary <- bind_rows(lapply(names(analysis_datasets), function(dataset_name) {
  summarize_hazard_bands(
    df = analysis_datasets[[dataset_name]],
    dataset_name = dataset_name,
    band_breaks = hazard_band_breaks_year,
    support_lookup = support_lookup
  )
})) %>%
  arrange(match(dataset, c("PNU", "SNU", "merged")), band_start_year)

supported_piecewise_rows <- piecewise_site_effect %>%
  filter(interval_label %in% c("0-1", "1-2")) %>%
  transmute(
    dataset = dataset,
    analysis_view = analysis_view,
    source_block = source_block,
    model_name = model_name,
    formula_rhs = formula_rhs,
    horizon_year = dplyr::case_when(interval_label == "0-1" ~ 1, interval_label == "1-2" ~ 2, TRUE ~ NA_real_),
    window_year = dplyr::case_when(interval_label == "0-1" ~ 1, interval_label == "1-2" ~ 2, TRUE ~ NA_real_),
    result_type = "hazard_ratio",
    site_counterfactual = NA_character_,
    contrast_direction = contrast_direction,
    interval_label = interval_label,
    estimate = estimate_hr,
    conf_low = conf_low_hr,
    conf_high = conf_high_hr,
    n_risk = NA_real_,
    max_followup_years = max_followup_years,
    available_within_observed_followup = available_within_observed_followup,
    bootstrap_success = NA_integer_,
    bootstrap_iterations = NA_integer_,
    bootstrap_success_rate = NA_real_,
    bootstrap_instability_flag = NA,
    bootstrap_note = NA_character_,
    estimable_flag = estimable_flag,
    availability_note = availability_note,
    support_issue_reason = support_issue_reason,
    analysis_subject_n_total = n_subjects,
    analysis_transition_events_total = n_events_total,
    at_risk_subject_n_total = at_risk_subject_n_total,
    transition_events_total = transition_events_total,
    pnu_event_n = pnu_event_n,
    snu_event_n = snu_event_n,
    support_label = support_label,
    support_priority = support_priority,
    primary_supported_flag = primary_supported_flag,
    common_restricted_window_flag = common_restricted_window_flag,
    support_basis = support_basis,
    reported_estimate_suppressed_flag = reported_estimate_suppressed_flag,
    output_label = "supported_short_horizon_contrast"
  )

supported_short_horizon_contrast <- bind_rows(
  observed_supported_rows,
  merged_site_free_context_rows,
  restricted_supported_rows,
  supported_piecewise_rows
) %>%
  arrange(
    factor(output_label, levels = c("supported_short_horizon_contrast")),
    factor(result_type, levels = c("observed_risk", "site_standardized_risk", "absolute_risk_difference", "hazard_ratio")),
    horizon_year,
    model_name,
    dataset
  ) %>%
  mutate(
    reporting_status = dplyr::case_when(
      reported_estimate_suppressed_flag ~ "suppressed",
      !is.na(bootstrap_instability_flag) & bootstrap_instability_flag ~ "reported_with_bootstrap_instability_flag",
      TRUE ~ "reported"
    ),
    plotting_status = dplyr::case_when(
      reported_estimate_suppressed_flag ~ "suppressed",
      !is.na(bootstrap_instability_flag) & bootstrap_instability_flag ~ "bootstrap_unstable",
      TRUE ~ "stable"
    ),
    plot_inclusion_flag = estimable_flag & !reported_estimate_suppressed_flag,
    reporting_note = dplyr::case_when(
      reporting_status == "reported_with_bootstrap_instability_flag" ~ coalesce_character(bootstrap_note, "bootstrap_instability_flagged"),
      reporting_status == "suppressed" ~ coalesce_character(support_issue_reason, availability_note),
      TRUE ~ availability_note
    )
  )

late_risk_rows <- risk_trajectory %>%
  filter(horizon_year %in% late_tail_horizons_year) %>%
  left_join(analysis_totals_lookup, by = "dataset") %>%
  mutate(
    at_risk_subject_n_total = ifelse(estimable_flag, n_risk, NA_real_),
    transition_events_total = ifelse(estimable_flag, cumulative_transition_n, NA_real_),
    pnu_event_n = case_when(
      dataset == "PNU" & estimable_flag ~ cumulative_transition_n,
      TRUE ~ NA_real_
    ),
    snu_event_n = case_when(
      dataset == "SNU" & estimable_flag ~ cumulative_transition_n,
      TRUE ~ NA_real_
    )
  ) %>%
  transmute(
    dataset = dataset,
    analysis_view = analysis_view,
    source_block = "observed_km",
    model_name = model_name,
    formula_rhs = NA_character_,
    horizon_year = horizon_year,
    window_year = NA_real_,
    result_type = "observed_risk",
    site_counterfactual = NA_character_,
    contrast_direction = NA_character_,
    interval_label = NA_character_,
    estimate = risk,
    conf_low = risk_conf_low,
    conf_high = risk_conf_high,
    n_risk = n_risk,
    max_followup_years = max_followup_years,
    available_within_observed_followup = available_within_observed_followup,
    bootstrap_success = NA_integer_,
    bootstrap_iterations = NA_integer_,
    bootstrap_success_rate = NA_real_,
    bootstrap_instability_flag = NA,
    bootstrap_note = NA_character_,
    estimable_flag = estimable_flag,
    availability_note = availability_note,
    support_issue_reason = case_when(
      estimable_flag ~ NA_character_,
      TRUE ~ "beyond_max_observed_followup"
    ),
    analysis_subject_n_total = analysis_subject_n_total,
    analysis_transition_events_total = analysis_transition_events_total,
    at_risk_subject_n_total = at_risk_subject_n_total,
    transition_events_total = transition_events_total,
    pnu_event_n = pnu_event_n,
    snu_event_n = snu_event_n,
    support_label = support_label,
    support_priority = support_priority,
    primary_supported_flag = primary_supported_flag,
    common_restricted_window_flag = common_restricted_window_flag,
    support_basis = support_basis,
    reported_estimate_suppressed_flag = !estimable_flag,
    output_label = "late_horizon_tail_contrast"
  )

tail_piecewise_rows <- piecewise_site_effect %>%
  filter(interval_label == ">2") %>%
  transmute(
    dataset = dataset,
    analysis_view = analysis_view,
    source_block = source_block,
    model_name = model_name,
    formula_rhs = formula_rhs,
    horizon_year = NA_real_,
    window_year = NA_real_,
    result_type = "hazard_ratio",
    site_counterfactual = NA_character_,
    contrast_direction = contrast_direction,
    interval_label = interval_label,
    estimate = estimate_hr,
    conf_low = conf_low_hr,
    conf_high = conf_high_hr,
    n_risk = NA_real_,
    max_followup_years = max_followup_years,
    available_within_observed_followup = available_within_observed_followup,
    bootstrap_success = NA_integer_,
    bootstrap_iterations = NA_integer_,
    bootstrap_success_rate = NA_real_,
    bootstrap_instability_flag = NA,
    bootstrap_note = NA_character_,
    estimable_flag = estimable_flag,
    availability_note = availability_note,
    support_issue_reason = support_issue_reason,
    analysis_subject_n_total = n_subjects,
    analysis_transition_events_total = n_events_total,
    at_risk_subject_n_total = at_risk_subject_n_total,
    transition_events_total = transition_events_total,
    pnu_event_n = pnu_event_n,
    snu_event_n = snu_event_n,
    support_label = support_label,
    support_priority = support_priority,
    primary_supported_flag = primary_supported_flag,
    common_restricted_window_flag = common_restricted_window_flag,
    support_basis = support_basis,
    reported_estimate_suppressed_flag = reported_estimate_suppressed_flag,
    output_label = "late_horizon_tail_contrast"
  )

late_horizon_tail_contrast <- bind_rows(
  late_risk_rows,
  tail_piecewise_rows
) %>%
  arrange(
    factor(result_type, levels = c("observed_risk", "hazard_ratio")),
    horizon_year,
    model_name,
    dataset
  ) %>%
  mutate(
    reporting_status = dplyr::case_when(
      reported_estimate_suppressed_flag ~ "suppressed",
      TRUE ~ "reported"
    ),
    plotting_status = dplyr::case_when(
      reported_estimate_suppressed_flag ~ "suppressed",
      TRUE ~ "stable"
    ),
    plot_inclusion_flag = estimable_flag & !reported_estimate_suppressed_flag,
    reporting_note = dplyr::case_when(
      reported_estimate_suppressed_flag ~ coalesce_character(support_issue_reason, availability_note),
      TRUE ~ availability_note
    )
  )

# 🔴 Check: count semantics and output alignment ===============================
validate_supported_short_counts(
  supported_df = supported_short_horizon_contrast,
  risk_trajectory = risk_trajectory,
  analysis_totals_lookup = analysis_totals_lookup,
  piecewise_interval_support = piecewise_interval_support
)

validate_late_tail_counts(
  late_df = late_horizon_tail_contrast,
  risk_trajectory = risk_trajectory,
  analysis_totals_lookup = analysis_totals_lookup,
  piecewise_interval_support = piecewise_interval_support
)

validate_reporting_flags(
  supported_df = supported_short_horizon_contrast,
  late_df = late_horizon_tail_contrast,
  piecewise_df = piecewise_site_effect
)

# 🔴 Assemble: metadata and manifest registries ===============================
stage4_metadata_registry <- tibble::tribble(
  ~metadata_group, ~metadata_name, ~metadata_value,
  "stage", "stage_name", "Stage 4 timing-difference",
  "stage", "stage_role", "Separate timing difference from cure interpretation before later cure blocks.",
  "stage", "canonical_framework", "5.Model Specification Framework_🇬🇧ENG.md",
  "stage", "site_effect_implementation", "Piecewise site effect implemented as 0-1, 1-2, and >2 years in merged Cox models.",
  "stage", "tail_visualization_rule", "Late-tail plot book page shows observed risk only; tail hazard-ratio estimates remain in CSV with estimability flags.",
  "inputs", "data_path", normalize_existing_path(data_path),
  "inputs", "export_path", normalize_existing_path(export_path),
  "inputs", "data_path_leaf", make_portable_file_reference(data_path),
  "inputs", "export_path_leaf", make_portable_file_reference(export_path),
  "inputs", "stage1_analysis_datasets_file", normalize_existing_path(stage1_analysis_datasets_file),
  "inputs", "stage1_dataset_registry_file", normalize_existing_path(stage1_dataset_registry_file),
  "inputs", "stage1_scaling_registry_file", normalize_existing_path(stage1_scaling_registry_file),
  "inputs", "stage1_metadata_registry_file", normalize_existing_path(stage1_metadata_registry_file),
  "inputs", "stage1_formula_registry_file", normalize_existing_path(stage1_formula_registry_file),
  "inputs", "stage1_horizon_registry_file", normalize_existing_path(stage1_horizon_registry_file),
  "inputs", "stage1_threshold_registry_file", normalize_existing_path(stage1_threshold_registry_file),
  "inputs", "stage1_analysis_datasets_filename", make_portable_file_reference(stage1_analysis_datasets_file),
  "inputs", "stage1_dataset_registry_filename", make_portable_file_reference(stage1_dataset_registry_file),
  "inputs", "stage1_scaling_registry_filename", make_portable_file_reference(stage1_scaling_registry_file),
  "inputs", "stage1_metadata_registry_filename", make_portable_file_reference(stage1_metadata_registry_file),
  "inputs", "stage1_formula_registry_filename", make_portable_file_reference(stage1_formula_registry_file),
  "inputs", "stage1_horizon_registry_filename", make_portable_file_reference(stage1_horizon_registry_file),
  "inputs", "stage1_threshold_registry_filename", make_portable_file_reference(stage1_threshold_registry_file),
  "time", "analysis_time_variable", "days_followup",
  "time", "reporting_time_variable", "time_year = days_followup / 365.25",
  "time", "horizon_vector", paste(common_horizons_year, collapse = ","),
  "time", "piecewise_cuts_year", paste(piecewise_cuts_year, collapse = ","),
  "event", "event_definition", "status_num == 1",
  "event", "main_censoring_definition", "status_num %in% c(0, 2)",
  "timing", "supported_short_horizons_year", paste(supported_short_horizons_year, collapse = ","),
  "timing", "late_tail_horizons_year", paste(late_tail_horizons_year, collapse = ","),
  "timing", "hazard_band_breaks_year", paste(hazard_band_breaks_year, collapse = ","),
  "timing", "threshold_vector", paste(format(risk_thresholds, trim = TRUE, scientific = FALSE), collapse = ","),
  "timing", "support_hierarchy_source", "Horizon-based support labels for PNU, SNU, and merged outputs are inherited from Stage 1 horizon_registry rather than redefined locally.",
  "timing", "observed_common_window_contrast_rule", "Separate-cohort observed PNU-vs-SNU contrasts use the more conservative common-window support label: 1 year primary_supported and 2 years sensitivity.",
  "timing", "merged_restricted_window_contrast_rule", "Merged restricted standardized contrasts inherit merged-horizon support from Stage 1, so 1-year and 2-year merged restricted rows remain primary_supported.",
  "bootstrap", "bootstrap_iterations", as.character(bootstrap_iterations),
  "bootstrap", "bootstrap_seed", as.character(bootstrap_seed),
  "bootstrap", "bootstrap_min_success_rate", format(bootstrap_min_success_rate, trim = TRUE),
  "bootstrap", "bootstrap_reporting_rule", "Rows with retained estimates but bootstrap_success_rate below the minimum stay reported, receive reporting_status = reported_with_bootstrap_instability_flag, and use a bootstrap_unstable plotting marker.",
  "piecewise", "piecewise_min_events_per_site", as.character(piecewise_min_events_per_site),
  "piecewise", "piecewise_min_subjects_per_site", as.character(piecewise_min_subjects_per_site),
  "piecewise", "piecewise_nonestimable_rule", "If one site has too few subjects, zero or sparse events, no person-time, or the interval coefficient is not estimable, suppress the reported piecewise HR and flag the row as non-estimable.",
  "piecewise", "piecewise_export_rule", "When piecewise HR rows are suppressed as non-estimable, exported raw HR diagnostics are blanked in the CSV; underlying fit objects remain in the RDS bundle.",
  "piecewise", "late_horizon_availability_rule", "Late-horizon observed-risk rows must retain explicit flags when the requested horizon exceeds observed follow-up.",
  "comparison", "short_output_label", "supported_short_horizon_contrast",
  "comparison", "tail_output_label", "late_horizon_tail_contrast",
  "comparison", "support_basis_rule", "Support labels and support_basis are exported together so readers can distinguish separate-cohort observed support, merged pooled observation, merged restricted-window standardization, merged piecewise interval support, and observed person-time band support.",
  "comparison", "short_output_block_rule", "The supported_short_horizon_contrast block contains all 1-year and 2-year common-window rows, while row-level support_label and primary_supported_flag continue to distinguish primary from sensitivity status.",
  "comparison", "analysis_subject_n_total_rule", "analysis_subject_n_total stores the full analysis cohort size for the row source, distinct from horizon- or interval-specific at_risk_subject_n_total.",
  "comparison", "analysis_transition_events_total_rule", "analysis_transition_events_total stores all observed transition events over full follow-up for the row source, distinct from horizon- or interval-specific transition_events_total.",
  "comparison", "at_risk_subject_n_total_rule", "at_risk_subject_n_total stores the true horizon or interval risk-set count; horizon rows beyond observed follow-up leave it blank.",
  "comparison", "transition_events_total_rule", "transition_events_total stores within-horizon or within-interval observed events only; horizon rows beyond observed follow-up leave it blank rather than reusing all-follow-up totals.",
  "comparison", "count_alignment_validation_rule", "Stage 4 stops if supported or late-tail comparison tables do not align their count fields with risk_trajectory or piecewise_interval_support.",
  "comparison", "reporting_alignment_validation_rule", "Stage 4 stops if bootstrap-instability, suppression, and plotting-status flags are not synchronized across supported, late-tail, and piecewise outputs.",
  "hazard_band", "partial_band_rule", "Hazard-band summaries retain nominal bands but flag and quantify partial observation when follow-up ends within a band.",
  "hazard_band", "partial_band_reporting_rule", "For partially observed hazard bands, crude hazard uses observed person-time only and conditional_interval_risk is reserved for fully observed bands; observed-portion risk is exported separately.",
  "outputs", "save_folder_rule", "Write all Stage 4 outputs into one export_path folder without creating extra subfolders.",
  "outputs", "site_interpretation_rule", "Interpret site as a proxy for broader treatment context, care pathway, selection, or follow-up structure rather than a clean treatment effect.",
  "outputs", "remission_note", "The main Stage 4 analysis treats remission as censoring; remission-sensitive competing-risk or multi-state analysis belongs to a later dedicated sensitivity stage.",
  "outputs", "manifest_portability_rule", "stage4_export_manifest.csv stores both absolute file_path and portable_file_reference so the manifest remains readable after path relocation.",
  "outputs", "stage4_revision_note", "This revised Stage 4 code keeps the original count-semantics fixes and additionally tags bootstrap-unstable short-horizon rows for reporting and plotting, while adding portable file references to the metadata and manifest outputs."
)

stage4_export_manifest <- tibble(
  file_name = c(
    "stage4_dataset_summary.csv",
    "stage4_risk_trajectory.csv",
    "stage4_hazard_pattern_summary.csv",
    "stage4_piecewise_interval_support.csv",
    "stage4_piecewise_site_effect.csv",
    "stage4_supported_short_horizon_contrast.csv",
    "stage4_late_horizon_tail_contrast.csv",
    "stage4_metadata_registry.csv",
    "stage4_model_objects.rds",
    "stage4_plot_book.pdf",
    "stage4_export_manifest.csv"
  ),
  object_name = c(
    "stage4_dataset_summary",
    "risk_trajectory",
    "hazard_pattern_summary",
    "piecewise_interval_support",
    "piecewise_site_effect",
    "supported_short_horizon_contrast",
    "late_horizon_tail_contrast",
    "stage4_metadata_registry",
    "stage4_model_objects",
    "stage4_plot_book",
    "stage4_export_manifest"
  ),
  description = c(
    "Dataset-level Stage 4 summary with Stage 1 aligned support counts",
    "Observed KM risk trajectory on the common 1-10 year grid with Stage 1 aligned support labels and support-basis annotation",
    "Hazard-pattern summary across pre-specified time bands with full/partial observation diagnostics and dataset-specific Stage 1 aligned support labels",
    "Piecewise interval-support table showing site-specific risk-set and event support",
    "Piecewise merged site-effect hazard ratio table with reported estimates only for estimable intervals",
    "Supported short-horizon contrast table with analysis totals separated from horizon-specific support counts, availability, bootstrap diagnostics, explicit reporting-status flags, and support-basis labels",
    "Late-horizon tail contrast table with analysis totals separated from horizon-specific support counts, availability, estimability, reporting-status flags, and support-basis diagnostics",
    "Stage 4 metadata registry",
    "Restricted and piecewise model objects bundled as RDS",
    "Single PDF plot book generated from exported data frames",
    "Manifest of all Stage 4 exported files"
  )
) %>%
  mutate(
    file_path = normalize_existing_path(file.path(export_path, file_name)),
    portable_file_reference = file_name
  )

stage4_model_objects <- list(
  restricted_models = restricted_model_outputs,
  piecewise_models = piecewise_model_objects,
  piecewise_interval_support = piecewise_interval_support,
  stage1_support_lookup = support_lookup,
  analysis_totals_lookup = analysis_totals_lookup,
  config = list(
    data_path = data_path,
    export_path = export_path,
    bootstrap_iterations = bootstrap_iterations,
    bootstrap_seed = bootstrap_seed,
    bootstrap_min_success_rate = bootstrap_min_success_rate,
    piecewise_cuts_year = piecewise_cuts_year,
    piecewise_min_events_per_site = piecewise_min_events_per_site,
    piecewise_min_subjects_per_site = piecewise_min_subjects_per_site,
    hazard_band_breaks_year = hazard_band_breaks_year,
    supported_short_horizons_year = supported_short_horizons_year,
    late_tail_horizons_year = late_tail_horizons_year,
    thresholds_from_stage1 = risk_thresholds,
    horizons_from_stage1 = common_horizons_year,
    stage1_horizon_registry_file = stage1_horizon_registry_file,
    support_hierarchy_source = "Stage 4 inherits horizon-based support labels from Stage 1 horizon_registry rather than redefining them locally.",
    observed_common_window_contrast_rule = "Separate-cohort observed PNU-vs-SNU contrasts use the more conservative common-window support label: 1 year primary_supported and 2 years sensitivity.",
    merged_restricted_window_contrast_rule = "Merged restricted standardized contrasts inherit merged-horizon support from Stage 1.",
    count_semantics_rule = "analysis_subject_n_total and analysis_transition_events_total store all-follow-up totals, while at_risk_subject_n_total and transition_events_total store horizon- or interval-specific support only.",
    late_nonestimable_rule = "Late-horizon rows beyond observed follow-up blank horizon-specific count fields rather than recycling all-follow-up totals.",
    partial_band_rule = "Hazard-band summaries retain nominal bands but flag and quantify partial observation when follow-up ends within a band.",
    bootstrap_reporting_rule = "Rows with retained estimates but bootstrap_success_rate below the minimum remain reported, receive reporting_status = reported_with_bootstrap_instability_flag, and use a bootstrap_unstable plotting marker.",
    manifest_portability_rule = "stage4_export_manifest.csv stores both absolute file_path and portable_file_reference so the manifest remains readable after path relocation."
  )
)

# 🔴 Export: tables, model bundle, and plot book ===============================
readr::write_csv(stage4_dataset_summary, file.path(export_path, "stage4_dataset_summary.csv"))
readr::write_csv(risk_trajectory, file.path(export_path, "stage4_risk_trajectory.csv"))
readr::write_csv(hazard_pattern_summary, file.path(export_path, "stage4_hazard_pattern_summary.csv"))
readr::write_csv(piecewise_interval_support, file.path(export_path, "stage4_piecewise_interval_support.csv"))
readr::write_csv(piecewise_site_effect, file.path(export_path, "stage4_piecewise_site_effect.csv"))
readr::write_csv(supported_short_horizon_contrast, file.path(export_path, "stage4_supported_short_horizon_contrast.csv"))
readr::write_csv(late_horizon_tail_contrast, file.path(export_path, "stage4_late_horizon_tail_contrast.csv"))
readr::write_csv(stage4_metadata_registry, file.path(export_path, "stage4_metadata_registry.csv"))
saveRDS(stage4_model_objects, file.path(export_path, "stage4_model_objects.rds"))

grDevices::pdf(file.path(export_path, "stage4_plot_book.pdf"), width = 11, height = 8.5, onefile = TRUE)
print(plot_risk_trajectory(risk_trajectory))
print(plot_supported_short(supported_short_horizon_contrast))
print(plot_piecewise_hr(piecewise_site_effect))
print(plot_piecewise_support(piecewise_interval_support))
print(plot_hazard_pattern(hazard_pattern_summary))
print(plot_late_tail(late_horizon_tail_contrast))
grDevices::dev.off()

readr::write_csv(stage4_export_manifest, file.path(export_path, "stage4_export_manifest.csv"))
