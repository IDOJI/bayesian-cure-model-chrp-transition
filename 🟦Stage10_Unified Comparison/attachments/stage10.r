rm(list = ls())

# 🔴 Configure: Stage-10 paths and comparison controls ===============================
data_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.데이터 분석 New/🟦Stage1/attachments'
export_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.데이터 분석 New/🟦Stage10_unified comparison/attachments'

stage3_export_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.데이터 분석 New/🟦Stage3/attachments'
stage5_export_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.데이터 분석 New/🟦Stage5_non-cure block/attachments'
stage7_export_path <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/STAGE7/stage7_run__2026-03-24_shard50_auto20_v1__00009fe20000'
stage8_export_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.데이터 분석 New/🟦Stage8_베이지안 방법/attachments'
stage9_export_path <- NA_character_
stage6_screening_file <- NA_character_

allow_partial_stage10_build <- TRUE
require_stage3_outputs <- TRUE
require_stage5_outputs <- TRUE
require_stage7_outputs <- FALSE
require_stage8_outputs <- FALSE
include_stage3_sensitivity_benchmarks <- FALSE
keep_only_stage8_admissible_fits <- TRUE
include_stage8_unsupportive_sensitivity_fits <- TRUE

risk_thresholds_override <- NULL
common_horizons_year <- 1:10
plot_horizons_year <- c(1L, 2L, 5L, 10L)

pnu_site_label <- 'PNU'
snu_site_label <- 'SNU'


# 🔴 Initialize: packages and runtime options ===============================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

dataset_order <- c('PNU', 'SNU', 'merged')
model_class_order <- c('km_benchmark', 'no_cure', 'frequentist_cure', 'bayesian_cure', 'remission_sensitivity')
family_order <- c('KM', 'exponential', 'weibull', 'lognormal', 'loglogistic', 'coxph', 'aft', 'other')

# 🔴 Define: scalar, path, and schema helpers ===============================
## 🟠 Define: low-level utilities ===============================
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

normalize_existing_path <- function(path) {
  normalizePath(path, winslash = '/', mustWork = FALSE)
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

safe_integer <- function(x) {
  suppressWarnings(as.integer(as.character(x)))
}

safe_character <- function(x) {
  as.character(x)
}

safe_divide <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

weighted_mean_safe <- function(x, w) {
  x <- safe_numeric(x)
  w <- safe_numeric(w)
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) {
    return(NA_real_)
  }
  sum(x[keep] * w[keep]) / sum(w[keep])
}

combine_nonempty_text <- function(...) {
  vals_list <- list(...)
  if (length(vals_list) == 0L) {
    return(NA_character_)
  }
  
  lens <- vapply(vals_list, length, integer(1))
  target_n <- max(1L, lens)
  
  vals_list <- lapply(vals_list, function(x) {
    x_chr <- safe_character(x)
    
    if (length(x_chr) == 0L) {
      return(rep(NA_character_, target_n))
    }
    
    if (length(x_chr) == 1L && target_n > 1L) {
      return(rep(x_chr, target_n))
    }
    
    if (length(x_chr) == target_n) {
      return(x_chr)
    }
    
    stop('combine_nonempty_text inputs must have length 1 or a common vector length.', call. = FALSE)
  })
  
  out <- vapply(seq_len(target_n), function(i) {
    vals <- vapply(vals_list, function(x) x[[i]], character(1))
    vals <- unique(vals[!is.na(vals) & nzchar(trimws(vals))])
    if (length(vals) == 0L) {
      NA_character_
    } else {
      paste(vals, collapse = '; ')
    }
  }, character(1))
  
  if (target_n == 1L) {
    out[[1]]
  } else {
    out
  }
}

normalize_dataset_label <- function(x) {
  x_chr <- trimws(safe_character(x))
  x_up <- toupper(x_chr)
  out <- ifelse(
    x_up == 'PNU',
    'PNU',
    ifelse(x_up == 'SNU', 'SNU', ifelse(x_up == 'MERGED', 'merged', x_chr))
  )
  out[is.na(x)] <- NA_character_
  out
}

first_existing_name <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) {
    return(NULL)
  }
  hit[[1]]
}

coalesce_numeric_columns <- function(df, candidates, default = NA_real_) {
  out <- rep(default, nrow(df))
  for (nm in candidates) {
    if (!(nm %in% names(df))) {
      next
    }
    vals <- safe_numeric(df[[nm]])
    idx <- is.na(out) & !is.na(vals)
    out[idx] <- vals[idx]
  }
  out
}

coalesce_character_columns <- function(df, candidates, default = NA_character_) {
  out <- rep(default, nrow(df))
  for (nm in candidates) {
    if (!(nm %in% names(df))) {
      next
    }
    vals <- safe_character(df[[nm]])
    idx <- (is.na(out) | out == '') & !is.na(vals) & vals != ''
    out[idx] <- vals[idx]
  }
  out
}

read_csv_if_exists <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

read_rds_if_exists <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  readRDS(path)
}

assert_required_file <- function(path, label, required = TRUE, allow_partial = TRUE) {
  if (file.exists(path)) {
    return(invisible(path))
  }
  if (isTRUE(required) && !isTRUE(allow_partial)) {
    stop(sprintf('%s not found: %s', label, path), call. = FALSE)
  }
  warning(sprintf('%s not found and will be skipped: %s', label, path), call. = FALSE)
  invisible(NULL)
}

as_named_factor <- function(x, levels) {
  factor(safe_character(x), levels = levels)
}

# 🔴 Define: support and formula harmonizers ===============================
## 🟠 Define: Stage-1 lookup builders ===============================
canonicalize_support_label <- function(x) {
  raw <- tolower(gsub('[^a-z0-9]+', '_', trimws(safe_character(x))))
  dplyr::case_when(
    is.na(raw) ~ NA_character_,
    raw %in% c('primary_supported', 'primary', 'primary_supported_horizon', 'primarysupported') ~ 'primary_supported',
    raw %in% c('secondary', 'secondary_supported') ~ 'secondary',
    raw %in% c('sensitivity', 'limited_support', 'limited_supported', 'sensitive') ~ 'sensitivity',
    raw %in% c('projection', 'tail_projection', 'projected') ~ 'projection',
    TRUE ~ raw
  )
}

support_priority_from_label <- function(x) {
  lbl <- canonicalize_support_label(x)
  dplyr::case_when(
    lbl == 'primary_supported' ~ 1L,
    lbl == 'secondary' ~ 2L,
    lbl == 'sensitivity' ~ 3L,
    lbl == 'projection' ~ 4L,
    TRUE ~ 5L
  )
}

normalize_horizon_registry <- function(df) {
  req <- c('dataset', 'horizon_year')
  missing_req <- setdiff(req, names(df))
  if (length(missing_req) > 0L) {
    stop(sprintf('Stage 1 horizon registry is missing required columns: %s', paste(missing_req, collapse = ', ')), call. = FALSE)
  }
  out <- tibble::as_tibble(df) %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      horizon_year = safe_integer(horizon_year),
      support_label = canonicalize_support_label(coalesce_character_columns(., c('interpretation_tier', 'horizon_support_label'))),
      support_priority = support_priority_from_label(support_label),
      primary_supported_flag = if ('primary_supported_flag' %in% names(.)) as.logical(primary_supported_flag) else support_label == 'primary_supported',
      interpretation_note = coalesce_character_columns(., c('interpretation_note'))
    ) %>%
    select(dataset, horizon_year, support_label, support_priority, primary_supported_flag, interpretation_note) %>%
    distinct()
  
  expected_horizons <- sort(unique(out$horizon_year))
  if (!identical(expected_horizons, as.integer(common_horizons_year))) {
    stop('Stage 1 horizon registry must remain locked to horizons 1:10.', call. = FALSE)
  }
  out
}

normalize_formula_registry <- function(df) {
  req <- c('dataset', 'formula_name', 'formula_label')
  missing_req <- setdiff(req, names(df))
  if (length(missing_req) > 0L) {
    stop(sprintf('Stage 1 formula registry is missing required columns: %s', paste(missing_req, collapse = ', ')), call. = FALSE)
  }
  tibble::as_tibble(df) %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      formula_name = safe_character(formula_name),
      formula_label = safe_character(formula_label),
      site_branch = coalesce_character_columns(., c('site_branch'), default = 'site_free'),
      interaction_branch = coalesce_character_columns(., c('interaction_branch'), default = 'no_age_sex_interaction')
    ) %>%
    select(dataset, formula_name, formula_label, site_branch, interaction_branch) %>%
    distinct()
}

normalize_threshold_registry <- function(df) {
  threshold_col <- first_existing_name(df, c('threshold', 'risk_threshold', 'threshold_value'))
  if (is.null(threshold_col)) {
    stop('Stage 1 threshold registry must contain a threshold column.', call. = FALSE)
  }
  sort(unique(safe_numeric(df[[threshold_col]])))
}

compute_riskset_lookup <- function(analysis_datasets, horizons_year) {
  bind_rows(lapply(names(analysis_datasets), function(dataset_name) {
    df <- tibble::as_tibble(analysis_datasets[[dataset_name]]) %>%
      mutate(
        dataset = dataset_name,
        time_year = safe_numeric(time_year),
        event_main = safe_integer(event_main)
      )
    
    bind_rows(lapply(horizons_year, function(h) {
      event_by_h <- df$event_main == 1L & df$time_year <= h
      known_nonevent_by_h <- !event_by_h & df$time_year >= h
      tibble(
        dataset = dataset_name,
        horizon_year = as.integer(h),
        risk_set_size = sum(df$time_year >= h, na.rm = TRUE),
        event_count_by_horizon = sum(event_by_h, na.rm = TRUE),
        known_nonevent_by_horizon = sum(known_nonevent_by_h, na.rm = TRUE),
        max_followup_year = max(df$time_year, na.rm = TRUE)
      )
    }))
  }))
}

join_stage1_support <- function(df, support_lookup, riskset_lookup) {
  out <- df %>%
    left_join(support_lookup, by = c('dataset', 'horizon_year')) %>%
    left_join(riskset_lookup, by = c('dataset', 'horizon_year'))
  
  if (!'support_label' %in% names(out)) {
    out$support_label <- NA_character_
  }
  if (!'support_priority' %in% names(out)) {
    out$support_priority <- NA_integer_
  }
  if (!'primary_supported_flag' %in% names(out)) {
    out$primary_supported_flag <- NA
  }
  if (!'interpretation_note' %in% names(out)) {
    out$interpretation_note <- NA_character_
  }
  
  out %>%
    mutate(
      instability_marker = dplyr::case_when(
        support_label == 'projection' ~ 'late_tail_projection',
        support_label == 'sensitivity' ~ 'limited_support_window',
        support_label == 'secondary' ~ 'secondary_supported_window',
        support_label == 'primary_supported' ~ 'primary_supported_window',
        TRUE ~ 'unclassified_horizon'
      )
    )
}

fill_formula_metadata <- function(df, formula_lookup) {
  out <- df
  if (!'formula_name' %in% names(out)) {
    out$formula_name <- NA_character_
  }
  if (!'site_branch' %in% names(out)) {
    out$site_branch <- NA_character_
  }
  if (!'interaction_branch' %in% names(out)) {
    out$interaction_branch <- NA_character_
  }
  if (!'formula_label' %in% names(out)) {
    out$formula_label <- NA_character_
  }
  
  specific <- out %>%
    filter(!is.na(formula_name) & formula_name != '' & toupper(formula_name) != 'ALL' & tolower(formula_name) != 'reference') %>%
    left_join(
      formula_lookup,
      by = c('dataset', 'formula_name'),
      suffix = c('', '.lookup')
    ) %>%
    mutate(
      formula_label = dplyr::coalesce(formula_label, formula_label.lookup),
      site_branch = dplyr::coalesce(site_branch, site_branch.lookup),
      interaction_branch = dplyr::coalesce(interaction_branch, interaction_branch.lookup)
    ) %>%
    select(-any_of(c('formula_label.lookup', 'site_branch.lookup', 'interaction_branch.lookup')))
  
  generic <- out %>%
    filter(is.na(formula_name) | formula_name == '' | toupper(formula_name) == 'ALL' | tolower(formula_name) == 'reference')
  
  expanded_generic <- tibble()
  if (nrow(generic) > 0L) {
    with_site <- generic %>% filter(!is.na(site_branch) & site_branch != '')
    no_site <- generic %>% filter(is.na(site_branch) | site_branch == '')
    
    if (nrow(with_site) > 0L) {
      expanded_generic <- bind_rows(
        expanded_generic,
        with_site %>%
          select(-any_of(c('formula_name', 'formula_label', 'interaction_branch'))) %>%
          inner_join(formula_lookup, by = c('dataset', 'site_branch'))
      )
    }
    
    if (nrow(no_site) > 0L) {
      expanded_generic <- bind_rows(
        expanded_generic,
        no_site %>%
          select(-any_of(c('formula_name', 'formula_label', 'site_branch', 'interaction_branch'))) %>%
          inner_join(formula_lookup, by = c('dataset'))
      )
    }
  }
  
  bind_rows(specific, expanded_generic) %>%
    mutate(
      formula_name = safe_character(formula_name),
      formula_label = dplyr::coalesce(formula_label, formula_name)
    )
}

make_model_class_label <- function(x) {
  dplyr::case_when(
    x == 'km_benchmark' ~ 'KM benchmark',
    x == 'no_cure' ~ 'No-cure',
    x == 'frequentist_cure' ~ 'Frequentist cure',
    x == 'bayesian_cure' ~ 'Bayesian cure',
    x == 'remission_sensitivity' ~ 'Remission sensitivity',
    TRUE ~ safe_character(x)
  )
}

make_model_display_label <- function(model_class, family, benchmark_label = NA_character_) {
  dplyr::case_when(
    model_class == 'km_benchmark' & !is.na(benchmark_label) & nzchar(benchmark_label) ~ paste0('KM: ', benchmark_label),
    model_class == 'km_benchmark' ~ 'KM benchmark',
    model_class == 'no_cure' ~ paste0('No-cure: ', family),
    model_class == 'frequentist_cure' ~ paste0('Cure MLE: ', family),
    model_class == 'bayesian_cure' ~ paste0('Bayesian cure: ', family),
    TRUE ~ paste(model_class, family)
  )
}

make_reporting_status <- function(x, fallback_reported = FALSE) {
  x_chr <- tolower(trimws(safe_character(x)))
  dplyr::case_when(
    x_chr %in% c('reported', 'ok', 'success') ~ 'reported',
    grepl('not_estimable', x_chr, fixed = TRUE) ~ 'suppressed',
    grepl('suppress', x_chr, fixed = TRUE) ~ 'suppressed',
    grepl('unstable', x_chr, fixed = TRUE) ~ 'reported_with_warning',
    grepl('warning', x_chr, fixed = TRUE) ~ 'reported_with_warning',
    is.na(x_chr) & fallback_reported ~ 'reported',
    TRUE ~ x_chr
  )
}

# 🔴 Define: stage-output readers ===============================
## 🟠 Define: staged file loaders ===============================
load_stage3_outputs <- function(export_dir, required = TRUE) {
  files <- list(
    benchmark_registry = file.path(export_dir, 'stage3_km_benchmark_registry.csv'),
    group_risk = file.path(export_dir, 'stage3_km_group_risk_table.csv'),
    classification = file.path(export_dir, 'stage3_km_benchmark_classification.csv')
  )
  
  missing_any <- any(!file.exists(unlist(files)))
  if (missing_any) {
    for (nm in names(files)) {
      assert_required_file(files[[nm]], sprintf('Stage 3 %s', nm), required = required, allow_partial = allow_partial_stage10_build)
    }
  }
  
  list(
    benchmark_registry = read_csv_if_exists(files$benchmark_registry),
    group_risk = read_csv_if_exists(files$group_risk),
    classification = read_csv_if_exists(files$classification)
  )
}

load_stage5_outputs <- function(export_dir, required = TRUE) {
  files <- list(
    model_registry = file.path(export_dir, 'stage5_model_registry.csv'),
    subject_horizon_risk_long = file.path(export_dir, 'stage5_subject_horizon_risk_long.csv'),
    model_performance_long = file.path(export_dir, 'stage5_model_performance_long.csv'),
    calibration_diagnostics_long = file.path(export_dir, 'stage5_calibration_diagnostics_long.csv')
  )
  
  missing_any <- any(!file.exists(unlist(files)))
  if (missing_any) {
    for (nm in names(files)) {
      assert_required_file(files[[nm]], sprintf('Stage 5 %s', nm), required = required, allow_partial = allow_partial_stage10_build)
    }
  }
  
  list(
    model_registry = read_csv_if_exists(files$model_registry),
    subject_horizon_risk_long = read_csv_if_exists(files$subject_horizon_risk_long),
    model_performance_long = read_csv_if_exists(files$model_performance_long),
    calibration_diagnostics_long = read_csv_if_exists(files$calibration_diagnostics_long)
  )
}

load_stage7_outputs <- function(export_dir, required = FALSE) {
  files <- list(
    fit_registry = file.path(export_dir, 'stage7_fit_registry.csv'),
    risk_summary = file.path(export_dir, 'stage7_risk_summary.csv'),
    delta_risk = file.path(export_dir, 'stage7_delta_risk.csv'),
    threshold_metrics = file.path(export_dir, 'stage7_threshold_metrics.csv')
  )
  
  missing_any <- any(!file.exists(unlist(files)))
  if (missing_any) {
    for (nm in names(files)) {
      assert_required_file(files[[nm]], sprintf('Stage 7 %s', nm), required = required, allow_partial = allow_partial_stage10_build)
    }
  }
  
  list(
    fit_registry = read_csv_if_exists(files$fit_registry),
    risk_summary = read_csv_if_exists(files$risk_summary),
    delta_risk = read_csv_if_exists(files$delta_risk),
    threshold_metrics = read_csv_if_exists(files$threshold_metrics)
  )
}

load_stage8_outputs <- function(export_dir, required = FALSE) {
  files <- list(
    model_registry = file.path(export_dir, 'bayes_stage8_model_registry.csv'),
    posterior_cohort_yearly = file.path(export_dir, 'bayes_stage8_posterior_cohort_yearly.csv'),
    posterior_classification = file.path(export_dir, 'bayes_stage8_posterior_classification.csv'),
    posterior_delta_vs_nocure = file.path(export_dir, 'bayes_stage8_posterior_delta_vs_nocure.csv'),
    carryforward_metadata = file.path(export_dir, 'bayes_stage8_carryforward_metadata.csv')
  )
  
  missing_any <- any(!file.exists(unlist(files)))
  if (missing_any) {
    for (nm in names(files)) {
      assert_required_file(files[[nm]], sprintf('Stage 8 %s', nm), required = required, allow_partial = allow_partial_stage10_build)
    }
  }
  
  list(
    model_registry = read_csv_if_exists(files$model_registry),
    posterior_cohort_yearly = read_csv_if_exists(files$posterior_cohort_yearly),
    posterior_classification = read_csv_if_exists(files$posterior_classification),
    posterior_delta_vs_nocure = read_csv_if_exists(files$posterior_delta_vs_nocure),
    carryforward_metadata = read_csv_if_exists(files$carryforward_metadata)
  )
}

load_stage6_outputs_optional <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext == 'csv') {
    return(read_csv_if_exists(path))
  }
  obj <- read_rds_if_exists(path)
  if (is.null(obj)) {
    return(NULL)
  }
  if (inherits(obj, 'data.frame')) {
    return(tibble::as_tibble(obj))
  }
  if (is.list(obj) && !is.null(obj$outputs) && !is.null(obj$outputs$carry_forward_flag_table)) {
    return(tibble::as_tibble(obj$outputs$carry_forward_flag_table))
  }
  NULL
}

# 🔴 Define: stage-specific standardizers ===============================
## 🟠 Define: KM benchmark normalizers ===============================
standardize_stage3_risk <- function(group_risk_tbl, formula_lookup, support_lookup, riskset_lookup) {
  if (is.null(group_risk_tbl) || nrow(group_risk_tbl) == 0L) {
    return(tibble())
  }
  
  out <- tibble::as_tibble(group_risk_tbl) %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      horizon_year = safe_integer(horizon_year),
      n_subjects_in_group = safe_numeric(n_subjects_in_group),
      km_risk = safe_numeric(km_risk),
      km_survival = safe_numeric(km_survival),
      benchmark_scope = safe_character(benchmark_scope),
      site_branch = safe_character(site_branch)
    )
  
  if (!isTRUE(include_stage3_sensitivity_benchmarks) && 'benchmark_scope' %in% names(out)) {
    out <- out %>% filter(.data$benchmark_scope %in% c('primary', 'primary_structural'))
  }
  
  out %>%
    group_by(dataset, benchmark_id, benchmark_label, benchmark_scope, site_branch, subgroup_kind, subgroup_variable, horizon_year) %>%
    summarise(
      risk = weighted_mean_safe(km_risk, n_subjects_in_group),
      survival = weighted_mean_safe(km_survival, n_subjects_in_group),
      n_subjects = sum(n_subjects_in_group, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      source_stage = 'Stage 3',
      model_class = 'km_benchmark',
      model_block = 'benchmark',
      family = 'KM',
      model_id = safe_character(benchmark_id),
      benchmark_id = safe_character(benchmark_id),
      formula_name = NA_character_,
      formula_label = NA_character_,
      interaction_branch = NA_character_,
      risk_lower = NA_real_,
      risk_upper = NA_real_,
      survival_lower = NA_real_,
      survival_upper = NA_real_,
      estimate_type = 'assigned_group_km',
      reporting_status = 'reported',
      reporting_note = benchmark_label
    ) %>%
    fill_formula_metadata(formula_lookup) %>%
    join_stage1_support(support_lookup, riskset_lookup) %>%
    mutate(
      model_class_label = make_model_class_label(model_class),
      model_display_label = make_model_display_label(model_class, family, benchmark_label = benchmark_label)
    )
}

standardize_stage3_threshold <- function(classification_tbl, formula_lookup, support_lookup, riskset_lookup) {
  if (is.null(classification_tbl) || nrow(classification_tbl) == 0L) {
    return(tibble())
  }
  
  out <- tibble::as_tibble(classification_tbl) %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      horizon_year = safe_integer(horizon_year),
      threshold = safe_numeric(threshold),
      benchmark_scope = safe_character(benchmark_scope),
      site_branch = safe_character(site_branch)
    )
  
  if (!isTRUE(include_stage3_sensitivity_benchmarks) && 'benchmark_scope' %in% names(out)) {
    out <- out %>% filter(.data$benchmark_scope %in% c('primary', 'primary_structural'))
  }
  
  out %>%
    transmute(
      dataset = dataset,
      source_stage = 'Stage 3',
      model_class = 'km_benchmark',
      model_block = 'benchmark',
      family = 'KM',
      model_id = safe_character(benchmark_id),
      benchmark_id = safe_character(benchmark_id),
      benchmark_label = safe_character(benchmark_label),
      formula_name = NA_character_,
      formula_label = NA_character_,
      site_branch = site_branch,
      interaction_branch = NA_character_,
      subgroup_kind = safe_character(subgroup_kind),
      subgroup_variable = safe_character(subgroup_variable),
      horizon_year = horizon_year,
      threshold = threshold,
      positive_classification_rate = safe_numeric(positive_classification_rate),
      false_positive_rate = coalesce_numeric_columns(., c('false_positive_rate', 'false_positive_burden_known_nonevent', 'false_positive_burden_nonevents')),
      false_positive_rate_lower = NA_real_,
      false_positive_rate_upper = NA_real_,
      false_positive_burden = coalesce_numeric_columns(., c('false_positive_burden_all', 'false_positive_burden')),
      false_positive_burden_lower = NA_real_,
      false_positive_burden_upper = NA_real_,
      false_positive_per_100 = coalesce_numeric_columns(., c('false_positive_per_100')),
      ppv = safe_numeric(ppv),
      tpr = coalesce_numeric_columns(., c('tpr', 'sensitivity')),
      net_benefit = NA_real_,
      net_benefit_lower = NA_real_,
      net_benefit_upper = NA_real_,
      net_reduction_unnecessary_per_100 = NA_real_,
      denom_case = safe_numeric(n_event_by_horizon),
      denom_control = safe_numeric(n_known_nonevent_by_horizon),
      reporting_status = ifelse(safe_numeric(n_evaluable_by_horizon) > 0, 'reported', 'suppressed'),
      reporting_note = safe_character(evaluation_rule)
    ) %>%
    fill_formula_metadata(formula_lookup) %>%
    join_stage1_support(support_lookup, riskset_lookup) %>%
    mutate(
      model_class_label = make_model_class_label(model_class),
      model_display_label = make_model_display_label(model_class, family, benchmark_label = benchmark_label)
    )
}

## 🟠 Define: no-cure normalizers ===============================
standardize_stage5_risk <- function(subject_risk_tbl, performance_tbl, formula_lookup, support_lookup, riskset_lookup) {
  if (!is.null(subject_risk_tbl) && nrow(subject_risk_tbl) > 0L) {
    out <- tibble::as_tibble(subject_risk_tbl) %>%
      mutate(
        dataset = normalize_dataset_label(dataset),
        horizon_year = safe_integer(horizon_year),
        risk_pred = safe_numeric(risk_pred),
        survival_pred = safe_numeric(survival_pred)
      ) %>%
      filter(model_class == 'non_cure') %>%
      group_by(dataset, model_id, model_family, formula_name, formula_label, site_branch, interaction_branch, horizon_year) %>%
      summarise(
        risk = mean(risk_pred, na.rm = TRUE),
        survival = mean(survival_pred, na.rm = TRUE),
        n_subjects = dplyr::n_distinct(unique_person_id),
        .groups = 'drop'
      ) %>%
      mutate(
        source_stage = 'Stage 5',
        model_class = 'no_cure',
        model_block = 'non_cure',
        family = safe_character(model_family),
        benchmark_id = NA_character_,
        risk_lower = NA_real_,
        risk_upper = NA_real_,
        survival_lower = NA_real_,
        survival_upper = NA_real_,
        estimate_type = 'mean_subject_predicted_risk',
        reporting_status = ifelse(is.na(risk), 'suppressed', 'reported'),
        reporting_note = 'Stage 5 mean subject-level predicted risk'
      )
    
    return(
      out %>%
        fill_formula_metadata(formula_lookup) %>%
        join_stage1_support(support_lookup, riskset_lookup) %>%
        mutate(
          model_class_label = make_model_class_label(model_class),
          model_display_label = make_model_display_label(model_class, family)
        )
    )
  }
  
  if (is.null(performance_tbl) || nrow(performance_tbl) == 0L) {
    return(tibble())
  }
  
  perf_wide <- tibble::as_tibble(performance_tbl) %>%
    filter(model_class == 'non_cure', metric_domain == 'horizon_summary', metric_name == 'mean_predicted_risk') %>%
    transmute(
      dataset = normalize_dataset_label(dataset),
      model_id = safe_character(model_id),
      family = safe_character(model_family),
      formula_name = safe_character(formula_name),
      formula_label = coalesce_character_columns(., c('formula_label'), default = NA_character_),
      site_branch = safe_character(site_branch),
      interaction_branch = safe_character(interaction_branch),
      horizon_year = safe_integer(horizon_year),
      risk = safe_numeric(metric_value)
    ) %>%
    mutate(
      source_stage = 'Stage 5',
      model_class = 'no_cure',
      model_block = 'non_cure',
      benchmark_id = NA_character_,
      survival = ifelse(!is.na(risk), 1 - risk, NA_real_),
      risk_lower = NA_real_,
      risk_upper = NA_real_,
      survival_lower = NA_real_,
      survival_upper = NA_real_,
      estimate_type = 'stage5_horizon_summary_mean_risk',
      reporting_status = ifelse(is.na(risk), 'suppressed', 'reported'),
      reporting_note = 'Fallback from Stage 5 horizon summary'
    )
  
  perf_wide %>%
    fill_formula_metadata(formula_lookup) %>%
    join_stage1_support(support_lookup, riskset_lookup) %>%
    mutate(
      model_class_label = make_model_class_label(model_class),
      model_display_label = make_model_display_label(model_class, family)
    )
}

standardize_stage5_threshold <- function(performance_tbl, formula_lookup, support_lookup, riskset_lookup) {
  if (is.null(performance_tbl) || nrow(performance_tbl) == 0L) {
    return(tibble())
  }
  
  metric_map <- c(
    positive_classification_rate = 'positive_classification_rate',
    false_positive_burden_nonevents = 'false_positive_rate',
    false_positive_burden_all = 'false_positive_burden',
    false_positive_per_100 = 'false_positive_per_100',
    ppv = 'ppv',
    sensitivity = 'tpr',
    net_benefit = 'net_benefit',
    net_reduction_unnecessary_per_100 = 'net_reduction_unnecessary_per_100',
    case_count_ipcw = 'denom_case',
    nonevent_count_ipcw = 'denom_control',
    binary_outcome_support_flag = 'binary_outcome_support_flag'
  )
  
  out <- tibble::as_tibble(performance_tbl) %>%
    filter(model_class == 'non_cure', metric_domain == 'threshold_summary', metric_name %in% names(metric_map)) %>%
    mutate(metric_name_std = unname(metric_map[metric_name])) %>%
    transmute(
      dataset = normalize_dataset_label(dataset),
      model_id = safe_character(model_id),
      family = safe_character(model_family),
      formula_name = safe_character(formula_name),
      formula_label = coalesce_character_columns(., c('formula_label'), default = NA_character_),
      site_branch = safe_character(site_branch),
      interaction_branch = safe_character(interaction_branch),
      horizon_year = safe_integer(horizon_year),
      threshold = safe_numeric(threshold),
      metric_name_std = metric_name_std,
      metric_value = safe_numeric(metric_value)
    ) %>%
    tidyr::pivot_wider(names_from = metric_name_std, values_from = metric_value) %>%
    mutate(
      source_stage = 'Stage 5',
      model_class = 'no_cure',
      model_block = 'non_cure',
      benchmark_id = NA_character_,
      false_positive_rate_lower = NA_real_,
      false_positive_rate_upper = NA_real_,
      false_positive_burden_lower = NA_real_,
      false_positive_burden_upper = NA_real_,
      net_benefit_lower = NA_real_,
      net_benefit_upper = NA_real_,
      reporting_status = ifelse(binary_outcome_support_flag == 1, 'reported', 'suppressed'),
      reporting_note = ifelse(binary_outcome_support_flag == 1, 'Stage 5 IPCW fixed-horizon threshold summary', 'Binary horizon outcome not supported at this horizon')
    ) %>%
    fill_formula_metadata(formula_lookup) %>%
    join_stage1_support(support_lookup, riskset_lookup) %>%
    mutate(
      model_class_label = make_model_class_label(model_class),
      model_display_label = make_model_display_label(model_class, family)
    )
  
  out
}

standardize_stage5_auxiliary <- function(performance_tbl, formula_lookup, support_lookup, riskset_lookup) {
  if (is.null(performance_tbl) || nrow(performance_tbl) == 0L) {
    return(tibble())
  }
  
  metric_keep <- c(
    'time_dependent_auc',
    'brier_score',
    'scaled_brier_score',
    'mean_calibration_difference',
    'calibration_intercept_offset',
    'calibration_intercept_free',
    'calibration_slope',
    'harrell_c',
    'uno_c',
    'restricted_integrated_brier_score'
  )
  
  tibble::as_tibble(performance_tbl) %>%
    filter(model_class == 'non_cure', metric_name %in% metric_keep) %>%
    transmute(
      dataset = normalize_dataset_label(dataset),
      source_stage = 'Stage 5',
      model_class = 'no_cure',
      model_block = 'non_cure',
      family = safe_character(model_family),
      model_id = safe_character(model_id),
      benchmark_id = NA_character_,
      formula_name = safe_character(formula_name),
      formula_label = coalesce_character_columns(., c('formula_label'), default = NA_character_),
      site_branch = safe_character(site_branch),
      interaction_branch = safe_character(interaction_branch),
      horizon_year = ifelse(is.na(safe_integer(horizon_year)), safe_integer(anchor_horizon_year), safe_integer(horizon_year)),
      threshold = NA_real_,
      metric_domain = safe_character(metric_domain),
      metric_name = safe_character(metric_name),
      metric_value = safe_numeric(metric_value),
      metric_lower = NA_real_,
      metric_upper = NA_real_,
      reporting_status = 'reported',
      reporting_note = paste('Stage 5', safe_character(metric_domain))
    ) %>%
    fill_formula_metadata(formula_lookup) %>%
    join_stage1_support(support_lookup, riskset_lookup) %>%
    mutate(
      model_class_label = make_model_class_label(model_class),
      model_display_label = make_model_display_label(model_class, family)
    )
}

## 🟠 Define: Stage 7 fit-registry compatibility helpers ===============================
stage7_status_is_success <- function(x) {
  x_chr <- tolower(trimws(safe_character(x)))
  startsWith(x_chr, "ok")
}

build_stage7_fit_info <- function(fit_registry_tbl) {
  if (is.null(fit_registry_tbl) || nrow(fit_registry_tbl) == 0L) {
    return(tibble())
  }
  
  tbl <- tibble::as_tibble(fit_registry_tbl)
  
  fit_status <- coalesce_character_columns(tbl, c("fit_status"))
  prediction_status <- coalesce_character_columns(tbl, c("prediction_status"))
  overall_status <- coalesce_character_columns(tbl, c("overall_status"))
  
  fit_success <- if ("fit_component_success" %in% names(tbl)) {
    as.logical(tbl$fit_component_success)
  } else {
    stage7_status_is_success(fit_status)
  }
  
  prediction_success <- if ("prediction_component_success" %in% names(tbl)) {
    as.logical(tbl$prediction_component_success)
  } else {
    stage7_status_is_success(prediction_status)
  }
  
  overall_success <- if ("overall_success" %in% names(tbl)) {
    as.logical(tbl$overall_success)
  } else if ("overall_status" %in% names(tbl)) {
    stage7_status_is_success(overall_status)
  } else {
    fit_success & prediction_success
  }
  
  failure_type <- if ("failure_type" %in% names(tbl)) {
    safe_character(tbl$failure_type)
  } else {
    dplyr::case_when(
      !dplyr::coalesce(fit_success, FALSE) & !dplyr::coalesce(prediction_success, FALSE) ~ "fit_and_prediction_failed",
      !dplyr::coalesce(fit_success, FALSE) ~ "fit_failed",
      dplyr::coalesce(fit_success, FALSE) & !dplyr::coalesce(prediction_success, FALSE) ~ "prediction_failed",
      !dplyr::coalesce(overall_success, FALSE) ~ "unknown_failure",
      TRUE ~ NA_character_
    )
  }
  
  tibble(
    dataset = normalize_dataset_label(tbl$dataset),
    formula_name = coalesce_character_columns(tbl, c("formula_variant", "formula_name")),
    model_id = safe_character(tbl$model_id),
    overall_success = dplyr::coalesce(as.logical(overall_success), TRUE),
    failure_type = failure_type,
    stage6_flag = coalesce_character_columns(tbl, c("stage6__cure_model_eligibility_flag", "stage6__final_decision_flag")),
    stage6_note = coalesce_character_columns(tbl, c("stage6__screening_note"))
  ) %>%
    distinct()
}

## 🟠 Define: frequentist cure normalizers ===============================
standardize_stage7_risk <- function(risk_tbl, fit_registry_tbl, formula_lookup, support_lookup, riskset_lookup) {
  if (is.null(risk_tbl) || nrow(risk_tbl) == 0L) {
    return(tibble())
  }
  
  fit_info <- build_stage7_fit_info(fit_registry_tbl)
  
  tibble::as_tibble(risk_tbl) %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      formula_name = coalesce_character_columns(., c("formula_variant", "formula_name")),
      formula_label = coalesce_character_columns(., c("formula_label"), default = NA_character_),
      horizon_year = safe_integer(horizon_year),
      site_branch = coalesce_character_columns(., c("site_branch"), default = NA_character_),
      interaction_branch = coalesce_character_columns(., c("interaction_branch"), default = NA_character_),
      family = coalesce_character_columns(., c("latency_type", "model_family"), default = "other"),
      risk = coalesce_numeric_columns(., c("mean_risk_overall", "overall_risk", "meanRisk", "mean_risk", "risk_mean")),
      risk_lower = coalesce_numeric_columns(., c("mean_risk_ci_lower", "risk_ci_lower", "ci_lower")),
      risk_upper = coalesce_numeric_columns(., c("mean_risk_ci_upper", "risk_ci_upper", "ci_upper")),
      survival = coalesce_numeric_columns(., c("mean_survival_overall", "overall_survival", "meanSurvival", "mean_survival", "survival_mean"))
    ) %>%
    mutate(
      survival = ifelse(is.na(survival) & !is.na(risk), 1 - risk, survival),
      survival_lower = ifelse(!is.na(risk_upper), 1 - risk_upper, NA_real_),
      survival_upper = ifelse(!is.na(risk_lower), 1 - risk_lower, NA_real_)
    ) %>%
    transmute(
      dataset = dataset,
      source_stage = "Stage 7",
      model_class = "frequentist_cure",
      model_block = coalesce_character_columns(., c("model_block"), default = "cure_mle"),
      family = family,
      model_id = safe_character(model_id),
      benchmark_id = NA_character_,
      formula_name = formula_name,
      formula_label = formula_label,
      site_branch = site_branch,
      interaction_branch = interaction_branch,
      latency_type = coalesce_character_columns(., c("latency_type")),
      horizon_year = horizon_year,
      risk = risk,
      risk_lower = risk_lower,
      risk_upper = risk_upper,
      survival = survival,
      survival_lower = survival_lower,
      survival_upper = survival_upper,
      estimate_type = "bootstrap_summary_mean_overall_risk",
      reporting_status = make_reporting_status(coalesce_character_columns(., c("cohort_reporting_status", "reporting_status")), fallback_reported = TRUE),
      reporting_note = coalesce_character_columns(., c("cohort_reporting_note", "reporting_note")),
      admissible_flag = NA,
      admissibility_reasons = NA_character_
    ) %>%
    left_join(fit_info, by = c("dataset", "formula_name", "model_id")) %>%
    mutate(
      overall_success = dplyr::coalesce(overall_success, TRUE),
      reporting_status = ifelse(!overall_success, "suppressed", reporting_status),
      reporting_note = combine_nonempty_text(
        reporting_note,
        ifelse(!overall_success, failure_type, NA_character_),
        stage6_note
      )
    ) %>%
    fill_formula_metadata(formula_lookup) %>%
    join_stage1_support(support_lookup, riskset_lookup) %>%
    mutate(
      screening_flag = stage6_flag,
      screening_note = stage6_note,
      model_class_label = make_model_class_label(model_class),
      model_display_label = make_model_display_label(model_class, family)
    ) %>%
    select(-any_of(c("overall_success", "failure_type", "stage6_flag", "stage6_note")))
}

standardize_stage7_threshold <- function(threshold_tbl, fit_registry_tbl, formula_lookup, support_lookup, riskset_lookup) {
  if (is.null(threshold_tbl) || nrow(threshold_tbl) == 0L) {
    return(tibble())
  }
  
  fit_info <- build_stage7_fit_info(fit_registry_tbl)
  
  tibble::as_tibble(threshold_tbl) %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      formula_name = coalesce_character_columns(., c("formula_variant", "formula_name")),
      formula_label = coalesce_character_columns(., c("formula_label"), default = NA_character_),
      horizon_year = safe_integer(horizon_year),
      threshold = safe_numeric(threshold),
      site_branch = coalesce_character_columns(., c("site_branch"), default = NA_character_),
      interaction_branch = coalesce_character_columns(., c("interaction_branch"), default = NA_character_),
      family = coalesce_character_columns(., c("latency_type", "model_family"), default = "other")
    ) %>%
    transmute(
      dataset = dataset,
      source_stage = "Stage 7",
      model_class = "frequentist_cure",
      model_block = coalesce_character_columns(., c("model_block"), default = "cure_mle"),
      family = family,
      model_id = safe_character(model_id),
      benchmark_id = NA_character_,
      formula_name = formula_name,
      formula_label = formula_label,
      site_branch = site_branch,
      interaction_branch = interaction_branch,
      latency_type = coalesce_character_columns(., c("latency_type")),
      horizon_year = horizon_year,
      threshold = threshold,
      positive_classification_rate = coalesce_numeric_columns(., c("positive_classification_rate")),
      false_positive_rate = {
        fp_rate <- coalesce_numeric_columns(., c("false_positive_rate", "FPR", "fpr", "false_positive_burden_primary", "false_positive_burden_nonevents"))
        ifelse(is.na(fp_rate) & "specificity" %in% names(.), 1 - safe_numeric(.$specificity), fp_rate)
      },
      false_positive_rate_lower = coalesce_numeric_columns(., c("false_positive_rate_ci_lower", "fpr_ci_lower", "false_positive_burden_primary_ci_lower", "false_positive_burden_ci_lower")),
      false_positive_rate_upper = coalesce_numeric_columns(., c("false_positive_rate_ci_upper", "fpr_ci_upper", "false_positive_burden_primary_ci_upper", "false_positive_burden_ci_upper")),
      false_positive_burden = coalesce_numeric_columns(., c("false_positive_burden_all", "false_positive_burden")),
      false_positive_burden_lower = coalesce_numeric_columns(., c("false_positive_burden_ci_lower", "false_positive_burden_all_ci_lower")),
      false_positive_burden_upper = coalesce_numeric_columns(., c("false_positive_burden_ci_upper", "false_positive_burden_all_ci_upper")),
      false_positive_per_100 = {
        fp100 <- coalesce_numeric_columns(., c("false_positive_per_100", "FP100"))
        ifelse(is.na(fp100) & "false_positive_weighted_per_n" %in% names(.), 100 * safe_numeric(.$false_positive_weighted_per_n), fp100)
      },
      ppv = coalesce_numeric_columns(., c("ppv", "PPV")),
      tpr = coalesce_numeric_columns(., c("tpr", "TPR", "sensitivity")),
      net_benefit = coalesce_numeric_columns(., c("net_benefit", "NB")),
      net_benefit_lower = coalesce_numeric_columns(., c("net_benefit_ci_lower", "NB_q025")),
      net_benefit_upper = coalesce_numeric_columns(., c("net_benefit_ci_upper", "NB_q975")),
      net_reduction_unnecessary_per_100 = coalesce_numeric_columns(., c("net_reduction_unnecessary_per_100", "unnecessary_high_risk_per_100_non_event", "unnecessary_high_risk_per_100_population")),
      denom_case = coalesce_numeric_columns(., c("case_count_ipcw", "denom_case", "n_event_by_horizon")),
      denom_control = coalesce_numeric_columns(., c("nonevent_count_ipcw", "denom_control", "n_known_nonevent_by_horizon")),
      reporting_status = make_reporting_status(coalesce_character_columns(., c("threshold_reporting_status", "classification_reporting_status", "reporting_status")), fallback_reported = FALSE),
      reporting_note = coalesce_character_columns(., c("threshold_reporting_note", "classification_reporting_note", "reporting_note")),
      threshold_projection_suppressed_flag = if ("threshold_projection_suppressed_flag" %in% names(.)) as.logical(.$threshold_projection_suppressed_flag) else FALSE
    ) %>%
    left_join(fit_info, by = c("dataset", "formula_name", "model_id")) %>%
    mutate(
      overall_success = dplyr::coalesce(overall_success, TRUE),
      threshold_projection_suppressed_flag = dplyr::coalesce(threshold_projection_suppressed_flag, FALSE),
      reporting_status = dplyr::case_when(
        threshold_projection_suppressed_flag ~ "suppressed",
        !overall_success ~ "suppressed",
        is.na(reporting_status) & !is.na(false_positive_rate) ~ "reported",
        is.na(reporting_status) ~ "suppressed",
        TRUE ~ reporting_status
      ),
      reporting_note = combine_nonempty_text(
        reporting_note,
        ifelse(threshold_projection_suppressed_flag, "Stage 7 threshold projection suppression", NA_character_),
        ifelse(!overall_success, failure_type, NA_character_),
        stage6_note
      )
    ) %>%
    fill_formula_metadata(formula_lookup) %>%
    join_stage1_support(support_lookup, riskset_lookup) %>%
    mutate(
      screening_flag = stage6_flag,
      screening_note = stage6_note,
      model_class_label = make_model_class_label(model_class),
      model_display_label = make_model_display_label(model_class, family)
    ) %>%
    select(-any_of(c("overall_success", "failure_type", "stage6_flag", "stage6_note")))
}

standardize_stage7_delta <- function(delta_tbl, formula_lookup, support_lookup, riskset_lookup) {
  if (is.null(delta_tbl) || nrow(delta_tbl) == 0L) {
    return(tibble())
  }
  
  tibble::as_tibble(delta_tbl) %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      formula_name = coalesce_character_columns(., c('formula_variant', 'formula_name')),
      formula_label = coalesce_character_columns(., c('formula_label'), default = NA_character_),
      horizon_year = safe_integer(horizon_year),
      family = coalesce_character_columns(., c('latency_type', 'model_family'), default = 'other')
    ) %>%
    transmute(
      dataset = dataset,
      source_stage = 'Stage 7',
      model_class = 'frequentist_cure',
      model_block = coalesce_character_columns(., c('model_block'), default = 'cure_mle'),
      family = family,
      model_id = safe_character(model_id),
      no_cure_model_id = coalesce_character_columns(., c('matched_nocure_model_id', 'no_cure_model_id')),
      formula_name = formula_name,
      formula_label = formula_label,
      site_branch = coalesce_character_columns(., c('site_branch'), default = NA_character_),
      interaction_branch = coalesce_character_columns(., c('interaction_branch'), default = NA_character_),
      horizon_year = horizon_year,
      threshold = NA_real_,
      metric_domain = 'delta_vs_nocure',
      metric_name = 'delta_risk_nc_minus_cure',
      metric_value = coalesce_numeric_columns(., c('mean_delta_risk_nc_minus_cure', 'delta_risk_nc_minus_cure', 'mean_delta_risk')),
      metric_lower = coalesce_numeric_columns(., c('mean_delta_risk_ci_lower', 'delta_risk_ci_lower')),
      metric_upper = coalesce_numeric_columns(., c('mean_delta_risk_ci_upper', 'delta_risk_ci_upper')),
      reporting_status = 'reported',
      reporting_note = 'Stage 7 cure versus no-cure risk shift'
    ) %>%
    fill_formula_metadata(formula_lookup) %>%
    join_stage1_support(support_lookup, riskset_lookup) %>%
    mutate(
      model_class_label = make_model_class_label(model_class),
      model_display_label = make_model_display_label(model_class, family)
    )
}

## 🟠 Define: Bayesian cure normalizers ===============================
standardize_stage8_risk <- function(cohort_tbl, formula_lookup, support_lookup, riskset_lookup) {
  if (is.null(cohort_tbl) || nrow(cohort_tbl) == 0L) {
    return(tibble())
  }
  
  out <- tibble::as_tibble(cohort_tbl) %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      formula_name = coalesce_character_columns(., c('formula_anchor', 'formula_name')),
      formula_label = coalesce_character_columns(., c('formula_label'), default = NA_character_),
      horizon_year = safe_integer(horizon_year),
      admissible_flag = if ('admissible_flag' %in% names(.)) as.logical(admissible_flag) else TRUE,
      screening_flag = coalesce_character_columns(., c('screening_flag')),
      screening_note = coalesce_character_columns(., c('screening_detail', 'cohort_reporting_note'))
    )
  
  if (isTRUE(keep_only_stage8_admissible_fits)) {
    out <- out %>% filter(dplyr::coalesce(admissible_flag, FALSE))
  }
  if (!isTRUE(include_stage8_unsupportive_sensitivity_fits)) {
    out <- out %>% filter(is.na(screening_flag) | screening_flag != 'unsupportive')
  }
  
  out %>%
    transmute(
      dataset = dataset,
      source_stage = 'Stage 8',
      model_class = 'bayesian_cure',
      model_block = 'bayesian_cure',
      family = coalesce_character_columns(., c('latency_family', 'family_code'), default = 'other'),
      model_id = safe_character(model_id),
      benchmark_id = NA_character_,
      formula_name = formula_name,
      formula_label = formula_label,
      site_branch = coalesce_character_columns(., c('site_branch'), default = NA_character_),
      interaction_branch = coalesce_character_columns(., c('interaction_branch'), default = NA_character_),
      horizon_year = horizon_year,
      risk = safe_numeric(meanRisk_Bayes_mean),
      risk_lower = safe_numeric(meanRisk_Bayes_q025),
      risk_upper = safe_numeric(meanRisk_Bayes_q975),
      survival = safe_numeric(meanSurvival_mean),
      survival_lower = safe_numeric(meanSurvival_q025),
      survival_upper = safe_numeric(meanSurvival_q975),
      estimate_type = 'posterior_mean_with_credible_interval',
      reporting_status = make_reporting_status(coalesce_character_columns(., c('cohort_reporting_status')), fallback_reported = TRUE),
      reporting_note = combine_nonempty_text(coalesce_character_columns(., c('cohort_reporting_note')), screening_note),
      admissible_flag = admissible_flag,
      admissibility_reasons = coalesce_character_columns(., c('admissibility_reasons')),
      screening_flag = screening_flag,
      screening_note = screening_note
    ) %>%
    fill_formula_metadata(formula_lookup) %>%
    join_stage1_support(support_lookup, riskset_lookup) %>%
    mutate(
      model_class_label = make_model_class_label(model_class),
      model_display_label = make_model_display_label(model_class, family)
    )
}

standardize_stage8_threshold <- function(class_tbl, formula_lookup, support_lookup, riskset_lookup) {
  if (is.null(class_tbl) || nrow(class_tbl) == 0L) {
    return(tibble())
  }
  
  out <- tibble::as_tibble(class_tbl) %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      formula_name = coalesce_character_columns(., c('formula_anchor', 'formula_name')),
      formula_label = coalesce_character_columns(., c('formula_label'), default = NA_character_),
      horizon_year = safe_integer(horizon_year),
      threshold = safe_numeric(threshold),
      admissible_flag = if ('admissible_flag' %in% names(.)) as.logical(admissible_flag) else TRUE,
      classification_estimable_flag = if ('classification_estimable_flag' %in% names(.)) as.logical(classification_estimable_flag) else !is.na(FPR_mean),
      screening_flag = coalesce_character_columns(., c('screening_flag')),
      screening_note = coalesce_character_columns(., c('screening_detail', 'classification_reporting_note'))
    )
  
  if (isTRUE(keep_only_stage8_admissible_fits)) {
    out <- out %>% filter(dplyr::coalesce(admissible_flag, FALSE))
  }
  if (!isTRUE(include_stage8_unsupportive_sensitivity_fits)) {
    out <- out %>% filter(is.na(screening_flag) | screening_flag != 'unsupportive')
  }
  
  out %>%
    transmute(
      dataset = dataset,
      source_stage = 'Stage 8',
      model_class = 'bayesian_cure',
      model_block = 'bayesian_cure',
      family = coalesce_character_columns(., c('latency_family', 'family_code'), default = 'other'),
      model_id = safe_character(model_id),
      benchmark_id = NA_character_,
      formula_name = formula_name,
      formula_label = formula_label,
      site_branch = coalesce_character_columns(., c('site_branch'), default = NA_character_),
      interaction_branch = coalesce_character_columns(., c('interaction_branch'), default = NA_character_),
      horizon_year = horizon_year,
      threshold = threshold,
      positive_classification_rate = safe_numeric(pos_rate_mean),
      false_positive_rate = safe_numeric(FPR_mean),
      false_positive_rate_lower = safe_numeric(FPR_q025),
      false_positive_rate_upper = safe_numeric(FPR_q975),
      false_positive_burden = safe_numeric(false_positive_burden_mean),
      false_positive_burden_lower = safe_numeric(false_positive_burden_q025),
      false_positive_burden_upper = safe_numeric(false_positive_burden_q975),
      false_positive_per_100 = safe_numeric(FP100_mean),
      ppv = safe_numeric(PPV_mean),
      tpr = safe_numeric(TPR_mean),
      net_benefit = safe_numeric(NB_mean),
      net_benefit_lower = safe_numeric(NB_q025),
      net_benefit_upper = safe_numeric(NB_q975),
      net_reduction_unnecessary_per_100 = NA_real_,
      denom_case = safe_numeric(denom_case),
      denom_control = safe_numeric(denom_control),
      reporting_status = make_reporting_status(coalesce_character_columns(., c('classification_reporting_status')), fallback_reported = FALSE),
      reporting_note = combine_nonempty_text(coalesce_character_columns(., c('classification_reporting_note')), screening_note),
      admissible_flag = admissible_flag,
      admissibility_reasons = coalesce_character_columns(., c('admissibility_reasons')),
      screening_flag = screening_flag,
      screening_note = screening_note,
      classification_estimable_flag = classification_estimable_flag
    ) %>%
    fill_formula_metadata(formula_lookup) %>%
    join_stage1_support(support_lookup, riskset_lookup) %>%
    mutate(
      reporting_status = ifelse(!classification_estimable_flag & reporting_status == 'reported', 'suppressed', reporting_status),
      model_class_label = make_model_class_label(model_class),
      model_display_label = make_model_display_label(model_class, family)
    )
}

standardize_stage8_delta <- function(delta_tbl, formula_lookup, support_lookup, riskset_lookup) {
  if (is.null(delta_tbl) || nrow(delta_tbl) == 0L) {
    return(tibble())
  }
  
  metric_map <- c(
    meanRisk = 'delta_risk_nc_minus_cure',
    FPR = 'delta_false_positive_rate',
    FP100 = 'delta_false_positive_per_100',
    NB = 'delta_net_benefit',
    false_positive_burden = 'delta_false_positive_burden'
  )
  
  tibble::as_tibble(delta_tbl) %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      formula_name = coalesce_character_columns(., c('formula_anchor', 'formula_name')),
      formula_label = coalesce_character_columns(., c('formula_label'), default = NA_character_),
      horizon_year = safe_integer(horizon_year),
      threshold = safe_numeric(threshold),
      metric = safe_character(metric)
    ) %>%
    filter(metric %in% names(metric_map)) %>%
    transmute(
      dataset = dataset,
      source_stage = 'Stage 8',
      model_class = 'bayesian_cure',
      model_block = 'bayesian_cure',
      family = coalesce_character_columns(., c('latency_family', 'family_code'), default = 'other'),
      model_id = safe_character(model_id),
      no_cure_model_id = safe_character(no_cure_model_id),
      formula_name = formula_name,
      formula_label = formula_label,
      site_branch = coalesce_character_columns(., c('site_branch'), default = NA_character_),
      interaction_branch = coalesce_character_columns(., c('interaction_branch'), default = NA_character_),
      horizon_year = horizon_year,
      threshold = threshold,
      metric_domain = 'delta_vs_nocure',
      metric_name = unname(metric_map[metric]),
      metric_value = safe_numeric(delta_mean),
      metric_lower = safe_numeric(delta_q025),
      metric_upper = safe_numeric(delta_q975),
      reporting_status = make_reporting_status(coalesce_character_columns(., c('reporting_status', 'classification_reporting_status')), fallback_reported = TRUE),
      reporting_note = coalesce_character_columns(., c('reporting_note', 'classification_reporting_note', 'screening_detail')),
      screening_flag = coalesce_character_columns(., c('screening_flag')),
      screening_note = coalesce_character_columns(., c('screening_detail'))
    ) %>%
    fill_formula_metadata(formula_lookup) %>%
    join_stage1_support(support_lookup, riskset_lookup) %>%
    mutate(
      model_class_label = make_model_class_label(model_class),
      model_display_label = make_model_display_label(model_class, family)
    )
}

# 🔴 Read: Stage-1 backbone and staged outputs ===============================
## 🟠 Read: shared registries, analysis datasets, and model blocks ===============================
stage1_analysis_datasets_file <- file.path(data_path, 'stage1_analysis_datasets.rds')
stage1_formula_registry_file <- file.path(data_path, 'stage1_formula_registry.csv')
stage1_horizon_registry_file <- file.path(data_path, 'stage1_horizon_registry.csv')
stage1_threshold_registry_file <- file.path(data_path, 'stage1_threshold_registry.csv')
stage1_dataset_registry_file <- file.path(data_path, 'stage1_dataset_registry.csv')
stage1_metadata_registry_file <- file.path(data_path, 'stage1_metadata_registry.csv')

assert_required_file(stage1_analysis_datasets_file, 'Stage 1 analysis datasets', required = TRUE, allow_partial = FALSE)
assert_required_file(stage1_formula_registry_file, 'Stage 1 formula registry', required = TRUE, allow_partial = FALSE)
assert_required_file(stage1_horizon_registry_file, 'Stage 1 horizon registry', required = TRUE, allow_partial = FALSE)
assert_required_file(stage1_threshold_registry_file, 'Stage 1 threshold registry', required = TRUE, allow_partial = FALSE)

analysis_datasets <- readRDS(stage1_analysis_datasets_file)
formula_lookup <- normalize_formula_registry(read_csv_if_exists(stage1_formula_registry_file))
support_lookup <- normalize_horizon_registry(read_csv_if_exists(stage1_horizon_registry_file))
stage1_thresholds <- normalize_threshold_registry(read_csv_if_exists(stage1_threshold_registry_file))
stage1_dataset_registry <- read_csv_if_exists(stage1_dataset_registry_file) %||% tibble()
stage1_metadata_registry <- read_csv_if_exists(stage1_metadata_registry_file) %||% tibble()

comparison_thresholds <- sort(unique(safe_numeric(risk_thresholds_override %||% stage1_thresholds)))
if (length(comparison_thresholds) == 0L || anyNA(comparison_thresholds)) {
  stop('No valid threshold values were recovered for Stage 10.', call. = FALSE)
}
if (!identical(sort(unique(safe_integer(common_horizons_year))), 1:10)) {
  stop('Stage 10 must preserve the common 1:10 horizon grid.', call. = FALSE)
}

riskset_lookup <- compute_riskset_lookup(analysis_datasets, as.integer(common_horizons_year))
stage6_optional <- load_stage6_outputs_optional(stage6_screening_file)

stage3_outputs <- load_stage3_outputs(stage3_export_path, required = require_stage3_outputs)
stage5_outputs <- load_stage5_outputs(stage5_export_path, required = require_stage5_outputs)
stage7_outputs <- load_stage7_outputs(stage7_export_path, required = require_stage7_outputs)
stage8_outputs <- load_stage8_outputs(stage8_export_path, required = require_stage8_outputs)

# 🔴 Build: canonical lookups and risk-set summaries ===============================
## 🟠 Build: stage availability and comparison metadata ===============================
stage_presence_registry <- tibble(
  source_stage = c('Stage 3', 'Stage 5', 'Stage 7', 'Stage 8', 'Stage 9'),
  configured_path = c(stage3_export_path, stage5_export_path, stage7_export_path, stage8_export_path, stage9_export_path),
  path_exists = c(dir.exists(stage3_export_path), dir.exists(stage5_export_path), dir.exists(stage7_export_path), dir.exists(stage8_export_path), !is.na(stage9_export_path) & dir.exists(stage9_export_path)),
  loaded_primary_output = c(
    !is.null(stage3_outputs$group_risk) && nrow(stage3_outputs$group_risk) > 0L,
    !is.null(stage5_outputs$subject_horizon_risk_long) && nrow(stage5_outputs$subject_horizon_risk_long) > 0L,
    !is.null(stage7_outputs$risk_summary) && nrow(stage7_outputs$risk_summary) > 0L,
    !is.null(stage8_outputs$posterior_cohort_yearly) && nrow(stage8_outputs$posterior_cohort_yearly) > 0L,
    FALSE
  ),
  stage_requirement = c(
    ifelse(require_stage3_outputs, 'required', 'optional'),
    ifelse(require_stage5_outputs, 'required', 'optional'),
    ifelse(require_stage7_outputs, 'required', 'optional'),
    ifelse(require_stage8_outputs, 'required', 'optional'),
    'optional'
  )
)

# 🔴 Standardize: Stage-3 KM benchmark outputs ===============================
## 🟠 Standardize: risk and threshold tables ===============================
stage3_risk_table <- standardize_stage3_risk(
  group_risk_tbl = stage3_outputs$group_risk,
  formula_lookup = formula_lookup,
  support_lookup = support_lookup,
  riskset_lookup = riskset_lookup
)

stage3_threshold_table <- standardize_stage3_threshold(
  classification_tbl = stage3_outputs$classification,
  formula_lookup = formula_lookup,
  support_lookup = support_lookup,
  riskset_lookup = riskset_lookup
) %>%
  filter(threshold %in% comparison_thresholds)

# 🔴 Standardize: Stage-5 no-cure outputs ===============================
## 🟠 Standardize: individualized comparator summaries ===============================
stage5_risk_table <- standardize_stage5_risk(
  subject_risk_tbl = stage5_outputs$subject_horizon_risk_long,
  performance_tbl = stage5_outputs$model_performance_long,
  formula_lookup = formula_lookup,
  support_lookup = support_lookup,
  riskset_lookup = riskset_lookup
)

stage5_threshold_table <- standardize_stage5_threshold(
  performance_tbl = stage5_outputs$model_performance_long,
  formula_lookup = formula_lookup,
  support_lookup = support_lookup,
  riskset_lookup = riskset_lookup
) %>%
  filter(threshold %in% comparison_thresholds)

stage5_auxiliary_metric_table <- standardize_stage5_auxiliary(
  performance_tbl = stage5_outputs$model_performance_long,
  formula_lookup = formula_lookup,
  support_lookup = support_lookup,
  riskset_lookup = riskset_lookup
)

# 🔴 Standardize: Stage-7 frequentist cure outputs ===============================
## 🟠 Standardize: cure-MLE summaries and risk shifts ===============================
stage7_risk_table <- standardize_stage7_risk(
  risk_tbl = stage7_outputs$risk_summary,
  fit_registry_tbl = stage7_outputs$fit_registry,
  formula_lookup = formula_lookup,
  support_lookup = support_lookup,
  riskset_lookup = riskset_lookup
)

stage7_threshold_table <- standardize_stage7_threshold(
  threshold_tbl = stage7_outputs$threshold_metrics,
  fit_registry_tbl = stage7_outputs$fit_registry,
  formula_lookup = formula_lookup,
  support_lookup = support_lookup,
  riskset_lookup = riskset_lookup
) %>%
  filter(threshold %in% comparison_thresholds)

stage7_delta_table <- standardize_stage7_delta(
  delta_tbl = stage7_outputs$delta_risk,
  formula_lookup = formula_lookup,
  support_lookup = support_lookup,
  riskset_lookup = riskset_lookup
)

# 🔴 Standardize: Stage-8 Bayesian cure outputs ===============================
## 🟠 Standardize: posterior summaries and posterior deltas ===============================
stage8_risk_table <- standardize_stage8_risk(
  cohort_tbl = stage8_outputs$posterior_cohort_yearly,
  formula_lookup = formula_lookup,
  support_lookup = support_lookup,
  riskset_lookup = riskset_lookup
)

stage8_threshold_table <- standardize_stage8_threshold(
  class_tbl = stage8_outputs$posterior_classification,
  formula_lookup = formula_lookup,
  support_lookup = support_lookup,
  riskset_lookup = riskset_lookup
) %>%
  filter(threshold %in% comparison_thresholds)

stage8_delta_table <- standardize_stage8_delta(
  delta_tbl = stage8_outputs$posterior_delta_vs_nocure,
  formula_lookup = formula_lookup,
  support_lookup = support_lookup,
  riskset_lookup = riskset_lookup
) %>%
  filter(is.na(threshold) | threshold %in% comparison_thresholds)

# 🔴 Assemble: unified comparison engine ===============================
## 🟠 Assemble: risk, false-positive, auxiliary, and delta blocks ===============================
risk_probability_by_horizon_long <- bind_rows(
  stage3_risk_table,
  stage5_risk_table,
  stage7_risk_table,
  stage8_risk_table
) %>%
  mutate(
    family = dplyr::case_when(
      is.na(family) | family == '' ~ 'other',
      TRUE ~ family
    ),
    model_class = factor(model_class, levels = model_class_order),
    family_for_order = factor(ifelse(family %in% family_order, family, 'other'), levels = family_order)
  ) %>%
  arrange(dataset_order[match(as.character(dataset), dataset_order)], formula_label, model_class, family_for_order, horizon_year) %>%
  mutate(model_class = as.character(model_class)) %>%
  select(-family_for_order)

false_positive_rate_by_horizon_long <- bind_rows(
  stage3_threshold_table,
  stage5_threshold_table,
  stage7_threshold_table,
  stage8_threshold_table
) %>%
  mutate(
    family = dplyr::case_when(
      is.na(family) | family == '' ~ 'other',
      TRUE ~ family
    ),
    model_class = factor(model_class, levels = model_class_order),
    family_for_order = factor(ifelse(family %in% family_order, family, 'other'), levels = family_order)
  ) %>%
  arrange(dataset_order[match(as.character(dataset), dataset_order)], formula_label, threshold, model_class, family_for_order, horizon_year) %>%
  mutate(model_class = as.character(model_class)) %>%
  select(-family_for_order)

auxiliary_metrics_long <- stage5_auxiliary_metric_table %>%
  mutate(
    family = dplyr::case_when(
      is.na(family) | family == '' ~ 'other',
      TRUE ~ family
    )
  ) %>%
  arrange(dataset_order[match(as.character(dataset), dataset_order)], formula_label, model_id, horizon_year, metric_domain, metric_name)

cure_vs_nocure_shift_long <- bind_rows(
  stage7_delta_table,
  stage8_delta_table
) %>%
  mutate(
    family = dplyr::case_when(
      is.na(family) | family == '' ~ 'other',
      TRUE ~ family
    )
  ) %>%
  arrange(dataset_order[match(as.character(dataset), dataset_order)], formula_label, model_id, horizon_year, threshold, metric_name)

## 🟠 Assemble: source-of-truth master long table ===============================
master_risk_long <- risk_probability_by_horizon_long %>%
  transmute(
    dataset,
    source_stage,
    model_class,
    model_class_label,
    model_block,
    family,
    model_id,
    benchmark_id,
    formula_name,
    formula_label,
    site_branch,
    interaction_branch,
    horizon_year,
    threshold = NA_real_,
    metric_domain = 'horizon_comparison',
    metric_name = 'risk',
    metric_value = risk,
    metric_lower = risk_lower,
    metric_upper = risk_upper,
    support_label,
    support_priority,
    primary_supported_flag,
    interpretation_note,
    reporting_status,
    reporting_note,
    risk_set_size,
    event_count_by_horizon,
    known_nonevent_by_horizon,
    instability_marker,
    model_display_label
  ) %>%
  bind_rows(
    risk_probability_by_horizon_long %>%
      transmute(
        dataset,
        source_stage,
        model_class,
        model_class_label,
        model_block,
        family,
        model_id,
        benchmark_id,
        formula_name,
        formula_label,
        site_branch,
        interaction_branch,
        horizon_year,
        threshold = NA_real_,
        metric_domain = 'horizon_comparison',
        metric_name = 'survival',
        metric_value = survival,
        metric_lower = survival_lower,
        metric_upper = survival_upper,
        support_label,
        support_priority,
        primary_supported_flag,
        interpretation_note,
        reporting_status,
        reporting_note,
        risk_set_size,
        event_count_by_horizon,
        known_nonevent_by_horizon,
        instability_marker,
        model_display_label
      )
  )

master_threshold_long <- bind_rows(
  false_positive_rate_by_horizon_long %>%
    transmute(
      dataset,
      source_stage,
      model_class,
      model_class_label,
      model_block,
      family,
      model_id,
      benchmark_id,
      formula_name,
      formula_label,
      site_branch,
      interaction_branch,
      horizon_year,
      threshold,
      metric_domain = 'threshold_comparison',
      metric_name = 'positive_classification_rate',
      metric_value = positive_classification_rate,
      metric_lower = NA_real_,
      metric_upper = NA_real_,
      support_label,
      support_priority,
      primary_supported_flag,
      interpretation_note,
      reporting_status,
      reporting_note,
      risk_set_size,
      event_count_by_horizon,
      known_nonevent_by_horizon,
      instability_marker,
      model_display_label
    ),
  false_positive_rate_by_horizon_long %>%
    transmute(
      dataset, source_stage, model_class, model_class_label, model_block, family, model_id, benchmark_id,
      formula_name, formula_label, site_branch, interaction_branch, horizon_year, threshold,
      metric_domain = 'threshold_comparison', metric_name = 'false_positive_rate',
      metric_value = false_positive_rate, metric_lower = false_positive_rate_lower, metric_upper = false_positive_rate_upper,
      support_label, support_priority, primary_supported_flag, interpretation_note, reporting_status, reporting_note,
      risk_set_size, event_count_by_horizon, known_nonevent_by_horizon, instability_marker, model_display_label
    ),
  false_positive_rate_by_horizon_long %>%
    transmute(
      dataset, source_stage, model_class, model_class_label, model_block, family, model_id, benchmark_id,
      formula_name, formula_label, site_branch, interaction_branch, horizon_year, threshold,
      metric_domain = 'threshold_comparison', metric_name = 'false_positive_burden',
      metric_value = false_positive_burden, metric_lower = false_positive_burden_lower, metric_upper = false_positive_burden_upper,
      support_label, support_priority, primary_supported_flag, interpretation_note, reporting_status, reporting_note,
      risk_set_size, event_count_by_horizon, known_nonevent_by_horizon, instability_marker, model_display_label
    ),
  false_positive_rate_by_horizon_long %>%
    transmute(
      dataset, source_stage, model_class, model_class_label, model_block, family, model_id, benchmark_id,
      formula_name, formula_label, site_branch, interaction_branch, horizon_year, threshold,
      metric_domain = 'threshold_comparison', metric_name = 'false_positive_per_100',
      metric_value = false_positive_per_100, metric_lower = NA_real_, metric_upper = NA_real_,
      support_label, support_priority, primary_supported_flag, interpretation_note, reporting_status, reporting_note,
      risk_set_size, event_count_by_horizon, known_nonevent_by_horizon, instability_marker, model_display_label
    ),
  false_positive_rate_by_horizon_long %>%
    transmute(
      dataset, source_stage, model_class, model_class_label, model_block, family, model_id, benchmark_id,
      formula_name, formula_label, site_branch, interaction_branch, horizon_year, threshold,
      metric_domain = 'threshold_comparison', metric_name = 'ppv',
      metric_value = ppv, metric_lower = NA_real_, metric_upper = NA_real_,
      support_label, support_priority, primary_supported_flag, interpretation_note, reporting_status, reporting_note,
      risk_set_size, event_count_by_horizon, known_nonevent_by_horizon, instability_marker, model_display_label
    ),
  false_positive_rate_by_horizon_long %>%
    transmute(
      dataset, source_stage, model_class, model_class_label, model_block, family, model_id, benchmark_id,
      formula_name, formula_label, site_branch, interaction_branch, horizon_year, threshold,
      metric_domain = 'threshold_comparison', metric_name = 'tpr',
      metric_value = tpr, metric_lower = NA_real_, metric_upper = NA_real_,
      support_label, support_priority, primary_supported_flag, interpretation_note, reporting_status, reporting_note,
      risk_set_size, event_count_by_horizon, known_nonevent_by_horizon, instability_marker, model_display_label
    ),
  false_positive_rate_by_horizon_long %>%
    transmute(
      dataset, source_stage, model_class, model_class_label, model_block, family, model_id, benchmark_id,
      formula_name, formula_label, site_branch, interaction_branch, horizon_year, threshold,
      metric_domain = 'threshold_comparison', metric_name = 'net_benefit',
      metric_value = net_benefit, metric_lower = net_benefit_lower, metric_upper = net_benefit_upper,
      support_label, support_priority, primary_supported_flag, interpretation_note, reporting_status, reporting_note,
      risk_set_size, event_count_by_horizon, known_nonevent_by_horizon, instability_marker, model_display_label
    ),
  false_positive_rate_by_horizon_long %>%
    transmute(
      dataset, source_stage, model_class, model_class_label, model_block, family, model_id, benchmark_id,
      formula_name, formula_label, site_branch, interaction_branch, horizon_year, threshold,
      metric_domain = 'threshold_comparison', metric_name = 'net_reduction_unnecessary_per_100',
      metric_value = net_reduction_unnecessary_per_100, metric_lower = NA_real_, metric_upper = NA_real_,
      support_label, support_priority, primary_supported_flag, interpretation_note, reporting_status, reporting_note,
      risk_set_size, event_count_by_horizon, known_nonevent_by_horizon, instability_marker, model_display_label
    )
)

master_auxiliary_long <- auxiliary_metrics_long %>%
  transmute(
    dataset,
    source_stage,
    model_class,
    model_class_label,
    model_block,
    family,
    model_id,
    benchmark_id,
    formula_name,
    formula_label,
    site_branch,
    interaction_branch,
    horizon_year,
    threshold,
    metric_domain,
    metric_name,
    metric_value,
    metric_lower,
    metric_upper,
    support_label,
    support_priority,
    primary_supported_flag,
    interpretation_note,
    reporting_status,
    reporting_note,
    risk_set_size,
    event_count_by_horizon,
    known_nonevent_by_horizon,
    instability_marker,
    model_display_label
  )

master_delta_long <- cure_vs_nocure_shift_long %>%
  transmute(
    dataset,
    source_stage,
    model_class,
    model_class_label,
    model_block,
    family,
    model_id,
    benchmark_id = no_cure_model_id,
    formula_name,
    formula_label,
    site_branch,
    interaction_branch,
    horizon_year,
    threshold,
    metric_domain,
    metric_name,
    metric_value,
    metric_lower,
    metric_upper,
    support_label,
    support_priority,
    primary_supported_flag,
    interpretation_note,
    reporting_status,
    reporting_note,
    risk_set_size,
    event_count_by_horizon,
    known_nonevent_by_horizon,
    instability_marker,
    model_display_label
  )

master_comparison_long <- bind_rows(
  master_risk_long,
  master_threshold_long,
  master_auxiliary_long,
  master_delta_long
) %>%
  filter(is.na(threshold) | threshold %in% comparison_thresholds) %>%
  arrange(
    factor(dataset, levels = dataset_order),
    factor(model_class, levels = model_class_order),
    formula_name,
    family,
    horizon_year,
    threshold,
    metric_domain,
    metric_name
  )

## 🟠 Assemble: convenience wide tables for manuscript use ===============================
risk_probability_by_horizon_wide <- risk_probability_by_horizon_long %>%
  select(
    dataset,
    formula_name,
    formula_label,
    site_branch,
    interaction_branch,
    source_stage,
    model_class,
    model_class_label,
    family,
    model_id,
    benchmark_id,
    model_display_label,
    horizon_year,
    risk,
    risk_lower,
    risk_upper,
    survival,
    survival_lower,
    survival_upper,
    support_label,
    reporting_status,
    risk_set_size
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = horizon_year,
    values_from = c(risk, risk_lower, risk_upper, survival, survival_lower, survival_upper, support_label, reporting_status, risk_set_size),
    names_glue = '{.value}_year_{horizon_year}'
  ) %>%
  arrange(factor(dataset, levels = dataset_order), formula_name, factor(model_class, levels = model_class_order), family)

false_positive_rate_by_horizon_wide <- false_positive_rate_by_horizon_long %>%
  select(
    dataset,
    formula_name,
    formula_label,
    site_branch,
    interaction_branch,
    source_stage,
    model_class,
    model_class_label,
    family,
    model_id,
    benchmark_id,
    model_display_label,
    threshold,
    horizon_year,
    false_positive_rate,
    false_positive_rate_lower,
    false_positive_rate_upper,
    false_positive_burden,
    false_positive_burden_lower,
    false_positive_burden_upper,
    false_positive_per_100,
    ppv,
    tpr,
    net_benefit,
    net_benefit_lower,
    net_benefit_upper,
    support_label,
    reporting_status,
    reporting_note
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = horizon_year,
    values_from = c(false_positive_rate, false_positive_rate_lower, false_positive_rate_upper, false_positive_burden, false_positive_burden_lower, false_positive_burden_upper, false_positive_per_100, ppv, tpr, net_benefit, net_benefit_lower, net_benefit_upper, support_label, reporting_status),
    names_glue = '{.value}_year_{horizon_year}'
  ) %>%
  arrange(factor(dataset, levels = dataset_order), formula_name, threshold, factor(model_class, levels = model_class_order), family)

## 🟠 Assemble: integrated model registry ===============================
stage10_model_registry <- bind_rows(
  risk_probability_by_horizon_long %>%
    transmute(
      dataset,
      source_stage,
      model_class,
      model_class_label,
      model_block,
      family,
      model_id,
      benchmark_id,
      formula_name,
      formula_label,
      site_branch,
      interaction_branch,
      model_display_label,
      has_risk_rows = TRUE,
      has_threshold_rows = FALSE,
      has_auxiliary_rows = FALSE,
      has_delta_rows = FALSE
    ),
  false_positive_rate_by_horizon_long %>%
    transmute(
      dataset,
      source_stage,
      model_class,
      model_class_label,
      model_block,
      family,
      model_id,
      benchmark_id,
      formula_name,
      formula_label,
      site_branch,
      interaction_branch,
      model_display_label,
      has_risk_rows = FALSE,
      has_threshold_rows = TRUE,
      has_auxiliary_rows = FALSE,
      has_delta_rows = FALSE
    ),
  auxiliary_metrics_long %>%
    transmute(
      dataset,
      source_stage,
      model_class,
      model_class_label,
      model_block,
      family,
      model_id,
      benchmark_id,
      formula_name,
      formula_label,
      site_branch,
      interaction_branch,
      model_display_label,
      has_risk_rows = FALSE,
      has_threshold_rows = FALSE,
      has_auxiliary_rows = TRUE,
      has_delta_rows = FALSE
    ),
  cure_vs_nocure_shift_long %>%
    transmute(
      dataset,
      source_stage,
      model_class,
      model_class_label,
      model_block,
      family,
      model_id,
      benchmark_id = no_cure_model_id,
      formula_name,
      formula_label,
      site_branch,
      interaction_branch,
      model_display_label,
      has_risk_rows = FALSE,
      has_threshold_rows = FALSE,
      has_auxiliary_rows = FALSE,
      has_delta_rows = TRUE
    )
) %>%
  group_by(
    dataset,
    source_stage,
    model_class,
    model_class_label,
    model_block,
    family,
    model_id,
    benchmark_id,
    formula_name,
    formula_label,
    site_branch,
    interaction_branch,
    model_display_label
  ) %>%
  summarise(
    has_risk_rows = any(has_risk_rows),
    has_threshold_rows = any(has_threshold_rows),
    has_auxiliary_rows = any(has_auxiliary_rows),
    has_delta_rows = any(has_delta_rows),
    n_supported_horizons = sum(risk_probability_by_horizon_long$model_id == first(model_id) & risk_probability_by_horizon_long$dataset == first(dataset), na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(factor(dataset, levels = dataset_order), formula_name, factor(model_class, levels = model_class_order), family)

# 🔴 Verify: unified-engine completeness and required tables ===============================
## 🟠 Verify: mandatory comparison outputs ===============================
if (nrow(risk_probability_by_horizon_long) == 0L) {
  stop('Stage 10 could not recover any risk-by-horizon rows. At least Stage 3 and Stage 5 outputs are required.', call. = FALSE)
}

if (nrow(false_positive_rate_by_horizon_long) == 0L) {
  warning('Stage 10 could not recover any threshold-based false-positive rows. The FPR comparison table will still be written, but it may be empty.', call. = FALSE)
}

required_risk_horizons <- as.integer(common_horizons_year)
missing_risk_horizons <- setdiff(required_risk_horizons, sort(unique(risk_probability_by_horizon_long$horizon_year)))
if (length(missing_risk_horizons) > 0L) {
  warning(
    sprintf('Some horizons are absent from the risk comparison table: %s', paste(missing_risk_horizons, collapse = ', ')),
    call. = FALSE
  )
}

# 🔴 Draw: time-wise risk and false-positive comparison figures ===============================
## 🟠 Draw: plot pages from the same CSV-ready data frames ===============================
summary_pdf_file <- file.path(export_path, 'stage10_summary_plots.pdf')
grDevices::pdf(summary_pdf_file, width = 13, height = 8.5, onefile = TRUE)

if (nrow(risk_probability_by_horizon_long) > 0L) {
  risk_plot_data <- risk_probability_by_horizon_long %>%
    filter(horizon_year %in% as.integer(common_horizons_year)) %>%
    mutate(
      dataset_formula_panel = factor(
        paste(dataset, formula_label, sep = ' | '),
        levels = unique(paste(dataset, formula_label, sep = ' | '))
      )
    )
  
  print(
    ggplot(risk_plot_data, aes(x = horizon_year, y = risk, color = model_display_label, fill = model_display_label, group = model_display_label)) +
      geom_ribbon(aes(ymin = risk_lower, ymax = risk_upper), alpha = 0.12, linewidth = 0, na.rm = TRUE, show.legend = FALSE) +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      geom_point(size = 1.4, na.rm = TRUE) +
      facet_wrap(~ dataset_formula_panel, scales = 'free_y') +
      scale_x_continuous(breaks = as.integer(common_horizons_year)) +
      labs(
        title = 'Stage 10: risk probability comparison across time horizons',
        x = 'Years after cohort entry',
        y = 'Predicted transition risk',
        color = 'Model'
      ) +
      theme_bw() +
      theme(legend.position = 'bottom')
  )
}

if (nrow(false_positive_rate_by_horizon_long) > 0L) {
  fpr_plot_data <- false_positive_rate_by_horizon_long %>%
    filter(horizon_year %in% as.integer(common_horizons_year)) %>%
    mutate(
      dataset_formula_panel = factor(
        paste(dataset, formula_label, sep = ' | '),
        levels = unique(paste(dataset, formula_label, sep = ' | '))
      )
    )
  
  print(
    ggplot(fpr_plot_data, aes(x = horizon_year, y = false_positive_rate, color = model_display_label, fill = model_display_label, group = model_display_label)) +
      geom_ribbon(aes(ymin = false_positive_rate_lower, ymax = false_positive_rate_upper), alpha = 0.12, linewidth = 0, na.rm = TRUE, show.legend = FALSE) +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      geom_point(size = 1.3, na.rm = TRUE) +
      facet_grid(dataset_formula_panel ~ threshold, scales = 'free_y') +
      scale_x_continuous(breaks = as.integer(common_horizons_year)) +
      labs(
        title = 'Stage 10: false-positive rate comparison across time horizons',
        x = 'Years after cohort entry',
        y = 'False-positive rate',
        color = 'Model'
      ) +
      theme_bw() +
      theme(legend.position = 'bottom')
  )
  
  print(
    ggplot(fpr_plot_data, aes(x = horizon_year, y = false_positive_burden, color = model_display_label, fill = model_display_label, group = model_display_label)) +
      geom_ribbon(aes(ymin = false_positive_burden_lower, ymax = false_positive_burden_upper), alpha = 0.12, linewidth = 0, na.rm = TRUE, show.legend = FALSE) +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      geom_point(size = 1.3, na.rm = TRUE) +
      facet_grid(dataset_formula_panel ~ threshold, scales = 'free_y') +
      scale_x_continuous(breaks = as.integer(common_horizons_year)) +
      labs(
        title = 'Stage 10: false-positive burden comparison across time horizons',
        x = 'Years after cohort entry',
        y = 'False-positive burden',
        color = 'Model'
      ) +
      theme_bw() +
      theme(legend.position = 'bottom')
  )
}

if (nrow(cure_vs_nocure_shift_long) > 0L) {
  shift_plot_data <- cure_vs_nocure_shift_long %>%
    filter(metric_name %in% c('delta_risk_nc_minus_cure', 'delta_false_positive_rate', 'delta_false_positive_burden', 'delta_net_benefit')) %>%
    mutate(
      dataset_formula_panel = factor(
        paste(dataset, formula_label, sep = ' | '),
        levels = unique(paste(dataset, formula_label, sep = ' | '))
      )
    )
  
  print(
    ggplot(shift_plot_data, aes(x = horizon_year, y = metric_value, color = model_display_label, fill = model_display_label, group = model_display_label)) +
      geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
      geom_ribbon(aes(ymin = metric_lower, ymax = metric_upper), alpha = 0.12, linewidth = 0, na.rm = TRUE, show.legend = FALSE) +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      geom_point(size = 1.3, na.rm = TRUE) +
      facet_grid(metric_name ~ dataset_formula_panel, scales = 'free_y') +
      scale_x_continuous(breaks = as.integer(common_horizons_year)) +
      labs(
        title = 'Stage 10: cure versus no-cure shifts on the common comparison grid',
        x = 'Years after cohort entry',
        y = 'Delta (no-cure minus cure)',
        color = 'Model'
      ) +
      theme_bw() +
      theme(legend.position = 'bottom')
  )
}

grDevices::dev.off()

# 🔴 Export: comparison tables and manuscript-ready artifacts ===============================
## 🟠 Export: CSV, PDF, and manifest outputs ===============================
stage10_metadata <- tibble(
  metadata_group = c(
    'stage', 'stage', 'stage', 'stage',
    'paths', 'paths', 'paths', 'paths', 'paths', 'paths',
    'inputs', 'inputs', 'inputs', 'inputs',
    'comparison', 'comparison', 'comparison', 'comparison',
    'notes'
  ),
  metadata_name = c(
    'stage_name', 'stage_role', 'allow_partial_stage10_build', 'stage6_optional_join',
    'data_path', 'stage3_export_path', 'stage5_export_path', 'stage7_export_path', 'stage8_export_path', 'export_path',
    'stage1_analysis_datasets_file', 'stage1_formula_registry_file', 'stage1_horizon_registry_file', 'stage1_threshold_registry_file',
    'common_horizon_vector', 'comparison_threshold_vector', 'primary_output_source_of_truth', 'plot_pdf_file',
    'main_scientific_priority'
  ),
  metadata_value = c(
    'Stage 10 unified comparison engine',
    'Combine KM benchmark, no-cure, frequentist cure, and Bayesian cure outputs on one common comparison grid.',
    as.character(allow_partial_stage10_build),
    ifelse(is.null(stage6_optional) || nrow(stage6_optional) == 0L, 'absent_or_not_joined', 'loaded_optional'),
    normalize_existing_path(data_path),
    normalize_existing_path(stage3_export_path),
    normalize_existing_path(stage5_export_path),
    normalize_existing_path(stage7_export_path),
    normalize_existing_path(stage8_export_path),
    normalize_existing_path(export_path),
    normalize_existing_path(stage1_analysis_datasets_file),
    normalize_existing_path(stage1_formula_registry_file),
    normalize_existing_path(stage1_horizon_registry_file),
    normalize_existing_path(stage1_threshold_registry_file),
    paste(as.integer(common_horizons_year), collapse = ','),
    paste(format(comparison_thresholds, trim = TRUE, scientific = FALSE), collapse = ','),
    'stage10_master_comparison_long.csv',
    basename(summary_pdf_file),
    'Time-wise false-positive rate comparison and time-wise risk probability comparison remain the primary Stage 10 tables.'
  )
)

risk_long_file <- file.path(export_path, 'stage10_risk_probability_by_horizon_long.csv')
risk_wide_file <- file.path(export_path, 'stage10_risk_probability_by_horizon_wide.csv')
fpr_long_file <- file.path(export_path, 'stage10_false_positive_rate_by_horizon_long.csv')
fpr_wide_file <- file.path(export_path, 'stage10_false_positive_rate_by_horizon_wide.csv')
aux_long_file <- file.path(export_path, 'stage10_auxiliary_metrics_long.csv')
delta_long_file <- file.path(export_path, 'stage10_cure_vs_nocure_shift_long.csv')
master_long_file <- file.path(export_path, 'stage10_master_comparison_long.csv')
model_registry_file <- file.path(export_path, 'stage10_model_registry.csv')
metadata_file <- file.path(export_path, 'stage10_stage_metadata.csv')
presence_registry_file <- file.path(export_path, 'stage10_stage_presence_registry.csv')
manifest_file <- file.path(export_path, 'stage10_output_manifest.csv')

readr::write_csv(risk_probability_by_horizon_long, risk_long_file)
readr::write_csv(risk_probability_by_horizon_wide, risk_wide_file)
readr::write_csv(false_positive_rate_by_horizon_long, fpr_long_file)
readr::write_csv(false_positive_rate_by_horizon_wide, fpr_wide_file)
readr::write_csv(auxiliary_metrics_long, aux_long_file)
readr::write_csv(cure_vs_nocure_shift_long, delta_long_file)
readr::write_csv(master_comparison_long, master_long_file)
readr::write_csv(stage10_model_registry, model_registry_file)
readr::write_csv(stage10_metadata, metadata_file)
readr::write_csv(stage_presence_registry, presence_registry_file)

stage10_output_manifest <- tibble(
  file_name = c(
    basename(risk_long_file),
    basename(risk_wide_file),
    basename(fpr_long_file),
    basename(fpr_wide_file),
    basename(aux_long_file),
    basename(delta_long_file),
    basename(master_long_file),
    basename(model_registry_file),
    basename(metadata_file),
    basename(presence_registry_file),
    basename(summary_pdf_file)
  ),
  file_role = c(
    'time-wise risk probability comparison long table',
    'time-wise risk probability comparison wide table',
    'time-wise false-positive rate comparison long table',
    'time-wise false-positive rate comparison wide table',
    'auxiliary discrimination/calibration/Brier metrics',
    'cure versus no-cure shift table',
    'source-of-truth master long comparison table',
    'integrated model registry',
    'Stage 10 metadata registry',
    'stage availability registry',
    'summary plot PDF generated from CSV-ready data frames'
  ),
  path = c(
    normalize_existing_path(risk_long_file),
    normalize_existing_path(risk_wide_file),
    normalize_existing_path(fpr_long_file),
    normalize_existing_path(fpr_wide_file),
    normalize_existing_path(aux_long_file),
    normalize_existing_path(delta_long_file),
    normalize_existing_path(master_long_file),
    normalize_existing_path(model_registry_file),
    normalize_existing_path(metadata_file),
    normalize_existing_path(presence_registry_file),
    normalize_existing_path(summary_pdf_file)
  )
)

readr::write_csv(stage10_output_manifest, manifest_file)

cat('\nStage 10 completed. Key outputs written to:\n')
print(stage10_output_manifest)