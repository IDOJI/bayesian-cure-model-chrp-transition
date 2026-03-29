# 🔴 Configure: stage-output paths and export target ===============================
data_path <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results'
export_path <- file.path(data_path, 'stage10_Unified comparison')

stage_paths <- list(
  stage1 = file.path(data_path, 'stage1_Backbone lock'),
  stage2 = file.path(data_path, 'stage2_Follow-up maturity'),
  stage3 = file.path(data_path, 'stage3_KM benchmark classifier'),
  stage4 = file.path(data_path, 'stage4_Timing difference'),
  stage5 = file.path(data_path, 'stage5_Noncure block'),
  stage6 = file.path(data_path, 'stage6_Cure appropriateness screening'),
  stage7 = file.path(data_path, 'stage7_Frequentist cure block'),
  stage8 = file.path(data_path, 'stage8_Bayesian cure block'),
  stage9 = file.path(data_path, 'stage9_Remission sensitivity')
)

file_overrides <- list(
  stage1 = list(backbone_bundle = NULL),
  stage2 = list(followup = NULL),
  stage3 = list(prediction = NULL, classification = NULL, performance = NULL, decision = NULL),
  stage4 = list(timing = NULL),
  stage5 = list(prediction = NULL, classification = NULL, performance = NULL, decision = NULL),
  stage6 = list(carryforward = NULL),
  stage7 = list(prediction = NULL, classification = NULL, performance = NULL, decision = NULL, decomposition = NULL, fit_contrast = NULL),
  stage8 = list(prediction = NULL, classification = NULL, performance = NULL, decision = NULL, decomposition = NULL, remission_delta = NULL, anchor_delta = NULL, incidence_update = NULL, diagnostics = NULL),
  stage9 = list(prediction = NULL, classification = NULL, performance = NULL, decision = NULL)
)

# 🔴 Initialize: packages and runtime options ===============================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(stringr)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# 🔴 Define: low-level readers and scalar helpers ===============================
## 🟠 Define: file readers ===============================
normalize_existing_path <- function(path) {
  normalizePath(path, winslash = '/', mustWork = FALSE)
}

read_table_any <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
    stop(sprintf('Table path does not exist: %s', path), call. = FALSE)
  }
  
  ext <- tolower(tools::file_ext(path))
  
  if (ext %in% c('csv', 'txt')) {
    return(
      readr::read_csv(
        file = path,
        show_col_types = FALSE,
        progress = FALSE,
        guess_max = 100000
      )
    )
  }
  
  if (ext == 'rds') {
    obj <- readRDS(path)
    
    if (inherits(obj, 'data.frame')) {
      return(tibble::as_tibble(obj))
    }
    
    if (is.list(obj)) {
      df_candidates <- purrr::keep(obj, ~ inherits(.x, 'data.frame'))
      if (length(df_candidates) == 1L) {
        return(tibble::as_tibble(df_candidates[[1L]]))
      }
    }
    
    stop(sprintf('RDS is not directly readable as a single table: %s', path), call. = FALSE)
  }
  
  stop(sprintf('Unsupported table extension for `%s`.', path), call. = FALSE)
}

safe_csv_header <- function(path) {
  tryCatch(
    names(readr::read_csv(path, n_max = 0, show_col_types = FALSE, progress = FALSE)),
    error = function(e) character(0)
  )
}

coerce_numeric_text <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

coalesce_cols <- function(df, candidates, default = NA) {
  found <- candidates[candidates %in% names(df)]
  
  if (length(found) == 0L) {
    return(rep(default, nrow(df)))
  }
  
  out <- df[[found[1L]]]
  if (length(found) >= 2L) {
    for (nm in found[-1L]) {
      out <- dplyr::coalesce(out, df[[nm]])
    }
  }
  
  out
}

first_non_missing_scalar <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(NA)
  }
  
  if (is.list(x)) {
    x <- unlist(x, recursive = FALSE, use.names = FALSE)
  }
  
  if (is.character(x)) {
    x <- trimws(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0L) return(NA_character_)
    return(x[[1L]])
  }
  
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA)
  x[[1L]]
}

trim_lower <- function(x) {
  tolower(trimws(as.character(x)))
}

header_has_any <- function(headers, candidates) {
  any(tolower(candidates) %in% tolower(headers))
}

name_has_any <- function(file_name, patterns) {
  any(vapply(patterns, function(ptn) stringr::str_detect(file_name, regex(ptn, ignore_case = TRUE)), logical(1)))
}

dataset_sort_num <- function(x, order_vec) {
  out <- match(as.character(x), order_vec)
  out[is.na(out)] <- length(order_vec) + 1L
  out
}

# 🔴 Define: backbone-compatible support helpers ===============================
## 🟠 Define: horizon-support rules ===============================
derive_support_tier <- function(dataset_name, horizon_year) {
  ds <- as.character(dataset_name)
  h0 <- as.integer(horizon_year)
  
  if (ds == 'PNU') {
    if (h0 == 1L) return('primary_supported')
    if (h0 == 2L) return('sensitivity')
    return('projection')
  }
  
  if (ds %in% c('SNU', 'merged')) {
    if (h0 %in% c(1L, 2L)) return('primary_supported')
    if (h0 <= 5L) return('secondary')
    return('projection')
  }
  
  NA_character_
}

derive_horizon_evidence_class <- function(dataset_name, horizon_year) {
  ds <- as.character(dataset_name)
  h0 <- as.integer(horizon_year)
  
  if (ds == 'PNU') {
    if (h0 == 1L) return('directly_observed_data_supported')
    if (h0 == 2L) return('partly_model_dependent')
    return('mostly_extrapolated')
  }
  
  if (ds %in% c('SNU', 'merged')) {
    if (h0 %in% c(1L, 2L)) return('directly_observed_data_supported')
    if (h0 <= 5L) return('partly_model_dependent')
    return('mostly_extrapolated')
  }
  
  NA_character_
}

derive_claim_restriction_flag <- function(horizon_evidence_class, support_tier, prior_tail_sensitive = FALSE) {
  if (is.na(horizon_evidence_class) || !nzchar(horizon_evidence_class)) {
    return(NA_character_)
  }
  
  if (horizon_evidence_class == 'directly_observed_data_supported') {
    return('primary_claim_allowed')
  }
  
  if (horizon_evidence_class == 'partly_model_dependent') {
    return('secondary_or_sensitivity_only')
  }
  
  if (horizon_evidence_class == 'mostly_extrapolated' && isTRUE(prior_tail_sensitive)) {
    return('projection_plus_prior_sensitive')
  }
  
  if (horizon_evidence_class == 'mostly_extrapolated') {
    return('projection_only')
  }
  
  'projection_only'
}

# 🔴 Define: value normalizers and row keys ===============================
## 🟠 Define: dataset and scale normalizers ===============================
normalize_dataset_key <- function(x) {
  x0 <- trimws(as.character(x))
  xl <- tolower(gsub('[^a-z0-9]+', '_', x0))
  
  out <- dplyr::case_when(
    xl %in% c('', 'na') ~ NA_character_,
    xl %in% c('pnu', 'site_pnu', 'pnu_site') ~ 'PNU',
    xl %in% c('snu', 'site_snu', 'snu_site') ~ 'SNU',
    xl %in% c('merged', 'overall_merged') ~ 'merged',
    xl %in% c('merged_site_free', 'merged__site_free', 'site_free', 'merged_free') ~ 'merged__site_free',
    xl %in% c('merged_site_adjusted', 'merged__site_adjusted', 'site_adjusted', 'merged_adjusted') ~ 'merged__site_adjusted',
    xl %in% c('overall_project', 'project_summary') ~ 'overall_project',
    TRUE ~ x0
  )
  
  out
}

dataset_parent_from_key <- function(dataset_key) {
  dplyr::case_when(
    dataset_key %in% c('PNU', 'SNU', 'merged') ~ dataset_key,
    dataset_key %in% c('merged__site_free', 'merged__site_adjusted') ~ 'merged',
    dataset_key == 'overall_project' ~ 'overall_project',
    TRUE ~ dataset_key
  )
}

derive_analysis_structure <- function(dataset_key) {
  dplyr::case_when(
    dataset_key %in% c('PNU', 'SNU') ~ 'single_site',
    dataset_key == 'merged__site_free' ~ 'site_free',
    dataset_key == 'merged__site_adjusted' ~ 'site_adjusted',
    dataset_key == 'merged' ~ 'merged_unspecified',
    dataset_key == 'overall_project' ~ 'project_summary',
    TRUE ~ 'unspecified'
  )
}

derive_reporting_priority <- function(dataset_key) {
  dplyr::case_when(
    dataset_key == 'merged__site_adjusted' ~ 'preferred_for_reporting',
    dataset_key == 'merged__site_free' ~ 'sensitivity_only',
    dataset_key == 'overall_project' ~ 'project_summary',
    TRUE ~ 'standard'
  )
}

normalize_risk_scale <- function(x, default = NA_character_) {
  x0 <- trimws(as.character(x))
  xl <- tolower(gsub('[^a-z0-9]+', '_', x0))
  
  out <- dplyr::case_when(
    xl %in% c('', 'na') ~ default,
    xl %in% c('transition_only_main', 'transition_only', 'main', 'main_scale') ~ 'transition_only_main',
    xl %in% c('transition_cif_competing', 'transition_cif', 'cif_competing', 'competing', 'competing_risk') ~ 'transition_cif_competing',
    TRUE ~ x0
  )
  
  out
}

normalize_prior_branch <- function(x) {
  x0 <- trimws(as.character(x))
  xl <- tolower(gsub('[^a-z0-9]+', '_', x0))
  
  out <- dplyr::case_when(
    xl %in% c('', 'na') ~ NA_character_,
    xl %in% c('anchor_informed', 'anchored', 'external_anchor', 'external_meta_anchor') ~ 'anchor_informed',
    xl %in% c('neutral_no_external_info', 'neutral_weakly_informative', 'neutral', 'no_external_info', 'neutral_no_external') ~ 'neutral_no_external_info',
    TRUE ~ x0
  )
  
  out
}

normalize_site_prior_family <- function(x) {
  x0 <- trimws(as.character(x))
  xl <- tolower(gsub('[^a-z0-9]+', '_', x0))
  
  out <- dplyr::case_when(
    xl %in% c('', 'na') ~ NA_character_,
    xl %in% c('normal_0_1_main', 'normal_main', 'gaussian_main') ~ 'normal_0_1_main',
    xl %in% c('student_t3_0_1_sensitivity', 'student_t3', 't3_sensitivity') ~ 'student_t3_0_1_sensitivity',
    TRUE ~ x0
  )
  
  out
}

normalize_model_family <- function(x) {
  x0 <- trimws(as.character(x))
  xl <- tolower(gsub('[^a-z0-9]+', '_', x0))
  
  out <- dplyr::case_when(
    xl %in% c('', 'na') ~ NA_character_,
    str_detect(xl, 'kaplan|^km$|benchmark_km') ~ 'KM',
    str_detect(xl, '^exp$|exponential') ~ 'exponential',
    str_detect(xl, 'weibull') ~ 'weibull',
    str_detect(xl, 'lognormal|log_normal') ~ 'lognormal',
    str_detect(xl, 'loglogistic|log_logistic') ~ 'loglogistic',
    str_detect(xl, 'cox') ~ 'cox',
    str_detect(xl, 'aft') ~ 'aft',
    str_detect(xl, 'aalen|\\baj\\b') ~ 'aalen_johansen',
    str_detect(xl, 'fine_gray|finegray') ~ 'fine_gray',
    TRUE ~ x0
  )
  
  out
}

normalize_branch <- function(x, default = NA_character_) {
  x0 <- trimws(as.character(x))
  xl <- tolower(gsub('[^a-z0-9]+', '_', x0))
  
  out <- dplyr::case_when(
    xl %in% c('', 'na') ~ default,
    str_detect(xl, '8a') & !str_detect(xl, '8b') ~ 'Stage8A',
    str_detect(xl, '8b|_cr|competing|remission') ~ 'Stage8B',
    str_detect(xl, 'stage9|^9$') ~ 'Stage9',
    TRUE ~ x0
  )
  
  out
}

infer_stage8_risk_scale_from_branch <- function(branch) {
  bl <- trim_lower(branch)
  
  dplyr::case_when(
    bl == 'stage8a' ~ 'transition_only_main',
    bl == 'stage8b' ~ 'transition_cif_competing',
    TRUE ~ NA_character_
  )
}

normalize_model_class <- function(x, default = NA_character_) {
  x0 <- trimws(as.character(x))
  xl <- tolower(gsub('[^a-z0-9]+', '_', x0))
  
  out <- dplyr::case_when(
    xl %in% c('', 'na') ~ default,
    str_detect(xl, 'km|benchmark') ~ 'KM_benchmark',
    str_detect(xl, 'no_cure|non_cure|noncure') ~ 'no_cure',
    str_detect(xl, 'frequentist_cure|mle_cure|stage7') ~ 'frequentist_cure',
    str_detect(xl, 'bayesian_cure|stage8|bayesian') ~ 'bayesian_cure',
    str_detect(xl, 'remission_sensitive_frequentist|stage9|competing_risk_frequentist') ~ 'remission_sensitive_frequentist',
    TRUE ~ x0
  )
  
  out
}

normalize_flag <- function(x) {
  if (is.logical(x)) return(x)
  
  xl <- trim_lower(x)
  
  out <- dplyr::case_when(
    xl %in% c('true', 't', '1', 'yes', 'y', 'admissible', 'retained', 'include', 'included', 'converged', 'ok', 'pass', 'supportive') ~ TRUE,
    xl %in% c('false', 'f', '0', 'no', 'n', 'non_admissible', 'excluded', 'exclude', 'dropped', 'fail', 'failed', 'not_retained', 'unsupportive') ~ FALSE,
    TRUE ~ NA
  )
  
  as.logical(out)
}

coerce_horizon_year <- function(df) {
  direct <- coalesce_cols(df, c('horizon_year', 'horizon', 'year', 'prediction_horizon', 'time_year'), default = NA)
  direct_num <- suppressWarnings(as.integer(round(coerce_numeric_text(direct))))
  
  if (!all(is.na(direct_num))) {
    return(direct_num)
  }
  
  horizon_id <- as.character(coalesce_cols(df, c('horizon_id'), default = NA_character_))
  suppressWarnings(as.integer(stringr::str_extract(horizon_id, '\\d+')))
}

make_fit_key <- function(df) {
  fields <- list(
    dataset_key = ifelse(is.na(df$dataset_key) | !nzchar(df$dataset_key), '__na__', df$dataset_key),
    model_class = ifelse(is.na(df$model_class) | !nzchar(df$model_class), '__na__', df$model_class),
    model_family = ifelse(is.na(df$model_family) | !nzchar(df$model_family), '__na__', df$model_family),
    branch = ifelse(is.na(df$branch) | !nzchar(df$branch), '__na__', df$branch),
    risk_scale = ifelse(is.na(df$risk_scale) | !nzchar(df$risk_scale), '__na__', df$risk_scale),
    prior_branch = ifelse(is.na(df$prior_branch) | !nzchar(df$prior_branch), '__na__', df$prior_branch),
    site_prior_family = ifelse(is.na(df$site_prior_family) | !nzchar(df$site_prior_family), '__na__', df$site_prior_family),
    retained_fit_id = ifelse(is.na(df$retained_fit_id) | !nzchar(df$retained_fit_id), '__na__', df$retained_fit_id)
  )
  
  do.call(paste, c(fields, sep = '||'))
}

make_row_key <- function(df) {
  fields <- list(
    fit_key = make_fit_key(df),
    horizon_year = ifelse(is.na(df$horizon_year), '__na__', as.character(df$horizon_year))
  )
  
  do.call(paste, c(fields, sep = '||'))
}

make_row_key_threshold <- function(df) {
  fields <- list(
    row_key = make_row_key(df),
    threshold = ifelse(is.na(df$threshold), '__na__', format(df$threshold, scientific = FALSE, trim = TRUE, digits = 8))
  )
  
  do.call(paste, c(fields, sep = '||'))
}

infer_model_class_from_stage <- function(stage_name) {
  dplyr::case_when(
    stage_name == 'Stage3' ~ 'KM_benchmark',
    stage_name == 'Stage5' ~ 'no_cure',
    stage_name == 'Stage7' ~ 'frequentist_cure',
    stage_name == 'Stage8' ~ 'bayesian_cure',
    stage_name == 'Stage9' ~ 'remission_sensitive_frequentist',
    TRUE ~ NA_character_
  )
}

infer_default_risk_scale_from_stage <- function(stage_name) {
  dplyr::case_when(
    stage_name %in% c('Stage3', 'Stage5', 'Stage7') ~ 'transition_only_main',
    stage_name == 'Stage9' ~ 'transition_cif_competing',
    TRUE ~ NA_character_
  )
}

infer_default_branch_from_stage <- function(stage_name) {
  dplyr::case_when(
    stage_name == 'Stage3' ~ 'Stage3KM',
    stage_name == 'Stage5' ~ 'Stage5',
    stage_name == 'Stage7' ~ 'Stage7',
    stage_name == 'Stage9' ~ 'Stage9',
    TRUE ~ NA_character_
  )
}

infer_significance_flag <- function(estimate, p_value, lower_ci, upper_ci, ratio_scale = FALSE) {
  if (!is.na(p_value)) {
    return(p_value < 0.05)
  }
  
  if (!is.na(lower_ci) && !is.na(upper_ci)) {
    if (ratio_scale) {
      return(lower_ci > 1 | upper_ci < 1)
    }
    return(lower_ci > 0 | upper_ci < 0)
  }
  
  if (!is.na(estimate)) {
    return(abs(estimate) > 0)
  }
  
  NA
}

# 🔴 Define: stage-directory scanners ===============================
## 🟠 Define: filename categorizer ===============================
categorize_csv <- function(file_name, headers) {
  fn <- tolower(basename(file_name))
  hd <- tolower(headers)
  
  if (header_has_any(hd, c('delta_risk_anchor_minus_neutral', 'delta_cure_fraction_anchor_minus_neutral', 'delta_fp100_anchor_minus_neutral'))) {
    return('anchor_delta')
  }
  
  if (header_has_any(hd, c('posterior_minus_prior_risk', 'prior_center_logit', 'posterior_mean_one_year_risk', 'age_sex_anchor_cell'))) {
    return('incidence_update')
  }
  
  if (header_has_any(hd, c('delta_risk_8b_minus_8a', 'delta_cure_fraction_8b_minus_8a', 'delta_fp100_8b_minus_8a'))) {
    return('remission_delta')
  }
  
  if (header_has_any(hd, c('lr_2delta_loglik', 'delta_aic_cure_minus_noncure', 'delta_bic_cure_minus_noncure', 'lrt_calibration_status'))) {
    return('fit_contrast')
  }
  
  if (header_has_any(hd, c('cure_model_eligibility_flag', 'receus_primary_class', 'followup_not_contradicted_flag', 'final_decision_flag'))) {
    return('carryforward')
  }
  
  if (header_has_any(hd, c('susceptible_fraction', 'cure_fraction', 'uncured_mean_support_flag', 'mstu', 'shape_class', 'merged_site_placement_label'))) {
    return('decomposition')
  }
  
  if (header_has_any(hd, c('percentage_method', 'cci', 'spt', 'late_horizon_instability_flag', 'risk_set_size')) &&
      name_has_any(fn, c('followup', 'maturity', 'support', 'completeness', 'instability'))) {
    return('followup')
  }
  
  if (header_has_any(hd, c('supported_short_horizon_contrast', 'intermediate_horizon_contrast', 'tail_diagnostic_contrast')) ||
      (name_has_any(fn, c('timing', 'contrast', 'restricted', 'interval')) &&
       header_has_any(hd, c('estimate', 'p_value', 'lower_ci', 'upper_ci', 'hazard_ratio')))) {
    return('timing')
  }
  
  if (header_has_any(hd, c('admissibility_flag', 'retained_flag', 'posterior_predictive_check_flag', 'convergence_status')) &&
      name_has_any(fn, c('diagnostic', 'diagnostics', 'admiss', 'retained', 'ppc', 'convergence'))) {
    return('diagnostics')
  }
  
  if (header_has_any(hd, c('net_benefit', 'net_reduction_in_unnecessary_intervention', 'net_reduction'))) {
    return('decision')
  }
  
  if (header_has_any(hd, c('positive_classification_rate', 'false_positive_count', 'false_positive_burden', 'false_positives_per100', 'fp100', 'ppv', 'tpr'))) {
    return('classification')
  }
  
  if (header_has_any(hd, c('auc_td', 'auc_t', 'auc', 'brier', 'ibs', 'calibration_slope', 'calibration_intercept', 'oe_ratio', 'c_index', 'uno_c', 'harrell_c'))) {
    return('performance')
  }
  
  if (header_has_any(hd, c('predicted_risk', 'pred_risk', 'risk', 'predicted_survival', 'pred_survival', 'survival', 'transition_cif')) &&
      header_has_any(hd, c('horizon_year', 'horizon', 'year', 'horizon_id'))) {
    return('prediction')
  }
  
  if (name_has_any(fn, c('prediction'))) return('prediction')
  if (name_has_any(fn, c('classification', 'false_positive', 'burden'))) return('classification')
  if (name_has_any(fn, c('performance', 'brier', 'auc', 'calibration'))) return('performance')
  if (name_has_any(fn, c('decision', 'net_benefit', 'dca'))) return('decision')
  if (name_has_any(fn, c('decomposition', 'uncured', 'susceptible', 'cure_fraction'))) return('decomposition')
  if (name_has_any(fn, c('anchor'))) return('anchor_delta')
  if (name_has_any(fn, c('incidence_update', 'prior_posterior'))) return('incidence_update')
  if (name_has_any(fn, c('delta', '8a', '8b', 'remission'))) return('remission_delta')
  if (name_has_any(fn, c('contrast', 'lrt', 'fit_gain'))) return('fit_contrast')
  if (name_has_any(fn, c('screen', 'carryforward', 'eligibility'))) return('carryforward')
  if (name_has_any(fn, c('followup', 'maturity', 'support'))) return('followup')
  if (name_has_any(fn, c('timing', 'restricted', 'interval'))) return('timing')
  if (name_has_any(fn, c('diagnostic', 'admiss'))) return('diagnostics')
  
  NA_character_
}

scan_stage_directory <- function(stage_dir, override_list = NULL) {
  if (!dir.exists(stage_dir)) {
    stop(sprintf('Stage directory does not exist: %s', stage_dir), call. = FALSE)
  }
  
  csv_files <- list.files(stage_dir, pattern = '\\.(csv|txt)$', full.names = TRUE, ignore.case = TRUE)
  
  if (length(csv_files) == 0L) {
    scanned <- tibble::tibble(
      path = character(),
      file_name = character(),
      headers = list(),
      category = character()
    )
  } else {
    scanned <- purrr::map_dfr(csv_files, function(path) {
      headers <- safe_csv_header(path)
      tibble::tibble(
        path = normalize_existing_path(path),
        file_name = basename(path),
        headers = list(headers),
        category = categorize_csv(path, headers)
      )
    })
  }
  
  if (!is.null(override_list)) {
    for (nm in names(override_list)) {
      override_path <- override_list[[nm]]
      if (!is.null(override_path) && nzchar(override_path)) {
        override_path <- normalize_existing_path(override_path)
        if (!file.exists(override_path)) {
          stop(sprintf('Override path does not exist for `%s`: %s', nm, override_path), call. = FALSE)
        }
        
        scanned <- scanned %>% filter(path != override_path)
        scanned <- bind_rows(
          scanned,
          tibble::tibble(
            path = override_path,
            file_name = basename(override_path),
            headers = list(safe_csv_header(override_path)),
            category = nm
          )
        )
      }
    }
  }
  
  scanned %>% arrange(category, file_name)
}

collect_category_paths <- function(scan_tbl, category_name) {
  if (nrow(scan_tbl) == 0L) return(character(0))
  scan_tbl %>%
    filter(category == category_name) %>%
    pull(path)
}

# 🔴 Define: canonical standardizers ===============================
## 🟠 Define: support metadata enricher ===============================
add_support_metadata <- function(df, horizon_registry) {
  if (nrow(df) == 0L) {
    return(df)
  }
  
  support_lookup <- horizon_registry %>%
    transmute(
      dataset_parent = normalize_dataset_key(dataset),
      horizon_year = as.integer(horizon_year),
      support_tier_registry = support_tier,
      horizon_evidence_class_registry = horizon_evidence_class,
      claim_restriction_flag_registry = claim_restriction_flag
    )
  
  df %>%
    left_join(support_lookup, by = c('dataset_parent', 'horizon_year')) %>%
    mutate(
      support_tier = dplyr::coalesce(support_tier, support_tier_registry),
      horizon_evidence_class = dplyr::coalesce(horizon_evidence_class, horizon_evidence_class_registry),
      claim_restriction_flag = dplyr::coalesce(
        claim_restriction_flag,
        ifelse(
          is.na(horizon_evidence_class),
          NA_character_,
          mapply(
            derive_claim_restriction_flag,
            horizon_evidence_class,
            support_tier,
            ifelse(is.na(prior_tail_sensitive), FALSE, prior_tail_sensitive),
            USE.NAMES = FALSE
          )
        ),
        claim_restriction_flag_registry
      )
    ) %>%
    select(-support_tier_registry, -horizon_evidence_class_registry, -claim_restriction_flag_registry)
}

standardize_base_fields <- function(df, stage_name, source_file) {
  branch_default_source <- if (stage_name == 'Stage8') basename(source_file) else infer_default_branch_from_stage(stage_name)
  
  out <- tibble::tibble(
    source_stage = stage_name,
    source_file = basename(source_file),
    dataset_key = normalize_dataset_key(coalesce_cols(df, c('dataset_key', 'dataset', 'dataset_name', 'analysis_dataset', 'cohort', 'analysis_scope'), default = NA_character_)),
    horizon_year = coerce_horizon_year(df),
    threshold = coerce_numeric_text(coalesce_cols(df, c('threshold', 'risk_threshold', 'threshold_value', 'c'), default = NA)),
    model_class = normalize_model_class(
      coalesce_cols(df, c('model_class', 'comparison_model_group', 'stage_model_group'), default = NA_character_),
      default = infer_model_class_from_stage(stage_name)
    ),
    model_family = normalize_model_family(coalesce_cols(df, c('model_family', 'family', 'latency_family', 'working_family', 'method', 'method_name'), default = NA_character_)),
    branch = normalize_branch(
      coalesce_cols(df, c('branch', 'stage_branch', 'analysis_branch', 'model_branch', 'risk_branch'), default = branch_default_source),
      default = infer_default_branch_from_stage(stage_name)
    ),
    risk_scale = normalize_risk_scale(
      coalesce_cols(df, c('risk_scale', 'scale', 'prediction_scale', 'risk_metric'), default = infer_default_risk_scale_from_stage(stage_name)),
      default = infer_default_risk_scale_from_stage(stage_name)
    ),
    prior_branch = normalize_prior_branch(coalesce_cols(df, c('prior_branch', 'prior_sensitivity_branch', 'anchor_branch'), default = NA_character_)),
    site_prior_family = normalize_site_prior_family(coalesce_cols(df, c('site_prior_family', 'site_prior', 'prior_family_site'), default = NA_character_)),
    retained_fit_id = as.character(coalesce_cols(df, c('retained_fit_id', 'fit_id', 'model_id', 'analysis_id', 'family_id'), default = NA_character_)),
    support_tier = as.character(coalesce_cols(df, c('support_tier', 'support_tier_standard'), default = NA_character_)),
    horizon_evidence_class = as.character(coalesce_cols(df, c('horizon_evidence_class'), default = NA_character_)),
    claim_restriction_flag = as.character(coalesce_cols(df, c('claim_restriction_flag'), default = NA_character_)),
    prior_tail_sensitive = normalize_flag(coalesce_cols(df, c('prior_tail_sensitive', 'tail_sensitive_flag'), default = NA)),
    admissibility_flag = normalize_flag(coalesce_cols(df, c('admissibility_flag', 'retained_flag', 'is_admissible'), default = NA))
  ) %>%
    mutate(
      dataset_key = dplyr::coalesce(dataset_key, ifelse(stage_name == 'Stage4', 'overall_project', NA_character_)),
      risk_scale = dplyr::coalesce(
        risk_scale,
        ifelse(stage_name == 'Stage8', infer_stage8_risk_scale_from_branch(branch), NA_character_)
      ),
      dataset_parent = dataset_parent_from_key(dataset_key),
      analysis_structure = derive_analysis_structure(dataset_key),
      reporting_priority = derive_reporting_priority(dataset_key)
    )
  
  out$fit_key <- make_fit_key(out)
  out$row_key <- make_row_key(out)
  out$row_key_threshold <- make_row_key_threshold(out)
  
  out
}

## 🟠 Define: table-specific standardizers ===============================
standardize_prediction_table <- function(df, stage_name, source_file, horizon_registry) {
  base <- standardize_base_fields(df, stage_name, source_file)
  
  id_raw <- as.character(coalesce_cols(df, c('person_id', 'unique_person_id', 'subject_id', 'id', 'uid'), default = NA_character_))
  site_raw <- as.character(coalesce_cols(df, c('site'), default = NA_character_))
  person_id <- ifelse(
    !is.na(id_raw) & nzchar(trimws(id_raw)) & !is.na(site_raw) & nzchar(trimws(site_raw)) & !str_detect(id_raw, '__'),
    paste(site_raw, id_raw, sep = '__'),
    id_raw
  )
  
  bind_cols(
    base,
    tibble::tibble(
      person_id = person_id,
      predicted_survival = coerce_numeric_text(coalesce_cols(df, c('predicted_survival', 'pred_survival', 'survival', 'posterior_survival'), default = NA)),
      predicted_risk = coerce_numeric_text(coalesce_cols(df, c('predicted_risk', 'pred_risk', 'risk', 'posterior_risk', 'km_risk'), default = NA)),
      transition_cif = coerce_numeric_text(coalesce_cols(df, c('transition_cif', 'cif', 'predicted_cif', 'posterior_cif'), default = NA)),
      lower_interval = coerce_numeric_text(coalesce_cols(df, c('lower_ci', 'lower_cl', 'ci_lower', 'lower_cri', 'credible_lower'), default = NA)),
      upper_interval = coerce_numeric_text(coalesce_cols(df, c('upper_ci', 'upper_cl', 'ci_upper', 'upper_cri', 'credible_upper'), default = NA)),
      interval_type = as.character(coalesce_cols(df, c('interval_type', 'ci_type'), default = NA_character_)),
      risk_set_size = coerce_numeric_text(coalesce_cols(df, c('risk_set_size', 'n_risk'), default = NA)),
      instability_marker = as.character(coalesce_cols(df, c('instability_marker', 'late_horizon_instability_flag', 'tail_instability_flag'), default = NA_character_))
    )
  ) %>%
    add_support_metadata(horizon_registry)
}

summarise_prediction_to_horizon <- function(pred_tbl) {
  if (nrow(pred_tbl) == 0L) {
    return(pred_tbl)
  }
  
  key_cols <- c(
    'source_stage', 'dataset_key', 'dataset_parent', 'analysis_structure', 'reporting_priority',
    'model_class', 'model_family', 'branch', 'risk_scale', 'prior_branch', 'site_prior_family',
    'retained_fit_id', 'support_tier', 'horizon_evidence_class', 'claim_restriction_flag',
    'prior_tail_sensitive', 'admissibility_flag', 'horizon_year', 'fit_key', 'row_key'
  )
  
  pred_tbl %>%
    group_by(across(all_of(key_cols))) %>%
    summarise(
      predicted_survival = if (all(is.na(predicted_survival))) NA_real_ else mean(predicted_survival, na.rm = TRUE),
      predicted_risk = if (all(is.na(predicted_risk))) NA_real_ else mean(predicted_risk, na.rm = TRUE),
      transition_cif = if (all(is.na(transition_cif))) NA_real_ else mean(transition_cif, na.rm = TRUE),
      lower_interval = if (all(is.na(lower_interval))) NA_real_ else mean(lower_interval, na.rm = TRUE),
      upper_interval = if (all(is.na(upper_interval))) NA_real_ else mean(upper_interval, na.rm = TRUE),
      risk_set_size = first_non_missing_scalar(risk_set_size),
      instability_marker = first_non_missing_scalar(instability_marker),
      prediction_rows = n(),
      .groups = 'drop'
    )
}

standardize_classification_table <- function(df, stage_name, source_file, horizon_registry) {
  base <- standardize_base_fields(df, stage_name, source_file)
  
  bind_cols(
    base,
    tibble::tibble(
      positive_classification_rate = coerce_numeric_text(coalesce_cols(df, c('positive_classification_rate', 'positive_rate', 'pos_rate'), default = NA)),
      false_positive_count = coerce_numeric_text(coalesce_cols(df, c('false_positive_count', 'fp_count'), default = NA)),
      false_positive_burden = coerce_numeric_text(coalesce_cols(df, c('false_positive_burden', 'fp_burden'), default = NA)),
      false_positives_per100 = coerce_numeric_text(coalesce_cols(df, c('false_positives_per100', 'fp100'), default = NA)),
      ppv = coerce_numeric_text(coalesce_cols(df, c('ppv', 'positive_predictive_value'), default = NA)),
      tpr = coerce_numeric_text(coalesce_cols(df, c('tpr', 'sensitivity'), default = NA)),
      risk_set_size = coerce_numeric_text(coalesce_cols(df, c('risk_set_size', 'n_risk'), default = NA)),
      instability_marker = as.character(coalesce_cols(df, c('instability_marker', 'late_horizon_instability_flag', 'tail_instability_flag'), default = NA_character_))
    )
  ) %>%
    add_support_metadata(horizon_registry)
}

standardize_decision_table <- function(df, stage_name, source_file, horizon_registry) {
  base <- standardize_base_fields(df, stage_name, source_file)
  
  bind_cols(
    base,
    tibble::tibble(
      net_benefit = coerce_numeric_text(coalesce_cols(df, c('net_benefit', 'nb'), default = NA)),
      net_reduction_in_unnecessary_intervention = coerce_numeric_text(coalesce_cols(df, c('net_reduction_in_unnecessary_intervention', 'net_reduction'), default = NA)),
      positive_classification_rate = coerce_numeric_text(coalesce_cols(df, c('positive_classification_rate', 'positive_rate', 'pos_rate'), default = NA)),
      risk_set_size = coerce_numeric_text(coalesce_cols(df, c('risk_set_size', 'n_risk'), default = NA)),
      instability_marker = as.character(coalesce_cols(df, c('instability_marker', 'late_horizon_instability_flag', 'tail_instability_flag'), default = NA_character_))
    )
  ) %>%
    add_support_metadata(horizon_registry)
}

standardize_performance_table <- function(df, stage_name, source_file, horizon_registry) {
  base <- standardize_base_fields(df, stage_name, source_file)
  
  bind_cols(
    base,
    tibble::tibble(
      predicted_survival = coerce_numeric_text(coalesce_cols(df, c('predicted_survival', 'survival'), default = NA)),
      predicted_risk = coerce_numeric_text(coalesce_cols(df, c('predicted_risk', 'risk', 'km_risk'), default = NA)),
      transition_cif = coerce_numeric_text(coalesce_cols(df, c('transition_cif', 'cif'), default = NA)),
      auc_td = coerce_numeric_text(coalesce_cols(df, c('auc_td', 'auc_t', 'auc'), default = NA)),
      c_index = coerce_numeric_text(coalesce_cols(df, c('c_index', 'uno_c', 'harrell_c'), default = NA)),
      calibration_intercept = coerce_numeric_text(coalesce_cols(df, c('calibration_intercept'), default = NA)),
      calibration_slope = coerce_numeric_text(coalesce_cols(df, c('calibration_slope'), default = NA)),
      oe_ratio = coerce_numeric_text(coalesce_cols(df, c('oe_ratio', 'o_e'), default = NA)),
      brier = coerce_numeric_text(coalesce_cols(df, c('brier', 'brier_score'), default = NA)),
      ibs = coerce_numeric_text(coalesce_cols(df, c('ibs', 'integrated_brier_score'), default = NA)),
      lower_interval = coerce_numeric_text(coalesce_cols(df, c('lower_ci', 'lower_cri'), default = NA)),
      upper_interval = coerce_numeric_text(coalesce_cols(df, c('upper_ci', 'upper_cri'), default = NA)),
      risk_set_size = coerce_numeric_text(coalesce_cols(df, c('risk_set_size', 'n_risk'), default = NA)),
      instability_marker = as.character(coalesce_cols(df, c('instability_marker', 'late_horizon_instability_flag', 'tail_instability_flag'), default = NA_character_))
    )
  ) %>%
    add_support_metadata(horizon_registry)
}

standardize_followup_table <- function(df, source_file, horizon_registry) {
  base <- standardize_base_fields(df, 'Stage2', source_file)
  
  bind_cols(
    base,
    tibble::tibble(
      risk_set_size = coerce_numeric_text(coalesce_cols(df, c('risk_set_size', 'n_risk'), default = NA)),
      cumulative_event_count = coerce_numeric_text(coalesce_cols(df, c('cumulative_event_count', 'cum_event_count'), default = NA)),
      cumulative_censoring_count = coerce_numeric_text(coalesce_cols(df, c('cumulative_censoring_count', 'cum_censor_count'), default = NA)),
      percentage_method = as.character(coalesce_cols(df, c('percentage_method'), default = NA_character_)),
      CCI = coerce_numeric_text(coalesce_cols(df, c('CCI', 'cci'), default = NA)),
      SPT = coerce_numeric_text(coalesce_cols(df, c('SPT', 'spt'), default = NA)),
      late_horizon_instability_flag = normalize_flag(coalesce_cols(df, c('late_horizon_instability_flag', 'tail_instability_flag'), default = NA)),
      site_dominance_warning = as.character(coalesce_cols(df, c('site_dominance_warning'), default = NA_character_))
    )
  ) %>%
    add_support_metadata(horizon_registry)
}

standardize_timing_table <- function(df, source_file) {
  base <- standardize_base_fields(df, 'Stage4', source_file)
  ratio_scale_flag <- header_has_any(names(df), c('hazard_ratio', 'hr', 'odds_ratio'))
  
  bind_cols(
    base,
    tibble::tibble(
      contrast_label = as.character(coalesce_cols(df, c('contrast_label', 'output_label', 'analysis_label', 'window_label', 'interval_label', 'contrast_type'), default = NA_character_)),
      estimate = coerce_numeric_text(coalesce_cols(df, c('estimate', 'difference', 'delta', 'beta', 'log_hr', 'hazard_ratio'), default = NA)),
      p_value = coerce_numeric_text(coalesce_cols(df, c('p_value', 'pvalue', 'p'), default = NA)),
      lower_ci = coerce_numeric_text(coalesce_cols(df, c('lower_ci', 'ci_lower', 'lcl'), default = NA)),
      upper_ci = coerce_numeric_text(coalesce_cols(df, c('upper_ci', 'ci_upper', 'ucl'), default = NA)),
      significance_flag = normalize_flag(coalesce_cols(df, c('significance_flag', 'is_significant'), default = NA))
    )
  ) %>%
    mutate(
      significance_flag = ifelse(
        is.na(significance_flag),
        mapply(infer_significance_flag, estimate, p_value, lower_ci, upper_ci, MoreArgs = list(ratio_scale = ratio_scale_flag)),
        significance_flag
      )
    )
}

standardize_carryforward_table <- function(df, source_file) {
  tibble::tibble(
    source_file = basename(source_file),
    dataset_key = normalize_dataset_key(coalesce_cols(df, c('dataset_key', 'dataset', 'dataset_name'), default = NA_character_)),
    cure_model_eligibility_flag = as.character(coalesce_cols(df, c('cure_model_eligibility_flag', 'final_decision_flag', 'stage6_final_class'), default = NA_character_)),
    primary_gate_method = as.character(coalesce_cols(df, c('primary_gate_method'), default = NA_character_)),
    primary_gate_flag = normalize_flag(coalesce_cols(df, c('primary_gate_flag'), default = NA)),
    receus_primary_class = as.character(coalesce_cols(df, c('receus_primary_class'), default = NA_character_)),
    presence_modifier_flag = as.character(coalesce_cols(df, c('presence_modifier_flag'), default = NA_character_)),
    cure_presence_support_flag = normalize_flag(coalesce_cols(df, c('cure_presence_support_flag'), default = NA)),
    followup_contradiction_flag = normalize_flag(coalesce_cols(df, c('followup_contradiction_flag'), default = NA)),
    followup_not_contradicted_flag = normalize_flag(coalesce_cols(df, c('followup_not_contradicted_flag'), default = NA)),
    descriptive_tail_summary_flag = as.character(coalesce_cols(df, c('descriptive_tail_summary_flag'), default = NA_character_)),
    supporting_methods = as.character(coalesce_cols(df, c('supporting_methods'), default = NA_character_)),
    contradicting_methods = as.character(coalesce_cols(df, c('contradicting_methods'), default = NA_character_)),
    screening_note = as.character(coalesce_cols(df, c('screening_note', 'screening_context_note'), default = NA_character_)),
    common_horizon_vector = as.character(coalesce_cols(df, c('common_horizon_vector'), default = NA_character_)),
    common_threshold_vector = as.character(coalesce_cols(df, c('common_threshold_vector'), default = NA_character_)),
    screening_context = as.character(coalesce_cols(df, c('screening_context'), default = NA_character_))
  ) %>%
    mutate(
      dataset_parent = dataset_parent_from_key(dataset_key),
      carry_forward_stage8 = ifelse(is.na(cure_model_eligibility_flag), NA, cure_model_eligibility_flag != 'unsupportive')
    )
}

standardize_decomposition_table <- function(df, stage_name, source_file, horizon_registry) {
  base <- standardize_base_fields(df, stage_name, source_file)
  
  extra_cols <- names(df)[str_detect(
    names(df),
    regex('^(uncured|susceptible)_(survival|risk)_(\\d+y|year_\\d+)$|^site_in_(neither|incidence_only|latency_only|both)$', ignore_case = TRUE)
  )]
  
  extra_tbl <- if (length(extra_cols) > 0L) {
    tibble::as_tibble(df[, extra_cols, drop = FALSE])
  } else {
    tibble::tibble()
  }
  
  bind_cols(
    base,
    tibble::tibble(
      incidence_coefficients = as.character(coalesce_cols(df, c('incidence_coefficients', 'incidence_coefs'), default = NA_character_)),
      latency_coefficients = as.character(coalesce_cols(df, c('latency_coefficients', 'latency_coefs'), default = NA_character_)),
      cure_fraction = coerce_numeric_text(coalesce_cols(df, c('cure_fraction', 'posterior_cure_fraction'), default = NA)),
      susceptible_fraction = coerce_numeric_text(coalesce_cols(df, c('susceptible_fraction'), default = NA)),
      uncured_mean_support_flag = normalize_flag(coalesce_cols(df, c('uncured_mean_support_flag'), default = NA)),
      MSTu = coerce_numeric_text(coalesce_cols(df, c('MSTu', 'mstu'), default = NA)),
      merged_site_placement_label = as.character(coalesce_cols(df, c('merged_site_placement_label', 'site_placement_label', 'site_structure_label'), default = NA_character_)),
      hazard_target = as.character(coalesce_cols(df, c('hazard_target'), default = NA_character_)),
      hazard_1y = coerce_numeric_text(coalesce_cols(df, c('hazard_1y'), default = NA)),
      hazard_2y = coerce_numeric_text(coalesce_cols(df, c('hazard_2y'), default = NA)),
      hazard_3y = coerce_numeric_text(coalesce_cols(df, c('hazard_3y'), default = NA)),
      hazard_4y = coerce_numeric_text(coalesce_cols(df, c('hazard_4y'), default = NA)),
      hazard_5y = coerce_numeric_text(coalesce_cols(df, c('hazard_5y'), default = NA)),
      hazard_6y = coerce_numeric_text(coalesce_cols(df, c('hazard_6y'), default = NA)),
      hazard_7y = coerce_numeric_text(coalesce_cols(df, c('hazard_7y'), default = NA)),
      hazard_8y = coerce_numeric_text(coalesce_cols(df, c('hazard_8y'), default = NA)),
      hazard_9y = coerce_numeric_text(coalesce_cols(df, c('hazard_9y'), default = NA)),
      hazard_10y = coerce_numeric_text(coalesce_cols(df, c('hazard_10y'), default = NA)),
      hazard_ratio_10y_vs_1y = coerce_numeric_text(coalesce_cols(df, c('hazard_ratio_10y_vs_1y'), default = NA)),
      shape_class = as.character(coalesce_cols(df, c('shape_class'), default = NA_character_))
    ),
    extra_tbl
  ) %>%
    mutate(
      risk_scale = dplyr::coalesce(
        risk_scale,
        normalize_risk_scale(coalesce_cols(df, c('risk_scale'), default = infer_default_risk_scale_from_stage(stage_name)), default = infer_default_risk_scale_from_stage(stage_name))
      )
    ) %>%
    add_support_metadata(horizon_registry)
}

standardize_fit_contrast_table <- function(df, source_file) {
  base <- standardize_base_fields(df, 'Stage7', source_file)
  
  bind_cols(
    base,
    tibble::tibble(
      family_pair = as.character(coalesce_cols(df, c('family_pair'), default = NA_character_)),
      logLik_noncure = coerce_numeric_text(coalesce_cols(df, c('logLik_noncure'), default = NA)),
      logLik_cure = coerce_numeric_text(coalesce_cols(df, c('logLik_cure'), default = NA)),
      delta_logLik_cure_minus_noncure = coerce_numeric_text(coalesce_cols(df, c('delta_logLik_cure_minus_noncure'), default = NA)),
      LR_2delta_logLik = coerce_numeric_text(coalesce_cols(df, c('LR_2delta_logLik'), default = NA)),
      delta_AIC_cure_minus_noncure = coerce_numeric_text(coalesce_cols(df, c('delta_AIC_cure_minus_noncure'), default = NA)),
      delta_BIC_cure_minus_noncure = coerce_numeric_text(coalesce_cols(df, c('delta_BIC_cure_minus_noncure'), default = NA)),
      lrt_calibration_status = as.character(coalesce_cols(df, c('lrt_calibration_status'), default = NA_character_)),
      lrt_pvalue_bootstrap = coerce_numeric_text(coalesce_cols(df, c('lrt_pvalue_bootstrap'), default = NA)),
      same_family_fit_gain_signal = as.character(coalesce_cols(df, c('same_family_fit_gain_signal'), default = NA_character_)),
      convergence_pair_flag = as.character(coalesce_cols(df, c('convergence_pair_flag'), default = NA_character_))
    )
  )
}

standardize_anchor_delta_table <- function(df, source_file, horizon_registry) {
  base <- standardize_base_fields(df, 'Stage8', source_file)
  
  bind_cols(
    base,
    tibble::tibble(
      delta_risk_anchor_minus_neutral = coerce_numeric_text(coalesce_cols(df, c('delta_risk_anchor_minus_neutral'), default = NA)),
      delta_cure_fraction_anchor_minus_neutral = coerce_numeric_text(coalesce_cols(df, c('delta_cure_fraction_anchor_minus_neutral'), default = NA)),
      delta_false_positive_burden_anchor_minus_neutral = coerce_numeric_text(coalesce_cols(df, c('delta_false_positive_burden_anchor_minus_neutral'), default = NA)),
      delta_FP100_anchor_minus_neutral = coerce_numeric_text(coalesce_cols(df, c('delta_FP100_anchor_minus_neutral'), default = NA)),
      delta_NB_anchor_minus_neutral = coerce_numeric_text(coalesce_cols(df, c('delta_NB_anchor_minus_neutral'), default = NA)),
      delta_PPV_anchor_minus_neutral = coerce_numeric_text(coalesce_cols(df, c('delta_PPV_anchor_minus_neutral'), default = NA)),
      delta_TPR_anchor_minus_neutral = coerce_numeric_text(coalesce_cols(df, c('delta_TPR_anchor_minus_neutral'), default = NA))
    )
  ) %>%
    add_support_metadata(horizon_registry)
}

standardize_incidence_update_table <- function(df, source_file) {
  base <- standardize_base_fields(df, 'Stage8', source_file)
  
  bind_cols(
    base,
    tibble::tibble(
      age_sex_anchor_cell = as.character(coalesce_cols(df, c('age_sex_anchor_cell'), default = NA_character_)),
      external_incidence_rate_per10k = coerce_numeric_text(coalesce_cols(df, c('external_incidence_rate_per10k'), default = NA)),
      external_one_year_risk = coerce_numeric_text(coalesce_cols(df, c('external_one_year_risk'), default = NA)),
      prior_center_logit = coerce_numeric_text(coalesce_cols(df, c('prior_center_logit'), default = NA)),
      posterior_mean_logit = coerce_numeric_text(coalesce_cols(df, c('posterior_mean_logit'), default = NA)),
      posterior_lower_logit = coerce_numeric_text(coalesce_cols(df, c('posterior_lower_logit'), default = NA)),
      posterior_upper_logit = coerce_numeric_text(coalesce_cols(df, c('posterior_upper_logit'), default = NA)),
      posterior_mean_one_year_risk = coerce_numeric_text(coalesce_cols(df, c('posterior_mean_one_year_risk'), default = NA)),
      posterior_lower_one_year_risk = coerce_numeric_text(coalesce_cols(df, c('posterior_lower_one_year_risk'), default = NA)),
      posterior_upper_one_year_risk = coerce_numeric_text(coalesce_cols(df, c('posterior_upper_one_year_risk'), default = NA)),
      posterior_minus_prior_logit = coerce_numeric_text(coalesce_cols(df, c('posterior_minus_prior_logit'), default = NA)),
      posterior_minus_prior_risk = coerce_numeric_text(coalesce_cols(df, c('posterior_minus_prior_risk'), default = NA))
    )
  )
}

standardize_remission_delta_table <- function(df, source_file, horizon_registry) {
  base <- standardize_base_fields(df, 'Stage8', source_file)
  
  bind_cols(
    base,
    tibble::tibble(
      delta_risk_8B_minus_8A = coerce_numeric_text(coalesce_cols(df, c('delta_risk_8B_minus_8A'), default = NA)),
      delta_cure_fraction_8B_minus_8A = coerce_numeric_text(coalesce_cols(df, c('delta_cure_fraction_8B_minus_8A'), default = NA)),
      delta_false_positive_burden_8B_minus_8A = coerce_numeric_text(coalesce_cols(df, c('delta_false_positive_burden_8B_minus_8A'), default = NA)),
      delta_FP100_8B_minus_8A = coerce_numeric_text(coalesce_cols(df, c('delta_FP100_8B_minus_8A'), default = NA)),
      delta_NB_8B_minus_8A = coerce_numeric_text(coalesce_cols(df, c('delta_NB_8B_minus_8A'), default = NA)),
      delta_PPV_8B_minus_8A = coerce_numeric_text(coalesce_cols(df, c('delta_PPV_8B_minus_8A'), default = NA)),
      delta_TPR_8B_minus_8A = coerce_numeric_text(coalesce_cols(df, c('delta_TPR_8B_minus_8A'), default = NA))
    )
  ) %>%
    add_support_metadata(horizon_registry)
}

standardize_diagnostics_table <- function(df, source_file) {
  base <- standardize_base_fields(df, 'Stage8', source_file)
  
  base %>%
    mutate(
      convergence_status = as.character(coalesce_cols(df, c('convergence_status'), default = NA_character_)),
      posterior_predictive_check_flag = normalize_flag(coalesce_cols(df, c('posterior_predictive_check_flag'), default = NA)),
      prior_sensitivity_ok_flag = normalize_flag(coalesce_cols(df, c('prior_sensitivity_ok_flag'), default = NA)),
      admissibility_flag = dplyr::coalesce(
        normalize_flag(coalesce_cols(df, c('admissibility_flag', 'retained_flag', 'is_admissible'), default = NA)),
        admissibility_flag
      )
    )
}

standardize_stage9_delta_from_main <- function(df, source_file, horizon_registry) {
  base <- standardize_base_fields(df, 'Stage9', source_file)
  
  bind_cols(
    base,
    tibble::tibble(
      delta_risk_9_minus_main = coerce_numeric_text(coalesce_cols(df, c('absolute_difference_from_main_risk', 'delta_risk_9_minus_main'), default = NA)),
      delta_false_positive_burden_9_minus_main = coerce_numeric_text(coalesce_cols(df, c('change_false_positive_burden', 'delta_false_positive_burden_9_minus_main'), default = NA)),
      delta_NB_9_minus_main = coerce_numeric_text(coalesce_cols(df, c('change_net_benefit', 'delta_NB_9_minus_main'), default = NA))
    )
  ) %>%
    add_support_metadata(horizon_registry)
}

ensure_columns <- function(df, col_defaults) {
  for (nm in names(col_defaults)) {
    if (!nm %in% names(df)) {
      val <- col_defaults[[nm]]
      
      if (length(val) == 1L) {
        df[[nm]] <- rep(val, nrow(df))
      } else if (length(val) == 0L) {
        template <- if (is.character(val)) {
          NA_character_
        } else if (is.integer(val)) {
          NA_integer_
        } else if (is.numeric(val)) {
          NA_real_
        } else if (is.logical(val)) {
          NA
        } else {
          NA
        }
        df[[nm]] <- rep(template, nrow(df))
      } else {
        df[[nm]] <- val
      }
    }
  }
  df
}

# 🔴 Define: stage assemblers and diagnostics join ===============================
## 🟠 Define: metric table assemblers ===============================
collapse_by_key <- function(df, key_col) {
  if (nrow(df) == 0L) {
    return(df)
  }
  
  df %>%
    group_by(.data[[key_col]]) %>%
    summarise(across(everything(), first_non_missing_scalar), .groups = 'drop')
}

filter_metric_rows <- function(df, metric_cols) {
  if (nrow(df) == 0L) {
    return(df)
  }
  
  use_cols <- intersect(metric_cols, names(df))
  if (length(use_cols) == 0L) {
    return(df[0, , drop = FALSE])
  }
  
  keep_flag <- Reduce(
    `|`,
    lapply(use_cols, function(nm) {
      x <- df[[nm]]
      if (is.character(x)) {
        !is.na(x) & nzchar(trimws(x))
      } else {
        !is.na(x)
      }
    })
  )
  
  df[keep_flag, , drop = FALSE]
}

ensure_threshold_col <- function(df) {
  if (!'threshold' %in% names(df)) {
    df$threshold <- NA_real_
  }
  df
}

assemble_metric_stage <- function(stage_name, scan_tbl, horizon_registry) {
  pred_paths <- collect_category_paths(scan_tbl, 'prediction')
  class_paths <- collect_category_paths(scan_tbl, 'classification')
  perf_paths <- collect_category_paths(scan_tbl, 'performance')
  dec_paths <- collect_category_paths(scan_tbl, 'decision')
  
  metric_paths <- unique(c(pred_paths, class_paths, perf_paths, dec_paths))
  if (length(metric_paths) == 0L) {
    return(
      list(
        prediction_table = tibble::tibble(),
        prediction_horizon_table = tibble::tibble(),
        classification_performance_table = tibble::tibble()
      )
    )
  }
  
  prediction_source_paths <- if (length(pred_paths) > 0L) unique(pred_paths) else metric_paths
  
  pred_tbl <- purrr::map_dfr(
    prediction_source_paths,
    ~ standardize_prediction_table(read_table_any(.x), stage_name, .x, horizon_registry)
  ) %>%
    filter_metric_rows(c('predicted_survival', 'predicted_risk', 'transition_cif', 'lower_interval', 'upper_interval'))
  
  pred_horizon_tbl <- summarise_prediction_to_horizon(pred_tbl)
  
  class_tbl <- purrr::map_dfr(
    metric_paths,
    ~ standardize_classification_table(read_table_any(.x), stage_name, .x, horizon_registry)
  ) %>%
    filter_metric_rows(c('positive_classification_rate', 'false_positive_count', 'false_positive_burden', 'false_positives_per100', 'ppv', 'tpr')) %>%
    ensure_threshold_col()
  
  perf_tbl <- purrr::map_dfr(
    metric_paths,
    ~ standardize_performance_table(read_table_any(.x), stage_name, .x, horizon_registry)
  ) %>%
    filter_metric_rows(c('predicted_survival', 'predicted_risk', 'transition_cif', 'auc_td', 'c_index', 'calibration_intercept', 'calibration_slope', 'oe_ratio', 'brier', 'ibs')) %>%
    ensure_threshold_col()
  
  dec_tbl <- purrr::map_dfr(
    metric_paths,
    ~ standardize_decision_table(read_table_any(.x), stage_name, .x, horizon_registry)
  ) %>%
    filter_metric_rows(c('net_benefit', 'net_reduction_in_unnecessary_intervention', 'positive_classification_rate')) %>%
    ensure_threshold_col()
  
  threshold_stack <- bind_rows(class_tbl, dec_tbl) %>%
    filter(!is.na(threshold))
  
  horizon_stack <- bind_rows(
    class_tbl %>% filter(is.na(threshold)),
    perf_tbl %>% filter(is.na(threshold)),
    dec_tbl %>% filter(is.na(threshold)),
    pred_horizon_tbl
  )
  
  threshold_merged <- collapse_by_key(threshold_stack, 'row_key_threshold')
  horizon_merged <- collapse_by_key(horizon_stack, 'row_key')
  
  if (nrow(threshold_merged) == 0L && nrow(horizon_merged) == 0L) {
    unified_tbl <- tibble::tibble()
  } else if (nrow(threshold_merged) > 0L && nrow(horizon_merged) > 0L) {
    horizon_suffix <- horizon_merged %>% select(-threshold, -row_key_threshold)
    unified_tbl <- threshold_merged %>%
      left_join(horizon_suffix, by = 'row_key', suffix = c('', '_h'))
  } else if (nrow(threshold_merged) > 0L) {
    unified_tbl <- threshold_merged
  } else {
    unified_tbl <- horizon_merged %>%
      mutate(
        threshold = NA_real_,
        row_key_threshold = make_row_key_threshold(.)
      )
  }
  
  list(
    prediction_table = pred_tbl,
    prediction_horizon_table = pred_horizon_tbl,
    classification_performance_table = unified_tbl
  )
}

apply_stage8_diagnostics <- function(df, diagnostics_tbl) {
  if (nrow(df) == 0L || nrow(diagnostics_tbl) == 0L) {
    return(df)
  }
  
  if (!'admissibility_flag' %in% names(df)) df$admissibility_flag <- NA
  if (!'convergence_status' %in% names(df)) df$convergence_status <- NA_character_
  if (!'posterior_predictive_check_flag' %in% names(df)) df$posterior_predictive_check_flag <- NA
  if (!'prior_sensitivity_ok_flag' %in% names(df)) df$prior_sensitivity_ok_flag <- NA
  
  diag_lookup <- diagnostics_tbl %>%
    select(fit_key, admissibility_flag, convergence_status, posterior_predictive_check_flag, prior_sensitivity_ok_flag) %>%
    distinct()
  
  df %>%
    left_join(diag_lookup, by = 'fit_key', suffix = c('', '_diag')) %>%
    mutate(
      admissibility_flag = dplyr::coalesce(admissibility_flag, admissibility_flag_diag),
      convergence_status = dplyr::coalesce(convergence_status, convergence_status_diag),
      posterior_predictive_check_flag = dplyr::coalesce(posterior_predictive_check_flag, posterior_predictive_check_flag_diag),
      prior_sensitivity_ok_flag = dplyr::coalesce(prior_sensitivity_ok_flag, prior_sensitivity_ok_flag_diag)
    ) %>%
    select(-admissibility_flag_diag, -convergence_status_diag, -posterior_predictive_check_flag_diag, -prior_sensitivity_ok_flag_diag)
}

# 🔴 Define: interpretation and triangulation helpers ===============================
## 🟠 Define: synthesized signal helpers ===============================
derive_timing_signal <- function(timing_tbl) {
  if (nrow(timing_tbl) == 0L) {
    return('timing_signal_unavailable')
  }
  
  lab <- trim_lower(timing_tbl$contrast_label)
  sig <- timing_tbl$significance_flag
  
  early_flag <- any(sig %in% TRUE & (str_detect(lab, 'supported_short|1_year|2_year|0_1|1_2') | str_detect(lab, 'restricted')), na.rm = TRUE)
  mid_flag <- any(sig %in% TRUE & str_detect(lab, 'intermediate|2_5'), na.rm = TRUE)
  tail_flag <- any(sig %in% TRUE & str_detect(lab, 'tail|>5|5_'), na.rm = TRUE)
  
  if (early_flag) return('early_supported_difference')
  if (!early_flag && (mid_flag || tail_flag)) return('late_only_difference')
  'weak_or_absent_early_difference'
}

derive_followup_signal <- function(dataset_key, followup_tbl, horizon_registry) {
  target_parent <- dataset_parent_from_key(dataset_key)
  
  target <- followup_tbl %>%
    filter(dataset_parent == target_parent, !is.na(horizon_year)) %>%
    arrange(horizon_year)
  
  if (nrow(target) == 0L) {
    target <- horizon_registry %>%
      transmute(
        dataset_parent = normalize_dataset_key(dataset),
        horizon_year = as.integer(horizon_year),
        support_tier = support_tier,
        horizon_evidence_class = horizon_evidence_class
      ) %>%
      filter(dataset_parent == target_parent)
  }
  
  direct_primary <- target %>%
    filter(support_tier == 'primary_supported', horizon_evidence_class == 'directly_observed_data_supported')
  
  if (nrow(direct_primary) > 0L) return('directly_observed_data_supported')
  
  partial_primary <- target %>%
    filter(support_tier %in% c('primary_supported', 'secondary', 'sensitivity'), horizon_evidence_class == 'partly_model_dependent')
  
  if (nrow(partial_primary) > 0L) return('partly_model_dependent')
  
  'mostly_extrapolated'
}

derive_anchor_dependence_signal <- function(anchor_tbl) {
  if (nrow(anchor_tbl) == 0L) {
    return('anchor_dependence_unavailable')
  }
  
  max_abs <- max(
    c(
      abs(anchor_tbl$delta_risk_anchor_minus_neutral),
      abs(anchor_tbl$delta_false_positive_burden_anchor_minus_neutral),
      abs(anchor_tbl$delta_FP100_anchor_minus_neutral),
      abs(anchor_tbl$delta_NB_anchor_minus_neutral)
    ),
    na.rm = TRUE
  )
  
  if (!is.finite(max_abs)) return('anchor_dependence_unavailable')
  if (max_abs >= 0.03) return('high_anchor_dependence')
  if (max_abs >= 0.01) return('moderate_anchor_dependence')
  'low_anchor_dependence'
}

derive_same_family_fit_gain_signal <- function(fit_contrast_tbl) {
  if (nrow(fit_contrast_tbl) == 0L) {
    return('fit_gain_unavailable')
  }
  
  if ('same_family_fit_gain_signal' %in% names(fit_contrast_tbl)) {
    explicit <- fit_contrast_tbl$same_family_fit_gain_signal
    explicit <- explicit[!is.na(explicit) & nzchar(trimws(explicit))]
    if (length(explicit) > 0L) return(explicit[[1L]])
  }
  
  supportive_n <- fit_contrast_tbl %>%
    filter(!is.na(delta_AIC_cure_minus_noncure), delta_AIC_cure_minus_noncure < 0) %>%
    nrow()
  
  negative_n <- fit_contrast_tbl %>%
    filter(!is.na(delta_AIC_cure_minus_noncure), delta_AIC_cure_minus_noncure > 0) %>%
    nrow()
  
  if (supportive_n >= 2L && negative_n == 0L) return('consistent_gain')
  if (supportive_n >= 1L && negative_n >= 1L) return('mixed_gain')
  if (supportive_n >= 1L) return('limited_gain')
  'no_clear_gain'
}

derive_uncured_only_support_signal <- function(decomp_tbl) {
  if (nrow(decomp_tbl) == 0L) {
    return('uncured_support_unavailable')
  }
  
  if ('uncured_mean_support_flag' %in% names(decomp_tbl) && any(decomp_tbl$uncured_mean_support_flag %in% TRUE, na.rm = TRUE)) {
    return('supported')
  }
  
  hazard_cols <- intersect(names(decomp_tbl), c('hazard_1y', 'hazard_2y', 'hazard_3y', 'hazard_4y', 'hazard_5y', 'hazard_6y', 'hazard_7y', 'hazard_8y', 'hazard_9y', 'hazard_10y'))
  grid_cols <- names(decomp_tbl)[str_detect(names(decomp_tbl), regex('^(uncured|susceptible)_(survival|risk)_(\\d+y|year_\\d+)$', ignore_case = TRUE))]
  use_cols <- unique(c(hazard_cols, grid_cols))
  
  if (length(use_cols) == 0L) {
    return('limited')
  }
  
  uncured_grid_present <- any(rowSums(!is.na(decomp_tbl[, use_cols, drop = FALSE])) > 0, na.rm = TRUE)
  
  if (uncured_grid_present) return('annual_grid_only')
  'limited'
}

derive_latency_family_signal <- function(main_tbl) {
  if (nrow(main_tbl) == 0L) {
    return('family_signal_unavailable')
  }
  
  fam_n <- main_tbl %>%
    filter(
      risk_scale == 'transition_only_main',
      support_tier %in% c('primary_supported', 'secondary', 'sensitivity'),
      model_class %in% c('no_cure', 'frequentist_cure', 'bayesian_cure'),
      !is.na(model_family),
      nzchar(model_family)
    ) %>%
    distinct(model_class, model_family) %>%
    nrow()
  
  if (fam_n >= 4L) return('stable_multi_family_signal')
  if (fam_n >= 2L) return('limited_multi_family_signal')
  if (fam_n >= 1L) return('single_family_only')
  'family_signal_unavailable'
}

derive_remission_distortion_signal <- function(delta_tbl) {
  if (nrow(delta_tbl) == 0L) {
    return('remission_signal_unavailable')
  }
  
  max_abs <- max(
    c(
      abs(delta_tbl$delta_risk_8B_minus_8A),
      abs(delta_tbl$delta_false_positive_burden_8B_minus_8A),
      abs(delta_tbl$delta_FP100_8B_minus_8A),
      abs(delta_tbl$delta_NB_8B_minus_8A),
      abs(delta_tbl$delta_risk_9_minus_main),
      abs(delta_tbl$delta_false_positive_burden_9_minus_main),
      abs(delta_tbl$delta_NB_9_minus_main)
    ),
    na.rm = TRUE
  )
  
  if (!is.finite(max_abs)) return('low')
  if (max_abs >= 0.03) return('high')
  if (max_abs >= 0.01) return('moderate')
  'low'
}

derive_interpretation_bucket <- function(axis_timing_difference, axis_followup_support, axis_cure_support, axis_remission_distortion) {
  if (axis_remission_distortion == 'high') {
    return('remission_distorted_apparent_cure')
  }
  
  if (axis_timing_difference == 'early_supported_difference' && axis_cure_support != 'supportive') {
    return('timing_difference_dominant')
  }
  
  if (axis_cure_support %in% c('supportive', 'equivocal') &&
      axis_followup_support != 'mostly_extrapolated' &&
      axis_remission_distortion != 'high') {
    return('cure_like_heterogeneity_plausible')
  }
  
  if (axis_followup_support == 'mostly_extrapolated') {
    return('tail_only_projection_dominant')
  }
  
  'mixed_or_indeterminate'
}

derive_future_recommendation_flag <- function(axis_cure_support, latency_family_signal, same_family_fit_gain_signal, axis_remission_distortion) {
  if (axis_cure_support == 'supportive' &&
      same_family_fit_gain_signal %in% c('consistent_gain', 'limited_gain') &&
      latency_family_signal %in% c('stable_multi_family_signal', 'limited_multi_family_signal') &&
      axis_remission_distortion != 'high') {
    return('consider_primary_cure_aware_modeling')
  }
  
  if (axis_cure_support == 'unsupportive') {
    return('no_primary_cure_modeling_recommendation')
  }
  
  'use_cure_modeling_as_sensitivity'
}

# 🔴 Load: stage-one backbone and registries ===============================
## 🟠 Read: backbone bundle ===============================
if (!dir.exists(stage_paths$stage1)) {
  stop('Stage 1 directory is required for Stage 10.', call. = FALSE)
}

stage1_bundle_file <- file_overrides$stage1$backbone_bundle
if (is.null(stage1_bundle_file) || !nzchar(stage1_bundle_file)) {
  candidate_rds <- list.files(stage_paths$stage1, pattern = 'stage1_.*backbone.*\\.rds$', full.names = TRUE, ignore.case = TRUE)
  if (length(candidate_rds) == 0L) {
    stop('Could not find `stage1_backbone_bundle.rds` automatically. Set `file_overrides$stage1$backbone_bundle`.', call. = FALSE)
  }
  stage1_bundle_file <- candidate_rds[[1L]]
}

backbone_bundle <- readRDS(stage1_bundle_file)

if (!is.list(backbone_bundle) || is.null(backbone_bundle$registries)) {
  stop('Stage 1 backbone bundle does not have the expected registry structure.', call. = FALSE)
}

horizon_registry <- backbone_bundle$registries$horizon_registry %>% tibble::as_tibble()
threshold_registry <- backbone_bundle$registries$threshold_registry %>% tibble::as_tibble()
metadata_registry_stage1 <- backbone_bundle$registries$metadata_registry %>% tibble::as_tibble()

dataset_order <- c('PNU', 'SNU', 'merged', 'merged__site_adjusted', 'merged__site_free', 'overall_project')

main_risk_scale <- first_non_missing_scalar(metadata_registry_stage1 %>% filter(metadata_name == 'main_risk_scale') %>% pull(metadata_value))
if (is.na(main_risk_scale) || !nzchar(main_risk_scale)) main_risk_scale <- 'transition_only_main'

supplementary_risk_scale <- first_non_missing_scalar(metadata_registry_stage1 %>% filter(metadata_name == 'supplementary_risk_scale') %>% pull(metadata_value))
if (is.na(supplementary_risk_scale) || !nzchar(supplementary_risk_scale)) supplementary_risk_scale <- 'transition_cif_competing'

# 🔴 Scan: stage directories for CSV outputs ===============================
## 🟠 Discover: available files ===============================
stage_scans <- list(
  Stage2 = scan_stage_directory(stage_paths$stage2, override_list = file_overrides$stage2),
  Stage3 = scan_stage_directory(stage_paths$stage3, override_list = file_overrides$stage3),
  Stage4 = scan_stage_directory(stage_paths$stage4, override_list = file_overrides$stage4),
  Stage5 = scan_stage_directory(stage_paths$stage5, override_list = file_overrides$stage5),
  Stage6 = scan_stage_directory(stage_paths$stage6, override_list = file_overrides$stage6),
  Stage7 = scan_stage_directory(stage_paths$stage7, override_list = file_overrides$stage7),
  Stage8 = scan_stage_directory(stage_paths$stage8, override_list = file_overrides$stage8),
  Stage9 = scan_stage_directory(stage_paths$stage9, override_list = file_overrides$stage9)
)

# 🔴 Assemble: stage-standardized tables ===============================
## 🟠 Build: stage-two follow-up support ===============================
stage2_followup_table <- purrr::map_dfr(
  collect_category_paths(stage_scans$Stage2, 'followup'),
  ~ standardize_followup_table(read_table_any(.x), .x, horizon_registry)
) %>%
  distinct()

## 🟠 Build: stage-four timing block ===============================
stage4_timing_table <- purrr::map_dfr(
  collect_category_paths(stage_scans$Stage4, 'timing'),
  ~ standardize_timing_table(read_table_any(.x), .x)
) %>%
  distinct()

## 🟠 Build: stage-six carry-forward flags ===============================
stage6_carryforward_table <- purrr::map_dfr(
  collect_category_paths(stage_scans$Stage6, 'carryforward'),
  ~ standardize_carryforward_table(read_table_any(.x), .x)
) %>%
  distinct(dataset_key, .keep_all = TRUE)

if (nrow(stage6_carryforward_table) == 0L) {
  stop('Stage 10 requires the Stage 6 carry-forward flag table.', call. = FALSE)
}

## 🟠 Build: stage-three, five, seven, eight, and nine metric bundles ===============================
stage3_bundle <- assemble_metric_stage('Stage3', stage_scans$Stage3, horizon_registry)
stage5_bundle <- assemble_metric_stage('Stage5', stage_scans$Stage5, horizon_registry)
stage7_bundle <- assemble_metric_stage('Stage7', stage_scans$Stage7, horizon_registry)
stage8_bundle <- assemble_metric_stage('Stage8', stage_scans$Stage8, horizon_registry)
stage9_bundle <- assemble_metric_stage('Stage9', stage_scans$Stage9, horizon_registry)

required_stage_rows <- c(
  Stage3 = nrow(stage3_bundle$classification_performance_table),
  Stage5 = nrow(stage5_bundle$classification_performance_table),
  Stage7 = nrow(stage7_bundle$classification_performance_table),
  Stage8 = nrow(stage8_bundle$classification_performance_table)
)

missing_required_stages <- names(required_stage_rows)[required_stage_rows == 0L]
if (length(missing_required_stages) > 0L) {
  stop(
    sprintf(
      'Stage 10 could not assemble the required classification/performance tables for: %s. Use `file_overrides` if auto-detection misses your filenames.',
      paste(missing_required_stages, collapse = ', ')
    ),
    call. = FALSE
  )
}

## 🟠 Build: stage-seven and stage-eight specialized tables ===============================
stage7_decomposition_table <- purrr::map_dfr(
  collect_category_paths(stage_scans$Stage7, 'decomposition'),
  ~ standardize_decomposition_table(read_table_any(.x), 'Stage7', .x, horizon_registry)
) %>%
  distinct() %>%
  ensure_columns(list(
    dataset_key = character(),
    dataset_parent = character(),
    source_stage = character(),
    model_family = character(),
    branch = character(),
    prior_branch = character(),
    horizon_year = integer(),
    risk_scale = character(),
    fit_key = character(),
    uncured_mean_support_flag = logical(),
    cure_fraction = numeric(),
    susceptible_fraction = numeric()
  ))

stage7_fit_contrast_table <- purrr::map_dfr(
  collect_category_paths(stage_scans$Stage7, 'fit_contrast'),
  ~ standardize_fit_contrast_table(read_table_any(.x), .x)
) %>%
  distinct() %>%
  ensure_columns(list(
    dataset_key = character(),
    dataset_parent = character(),
    delta_AIC_cure_minus_noncure = numeric(),
    same_family_fit_gain_signal = character()
  ))

stage8_decomposition_table <- purrr::map_dfr(
  collect_category_paths(stage_scans$Stage8, 'decomposition'),
  ~ standardize_decomposition_table(read_table_any(.x), 'Stage8', .x, horizon_registry)
) %>%
  distinct() %>%
  ensure_columns(list(
    dataset_key = character(),
    dataset_parent = character(),
    source_stage = character(),
    model_family = character(),
    branch = character(),
    prior_branch = character(),
    horizon_year = integer(),
    risk_scale = character(),
    fit_key = character(),
    uncured_mean_support_flag = logical(),
    cure_fraction = numeric(),
    susceptible_fraction = numeric()
  ))

stage8_anchor_delta_table <- purrr::map_dfr(
  collect_category_paths(stage_scans$Stage8, 'anchor_delta'),
  ~ standardize_anchor_delta_table(read_table_any(.x), .x, horizon_registry)
) %>%
  distinct() %>%
  ensure_columns(list(
    dataset_key = character(),
    dataset_parent = character(),
    branch = character(),
    prior_branch = character(),
    horizon_year = integer(),
    threshold = numeric(),
    fit_key = character(),
    delta_risk_anchor_minus_neutral = numeric(),
    delta_false_positive_burden_anchor_minus_neutral = numeric(),
    delta_FP100_anchor_minus_neutral = numeric(),
    delta_NB_anchor_minus_neutral = numeric()
  ))

stage8_incidence_update_table <- purrr::map_dfr(
  collect_category_paths(stage_scans$Stage8, 'incidence_update'),
  ~ standardize_incidence_update_table(read_table_any(.x), .x)
) %>%
  distinct() %>%
  ensure_columns(list(
    dataset_key = character(),
    dataset_parent = character(),
    branch = character(),
    prior_branch = character(),
    fit_key = character(),
    age_sex_anchor_cell = character(),
    posterior_minus_prior_risk = numeric()
  ))

stage8_remission_delta_table <- purrr::map_dfr(
  collect_category_paths(stage_scans$Stage8, 'remission_delta'),
  ~ standardize_remission_delta_table(read_table_any(.x), .x, horizon_registry)
) %>%
  distinct() %>%
  ensure_columns(list(
    dataset_key = character(),
    dataset_parent = character(),
    branch = character(),
    prior_branch = character(),
    horizon_year = integer(),
    threshold = numeric(),
    fit_key = character(),
    delta_risk_8B_minus_8A = numeric(),
    delta_false_positive_burden_8B_minus_8A = numeric(),
    delta_FP100_8B_minus_8A = numeric(),
    delta_NB_8B_minus_8A = numeric()
  ))

stage8_diagnostics_table <- purrr::map_dfr(
  collect_category_paths(stage_scans$Stage8, 'diagnostics'),
  ~ standardize_diagnostics_table(read_table_any(.x), .x)
) %>%
  distinct() %>%
  ensure_columns(list(
    fit_key = character(),
    admissibility_flag = logical(),
    convergence_status = character(),
    posterior_predictive_check_flag = logical(),
    prior_sensitivity_ok_flag = logical()
  ))

## 🟠 Apply: stage-eight admissibility enrichment ===============================
stage8_bundle$prediction_table <- apply_stage8_diagnostics(stage8_bundle$prediction_table, stage8_diagnostics_table)
stage8_bundle$classification_performance_table <- apply_stage8_diagnostics(stage8_bundle$classification_performance_table, stage8_diagnostics_table)
stage8_decomposition_table <- apply_stage8_diagnostics(stage8_decomposition_table, stage8_diagnostics_table)
stage8_anchor_delta_table <- apply_stage8_diagnostics(stage8_anchor_delta_table, stage8_diagnostics_table)
stage8_incidence_update_table <- apply_stage8_diagnostics(stage8_incidence_update_table, stage8_diagnostics_table)
stage8_remission_delta_table <- apply_stage8_diagnostics(stage8_remission_delta_table, stage8_diagnostics_table)

stage8_bundle$prediction_table <- stage8_bundle$prediction_table %>% filter(is.na(admissibility_flag) | admissibility_flag)
stage8_bundle$classification_performance_table <- stage8_bundle$classification_performance_table %>% filter(is.na(admissibility_flag) | admissibility_flag)
stage8_decomposition_table <- stage8_decomposition_table %>% filter(is.na(admissibility_flag) | admissibility_flag)
stage8_anchor_delta_table <- stage8_anchor_delta_table %>% filter(is.na(admissibility_flag) | admissibility_flag)
stage8_incidence_update_table <- stage8_incidence_update_table %>% filter(is.na(admissibility_flag) | admissibility_flag)
stage8_remission_delta_table <- stage8_remission_delta_table %>% filter(is.na(admissibility_flag) | admissibility_flag)

## 🟠 Build: stage-nine remission-versus-main deltas ===============================
stage9_delta_from_main_table <- purrr::map_dfr(
  unique(c(
    collect_category_paths(stage_scans$Stage9, 'classification'),
    collect_category_paths(stage_scans$Stage9, 'decision'),
    collect_category_paths(stage_scans$Stage9, 'performance')
  )),
  ~ standardize_stage9_delta_from_main(read_table_any(.x), .x, horizon_registry)
) %>%
  filter(!(is.na(delta_risk_9_minus_main) & is.na(delta_false_positive_burden_9_minus_main) & is.na(delta_NB_9_minus_main))) %>%
  distinct() %>%
  ensure_columns(list(
    dataset_key = character(),
    dataset_parent = character(),
    branch = character(),
    prior_branch = character(),
    horizon_year = integer(),
    threshold = numeric(),
    delta_risk_9_minus_main = numeric(),
    delta_false_positive_burden_9_minus_main = numeric(),
    delta_NB_9_minus_main = numeric()
  ))

# 🔴 Construct: unified long-format source tables ===============================
## 🟠 Combine: prediction source-of-truth table ===============================
long_prediction_table <- bind_rows(
  stage3_bundle$prediction_table,
  stage5_bundle$prediction_table,
  stage7_bundle$prediction_table,
  stage8_bundle$prediction_table,
  stage9_bundle$prediction_table
) %>%
  left_join(
    stage6_carryforward_table %>% select(-dataset_parent, -carry_forward_stage8),
    by = 'dataset_key'
  ) %>%
  mutate(
    risk_scale = normalize_risk_scale(risk_scale),
    model_display_group = factor(
      model_class,
      levels = c('KM_benchmark', 'no_cure', 'frequentist_cure', 'bayesian_cure', 'remission_sensitive_frequentist')
    ),
    prior_branch = normalize_prior_branch(prior_branch),
    site_prior_family = normalize_site_prior_family(site_prior_family)
  ) %>%
  arrange(dataset_sort_num(dataset_key, dataset_order), horizon_year, model_display_group, model_family, branch, prior_branch)

if (nrow(long_prediction_table) == 0L) {
  stop('Stage 10 requires at least one prediction source table. If your filenames are unusual, set `file_overrides`.', call. = FALSE)
}

## 🟠 Combine: classification and performance source-of-truth table ===============================
long_classification_performance_table <- bind_rows(
  stage3_bundle$classification_performance_table,
  stage5_bundle$classification_performance_table,
  stage7_bundle$classification_performance_table,
  stage8_bundle$classification_performance_table,
  stage9_bundle$classification_performance_table
) %>%
  left_join(
    stage6_carryforward_table %>% select(-dataset_parent, -carry_forward_stage8),
    by = 'dataset_key'
  ) %>%
  mutate(
    risk_scale = normalize_risk_scale(risk_scale),
    model_display_group = factor(
      model_class,
      levels = c('KM_benchmark', 'no_cure', 'frequentist_cure', 'bayesian_cure', 'remission_sensitive_frequentist')
    ),
    prior_branch = normalize_prior_branch(prior_branch),
    site_prior_family = normalize_site_prior_family(site_prior_family)
  ) %>%
  arrange(dataset_sort_num(dataset_key, dataset_order), risk_scale, horizon_year, threshold, model_display_group, model_family, branch, prior_branch)

## 🟠 Derive: reporting metadata table ===============================
reporting_metadata_table <- long_classification_performance_table %>%
  transmute(
    dataset_key,
    dataset_parent,
    analysis_structure,
    reporting_priority,
    risk_scale,
    horizon_year,
    threshold,
    support_tier,
    horizon_evidence_class,
    claim_restriction_flag,
    prior_tail_sensitive = ifelse(is.na(prior_tail_sensitive), FALSE, prior_tail_sensitive),
    prior_branch,
    site_prior_family,
    branch,
    model_class,
    model_family,
    retained_fit_id,
    risk_set_size,
    instability_marker
  ) %>%
  distinct() %>%
  arrange(dataset_sort_num(dataset_key, dataset_order), risk_scale, horizon_year, threshold, model_class, model_family, branch, prior_branch)

# 🔴 Construct: five coordinated reporting objects ===============================
## 🟠 Create: main common-scale comparison table ===============================
main_common_scale_comparison_table <- long_classification_performance_table %>%
  filter(risk_scale == main_risk_scale) %>%
  filter(model_class %in% c('KM_benchmark', 'no_cure', 'frequentist_cure', 'bayesian_cure')) %>%
  mutate(
    remission_reporting_status = 'primary_common_scale',
    merged_site_interpretation_note = ifelse(
      grepl('merged', as.character(dataset_key), ignore.case = TRUE),
      'Merged site terms are structural context proxies, not causal treatment effects.',
      NA_character_
    )
  ) %>%
  arrange(dataset_sort_num(dataset_key, dataset_order), horizon_year, threshold, model_display_group, model_family, branch, prior_branch)

## 🟠 Create: remission-aware comparison table ===============================
remission_aware_comparison_table <- long_classification_performance_table %>%
  filter(risk_scale == supplementary_risk_scale) %>%
  filter(model_class %in% c('bayesian_cure', 'remission_sensitive_frequentist')) %>%
  mutate(
    remission_reporting_status = 'supplementary_competing_risk_scale',
    cross_scale_ranking_forbidden_note = 'Cross-scale Stage8A-versus-Stage8B differences are remission-aware change summaries, not fair within-scale rankings.'
  ) %>%
  arrange(dataset_sort_num(dataset_key, dataset_order), horizon_year, threshold, model_display_group, model_family, branch, prior_branch)

## 🟠 Create: cure-model-only supporting decomposition block ===============================
cure_model_only_supporting_decomposition_table <- bind_rows(stage7_decomposition_table, stage8_decomposition_table) %>%
  left_join(
    stage6_carryforward_table %>% select(-dataset_parent, -carry_forward_stage8),
    by = 'dataset_key'
  ) %>%
  mutate(
    decomposition_separation_rule = 'Keep separate from the main common-scale comparison table.',
    merged_site_proxy_warning = ifelse(
      grepl('merged', dataset_key, ignore.case = TRUE),
      'Merged site-placement terms are structural decomposition labels, not causal treatment effects.',
      NA_character_
    )
  ) %>%
  arrange(dataset_sort_num(dataset_key, dataset_order), source_stage, model_family, branch, prior_branch, horizon_year)

## 🟠 Create: remission-aware delta block ===============================
stage8_remission_delta_table <- stage8_remission_delta_table %>%
  mutate(
    delta_source = 'Stage8_exported_8A_vs_8B',
    delta_risk_9_minus_main = NA_real_,
    delta_false_positive_burden_9_minus_main = NA_real_,
    delta_NB_9_minus_main = NA_real_
  )

stage9_delta_from_main_table <- stage9_delta_from_main_table %>%
  mutate(
    delta_source = 'Stage9_exported_vs_main',
    delta_risk_8B_minus_8A = NA_real_,
    delta_cure_fraction_8B_minus_8A = NA_real_,
    delta_false_positive_burden_8B_minus_8A = delta_false_positive_burden_9_minus_main,
    delta_FP100_8B_minus_8A = NA_real_,
    delta_NB_8B_minus_8A = delta_NB_9_minus_main,
    delta_PPV_8B_minus_8A = NA_real_,
    delta_TPR_8B_minus_8A = NA_real_
  )

remission_aware_delta_table <- bind_rows(stage8_remission_delta_table, stage9_delta_from_main_table) %>%
  mutate(
    remission_delta_interpretation_note = 'Use as remission-aware change summaries only; do not treat as a fair within-scale model-ranking contest.'
  ) %>%
  arrange(dataset_sort_num(dataset_key, dataset_order), delta_source, horizon_year, threshold, branch, prior_branch)

## 🟠 Create: external-anchor value block ===============================
external_anchor_value_table <- stage8_anchor_delta_table %>%
  mutate(
    anchor_interpretation_note = 'Quantifies what changes when the external incidence anchor is removed; anchored fits are not automatically correct.',
    incidence_anchor_update_available_flag = fit_key %in% stage8_incidence_update_table$fit_key
  ) %>%
  arrange(dataset_sort_num(dataset_key, dataset_order), branch, prior_branch, horizon_year, threshold)

# 🔴 Derive: heterogeneity triangulation and interpretation axes ===============================
## 🟠 Prepare: dataset scopes ===============================
heterogeneity_scope_tbl <- tibble::tibble(
  dataset_key = unique(c(
    'overall_project',
    as.character(long_classification_performance_table$dataset_key),
    as.character(stage6_carryforward_table$dataset_key),
    as.character(stage2_followup_table$dataset_key),
    as.character(cure_model_only_supporting_decomposition_table$dataset_key)
  ))
) %>%
  filter(!is.na(dataset_key) & nzchar(dataset_key)) %>%
  mutate(
    dataset_parent = dataset_parent_from_key(dataset_key),
    analysis_structure = derive_analysis_structure(dataset_key),
    reporting_priority = derive_reporting_priority(dataset_key)
  ) %>%
  distinct() %>%
  arrange(dataset_sort_num(dataset_key, dataset_order))

## 🟠 Compute: dataset-level triangulation rows ===============================
heterogeneity_triangulation_table <- heterogeneity_scope_tbl %>%
  rowwise() %>%
  mutate(
    timing_signal = {
      this_dataset <- dataset_key
      if (this_dataset %in% c('merged', 'merged__site_free', 'merged__site_adjusted', 'overall_project')) {
        derive_timing_signal(stage4_timing_table)
      } else {
        'single_site_no_within_dataset_site_contrast'
      }
    },
    followup_signal = derive_followup_signal(dataset_key, stage2_followup_table, horizon_registry),
    screening_signal = {
      picked <- stage6_carryforward_table %>%
        filter(dataset_key == !!dataset_key | dataset_parent == !!dataset_parent) %>%
        arrange(dataset_sort_num(dataset_key, dataset_order))
      if (nrow(picked) == 0L) 'screening_unavailable' else first_non_missing_scalar(picked$cure_model_eligibility_flag)
    },
    latency_family_signal = derive_latency_family_signal(
      main_common_scale_comparison_table %>% filter(as.character(dataset_key) == !!dataset_key | dataset_parent == !!dataset_parent)
    ),
    remission_signal = derive_remission_distortion_signal(
      remission_aware_delta_table %>% filter(as.character(dataset_key) == !!dataset_key | dataset_parent == !!dataset_parent)
    ),
    anchor_dependence_signal = derive_anchor_dependence_signal(
      external_anchor_value_table %>% filter(as.character(dataset_key) == !!dataset_key | dataset_parent == !!dataset_parent)
    ),
    same_family_fit_gain_signal = derive_same_family_fit_gain_signal(
      stage7_fit_contrast_table %>% filter(as.character(dataset_key) == !!dataset_key | dataset_parent == !!dataset_parent)
    ),
    uncured_only_support_signal = derive_uncured_only_support_signal(
      cure_model_only_supporting_decomposition_table %>% filter(as.character(dataset_key) == !!dataset_key | dataset_parent == !!dataset_parent)
    ),
    axis_timing_difference = dplyr::case_when(
      timing_signal %in% c('early_supported_difference', 'late_only_difference', 'weak_or_absent_early_difference') ~ timing_signal,
      TRUE ~ 'weak_or_absent_early_difference'
    ),
    axis_followup_support = dplyr::case_when(
      followup_signal %in% c('directly_observed_data_supported', 'partly_model_dependent', 'mostly_extrapolated') ~ followup_signal,
      TRUE ~ 'mostly_extrapolated'
    ),
    axis_cure_support = dplyr::case_when(
      screening_signal %in% c('supportive', 'equivocal', 'unsupportive') ~ screening_signal,
      TRUE ~ 'equivocal'
    ),
    axis_remission_distortion = dplyr::case_when(
      remission_signal %in% c('low', 'moderate', 'high') ~ remission_signal,
      TRUE ~ 'moderate'
    ),
    interpretation_bucket = derive_interpretation_bucket(axis_timing_difference, axis_followup_support, axis_cure_support, axis_remission_distortion),
    triangulated_heterogeneity_bucket = interpretation_bucket,
    future_cure_modeling_recommendation_flag = derive_future_recommendation_flag(axis_cure_support, latency_family_signal, same_family_fit_gain_signal, axis_remission_distortion)
  ) %>%
  ungroup() %>%
  arrange(dataset_sort_num(dataset_key, dataset_order))

# 🔴 Assemble: mandatory figure-ready source tables ===============================
## 🟠 Create: horizon support panel source ===============================
risk_scale_levels <- unique(c(main_risk_scale, supplementary_risk_scale, long_classification_performance_table$risk_scale))
risk_scale_levels <- risk_scale_levels[!is.na(risk_scale_levels) & nzchar(risk_scale_levels)]

if (nrow(stage2_followup_table) > 0L) {
  horizon_support_seed <- stage2_followup_table %>%
    select(dataset_key, dataset_parent, analysis_structure, reporting_priority, horizon_year, risk_set_size, support_tier, horizon_evidence_class, claim_restriction_flag, late_horizon_instability_flag) %>%
    distinct()
} else {
  horizon_support_seed <- horizon_registry %>%
    transmute(
      dataset_key = normalize_dataset_key(dataset),
      dataset_parent = normalize_dataset_key(dataset),
      analysis_structure = derive_analysis_structure(dataset_key),
      reporting_priority = derive_reporting_priority(dataset_key),
      horizon_year = as.integer(horizon_year),
      risk_set_size = NA_real_,
      support_tier = support_tier,
      horizon_evidence_class = horizon_evidence_class,
      claim_restriction_flag = claim_restriction_flag,
      late_horizon_instability_flag = NA
    )
}

horizon_support_panel_source <- horizon_support_seed %>%
  tidyr::expand_grid(risk_scale = risk_scale_levels) %>%
  left_join(
    reporting_metadata_table %>%
      group_by(dataset_key, horizon_year, risk_scale) %>%
      summarise(prior_tail_sensitive = any(prior_tail_sensitive %in% TRUE, na.rm = TRUE), .groups = 'drop'),
    by = c('dataset_key', 'horizon_year', 'risk_scale')
  ) %>%
  mutate(prior_tail_sensitive = ifelse(is.na(prior_tail_sensitive), FALSE, prior_tail_sensitive)) %>%
  arrange(dataset_sort_num(dataset_key, dataset_order), risk_scale, horizon_year)

## 🟠 Create: 8A-vs-8B delta panel source ===============================
stage8A_vs_stage8B_delta_panel_source <- remission_aware_delta_table %>%
  filter(delta_source == 'Stage8_exported_8A_vs_8B') %>%
  arrange(dataset_sort_num(dataset_key, dataset_order), horizon_year, threshold, branch, prior_branch)

## 🟠 Create: anchor-versus-neutral delta panel source ===============================
anchor_vs_neutral_delta_panel_source <- external_anchor_value_table %>%
  arrange(dataset_sort_num(dataset_key, dataset_order), horizon_year, threshold, branch, prior_branch)

## 🟠 Create: incidence anchor-update panel source ===============================
incidence_anchor_update_panel_source <- stage8_incidence_update_table %>%
  arrange(dataset_sort_num(dataset_key, dataset_order), branch, prior_branch, age_sex_anchor_cell)

## 🟠 Create: heterogeneity triangulation panel source ===============================
heterogeneity_triangulation_panel_source <- heterogeneity_triangulation_table %>%
  arrange(dataset_sort_num(dataset_key, dataset_order))

# 🔴 Export: stage-ten coordinated outputs ===============================
## 🟠 Save: CSV source tables ===============================
readr::write_csv(long_prediction_table, file.path(export_path, 'stage10_long_prediction_table.csv'))
readr::write_csv(long_classification_performance_table, file.path(export_path, 'stage10_long_classification_performance_table.csv'))
readr::write_csv(reporting_metadata_table, file.path(export_path, 'stage10_reporting_metadata_table.csv'))
readr::write_csv(main_common_scale_comparison_table, file.path(export_path, 'stage10_main_common_scale_comparison_table.csv'))
readr::write_csv(remission_aware_comparison_table, file.path(export_path, 'stage10_remission_aware_comparison_table.csv'))
readr::write_csv(cure_model_only_supporting_decomposition_table, file.path(export_path, 'stage10_cure_model_only_supporting_decomposition_table.csv'))
readr::write_csv(remission_aware_delta_table, file.path(export_path, 'stage10_remission_aware_delta_table.csv'))
readr::write_csv(external_anchor_value_table, file.path(export_path, 'stage10_external_anchor_value_table.csv'))
readr::write_csv(heterogeneity_triangulation_table, file.path(export_path, 'stage10_heterogeneity_triangulation_table.csv'))
readr::write_csv(horizon_support_panel_source, file.path(export_path, 'stage10_horizon_support_panel_source.csv'))
readr::write_csv(stage8A_vs_stage8B_delta_panel_source, file.path(export_path, 'stage10_8A_vs_8B_delta_panel_source.csv'))
readr::write_csv(anchor_vs_neutral_delta_panel_source, file.path(export_path, 'stage10_anchor_vs_neutral_delta_panel_source.csv'))
readr::write_csv(incidence_anchor_update_panel_source, file.path(export_path, 'stage10_incidence_anchor_update_panel_source.csv'))
readr::write_csv(heterogeneity_triangulation_panel_source, file.path(export_path, 'stage10_heterogeneity_triangulation_panel_source.csv'))

## 🟠 Save: reusable bundle ===============================
stage10_bundle <- list(
  stage = 'Stage 10',
  created_at = as.character(Sys.time()),
  session_info = utils::sessionInfo(),
  config = list(
    data_path = data_path,
    export_path = export_path,
    stage_paths = stage_paths,
    file_overrides = file_overrides,
    main_risk_scale = main_risk_scale,
    supplementary_risk_scale = supplementary_risk_scale
  ),
  stage_file_scan = stage_scans,
  stage_inputs = list(
    horizon_registry = horizon_registry,
    threshold_registry = threshold_registry,
    stage2_followup_table = stage2_followup_table,
    stage4_timing_table = stage4_timing_table,
    stage6_carryforward_table = stage6_carryforward_table,
    stage7_fit_contrast_table = stage7_fit_contrast_table,
    stage8_diagnostics_table = stage8_diagnostics_table,
    stage8_incidence_update_table = stage8_incidence_update_table
  ),
  outputs = list(
    long_prediction_table = long_prediction_table,
    long_classification_performance_table = long_classification_performance_table,
    reporting_metadata_table = reporting_metadata_table,
    main_common_scale_comparison_table = main_common_scale_comparison_table,
    remission_aware_comparison_table = remission_aware_comparison_table,
    cure_model_only_supporting_decomposition_table = cure_model_only_supporting_decomposition_table,
    remission_aware_delta_table = remission_aware_delta_table,
    external_anchor_value_table = external_anchor_value_table,
    heterogeneity_triangulation_table = heterogeneity_triangulation_table,
    horizon_support_panel_source = horizon_support_panel_source,
    stage8A_vs_stage8B_delta_panel_source = stage8A_vs_stage8B_delta_panel_source,
    anchor_vs_neutral_delta_panel_source = anchor_vs_neutral_delta_panel_source,
    incidence_anchor_update_panel_source = incidence_anchor_update_panel_source,
    heterogeneity_triangulation_panel_source = heterogeneity_triangulation_panel_source
  )
)

saveRDS(stage10_bundle, file.path(export_path, 'stage10_unified_comparison_bundle.rds'))

## 🟠 Print: completion summary ===============================
cat('\nStage 10 completed. Exported files:\n')
cat(' - stage10_long_prediction_table.csv\n')
cat(' - stage10_long_classification_performance_table.csv\n')
cat(' - stage10_reporting_metadata_table.csv\n')
cat(' - stage10_main_common_scale_comparison_table.csv\n')
cat(' - stage10_remission_aware_comparison_table.csv\n')
cat(' - stage10_cure_model_only_supporting_decomposition_table.csv\n')
cat(' - stage10_remission_aware_delta_table.csv\n')
cat(' - stage10_external_anchor_value_table.csv\n')
cat(' - stage10_heterogeneity_triangulation_table.csv\n')
cat(' - stage10_horizon_support_panel_source.csv\n')
cat(' - stage10_8A_vs_8B_delta_panel_source.csv\n')
cat(' - stage10_anchor_vs_neutral_delta_panel_source.csv\n')
cat(' - stage10_incidence_anchor_update_panel_source.csv\n')
cat(' - stage10_heterogeneity_triangulation_panel_source.csv\n')
cat(' - stage10_unified_comparison_bundle.rds\n')