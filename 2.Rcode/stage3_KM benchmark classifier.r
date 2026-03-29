# 🔴 Configure: stage-three paths and reusable inputs ===============================
data_path <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/stage1_Backbone lock'
export_path <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/stage3_KM benchmark classifier'

stage1_bundle_file <- file.path(data_path, "stage1_backbone_bundle.rds")
stage1_analysis_datasets_file <- file.path(data_path, "stage1_analysis_datasets.rds")
stage1_manifest_file <- file.path(data_path, "stage1_export_manifest.csv")
stage1_dataset_registry_file <- file.path(data_path, "stage1_dataset_registry.csv")
stage1_horizon_registry_file <- file.path(data_path, "stage1_horizon_registry.csv")
stage1_threshold_registry_file <- file.path(data_path, "stage1_threshold_registry.csv")
stage1_metadata_registry_file <- file.path(data_path, "stage1_metadata_registry.csv")
stage1_scaling_registry_file <- file.path(data_path, "stage1_scaling_registry.csv")
stage1_modeling_registry_file <- file.path(data_path, "stage1_modeling_registry.csv")

run_subgroup_km <- FALSE
subgroup_variables <- c("sex_label")
subgroup_min_n <- 30L

# 🔴 Initialize: package state and global options ===============================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(survival)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

dataset_order <- c("PNU", "SNU", "merged")

# 🔴 Define: backbone-safe helper functions ===============================
## 🟠 Define: file and scalar utilities ===============================
assert_existing_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

read_stage1_csv <- function(path) {
  readr::read_csv(
    file = path,
    show_col_types = FALSE,
    progress = FALSE
  )
}

normalize_existing_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

null_to_na_character <- function(x) {
  if (is.null(x)) {
    return(NA_character_)
  }
  as.character(x)
}

safe_divide <- function(num, den) {
  ifelse(is.na(den) | den <= 0, NA_real_, num / den)
}

collapse_unique_chr <- function(x) {
  paste(sort(unique(as.character(x))), collapse = "|")
}

apply_threshold_rule <- function(risk_value, threshold_value) {
  as.logical(risk_value >= threshold_value)
}

make_stage3_dataset_key <- function(dataset_name, site_branch, benchmark_scope, subgroup_variable = NULL) {
  if (identical(benchmark_scope, "sensitivity") && !is.null(subgroup_variable) && !is.na(subgroup_variable)) {
    return(paste(dataset_name, paste0("subgroup_", subgroup_variable), sep = "__"))
  }
  if (dataset_name == "merged" && site_branch == "site_adjusted") {
    return("merged__site_adjusted")
  }
  if (dataset_name == "merged") {
    return("merged__site_free")
  }
  dataset_name
}

make_stage3_analysis_structure <- function(dataset_name, site_branch, benchmark_scope) {
  if (identical(benchmark_scope, "sensitivity")) {
    return("subgroup_sensitivity")
  }
  if (dataset_name == "merged" && site_branch == "site_adjusted") {
    return("site_adjusted_structural")
  }
  if (dataset_name == "merged") {
    return("site_free_backbone")
  }
  "single_site"
}

# 🔴 Define: Stage-1 compatibility checks ===============================
## 🟠 Define: inherited backbone validator ===============================
validate_stage1_backbone <- function(
    bundle,
    analysis_datasets,
    dataset_registry,
    horizon_registry,
    threshold_registry,
    metadata_registry,
    scaling_registry,
    modeling_registry
) {
  required_dataset_names <- c("PNU", "SNU", "merged")
  required_dataset_cols <- c(
    "unique_person_id", "site", "sex_num", "sex_label",
    "age_exact_entry", "age_s", "days_followup", "time_year",
    "status_num", "event_main", "right_censor_flag", "remission_flag", "censor_main"
  )
  required_horizon_cols <- c(
    "dataset", "dataset_key", "horizon_id", "horizon_year", "horizon_days",
    "support_tier", "horizon_evidence_class", "claim_restriction_flag"
  )
  required_threshold_cols <- c("threshold_id", "threshold", "threshold_label", "positive_rule")

  if (!is.list(bundle)) {
    stop("Stage 1 backbone bundle must be a list.", call. = FALSE)
  }

  if (!identical(bundle$stage, "Stage 1")) {
    stop("Stage 1 backbone bundle must declare `stage = \"Stage 1\"`.", call. = FALSE)
  }

  if (!is.list(bundle$datasets) || !all(required_dataset_names %in% names(bundle$datasets))) {
    stop("Stage 1 backbone bundle must contain datasets for PNU, SNU, and merged.", call. = FALSE)
  }

  if (!is.list(analysis_datasets) || !all(required_dataset_names %in% names(analysis_datasets))) {
    stop("Stage 1 analysis datasets RDS must contain PNU, SNU, and merged.", call. = FALSE)
  }

  if (nrow(dataset_registry) == 0 || nrow(horizon_registry) == 0 || nrow(threshold_registry) == 0 ||
      nrow(metadata_registry) == 0 || nrow(scaling_registry) == 0 || nrow(modeling_registry) == 0) {
    stop("One or more required Stage 1 registry/metadata/scaling/modeling files are empty.", call. = FALSE)
  }

  missing_horizon_cols <- setdiff(required_horizon_cols, names(horizon_registry))
  if (length(missing_horizon_cols) > 0) {
    stop(
      sprintf("Stage 1 horizon registry is missing required columns: %s", paste(missing_horizon_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  missing_threshold_cols <- setdiff(required_threshold_cols, names(threshold_registry))
  if (length(missing_threshold_cols) > 0) {
    stop(
      sprintf("Stage 1 threshold registry is missing required columns: %s", paste(missing_threshold_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  observed_horizons <- sort(unique(as.integer(horizon_registry$horizon_year)))
  if (!identical(observed_horizons, 1:10)) {
    stop("Stage 1 horizon registry must contain the common 1:10 year horizon grid.", call. = FALSE)
  }

  bundle_thresholds <- sort(unique(as.numeric(bundle$config$risk_thresholds)))
  registry_thresholds <- sort(unique(as.numeric(threshold_registry$threshold)))
  if (!identical(bundle_thresholds, registry_thresholds)) {
    stop("Stage 1 threshold registry does not match Stage 1 bundle risk_thresholds.", call. = FALSE)
  }

  if (!identical(as.character(bundle$config$main_risk_scale), "transition_only_main")) {
    stop("Stage 1 bundle main risk scale must be `transition_only_main`.", call. = FALSE)
  }

  if (!all(required_dataset_names %in% dataset_registry$dataset)) {
    stop("Stage 1 dataset registry must include PNU, SNU, and merged.", call. = FALSE)
  }

  if (!all(required_dataset_names %in% scaling_registry$dataset)) {
    stop("Stage 1 scaling registry must include PNU, SNU, and merged.", call. = FALSE)
  }

  if (!all(required_dataset_names %in% horizon_registry$dataset)) {
    stop("Stage 1 horizon registry must include PNU, SNU, and merged.", call. = FALSE)
  }

  for (dataset_name in required_dataset_names) {
    df <- analysis_datasets[[dataset_name]]

    missing_cols <- setdiff(required_dataset_cols, names(df))
    if (length(missing_cols) > 0) {
      stop(
        sprintf("[%s] Stage 1 analysis dataset is missing required columns: %s",
                dataset_name, paste(missing_cols, collapse = ", ")),
        call. = FALSE
      )
    }

    if (nrow(df) != dplyr::n_distinct(df$unique_person_id)) {
      stop(sprintf("[%s] `unique_person_id` must be unique.", dataset_name), call. = FALSE)
    }

    if (anyNA(df[, required_dataset_cols])) {
      stop(sprintf("[%s] Missing values detected in required Stage 1 prepared columns.", dataset_name), call. = FALSE)
    }

    if (any(df$days_followup < 0, na.rm = TRUE)) {
      stop(sprintf("[%s] Negative `days_followup` detected.", dataset_name), call. = FALSE)
    }

    if (any(!df$sex_num %in% c(0L, 1L), na.rm = TRUE)) {
      stop(sprintf("[%s] `sex_num` must remain coded as 0/1.", dataset_name), call. = FALSE)
    }

    if (any(!df$status_num %in% c(0L, 1L, 2L), na.rm = TRUE)) {
      stop(sprintf("[%s] `status_num` must remain coded as 0/1/2.", dataset_name), call. = FALSE)
    }

    if (dataset_name %in% c("PNU", "SNU") && dplyr::n_distinct(df$site) != 1L) {
      stop(sprintf("[%s] Single-cohort Stage 1 dataset must contain exactly one site value.", dataset_name), call. = FALSE)
    }

    if (dataset_name == "merged" && dplyr::n_distinct(df$site) < 2L) {
      stop("[merged] Merged Stage 1 dataset must contain at least two site values.", call. = FALSE)
    }
  }

  invisible(TRUE)
}

# 🔴 Define: Kaplan-Meier benchmark constructors ===============================
## 🟠 Define: survfit extraction helpers ===============================
extract_survfit_lookup <- function(fit, horizons_tbl, strata_variable = NULL) {
  fit_summary <- summary(fit, times = horizons_tbl$horizon_year, extend = TRUE)

  if (is.null(fit_summary$strata)) {
    out <- tibble(
      assignment_value = "overall",
      horizon_id = horizons_tbl$horizon_id,
      horizon_year = as.integer(horizons_tbl$horizon_year),
      horizon_days = as.numeric(horizons_tbl$horizon_days),
      km_survival = as.numeric(fit_summary$surv)
    )
  } else {
    if (is.null(strata_variable) || is.na(strata_variable) || identical(strata_variable, "")) {
      stop("`strata_variable` must be provided for stratified KM fits.", call. = FALSE)
    }

    strata_prefix <- paste0("^", strata_variable, "=")
    n_strata <- length(fit$strata)

    out <- tibble(
      strata_label = as.character(fit_summary$strata),
      horizon_id = rep(horizons_tbl$horizon_id, times = n_strata),
      horizon_year = rep(as.integer(horizons_tbl$horizon_year), times = n_strata),
      horizon_days = rep(as.numeric(horizons_tbl$horizon_days), times = n_strata),
      km_survival = as.numeric(fit_summary$surv)
    ) %>%
      mutate(
        assignment_value = sub(strata_prefix, "", strata_label)
      ) %>%
      select(assignment_value, horizon_id, horizon_year, horizon_days, km_survival)
  }

  out %>%
    mutate(
      km_survival = pmin(pmax(km_survival, 0), 1),
      km_risk = 1 - km_survival
    )
}

build_benchmark_registry_row <- function(df, risk_lookup, meta) {
  assignment_variable_export <- if (identical(meta$assignment_variable, "__overall__")) {
    "overall_constant"
  } else {
    null_to_na_character(meta$assignment_variable)
  }

  tibble(
    dataset = meta$dataset,
    dataset_key = meta$dataset_key,
    base_dataset_key = meta$base_dataset_key,
    analysis_structure = meta$analysis_structure,
    risk_scale = meta$risk_scale,
    model_class = "KM_benchmark",
    model_family = "KM",
    benchmark_id = meta$benchmark_id,
    benchmark_label = meta$benchmark_label,
    benchmark_scope = meta$benchmark_scope,
    benchmark_structure = meta$benchmark_structure,
    site_branch = meta$site_branch,
    subgroup_kind = meta$subgroup_kind,
    subgroup_variable = null_to_na_character(meta$subgroup_variable),
    assignment_variable = assignment_variable_export,
    risk_assignment_rule = meta$risk_assignment_rule,
    primary_benchmark_flag = identical(meta$benchmark_scope, "primary"),
    structural_site_adjusted_flag = identical(meta$benchmark_scope, "primary_structural"),
    sensitivity_benchmark_flag = identical(meta$benchmark_scope, "sensitivity"),
    n_subjects = nrow(df),
    n_groups = dplyr::n_distinct(risk_lookup$assignment_value),
    group_values = collapse_unique_chr(risk_lookup$assignment_value),
    threshold_source = "Stage 1 threshold registry",
    horizon_source = "Stage 1 horizon registry",
    event_definition = "status_num == 1",
    censoring_definition = "status_num %in% c(0, 2)",
    site_effect_interpretation = meta$site_effect_interpretation,
    prediction_type = "group-level KM risk assigned to subjects according to benchmark grouping"
  )
}

build_subject_predictions <- function(df, horizons_tbl, risk_lookup, meta) {
  assignment_variable <- meta$assignment_variable

  subject_base <- df %>%
    transmute(
      unique_person_id,
      site,
      sex_num,
      sex_label = as.character(sex_label),
      days_followup,
      time_year,
      status_num,
      event_main,
      right_censor_flag,
      remission_flag,
      censor_main,
      subgroup_value_assigned = if (identical(assignment_variable, "__overall__")) {
        "overall"
      } else {
        as.character(.data[[assignment_variable]])
      }
    )

  out <- tidyr::crossing(
    subject_base,
    horizons_tbl %>%
      select(
        horizon_id,
        horizon_year,
        horizon_days,
        support_tier,
        support_tier_standard,
        interpretation_tier,
        primary_supported_flag,
        horizon_evidence_class,
        claim_restriction_flag,
        interpretation_note
      )
  ) %>%
    left_join(
      risk_lookup %>%
        distinct(assignment_value, horizon_id, horizon_year, horizon_days, km_survival, km_risk),
      by = c("subgroup_value_assigned" = "assignment_value", "horizon_id", "horizon_year", "horizon_days")
    ) %>%
    mutate(
      dataset = meta$dataset,
      dataset_key = meta$dataset_key,
      base_dataset_key = meta$base_dataset_key,
      analysis_structure = meta$analysis_structure,
      risk_scale = meta$risk_scale,
      model_class = "KM_benchmark",
      model_family = "KM",
      benchmark_id = meta$benchmark_id,
      benchmark_label = meta$benchmark_label,
      benchmark_scope = meta$benchmark_scope,
      benchmark_structure = meta$benchmark_structure,
      site_branch = meta$site_branch,
      subgroup_kind = meta$subgroup_kind,
      subgroup_variable = null_to_na_character(meta$subgroup_variable),
      risk_assignment_rule = meta$risk_assignment_rule,
      source_stage = "Stage 3"
    ) %>%
    select(
      dataset,
      dataset_key,
      base_dataset_key,
      analysis_structure,
      risk_scale,
      model_class,
      model_family,
      benchmark_id,
      benchmark_label,
      benchmark_scope,
      benchmark_structure,
      site_branch,
      subgroup_kind,
      subgroup_variable,
      subgroup_value_assigned,
      unique_person_id,
      site,
      sex_num,
      sex_label,
      days_followup,
      time_year,
      status_num,
      event_main,
      right_censor_flag,
      remission_flag,
      censor_main,
      horizon_id,
      horizon_year,
      horizon_days,
      support_tier,
      support_tier_standard,
      interpretation_tier,
      primary_supported_flag,
      horizon_evidence_class,
      claim_restriction_flag,
      interpretation_note,
      km_survival,
      km_risk,
      risk_assignment_rule,
      source_stage
    )

  if (anyNA(out$km_survival) || anyNA(out$km_risk)) {
    stop(
      sprintf("[%s - %s] Missing KM survival/risk values after benchmark assignment.",
              meta$dataset_key, meta$benchmark_id),
      call. = FALSE
    )
  }

  out
}

## 🟠 Define: primary and optional KM builders ===============================
make_overall_benchmark <- function(df, dataset_name, horizons_tbl, risk_scale_label) {
  fit <- survival::survfit(
    survival::Surv(time_year, event_main) ~ 1,
    data = df
  )

  risk_lookup <- extract_survfit_lookup(
    fit = fit,
    horizons_tbl = horizons_tbl,
    strata_variable = NULL
  )

  meta <- list(
    dataset = dataset_name,
    dataset_key = make_stage3_dataset_key(dataset_name, site_branch = "site_free", benchmark_scope = "primary"),
    base_dataset_key = dataset_name,
    analysis_structure = make_stage3_analysis_structure(dataset_name, site_branch = "site_free", benchmark_scope = "primary"),
    risk_scale = risk_scale_label,
    benchmark_id = "km_overall",
    benchmark_label = if (dataset_name == "merged") "Overall KM benchmark (merged site-free)" else "Overall KM benchmark",
    benchmark_scope = "primary",
    benchmark_structure = "overall",
    site_branch = "site_free",
    subgroup_kind = "overall",
    subgroup_variable = NULL,
    assignment_variable = "__overall__",
    site_effect_interpretation = if (dataset_name == "merged") {
      "merged_site_free_benchmark_no_site_specific_structural_assignment"
    } else {
      "single_site_dataset"
    },
    risk_assignment_rule = if (dataset_name == "merged") {
      "overall merged KM risk applied to all subjects in the site-free merged benchmark"
    } else {
      "overall cohort-level KM risk applied to all subjects in the dataset"
    }
  )

  list(
    fit = fit,
    registry = build_benchmark_registry_row(df, risk_lookup, meta),
    subject_predictions = build_subject_predictions(df, horizons_tbl, risk_lookup, meta)
  )
}

make_site_adjusted_benchmark <- function(df, dataset_name, horizons_tbl, risk_scale_label) {
  fit <- survival::survfit(
    survival::Surv(time_year, event_main) ~ site,
    data = df
  )

  risk_lookup <- extract_survfit_lookup(
    fit = fit,
    horizons_tbl = horizons_tbl,
    strata_variable = "site"
  )

  meta <- list(
    dataset = dataset_name,
    dataset_key = make_stage3_dataset_key(dataset_name, site_branch = "site_adjusted", benchmark_scope = "primary_structural"),
    base_dataset_key = dataset_name,
    analysis_structure = make_stage3_analysis_structure(dataset_name, site_branch = "site_adjusted", benchmark_scope = "primary_structural"),
    risk_scale = risk_scale_label,
    benchmark_id = "km_site_adjusted",
    benchmark_label = "Site-adjusted KM benchmark",
    benchmark_scope = "primary_structural",
    benchmark_structure = "site_adjusted",
    site_branch = "site_adjusted",
    subgroup_kind = "site_adjusted",
    subgroup_variable = "site",
    assignment_variable = "site",
    site_effect_interpretation = "site_specific_structural_context_proxy_not_causal_treatment_effect",
    risk_assignment_rule = "site-specific KM risk applied by each subject's site within the merged dataset"
  )

  list(
    fit = fit,
    registry = build_benchmark_registry_row(df, risk_lookup, meta),
    subject_predictions = build_subject_predictions(df, horizons_tbl, risk_lookup, meta)
  )
}

make_optional_subgroup_benchmark <- function(df, dataset_name, horizons_tbl, subgroup_variable, subgroup_min_n_value, risk_scale_label) {
  if (!subgroup_variable %in% names(df)) {
    return(NULL)
  }

  subgroup_counts <- df %>%
    mutate(subgroup_value = as.character(.data[[subgroup_variable]])) %>%
    count(subgroup_value, name = "n_group")

  if (nrow(subgroup_counts) < 2L) {
    return(NULL)
  }

  if (any(subgroup_counts$n_group < subgroup_min_n_value)) {
    return(NULL)
  }

  fit_formula <- stats::as.formula(
    paste0("survival::Surv(time_year, event_main) ~ `", subgroup_variable, "`")
  )

  fit <- survival::survfit(
    fit_formula,
    data = df
  )

  risk_lookup <- extract_survfit_lookup(
    fit = fit,
    horizons_tbl = horizons_tbl,
    strata_variable = subgroup_variable
  )

  meta <- list(
    dataset = dataset_name,
    dataset_key = make_stage3_dataset_key(dataset_name, site_branch = "site_free", benchmark_scope = "sensitivity", subgroup_variable = subgroup_variable),
    base_dataset_key = dataset_name,
    analysis_structure = make_stage3_analysis_structure(dataset_name, site_branch = "site_free", benchmark_scope = "sensitivity"),
    risk_scale = risk_scale_label,
    benchmark_id = paste0("km_subgroup_", subgroup_variable),
    benchmark_label = paste("Subgroup KM benchmark:", subgroup_variable),
    benchmark_scope = "sensitivity",
    benchmark_structure = "subgroup",
    site_branch = "site_free",
    subgroup_kind = "subgroup",
    subgroup_variable = subgroup_variable,
    assignment_variable = subgroup_variable,
    site_effect_interpretation = if (dataset_name == "merged") {
      "merged_subgroup_sensitivity_not_causal_treatment_effect"
    } else {
      "single_site_subgroup_sensitivity"
    },
    risk_assignment_rule = paste0("subgroup-specific KM risk applied by `", subgroup_variable, "` within dataset")
  )

  list(
    fit = fit,
    registry = build_benchmark_registry_row(df, risk_lookup, meta),
    subject_predictions = build_subject_predictions(df, horizons_tbl, risk_lookup, meta)
  )
}

# 🔴 Define: standardized Stage-3 summary builders ===============================
## 🟠 Define: KM group-risk table builder ===============================
make_group_risk_table <- function(subject_predictions) {
  subject_predictions %>%
    group_by(
      dataset,
      dataset_key,
      base_dataset_key,
      analysis_structure,
      risk_scale,
      model_class,
      model_family,
      benchmark_id,
      benchmark_label,
      benchmark_scope,
      benchmark_structure,
      site_branch,
      subgroup_kind,
      subgroup_variable,
      subgroup_value_assigned,
      horizon_id,
      horizon_year,
      horizon_days,
      support_tier,
      support_tier_standard,
      interpretation_tier,
      primary_supported_flag,
      horizon_evidence_class,
      claim_restriction_flag,
      interpretation_note,
      km_survival,
      km_risk,
      risk_assignment_rule
    ) %>%
    summarise(
      n_subjects_in_group = n_distinct(unique_person_id),
      .groups = "drop"
    ) %>%
    arrange(
      match(dataset, dataset_order),
      dataset_key,
      benchmark_scope,
      benchmark_id,
      subgroup_value_assigned,
      horizon_year
    )
}

## 🟠 Define: threshold-based benchmark classification builder ===============================
make_classification_table <- function(subject_predictions, threshold_registry) {
  subject_predictions %>%
    tidyr::crossing(
      threshold_registry %>%
        select(
          threshold_id,
          threshold,
          threshold_label,
          positive_rule
        )
    ) %>%
    mutate(
      predicted_positive = apply_threshold_rule(km_risk, threshold),
      event_by_horizon = event_main == 1L & time_year <= horizon_year,
      known_nonevent_by_horizon = !event_by_horizon & time_year >= horizon_year,
      evaluable_by_horizon = event_by_horizon | known_nonevent_by_horizon,
      unknown_due_to_early_censoring = !evaluable_by_horizon
    ) %>%
    group_by(
      dataset,
      dataset_key,
      base_dataset_key,
      analysis_structure,
      risk_scale,
      model_class,
      model_family,
      benchmark_id,
      benchmark_label,
      benchmark_scope,
      benchmark_structure,
      site_branch,
      subgroup_kind,
      subgroup_variable,
      threshold_id,
      threshold,
      threshold_label,
      positive_rule,
      horizon_id,
      horizon_year,
      horizon_days,
      support_tier,
      support_tier_standard,
      interpretation_tier,
      primary_supported_flag,
      horizon_evidence_class,
      claim_restriction_flag,
      interpretation_note
    ) %>%
    summarise(
      group_values_used = collapse_unique_chr(subgroup_value_assigned),
      n_groups_applied = dplyr::n_distinct(subgroup_value_assigned),
      n_subjects = n(),
      n_event_by_horizon = sum(event_by_horizon),
      n_known_nonevent_by_horizon = sum(known_nonevent_by_horizon),
      n_evaluable_by_horizon = sum(evaluable_by_horizon),
      n_unknown_due_to_early_censoring = sum(unknown_due_to_early_censoring),
      n_predicted_positive = sum(predicted_positive),
      n_predicted_negative = sum(!predicted_positive),
      n_predicted_positive_evaluable = sum(predicted_positive & evaluable_by_horizon),
      n_predicted_positive_unknown = sum(predicted_positive & unknown_due_to_early_censoring),
      tp_count = sum(predicted_positive & event_by_horizon),
      fp_count = sum(predicted_positive & known_nonevent_by_horizon),
      tn_count = sum(!predicted_positive & known_nonevent_by_horizon),
      fn_count = sum(!predicted_positive & event_by_horizon),
      km_survival_value = if (dplyr::n_distinct(km_survival) == 1L) dplyr::first(km_survival) else NA_real_,
      km_risk_value = if (dplyr::n_distinct(km_risk) == 1L) dplyr::first(km_risk) else NA_real_,
      km_survival_mean = mean(km_survival),
      km_survival_min = min(km_survival),
      km_survival_max = max(km_survival),
      km_risk_mean = mean(km_risk),
      km_risk_min = min(km_risk),
      km_risk_max = max(km_risk),
      .groups = "drop"
    ) %>%
    mutate(
      positive_classification_rate = safe_divide(n_predicted_positive, n_subjects),
      km_positive_classification_rate = positive_classification_rate,
      positive_classification_rate_evaluable = safe_divide(n_predicted_positive_evaluable, n_evaluable_by_horizon),
      false_positive_rate = safe_divide(fp_count, n_known_nonevent_by_horizon),
      false_positive_burden = safe_divide(fp_count, n_subjects),
      false_positive_burden_all = false_positive_burden,
      false_positive_burden_known_nonevent = safe_divide(fp_count, n_known_nonevent_by_horizon),
      false_positives_per_100 = 100 * false_positive_burden,
      FP100 = false_positives_per_100,
      specificity = safe_divide(tn_count, tn_count + fp_count),
      ppv = safe_divide(tp_count, tp_count + fp_count),
      tpr = safe_divide(tp_count, tp_count + fn_count),
      fdp = ifelse(is.na(ppv), NA_real_, 1 - ppv),
      benchmark_positive_pattern = case_when(
        is.na(positive_classification_rate) ~ NA_character_,
        positive_classification_rate == 1 ~ "all_high_risk",
        positive_classification_rate == 0 ~ "all_low_risk",
        TRUE ~ "mixed_group_assigned_risk"
      ),
      evaluation_rule = "event by horizon OR observed event-free through horizon; early censoring before horizon excluded from outcome-based metrics"
    ) %>%
    arrange(
      match(dataset, dataset_order),
      dataset_key,
      benchmark_scope,
      benchmark_id,
      horizon_year,
      threshold
    )
}

## 🟠 Define: Stage-3 metadata builder ===============================
make_stage3_metadata <- function(
    stage1_files,
    stage1_bundle,
    threshold_registry,
    subgroup_flag,
    subgroup_vars,
    subgroup_min_n_value,
    prediction_file_name,
    classification_file_name
) {
  shared_master_spec <- if (!is.null(stage1_bundle$config$shared_master_spec_file)) stage1_bundle$config$shared_master_spec_file else "Integrated_Modeling_Master_Specification_English_REVISED_v5.md"
  code_rules_spec <- if (!is.null(stage1_bundle$config$code_rules_file)) stage1_bundle$config$code_rules_file else "Rules Before Generating R Code_🇬🇧ENG.md"
  data_dictionary_spec <- if (!is.null(stage1_bundle$config$data_dictionary_file)) stage1_bundle$config$data_dictionary_file else "3.Data Dictionary_🇬🇧ENG.md"

  tibble(
    metadata_group = c(
      "stage", "stage", "stage", "stage",
      "documents", "documents", "documents",
      "inputs", "inputs", "inputs", "inputs", "inputs", "inputs", "inputs", "inputs", "inputs",
      "event", "event", "event",
      "horizons", "thresholds",
      "benchmark", "benchmark", "benchmark", "benchmark",
      "subgroup", "subgroup", "subgroup",
      "outputs", "outputs",
      "stage1"
    ),
    metadata_name = c(
      "stage_name", "stage_role", "model_fitting_allowed", "risk_scale",
      "canonical_common_rules", "governing_code_rules", "governing_data_dictionary",
      "data_path", "export_path", "stage1_bundle_file", "stage1_analysis_datasets_file",
      "stage1_manifest_file", "stage1_dataset_registry_file", "stage1_horizon_registry_file",
      "stage1_threshold_registry_file", "stage1_modeling_registry_file",
      "event_definition", "main_censoring_definition", "threshold_positive_rule",
      "horizon_vector", "threshold_vector",
      "primary_benchmark_rule", "merged_site_structure_rule", "benchmark_output_source_of_truth", "classification_output_file",
      "run_subgroup_km", "subgroup_variables", "subgroup_min_n",
      "subject_prediction_file", "group_risk_file",
      "stage1_created_at"
    ),
    metadata_value = c(
      "Stage 3 KM benchmark classifier",
      "Operationalize overall KM as the formal horizon-specific benchmark classifier and retain merged site-adjusted structural KM as a comparator.",
      "TRUE_KM_only",
      as.character(stage1_bundle$config$main_risk_scale),
      shared_master_spec,
      code_rules_spec,
      data_dictionary_spec,
      normalize_existing_path(data_path),
      normalize_existing_path(export_path),
      stage1_files$bundle,
      stage1_files$analysis_datasets,
      stage1_files$manifest,
      stage1_files$dataset_registry,
      stage1_files$horizon_registry,
      stage1_files$threshold_registry,
      stage1_files$modeling_registry,
      "status_num == 1",
      "status_num %in% c(0, 2)",
      "Classify as high risk when predicted risk >= threshold",
      paste(sort(unique(as.integer(stage1_bundle$config$common_horizons_year))), collapse = ","),
      paste(format(sort(unique(as.numeric(threshold_registry$threshold))), trim = TRUE, scientific = FALSE), collapse = ","),
      "Overall KM is the primary benchmark classifier on the transition-only main risk scale.",
      "Merged outputs retain both a site-free overall benchmark and a site-adjusted structural benchmark.",
      prediction_file_name,
      classification_file_name,
      as.character(subgroup_flag),
      paste(subgroup_vars, collapse = "|"),
      as.character(subgroup_min_n_value),
      prediction_file_name,
      "stage3_km_group_risk_table.csv",
      as.character(stage1_bundle$created_at)
    )
  )
}

# 🔴 Read: Stage-1 backbone artifacts ===============================
## 🟠 Read: required CSV and RDS inputs ===============================
required_stage1_files <- list(
  bundle = stage1_bundle_file,
  analysis_datasets = stage1_analysis_datasets_file,
  manifest = stage1_manifest_file,
  dataset_registry = stage1_dataset_registry_file,
  horizon_registry = stage1_horizon_registry_file,
  threshold_registry = stage1_threshold_registry_file,
  metadata_registry = stage1_metadata_registry_file,
  scaling_registry = stage1_scaling_registry_file,
  modeling_registry = stage1_modeling_registry_file
)

invisible(purrr::imap(required_stage1_files, ~ assert_existing_file(.x, paste("Required Stage 1 file", .y))))

stage1_bundle <- readRDS(stage1_bundle_file)
analysis_datasets <- readRDS(stage1_analysis_datasets_file)
stage1_manifest <- read_stage1_csv(stage1_manifest_file)
stage1_dataset_registry <- read_stage1_csv(stage1_dataset_registry_file)
stage1_horizon_registry <- read_stage1_csv(stage1_horizon_registry_file)
stage1_threshold_registry <- read_stage1_csv(stage1_threshold_registry_file)
stage1_metadata_registry <- read_stage1_csv(stage1_metadata_registry_file)
stage1_scaling_registry <- read_stage1_csv(stage1_scaling_registry_file)
stage1_modeling_registry <- read_stage1_csv(stage1_modeling_registry_file)

# 🔴 Validate: Stage-1 inherited common backbone ===============================
## 🟠 Check: comparability prerequisites ===============================
validate_stage1_backbone(
  bundle = stage1_bundle,
  analysis_datasets = analysis_datasets,
  dataset_registry = stage1_dataset_registry,
  horizon_registry = stage1_horizon_registry,
  threshold_registry = stage1_threshold_registry,
  metadata_registry = stage1_metadata_registry,
  scaling_registry = stage1_scaling_registry,
  modeling_registry = stage1_modeling_registry
)

horizon_registry <- stage1_horizon_registry %>%
  mutate(
    horizon_year = as.integer(horizon_year),
    horizon_days = as.numeric(horizon_days),
    support_tier = as.character(support_tier),
    horizon_evidence_class = as.character(horizon_evidence_class),
    claim_restriction_flag = as.character(claim_restriction_flag),
    interpretation_note = as.character(interpretation_note)
  )

if ("support_tier_standard" %in% names(horizon_registry)) {
  horizon_registry <- horizon_registry %>% mutate(support_tier_standard = as.character(support_tier_standard))
} else {
  horizon_registry <- horizon_registry %>% mutate(support_tier_standard = support_tier)
}

if ("interpretation_tier" %in% names(horizon_registry)) {
  horizon_registry <- horizon_registry %>% mutate(interpretation_tier = as.character(interpretation_tier))
} else {
  horizon_registry <- horizon_registry %>% mutate(
    interpretation_tier = if_else(support_tier == "primary_supported", "primary-supported", support_tier)
  )
}

if ("primary_supported_flag" %in% names(horizon_registry)) {
  horizon_registry <- horizon_registry %>% mutate(primary_supported_flag = as.logical(primary_supported_flag))
} else {
  horizon_registry <- horizon_registry %>% mutate(primary_supported_flag = support_tier == "primary_supported")
}

horizon_registry <- horizon_registry %>%
  arrange(match(dataset, dataset_order), horizon_year)

threshold_registry <- stage1_threshold_registry %>%
  mutate(
    threshold = as.numeric(threshold)
  ) %>%
  arrange(threshold)

main_risk_scale <- as.character(stage1_bundle$config$main_risk_scale)

# 🔴 Fit: Stage-3 KM benchmark variants ===============================
## 🟠 Fit: overall KM benchmarks for all datasets ===============================
overall_results <- purrr::map(
  dataset_order,
  function(dataset_name) {
    make_overall_benchmark(
      df = analysis_datasets[[dataset_name]],
      dataset_name = dataset_name,
      horizons_tbl = horizon_registry %>% filter(dataset == dataset_name),
      risk_scale_label = main_risk_scale
    )
  }
)
names(overall_results) <- dataset_order

## 🟠 Fit: merged site-adjusted KM benchmark ===============================
site_adjusted_result <- make_site_adjusted_benchmark(
  df = analysis_datasets[["merged"]],
  dataset_name = "merged",
  horizons_tbl = horizon_registry %>% filter(dataset == "merged"),
  risk_scale_label = main_risk_scale
)

## 🟠 Fit: optional subgroup KM sensitivity benchmarks ===============================
subgroup_results <- list()

if (isTRUE(run_subgroup_km)) {
  for (dataset_name in dataset_order) {
    for (subgroup_var in subgroup_variables) {
      subgroup_candidate <- make_optional_subgroup_benchmark(
        df = analysis_datasets[[dataset_name]],
        dataset_name = dataset_name,
        horizons_tbl = horizon_registry %>% filter(dataset == dataset_name),
        subgroup_variable = subgroup_var,
        subgroup_min_n_value = subgroup_min_n,
        risk_scale_label = main_risk_scale
      )

      if (!is.null(subgroup_candidate)) {
        subgroup_results[[length(subgroup_results) + 1L]] <- subgroup_candidate
      }
    }
  }
}

# 🔴 Summarize: threshold-based benchmark burden outputs ===============================
## 🟠 Collect: comparison-ready Stage-3 long objects ===============================
all_benchmark_results <- purrr::compact(c(
  overall_results,
  list(site_adjusted_result),
  subgroup_results
))

benchmark_fit_names <- purrr::map_chr(
  all_benchmark_results,
  ~ paste0(.x$registry$dataset_key[[1]], "__", .x$registry$benchmark_id[[1]])
)

km_fit_objects <- stats::setNames(
  object = purrr::map(all_benchmark_results, "fit"),
  nm = benchmark_fit_names
)

km_benchmark_registry <- bind_rows(purrr::map(all_benchmark_results, "registry")) %>%
  arrange(
    match(dataset, dataset_order),
    dataset_key,
    benchmark_scope,
    benchmark_id
  )

km_subject_horizon_predictions <- bind_rows(purrr::map(all_benchmark_results, "subject_predictions")) %>%
  arrange(
    match(dataset, dataset_order),
    dataset_key,
    benchmark_scope,
    benchmark_id,
    unique_person_id,
    horizon_year
  )

km_group_risk_table <- make_group_risk_table(km_subject_horizon_predictions)

km_benchmark_classification <- make_classification_table(
  subject_predictions = km_subject_horizon_predictions,
  threshold_registry = threshold_registry
)

# 🔴 Export: Stage-3 CSV and RDS deliverables ===============================
## 🟠 Export: compact file set in one folder ===============================
km_registry_file <- "stage3_km_benchmark_registry.csv"
km_subject_predictions_file <- "stage3_km_subject_horizon_predictions.csv"
km_group_risk_file <- "stage3_km_group_risk_table.csv"
km_classification_file <- "stage3_km_benchmark_classification.csv"
stage3_metadata_file <- "stage3_km_stage_metadata.csv"
km_fit_rds_file <- "stage3_km_fit_objects.rds"
km_bundle_rds_file <- "stage3_km_output_bundle.rds"
stage3_manifest_file <- "stage3_export_manifest.csv"

stage3_metadata <- make_stage3_metadata(
  stage1_files = list(
    bundle = normalize_existing_path(stage1_bundle_file),
    analysis_datasets = normalize_existing_path(stage1_analysis_datasets_file),
    manifest = normalize_existing_path(stage1_manifest_file),
    dataset_registry = normalize_existing_path(stage1_dataset_registry_file),
    horizon_registry = normalize_existing_path(stage1_horizon_registry_file),
    threshold_registry = normalize_existing_path(stage1_threshold_registry_file),
    modeling_registry = normalize_existing_path(stage1_modeling_registry_file)
  ),
  stage1_bundle = stage1_bundle,
  threshold_registry = threshold_registry,
  subgroup_flag = run_subgroup_km,
  subgroup_vars = subgroup_variables,
  subgroup_min_n_value = subgroup_min_n,
  prediction_file_name = km_subject_predictions_file,
  classification_file_name = km_classification_file
)

readr::write_csv(km_benchmark_registry, file.path(export_path, km_registry_file))
readr::write_csv(km_subject_horizon_predictions, file.path(export_path, km_subject_predictions_file))
readr::write_csv(km_group_risk_table, file.path(export_path, km_group_risk_file))
readr::write_csv(km_benchmark_classification, file.path(export_path, km_classification_file))
readr::write_csv(stage3_metadata, file.path(export_path, stage3_metadata_file))

saveRDS(km_fit_objects, file.path(export_path, km_fit_rds_file))

stage3_output_bundle <- list(
  stage = "Stage 3",
  created_at = as.character(Sys.time()),
  session_info = utils::sessionInfo(),
  config = list(
    data_path = normalize_existing_path(data_path),
    export_path = normalize_existing_path(export_path),
    main_risk_scale = main_risk_scale,
    run_subgroup_km = run_subgroup_km,
    subgroup_variables = subgroup_variables,
    subgroup_min_n = subgroup_min_n
  ),
  stage1_input_files = required_stage1_files,
  stage1_horizon_registry = horizon_registry,
  stage1_threshold_registry = threshold_registry,
  benchmark_registry = km_benchmark_registry,
  subject_horizon_predictions = km_subject_horizon_predictions,
  group_risk_table = km_group_risk_table,
  benchmark_classification = km_benchmark_classification,
  stage3_metadata = stage3_metadata
)

saveRDS(stage3_output_bundle, file.path(export_path, km_bundle_rds_file))

stage3_export_manifest <- tibble(
  file_name = c(
    km_registry_file,
    km_subject_predictions_file,
    km_group_risk_file,
    km_classification_file,
    stage3_metadata_file,
    km_fit_rds_file,
    km_bundle_rds_file,
    stage3_manifest_file
  ),
  object_name = c(
    "km_benchmark_registry",
    "km_subject_horizon_predictions",
    "km_group_risk_table",
    "km_benchmark_classification",
    "stage3_metadata",
    "km_fit_objects",
    "stage3_output_bundle",
    "stage3_export_manifest"
  ),
  description = c(
    "Registry of Stage 3 KM benchmark variants, including merged site-free and site-adjusted structures",
    "Subject-by-horizon KM prediction table; main source of truth for Stage 3 benchmark outputs",
    "Group-level KM survival and risk table derived from the subject-horizon source object",
    "Long-format KM benchmark classification table aligned to later non-cure, cure, and Bayesian outputs",
    "Stage 3 metadata and interpretation-governance notes",
    "Saved survfit objects for KM benchmark variants",
    "Reusable Stage 3 KM bundle containing comparison-ready outputs",
    "Manifest of all Stage 3 exported files"
  ),
  file_path = file.path(
    export_path,
    c(
      km_registry_file,
      km_subject_predictions_file,
      km_group_risk_file,
      km_classification_file,
      stage3_metadata_file,
      km_fit_rds_file,
      km_bundle_rds_file,
      stage3_manifest_file
    )
  )
)

readr::write_csv(stage3_export_manifest, file.path(export_path, stage3_manifest_file))
