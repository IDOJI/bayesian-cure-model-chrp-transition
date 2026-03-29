# 🔴 Configure: paths and cache controls ===============================
rm(list = ls())

data_path <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/stage1_Backbone lock'
export_path <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/stage5_Individualized no-cure comparator'

stage1_bundle_file <- file.path(data_path, "stage1_backbone_bundle.rds")
stage1_datasets_file <- file.path(data_path, "stage1_analysis_datasets.rds")
stage1_dataset_registry_file <- file.path(data_path, "stage1_dataset_registry.csv")
stage1_formula_registry_file <- file.path(data_path, "stage1_formula_registry.csv")
stage1_horizon_registry_file <- file.path(data_path, "stage1_horizon_registry.csv")
stage1_threshold_registry_file <- file.path(data_path, "stage1_threshold_registry.csv")
stage1_scaling_registry_file <- file.path(data_path, "stage1_scaling_registry.csv")
stage1_metadata_registry_file <- file.path(data_path, "stage1_metadata_registry.csv")

reuse_existing_stage5_core <- TRUE
force_recompute_stage5_core <- FALSE
force_regenerate_plots <- FALSE

calibration_group_count <- 10L
prediction_horizons_for_plots <- c(1L, 2L, 5L, 10L)
probability_clip_epsilon <- 1e-06
time_origin_epsilon_year <- 1e-08
qc_tolerance <- 1e-08
time_range_score_min_sd <- 1e-10
uno_harrell_discordance_tolerance <- 0.20
validation_level_default <- "apparent"
optimism_correction_method_default <- "none"
calibration_target_label <- "IPCW fixed-horizon binary outcome"
mean_calibration_reference_label <- "KM observed risk at horizon"
late_horizon_warning_start_year <- 8L
late_horizon_known_count_warning_min <- 55L

stage_name <- "Stage 5"
stage_role <- "individualized_no_cure_comparator"
stage_spec_version <- "stage5_noncure_revised_v7"

plot_width_in <- 12
plot_height_in <- 8
plot_dpi <- 320

model_registry_file <- file.path(export_path, "stage5_model_registry.csv")
risk_long_file <- file.path(export_path, "stage5_subject_horizon_risk_long.csv")
performance_long_file <- file.path(export_path, "stage5_model_performance_long.csv")
calibration_bins_file <- file.path(export_path, "stage5_calibration_bins_long.csv")
calibration_diagnostics_file <- file.path(export_path, "stage5_calibration_diagnostics_long.csv")
qc_summary_file <- file.path(export_path, "stage5_qc_summary.csv")
family_fit_summary_file <- file.path(export_path, "stage5_family_fit_summary.csv")
threshold_classification_file <- file.path(export_path, "stage5_threshold_classification_long.csv")
false_positive_burden_file <- file.path(export_path, "stage5_false_positive_burden_long.csv")
decision_curve_file <- file.path(export_path, "stage5_decision_curve_long.csv")

plot_source_risk_file <- file.path(export_path, "stage5_plot_source_risk_summary.csv")
plot_source_auc_file <- file.path(export_path, "stage5_plot_source_auc_summary.csv")
plot_source_brier_file <- file.path(export_path, "stage5_plot_source_brier_summary.csv")
plot_source_calibration_file <- file.path(export_path, "stage5_plot_source_calibration_selected_horizons.csv")
plot_source_net_benefit_file <- file.path(export_path, "stage5_plot_source_net_benefit_selected_horizons.csv")

pdf_file <- file.path(export_path, "stage5_summary_plots.pdf")
plot_file_png_risk <- file.path(export_path, "stage5_plot_risk_summary.png")
plot_file_png_auc <- file.path(export_path, "stage5_plot_auc_summary.png")
plot_file_png_brier <- file.path(export_path, "stage5_plot_brier_summary.png")
plot_file_png_calibration <- file.path(export_path, "stage5_plot_calibration_selected_horizons.png")
plot_file_png_net_benefit <- file.path(export_path, "stage5_plot_net_benefit_selected_horizons.png")

fitted_models_file <- file.path(export_path, "stage5_fitted_models.rds")
results_bundle_file <- file.path(export_path, "stage5_results_bundle.rds")
export_manifest_file <- file.path(export_path, "stage5_export_manifest.csv")

formula_plot_levels <- c(
  "Base",
  "Interaction",
  "Site-added",
  "Site + interaction",
  "Reference"
)

model_family_plot_levels <- c(
  "exponential",
  "weibull",
  "lognormal",
  "loglogistic",
  "coxph (no tail extrapolation)",
  "treat_all",
  "treat_none",
  "observed_km"
)

threshold_outcome_metric_names <- c(
  "true_positive_count_ipcw",
  "false_positive_count_ipcw",
  "false_negative_count_ipcw",
  "true_negative_count_ipcw",
  "false_positive_burden_all",
  "false_positive_burden_nonevents",
  "false_positive_per_100",
  "sensitivity",
  "specificity",
  "ppv",
  "npv",
  "false_discovery_proportion",
  "event_prevalence_ipcw",
  "event_rate_per_subject_ipcw",
  "net_benefit",
  "net_benefit_treat_all",
  "net_reduction_unnecessary_per_100"
)

# 🔴 Attach: packages and options ===============================
required_packages <- c(
  "dplyr",
  "tidyr",
  "purrr",
  "tibble",
  "readr",
  "ggplot2",
  "survival",
  "flexsurv"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    sprintf(
      "Install required packages before running this script: %s",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(survival)
  library(flexsurv)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

dataset_order <- c("PNU", "SNU", "merged")
has_survAUC <- FALSE

# 🔴 Define: scalar and text helpers ===============================
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

clip_prob <- function(x, eps = probability_clip_epsilon) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

safe_char <- function(x) {
  as.character(x %||% NA_character_)
}

safe_num <- function(x) {
  as.numeric(x %||% NA_real_)
}

safe_max_abs <- function(x) {
  x <- abs(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  max(x)
}

safe_sd <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) <= 1L) {
    return(NA_real_)
  }
  stats::sd(x)
}

safe_divide <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

weighted_mean_safe <- function(x, w) {
  valid <- !is.na(x) & !is.na(w) & is.finite(x) & is.finite(w) & w > 0
  if (!any(valid)) {
    return(NA_real_)
  }
  sum(x[valid] * w[valid]) / sum(w[valid])
}

collapse_nonempty_chr <- function(x, sep = "; ") {
  x <- trimws(as.character(x))
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x) == 0L) {
    return(NA_character_)
  }
  paste(x, collapse = sep)
}

normalize_existing_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
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

strip_formula_rhs <- function(x) {
  x_chr <- trimws(as.character(x))
  gsub("^~\\s*", "", x_chr)
}

bind_rows_from_df_lists <- function(..., object_name = "input_objects") {
  raw_inputs <- list(...)
  out <- list()

  for (obj in raw_inputs) {
    if (is.null(obj)) {
      next
    }

    if (inherits(obj, "data.frame")) {
      out[[length(out) + 1L]] <- obj
      next
    }

    if (is.list(obj)) {
      obj_flat <- purrr::compact(obj)
      out <- c(out, obj_flat)
      next
    }

    out[[length(out) + 1L]] <- obj
  }

  if (length(out) == 0L) {
    return(tibble())
  }

  is_df <- vapply(out, inherits, logical(1), what = "data.frame")

  if (!all(is_df)) {
    bad_positions <- which(!is_df)
    stop(
      sprintf(
        "`%s` contains at least one non-data.frame object at positions: %s",
        object_name,
        paste(bad_positions, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  dplyr::bind_rows(out)
}

collapse_dataset_key_suffixes <- function(df) {
  if (is.null(df) || !inherits(df, "data.frame")) {
    return(df)
  }

  out <- tibble::as_tibble(df)

  if (!"dataset_key" %in% names(out)) {
    out$dataset_key <- NA_character_
  }

  if ("dataset_key.x" %in% names(out)) {
    out$dataset_key <- dplyr::coalesce(out$dataset_key, as.character(out[["dataset_key.x"]]))
  }

  if ("dataset_key.y" %in% names(out)) {
    out$dataset_key <- dplyr::coalesce(out$dataset_key, as.character(out[["dataset_key.y"]]))
  }

  out %>%
    select(-any_of(c("dataset_key.x", "dataset_key.y")))
}

assert_no_dataset_key_suffixes <- function(df, object_name) {
  bad_cols <- intersect(c("dataset_key.x", "dataset_key.y"), names(df))
  if (length(bad_cols) > 0L) {
    stop(
      sprintf(
        "`%s` still contains unresolved dataset_key suffix columns: %s",
        object_name,
        paste(bad_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

safe_read_csv_file <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(
    collapse_dataset_key_suffixes(
      readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
    ),
    error = function(e) NULL
  )
}

safe_read_rds_file <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(readRDS(path), error = function(e) NULL)
}

make_model_family_plot_label <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x == "coxph" ~ "coxph (no tail extrapolation)",
    TRUE ~ x
  )
}

apply_model_family_plot_label <- function(df) {
  df %>%
    mutate(
      model_family_plot_label = factor(
        make_model_family_plot_label(model_family),
        levels = model_family_plot_levels
      )
    )
}

make_dataset_formula_panel <- function(dataset, formula_label) {
  paste0(dataset, " | ", formula_label)
}

apply_dataset_formula_panel <- function(df) {
  panel_levels <- df %>%
    distinct(dataset, formula_label) %>%
    mutate(formula_label = if_else(is.na(formula_label), "Reference", formula_label)) %>%
    arrange(
      match(dataset, dataset_order),
      match(formula_label, formula_plot_levels)
    ) %>%
    mutate(dataset_formula_panel = make_dataset_formula_panel(dataset, formula_label)) %>%
    pull(dataset_formula_panel)

  df %>%
    mutate(
      formula_label = if_else(is.na(formula_label), "Reference", formula_label),
      dataset_formula_panel = factor(
        make_dataset_formula_panel(dataset, formula_label),
        levels = unique(panel_levels)
      )
    )
}

make_extrapolation_role <- function(model_family) {
  dplyr::case_when(
    model_family == "coxph" ~ "observed_window_comparator",
    model_family %in% c("exponential", "weibull", "lognormal", "loglogistic") ~ "parametric_tail_extrapolator",
    model_family %in% c("observed_km", "treat_all", "treat_none") ~ "not_applicable",
    TRUE ~ "not_applicable"
  )
}

make_tail_extrapolation_capable_flag <- function(model_family) {
  as.integer(model_family %in% c("exponential", "weibull", "lognormal", "loglogistic"))
}

make_likelihood_basis <- function(model_family) {
  dplyr::case_when(
    model_family == "coxph" ~ "partial_likelihood",
    model_family %in% c("exponential", "weibull", "lognormal", "loglogistic") ~ "full_likelihood",
    TRUE ~ "not_applicable"
  )
}

make_information_criterion_comparable_flag <- function(model_family) {
  as.integer(model_family %in% c("exponential", "weibull", "lognormal", "loglogistic"))
}

binary_outcome_support_from_ipcw <- function(ipcw_info) {
  as.integer(isTRUE(ipcw_info$binary_outcome_support_flag))
}

choose_time_range_score <- function(pred_result, anchor_horizon_index, anchor_horizon_year) {
  anchor_vec <- as.numeric(pred_result$risk[, anchor_horizon_index])

  list(
    score = anchor_vec,
    score_type = paste0("anchor_risk_year_", anchor_horizon_year)
  )
}

# 🔴 Define: Stage-1 bundle readers ===============================
read_stage1_bundle <- function(bundle_path, datasets_path = NULL) {
  if (!file.exists(bundle_path)) {
    stop(sprintf("Stage 1 backbone bundle not found: %s", bundle_path), call. = FALSE)
  }

  bundle <- readRDS(bundle_path)

  if (!is.list(bundle)) {
    stop("`stage1_backbone_bundle.rds` must be an R list.", call. = FALSE)
  }

  if (is.null(bundle$datasets) && !is.null(datasets_path) && file.exists(datasets_path)) {
    bundle$datasets <- readRDS(datasets_path)
  }

  bundle
}

extract_stage1_piece <- function(bundle, piece_name) {
  find_first_data_frame <- function(x) {
    if (inherits(x, "data.frame")) {
      return(tibble::as_tibble(x))
    }

    if (!is.list(x) || length(x) == 0L) {
      return(NULL)
    }

    preferred_names <- c("data", "df", "table", "registry", "value")

    for (nm in preferred_names) {
      if (!is.null(x[[nm]])) {
        out <- find_first_data_frame(x[[nm]])
        if (!is.null(out)) {
          return(out)
        }
      }
    }

    for (idx in seq_along(x)) {
      out <- find_first_data_frame(x[[idx]])
      if (!is.null(out)) {
        return(out)
      }
    }

    NULL
  }

  csv_lookup <- list(
    dataset_registry = stage1_dataset_registry_file,
    formula_registry = stage1_formula_registry_file,
    horizon_registry = stage1_horizon_registry_file,
    threshold_registry = stage1_threshold_registry_file,
    scaling_registry = stage1_scaling_registry_file,
    metadata_registry = stage1_metadata_registry_file
  )

  csv_path <- csv_lookup[[piece_name]] %||% NULL

  if (!is.null(csv_path) && file.exists(csv_path)) {
    piece <- readr::read_csv(csv_path, show_col_types = FALSE)
  } else {
    piece <- bundle$registries[[piece_name]] %||% NULL

    if (is.null(piece)) {
      stop(sprintf("Stage 1 bundle is missing registry `%s`.", piece_name), call. = FALSE)
    }

    piece <- find_first_data_frame(piece)
  }

  if (!inherits(piece, "data.frame")) {
    stop(
      sprintf(
        "Stage 1 registry `%s` could not be recovered as a rectangular data frame.",
        piece_name
      ),
      call. = FALSE
    )
  }

  list_cols <- names(piece)[vapply(piece, is.list, logical(1))]

  if (length(list_cols) > 0L) {
    stop(
      sprintf(
        "Stage 1 registry `%s` contains list-columns that cannot be grouped safely: %s",
        piece_name,
        paste(list_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  tibble::as_tibble(piece)
}

extract_stage1_datasets <- function(bundle) {
  datasets <- bundle$datasets %||% NULL

  if (is.null(datasets)) {
    stop("Stage 1 bundle does not contain `datasets`.", call. = FALSE)
  }

  expected_names <- c("PNU", "SNU", "merged")

  if (!all(expected_names %in% names(datasets))) {
    stop("Stage 1 datasets must contain `PNU`, `SNU`, and `merged`.", call. = FALSE)
  }

  datasets
}

extract_stage1_thresholds <- function(bundle, threshold_registry) {
  out <- bundle$config$risk_thresholds %||% bundle$config$risk_thresholds %||% NULL

  if (is.null(out)) {
    out <- threshold_registry$threshold
  }

  out <- sort(unique(as.numeric(out)))

  if (length(out) == 0L || anyNA(out) || any(out <= 0 | out >= 1)) {
    stop("Stage 5 requires a valid Stage 1 threshold vector.", call. = FALSE)
  }

  out
}

extract_stage1_horizons <- function(bundle, horizon_registry) {
  out <- bundle$config$common_horizons_year %||% NULL

  if (is.null(out)) {
    out <- sort(unique(as.integer(horizon_registry$horizon_year)))
  }

  out <- as.integer(out)

  if (!identical(out, 1:10)) {
    stop("Stage 5 requires the Stage 1 common horizon grid to be exactly 1:10 years.", call. = FALSE)
  }

  out
}

# 🔴 Define: registry harmonizers ===============================
required_formula_metadata_cols <- c(
  "dataset",
  "dataset_key",
  "formula_id",
  "formula_name",
  "formula_label",
  "formula_rhs",
  "formula_full",
  "uses_site",
  "uses_age_sex_interaction",
  "site_branch",
  "interaction_branch",
  "risk_scale",
  "formula_scope",
  "site_term_interpretation"
)

required_horizon_metadata_cols <- c(
  "dataset",
  "dataset_key",
  "horizon_year",
  "horizon_id",
  "support_tier",
  "support_tier_standard",
  "interpretation_tier",
  "primary_supported_flag",
  "horizon_evidence_class",
  "claim_restriction_flag",
  "interpretation_note"
)

prepare_formula_registry_core <- function(formula_registry) {
  formula_registry %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      dataset_key = as.character(dataset_key %||% dataset),
      formula_full = as.character(formula_full %||% paste("~", formula_rhs)),
      uses_site = as.logical(uses_site),
      uses_age_sex_interaction = as.logical(uses_age_sex_interaction),
      site_branch = as.character(site_branch %||% ifelse(uses_site, "site_adjusted", "site_free")),
      interaction_branch = as.character(interaction_branch %||% ifelse(uses_age_sex_interaction, "age_sex_interaction", "no_age_sex_interaction")),
      risk_scale = as.character(risk_scale %||% "transition_only_main"),
      formula_scope = as.character(formula_scope %||% "main_transition_only_scale"),
      site_term_interpretation = as.character(site_term_interpretation %||% ifelse(uses_site, "structural_context_proxy_not_causal_treatment_effect", "not_applicable"))
    ) %>%
    select(all_of(required_formula_metadata_cols)) %>%
    distinct()
}

prepare_horizon_registry <- function(horizon_registry) {
  horizon_registry %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      dataset_key = as.character(dataset_key %||% dataset),
      horizon_year = as.integer(horizon_year),
      support_tier = as.character(support_tier),
      support_tier_standard = as.character(support_tier_standard %||% support_tier),
      interpretation_tier = as.character(interpretation_tier),
      primary_supported_flag = as.logical(primary_supported_flag),
      horizon_evidence_class = as.character(horizon_evidence_class),
      claim_restriction_flag = as.character(claim_restriction_flag),
      interpretation_note = as.character(interpretation_note)
    ) %>%
    select(all_of(required_horizon_metadata_cols)) %>%
    distinct()
}

drop_horizon_metadata_columns <- function(df) {
  drop_cols <- intersect(
    names(df),
    c(
      "horizon_id",
      "horizon_label",
      "support_tier",
      "support_tier_standard",
      "interpretation_tier",
      "primary_supported_flag",
      "horizon_evidence_class",
      "claim_restriction_flag",
      "interpretation_note"
    )
  )
  if (length(drop_cols) == 0L) {
    return(df)
  }
  df %>% select(-all_of(drop_cols))
}

join_horizon_metadata <- function(df, horizon_registry) {
  if (!all(c("dataset", "horizon_year") %in% names(df))) {
    return(df)
  }

  horizon_meta <- horizon_registry %>%
    mutate(horizon_label = paste0("Year ", horizon_year))

  overlap_cols <- intersect(
    setdiff(names(horizon_meta), c("dataset", "horizon_year")),
    names(df)
  )

  df %>%
    collapse_dataset_key_suffixes() %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      horizon_year = as.integer(horizon_year)
    ) %>%
    drop_horizon_metadata_columns() %>%
    select(-any_of(overlap_cols)) %>%
    left_join(horizon_meta, by = c("dataset", "horizon_year")) %>%
    collapse_dataset_key_suffixes()
}

join_formula_metadata <- function(df, formula_registry_core) {
  if (!all(c("dataset", "formula_id") %in% names(df))) {
    return(df)
  }

  formula_meta <- formula_registry_core %>%
    select(all_of(required_formula_metadata_cols))

  overlap_cols <- intersect(
    setdiff(names(formula_meta), c("dataset", "formula_id")),
    names(df)
  )

  df_out <- df %>%
    mutate(dataset = normalize_dataset_label(dataset)) %>%
    select(-all_of(overlap_cols)) %>%
    left_join(formula_meta, by = c("dataset", "formula_id")) %>%
    mutate(
      dataset_key = if_else(is.na(dataset_key), dataset, dataset_key),
      formula_name = if_else(is.na(formula_name) & formula_id == "reference", "reference", formula_name),
      formula_label = if_else(is.na(formula_label) & formula_id == "reference", "Reference", formula_label),
      formula_scope = if_else(is.na(formula_scope) & formula_id == "reference", "benchmark_reference", formula_scope),
      formula_rhs = if_else(is.na(formula_rhs) & formula_id == "reference", "1", formula_rhs),
      formula_full = if_else(is.na(formula_full) & formula_id == "reference", "~ 1", formula_full),
      uses_site = dplyr::coalesce(uses_site, FALSE),
      uses_age_sex_interaction = dplyr::coalesce(uses_age_sex_interaction, FALSE),
      site_branch = if_else(is.na(site_branch), "site_free", site_branch),
      interaction_branch = if_else(is.na(interaction_branch), "no_age_sex_interaction", interaction_branch),
      risk_scale = if_else(is.na(risk_scale), "transition_only_main", risk_scale),
      site_term_interpretation = if_else(is.na(site_term_interpretation), "not_applicable", site_term_interpretation)
    )

  df_out
}

horizon_label_lookup <- function(dataset_name, horizon_year, horizon_registry) {
  out <- horizon_registry %>%
    filter(dataset == dataset_name, horizon_year == !!horizon_year) %>%
    distinct(
      dataset,
      dataset_key,
      horizon_id,
      support_tier,
      support_tier_standard,
      interpretation_tier,
      primary_supported_flag,
      horizon_evidence_class,
      claim_restriction_flag,
      interpretation_note,
      .keep_all = FALSE
    ) %>%
    slice(1)

  if (nrow(out) == 0L) {
    stop(
      sprintf(
        "Missing horizon registry row for dataset `%s` and horizon `%s`.",
        dataset_name,
        horizon_year
      ),
      call. = FALSE
    )
  }

  out %>%
    mutate(
      horizon_year = as.integer(horizon_year),
      horizon_label = paste0("Year ", horizon_year)
    )
}

augment_anchor_horizon_metadata <- function(model_registry, horizon_registry) {
  if (!all(c("dataset", "anchor_horizon_year") %in% names(model_registry))) {
    return(model_registry)
  }

  anchor_meta <- horizon_registry %>%
    select(dataset, horizon_year, support_tier, horizon_evidence_class, claim_restriction_flag) %>%
    rename(
      anchor_horizon_year = horizon_year,
      time_range_support_tier = support_tier,
      time_range_horizon_evidence_class = horizon_evidence_class,
      time_range_claim_restriction_flag = claim_restriction_flag
    )

  model_registry %>%
    select(-any_of(c("time_range_support_tier", "time_range_horizon_evidence_class", "time_range_claim_restriction_flag"))) %>%
    left_join(anchor_meta, by = c("dataset", "anchor_horizon_year"))
}

# 🔴 Define: dataset harmonizers ===============================
normalize_dataset_for_stage5 <- function(df, dataset_name) {
  required_cols <- c(
    "unique_person_id",
    "site",
    "id",
    "sex_num",
    "age_exact_entry",
    "age_s",
    "days_followup",
    "time_year",
    "event_main",
    "censor_main",
    "status_num",
    "right_censor_flag",
    "remission_flag"
  )

  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "[%s] Stage 1 dataset is missing required columns: %s",
        dataset_name,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  out <- df %>%
    mutate(
      unique_person_id = as.character(unique_person_id),
      site = factor(as.character(site)),
      id = as.character(id),
      sex_num = as.integer(sex_num),
      age_exact_entry = as.numeric(age_exact_entry),
      age_s = as.numeric(age_s),
      days_followup = as.numeric(days_followup),
      time_year = as.numeric(time_year),
      time_year_model = pmax(as.numeric(time_year), time_origin_epsilon_year),
      event_main = as.integer(event_main),
      censor_main = as.integer(censor_main),
      status_num = as.integer(status_num),
      right_censor_flag = as.integer(right_censor_flag),
      remission_flag = as.integer(remission_flag)
    )

  if (anyNA(out$unique_person_id) || any(out$unique_person_id == "")) {
    stop(sprintf("[%s] Invalid `unique_person_id` values detected.", dataset_name), call. = FALSE)
  }

  if (nrow(out) != dplyr::n_distinct(out$unique_person_id)) {
    stop(sprintf("[%s] `unique_person_id` must be unique.", dataset_name), call. = FALSE)
  }

  if (anyNA(out$site)) {
    stop(sprintf("[%s] Missing `site` values detected.", dataset_name), call. = FALSE)
  }

  if (anyNA(out$sex_num) || any(!out$sex_num %in% c(0L, 1L))) {
    stop(sprintf("[%s] `sex_num` must be 0/1 only.", dataset_name), call. = FALSE)
  }

  if (anyNA(out$event_main) || any(!out$event_main %in% c(0L, 1L))) {
    stop(sprintf("[%s] `event_main` must be 0/1 only.", dataset_name), call. = FALSE)
  }

  if (anyNA(out$censor_main) || any(!out$censor_main %in% c(0L, 1L))) {
    stop(sprintf("[%s] `censor_main` must be 0/1 only.", dataset_name), call. = FALSE)
  }

  if (any(out$event_main + out$censor_main != 1L, na.rm = TRUE)) {
    stop(sprintf("[%s] `event_main + censor_main` must equal 1 for every row.", dataset_name), call. = FALSE)
  }

  if (anyNA(out$status_num) || any(!out$status_num %in% c(0L, 1L, 2L))) {
    stop(sprintf("[%s] `status_num` must be coded as 0/1/2 only.", dataset_name), call. = FALSE)
  }

  if (any(out$days_followup < 0, na.rm = TRUE) || any(out$time_year < 0, na.rm = TRUE)) {
    stop(sprintf("[%s] Negative follow-up time detected.", dataset_name), call. = FALSE)
  }

  out
}

normalize_stage1_datasets <- function(datasets) {
  expected_names <- c("PNU", "SNU", "merged")
  missing_names <- setdiff(expected_names, names(datasets))

  if (length(missing_names) > 0L) {
    stop(
      sprintf(
        "Stage 1 datasets are missing required names: %s",
        paste(missing_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  out <- stats::setNames(vector("list", length(expected_names)), expected_names)

  for (dataset_name in expected_names) {
    out[[dataset_name]] <- normalize_dataset_for_stage5(datasets[[dataset_name]], dataset_name)
  }

  out
}

assign_anchor_horizon <- function(dataset_name) {
  if (dataset_name == "PNU") {
    return(1L)
  }

  if (dataset_name %in% c("SNU", "merged")) {
    return(2L)
  }

  stop(sprintf("Unknown dataset for anchor horizon: %s", dataset_name), call. = FALSE)
}

# 🔴 Define: IPCW and prediction helpers ===============================
evaluate_survfit <- function(fit, times) {
  times <- as.numeric(times)

  if (length(times) == 0L) {
    return(numeric(0))
  }

  if (is.null(fit) || length(fit$time) == 0L) {
    return(rep(1, length(times)))
  }

  summary_obj <- summary(fit, times = times, extend = TRUE)
  surv_vec <- as.numeric(summary_obj$surv)

  if (length(surv_vec) != length(times)) {
    stop("Unexpected length mismatch while evaluating a survfit object.", call. = FALSE)
  }

  surv_vec
}

fit_event_survfit <- function(time_year, event_main) {
  survival::survfit(survival::Surv(time_year, event_main) ~ 1)
}

fit_censor_survfit <- function(time_year, censor_main) {
  survival::survfit(survival::Surv(time_year, censor_main) ~ 1)
}

build_ipcw_cache <- function(data, horizons_year) {
  censor_fit <- fit_censor_survfit(data$time_year, data$censor_main)
  event_fit <- fit_event_survfit(data$time_year, data$event_main)

  cache <- vector("list", length(horizons_year))
  names(cache) <- paste0("year_", horizons_year)

  for (idx in seq_along(horizons_year)) {
    horizon <- as.numeric(horizons_year[idx])

    g_horizon <- evaluate_survfit(censor_fit, horizon)
    g_left <- evaluate_survfit(censor_fit, pmax(data$time_year - 1e-08, 0))
    observed_km_risk <- 1 - evaluate_survfit(event_fit, horizon)

    case_ind <- data$event_main == 1L & data$time_year <= horizon
    control_ind <- data$time_year > horizon

    w_case <- ifelse(case_ind, 1 / pmax(g_left, time_origin_epsilon_year), 0)
    w_control <- ifelse(control_ind, 1 / pmax(g_horizon, time_origin_epsilon_year), 0)

    y_horizon <- ifelse(case_ind, 1, ifelse(control_ind, 0, NA_real_))
    w_ipcw <- ifelse(case_ind, w_case, ifelse(control_ind, w_control, NA_real_))

    case_count_ipcw <- sum(w_case, na.rm = TRUE)
    nonevent_count_ipcw <- sum(w_control, na.rm = TRUE)
    case_support_flag <- is.finite(case_count_ipcw) && case_count_ipcw > 0
    nonevent_support_flag <- is.finite(nonevent_count_ipcw) && nonevent_count_ipcw > 0
    binary_outcome_support_flag <- case_support_flag && nonevent_support_flag

    observed_ipcw_risk_raw <- if ((case_count_ipcw + nonevent_count_ipcw) > 0) {
      case_count_ipcw / (case_count_ipcw + nonevent_count_ipcw)
    } else {
      NA_real_
    }

    observed_ipcw_risk <- if (binary_outcome_support_flag) observed_ipcw_risk_raw else NA_real_

    cache[[idx]] <- list(
      horizon_year = horizon,
      censor_fit = censor_fit,
      event_fit = event_fit,
      censor_survival_at_horizon = as.numeric(g_horizon),
      observed_km_risk = as.numeric(observed_km_risk),
      observed_ipcw_risk = as.numeric(observed_ipcw_risk),
      observed_ipcw_risk_raw = as.numeric(observed_ipcw_risk_raw),
      case_ind = case_ind,
      control_ind = control_ind,
      y_horizon = y_horizon,
      w_case = w_case,
      w_control = w_control,
      w_ipcw = w_ipcw,
      case_count_ipcw = case_count_ipcw,
      nonevent_count_ipcw = nonevent_count_ipcw,
      known_count = sum(case_ind | control_ind),
      case_support_flag = as.logical(case_support_flag),
      nonevent_support_flag = as.logical(nonevent_support_flag),
      binary_outcome_support_flag = as.logical(binary_outcome_support_flag)
    )
  }

  list(
    censor_fit = censor_fit,
    event_fit = event_fit,
    horizon_cache = cache
  )
}

make_horizon_support_info <- function(ipcw_info, horizon_year, last_event_time_year, last_followup_time_year) {
  known_count_at_horizon <- as.integer(ipcw_info$known_count %||% NA_integer_)
  case_support_flag <- as.integer(isTRUE(ipcw_info$case_support_flag))
  nonevent_support_flag <- as.integer(isTRUE(ipcw_info$nonevent_support_flag))
  binary_outcome_support_flag <- as.integer(isTRUE(ipcw_info$binary_outcome_support_flag))
  observed_horizon_information_flag <- as.integer((ipcw_info$known_count %||% 0) > 0)
  partial_binary_outcome_support_flag <- as.integer(
    observed_horizon_information_flag == 1L &&
      binary_outcome_support_flag == 0L
  )
  beyond_last_event_flag <- as.integer(
    is.finite(last_event_time_year) &&
      horizon_year > last_event_time_year
  )
  beyond_last_followup_flag <- as.integer(
    is.finite(last_followup_time_year) &&
      horizon_year > last_followup_time_year
  )
  descriptive_projection_flag <- as.integer(beyond_last_followup_flag == 1L)
  low_known_count_warning_flag <- as.integer(
    is.finite(known_count_at_horizon) &&
      known_count_at_horizon < late_horizon_known_count_warning_min
  )
  late_horizon_instability_flag <- as.integer(
    binary_outcome_support_flag == 1L &&
      horizon_year >= late_horizon_warning_start_year &&
      is.finite(known_count_at_horizon) &&
      known_count_at_horizon < late_horizon_known_count_warning_min
  )

  horizon_data_status <- dplyr::case_when(
    binary_outcome_support_flag == 1L ~ "binary_outcome_supported",
    observed_horizon_information_flag == 1L &&
      case_support_flag == 1L &&
      nonevent_support_flag == 0L ~ "event_only_partial_support",
    observed_horizon_information_flag == 1L &&
      case_support_flag == 0L &&
      nonevent_support_flag == 1L ~ "nonevent_only_partial_support",
    observed_horizon_information_flag == 0L ~ "no_observed_horizon_information",
    TRUE ~ "partial_support_other"
  )

  descriptive_value_status <- dplyr::case_when(
    beyond_last_followup_flag == 1L ~ "projection_beyond_last_followup",
    binary_outcome_support_flag == 1L ~ "supported_fixed_horizon_estimation",
    observed_horizon_information_flag == 1L ~ "descriptive_only_partial_outcome_support",
    TRUE ~ "model_only_no_observed_horizon_information"
  )

  late_horizon_warning_reason <- dplyr::case_when(
    late_horizon_instability_flag == 1L ~ sprintf(
      "supported_horizon_year_gte_%s_with_known_count_below_%s",
      late_horizon_warning_start_year,
      late_horizon_known_count_warning_min
    ),
    low_known_count_warning_flag == 1L ~ sprintf(
      "known_count_below_%s",
      late_horizon_known_count_warning_min
    ),
    TRUE ~ "none"
  )

  list(
    known_count_at_horizon = known_count_at_horizon,
    case_support_flag = case_support_flag,
    nonevent_support_flag = nonevent_support_flag,
    observed_horizon_information_flag = observed_horizon_information_flag,
    partial_binary_outcome_support_flag = partial_binary_outcome_support_flag,
    binary_outcome_support_flag = binary_outcome_support_flag,
    beyond_last_event_flag = beyond_last_event_flag,
    beyond_last_followup_flag = beyond_last_followup_flag,
    descriptive_projection_flag = descriptive_projection_flag,
    low_known_count_warning_flag = low_known_count_warning_flag,
    late_horizon_instability_flag = late_horizon_instability_flag,
    late_horizon_warning_reason = late_horizon_warning_reason,
    horizon_data_status = horizon_data_status,
    descriptive_value_status = descriptive_value_status
  )
}

extract_flexsurv_matrix <- function(summary_list, target_times) {
  if (!is.list(summary_list)) {
    summary_list <- list(summary_list)
  }

  target_times <- as.numeric(target_times)

  pred_mat <- t(
    vapply(
      summary_list,
      FUN.VALUE = numeric(length(target_times)),
      FUN = function(x) {
        xx <- as.data.frame(x)

        if ("time" %in% names(xx)) {
          time_vec <- as.numeric(xx$time)

          if (nrow(xx) == length(target_times) && all(abs(time_vec - target_times) < 1e-08)) {
            row_index <- seq_along(target_times)
          } else {
            row_index <- match(round(target_times, 8), round(time_vec, 8))
          }

          xx <- xx[row_index, , drop = FALSE]
        }

        if ("est" %in% names(xx)) {
          return(as.numeric(xx$est))
        }

        as.numeric(xx[[ncol(xx)]])
      }
    )
  )

  pred_mat
}

extract_flexsurv_vector <- function(summary_list) {
  if (!is.list(summary_list)) {
    summary_list <- list(summary_list)
  }

  vapply(
    summary_list,
    FUN.VALUE = numeric(1),
    FUN = function(x) {
      xx <- as.data.frame(x)

      if ("est" %in% names(xx)) {
        return(as.numeric(xx$est[1]))
      }

      as.numeric(xx[[ncol(xx)]][1])
    }
  )
}

predict_parametric_survival <- function(fit, newdata, horizons_year) {
  surv_list <- summary(
    fit,
    newdata = newdata,
    t = horizons_year,
    type = "survival",
    ci = FALSE,
    se = FALSE
  )

  median_list <- tryCatch(
    summary(
      fit,
      newdata = newdata,
      type = "quantile",
      quantiles = 0.5,
      ci = FALSE,
      se = FALSE
    ),
    error = function(e) NULL
  )

  surv_mat <- extract_flexsurv_matrix(surv_list, horizons_year)
  median_vec <- if (is.null(median_list)) {
    rep(NA_real_, nrow(newdata))
  } else {
    extract_flexsurv_vector(median_list)
  }

  list(
    survival = surv_mat,
    risk = 1 - surv_mat,
    median_survival = median_vec,
    score = NULL
  )
}

predict_cox_survival <- function(fit, newdata, horizons_year) {
  linear_predictor <- as.numeric(stats::predict(fit, newdata = newdata, type = "lp"))
  baseline_hazard <- survival::basehaz(fit, centered = FALSE)

  if (nrow(baseline_hazard) == 0L) {
    cumulative_baseline <- rep(0, length(horizons_year))
  } else {
    baseline_hazard <- baseline_hazard[order(baseline_hazard$time), , drop = FALSE]
    cumulative_fun <- stats::stepfun(
      baseline_hazard$time,
      c(0, baseline_hazard$hazard),
      right = TRUE
    )
    cumulative_baseline <- pmax(cumulative_fun(horizons_year), 0)
  }

  survival_mat <- exp(-outer(exp(linear_predictor), cumulative_baseline, `*`))

  list(
    survival = survival_mat,
    risk = 1 - survival_mat,
    median_survival = rep(NA_real_, nrow(newdata)),
    score = linear_predictor
  )
}

# 🔴 Define: discrimination, calibration, and threshold scorers ===============================
compute_weighted_auc_cd <- function(score, w_case, w_control) {
  valid <- !is.na(score) & ((w_case > 0) | (w_control > 0))

  if (sum(w_case[valid] > 0) == 0L || sum(w_control[valid] > 0) == 0L) {
    return(NA_real_)
  }

  auc_df <- tibble(
    score = as.numeric(score[valid]),
    w_case = as.numeric(w_case[valid]),
    w_control = as.numeric(w_control[valid])
  ) %>%
    group_by(score) %>%
    summarise(
      w_case = sum(w_case, na.rm = TRUE),
      w_control = sum(w_control, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(score)

  cum_control_before <- c(0, head(cumsum(auc_df$w_control), -1))
  concordant_weight <- sum(auc_df$w_case * cum_control_before, na.rm = TRUE) +
    0.5 * sum(auc_df$w_case * auc_df$w_control, na.rm = TRUE)

  denominator <- sum(auc_df$w_case, na.rm = TRUE) * sum(auc_df$w_control, na.rm = TRUE)

  if (denominator <= 0) {
    return(NA_real_)
  }

  concordant_weight / denominator
}

compute_concordance_metric <- function(time_year, event_main, risk_score, ymax = NULL, timewt = c("n", "n/G2")) {
  timewt <- match.arg(timewt)

  concordance_df <- tibble(
    time_year = as.numeric(time_year),
    event_main = as.integer(event_main),
    risk_score = as.numeric(risk_score)
  ) %>%
    filter(
      !is.na(time_year),
      is.finite(time_year),
      !is.na(event_main),
      event_main %in% c(0L, 1L),
      !is.na(risk_score),
      is.finite(risk_score)
    )

  if (nrow(concordance_df) < 2L) {
    return(list(
      estimate = NA_real_,
      se = NA_real_,
      message = "Insufficient non-missing rows for concordance."
    ))
  }

  if (sum(concordance_df$event_main == 1L, na.rm = TRUE) == 0L) {
    return(list(
      estimate = NA_real_,
      se = NA_real_,
      message = "No observed events available for concordance."
    ))
  }

  if (length(unique(concordance_df$risk_score)) < 2L) {
    return(list(
      estimate = NA_real_,
      se = NA_real_,
      message = "Risk score has fewer than two unique values."
    ))
  }

  result <- tryCatch(
    {
      if (is.null(ymax) || !is.finite(ymax)) {
        survival::concordance(
          survival::Surv(time_year, event_main) ~ risk_score,
          data = concordance_df,
          reverse = TRUE,
          timewt = timewt
        )
      } else {
        survival::concordance(
          survival::Surv(time_year, event_main) ~ risk_score,
          data = concordance_df,
          reverse = TRUE,
          timewt = timewt,
          ymin = 0,
          ymax = as.numeric(ymax)
        )
      }
    },
    error = function(e) e
  )

  if (inherits(result, "error")) {
    return(list(
      estimate = NA_real_,
      se = NA_real_,
      message = conditionMessage(result)
    ))
  }

  estimate_val <- tryCatch(
    as.numeric(result$concordance)[1],
    error = function(e) NA_real_
  )

  var_val <- tryCatch(
    as.numeric(result$var)[1],
    error = function(e) NA_real_
  )

  se_val <- if (is.finite(var_val) && var_val >= 0) sqrt(var_val) else NA_real_

  list(
    estimate = estimate_val,
    se = se_val,
    message = NA_character_
  )
}

compute_harrell_c <- function(time_year, event_main, risk_score, ymax = NULL) {
  compute_concordance_metric(
    time_year = time_year,
    event_main = event_main,
    risk_score = risk_score,
    ymax = ymax,
    timewt = "n"
  )
}

compute_uno_c <- function(time_year, event_main, risk_score, ymax = NULL) {
  compute_concordance_metric(
    time_year = time_year,
    event_main = event_main,
    risk_score = risk_score,
    ymax = ymax,
    timewt = "n/G2"
  )
}

weighted_binomial_nll <- function(eta, y, w) {
  p <- clip_prob(stats::plogis(eta))
  -sum(w * (y * log(p) + (1 - y) * log1p(-p)), na.rm = TRUE)
}

safe_calibration_optimizer <- function(par, objective_fn, lower, upper) {
  par <- as.numeric(par)
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)

  tryCatch(
    stats::optim(
      par = pmin(pmax(par, lower + 1e-06), upper - 1e-06),
      fn = objective_fn,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = 1000)
    ),
    error = function(e) e
  )
}

fit_calibration_offset_model <- function(y, lp, w) {
  prev <- weighted_mean_safe(y, w)

  if (!is.finite(prev) || prev <= 0 || prev >= 1) {
    return(list(
      estimate = NA_real_,
      converged = FALSE,
      boundary = FALSE,
      message = "Weighted event prevalence is outside (0, 1)."
    ))
  }

  start_alpha <- qlogis(prev) - weighted_mean_safe(lp, w)
  lower <- -25
  upper <- 25

  opt <- safe_calibration_optimizer(
    par = start_alpha,
    objective_fn = function(par) weighted_binomial_nll(par[1] + lp, y, w),
    lower = lower,
    upper = upper
  )

  if (inherits(opt, "error")) {
    return(list(
      estimate = NA_real_,
      converged = FALSE,
      boundary = FALSE,
      message = conditionMessage(opt)
    ))
  }

  est <- as.numeric(opt$par[1])
  boundary_flag <- is.finite(est) && (est <= lower + 1e-04 || est >= upper - 1e-04)

  list(
    estimate = est,
    converged = isTRUE(opt$convergence == 0L),
    boundary = boundary_flag,
    message = safe_char(opt$message)
  )
}

fit_calibration_slope_model <- function(y, lp, w) {
  prev <- weighted_mean_safe(y, w)
  lp_sd <- stats::sd(lp)

  if (!is.finite(prev) || prev <= 0 || prev >= 1) {
    return(list(
      intercept = NA_real_,
      slope = NA_real_,
      converged = FALSE,
      boundary = FALSE,
      message = "Weighted event prevalence is outside (0, 1)."
    ))
  }

  if (!is.finite(lp_sd) || lp_sd < 1e-10) {
    return(list(
      intercept = NA_real_,
      slope = NA_real_,
      converged = FALSE,
      boundary = FALSE,
      message = "Logit prediction has near-zero variance."
    ))
  }

  start_par <- c(qlogis(prev), 1)
  lower <- c(-25, -10)
  upper <- c(25, 10)

  opt <- safe_calibration_optimizer(
    par = start_par,
    objective_fn = function(par) weighted_binomial_nll(par[1] + par[2] * lp, y, w),
    lower = lower,
    upper = upper
  )

  if (inherits(opt, "error")) {
    return(list(
      intercept = NA_real_,
      slope = NA_real_,
      converged = FALSE,
      boundary = FALSE,
      message = conditionMessage(opt)
    ))
  }

  est <- as.numeric(opt$par)
  boundary_flag <- any(est <= lower + 1e-04 | est >= upper - 1e-04)

  list(
    intercept = est[1],
    slope = est[2],
    converged = isTRUE(opt$convergence == 0L),
    boundary = boundary_flag,
    message = safe_char(opt$message)
  )
}

finalize_calibration_statistics <- function(out) {
  calibration_fit_success_flag <- as.integer(
    out$calibration_support_flag == 1L &&
      out$calibration_offset_converged_flag == 1L &&
      out$calibration_slope_converged_flag == 1L &&
      out$calibration_offset_boundary_flag == 0L &&
      out$calibration_slope_boundary_flag == 0L &&
      is.finite(out$calibration_intercept_offset) &&
      is.finite(out$calibration_intercept_free) &&
      is.finite(out$calibration_slope)
  )

  out$calibration_regression_suppressed_flag <- as.integer(
    out$calibration_support_flag == 1L &&
      calibration_fit_success_flag == 0L
  )

  out$calibration_fit_success_flag <- calibration_fit_success_flag

  if (calibration_fit_success_flag != 1L) {
    out$calibration_intercept_offset <- NA_real_
    out$calibration_intercept_free <- NA_real_
    out$calibration_slope <- NA_real_
  }

  out
}

compute_calibration_statistics <- function(pred_risk, ipcw_info, observed_reference_risk) {
  pred_risk <- clip_prob(pred_risk)

  known_df <- tibble(
    y = ipcw_info$y_horizon,
    lp = qlogis(pred_risk),
    w = ipcw_info$w_ipcw
  ) %>%
    filter(
      !is.na(y),
      !is.na(lp),
      is.finite(lp),
      !is.na(w),
      is.finite(w),
      w > 0
    )

  calibration_case_weight <- sum(known_df$w * known_df$y, na.rm = TRUE)
  calibration_nonevent_weight <- sum(known_df$w * (1 - known_df$y), na.rm = TRUE)
  binary_outcome_support_flag <- as.integer(
    is.finite(calibration_case_weight) &&
      is.finite(calibration_nonevent_weight) &&
      calibration_case_weight > 0 &&
      calibration_nonevent_weight > 0
  )
  calibration_minimum_n_flag <- as.integer(nrow(known_df) >= 10)
  calibration_support_flag <- as.integer(
    binary_outcome_support_flag == 1L &&
      calibration_minimum_n_flag == 1L
  )

  out <- list(
    mean_predicted_risk = mean(pred_risk, na.rm = TRUE),
    mean_calibration_difference = observed_reference_risk - mean(pred_risk, na.rm = TRUE),
    calibration_intercept_offset = NA_real_,
    calibration_intercept_free = NA_real_,
    calibration_slope = NA_real_,
    calibration_known_n = nrow(known_df),
    calibration_case_weight = calibration_case_weight,
    calibration_nonevent_weight = calibration_nonevent_weight,
    binary_outcome_support_flag = binary_outcome_support_flag,
    calibration_minimum_n_flag = calibration_minimum_n_flag,
    calibration_support_flag = calibration_support_flag,
    calibration_fit_success_flag = 0L,
    calibration_offset_converged_flag = 0L,
    calibration_slope_converged_flag = 0L,
    calibration_offset_boundary_flag = 0L,
    calibration_slope_boundary_flag = 0L,
    calibration_regression_suppressed_flag = 0L,
    calibration_offset_message = NA_character_,
    calibration_slope_message = NA_character_
  )

  if (calibration_support_flag != 1L) {
    return(out)
  }

  offset_fit <- fit_calibration_offset_model(
    y = known_df$y,
    lp = known_df$lp,
    w = known_df$w
  )

  slope_fit <- fit_calibration_slope_model(
    y = known_df$y,
    lp = known_df$lp,
    w = known_df$w
  )

  out$calibration_intercept_offset <- offset_fit$estimate
  out$calibration_intercept_free <- slope_fit$intercept
  out$calibration_slope <- slope_fit$slope
  out$calibration_offset_converged_flag <- as.integer(isTRUE(offset_fit$converged))
  out$calibration_slope_converged_flag <- as.integer(isTRUE(slope_fit$converged))
  out$calibration_offset_boundary_flag <- as.integer(isTRUE(offset_fit$boundary))
  out$calibration_slope_boundary_flag <- as.integer(isTRUE(slope_fit$boundary))
  out$calibration_offset_message <- safe_char(offset_fit$message)
  out$calibration_slope_message <- safe_char(slope_fit$message)

  finalize_calibration_statistics(out)
}

derive_calibration_failure_reason <- function(calibration_stats, horizon_support_info = NULL) {
  support_info <- horizon_support_info %||% list()

  if (isTRUE(calibration_stats$calibration_fit_success_flag == 1L)) {
    return("none")
  }

  reasons <- character(0)

  if (isTRUE(support_info$binary_outcome_support_flag == 0L)) {
    reasons <- c(reasons, "binary_outcome_not_supported")
  }

  if (isTRUE(support_info$beyond_last_followup_flag == 1L)) {
    reasons <- c(reasons, "projection_beyond_last_followup")
  }

  if (isTRUE(calibration_stats$calibration_minimum_n_flag == 0L)) {
    reasons <- c(reasons, "known_n_below_minimum")
  }

  if (isTRUE(calibration_stats$calibration_offset_converged_flag == 0L)) {
    reasons <- c(reasons, "offset_model_not_converged")
  }

  if (isTRUE(calibration_stats$calibration_slope_converged_flag == 0L)) {
    reasons <- c(reasons, "slope_model_not_converged")
  }

  if (isTRUE(calibration_stats$calibration_offset_boundary_flag == 1L)) {
    reasons <- c(reasons, "offset_boundary")
  }

  if (isTRUE(calibration_stats$calibration_slope_boundary_flag == 1L)) {
    reasons <- c(reasons, "slope_boundary")
  }

  if (length(reasons) == 0L && isTRUE(calibration_stats$calibration_support_flag == 0L)) {
    reasons <- c(reasons, "calibration_support_not_met")
  }

  if (length(reasons) == 0L && isTRUE(calibration_stats$calibration_fit_success_flag == 0L)) {
    reasons <- c(reasons, "unspecified_failure")
  }

  collapse_nonempty_chr(reasons)
}

compute_calibration_bins <- function(pred_risk, ipcw_info, group_count = calibration_group_count) {
  pred_risk <- clip_prob(pred_risk)
  unique_pred_n <- dplyr::n_distinct(pred_risk)

  overall_binary_support_flag <- if (!is.null(ipcw_info$binary_outcome_support_flag)) {
    isTRUE(ipcw_info$binary_outcome_support_flag)
  } else {
    sum(ipcw_info$w_case, na.rm = TRUE) > 0 &&
      sum(ipcw_info$w_control, na.rm = TRUE) > 0
  }

  if (unique_pred_n <= 1L) {
    bin_factor <- factor(rep("All", length(pred_risk)))
  } else {
    probs <- seq(0, 1, length.out = min(group_count, unique_pred_n) + 1)
    breaks_vec <- unique(stats::quantile(pred_risk, probs = probs, na.rm = TRUE, type = 2))

    if (length(breaks_vec) <= 2L) {
      bin_factor <- factor(rep("All", length(pred_risk)))
    } else {
      bin_factor <- cut(
        pred_risk,
        breaks = breaks_vec,
        include.lowest = TRUE,
        dig.lab = 12
      )
    }
  }

  tibble(
    bin_label = bin_factor,
    predicted_risk = pred_risk,
    w_case = ipcw_info$w_case,
    w_control = ipcw_info$w_control
  ) %>%
    group_by(bin_label) %>%
    summarise(
      n_subject = n(),
      mean_predicted_risk = mean(predicted_risk, na.rm = TRUE),
      case_count_ipcw = sum(w_case, na.rm = TRUE),
      nonevent_count_ipcw = sum(w_control, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      bin_index = row_number(),
      case_support_flag = as.integer(case_count_ipcw > 0),
      nonevent_support_flag = as.integer(nonevent_count_ipcw > 0),
      bin_binary_support_flag = as.integer(case_count_ipcw > 0 & nonevent_count_ipcw > 0),
      binary_outcome_support_flag = as.integer(overall_binary_support_flag),
      observed_risk_ipcw = dplyr::if_else(
        bin_binary_support_flag == 1L,
        safe_divide(case_count_ipcw, case_count_ipcw + nonevent_count_ipcw),
        NA_real_
      )
    )
}

compute_threshold_vector <- function(pred_positive, ipcw_info, threshold, n_total) {
  pred_positive <- as.logical(pred_positive)
  binary_support_flag <- isTRUE(ipcw_info$binary_outcome_support_flag)

  positive_classification_rate <- mean(pred_positive)

  base_out <- c(
    positive_classification_rate = positive_classification_rate,
    true_positive_count_ipcw = NA_real_,
    false_positive_count_ipcw = NA_real_,
    false_negative_count_ipcw = NA_real_,
    true_negative_count_ipcw = NA_real_,
    false_positive_burden_all = NA_real_,
    false_positive_burden_nonevents = NA_real_,
    false_positive_per_100 = NA_real_,
    sensitivity = NA_real_,
    specificity = NA_real_,
    ppv = NA_real_,
    npv = NA_real_,
    false_discovery_proportion = NA_real_,
    event_prevalence_ipcw = NA_real_,
    event_rate_per_subject_ipcw = NA_real_,
    net_benefit = NA_real_,
    net_benefit_treat_all = NA_real_,
    net_reduction_unnecessary_per_100 = NA_real_,
    case_count_ipcw = ipcw_info$case_count_ipcw,
    nonevent_count_ipcw = ipcw_info$nonevent_count_ipcw,
    known_count = ipcw_info$known_count,
    binary_outcome_support_flag = as.integer(binary_support_flag)
  )

  if (!binary_support_flag) {
    return(base_out)
  }

  tp <- sum(ipcw_info$w_case[pred_positive], na.rm = TRUE)
  fp <- sum(ipcw_info$w_control[pred_positive], na.rm = TRUE)
  fn <- sum(ipcw_info$w_case[!pred_positive], na.rm = TRUE)
  tn <- sum(ipcw_info$w_control[!pred_positive], na.rm = TRUE)

  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  ppv <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  npv <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
  fdp <- if ((tp + fp) > 0) fp / (tp + fp) else NA_real_

  total_case <- sum(ipcw_info$w_case, na.rm = TRUE)
  total_nonevent <- sum(ipcw_info$w_control, na.rm = TRUE)
  prevalence_ipcw <- if ((total_case + total_nonevent) > 0) {
    total_case / (total_case + total_nonevent)
  } else {
    NA_real_
  }

  threshold_odds <- threshold / (1 - threshold)
  event_rate_per_subject_ipcw <- total_case / n_total
  net_benefit <- (tp / n_total) - (fp / n_total) * threshold_odds
  net_benefit_all <- (total_case / n_total) - (total_nonevent / n_total) * threshold_odds
  net_reduction_unnecessary_per_100 <- if (threshold_odds > 0) {
    100 * (net_benefit - net_benefit_all) / threshold_odds
  } else {
    NA_real_
  }

  c(
    positive_classification_rate = positive_classification_rate,
    true_positive_count_ipcw = tp,
    false_positive_count_ipcw = fp,
    false_negative_count_ipcw = fn,
    true_negative_count_ipcw = tn,
    false_positive_burden_all = fp / n_total,
    false_positive_burden_nonevents = if ((tn + fp) > 0) fp / (tn + fp) else NA_real_,
    false_positive_per_100 = 100 * fp / n_total,
    sensitivity = sensitivity,
    specificity = specificity,
    ppv = ppv,
    npv = npv,
    false_discovery_proportion = fdp,
    event_prevalence_ipcw = prevalence_ipcw,
    event_rate_per_subject_ipcw = event_rate_per_subject_ipcw,
    net_benefit = net_benefit,
    net_benefit_treat_all = net_benefit_all,
    net_reduction_unnecessary_per_100 = net_reduction_unnecessary_per_100,
    case_count_ipcw = ipcw_info$case_count_ipcw,
    nonevent_count_ipcw = ipcw_info$nonevent_count_ipcw,
    known_count = ipcw_info$known_count,
    binary_outcome_support_flag = as.integer(binary_support_flag)
  )
}

# 🔴 Define: output row constructors and finalizers ===============================
make_metric_rows <- function(meta_row, horizon_info, threshold_value, metric_domain, metrics_named, horizon_support_info = NULL) {
  support_info <- horizon_support_info %||% list()

  tibble(
    stage = stage_name,
    stage_role = stage_role,
    stage_spec_version = stage_spec_version,
    validation_level = meta_row$validation_level %||% validation_level_default,
    optimism_correction_method = meta_row$optimism_correction_method %||% optimism_correction_method_default,
    dataset = meta_row$dataset %||% NA_character_,
    dataset_key = meta_row$dataset_key %||% meta_row$dataset %||% NA_character_,
    risk_scale = meta_row$risk_scale %||% "transition_only_main",
    model_class = meta_row$model_class %||% NA_character_,
    model_family = meta_row$model_family %||% NA_character_,
    model_id = meta_row$model_id %||% NA_character_,
    strategy_family = meta_row$strategy_family %||% NA_character_,
    formula_id = meta_row$formula_id %||% NA_character_,
    formula_name = meta_row$formula_name %||% NA_character_,
    formula_label = meta_row$formula_label %||% NA_character_,
    formula_scope = meta_row$formula_scope %||% NA_character_,
    site_branch = meta_row$site_branch %||% NA_character_,
    interaction_branch = meta_row$interaction_branch %||% NA_character_,
    anchor_horizon_year = meta_row$anchor_horizon_year %||% NA_integer_,
    tail_extrapolation_capable_flag = meta_row$tail_extrapolation_capable_flag %||% NA_integer_,
    extrapolation_role = meta_row$extrapolation_role %||% NA_character_,
    horizon_id = horizon_info$horizon_id %||% NA_character_,
    horizon_year = horizon_info$horizon_year %||% NA_integer_,
    horizon_label = horizon_info$horizon_label %||% NA_character_,
    support_tier = horizon_info$support_tier %||% NA_character_,
    support_tier_standard = horizon_info$support_tier_standard %||% horizon_info$support_tier %||% NA_character_,
    interpretation_tier = horizon_info$interpretation_tier %||% NA_character_,
    primary_supported_flag = horizon_info$primary_supported_flag %||% NA,
    horizon_evidence_class = horizon_info$horizon_evidence_class %||% NA_character_,
    claim_restriction_flag = horizon_info$claim_restriction_flag %||% NA_character_,
    interpretation_note = horizon_info$interpretation_note %||% NA_character_,
    observed_case_support_flag = as.integer(support_info$case_support_flag %||% NA_integer_),
    observed_nonevent_support_flag = as.integer(support_info$nonevent_support_flag %||% NA_integer_),
    observed_horizon_information_flag = as.integer(support_info$observed_horizon_information_flag %||% NA_integer_),
    partial_binary_outcome_support_flag = as.integer(support_info$partial_binary_outcome_support_flag %||% NA_integer_),
    beyond_last_event_flag = as.integer(support_info$beyond_last_event_flag %||% NA_integer_),
    beyond_last_followup_flag = as.integer(support_info$beyond_last_followup_flag %||% NA_integer_),
    descriptive_projection_flag = as.integer(support_info$descriptive_projection_flag %||% NA_integer_),
    known_count_at_horizon = as.integer(support_info$known_count_at_horizon %||% NA_integer_),
    low_known_count_warning_flag = as.integer(support_info$low_known_count_warning_flag %||% NA_integer_),
    late_horizon_instability_flag = as.integer(support_info$late_horizon_instability_flag %||% NA_integer_),
    late_horizon_warning_reason = safe_char(support_info$late_horizon_warning_reason),
    horizon_data_status = safe_char(support_info$horizon_data_status),
    descriptive_value_status = safe_char(support_info$descriptive_value_status),
    threshold = threshold_value,
    metric_domain = metric_domain,
    metric_name = names(metrics_named),
    metric_value = as.numeric(unname(metrics_named))
  )
}

make_calibration_diagnostic_row <- function(meta_row, horizon_info, calibration_stats, horizon_support_info = NULL) {
  support_info <- horizon_support_info %||% list()

  tibble(
    stage = stage_name,
    stage_role = stage_role,
    stage_spec_version = stage_spec_version,
    validation_level = meta_row$validation_level %||% validation_level_default,
    optimism_correction_method = meta_row$optimism_correction_method %||% optimism_correction_method_default,
    dataset = meta_row$dataset %||% NA_character_,
    dataset_key = meta_row$dataset_key %||% meta_row$dataset %||% NA_character_,
    risk_scale = meta_row$risk_scale %||% "transition_only_main",
    model_class = meta_row$model_class %||% NA_character_,
    model_family = meta_row$model_family %||% NA_character_,
    model_id = meta_row$model_id %||% NA_character_,
    strategy_family = meta_row$strategy_family %||% NA_character_,
    formula_id = meta_row$formula_id %||% NA_character_,
    formula_name = meta_row$formula_name %||% NA_character_,
    formula_label = meta_row$formula_label %||% NA_character_,
    formula_scope = meta_row$formula_scope %||% NA_character_,
    site_branch = meta_row$site_branch %||% NA_character_,
    interaction_branch = meta_row$interaction_branch %||% NA_character_,
    anchor_horizon_year = meta_row$anchor_horizon_year %||% NA_integer_,
    tail_extrapolation_capable_flag = meta_row$tail_extrapolation_capable_flag %||% NA_integer_,
    extrapolation_role = meta_row$extrapolation_role %||% NA_character_,
    horizon_id = horizon_info$horizon_id %||% NA_character_,
    horizon_year = horizon_info$horizon_year %||% NA_integer_,
    horizon_label = horizon_info$horizon_label %||% NA_character_,
    support_tier = horizon_info$support_tier %||% NA_character_,
    support_tier_standard = horizon_info$support_tier_standard %||% horizon_info$support_tier %||% NA_character_,
    interpretation_tier = horizon_info$interpretation_tier %||% NA_character_,
    primary_supported_flag = horizon_info$primary_supported_flag %||% NA,
    horizon_evidence_class = horizon_info$horizon_evidence_class %||% NA_character_,
    claim_restriction_flag = horizon_info$claim_restriction_flag %||% NA_character_,
    interpretation_note = horizon_info$interpretation_note %||% NA_character_,
    calibration_target = calibration_target_label,
    mean_calibration_reference = mean_calibration_reference_label,
    mean_predicted_risk = safe_num(calibration_stats$mean_predicted_risk),
    mean_calibration_difference = safe_num(calibration_stats$mean_calibration_difference),
    calibration_intercept_offset = safe_num(calibration_stats$calibration_intercept_offset),
    calibration_intercept_free = safe_num(calibration_stats$calibration_intercept_free),
    calibration_slope = safe_num(calibration_stats$calibration_slope),
    calibration_known_n = as.integer(calibration_stats$calibration_known_n %||% NA_integer_),
    calibration_case_weight = safe_num(calibration_stats$calibration_case_weight),
    calibration_nonevent_weight = safe_num(calibration_stats$calibration_nonevent_weight),
    binary_outcome_support_flag = as.integer(calibration_stats$binary_outcome_support_flag %||% NA_integer_),
    calibration_minimum_n_flag = as.integer(calibration_stats$calibration_minimum_n_flag %||% NA_integer_),
    calibration_support_flag = as.integer(calibration_stats$calibration_support_flag %||% NA_integer_),
    calibration_fit_success_flag = as.integer(calibration_stats$calibration_fit_success_flag %||% NA_integer_),
    calibration_offset_converged_flag = as.integer(calibration_stats$calibration_offset_converged_flag %||% NA_integer_),
    calibration_slope_converged_flag = as.integer(calibration_stats$calibration_slope_converged_flag %||% NA_integer_),
    calibration_offset_boundary_flag = as.integer(calibration_stats$calibration_offset_boundary_flag %||% NA_integer_),
    calibration_slope_boundary_flag = as.integer(calibration_stats$calibration_slope_boundary_flag %||% NA_integer_),
    calibration_regression_suppressed_flag = as.integer(calibration_stats$calibration_regression_suppressed_flag %||% NA_integer_),
    observed_case_support_flag = as.integer(support_info$case_support_flag %||% NA_integer_),
    observed_nonevent_support_flag = as.integer(support_info$nonevent_support_flag %||% NA_integer_),
    observed_horizon_information_flag = as.integer(support_info$observed_horizon_information_flag %||% NA_integer_),
    partial_binary_outcome_support_flag = as.integer(support_info$partial_binary_outcome_support_flag %||% NA_integer_),
    beyond_last_event_flag = as.integer(support_info$beyond_last_event_flag %||% NA_integer_),
    beyond_last_followup_flag = as.integer(support_info$beyond_last_followup_flag %||% NA_integer_),
    descriptive_projection_flag = as.integer(support_info$descriptive_projection_flag %||% NA_integer_),
    known_count_at_horizon = as.integer(support_info$known_count_at_horizon %||% NA_integer_),
    low_known_count_warning_flag = as.integer(support_info$low_known_count_warning_flag %||% NA_integer_),
    late_horizon_instability_flag = as.integer(support_info$late_horizon_instability_flag %||% NA_integer_),
    late_horizon_warning_reason = safe_char(support_info$late_horizon_warning_reason),
    horizon_data_status = safe_char(support_info$horizon_data_status),
    descriptive_value_status = safe_char(support_info$descriptive_value_status),
    calibration_failure_reason = derive_calibration_failure_reason(calibration_stats, support_info),
    calibration_offset_message = safe_char(calibration_stats$calibration_offset_message),
    calibration_slope_message = safe_char(calibration_stats$calibration_slope_message)
  )
}


ensure_default_columns <- function(df, defaults) {
  out <- tibble::as_tibble(df)
  for (nm in names(defaults)) {
    if (!nm %in% names(out)) {
      value <- defaults[[nm]]
      if (length(value) %in% c(0L, 1L)) {
        out[[nm]] <- rep(value, nrow(out))
      } else if (length(value) == nrow(out)) {
        out[[nm]] <- value
      } else {
        stop(sprintf("Default column `%s` has incompatible length.", nm), call. = FALSE)
      }
    }
  }
  out
}

finalize_model_registry <- function(model_registry, formula_registry_core, horizon_registry) {
  out <- model_registry %>%
    ensure_default_columns(
      list(
        validation_level = validation_level_default,
        optimism_correction_method = optimism_correction_method_default,
        analysis_time_variable = "days_followup",
        reporting_time_variable = "time_year",
        event_definition = "status_num == 1",
        censoring_definition = "status_num %in% c(0, 2)",
        risk_scale = "transition_only_main",
        model_class = "non_cure",
        effective_parameter_count = NA_integer_,
        n_parameter = NA_integer_,
        likelihood_basis = NA_character_,
        information_criterion_comparable_across_families_flag = NA_integer_,
        tail_extrapolation_capable_flag = NA_integer_,
        extrapolation_role = NA_character_
      )
    ) %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      stage = stage_name,
      stage_role = stage_role,
      stage_spec_version = stage_spec_version,
      validation_level = as.character(validation_level),
      optimism_correction_method = as.character(optimism_correction_method),
      analysis_time_variable = as.character(analysis_time_variable),
      reporting_time_variable = as.character(reporting_time_variable),
      event_definition = as.character(event_definition),
      censoring_definition = as.character(censoring_definition),
      risk_scale = as.character(risk_scale),
      model_class = as.character(model_class),
      effective_parameter_count = as.integer(dplyr::coalesce(effective_parameter_count, n_parameter)),
      fit_engine = dplyr::case_when(
        model_family == "coxph" ~ "survival::coxph",
        model_family %in% c("exponential", "weibull", "lognormal", "loglogistic") ~ "flexsurv::flexsurvreg",
        TRUE ~ NA_character_
      ),
      likelihood_basis = dplyr::coalesce(likelihood_basis, make_likelihood_basis(model_family)),
      information_criterion_comparable_across_families_flag = dplyr::coalesce(
        information_criterion_comparable_across_families_flag,
        make_information_criterion_comparable_flag(model_family)
      ),
      tail_extrapolation_capable_flag = dplyr::coalesce(
        tail_extrapolation_capable_flag,
        make_tail_extrapolation_capable_flag(model_family)
      ),
      extrapolation_role = dplyr::coalesce(
        extrapolation_role,
        make_extrapolation_role(model_family)
      )
    ) %>%
    join_formula_metadata(formula_registry_core) %>%
    augment_anchor_horizon_metadata(horizon_registry) %>%
    collapse_dataset_key_suffixes() %>%
    arrange(match(dataset, dataset_order), formula_name, model_family)

  out
}

finalize_subject_horizon_risk_long <- function(subject_horizon_risk_long, formula_registry_core, horizon_registry) {
  subject_horizon_risk_long %>%
    ensure_default_columns(
      list(
        validation_level = validation_level_default,
        optimism_correction_method = optimism_correction_method_default,
        risk_scale = "transition_only_main"
      )
    ) %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      stage = stage_name,
      stage_role = stage_role,
      stage_spec_version = stage_spec_version,
      validation_level = as.character(validation_level),
      optimism_correction_method = as.character(optimism_correction_method),
      risk_scale = as.character(risk_scale)
    ) %>%
    join_formula_metadata(formula_registry_core) %>%
    join_horizon_metadata(horizon_registry) %>%
    collapse_dataset_key_suffixes() %>%
    arrange(match(dataset, dataset_order), formula_name, model_family, unique_person_id, horizon_year)
}

finalize_model_performance_long <- function(model_performance_long, formula_registry_core, horizon_registry) {
  model_performance_long %>%
    ensure_default_columns(
      list(
        validation_level = validation_level_default,
        optimism_correction_method = optimism_correction_method_default,
        risk_scale = "transition_only_main"
      )
    ) %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      stage = stage_name,
      stage_role = stage_role,
      stage_spec_version = stage_spec_version,
      validation_level = as.character(validation_level),
      optimism_correction_method = as.character(optimism_correction_method),
      risk_scale = as.character(risk_scale)
    ) %>%
    join_formula_metadata(formula_registry_core) %>%
    join_horizon_metadata(horizon_registry) %>%
    collapse_dataset_key_suffixes() %>%
    arrange(match(dataset, dataset_order), formula_name, model_family, horizon_year, threshold, metric_domain, metric_name)
}

finalize_calibration_bins_long <- function(calibration_bins_long, formula_registry_core, horizon_registry) {
  calibration_bins_long %>%
    ensure_default_columns(
      list(
        validation_level = validation_level_default,
        optimism_correction_method = optimism_correction_method_default,
        risk_scale = "transition_only_main"
      )
    ) %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      stage = stage_name,
      stage_role = stage_role,
      stage_spec_version = stage_spec_version,
      validation_level = as.character(validation_level),
      optimism_correction_method = as.character(optimism_correction_method),
      risk_scale = as.character(risk_scale)
    ) %>%
    join_formula_metadata(formula_registry_core) %>%
    join_horizon_metadata(horizon_registry) %>%
    collapse_dataset_key_suffixes() %>%
    arrange(match(dataset, dataset_order), formula_name, model_family, horizon_year, bin_index)
}

finalize_calibration_diagnostics_long <- function(calibration_diagnostics_long, formula_registry_core, horizon_registry) {
  calibration_diagnostics_long %>%
    ensure_default_columns(
      list(
        validation_level = validation_level_default,
        optimism_correction_method = optimism_correction_method_default,
        risk_scale = "transition_only_main"
      )
    ) %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      stage = stage_name,
      stage_role = stage_role,
      stage_spec_version = stage_spec_version,
      validation_level = as.character(validation_level),
      optimism_correction_method = as.character(optimism_correction_method),
      risk_scale = as.character(risk_scale)
    ) %>%
    join_formula_metadata(formula_registry_core) %>%
    join_horizon_metadata(horizon_registry) %>%
    collapse_dataset_key_suffixes() %>%
    arrange(match(dataset, dataset_order), formula_name, model_family, horizon_year)
}

# 🔴 Define: model fitting helpers ===============================
fit_single_noncure_model <- function(data, formula_rhs, model_family, flexsurv_dist = NA_character_) {
  surv_formula <- stats::as.formula(
    paste0("survival::Surv(time_year_model, event_main) ~ ", strip_formula_rhs(formula_rhs))
  )

  if (model_family == "coxph") {
    fit_object <- tryCatch(
      survival::coxph(
        formula = surv_formula,
        data = data,
        ties = "breslow",
        x = TRUE,
        y = TRUE,
        model = TRUE,
        singular.ok = TRUE
      ),
      error = function(e) e
    )

    if (inherits(fit_object, "error")) {
      return(list(
        fit = NULL,
        converged = FALSE,
        convergence_code = NA_integer_,
        error_message = conditionMessage(fit_object)
      ))
    }

    return(list(
      fit = fit_object,
      converged = is.null(fit_object$fail),
      convergence_code = if (is.null(fit_object$fail)) 0L else 1L,
      error_message = safe_char(fit_object$fail)
    ))
  }

  fit_object <- tryCatch(
    flexsurv::flexsurvreg(
      formula = surv_formula,
      data = data,
      dist = flexsurv_dist
    ),
    error = function(e) e
  )

  if (inherits(fit_object, "error")) {
    return(list(
      fit = NULL,
      converged = FALSE,
      convergence_code = NA_integer_,
      error_message = conditionMessage(fit_object)
    ))
  }

  conv_code <- fit_object$opt$convergence %||% fit_object$optim$convergence %||% NA_integer_
  conv_flag <- if (is.na(conv_code)) TRUE else conv_code == 0L

  list(
    fit = fit_object,
    converged = conv_flag,
    convergence_code = conv_code,
    error_message = NA_character_
  )
}

extract_fit_statistics <- function(fit_object, data, model_family, convergence_code, converged, error_message) {
  last_event_time_year <- if (sum(data$event_main, na.rm = TRUE) > 0) {
    max(data$time_year[data$event_main == 1], na.rm = TRUE)
  } else {
    NA_real_
  }

  last_followup_time_year <- max(data$time_year, na.rm = TRUE)
  likelihood_basis <- make_likelihood_basis(model_family)
  information_criterion_comparable_across_families_flag <- make_information_criterion_comparable_flag(model_family)
  information_criterion_note <- dplyr::case_when(
    model_family == "coxph" ~ "cox_partial_likelihood_not_directly_comparable_to_parametric_full_likelihood",
    model_family %in% c("exponential", "weibull", "lognormal", "loglogistic") ~ "full_likelihood_information_criteria_comparable_across_parametric_families",
    TRUE ~ "not_applicable"
  )

  if (is.null(fit_object)) {
    return(tibble(
      converged = FALSE,
      convergence_code = convergence_code,
      error_message = error_message,
      n_obs = nrow(data),
      n_event = sum(data$event_main, na.rm = TRUE),
      n_parameter = NA_integer_,
      effective_parameter_count = NA_integer_,
      likelihood_basis = likelihood_basis,
      information_criterion_comparable_across_families_flag = information_criterion_comparable_across_families_flag,
      information_criterion_note = information_criterion_note,
      logLik = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      partial_logLik = NA_real_,
      partial_AIC = NA_real_,
      partial_BIC = NA_real_,
      last_event_time_year = last_event_time_year,
      last_followup_time_year = last_followup_time_year,
      extrapolation_rule = ifelse(
        model_family == "coxph",
        "baseline_carried_forward_after_last_event",
        "parametric_extrapolation"
      )
    ))
  }

  n_parameter <- length(stats::coef(fit_object))
  raw_loglik_value <- tryCatch(as.numeric(stats::logLik(fit_object))[1], error = function(e) NA_real_)
  raw_aic_value <- tryCatch(as.numeric(stats::AIC(fit_object))[1], error = function(e) NA_real_)
  raw_bic_value <- if (!is.na(raw_loglik_value)) {
    -2 * raw_loglik_value + log(nrow(data)) * n_parameter
  } else {
    NA_real_
  }

  full_loglik_value <- if (model_family == "coxph") NA_real_ else raw_loglik_value
  full_aic_value <- if (model_family == "coxph") NA_real_ else raw_aic_value
  full_bic_value <- if (model_family == "coxph") NA_real_ else raw_bic_value
  partial_loglik_value <- if (model_family == "coxph") raw_loglik_value else NA_real_
  partial_aic_value <- if (model_family == "coxph") raw_aic_value else NA_real_
  partial_bic_value <- if (model_family == "coxph") raw_bic_value else NA_real_

  tibble(
    converged = as.logical(converged),
    convergence_code = as.integer(convergence_code),
    error_message = safe_char(error_message),
    n_obs = nrow(data),
    n_event = sum(data$event_main, na.rm = TRUE),
    n_parameter = as.integer(n_parameter),
    effective_parameter_count = as.integer(n_parameter),
    likelihood_basis = likelihood_basis,
    information_criterion_comparable_across_families_flag = information_criterion_comparable_across_families_flag,
    information_criterion_note = information_criterion_note,
    logLik = full_loglik_value,
    AIC = full_aic_value,
    BIC = full_bic_value,
    partial_logLik = partial_loglik_value,
    partial_AIC = partial_aic_value,
    partial_BIC = partial_bic_value,
    last_event_time_year = last_event_time_year,
    last_followup_time_year = last_followup_time_year,
    extrapolation_rule = ifelse(
      model_family == "coxph",
      "baseline_carried_forward_after_last_event",
      "parametric_extrapolation"
    )
  )
}

# 🔴 Define: cache readers and reuse validators ===============================
load_existing_stage5_core <- function() {
  list(
    model_registry = safe_read_csv_file(model_registry_file),
    subject_horizon_risk_long = safe_read_csv_file(risk_long_file),
    model_performance_long = safe_read_csv_file(performance_long_file),
    calibration_bins_long = safe_read_csv_file(calibration_bins_file),
    calibration_diagnostics_long = safe_read_csv_file(calibration_diagnostics_file),
    fitted_models = safe_read_rds_file(fitted_models_file)
  )
}

validate_existing_model_registry <- function(model_registry, model_plan) {
  if (is.null(model_registry)) {
    return(FALSE)
  }

  required_cols <- c("dataset", "model_id", "model_family", "formula_id", "converged")
  if (!all(required_cols %in% names(model_registry))) {
    return(FALSE)
  }

  model_registry_key <- model_registry %>%
    transmute(
      dataset = normalize_dataset_label(dataset),
      model_id = as.character(model_id),
      model_family = as.character(model_family),
      formula_id = as.character(formula_id)
    ) %>%
    distinct()

  plan_key <- model_plan %>%
    transmute(
      dataset = normalize_dataset_label(dataset),
      model_id = as.character(model_id),
      model_family = as.character(model_family),
      formula_id = as.character(formula_id)
    ) %>%
    distinct()

  nrow(model_registry_key) == nrow(plan_key) &&
    nrow(dplyr::anti_join(plan_key, model_registry_key, by = c("dataset", "model_id", "model_family", "formula_id"))) == 0L &&
    nrow(dplyr::anti_join(model_registry_key, plan_key, by = c("dataset", "model_id", "model_family", "formula_id"))) == 0L
}

validate_existing_subject_horizon_risk_long <- function(subject_horizon_risk_long, model_registry, analysis_datasets, common_horizons_year) {
  if (is.null(subject_horizon_risk_long) || is.null(model_registry)) {
    return(FALSE)
  }

  required_cols <- c("dataset", "model_id", "unique_person_id", "horizon_year", "survival_pred", "risk_pred", "formula_id")
  if (!all(required_cols %in% names(subject_horizon_risk_long))) {
    return(FALSE)
  }

  risk_long_key_count <- subject_horizon_risk_long %>%
    transmute(
      dataset = normalize_dataset_label(dataset),
      model_id = as.character(model_id),
      unique_person_id = as.character(unique_person_id),
      horizon_year = as.integer(horizon_year)
    ) %>%
    distinct() %>%
    nrow()

  if (risk_long_key_count != nrow(subject_horizon_risk_long)) {
    return(FALSE)
  }

  converged_registry <- model_registry %>%
    mutate(dataset = normalize_dataset_label(dataset)) %>%
    filter(as.logical(converged))

  expected_prediction_rows <- if (nrow(converged_registry) > 0) {
    sum(
      vapply(
        converged_registry$dataset,
        function(dataset_name) {
          nrow(analysis_datasets[[dataset_name]]) * length(common_horizons_year)
        },
        numeric(1)
      )
    )
  } else {
    0
  }

  actual_prediction_rows <- subject_horizon_risk_long %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      horizon_year = as.integer(horizon_year)
    ) %>%
    filter(model_id %in% converged_registry$model_id, horizon_year %in% common_horizons_year) %>%
    nrow()

  expected_prediction_rows == actual_prediction_rows
}

validate_existing_model_performance_long <- function(model_performance_long, risk_thresholds) {
  if (is.null(model_performance_long)) {
    return(FALSE)
  }

  required_cols <- c("dataset", "model_id", "metric_domain", "metric_name", "metric_value")
  if (!all(required_cols %in% names(model_performance_long))) {
    return(FALSE)
  }

  if ("threshold" %in% names(model_performance_long)) {
    thresholds_present <- model_performance_long %>%
      filter(!is.na(threshold)) %>%
      pull(threshold) %>%
      unique() %>%
      as.numeric()

    if (length(thresholds_present) > 0L && !setequal(thresholds_present, risk_thresholds)) {
      return(FALSE)
    }
  }

  TRUE
}

validate_existing_calibration_tables <- function(calibration_bins_long, calibration_diagnostics_long) {
  if (is.null(calibration_bins_long) || is.null(calibration_diagnostics_long)) {
    return(FALSE)
  }

  bins_required <- c("dataset", "model_id", "horizon_year", "bin_index", "mean_predicted_risk")
  diag_required <- c("dataset", "model_id", "horizon_year", "binary_outcome_support_flag", "calibration_support_flag")

  all(bins_required %in% names(calibration_bins_long)) &&
    all(diag_required %in% names(calibration_diagnostics_long))
}

can_reuse_existing_stage5_core <- function(existing_core, model_plan, analysis_datasets, common_horizons_year, risk_thresholds) {
  isTRUE(
    validate_existing_model_registry(existing_core$model_registry, model_plan) &&
      validate_existing_subject_horizon_risk_long(existing_core$subject_horizon_risk_long, existing_core$model_registry, analysis_datasets, common_horizons_year) &&
      validate_existing_model_performance_long(existing_core$model_performance_long, risk_thresholds) &&
      validate_existing_calibration_tables(existing_core$calibration_bins_long, existing_core$calibration_diagnostics_long)
  )
}

# 🔴 Define: QC builders and validators ===============================
build_stage5_qc_summary <- function(
    model_plan,
    model_registry,
    subject_horizon_risk_long,
    model_performance_long,
    calibration_diagnostics_long,
    calibration_bins_long,
    analysis_datasets,
    common_horizons_year
) {
  converged_registry <- model_registry %>%
    filter(converged == TRUE)

  time_range_score_rule_violation_count <- if (nrow(converged_registry) > 0) {
    sum(
      is.na(converged_registry$time_range_score_type) |
        converged_registry$time_range_score_type != paste0("anchor_risk_year_", converged_registry$anchor_horizon_year),
      na.rm = TRUE
    )
  } else {
    0L
  }

  expected_prediction_rows <- if (nrow(converged_registry) > 0) {
    sum(
      vapply(
        converged_registry$dataset,
        function(dataset_name) {
          nrow(analysis_datasets[[dataset_name]]) * length(common_horizons_year)
        },
        numeric(1)
      )
    )
  } else {
    0
  }

  risk_range_violation_count <- if (nrow(subject_horizon_risk_long) > 0) {
    sum(
      subject_horizon_risk_long$risk_pred < -qc_tolerance |
        subject_horizon_risk_long$risk_pred > 1 + qc_tolerance,
      na.rm = TRUE
    )
  } else {
    0L
  }

  survival_range_violation_count <- if (nrow(subject_horizon_risk_long) > 0) {
    sum(
      subject_horizon_risk_long$survival_pred < -qc_tolerance |
        subject_horizon_risk_long$survival_pred > 1 + qc_tolerance,
      na.rm = TRUE
    )
  } else {
    0L
  }

  max_survival_plus_risk_error <- if (nrow(subject_horizon_risk_long) > 0) {
    safe_max_abs(subject_horizon_risk_long$survival_pred + subject_horizon_risk_long$risk_pred - 1)
  } else {
    NA_real_
  }

  monotone_risk_violation_count <- if (nrow(subject_horizon_risk_long) > 0) {
    monotone_check <- subject_horizon_risk_long %>%
      arrange(model_id, unique_person_id, horizon_year) %>%
      group_by(model_id, unique_person_id) %>%
      summarise(
        monotone_ok = {
          d <- diff(risk_pred)
          length(d) == 0L || all(is.finite(d) & d >= -qc_tolerance)
        },
        .groups = "drop"
      )

    sum(!monotone_check$monotone_ok, na.rm = TRUE)
  } else {
    0L
  }

  treat_all_consistency <- model_performance_long %>%
    filter(
      model_family == "treat_all",
      metric_name %in% c("net_benefit", "net_benefit_treat_all")
    ) %>%
    select(dataset, horizon_year, threshold, metric_name, metric_value) %>%
    distinct() %>%
    tidyr::pivot_wider(
      names_from = metric_name,
      values_from = metric_value
    )

  treat_all_net_benefit_max_abs_diff <- if (nrow(treat_all_consistency) > 0) {
    safe_max_abs(treat_all_consistency$net_benefit - treat_all_consistency$net_benefit_treat_all)
  } else {
    NA_real_
  }

  horizon_metric_qc <- model_performance_long %>%
    filter(
      model_class == "non_cure",
      metric_domain == "horizon_summary",
      metric_name %in% c("time_dependent_auc", "brier_score", "binary_outcome_support_flag", "observed_ipcw_risk")
    ) %>%
    select(dataset, model_id, horizon_year, metric_name, metric_value) %>%
    distinct() %>%
    tidyr::pivot_wider(
      names_from = metric_name,
      values_from = metric_value
    )

  auc_missing_with_support_count <- if (nrow(horizon_metric_qc) > 0) {
    sum(
      is.na(horizon_metric_qc$time_dependent_auc) &
        horizon_metric_qc$binary_outcome_support_flag == 1,
      na.rm = TRUE
    )
  } else {
    0L
  }

  brier_missing_with_support_count <- if (nrow(horizon_metric_qc) > 0) {
    sum(
      is.na(horizon_metric_qc$brier_score) &
        horizon_metric_qc$binary_outcome_support_flag == 1,
      na.rm = TRUE
    )
  } else {
    0L
  }

  unsupported_auc_nonmissing_count <- if (nrow(horizon_metric_qc) > 0) {
    sum(
      !is.na(horizon_metric_qc$time_dependent_auc) &
        horizon_metric_qc$binary_outcome_support_flag == 0,
      na.rm = TRUE
    )
  } else {
    0L
  }

  unsupported_brier_nonmissing_count <- if (nrow(horizon_metric_qc) > 0) {
    sum(
      !is.na(horizon_metric_qc$brier_score) &
        horizon_metric_qc$binary_outcome_support_flag == 0,
      na.rm = TRUE
    )
  } else {
    0L
  }

  supported_observed_ipcw_missing_count <- if (nrow(horizon_metric_qc) > 0) {
    sum(
      is.na(horizon_metric_qc$observed_ipcw_risk) &
        horizon_metric_qc$binary_outcome_support_flag == 1,
      na.rm = TRUE
    )
  } else {
    0L
  }

  unsupported_observed_ipcw_nonmissing_count <- if (nrow(horizon_metric_qc) > 0) {
    sum(
      !is.na(horizon_metric_qc$observed_ipcw_risk) &
        horizon_metric_qc$binary_outcome_support_flag == 0,
      na.rm = TRUE
    )
  } else {
    0L
  }

  threshold_metric_qc <- model_performance_long %>%
    filter(
      metric_domain %in% c("threshold_summary", "threshold_reference"),
      metric_name %in% c("binary_outcome_support_flag", threshold_outcome_metric_names)
    ) %>%
    select(dataset, model_id, model_family, horizon_year, threshold, metric_domain, metric_name, metric_value) %>%
    distinct() %>%
    tidyr::pivot_wider(
      names_from = metric_name,
      values_from = metric_value
    )

  unsupported_threshold_outcome_row_count <- if (nrow(threshold_metric_qc) > 0) {
    threshold_outcome_cols <- intersect(threshold_outcome_metric_names, names(threshold_metric_qc))
    outcome_any_nonmissing <- if (length(threshold_outcome_cols) > 0) {
      apply(
        X = as.matrix(threshold_metric_qc[, threshold_outcome_cols, drop = FALSE]),
        MARGIN = 1,
        FUN = function(x) any(!is.na(x))
      )
    } else {
      rep(FALSE, nrow(threshold_metric_qc))
    }

    sum(
      threshold_metric_qc$binary_outcome_support_flag == 0 &
        outcome_any_nonmissing,
      na.rm = TRUE
    )
  } else {
    0L
  }

  observed_ipcw_metric_qc <- model_performance_long %>%
    filter(
      metric_domain %in% c("horizon_reference", "horizon_summary"),
      metric_name %in% c("binary_outcome_support_flag", "observed_ipcw_risk")
    ) %>%
    select(dataset, model_id, horizon_year, metric_domain, metric_name, metric_value) %>%
    distinct() %>%
    tidyr::pivot_wider(
      names_from = metric_name,
      values_from = metric_value
    )

  observed_ipcw_support_missing_count <- if (nrow(observed_ipcw_metric_qc) > 0) {
    sum(
      is.na(observed_ipcw_metric_qc$observed_ipcw_risk) &
        observed_ipcw_metric_qc$binary_outcome_support_flag == 1,
      na.rm = TRUE
    )
  } else {
    0L
  }

  observed_ipcw_unsupported_nonmissing_count <- if (nrow(observed_ipcw_metric_qc) > 0) {
    sum(
      !is.na(observed_ipcw_metric_qc$observed_ipcw_risk) &
        observed_ipcw_metric_qc$binary_outcome_support_flag == 0,
      na.rm = TRUE
    )
  } else {
    0L
  }

  calibration_bin_qc <- calibration_bins_long %>%
    transmute(
      binary_outcome_support_flag = as.integer(binary_outcome_support_flag),
      bin_binary_support_flag = as.integer(bin_binary_support_flag),
      observed_risk_ipcw = as.numeric(observed_risk_ipcw)
    )

  unsupported_bin_observed_risk_nonmissing_count <- if (nrow(calibration_bin_qc) > 0) {
    sum(
      !is.na(calibration_bin_qc$observed_risk_ipcw) &
        (calibration_bin_qc$binary_outcome_support_flag == 0 | calibration_bin_qc$bin_binary_support_flag == 0),
      na.rm = TRUE
    )
  } else {
    0L
  }

  supported_bin_observed_risk_missing_count <- if (nrow(calibration_bin_qc) > 0) {
    sum(
      is.na(calibration_bin_qc$observed_risk_ipcw) &
        calibration_bin_qc$binary_outcome_support_flag == 1 &
        calibration_bin_qc$bin_binary_support_flag == 1,
      na.rm = TRUE
    )
  } else {
    0L
  }

  supported_calibration_rows <- calibration_diagnostics_long %>%
    filter(calibration_support_flag == 1)

  successful_calibration_row_count <- if (nrow(supported_calibration_rows) > 0) {
    sum(supported_calibration_rows$calibration_fit_success_flag == 1, na.rm = TRUE)
  } else {
    0L
  }

  calibration_boundary_row_count <- if (nrow(supported_calibration_rows) > 0) {
    sum(
      supported_calibration_rows$calibration_offset_boundary_flag == 1 |
        supported_calibration_rows$calibration_slope_boundary_flag == 1,
      na.rm = TRUE
    )
  } else {
    0L
  }

  calibration_nonconverged_row_count <- if (nrow(supported_calibration_rows) > 0) {
    sum(
      !(supported_calibration_rows$calibration_offset_converged_flag == 1 &
          supported_calibration_rows$calibration_slope_converged_flag == 1),
      na.rm = TRUE
    )
  } else {
    0L
  }

  calibration_suppressed_row_count <- if (nrow(supported_calibration_rows) > 0) {
    sum(
      supported_calibration_rows$calibration_regression_suppressed_flag == 1,
      na.rm = TRUE
    )
  } else {
    0L
  }

  unstable_calibration_estimate_nonmissing_count <- if (nrow(supported_calibration_rows) > 0) {
    sum(
      supported_calibration_rows$calibration_fit_success_flag == 0 &
        (
          !is.na(supported_calibration_rows$calibration_intercept_offset) |
            !is.na(supported_calibration_rows$calibration_intercept_free) |
            !is.na(supported_calibration_rows$calibration_slope)
        ),
      na.rm = TRUE
    )
  } else {
    0L
  }

  time_range_metric_qc <- model_performance_long %>%
    filter(
      metric_domain == "time_range_discrimination",
      metric_name %in% c("harrell_c", "uno_c")
    ) %>%
    select(dataset, model_id, metric_name, metric_value) %>%
    distinct() %>%
    tidyr::pivot_wider(
      names_from = metric_name,
      values_from = metric_value
    )

  uno_harrell_large_gap_count <- if (nrow(time_range_metric_qc) > 0) {
    sum(
      is.finite(time_range_metric_qc$harrell_c) &
        is.finite(time_range_metric_qc$uno_c) &
        abs(time_range_metric_qc$harrell_c - time_range_metric_qc$uno_c) > uno_harrell_discordance_tolerance,
      na.rm = TRUE
    )
  } else {
    0L
  }

  max_uno_harrell_abs_diff <- if (nrow(time_range_metric_qc) > 0) {
    safe_max_abs(time_range_metric_qc$harrell_c - time_range_metric_qc$uno_c)
  } else {
    NA_real_
  }

  cox_model_registry_qc <- model_registry %>%
    filter(model_family == "coxph")

  parametric_model_registry_qc <- model_registry %>%
    filter(model_family %in% c("exponential", "weibull", "lognormal", "loglogistic"))

  cox_full_ic_nonmissing_count <- if (nrow(cox_model_registry_qc) > 0) {
    sum(
      !is.na(cox_model_registry_qc$logLik) |
        !is.na(cox_model_registry_qc$AIC) |
        !is.na(cox_model_registry_qc$BIC),
      na.rm = TRUE
    )
  } else {
    0L
  }

  cox_partial_ic_missing_count <- if (nrow(cox_model_registry_qc) > 0) {
    sum(
      is.na(cox_model_registry_qc$partial_logLik) |
        is.na(cox_model_registry_qc$partial_AIC) |
        is.na(cox_model_registry_qc$partial_BIC),
      na.rm = TRUE
    )
  } else {
    0L
  }

  parametric_partial_ic_nonmissing_count <- if (nrow(parametric_model_registry_qc) > 0) {
    sum(
      !is.na(parametric_model_registry_qc$partial_logLik) |
        !is.na(parametric_model_registry_qc$partial_AIC) |
        !is.na(parametric_model_registry_qc$partial_BIC),
      na.rm = TRUE
    )
  } else {
    0L
  }

  information_criterion_flag_mismatch_count <- if (nrow(model_registry) > 0) {
    sum(
      (
        model_registry$model_family == "coxph" &
          model_registry$information_criterion_comparable_across_families_flag != 0L
      ) |
        (
          model_registry$model_family %in% c("exponential", "weibull", "lognormal", "loglogistic") &
            model_registry$information_criterion_comparable_across_families_flag != 1L
        ),
      na.rm = TRUE
    )
  } else {
    0L
  }

  low_known_count_warning_row_count <- if (nrow(calibration_diagnostics_long) > 0) {
    sum(calibration_diagnostics_long$low_known_count_warning_flag == 1L, na.rm = TRUE)
  } else {
    0L
  }

  late_horizon_instability_row_count <- if (nrow(calibration_diagnostics_long) > 0) {
    sum(calibration_diagnostics_long$late_horizon_instability_flag == 1L, na.rm = TRUE)
  } else {
    0L
  }

  risk_scale_mismatch_count <- if (nrow(model_registry) > 0) {
    sum(model_registry$risk_scale != "transition_only_main", na.rm = TRUE)
  } else {
    0L
  }

  effective_parameter_missing_count <- if (nrow(model_registry) > 0) {
    sum(is.na(model_registry$effective_parameter_count), na.rm = TRUE)
  } else {
    0L
  }

  horizon_metadata_missing_count <- if (nrow(model_performance_long) > 0) {
    sum(
      !is.na(model_performance_long$horizon_year) &
        (
          is.na(model_performance_long$support_tier) |
            is.na(model_performance_long$horizon_evidence_class) |
            is.na(model_performance_long$claim_restriction_flag)
        ),
      na.rm = TRUE
    )
  } else {
    0L
  }

  projection_primary_claim_mismatch_count <- if (nrow(model_performance_long) > 0) {
    sum(
      model_performance_long$horizon_evidence_class == "mostly_extrapolated" &
        model_performance_long$claim_restriction_flag == "primary_claim_allowed",
      na.rm = TRUE
    )
  } else {
    0L
  }

  tibble(
    stage = stage_name,
    stage_role = stage_role,
    stage_spec_version = stage_spec_version,
    planned_model_count = nrow(model_plan),
    fitted_model_registry_rows = nrow(model_registry),
    converged_model_count = sum(model_registry$converged == TRUE, na.rm = TRUE),
    time_range_score_rule_violation_count = as.integer(time_range_score_rule_violation_count),
    expected_prediction_rows = expected_prediction_rows,
    actual_prediction_rows = nrow(subject_horizon_risk_long),
    risk_range_violation_count = as.integer(risk_range_violation_count),
    survival_range_violation_count = as.integer(survival_range_violation_count),
    max_survival_plus_risk_error = max_survival_plus_risk_error,
    monotone_risk_violation_count = as.integer(monotone_risk_violation_count),
    auc_missing_with_support_count = as.integer(auc_missing_with_support_count),
    brier_missing_with_support_count = as.integer(brier_missing_with_support_count),
    unsupported_auc_nonmissing_count = as.integer(unsupported_auc_nonmissing_count),
    unsupported_brier_nonmissing_count = as.integer(unsupported_brier_nonmissing_count),
    supported_observed_ipcw_missing_count = as.integer(supported_observed_ipcw_missing_count),
    unsupported_observed_ipcw_nonmissing_count = as.integer(unsupported_observed_ipcw_nonmissing_count),
    unsupported_threshold_outcome_row_count = as.integer(unsupported_threshold_outcome_row_count),
    observed_ipcw_support_missing_count = as.integer(observed_ipcw_support_missing_count),
    observed_ipcw_unsupported_nonmissing_count = as.integer(observed_ipcw_unsupported_nonmissing_count),
    supported_bin_observed_risk_missing_count = as.integer(supported_bin_observed_risk_missing_count),
    unsupported_bin_observed_risk_nonmissing_count = as.integer(unsupported_bin_observed_risk_nonmissing_count),
    supported_calibration_row_count = nrow(supported_calibration_rows),
    successful_calibration_row_count = as.integer(successful_calibration_row_count),
    calibration_boundary_row_count = as.integer(calibration_boundary_row_count),
    calibration_nonconverged_row_count = as.integer(calibration_nonconverged_row_count),
    calibration_suppressed_row_count = as.integer(calibration_suppressed_row_count),
    unstable_calibration_estimate_nonmissing_count = as.integer(unstable_calibration_estimate_nonmissing_count),
    uno_harrell_large_gap_count = as.integer(uno_harrell_large_gap_count),
    max_uno_harrell_abs_diff = max_uno_harrell_abs_diff,
    cox_full_ic_nonmissing_count = as.integer(cox_full_ic_nonmissing_count),
    cox_partial_ic_missing_count = as.integer(cox_partial_ic_missing_count),
    parametric_partial_ic_nonmissing_count = as.integer(parametric_partial_ic_nonmissing_count),
    information_criterion_flag_mismatch_count = as.integer(information_criterion_flag_mismatch_count),
    low_known_count_warning_row_count = as.integer(low_known_count_warning_row_count),
    late_horizon_instability_row_count = as.integer(late_horizon_instability_row_count),
    risk_scale_mismatch_count = as.integer(risk_scale_mismatch_count),
    effective_parameter_missing_count = as.integer(effective_parameter_missing_count),
    horizon_metadata_missing_count = as.integer(horizon_metadata_missing_count),
    projection_primary_claim_mismatch_count = as.integer(projection_primary_claim_mismatch_count),
    treat_all_net_benefit_max_abs_diff = treat_all_net_benefit_max_abs_diff
  )
}

validate_stage5_outputs <- function(qc_summary) {
  qc <- qc_summary[1, , drop = FALSE]

  if (qc$converged_model_count <= 0) {
    stop("Stage 5 produced zero converged models.", call. = FALSE)
  }

  if (qc$planned_model_count != qc$fitted_model_registry_rows) {
    stop("Stage 5 model registry row count does not match the planned model count.", call. = FALSE)
  }

  if (qc$time_range_score_rule_violation_count > 0) {
    stop("Stage 5 time-range discrimination score type no longer follows the unified anchor-risk rule.", call. = FALSE)
  }

  if (qc$expected_prediction_rows != qc$actual_prediction_rows) {
    stop("Stage 5 prediction row count does not match the expected count.", call. = FALSE)
  }

  if (qc$risk_range_violation_count > 0) {
    stop("Stage 5 predicted risks left the [0, 1] range.", call. = FALSE)
  }

  if (qc$survival_range_violation_count > 0) {
    stop("Stage 5 predicted survivals left the [0, 1] range.", call. = FALSE)
  }

  if (is.finite(qc$max_survival_plus_risk_error) && qc$max_survival_plus_risk_error > qc_tolerance) {
    stop("Stage 5 survival and risk predictions no longer sum to 1 within tolerance.", call. = FALSE)
  }

  if (qc$monotone_risk_violation_count > 0) {
    stop("Stage 5 predicted risk is not monotone non-decreasing across horizons for at least one subject.", call. = FALSE)
  }

  if (qc$auc_missing_with_support_count > 0) {
    stop("Stage 5 AUC is missing despite binary outcome support being available.", call. = FALSE)
  }

  if (qc$brier_missing_with_support_count > 0) {
    stop("Stage 5 Brier score is missing despite binary outcome support being available.", call. = FALSE)
  }

  if (qc$unsupported_auc_nonmissing_count > 0) {
    stop("Stage 5 AUC remained populated despite binary outcome support being unavailable.", call. = FALSE)
  }

  if (qc$unsupported_brier_nonmissing_count > 0) {
    stop("Stage 5 Brier score remained populated despite binary outcome support being unavailable.", call. = FALSE)
  }

  if (qc$supported_observed_ipcw_missing_count > 0 || qc$observed_ipcw_support_missing_count > 0) {
    stop("Stage 5 observed IPCW risk is missing despite binary outcome support being available.", call. = FALSE)
  }

  if (qc$unsupported_observed_ipcw_nonmissing_count > 0 || qc$observed_ipcw_unsupported_nonmissing_count > 0) {
    stop("Stage 5 observed IPCW risk remained populated despite binary outcome support being unavailable.", call. = FALSE)
  }

  if (qc$unsupported_threshold_outcome_row_count > 0) {
    stop("Stage 5 threshold-based outcome metrics remained populated despite binary outcome support being unavailable.", call. = FALSE)
  }

  if (qc$supported_bin_observed_risk_missing_count > 0) {
    stop("Stage 5 calibration-bin observed risk is missing despite bin-level binary outcome support being available.", call. = FALSE)
  }

  if (qc$unsupported_bin_observed_risk_nonmissing_count > 0) {
    stop("Stage 5 calibration-bin observed risk remained populated despite bin-level binary outcome support being unavailable.", call. = FALSE)
  }

  if (qc$cox_full_ic_nonmissing_count > 0) {
    stop("Stage 5 coxph full-likelihood-comparable information-criterion fields should be NA.", call. = FALSE)
  }

  if (qc$cox_partial_ic_missing_count > 0) {
    stop("Stage 5 coxph partial-likelihood information-criterion fields are missing.", call. = FALSE)
  }

  if (qc$parametric_partial_ic_nonmissing_count > 0) {
    stop("Stage 5 parametric models unexpectedly populated partial-likelihood information-criterion fields.", call. = FALSE)
  }

  if (qc$information_criterion_flag_mismatch_count > 0) {
    stop("Stage 5 information-criterion comparability flags are inconsistent with the model family.", call. = FALSE)
  }

  if (qc$risk_scale_mismatch_count > 0) {
    stop("Stage 5 exported a non-main risk scale inside the non-cure block.", call. = FALSE)
  }

  if (qc$effective_parameter_missing_count > 0) {
    stop("Stage 5 effective parameter counts are missing.", call. = FALSE)
  }

  if (qc$horizon_metadata_missing_count > 0) {
    stop("Stage 5 horizon metadata fields are missing for one or more horizon-specific rows.", call. = FALSE)
  }

  if (qc$projection_primary_claim_mismatch_count > 0) {
    stop("Stage 5 exported a projection-dominant row with a primary-claim flag.", call. = FALSE)
  }

  if (is.finite(qc$treat_all_net_benefit_max_abs_diff) && qc$treat_all_net_benefit_max_abs_diff > qc_tolerance) {
    stop("Stage 5 treat-all net benefit self-consistency check failed.", call. = FALSE)
  }

  if (qc$supported_calibration_row_count > 0 && qc$successful_calibration_row_count == 0) {
    stop("Stage 5 calibration fitting failed for every supported horizon.", call. = FALSE)
  }

  if (qc$unstable_calibration_estimate_nonmissing_count > 0) {
    stop("Stage 5 unstable calibration regression estimates were not suppressed to NA.", call. = FALSE)
  }

  if (qc$supported_calibration_row_count > 0 &&
      (qc$calibration_boundary_row_count > 0 || qc$calibration_nonconverged_row_count > 0)) {
    warning(
      "Stage 5 retained supported horizons with calibration boundary or nonconverged fits; inspect `stage5_calibration_diagnostics_long.csv`.",
      call. = FALSE
    )
  }

  if (qc$late_horizon_instability_row_count > 0) {
    warning(
      sprintf(
        "Stage 5 retained %s supported late-horizon rows with low known_count; review `late_horizon_instability_flag` in the exported tables.",
        qc$late_horizon_instability_row_count
      ),
      call. = FALSE
    )
  }

  if (is.finite(qc$max_uno_harrell_abs_diff) && qc$max_uno_harrell_abs_diff > uno_harrell_discordance_tolerance) {
    warning(
      "Stage 5 Harrell/Uno time-range concordance gap exceeded the review tolerance; inspect score choice and censoring weighting.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# 🔴 Define: plot builders and export helpers ===============================
make_empty_plot <- function(title, subtitle = NULL) {
  ggplot(data.frame(x = 0, y = 0), aes(x = x, y = y)) +
    geom_blank() +
    annotate("text", x = 0, y = 0, label = "No plottable data available.") +
    xlim(-1, 1) +
    ylim(-1, 1) +
    labs(title = title, subtitle = subtitle) +
    theme_void()
}

save_plot_png <- function(plot_obj, file_path, width = plot_width_in, height = plot_height_in, dpi = plot_dpi) {
  ggplot2::ggsave(
    filename = file_path,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    units = "in"
  )
}

save_plot_pdf_bundle <- function(plot_list, pdf_path, width = plot_width_in, height = plot_height_in) {
  grDevices::pdf(pdf_path, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (plot_obj in plot_list) {
    print(plot_obj)
  }
}

prepare_stage5_plot_source_backbone <- function(df) {
  if (is.null(df)) {
    return(tibble())
  }

  out <- collapse_dataset_key_suffixes(df)

  if ("dataset" %in% names(out)) {
    out <- out %>%
      mutate(dataset = normalize_dataset_label(dataset))
  }

  if (!"dataset_key" %in% names(out) && "dataset" %in% names(out)) {
    out <- out %>%
      mutate(dataset_key = dataset)
  }

  if (!"formula_name" %in% names(out)) {
    if ("formula_label" %in% names(out)) {
      out <- out %>%
        mutate(formula_name = as.character(formula_label))
    } else if ("model_family" %in% names(out)) {
      out <- out %>%
        mutate(formula_name = if_else(model_family %in% c("observed_km", "treat_all", "treat_none"), "reference", "unknown_formula_name"))
    } else {
      out <- out %>%
        mutate(formula_name = "unknown_formula_name")
    }
  }

  if (!"formula_label" %in% names(out)) {
    out <- out %>%
      mutate(formula_label = if_else(formula_name == "reference", "Reference", as.character(formula_name)))
  }

  if (!"formula_id" %in% names(out)) {
    out <- out %>%
      mutate(formula_id = if_else(formula_name == "reference", "reference", paste0("derived_", make.names(formula_name))))
  }

  if (!"formula_scope" %in% names(out)) {
    out <- out %>%
      mutate(formula_scope = if_else(formula_name == "reference", "benchmark_reference", "main_transition_only_scale"))
  }

  if ("horizon_year" %in% names(out)) {
    out <- out %>%
      mutate(horizon_year = as.integer(horizon_year))
  }

  if ("threshold" %in% names(out)) {
    out <- out %>%
      mutate(threshold = as.numeric(threshold))
  }

  if ("metric_value" %in% names(out)) {
    out <- out %>%
      mutate(metric_value = as.numeric(metric_value))
  }

  out
}

build_stage5_plot_sources <- function(
    subject_horizon_risk_long,
    model_performance_long,
    calibration_bins_long,
    common_horizons_year,
    prediction_horizons_for_plots
) {
  subject_horizon_risk_long <- prepare_stage5_plot_source_backbone(subject_horizon_risk_long)
  model_performance_long <- prepare_stage5_plot_source_backbone(model_performance_long)
  calibration_bins_long <- prepare_stage5_plot_source_backbone(calibration_bins_long)

  risk_plot_source <- subject_horizon_risk_long %>%
    filter(model_class == "non_cure") %>%
    group_by(dataset, dataset_key, formula_id, formula_name, formula_label, formula_scope, model_family, horizon_year) %>%
    summarise(
      median_risk = median(risk_pred, na.rm = TRUE),
      risk_p25 = stats::quantile(risk_pred, probs = 0.25, na.rm = TRUE, names = FALSE),
      risk_p75 = stats::quantile(risk_pred, probs = 0.75, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    ) %>%
    apply_model_family_plot_label() %>%
    apply_dataset_formula_panel()

  auc_plot_source <- model_performance_long %>%
    filter(
      metric_name == "time_dependent_auc",
      model_class == "non_cure",
      !is.na(metric_value)
    ) %>%
    arrange(dataset, formula_label, model_family, horizon_year) %>%
    group_by(dataset, formula_label, model_family) %>%
    mutate(plot_segment_id = cumsum(c(TRUE, diff(horizon_year) != 1L))) %>%
    ungroup() %>%
    apply_model_family_plot_label() %>%
    apply_dataset_formula_panel()

  brier_plot_source <- model_performance_long %>%
    filter(
      metric_name == "brier_score",
      model_class == "non_cure",
      !is.na(metric_value)
    ) %>%
    arrange(dataset, formula_label, model_family, horizon_year) %>%
    group_by(dataset, formula_label, model_family) %>%
    mutate(plot_segment_id = cumsum(c(TRUE, diff(horizon_year) != 1L))) %>%
    ungroup() %>%
    apply_model_family_plot_label() %>%
    apply_dataset_formula_panel()

  calibration_plot_source <- calibration_bins_long %>%
    filter(
      horizon_year %in% prediction_horizons_for_plots,
      binary_outcome_support_flag == 1,
      bin_binary_support_flag == 1,
      !is.na(mean_predicted_risk),
      !is.na(observed_risk_ipcw)
    ) %>%
    arrange(dataset, formula_label, model_family, horizon_year, bin_index) %>%
    group_by(dataset, formula_label, model_family, horizon_year) %>%
    mutate(calibration_line_segment_id = cumsum(c(TRUE, diff(bin_index) != 1L))) %>%
    ungroup() %>%
    apply_model_family_plot_label()

  net_benefit_plot_source <- model_performance_long %>%
    filter(
      metric_name == "net_benefit",
      horizon_year %in% prediction_horizons_for_plots,
      model_family %in% c(
        "exponential",
        "weibull",
        "lognormal",
        "loglogistic",
        "coxph",
        "treat_all",
        "treat_none"
      ),
      !is.na(metric_value)
    ) %>%
    mutate(
      formula_label = if_else(is.na(formula_label), "Reference", formula_label)
    ) %>%
    apply_model_family_plot_label()

  list(
    risk = risk_plot_source,
    auc = auc_plot_source,
    brier = brier_plot_source,
    calibration = calibration_plot_source,
    net_benefit = net_benefit_plot_source
  )
}

build_stage5_plot_objects <- function(plot_sources, common_horizons_year, prediction_horizons_for_plots) {
  selected_plot_horizon_caption <- paste(prediction_horizons_for_plots, collapse = ", ")

  if (nrow(plot_sources$risk) > 0) {
    risk_plot <- ggplot(
      plot_sources$risk,
      aes(x = horizon_year, y = median_risk, color = model_family_plot_label, fill = model_family_plot_label)
    ) +
      geom_line(linewidth = 0.8) +
      geom_ribbon(aes(ymin = risk_p25, ymax = risk_p75), alpha = 0.15, linewidth = 0) +
      facet_wrap(~ dataset_formula_panel, scales = "free_y", ncol = 4) +
      scale_x_continuous(breaks = common_horizons_year) +
      labs(
        title = "Stage 5 non-cure predicted risk summary",
        subtitle = "Panels follow the current Stage 1 formula registry and the transition-only main risk scale.",
        x = "Prediction horizon (years)",
        y = "Median predicted risk",
        color = "Model family",
        fill = "Model family"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
  } else {
    risk_plot <- make_empty_plot(
      title = "Stage 5 non-cure predicted risk summary",
      subtitle = "No plottable risk-summary rows were available."
    )
  }

  auc_line_source <- plot_sources$auc %>%
    add_count(dataset, formula_label, model_family, plot_segment_id, name = "segment_n") %>%
    filter(segment_n >= 2L)

  if (nrow(plot_sources$auc) > 0) {
    auc_plot <- ggplot(plot_sources$auc, aes(x = horizon_year, y = metric_value, color = model_family_plot_label)) +
      geom_line(
        data = auc_line_source,
        aes(group = interaction(model_family_plot_label, plot_segment_id)),
        linewidth = 0.8
      ) +
      geom_point(size = 1.5) +
      facet_wrap(~ dataset_formula_panel, ncol = 4) +
      scale_x_continuous(breaks = common_horizons_year) +
      labs(
        title = "Stage 5 cumulative/dynamic time-dependent AUC",
        subtitle = "Shown only where the binary IPCW horizon outcome is estimable; late supported horizons with low known_count remain flagged in the source tables.",
        x = "Prediction horizon (years)",
        y = "AUC",
        color = "Model family"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
  } else {
    auc_plot <- make_empty_plot(
      title = "Stage 5 cumulative/dynamic time-dependent AUC",
      subtitle = "No plottable AUC rows were available."
    )
  }

  brier_line_source <- plot_sources$brier %>%
    add_count(dataset, formula_label, model_family, plot_segment_id, name = "segment_n") %>%
    filter(segment_n >= 2L)

  if (nrow(plot_sources$brier) > 0) {
    brier_plot <- ggplot(plot_sources$brier, aes(x = horizon_year, y = metric_value, color = model_family_plot_label)) +
      geom_line(
        data = brier_line_source,
        aes(group = interaction(model_family_plot_label, plot_segment_id)),
        linewidth = 0.8
      ) +
      geom_point(size = 1.5) +
      facet_wrap(~ dataset_formula_panel, ncol = 4) +
      scale_x_continuous(breaks = common_horizons_year) +
      labs(
        title = "Stage 5 horizon-specific Brier score",
        subtitle = "Shown only where the binary IPCW horizon outcome is estimable; late supported horizons with low known_count remain flagged in the source tables.",
        x = "Prediction horizon (years)",
        y = "Brier score",
        color = "Model family"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
  } else {
    brier_plot <- make_empty_plot(
      title = "Stage 5 horizon-specific Brier score",
      subtitle = "No plottable Brier-score rows were available."
    )
  }

  calibration_line_source <- plot_sources$calibration %>%
    add_count(
      dataset, formula_label, model_family, horizon_year, calibration_line_segment_id,
      name = "segment_n"
    ) %>%
    filter(segment_n >= 2L)

  if (nrow(plot_sources$calibration) > 0) {
    calibration_plot <- ggplot(
      plot_sources$calibration,
      aes(x = mean_predicted_risk, y = observed_risk_ipcw, color = model_family_plot_label)
    ) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      geom_point(size = 1.5) +
      geom_line(
        data = calibration_line_source,
        aes(group = interaction(model_family_plot_label, calibration_line_segment_id)),
        linewidth = 0.8
      ) +
      facet_grid(dataset + formula_label ~ horizon_year, scales = "free") +
      labs(
        title = sprintf("Stage 5 grouped calibration at selected horizons (%s years)", selected_plot_horizon_caption),
        subtitle = "Selected snapshot panels may include primary-supported, secondary, sensitivity, or projection horizons depending on dataset; supported late horizons with low known_count are flagged in the source tables.",
        x = "Mean predicted risk",
        y = "Observed risk (IPCW)",
        color = "Model family"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
  } else {
    calibration_plot <- make_empty_plot(
      title = sprintf("Stage 5 grouped calibration at selected horizons (%s years)", selected_plot_horizon_caption),
      subtitle = "No plottable grouped calibration rows were available."
    )
  }

  if (nrow(plot_sources$net_benefit) > 0) {
    net_benefit_plot <- ggplot(
      plot_sources$net_benefit,
      aes(x = threshold, y = metric_value, color = model_family_plot_label)
    ) +
      geom_line(linewidth = 0.8) +
      facet_grid(dataset + formula_label ~ horizon_year, scales = "free_y") +
      labs(
        title = sprintf("Stage 5 decision-curve net benefit at selected horizons (%s years)", selected_plot_horizon_caption),
        subtitle = "Selected snapshot panels may include primary-supported, secondary, sensitivity, or projection horizons depending on dataset; late supported horizons with low known_count remain flagged in the source tables.",
        x = "Risk threshold",
        y = "Net benefit",
        color = "Model family"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
  } else {
    net_benefit_plot <- make_empty_plot(
      title = sprintf("Stage 5 decision-curve net benefit at selected horizons (%s years)", selected_plot_horizon_caption),
      subtitle = "No plottable net-benefit rows were available."
    )
  }

  list(
    risk = risk_plot,
    auc = auc_plot,
    brier = brier_plot,
    calibration = calibration_plot,
    net_benefit = net_benefit_plot
  )
}

# 🔴 Load: Stage 1 backbone and frozen inputs ===============================
stage1_bundle <- read_stage1_bundle(
  bundle_path = stage1_bundle_file,
  datasets_path = stage1_datasets_file
)

dataset_registry <- extract_stage1_piece(stage1_bundle, "dataset_registry") %>%
  mutate(dataset = normalize_dataset_label(dataset))

formula_registry <- extract_stage1_piece(stage1_bundle, "formula_registry") %>%
  mutate(dataset = normalize_dataset_label(dataset))

horizon_registry <- extract_stage1_piece(stage1_bundle, "horizon_registry") %>%
  mutate(dataset = normalize_dataset_label(dataset))

threshold_registry <- extract_stage1_piece(stage1_bundle, "threshold_registry")
scaling_registry <- extract_stage1_piece(stage1_bundle, "scaling_registry")
metadata_registry <- extract_stage1_piece(stage1_bundle, "metadata_registry")
analysis_datasets <- normalize_stage1_datasets(extract_stage1_datasets(stage1_bundle))

formula_registry_core <- prepare_formula_registry_core(formula_registry)
horizon_registry <- prepare_horizon_registry(horizon_registry)

common_horizons_year <- extract_stage1_horizons(stage1_bundle, horizon_registry)
risk_thresholds <- extract_stage1_thresholds(stage1_bundle, threshold_registry)
prediction_horizons_for_plots <- intersect(common_horizons_year, prediction_horizons_for_plots)
selected_plot_horizon_caption <- paste(prediction_horizons_for_plots, collapse = ", ")

cat("stage1_bundle_file:\n", normalize_existing_path(stage1_bundle_file), "\n\n")
cat("stage1_datasets_file:\n", normalize_existing_path(stage1_datasets_file), "\n\n")

required_datasets <- c("PNU", "SNU", "merged")

if (!all(required_datasets %in% dataset_registry$dataset)) {
  stop("Stage 1 dataset registry must contain PNU, SNU, and merged.", call. = FALSE)
}

if (!setequal(unique(formula_registry_core$dataset), required_datasets)) {
  stop(
    sprintf(
      "Stage 1 formula registry datasets mismatch. observed = [%s]",
      paste(sort(unique(formula_registry_core$dataset)), collapse = ", ")
    ),
    call. = FALSE
  )
}

if (anyDuplicated(formula_registry_core$formula_id) > 0) {
  stop("Stage 1 formula registry contains duplicated `formula_id` values.", call. = FALSE)
}

if (anyNA(formula_registry_core$formula_rhs) || any(trimws(formula_registry_core$formula_rhs) == "")) {
  stop("Stage 1 formula registry contains missing or empty `formula_rhs` values.", call. = FALSE)
}

expected_formula_counts <- c(PNU = 2L, SNU = 2L, merged = 4L)
observed_formula_counts <- table(formula_registry_core$dataset)

if (!all(names(expected_formula_counts) %in% names(observed_formula_counts))) {
  stop("Stage 1 formula registry is missing one or more required dataset blocks.", call. = FALSE)
}

if (!identical(
  as.integer(observed_formula_counts[names(expected_formula_counts)]),
  as.integer(expected_formula_counts)
)) {
  stop(
    sprintf(
      "Stage 1 formula registry distinct formula counts mismatch. observed = [%s]",
      paste(
        sprintf("%s:%s", names(observed_formula_counts), as.integer(observed_formula_counts)),
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

if (anyNA(sort(unique(as.integer(horizon_registry$horizon_year)))) || !identical(sort(unique(as.integer(horizon_registry$horizon_year))), 1:10)) {
  stop("Stage 1 horizon registry is not fixed to 1:10 years.", call. = FALSE)
}

if (!all(risk_thresholds %in% threshold_registry$threshold)) {
  stop("Stage 1 threshold registry and extracted threshold vector are inconsistent.", call. = FALSE)
}

scaling_check <- scaling_registry %>%
  filter(variable == "age_exact_entry", scaled_variable == "age_s")

if (nrow(scaling_check) < 3L) {
  stop("Stage 1 scaling registry must contain dataset-specific age_s scaling rows.", call. = FALSE)
}

cat("formula count by dataset:\n")
print(as.data.frame(table(formula_registry_core$dataset, useNA = "ifany")))
cat("\ncommon horizons:\n")
print(common_horizons_year)
cat("\nrisk thresholds:\n")
print(risk_thresholds)

# 🔴 Plan: non-cure model menu ===============================
family_registry <- tibble(
  model_family = c("exponential", "weibull", "lognormal", "loglogistic", "coxph"),
  fit_engine = c("flexsurv::flexsurvreg", "flexsurv::flexsurvreg", "flexsurv::flexsurvreg", "flexsurv::flexsurvreg", "survival::coxph"),
  flexsurv_dist = c("exp", "weibull", "lnorm", "llogis", NA_character_),
  model_class = "non_cure"
)

model_plan <- tidyr::crossing(
  formula_registry_core,
  family_registry
) %>%
  mutate(
    stage = stage_name,
    stage_role = stage_role,
    stage_spec_version = stage_spec_version,
    model_id = paste(dataset, formula_name, model_family, sep = "__"),
    strategy_family = model_family,
    anchor_horizon_year = vapply(dataset, assign_anchor_horizon, integer(1)),
    tail_extrapolation_capable_flag = make_tail_extrapolation_capable_flag(model_family),
    extrapolation_role = make_extrapolation_role(model_family),
    validation_level = validation_level_default,
    optimism_correction_method = optimism_correction_method_default,
    analysis_time_variable = "days_followup",
    reporting_time_variable = "time_year",
    event_definition = "status_num == 1",
    censoring_definition = "status_num %in% c(0, 2)"
  ) %>%
  arrange(match(dataset, dataset_order), formula_name, model_family)

if (nrow(model_plan) == 0L) {
  stop("No Stage 5 non-cure models were created in the model plan.", call. = FALSE)
}

if (anyDuplicated(model_plan$model_id) > 0) {
  stop("Stage 5 model plan contains duplicated `model_id` values.", call. = FALSE)
}

# 🔴 Build: core results by reuse-or-recompute logic ===============================
existing_core <- load_existing_stage5_core()
fitted_models <- existing_core$fitted_models %||% list()
reuse_existing_core_candidate <- isTRUE(reuse_existing_stage5_core) &&
  !isTRUE(force_recompute_stage5_core) &&
  can_reuse_existing_stage5_core(
    existing_core = existing_core,
    model_plan = model_plan,
    analysis_datasets = analysis_datasets,
    common_horizons_year = common_horizons_year,
    risk_thresholds = risk_thresholds
  )

core_build_mode <- NA_character_

if (reuse_existing_core_candidate) {
  cat("\nReusing existing Stage 5 core tables because the stored results remain structurally compatible.\n")

  model_registry <- finalize_model_registry(
    existing_core$model_registry,
    formula_registry_core = formula_registry_core,
    horizon_registry = horizon_registry
  )

  subject_horizon_risk_long <- finalize_subject_horizon_risk_long(
    existing_core$subject_horizon_risk_long,
    formula_registry_core = formula_registry_core,
    horizon_registry = horizon_registry
  )

  model_performance_long <- finalize_model_performance_long(
    existing_core$model_performance_long,
    formula_registry_core = formula_registry_core,
    horizon_registry = horizon_registry
  )

  calibration_bins_long <- finalize_calibration_bins_long(
    existing_core$calibration_bins_long,
    formula_registry_core = formula_registry_core,
    horizon_registry = horizon_registry
  )

  calibration_diagnostics_long <- finalize_calibration_diagnostics_long(
    existing_core$calibration_diagnostics_long,
    formula_registry_core = formula_registry_core,
    horizon_registry = horizon_registry
  )

  assert_no_dataset_key_suffixes(subject_horizon_risk_long, "subject_horizon_risk_long")
  assert_no_dataset_key_suffixes(model_performance_long, "model_performance_long")
  assert_no_dataset_key_suffixes(calibration_bins_long, "calibration_bins_long")
  assert_no_dataset_key_suffixes(calibration_diagnostics_long, "calibration_diagnostics_long")

  qc_summary <- build_stage5_qc_summary(
    model_plan = model_plan,
    model_registry = model_registry,
    subject_horizon_risk_long = subject_horizon_risk_long,
    model_performance_long = model_performance_long,
    calibration_diagnostics_long = calibration_diagnostics_long,
    calibration_bins_long = calibration_bins_long,
    analysis_datasets = analysis_datasets,
    common_horizons_year = common_horizons_year
  )

  reuse_qc_ok <- TRUE
  tryCatch(
    validate_stage5_outputs(qc_summary),
    error = function(e) {
      reuse_qc_ok <<- FALSE
      message("Existing Stage 5 core tables failed current QC and will be recomputed from scratch: ", conditionMessage(e))
    }
  )

  if (reuse_qc_ok) {
    core_build_mode <- "reused_existing_core"
  }
}

if (is.na(core_build_mode)) {
  cat("\nComputing Stage 5 core results from scratch under the current Stage 1 backbone and revised specification.\n")

  dataset_reference_rows <- list()
  dataset_threshold_reference_rows <- list()
  dataset_ipcw_registry <- list()

  for (dataset_name in dataset_order) {
    dataset_data <- analysis_datasets[[dataset_name]]
    dataset_ipcw <- build_ipcw_cache(dataset_data, common_horizons_year)
    dataset_ipcw_registry[[dataset_name]] <- dataset_ipcw

    dataset_last_event_time_year <- if (sum(dataset_data$event_main, na.rm = TRUE) > 0) {
      max(dataset_data$time_year[dataset_data$event_main == 1L], na.rm = TRUE)
    } else {
      NA_real_
    }
    dataset_last_followup_time_year <- max(dataset_data$time_year, na.rm = TRUE)

    observed_horizon_rows <- list()
    observed_threshold_rows <- list()

    for (horizon_year in common_horizons_year) {
      horizon_key <- paste0("year_", horizon_year)
      ipcw_info <- dataset_ipcw$horizon_cache[[horizon_key]]
      horizon_info <- horizon_label_lookup(dataset_name, horizon_year, horizon_registry)
      horizon_support_info <- make_horizon_support_info(
        ipcw_info = ipcw_info,
        horizon_year = horizon_year,
        last_event_time_year = dataset_last_event_time_year,
        last_followup_time_year = dataset_last_followup_time_year
      )

      observed_meta <- list(
        dataset = dataset_name,
        dataset_key = dataset_name,
        risk_scale = "transition_only_main",
        model_class = "benchmark",
        model_family = "observed_km",
        model_id = paste(dataset_name, "observed_km", sep = "__"),
        strategy_family = "observed_km",
        formula_id = "reference",
        formula_name = "reference",
        formula_label = "Reference",
        formula_scope = "benchmark_reference",
        site_branch = "site_free",
        interaction_branch = "no_age_sex_interaction",
        anchor_horizon_year = assign_anchor_horizon(dataset_name),
        tail_extrapolation_capable_flag = make_tail_extrapolation_capable_flag("observed_km"),
        extrapolation_role = make_extrapolation_role("observed_km"),
        validation_level = validation_level_default,
        optimism_correction_method = optimism_correction_method_default
      )

      observed_metrics <- c(
        observed_km_risk = ipcw_info$observed_km_risk,
        observed_ipcw_risk = ipcw_info$observed_ipcw_risk,
        case_count_ipcw = ipcw_info$case_count_ipcw,
        nonevent_count_ipcw = ipcw_info$nonevent_count_ipcw,
        known_count = ipcw_info$known_count,
        cohort_n = nrow(dataset_data),
        binary_outcome_support_flag = binary_outcome_support_from_ipcw(ipcw_info)
      )

      observed_horizon_rows[[horizon_key]] <- make_metric_rows(
        meta_row = observed_meta,
        horizon_info = horizon_info,
        threshold_value = NA_real_,
        metric_domain = "horizon_reference",
        metrics_named = observed_metrics,
        horizon_support_info = horizon_support_info
      )

      for (threshold_value in risk_thresholds) {
        treat_all_meta <- observed_meta
        treat_all_meta$model_class <- "reference"
        treat_all_meta$model_family <- "treat_all"
        treat_all_meta$model_id <- paste(dataset_name, "treat_all", sep = "__")
        treat_all_meta$strategy_family <- "treat_all"
        treat_all_meta$tail_extrapolation_capable_flag <- make_tail_extrapolation_capable_flag("treat_all")
        treat_all_meta$extrapolation_role <- make_extrapolation_role("treat_all")

        treat_none_meta <- observed_meta
        treat_none_meta$model_class <- "reference"
        treat_none_meta$model_family <- "treat_none"
        treat_none_meta$model_id <- paste(dataset_name, "treat_none", sep = "__")
        treat_none_meta$strategy_family <- "treat_none"
        treat_none_meta$tail_extrapolation_capable_flag <- make_tail_extrapolation_capable_flag("treat_none")
        treat_none_meta$extrapolation_role <- make_extrapolation_role("treat_none")

        treat_all_metrics <- compute_threshold_vector(
          pred_positive = rep(TRUE, nrow(dataset_data)),
          ipcw_info = ipcw_info,
          threshold = threshold_value,
          n_total = nrow(dataset_data)
        )

        treat_none_metrics <- compute_threshold_vector(
          pred_positive = rep(FALSE, nrow(dataset_data)),
          ipcw_info = ipcw_info,
          threshold = threshold_value,
          n_total = nrow(dataset_data)
        )

        observed_threshold_rows[[paste0(horizon_key, "_all_", threshold_value)]] <- make_metric_rows(
          meta_row = treat_all_meta,
          horizon_info = horizon_info,
          threshold_value = threshold_value,
          metric_domain = "threshold_reference",
          metrics_named = treat_all_metrics,
          horizon_support_info = horizon_support_info
        )

        observed_threshold_rows[[paste0(horizon_key, "_none_", threshold_value)]] <- make_metric_rows(
          meta_row = treat_none_meta,
          horizon_info = horizon_info,
          threshold_value = threshold_value,
          metric_domain = "threshold_reference",
          metrics_named = treat_none_metrics,
          horizon_support_info = horizon_support_info
        )
      }
    }

    dataset_reference_rows[[dataset_name]] <- bind_rows_from_df_lists(
      observed_horizon_rows,
      object_name = paste0("dataset_reference_rows__", dataset_name)
    )

    dataset_threshold_reference_rows[[dataset_name]] <- bind_rows_from_df_lists(
      observed_threshold_rows,
      object_name = paste0("dataset_threshold_reference_rows__", dataset_name)
    )
  }

  model_registry_rows <- list()
  risk_long_rows <- list()
  performance_rows <- list()
  calibration_bin_rows <- list()
  calibration_diagnostic_rows <- list()
  fitted_models <- list()

  for (plan_idx in seq_len(nrow(model_plan))) {
    plan_row <- model_plan[plan_idx, , drop = FALSE]

    dataset_name <- plan_row$dataset[[1]]
    dataset_data <- analysis_datasets[[dataset_name]]
    dataset_ipcw <- dataset_ipcw_registry[[dataset_name]]
    formula_rhs <- plan_row$formula_rhs[[1]]
    model_family <- plan_row$model_family[[1]]
    flexsurv_dist <- plan_row$flexsurv_dist[[1]]
    anchor_horizon_year <- plan_row$anchor_horizon_year[[1]]
    anchor_horizon_index <- match(anchor_horizon_year, common_horizons_year)

    dataset_last_event_time_year <- if (sum(dataset_data$event_main, na.rm = TRUE) > 0) {
      max(dataset_data$time_year[dataset_data$event_main == 1L], na.rm = TRUE)
    } else {
      NA_real_
    }
    dataset_last_followup_time_year <- max(dataset_data$time_year, na.rm = TRUE)

    fit_result <- fit_single_noncure_model(
      data = dataset_data,
      formula_rhs = formula_rhs,
      model_family = model_family,
      flexsurv_dist = flexsurv_dist
    )

    fit_stats <- extract_fit_statistics(
      fit_object = fit_result$fit,
      data = dataset_data,
      model_family = model_family,
      convergence_code = fit_result$convergence_code,
      converged = fit_result$converged,
      error_message = fit_result$error_message
    )

    registry_row <- plan_row %>%
      mutate(
        stage = stage_name,
        stage_role = stage_role,
        stage_spec_version = stage_spec_version,
        time_range_score_type = NA_character_,
        time_range_discrimination_upper_time_year = anchor_horizon_year
      ) %>%
      bind_cols(fit_stats)

    model_registry_rows[[plan_row$model_id[[1]]]] <- registry_row

    if (!isTRUE(fit_result$converged) || is.null(fit_result$fit)) {
      next
    }

    pred_result <- if (model_family == "coxph") {
      predict_cox_survival(
        fit = fit_result$fit,
        newdata = dataset_data,
        horizons_year = common_horizons_year
      )
    } else {
      tryCatch(
        predict_parametric_survival(
          fit = fit_result$fit,
          newdata = dataset_data,
          horizons_year = common_horizons_year
        ),
        error = function(e) NULL
      )
    }

    if (is.null(pred_result)) {
      model_registry_rows[[plan_row$model_id[[1]]]] <- model_registry_rows[[plan_row$model_id[[1]]]] %>%
        mutate(
          converged = FALSE,
          error_message = "Prediction extraction failed after successful model fit."
        )
      next
    }

    fitted_models[[plan_row$model_id[[1]]]] <- fit_result$fit

    anchor_risk_score_vec <- pred_result$risk[, anchor_horizon_index]
    time_range_score_info <- choose_time_range_score(
      pred_result = pred_result,
      anchor_horizon_index = anchor_horizon_index,
      anchor_horizon_year = anchor_horizon_year
    )
    time_range_score_vec <- as.numeric(time_range_score_info$score)

    model_registry_rows[[plan_row$model_id[[1]]]] <- model_registry_rows[[plan_row$model_id[[1]]]] %>%
      mutate(
        time_range_score_type = time_range_score_info$score_type,
        time_range_discrimination_upper_time_year = anchor_horizon_year
      )

    prediction_long <- tibble(
      stage = stage_name,
      stage_role = stage_role,
      stage_spec_version = stage_spec_version,
      validation_level = validation_level_default,
      optimism_correction_method = optimism_correction_method_default,
      dataset = dataset_name,
      dataset_key = dataset_name,
      risk_scale = "transition_only_main",
      model_class = "non_cure",
      model_family = model_family,
      model_id = plan_row$model_id[[1]],
      strategy_family = plan_row$strategy_family[[1]],
      formula_id = plan_row$formula_id[[1]],
      formula_name = plan_row$formula_name[[1]],
      formula_label = plan_row$formula_label[[1]],
      formula_scope = plan_row$formula_scope[[1]],
      site_branch = plan_row$site_branch[[1]],
      interaction_branch = plan_row$interaction_branch[[1]],
      tail_extrapolation_capable_flag = plan_row$tail_extrapolation_capable_flag[[1]],
      extrapolation_role = plan_row$extrapolation_role[[1]],
      unique_person_id = rep(dataset_data$unique_person_id, each = length(common_horizons_year)),
      site = rep(as.character(dataset_data$site), each = length(common_horizons_year)),
      id = rep(dataset_data$id, each = length(common_horizons_year)),
      sex_num = rep(dataset_data$sex_num, each = length(common_horizons_year)),
      age_exact_entry = rep(dataset_data$age_exact_entry, each = length(common_horizons_year)),
      age_s = rep(dataset_data$age_s, each = length(common_horizons_year)),
      horizon_year = rep(common_horizons_year, times = nrow(dataset_data)),
      survival_pred = as.vector(t(pred_result$survival)),
      risk_pred = as.vector(t(pred_result$risk)),
      median_survival_pred = rep(pred_result$median_survival, each = length(common_horizons_year)),
      anchor_horizon_year = anchor_horizon_year,
      anchor_risk_score = rep(anchor_risk_score_vec, each = length(common_horizons_year)),
      c_index_score = rep(time_range_score_vec, each = length(common_horizons_year)),
      c_index_score_type = time_range_score_info$score_type
    )

    risk_long_rows[[plan_row$model_id[[1]]]] <- prediction_long

    cindex_meta <- list(
      dataset = dataset_name,
      dataset_key = dataset_name,
      risk_scale = "transition_only_main",
      model_class = "non_cure",
      model_family = model_family,
      model_id = plan_row$model_id[[1]],
      strategy_family = plan_row$strategy_family[[1]],
      formula_id = plan_row$formula_id[[1]],
      formula_name = plan_row$formula_name[[1]],
      formula_label = plan_row$formula_label[[1]],
      formula_scope = plan_row$formula_scope[[1]],
      site_branch = plan_row$site_branch[[1]],
      interaction_branch = plan_row$interaction_branch[[1]],
      anchor_horizon_year = anchor_horizon_year,
      tail_extrapolation_capable_flag = plan_row$tail_extrapolation_capable_flag[[1]],
      extrapolation_role = plan_row$extrapolation_role[[1]],
      validation_level = validation_level_default,
      optimism_correction_method = optimism_correction_method_default
    )

    harrell_c <- compute_harrell_c(
      time_year = dataset_data$time_year,
      event_main = dataset_data$event_main,
      risk_score = time_range_score_vec,
      ymax = anchor_horizon_year
    )

    uno_c <- compute_uno_c(
      time_year = dataset_data$time_year,
      event_main = dataset_data$event_main,
      risk_score = time_range_score_vec,
      ymax = anchor_horizon_year
    )

    cindex_horizon_info <- horizon_label_lookup(dataset_name, anchor_horizon_year, horizon_registry)

    performance_rows[[paste0(plan_row$model_id[[1]], "_cindex")]] <- make_metric_rows(
      meta_row = cindex_meta,
      horizon_info = cindex_horizon_info,
      threshold_value = NA_real_,
      metric_domain = "time_range_discrimination",
      metrics_named = c(
        harrell_c = harrell_c$estimate,
        harrell_c_se = harrell_c$se,
        uno_c = uno_c$estimate,
        uno_c_se = uno_c$se,
        time_range_upper_year = anchor_horizon_year
      )
    )

    brier_by_horizon <- rep(NA_real_, length(common_horizons_year))

    for (h_idx in seq_along(common_horizons_year)) {
      horizon_year <- common_horizons_year[h_idx]
      horizon_key <- paste0("year_", horizon_year)

      horizon_info <- horizon_label_lookup(dataset_name, horizon_year, horizon_registry)
      ipcw_info <- dataset_ipcw$horizon_cache[[horizon_key]]
      horizon_support_info <- make_horizon_support_info(
        ipcw_info = ipcw_info,
        horizon_year = horizon_year,
        last_event_time_year = dataset_last_event_time_year,
        last_followup_time_year = dataset_last_followup_time_year
      )
      pred_risk_h <- pred_result$risk[, h_idx]
      binary_support_flag <- isTRUE(horizon_support_info$binary_outcome_support_flag == 1L)

      auc_est <- compute_weighted_auc_cd(
        score = pred_risk_h,
        w_case = ipcw_info$w_case,
        w_control = ipcw_info$w_control
      )

      y_known <- ipcw_info$y_horizon
      w_known <- ipcw_info$w_ipcw
      known_index <- !is.na(y_known) & !is.na(w_known) & is.finite(w_known) & w_known > 0

      if (binary_support_flag && any(known_index)) {
        brier_score <- sum(
          w_known[known_index] * (y_known[known_index] - pred_risk_h[known_index])^2,
          na.rm = TRUE
        ) / nrow(dataset_data)

        brier_null <- sum(
          w_known[known_index] * (y_known[known_index] - ipcw_info$observed_km_risk)^2,
          na.rm = TRUE
        ) / nrow(dataset_data)

        scaled_brier <- if (!is.na(brier_null) && brier_null > 0) {
          1 - (brier_score / brier_null)
        } else {
          NA_real_
        }
      } else {
        brier_score <- NA_real_
        scaled_brier <- NA_real_
      }

      brier_by_horizon[h_idx] <- brier_score

      calibration_stats <- compute_calibration_statistics(
        pred_risk = pred_risk_h,
        ipcw_info = ipcw_info,
        observed_reference_risk = ipcw_info$observed_km_risk
      )

      horizon_metric_values <- c(
        mean_predicted_risk = calibration_stats$mean_predicted_risk,
        observed_km_risk = ipcw_info$observed_km_risk,
        observed_ipcw_risk = ipcw_info$observed_ipcw_risk,
        mean_calibration_difference = calibration_stats$mean_calibration_difference,
        calibration_intercept_offset = calibration_stats$calibration_intercept_offset,
        calibration_intercept_free = calibration_stats$calibration_intercept_free,
        calibration_slope = calibration_stats$calibration_slope,
        calibration_known_n = calibration_stats$calibration_known_n,
        calibration_case_weight = calibration_stats$calibration_case_weight,
        calibration_nonevent_weight = calibration_stats$calibration_nonevent_weight,
        binary_outcome_support_flag = calibration_stats$binary_outcome_support_flag,
        calibration_minimum_n_flag = calibration_stats$calibration_minimum_n_flag,
        calibration_support_flag = calibration_stats$calibration_support_flag,
        calibration_fit_success_flag = calibration_stats$calibration_fit_success_flag,
        calibration_offset_converged_flag = calibration_stats$calibration_offset_converged_flag,
        calibration_slope_converged_flag = calibration_stats$calibration_slope_converged_flag,
        calibration_offset_boundary_flag = calibration_stats$calibration_offset_boundary_flag,
        calibration_slope_boundary_flag = calibration_stats$calibration_slope_boundary_flag,
        calibration_regression_suppressed_flag = calibration_stats$calibration_regression_suppressed_flag,
        time_dependent_auc = if (binary_support_flag) auc_est else NA_real_,
        brier_score = brier_score,
        scaled_brier_score = scaled_brier,
        case_count_ipcw = ipcw_info$case_count_ipcw,
        nonevent_count_ipcw = ipcw_info$nonevent_count_ipcw,
        known_count = ipcw_info$known_count
      )

      performance_rows[[paste0(plan_row$model_id[[1]], "_horizon_", horizon_year)]] <- make_metric_rows(
        meta_row = cindex_meta,
        horizon_info = horizon_info,
        threshold_value = NA_real_,
        metric_domain = "horizon_summary",
        metrics_named = horizon_metric_values,
        horizon_support_info = horizon_support_info
      )

      calibration_diagnostic_rows[[paste0(plan_row$model_id[[1]], "_caldiag_", horizon_year)]] <- make_calibration_diagnostic_row(
        meta_row = cindex_meta,
        horizon_info = horizon_info,
        calibration_stats = calibration_stats,
        horizon_support_info = horizon_support_info
      )

      calibration_bins <- compute_calibration_bins(
        pred_risk = pred_risk_h,
        ipcw_info = ipcw_info,
        group_count = calibration_group_count
      ) %>%
        mutate(
          stage = stage_name,
          stage_role = stage_role,
          stage_spec_version = stage_spec_version,
          validation_level = validation_level_default,
          optimism_correction_method = optimism_correction_method_default,
          dataset = dataset_name,
          dataset_key = dataset_name,
          risk_scale = "transition_only_main",
          model_class = "non_cure",
          model_family = model_family,
          model_id = plan_row$model_id[[1]],
          strategy_family = plan_row$strategy_family[[1]],
          formula_id = plan_row$formula_id[[1]],
          formula_name = plan_row$formula_name[[1]],
          formula_label = plan_row$formula_label[[1]],
          formula_scope = plan_row$formula_scope[[1]],
          site_branch = plan_row$site_branch[[1]],
          interaction_branch = plan_row$interaction_branch[[1]],
          anchor_horizon_year = anchor_horizon_year,
          tail_extrapolation_capable_flag = plan_row$tail_extrapolation_capable_flag[[1]],
          extrapolation_role = plan_row$extrapolation_role[[1]],
          horizon_year = horizon_year,
          horizon_id = horizon_info$horizon_id,
          horizon_label = horizon_info$horizon_label,
          support_tier = horizon_info$support_tier,
          support_tier_standard = horizon_info$support_tier_standard,
          interpretation_tier = horizon_info$interpretation_tier,
          primary_supported_flag = horizon_info$primary_supported_flag,
          horizon_evidence_class = horizon_info$horizon_evidence_class,
          claim_restriction_flag = horizon_info$claim_restriction_flag,
          interpretation_note = horizon_info$interpretation_note,
          observed_case_support_flag = horizon_support_info$case_support_flag,
          observed_nonevent_support_flag = horizon_support_info$nonevent_support_flag,
          observed_horizon_information_flag = horizon_support_info$observed_horizon_information_flag,
          partial_binary_outcome_support_flag = horizon_support_info$partial_binary_outcome_support_flag,
          beyond_last_event_flag = horizon_support_info$beyond_last_event_flag,
          beyond_last_followup_flag = horizon_support_info$beyond_last_followup_flag,
          descriptive_projection_flag = horizon_support_info$descriptive_projection_flag,
          known_count_at_horizon = horizon_support_info$known_count_at_horizon,
          low_known_count_warning_flag = horizon_support_info$low_known_count_warning_flag,
          late_horizon_instability_flag = horizon_support_info$late_horizon_instability_flag,
          late_horizon_warning_reason = horizon_support_info$late_horizon_warning_reason,
          horizon_data_status = horizon_support_info$horizon_data_status,
          descriptive_value_status = horizon_support_info$descriptive_value_status
        )

      calibration_bin_rows[[paste0(plan_row$model_id[[1]], "_calbin_", horizon_year)]] <- calibration_bins

      for (threshold_value in risk_thresholds) {
        threshold_metrics <- compute_threshold_vector(
          pred_positive = pred_risk_h >= threshold_value,
          ipcw_info = ipcw_info,
          threshold = threshold_value,
          n_total = nrow(dataset_data)
        )

        performance_rows[[paste0(plan_row$model_id[[1]], "_threshold_", horizon_year, "_", threshold_value)]] <- make_metric_rows(
          meta_row = cindex_meta,
          horizon_info = horizon_info,
          threshold_value = threshold_value,
          metric_domain = "threshold_summary",
          metrics_named = threshold_metrics,
          horizon_support_info = horizon_support_info
        )
      }
    }

    ibs_valid <- !is.na(brier_by_horizon)

    ibs_value <- if (sum(ibs_valid) >= 2) {
      valid_horizons <- common_horizons_year[ibs_valid]
      valid_brier <- brier_by_horizon[ibs_valid]

      sum(
        diff(valid_horizons) * (head(valid_brier, -1) + tail(valid_brier, -1)) / 2
      ) / (max(valid_horizons) - min(valid_horizons))
    } else {
      NA_real_
    }

    ibs_horizon_info <- horizon_label_lookup(dataset_name, anchor_horizon_year, horizon_registry)
    supported_horizons_for_ibs <- common_horizons_year[ibs_valid]

    performance_rows[[paste0(plan_row$model_id[[1]], "_ibs")]] <- make_metric_rows(
      meta_row = cindex_meta,
      horizon_info = ibs_horizon_info,
      threshold_value = NA_real_,
      metric_domain = "integrated_summary",
      metrics_named = c(
        restricted_integrated_brier_score = ibs_value,
        restricted_ibs_supported_horizon_count = sum(ibs_valid),
        restricted_ibs_min_supported_horizon_year = if (length(supported_horizons_for_ibs) > 0) min(supported_horizons_for_ibs) else NA_real_,
        restricted_ibs_max_supported_horizon_year = if (length(supported_horizons_for_ibs) > 0) max(supported_horizons_for_ibs) else NA_real_,
        restricted_ibs_full_grid_flag = as.integer(length(supported_horizons_for_ibs) == length(common_horizons_year))
      )
    )
  }

  model_registry <- bind_rows_from_df_lists(
    model_registry_rows,
    object_name = "model_registry_rows"
  )

  subject_horizon_risk_long <- bind_rows_from_df_lists(
    risk_long_rows,
    object_name = "risk_long_rows"
  )

  all_performance_rows <- c(
    dataset_reference_rows,
    dataset_threshold_reference_rows,
    performance_rows
  )

  all_performance_rows <- purrr::compact(all_performance_rows)

  if (!all(vapply(all_performance_rows, inherits, logical(1), what = "data.frame"))) {
    stop("`all_performance_rows` contains at least one non-data.frame object.", call. = FALSE)
  }

  model_performance_long <- bind_rows(all_performance_rows)

  calibration_bins_long <- bind_rows_from_df_lists(
    calibration_bin_rows,
    object_name = "calibration_bin_rows"
  )

  calibration_diagnostics_long <- bind_rows_from_df_lists(
    calibration_diagnostic_rows,
    object_name = "calibration_diagnostic_rows"
  )

  model_registry <- finalize_model_registry(
    model_registry,
    formula_registry_core = formula_registry_core,
    horizon_registry = horizon_registry
  )

  subject_horizon_risk_long <- finalize_subject_horizon_risk_long(
    subject_horizon_risk_long,
    formula_registry_core = formula_registry_core,
    horizon_registry = horizon_registry
  )

  model_performance_long <- finalize_model_performance_long(
    model_performance_long,
    formula_registry_core = formula_registry_core,
    horizon_registry = horizon_registry
  )

  calibration_bins_long <- finalize_calibration_bins_long(
    calibration_bins_long,
    formula_registry_core = formula_registry_core,
    horizon_registry = horizon_registry
  )

  calibration_diagnostics_long <- finalize_calibration_diagnostics_long(
    calibration_diagnostics_long,
    formula_registry_core = formula_registry_core,
    horizon_registry = horizon_registry
  )

  assert_no_dataset_key_suffixes(subject_horizon_risk_long, "subject_horizon_risk_long")
  assert_no_dataset_key_suffixes(model_performance_long, "model_performance_long")
  assert_no_dataset_key_suffixes(calibration_bins_long, "calibration_bins_long")
  assert_no_dataset_key_suffixes(calibration_diagnostics_long, "calibration_diagnostics_long")

  qc_summary <- build_stage5_qc_summary(
    model_plan = model_plan,
    model_registry = model_registry,
    subject_horizon_risk_long = subject_horizon_risk_long,
    model_performance_long = model_performance_long,
    calibration_diagnostics_long = calibration_diagnostics_long,
    calibration_bins_long = calibration_bins_long,
    analysis_datasets = analysis_datasets,
    common_horizons_year = common_horizons_year
  )

  validate_stage5_outputs(qc_summary)

  core_build_mode <- "recomputed_current_spec"
}

# 🔴 Derive: summary tables and plot sources ===============================
family_fit_summary <- model_registry %>%
  transmute(
    stage,
    stage_role,
    stage_spec_version,
    dataset,
    dataset_key,
    risk_scale,
    model_class,
    model_family,
    model_id,
    formula_id,
    formula_name,
    formula_label,
    formula_scope,
    site_branch,
    interaction_branch,
    fit_engine,
    converged,
    convergence_code,
    error_message,
    n_obs,
    n_event,
    effective_parameter_count,
    likelihood_basis,
    information_criterion_comparable_across_families_flag,
    information_criterion_note,
    logLik,
    AIC,
    BIC,
    partial_logLik,
    partial_AIC,
    partial_BIC,
    last_event_time_year,
    last_followup_time_year,
    anchor_horizon_year,
    time_range_support_tier,
    time_range_horizon_evidence_class,
    time_range_claim_restriction_flag,
    extrapolation_rule,
    extrapolation_role,
    tail_extrapolation_capable_flag
  ) %>%
  arrange(match(dataset, dataset_order), formula_name, model_family)

threshold_classification_long <- model_performance_long %>%
  filter(metric_domain %in% c("threshold_summary", "threshold_reference")) %>%
  arrange(match(dataset, dataset_order), formula_name, model_family, horizon_year, threshold, metric_name)

false_positive_burden_long <- threshold_classification_long %>%
  filter(metric_name %in% c(
    "binary_outcome_support_flag",
    "case_count_ipcw",
    "nonevent_count_ipcw",
    "known_count",
    "positive_classification_rate",
    "false_positive_count_ipcw",
    "false_positive_burden_all",
    "false_positive_burden_nonevents",
    "false_positive_per_100",
    "ppv",
    "false_discovery_proportion"
  )) %>%
  arrange(match(dataset, dataset_order), formula_name, model_family, horizon_year, threshold, metric_name)

decision_curve_long <- threshold_classification_long %>%
  filter(metric_name %in% c(
    "binary_outcome_support_flag",
    "event_prevalence_ipcw",
    "event_rate_per_subject_ipcw",
    "net_benefit",
    "net_benefit_treat_all",
    "net_reduction_unnecessary_per_100"
  )) %>%
  arrange(match(dataset, dataset_order), formula_name, model_family, horizon_year, threshold, metric_name)

plot_sources <- build_stage5_plot_sources(
  subject_horizon_risk_long = subject_horizon_risk_long,
  model_performance_long = model_performance_long,
  calibration_bins_long = calibration_bins_long,
  common_horizons_year = common_horizons_year,
  prediction_horizons_for_plots = prediction_horizons_for_plots
)

plot_objects <- build_stage5_plot_objects(
  plot_sources = plot_sources,
  common_horizons_year = common_horizons_year,
  prediction_horizons_for_plots = prediction_horizons_for_plots
)

# 🔴 Export: tables, plots, and reusable bundles ===============================
readr::write_csv(model_registry, model_registry_file)
readr::write_csv(subject_horizon_risk_long, risk_long_file)
readr::write_csv(model_performance_long, performance_long_file)
readr::write_csv(calibration_bins_long, calibration_bins_file)
readr::write_csv(calibration_diagnostics_long, calibration_diagnostics_file)
readr::write_csv(qc_summary, qc_summary_file)
readr::write_csv(family_fit_summary, family_fit_summary_file)
readr::write_csv(threshold_classification_long, threshold_classification_file)
readr::write_csv(false_positive_burden_long, false_positive_burden_file)
readr::write_csv(decision_curve_long, decision_curve_file)

readr::write_csv(plot_sources$risk, plot_source_risk_file)
readr::write_csv(plot_sources$auc, plot_source_auc_file)
readr::write_csv(plot_sources$brier, plot_source_brier_file)
readr::write_csv(plot_sources$calibration, plot_source_calibration_file)
readr::write_csv(plot_sources$net_benefit, plot_source_net_benefit_file)

save_plot_pdf_bundle(
  plot_list = list(
    plot_objects$risk,
    plot_objects$auc,
    plot_objects$brier,
    plot_objects$calibration,
    plot_objects$net_benefit
  ),
  pdf_path = pdf_file
)

save_plot_png(plot_objects$risk, plot_file_png_risk)
save_plot_png(plot_objects$auc, plot_file_png_auc)
save_plot_png(plot_objects$brier, plot_file_png_brier)
save_plot_png(plot_objects$calibration, plot_file_png_calibration)
save_plot_png(plot_objects$net_benefit, plot_file_png_net_benefit)

if (!is.null(fitted_models)) {
  saveRDS(fitted_models, fitted_models_file)
}

package_versions_out <- c(
  vapply(
    required_packages,
    function(pkg) as.character(utils::packageVersion(pkg)),
    character(1)
  ),
  survAUC = NA_character_
)

stage5_results_bundle <- list(
  stage = stage_name,
  stage_role = stage_role,
  stage_spec_version = stage_spec_version,
  created_at = as.character(Sys.time()),
  data_path = normalize_existing_path(data_path),
  export_path = normalize_existing_path(export_path),
  stage1_bundle_file = normalize_existing_path(stage1_bundle_file),
  stage1_datasets_file = normalize_existing_path(stage1_datasets_file),
  validation_level = validation_level_default,
  optimism_correction_method = optimism_correction_method_default,
  calibration_target = calibration_target_label,
  mean_calibration_reference = mean_calibration_reference_label,
  common_horizons_year = common_horizons_year,
  prediction_horizons_for_plots = prediction_horizons_for_plots,
  risk_thresholds = risk_thresholds,
  calibration_group_count = calibration_group_count,
  late_horizon_warning_start_year = late_horizon_warning_start_year,
  late_horizon_known_count_warning_min = late_horizon_known_count_warning_min,
  has_survAUC = has_survAUC,
  uno_c_backend = "survival::concordance(timewt = 'n/G2')",
  time_range_discrimination_rule = "restricted_to_anchor_horizon_year_using_anchor_risk_score",
  risk_scale = "transition_only_main",
  core_build_mode = core_build_mode,
  reuse_existing_stage5_core = reuse_existing_stage5_core,
  force_recompute_stage5_core = force_recompute_stage5_core,
  force_regenerate_plots = force_regenerate_plots,
  package_versions = package_versions_out,
  model_registry = model_registry,
  subject_horizon_risk_long = subject_horizon_risk_long,
  model_performance_long = model_performance_long,
  calibration_bins_long = calibration_bins_long,
  calibration_diagnostics_long = calibration_diagnostics_long,
  qc_summary = qc_summary,
  family_fit_summary = family_fit_summary,
  threshold_classification_long = threshold_classification_long,
  false_positive_burden_long = false_positive_burden_long,
  decision_curve_long = decision_curve_long,
  plot_sources = plot_sources
)

saveRDS(stage5_results_bundle, results_bundle_file)

export_manifest <- tibble(
  file_name = c(
    "stage5_model_registry.csv",
    "stage5_subject_horizon_risk_long.csv",
    "stage5_model_performance_long.csv",
    "stage5_calibration_bins_long.csv",
    "stage5_calibration_diagnostics_long.csv",
    "stage5_qc_summary.csv",
    "stage5_family_fit_summary.csv",
    "stage5_threshold_classification_long.csv",
    "stage5_false_positive_burden_long.csv",
    "stage5_decision_curve_long.csv",
    "stage5_plot_source_risk_summary.csv",
    "stage5_plot_source_auc_summary.csv",
    "stage5_plot_source_brier_summary.csv",
    "stage5_plot_source_calibration_selected_horizons.csv",
    "stage5_plot_source_net_benefit_selected_horizons.csv",
    "stage5_summary_plots.pdf",
    "stage5_fitted_models.rds",
    "stage5_results_bundle.rds",
    "stage5_export_manifest.csv",
    "stage5_plot_risk_summary.png",
    "stage5_plot_auc_summary.png",
    "stage5_plot_brier_summary.png",
    "stage5_plot_calibration_selected_horizons.png",
    "stage5_plot_net_benefit_selected_horizons.png"
  ),
  object_name = c(
    "model_registry",
    "subject_horizon_risk_long",
    "model_performance_long",
    "calibration_bins_long",
    "calibration_diagnostics_long",
    "qc_summary",
    "family_fit_summary",
    "threshold_classification_long",
    "false_positive_burden_long",
    "decision_curve_long",
    "plot_source_risk_summary",
    "plot_source_auc_summary",
    "plot_source_brier_summary",
    "plot_source_calibration_selected_horizons",
    "plot_source_net_benefit_selected_horizons",
    "summary_plots_pdf",
    "fitted_models",
    "stage5_results_bundle",
    "export_manifest",
    "plot_risk_png",
    "plot_auc_png",
    "plot_brier_png",
    "plot_calibration_png",
    "plot_net_benefit_png"
  ),
  description = c(
    "Model-level fit registry for all Stage 5 non-cure models, including current Stage 1 backbone metadata, effective parameter count, and anchor-horizon support labels.",
    "Subject-by-horizon predicted survival and risk long table, aligned to the Stage 1 formula registry and current horizon-support metadata.",
    "Long-format performance table containing horizon summaries, support-gated threshold summaries, benchmark references, and restricted IBS summaries.",
    "Grouped calibration source-of-truth table with support-tier, evidence-class, and claim-restriction metadata for every plotted bin.",
    "Model-by-horizon calibration fit diagnostics including support, projection status, convergence, boundary, and suppression flags.",
    "Automated Stage 5 QC summary for code-output consistency checks under the revised specification.",
    "Family-level fit summary table with convergence status, effective parameter count, log-likelihood, AIC/BIC, and anchor-horizon support labels.",
    "Threshold-based classification source-of-truth table for all non-cure models and treat-all/treat-none references.",
    "False-positive-burden view of the threshold table, exported separately for direct interpretation and later reporting.",
    "Decision-curve source-of-truth table including net benefit and net reduction in unnecessary high-risk classification.",
    "Risk-summary plot source table used for both CSV export and the risk PNG/PDF figures.",
    "AUC plot source table used for both CSV export and the AUC PNG/PDF figures.",
    "Brier-score plot source table used for both CSV export and the Brier PNG/PDF figures.",
    "Grouped calibration plot source table used for both CSV export and the calibration PNG/PDF figures.",
    "Net-benefit plot source table used for both CSV export and the decision-curve PNG/PDF figures.",
    "Combined multi-page PDF containing the five main Stage 5 figures.",
    "Reusable fitted non-cure model objects as an RDS list.",
    "Reusable Stage 5 results bundle with exported tables, QC summary, and plot-source objects.",
    "Manifest of all Stage 5 exported files, including PNG figures.",
    "Standalone PNG for the risk-summary figure.",
    "Standalone PNG for the AUC figure.",
    "Standalone PNG for the Brier-score figure.",
    "Standalone PNG for the grouped calibration figure.",
    "Standalone PNG for the decision-curve net-benefit figure."
  ),
  file_path = c(
    model_registry_file,
    risk_long_file,
    performance_long_file,
    calibration_bins_file,
    calibration_diagnostics_file,
    qc_summary_file,
    family_fit_summary_file,
    threshold_classification_file,
    false_positive_burden_file,
    decision_curve_file,
    plot_source_risk_file,
    plot_source_auc_file,
    plot_source_brier_file,
    plot_source_calibration_file,
    plot_source_net_benefit_file,
    pdf_file,
    fitted_models_file,
    results_bundle_file,
    export_manifest_file,
    plot_file_png_risk,
    plot_file_png_auc,
    plot_file_png_brier,
    plot_file_png_calibration,
    plot_file_png_net_benefit
  )
)

readr::write_csv(export_manifest, export_manifest_file)

cat("\nStage 5 exports completed:\n")
print(export_manifest)
cat("\nCore build mode:\n")
print(core_build_mode)
