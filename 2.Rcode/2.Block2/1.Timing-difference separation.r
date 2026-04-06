# 🔴 Configure: paths and Block 2 constants ===============================
dropbox_project_root <- Sys.getenv(
  "CHR_COHORTS_DROPBOX_ROOT",
  unset = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model"
)
data_path <- Sys.getenv(
  "BLOCK2_STAGE1_DATA_PATH",
  unset = file.path(dropbox_project_root, "old", "stage1_Backbone lock")
)
export_path <- Sys.getenv(
  "BLOCK2_EXPORT_PATH",
  unset = file.path(dropbox_project_root, "2.Block2")
)
supporting_output_dir <- file.path(export_path, "sub_supporting")
legacy_stage4_export_path <- Sys.getenv(
  "BLOCK2_LEGACY_STAGE4_EXPORT_PATH",
  unset = file.path(dropbox_project_root, "old", "stage4_Timing-difference separation")
)

stage1_analysis_datasets_file <- file.path(data_path, "stage1_analysis_datasets.rds")
stage1_backbone_bundle_file <- file.path(data_path, "stage1_backbone_bundle.rds")
stage1_dataset_registry_file <- file.path(data_path, "stage1_dataset_registry.csv")
stage1_scaling_registry_file <- file.path(data_path, "stage1_scaling_registry.csv")
stage1_metadata_registry_file <- file.path(data_path, "stage1_metadata_registry.csv")
stage1_formula_registry_file <- file.path(data_path, "stage1_formula_registry.csv")
stage1_horizon_registry_file <- file.path(data_path, "stage1_horizon_registry.csv")
stage1_threshold_registry_file <- file.path(data_path, "stage1_threshold_registry.csv")

bootstrap_iterations <- 1000L
bootstrap_seed <- 20260322L
bootstrap_min_success_rate <- 0.90

time_zero_epsilon <- 1e-08
piecewise_cuts_year <- c(1, 2, 5)
hazard_band_breaks_year <- c(0, 1, 2, 5, 10)
supported_short_horizons_year <- c(1, 2)
intermediate_interval_label <- "2-5"
tail_interval_label <- ">5"

piecewise_min_events_per_site <- 1L
piecewise_min_subjects_per_site <- 5L

pnu_site_label <- "PNU"
snu_site_label <- "SNU"

block2_spec_version <- "block2_timing_difference_v1_20260402"
block2_caching_policy <- "prefer_block2_cache_then_reuse_compatible_legacy_stage4_outputs_before_recomputation"

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
dir.create(supporting_output_dir, recursive = TRUE, showWarnings = FALSE)

dataset_order <- c("PNU", "SNU", "merged")
png_files_do_not_count_toward_export_limit <- TRUE

# 🔴 Define: reusable Block 2 helpers ===============================
## 🟠 Define: scalar and path utilities ===============================
normalize_existing_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Required CSV file does not exist: %s", path), call. = FALSE)
  }
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

safe_read_optional_csv <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

safe_read_optional_rds <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  readRDS(path)
}

resolve_existing_output_file <- function(primary_dir, primary_file_name, secondary_dir = NULL, secondary_file_name = primary_file_name, legacy_dir = NULL, legacy_file_name = primary_file_name) {
  primary_path <- file.path(primary_dir, primary_file_name)
  if (file.exists(primary_path)) {
    return(primary_path)
  }

  if (!is.null(secondary_dir) && nzchar(secondary_dir)) {
    secondary_path <- file.path(secondary_dir, secondary_file_name)
    if (file.exists(secondary_path)) {
      return(secondary_path)
    }
  }

  if (!is.null(legacy_dir) && nzchar(legacy_dir)) {
    legacy_path <- file.path(legacy_dir, legacy_file_name)
    if (file.exists(legacy_path)) {
      return(legacy_path)
    }
  }

  primary_path
}

safe_divide <- function(num, den) {
  ifelse(is.na(den) | den <= 0, NA_real_, num / den)
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
    "dataset",
    "dataset_label",
    "source_dataset",
    "dataset_scope",
    "horizon_year",
    "time_year",
    "interval_label",
    "result_type",
    "site_counterfactual",
    "contrast_label",
    "analysis_group",
    "reporting_priority",
    "support_tier",
    "support_display_label"
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

make_relative_export_reference <- function(path) {
  normalized_path <- normalize_existing_path(path)
  normalized_export_path <- normalize_existing_path(export_path)
  export_prefix <- paste0(normalized_export_path, "/")

  if (startsWith(normalized_path, export_prefix)) {
    return(substring(normalized_path, nchar(export_prefix) + 1L))
  }

  make_portable_file_reference(normalized_path)
}

format_manifest_bullets <- function(df) {
  if (nrow(df) == 0L) {
    return("- none")
  }

  vapply(
    seq_len(nrow(df)),
    function(i) sprintf("- `%s`: %s", df$file_name[i], df$description[i]),
    character(1)
  )
}

write_block2_readme <- function(readme_path, manifest_df, metadata_df) {
  top_level_df <- manifest_df %>%
    filter(!startsWith(file_name, "sub_supporting/"))
  navigation_df <- top_level_df %>%
    filter(file_name == "README.md")
  interpretation_df <- top_level_df %>%
    filter(
      file_name != "README.md",
      !startsWith(file_name, "stage1_")
    )
  audit_df <- top_level_df %>%
    filter(startsWith(file_name, "stage1_"))
  supporting_df <- manifest_df %>%
    filter(startsWith(file_name, "sub_supporting/"))

  metadata_lookup <- setNames(
    metadata_df$metadata_value,
    paste(metadata_df$metadata_group, metadata_df$metadata_name, sep = "::")
  )

  lines <- c(
    "# Block 2 Output README",
    "",
    "## Overview",
    "This folder contains Block 2 timing-difference separation outputs.",
    "",
    sprintf("- Main specification: `%s`", metadata_lookup[["stage::canonical_framework"]]),
    sprintf("- Export root: `%s`", metadata_lookup[["inputs::export_path"]]),
    sprintf("- Supporting subfolder: `%s`", metadata_lookup[["inputs::supporting_output_dir"]]),
    "- Main event definition: `status_num == 1`",
    "- Main censoring definition: `status_num %in% c(0, 2)`",
    "- Time scale: `time_year = days_followup / 365.25`",
    "",
    "Top-level files are intended for interpretation and audit. The `sub_supporting/` folder contains supporting tables, logs, metadata, model objects, and standalone PNG files.",
    "",
    "## Recommended Reading Order",
    "1. `README.md`",
    "2. `block2_plot_book.pdf`",
    "3. `block2_supported_short_horizon_contrast.csv`",
    "4. `block2_intermediate_horizon_contrast.csv`",
    "5. `block2_tail_diagnostic_contrast.csv`",
    "6. `stage1_horizon_registry.csv`",
    "7. `stage1_dataset_registry.csv`",
    "8. `stage1_analysis_datasets.rds`",
    "9. `sub_supporting/block2_export_manifest.csv`",
    "",
    "## Top-Level Navigation File",
    format_manifest_bullets(navigation_df),
    "",
    "## Top-Level Interpretation Files",
    format_manifest_bullets(interpretation_df),
    "",
    "## Top-Level Audit Files",
    format_manifest_bullets(audit_df),
    "",
    "## Supporting Files In `sub_supporting/`",
    format_manifest_bullets(supporting_df),
    "",
    "## Notes",
    "- `sub_supporting/block2_export_manifest.csv` is the machine-readable registry for the full export bundle.",
    "- `sub_supporting/block2_metadata_registry.csv` records the specification, paths, and output rules used for this export.",
    "- `2-5` and `>5` contrasts can appear as suppressed when interval support is inadequate on the PNU side; this is expected under the Block 2 support rules.",
    "- Standalone PNG files are supporting copies of the figures already collected in `block2_plot_book.pdf`."
  )

  writeLines(lines, con = readme_path)
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

## 🟠 Define: build-log helpers ===============================
new_build_log <- function() {
  tibble(
    component = character(),
    status = character(),
    source_file = character(),
    detail = character()
  )
}

append_build_log <- function(build_log, component, status, source_file, detail) {
  bind_rows(
    build_log,
    tibble(
      component = component,
      status = status,
      source_file = source_file,
      detail = detail
    )
  )
}

## 🟠 Define: Stage 1 readers and harmonizers ===============================
read_stage1_inputs <- function() {
  backbone_bundle <- safe_read_optional_rds(stage1_backbone_bundle_file)
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

normalize_interpretation_tier <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x == "primary-supported" ~ "primary_supported",
    TRUE ~ as.character(x)
  )
}

make_support_priority <- function(support_tier) {
  dplyr::case_when(
    is.na(support_tier) ~ NA_character_,
    support_tier == "primary_supported" ~ "primary",
    support_tier == "sensitivity" ~ "sensitivity",
    support_tier == "secondary" ~ "secondary",
    support_tier == "projection" ~ "projection",
    TRUE ~ "context"
  )
}

derive_support_tier_from_dataset_horizon <- function(dataset_name, horizon_year) {
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
  stop(sprintf("Unknown dataset `%s` while deriving support tier.", dataset_name), call. = FALSE)
}

derive_horizon_evidence_class_from_dataset_horizon <- function(dataset_name, horizon_year) {
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
  stop(sprintf("Unknown dataset `%s` while deriving evidence class.", dataset_name), call. = FALSE)
}

derive_claim_restriction_flag_from_evidence <- function(horizon_evidence_class) {
  dplyr::case_when(
    is.na(horizon_evidence_class) ~ NA_character_,
    horizon_evidence_class == "directly_observed_data_supported" ~ "primary_claim_allowed",
    horizon_evidence_class == "partly_model_dependent" ~ "secondary_or_sensitivity_only",
    horizon_evidence_class == "mostly_extrapolated" ~ "projection_only",
    TRUE ~ "projection_only"
  )
}

make_stage1_support_lookup <- function(horizon_registry) {
  n_rows <- nrow(horizon_registry)
  out <- tibble(
    dataset = as.character(horizon_registry$dataset),
    horizon_year = as.integer(horizon_registry$horizon_year),
    interpretation_tier = if ("interpretation_tier" %in% names(horizon_registry)) as.character(horizon_registry$interpretation_tier) else rep(NA_character_, n_rows),
    interpretation_note = if ("interpretation_note" %in% names(horizon_registry)) as.character(horizon_registry$interpretation_note) else rep(NA_character_, n_rows),
    primary_supported_flag = if ("primary_supported_flag" %in% names(horizon_registry)) as.logical(horizon_registry$primary_supported_flag) else rep(NA, n_rows),
    support_tier = if ("support_tier" %in% names(horizon_registry)) as.character(horizon_registry$support_tier) else rep(NA_character_, n_rows),
    horizon_evidence_class = if ("horizon_evidence_class" %in% names(horizon_registry)) as.character(horizon_registry$horizon_evidence_class) else rep(NA_character_, n_rows),
    claim_restriction_flag = if ("claim_restriction_flag" %in% names(horizon_registry)) as.character(horizon_registry$claim_restriction_flag) else rep(NA_character_, n_rows)
  )
  
  out <- out %>%
    mutate(
      support_tier = dplyr::case_when(
        !is.na(support_tier) ~ support_tier,
        !is.na(interpretation_tier) ~ normalize_interpretation_tier(interpretation_tier),
        TRUE ~ mapply(derive_support_tier_from_dataset_horizon, dataset, horizon_year, USE.NAMES = FALSE)
      ),
      support_label = support_tier,
      support_priority = make_support_priority(support_tier),
      primary_supported_flag = dplyr::case_when(
        !is.na(primary_supported_flag) ~ primary_supported_flag,
        TRUE ~ support_tier == "primary_supported"
      ),
      horizon_evidence_class = dplyr::case_when(
        !is.na(horizon_evidence_class) ~ horizon_evidence_class,
        TRUE ~ mapply(derive_horizon_evidence_class_from_dataset_horizon, dataset, horizon_year, USE.NAMES = FALSE)
      ),
      claim_restriction_flag = dplyr::case_when(
        !is.na(claim_restriction_flag) ~ claim_restriction_flag,
        TRUE ~ derive_claim_restriction_flag_from_evidence(horizon_evidence_class)
      ),
      interpretation_tier = dplyr::case_when(
        !is.na(interpretation_tier) ~ interpretation_tier,
        support_tier == "primary_supported" ~ "primary-supported",
        TRUE ~ support_tier
      ),
      interpretation_note = dplyr::case_when(
        !is.na(interpretation_note) ~ interpretation_note,
        TRUE ~ "Inherited from Stage 1 support hierarchy."
      )
    ) %>%
    arrange(match(dataset, dataset_order), horizon_year)
  
  out
}

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
  
  if (!all(c("dataset", "horizon_year") %in% names(stage1_inputs$horizon_registry))) {
    stop("Stage 1 horizon registry is missing required dataset/horizon columns.", call. = FALSE)
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
  
  if (!all(c("metadata_name", "metadata_value") %in% names(stage1_inputs$metadata_registry))) {
    stop("Stage 1 metadata registry is missing required columns.", call. = FALSE)
  }
  
  invisible(TRUE)
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
      support_tier = .data$support_tier,
      support_label = .data$support_label,
      support_priority = .data$support_priority,
      primary_supported_flag = .data$primary_supported_flag,
      horizon_evidence_class = .data$horizon_evidence_class,
      claim_restriction_flag = .data$claim_restriction_flag,
      interpretation_tier = .data$interpretation_tier,
      interpretation_note = .data$interpretation_note
    )
}

make_conservative_observed_contrast_support <- function(horizon_year) {
  h0 <- as.integer(horizon_year)
  support_tier <- dplyr::case_when(
    h0 == 1L ~ "primary_supported",
    h0 == 2L ~ "sensitivity",
    h0 >= 3L ~ "projection",
    TRUE ~ "projection"
  )
  horizon_evidence_class <- dplyr::case_when(
    h0 == 1L ~ "directly_observed_data_supported",
    h0 == 2L ~ "partly_model_dependent",
    TRUE ~ "mostly_extrapolated"
  )
  tibble(
    support_tier = support_tier,
    support_label = support_tier,
    support_priority = make_support_priority(support_tier),
    primary_supported_flag = support_tier == "primary_supported",
    horizon_evidence_class = horizon_evidence_class,
    claim_restriction_flag = derive_claim_restriction_flag_from_evidence(horizon_evidence_class),
    interpretation_tier = ifelse(support_tier == "primary_supported", "primary-supported", support_tier),
    interpretation_note = dplyr::case_when(
      h0 == 1L ~ "Common-window observed contrast at 1 year remains primary-supported.",
      h0 == 2L ~ "Common-window observed contrast at 2 years is treated conservatively as sensitivity only.",
      TRUE ~ "Common-window observed contrast beyond 2 years is projection-dominant only."
    )
  )
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

make_common_restricted_window_flag <- function(horizon_year) {
  horizon_year %in% supported_short_horizons_year
}

make_interval_support_fields <- function(interval_label) {
  out <- dplyr::case_when(
    interval_label == "0-1" ~ "primary_supported",
    interval_label == "1-2" ~ "primary_supported",
    interval_label == "2-5" ~ "secondary",
    interval_label == ">5" ~ "projection",
    TRUE ~ "projection"
  )
  evidence <- dplyr::case_when(
    interval_label %in% c("0-1", "1-2") ~ "directly_observed_data_supported",
    interval_label == "2-5" ~ "partly_model_dependent",
    TRUE ~ "mostly_extrapolated"
  )
  label <- dplyr::case_when(
    interval_label %in% c("0-1", "1-2") ~ "supported_short_horizon_contrast",
    interval_label == "2-5" ~ "intermediate_horizon_contrast",
    TRUE ~ "tail_diagnostic_contrast"
  )
  tibble(
    support_tier = out,
    support_label = out,
    support_priority = make_support_priority(out),
    primary_supported_flag = out == "primary_supported",
    horizon_evidence_class = evidence,
    claim_restriction_flag = derive_claim_restriction_flag_from_evidence(evidence),
    common_restricted_window_flag = interval_label %in% c("0-1", "1-2"),
    output_label = label
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
    arrange(match(dataset, dataset_order))
}

make_block2_dataset_summary <- function(analysis_datasets, stage1_dataset_registry, support_lookup) {
  support_counts <- support_lookup %>%
    group_by(dataset) %>%
    summarise(
      n_primary_supported_horizons = sum(primary_supported_flag, na.rm = TRUE),
      .groups = "drop"
    )
  
  source_block <- if (all(c("dataset", "source_description", "site_values") %in% names(stage1_dataset_registry))) {
    stage1_dataset_registry %>%
      transmute(
        dataset = as.character(dataset),
        source_description = as.character(source_description),
        site_values = as.character(site_values)
      )
  } else {
    tibble(dataset = names(analysis_datasets), source_description = NA_character_, site_values = NA_character_)
  }
  
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
      block2_role = dplyr::case_when(
        dataset_name == "merged" ~ "site_free_descriptive_plus_site_adjusted_contrast",
        TRUE ~ "separate_cohort_descriptive_and_contrast"
      )
    )
  })) %>%
    left_join(support_counts, by = "dataset") %>%
    left_join(source_block, by = "dataset") %>%
    arrange(match(dataset, dataset_order))
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
      formula_rhs = NA_character_,
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
          support_tier,
          support_label,
          support_priority,
          primary_supported_flag,
          horizon_evidence_class,
          claim_restriction_flag
        ),
      by = "horizon_year"
    ) %>%
    mutate(
      common_restricted_window_flag = make_common_restricted_window_flag(horizon_year),
      support_basis = dplyr::case_when(
        dataset == "merged" ~ make_support_basis("merged_pooled_observed_support"),
        TRUE ~ make_support_basis("site_specific_observed_followup")
      ),
      output_label = dplyr::case_when(
        horizon_year %in% supported_short_horizons_year ~ "supported_short_horizon_contrast",
        horizon_year <= 5 ~ "intermediate_horizon_contrast",
        TRUE ~ "tail_diagnostic_contrast"
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
  
  coef_vec <- tryCatch(stats::coef(fit_object), error = function(e) NULL)
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
    support_tier = rep(merged_support$support_tier, 3L),
    support_label = rep(merged_support$support_label, 3L),
    support_priority = rep(merged_support$support_priority, 3L),
    primary_supported_flag = rep(merged_support$primary_supported_flag, 3L),
    horizon_evidence_class = rep(merged_support$horizon_evidence_class, 3L),
    claim_restriction_flag = rep(merged_support$claim_restriction_flag, 3L),
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
      support_tier = merged_support$support_tier,
      support_label = merged_support$support_label,
      support_priority = merged_support$support_priority,
      primary_supported_flag = merged_support$primary_supported_flag,
      horizon_evidence_class = merged_support$horizon_evidence_class,
      claim_restriction_flag = merged_support$claim_restriction_flag,
      common_restricted_window_flag = make_common_restricted_window_flag(h0),
      support_basis = make_support_basis("merged_restricted_window_model_support"),
      reported_estimate_suppressed_flag = is.na(risk_estimate),
      output_label = "supported_short_horizon_contrast"
    )
  }))
  
  pnu_est <- site_rows %>% filter(site_counterfactual == pnu_site_label) %>% pull(estimate)
  pnu_est <- if (length(pnu_est) >= 1L) pnu_est[1] else NA_real_
  snu_est <- site_rows %>% filter(site_counterfactual == snu_site_label) %>% pull(estimate)
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
    support_tier = merged_support$support_tier,
    support_label = merged_support$support_label,
    support_priority = merged_support$support_priority,
    primary_supported_flag = merged_support$primary_supported_flag,
    horizon_evidence_class = merged_support$horizon_evidence_class,
    claim_restriction_flag = merged_support$claim_restriction_flag,
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
    filter(result_type == "site_standardized_risk", as.character(site_counterfactual) == pnu_site_label) %>%
    pull(estimate)
  snu_val <- rows %>%
    filter(result_type == "site_standardized_risk", as.character(site_counterfactual) == snu_site_label) %>%
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
  
  reference_df <- tryCatch(build_restricted_dataset(df, h0), error = function(e) NULL)
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
      boot_fit <- tryCatch(fit_restricted_cox(boot_df, h0, formula_rhs, model_name), error = function(e) NULL)
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
  
  support_tier_vec <- c(pnu_support$support_tier, snu_support$support_tier, contrast_support$support_tier)
  support_label_vec <- c(pnu_support$support_label, snu_support$support_label, contrast_support$support_label)
  support_priority_vec <- c(pnu_support$support_priority, snu_support$support_priority, contrast_support$support_priority)
  primary_flag_vec <- c(pnu_support$primary_supported_flag, snu_support$primary_supported_flag, contrast_support$primary_supported_flag)
  evidence_vec <- c(pnu_support$horizon_evidence_class, snu_support$horizon_evidence_class, contrast_support$horizon_evidence_class)
  claim_vec <- c(pnu_support$claim_restriction_flag, snu_support$claim_restriction_flag, contrast_support$claim_restriction_flag)
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
    support_tier = support_tier_vec,
    support_label = support_label_vec,
    support_priority = support_priority_vec,
    primary_supported_flag = primary_flag_vec,
    horizon_evidence_class = evidence_vec,
    claim_restriction_flag = claim_vec,
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
      pnu_risk <- tryCatch(compute_km_risk_value(boot_pnu, h0), error = function(e) NA_real_)
      snu_risk <- tryCatch(compute_km_risk_value(boot_snu, h0), error = function(e) NA_real_)
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

compute_supported_short_nonpiecewise_fresh <- function(pnu_data, snu_data, merged_data, support_lookup, merged_site_formula_registry, build_log) {
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
      support_tier = support_tier,
      support_label = support_label,
      support_priority = support_priority,
      primary_supported_flag = primary_supported_flag,
      horizon_evidence_class = horizon_evidence_class,
      claim_restriction_flag = claim_restriction_flag,
      common_restricted_window_flag = common_restricted_window_flag,
      support_basis = make_support_basis("merged_pooled_observed_support"),
      reported_estimate_suppressed_flag = !estimable_flag,
      output_label = "supported_short_horizon_contrast"
    )
  
  restricted_model_outputs_local <- list()
  restricted_supported_rows_list <- list()
  restricted_index <- 1L
  
  for (i in seq_len(nrow(merged_site_formula_registry))) {
    model_name <- paste0("restricted_", merged_site_formula_registry$formula_name[i])
    formula_rhs <- merged_site_formula_registry$formula_rhs[i]
    for (h in supported_short_horizons_year) {
      restricted_model_outputs_local[[paste(model_name, h, sep = "__")]] <- fit_restricted_cox(
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
  
  build_log <<- append_build_log(
    build_log,
    component = "supported_short_nonpiecewise",
    status = "recomputed",
    source_file = NA_character_,
    detail = "Observed 1y/2y and restricted standardized short-horizon rows were recomputed."
  )
  
  list(
    rows = bind_rows(observed_supported_rows, merged_site_free_context_rows, bind_rows(restricted_supported_rows_list)),
    restricted_models = restricted_model_outputs_local
  )
}

validate_existing_supported_short_nonpiecewise_cache <- function(df, merged_site_formula_registry) {
  required_cols <- c(
    "dataset", "source_block", "model_name", "formula_rhs", "horizon_year", "result_type",
    "estimate", "conf_low", "conf_high", "analysis_subject_n_total", "analysis_transition_events_total",
    "at_risk_subject_n_total", "transition_events_total", "pnu_event_n", "snu_event_n",
    "estimable_flag", "reported_estimate_suppressed_flag"
  )
  if (is.null(df) || !all(required_cols %in% names(df))) {
    return(FALSE)
  }
  nonpiecewise <- df %>%
    filter(source_block %in% c("observed_km_restricted", "observed_km_pooled_merged", "merged_restricted_cox_standardized"))
  if (nrow(nonpiecewise) == 0L) {
    return(FALSE)
  }
  if (!identical(sort(unique(as.integer(nonpiecewise$horizon_year))), supported_short_horizons_year)) {
    return(FALSE)
  }
  expected_restricted_models <- paste0("restricted_", merged_site_formula_registry$formula_name)
  restricted_rows <- nonpiecewise %>% filter(source_block == "merged_restricted_cox_standardized")
  if (!all(expected_restricted_models %in% unique(restricted_rows$model_name))) {
    return(FALSE)
  }
  formula_map <- setNames(merged_site_formula_registry$formula_rhs, expected_restricted_models)
  formula_ok <- restricted_rows %>%
    mutate(expected_formula_rhs = unname(formula_map[model_name])) %>%
    summarise(ok = all(expected_formula_rhs == formula_rhs)) %>%
    pull(ok)
  if (!isTRUE(formula_ok)) {
    return(FALSE)
  }
  TRUE
}

harmonize_supported_short_nonpiecewise_cache <- function(df, support_lookup) {
  nonpiecewise <- df %>%
    filter(source_block %in% c("observed_km_restricted", "observed_km_pooled_merged", "merged_restricted_cox_standardized"))
  
  if (!"support_tier" %in% names(nonpiecewise)) nonpiecewise$support_tier <- NA_character_
  if (!"support_label" %in% names(nonpiecewise)) nonpiecewise$support_label <- nonpiecewise$support_tier
  if (!"support_priority" %in% names(nonpiecewise)) nonpiecewise$support_priority <- make_support_priority(nonpiecewise$support_tier)
  if (!"primary_supported_flag" %in% names(nonpiecewise)) nonpiecewise$primary_supported_flag <- nonpiecewise$support_tier == "primary_supported"
  if (!"horizon_evidence_class" %in% names(nonpiecewise)) nonpiecewise$horizon_evidence_class <- NA_character_
  if (!"claim_restriction_flag" %in% names(nonpiecewise)) nonpiecewise$claim_restriction_flag <- NA_character_
  if (!"common_restricted_window_flag" %in% names(nonpiecewise)) nonpiecewise$common_restricted_window_flag <- nonpiecewise$horizon_year %in% supported_short_horizons_year
  if (!"support_basis" %in% names(nonpiecewise)) nonpiecewise$support_basis <- NA_character_
  if (!"bootstrap_success" %in% names(nonpiecewise)) nonpiecewise$bootstrap_success <- NA_integer_
  if (!"bootstrap_iterations" %in% names(nonpiecewise)) nonpiecewise$bootstrap_iterations <- NA_integer_
  if (!"bootstrap_success_rate" %in% names(nonpiecewise)) nonpiecewise$bootstrap_success_rate <- NA_real_
  if (!"bootstrap_instability_flag" %in% names(nonpiecewise)) nonpiecewise$bootstrap_instability_flag <- NA
  if (!"bootstrap_note" %in% names(nonpiecewise)) nonpiecewise$bootstrap_note <- NA_character_
  if (!"max_followup_years" %in% names(nonpiecewise)) nonpiecewise$max_followup_years <- NA_real_
  if (!"available_within_observed_followup" %in% names(nonpiecewise)) nonpiecewise$available_within_observed_followup <- nonpiecewise$estimable_flag
  if (!"availability_note" %in% names(nonpiecewise)) nonpiecewise$availability_note <- ifelse(nonpiecewise$estimable_flag, "inherited_existing_estimate", "not_estimable")
  if (!"support_issue_reason" %in% names(nonpiecewise)) nonpiecewise$support_issue_reason <- NA_character_
  if (!"site_counterfactual" %in% names(nonpiecewise)) nonpiecewise$site_counterfactual <- NA_character_
  if (!"contrast_direction" %in% names(nonpiecewise)) nonpiecewise$contrast_direction <- NA_character_
  if (!"interval_label" %in% names(nonpiecewise)) nonpiecewise$interval_label <- NA_character_
  if (!"window_year" %in% names(nonpiecewise)) nonpiecewise$window_year <- nonpiecewise$horizon_year
  if (!"output_label" %in% names(nonpiecewise)) nonpiecewise$output_label <- "supported_short_horizon_contrast"
  
  join_from_lookup <- function(df_in, dataset_value) {
    df_in %>%
      left_join(
        support_lookup %>%
          filter(dataset == dataset_value) %>%
          select(horizon_year, support_tier, support_label, support_priority, primary_supported_flag, horizon_evidence_class, claim_restriction_flag),
        by = "horizon_year",
        suffix = c("", "__new")
      ) %>%
      mutate(
        support_tier = coalesce_character(support_tier, support_tier__new),
        support_label = coalesce_character(support_label, support_label__new),
        support_priority = coalesce_character(support_priority, support_priority__new),
        primary_supported_flag = dplyr::coalesce(primary_supported_flag, primary_supported_flag__new),
        horizon_evidence_class = coalesce_character(horizon_evidence_class, horizon_evidence_class__new),
        claim_restriction_flag = coalesce_character(claim_restriction_flag, claim_restriction_flag__new)
      ) %>%
      select(-ends_with("__new"))
  }
  
  observed_by_site <- nonpiecewise %>% filter(source_block == "observed_km_restricted", dataset %in% c("PNU", "SNU"))
  observed_by_site <- bind_rows(
    join_from_lookup(observed_by_site %>% filter(dataset == "PNU"), "PNU"),
    join_from_lookup(observed_by_site %>% filter(dataset == "SNU"), "SNU")
  )
  
  observed_contrast <- nonpiecewise %>%
    filter(source_block == "observed_km_restricted", dataset == "PNU_vs_SNU") %>%
    rowwise() %>%
    mutate(
      support_stub = list(make_conservative_observed_contrast_support(horizon_year)),
      support_tier = coalesce_character(support_tier, support_stub[[1]]$support_tier),
      support_label = coalesce_character(support_label, support_stub[[1]]$support_label),
      support_priority = coalesce_character(support_priority, support_stub[[1]]$support_priority),
      primary_supported_flag = dplyr::coalesce(primary_supported_flag, support_stub[[1]]$primary_supported_flag),
      horizon_evidence_class = coalesce_character(horizon_evidence_class, support_stub[[1]]$horizon_evidence_class),
      claim_restriction_flag = coalesce_character(claim_restriction_flag, support_stub[[1]]$claim_restriction_flag)
    ) %>%
    ungroup() %>%
    select(-support_stub)
  
  merged_rows <- nonpiecewise %>%
    filter(source_block %in% c("observed_km_pooled_merged", "merged_restricted_cox_standardized")) %>%
    left_join(
      support_lookup %>%
        filter(dataset == "merged") %>%
        select(horizon_year, support_tier, support_label, support_priority, primary_supported_flag, horizon_evidence_class, claim_restriction_flag),
      by = "horizon_year",
      suffix = c("", "__new")
    ) %>%
    mutate(
      support_tier = coalesce_character(support_tier, support_tier__new),
      support_label = coalesce_character(support_label, support_label__new),
      support_priority = coalesce_character(support_priority, support_priority__new),
      primary_supported_flag = dplyr::coalesce(primary_supported_flag, primary_supported_flag__new),
      horizon_evidence_class = coalesce_character(horizon_evidence_class, horizon_evidence_class__new),
      claim_restriction_flag = coalesce_character(claim_restriction_flag, claim_restriction_flag__new)
    ) %>%
    select(-ends_with("__new"))
  
  bind_rows(observed_by_site, observed_contrast, merged_rows) %>%
    mutate(
      support_label = coalesce_character(support_label, support_tier),
      support_priority = coalesce_character(support_priority, make_support_priority(support_tier)),
      common_restricted_window_flag = horizon_year %in% supported_short_horizons_year,
      output_label = "supported_short_horizon_contrast"
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
          interval_id == 3L ~ "2-5",
          TRUE ~ ">5"
        ),
        levels = c("0-1", "1-2", "2-5", ">5")
      ),
      site_snu = as.integer(as.character(site) == snu_site_label),
      site_snu_int_0_1 = as.integer(interval_label == "0-1") * site_snu,
      site_snu_int_1_2 = as.integer(interval_label == "1-2") * site_snu,
      site_snu_int_2_5 = as.integer(interval_label == "2-5") * site_snu,
      site_snu_int_gt_5 = as.integer(interval_label == ">5") * site_snu
    )
}

summarize_piecewise_interval_support <- function(split_df, min_events_per_site, min_subjects_per_site) {
  site_support <- split_df %>%
    mutate(interval_label = as.character(interval_label), site = as.character(site)) %>%
    group_by(interval_label, site) %>%
    summarise(
      site_subject_n = dplyr::n_distinct(unique_person_id),
      site_event_n = sum(event_piece == 1L),
      site_person_years = sum(pmax(tstop_year - tstart_year, 0)),
      .groups = "drop"
    )
  
  site_template <- tidyr::expand_grid(
    interval_label = c("0-1", "1-2", "2-5", ">5"),
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
      interval_label = factor(interval_label, levels = c("0-1", "1-2", "2-5", ">5")),
      interval_start_year = case_when(
        interval_label == "0-1" ~ 0,
        interval_label == "1-2" ~ 1,
        interval_label == "2-5" ~ 2,
        TRUE ~ 5
      ),
      interval_end_year = case_when(
        interval_label == "0-1" ~ 1,
        interval_label == "1-2" ~ 2,
        interval_label == "2-5" ~ 5,
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
      )
    )
  
  support_stub <- bind_rows(lapply(as.character(interval_support$interval_label), function(lbl) {
    bind_cols(tibble(interval_label = lbl), make_interval_support_fields(lbl))
  }))
  
  interval_support %>%
    left_join(support_stub, by = "interval_label") %>%
    mutate(
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
      support_tier, support_label, support_priority, primary_supported_flag,
      horizon_evidence_class, claim_restriction_flag,
      common_restricted_window_flag, support_basis, output_label
    ) %>%
    arrange(match(interval_label, c("0-1", "1-2", "2-5", ">5")))
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
    "site_snu_int_2_5",
    "site_snu_int_gt_5",
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
    term = c("site_snu_int_0_1", "site_snu_int_1_2", "site_snu_int_2_5", "site_snu_int_gt_5"),
    interval_label = c("0-1", "1-2", "2-5", ">5"),
    interval_start_year = c(0, 1, 2, 5),
    interval_end_year = c(1, 2, 5, NA_real_)
  )
  
  interval_support_join <- interval_support_df %>%
    select(
      interval_label,
      at_risk_subject_n_total, transition_events_total, person_years_total,
      pnu_subject_n, snu_subject_n, pnu_event_n, snu_event_n,
      pnu_person_years, snu_person_years,
      estimable_flag, support_issue_reason,
      support_tier, support_label, support_priority, primary_supported_flag,
      horizon_evidence_class, claim_restriction_flag,
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
      support_tier, support_label, support_priority, primary_supported_flag, horizon_evidence_class, claim_restriction_flag,
      common_restricted_window_flag, support_basis, output_label
    )
  
  list(
    fit = fit_object,
    summary_table = summary_table,
    model_name = model_name,
    formula_rhs = ifelse(nzchar(trimws(covariate_rhs)), covariate_rhs, NA_character_)
  )
}

validate_existing_piecewise_cache <- function(piecewise_site_effect_df, piecewise_interval_support_df) {
  if (is.null(piecewise_site_effect_df) || is.null(piecewise_interval_support_df)) {
    return(FALSE)
  }
  required_site_cols <- c("interval_label", "model_name", "estimate_hr", "conf_low_hr", "conf_high_hr", "output_label")
  required_support_cols <- c("interval_label", "pnu_event_n", "snu_event_n", "output_label")
  if (!all(required_site_cols %in% names(piecewise_site_effect_df))) return(FALSE)
  if (!all(required_support_cols %in% names(piecewise_interval_support_df))) return(FALSE)
  if (!identical(sort(unique(as.character(piecewise_site_effect_df$interval_label))), sort(c("0-1", "1-2", "2-5", ">5")))) return(FALSE)
  if (!identical(sort(unique(as.character(piecewise_interval_support_df$interval_label))), sort(c("0-1", "1-2", "2-5", ">5")))) return(FALSE)
  expected_models <- c("piecewise_unadjusted", "piecewise_site_added", "piecewise_site_interaction")
  if (!all(expected_models %in% unique(as.character(piecewise_site_effect_df$model_name)))) return(FALSE)
  TRUE
}

compute_piecewise_components_fresh <- function(merged_data, build_log) {
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
    arrange(match(model_name, piecewise_model_specs$model_name), match(interval_label, c("0-1", "1-2", "2-5", ">5"))) %>%
    mutate(
      reporting_status = make_reporting_status(estimable_flag, reported_estimate_suppressed_flag),
      plotting_status = make_plotting_status(estimable_flag, reported_estimate_suppressed_flag),
      plot_inclusion_flag = estimable_flag & !reported_estimate_suppressed_flag,
      reporting_note = make_reporting_note(
        reporting_status = reporting_status,
        availability_note = availability_note,
        support_issue_reason = support_issue_reason
      )
    )
  
  build_log <<- append_build_log(
    build_log,
    component = "piecewise_backbone",
    status = "recomputed",
    source_file = NA_character_,
    detail = "Piecewise interval support and piecewise site-effect models were recomputed under the 0-1 / 1-2 / 2-5 / >5 backbone."
  )
  
  list(
    piecewise_interval_support = piecewise_interval_support,
    piecewise_site_effect = piecewise_site_effect,
    piecewise_model_objects = piecewise_model_objects
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
      output_label = case_when(
        band_end_year <= 2 ~ "supported_short_horizon_contrast",
        band_end_year <= 5 ~ "intermediate_horizon_contrast",
        TRUE ~ "tail_diagnostic_contrast"
      )
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
        select(horizon_year, support_tier, support_label, support_priority, primary_supported_flag, horizon_evidence_class, claim_restriction_flag),
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
      support_tier, support_label, support_priority, primary_supported_flag,
      horizon_evidence_class, claim_restriction_flag,
      common_restricted_window_flag, support_basis, output_label
    )
}

## 🟠 Define: interval contrast extractors ===============================
convert_piecewise_rows_to_contrast <- function(piecewise_site_effect, interval_labels, output_label_value, horizon_year_value = NA_real_) {
  piecewise_site_effect %>%
    filter(interval_label %in% interval_labels) %>%
    transmute(
      dataset = dataset,
      analysis_view = analysis_view,
      source_block = source_block,
      model_name = model_name,
      formula_rhs = formula_rhs,
      horizon_year = horizon_year_value,
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
      support_tier = support_tier,
      support_label = support_label,
      support_priority = support_priority,
      primary_supported_flag = primary_supported_flag,
      horizon_evidence_class = horizon_evidence_class,
      claim_restriction_flag = claim_restriction_flag,
      common_restricted_window_flag = common_restricted_window_flag,
      support_basis = support_basis,
      reported_estimate_suppressed_flag = reported_estimate_suppressed_flag,
      output_label = output_label_value
    )
}

## 🟠 Define: validation helpers ===============================
validate_supported_short_counts <- function(supported_df, risk_trajectory, analysis_totals_lookup, piecewise_interval_support) {
  analysis_lookup <- analysis_totals_lookup %>%
    rename(
      expected_analysis_subject_n_total = analysis_subject_n_total,
      expected_analysis_transition_events_total = analysis_transition_events_total
    )
  
  observed_noncontrast_check <- supported_df %>%
    filter(result_type == "observed_risk", source_block %in% c("observed_km_restricted", "observed_km_pooled_merged")) %>%
    select(dataset, horizon_year, analysis_subject_n_total, analysis_transition_events_total, at_risk_subject_n_total, transition_events_total) %>%
    left_join(
      risk_trajectory %>% select(dataset, horizon_year, n_risk, cumulative_transition_n),
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
    filter(dataset == "PNU_vs_SNU", source_block == "observed_km_restricted", result_type == "absolute_risk_difference") %>%
    select(horizon_year, analysis_subject_n_total, analysis_transition_events_total, at_risk_subject_n_total, transition_events_total, pnu_event_n, snu_event_n) %>%
    left_join(
      risk_trajectory %>% filter(dataset == "PNU") %>% transmute(horizon_year, expected_pnu_n_risk = n_risk, expected_pnu_transition_events = cumulative_transition_n),
      by = "horizon_year"
    ) %>%
    left_join(
      risk_trajectory %>% filter(dataset == "SNU") %>% transmute(horizon_year, expected_snu_n_risk = n_risk, expected_snu_transition_events = cumulative_transition_n),
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
    filter(source_block == "merged_restricted_cox_standardized", result_type %in% c("site_standardized_risk", "absolute_risk_difference")) %>%
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
    select(model_name, interval_label, analysis_subject_n_total, analysis_transition_events_total, at_risk_subject_n_total, transition_events_total, pnu_event_n, snu_event_n) %>%
    left_join(
      piecewise_interval_support %>%
        select(interval_label, expected_at_risk_subject_n_total = at_risk_subject_n_total, expected_transition_events_total = transition_events_total, expected_pnu_event_n = pnu_event_n, expected_snu_event_n = snu_event_n),
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
    if (any(!observed_noncontrast_check$passed)) sprintf("supported observed rows: %s", format_check_locations(observed_noncontrast_check %>% filter(!passed), c("dataset", "horizon_year"))) else NULL,
    if (any(!contrast_check$passed)) sprintf("supported observed contrast rows: %s", format_check_locations(contrast_check %>% filter(!passed), c("horizon_year"))) else NULL,
    if (any(!standardized_check$passed)) sprintf("supported standardized rows: %s", format_check_locations(standardized_check %>% filter(!passed), c("model_name", "horizon_year", "result_type"))) else NULL,
    if (any(!piecewise_short_check$passed)) sprintf("supported piecewise rows: %s", format_check_locations(piecewise_short_check %>% filter(!passed), c("model_name", "interval_label"))) else NULL
  )
  
  if (length(failed_ids) > 0L) {
    stop(sprintf("Supported short-horizon count alignment failed: %s", paste(failed_ids, collapse = " | ")), call. = FALSE)
  }
  invisible(TRUE)
}

validate_piecewise_interval_counts <- function(interval_df, piecewise_interval_support, analysis_totals_lookup, interval_labels, object_label) {
  if (nrow(interval_df) == 0L) {
    return(invisible(TRUE))
  }
  analysis_lookup <- analysis_totals_lookup %>%
    rename(
      expected_analysis_subject_n_total = analysis_subject_n_total,
      expected_analysis_transition_events_total = analysis_transition_events_total
    )
  check_df <- interval_df %>%
    filter(result_type == "hazard_ratio", interval_label %in% interval_labels) %>%
    select(model_name, interval_label, analysis_subject_n_total, analysis_transition_events_total, at_risk_subject_n_total, transition_events_total, pnu_event_n, snu_event_n) %>%
    left_join(
      piecewise_interval_support %>%
        select(interval_label, expected_at_risk_subject_n_total = at_risk_subject_n_total, expected_transition_events_total = transition_events_total, expected_pnu_event_n = pnu_event_n, expected_snu_event_n = snu_event_n),
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
  if (any(!check_df$passed)) {
    stop(sprintf("%s count alignment failed: %s", object_label, format_check_locations(check_df %>% filter(!passed), c("model_name", "interval_label"))), call. = FALSE)
  }
  invisible(TRUE)
}

validate_reporting_flags <- function(supported_df, intermediate_df, tail_df, piecewise_df) {
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
  
  intermediate_suppressed_check <- intermediate_df %>%
    filter(reported_estimate_suppressed_flag) %>%
    mutate(
      passed = reporting_status == "suppressed" &
        plotting_status == "suppressed" &
        !plot_inclusion_flag
    )
  
  tail_suppressed_check <- tail_df %>%
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
    if (nrow(supported_bootstrap_check) > 0L && any(!supported_bootstrap_check$passed)) sprintf("supported bootstrap-flag rows: %s", format_check_locations(supported_bootstrap_check %>% filter(!passed), c("dataset", "horizon_year", "model_name", "result_type"))) else NULL,
    if (nrow(supported_suppressed_check) > 0L && any(!supported_suppressed_check$passed)) sprintf("supported suppressed rows: %s", format_check_locations(supported_suppressed_check %>% filter(!passed), c("dataset", "horizon_year", "model_name", "result_type"))) else NULL,
    if (nrow(intermediate_suppressed_check) > 0L && any(!intermediate_suppressed_check$passed)) sprintf("intermediate suppressed rows: %s", format_check_locations(intermediate_suppressed_check %>% filter(!passed), c("model_name", "interval_label"))) else NULL,
    if (nrow(tail_suppressed_check) > 0L && any(!tail_suppressed_check$passed)) sprintf("tail suppressed rows: %s", format_check_locations(tail_suppressed_check %>% filter(!passed), c("model_name", "interval_label"))) else NULL,
    if (nrow(piecewise_suppressed_check) > 0L && any(!piecewise_suppressed_check$passed)) sprintf("piecewise suppressed rows: %s", format_check_locations(piecewise_suppressed_check %>% filter(!passed), c("model_name", "interval_label"))) else NULL
  )
  
  if (length(failed_ids) > 0L) {
    stop(sprintf("Reporting-flag alignment failed: %s", paste(failed_ids, collapse = " | ")), call. = FALSE)
  }
  invisible(TRUE)
}

## 🟠 Define: plot constructors ===============================
make_empty_plot <- function(title_text, subtitle_text, body_text) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = body_text, size = 5) +
    xlim(0, 1) + ylim(0, 1) +
    labs(title = title_text, subtitle = subtitle_text, x = NULL, y = NULL) +
    theme_void()
}


plot_interval_support_fallback <- function(interval_df, piecewise_support_df, interval_label_value, title_text, subtitle_text) {
  support_row <- piecewise_support_df %>%
    filter(as.character(interval_label) == interval_label_value)
  
  if (nrow(support_row) == 0L) {
    return(make_empty_plot(
      title_text = title_text,
      subtitle_text = subtitle_text,
      body_text = paste("No support row found for interval", interval_label_value)
    ))
  }
  
  model_rows <- interval_df %>%
    filter(result_type == "hazard_ratio", as.character(interval_label) == interval_label_value) %>%
    mutate(
      support_or_fit_reason = coalesce_character(support_issue_reason, availability_note),
      support_or_fit_reason = coalesce_character(support_or_fit_reason, "suppressed_without_explicit_reason"),
      display_reason = paste0(model_name, ": ", support_or_fit_reason)
    )
  
  support_df <- tibble(
    site_name = factor(c(pnu_site_label, snu_site_label), levels = c(pnu_site_label, snu_site_label)),
    event_n = c(support_row$pnu_event_n[1], support_row$snu_event_n[1]),
    subject_n = c(support_row$pnu_subject_n[1], support_row$snu_subject_n[1]),
    person_years = c(support_row$pnu_person_years[1], support_row$snu_person_years[1])
  ) %>%
    mutate(
      event_n = dplyr::coalesce(as.numeric(event_n), 0),
      subject_n = dplyr::coalesce(as.numeric(subject_n), 0),
      person_years = dplyr::coalesce(as.numeric(person_years), 0),
      support_label_text = sprintf(
        "events=%s\nsubjects=%s\nPY=%s",
        event_n,
        subject_n,
        format(round(person_years, 2), nsmall = 2)
      )
    )
  
  suppression_text <- if (nrow(model_rows) > 0L) {
    paste(model_rows$display_reason, collapse = "\n")
  } else {
    "No model rows were available for this interval."
  }
  
  max_y <- suppressWarnings(max(support_df$event_n, na.rm = TRUE))
  if (!is.finite(max_y) || max_y < 1) {
    max_y <- 1
  }
  
  interval_support_reason <- coalesce_character(support_row$support_issue_reason[1], "none_recorded")
  
  ggplot(support_df, aes(x = site_name, y = event_n)) +
    geom_col(width = 0.6, na.rm = TRUE) +
    geom_text(aes(label = support_label_text), vjust = -0.15, size = 3, lineheight = 0.95, na.rm = TRUE) +
    annotate(
      geom = "label",
      x = 1.5,
      y = max_y * 1.65,
      label = paste(
        "Diagnostic fallback: no reportable hazard-ratio rows",
        suppression_text,
        sep = "\n"
      ),
      size = 3,
      label.size = 0.25,
      hjust = 0.5,
      vjust = 1
    ) +
    coord_cartesian(ylim = c(0, max_y * 1.95), clip = "off") +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Site",
      y = "Transition events in interval",
      caption = paste0("Interval support reason: ", interval_support_reason)
    ) +
    theme_bw()
}

plot_risk_trajectory <- function(risk_trajectory_df) {
  ggplot(risk_trajectory_df, aes(x = horizon_year, y = risk, color = dataset, group = dataset)) +
    geom_line(na.rm = TRUE, linewidth = 0.8) +
    geom_point(na.rm = TRUE, size = 2) +
    scale_x_continuous(breaks = 1:10) +
    labs(
      title = "Block 2 risk trajectory: observed KM risk by dataset",
      subtitle = "Annual KM risks are retained for context, but main timing separation uses interval-specific contrasts (0-1, 1-2, 2-5, >5).",
      x = "Horizon (years)",
      y = "Risk = 1 - KM survival",
      color = "Dataset"
    ) +
    theme_bw()
}

plot_supported_short <- function(supported_df) {
  plot_df <- supported_df %>%
    filter(result_type == "absolute_risk_difference") %>%
    mutate(
      model_name = factor(model_name, levels = unique(model_name)),
      plotting_status = factor(plotting_status, levels = c("stable", "bootstrap_unstable", "suppressed"))
    )
  if (nrow(plot_df) == 0L) {
    return(make_empty_plot(
      title_text = "Supported short-horizon contrast",
      subtitle_text = "No estimable short-horizon contrast rows were available for plotting.",
      body_text = "No plotted short-horizon rows"
    ))
  }
  dodge_width <- 0.35
  ggplot(plot_df, aes(x = factor(horizon_year), y = estimate, color = model_name, group = model_name)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_linerange(aes(ymin = conf_low, ymax = conf_high), position = position_dodge(width = dodge_width), na.rm = TRUE) +
    geom_point(aes(shape = plotting_status), position = position_dodge(width = dodge_width), size = 2.5, na.rm = TRUE) +
    geom_text(
      data = plot_df %>% filter(plotting_status == "bootstrap_unstable"),
      aes(label = "unstable"),
      position = position_dodge(width = dodge_width),
      vjust = -0.8,
      size = 3,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    scale_shape_manual(values = c(stable = 16, bootstrap_unstable = 1, suppressed = 4), drop = FALSE) +
    labs(
      title = "Supported short-horizon contrast: PNU minus SNU absolute risk difference",
      subtitle = "Only 1-year and 2-year common-window contrasts are treated as the short-horizon supported block.",
      x = "Restricted horizon (years)",
      y = "Absolute risk difference",
      color = "Source",
      shape = "Reporting status"
    ) +
    theme_bw()
}

plot_intermediate_interval <- function(intermediate_df, piecewise_support_df) {
  plot_df <- intermediate_df %>%
    filter(result_type == "hazard_ratio", plot_inclusion_flag)
  if (nrow(plot_df) == 0L) {
    return(plot_interval_support_fallback(
      interval_df = intermediate_df,
      piecewise_support_df = piecewise_support_df,
      interval_label_value = intermediate_interval_label,
      title_text = "Intermediate-horizon contrast: 2-5 year site-effect hazard ratio",
      subtitle_text = "No reportable hazard-ratio rows survived suppression; fallback diagnostics show site-specific interval support and suppression reasons."
    ))
  }
  ggplot(plot_df, aes(x = interval_label, y = estimate, color = model_name, group = model_name)) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
    geom_pointrange(aes(ymin = conf_low, ymax = conf_high), position = position_dodge(width = 0.3), na.rm = TRUE) +
    scale_y_log10() +
    labs(
      title = "Intermediate-horizon contrast: 2-5 year site-effect hazard ratio",
      subtitle = "This interval is secondary and partly model-dependent under the current support hierarchy.",
      x = "Time interval (years)",
      y = "Hazard ratio (log scale)",
      color = "Model"
    ) +
    theme_bw()
}

plot_tail_interval <- function(tail_df, piecewise_support_df) {
  plot_df <- tail_df %>%
    filter(result_type == "hazard_ratio", plot_inclusion_flag)
  if (nrow(plot_df) == 0L) {
    return(plot_interval_support_fallback(
      interval_df = tail_df,
      piecewise_support_df = piecewise_support_df,
      interval_label_value = tail_interval_label,
      title_text = "Tail diagnostic contrast: >5 year site-effect hazard ratio",
      subtitle_text = "No reportable tail hazard-ratio rows survived suppression; fallback diagnostics show site-specific interval support and suppression reasons."
    ))
  }
  ggplot(plot_df, aes(x = interval_label, y = estimate, color = model_name, group = model_name)) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
    geom_pointrange(aes(ymin = conf_low, ymax = conf_high), position = position_dodge(width = 0.3), na.rm = TRUE) +
    scale_y_log10() +
    labs(
      title = "Tail diagnostic contrast: >5 year site-effect hazard ratio",
      subtitle = "This tail segment is a diagnostic extension only and must not drive primary timing interpretation.",
      x = "Time interval (years)",
      y = "Hazard ratio (log scale)",
      color = "Model"
    ) +
    theme_bw()
}

plot_piecewise_support <- function(piecewise_support_df) {
  plot_df <- piecewise_support_df %>%
    select(interval_label, pnu_event_n, snu_event_n) %>%
    pivot_longer(cols = c(pnu_event_n, snu_event_n), names_to = "site_name", values_to = "event_n") %>%
    mutate(
      interval_label = factor(interval_label, levels = c("0-1", "1-2", "2-5", ">5")),
      site_name = factor(
        case_when(site_name == "pnu_event_n" ~ pnu_site_label, TRUE ~ snu_site_label),
        levels = c(pnu_site_label, snu_site_label)
      )
    )
  ggplot(plot_df, aes(x = interval_label, y = event_n, fill = site_name)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, na.rm = TRUE) +
    labs(
      title = "Piecewise interval support: site-specific transition events",
      subtitle = "The 2-5 year interval is the intermediate backbone, whereas >5 years is treated only as a tail diagnostic extension.",
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
      subtitle = "Hazard bands retain 0-1, 1-2, 2-5, and 5-10 year structure; the 5-10 band is used only as a tail-context diagnostic.",
      x = "Time band (years)",
      y = "Crude hazard per 100 person-years",
      color = "Dataset"
    ) +
    theme_bw()
}

save_plot_png <- function(plot_object, filename, width = 11, height = 8.5, dpi = 300) {
  ggplot2::ggsave(
    filename = file.path(supporting_output_dir, filename),
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    units = "in"
  )
}

# 🔴 Read: Stage 1 backbone inputs ===============================
stage1_inputs <- read_stage1_inputs()

# 🔴 Validate: inherited registries and prepared datasets ===============================
validate_stage1_inputs(stage1_inputs)

common_horizons_year <- sort(unique(as.integer(stage1_inputs$horizon_registry$horizon_year)))
risk_thresholds <- sort(unique(as.numeric(stage1_inputs$threshold_registry$threshold)))
support_lookup <- make_stage1_support_lookup(stage1_inputs$horizon_registry)
analysis_datasets <- lapply(stage1_inputs$analysis_datasets, add_model_time_year)
analysis_totals_lookup <- make_analysis_totals_lookup(analysis_datasets)
block2_dataset_summary <- make_block2_dataset_summary(
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

pnu_data <- analysis_datasets[["PNU"]]
snu_data <- analysis_datasets[["SNU"]]
merged_data <- analysis_datasets[["merged"]]

build_log <- new_build_log()

# 🔴 Compute: observed trajectories and hazard summaries ===============================
risk_trajectory <- bind_rows(lapply(names(analysis_datasets), function(dataset_name) {
  compute_km_risk_trajectory(
    df = analysis_datasets[[dataset_name]],
    dataset_name = dataset_name,
    support_lookup = support_lookup
  )
})) %>%
  arrange(match(dataset, dataset_order), horizon_year)

hazard_pattern_summary <- bind_rows(lapply(names(analysis_datasets), function(dataset_name) {
  summarize_hazard_bands(
    df = analysis_datasets[[dataset_name]],
    dataset_name = dataset_name,
    band_breaks = hazard_band_breaks_year,
    support_lookup = support_lookup
  )
})) %>%
  arrange(match(dataset, dataset_order), band_start_year)

build_log <- append_build_log(
  build_log,
  component = "risk_and_hazard_context",
  status = "recomputed",
  source_file = NA_character_,
  detail = "Annual risk trajectories and hazard-band summaries were recomputed from Stage 1 analysis datasets."
)

# 🔴 Recover: reusable block components when compatible ===============================
existing_supported_short_file <- resolve_existing_output_file(
  export_path,
  "block2_supported_short_horizon_contrast.csv",
  legacy_stage4_export_path,
  "stage4_supported_short_horizon_contrast.csv"
)
existing_piecewise_support_file <- resolve_existing_output_file(
  supporting_output_dir,
  "block2_piecewise_interval_support.csv",
  export_path,
  "block2_piecewise_interval_support.csv",
  legacy_stage4_export_path,
  "stage4_piecewise_interval_support.csv"
)
existing_piecewise_site_effect_file <- resolve_existing_output_file(
  supporting_output_dir,
  "block2_piecewise_site_effect.csv",
  export_path,
  "block2_piecewise_site_effect.csv",
  legacy_stage4_export_path,
  "stage4_piecewise_site_effect.csv"
)
existing_model_objects_file <- resolve_existing_output_file(
  supporting_output_dir,
  "block2_model_objects.rds",
  export_path,
  "block2_model_objects.rds",
  legacy_stage4_export_path,
  "stage4_model_objects.rds"
)

existing_supported_short <- safe_read_optional_csv(existing_supported_short_file)
existing_piecewise_support <- safe_read_optional_csv(existing_piecewise_support_file)
existing_piecewise_site_effect <- safe_read_optional_csv(existing_piecewise_site_effect_file)
existing_model_objects <- safe_read_optional_rds(existing_model_objects_file)

use_existing_supported_short_nonpiecewise <- validate_existing_supported_short_nonpiecewise_cache(existing_supported_short, merged_site_formula_registry)
use_existing_piecewise_cache <- validate_existing_piecewise_cache(existing_piecewise_site_effect, existing_piecewise_support)

# 🔴 Build: piecewise interval backbone and interval contrasts ===============================
if (use_existing_piecewise_cache) {
  piecewise_interval_support <- existing_piecewise_support
  piecewise_site_effect <- existing_piecewise_site_effect
  piecewise_model_objects <- if (!is.null(existing_model_objects) && !is.null(existing_model_objects$piecewise_models)) {
    existing_model_objects$piecewise_models
  } else {
    list()
  }
  build_log <- append_build_log(
    build_log,
    component = "piecewise_backbone",
    status = "reused",
    source_file = existing_piecewise_site_effect_file,
    detail = "Existing piecewise interval support and site-effect outputs matched the new 0-1 / 1-2 / 2-5 / >5 schema and were reused."
  )
} else {
  piecewise_components <- compute_piecewise_components_fresh(merged_data, build_log)
  piecewise_interval_support <- piecewise_components$piecewise_interval_support
  piecewise_site_effect <- piecewise_components$piecewise_site_effect
  piecewise_model_objects <- piecewise_components$piecewise_model_objects
}

supported_piecewise_rows <- convert_piecewise_rows_to_contrast(
  piecewise_site_effect = piecewise_site_effect,
  interval_labels = c("0-1", "1-2"),
  output_label_value = "supported_short_horizon_contrast",
  horizon_year_value = NA_real_
)

if (nrow(supported_piecewise_rows) > 0L) {
  supported_piecewise_rows <- supported_piecewise_rows %>%
    mutate(
      horizon_year = case_when(
        interval_label == "0-1" ~ 1,
        interval_label == "1-2" ~ 2,
        TRUE ~ horizon_year
      ),
      window_year = horizon_year
    )
}

intermediate_horizon_contrast <- convert_piecewise_rows_to_contrast(
  piecewise_site_effect = piecewise_site_effect,
  interval_labels = intermediate_interval_label,
  output_label_value = "intermediate_horizon_contrast",
  horizon_year_value = 5
) %>%
  mutate(
    reporting_status = make_reporting_status(estimable_flag, reported_estimate_suppressed_flag),
    plotting_status = make_plotting_status(estimable_flag, reported_estimate_suppressed_flag),
    plot_inclusion_flag = estimable_flag & !reported_estimate_suppressed_flag,
    reporting_note = make_reporting_note(
      reporting_status = reporting_status,
      availability_note = availability_note,
      support_issue_reason = support_issue_reason
    )
  )

tail_diagnostic_contrast <- convert_piecewise_rows_to_contrast(
  piecewise_site_effect = piecewise_site_effect,
  interval_labels = tail_interval_label,
  output_label_value = "tail_diagnostic_contrast",
  horizon_year_value = NA_real_
) %>%
  mutate(
    reporting_status = make_reporting_status(estimable_flag, reported_estimate_suppressed_flag),
    plotting_status = make_plotting_status(estimable_flag, reported_estimate_suppressed_flag),
    plot_inclusion_flag = estimable_flag & !reported_estimate_suppressed_flag,
    reporting_note = make_reporting_note(
      reporting_status = reporting_status,
      availability_note = availability_note,
      support_issue_reason = support_issue_reason
    )
  )

# 🔴 Assemble: supported short-horizon contrast block ===============================
restricted_model_outputs <- list()
reused_nonpiecewise_short <- FALSE
supported_short_nonpiecewise_rows <- NULL

if (use_existing_supported_short_nonpiecewise) {
  supported_short_nonpiecewise_rows <- harmonize_supported_short_nonpiecewise_cache(existing_supported_short, support_lookup)
  reused_nonpiecewise_short <- TRUE
  build_log <- append_build_log(
    build_log,
    component = "supported_short_nonpiecewise",
    status = "reused",
    source_file = existing_supported_short_file,
    detail = "Legacy short-horizon observed and restricted standardized rows were reused because the 1y/2y estimands remain valid under the revised specification."
  )
}

if (is.null(supported_short_nonpiecewise_rows)) {
  short_nonpiecewise_obj <- compute_supported_short_nonpiecewise_fresh(
    pnu_data = pnu_data,
    snu_data = snu_data,
    merged_data = merged_data,
    support_lookup = support_lookup,
    merged_site_formula_registry = merged_site_formula_registry,
    build_log = build_log
  )
  supported_short_nonpiecewise_rows <- short_nonpiecewise_obj$rows
  restricted_model_outputs <- short_nonpiecewise_obj$restricted_models
} else {
  restricted_model_outputs <- if (!is.null(existing_model_objects) && !is.null(existing_model_objects$restricted_models)) {
    existing_model_objects$restricted_models
  } else {
    list()
  }
}

supported_short_horizon_contrast <- bind_rows(
  supported_short_nonpiecewise_rows,
  supported_piecewise_rows
) %>%
  arrange(
    factor(result_type, levels = c("observed_risk", "site_standardized_risk", "absolute_risk_difference", "hazard_ratio")),
    horizon_year,
    model_name,
    dataset,
    interval_label
  ) %>%
  mutate(
    support_label = coalesce_character(support_label, support_tier),
    support_priority = coalesce_character(support_priority, make_support_priority(support_tier)),
    common_restricted_window_flag = dplyr::coalesce(common_restricted_window_flag, TRUE),
    output_label = "supported_short_horizon_contrast",
    reporting_status = make_reporting_status(estimable_flag, reported_estimate_suppressed_flag, bootstrap_instability_flag),
    plotting_status = make_plotting_status(estimable_flag, reported_estimate_suppressed_flag, bootstrap_instability_flag),
    plot_inclusion_flag = estimable_flag & !reported_estimate_suppressed_flag,
    reporting_note = make_reporting_note(
      reporting_status = reporting_status,
      availability_note = availability_note,
      bootstrap_note = bootstrap_note,
      support_issue_reason = support_issue_reason
    )
  )

short_validation_ok <- tryCatch({
  validate_supported_short_counts(
    supported_df = supported_short_horizon_contrast,
    risk_trajectory = risk_trajectory,
    analysis_totals_lookup = analysis_totals_lookup,
    piecewise_interval_support = piecewise_interval_support
  )
  TRUE
}, error = function(e) {
  if (!reused_nonpiecewise_short) stop(e)
  FALSE
})

if (!short_validation_ok) {
  short_nonpiecewise_obj <- compute_supported_short_nonpiecewise_fresh(
    pnu_data = pnu_data,
    snu_data = snu_data,
    merged_data = merged_data,
    support_lookup = support_lookup,
    merged_site_formula_registry = merged_site_formula_registry,
    build_log = build_log
  )
  supported_short_horizon_contrast <- bind_rows(
    short_nonpiecewise_obj$rows,
    supported_piecewise_rows
  ) %>%
    arrange(
      factor(result_type, levels = c("observed_risk", "site_standardized_risk", "absolute_risk_difference", "hazard_ratio")),
      horizon_year,
      model_name,
      dataset,
      interval_label
    ) %>%
    mutate(
      support_label = coalesce_character(support_label, support_tier),
      support_priority = coalesce_character(support_priority, make_support_priority(support_tier)),
      common_restricted_window_flag = dplyr::coalesce(common_restricted_window_flag, TRUE),
      output_label = "supported_short_horizon_contrast",
      reporting_status = make_reporting_status(estimable_flag, reported_estimate_suppressed_flag, bootstrap_instability_flag),
      plotting_status = make_plotting_status(estimable_flag, reported_estimate_suppressed_flag, bootstrap_instability_flag),
      plot_inclusion_flag = estimable_flag & !reported_estimate_suppressed_flag,
      reporting_note = make_reporting_note(
        reporting_status = reporting_status,
        availability_note = availability_note,
        bootstrap_note = bootstrap_note,
        support_issue_reason = support_issue_reason
      )
    )
  restricted_model_outputs <- short_nonpiecewise_obj$restricted_models
  validate_supported_short_counts(
    supported_df = supported_short_horizon_contrast,
    risk_trajectory = risk_trajectory,
    analysis_totals_lookup = analysis_totals_lookup,
    piecewise_interval_support = piecewise_interval_support
  )
}

validate_piecewise_interval_counts(
  interval_df = intermediate_horizon_contrast,
  piecewise_interval_support = piecewise_interval_support,
  analysis_totals_lookup = analysis_totals_lookup,
  interval_labels = intermediate_interval_label,
  object_label = "Intermediate interval"
)

validate_piecewise_interval_counts(
  interval_df = tail_diagnostic_contrast,
  piecewise_interval_support = piecewise_interval_support,
  analysis_totals_lookup = analysis_totals_lookup,
  interval_labels = tail_interval_label,
  object_label = "Tail diagnostic interval"
)

validate_reporting_flags(
  supported_df = supported_short_horizon_contrast,
  intermediate_df = intermediate_horizon_contrast,
  tail_df = tail_diagnostic_contrast,
  piecewise_df = piecewise_site_effect
)

build_log <- append_build_log(
  build_log,
  component = "interval_contrast_tables",
  status = "assembled",
  source_file = NA_character_,
  detail = "supported_short_horizon_contrast, intermediate_horizon_contrast, and tail_diagnostic_contrast were assembled under the revised Block 2 labeling contract."
)

# 🔴 Assemble: metadata and manifest registries ===============================
block2_metadata_registry <- tibble::tribble(
  ~metadata_group, ~metadata_name, ~metadata_value,
  "stage", "stage_name", "Block 2 timing-difference separation",
  "stage", "stage_role", "Separate timing difference from cure interpretation before later cure blocks.",
  "stage", "canonical_framework", "1.Model specifciation/spec.md",
  "stage", "block2_spec_version", block2_spec_version,
  "stage", "block2_caching_policy", block2_caching_policy,
  "stage", "block2_revision_note", "Revised Block 2 uses the interval backbone 0-1, 1-2, 2-5, and >5, preserves 1y/2y restricted-window contrasts, and saves both a plot-book PDF and separate PNG files.",
  "inputs", "dropbox_project_root", normalize_existing_path(dropbox_project_root),
  "stage", "main_interval_backbone", "0-1,1-2,2-5,>5",
  "stage", "tail_interval_rule", ">5 years is a projection-dominant tail diagnostic only and must not drive the primary timing interpretation.",
  "stage", "intermediate_interval_rule", "2-5 years is the main intermediate timing-separation interval and is secondary / partly model-dependent.",
  "inputs", "data_path", normalize_existing_path(data_path),
  "inputs", "export_path", normalize_existing_path(export_path),
  "inputs", "supporting_output_dir", normalize_existing_path(supporting_output_dir),
  "inputs", "legacy_stage4_export_path", normalize_existing_path(legacy_stage4_export_path),
  "inputs", "data_path_leaf", make_portable_file_reference(data_path),
  "inputs", "export_path_leaf", make_portable_file_reference(export_path),
  "inputs", "stage1_analysis_datasets_file", normalize_existing_path(stage1_analysis_datasets_file),
  "inputs", "stage1_dataset_registry_file", normalize_existing_path(stage1_dataset_registry_file),
  "inputs", "stage1_scaling_registry_file", normalize_existing_path(stage1_scaling_registry_file),
  "inputs", "stage1_metadata_registry_file", normalize_existing_path(stage1_metadata_registry_file),
  "inputs", "stage1_formula_registry_file", normalize_existing_path(stage1_formula_registry_file),
  "inputs", "stage1_horizon_registry_file", normalize_existing_path(stage1_horizon_registry_file),
  "inputs", "stage1_threshold_registry_file", normalize_existing_path(stage1_threshold_registry_file),
  "time", "analysis_time_variable", "days_followup",
  "time", "reporting_time_variable", "time_year = days_followup / 365.25",
  "time", "horizon_vector", paste(common_horizons_year, collapse = ","),
  "time", "piecewise_cuts_year", paste(piecewise_cuts_year, collapse = ","),
  "event", "event_definition", "status_num == 1",
  "event", "main_censoring_definition", "status_num %in% c(0, 2)",
  "timing", "supported_short_horizons_year", paste(supported_short_horizons_year, collapse = ","),
  "timing", "intermediate_interval_label", intermediate_interval_label,
  "timing", "tail_interval_label", tail_interval_label,
  "timing", "hazard_band_breaks_year", paste(hazard_band_breaks_year, collapse = ","),
  "timing", "threshold_vector", paste(format(risk_thresholds, trim = TRUE, scientific = FALSE), collapse = ","),
  "timing", "support_hierarchy_source", "Horizon-based support labels are inherited from Stage 1 for horizon-specific rows and mapped conservatively for interval-specific contrasts.",
  "bootstrap", "bootstrap_iterations", as.character(bootstrap_iterations),
  "bootstrap", "bootstrap_seed", as.character(bootstrap_seed),
  "bootstrap", "bootstrap_min_success_rate", format(bootstrap_min_success_rate, trim = TRUE),
  "bootstrap", "bootstrap_reporting_rule", "Rows with retained estimates but bootstrap_success_rate below the minimum stay reported and are flagged as bootstrap_unstable for plotting and interpretation.",
  "piecewise", "piecewise_min_events_per_site", as.character(piecewise_min_events_per_site),
  "piecewise", "piecewise_min_subjects_per_site", as.character(piecewise_min_subjects_per_site),
  "piecewise", "piecewise_nonestimable_rule", "If one site has too few subjects, sparse events, or no person-time in an interval, the interval-specific hazard ratio is suppressed.",
  "piecewise", "piecewise_export_rule", "Suppressed piecewise hazard ratios retain explicit support diagnostics in CSV outputs and are omitted from plotted point ranges.",
  "comparison", "short_output_label", "supported_short_horizon_contrast",
  "comparison", "intermediate_output_label", "intermediate_horizon_contrast",
  "comparison", "tail_output_label", "tail_diagnostic_contrast",
  "comparison", "analysis_subject_n_total_rule", "analysis_subject_n_total stores the full source-cohort size and is distinct from interval-specific at_risk_subject_n_total.",
  "comparison", "count_alignment_validation_rule", "Block 2 stops or falls back to recomputation if supported/intermediate/tail contrast tables do not align with risk_trajectory or piecewise_interval_support.",
  "comparison", "reuse_rule", "Existing Block 2 files are reused first; if absent, compatible legacy Stage 4 outputs are loaded from old/stage4_Timing-difference separation; otherwise the corresponding component is recomputed.",
  "outputs", "save_folder_rule", "Keep interpretive outputs and requested Stage 1 audit inputs at the top level and write auxiliary tables, logs, metadata, model objects, and standalone PNG files into sub_supporting.",
  "outputs", "top_level_navigation_files", "README.md",
  "outputs", "top_level_interpretation_files", "block2_supported_short_horizon_contrast.csv;block2_intermediate_horizon_contrast.csv;block2_tail_diagnostic_contrast.csv;block2_plot_book.pdf",
  "outputs", "top_level_audit_files", "stage1_analysis_datasets.rds;stage1_horizon_registry.csv;stage1_dataset_registry.csv",
  "outputs", "plot_png_rule", "Standalone PNG files are written into sub_supporting, while the combined PDF plot-book stays at the top level.",
  "outputs", "readme_rule", "Generate a top-level README.md that describes every exported file in both the top level and sub_supporting.",
  "outputs", "manifest_portability_rule", "block2_export_manifest.csv stores relative file_name / portable_file_reference paths from export_path so moved auxiliary files remain traceable.",
  "outputs", "site_interpretation_rule", "Interpret site as a proxy for broader treatment context, care pathway, selection, or follow-up structure rather than a clean treatment effect."
)

plot_filenames <- c(
  "block2_plot_risk_trajectory.png",
  "block2_plot_supported_short_horizon_contrast.png",
  "block2_plot_intermediate_horizon_contrast.png",
  "block2_plot_tail_diagnostic_contrast.png",
  "block2_plot_piecewise_interval_support.png",
  "block2_plot_hazard_pattern_summary.png"
)

block2_output_paths <- list(
  readme = file.path(export_path, "README.md"),
  dataset_summary = file.path(supporting_output_dir, "block2_dataset_summary.csv"),
  risk_trajectory = file.path(supporting_output_dir, "block2_risk_trajectory.csv"),
  hazard_pattern_summary = file.path(supporting_output_dir, "block2_hazard_pattern_summary.csv"),
  piecewise_interval_support = file.path(supporting_output_dir, "block2_piecewise_interval_support.csv"),
  piecewise_site_effect = file.path(supporting_output_dir, "block2_piecewise_site_effect.csv"),
  supported_short_horizon_contrast = file.path(export_path, "block2_supported_short_horizon_contrast.csv"),
  intermediate_horizon_contrast = file.path(export_path, "block2_intermediate_horizon_contrast.csv"),
  tail_diagnostic_contrast = file.path(export_path, "block2_tail_diagnostic_contrast.csv"),
  standard_error_registry = file.path(supporting_output_dir, "block2_standard_error_table_registry.csv"),
  standard_error_long = file.path(supporting_output_dir, "block2_standard_error_long.csv"),
  audit_stage1_analysis_datasets = file.path(export_path, "stage1_analysis_datasets.rds"),
  audit_stage1_horizon_registry = file.path(export_path, "stage1_horizon_registry.csv"),
  audit_stage1_dataset_registry = file.path(export_path, "stage1_dataset_registry.csv"),
  build_log = file.path(supporting_output_dir, "block2_build_log.csv"),
  metadata_registry = file.path(supporting_output_dir, "block2_metadata_registry.csv"),
  model_objects = file.path(supporting_output_dir, "block2_model_objects.rds"),
  plot_book = file.path(export_path, "block2_plot_book.pdf"),
  plot_risk_trajectory = file.path(supporting_output_dir, "block2_plot_risk_trajectory.png"),
  plot_supported_short = file.path(supporting_output_dir, "block2_plot_supported_short_horizon_contrast.png"),
  plot_intermediate = file.path(supporting_output_dir, "block2_plot_intermediate_horizon_contrast.png"),
  plot_tail = file.path(supporting_output_dir, "block2_plot_tail_diagnostic_contrast.png"),
  plot_piecewise_support = file.path(supporting_output_dir, "block2_plot_piecewise_interval_support.png"),
  plot_hazard_pattern = file.path(supporting_output_dir, "block2_plot_hazard_pattern_summary.png"),
  export_manifest = file.path(supporting_output_dir, "block2_export_manifest.csv")
)

block2_export_manifest <- tibble(
  object_name = c(
    "block2_readme",
    "block2_dataset_summary",
    "risk_trajectory",
    "hazard_pattern_summary",
    "piecewise_interval_support",
    "piecewise_site_effect",
    "supported_short_horizon_contrast",
    "intermediate_horizon_contrast",
    "tail_diagnostic_contrast",
    "block2_standard_error_registry",
    "block2_standard_error_long",
    "stage1_analysis_datasets_audit_copy",
    "stage1_horizon_registry_audit_copy",
    "stage1_dataset_registry_audit_copy",
    "build_log",
    "block2_metadata_registry",
    "block2_model_objects",
    "block2_plot_book",
    rep("block2_plot_png", length(plot_filenames)),
    "block2_export_manifest"
  ),
  description = c(
    "Human-readable guide describing every exported Block 2 file",
    "Dataset-level Block 2 summary aligned to the current Stage 1 backbone",
    "Observed KM risk trajectory on the annual 1-10 year grid with inherited support labels",
    "Hazard-pattern summary across 0-1, 1-2, 2-5, and 5-10 year bands",
    "Piecewise interval-support table for 0-1, 1-2, 2-5, and >5 years",
    "Piecewise merged site-effect hazard ratio table with interval-specific support diagnostics",
    "Supported short-horizon contrast block combining 1-year and 2-year common-window comparisons",
    "Intermediate-horizon contrast block for the 2-5 year interval",
    "Projection-dominant tail diagnostic contrast block for the >5 year interval",
    "Registry of Block 2 tables that contain standard-error or uncertainty columns",
    "Long-form Block 2 standard-error table assembled from the exported Block 2 CSV sources",
    "Top-level audit copy of the Stage 1 analysis datasets object used to build Block 2",
    "Top-level audit copy of the Stage 1 horizon registry used to inherit support hierarchy",
    "Top-level audit copy of the Stage 1 dataset registry used to trace dataset provenance",
    "Component-level build and reuse log",
    "Block 2 metadata registry",
    "Restricted and piecewise model objects bundled as RDS",
    "Single PDF plot book generated from exported data frames",
    rep("Standalone PNG version of a Block 2 plot generated from the same data-frame source as the PDF page", length(plot_filenames)),
    "Manifest of all Block 2 exported files"
  ),
  file_path = c(
    block2_output_paths$readme,
    block2_output_paths$dataset_summary,
    block2_output_paths$risk_trajectory,
    block2_output_paths$hazard_pattern_summary,
    block2_output_paths$piecewise_interval_support,
    block2_output_paths$piecewise_site_effect,
    block2_output_paths$supported_short_horizon_contrast,
    block2_output_paths$intermediate_horizon_contrast,
    block2_output_paths$tail_diagnostic_contrast,
    block2_output_paths$standard_error_registry,
    block2_output_paths$standard_error_long,
    block2_output_paths$audit_stage1_analysis_datasets,
    block2_output_paths$audit_stage1_horizon_registry,
    block2_output_paths$audit_stage1_dataset_registry,
    block2_output_paths$build_log,
    block2_output_paths$metadata_registry,
    block2_output_paths$model_objects,
    block2_output_paths$plot_book,
    block2_output_paths$plot_risk_trajectory,
    block2_output_paths$plot_supported_short,
    block2_output_paths$plot_intermediate,
    block2_output_paths$plot_tail,
    block2_output_paths$plot_piecewise_support,
    block2_output_paths$plot_hazard_pattern,
    block2_output_paths$export_manifest
  )
) %>%
  mutate(
    file_path = normalize_existing_path(file_path),
    file_name = vapply(file_path, make_relative_export_reference, character(1)),
    portable_file_reference = file_name
  ) %>%
  select(file_name, object_name, description, file_path, portable_file_reference)

block2_model_objects <- list(
  restricted_models = restricted_model_outputs,
  piecewise_models = piecewise_model_objects,
  piecewise_interval_support = piecewise_interval_support,
  stage1_support_lookup = support_lookup,
  analysis_totals_lookup = analysis_totals_lookup,
  config = list(
    data_path = data_path,
    export_path = export_path,
    supporting_output_dir = supporting_output_dir,
    bootstrap_iterations = bootstrap_iterations,
    bootstrap_seed = bootstrap_seed,
    bootstrap_min_success_rate = bootstrap_min_success_rate,
    piecewise_cuts_year = piecewise_cuts_year,
    piecewise_min_events_per_site = piecewise_min_events_per_site,
    piecewise_min_subjects_per_site = piecewise_min_subjects_per_site,
    hazard_band_breaks_year = hazard_band_breaks_year,
    supported_short_horizons_year = supported_short_horizons_year,
    intermediate_interval_label = intermediate_interval_label,
    tail_interval_label = tail_interval_label,
    thresholds_from_stage1 = risk_thresholds,
    horizons_from_stage1 = common_horizons_year,
    stage1_horizon_registry_file = stage1_horizon_registry_file,
    block2_spec_version = block2_spec_version,
    block2_caching_policy = block2_caching_policy
  )
)

block2_standard_error_bundle <- build_standard_error_export_bundle(list(
  dataset_summary = block2_dataset_summary,
  risk_trajectory = risk_trajectory,
  hazard_pattern_summary = hazard_pattern_summary,
  piecewise_interval_support = piecewise_interval_support,
  piecewise_site_effect = piecewise_site_effect,
  supported_short_horizon_contrast = supported_short_horizon_contrast,
  intermediate_horizon_contrast = intermediate_horizon_contrast,
  tail_diagnostic_contrast = tail_diagnostic_contrast
))

# 🔴 Export: tables, model bundle, and plot files ===============================
readr::write_csv(block2_dataset_summary, block2_output_paths$dataset_summary)
readr::write_csv(risk_trajectory, block2_output_paths$risk_trajectory)
readr::write_csv(hazard_pattern_summary, block2_output_paths$hazard_pattern_summary)
readr::write_csv(piecewise_interval_support, block2_output_paths$piecewise_interval_support)
readr::write_csv(piecewise_site_effect, block2_output_paths$piecewise_site_effect)
readr::write_csv(supported_short_horizon_contrast, block2_output_paths$supported_short_horizon_contrast)
readr::write_csv(intermediate_horizon_contrast, block2_output_paths$intermediate_horizon_contrast)
readr::write_csv(tail_diagnostic_contrast, block2_output_paths$tail_diagnostic_contrast)
readr::write_csv(block2_standard_error_bundle$registry, block2_output_paths$standard_error_registry)
readr::write_csv(block2_standard_error_bundle$long, block2_output_paths$standard_error_long)
if (file.exists(stage1_analysis_datasets_file)) {
  ok <- file.copy(stage1_analysis_datasets_file, block2_output_paths$audit_stage1_analysis_datasets, overwrite = TRUE)
  if (!isTRUE(ok)) {
    stop("Failed to copy stage1_analysis_datasets.rds into the Block 2 top-level audit bundle.", call. = FALSE)
  }
} else {
  saveRDS(stage1_inputs$analysis_datasets, block2_output_paths$audit_stage1_analysis_datasets)
}
if (file.exists(stage1_horizon_registry_file)) {
  ok <- file.copy(stage1_horizon_registry_file, block2_output_paths$audit_stage1_horizon_registry, overwrite = TRUE)
  if (!isTRUE(ok)) {
    stop("Failed to copy stage1_horizon_registry.csv into the Block 2 top-level audit bundle.", call. = FALSE)
  }
} else {
  readr::write_csv(stage1_inputs$horizon_registry, block2_output_paths$audit_stage1_horizon_registry)
}
if (file.exists(stage1_dataset_registry_file)) {
  ok <- file.copy(stage1_dataset_registry_file, block2_output_paths$audit_stage1_dataset_registry, overwrite = TRUE)
  if (!isTRUE(ok)) {
    stop("Failed to copy stage1_dataset_registry.csv into the Block 2 top-level audit bundle.", call. = FALSE)
  }
} else {
  readr::write_csv(stage1_inputs$dataset_registry, block2_output_paths$audit_stage1_dataset_registry)
}
readr::write_csv(build_log, block2_output_paths$build_log)
readr::write_csv(block2_metadata_registry, block2_output_paths$metadata_registry)
saveRDS(block2_model_objects, block2_output_paths$model_objects)

plot_risk_trajectory_obj <- plot_risk_trajectory(risk_trajectory)
plot_supported_short_obj <- plot_supported_short(supported_short_horizon_contrast)
plot_intermediate_obj <- plot_intermediate_interval(intermediate_horizon_contrast, piecewise_interval_support)
plot_tail_obj <- plot_tail_interval(tail_diagnostic_contrast, piecewise_interval_support)
plot_piecewise_support_obj <- plot_piecewise_support(piecewise_interval_support)
plot_hazard_pattern_obj <- plot_hazard_pattern(hazard_pattern_summary)

grDevices::pdf(block2_output_paths$plot_book, width = 11, height = 8.5, onefile = TRUE)
print(plot_risk_trajectory_obj)
print(plot_supported_short_obj)
print(plot_intermediate_obj)
print(plot_tail_obj)
print(plot_piecewise_support_obj)
print(plot_hazard_pattern_obj)
grDevices::dev.off()

save_plot_png(plot_risk_trajectory_obj, "block2_plot_risk_trajectory.png")
save_plot_png(plot_supported_short_obj, "block2_plot_supported_short_horizon_contrast.png")
save_plot_png(plot_intermediate_obj, "block2_plot_intermediate_horizon_contrast.png")
save_plot_png(plot_tail_obj, "block2_plot_tail_diagnostic_contrast.png")
save_plot_png(plot_piecewise_support_obj, "block2_plot_piecewise_interval_support.png")
save_plot_png(plot_hazard_pattern_obj, "block2_plot_hazard_pattern_summary.png")

write_block2_readme(block2_output_paths$readme, block2_export_manifest, block2_metadata_registry)
readr::write_csv(block2_export_manifest, block2_output_paths$export_manifest)
