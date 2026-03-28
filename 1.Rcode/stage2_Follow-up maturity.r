# 🔴 Configure: stage-two paths and maturity controls ===============================
# Edit the two paths below before running Stage 2.
data_path <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/stage1_Backbone lock'
export_path <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/stage2_Follow-up maturity'

reuse_existing_core_results <- TRUE
force_recompute_core_results <- FALSE
rebuild_plots_from_source_tables <- TRUE

plot_time_step_year <- 0.25
plot_max_year <- 10
risk_table_times_year <- c(0, 1:10)

completeness_sparse_risk_n <- 10L
completeness_very_sparse_risk_n <- 5L
completeness_low_fraction_cutoff <- 0.10

export_sex_panel_plots <- TRUE
site_dominance_warning_fraction <- 0.80
site_dominance_warning_min_horizon_year <- 3L
site_dominance_warning_exempt_primary_supported <- TRUE

current_patch_level <- 'stage2_followup_maturity_reuse_cache_png_v8'
scientific_compatibility_signature <- 'stage2_followup_maturity_science_v7'
compatible_reuse_patch_levels <- c(
  'patched_user_specified_paths_preferred_site_adjusted_reporting_v7',
  'stage2_followup_maturity_reuse_cache_png_v8'
)
compatible_science_signatures <- c(
  'stage2_followup_maturity_science_v7'
)

# 🔴 Initialize: packages and runtime state ===============================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(survival)
  library(scales)
})

options(stringsAsFactors = FALSE, scipen = 999)

# 🔴 Define: basic helpers ===============================
`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}

normalize_existing_directory <- function(path_label, path_value) {
  path_chr <- trimws(as.character(path_value))
  if (!nzchar(path_chr)) {
    stop(sprintf('`%s` must be a non-empty directory path.', path_label), call. = FALSE)
  }
  if (!dir.exists(path_chr)) {
    stop(sprintf('Directory does not exist for `%s`: %s', path_label, path_chr), call. = FALSE)
  }
  normalizePath(path_chr, winslash = '/', mustWork = TRUE)
}

normalize_output_directory <- function(path_label, path_value) {
  path_chr <- trimws(as.character(path_value))
  if (!nzchar(path_chr)) {
    stop(sprintf('`%s` must be a non-empty directory path.', path_label), call. = FALSE)
  }
  normalizePath(path_chr, winslash = '/', mustWork = FALSE)
}

collect_existing_files <- function(paths) {
  path_chr <- unique(as.character(paths))
  path_chr <- path_chr[!is.na(path_chr) & nzchar(path_chr) & file.exists(path_chr)]
  if (length(path_chr) == 0L) {
    return(character())
  }
  normalizePath(path_chr, winslash = '/', mustWork = TRUE)
}

compute_file_signature <- function(paths) {
  existing_paths <- collect_existing_files(paths)
  if (length(existing_paths) == 0L) {
    return(NA_character_)
  }

  md5_vec <- unname(tools::md5sum(existing_paths))
  paste(paste(basename(existing_paths), md5_vec, sep = '='), collapse = '|')
}


data_path <- normalize_existing_directory('data_path', data_path)
export_path <- normalize_output_directory('export_path', export_path)

dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(export_path)) {
  stop(sprintf('Failed to create the Stage 2 output directory: %s', export_path), call. = FALSE)
}

write_probe_file <- file.path(export_path, '.stage2_write_probe.tmp')
write_probe_ok <- tryCatch({
  probe_conn <- file(write_probe_file, open = 'wt')
  writeLines('stage2 write probe', probe_conn)
  close(probe_conn)
  TRUE
}, error = function(e) {
  FALSE
})
if (!isTRUE(write_probe_ok) || !file.exists(write_probe_file)) {
  stop(sprintf('The Stage 2 output directory is not writable: %s', export_path), call. = FALSE)
}
unlink(write_probe_file, force = TRUE)

message(sprintf('[Stage 2] Input directory   : %s', data_path))
message(sprintf('[Stage 2] Output directory  : %s', export_path))

dataset_order <- c('PNU', 'SNU', 'merged')

stage1_analysis_datasets_file <- file.path(data_path, 'stage1_analysis_datasets.rds')
stage1_backbone_bundle_file <- file.path(data_path, 'stage1_backbone_bundle.rds')
stage1_scaling_registry_file <- file.path(data_path, 'stage1_scaling_registry.csv')
stage1_horizon_registry_file <- file.path(data_path, 'stage1_horizon_registry.csv')
stage1_threshold_registry_file <- file.path(data_path, 'stage1_threshold_registry.csv')
stage1_metadata_registry_file <- file.path(data_path, 'stage1_metadata_registry.csv')
stage1_dataset_registry_file <- file.path(data_path, 'stage1_dataset_registry.csv')
stage1_formula_registry_file <- file.path(data_path, 'stage1_formula_registry.csv')
stage1_export_manifest_file <- file.path(data_path, 'stage1_export_manifest.csv')

current_stage1_source_signature <- compute_file_signature(c(
  stage1_analysis_datasets_file,
  stage1_backbone_bundle_file,
  stage1_scaling_registry_file,
  stage1_horizon_registry_file,
  stage1_threshold_registry_file,
  stage1_metadata_registry_file,
  stage1_dataset_registry_file,
  stage1_formula_registry_file,
  stage1_export_manifest_file
))

stage2_core_paths <- list(
  site_admin_lookup = file.path(export_path, 'stage2_site_admin_end_dates.csv'),
  reverse_km_summary = file.path(export_path, 'stage2_reverse_km_summary.csv'),
  numbers_at_risk = file.path(export_path, 'stage2_numbers_at_risk.csv'),
  composition = file.path(export_path, 'stage2_followup_composition.csv'),
  site_contribution = file.path(export_path, 'stage2_merged_site_contribution.csv'),
  horizon_summary = file.path(export_path, 'stage2_followup_horizon_summary.csv'),
  curve_plot_data = file.path(export_path, 'stage2_followup_curve_data.csv'),
  validation_summary = file.path(export_path, 'stage2_validation_summary.csv'),
  validation_issue_rows = file.path(export_path, 'stage2_validation_issue_rows.csv'),
  metadata_registry = file.path(export_path, 'stage2_metadata_registry.csv'),
  quality_flags = file.path(export_path, 'stage2_quality_flags.csv'),
  bundle = file.path(export_path, 'stage2_followup_bundle.rds')
)

assert_file_exists <- function(path_label, path_value) {
  if (!file.exists(path_value)) {
    stop(sprintf('Required input file is missing for `%s`: %s', path_label, path_value), call. = FALSE)
  }
  invisible(path_value)
}

read_csv_required <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

safe_read_csv <- function(path) {
  if (file.exists(path)) {
    read_csv_required(path)
  } else {
    NULL
  }
}

read_input_dataset <- function(path) {
  assert_file_exists('input_dataset', path)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c('csv', 'txt')) {
    return(
      readr::read_csv(
        path,
        col_types = readr::cols(.default = readr::col_character()),
        show_col_types = FALSE,
        progress = FALSE
      )
    )
  }
  if (ext == 'rds') {
    return(readRDS(path))
  }
  stop(sprintf('Unsupported input extension for `%s`.', path), call. = FALSE)
}

parse_date_flex <- function(x) {
  if (inherits(x, 'Date')) {
    return(x)
  }

  x_chr <- as.character(x)
  out <- as.Date(x_chr)

  if (anyNA(out)) {
    try_formats <- c('%Y-%m-%d', '%Y/%m/%d', '%m/%d/%Y', '%d/%m/%Y', '%Y%m%d')
    for (fmt in try_formats) {
      idx <- is.na(out) & !is.na(x_chr) & nzchar(x_chr)
      if (!any(idx)) {
        break
      }
      out[idx] <- as.Date(x_chr[idx], format = fmt)
    }
  }

  out
}

safe_divide <- function(num, den) {
  common_len <- max(length(num), length(den))
  num_vec <- rep_len(as.numeric(num), common_len)
  den_vec <- rep_len(as.numeric(den), common_len)

  out <- num_vec / den_vec
  out[is.na(den_vec) | den_vec == 0] <- NA_real_

  out
}

collapse_unique_chr <- function(x) {
  paste(sort(unique(as.character(x))), collapse = '|')
}

coalesce_character <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (length(v) > 0L && any(!is.na(v) & nzchar(as.character(v)))) {
      first_ok <- which(!is.na(v) & nzchar(as.character(v)))[1]
      return(as.character(v[[first_ok]]))
    }
  }
  NA_character_
}

first_time_below_survival <- function(time_vec, surv_vec, threshold = 0.5) {
  idx <- which(!is.na(surv_vec) & surv_vec <= threshold)
  if (length(idx) == 0L) {
    return(NA_real_)
  }
  as.numeric(time_vec[min(idx)])
}

vec_to_len <- function(x, n, fill = NA_real_) {
  out <- rep(fill, n)
  if (is.null(x)) {
    return(out)
  }
  x <- as.numeric(x)
  if (length(x) == 0L) {
    return(out)
  }
  copy_n <- min(length(x), n)
  out[seq_len(copy_n)] <- x[seq_len(copy_n)]
  out
}

get_named_numeric <- function(x, name) {
  if (is.null(x)) {
    return(NA_real_)
  }

  if (is.matrix(x) || is.data.frame(x)) {
    if (!is.null(colnames(x)) && name %in% colnames(x)) {
      return(as.numeric(x[1, name, drop = TRUE]))
    }
    return(NA_real_)
  }

  if (!is.null(names(x)) && name %in% names(x)) {
    return(as.numeric(x[[name]]))
  }

  NA_real_
}

factor_panel_labels <- function(x, panel_order) {
  factor(as.character(x), levels = unique(as.character(panel_order)))
}

safe_max_abs_diff <- function(x, y) {
  diff_vec <- abs(as.numeric(x) - as.numeric(y))
  if (length(diff_vec) == 0L || all(is.na(diff_vec))) {
    return(NA_real_)
  }
  max(diff_vec, na.rm = TRUE)
}

pull_first_or_na <- function(data, col_name) {
  if (nrow(data) == 0L || !col_name %in% names(data)) {
    return(NA_real_)
  }
  as.numeric(data[[col_name]][[1]])
}

empty_site_contribution_tbl <- function() {
  tibble(
    dataset = character(),
    subgroup_kind = character(),
    subgroup_value = character(),
    group_id = character(),
    panel_label = character(),
    analysis_structure = character(),
    horizon_year = integer(),
    interpretation_tier = character(),
    interpretation_note = character(),
    support_tier_standard = character(),
    support_tier = character(),
    support_display_group = character(),
    support_display_label = character(),
    horizon_evidence_class = character(),
    claim_restriction_flag = character(),
    site = character(),
    site_eligible_n = integer(),
    total_eligible_n = integer(),
    eligible_site_fraction = double(),
    eligible_site_count = integer(),
    dominant_site = character(),
    dominant_site_n = integer(),
    dominant_site_fraction = double(),
    site_dominance_status = character(),
    site_dominance_note = character()
  )
}

empty_quality_flags_tbl <- function() {
  tibble(
    flag_domain = character(),
    flag_severity = character(),
    flag_code = character(),
    dataset = character(),
    subgroup_kind = character(),
    subgroup_value = character(),
    group_id = character(),
    panel_label = character(),
    analysis_structure = character(),
    summary_view = character(),
    horizon_year = integer(),
    flag_text_value = character(),
    flag_numeric_value = double(),
    flag_message = character()
  )
}

# 🔴 Define: stage-one bundle readers ===============================
read_stage1_inputs <- function() {
  backbone_bundle <- if (file.exists(stage1_backbone_bundle_file)) readRDS(stage1_backbone_bundle_file) else NULL
  analysis_datasets <- if (file.exists(stage1_analysis_datasets_file)) {
    readRDS(stage1_analysis_datasets_file)
  } else if (!is.null(backbone_bundle) && !is.null(backbone_bundle$datasets)) {
    backbone_bundle$datasets
  } else {
    stop('Could not find `stage1_analysis_datasets.rds` or datasets inside `stage1_backbone_bundle.rds`.', call. = FALSE)
  }

  list(
    backbone_bundle = backbone_bundle,
    analysis_datasets = analysis_datasets,
    dataset_registry = safe_read_csv(stage1_dataset_registry_file),
    scaling_registry = safe_read_csv(stage1_scaling_registry_file),
    metadata_registry = safe_read_csv(stage1_metadata_registry_file),
    formula_registry = safe_read_csv(stage1_formula_registry_file),
    horizon_registry = safe_read_csv(stage1_horizon_registry_file),
    threshold_registry = safe_read_csv(stage1_threshold_registry_file),
    export_manifest = safe_read_csv(stage1_export_manifest_file)
  )
}

validate_stage1_inputs <- function(stage1_inputs) {
  required_datasets <- c('PNU', 'SNU', 'merged')

  if (!is.list(stage1_inputs$analysis_datasets) || !all(required_datasets %in% names(stage1_inputs$analysis_datasets))) {
    stop('Stage 1 analysis datasets must be an R list containing PNU, SNU, and merged.', call. = FALSE)
  }

  required_cols <- c(
    'id', 'site', 'unique_person_id', 'sex_num', 'age_exact_entry', 'age_s',
    'days_followup', 'time_year', 'status_num', 'event_main', 'censor_main',
    'right_censor_flag', 'remission_flag'
  )

  for (dataset_name in required_datasets) {
    df <- stage1_inputs$analysis_datasets[[dataset_name]]

    if (!all(required_cols %in% names(df))) {
      missing_cols <- setdiff(required_cols, names(df))
      stop(sprintf('[%s] Missing required Stage 1 columns: %s', dataset_name, paste(missing_cols, collapse = ', ')), call. = FALSE)
    }

    if (nrow(df) == 0L) {
      stop(sprintf('[%s] Stage 1 dataset contains zero rows.', dataset_name), call. = FALSE)
    }

    if (nrow(df) != dplyr::n_distinct(df$unique_person_id)) {
      stop(sprintf('[%s] `unique_person_id` must be unique.', dataset_name), call. = FALSE)
    }

    if (anyNA(df[, required_cols])) {
      stop(sprintf('[%s] Missing values detected in required Stage 1 columns.', dataset_name), call. = FALSE)
    }
  }

  if (is.null(stage1_inputs$horizon_registry) || !all(c('dataset', 'horizon_year') %in% names(stage1_inputs$horizon_registry))) {
    stop('Stage 1 horizon registry is missing required columns `dataset` and `horizon_year`.', call. = FALSE)
  }

  invisible(stage1_inputs)
}

infer_site_labels_from_stage1 <- function(stage1_inputs) {
  bundle_cfg <- stage1_inputs$backbone_bundle$config %||% list()
  pnu_site_label <- as.character(bundle_cfg$pnu_site_label %||% 'PNU')
  snu_site_label <- as.character(bundle_cfg$snu_site_label %||% 'SNU')
  c(pnu_site_label = pnu_site_label, snu_site_label = snu_site_label)
}

standardize_known_site_labels <- function(df, pnu_label, snu_label) {
  if (!'site' %in% names(df)) {
    stop('Input must contain column `site`.', call. = FALSE)
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

get_stage1_merged_source_file <- function(stage1_inputs) {
  bundle_cfg <- stage1_inputs$backbone_bundle$config %||% list()
  source_lookup <- stage1_inputs$backbone_bundle$source_lookup %||% list()

  candidate_paths <- c(
    bundle_cfg$merged_file %||% NA_character_,
    source_lookup$merged_file %||% NA_character_,
    source_lookup$merged_dataset_file %||% NA_character_
  )

  candidate_paths <- unique(candidate_paths[!is.na(candidate_paths) & nzchar(candidate_paths)])

  for (path in candidate_paths) {
    if (file.exists(path)) {
      return(normalizePath(path, winslash = '/', mustWork = TRUE))
    }
  }

  NA_character_
}

build_date_entry_lookup_from_raw <- function(stage1_inputs, pnu_label, snu_label) {
  raw_merged_file <- get_stage1_merged_source_file(stage1_inputs)

  if (is.na(raw_merged_file)) {
    stop(
      'Stage 1 datasets do not contain `date_entry`, and no valid `merged_file` path was found in `stage1_backbone_bundle.rds`.',
      call. = FALSE
    )
  }

  raw_df <- read_input_dataset(raw_merged_file)
  raw_df <- standardize_known_site_labels(raw_df, pnu_label = pnu_label, snu_label = snu_label)

  if (!all(c('id', 'site') %in% names(raw_df))) {
    stop('Raw merged file must contain `id` and `site` for date-entry recovery.', call. = FALSE)
  }

  date_entry_raw <- if ('date_entry' %in% names(raw_df)) {
    raw_df[['date_entry']]
  } else {
    out <- rep(NA_character_, nrow(raw_df))
    if ('mri' %in% names(raw_df)) {
      out[toupper(raw_df$site) == toupper(pnu_label)] <- as.character(raw_df[['mri']][toupper(raw_df$site) == toupper(pnu_label)])
    }
    if ('enter' %in% names(raw_df)) {
      out[toupper(raw_df$site) == toupper(snu_label)] <- as.character(raw_df[['enter']][toupper(raw_df$site) == toupper(snu_label)])
    }
    out
  }

  lookup <- raw_df %>%
    transmute(
      unique_person_id = paste(trimws(as.character(site)), trimws(as.character(id)), sep = '_'),
      date_entry = parse_date_flex(date_entry_raw)
    ) %>%
    distinct()

  if (nrow(lookup) != dplyr::n_distinct(lookup$unique_person_id)) {
    stop('Recovered date-entry lookup must be unique by `unique_person_id`.', call. = FALSE)
  }

  if (anyNA(lookup$date_entry)) {
    stop(
      sprintf(
        'Unable to recover `date_entry` for all subjects from raw merged file: %s',
        raw_merged_file
      ),
      call. = FALSE
    )
  }

  list(
    date_lookup = lookup,
    raw_merged_file = raw_merged_file
  )
}

ensure_date_entry_in_dataset <- function(df, dataset_name, date_lookup) {
  out <- tibble::as_tibble(df)

  if ('date_entry' %in% names(out)) {
    out <- out %>% mutate(date_entry = parse_date_flex(date_entry))
  } else {
    out <- out %>%
      left_join(date_lookup, by = 'unique_person_id')
  }

  if (anyNA(out$date_entry)) {
    stop(sprintf('[%s] Missing `date_entry` remained after Stage 1 recovery.', dataset_name), call. = FALSE)
  }

  out %>%
    mutate(
      date_entry = parse_date_flex(date_entry),
      followup_date = date_entry + round(as.numeric(days_followup))
    )
}

# 🔴 Define: dataset preparation and subgroup helpers ===============================
prepare_stage2_dataset <- function(df, dataset_name, site_mode = c('single', 'merged')) {
  site_mode <- match.arg(site_mode)

  required_cols <- c('site', 'id', 'sex_num', 'age_exact_entry', 'days_followup', 'status_num', 'date_entry')
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0L) {
    stop(
      sprintf('[%s] Missing required columns from Stage 1 dataset: %s', dataset_name, paste(missing_cols, collapse = ', ')),
      call. = FALSE
    )
  }

  out <- df %>%
    mutate(
      dataset = dataset_name,
      site = trimws(as.character(site)),
      id = trimws(as.character(id)),
      sex_num = as.integer(as.numeric(sex_num)),
      age_exact_entry = as.numeric(age_exact_entry),
      days_followup = as.numeric(days_followup),
      status_num = as.integer(as.numeric(status_num)),
      date_entry = parse_date_flex(date_entry),
      followup_date = date_entry + round(days_followup),
      unique_person_id = paste(site, id, sep = '_')
    )

  if (nrow(out) == 0L) {
    stop(sprintf('[%s] Dataset has zero rows after loading.', dataset_name), call. = FALSE)
  }

  if (any(vapply(out[required_cols], function(col) anyNA(col), logical(1)))) {
    stop(sprintf('[%s] Missing values detected in required Stage 2 columns.', dataset_name), call. = FALSE)
  }

  if (any(out$id == '')) {
    stop(sprintf('[%s] Blank `id` values detected.', dataset_name), call. = FALSE)
  }

  if (any(out$site == '')) {
    stop(sprintf('[%s] Blank `site` values detected.', dataset_name), call. = FALSE)
  }

  if (any(!out$sex_num %in% c(0L, 1L))) {
    stop(sprintf('[%s] `sex_num` must be coded as 0/1.', dataset_name), call. = FALSE)
  }

  if (any(!out$status_num %in% c(0L, 1L, 2L))) {
    stop(sprintf('[%s] `status_num` must be coded as 0/1/2.', dataset_name), call. = FALSE)
  }

  if (any(out$days_followup < 0)) {
    stop(sprintf('[%s] Negative `days_followup` values detected.', dataset_name), call. = FALSE)
  }

  if (anyNA(out$date_entry)) {
    stop(sprintf('[%s] Unable to parse all `date_entry` values as Date.', dataset_name), call. = FALSE)
  }

  if (anyNA(out$followup_date)) {
    stop(sprintf('[%s] Unable to derive `followup_date` for all subjects.', dataset_name), call. = FALSE)
  }

  if (nrow(out) != dplyr::n_distinct(out$unique_person_id)) {
    stop(sprintf('[%s] `site + id` must remain unique in Stage 2.', dataset_name), call. = FALSE)
  }

  n_site_levels <- dplyr::n_distinct(out$site)

  if (site_mode == 'single' && n_site_levels != 1L) {
    stop(sprintf('[%s] Single-cohort input must contain exactly one site level.', dataset_name), call. = FALSE)
  }

  if (site_mode == 'merged' && n_site_levels < 2L) {
    stop(sprintf('[%s] Merged input must contain at least two site levels.', dataset_name), call. = FALSE)
  }

  age_mean <- mean(out$age_exact_entry)
  age_sd <- stats::sd(out$age_exact_entry)

  if (is.na(age_sd) || age_sd <= 0) {
    stop(sprintf('[%s] `age_exact_entry` must have positive standard deviation.', dataset_name), call. = FALSE)
  }

  out <- out %>%
    mutate(
      time_year = days_followup / 365.25,
      event_main = as.integer(status_num == 1L),
      censor_main = as.integer(status_num %in% c(0L, 2L)),
      remission_flag = as.integer(status_num == 2L),
      right_censor_flag = as.integer(status_num == 0L),
      sex_label = factor(if_else(sex_num == 0L, 'Female', 'Male'), levels = c('Female', 'Male')),
      status_label = factor(
        case_when(
          status_num == 0L ~ 'right_censoring',
          status_num == 2L ~ 'remission',
          TRUE ~ 'transition'
        ),
        levels = c('right_censoring', 'remission', 'transition')
      ),
      age_s = (age_exact_entry - age_mean) / (2 * age_sd)
    )

  scaling_row <- tibble(
    dataset = dataset_name,
    variable = 'age_exact_entry',
    scaled_variable = 'age_s',
    center_mean = age_mean,
    scale_sd = age_sd,
    scale_two_sd = 2 * age_sd,
    scaling_rule = '(age_exact_entry - mean(age_exact_entry)) / (2 * sd(age_exact_entry))'
  )

  list(data = out, scaling = scaling_row)
}

compute_site_admin_lookup <- function(merged_df) {
  merged_df %>%
    group_by(site) %>%
    summarise(
      site_admin_end_date = max(followup_date, na.rm = TRUE),
      site_admin_end_days_followup_max = max(days_followup, na.rm = TRUE),
      site_admin_end_years_followup_max = max(time_year, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(site)
}

attach_admin_dates <- function(df, dataset_name, site_admin_lookup, merged_admin_end_date) {
  dataset_admin_end_date <- if (dataset_name == 'merged') {
    merged_admin_end_date
  } else {
    max(df$followup_date, na.rm = TRUE)
  }

  out <- df %>%
    left_join(site_admin_lookup, by = 'site') %>%
    mutate(
      site_admin_end_date = as.Date(site_admin_end_date),
      global_admin_end_date = as.Date(dataset_admin_end_date),
      site_specific_potential_followup_days = pmax(as.numeric(site_admin_end_date - date_entry), 0),
      site_specific_potential_followup_year = site_specific_potential_followup_days / 365.25,
      global_potential_followup_days = pmax(as.numeric(global_admin_end_date - date_entry), 0),
      global_potential_followup_year = global_potential_followup_days / 365.25
    )

  if (anyNA(out$site_admin_end_date)) {
    stop(sprintf('[%s] Every row must receive a site-specific administrative end date.', dataset_name), call. = FALSE)
  }

  out
}

build_site_support_dataset_lookup <- function(analysis_datasets_stage2) {
  single_dataset_names <- intersect(c('PNU', 'SNU'), names(analysis_datasets_stage2))

  lookup_rows <- lapply(single_dataset_names, function(dataset_name) {
    site_values <- analysis_datasets_stage2[[dataset_name]] %>%
      distinct(site) %>%
      pull(site) %>%
      as.character() %>%
      unique() %>%
      sort()

    if (length(site_values) != 1L) {
      stop(
        sprintf('[%s] Single-site dataset must contribute exactly one site label to the support lookup.', dataset_name),
        call. = FALSE
      )
    }

    tibble(
      site = site_values,
      support_dataset = dataset_name
    )
  })

  out <- bind_rows(lookup_rows) %>%
    distinct()

  if (nrow(out) != length(single_dataset_names)) {
    stop('Stage 2 site-support lookup must contain one unique mapping row for each single-site dataset.', call. = FALSE)
  }

  if (anyDuplicated(out$site)) {
    dup_sites <- unique(out$site[duplicated(out$site)])
    stop(
      sprintf('Stage 2 site-support lookup contains duplicated site labels: %s', paste(dup_sites, collapse = ', ')),
      call. = FALSE
    )
  }

  out %>%
    arrange(site)
}

resolve_support_dataset_for_group <- function(dataset_name, subgroup_kind, subgroup_value, site_support_dataset_lookup) {
  if (dataset_name == 'merged' && subgroup_kind == 'site') {
    matched_dataset <- site_support_dataset_lookup %>%
      filter(site == subgroup_value) %>%
      pull(support_dataset) %>%
      unique()

    if (length(matched_dataset) == 1L && !is.na(matched_dataset)) {
      return(as.character(matched_dataset))
    }

    if (length(matched_dataset) == 0L) {
      stop(
        sprintf('Merged site subgroup `%s` could not be matched to a single-site support dataset.', subgroup_value),
        call. = FALSE
      )
    }

    stop(
      sprintf(
        'Merged site subgroup `%s` matched multiple support datasets: %s',
        subgroup_value,
        paste(matched_dataset, collapse = ', ')
      ),
      call. = FALSE
    )
  }

  dataset_name
}

get_horizon_registry_for_group <- function(dataset_name, subgroup_kind, subgroup_value, horizon_registry, site_support_dataset_lookup) {
  support_dataset <- resolve_support_dataset_for_group(
    dataset_name = dataset_name,
    subgroup_kind = subgroup_kind,
    subgroup_value = subgroup_value,
    site_support_dataset_lookup = site_support_dataset_lookup
  )

  out <- horizon_registry %>%
    filter(dataset == support_dataset) %>%
    arrange(horizon_year)

  if (nrow(out) == 0L) {
    stop(
      sprintf(
        'No horizon registry rows were found for support dataset `%s` (requested by dataset=%s, subgroup_kind=%s, subgroup_value=%s).',
        support_dataset, dataset_name, subgroup_kind, subgroup_value
      ),
      call. = FALSE
    )
  }

  out
}

build_group_registry <- function(analysis_datasets) {
  group_rows <- list()

  for (dataset_name in names(analysis_datasets)) {
    df <- analysis_datasets[[dataset_name]]

    group_rows[[length(group_rows) + 1L]] <- tibble(
      dataset = dataset_name,
      subgroup_kind = 'overall',
      subgroup_value = 'Overall',
      group_id = paste(dataset_name, 'overall', 'Overall', sep = '__'),
      panel_label = paste0(dataset_name, ' | Overall')
    )

    sex_values <- df %>%
      distinct(sex_label) %>%
      pull(sex_label) %>%
      as.character() %>%
      sort()

    for (sex_value in sex_values) {
      group_rows[[length(group_rows) + 1L]] <- tibble(
        dataset = dataset_name,
        subgroup_kind = 'sex',
        subgroup_value = sex_value,
        group_id = paste(dataset_name, 'sex', sex_value, sep = '__'),
        panel_label = paste0(dataset_name, ' | Sex = ', sex_value)
      )
    }

    if (dataset_name == 'merged') {
      site_values <- df %>%
        distinct(site) %>%
        pull(site) %>%
        sort()

      for (site_value in site_values) {
        group_rows[[length(group_rows) + 1L]] <- tibble(
          dataset = dataset_name,
          subgroup_kind = 'site',
          subgroup_value = site_value,
          group_id = paste(dataset_name, 'site', site_value, sep = '__'),
          panel_label = paste0(dataset_name, ' | Site = ', site_value)
        )
      }
    }
  }

  bind_rows(group_rows) %>%
    mutate(dataset = factor(dataset, levels = dataset_order) %>% as.character())
}

filter_group_data <- function(df, subgroup_kind, subgroup_value) {
  if (subgroup_kind == 'overall') {
    return(df)
  }

  if (subgroup_kind == 'sex') {
    return(df %>% filter(as.character(sex_label) == subgroup_value))
  }

  if (subgroup_kind == 'site') {
    return(df %>% filter(site == subgroup_value))
  }

  stop(sprintf('Unsupported subgroup_kind: %s', subgroup_kind), call. = FALSE)
}

get_plot_panel_order <- function(group_registry, subgroup_kind_order) {
  group_registry %>%
    filter(subgroup_kind %in% subgroup_kind_order) %>%
    arrange(
      factor(dataset, levels = dataset_order),
      match(subgroup_kind, subgroup_kind_order),
      subgroup_value
    ) %>%
    pull(panel_label) %>%
    unique()
}

get_completeness_structures <- function(dataset_name, subgroup_kind) {
  if (dataset_name %in% c('PNU', 'SNU')) {
    return(
      tibble(
        analysis_structure = 'single_site',
        potential_followup_var = 'site_specific_potential_followup_year'
      )
    )
  }

  if (dataset_name == 'merged' && subgroup_kind %in% c('overall', 'sex')) {
    return(
      tibble(
        analysis_structure = c('site_free', 'site_adjusted'),
        potential_followup_var = c('global_potential_followup_year', 'site_specific_potential_followup_year')
      )
    )
  }

  if (dataset_name == 'merged' && subgroup_kind == 'site') {
    return(
      tibble(
        analysis_structure = 'site_adjusted',
        potential_followup_var = 'site_specific_potential_followup_year'
      )
    )
  }

  stop(
    sprintf('Unsupported completeness structure request for dataset=%s subgroup_kind=%s', dataset_name, subgroup_kind),
    call. = FALSE
  )
}

# 🔴 Define: horizon-tier standardization helpers ===============================
standard_to_stage_label <- function(x) {
  x_chr <- as.character(x)
  dplyr::case_when(
    x_chr == 'primary_supported' ~ 'primary-supported',
    x_chr == 'secondary' ~ 'secondary',
    x_chr == 'sensitivity' ~ 'sensitivity',
    x_chr == 'projection' ~ 'projection',
    TRUE ~ x_chr
  )
}

derive_default_support_tier <- function(dataset, horizon_year) {
  dataset_chr <- as.character(dataset)
  horizon_int <- as.integer(horizon_year)

  dplyr::case_when(
    dataset_chr == 'PNU' & horizon_int == 1L ~ 'primary_supported',
    dataset_chr == 'PNU' & horizon_int == 2L ~ 'sensitivity',
    dataset_chr == 'PNU' & horizon_int >= 3L ~ 'projection',
    dataset_chr %in% c('SNU', 'merged') & horizon_int %in% c(1L, 2L) ~ 'primary_supported',
    dataset_chr %in% c('SNU', 'merged') & horizon_int %in% c(3L, 4L, 5L) ~ 'secondary',
    dataset_chr %in% c('SNU', 'merged') & horizon_int >= 6L ~ 'projection',
    TRUE ~ 'sensitivity'
  )
}

derive_horizon_evidence_class <- function(dataset, horizon_year) {
  dataset_chr <- as.character(dataset)
  horizon_int <- as.integer(horizon_year)

  dplyr::case_when(
    dataset_chr == 'PNU' & horizon_int == 1L ~ 'directly_observed_data_supported',
    dataset_chr == 'PNU' & horizon_int == 2L ~ 'partly_model_dependent',
    dataset_chr == 'PNU' & horizon_int >= 3L ~ 'mostly_extrapolated',
    dataset_chr %in% c('SNU', 'merged') & horizon_int %in% c(1L, 2L) ~ 'directly_observed_data_supported',
    dataset_chr %in% c('SNU', 'merged') & horizon_int %in% c(3L, 4L, 5L) ~ 'partly_model_dependent',
    dataset_chr %in% c('SNU', 'merged') & horizon_int >= 6L ~ 'mostly_extrapolated',
    TRUE ~ 'partly_model_dependent'
  )
}

derive_claim_restriction_flag <- function(horizon_evidence_class, support_tier_standard) {
  dplyr::case_when(
    horizon_evidence_class == 'directly_observed_data_supported' ~ 'primary_claim_allowed',
    horizon_evidence_class == 'partly_model_dependent' ~ 'secondary_or_sensitivity_only',
    horizon_evidence_class == 'mostly_extrapolated' ~ 'projection_only',
    TRUE ~ 'secondary_or_sensitivity_only'
  )
}

canonicalize_support_tier <- function(x, dataset, horizon_year) {
  raw_chr <- trimws(as.character(x))
  normalized_chr <- tolower(gsub('[^a-z0-9]+', '_', raw_chr))
  default_tier <- derive_default_support_tier(dataset = dataset, horizon_year = horizon_year)

  out <- dplyr::case_when(
    is.na(raw_chr) | raw_chr == '' ~ default_tier,
    normalized_chr %in% c('primary', 'primary_support', 'primary_supported', 'primarysupported') ~ 'primary_supported',
    normalized_chr %in% c('secondary', 'secondary_support', 'secondary_supported', 'secondarysupported') ~ 'secondary',
    normalized_chr %in% c('sensitivity', 'sensitive', 'sensitivity_window', 'limited_support', 'limited_supported', 'intermediate') ~ 'sensitivity',
    normalized_chr %in% c('projection', 'projected', 'tail', 'tail_projection') ~ 'projection',
    TRUE ~ default_tier
  )

  as.character(out)
}

make_support_note <- function(support_tier_standard, horizon_evidence_class = NULL) {
  tier_chr <- as.character(support_tier_standard)
  evidence_chr <- if (is.null(horizon_evidence_class)) rep(NA_character_, length(tier_chr)) else as.character(horizon_evidence_class)

  dplyr::case_when(
    tier_chr == 'primary_supported' & evidence_chr == 'directly_observed_data_supported' ~ 'Primary supported horizon with comparatively direct follow-up support',
    tier_chr == 'secondary' ~ 'Secondary horizon with growing tail uncertainty',
    tier_chr == 'sensitivity' ~ 'Sensitivity or limited-support horizon; report cautiously',
    tier_chr == 'projection' ~ 'Projection horizon with substantial tail uncertainty',
    TRUE ~ 'Horizon-support classification not explicitly defined'
  )
}

make_support_display_group <- function(support_tier_standard, time_value = NULL) {
  tier_chr <- as.character(support_tier_standard)

  if (is.null(time_value)) {
    time_value <- rep(NA_real_, length(tier_chr))
  }

  dplyr::case_when(
    !is.na(time_value) & as.numeric(time_value) == 0 ~ 'baseline',
    tier_chr == 'primary_supported' ~ 'primary_supported',
    tier_chr %in% c('secondary', 'sensitivity') ~ 'secondary_or_sensitivity',
    tier_chr == 'projection' ~ 'projection',
    TRUE ~ 'other'
  )
}

make_support_display_label <- function(support_display_group) {
  group_chr <- as.character(support_display_group)

  dplyr::case_when(
    group_chr == 'baseline' ~ 'Baseline',
    group_chr == 'primary_supported' ~ 'Primary supported',
    group_chr == 'secondary_or_sensitivity' ~ 'Secondary / sensitivity',
    group_chr == 'projection' ~ 'Projection',
    TRUE ~ 'Other'
  )
}

normalize_horizon_registry <- function(horizon_registry) {
  required_cols <- c('dataset', 'horizon_year')
  missing_cols <- setdiff(required_cols, names(horizon_registry))

  if (length(missing_cols) > 0L) {
    stop(
      sprintf('Stage 1 horizon registry is missing required columns: %s', paste(missing_cols, collapse = ', ')),
      call. = FALSE
    )
  }

  out <- horizon_registry %>%
    mutate(
      dataset = as.character(dataset),
      horizon_year = as.integer(horizon_year)
    )

  raw_tier <- rep(NA_character_, nrow(out))
  if ('support_tier_standard' %in% names(out)) {
    raw_tier <- as.character(out$support_tier_standard)
  } else if ('support_tier' %in% names(out)) {
    raw_tier <- as.character(out$support_tier)
  } else if ('interpretation_tier' %in% names(out)) {
    raw_tier <- as.character(out$interpretation_tier)
  }

  support_tier_standard <- canonicalize_support_tier(
    x = raw_tier,
    dataset = out$dataset,
    horizon_year = out$horizon_year
  )

  evidence_class <- if ('horizon_evidence_class' %in% names(out)) {
    as.character(out$horizon_evidence_class)
  } else {
    derive_horizon_evidence_class(out$dataset, out$horizon_year)
  }

  claim_restriction_flag <- if ('claim_restriction_flag' %in% names(out)) {
    as.character(out$claim_restriction_flag)
  } else {
    derive_claim_restriction_flag(evidence_class, support_tier_standard)
  }

  out <- out %>%
    mutate(
      support_tier_standard = support_tier_standard,
      support_tier = support_tier_standard,
      interpretation_tier = standard_to_stage_label(support_tier_standard),
      horizon_evidence_class = evidence_class,
      claim_restriction_flag = claim_restriction_flag
    )

  if (!'interpretation_note' %in% names(out)) {
    out <- out %>%
      mutate(
        interpretation_note = make_support_note(support_tier_standard, horizon_evidence_class)
      )
    }
  out %>%
    mutate(
      support_display_group = make_support_display_group(
        support_tier_standard = support_tier_standard,
        time_value = horizon_year
      ),
      support_display_label = make_support_display_label(support_display_group),
      primary_supported_flag = support_tier_standard == 'primary_supported',
      projection_flag = support_tier_standard == 'projection'
    ) %>%
    select(
      dataset, horizon_year,
      interpretation_tier, interpretation_note,
      support_tier_standard, support_tier,
      support_display_group, support_display_label,
      horizon_evidence_class, claim_restriction_flag,
      primary_supported_flag, projection_flag
    ) %>%
    distinct() %>%
    arrange(factor(dataset, levels = dataset_order), horizon_year)
}

extract_threshold_vector <- function(threshold_registry) {
  threshold_col <- intersect(c('threshold', 'risk_threshold', 'threshold_value'), names(threshold_registry))
  if (length(threshold_col) == 0L) {
    stop('Stage 1 threshold registry must contain one of: threshold, risk_threshold, threshold_value.', call. = FALSE)
  }

  threshold_vector <- threshold_registry[[threshold_col[1]]] %>%
    as.numeric() %>%
    unique() %>%
    sort()

  if (length(threshold_vector) == 0L || anyNA(threshold_vector)) {
    stop('Threshold registry contains invalid threshold values.', call. = FALSE)
  }

  threshold_vector
}

make_reporting_preference_note <- function(dataset, subgroup_kind, analysis_structure) {
  dataset_chr <- as.character(dataset)
  subgroup_chr <- as.character(subgroup_kind)
  structure_chr <- as.character(analysis_structure)

  dplyr::case_when(
    dataset_chr == 'merged' & subgroup_chr %in% c('overall', 'sex') & structure_chr == 'site_adjusted' ~
      'Preferred reporting view for merged overall/sex completeness because it uses site-specific administrative end dates.',
    dataset_chr == 'merged' & subgroup_chr %in% c('overall', 'sex') & structure_chr == 'site_free' ~
      'Sensitivity-only completeness view for merged overall/sex summaries; do not cite as the main merged completeness result.',
    TRUE ~ 'Only available completeness structure for this panel; acceptable for reporting.'
  )
}

add_reporting_preference_fields <- function(horizon_summary) {
  horizon_summary %>%
    mutate(
      preferred_for_reporting = dplyr::case_when(
        dataset == 'merged' & subgroup_kind %in% c('overall', 'sex') & analysis_structure == 'site_adjusted' ~ TRUE,
        dataset == 'merged' & subgroup_kind %in% c('overall', 'sex') & analysis_structure == 'site_free' ~ FALSE,
        TRUE ~ TRUE
      ),
      reporting_priority = dplyr::case_when(
        preferred_for_reporting ~ 1L,
        TRUE ~ 2L
      ),
      reporting_preference_note = make_reporting_preference_note(
        dataset = dataset,
        subgroup_kind = subgroup_kind,
        analysis_structure = analysis_structure
      )
    )
}


# 🔴 Define: survival-curve summarizers ===============================
safe_survfit <- function(time_year, event_indicator) {
  if (length(time_year) == 0L) {
    stop('safe_survfit received an empty time vector.', call. = FALSE)
  }

  if (anyNA(time_year) || anyNA(event_indicator)) {
    stop('safe_survfit requires non-missing time and event vectors.', call. = FALSE)
  }

  if (any(time_year < 0)) {
    stop('safe_survfit requires non-negative follow-up time.', call. = FALSE)
  }

  if (any(!event_indicator %in% c(0, 1))) {
    stop('safe_survfit requires a 0/1 event indicator.', call. = FALSE)
  }

  survival::survfit(survival::Surv(time_year, event_indicator) ~ 1, conf.type = 'log-log')
}

summarise_survfit_at_times <- function(fit, time_grid_year, curve_role, summary_view) {
  sm <- summary(fit, times = time_grid_year, extend = TRUE)
  n_time <- length(time_grid_year)

  tibble(
    time_year = as.numeric(time_grid_year),
    estimate = vec_to_len(sm$surv, n_time),
    conf_low = vec_to_len(sm$lower, n_time),
    conf_high = vec_to_len(sm$upper, n_time),
    std_error = vec_to_len(sm$std.err, n_time),
    n_risk = vec_to_len(sm$n.risk, n_time),
    n_event = vec_to_len(sm$n.event, n_time),
    n_censor = vec_to_len(sm$n.censor, n_time),
    curve_role = curve_role,
    summary_view = summary_view
  )
}

extract_reverse_km_table <- function(fit, summary_view) {
  fit_summary <- tryCatch(summary(fit), error = function(e) NULL)
  fit_table <- if (!is.null(fit_summary) && !is.null(fit_summary$table)) fit_summary$table else NULL

  records <- get_named_numeric(fit_table, 'records')
  if (is.na(records) && !is.null(fit$n)) {
    records <- as.numeric(fit$n)
  }

  events <- get_named_numeric(fit_table, 'events')
  if (is.na(events) && !is.null(fit$n.event)) {
    events <- sum(as.numeric(fit$n.event), na.rm = TRUE)
  }

  median_followup <- get_named_numeric(fit_table, 'median')
  if (is.na(median_followup)) {
    median_followup <- first_time_below_survival(fit$time, fit$surv, threshold = 0.5)
  }

  median_lcl <- get_named_numeric(fit_table, '0.95LCL')
  if (is.na(median_lcl) && !is.null(fit$lower)) {
    median_lcl <- first_time_below_survival(fit$time, fit$lower, threshold = 0.5)
  }

  median_ucl <- get_named_numeric(fit_table, '0.95UCL')
  if (is.na(median_ucl) && !is.null(fit$upper)) {
    median_ucl <- first_time_below_survival(fit$time, fit$upper, threshold = 0.5)
  }

  tibble(
    summary_view = summary_view,
    reverse_km_records = records,
    reverse_km_events = events,
    reverse_km_median_followup_year = median_followup,
    reverse_km_median_followup_lcl_year = median_lcl,
    reverse_km_median_followup_ucl_year = median_ucl
  )
}

compute_curve_bundle <- function(df, dataset_name, subgroup_kind, subgroup_value, group_id, panel_label, plot_time_grid_year) {
  transition_event <- as.integer(df$status_num == 1L)
  reverse_event_transition_only <- as.integer(df$status_num %in% c(0L, 2L))
  event_any_terminal <- as.integer(df$status_num %in% c(1L, 2L))
  reverse_event_remission_sensitive <- as.integer(df$status_num == 0L)

  km_transition_fit <- safe_survfit(df$time_year, transition_event)
  km_any_terminal_fit <- safe_survfit(df$time_year, event_any_terminal)
  rev_transition_fit <- safe_survfit(df$time_year, reverse_event_transition_only)
  rev_remission_fit <- safe_survfit(df$time_year, reverse_event_remission_sensitive)

  curve_df <- bind_rows(
    summarise_survfit_at_times(km_transition_fit, plot_time_grid_year, 'km_curve', 'transition_only'),
    summarise_survfit_at_times(km_any_terminal_fit, plot_time_grid_year, 'km_curve', 'remission_sensitive'),
    summarise_survfit_at_times(rev_transition_fit, plot_time_grid_year, 'reverse_km_curve', 'transition_only'),
    summarise_survfit_at_times(rev_remission_fit, plot_time_grid_year, 'reverse_km_curve', 'remission_sensitive')
  ) %>%
    mutate(
      dataset = dataset_name,
      subgroup_kind = subgroup_kind,
      subgroup_value = subgroup_value,
      group_id = group_id,
      panel_label = panel_label,
      n_total = nrow(df)
    ) %>%
    relocate(dataset, subgroup_kind, subgroup_value, group_id, panel_label, summary_view, curve_role)

  reverse_summary_df <- bind_rows(
    extract_reverse_km_table(rev_transition_fit, 'transition_only'),
    extract_reverse_km_table(rev_remission_fit, 'remission_sensitive')
  ) %>%
    mutate(
      dataset = dataset_name,
      subgroup_kind = subgroup_kind,
      subgroup_value = subgroup_value,
      group_id = group_id,
      panel_label = panel_label,
      n_total = nrow(df),
      max_observed_followup_year = max(df$time_year, na.rm = TRUE),
      max_site_specific_potential_followup_year = max(df$site_specific_potential_followup_year, na.rm = TRUE),
      max_global_potential_followup_year = max(df$global_potential_followup_year, na.rm = TRUE)
    ) %>%
    relocate(dataset, subgroup_kind, subgroup_value, group_id, panel_label, summary_view)

  list(curves = curve_df, reverse_summary = reverse_summary_df)
}

# 🔴 Define: horizon-summary builders ===============================
compute_horizon_row <- function(
    df,
    dataset_name,
    subgroup_kind,
    subgroup_value,
    group_id,
    panel_label,
    horizon_year,
    interpretation_tier,
    interpretation_note,
    support_tier_standard,
    support_tier,
    support_display_group,
    support_display_label,
    horizon_evidence_class,
    claim_restriction_flag,
    analysis_structure,
    potential_followup_var,
    summary_view
) {
  tol <- 1e-10
  n_total <- nrow(df)
  risk_set_n <- sum(df$time_year >= (horizon_year - tol), na.rm = TRUE)

  cumulative_transition_n <- sum(df$status_num == 1L & df$time_year <= (horizon_year + tol), na.rm = TRUE)
  cumulative_remission_n <- sum(df$status_num == 2L & df$time_year <= (horizon_year + tol), na.rm = TRUE)
  cumulative_right_censoring_n <- sum(df$status_num == 0L & df$time_year <= (horizon_year + tol), na.rm = TRUE)

  eligible_df <- df %>%
    filter(.data[[potential_followup_var]] >= (horizon_year - tol))

  eligible_n <- nrow(eligible_df)

  if (eligible_n == 0L) {
    return(
      tibble(
        dataset = dataset_name,
        subgroup_kind = subgroup_kind,
        subgroup_value = subgroup_value,
        group_id = group_id,
        panel_label = panel_label,
        analysis_structure = analysis_structure,
        summary_view = summary_view,
        horizon_year = horizon_year,
        interpretation_tier = interpretation_tier,
        interpretation_note = interpretation_note,
        support_tier_standard = support_tier_standard,
        support_tier = support_tier,
        support_display_group = support_display_group,
        support_display_label = support_display_label,
        horizon_evidence_class = horizon_evidence_class,
        claim_restriction_flag = claim_restriction_flag,
        n_total = n_total,
        eligible_n = 0L,
        eligibility_fraction = 0,
        risk_set_n = risk_set_n,
        risk_set_fraction = safe_divide(risk_set_n, n_total),
        cumulative_transition_n = cumulative_transition_n,
        cumulative_remission_n = cumulative_remission_n,
        cumulative_right_censoring_n = cumulative_right_censoring_n,
        cumulative_censoring_n = if (summary_view == 'transition_only') cumulative_remission_n + cumulative_right_censoring_n else cumulative_right_censoring_n,
        cumulative_eventlike_n = if (summary_view == 'transition_only') cumulative_transition_n else cumulative_transition_n + cumulative_remission_n,
        observed_person_time_year = NA_real_,
        percentage_method = NA_real_,
        CCI = NA_real_,
        SPT = NA_real_,
        cumulative_transition_prop = safe_divide(cumulative_transition_n, n_total),
        cumulative_remission_prop = safe_divide(cumulative_remission_n, n_total),
        cumulative_right_censoring_prop = safe_divide(cumulative_right_censoring_n, n_total),
        cumulative_censoring_prop = safe_divide(if (summary_view == 'transition_only') cumulative_remission_n + cumulative_right_censoring_n else cumulative_right_censoring_n, n_total),
        cumulative_eventlike_prop = safe_divide(if (summary_view == 'transition_only') cumulative_transition_n else cumulative_transition_n + cumulative_remission_n, n_total),
        instability_marker = 'no_eligible_subjects'
      )
    )
  }

  obs <- eligible_df$time_year

  if (summary_view == 'transition_only') {
    event_obs <- eligible_df$status_num == 1L & obs <= (horizon_year + tol)
    early_censor <- eligible_df$status_num %in% c(0L, 2L) & obs < (horizon_year - tol)
  } else {
    event_obs <- eligible_df$status_num %in% c(1L, 2L) & obs <= (horizon_year + tol)
    early_censor <- eligible_df$status_num == 0L & obs < (horizon_year - tol)
  }

  observed_complete <- event_obs | obs >= (horizon_year - tol)

  observed_person_time_year <- sum(pmin(obs, horizon_year), na.rm = TRUE)
  percentage_method <- safe_divide(sum(observed_complete, na.rm = TRUE), eligible_n)
  cci_denom <- sum(ifelse(event_obs, pmin(obs, horizon_year), horizon_year), na.rm = TRUE)
  spt_numer <- sum(ifelse(early_censor, pmin(obs, horizon_year), horizon_year), na.rm = TRUE)

  instability_marker <- dplyr::case_when(
    risk_set_n == 0L ~ 'no_risk_set',
    risk_set_n < completeness_very_sparse_risk_n ~ 'very_sparse',
    risk_set_n < completeness_sparse_risk_n ~ 'sparse',
    risk_set_n / n_total < completeness_low_fraction_cutoff ~ 'low_fraction',
    TRUE ~ 'stable_enough'
  )

  tibble(
    dataset = dataset_name,
    subgroup_kind = subgroup_kind,
    subgroup_value = subgroup_value,
    group_id = group_id,
    panel_label = panel_label,
    analysis_structure = analysis_structure,
    summary_view = summary_view,
    horizon_year = horizon_year,
    interpretation_tier = interpretation_tier,
    interpretation_note = interpretation_note,
    support_tier_standard = support_tier_standard,
    support_tier = support_tier,
    support_display_group = support_display_group,
    support_display_label = support_display_label,
    horizon_evidence_class = horizon_evidence_class,
    claim_restriction_flag = claim_restriction_flag,
    n_total = n_total,
    eligible_n = eligible_n,
    eligibility_fraction = safe_divide(eligible_n, n_total),
    risk_set_n = risk_set_n,
    risk_set_fraction = safe_divide(risk_set_n, n_total),
    cumulative_transition_n = cumulative_transition_n,
    cumulative_remission_n = cumulative_remission_n,
    cumulative_right_censoring_n = cumulative_right_censoring_n,
    cumulative_censoring_n = if (summary_view == 'transition_only') cumulative_remission_n + cumulative_right_censoring_n else cumulative_right_censoring_n,
    cumulative_eventlike_n = if (summary_view == 'transition_only') cumulative_transition_n else cumulative_transition_n + cumulative_remission_n,
    observed_person_time_year = observed_person_time_year,
    percentage_method = percentage_method,
    CCI = safe_divide(observed_person_time_year, cci_denom),
    SPT = safe_divide(spt_numer, eligible_n * horizon_year),
    cumulative_transition_prop = safe_divide(cumulative_transition_n, n_total),
    cumulative_remission_prop = safe_divide(cumulative_remission_n, n_total),
    cumulative_right_censoring_prop = safe_divide(cumulative_right_censoring_n, n_total),
    cumulative_censoring_prop = safe_divide(if (summary_view == 'transition_only') cumulative_remission_n + cumulative_right_censoring_n else cumulative_right_censoring_n, n_total),
    cumulative_eventlike_prop = safe_divide(if (summary_view == 'transition_only') cumulative_transition_n else cumulative_transition_n + cumulative_remission_n, n_total),
    instability_marker = instability_marker
  )
}

compute_horizon_summary <- function(df, dataset_name, subgroup_kind, subgroup_value, group_id, panel_label, horizon_registry_dataset) {
  if (nrow(horizon_registry_dataset) == 0L) {
    stop(sprintf('[%s] Horizon registry is empty for this dataset.', dataset_name), call. = FALSE)
  }

  views <- c('transition_only', 'remission_sensitive')
  structure_registry <- get_completeness_structures(dataset_name = dataset_name, subgroup_kind = subgroup_kind)

  bind_rows(lapply(views, function(view_name) {
    bind_rows(lapply(seq_len(nrow(structure_registry)), function(s) {
      bind_rows(lapply(seq_len(nrow(horizon_registry_dataset)), function(i) {
        compute_horizon_row(
          df = df,
          dataset_name = dataset_name,
          subgroup_kind = subgroup_kind,
          subgroup_value = subgroup_value,
          group_id = group_id,
          panel_label = panel_label,
          horizon_year = horizon_registry_dataset$horizon_year[i],
          interpretation_tier = horizon_registry_dataset$interpretation_tier[i],
          interpretation_note = horizon_registry_dataset$interpretation_note[i],
          support_tier_standard = horizon_registry_dataset$support_tier_standard[i],
          support_tier = horizon_registry_dataset$support_tier[i],
          support_display_group = horizon_registry_dataset$support_display_group[i],
          support_display_label = horizon_registry_dataset$support_display_label[i],
          horizon_evidence_class = horizon_registry_dataset$horizon_evidence_class[i],
          claim_restriction_flag = horizon_registry_dataset$claim_restriction_flag[i],
          analysis_structure = structure_registry$analysis_structure[s],
          potential_followup_var = structure_registry$potential_followup_var[s],
          summary_view = view_name
        )
      }))
    }))
  }))
}

add_completeness_display_fields <- function(horizon_summary) {
  horizon_summary %>%
    mutate(
      completeness_plot_status = case_when(
        eligible_n <= 0 ~ 'masked_no_eligible',
        risk_set_n <= 0 ~ 'masked_no_risk_set',
        risk_set_n < completeness_very_sparse_risk_n ~ 'very_sparse',
        risk_set_n < completeness_sparse_risk_n ~ 'sparse',
        risk_set_fraction < completeness_low_fraction_cutoff ~ 'low_fraction',
        TRUE ~ 'stable'
      ),
      completeness_plot_note = case_when(
        completeness_plot_status == 'masked_no_eligible' ~ 'Masked in completeness plot because no eligible subjects remain for this horizon.',
        completeness_plot_status == 'masked_no_risk_set' ~ 'Masked in completeness plot because the risk set is zero at this horizon.',
        completeness_plot_status == 'very_sparse' ~ 'Displayed with caution: very sparse risk set.',
        completeness_plot_status == 'sparse' ~ 'Displayed with caution: sparse risk set.',
        completeness_plot_status == 'low_fraction' ~ 'Displayed with caution: low remaining risk-set fraction.',
        TRUE ~ 'Displayed without masking.'
      ),
      completeness_plot_alpha = case_when(
        completeness_plot_status == 'stable' ~ 1.00,
        completeness_plot_status == 'low_fraction' ~ 0.90,
        completeness_plot_status == 'sparse' ~ 0.80,
        completeness_plot_status == 'very_sparse' ~ 0.65,
        TRUE ~ 0.00
      ),
      completeness_plot_masked = completeness_plot_status %in% c('masked_no_eligible', 'masked_no_risk_set'),
      percentage_method_plot = if_else(completeness_plot_masked, NA_real_, percentage_method),
      CCI_plot = if_else(completeness_plot_masked, NA_real_, CCI),
      SPT_plot = if_else(completeness_plot_masked, NA_real_, SPT)
    )
}


# 🔴 Define: composition, risk-table, and site-contribution summaries ===============================
compute_composition_over_time <- function(df, dataset_name, subgroup_kind, subgroup_value, group_id, panel_label, time_grid_year) {
  n_total <- nrow(df)

  if (n_total <= 0L) {
    stop(sprintf('[%s] Cannot compute composition for an empty dataset.', dataset_name), call. = FALSE)
  }

  views <- c('transition_only', 'remission_sensitive')
  all_state_levels <- c('transition', 'remission', 'censoring', 'right_censoring', 'still_under_followup')

  build_one_time <- function(t_year, view_name, state_levels) {
    if (view_name == 'transition_only') {
      state_vec <- dplyr::case_when(
        df$status_num == 1L & df$time_year <= t_year ~ 'transition',
        df$status_num %in% c(0L, 2L) & df$time_year <= t_year ~ 'censoring',
        TRUE ~ 'still_under_followup'
      )
    } else {
      state_vec <- dplyr::case_when(
        df$status_num == 1L & df$time_year <= t_year ~ 'transition',
        df$status_num == 2L & df$time_year <= t_year ~ 'remission',
        df$status_num == 0L & df$time_year <= t_year ~ 'right_censoring',
        TRUE ~ 'still_under_followup'
      )
    }

    state_vec <- as.character(state_vec)
    state_tab <- table(factor(state_vec, levels = state_levels))

    tibble(
      dataset = dataset_name,
      subgroup_kind = subgroup_kind,
      subgroup_value = subgroup_value,
      group_id = group_id,
      panel_label = panel_label,
      summary_view = view_name,
      time_year = rep(as.numeric(t_year), length(state_levels)),
      n_total = rep(n_total, length(state_levels)),
      state = factor(state_levels, levels = all_state_levels),
      state_n = as.integer(state_tab),
      state_prop = as.numeric(state_tab) / n_total
    )
  }

  bind_rows(lapply(views, function(view_name) {
    state_levels <- if (view_name == 'transition_only') {
      c('transition', 'censoring', 'still_under_followup')
    } else {
      c('transition', 'remission', 'right_censoring', 'still_under_followup')
    }

    bind_rows(lapply(time_grid_year, function(t_year) {
      build_one_time(
        t_year = t_year,
        view_name = view_name,
        state_levels = state_levels
      )
    }))
  })) %>%
    arrange(summary_view, time_year, state)
}

compute_numbers_at_risk <- function(df, dataset_name, subgroup_kind, subgroup_value, group_id, panel_label, risk_table_times_year, horizon_registry_dataset) {
  if (nrow(horizon_registry_dataset) == 0L) {
    stop(sprintf('[%s] Horizon registry is empty for risk-table computation.', dataset_name), call. = FALSE)
  }

  matched_idx <- match(as.integer(risk_table_times_year), horizon_registry_dataset$horizon_year)
  matched_tier <- horizon_registry_dataset$interpretation_tier[matched_idx]
  matched_standard <- horizon_registry_dataset$support_tier_standard[matched_idx]
  matched_support <- horizon_registry_dataset$support_tier[matched_idx]
  matched_group <- horizon_registry_dataset$support_display_group[matched_idx]
  matched_label <- horizon_registry_dataset$support_display_label[matched_idx]
  matched_evidence <- horizon_registry_dataset$horizon_evidence_class[matched_idx]
  matched_claim <- horizon_registry_dataset$claim_restriction_flag[matched_idx]

  tibble(
    dataset = dataset_name,
    subgroup_kind = subgroup_kind,
    subgroup_value = subgroup_value,
    group_id = group_id,
    panel_label = panel_label,
    time_year = as.numeric(risk_table_times_year)
  ) %>%
    mutate(
      n_risk = purrr::map_int(time_year, ~ sum(df$time_year >= (.x - 1e-10), na.rm = TRUE)),
      n_total = nrow(df),
      n_risk_fraction = safe_divide(n_risk, n_total),
      interpretation_tier = if_else(time_year == 0, 'baseline', dplyr::coalesce(matched_tier, 'unregistered')),
      support_tier_standard = if_else(time_year == 0, 'baseline', dplyr::coalesce(matched_standard, 'unregistered')),
      support_tier = if_else(time_year == 0, 'baseline', dplyr::coalesce(matched_support, 'unregistered')),
      support_display_group = if_else(time_year == 0, 'baseline', dplyr::coalesce(matched_group, 'other')),
      support_display_label = if_else(time_year == 0, 'Baseline', dplyr::coalesce(matched_label, 'Other')),
      horizon_evidence_class = if_else(time_year == 0, 'baseline', dplyr::coalesce(matched_evidence, 'partly_model_dependent')),
      claim_restriction_flag = if_else(time_year == 0, 'primary_claim_allowed', dplyr::coalesce(matched_claim, 'secondary_or_sensitivity_only')),
      risk_table_label = format(as.integer(n_risk), trim = TRUE, scientific = FALSE)
    )
}

should_raise_site_dominance_warning <- function(site_dominance_status, support_tier_standard, horizon_year) {
  status_chr <- as.character(site_dominance_status)
  support_chr <- as.character(support_tier_standard)
  horizon_int <- as.integer(horizon_year)

  dplyr::case_when(
    is.na(status_chr) ~ FALSE,
    status_chr %in% c('not_applicable', 'not_computed', 'multi_site_balanced', 'no_eligible_subjects') ~ FALSE,
    status_chr == 'single_site_only' ~ TRUE,
    status_chr == 'single_site_dominant' &
      !is.na(horizon_int) &
      horizon_int >= site_dominance_warning_min_horizon_year &
      !(isTRUE(site_dominance_warning_exempt_primary_supported) & support_chr == 'primary_supported') ~ TRUE,
    TRUE ~ FALSE
  )
}

make_site_dominance_warning_note <- function(site_dominance_status, support_tier_standard, horizon_year) {
  status_chr <- as.character(site_dominance_status)
  support_chr <- as.character(support_tier_standard)
  horizon_int <- as.integer(horizon_year)

  dplyr::case_when(
    is.na(status_chr) ~ 'Merged site-dominance warning was not evaluated.',
    status_chr == 'not_applicable' ~ 'Merged site-dominance warning is not applicable outside merged overall/sex summaries.',
    status_chr == 'not_computed' ~ 'Merged site-dominance warning could not be evaluated because site contribution totals were unavailable.',
    status_chr == 'no_eligible_subjects' ~ 'No merged site-dominance warning because no eligible subjects remain for this horizon.',
    status_chr == 'multi_site_balanced' ~ 'No merged site-dominance warning because eligible subjects remain distributed across multiple sites.',
    status_chr == 'single_site_only' ~ 'Merged site-dominance warning retained because only one site contributes eligible subjects at this horizon.',
    status_chr == 'single_site_dominant' & !is.na(horizon_int) & horizon_int < site_dominance_warning_min_horizon_year ~
      sprintf('Descriptive only: dominant-site share is not escalated to a warning before year %d.', site_dominance_warning_min_horizon_year),
    status_chr == 'single_site_dominant' & isTRUE(site_dominance_warning_exempt_primary_supported) & support_chr == 'primary_supported' ~
      'Descriptive only: dominant-site share is not escalated to a warning on primary-supported horizons.',
    status_chr == 'single_site_dominant' ~
      'Merged site-dominance warning retained because one site dominates the eligible set at a late, non-primary-supported horizon.',
    TRUE ~ 'Merged site-dominance warning was not triggered.'
  )
}

compute_site_contribution_profile <- function(
    df,
    dataset_name,
    subgroup_kind,
    subgroup_value,
    group_id,
    panel_label,
    horizon_registry_dataset
) {
  if (!(dataset_name == 'merged' && subgroup_kind %in% c('overall', 'sex'))) {
    return(empty_site_contribution_tbl())
  }

  structure_registry <- get_completeness_structures(dataset_name = dataset_name, subgroup_kind = subgroup_kind)
  site_levels <- sort(unique(as.character(df$site)))
  tol <- 1e-10

  bind_rows(lapply(seq_len(nrow(structure_registry)), function(s) {
    analysis_structure <- structure_registry$analysis_structure[s]
    potential_followup_var <- structure_registry$potential_followup_var[s]

    bind_rows(lapply(seq_len(nrow(horizon_registry_dataset)), function(i) {
      horizon_year <- horizon_registry_dataset$horizon_year[i]

      site_counts <- df %>%
        mutate(is_eligible = .data[[potential_followup_var]] >= (horizon_year - tol)) %>%
        group_by(site) %>%
        summarise(site_eligible_n = sum(is_eligible, na.rm = TRUE), .groups = 'drop') %>%
        right_join(tibble(site = site_levels), by = 'site') %>%
        mutate(site_eligible_n = dplyr::coalesce(as.integer(site_eligible_n), 0L)) %>%
        arrange(site)

      total_eligible_n <- sum(site_counts$site_eligible_n, na.rm = TRUE)
      eligible_site_count <- sum(site_counts$site_eligible_n > 0L, na.rm = TRUE)

      if (total_eligible_n > 0L) {
        dominant_row <- site_counts %>%
          arrange(desc(site_eligible_n), site) %>%
          slice(1L)
        dominant_site <- dominant_row$site[[1]]
        dominant_site_n <- as.integer(dominant_row$site_eligible_n[[1]])
      } else {
        dominant_site <- NA_character_
        dominant_site_n <- NA_integer_
      }

      dominant_site_fraction <- safe_divide(dominant_site_n, total_eligible_n)

      site_dominance_status <- dplyr::case_when(
        total_eligible_n <= 0L ~ 'no_eligible_subjects',
        eligible_site_count <= 1L ~ 'single_site_only',
        !is.na(dominant_site_fraction) & dominant_site_fraction >= site_dominance_warning_fraction ~ 'single_site_dominant',
        TRUE ~ 'multi_site_balanced'
      )

      site_dominance_note <- dplyr::case_when(
        site_dominance_status == 'no_eligible_subjects' ~ 'No eligible subjects remain for this merged completeness horizon.',
        site_dominance_status == 'single_site_only' ~ 'Eligible subjects come from only one site at this merged horizon.',
        site_dominance_status == 'single_site_dominant' ~ sprintf('Eligible subjects are dominated by one site (>= %.0f%% of the eligible set).', 100 * site_dominance_warning_fraction),
        TRUE ~ 'Eligible subjects remain distributed across multiple sites without a dominant-site warning.'
      )

      site_counts %>%
        mutate(
          dataset = dataset_name,
          subgroup_kind = subgroup_kind,
          subgroup_value = subgroup_value,
          group_id = group_id,
          panel_label = panel_label,
          analysis_structure = analysis_structure,
          horizon_year = horizon_year,
          interpretation_tier = horizon_registry_dataset$interpretation_tier[i],
          interpretation_note = horizon_registry_dataset$interpretation_note[i],
          support_tier_standard = horizon_registry_dataset$support_tier_standard[i],
          support_tier = horizon_registry_dataset$support_tier[i],
          support_display_group = horizon_registry_dataset$support_display_group[i],
          support_display_label = horizon_registry_dataset$support_display_label[i],
          horizon_evidence_class = horizon_registry_dataset$horizon_evidence_class[i],
          claim_restriction_flag = horizon_registry_dataset$claim_restriction_flag[i],
          total_eligible_n = as.integer(total_eligible_n),
          eligible_site_fraction = safe_divide(site_eligible_n, total_eligible_n),
          eligible_site_count = as.integer(eligible_site_count),
          dominant_site = dominant_site,
          dominant_site_n = dominant_site_n,
          dominant_site_fraction = dominant_site_fraction,
          site_dominance_status = site_dominance_status,
          site_dominance_note = site_dominance_note
        ) %>%
        select(
          dataset, subgroup_kind, subgroup_value, group_id, panel_label,
          analysis_structure, horizon_year,
          interpretation_tier, interpretation_note,
          support_tier_standard, support_tier,
          support_display_group, support_display_label,
          horizon_evidence_class, claim_restriction_flag,
          site, site_eligible_n, total_eligible_n, eligible_site_fraction,
          eligible_site_count, dominant_site, dominant_site_n, dominant_site_fraction,
          site_dominance_status, site_dominance_note
        )
    }))
  }))
}

summarise_site_contribution_profile <- function(site_contribution_data) {
  if (nrow(site_contribution_data) == 0L) {
    return(
      tibble(
        dataset = character(),
        subgroup_kind = character(),
        subgroup_value = character(),
        group_id = character(),
        panel_label = character(),
        analysis_structure = character(),
        horizon_year = integer(),
        total_eligible_n = integer(),
        eligible_site_count = integer(),
        dominant_site = character(),
        dominant_site_n = integer(),
        dominant_site_fraction = double(),
        site_dominance_status = character(),
        site_dominance_note = character()
      )
    )
  }

  site_contribution_data %>%
    group_by(
      dataset, subgroup_kind, subgroup_value, group_id, panel_label,
      analysis_structure, horizon_year
    ) %>%
    summarise(
      total_eligible_n = first(total_eligible_n),
      eligible_site_count = first(eligible_site_count),
      dominant_site = first(dominant_site),
      dominant_site_n = first(dominant_site_n),
      dominant_site_fraction = first(dominant_site_fraction),
      site_dominance_status = first(site_dominance_status),
      site_dominance_note = first(site_dominance_note),
      .groups = 'drop'
    )
}

attach_site_contribution_to_horizon <- function(horizon_summary, site_contribution_summary) {
  out <- horizon_summary %>%
    left_join(
      site_contribution_summary,
      by = c(
        'dataset', 'subgroup_kind', 'subgroup_value', 'group_id', 'panel_label',
        'analysis_structure', 'horizon_year'
      )
    ) %>%
    mutate(
      site_dominance_status = case_when(
        dataset == 'merged' & subgroup_kind %in% c('overall', 'sex') & analysis_structure %in% c('site_free', 'site_adjusted') ~ dplyr::coalesce(site_dominance_status, 'not_computed'),
        TRUE ~ 'not_applicable'
      ),
      site_dominance_note = case_when(
        dataset == 'merged' & subgroup_kind %in% c('overall', 'sex') & analysis_structure %in% c('site_free', 'site_adjusted') ~ dplyr::coalesce(site_dominance_note, 'Merged site contribution summary was not available for this horizon.'),
        TRUE ~ 'Not applicable outside merged overall/sex completeness summaries.'
      ),
      site_dominance_warning_flag = should_raise_site_dominance_warning(
        site_dominance_status = site_dominance_status,
        support_tier_standard = support_tier_standard,
        horizon_year = horizon_year
      ),
      site_dominance_warning_note = make_site_dominance_warning_note(
        site_dominance_status = site_dominance_status,
        support_tier_standard = support_tier_standard,
        horizon_year = horizon_year
      )
    )

  out %>%
    mutate(
      eligible_site_count = if_else(site_dominance_status == 'not_applicable', NA_integer_, as.integer(eligible_site_count)),
      dominant_site_n = if_else(site_dominance_status == 'not_applicable', NA_integer_, as.integer(dominant_site_n)),
      dominant_site_fraction = if_else(site_dominance_status == 'not_applicable', NA_real_, as.numeric(dominant_site_fraction))
    )
}

# 🔴 Define: reshaping helpers ===============================
make_horizon_long <- function(horizon_summary) {
  horizon_summary %>%
    pivot_longer(
      cols = c(
        eligible_n, eligibility_fraction,
        risk_set_n, risk_set_fraction,
        cumulative_transition_n, cumulative_remission_n, cumulative_right_censoring_n,
        cumulative_censoring_n, cumulative_eventlike_n,
        observed_person_time_year,
        percentage_method, CCI, SPT,
        cumulative_transition_prop, cumulative_remission_prop, cumulative_right_censoring_prop,
        cumulative_censoring_prop, cumulative_eventlike_prop
      ),
      names_to = 'metric_name',
      values_to = 'metric_value'
    ) %>%
    mutate(metric_scope = 'horizon_level') %>%
    arrange(dataset, subgroup_kind, subgroup_value, analysis_structure, summary_view, horizon_year, metric_name)
}

make_side_by_side_summary <- function(horizon_summary) {
  horizon_summary %>%
    select(
      dataset, subgroup_kind, subgroup_value, group_id, panel_label,
      analysis_structure, preferred_for_reporting, reporting_priority, reporting_preference_note,
      horizon_year,
      interpretation_tier, interpretation_note,
      support_tier_standard, support_tier,
      support_display_group, support_display_label,
      horizon_evidence_class, claim_restriction_flag,
      completeness_plot_status, completeness_plot_note,
      eligible_site_count, dominant_site, dominant_site_n, dominant_site_fraction,
      site_dominance_status, site_dominance_note, site_dominance_warning_flag, site_dominance_warning_note,
      summary_view,
      eligible_n, risk_set_n,
      cumulative_transition_n, cumulative_remission_n, cumulative_right_censoring_n,
      cumulative_censoring_n, cumulative_eventlike_n,
      percentage_method, CCI, SPT, instability_marker
    ) %>%
    pivot_wider(
      names_from = summary_view,
      values_from = c(
        eligible_n, risk_set_n,
        cumulative_transition_n, cumulative_remission_n, cumulative_right_censoring_n,
        cumulative_censoring_n, cumulative_eventlike_n,
        percentage_method, CCI, SPT, instability_marker
      ),
      names_sep = '__'
    ) %>%
    arrange(factor(dataset, levels = dataset_order), subgroup_kind, subgroup_value, reporting_priority, analysis_structure, horizon_year)
}

build_completeness_plot_data <- function(horizon_summary) {
  horizon_summary %>%
    select(
      dataset, subgroup_kind, subgroup_value, group_id, panel_label,
      analysis_structure, preferred_for_reporting, reporting_priority, reporting_preference_note,
      summary_view, horizon_year,
      interpretation_tier, interpretation_note,
      support_tier_standard, support_tier,
      support_display_group, support_display_label,
      horizon_evidence_class, claim_restriction_flag,
      completeness_plot_status, completeness_plot_note, completeness_plot_alpha,
      eligible_n, risk_set_n, risk_set_fraction,
      percentage_method_plot, CCI_plot, SPT_plot,
      site_dominance_status, site_dominance_warning_flag, site_dominance_warning_note
    ) %>%
    pivot_longer(
      cols = c(percentage_method_plot, CCI_plot, SPT_plot),
      names_to = 'completeness_metric',
      values_to = 'completeness_value'
    ) %>%
    mutate(
      completeness_metric = recode(
        completeness_metric,
        percentage_method_plot = 'percentage_method',
        CCI_plot = 'CCI',
        SPT_plot = 'SPT'
      )
    ) %>%
    arrange(factor(dataset, levels = dataset_order), subgroup_kind, subgroup_value, reporting_priority, analysis_structure, summary_view, horizon_year, completeness_metric)
}

# 🔴 Define: validation audits and quality flags ===============================
audit_numbers_at_risk_monotonic <- function(numbers_at_risk_data) {
  detail <- numbers_at_risk_data %>%
    arrange(dataset, subgroup_kind, subgroup_value, time_year) %>%
    group_by(dataset, subgroup_kind, subgroup_value) %>%
    summarise(
      monotone_nonincreasing = all(diff(n_risk) <= 0),
      .groups = 'drop'
    ) %>%
    filter(!monotone_nonincreasing) %>%
    mutate(check_name = 'numbers_at_risk_monotonic')

  summary <- tibble(
    validation_group = 'internal_consistency',
    check_name = 'numbers_at_risk_monotonic',
    status = if (nrow(detail) == 0L) 'pass' else 'fail',
    n_problem_rows = nrow(detail),
    check_note = 'Numbers-at-risk must be monotone non-increasing within each dataset/subgroup series.'
  )

  list(summary = summary, detail = detail)
}

audit_composition_totals <- function(composition_data) {
  detail <- composition_data %>%
    group_by(dataset, subgroup_kind, subgroup_value, summary_view, time_year) %>%
    summarise(state_prop_sum = sum(state_prop), .groups = 'drop') %>%
    filter(abs(state_prop_sum - 1) > 1e-8) %>%
    mutate(check_name = 'composition_totals')

  summary <- tibble(
    validation_group = 'internal_consistency',
    check_name = 'composition_totals',
    status = if (nrow(detail) == 0L) 'pass' else 'fail',
    n_problem_rows = nrow(detail),
    check_note = 'Composition proportions must sum to 1 at every dataset/subgroup/view/time combination.'
  )

  list(summary = summary, detail = detail)
}

audit_horizon_monotonicity <- function(horizon_summary) {
  detail <- horizon_summary %>%
    arrange(dataset, subgroup_kind, subgroup_value, analysis_structure, summary_view, horizon_year) %>%
    group_by(dataset, subgroup_kind, subgroup_value, analysis_structure, summary_view) %>%
    summarise(
      risk_set_nonincreasing = all(diff(risk_set_n) <= 0),
      eligible_nonincreasing = all(diff(eligible_n) <= 0),
      transition_nondecreasing = all(diff(cumulative_transition_n) >= 0),
      remission_nondecreasing = all(diff(cumulative_remission_n) >= 0),
      censor_nondecreasing = all(diff(cumulative_censoring_n) >= 0),
      eventlike_nondecreasing = all(diff(cumulative_eventlike_n) >= 0),
      .groups = 'drop'
    ) %>%
    filter(
      !risk_set_nonincreasing |
        !eligible_nonincreasing |
        !transition_nondecreasing |
        !remission_nondecreasing |
        !censor_nondecreasing |
        !eventlike_nondecreasing
    ) %>%
    mutate(check_name = 'horizon_monotonicity')

  summary <- tibble(
    validation_group = 'internal_consistency',
    check_name = 'horizon_monotonicity',
    status = if (nrow(detail) == 0L) 'pass' else 'fail',
    n_problem_rows = nrow(detail),
    check_note = 'Risk-set and eligibility counts must decline with horizon; cumulative event/censoring counts must not decrease.'
  )

  list(summary = summary, detail = detail)
}

audit_required_stage2_views <- function(horizon_summary) {
  detail <- horizon_summary %>%
    distinct(dataset, subgroup_kind, subgroup_value, analysis_structure, summary_view) %>%
    group_by(dataset, subgroup_kind, subgroup_value, analysis_structure) %>%
    summarise(
      n_views = n_distinct(summary_view),
      .groups = 'drop'
    ) %>%
    filter(n_views != 2L) %>%
    mutate(check_name = 'required_stage2_views')

  summary <- tibble(
    validation_group = 'internal_consistency',
    check_name = 'required_stage2_views',
    status = if (nrow(detail) == 0L) 'pass' else 'fail',
    n_problem_rows = nrow(detail),
    check_note = 'Every analysis_structure must include both transition_only and remission_sensitive summaries.'
  )

  list(summary = summary, detail = detail)
}

audit_merged_site_consistency <- function(horizon_summary, reverse_km_summary) {
  site_map <- c('PNU', 'SNU')

  compare_cols_horizon <- c(
    'n_total', 'eligible_n', 'risk_set_n',
    'cumulative_transition_n', 'cumulative_remission_n', 'cumulative_right_censoring_n',
    'cumulative_censoring_n', 'cumulative_eventlike_n',
    'observed_person_time_year', 'percentage_method', 'CCI', 'SPT'
  )

  compare_cols_reverse <- c(
    'reverse_km_records', 'reverse_km_events',
    'reverse_km_median_followup_year', 'reverse_km_median_followup_lcl_year',
    'reverse_km_median_followup_ucl_year', 'n_total',
    'max_observed_followup_year',
    'max_site_specific_potential_followup_year'
  )

  detail_horizon <- bind_rows(lapply(site_map, function(site_name) {
    single_df <- horizon_summary %>%
      filter(
        dataset == site_name,
        subgroup_kind == 'overall',
        subgroup_value == 'Overall',
        analysis_structure == 'single_site'
      ) %>%
      arrange(summary_view, horizon_year)

    merged_df <- horizon_summary %>%
      filter(
        dataset == 'merged',
        subgroup_kind == 'site',
        subgroup_value == site_name,
        analysis_structure == 'site_adjusted'
      ) %>%
      arrange(summary_view, horizon_year)

    if (nrow(single_df) != nrow(merged_df) || nrow(single_df) == 0L) {
      return(tibble(site = site_name, compare_block = 'horizon_summary', compare_column = 'row_alignment', max_abs_diff = NA_real_, issue_type = 'row_mismatch'))
    }

    tibble(
      site = site_name,
      compare_block = 'horizon_summary',
      compare_column = compare_cols_horizon,
      max_abs_diff = map_dbl(compare_cols_horizon, function(col_name) safe_max_abs_diff(single_df[[col_name]], merged_df[[col_name]])),
      issue_type = 'value_difference'
    ) %>% filter(is.na(max_abs_diff) | max_abs_diff > 1e-8)
  })) %>% mutate(check_name = 'merged_site_consistency_horizon')

  detail_reverse <- bind_rows(lapply(site_map, function(site_name) {
    single_rev <- reverse_km_summary %>%
      filter(dataset == site_name, subgroup_kind == 'overall', subgroup_value == 'Overall') %>%
      arrange(summary_view)

    merged_rev <- reverse_km_summary %>%
      filter(dataset == 'merged', subgroup_kind == 'site', subgroup_value == site_name) %>%
      arrange(summary_view)

    if (nrow(single_rev) != nrow(merged_rev) || nrow(single_rev) == 0L) {
      return(tibble(site = site_name, compare_block = 'reverse_km_summary', compare_column = 'row_alignment', max_abs_diff = NA_real_, issue_type = 'row_mismatch'))
    }

    tibble(
      site = site_name,
      compare_block = 'reverse_km_summary',
      compare_column = compare_cols_reverse,
      max_abs_diff = map_dbl(compare_cols_reverse, function(col_name) safe_max_abs_diff(single_rev[[col_name]], merged_rev[[col_name]])),
      issue_type = 'value_difference'
    ) %>% filter(is.na(max_abs_diff) | max_abs_diff > 1e-8)
  })) %>% mutate(check_name = 'merged_site_consistency_reverse_km')

  summary <- bind_rows(
    tibble(validation_group = 'cross_dataset_consistency', check_name = 'merged_site_consistency_horizon', status = if (nrow(detail_horizon) == 0L) 'pass' else 'fail', n_problem_rows = nrow(detail_horizon), check_note = 'Merged site_adjusted site-specific horizon summaries must match the corresponding single-site overall summaries.'),
    tibble(validation_group = 'cross_dataset_consistency', check_name = 'merged_site_consistency_reverse_km', status = if (nrow(detail_reverse) == 0L) 'pass' else 'fail', n_problem_rows = nrow(detail_reverse), check_note = 'Merged site-specific reverse-KM summaries must match the corresponding single-site overall reverse-KM summaries.')
  )

  list(summary = summary, detail = bind_rows(detail_horizon, detail_reverse))
}

audit_site_contribution_alignment <- function(horizon_summary, site_contribution_summary) {
  target_rows <- horizon_summary %>%
    filter(dataset == 'merged', subgroup_kind %in% c('overall', 'sex'), analysis_structure %in% c('site_free', 'site_adjusted')) %>%
    distinct(dataset, subgroup_kind, subgroup_value, group_id, panel_label, analysis_structure, horizon_year, eligible_n)

  detail <- target_rows %>%
    left_join(
      site_contribution_summary %>%
        select(dataset, subgroup_kind, subgroup_value, group_id, panel_label, analysis_structure, horizon_year, total_eligible_n),
      by = c('dataset', 'subgroup_kind', 'subgroup_value', 'group_id', 'panel_label', 'analysis_structure', 'horizon_year')
    ) %>%
    mutate(eligible_difference = abs(as.numeric(eligible_n) - as.numeric(total_eligible_n))) %>%
    filter(is.na(total_eligible_n) | is.na(eligible_n) | eligible_difference > 1e-8) %>%
    mutate(check_name = 'site_contribution_alignment')

  summary <- tibble(
    validation_group = 'internal_consistency',
    check_name = 'site_contribution_alignment',
    status = if (nrow(detail) == 0L) 'pass' else 'fail',
    n_problem_rows = nrow(detail),
    check_note = 'Merged site-contribution totals must reproduce eligible_n in the merged overall/sex completeness summaries.'
  )

  list(summary = summary, detail = detail)
}

compare_support_label_frames <- function(reference_df, comparison_df, key_col, compare_cols, site_value, support_dataset, compare_block) {
  ref_key <- as.numeric(reference_df[[key_col]])
  cmp_key <- as.numeric(comparison_df[[key_col]])

  if (nrow(reference_df) != nrow(comparison_df) || nrow(reference_df) == 0L || !identical(ref_key, cmp_key)) {
    return(tibble(site = site_value, support_dataset = support_dataset, compare_block = compare_block, grid_value = NA_real_, compare_column = 'row_alignment', reference_value = NA_character_, merged_value = NA_character_, issue_type = 'row_mismatch'))
  }

  bind_rows(lapply(compare_cols, function(col_name) {
    ref_chr <- as.character(reference_df[[col_name]])
    cmp_chr <- as.character(comparison_df[[col_name]])
    ref_cmp <- ifelse(is.na(ref_chr), '__NA__', ref_chr)
    cmp_cmp <- ifelse(is.na(cmp_chr), '__NA__', cmp_chr)
    mismatch_idx <- which(ref_cmp != cmp_cmp)

    if (length(mismatch_idx) == 0L) {
      return(tibble())
    }

    tibble(
      site = site_value,
      support_dataset = support_dataset,
      compare_block = compare_block,
      grid_value = as.numeric(reference_df[[key_col]][mismatch_idx]),
      compare_column = col_name,
      reference_value = ref_chr[mismatch_idx],
      merged_value = cmp_chr[mismatch_idx],
      issue_type = 'label_mismatch'
    )
  }))
}

audit_site_support_label_alignment <- function(numbers_at_risk_data, horizon_summary, site_support_dataset_lookup) {
  if (nrow(site_support_dataset_lookup) == 0L) {
    summary <- tibble(
      validation_group = 'cross_dataset_consistency',
      check_name = 'site_support_label_alignment',
      status = 'pass',
      n_problem_rows = 0L,
      check_note = 'No site-support lookup rows were available, so no merged site support-label alignment check was needed.'
    )
    return(list(summary = summary, detail = tibble()))
  }

  compare_cols <- c(
    'interpretation_tier', 'support_tier_standard', 'support_tier', 'support_display_group',
    'support_display_label', 'horizon_evidence_class', 'claim_restriction_flag'
  )

  detail <- bind_rows(lapply(seq_len(nrow(site_support_dataset_lookup)), function(i) {
    site_value <- as.character(site_support_dataset_lookup$site[i])
    support_dataset <- as.character(site_support_dataset_lookup$support_dataset[i])

    reference_horizon <- horizon_summary %>%
      filter(dataset == support_dataset, subgroup_kind == 'overall', subgroup_value == 'Overall', analysis_structure == 'single_site', summary_view == 'transition_only') %>%
      arrange(horizon_year) %>%
      select(horizon_year, all_of(compare_cols))

    merged_horizon <- horizon_summary %>%
      filter(dataset == 'merged', subgroup_kind == 'site', subgroup_value == site_value, analysis_structure == 'site_adjusted', summary_view == 'transition_only') %>%
      arrange(horizon_year) %>%
      select(horizon_year, all_of(compare_cols))

    reference_risk <- numbers_at_risk_data %>%
      filter(dataset == support_dataset, subgroup_kind == 'overall', subgroup_value == 'Overall') %>%
      arrange(time_year) %>%
      select(time_year, all_of(compare_cols))

    merged_risk <- numbers_at_risk_data %>%
      filter(dataset == 'merged', subgroup_kind == 'site', subgroup_value == site_value) %>%
      arrange(time_year) %>%
      select(time_year, all_of(compare_cols))

    bind_rows(
      compare_support_label_frames(reference_horizon, merged_horizon, 'horizon_year', compare_cols, site_value, support_dataset, 'horizon_summary'),
      compare_support_label_frames(reference_risk, merged_risk, 'time_year', compare_cols, site_value, support_dataset, 'numbers_at_risk')
    )
  })) %>% mutate(check_name = 'site_support_label_alignment')

  summary <- tibble(
    validation_group = 'cross_dataset_consistency',
    check_name = 'site_support_label_alignment',
    status = if (nrow(detail) == 0L) 'pass' else 'fail',
    n_problem_rows = nrow(detail),
    check_note = 'Merged site panels must inherit support-tier labels from the corresponding single-site cohort in both horizon summaries and numbers-at-risk tables.'
  )

  list(summary = summary, detail = detail)
}

audit_reporting_preference_rule <- function(horizon_summary) {
  detail <- horizon_summary %>%
    distinct(dataset, subgroup_kind, subgroup_value, group_id, panel_label, analysis_structure, preferred_for_reporting, reporting_priority, reporting_preference_note) %>%
    mutate(
      expected_preferred_for_reporting = dplyr::case_when(
        dataset == 'merged' & subgroup_kind %in% c('overall', 'sex') ~ analysis_structure == 'site_adjusted',
        TRUE ~ TRUE
      ),
      expected_reporting_priority = dplyr::case_when(expected_preferred_for_reporting ~ 1L, TRUE ~ 2L)
    ) %>%
    filter(preferred_for_reporting != expected_preferred_for_reporting | reporting_priority != expected_reporting_priority) %>%
    mutate(check_name = 'reporting_preference_rule')

  summary <- tibble(
    validation_group = 'internal_consistency',
    check_name = 'reporting_preference_rule',
    status = if (nrow(detail) == 0L) 'pass' else 'fail',
    n_problem_rows = nrow(detail),
    check_note = 'Preferred-reporting fields must mark merged overall/sex site_adjusted completeness as primary and site_free as sensitivity only.'
  )

  list(summary = summary, detail = detail)
}

audit_group_size_alignment_to_raw <- function(analysis_datasets_stage2, group_registry, numbers_at_risk_data, reverse_km_summary, horizon_summary) {
  detail <- bind_rows(lapply(seq_len(nrow(group_registry)), function(i) {
    group_row <- group_registry[i, ]
    dataset_name <- as.character(group_row$dataset)
    subgroup_kind_value <- as.character(group_row$subgroup_kind)
    subgroup_value_value <- as.character(group_row$subgroup_value)
    group_id_value <- as.character(group_row$group_id)
    panel_label_value <- as.character(group_row$panel_label)

    raw_df <- filter_group_data(analysis_datasets_stage2[[dataset_name]], subgroup_kind_value, subgroup_value_value)
    raw_n_total <- nrow(raw_df)

    output_checks <- tibble(
      compare_block = c('numbers_at_risk_n_total', 'numbers_at_risk_time0', 'reverse_km_n_total', 'reverse_km_records', 'horizon_summary_n_total'),
      expected_value = c(raw_n_total, raw_n_total, raw_n_total, raw_n_total, raw_n_total),
      observed_value = c(
        pull_first_or_na(numbers_at_risk_data %>% filter(.data$group_id == group_id_value), 'n_total'),
        pull_first_or_na(numbers_at_risk_data %>% filter(.data$group_id == group_id_value, time_year == 0), 'n_risk'),
        pull_first_or_na(reverse_km_summary %>% filter(.data$group_id == group_id_value), 'n_total'),
        pull_first_or_na(reverse_km_summary %>% filter(.data$group_id == group_id_value), 'reverse_km_records'),
        pull_first_or_na(horizon_summary %>% filter(.data$group_id == group_id_value), 'n_total')
      )
    ) %>%
      mutate(
        dataset = dataset_name,
        subgroup_kind = subgroup_kind_value,
        subgroup_value = subgroup_value_value,
        group_id = group_id_value,
        panel_label = panel_label_value,
        issue_type = case_when(
          is.na(observed_value) ~ 'missing_output_value',
          abs(as.numeric(observed_value) - as.numeric(expected_value)) > 1e-8 ~ 'value_mismatch',
          TRUE ~ 'match'
        )
      ) %>%
      filter(issue_type != 'match')

    output_checks
  })) %>%
    mutate(check_name = 'group_size_alignment_to_raw')

  summary <- tibble(
    validation_group = 'raw_reconstruction',
    check_name = 'group_size_alignment_to_raw',
    status = if (nrow(detail) == 0L) 'pass' else 'fail',
    n_problem_rows = nrow(detail),
    check_note = 'Raw group sizes must match Stage 2 n_total fields, time-0 numbers at risk, and reverse-KM record counts.'
  )

  list(summary = summary, detail = detail)
}

audit_site_admin_lookup_reconstruction <- function(analysis_datasets_stage2, site_admin_lookup, merged_admin_end_date) {
  expected_lookup <- analysis_datasets_stage2[['merged']] %>%
    group_by(site) %>%
    summarise(
      expected_site_admin_end_date = max(followup_date, na.rm = TRUE),
      expected_site_admin_end_days_followup_max = max(days_followup, na.rm = TRUE),
      expected_site_admin_end_years_followup_max = max(time_year, na.rm = TRUE),
      .groups = 'drop'
    )

  detail_lookup <- expected_lookup %>%
    full_join(site_admin_lookup, by = 'site') %>%
    mutate(
      date_diff_days = abs(as.numeric(expected_site_admin_end_date - site_admin_end_date)),
      days_followup_diff = abs(as.numeric(expected_site_admin_end_days_followup_max) - as.numeric(site_admin_end_days_followup_max)),
      years_followup_diff = abs(as.numeric(expected_site_admin_end_years_followup_max) - as.numeric(site_admin_end_years_followup_max)),
      issue_type = case_when(
        is.na(site_admin_end_date) | is.na(expected_site_admin_end_date) ~ 'missing_site_row',
        date_diff_days > 0 | days_followup_diff > 1e-8 | years_followup_diff > 1e-8 ~ 'value_mismatch',
        TRUE ~ 'match'
      )
    ) %>%
    filter(issue_type != 'match') %>%
    mutate(check_name = 'site_admin_lookup_reconstruction')

  expected_global_admin_end_date <- max(analysis_datasets_stage2[['merged']]$followup_date, na.rm = TRUE)
  detail_global <- tibble(
    site = NA_character_,
    expected_site_admin_end_date = as.Date(NA),
    expected_site_admin_end_days_followup_max = NA_real_,
    expected_site_admin_end_years_followup_max = NA_real_,
    site_admin_end_date = as.Date(NA),
    site_admin_end_days_followup_max = NA_real_,
    site_admin_end_years_followup_max = NA_real_,
    date_diff_days = abs(as.numeric(expected_global_admin_end_date - merged_admin_end_date)),
    days_followup_diff = NA_real_,
    years_followup_diff = NA_real_,
    issue_type = ifelse(abs(as.numeric(expected_global_admin_end_date - merged_admin_end_date)) > 0, 'global_admin_mismatch', 'match'),
    check_name = 'site_admin_lookup_reconstruction'
  ) %>% filter(issue_type != 'match')

  detail <- bind_rows(detail_lookup, detail_global)

  summary <- tibble(
    validation_group = 'raw_reconstruction',
    check_name = 'site_admin_lookup_reconstruction',
    status = if (nrow(detail) == 0L) 'pass' else 'fail',
    n_problem_rows = nrow(detail),
    check_note = 'Site-level and merged administrative end dates must be reproducible directly from the merged Stage 2 raw data.'
  )

  list(summary = summary, detail = detail)
}

run_stage2_validation_audits <- function(
    analysis_datasets_stage2,
    group_registry,
    site_admin_lookup,
    merged_admin_end_date,
    numbers_at_risk_data,
    composition_data,
    horizon_summary,
    reverse_km_summary,
    site_contribution_summary,
    site_support_dataset_lookup
) {
  audit_results <- list(
    audit_numbers_at_risk_monotonic(numbers_at_risk_data),
    audit_composition_totals(composition_data),
    audit_horizon_monotonicity(horizon_summary),
    audit_required_stage2_views(horizon_summary),
    audit_merged_site_consistency(horizon_summary, reverse_km_summary),
    audit_site_contribution_alignment(horizon_summary, site_contribution_summary),
    audit_site_support_label_alignment(numbers_at_risk_data, horizon_summary, site_support_dataset_lookup),
    audit_reporting_preference_rule(horizon_summary),
    audit_group_size_alignment_to_raw(analysis_datasets_stage2, group_registry, numbers_at_risk_data, reverse_km_summary, horizon_summary),
    audit_site_admin_lookup_reconstruction(analysis_datasets_stage2, site_admin_lookup, merged_admin_end_date)
  )

  validation_summary <- bind_rows(lapply(audit_results, `[[`, 'summary')) %>%
    arrange(validation_group, check_name)

  validation_issue_rows <- bind_rows(lapply(audit_results, `[[`, 'detail')) %>%
    arrange(check_name)

  list(
    validation_summary = validation_summary,
    validation_issue_rows = validation_issue_rows
  )
}

make_quality_flags <- function(horizon_summary, reverse_km_summary) {
  masked_flags <- horizon_summary %>%
    filter(completeness_plot_status %in% c('masked_no_eligible', 'masked_no_risk_set')) %>%
    distinct(dataset, subgroup_kind, subgroup_value, group_id, panel_label, analysis_structure, summary_view, horizon_year, completeness_plot_status, completeness_plot_note) %>%
    transmute(
      flag_domain = 'completeness',
      flag_severity = 'warning',
      flag_code = paste0('completeness_', completeness_plot_status),
      dataset, subgroup_kind, subgroup_value, group_id, panel_label,
      analysis_structure, summary_view, horizon_year,
      flag_text_value = completeness_plot_status,
      flag_numeric_value = NA_real_,
      flag_message = completeness_plot_note
    )

  sparse_flags <- horizon_summary %>%
    filter(completeness_plot_status %in% c('very_sparse', 'sparse', 'low_fraction')) %>%
    distinct(dataset, subgroup_kind, subgroup_value, group_id, panel_label, analysis_structure, summary_view, horizon_year, completeness_plot_status, completeness_plot_note, risk_set_n, risk_set_fraction) %>%
    transmute(
      flag_domain = 'completeness',
      flag_severity = case_when(completeness_plot_status == 'very_sparse' ~ 'warning', TRUE ~ 'caution'),
      flag_code = paste0('completeness_', completeness_plot_status),
      dataset, subgroup_kind, subgroup_value, group_id, panel_label, analysis_structure, summary_view, horizon_year,
      flag_text_value = completeness_plot_status,
      flag_numeric_value = case_when(
        completeness_plot_status %in% c('very_sparse', 'sparse') ~ as.numeric(risk_set_n),
        completeness_plot_status == 'low_fraction' ~ as.numeric(risk_set_fraction),
        TRUE ~ NA_real_
      ),
      flag_message = completeness_plot_note
    )

  dominance_flags <- horizon_summary %>%
    filter(site_dominance_warning_flag) %>%
    distinct(dataset, subgroup_kind, subgroup_value, group_id, panel_label, analysis_structure, horizon_year, dominant_site, dominant_site_fraction, site_dominance_status, site_dominance_warning_note) %>%
    transmute(
      flag_domain = 'merged_site_contribution',
      flag_severity = case_when(site_dominance_status == 'single_site_only' ~ 'warning', TRUE ~ 'caution'),
      flag_code = paste0('site_dominance_', site_dominance_status),
      dataset, subgroup_kind, subgroup_value, group_id, panel_label, analysis_structure,
      summary_view = NA_character_,
      horizon_year,
      flag_text_value = dominant_site,
      flag_numeric_value = dominant_site_fraction,
      flag_message = site_dominance_warning_note
    )

  reverse_flags <- reverse_km_summary %>%
    filter(!is.na(reverse_km_median_followup_year) & (is.na(reverse_km_median_followup_lcl_year) | is.na(reverse_km_median_followup_ucl_year))) %>%
    transmute(
      flag_domain = 'reverse_km',
      flag_severity = 'caution',
      flag_code = 'reverse_km_median_ci_incomplete',
      dataset, subgroup_kind, subgroup_value, group_id, panel_label,
      analysis_structure = NA_character_,
      summary_view,
      horizon_year = NA_integer_,
      flag_text_value = NA_character_,
      flag_numeric_value = reverse_km_median_followup_year,
      flag_message = 'Reverse-KM median follow-up was estimable but at least one confidence bound was missing; this usually reflects sparse subgroup tail support.'
    )

  bind_rows(empty_quality_flags_tbl(), masked_flags, sparse_flags, dominance_flags, reverse_flags) %>%
    distinct() %>%
    arrange(flag_domain, dataset, subgroup_kind, subgroup_value, analysis_structure, summary_view, horizon_year, flag_code)
}


# 🔴 Define: plot builders and exporters ===============================
theme_stage2 <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = 'bold'),
      legend.position = 'bottom',
      legend.title = element_text(face = 'bold')
    )
}

build_stage2_plot_objects <- function(
    curve_plot_data,
    numbers_at_risk_data,
    completeness_plot_data,
    panel_order,
    subgroup_kind_filter,
    panel_family_label
) {
  plot_group_filter <- curve_plot_data %>%
    distinct(dataset, subgroup_kind, subgroup_value, panel_label) %>%
    filter(subgroup_kind %in% subgroup_kind_filter)

  if (nrow(plot_group_filter) == 0L) {
    warning(sprintf('No plot-eligible groups were found for panel family `%s`.', panel_family_label), call. = FALSE)
    return(NULL)
  }

  curve_plot_subset <- curve_plot_data %>%
    semi_join(plot_group_filter, by = c('dataset', 'subgroup_kind', 'subgroup_value', 'panel_label')) %>%
    mutate(panel_label = factor_panel_labels(panel_label, panel_order))

  risk_plot_subset <- numbers_at_risk_data %>%
    semi_join(plot_group_filter, by = c('dataset', 'subgroup_kind', 'subgroup_value', 'panel_label')) %>%
    mutate(
      panel_label = factor_panel_labels(panel_label, panel_order),
      time_year_factor = factor(time_year, levels = risk_table_times_year)
    )

  completeness_plot_subset <- completeness_plot_data %>%
    semi_join(plot_group_filter, by = c('dataset', 'subgroup_kind', 'subgroup_value', 'panel_label')) %>%
    mutate(panel_label = factor_panel_labels(panel_label, panel_order))

  composition_plot_data <- curve_plot_subset %>%
    filter(curve_role == 'composition_state') %>%
    filter(!is.na(state), is.finite(time_year), is.finite(state_prop)) %>%
    arrange(summary_view, panel_label, state, time_year)

  km_curve_plot <- curve_plot_subset %>%
    filter(curve_role == 'km_curve') %>%
    ggplot(aes(x = time_year, y = estimate, color = summary_view, fill = summary_view)) +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.15, color = NA, na.rm = TRUE) +
    geom_step(linewidth = 0.7, na.rm = TRUE) +
    facet_wrap(~ panel_label, ncol = 2) +
    scale_x_continuous(breaks = 0:10) +
    scale_y_continuous(labels = label_percent(accuracy = 1)) +
    coord_cartesian(xlim = c(0, plot_max_year), ylim = c(0, 1)) +
    labs(
      title = sprintf('Stage 2 Kaplan-Meier curves (%s)', panel_family_label),
      subtitle = 'transition_only vs remission_sensitive',
      x = 'Years from cohort entry',
      y = 'Survival probability',
      color = 'Summary view',
      fill = 'Summary view'
    ) +
    theme_stage2()

  reverse_km_plot <- curve_plot_subset %>%
    filter(curve_role == 'reverse_km_curve') %>%
    ggplot(aes(x = time_year, y = estimate, color = summary_view, fill = summary_view)) +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), alpha = 0.15, color = NA, na.rm = TRUE) +
    geom_step(linewidth = 0.7, na.rm = TRUE) +
    facet_wrap(~ panel_label, ncol = 2) +
    scale_x_continuous(breaks = 0:10) +
    scale_y_continuous(labels = label_percent(accuracy = 1)) +
    coord_cartesian(xlim = c(0, plot_max_year), ylim = c(0, 1)) +
    labs(
      title = sprintf('Stage 2 reverse Kaplan-Meier curves (%s)', panel_family_label),
      subtitle = 'Descriptive follow-up maturity curves',
      x = 'Years from cohort entry',
      y = 'Reverse-KM survival probability',
      color = 'Summary view',
      fill = 'Summary view'
    ) +
    theme_stage2()

  n_risk_plot <- risk_plot_subset %>%
    ggplot(aes(x = time_year_factor, y = 'At risk')) +
    geom_tile(
      aes(fill = support_display_group),
      width = 0.95,
      height = 0.75,
      color = 'white',
      linewidth = 0.35
    ) +
    geom_text(
      aes(label = risk_table_label),
      size = 3.4,
      fontface = 'bold',
      na.rm = TRUE
    ) +
    facet_wrap(~ panel_label, ncol = 2) +
    scale_fill_manual(
      values = c(
        baseline = '#d9d9d9',
        primary_supported = '#cfe8cf',
        secondary_or_sensitivity = '#f5deb3',
        projection = '#e4d9ff',
        other = '#f0f0f0'
      ),
      breaks = c('baseline', 'primary_supported', 'secondary_or_sensitivity', 'projection'),
      labels = c('Baseline', 'Primary supported', 'Secondary / sensitivity', 'Projection')
    ) +
    labs(
      title = sprintf('Stage 2 numbers at risk (%s)', panel_family_label),
      subtitle = 'Table-like display at 0, 1, ..., 10 years',
      x = 'Years from cohort entry',
      y = NULL,
      fill = 'Support display',
      caption = 'Cell fill shows simplified horizon-support group. Secondary and sensitivity horizons are collapsed for display only.'
    ) +
    theme_stage2() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(face = 'bold'),
      axis.ticks = element_blank()
    )

  composition_plot <- composition_plot_data %>%
    ggplot(aes(x = time_year, y = state_prop, fill = state, group = state)) +
    geom_area(
      stat = 'identity',
      position = 'stack',
      alpha = 0.9,
      color = 'white',
      linewidth = 0.15,
      na.rm = TRUE
    ) +
    facet_grid(summary_view ~ panel_label) +
    scale_x_continuous(breaks = 0:10) +
    scale_y_continuous(labels = label_percent(accuracy = 1)) +
    coord_cartesian(xlim = c(0, plot_max_year), ylim = c(0, 1)) +
    labs(
      title = sprintf('Stage 2 event / censoring / remission composition (%s)', panel_family_label),
      subtitle = 'Stacked composition by summary view',
      x = 'Years from cohort entry',
      y = 'Composition proportion',
      fill = 'State'
    ) +
    theme_stage2()

  completeness_plot_data_stable <- completeness_plot_subset %>%
    filter(completeness_plot_status == 'stable')

  completeness_plot_data_sparse <- completeness_plot_subset %>%
    filter(completeness_plot_status %in% c('low_fraction', 'sparse'))

  completeness_plot_data_very_sparse <- completeness_plot_subset %>%
    filter(completeness_plot_status == 'very_sparse')

  completeness_plot <- completeness_plot_subset %>%
    ggplot(aes(
      x = horizon_year,
      y = completeness_value,
      color = completeness_metric,
      linetype = analysis_structure,
      group = interaction(completeness_metric, analysis_structure)
    )) +
    geom_line(linewidth = 0.7, na.rm = TRUE) +
    geom_point(data = completeness_plot_data_stable, size = 1.7, na.rm = TRUE) +
    geom_point(data = completeness_plot_data_sparse, size = 2.0, shape = 1, stroke = 0.8, na.rm = TRUE) +
    geom_point(data = completeness_plot_data_very_sparse, size = 2.1, shape = 4, stroke = 0.9, na.rm = TRUE) +
    facet_grid(summary_view ~ panel_label) +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(labels = label_percent(accuracy = 1)) +
    coord_cartesian(xlim = c(1, 10), ylim = c(0, 1)) +
    labs(
      title = sprintf('Stage 2 follow-up completeness summaries (%s)', panel_family_label),
      subtitle = 'percentage_method, CCI, and SPT on the common 1-10 year grid',
      x = 'Horizon (years)',
      y = 'Completeness summary',
      color = 'Completeness metric',
      linetype = 'Analysis structure',
      caption = paste(
        'Source-of-truth CSV: stage2_followup_completeness_plot_data.csv.',
        'For merged overall/sex panels, site_adjusted is the preferred reporting view and site_free is retained as sensitivity.',
        'Horizons with zero eligible subjects or zero risk set are masked. Open circles and × markers indicate sparse late-horizon support.'
      )
    ) +
    theme_stage2()

  list(
    km_curve = km_curve_plot,
    reverse_km = reverse_km_plot,
    numbers_at_risk = n_risk_plot,
    composition = composition_plot,
    completeness = completeness_plot
  )
}

save_stage2_plot_collection <- function(plot_objects, output_pdf_path, png_prefix) {
  if (is.null(plot_objects) || length(plot_objects) == 0L) {
    return(character())
  }

  grDevices::pdf(output_pdf_path, width = 14, height = 10, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (plot_name in names(plot_objects)) {
    print(plot_objects[[plot_name]])
  }

  png_paths <- character()

  for (plot_name in names(plot_objects)) {
    png_file <- file.path(export_path, sprintf('%s_%s.png', png_prefix, plot_name))
    ggplot2::ggsave(
      filename = png_file,
      plot = plot_objects[[plot_name]],
      width = 14,
      height = 10,
      units = 'in',
      dpi = 300
    )
    png_paths <- c(png_paths, png_file)
  }

  png_paths
}

# 🔴 Define: metadata, cache, and manifest helpers ===============================
make_stage2_metadata <- function(
    stage1_inputs,
    site_admin_lookup,
    plot_time_grid_year,
    threshold_vector,
    merged_admin_end_date,
    main_plot_panel_order,
    sex_plot_panel_order,
    data_path,
    export_path,
    reused_existing_core_results,
    core_reuse_mode,
    raw_merged_file,
    stage1_source_signature
) {
  bundle_cfg <- stage1_inputs$backbone_bundle$config %||% list()

  tibble::tribble(
    ~metadata_group, ~metadata_name, ~metadata_value,
    'stage', 'stage_name', 'Stage 2 follow-up maturity',
    'stage', 'stage_role', 'Describe follow-up maturity before any cure claim',
    'stage', 'model_fitting_allowed', 'FALSE',
    'stage', 'patch_level', current_patch_level,
    'stage', 'scientific_compatibility_signature', scientific_compatibility_signature,
    'stage', 'compatible_reuse_patch_levels', paste(compatible_reuse_patch_levels, collapse = '|'),
    'stage', 'compatible_science_signatures', paste(compatible_science_signatures, collapse = '|'),
    'stage', 'reused_existing_core_results', if (isTRUE(reused_existing_core_results)) 'TRUE' else 'FALSE',
    'stage', 'core_reuse_mode', as.character(core_reuse_mode),
    'inputs', 'stage1_source_signature', as.character(stage1_source_signature),
    'inputs', 'data_path_resolved', data_path,
    'inputs', 'export_path_resolved', export_path,
    'inputs', 'stage1_analysis_datasets_file', stage1_analysis_datasets_file,
    'inputs', 'stage1_backbone_bundle_file', stage1_backbone_bundle_file,
    'inputs', 'stage1_scaling_registry_file', stage1_scaling_registry_file,
    'inputs', 'stage1_horizon_registry_file', stage1_horizon_registry_file,
    'inputs', 'stage1_threshold_registry_file', stage1_threshold_registry_file,
    'inputs', 'stage1_metadata_registry_file', stage1_metadata_registry_file,
    'inputs', 'stage1_dataset_registry_file', stage1_dataset_registry_file,
    'inputs', 'stage1_formula_registry_file', stage1_formula_registry_file,
    'inputs', 'raw_merged_file_for_date_recovery', coalesce_character(raw_merged_file, bundle_cfg$merged_file %||% NA_character_),
    'time', 'analysis_time_variable', 'days_followup / time_year',
    'time', 'horizon_vector_year', paste(sort(unique(stage1_inputs$horizon_registry$horizon_year)), collapse = ','),
    'time', 'plot_grid_year', paste(plot_time_grid_year, collapse = ','),
    'time', 'threshold_vector_from_stage1', paste(format(threshold_vector, trim = TRUE, scientific = FALSE), collapse = ','),
    'summaries', 'main_transition_only_rule', 'transition_only = status_num == 1 as event; status_num %in% c(0,2) as censoring',
    'summaries', 'remission_sensitive_rule', 'remission_sensitive = status_num %in% c(1,2) as informative terminal outcomes; status_num == 0 as censoring',
    'summaries', 'reverse_km_note', 'Reverse KM is descriptive of follow-up maturity and should not be over-interpreted as a direct follow-up-rate estimator.',
    'summaries', 'followup_completeness_note', 'percentage_method, CCI, and SPT are retained as raw maturity indicators; plotting masks horizons with zero eligible subjects or zero risk set.',
    'summaries', 'analysis_structure_note', 'Merged overall and merged sex summaries retain both site_free and site_adjusted completeness structures; merged site-stratified summaries use site_adjusted completeness.',
    'summaries', 'preferred_reporting_rule', 'For merged overall and merged sex completeness outputs, site_adjusted is the preferred reporting view; site_free is retained only as sensitivity.',
    'summaries', 'horizon_summary_scope_note', 'stage2_followup_horizon_summary.csv contains horizon-level metrics plus merged site-contribution descriptors, late-horizon site-dominance warning fields, plot-masking flags, preferred_for_reporting fields, and horizon support/evidence/claim-restriction fields.',
    'summaries', 'completeness_plot_source_note', 'stage2_followup_completeness_plot_data.csv is the explicit source-of-truth for the completeness plot and is derived directly from stage2_followup_horizon_summary.csv.',
    'summaries', 'png_export_note', 'Each Stage 2 plot is exported both within bundled PDF files and as independent PNG files generated from the same plot-source data frames used for the CSV exports.',
    'merged', 'site_specific_admin_end_date_sites', collapse_unique_chr(site_admin_lookup$site),
    'merged', 'merged_global_admin_end_date', as.character(merged_admin_end_date),
    'merged', 'merged_completeness_rule', 'Use site-specific administrative end date for merged site_adjusted completeness; use pooled merged administrative end date for merged site_free completeness.',
    'merged', 'site_dominance_rule', sprintf('Describe site composition for all merged overall/sex horizons, but escalate warnings only when one site fully remains eligible or when one site exceeds %.0f%% of eligible subjects at horizons >= %d years outside primary-supported support tiers.', 100 * site_dominance_warning_fraction, site_dominance_warning_min_horizon_year),
    'merged', 'site_dominance_primary_supported_exemption', if (isTRUE(site_dominance_warning_exempt_primary_supported)) 'TRUE' else 'FALSE',
    'plots', 'numbers_at_risk_plot_note', 'Numbers-at-risk page uses table-like tiles; fill encodes simplified horizon-support groups and avoids the prior floating-text layout.',
    'plots', 'completeness_plot_masking_rule', sprintf('Mask completeness plot points when eligible_n <= 0 or risk_set_n <= 0; keep sparse points visible with caution thresholds risk_set_n < %d and risk_set_fraction < %.2f.', completeness_sparse_risk_n, completeness_low_fraction_cutoff),
    'plots', 'completeness_plot_reporting_note', 'For merged overall/sex panels, site_adjusted is the preferred reporting view and site_free is retained as sensitivity.',
    'plots', 'support_display_group_note', 'support_display_group collapses secondary and sensitivity horizons for display only; merged site panels inherit support tiers from the matched single-site cohort when available.',
    'plots', 'main_plot_panels', paste(main_plot_panel_order, collapse = ' | '),
    'plots', 'sex_plot_panels', if (length(sex_plot_panel_order) > 0L) paste(sex_plot_panel_order, collapse = ' | ') else 'not_exported',
    'validation', 'validation_log_note', 'Write stage2_validation_summary.csv and stage2_validation_issue_rows.csv before any failure stop so audit traces remain inspectable.',
    'validation', 'quality_flag_note', 'stage2_quality_flags.csv centralizes non-fatal sparse-tail, masking, reverse-KM CI, and late-horizon merged site-dominance cautions after the primary-supported exemption gate is applied.',
    'validation', 'raw_reconstruction_note', 'Validation includes direct raw-data reconstructions of group sizes, numbers at risk, merged site administrative end dates, and preferred-reporting flags.',
    'validation', 'sex_plot_exported', if (isTRUE(export_sex_panel_plots)) 'TRUE' else 'FALSE'
  )
}

detect_reusable_stage2_state <- function(export_path, expected_stage1_source_signature) {
  metadata_path <- file.path(export_path, 'stage2_metadata_registry.csv')
  bundle_path <- file.path(export_path, 'stage2_followup_bundle.rds')

  if (isTRUE(reuse_existing_core_results) && !isTRUE(force_recompute_core_results) && file.exists(bundle_path)) {
    bundle <- tryCatch(readRDS(bundle_path), error = function(e) NULL)
    bundle_cfg <- bundle$config %||% list()
    bundle_stage1_signature <- as.character(bundle_cfg$stage1_source_signature %||% NA_character_)
    bundle_patch <- as.character(bundle_cfg$patch_level %||% NA_character_)
    bundle_science <- as.character(bundle_cfg$scientific_compatibility_signature %||% NA_character_)

    if (is.list(bundle) &&
        all(c('site_admin_lookup', 'reverse_km_summary', 'numbers_at_risk_data', 'composition_data',
              'site_contribution_data', 'horizon_summary', 'curve_plot_data') %in% names(bundle)) &&
        !is.na(bundle_stage1_signature) && identical(bundle_stage1_signature, expected_stage1_source_signature) &&
        bundle_patch %in% compatible_reuse_patch_levels &&
        bundle_science %in% compatible_science_signatures) {
      return(list(mode = 'bundle', bundle = bundle))
    }
  }

  minimal_files <- c(
    stage2_core_paths$site_admin_lookup,
    stage2_core_paths$reverse_km_summary,
    stage2_core_paths$numbers_at_risk,
    stage2_core_paths$composition,
    stage2_core_paths$site_contribution,
    stage2_core_paths$horizon_summary,
    stage2_core_paths$curve_plot_data,
    metadata_path
  )

  if (isTRUE(reuse_existing_core_results) && !isTRUE(force_recompute_core_results) && all(file.exists(minimal_files))) {
    metadata_tbl <- tryCatch(read_csv_required(metadata_path), error = function(e) NULL)

    if (!is.null(metadata_tbl) && nrow(metadata_tbl) > 0L) {
      stage1_signature_csv <- metadata_tbl %>%
        filter(metadata_name == 'stage1_source_signature') %>%
        pull(metadata_value) %>%
        unique() %>%
        as.character()
      patch_csv <- metadata_tbl %>%
        filter(metadata_name == 'patch_level') %>%
        pull(metadata_value) %>%
        unique() %>%
        as.character()
      science_csv <- metadata_tbl %>%
        filter(metadata_name == 'scientific_compatibility_signature') %>%
        pull(metadata_value) %>%
        unique() %>%
        as.character()

      if (length(stage1_signature_csv) == 1L && !is.na(stage1_signature_csv) && identical(stage1_signature_csv, expected_stage1_source_signature) &&
          length(patch_csv) == 1L && patch_csv %in% compatible_reuse_patch_levels &&
          length(science_csv) == 1L && science_csv %in% compatible_science_signatures) {
        return(list(mode = 'csv', bundle = NULL))
      }
    }
  }

  list(mode = 'none', bundle = NULL)
}

augment_horizon_support_fields <- function(df, horizon_registry, key_time_col) {
  if (nrow(df) == 0L) {
    return(df)
  }

  if (!'dataset' %in% names(df) || !key_time_col %in% names(df)) {
    return(df)
  }

  join_df <- horizon_registry %>%
    select(dataset, horizon_year, interpretation_tier, interpretation_note, support_tier_standard, support_tier, support_display_group, support_display_label, horizon_evidence_class, claim_restriction_flag)

  time_values <- as.numeric(df[[key_time_col]])
  out <- df

  if (all(time_values %in% c(0, 1:10))) {
    out <- out %>%
      mutate(.join_horizon_year = if_else(as.numeric(.data[[key_time_col]]) == 0, NA_integer_, as.integer(as.numeric(.data[[key_time_col]])))) %>%
      left_join(join_df, by = c('dataset', '.join_horizon_year' = 'horizon_year')) %>%
      mutate(
        interpretation_tier = if_else(as.numeric(.data[[key_time_col]]) == 0, 'baseline', interpretation_tier),
        support_tier_standard = if_else(as.numeric(.data[[key_time_col]]) == 0, 'baseline', support_tier_standard),
        support_tier = if_else(as.numeric(.data[[key_time_col]]) == 0, 'baseline', support_tier),
        support_display_group = if_else(as.numeric(.data[[key_time_col]]) == 0, 'baseline', support_display_group),
        support_display_label = if_else(as.numeric(.data[[key_time_col]]) == 0, 'Baseline', support_display_label),
        horizon_evidence_class = if_else(as.numeric(.data[[key_time_col]]) == 0, 'baseline', horizon_evidence_class),
        claim_restriction_flag = if_else(as.numeric(.data[[key_time_col]]) == 0, 'primary_claim_allowed', claim_restriction_flag)
      ) %>%
      select(-.join_horizon_year)
  }

  out
}

upgrade_stage2_outputs_for_current_spec <- function(outputs, horizon_registry) {
  outputs$horizon_summary <- outputs$horizon_summary %>%
    left_join(
      horizon_registry %>%
        select(dataset, horizon_year, interpretation_tier, interpretation_note, support_tier_standard, support_tier, support_display_group, support_display_label, horizon_evidence_class, claim_restriction_flag),
      by = c('dataset', 'horizon_year'),
      suffix = c('', '.hr')
    ) %>%
    mutate(
      interpretation_tier = coalesce(interpretation_tier, interpretation_tier.hr),
      interpretation_note = coalesce(interpretation_note, interpretation_note.hr),
      support_tier_standard = coalesce(support_tier_standard, support_tier_standard.hr),
      support_tier = coalesce(support_tier, support_tier.hr, support_tier_standard),
      support_display_group = coalesce(support_display_group, support_display_group.hr),
      support_display_label = coalesce(support_display_label, support_display_label.hr),
      horizon_evidence_class = coalesce(horizon_evidence_class, horizon_evidence_class.hr),
      claim_restriction_flag = coalesce(claim_restriction_flag, claim_restriction_flag.hr)
    ) %>%
    select(-ends_with('.hr'))

  outputs$numbers_at_risk_data <- augment_horizon_support_fields(outputs$numbers_at_risk_data, horizon_registry, 'time_year')

  outputs$site_contribution_data <- outputs$site_contribution_data %>%
    left_join(
      horizon_registry %>%
        select(dataset, horizon_year, interpretation_tier, interpretation_note, support_tier_standard, support_tier, support_display_group, support_display_label, horizon_evidence_class, claim_restriction_flag),
      by = c('dataset', 'horizon_year'),
      suffix = c('', '.hr')
    ) %>%
    mutate(
      interpretation_tier = coalesce(interpretation_tier, interpretation_tier.hr),
      interpretation_note = coalesce(interpretation_note, interpretation_note.hr),
      support_tier_standard = coalesce(support_tier_standard, support_tier_standard.hr),
      support_tier = coalesce(support_tier, support_tier.hr, support_tier_standard),
      support_display_group = coalesce(support_display_group, support_display_group.hr),
      support_display_label = coalesce(support_display_label, support_display_label.hr),
      horizon_evidence_class = coalesce(horizon_evidence_class, horizon_evidence_class.hr),
      claim_restriction_flag = coalesce(claim_restriction_flag, claim_restriction_flag.hr)
    ) %>%
    select(-ends_with('.hr'))

  outputs
}

confirm_stage2_exports <- function(export_manifest) {
  expected_paths <- export_manifest$file_path
  missing_paths <- expected_paths[!file.exists(expected_paths)]

  if (length(missing_paths) > 0L) {
    stop(sprintf('Stage 2 completed, but some expected export files are missing: %s', paste(missing_paths, collapse = ', ')), call. = FALSE)
  }

  invisible(expected_paths)
}

# 🔴 Load: stage-one reusable backbone ===============================
stage1_inputs <- read_stage1_inputs()
validate_stage1_inputs(stage1_inputs)

site_labels <- infer_site_labels_from_stage1(stage1_inputs)
pnu_site_label <- site_labels[['pnu_site_label']]
snu_site_label <- site_labels[['snu_site_label']]

analysis_datasets_stage1 <- stage1_inputs$analysis_datasets[dataset_order]
horizon_registry <- normalize_horizon_registry(stage1_inputs$horizon_registry)
threshold_vector <- extract_threshold_vector(stage1_inputs$threshold_registry)
scaling_registry_stage1 <- stage1_inputs$scaling_registry
dataset_registry <- stage1_inputs$dataset_registry

need_date_recovery <- !all(vapply(analysis_datasets_stage1, function(df) 'date_entry' %in% names(df), logical(1)))

if (isTRUE(need_date_recovery)) {
  date_lookup_info <- build_date_entry_lookup_from_raw(stage1_inputs, pnu_label = pnu_site_label, snu_label = snu_site_label)
  date_lookup <- date_lookup_info$date_lookup
  raw_merged_file_used <- date_lookup_info$raw_merged_file
} else {
  date_lookup <- tibble(unique_person_id = character(), date_entry = as.Date(character()))
  raw_merged_file_used <- NA_character_
}

analysis_datasets_with_dates <- list(
  PNU = ensure_date_entry_in_dataset(analysis_datasets_stage1[['PNU']], 'PNU', date_lookup),
  SNU = ensure_date_entry_in_dataset(analysis_datasets_stage1[['SNU']], 'SNU', date_lookup),
  merged = ensure_date_entry_in_dataset(analysis_datasets_stage1[['merged']], 'merged', date_lookup)
)

analysis_datasets_prepped <- list(
  PNU = prepare_stage2_dataset(analysis_datasets_with_dates[['PNU']], 'PNU', site_mode = 'single'),
  SNU = prepare_stage2_dataset(analysis_datasets_with_dates[['SNU']], 'SNU', site_mode = 'single'),
  merged = prepare_stage2_dataset(analysis_datasets_with_dates[['merged']], 'merged', site_mode = 'merged')
)

analysis_datasets_prepped_df <- purrr::map(analysis_datasets_prepped, 'data')

stage2_scaling_registry <- bind_rows(
  scaling_registry_stage1,
  purrr::map(analysis_datasets_prepped, 'scaling')
) %>%
  distinct()

site_admin_lookup_current <- compute_site_admin_lookup(analysis_datasets_prepped_df[['merged']])
merged_admin_end_date_current <- max(analysis_datasets_prepped_df[['merged']]$followup_date, na.rm = TRUE)

analysis_datasets_stage2 <- list(
  PNU = attach_admin_dates(analysis_datasets_prepped_df[['PNU']], 'PNU', site_admin_lookup_current, merged_admin_end_date_current),
  SNU = attach_admin_dates(analysis_datasets_prepped_df[['SNU']], 'SNU', site_admin_lookup_current, merged_admin_end_date_current),
  merged = attach_admin_dates(analysis_datasets_prepped_df[['merged']], 'merged', site_admin_lookup_current, merged_admin_end_date_current)
)

site_support_dataset_lookup <- build_site_support_dataset_lookup(analysis_datasets_stage2)
group_registry <- build_group_registry(analysis_datasets_stage2)
main_plot_panel_order <- get_plot_panel_order(group_registry, subgroup_kind_order = c('overall', 'site'))
sex_plot_panel_order <- get_plot_panel_order(group_registry, subgroup_kind_order = c('sex'))

plot_time_grid_year <- sort(unique(c(seq(0, plot_max_year, by = plot_time_step_year), 1:10)))
risk_table_times_year <- sort(unique(as.numeric(risk_table_times_year)))

# 🔴 Reuse: compatible Stage 2 core outputs when possible ===============================
reuse_state <- detect_reusable_stage2_state(export_path, expected_stage1_source_signature = current_stage1_source_signature)
reused_existing_core_results <- FALSE
core_reuse_mode <- 'none'

if (reuse_state$mode == 'bundle') {
  reused_existing_core_results <- TRUE
  core_reuse_mode <- 'bundle'

  existing_bundle <- reuse_state$bundle

  site_admin_lookup <- tibble::as_tibble(existing_bundle$site_admin_lookup)
  reverse_km_summary <- tibble::as_tibble(existing_bundle$reverse_km_summary)
  numbers_at_risk_data <- tibble::as_tibble(existing_bundle$numbers_at_risk_data)
  composition_data <- tibble::as_tibble(existing_bundle$composition_data)
  site_contribution_data <- tibble::as_tibble(existing_bundle$site_contribution_data)
  horizon_summary <- tibble::as_tibble(existing_bundle$horizon_summary)
  curve_plot_data <- tibble::as_tibble(existing_bundle$curve_plot_data)

  outputs_upgrade <- list(
    site_admin_lookup = site_admin_lookup,
    reverse_km_summary = reverse_km_summary,
    numbers_at_risk_data = numbers_at_risk_data,
    composition_data = composition_data,
    site_contribution_data = site_contribution_data,
    horizon_summary = horizon_summary,
    curve_plot_data = curve_plot_data
  )
  outputs_upgrade <- upgrade_stage2_outputs_for_current_spec(outputs_upgrade, horizon_registry)

  site_admin_lookup <- outputs_upgrade$site_admin_lookup
  reverse_km_summary <- outputs_upgrade$reverse_km_summary
  numbers_at_risk_data <- outputs_upgrade$numbers_at_risk_data
  composition_data <- outputs_upgrade$composition_data
  site_contribution_data <- outputs_upgrade$site_contribution_data
  horizon_summary <- outputs_upgrade$horizon_summary
  curve_plot_data <- outputs_upgrade$curve_plot_data
} else if (reuse_state$mode == 'csv') {
  reused_existing_core_results <- TRUE
  core_reuse_mode <- 'csv'

  site_admin_lookup <- read_csv_required(stage2_core_paths$site_admin_lookup)
  reverse_km_summary <- read_csv_required(stage2_core_paths$reverse_km_summary)
  numbers_at_risk_data <- read_csv_required(stage2_core_paths$numbers_at_risk)
  composition_data <- read_csv_required(stage2_core_paths$composition)
  site_contribution_data <- read_csv_required(stage2_core_paths$site_contribution)
  horizon_summary <- read_csv_required(stage2_core_paths$horizon_summary)
  curve_plot_data <- read_csv_required(stage2_core_paths$curve_plot_data)

  outputs_upgrade <- list(
    site_admin_lookup = site_admin_lookup,
    reverse_km_summary = reverse_km_summary,
    numbers_at_risk_data = numbers_at_risk_data,
    composition_data = composition_data,
    site_contribution_data = site_contribution_data,
    horizon_summary = horizon_summary,
    curve_plot_data = curve_plot_data
  )
  outputs_upgrade <- upgrade_stage2_outputs_for_current_spec(outputs_upgrade, horizon_registry)

  site_admin_lookup <- outputs_upgrade$site_admin_lookup
  reverse_km_summary <- outputs_upgrade$reverse_km_summary
  numbers_at_risk_data <- outputs_upgrade$numbers_at_risk_data
  composition_data <- outputs_upgrade$composition_data
  site_contribution_data <- outputs_upgrade$site_contribution_data
  horizon_summary <- outputs_upgrade$horizon_summary
  curve_plot_data <- outputs_upgrade$curve_plot_data
}

# 🔴 Compute: follow-up maturity objects when cache is not reused ===============================
if (!isTRUE(reused_existing_core_results)) {
  group_output_list <- lapply(seq_len(nrow(group_registry)), function(i) {
    group_row <- group_registry[i, ]
    dataset_name <- group_row$dataset
    subgroup_kind <- group_row$subgroup_kind
    subgroup_value <- group_row$subgroup_value
    group_id <- group_row$group_id
    panel_label <- group_row$panel_label

    df_group <- filter_group_data(
      df = analysis_datasets_stage2[[dataset_name]],
      subgroup_kind = subgroup_kind,
      subgroup_value = subgroup_value
    )

    if (nrow(df_group) == 0L) {
      return(NULL)
    }

    horizon_registry_dataset <- get_horizon_registry_for_group(
      dataset_name = dataset_name,
      subgroup_kind = subgroup_kind,
      subgroup_value = subgroup_value,
      horizon_registry = horizon_registry,
      site_support_dataset_lookup = site_support_dataset_lookup
    )

    curve_bundle <- compute_curve_bundle(
      df = df_group,
      dataset_name = dataset_name,
      subgroup_kind = subgroup_kind,
      subgroup_value = subgroup_value,
      group_id = group_id,
      panel_label = panel_label,
      plot_time_grid_year = plot_time_grid_year
    )

    numbers_at_risk_df <- compute_numbers_at_risk(
      df = df_group,
      dataset_name = dataset_name,
      subgroup_kind = subgroup_kind,
      subgroup_value = subgroup_value,
      group_id = group_id,
      panel_label = panel_label,
      risk_table_times_year = risk_table_times_year,
      horizon_registry_dataset = horizon_registry_dataset
    )

    composition_df <- compute_composition_over_time(
      df = df_group,
      dataset_name = dataset_name,
      subgroup_kind = subgroup_kind,
      subgroup_value = subgroup_value,
      group_id = group_id,
      panel_label = panel_label,
      time_grid_year = plot_time_grid_year
    )

    horizon_summary_df <- compute_horizon_summary(
      df = df_group,
      dataset_name = dataset_name,
      subgroup_kind = subgroup_kind,
      subgroup_value = subgroup_value,
      group_id = group_id,
      panel_label = panel_label,
      horizon_registry_dataset = horizon_registry_dataset
    )

    site_contribution_df <- compute_site_contribution_profile(
      df = df_group,
      dataset_name = dataset_name,
      subgroup_kind = subgroup_kind,
      subgroup_value = subgroup_value,
      group_id = group_id,
      panel_label = panel_label,
      horizon_registry_dataset = horizon_registry_dataset
    )

    list(
      curves = curve_bundle$curves,
      reverse_summary = curve_bundle$reverse_summary,
      numbers_at_risk = numbers_at_risk_df,
      composition = composition_df,
      horizon_summary = horizon_summary_df,
      site_contribution = site_contribution_df
    )
  })

  group_output_list <- purrr::compact(group_output_list)

  if (length(group_output_list) == 0L) {
    stop('No non-empty Stage 2 group outputs were generated.', call. = FALSE)
  }

  survival_curve_data <- bind_rows(lapply(group_output_list, `[[`, 'curves')) %>%
    arrange(factor(dataset, levels = dataset_order), subgroup_kind, subgroup_value, summary_view, curve_role, time_year)

  reverse_km_summary <- bind_rows(lapply(group_output_list, `[[`, 'reverse_summary')) %>%
    distinct() %>%
    arrange(factor(dataset, levels = dataset_order), subgroup_kind, subgroup_value, summary_view)

  numbers_at_risk_data <- bind_rows(lapply(group_output_list, `[[`, 'numbers_at_risk')) %>%
    arrange(factor(dataset, levels = dataset_order), subgroup_kind, subgroup_value, time_year)

  composition_data <- bind_rows(lapply(group_output_list, `[[`, 'composition')) %>%
    arrange(factor(dataset, levels = dataset_order), subgroup_kind, subgroup_value, summary_view, time_year, state)

  site_contribution_data <- bind_rows(lapply(group_output_list, `[[`, 'site_contribution')) %>%
    arrange(factor(dataset, levels = dataset_order), subgroup_kind, subgroup_value, analysis_structure, horizon_year, site)

  site_contribution_summary <- summarise_site_contribution_profile(site_contribution_data)

  horizon_summary <- bind_rows(lapply(group_output_list, `[[`, 'horizon_summary')) %>%
    arrange(factor(dataset, levels = dataset_order), subgroup_kind, subgroup_value, analysis_structure, summary_view, horizon_year)

  curve_data_numbers <- numbers_at_risk_data %>%
    mutate(
      curve_role = 'numbers_at_risk',
      summary_view = NA_character_,
      estimate = as.numeric(n_risk),
      conf_low = NA_real_,
      conf_high = NA_real_,
      std_error = NA_real_,
      n_event = NA_real_,
      n_censor = NA_real_,
      state = NA_character_,
      state_n = NA_real_,
      state_prop = NA_real_
    ) %>%
    select(
      dataset, subgroup_kind, subgroup_value, group_id, panel_label,
      summary_view, curve_role, time_year,
      estimate, conf_low, conf_high, std_error, n_risk, n_event, n_censor,
      n_total, state, state_n, state_prop
    )

  curve_data_composition <- composition_data %>%
    mutate(
      curve_role = 'composition_state',
      estimate = state_prop,
      conf_low = NA_real_,
      conf_high = NA_real_,
      std_error = NA_real_,
      n_risk = NA_real_,
      n_event = NA_real_,
      n_censor = NA_real_
    ) %>%
    select(
      dataset, subgroup_kind, subgroup_value, group_id, panel_label,
      summary_view, curve_role, time_year,
      estimate, conf_low, conf_high, std_error, n_risk, n_event, n_censor,
      n_total, state, state_n, state_prop
    )

  curve_plot_data <- bind_rows(
    survival_curve_data %>% mutate(state = NA_character_, state_n = NA_real_, state_prop = NA_real_),
    curve_data_numbers,
    curve_data_composition
  ) %>%
    arrange(factor(dataset, levels = dataset_order), subgroup_kind, subgroup_value, curve_role, summary_view, time_year)

  site_admin_lookup <- site_admin_lookup_current
}

site_admin_lookup <- site_admin_lookup %>%
  mutate(
    site_admin_end_date = parse_date_flex(site_admin_end_date),
    site_admin_end_days_followup_max = as.numeric(site_admin_end_days_followup_max),
    site_admin_end_years_followup_max = as.numeric(site_admin_end_years_followup_max)
  )

horizon_summary <- horizon_summary %>%
  select(-any_of(c(
    'total_eligible_n', 'eligible_site_count', 'dominant_site', 'dominant_site_n', 'dominant_site_fraction',
    'site_dominance_status', 'site_dominance_note', 'site_dominance_warning_flag', 'site_dominance_warning_note',
    'preferred_for_reporting', 'reporting_priority', 'reporting_preference_note',
    'completeness_plot_status', 'completeness_plot_note', 'completeness_plot_alpha',
    'completeness_plot_masked', 'percentage_method_plot', 'CCI_plot', 'SPT_plot'
  )))

# 🔴 Derive: reusable long tables and validation products ===============================
site_contribution_summary <- summarise_site_contribution_profile(site_contribution_data)

horizon_summary <- horizon_summary %>%
  attach_site_contribution_to_horizon(site_contribution_summary) %>%
  add_completeness_display_fields() %>%
  add_reporting_preference_fields()

followup_maturity_long <- make_horizon_long(horizon_summary)
followup_side_by_side <- make_side_by_side_summary(horizon_summary)
completeness_plot_data <- build_completeness_plot_data(horizon_summary)

validation_outputs <- run_stage2_validation_audits(
  analysis_datasets_stage2 = analysis_datasets_stage2,
  group_registry = group_registry,
  site_admin_lookup = site_admin_lookup,
  merged_admin_end_date = merged_admin_end_date_current,
  numbers_at_risk_data = numbers_at_risk_data,
  composition_data = composition_data,
  horizon_summary = horizon_summary,
  reverse_km_summary = reverse_km_summary,
  site_contribution_summary = site_contribution_summary,
  site_support_dataset_lookup = site_support_dataset_lookup
)

validation_summary <- validation_outputs$validation_summary
validation_issue_rows <- validation_outputs$validation_issue_rows

readr::write_csv(validation_summary, file.path(export_path, 'stage2_validation_summary.csv'))
readr::write_csv(validation_issue_rows, file.path(export_path, 'stage2_validation_issue_rows.csv'))

validation_failures <- validation_summary %>%
  filter(status != 'pass')

if (nrow(validation_failures) > 0L) {
  stop(
    paste('Stage 2 validation failed for the following checks:', paste(validation_failures$check_name, collapse = ', ')),
    call. = FALSE
  )
}

quality_flags <- make_quality_flags(horizon_summary = horizon_summary, reverse_km_summary = reverse_km_summary)

stage2_metadata_registry <- make_stage2_metadata(
  stage1_inputs = stage1_inputs,
  site_admin_lookup = site_admin_lookup,
  plot_time_grid_year = plot_time_grid_year,
  threshold_vector = threshold_vector,
  merged_admin_end_date = merged_admin_end_date_current,
  main_plot_panel_order = main_plot_panel_order,
  sex_plot_panel_order = sex_plot_panel_order,
  data_path = data_path,
  export_path = export_path,
  reused_existing_core_results = reused_existing_core_results,
  core_reuse_mode = core_reuse_mode,
  raw_merged_file = raw_merged_file_used,
  stage1_source_signature = current_stage1_source_signature
)

# 🔴 Save: Stage 2 core outputs and bundle ===============================
readr::write_csv(site_admin_lookup, file.path(export_path, 'stage2_site_admin_end_dates.csv'))
readr::write_csv(reverse_km_summary, file.path(export_path, 'stage2_reverse_km_summary.csv'))
readr::write_csv(numbers_at_risk_data, file.path(export_path, 'stage2_numbers_at_risk.csv'))
readr::write_csv(composition_data, file.path(export_path, 'stage2_followup_composition.csv'))
readr::write_csv(site_contribution_data, file.path(export_path, 'stage2_merged_site_contribution.csv'))
readr::write_csv(horizon_summary, file.path(export_path, 'stage2_followup_horizon_summary.csv'))
readr::write_csv(followup_maturity_long, file.path(export_path, 'stage2_followup_maturity_long.csv'))
readr::write_csv(followup_side_by_side, file.path(export_path, 'stage2_followup_side_by_side.csv'))
readr::write_csv(curve_plot_data, file.path(export_path, 'stage2_followup_curve_data.csv'))
readr::write_csv(completeness_plot_data, file.path(export_path, 'stage2_followup_completeness_plot_data.csv'))
readr::write_csv(stage2_metadata_registry, file.path(export_path, 'stage2_metadata_registry.csv'))
readr::write_csv(quality_flags, file.path(export_path, 'stage2_quality_flags.csv'))

stage2_followup_bundle <- list(
  stage = 'Stage 2',
  created_at = as.character(Sys.time()),
  session_info = utils::sessionInfo(),
  config = list(
    data_path = data_path,
    export_path = export_path,
    stage1_analysis_datasets_file = stage1_analysis_datasets_file,
    stage1_backbone_bundle_file = stage1_backbone_bundle_file,
    plot_time_step_year = plot_time_step_year,
    plot_max_year = plot_max_year,
    risk_table_times_year = risk_table_times_year,
    threshold_vector = threshold_vector,
    completeness_sparse_risk_n = completeness_sparse_risk_n,
    completeness_very_sparse_risk_n = completeness_very_sparse_risk_n,
    completeness_low_fraction_cutoff = completeness_low_fraction_cutoff,
    export_sex_panel_plots = export_sex_panel_plots,
    site_dominance_warning_fraction = site_dominance_warning_fraction,
    site_dominance_warning_min_horizon_year = site_dominance_warning_min_horizon_year,
    site_dominance_warning_exempt_primary_supported = site_dominance_warning_exempt_primary_supported,
    patch_level = current_patch_level,
    scientific_compatibility_signature = scientific_compatibility_signature,
    stage1_source_signature = current_stage1_source_signature,
    reused_existing_core_results = reused_existing_core_results,
    core_reuse_mode = core_reuse_mode,
    raw_merged_file_for_date_recovery = raw_merged_file_used
  ),
  stage1_inputs = list(
    scaling_registry = scaling_registry_stage1,
    horizon_registry = horizon_registry,
    threshold_registry = stage1_inputs$threshold_registry,
    metadata_registry = stage1_inputs$metadata_registry,
    dataset_registry = dataset_registry,
    formula_registry = stage1_inputs$formula_registry
  ),
  stage2_scaling_registry = stage2_scaling_registry,
  site_admin_lookup = site_admin_lookup,
  merged_admin_end_date = merged_admin_end_date_current,
  site_support_dataset_lookup = site_support_dataset_lookup,
  group_registry = group_registry,
  main_plot_panel_order = main_plot_panel_order,
  sex_plot_panel_order = sex_plot_panel_order,
  reverse_km_summary = reverse_km_summary,
  numbers_at_risk_data = numbers_at_risk_data,
  composition_data = composition_data,
  site_contribution_data = site_contribution_data,
  site_contribution_summary = site_contribution_summary,
  horizon_summary = horizon_summary,
  followup_maturity_long = followup_maturity_long,
  followup_side_by_side = followup_side_by_side,
  curve_plot_data = curve_plot_data,
  completeness_plot_data = completeness_plot_data,
  validation_summary = validation_summary,
  validation_issue_rows = validation_issue_rows,
  quality_flags = quality_flags,
  metadata_registry = stage2_metadata_registry
)

saveRDS(stage2_followup_bundle, file.path(export_path, 'stage2_followup_bundle.rds'))

# 🔴 Render: bundled PDFs and individual PNG plots ===============================
main_plot_objects <- build_stage2_plot_objects(
  curve_plot_data = curve_plot_data,
  numbers_at_risk_data = numbers_at_risk_data,
  completeness_plot_data = completeness_plot_data,
  panel_order = main_plot_panel_order,
  subgroup_kind_filter = c('overall', 'site'),
  panel_family_label = 'overall and site panels'
)

main_plot_pdf <- file.path(export_path, 'stage2_followup_plots.pdf')
main_png_paths <- save_stage2_plot_collection(
  plot_objects = main_plot_objects,
  output_pdf_path = main_plot_pdf,
  png_prefix = 'stage2_plot_main'
)

sex_plot_pdf <- file.path(export_path, 'stage2_followup_plots_sex.pdf')
sex_png_paths <- character()

if (isTRUE(export_sex_panel_plots)) {
  sex_plot_objects <- build_stage2_plot_objects(
    curve_plot_data = curve_plot_data,
    numbers_at_risk_data = numbers_at_risk_data,
    completeness_plot_data = completeness_plot_data,
    panel_order = sex_plot_panel_order,
    subgroup_kind_filter = c('sex'),
    panel_family_label = 'sex panels'
  )

  sex_png_paths <- save_stage2_plot_collection(
    plot_objects = sex_plot_objects,
    output_pdf_path = sex_plot_pdf,
    png_prefix = 'stage2_plot_sex'
  )
}

# 🔴 Export: manifest and run log ===============================
export_manifest <- tibble(
  file_name = c(
    'stage2_site_admin_end_dates.csv',
    'stage2_reverse_km_summary.csv',
    'stage2_numbers_at_risk.csv',
    'stage2_followup_composition.csv',
    'stage2_merged_site_contribution.csv',
    'stage2_followup_horizon_summary.csv',
    'stage2_followup_maturity_long.csv',
    'stage2_followup_side_by_side.csv',
    'stage2_followup_curve_data.csv',
    'stage2_followup_completeness_plot_data.csv',
    'stage2_metadata_registry.csv',
    'stage2_validation_summary.csv',
    'stage2_validation_issue_rows.csv',
    'stage2_quality_flags.csv',
    'stage2_followup_bundle.rds',
    'stage2_followup_plots.pdf',
    'stage2_export_manifest.csv',
    'stage2_run_log.txt'
  ),
  object_name = c(
    'site_admin_lookup',
    'reverse_km_summary',
    'numbers_at_risk_data',
    'composition_data',
    'site_contribution_data',
    'horizon_summary',
    'followup_maturity_long',
    'followup_side_by_side',
    'curve_plot_data',
    'completeness_plot_data',
    'stage2_metadata_registry',
    'validation_summary',
    'validation_issue_rows',
    'quality_flags',
    'stage2_followup_bundle',
    'stage2_followup_plots',
    'export_manifest',
    'stage2_run_log'
  ),
  description = c(
    'Site-specific administrative end dates used for merged completeness recalculation',
    'Group-level reverse Kaplan-Meier follow-up summaries with max observed and potential follow-up windows',
    'Numbers-at-risk table at 0, 1, ..., 10 years with standardized support-display fields',
    'Event/censoring/remission composition over time',
    'Merged overall/sex horizon-level site contribution table showing eligible counts by site and descriptive dominance status',
    'Horizon-specific follow-up maturity table with support-tier, evidence, claim-restriction, site-dominance, plot-masking, and preferred-reporting fields',
    'Long-format horizon-level follow-up maturity table',
    'Transition-only and remission-sensitive horizon summaries saved side by side',
    'Curve-level source-of-truth data for KM, reverse KM, numbers-at-risk, and event/censoring/remission composition plots',
    'Explicit source-of-truth data frame for the completeness plot, derived from the horizon summary',
    'Stage 2 metadata registry',
    'Pass/fail summary of Stage 2 validation checks',
    'Issue-level rows for Stage 2 validation checks; empty when all validations pass',
    'Non-fatal sparse-tail, masking, reverse-KM CI, and merged site-dominance caution flags',
    'Reusable Stage 2 follow-up maturity bundle',
    'Stage 2 follow-up maturity plots for overall and site panels',
    'Manifest of all Stage 2 exported files',
    'Plain-text run log reporting the input/output paths, cache mode, and export confirmation'
  ),
  file_path = file.path(export_path, c(
    'stage2_site_admin_end_dates.csv',
    'stage2_reverse_km_summary.csv',
    'stage2_numbers_at_risk.csv',
    'stage2_followup_composition.csv',
    'stage2_merged_site_contribution.csv',
    'stage2_followup_horizon_summary.csv',
    'stage2_followup_maturity_long.csv',
    'stage2_followup_side_by_side.csv',
    'stage2_followup_curve_data.csv',
    'stage2_followup_completeness_plot_data.csv',
    'stage2_metadata_registry.csv',
    'stage2_validation_summary.csv',
    'stage2_validation_issue_rows.csv',
    'stage2_quality_flags.csv',
    'stage2_followup_bundle.rds',
    'stage2_followup_plots.pdf',
    'stage2_export_manifest.csv',
    'stage2_run_log.txt'
  )),
  file_path_relative = c(
    'stage2_site_admin_end_dates.csv',
    'stage2_reverse_km_summary.csv',
    'stage2_numbers_at_risk.csv',
    'stage2_followup_composition.csv',
    'stage2_merged_site_contribution.csv',
    'stage2_followup_horizon_summary.csv',
    'stage2_followup_maturity_long.csv',
    'stage2_followup_side_by_side.csv',
    'stage2_followup_curve_data.csv',
    'stage2_followup_completeness_plot_data.csv',
    'stage2_metadata_registry.csv',
    'stage2_validation_summary.csv',
    'stage2_validation_issue_rows.csv',
    'stage2_quality_flags.csv',
    'stage2_followup_bundle.rds',
    'stage2_followup_plots.pdf',
    'stage2_export_manifest.csv',
    'stage2_run_log.txt'
  ),
  counts_against_limit = TRUE
)

if (isTRUE(export_sex_panel_plots)) {
  export_manifest <- bind_rows(
    export_manifest,
    tibble(
      file_name = 'stage2_followup_plots_sex.pdf',
      object_name = 'stage2_followup_plots_sex',
      description = 'Stage 2 follow-up maturity plots for sex panels',
      file_path = file.path(export_path, 'stage2_followup_plots_sex.pdf'),
      file_path_relative = 'stage2_followup_plots_sex.pdf',
      counts_against_limit = TRUE
    )
  )
}

readr::write_csv(export_manifest, file.path(export_path, 'stage2_export_manifest.csv'))

run_log_lines <- c(
  sprintf('[Stage 2] Completed at: %s', as.character(Sys.time())),
  sprintf('[Stage 2] Input directory   : %s', data_path),
  sprintf('[Stage 2] Output directory  : %s', export_path),
  sprintf('[Stage 2] Patch level       : %s', current_patch_level),
  sprintf('[Stage 2] Scientific compat : %s', scientific_compatibility_signature),
  sprintf('[Stage 2] Reused core       : %s', if (isTRUE(reused_existing_core_results)) 'TRUE' else 'FALSE'),
  sprintf('[Stage 2] Reuse mode        : %s', core_reuse_mode),
  sprintf('[Stage 2] Raw merged file   : %s', raw_merged_file_used),
  sprintf('[Stage 2] Exported counted files: %d', sum(export_manifest$counts_against_limit, na.rm = TRUE)),
  sprintf('[Stage 2] Main PNG files    : %d', length(main_png_paths)),
  sprintf('[Stage 2] Sex PNG files     : %d', length(sex_png_paths)),
  paste(sprintf(' - %s', export_manifest$file_name), collapse = '\n'),
  if (length(c(main_png_paths, sex_png_paths)) > 0L) {
    paste(c('PNG exports (not counted toward the file limit):', paste(sprintf(' - %s', basename(c(main_png_paths, sex_png_paths))), collapse = '\n')), collapse = '\n')
  } else {
    'PNG exports (not counted toward the file limit): none'
  }
)
writeLines(run_log_lines, con = file.path(export_path, 'stage2_run_log.txt'))

confirm_stage2_exports(export_manifest = export_manifest)

message(sprintf('[Stage 2] Saved %d counted files to: %s', sum(export_manifest$counts_against_limit, na.rm = TRUE), export_path))
message(sprintf('[Stage 2] Saved %d PNG files (not counted toward the file limit).', length(c(main_png_paths, sex_png_paths))))
