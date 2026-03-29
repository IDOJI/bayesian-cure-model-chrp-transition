
# 🔴 Configure: paths, reuse policy, and spec guards ===============================
run_root_dir <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/stage6_Cure-appropriateness screening'
# If the completed Stage 6 run already lives in the folder you want to export into,
# you usually only need to change `export_path`. The script will auto-resolve
# `run_root_dir` from `export_path` whenever possible.
# Set `run_root_dir` explicitly only when the source completed run and the export
# folder are different.
export_path <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/stage6_Cure-appropriateness screening'

data_path <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/stage1_Backbone lock'

shared_master_spec_file <- "Integrated_Modeling_Master_Specification_English_REVISED_v5.md"
bayesian_companion_spec_file <- "Bayesian_Modeling_Specification_Stage8_REVISED_v5.md"
code_rules_file_candidates <- c(
  "Rules Before Generating R Code_🇬🇧ENG.md",
  "Rules Before Generating R Code_🇬🇧ENG_.md"
)
data_dictionary_file <- "3.Data Dictionary_🇬🇧ENG.md"

main_risk_scale <- "transition_only_main"
supplementary_risk_scale <- "transition_cif_competing"
screening_context <- "observational_cohort"
decision_rule_version <- "stage6_receus_aic_primary_hsu_familyset_xie_contradiction_v5"
alpha_screening <- 0.05
candidate_receus_families <- c("exponential", "weibull", "gamma", "loglogistic")
candidate_hsu_families <- c("weibull", "lognormal", "loglogistic")
hsu_familyset_method_name <- "hsu_supscore_familyset_approx"
hsu_family_set_name <- paste(candidate_hsu_families, collapse = "|")
receus_pi_threshold <- 0.025
receus_r_threshold <- 0.05
receus_pi_equivocal <- 0.010
receus_r_equivocal <- 0.10

refresh_stage6_exports <- TRUE
refresh_bootstrap_exports <- FALSE
reuse_existing_bootstrap_exports_if_present <- TRUE
rebuild_visual_summary_pdf <- TRUE
save_fit_objects <- TRUE
reuse_existing_completed_run_only <- TRUE

bootstrap_hist_bins <- 35L
bootstrap_plot_width <- 11
bootstrap_plot_height <- 8.5
bootstrap_density_adjust <- 1
core_plot_width <- 13
core_plot_height <- 8
core_plot_height_tall <- 10
plot_dpi <- 320

# 🔴 Initialize: packages and runtime options ===============================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

options(stringsAsFactors = FALSE, scipen = 999)

# 🔴 Define: file-system helpers ===============================
assert_exists <- function(path, label = basename(path)) {
  if (!file.exists(path) && !dir.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

make_temp_export_path <- function(path) {
  ext <- tools::file_ext(path)
  tempfile(
    pattern = paste0(tools::file_path_sans_ext(basename(path)), "_"),
    tmpdir = dirname(path),
    fileext = if (nzchar(ext)) paste0(".", ext) else ""
  )
}

replace_file_atomically <- function(temp_path, final_path) {
  if (file.exists(final_path)) {
    unlink(final_path)
  }
  renamed <- file.rename(temp_path, final_path)
  if (isTRUE(renamed)) {
    return(TRUE)
  }
  copied <- file.copy(temp_path, final_path, overwrite = TRUE)
  if (isTRUE(copied)) {
    unlink(temp_path)
    return(TRUE)
  }
  FALSE
}

safe_write_csv_atomic <- function(df, path) {
  temp_path <- make_temp_export_path(path)
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  readr::write_csv(df, temp_path)
  ok <- replace_file_atomically(temp_path, path)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to write CSV atomically: %s", path), call. = FALSE)
  }
  invisible(path)
}

safe_save_rds_atomic <- function(object, path) {
  temp_path <- make_temp_export_path(path)
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  saveRDS(object, temp_path)
  ok <- replace_file_atomically(temp_path, path)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to write RDS atomically: %s", path), call. = FALSE)
  }
  invisible(path)
}

safe_save_pdf_atomic <- function(path, width, height, plot_fun) {
  temp_path <- make_temp_export_path(path)
  pdf_open <- FALSE
  on.exit({
    if (pdf_open) {
      try(grDevices::dev.off(), silent = TRUE)
    }
    if (file.exists(temp_path)) {
      unlink(temp_path)
    }
  }, add = TRUE)

  grDevices::pdf(temp_path, width = width, height = height, onefile = TRUE)
  pdf_open <- TRUE
  plot_fun()
  grDevices::dev.off()
  pdf_open <- FALSE

  ok <- replace_file_atomically(temp_path, path)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to write PDF atomically: %s", path), call. = FALSE)
  }
  invisible(path)
}

safe_save_plot <- function(plot_object, path, width, height, dpi = plot_dpi, bg = "white") {
  temp_path <- make_temp_export_path(path)
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  ggplot2::ggsave(
    filename = temp_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = bg,
    limitsize = FALSE
  )
  ok <- replace_file_atomically(temp_path, path)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to write plot atomically: %s", path), call. = FALSE)
  }
  invisible(path)
}

safe_read_csv_if_exists <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
}

normalize_existing_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

normalize_blank_na_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

normalize_candidate_paths <- function(paths) {
  paths <- normalize_blank_na_chr(paths)
  paths <- paths[!is.na(paths) & nzchar(paths)]
  if (length(paths) == 0L) {
    return(character())
  }
  unique(normalize_existing_path(paths))
}

path_has_stage6_reuse_markers <- function(path) {
  path <- normalize_blank_na_chr(path)
  if (is.na(path) || !nzchar(path) || !dir.exists(path)) {
    return(FALSE)
  }

  observed_dir <- file.path(path, "observed")
  shards_dir <- file.path(path, "shards")

  has_observed_dir <- dir.exists(observed_dir)
  has_shards_dir <- dir.exists(shards_dir)
  has_bundle <- file.exists(file.path(path, "stage6_screening_bundle.rds")) ||
    file.exists(file.path(path, "state", "stage6_state_bundle.rds"))
  has_registry <- file.exists(file.path(path, "stage6_shard_registry.csv"))

  has_observed_rds <- has_observed_dir && length(list.files(observed_dir, pattern = "\\.rds$", full.names = TRUE)) > 0L
  has_shard_rds <- has_shards_dir && length(list.files(shards_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)) > 0L

  has_observed_rds && has_shard_rds && (has_bundle || has_registry)
}

resolve_stage6_run_root <- function(run_root_dir = NA_character_, export_path = NA_character_) {
  candidates <- normalize_candidate_paths(c(
    export_path,
    run_root_dir,
    if (!is.na(export_path) && nzchar(export_path)) dirname(export_path) else NA_character_,
    if (!is.na(run_root_dir) && nzchar(run_root_dir)) dirname(run_root_dir) else NA_character_
  ))

  valid_idx <- which(vapply(candidates, path_has_stage6_reuse_markers, logical(1)))
  if (length(valid_idx) == 0L) {
    stop(
      paste0(
        "Could not locate an existing completed Stage 6 run root. Checked candidates: ",
        paste(candidates, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  candidates[[valid_idx[[1]]]]
}

resolve_existing_file_name <- function(candidates, directory) {
  for (nm in unique(candidates)) {
    candidate_path <- file.path(directory, nm)
    if (file.exists(candidate_path)) {
      return(nm)
    }
  }
  NA_character_
}

# 🔴 Define: scalar and parsing helpers ===============================
as_integer_or_na <- function(x) {
  suppressWarnings(as.integer(as.character(x)))
}

as_numeric_or_na <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

coalesce_chr_scalar <- function(...) {
  xs <- list(...)
  for (ii in seq_along(xs)) {
    val <- xs[[ii]]
    if (length(val) == 0L || is.null(val)) {
      next
    }
    val <- as.character(val[[1]])
    if (!is.na(val) && nzchar(val)) {
      return(val)
    }
  }
  NA_character_
}

coalesce_num_scalar <- function(...) {
  xs <- list(...)
  for (ii in seq_along(xs)) {
    val <- xs[[ii]]
    if (length(val) == 0L || is.null(val)) {
      next
    }
    val <- suppressWarnings(as.numeric(val[[1]]))
    if (is.finite(val)) {
      return(val)
    }
  }
  NA_real_
}

coalesce_lgl_scalar <- function(...) {
  xs <- list(...)
  for (ii in seq_along(xs)) {
    val <- xs[[ii]]
    if (length(val) == 0L || is.null(val)) {
      next
    }
    if (is.logical(val)) {
      out <- val[[1]]
    } else {
      txt <- tolower(trimws(as.character(val[[1]])))
      out <- dplyr::case_when(
        txt %in% c("true", "t", "1", "yes", "y") ~ TRUE,
        txt %in% c("false", "f", "0", "no", "n") ~ FALSE,
        TRUE ~ NA
      )
    }
    if (!is.na(out)) {
      return(out)
    }
  }
  NA
}

safe_divide <- function(numerator, denominator) {
  ifelse(is.na(denominator) | denominator == 0, NA_real_, numerator / denominator)
}

extract_list_scalar <- function(x, name, default = NA) {
  if (is.null(x) || is.null(x[[name]]) || length(x[[name]]) == 0L) {
    return(default)
  }
  x[[name]][[1]]
}

extract_frame_cell <- function(df, column, default = NA) {
  if (is.null(df) || nrow(df) == 0L || !column %in% names(df)) {
    return(default)
  }
  df[[column]][[1]]
}

strip_formula_rhs <- function(formula_text) {
  if (is.null(formula_text) || length(formula_text) == 0L) {
    return(NA_character_)
  }
  formula_text <- as.character(formula_text[[1]])
  if (is.na(formula_text) || !nzchar(formula_text)) {
    return(NA_character_)
  }
  gsub("^\\s*~\\s*", "", formula_text)
}

safe_quantile <- function(x, probs) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(rep(NA_real_, length(probs)))
  }
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, type = 8, names = FALSE))
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  mean(x)
}

safe_sd <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) <= 1L) {
    return(NA_real_)
  }
  stats::sd(x)
}

safe_min <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  min(x)
}

safe_max <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  max(x)
}

append_note <- function(existing_note, addition) {
  existing_note <- as.character(existing_note)
  addition <- as.character(addition)
  if (is.na(addition) || !nzchar(addition)) {
    return(existing_note)
  }
  if (is.na(existing_note) || !nzchar(existing_note)) {
    return(addition)
  }
  paste(existing_note, addition, sep = " | ")
}

compact_character_set <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) {
    return(NA_character_)
  }
  paste(sort(x), collapse = "|")
}

normalize_flag_chr <- function(x, allowed = c("supportive", "equivocal", "unsupportive")) {
  out <- as.character(x)
  out[!out %in% allowed] <- NA_character_
  out
}

normalize_method_flag <- function(df) {
  df$method_flag <- normalize_flag_chr(df$method_flag)
  df
}

to_logical_column <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  txt <- tolower(trimws(as.character(x)))
  out <- dplyr::case_when(
    txt %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    txt %in% c("false", "f", "0", "no", "n") ~ FALSE,
    TRUE ~ NA
  )
  as.logical(out)
}

# 🔴 Define: stage-1 ingest and compatibility checks ===============================
read_stage1_inputs <- function(data_path) {
  stage1_bundle_file <- file.path(data_path, "stage1_backbone_bundle.rds")
  stage1_datasets_file <- file.path(data_path, "stage1_analysis_datasets.rds")
  stage1_dataset_registry_file <- file.path(data_path, "stage1_dataset_registry.csv")
  stage1_scaling_registry_file <- file.path(data_path, "stage1_scaling_registry.csv")
  stage1_formula_registry_file <- file.path(data_path, "stage1_formula_registry.csv")
  stage1_modeling_registry_file <- file.path(data_path, "stage1_modeling_registry.csv")
  stage1_horizon_registry_file <- file.path(data_path, "stage1_horizon_registry.csv")
  stage1_threshold_registry_file <- file.path(data_path, "stage1_threshold_registry.csv")
  stage1_metadata_registry_file <- file.path(data_path, "stage1_metadata_registry.csv")

  required_files <- c(
    stage1_bundle_file,
    stage1_datasets_file,
    stage1_dataset_registry_file,
    stage1_scaling_registry_file,
    stage1_formula_registry_file,
    stage1_modeling_registry_file,
    stage1_horizon_registry_file,
    stage1_threshold_registry_file,
    stage1_metadata_registry_file
  )
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0L) {
    stop(sprintf("Missing Stage 1 input file(s): %s", paste(missing_files, collapse = ", ")), call. = FALSE)
  }

  inputs <- list(
    backbone_bundle = readRDS(stage1_bundle_file),
    analysis_datasets = readRDS(stage1_datasets_file),
    dataset_registry = readr::read_csv(stage1_dataset_registry_file, show_col_types = FALSE, progress = FALSE),
    scaling_registry = readr::read_csv(stage1_scaling_registry_file, show_col_types = FALSE, progress = FALSE),
    formula_registry = readr::read_csv(stage1_formula_registry_file, show_col_types = FALSE, progress = FALSE),
    modeling_registry = readr::read_csv(stage1_modeling_registry_file, show_col_types = FALSE, progress = FALSE),
    horizon_registry = readr::read_csv(stage1_horizon_registry_file, show_col_types = FALSE, progress = FALSE),
    threshold_registry = readr::read_csv(stage1_threshold_registry_file, show_col_types = FALSE, progress = FALSE),
    metadata_registry = readr::read_csv(stage1_metadata_registry_file, show_col_types = FALSE, progress = FALSE),
    code_rules_file = resolve_existing_file_name(code_rules_file_candidates, dirname(stage1_metadata_registry_file))
  )

  metadata <- inputs$metadata_registry
  canonical_common_rules <- metadata %>%
    filter(metadata_name == "canonical_common_rules") %>%
    pull(metadata_value)
  governing_shared <- metadata %>%
    filter(metadata_name == "governing_shared_master_spec") %>%
    pull(metadata_value)

  master_candidates <- unique(c(canonical_common_rules, governing_shared))
  if (!(shared_master_spec_file %in% master_candidates)) {
    stop(
      sprintf(
        "Current Stage 1 outputs are not locked to `%s`. Found: %s",
        shared_master_spec_file,
        paste(master_candidates, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  horizon_vector <- sort(unique(as.integer(inputs$horizon_registry$horizon_year)))
  if (!identical(horizon_vector, 1:10)) {
    stop(sprintf("Stage 1 horizon registry is not the required 1:10 grid. Found: %s", paste(horizon_vector, collapse = ", ")), call. = FALSE)
  }

  inputs
}

check_stage6_reuse_compatibility <- function(stage1_inputs, screening_outputs) {
  dataset_registry <- stage1_inputs$dataset_registry
  followup_df <- screening_outputs$followup_sufficiency_summary
  if (is.null(followup_df) || nrow(followup_df) == 0L) {
    stop("Existing Stage 6 outputs do not contain follow-up summary data.", call. = FALSE)
  }

  expected_counts <- dataset_registry %>%
    transmute(dataset = as.character(dataset), n_rows = as.integer(n_rows))

  compare_tbl <- tibble::tibble(
    dataset = c("PNU", "SNU", "merged__site_free", "merged__site_adjusted"),
    source_dataset = c("PNU", "SNU", "merged", "merged"),
    observed_n_total = c(
      extract_frame_cell(followup_df[followup_df$dataset_key == "PNU", , drop = FALSE], "n_total", NA_real_),
      extract_frame_cell(followup_df[followup_df$dataset_key == "SNU", , drop = FALSE], "n_total", NA_real_),
      extract_frame_cell(followup_df[followup_df$dataset_key == "merged__site_free", , drop = FALSE], "n_total", NA_real_),
      extract_frame_cell(followup_df[followup_df$dataset_key == "merged__site_adjusted", , drop = FALSE], "n_total", NA_real_)
    )
  ) %>%
    left_join(expected_counts, by = c("source_dataset" = "dataset")) %>%
    mutate(matches_stage1 = as.integer(observed_n_total) == as.integer(n_rows))

  if (any(is.na(compare_tbl$matches_stage1)) || any(!compare_tbl$matches_stage1)) {
    stop(
      paste0(
        "Existing Stage 6 completed-run outputs are not compatible with the current Stage 1 backbone. ",
        paste(
          sprintf(
            "%s: observed_n_total=%s vs stage1_n_rows=%s",
            compare_tbl$dataset,
            compare_tbl$observed_n_total,
            compare_tbl$n_rows
          ),
          collapse = "; "
        )
      ),
      call. = FALSE
    )
  }

  invisible(compare_tbl)
}

# 🔴 Define: Stage 6 bundle readers ===============================
read_stage6_bundle <- function(run_root_dir) {
  screening_bundle_file <- file.path(run_root_dir, "stage6_screening_bundle.rds")
  stage6_fitted_objects_file <- file.path(run_root_dir, "stage6_fitted_objects.rds")
  stage6_shard_registry_file <- file.path(run_root_dir, "stage6_shard_registry.csv")
  state_bundle_file <- file.path(run_root_dir, "state", "stage6_state_bundle.rds")
  observed_dir <- file.path(run_root_dir, "observed")
  shards_dir <- file.path(run_root_dir, "shards")

  assert_exists(run_root_dir, "run_root_dir")
  assert_exists(observed_dir, "observed_dir")
  assert_exists(shards_dir, "shards_dir")

  screening_bundle <- if (file.exists(screening_bundle_file)) readRDS(screening_bundle_file) else NULL
  state_bundle <- if (file.exists(state_bundle_file)) readRDS(state_bundle_file) else screening_bundle

  if (is.null(screening_bundle) && is.null(state_bundle)) {
    stop(
      paste0(
        "No completed Stage 6 bundle could be found. This revised script is designed to reuse an existing completed run. ",
        "Because `reuse_existing_completed_run_only = TRUE`, a brand-new full rerun is intentionally not attempted here."
      ),
      call. = FALSE
    )
  }

  shard_registry <- if (file.exists(stage6_shard_registry_file)) {
    readr::read_csv(stage6_shard_registry_file, show_col_types = FALSE, progress = FALSE)
  } else {
    tibble::tibble()
  }
  fit_objects <- if (file.exists(stage6_fitted_objects_file)) readRDS(stage6_fitted_objects_file) else NULL

  list(
    screening_bundle = screening_bundle,
    state_bundle = state_bundle,
    fit_objects = fit_objects,
    shard_registry = shard_registry,
    observed_dir = observed_dir,
    shards_dir = shards_dir,
    screening_bundle_file = screening_bundle_file,
    stage6_fitted_objects_file = stage6_fitted_objects_file,
    stage6_shard_registry_file = stage6_shard_registry_file
  )
}

extract_stage6_outputs <- function(stage6_objects) {
  screening_bundle <- stage6_objects$screening_bundle
  state_bundle <- stage6_objects$state_bundle

  variant_registry <- NULL
  metadata_registry <- NULL
  followup_sufficiency_summary <- NULL
  screening_method_results <- NULL
  receus_candidate_fits <- NULL
  screening_summary <- NULL
  carry_forward_flag_table <- NULL

  if (!is.null(screening_bundle)) {
    if (!is.null(screening_bundle$registries)) {
      variant_registry <- screening_bundle$registries$variant_registry
      metadata_registry <- screening_bundle$registries$stage6_metadata_registry
      if (is.null(metadata_registry)) {
        metadata_registry <- screening_bundle$registries$metadata_registry
      }
    }
    if (!is.null(screening_bundle$outputs)) {
      followup_sufficiency_summary <- screening_bundle$outputs$followup_sufficiency_summary
      screening_method_results <- screening_bundle$outputs$screening_method_results
      receus_candidate_fits <- screening_bundle$outputs$receus_candidate_fits
      screening_summary <- screening_bundle$outputs$screening_summary
      carry_forward_flag_table <- screening_bundle$outputs$carry_forward_flag_table
    }
  }

  if (is.null(variant_registry) && !is.null(state_bundle) && !is.null(state_bundle$registries)) {
    variant_registry <- state_bundle$registries$variant_registry
  }

  if (is.null(followup_sufficiency_summary) && !is.null(state_bundle) && !is.null(state_bundle$observed_outputs)) {
    followup_sufficiency_summary <- state_bundle$observed_outputs$followup_summary_base
    receus_candidate_fits <- state_bundle$observed_outputs$receus_candidate_base
  }

  list(
    variant_registry = tibble::as_tibble(variant_registry),
    metadata_registry = tibble::as_tibble(metadata_registry),
    followup_sufficiency_summary = tibble::as_tibble(followup_sufficiency_summary),
    screening_method_results = tibble::as_tibble(screening_method_results),
    receus_candidate_fits = tibble::as_tibble(receus_candidate_fits),
    screening_summary = tibble::as_tibble(screening_summary),
    carry_forward_flag_table = tibble::as_tibble(carry_forward_flag_table)
  )
}

# 🔴 Define: Stage 6 table patchers ===============================
dataset_order_full <- c("PNU", "SNU", "merged__site_free", "merged__site_adjusted")
method_order <- c(
  "maller_zhou_dn",
  "maller_zhou_alpha_n",
  "xie_sufficient_followup",
  "receus_weibull",
  "receus_aic",
  "hsu_supscore_weibull_approx",
  "hsu_supscore_lognormal_approx",
  "hsu_supscore_loglogistic_approx",
  hsu_familyset_method_name
)

canonicalize_variant_registry <- function(df, stage1_inputs) {
  if (nrow(df) == 0L) {
    source_lookup <- stage1_inputs$backbone_bundle$source_lookup
    formula_registry <- stage1_inputs$formula_registry
    get_hsu_formula_spec <- function(dataset_key) {
      if (dataset_key %in% c("PNU", "SNU")) {
        dataset_name <- dataset_key
        latency_rhs <- formula_registry %>%
          filter(dataset == dataset_name, formula_name == "base") %>%
          pull(formula_rhs) %>%
          unique()
        incidence_rhs_set <- formula_registry %>%
          filter(dataset == dataset_name) %>%
          pull(formula_rhs) %>%
          unique()
        tibble::tibble(
          dataset_key = dataset_key,
          source_dataset = dataset_name,
          source_description = unname(source_lookup[[dataset_name]]),
          analysis_variant = "single_cohort",
          screening_context = screening_context,
          stage1_dataset_name = dataset_name,
          tail_metric_source = "self",
          hsu_formula_branch = "site_free",
          hsu_latency_formula_rhs = latency_rhs[[1]],
          hsu_incidence_formula_rhs_set = paste(sort(unique(c("age_s + sex_num", incidence_rhs_set))), collapse = "|"),
          hsu_incidence_formula_count = length(unique(c("age_s + sex_num", incidence_rhs_set))),
          hsu_working_family_set = hsu_family_set_name,
          screening_variant_note = ifelse(dataset_key == "PNU", "Direct single-cohort screening on Stage 1 PNU data.", "Direct single-cohort screening on Stage 1 SNU data.")
        )
      } else {
        latency_rhs <- formula_registry %>%
          filter(dataset == "merged", formula_name == ifelse(dataset_key == "merged__site_adjusted", "site_added", "base")) %>%
          pull(formula_rhs) %>%
          unique()
        incidence_rhs_set <- formula_registry %>%
          filter(dataset == "merged", site_branch == ifelse(dataset_key == "merged__site_adjusted", "site_adjusted", "site_free")) %>%
          pull(formula_rhs) %>%
          unique()
        tibble::tibble(
          dataset_key = dataset_key,
          source_dataset = "merged",
          source_description = unname(source_lookup[["merged"]]),
          analysis_variant = ifelse(dataset_key == "merged__site_adjusted", "site_adjusted", "site_free"),
          screening_context = screening_context,
          stage1_dataset_name = "merged",
          tail_metric_source = ifelse(dataset_key == "merged__site_adjusted", "merged__site_free", "self"),
          hsu_formula_branch = ifelse(dataset_key == "merged__site_adjusted", "site_adjusted", "site_free"),
          hsu_latency_formula_rhs = latency_rhs[[1]],
          hsu_incidence_formula_rhs_set = paste(sort(unique(c(ifelse(dataset_key == "merged__site_adjusted", "age_s + sex_num + site", "age_s + sex_num"), incidence_rhs_set))), collapse = "|"),
          hsu_incidence_formula_count = length(unique(incidence_rhs_set)),
          hsu_working_family_set = hsu_family_set_name,
          screening_variant_note = ifelse(
            dataset_key == "merged__site_adjusted",
            "Raw-tail metrics inherit the merged site-free data; only the HSU working-model branch is rerun with site adjustment.",
            "Direct merged-data screening for the site-free structural branch."
          )
        )
      }
    }

    out <- bind_rows(lapply(dataset_order_full, get_hsu_formula_spec))
    return(out %>% arrange(match(dataset_key, dataset_order_full)))
  }

  for (nm in setdiff(c("screening_context", "dataset_key", "source_dataset", "analysis_variant"), names(df))) {
    df[[nm]] <- NA
  }

  df %>%
    mutate(
      screening_context = coalesce(as.character(screening_context), .env$screening_context),
      dataset_key = as.character(dataset_key),
      source_dataset = as.character(source_dataset),
      analysis_variant = as.character(analysis_variant)
    ) %>%
    arrange(match(dataset_key, dataset_order_full))
}

patch_followup_sufficiency_summary <- function(df) {
  out <- tibble::as_tibble(df)
  if (nrow(out) == 0L) {
    stop("Stage 6 follow-up summary is missing and cannot be reconstructed from the available completed-run objects.", call. = FALSE)
  }
  if (!"summary_note" %in% names(out)) {
    out$summary_note <- NA_character_
  }
  out %>%
    mutate(
      dataset_key = as.character(dataset_key),
      source_dataset = as.character(source_dataset),
      analysis_variant = as.character(analysis_variant),
      screening_context = .env$screening_context
    ) %>%
    arrange(match(dataset_key, dataset_order_full))
}

patch_screening_method_results <- function(df) {
  out <- tibble::as_tibble(df)
  if (nrow(out) == 0L) {
    stop("Stage 6 screening method results are missing and cannot be reconstructed from the available completed-run objects.", call. = FALSE)
  }

  required_columns <- c(
    "dataset_key", "source_dataset", "source_description", "analysis_variant", "screening_context",
    "primary_gate_method", "primary_gate_flag", "receus_aic_flag", "cure_model_eligibility_flag",
    "final_decision_flag", "receus_primary_class", "presence_modifier_flag", "cure_presence_support_flag",
    "presence_support_flag", "followup_contradiction_flag", "followup_not_contradicted_flag",
    "xie_centered_bootstrap_p_value", "descriptive_tail_summary_flag", "supporting_methods",
    "contradicting_methods", "screening_note", "common_horizon_vector", "common_threshold_vector",
    "stage6_final_class", "carry_forward_stage8"
  )
  for (nm in setdiff(required_columns, names(out))) {
    out[[nm]] <- NA
  }

  out <- out %>%
    mutate(
      dataset_key = as.character(dataset_key),
      source_dataset = as.character(source_dataset),
      analysis_variant = as.character(analysis_variant),
      method_name = as.character(method_name),
      method_role = as.character(method_role),
      model_family = as.character(model_family),
      selected_family = as.character(selected_family),
      family_set_name = as.character(family_set_name),
      latency_formula_rhs = as.character(latency_formula_rhs),
      incidence_formula_rhs = as.character(incidence_formula_rhs),
      implementation_label = as.character(implementation_label),
      note = as.character(note)
    )

  out <- normalize_method_flag(out)

  required_methods <- method_order
  missing_methods <- setdiff(required_methods, unique(out$method_name))
  if (length(missing_methods) > 0L) {
    stop(
      sprintf(
        "Existing Stage 6 method results are missing required method rows under the revised specification: %s",
        paste(missing_methods, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  out %>% arrange(match(dataset_key, dataset_order_full), match(method_name, method_order))
}

patch_receus_candidate_fits <- function(df) {
  out <- tibble::as_tibble(df)
  if (nrow(out) == 0L) {
    stop("Stage 6 RECeUS candidate-family summary is missing and cannot be reconstructed from the available completed-run objects.", call. = FALSE)
  }

  if (!"selected_by_aic" %in% names(out)) {
    out$selected_by_aic <- FALSE
  }

  out %>%
    mutate(
      dataset_key = as.character(dataset_key),
      source_dataset = as.character(source_dataset),
      analysis_variant = as.character(analysis_variant),
      family = as.character(family),
      receus_flag = normalize_flag_chr(receus_flag),
      selected_by_aic = to_logical_column(selected_by_aic),
      note = as.character(note)
    ) %>%
    arrange(match(dataset_key, dataset_order_full), match(family, candidate_receus_families))
}

patch_carry_forward_flag_table <- function(df, stage1_inputs, variant_registry, screening_method_results, followup_sufficiency_summary) {
  out <- tibble::as_tibble(df)

  if (nrow(out) == 0L) {
    stop("Stage 6 carry-forward table is missing and cannot be reconstructed from the available completed-run objects.", call. = FALSE)
  }

  required_columns <- c(
    "dataset_key", "source_dataset", "source_description", "analysis_variant", "screening_context",
    "decision_rule_version",
    "primary_gate_method", "primary_gate_flag", "receus_aic_flag", "cure_model_eligibility_flag",
    "final_decision_flag", "receus_primary_class", "presence_modifier_flag", "cure_presence_support_flag",
    "presence_support_flag", "followup_contradiction_flag", "followup_not_contradicted_flag",
    "xie_centered_bootstrap_p_value", "descriptive_tail_summary_flag", "supporting_methods",
    "contradicting_methods", "screening_note", "common_horizon_vector", "common_threshold_vector",
    "stage6_final_class", "carry_forward_stage8"
  )
  for (nm in setdiff(required_columns, names(out))) {
    out[[nm]] <- NA
  }

  contradiction_from_xie <- if ("xie_centered_bootstrap_p_value" %in% names(out)) {
    to_logical_column(as_numeric_or_na(out$xie_centered_bootstrap_p_value) < alpha_screening)
  } else {
    rep(NA, nrow(out))
  }

  out <- out %>%
    mutate(
      dataset_key = as.character(dataset_key),
      source_dataset = as.character(source_dataset),
      source_description = as.character(source_description),
      analysis_variant = as.character(analysis_variant),
      screening_context = coalesce(as.character(screening_context), .env$screening_context),
      decision_rule_version = .env$decision_rule_version,
      primary_gate_method = coalesce(as.character(primary_gate_method), "RECeUS-AIC"),
      primary_gate_flag = normalize_flag_chr(coalesce(as.character(primary_gate_flag), as.character(receus_aic_flag))),
      cure_model_eligibility_flag = normalize_flag_chr(coalesce(as.character(cure_model_eligibility_flag), as.character(final_decision_flag))),
      receus_primary_class = normalize_flag_chr(coalesce(as.character(receus_primary_class), as.character(receus_aic_flag), as.character(primary_gate_flag))),
      presence_modifier_flag = coalesce(to_logical_column(presence_modifier_flag), to_logical_column(presence_support_flag)),
      cure_presence_support_flag = coalesce(to_logical_column(cure_presence_support_flag), to_logical_column(presence_support_flag)),
      followup_contradiction_flag = coalesce(to_logical_column(followup_contradiction_flag), contradiction_from_xie),
      followup_not_contradicted_flag = coalesce(to_logical_column(followup_not_contradicted_flag), ifelse(is.na(followup_contradiction_flag), NA, !to_logical_column(followup_contradiction_flag))),
      descriptive_tail_summary_flag = coalesce(as.character(descriptive_tail_summary_flag), "descriptive_only"),
      supporting_methods = as.character(supporting_methods),
      contradicting_methods = as.character(contradicting_methods),
      screening_note = as.character(screening_note),
      common_horizon_vector = coalesce(as.character(common_horizon_vector), paste(sort(unique(as.integer(stage1_inputs$horizon_registry$horizon_year))), collapse = "|")),
      common_threshold_vector = coalesce(as.character(common_threshold_vector), paste(format(sort(unique(as.numeric(stage1_inputs$threshold_registry$threshold))), trim = TRUE, scientific = FALSE), collapse = "|")),
      final_decision_flag = cure_model_eligibility_flag,
      stage6_final_class = cure_model_eligibility_flag,
      carry_forward_stage8 = cure_model_eligibility_flag != "unsupportive"
    )

  if (!"presence_support_flag" %in% names(out)) {
    out$presence_support_flag <- out$cure_presence_support_flag
  } else {
    out$presence_support_flag <- coalesce(to_logical_column(out$presence_support_flag), out$cure_presence_support_flag)
  }

  if (!"followup_support_flag" %in% names(out)) {
    out$followup_support_flag <- dplyr::case_when(
      out$primary_gate_flag == "supportive" & !is.na(out$followup_not_contradicted_flag) & out$followup_not_contradicted_flag ~ TRUE,
      out$primary_gate_flag == "supportive" & !is.na(out$followup_not_contradicted_flag) & !out$followup_not_contradicted_flag ~ FALSE,
      out$primary_gate_flag == "unsupportive" ~ FALSE,
      TRUE ~ NA
    )
  }

  out %>% arrange(match(dataset_key, dataset_order_full))
}

patch_screening_summary <- function(df, carry_forward_flag_table, screening_method_results, followup_sufficiency_summary, variant_registry) {
  out <- tibble::as_tibble(df)
  if (nrow(out) == 0L) {
    build_row <- function(dataset_key) {
      carry_row <- carry_forward_flag_table[carry_forward_flag_table$dataset_key == dataset_key, , drop = FALSE]
      follow_row <- followup_sufficiency_summary[followup_sufficiency_summary$dataset_key == dataset_key, , drop = FALSE]
      method_subset <- screening_method_results[screening_method_results$dataset_key == dataset_key, , drop = FALSE]
      tibble::tibble(
        dataset_key = dataset_key,
        source_dataset = extract_frame_cell(carry_row, "source_dataset", NA_character_),
        source_description = extract_frame_cell(carry_row, "source_description", NA_character_),
        analysis_variant = extract_frame_cell(carry_row, "analysis_variant", NA_character_),
        screening_context = screening_context,
        tail_metric_source = extract_frame_cell(carry_row, "tail_metric_source", NA_character_),
        hsu_formula_branch = extract_frame_cell(carry_row, "hsu_formula_branch", NA_character_),
        decision_rule_version = decision_rule_version,
        primary_gate_method = extract_frame_cell(carry_row, "primary_gate_method", NA_character_),
        primary_gate_flag = extract_frame_cell(carry_row, "primary_gate_flag", NA_character_),
        presence_modifier_flag = extract_frame_cell(carry_row, "presence_modifier_flag", NA),
        followup_contradiction_flag = extract_frame_cell(carry_row, "followup_contradiction_flag", NA),
        descriptive_tail_summary_flag = extract_frame_cell(carry_row, "descriptive_tail_summary_flag", NA_character_),
        modifier_evidence_label = extract_frame_cell(carry_row, "modifier_evidence_label", NA_character_),
        n_total = extract_frame_cell(follow_row, "n_total", NA_real_),
        n_event = extract_frame_cell(follow_row, "n_event", NA_real_),
        n_censor_main = extract_frame_cell(follow_row, "n_censor_main", NA_real_),
        n_right_censor = extract_frame_cell(follow_row, "n_right_censor", NA_real_),
        n_remission = extract_frame_cell(follow_row, "n_remission", NA_real_),
        censor_rate_main = extract_frame_cell(follow_row, "censor_rate_main", NA_real_),
        person_time_year = extract_frame_cell(follow_row, "person_time_year", NA_real_),
        tau_year = extract_frame_cell(follow_row, "tau_year", NA_real_),
        largest_event_time_year = extract_frame_cell(follow_row, "largest_event_time_year", NA_real_),
        plateau_length_year = extract_frame_cell(follow_row, "plateau_length_year", NA_real_),
        n_censored_after_last_event = extract_frame_cell(follow_row, "n_censored_after_last_event", NA_real_),
        km_tail_survival = extract_frame_cell(follow_row, "km_tail_survival", NA_real_),
        km_tail_noncure = extract_frame_cell(follow_row, "km_tail_noncure", NA_real_),
        mz_interval_lower_year = extract_frame_cell(follow_row, "mz_interval_lower_year", NA_real_),
        mz_interval_event_count = extract_frame_cell(follow_row, "mz_interval_event_count", NA_real_),
        maller_zhou_alpha_n = extract_frame_cell(method_subset[method_subset$method_name == "maller_zhou_alpha_n", , drop = FALSE], "statistic_value", NA_real_),
        maller_zhou_alpha_flag = extract_frame_cell(method_subset[method_subset$method_name == "maller_zhou_alpha_n", , drop = FALSE], "method_flag", NA_character_),
        maller_zhou_q_n = extract_frame_cell(follow_row, "mz_q_n", NA_real_),
        maller_zhou_dn_stat = extract_frame_cell(method_subset[method_subset$method_name == "maller_zhou_dn", , drop = FALSE], "statistic_value", NA_real_),
        maller_zhou_dn_p_value = extract_frame_cell(method_subset[method_subset$method_name == "maller_zhou_dn", , drop = FALSE], "p_value", NA_real_),
        maller_zhou_dn_flag = extract_frame_cell(method_subset[method_subset$method_name == "maller_zhou_dn", , drop = FALSE], "method_flag", NA_character_),
        xie_epsilon = extract_frame_cell(follow_row, "xie_epsilon", NA_real_),
        xie_p_hat_n = extract_frame_cell(follow_row, "xie_p_hat_n", NA_real_),
        xie_p_hat_G = extract_frame_cell(follow_row, "xie_p_hat_G", NA_real_),
        xie_T_n = extract_frame_cell(method_subset[method_subset$method_name == "xie_sufficient_followup", , drop = FALSE], "statistic_value", NA_real_),
        xie_centered_bootstrap_p_value = extract_frame_cell(method_subset[method_subset$method_name == "xie_sufficient_followup", , drop = FALSE], "p_value", NA_real_),
        xie_flag = extract_frame_cell(method_subset[method_subset$method_name == "xie_sufficient_followup", , drop = FALSE], "method_flag", NA_character_),
        receus_weibull_family = extract_frame_cell(method_subset[method_subset$method_name == "receus_weibull", , drop = FALSE], "selected_family", NA_character_),
        receus_weibull_aic = extract_frame_cell(method_subset[method_subset$method_name == "receus_weibull", , drop = FALSE], "aic", NA_real_),
        receus_weibull_cure_fraction = extract_frame_cell(method_subset[method_subset$method_name == "receus_weibull", , drop = FALSE], "cure_fraction_hat", NA_real_),
        receus_weibull_ratio = extract_frame_cell(method_subset[method_subset$method_name == "receus_weibull", , drop = FALSE], "receus_ratio_hat", NA_real_),
        receus_weibull_flag = extract_frame_cell(method_subset[method_subset$method_name == "receus_weibull", , drop = FALSE], "method_flag", NA_character_),
        receus_aic_selected_family = extract_frame_cell(method_subset[method_subset$method_name == "receus_aic", , drop = FALSE], "selected_family", NA_character_),
        receus_aic_aic = extract_frame_cell(method_subset[method_subset$method_name == "receus_aic", , drop = FALSE], "aic", NA_real_),
        receus_aic_cure_fraction = extract_frame_cell(method_subset[method_subset$method_name == "receus_aic", , drop = FALSE], "cure_fraction_hat", NA_real_),
        receus_aic_ratio = extract_frame_cell(method_subset[method_subset$method_name == "receus_aic", , drop = FALSE], "receus_ratio_hat", NA_real_),
        receus_aic_flag = extract_frame_cell(method_subset[method_subset$method_name == "receus_aic", , drop = FALSE], "method_flag", NA_character_),
        hsu_weibull_stat = extract_frame_cell(method_subset[method_subset$method_name == "hsu_supscore_weibull_approx", , drop = FALSE], "statistic_value", NA_real_),
        hsu_weibull_p_value = extract_frame_cell(method_subset[method_subset$method_name == "hsu_supscore_weibull_approx", , drop = FALSE], "p_value", NA_real_),
        hsu_weibull_adjusted_p_value = extract_frame_cell(method_subset[method_subset$method_name == "hsu_supscore_weibull_approx", , drop = FALSE], "adjusted_p_value", NA_real_),
        hsu_weibull_flag = extract_frame_cell(method_subset[method_subset$method_name == "hsu_supscore_weibull_approx", , drop = FALSE], "method_flag", NA_character_),
        hsu_lognormal_stat = extract_frame_cell(method_subset[method_subset$method_name == "hsu_supscore_lognormal_approx", , drop = FALSE], "statistic_value", NA_real_),
        hsu_lognormal_p_value = extract_frame_cell(method_subset[method_subset$method_name == "hsu_supscore_lognormal_approx", , drop = FALSE], "p_value", NA_real_),
        hsu_lognormal_adjusted_p_value = extract_frame_cell(method_subset[method_subset$method_name == "hsu_supscore_lognormal_approx", , drop = FALSE], "adjusted_p_value", NA_real_),
        hsu_lognormal_flag = extract_frame_cell(method_subset[method_subset$method_name == "hsu_supscore_lognormal_approx", , drop = FALSE], "method_flag", NA_character_),
        hsu_loglogistic_stat = extract_frame_cell(method_subset[method_subset$method_name == "hsu_supscore_loglogistic_approx", , drop = FALSE], "statistic_value", NA_real_),
        hsu_loglogistic_p_value = extract_frame_cell(method_subset[method_subset$method_name == "hsu_supscore_loglogistic_approx", , drop = FALSE], "p_value", NA_real_),
        hsu_loglogistic_adjusted_p_value = extract_frame_cell(method_subset[method_subset$method_name == "hsu_supscore_loglogistic_approx", , drop = FALSE], "adjusted_p_value", NA_real_),
        hsu_loglogistic_flag = extract_frame_cell(method_subset[method_subset$method_name == "hsu_supscore_loglogistic_approx", , drop = FALSE], "method_flag", NA_character_),
        hsu_familyset_stat = extract_frame_cell(method_subset[method_subset$method_name == hsu_familyset_method_name, , drop = FALSE], "statistic_value", NA_real_),
        hsu_familyset_p_value = extract_frame_cell(method_subset[method_subset$method_name == hsu_familyset_method_name, , drop = FALSE], "p_value", NA_real_),
        hsu_familyset_flag = extract_frame_cell(method_subset[method_subset$method_name == hsu_familyset_method_name, , drop = FALSE], "method_flag", NA_character_),
        hsu_familyset_name = extract_frame_cell(method_subset[method_subset$method_name == hsu_familyset_method_name, , drop = FALSE], "family_set_name", NA_character_),
        hsu_familyset_selected_family = extract_frame_cell(method_subset[method_subset$method_name == hsu_familyset_method_name, , drop = FALSE], "selected_family", NA_character_),
        hsu_latency_formula_rhs = extract_frame_cell(method_subset[method_subset$method_name == hsu_familyset_method_name, , drop = FALSE], "latency_formula_rhs", NA_character_),
        hsu_selected_incidence_formula_rhs = extract_frame_cell(method_subset[method_subset$method_name == hsu_familyset_method_name, , drop = FALSE], "incidence_formula_rhs", NA_character_),
        cure_model_eligibility_flag = extract_frame_cell(carry_row, "cure_model_eligibility_flag", NA_character_),
        final_decision_flag = extract_frame_cell(carry_row, "cure_model_eligibility_flag", NA_character_),
        presence_evidence_label = extract_frame_cell(carry_row, "presence_evidence_label", NA_character_),
        followup_evidence_label = extract_frame_cell(carry_row, "followup_evidence_label", NA_character_),
        supporting_methods = extract_frame_cell(carry_row, "supporting_methods", NA_character_),
        contradicting_methods = extract_frame_cell(carry_row, "contradicting_methods", NA_character_),
        followup_summary_note = extract_frame_cell(follow_row, "summary_note", NA_character_),
        screening_note = extract_frame_cell(carry_row, "screening_note", NA_character_)
      )
    }
    out <- bind_rows(lapply(dataset_order_full, build_row))
  }

  for (nm in setdiff(c("screening_context", "decision_rule_version", "dataset_key", "source_dataset", "source_description", "analysis_variant"), names(out))) {
    out[[nm]] <- NA
  }

  out %>%
    mutate(
      dataset_key = as.character(dataset_key),
      source_dataset = as.character(source_dataset),
      source_description = as.character(source_description),
      analysis_variant = as.character(analysis_variant),
      screening_context = coalesce(as.character(screening_context), .env$screening_context),
      decision_rule_version = .env$decision_rule_version
    ) %>%
    arrange(match(dataset_key, dataset_order_full))
}

patch_stage6_metadata_registry <- function(df, stage1_inputs, reuse_completed_run = TRUE) {
  out <- tibble::as_tibble(df)
  if (nrow(out) == 0L) {
    out <- tibble::tibble(metadata_group = character(), metadata_name = character(), metadata_value = character())
  }

  code_rules_file <- coalesce_chr_scalar(stage1_inputs$code_rules_file, code_rules_file_candidates[[1]])

  patch_tbl <- tibble::tribble(
    ~metadata_group, ~metadata_name, ~metadata_value,
    "stage", "stage_name", "Stage 6 cure appropriateness screening",
    "stage", "stage_role", "Frequentist-only cure appropriateness screening and carry-forward eligibility gate",
    "stage", "screening_repeated_in_stage8", "FALSE",
    "stage", "spec_patch_level", "revised_v5_stage1_aligned_completed_run_reuse_png_split",
    "documents", "canonical_common_rules", shared_master_spec_file,
    "documents", "canonical_framework", shared_master_spec_file,
    "documents", "governing_shared_master_spec", shared_master_spec_file,
    "documents", "governing_stage8_bayesian_spec", bayesian_companion_spec_file,
    "documents", "governing_code_rules", code_rules_file,
    "documents", "governing_data_dictionary", data_dictionary_file,
    "documents", "stage6_specification_file", shared_master_spec_file,
    "documents", "stage6_specification_section", "6.6",
    "screening", "alpha_screening", as.character(alpha_screening),
    "screening", "xie_bootstrap_reps", as.character(NA),
    "screening", "hsu_bootstrap_reps", as.character(NA),
    "screening", "receus_candidate_families", paste(candidate_receus_families, collapse = "|"),
    "screening", "hsu_candidate_families", paste(candidate_hsu_families, collapse = "|"),
    "screening", "hsu_family_set_name", hsu_family_set_name,
    "screening", "decision_rule_version", decision_rule_version,
    "screening", "primary_gate_method", "RECeUS-AIC",
    "screening", "screening_context", screening_context,
    "screening", "reuse_existing_completed_run_only", as.character(reuse_completed_run),
    "screening", "reuse_existing_bootstrap_exports_if_present", as.character(reuse_existing_bootstrap_exports_if_present),
    "screening", "results_reuse_rule", "Reuse existing completed-run Stage 6 artifacts when compatible with the current Stage 1 backbone and revised Stage 6 specification; do not rerun heavy screening calculations unless compatibility fails.",
    "screening", "bootstrap_export_reuse_rule", "If refresh_bootstrap_exports = FALSE and existing bootstrap CSV/PDF exports are already present in export_path, reuse them and skip bootstrap reconstruction/regeneration.",
    "screening", "bootstrap_histogram_png_rule", "Each bootstrap histogram page is additionally saved as a standalone PNG that does not count toward the export file-number limit.",
    "outputs", "save_folder_rule", "Write all outputs into the single user-specified run_root_dir without creating extra subfolders.",
    "inputs", "data_path", normalize_existing_path(data_path),
    "inputs", "export_path", normalize_existing_path(export_path)
  )

  out <- out %>%
    filter(!metadata_name %in% patch_tbl$metadata_name) %>%
    bind_rows(patch_tbl) %>%
    distinct(metadata_group, metadata_name, .keep_all = TRUE)

  out %>% arrange(metadata_group, metadata_name)
}

# 🔴 Define: plot builders for Stage 6 summary figures ===============================
make_followup_plot <- function(followup_sufficiency_summary) {
  followup_metric_lookup <- tibble::tibble(
    metric = c(
      "tau_year",
      "largest_event_time_year",
      "plateau_length_year",
      "censor_rate_main",
      "km_tail_survival",
      "n_censored_after_last_event"
    ),
    metric_label = c(
      "Max follow-up (years)",
      "Largest event time (years)",
      "Plateau length (years)",
      "Main censoring rate",
      "KM tail survival",
      "Censored after last event"
    ),
    label_type = c("number", "number", "number", "percent", "percent", "count")
  )

  followup_plot_data <- bind_rows(lapply(followup_metric_lookup$metric, function(metric_name) {
    followup_sufficiency_summary %>%
      transmute(
        dataset_key = dataset_key,
        analysis_variant = analysis_variant,
        metric = metric_name,
        metric_value = .data[[metric_name]]
      )
  })) %>%
    left_join(followup_metric_lookup, by = "metric") %>%
    mutate(
      dataset_key = factor(dataset_key, levels = dataset_order_full),
      metric_label = factor(metric_label, levels = followup_metric_lookup$metric_label),
      metric_value_label = case_when(
        label_type == "percent" ~ scales::label_percent(accuracy = 0.1)(metric_value),
        label_type == "count" ~ ifelse(is.na(metric_value), NA_character_, format(round(metric_value), big.mark = ",", trim = TRUE)),
        TRUE ~ ifelse(is.na(metric_value), NA_character_, format(round(metric_value, 2), nsmall = 2, trim = TRUE))
      )
    )

  followup_plot <- ggplot2::ggplot(followup_plot_data, ggplot2::aes(x = metric_value, y = dataset_key)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = metric_value, yend = dataset_key), linewidth = 0.45, colour = "#BDBDBD", na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(colour = analysis_variant), size = 2.8, na.rm = TRUE) +
    ggplot2::geom_text(ggplot2::aes(label = metric_value_label), hjust = -0.10, size = 3.1, na.rm = TRUE) +
    ggplot2::facet_wrap(~metric_label, scales = "free_x", ncol = 2) +
    ggplot2::scale_colour_manual(values = c(single_cohort = "#1F78B4", site_free = "#33A02C", site_adjusted = "#E31A1C"), drop = FALSE) +
    ggplot2::labs(
      title = "Stage 6 follow-up sufficiency summary",
      subtitle = "Observed raw-tail summaries carried forward under the revised Stage 6 specification",
      x = NULL,
      y = NULL,
      colour = "Analysis variant"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::coord_cartesian(clip = "off")

  list(data = followup_plot_data, plot = followup_plot)
}

make_screening_method_plot <- function(screening_method_results) {
  screening_method_plot_data <- screening_method_results %>%
    mutate(
      dataset_key = factor(dataset_key, levels = dataset_order_full),
      method_name = factor(method_name, levels = method_order),
      method_label = factor(stringr::str_replace_all(as.character(method_name), "_", "\n"), levels = stringr::str_replace_all(method_order, "_", "\n")),
      method_flag = factor(method_flag, levels = c("supportive", "equivocal", "unsupportive")),
      plot_label = case_when(
        !is.na(p_value) ~ paste0("p=", formatC(p_value, format = "f", digits = 3)),
        !is.na(receus_ratio_hat) & !is.na(cure_fraction_hat) ~ paste0("r=", formatC(receus_ratio_hat, format = "f", digits = 2), "\ncf=", formatC(cure_fraction_hat, format = "f", digits = 2)),
        !is.na(statistic_value) ~ formatC(statistic_value, format = "f", digits = 2),
        TRUE ~ ""
      )
    )

  screening_method_plot <- ggplot2::ggplot(screening_method_plot_data, ggplot2::aes(x = method_label, y = dataset_key, fill = method_flag)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = plot_label), size = 2.6, lineheight = 0.95, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = c(supportive = "#1B9E77", equivocal = "#E6AB02", unsupportive = "#D95F02"), drop = FALSE, na.value = "#D9D9D9") +
    ggplot2::labs(
      title = "Stage 6 screening method flags",
      subtitle = "RECeUS-AIC primary gate, family-specific HSU rows, and HSU family-set row",
      x = NULL,
      y = NULL,
      fill = "Method flag"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 8)
    )

  list(data = screening_method_plot_data, plot = screening_method_plot)
}

make_receus_plot <- function(receus_candidate_fits) {
  receus_plot_data <- receus_candidate_fits %>%
    mutate(
      dataset_key = factor(dataset_key, levels = dataset_order_full),
      family = factor(family, levels = candidate_receus_families),
      selected_status = if_else(selected_by_aic, "AIC selected", "Candidate"),
      selected_status = factor(selected_status, levels = c("Candidate", "AIC selected"))
    )

  receus_plot <- ggplot2::ggplot(receus_plot_data, ggplot2::aes(x = receus_ratio_hat, y = cure_fraction_hat)) +
    ggplot2::geom_vline(xintercept = receus_r_threshold, linetype = "dashed", colour = "#636363") +
    ggplot2::geom_hline(yintercept = receus_pi_threshold, linetype = "dashed", colour = "#636363") +
    ggplot2::geom_point(ggplot2::aes(shape = family, colour = selected_status), size = 3.0, alpha = 0.9, na.rm = TRUE) +
    ggplot2::facet_wrap(~dataset_key, ncol = 2) +
    ggplot2::scale_colour_manual(values = c("Candidate" = "#6A3D9A", "AIC selected" = "#E31A1C"), drop = FALSE) +
    ggplot2::labs(
      title = "Stage 6 RECeUS candidate-family map",
      subtitle = "RECeUS-AIC remains the primary eligibility gate",
      x = "RECeUS ratio estimate",
      y = "Estimated cure fraction",
      colour = "Selection status",
      shape = "Latency family"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )

  list(data = receus_plot_data, plot = receus_plot)
}

make_carry_forward_plot <- function(carry_forward_flag_table) {
  carry_forward_base <- carry_forward_flag_table %>%
    mutate(
      presence_support_flag = coalesce(to_logical_column(cure_presence_support_flag), to_logical_column(presence_support_flag)),
      followup_not_contradicted_flag = coalesce(to_logical_column(followup_not_contradicted_flag), ifelse(is.na(followup_contradiction_flag), NA, !to_logical_column(followup_contradiction_flag))),
      dataset_key = factor(dataset_key, levels = dataset_order_full),
      cure_model_eligibility_flag = factor(cure_model_eligibility_flag, levels = c("supportive", "equivocal", "unsupportive")),
      presence_component = case_when(
        !is.na(presence_support_flag) & presence_support_flag ~ "supportive",
        !is.na(presence_support_flag) & !presence_support_flag ~ "unsupportive",
        TRUE ~ "equivocal"
      ),
      followup_component = case_when(
        primary_gate_flag == "supportive" & !is.na(followup_not_contradicted_flag) & followup_not_contradicted_flag ~ "supportive",
        primary_gate_flag == "unsupportive" ~ "unsupportive",
        primary_gate_flag == "supportive" & !is.na(followup_not_contradicted_flag) & !followup_not_contradicted_flag ~ "unsupportive",
        TRUE ~ "equivocal"
      )
    )

  carry_forward_plot_data <- bind_rows(
    carry_forward_base %>%
      transmute(
        dataset_key = dataset_key,
        component = "Presence evidence",
        decision_flag = presence_component,
        component_label = presence_evidence_label
      ),
    carry_forward_base %>%
      transmute(
        dataset_key = dataset_key,
        component = "Primary gate / Xie check",
        decision_flag = followup_component,
        component_label = followup_evidence_label
      ),
    carry_forward_base %>%
      transmute(
        dataset_key = dataset_key,
        component = "Eligibility flag",
        decision_flag = as.character(cure_model_eligibility_flag),
        component_label = as.character(cure_model_eligibility_flag)
      )
  ) %>%
    mutate(
      component = factor(component, levels = c("Presence evidence", "Primary gate / Xie check", "Eligibility flag")),
      decision_flag = factor(decision_flag, levels = c("supportive", "equivocal", "unsupportive")),
      component_label = as.character(component_label)
    )

  carry_forward_plot <- ggplot2::ggplot(carry_forward_plot_data, ggplot2::aes(x = component, y = dataset_key, fill = decision_flag)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = stringr::str_replace_all(component_label, "_", "\n")), size = 2.8, lineheight = 0.92, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = c(supportive = "#1B9E77", equivocal = "#7570B3", unsupportive = "#D95F02"), drop = FALSE, na.value = "#D9D9D9") +
    ggplot2::labs(
      title = "Stage 6 carry-forward eligibility map",
      subtitle = "Decision rule: RECeUS-AIC primary gate, presence-only modifier, Xie contradiction downgrade",
      x = NULL,
      y = NULL,
      fill = "Decision flag"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank()
    )

  list(data = carry_forward_plot_data, plot = carry_forward_plot)
}

# 🔴 Define: bootstrap reconstruction helpers ===============================
source_dataset_from_key <- function(dataset_key) {
  dplyr::case_when(
    dataset_key == "PNU" ~ "PNU",
    dataset_key == "SNU" ~ "SNU",
    dataset_key %in% c("merged__site_free", "merged__site_adjusted") ~ "merged",
    TRUE ~ NA_character_
  )
}

analysis_variant_from_key <- function(dataset_key) {
  dplyr::case_when(
    dataset_key %in% c("PNU", "SNU") ~ "single_cohort",
    dataset_key == "merged__site_free" ~ "site_free",
    dataset_key == "merged__site_adjusted" ~ "site_adjusted",
    TRUE ~ NA_character_
  )
}

hsu_method_name <- function(family) {
  switch(
    as.character(family),
    weibull = "hsu_supscore_weibull_approx",
    lognormal = "hsu_supscore_lognormal_approx",
    loglogistic = "hsu_supscore_loglogistic_approx",
    stop(sprintf("Unsupported HSU family for method naming: %s", family), call. = FALSE)
  )
}

stage6_shard_registry_columns <- c(
  "task_id", "task_type", "dataset_key", "source_dataset", "analysis_variant",
  "working_family", "shard_id", "n_shards", "rep_id_start", "rep_id_end", "n_reps",
  "status", "worker_id", "attempt_count", "created_at_utc", "started_at_utc",
  "heartbeat_at_utc", "finished_at_utc", "n_valid_stats", "shard_file", "observed_file", "note"
)

normalize_stage6_shard_registry <- function(df) {
  out <- tibble::as_tibble(df)
  for (nm in setdiff(stage6_shard_registry_columns, names(out))) {
    out[[nm]] <- NA_character_
  }
  out <- out[, stage6_shard_registry_columns, drop = FALSE] %>%
    mutate(across(everything(), as.character)) %>%
    mutate(
      task_type = normalize_blank_na_chr(task_type),
      dataset_key = normalize_blank_na_chr(dataset_key),
      working_family = normalize_blank_na_chr(working_family),
      status = tolower(trimws(as.character(status))),
      status = ifelse(is.na(status) | !nzchar(status), NA_character_, status),
      shard_file = normalize_blank_na_chr(shard_file),
      observed_file = normalize_blank_na_chr(observed_file)
    )
  out
}

parse_shard_file_metadata <- function(path) {
  nm <- basename(path)
  out <- list(
    task_type = NA_character_,
    dataset_key = NA_character_,
    working_family = NA_character_,
    shard_id = NA_integer_,
    n_shards = NA_integer_
  )

  shard_match <- regexec("__shard_([0-9]+)_of_([0-9]+)\\.rds$", nm)
  shard_parts <- regmatches(nm, shard_match)[[1]]
  if (length(shard_parts) == 3L) {
    out$shard_id <- as.integer(shard_parts[[2]])
    out$n_shards <- as.integer(shard_parts[[3]])
  }

  if (startsWith(nm, "xie__")) {
    core <- sub("^xie__", "", nm)
    core <- sub("__shard_[0-9]+_of_[0-9]+\\.rds$", "", core)
    out$task_type <- "xie"
    out$dataset_key <- core
    return(out)
  }

  if (startsWith(nm, "hsu__")) {
    core <- sub("^hsu__", "", nm)
    core <- sub("__shard_[0-9]+_of_[0-9]+\\.rds$", "", core)
    family_pattern <- paste(candidate_hsu_families, collapse = "|")
    if (grepl(paste0("__(", family_pattern, ")$"), core)) {
      out$task_type <- "hsu"
      out$working_family <- sub(paste0("^.*__(", family_pattern, ")$"), "\\1", core)
      out$dataset_key <- sub(paste0("__(", family_pattern, ")$"), "", core)
      return(out)
    }
  }

  out
}

read_one_existing_shard <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }

  finfo <- file.info(path)
  if (nrow(finfo) == 0L || is.na(finfo$size[[1]]) || finfo$size[[1]] <= 0L) {
    return(NULL)
  }

  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(obj) || !is.list(obj)) {
    return(NULL)
  }

  meta <- parse_shard_file_metadata(path)

  task_type <- normalize_blank_na_chr(extract_list_scalar(obj, "task_type", meta$task_type))
  dataset_key <- normalize_blank_na_chr(extract_list_scalar(obj, "dataset_key", meta$dataset_key))
  working_family <- normalize_blank_na_chr(extract_list_scalar(obj, "working_family", meta$working_family))
  rep_ids <- as_integer_or_na(obj$rep_ids)
  stats <- as_numeric_or_na(obj$statistics)

  if (length(rep_ids) == 0L || length(stats) == 0L || is.na(task_type) || is.na(dataset_key)) {
    return(NULL)
  }

  n_use <- min(length(rep_ids), length(stats))
  rep_ids <- rep_ids[seq_len(n_use)]
  stats <- stats[seq_len(n_use)]
  keep_idx <- which(!is.na(rep_ids))
  if (length(keep_idx) == 0L) {
    return(NULL)
  }

  tibble(
    task_type = task_type,
    dataset_key = dataset_key,
    working_family = ifelse(identical(task_type, "xie"), NA_character_, working_family),
    rep_id = rep_ids[keep_idx],
    statistic = stats[keep_idx],
    shard_file = normalize_existing_path(path),
    shard_id = meta$shard_id,
    n_shards = meta$n_shards,
    file_size_bytes = as.numeric(finfo$size[[1]]),
    file_mtime_utc = format(finfo$mtime[[1]], "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )
}

build_scanned_shard_registry <- function(shards_dir, observed_dir) {
  shard_files <- c(
    list.files(file.path(shards_dir, "xie"), pattern = "\\.rds$", full.names = TRUE),
    list.files(file.path(shards_dir, "hsu"), pattern = "\\.rds$", full.names = TRUE)
  )

  if (length(shard_files) == 0L) {
    return(tibble::tibble())
  }

  long_df <- bind_rows(lapply(shard_files, read_one_existing_shard))
  if (nrow(long_df) == 0L) {
    return(tibble::tibble())
  }

  long_df %>%
    group_by(task_type, dataset_key, working_family, shard_file, shard_id, n_shards, file_mtime_utc) %>%
    summarise(
      rep_id_start = min(rep_id, na.rm = TRUE),
      rep_id_end = max(rep_id, na.rm = TRUE),
      n_reps = n(),
      n_valid_stats = sum(is.finite(statistic)),
      .groups = "drop"
    ) %>%
    mutate(
      source_dataset = source_dataset_from_key(dataset_key),
      analysis_variant = analysis_variant_from_key(dataset_key),
      task_id = ifelse(
        task_type == "xie",
        sprintf("xie__%s__%03d", dataset_key, as_integer_or_na(shard_id)),
        sprintf("hsu__%s__%s__%03d", dataset_key, working_family, as_integer_or_na(shard_id))
      ),
      status = "completed",
      worker_id = NA_character_,
      attempt_count = NA_character_,
      created_at_utc = NA_character_,
      started_at_utc = NA_character_,
      heartbeat_at_utc = file_mtime_utc,
      finished_at_utc = file_mtime_utc,
      observed_file = ifelse(
        task_type == "xie",
        file.path(observed_dir, paste0("xie__", dataset_key, "__observed.rds")),
        file.path(observed_dir, paste0("hsu__", dataset_key, "__", working_family, "__observed.rds"))
      ),
      note = "Reconstructed from existing shard RDS files."
    ) %>%
    normalize_stage6_shard_registry()
}

choose_best_stage6_shard_registry <- function(existing_registry, scanned_registry) {
  existing_norm <- normalize_stage6_shard_registry(existing_registry)
  scanned_norm <- normalize_stage6_shard_registry(scanned_registry)

  existing_completed <- sum(existing_norm$status == "completed", na.rm = TRUE)
  scanned_completed <- nrow(scanned_norm)

  if (scanned_completed > existing_completed) {
    return(scanned_norm)
  }

  if (nrow(existing_norm) == 0L && scanned_completed > 0L) {
    return(scanned_norm)
  }

  existing_norm
}

collect_completed_bootstrap_stats_long_from_shards <- function(shards_dir, task_type, dataset_key, working_family = NA_character_) {
  task_type <- normalize_blank_na_chr(task_type)
  dataset_key <- normalize_blank_na_chr(dataset_key)
  working_family_filter <- normalize_blank_na_chr(working_family)

  subdir <- file.path(shards_dir, task_type)
  if (!dir.exists(subdir)) {
    return(tibble(rep_id = integer(), statistic = numeric()))
  }

  shard_files <- list.files(subdir, pattern = "\\.rds$", full.names = TRUE)
  if (length(shard_files) == 0L) {
    return(tibble(rep_id = integer(), statistic = numeric()))
  }

  if (identical(task_type, "xie")) {
    shard_files <- shard_files[startsWith(basename(shard_files), paste0("xie__", dataset_key, "__shard_"))]
  } else {
    if (is.na(working_family_filter)) {
      return(tibble(rep_id = integer(), statistic = numeric()))
    }
    shard_files <- shard_files[startsWith(basename(shard_files), paste0("hsu__", dataset_key, "__", working_family_filter, "__shard_"))]
  }

  if (length(shard_files) == 0L) {
    return(tibble(rep_id = integer(), statistic = numeric()))
  }

  long_df <- bind_rows(lapply(shard_files, read_one_existing_shard))
  if (nrow(long_df) == 0L) {
    return(tibble(rep_id = integer(), statistic = numeric()))
  }

  long_df <- long_df %>%
    filter(task_type == !!task_type, dataset_key == !!dataset_key)

  if (!is.na(working_family_filter)) {
    long_df <- long_df %>% filter(working_family == !!working_family_filter)
  }

  if (nrow(long_df) == 0L) {
    return(tibble(rep_id = integer(), statistic = numeric()))
  }

  long_df %>%
    arrange(rep_id, file_mtime_utc, shard_file) %>%
    group_by(rep_id) %>%
    slice_tail(n = 1L) %>%
    ungroup() %>%
    arrange(rep_id) %>%
    select(rep_id, statistic)
}

collect_completed_bootstrap_stats_long <- function(registry_df, task_type, dataset_key, working_family = NA_character_) {
  subset_df <- registry_df %>%
    dplyr::filter(task_type == !!task_type, dataset_key == !!dataset_key, status == "completed")

  if (!is.na(working_family)) {
    subset_df <- subset_df %>% dplyr::filter(working_family == !!working_family)
  }

  if (nrow(subset_df) == 0L) {
    return(tibble(rep_id = integer(), statistic = numeric()))
  }

  shard_rows <- list()
  for (ii in seq_len(nrow(subset_df))) {
    shard_file <- subset_df$shard_file[[ii]]
    if (!file.exists(shard_file)) {
      next
    }

    obj <- tryCatch(readRDS(shard_file), error = function(e) NULL)
    if (is.null(obj)) {
      next
    }

    rep_ids <- as_integer_or_na(obj$rep_ids)
    stats <- as_numeric_or_na(obj$statistics)
    if (length(rep_ids) == 0L || length(stats) == 0L) {
      next
    }

    n_use <- min(length(rep_ids), length(stats))
    rep_ids <- rep_ids[seq_len(n_use)]
    stats <- stats[seq_len(n_use)]
    keep_idx <- which(!is.na(rep_ids))
    if (length(keep_idx) == 0L) {
      next
    }

    shard_rows[[length(shard_rows) + 1L]] <- tibble(
      rep_id = rep_ids[keep_idx],
      statistic = stats[keep_idx]
    )
  }

  if (length(shard_rows) == 0L) {
    return(tibble(rep_id = integer(), statistic = numeric()))
  }

  bind_rows(shard_rows) %>%
    group_by(rep_id) %>%
    slice_tail(n = 1L) %>%
    ungroup() %>%
    arrange(rep_id)
}

finalize_xie_from_bootstrap_stats <- function(observed_object, bootstrap_stats, alpha_screening) {
  observed <- observed_object$observed
  bootstrap_stats <- as.numeric(bootstrap_stats)

  if (is.na(observed$T_n)) {
    return(list(observed = observed, p_centered = NA_real_, p_raw = NA_real_, critical_centered = NA_real_))
  }

  if (isTRUE(observed$degenerate) && observed$T_n == 0) {
    return(list(observed = observed, p_centered = 1, p_raw = 1, critical_centered = 0))
  }

  centered_stats <- bootstrap_stats - observed$T_n
  n_valid_centered <- sum(is.finite(centered_stats))
  n_valid_raw <- sum(is.finite(bootstrap_stats))

  critical_centered <- if (n_valid_centered > 0L) {
    as.numeric(stats::quantile(centered_stats, probs = 1 - alpha_screening, na.rm = TRUE, type = 8, names = FALSE))
  } else {
    NA_real_
  }

  p_centered <- if (n_valid_centered > 0L) {
    (sum(centered_stats >= observed$T_n, na.rm = TRUE) + 1) / (n_valid_centered + 1)
  } else {
    NA_real_
  }

  p_raw <- if (n_valid_raw > 0L) {
    (sum(bootstrap_stats >= observed$T_n, na.rm = TRUE) + 1) / (n_valid_raw + 1)
  } else {
    NA_real_
  }

  list(observed = observed, p_centered = p_centered, p_raw = p_raw, critical_centered = critical_centered)
}

finalize_hsu_family_from_bootstrap_stats <- function(observed_object, bootstrap_stats, alpha_screening) {
  observed <- observed_object$result
  bootstrap_stats <- as.numeric(bootstrap_stats)
  working_family <- observed_object$working_family

  if (sum(is.finite(bootstrap_stats)) == 0L || is.null(observed$null_fit) || !isTRUE(observed$null_fit$converged) || is.na(observed$statistic)) {
    return(list(
      statistic = observed$statistic,
      p_value = NA_real_,
      critical_value = NA_real_,
      working_family = working_family,
      selected_formula = NA_character_,
      latency_formula = NA_character_,
      note = paste0("HSU bootstrap export could not finalize for family `", working_family, "`.")
    ))
  }

  n_valid <- sum(is.finite(bootstrap_stats))
  p_value <- if (n_valid > 0L) {
    (sum(bootstrap_stats >= observed$statistic, na.rm = TRUE) + 1) / (n_valid + 1)
  } else {
    NA_real_
  }

  critical_value <- if (n_valid > 0L) {
    as.numeric(stats::quantile(bootstrap_stats, probs = 1 - alpha_screening, na.rm = TRUE, type = 8, names = FALSE))
  } else {
    NA_real_
  }

  list(
    statistic = observed$statistic,
    p_value = p_value,
    critical_value = critical_value,
    working_family = working_family,
    selected_formula = extract_list_scalar(observed, "selected_formula", NA_character_),
    latency_formula = extract_list_scalar(observed, "latency_formula", NA_character_),
    note = paste0(
      "Self-contained pragmatic approximation to the Hsu sup-score test using a ",
      working_family,
      " AFT working latency model and sharded parametric bootstrap under the no-cure null."
    )
  )
}

build_bootstrap_distribution_table <- function(shards_dir, observed_dir, alpha_screening, candidate_hsu_families) {
  xie_dataset_keys <- c("PNU", "SNU", "merged__site_free")
  out_rows <- list()

  for (dataset_key in xie_dataset_keys) {
    observed_file <- file.path(observed_dir, paste0("xie__", dataset_key, "__observed.rds"))
    if (!file.exists(observed_file)) {
      next
    }

    observed_object <- readRDS(observed_file)
    shard_tbl <- collect_completed_bootstrap_stats_long_from_shards(shards_dir = shards_dir, task_type = "xie", dataset_key = dataset_key)
    if (nrow(shard_tbl) == 0L) {
      next
    }

    finalized <- finalize_xie_from_bootstrap_stats(observed_object, shard_tbl$statistic, alpha_screening = alpha_screening)
    out_rows[[length(out_rows) + 1L]] <- tibble(
      record_type = "distribution",
      task_type = "xie",
      dataset_key = dataset_key,
      source_dataset = source_dataset_from_key(dataset_key),
      analysis_variant = analysis_variant_from_key(dataset_key),
      working_family = NA_character_,
      method_name = "xie_sufficient_followup",
      statistic_name = "T_n",
      rep_id = shard_tbl$rep_id,
      bootstrap_statistic = as.numeric(shard_tbl$statistic),
      centered_statistic = as.numeric(shard_tbl$statistic) - finalized$observed$T_n,
      observed_statistic = finalized$observed$T_n,
      p_value = finalized$p_centered,
      critical_value = finalized$critical_centered,
      critical_scale = "centered",
      plot_metric_name = "centered_bootstrap_statistic",
      plot_statistic_value = as.numeric(shard_tbl$statistic) - finalized$observed$T_n,
      plot_observed_line = finalized$observed$T_n,
      plot_critical_line = finalized$critical_centered,
      valid_statistic = is.finite(as.numeric(shard_tbl$statistic)),
      tau_year = as.numeric(observed_object$tau_year),
      selected_formula_rhs = NA_character_,
      latency_formula_rhs = NA_character_,
      note = "Centered Xie bootstrap distribution: histogram is drawn for T*_n - T_n and compared against the observed T_n and centered critical value."
    )
  }

  for (dataset_key in dataset_order_full) {
    for (family_name in candidate_hsu_families) {
      observed_file <- file.path(observed_dir, paste0("hsu__", dataset_key, "__", family_name, "__observed.rds"))
      if (!file.exists(observed_file)) {
        next
      }

      observed_object <- readRDS(observed_file)
      shard_tbl <- collect_completed_bootstrap_stats_long_from_shards(
        shards_dir = shards_dir,
        task_type = "hsu",
        dataset_key = dataset_key,
        working_family = family_name
      )
      if (nrow(shard_tbl) == 0L) {
        next
      }

      finalized <- finalize_hsu_family_from_bootstrap_stats(observed_object, shard_tbl$statistic, alpha_screening = alpha_screening)
      out_rows[[length(out_rows) + 1L]] <- tibble(
        record_type = "distribution",
        task_type = "hsu",
        dataset_key = dataset_key,
        source_dataset = source_dataset_from_key(dataset_key),
        analysis_variant = analysis_variant_from_key(dataset_key),
        working_family = family_name,
        method_name = hsu_method_name(family_name),
        statistic_name = "sup_lr_boot",
        rep_id = shard_tbl$rep_id,
        bootstrap_statistic = as.numeric(shard_tbl$statistic),
        centered_statistic = as.numeric(shard_tbl$statistic) - finalized$statistic,
        observed_statistic = finalized$statistic,
        p_value = finalized$p_value,
        critical_value = finalized$critical_value,
        critical_scale = "raw",
        plot_metric_name = "bootstrap_statistic",
        plot_statistic_value = as.numeric(shard_tbl$statistic),
        plot_observed_line = finalized$statistic,
        plot_critical_line = finalized$critical_value,
        valid_statistic = is.finite(as.numeric(shard_tbl$statistic)),
        tau_year = as.numeric(observed_object$tau_year),
        selected_formula_rhs = strip_formula_rhs(finalized$selected_formula),
        latency_formula_rhs = strip_formula_rhs(finalized$latency_formula),
        note = extract_list_scalar(finalized, "note", NA_character_)
      )
    }
  }

  dist_df <- bind_rows(out_rows)
  if (nrow(dist_df) == 0L) {
    return(tibble())
  }

  summary_df <- dist_df %>%
    group_by(
      task_type, dataset_key, source_dataset, analysis_variant, working_family,
      method_name, statistic_name, observed_statistic, p_value, critical_value,
      critical_scale, plot_metric_name, plot_observed_line, plot_critical_line,
      tau_year, selected_formula_rhs, latency_formula_rhs, note
    ) %>%
    summarise(
      n_bootstrap = n(),
      n_valid = sum(valid_statistic),
      valid_fraction = safe_divide(n_valid, n_bootstrap),
      bootstrap_mean = safe_mean(bootstrap_statistic),
      bootstrap_sd = safe_sd(bootstrap_statistic),
      bootstrap_min = safe_min(bootstrap_statistic),
      bootstrap_q025 = safe_quantile(bootstrap_statistic, probs = 0.025),
      bootstrap_q500 = safe_quantile(bootstrap_statistic, probs = 0.500),
      bootstrap_q975 = safe_quantile(bootstrap_statistic, probs = 0.975),
      centered_mean = safe_mean(centered_statistic),
      centered_sd = safe_sd(centered_statistic),
      centered_q025 = safe_quantile(centered_statistic, probs = 0.025),
      centered_q500 = safe_quantile(centered_statistic, probs = 0.500),
      centered_q975 = safe_quantile(centered_statistic, probs = 0.975),
      plot_min = safe_min(plot_statistic_value),
      plot_max = safe_max(plot_statistic_value),
      .groups = "drop"
    )

  dist_df %>%
    left_join(
      summary_df,
      by = c(
        "task_type", "dataset_key", "source_dataset", "analysis_variant", "working_family",
        "method_name", "statistic_name", "observed_statistic", "p_value", "critical_value",
        "critical_scale", "plot_metric_name", "plot_observed_line", "plot_critical_line",
        "tau_year", "selected_formula_rhs", "latency_formula_rhs", "note"
      )
    ) %>%
    arrange(
      factor(dataset_key, levels = dataset_order_full),
      task_type,
      working_family,
      rep_id
    )
}

build_bootstrap_plot_index <- function(bootstrap_df) {
  family_order <- c("weibull", "lognormal", "loglogistic")
  if (nrow(bootstrap_df) == 0L) {
    return(tibble())
  }
  bootstrap_df %>%
    mutate(
      plot_group_id = case_when(
        task_type == "xie" ~ paste0("xie__", dataset_key),
        TRUE ~ paste0("hsu__", dataset_key, "__", working_family)
      ),
      png_file_name = case_when(
        task_type == "xie" ~ paste0("stage6_bootstrap_histogram__", plot_group_id, ".png"),
        TRUE ~ paste0("stage6_bootstrap_histogram__", plot_group_id, ".png")
      ),
      png_file_path = file.path(export_path, png_file_name)
    ) %>%
    distinct(
      plot_group_id, task_type, dataset_key, source_dataset, analysis_variant,
      working_family, method_name, statistic_name, plot_metric_name,
      observed_statistic, p_value, critical_value, plot_observed_line,
      plot_critical_line, n_bootstrap, n_valid, valid_fraction,
      selected_formula_rhs, latency_formula_rhs, tau_year, png_file_name, png_file_path
    ) %>%
    arrange(
      factor(dataset_key, levels = dataset_order_full),
      task_type,
      factor(working_family, levels = family_order)
    )
}

generate_bootstrap_histograms <- function(bootstrap_df, plot_index_df, pdf_path) {
  if (nrow(plot_index_df) == 0L) {
    safe_save_pdf_atomic(pdf_path, width = bootstrap_plot_width, height = bootstrap_plot_height, plot_fun = function() {
      p <- ggplot() +
        annotate("text", x = 0, y = 0, label = "No completed bootstrap shard results were found.", size = 5) +
        theme_void()
      print(p)
    })
    return(invisible(TRUE))
  }

  build_one_plot <- function(current_group_id) {
    current_index <- plot_index_df %>% filter(plot_group_id == current_group_id)
    current_df <- bootstrap_df %>%
      mutate(plot_group_id = case_when(
        task_type == "xie" ~ paste0("xie__", dataset_key),
        TRUE ~ paste0("hsu__", dataset_key, "__", working_family)
      )) %>%
      filter(plot_group_id == current_group_id)

    title_text <- if (current_index$task_type[[1]] == "xie") {
      paste0("Stage 6 bootstrap distribution: Xie | ", current_index$dataset_key[[1]])
    } else {
      paste0(
        "Stage 6 bootstrap distribution: HSU-",
        current_index$working_family[[1]],
        " | ",
        current_index$dataset_key[[1]]
      )
    }

    subtitle_text <- paste0(
      "Observed = ", formatC(current_index$plot_observed_line[[1]], format = "f", digits = 4),
      "; Critical = ", formatC(current_index$plot_critical_line[[1]], format = "f", digits = 4),
      "; p = ", formatC(current_index$p_value[[1]], format = "f", digits = 4),
      "; valid = ", current_index$n_valid[[1]], "/", current_index$n_bootstrap[[1]]
    )

    ggplot(current_df, aes(x = plot_statistic_value)) +
      geom_histogram(bins = bootstrap_hist_bins, fill = "#4C78A8", colour = "white", linewidth = 0.3, na.rm = TRUE) +
      geom_density(adjust = bootstrap_density_adjust, linewidth = 0.6, colour = "#1F1F1F", na.rm = TRUE) +
      geom_vline(aes(xintercept = plot_observed_line), colour = "#D95F02", linewidth = 0.7) +
      geom_vline(aes(xintercept = plot_critical_line), colour = "#7570B3", linewidth = 0.7, linetype = "dashed") +
      labs(
        title = title_text,
        subtitle = subtitle_text,
        x = current_index$plot_metric_name[[1]],
        y = "Count"
      ) +
      theme_minimal(base_size = 11) +
      theme(panel.grid.minor = element_blank())
  }

  safe_save_pdf_atomic(pdf_path, width = bootstrap_plot_width, height = bootstrap_plot_height, plot_fun = function() {
    for (ii in seq_len(nrow(plot_index_df))) {
      print(build_one_plot(plot_index_df$plot_group_id[[ii]]))
    }
  })

  for (ii in seq_len(nrow(plot_index_df))) {
    p <- build_one_plot(plot_index_df$plot_group_id[[ii]])
    safe_save_plot(
      p,
      file.path(export_path, plot_index_df$png_file_name[[ii]]),
      width = bootstrap_plot_width,
      height = bootstrap_plot_height,
      dpi = plot_dpi
    )
  }

  invisible(TRUE)
}

# 🔴 Load: completed Stage 6 run and Stage 1 backbone ===============================
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)
export_path <- normalize_existing_path(export_path)
data_path <- normalize_existing_path(data_path)
run_root_dir <- resolve_stage6_run_root(run_root_dir = run_root_dir, export_path = export_path)
run_root_dir <- normalize_existing_path(run_root_dir)

message("Using existing Stage 6 run root: ", run_root_dir)

stage1_inputs <- read_stage1_inputs(data_path)
stage6_objects <- read_stage6_bundle(run_root_dir)
stage6_outputs <- extract_stage6_outputs(stage6_objects)

followup_sufficiency_summary <- patch_followup_sufficiency_summary(stage6_outputs$followup_sufficiency_summary)
check_stage6_reuse_compatibility(stage1_inputs, list(followup_sufficiency_summary = followup_sufficiency_summary))

variant_registry <- canonicalize_variant_registry(stage6_outputs$variant_registry, stage1_inputs)
screening_method_results <- patch_screening_method_results(stage6_outputs$screening_method_results)
receus_candidate_fits <- patch_receus_candidate_fits(stage6_outputs$receus_candidate_fits)
carry_forward_flag_table <- patch_carry_forward_flag_table(
  stage6_outputs$carry_forward_flag_table,
  stage1_inputs,
  variant_registry,
  screening_method_results,
  followup_sufficiency_summary
)
screening_summary <- patch_screening_summary(
  stage6_outputs$screening_summary,
  carry_forward_flag_table,
  screening_method_results,
  followup_sufficiency_summary,
  variant_registry
)
stage6_metadata_registry <- patch_stage6_metadata_registry(stage6_outputs$metadata_registry, stage1_inputs, reuse_completed_run = reuse_existing_completed_run_only)
scanned_stage6_shard_registry <- build_scanned_shard_registry(
  shards_dir = stage6_objects$shards_dir,
  observed_dir = stage6_objects$observed_dir
)
stage6_shard_registry <- choose_best_stage6_shard_registry(
  existing_registry = stage6_objects$shard_registry,
  scanned_registry = scanned_stage6_shard_registry
)
fit_objects_final <- stage6_objects$fit_objects

# 🔴 Transform: visual outputs and bootstrap exports ===============================
followup_plot_bundle <- make_followup_plot(followup_sufficiency_summary)
screening_method_plot_bundle <- make_screening_method_plot(screening_method_results)
receus_plot_bundle <- make_receus_plot(receus_candidate_fits)
carry_forward_plot_bundle <- make_carry_forward_plot(carry_forward_flag_table)

visual_plot_data <- bind_rows(
  followup_plot_bundle$data %>% mutate(plot_name = "followup_sufficiency_summary"),
  screening_method_plot_bundle$data %>% mutate(plot_name = "screening_method_results"),
  receus_plot_bundle$data %>% mutate(plot_name = "receus_candidate_fits"),
  carry_forward_plot_bundle$data %>% mutate(plot_name = "carry_forward_flag_table")
)

existing_bootstrap_results_file <- file.path(export_path, "stage6_bootstrap_results.csv")
existing_bootstrap_plot_index_file <- file.path(export_path, "stage6_bootstrap_plot_index.csv")
existing_bootstrap_histograms_pdf_file <- file.path(export_path, "stage6_bootstrap_histograms.pdf")

reuse_existing_bootstrap_tables <- FALSE
reuse_existing_bootstrap_histograms <- FALSE

if (isTRUE(reuse_existing_bootstrap_exports_if_present) && !isTRUE(refresh_bootstrap_exports)) {
  bootstrap_results_existing <- safe_read_csv_if_exists(existing_bootstrap_results_file)
  bootstrap_plot_index_existing <- safe_read_csv_if_exists(existing_bootstrap_plot_index_file)

  if (!is.null(bootstrap_results_existing) && !is.null(bootstrap_plot_index_existing)) {
    bootstrap_results <- tibble::as_tibble(bootstrap_results_existing)
    bootstrap_plot_index <- tibble::as_tibble(bootstrap_plot_index_existing)
    reuse_existing_bootstrap_tables <- TRUE
    message("Reusing existing bootstrap CSV exports from export_path: ", export_path)
  }
}

if (!isTRUE(reuse_existing_bootstrap_tables)) {
  bootstrap_results <- build_bootstrap_distribution_table(
    shards_dir = stage6_objects$shards_dir,
    observed_dir = stage6_objects$observed_dir,
    alpha_screening = alpha_screening,
    candidate_hsu_families = candidate_hsu_families
  )
  bootstrap_plot_index <- build_bootstrap_plot_index(bootstrap_results)
}

reuse_existing_bootstrap_histograms <- isTRUE(reuse_existing_bootstrap_tables) &&
  isTRUE(reuse_existing_bootstrap_exports_if_present) &&
  !isTRUE(refresh_bootstrap_exports) &&
  file.exists(existing_bootstrap_histograms_pdf_file)

if (nrow(bootstrap_results) == 0L) {
  warning(
    paste0(
      "No completed bootstrap shard results were recovered from the shard files under: ",
      stage6_objects$shards_dir,
      ". Check that `export_path` or `run_root_dir` points to the completed Stage 6 run folder."
    ),
    call. = FALSE
  )
}

# 🔴 Export: revised Stage 6 artifacts ===============================
file_variants <- list(
  stage6_variant_registry = file.path(export_path, "stage6_variant_registry.csv"),
  stage6_metadata_registry = file.path(export_path, "stage6_metadata_registry.csv"),
  stage6_followup_sufficiency_summary = file.path(export_path, "stage6_followup_sufficiency_summary.csv"),
  stage6_screening_method_results = file.path(export_path, "stage6_screening_method_results.csv"),
  stage6_receus_candidate_fits = file.path(export_path, "stage6_receus_candidate_fits.csv"),
  stage6_screening_summary = file.path(export_path, "stage6_screening_summary.csv"),
  stage6_carry_forward_flag_table = file.path(export_path, "stage6_carry_forward_flag_table.csv"),
  stage6_shard_registry = file.path(export_path, "stage6_shard_registry.csv"),
  stage6_visual_plot_data = file.path(export_path, "stage6_visual_plot_data.csv"),
  stage6_bootstrap_results = file.path(export_path, "stage6_bootstrap_results.csv"),
  stage6_bootstrap_plot_index = file.path(export_path, "stage6_bootstrap_plot_index.csv"),
  stage6_screening_bundle = file.path(export_path, "stage6_screening_bundle.rds"),
  stage6_fitted_objects = file.path(export_path, "stage6_fitted_objects.rds"),
  stage6_visual_summary_pdf = file.path(export_path, "stage6_visual_summary.pdf"),
  stage6_bootstrap_histograms_pdf = file.path(export_path, "stage6_bootstrap_histograms.pdf"),
  stage6_export_manifest = file.path(export_path, "stage6_export_manifest.csv")
)

safe_write_csv_atomic(variant_registry, file_variants$stage6_variant_registry)
safe_write_csv_atomic(stage6_metadata_registry, file_variants$stage6_metadata_registry)
safe_write_csv_atomic(followup_sufficiency_summary, file_variants$stage6_followup_sufficiency_summary)
safe_write_csv_atomic(screening_method_results, file_variants$stage6_screening_method_results)
safe_write_csv_atomic(receus_candidate_fits, file_variants$stage6_receus_candidate_fits)
safe_write_csv_atomic(screening_summary, file_variants$stage6_screening_summary)
safe_write_csv_atomic(carry_forward_flag_table, file_variants$stage6_carry_forward_flag_table)
safe_write_csv_atomic(stage6_shard_registry, file_variants$stage6_shard_registry)
safe_write_csv_atomic(visual_plot_data, file_variants$stage6_visual_plot_data)

if (!isTRUE(reuse_existing_bootstrap_tables) || !file.exists(file_variants$stage6_bootstrap_results)) {
  safe_write_csv_atomic(bootstrap_results, file_variants$stage6_bootstrap_results)
}
if (!isTRUE(reuse_existing_bootstrap_tables) || !file.exists(file_variants$stage6_bootstrap_plot_index)) {
  safe_write_csv_atomic(bootstrap_plot_index, file_variants$stage6_bootstrap_plot_index)
}

screening_bundle_revised <- list(
  stage = "Stage 6",
  created_at = as.character(Sys.time()),
  session_info = utils::sessionInfo(),
  config = list(
    run_root_dir = run_root_dir,
    export_path = export_path,
    stage6_specification_file = shared_master_spec_file,
    stage6_specification_section = "6.6",
    alpha_screening = alpha_screening,
    candidate_receus_families = candidate_receus_families,
    candidate_hsu_families = candidate_hsu_families,
    hsu_family_set_name = hsu_family_set_name,
    decision_rule_version = decision_rule_version,
    screening_context = screening_context,
    reuse_existing_completed_run_only = reuse_existing_completed_run_only,
    refresh_bootstrap_exports = refresh_bootstrap_exports,
    reuse_existing_bootstrap_exports_if_present = reuse_existing_bootstrap_exports_if_present
  ),
  source_lookup = stage1_inputs$backbone_bundle$source_lookup,
  registries = list(
    variant_registry = variant_registry,
    stage6_metadata_registry = stage6_metadata_registry,
    shard_registry = stage6_shard_registry
  ),
  outputs = list(
    followup_sufficiency_summary = followup_sufficiency_summary,
    screening_method_results = screening_method_results,
    receus_candidate_fits = receus_candidate_fits,
    screening_summary = screening_summary,
    carry_forward_flag_table = carry_forward_flag_table,
    visual_plot_data = visual_plot_data,
    bootstrap_results = bootstrap_results,
    bootstrap_plot_index = bootstrap_plot_index
  )
)
safe_save_rds_atomic(screening_bundle_revised, file_variants$stage6_screening_bundle)

if (isTRUE(save_fit_objects) && !is.null(fit_objects_final)) {
  safe_save_rds_atomic(fit_objects_final, file_variants$stage6_fitted_objects)
}

safe_save_plot(
  followup_plot_bundle$plot,
  file.path(export_path, "stage6_followup_sufficiency_summary.png"),
  width = core_plot_width,
  height = core_plot_height_tall,
  dpi = plot_dpi
)
safe_save_plot(
  screening_method_plot_bundle$plot,
  file.path(export_path, "stage6_screening_method_results.png"),
  width = core_plot_width,
  height = core_plot_height,
  dpi = plot_dpi
)
safe_save_plot(
  receus_plot_bundle$plot,
  file.path(export_path, "stage6_receus_candidate_fits.png"),
  width = core_plot_width,
  height = core_plot_height,
  dpi = plot_dpi
)
safe_save_plot(
  carry_forward_plot_bundle$plot,
  file.path(export_path, "stage6_carry_forward_flag_table.png"),
  width = core_plot_width,
  height = core_plot_height,
  dpi = plot_dpi
)

if (isTRUE(rebuild_visual_summary_pdf)) {
  safe_save_pdf_atomic(
    file_variants$stage6_visual_summary_pdf,
    width = core_plot_width,
    height = core_plot_height,
    plot_fun = function() {
      print(followup_plot_bundle$plot)
      print(screening_method_plot_bundle$plot)
      print(receus_plot_bundle$plot)
      print(carry_forward_plot_bundle$plot)
    }
  )
}

if (!isTRUE(reuse_existing_bootstrap_histograms) || !file.exists(file_variants$stage6_bootstrap_histograms_pdf)) {
  generate_bootstrap_histograms(bootstrap_results, bootstrap_plot_index, file_variants$stage6_bootstrap_histograms_pdf)
} else {
  message("Reusing existing bootstrap histogram PDF and standalone PNG files from export_path: ", export_path)
}

export_manifest <- tibble::tibble(
  file_name = c(
    "stage6_variant_registry.csv",
    "stage6_metadata_registry.csv",
    "stage6_followup_sufficiency_summary.csv",
    "stage6_screening_method_results.csv",
    "stage6_receus_candidate_fits.csv",
    "stage6_screening_summary.csv",
    "stage6_carry_forward_flag_table.csv",
    "stage6_shard_registry.csv",
    "stage6_visual_plot_data.csv",
    "stage6_bootstrap_results.csv",
    "stage6_bootstrap_plot_index.csv",
    "stage6_screening_bundle.rds",
    if (isTRUE(save_fit_objects) && !is.null(fit_objects_final)) "stage6_fitted_objects.rds" else NA_character_,
    "stage6_visual_summary.pdf",
    "stage6_bootstrap_histograms.pdf"
  ),
  object_name = c(
    "variant_registry",
    "stage6_metadata_registry",
    "followup_sufficiency_summary",
    "screening_method_results",
    "receus_candidate_fits",
    "screening_summary",
    "carry_forward_flag_table",
    "stage6_shard_registry",
    "visual_plot_data",
    "bootstrap_results",
    "bootstrap_plot_index",
    "screening_bundle",
    if (isTRUE(save_fit_objects) && !is.null(fit_objects_final)) "fit_objects" else NA_character_,
    "visual_summary_pdf",
    "bootstrap_histogram_pdf"
  ),
  description = c(
    "Stage 6 structural variant registry keyed by dataset_key.",
    "Stage 6 metadata and carry-forward rules updated to the revised integrated specification.",
    "Dataset-level follow-up sufficiency summary including extreme-tail quantities.",
    "Long-format Stage 6 screening method results including family-specific HSU rows and HSU family-set row.",
    "RECeUS and RECeUS-AIC candidate-family fit summary table.",
    "Wide Stage 6 screening summary by dataset_key.",
    "Carry-forward eligibility flag table with canonical Stage 7/8 fields patched under the revised specification.",
    "Shard registry tracking completed shard tasks.",
    "Unified source-of-truth table used to draw the main Stage 6 summary plots.",
    "Long-format bootstrap source-of-truth table reconstructed from completed Stage 6 shard files or reused from an existing export when present.",
    "Index of bootstrap histogram groups and their standalone PNG filenames, reused when existing export files are present and refresh_bootstrap_exports = FALSE.",
    "Reusable Stage 6 bundle with patched tables and bootstrap source-of-truth objects.",
    if (isTRUE(save_fit_objects) && !is.null(fit_objects_final)) "Observed objects and reduced bootstrap/fitted-object bundle reused from the completed Stage 6 run." else NA_character_,
    "Multi-page PDF containing the four main Stage 6 summary figures.",
    "Multi-page PDF of bootstrap histograms, with standalone PNGs exported for each page unless existing histogram exports are reused."
  ),
  file_path = c(
    file_variants$stage6_variant_registry,
    file_variants$stage6_metadata_registry,
    file_variants$stage6_followup_sufficiency_summary,
    file_variants$stage6_screening_method_results,
    file_variants$stage6_receus_candidate_fits,
    file_variants$stage6_screening_summary,
    file_variants$stage6_carry_forward_flag_table,
    file_variants$stage6_shard_registry,
    file_variants$stage6_visual_plot_data,
    file_variants$stage6_bootstrap_results,
    file_variants$stage6_bootstrap_plot_index,
    file_variants$stage6_screening_bundle,
    if (isTRUE(save_fit_objects) && !is.null(fit_objects_final)) file_variants$stage6_fitted_objects else NA_character_,
    file_variants$stage6_visual_summary_pdf,
    file_variants$stage6_bootstrap_histograms_pdf
  )
) %>%
  filter(!is.na(file_name))

safe_write_csv_atomic(export_manifest, file_variants$stage6_export_manifest)

message("Stage 6 revised exports refreshed from the compatible completed run:")
message("  - ", file_variants$stage6_variant_registry)
message("  - ", file_variants$stage6_screening_method_results)
message("  - ", file_variants$stage6_carry_forward_flag_table)
message("  - ", file_variants$stage6_bootstrap_results, if (isTRUE(reuse_existing_bootstrap_tables)) " (reused existing CSV)" else " (refreshed from shard files)")
message("  - ", file_variants$stage6_bootstrap_histograms_pdf, if (isTRUE(reuse_existing_bootstrap_histograms)) " (reused existing PDF/PNG set)" else " (refreshed)")
message("  - individual PNG files for every main figure and bootstrap histogram page were also written to export_path unless existing bootstrap histogram exports were reused")
