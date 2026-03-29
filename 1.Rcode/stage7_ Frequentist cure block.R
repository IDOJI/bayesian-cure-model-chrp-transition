# 🔴 Configure: Stage-7 paths, controls, and reporting thresholds ===============================
# 🟧 Configure: OS-aware project root ===============================
sys_name <- Sys.info()[["sysname"]]

project_root <- switch(
  sys_name,
  "Darwin"  = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  "Windows" = "C:/Users/clair/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  stop("Unsupported OS: ", sys_name)
)

# 🟧 Configure: Stage folder paths ===============================
run_root_dir <- file.path(project_root, "stage7_Frequentist cure block")
export_path <- file.path(project_root, "stage7_Frequentist cure block")
stage1_data_path <- file.path(project_root, "stage1_Backbone lock")
stage6_screening_file <- file.path(project_root, "stage6_Cure-appropriateness screening", "stage6_carry_forward_flag_table.csv")
stage5_root_dir <- file.path(project_root, "stage5_Individualized no-cure comparator")

stage7_refresh_export_suffix <- ""

skip_stage7_core_refit_if_outputs_exist <- TRUE
force_stage7_core_refit <- TRUE
save_stage7_fit_rds <- FALSE
stage7_model_rds_prefix <- "stage7_fit_object"
stage7_expected_bootstrap_reps <- 0L

refresh_patched_exports_even_when_outputs_exist <- TRUE
refresh_ipcw_registry_even_when_output_exists <- FALSE
refresh_new_spec_tables_even_when_outputs_exist <- TRUE
regenerate_visuals_even_when_outputs_exist <- TRUE
run_bootstrap_extraction <- TRUE
refresh_bootstrap_extraction_even_when_outputs_exist <- FALSE
write_full_metric_values_long_csv <- FALSE
write_individual_plot_pngs <- TRUE
include_pngs_in_manifest <- TRUE

main_risk_scale <- "transition_only_main"
common_horizons_year <- 1:10
threshold_plot_horizons_year <- c(1L, 2L, 5L, 10L)
pnu_threshold_suppress_from_year <- 3L
plot_width_in <- 12
plot_height_in <- 8
plot_dpi <- 320

top_n_variable_metric_cells <- 20L
top_n_failure_patterns <- 30L
top_n_failure_models <- 30L

# 🔴 Initialize: packages and global options ===============================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(survival)
  library(stringr)
})

options(stringsAsFactors = FALSE, scipen = 999)
requested_run_root_dir <- run_root_dir
requested_export_path <- export_path

# 🔴 Define: reusable file helpers ===============================
assert_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

normalize_existing_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

safe_read_csv <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
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

resolve_stage7_root_dir <- function(root_dir) {
  root_dir <- normalizePath(root_dir, winslash = "/", mustWork = FALSE)
  
  if (!dir.exists(root_dir)) {
    stop(sprintf("Supplied run_root_dir does not exist: %s", root_dir), call. = FALSE)
  }
  
  if (file.exists(file.path(root_dir, "stage7_fit_registry.csv"))) {
    return(root_dir)
  }
  
  subdirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
  candidate_dirs <- subdirs[file.exists(file.path(subdirs, "stage7_fit_registry.csv"))]
  
  if (length(candidate_dirs) == 1L) {
    message(sprintf("ℹ️  Stage 7 root auto-resolved to nested run directory: %s", candidate_dirs[[1]]))
    return(normalizePath(candidate_dirs[[1]], winslash = "/", mustWork = FALSE))
  }
  
  if (length(candidate_dirs) > 1L) {
    stop(
      paste0(
        "Multiple Stage 7 run directories were found under the supplied run_root_dir. ",
        "Please point run_root_dir to the exact folder containing stage7_fit_registry.csv.\nCandidates:\n- ",
        paste(candidate_dirs, collapse = "\n- ")
      ),
      call. = FALSE
    )
  }
  
  stop(
    paste0(
      "Could not find stage7_fit_registry.csv in the supplied run_root_dir or its immediate subdirectories: ",
      root_dir
    ),
    call. = FALSE
  )
}

resolve_stage7_export_path <- function(export_path, resolved_run_root_dir, suffix = "__stage1_backbone_lock_refresh") {
  export_missing <- is.null(export_path) || length(export_path) == 0L || is.na(export_path) || !nzchar(export_path)
  if (isTRUE(export_missing)) {
    if (!is.na(suffix) && nzchar(suffix)) {
      export_path <- file.path(dirname(resolved_run_root_dir), paste0(basename(resolved_run_root_dir), suffix))
    } else {
      export_path <- resolved_run_root_dir
    }
  }
  
  normalizePath(export_path, winslash = "/", mustWork = FALSE)
}

safe_write_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(df, path)
  normalize_existing_path(path)
}

ensure_required_cols_present <- function(df, required_cols = character()) {
  df <- tibble::as_tibble(df)
  if (length(required_cols) == 0L) {
    return(df)
  }
  n0 <- nrow(df)
  for (nm in required_cols) {
    if (!nm %in% names(df)) {
      df[[nm]] <- if (n0 == 0L) character(0) else rep(NA_character_, n0)
    }
  }
  df
}

safe_save_rds <- function(object, path) {
  saveRDS(object, file = path)
  normalize_existing_path(path)
}

make_output_path <- function(file_name) {
  file.path(export_path, file_name)
}

file_has_required_columns <- function(path, required_cols) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  if (length(required_cols) == 0L) {
    return(TRUE)
  }
  existing <- tryCatch(safe_read_csv(path), error = function(e) NULL)
  if (is.null(existing)) {
    return(FALSE)
  }
  all(required_cols %in% names(existing))
}

should_rebuild_output <- function(path, force_rebuild = FALSE, required_cols = character()) {
  if (!file.exists(path)) {
    return(TRUE)
  }
  if (isTRUE(force_rebuild)) {
    return(TRUE)
  }
  !file_has_required_columns(path, required_cols)
}

load_or_rebuild_csv_output <- function(path, builder_fn, force_rebuild = FALSE, required_cols = character()) {
  if (!should_rebuild_output(path, force_rebuild = force_rebuild, required_cols = required_cols)) {
    return(safe_read_csv(path))
  }
  out <- builder_fn()
  out <- ensure_required_cols_present(out, required_cols = required_cols)
  safe_write_csv(out, path)
  out
}

# 🔴 Define: scalar and table helpers ===============================
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

safe_divide <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

coalesce_character <- function(...) {
  xs <- list(...)
  if (length(xs) == 0L) {
    return(character(0))
  }
  out <- as.character(xs[[1]])
  if (length(xs) > 1L) {
    for (ii in 2:length(xs)) {
      cand <- as.character(xs[[ii]])
      idx <- is.na(out) | out == ""
      out[idx] <- cand[idx]
    }
  }
  out
}

col_or_chr <- function(df, nm, default = NA_character_) {
  if (!nm %in% names(df)) return(rep(default, nrow(df)))
  as.character(df[[nm]])
}

col_or_num <- function(df, nm, default = NA_real_) {
  if (!nm %in% names(df)) return(rep(default, nrow(df)))
  suppressWarnings(as.numeric(df[[nm]]))
}

col_or_lgl <- function(df, nm, default = NA) {
  if (!nm %in% names(df)) return(rep(default, nrow(df)))
  as.logical(df[[nm]])
}

safe_quantile <- function(x, prob) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, names = FALSE, na.rm = TRUE))
}

collapse_unique_sorted <- function(x, empty_value = "") {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) return(empty_value)
  paste(sort(x), collapse = "|")
}

is_success_status <- function(x) {
  stringr::str_starts(dplyr::coalesce(as.character(x), ""), "ok")
}

make_metric_key <- function(df) {
  paste0(
    df$table_name,
    " | ", ifelse(is.na(df$dataset) | df$dataset == "", "<no_dataset>", df$dataset),
    " | ", ifelse(is.na(df$formula_variant) | df$formula_variant == "", "<no_formula>", df$formula_variant),
    " | ", ifelse(is.na(df$model_id) | df$model_id == "", "<no_model>", df$model_id),
    ifelse(!is.na(df$component) & nzchar(df$component), paste0(" | ", df$component), ""),
    ifelse(!is.na(df$term) & nzchar(df$term), paste0(" | ", df$term), ""),
    ifelse(!is.na(df$horizon_year), paste0(" | h=", df$horizon_year), ""),
    ifelse(!is.na(df$threshold), paste0(" | thr=", format(df$threshold, trim = TRUE, scientific = FALSE)), ""),
    " | ", df$metric_name
  )
}

first_existing_value <- function(df, candidates, default = NA_real_) {
  if (nrow(df) == 0L) {
    return(numeric(0))
  }
  out <- rep(default, nrow(df))
  for (nm in candidates) {
    if (!nm %in% names(df)) next
    cand <- suppressWarnings(as.numeric(df[[nm]]))
    idx <- is.na(out) & !is.na(cand)
    out[idx] <- cand[idx]
  }
  out
}

first_existing_character <- function(df, candidates, default = NA_character_) {
  if (nrow(df) == 0L) {
    return(character(0))
  }
  out <- rep(default, nrow(df))
  for (nm in candidates) {
    if (!nm %in% names(df)) next
    cand <- as.character(df[[nm]])
    idx <- (is.na(out) | out == "") & !is.na(cand) & cand != ""
    out[idx] <- cand[idx]
  }
  out
}

first_existing_column_name <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0L) return(NA_character_)
  hit[[1]]
}

ensure_columns_exist <- function(df, col_defaults) {
  n0 <- nrow(df)
  for (nm in names(col_defaults)) {
    if (!nm %in% names(df)) {
      default_value <- col_defaults[[nm]]
      if (length(default_value) == n0) {
        df[[nm]] <- default_value
      } else if (length(default_value) == 1L) {
        df[[nm]] <- rep(default_value, n0)
      } else if (n0 == 0L) {
        df[[nm]] <- default_value[0]
      } else {
        df[[nm]] <- rep_len(default_value, n0)
      }
    }
  }
  df
}

empty_like_n <- function(example, n) {
  if (is.logical(example)) return(rep(NA, n))
  if (is.integer(example)) return(rep(NA_integer_, n))
  if (is.numeric(example)) return(rep(NA_real_, n))
  if (inherits(example, "Date")) return(as.Date(rep(NA_real_, n), origin = "1970-01-01"))
  rep(NA_character_, n)
}

# 🔴 Define: Stage 1 backbone readers ===============================
load_stage1_backbone <- function(stage1_dir) {
  required_files <- c(
    "stage1_backbone_bundle.rds",
    "stage1_analysis_datasets.rds",
    "stage1_formula_registry.csv",
    "stage1_modeling_registry.csv",
    "stage1_horizon_registry.csv",
    "stage1_threshold_registry.csv"
  )
  
  for (file_name in required_files) {
    assert_file_exists(file.path(stage1_dir, file_name), paste("Stage 1", file_name))
  }
  
  bundle <- readRDS(file.path(stage1_dir, "stage1_backbone_bundle.rds"))
  analysis_datasets <- readRDS(file.path(stage1_dir, "stage1_analysis_datasets.rds"))
  formula_registry <- safe_read_csv(file.path(stage1_dir, "stage1_formula_registry.csv"))
  modeling_registry <- safe_read_csv(file.path(stage1_dir, "stage1_modeling_registry.csv"))
  horizon_registry <- safe_read_csv(file.path(stage1_dir, "stage1_horizon_registry.csv"))
  threshold_registry <- safe_read_csv(file.path(stage1_dir, "stage1_threshold_registry.csv"))
  
  expected_horizons <- sort(unique(as.integer(common_horizons_year)))
  stage1_horizons <- sort(unique(as.integer(horizon_registry$horizon_year)))
  
  if (!identical(stage1_horizons, expected_horizons)) {
    stop(
      sprintf(
        "Stage 1 horizon registry mismatch. Expected horizons %s but found %s.",
        paste(expected_horizons, collapse = ","),
        paste(stage1_horizons, collapse = ",")
      ),
      call. = FALSE
    )
  }
  
  list(
    bundle = bundle,
    analysis_datasets = analysis_datasets,
    formula_registry = formula_registry,
    modeling_registry = modeling_registry,
    horizon_registry = horizon_registry,
    threshold_registry = threshold_registry,
    main_risk_scale = bundle$config$main_risk_scale %||% main_risk_scale,
    common_horizons_year = as.integer(bundle$config$common_horizons_year %||% common_horizons_year)
  )
}

prepare_formula_lookup <- function(formula_registry) {
  formula_registry %>%
    mutate(
      formula_variant = as.character(formula_name),
      formula_label = as.character(formula_label),
      site_branch = as.character(site_branch),
      interaction_branch = as.character(interaction_branch),
      uses_site = as.logical(uses_site),
      uses_age_sex_interaction = as.logical(uses_age_sex_interaction),
      risk_scale = as.character(risk_scale)
    ) %>%
    select(
      dataset,
      formula_variant,
      formula_label,
      formula_rhs,
      formula_full,
      uses_site,
      uses_age_sex_interaction,
      site_branch,
      interaction_branch,
      risk_scale,
      formula_scope,
      site_term_interpretation
    )
}

prepare_horizon_lookup <- function(horizon_registry) {
  keep_cols <- intersect(
    c(
      "dataset",
      "dataset_key",
      "horizon_year",
      "support_tier",
      "interpretation_tier",
      "horizon_evidence_class",
      "claim_restriction_flag",
      "interpretation_note"
    ),
    names(horizon_registry)
  )
  
  horizon_registry %>%
    mutate(
      dataset = as.character(dataset),
      horizon_year = as.integer(horizon_year)
    ) %>%
    select(all_of(keep_cols))
}

# 🔴 Define: Stage 6 carry-forward readers ===============================
normalize_stage6_flag_text <- function(x) {
  out <- tolower(trimws(as.character(x)))
  out[out %in% c("", "na", "nan", "null")] <- NA_character_
  out
}

coerce_logical_like <- function(x) {
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

make_empty_stage6_carry_forward_tbl <- function() {
  tibble(
    dataset_key = character(),
    source_dataset = character(),
    source_description = character(),
    analysis_variant = character(),
    screening_context = character(),
    primary_gate_method = character(),
    primary_gate_flag = character(),
    receus_aic_flag = character(),
    cure_model_eligibility_flag = character(),
    final_decision_flag = character(),
    receus_primary_class = character(),
    presence_modifier_flag = logical(),
    cure_presence_support_flag = logical(),
    presence_support_flag = logical(),
    followup_contradiction_flag = logical(),
    followup_not_contradicted_flag = logical(),
    followup_support_flag = logical(),
    xie_centered_bootstrap_p_value = double(),
    descriptive_tail_summary_flag = character(),
    supporting_methods = character(),
    contradicting_methods = character(),
    screening_note = character(),
    common_horizon_vector = character(),
    common_threshold_vector = character(),
    stage6_final_class = character(),
    carry_forward_stage7 = logical(),
    carry_forward_stage8 = logical(),
    screening_repeated_in_stage8 = logical()
  )
}

normalize_stage6_carry_forward <- function(tbl) {
  tbl <- tibble::as_tibble(tbl)
  if (nrow(tbl) == 0L) {
    return(make_empty_stage6_carry_forward_tbl())
  }
  
  tbl <- ensure_columns_exist(
    tbl,
    list(
      source_dataset = NA_character_,
      source_description = NA_character_,
      analysis_variant = NA_character_,
      dataset = NA_character_,
      formula_variant = NA_character_,
      site_branch = NA_character_
    )
  )
  
  if (!"dataset_key" %in% names(tbl)) {
    dataset_like <- first_existing_character(tbl, c("dataset", "source_dataset"), default = NA_character_)
    analysis_like <- first_existing_character(tbl, c("analysis_variant", "site_branch", "formula_variant"), default = NA_character_)
    tbl$dataset_key <- dplyr::case_when(
      dataset_like %in% c("PNU", "SNU") ~ dataset_like,
      dataset_like == "merged" & analysis_like %in% c("site_adjusted", "site-added", "site_interaction", "site_added") ~ "merged__site_adjusted",
      dataset_like == "merged" ~ "merged__site_free",
      TRUE ~ dataset_like
    )
  }
  
  out <- tbl %>%
    mutate(
      dataset_key = as.character(dataset_key),
      source_dataset = dplyr::coalesce(
        first_existing_character(tbl, c("source_dataset", "dataset"), default = NA_character_),
        dplyr::case_when(
          dataset_key %in% c("PNU", "SNU") ~ dataset_key,
          dataset_key %in% c("merged__site_free", "merged__site_adjusted") ~ "merged",
          TRUE ~ NA_character_
        )
      ),
      source_description = first_existing_character(tbl, c("source_description"), default = NA_character_),
      analysis_variant = dplyr::coalesce(
        first_existing_character(tbl, c("analysis_variant", "site_branch", "formula_variant"), default = NA_character_),
        dplyr::case_when(
          dataset_key %in% c("PNU", "SNU") ~ "single_cohort",
          dataset_key == "merged__site_adjusted" ~ "site_adjusted",
          dataset_key == "merged__site_free" ~ "site_free",
          TRUE ~ NA_character_
        )
      ),
      screening_context = first_existing_character(tbl, c("screening_context"), default = NA_character_),
      primary_gate_method = first_existing_character(tbl, c("primary_gate_method"), default = NA_character_),
      primary_gate_flag = normalize_stage6_flag_text(first_existing_character(tbl, c("primary_gate_flag", "receus_aic_flag"), default = NA_character_)),
      receus_aic_flag = normalize_stage6_flag_text(first_existing_character(tbl, c("receus_aic_flag", "primary_gate_flag"), default = NA_character_)),
      cure_model_eligibility_flag = normalize_stage6_flag_text(first_existing_character(tbl, c("cure_model_eligibility_flag", "final_decision_flag", "stage6_final_class"), default = NA_character_)),
      final_decision_flag = normalize_stage6_flag_text(first_existing_character(tbl, c("final_decision_flag", "cure_model_eligibility_flag", "stage6_final_class"), default = NA_character_)),
      receus_primary_class = normalize_stage6_flag_text(first_existing_character(tbl, c("receus_primary_class", "receus_aic_flag", "primary_gate_flag"), default = NA_character_)),
      presence_modifier_flag = dplyr::coalesce(col_or_lgl(tbl, "presence_modifier_flag"), col_or_lgl(tbl, "presence_support_flag")),
      cure_presence_support_flag = dplyr::coalesce(
        col_or_lgl(tbl, "cure_presence_support_flag"),
        col_or_lgl(tbl, "presence_support_flag"),
        ifelse(final_decision_flag == "supportive", TRUE, ifelse(final_decision_flag == "unsupportive", FALSE, NA))
      ),
      presence_support_flag = dplyr::coalesce(
        col_or_lgl(tbl, "presence_support_flag"),
        col_or_lgl(tbl, "cure_presence_support_flag"),
        cure_presence_support_flag
      ),
      followup_contradiction_flag = dplyr::coalesce(
        col_or_lgl(tbl, "followup_contradiction_flag"),
        ifelse(!is.na(col_or_lgl(tbl, "followup_not_contradicted_flag")), !col_or_lgl(tbl, "followup_not_contradicted_flag"), NA)
      ),
      followup_not_contradicted_flag = dplyr::coalesce(
        col_or_lgl(tbl, "followup_not_contradicted_flag"),
        ifelse(!is.na(followup_contradiction_flag), !followup_contradiction_flag, NA)
      ),
      followup_support_flag = dplyr::coalesce(
        col_or_lgl(tbl, "followup_support_flag"),
        dplyr::case_when(
          primary_gate_flag == "supportive" & !is.na(followup_not_contradicted_flag) & followup_not_contradicted_flag ~ TRUE,
          primary_gate_flag == "unsupportive" ~ FALSE,
          primary_gate_flag == "supportive" & !is.na(followup_not_contradicted_flag) & !followup_not_contradicted_flag ~ FALSE,
          TRUE ~ NA
        )
      ),
      xie_centered_bootstrap_p_value = col_or_num(tbl, "xie_centered_bootstrap_p_value"),
      descriptive_tail_summary_flag = first_existing_character(tbl, c("descriptive_tail_summary_flag"), default = NA_character_),
      supporting_methods = first_existing_character(tbl, c("supporting_methods"), default = NA_character_),
      contradicting_methods = first_existing_character(tbl, c("contradicting_methods"), default = NA_character_),
      screening_note = first_existing_character(tbl, c("screening_note"), default = NA_character_),
      common_horizon_vector = first_existing_character(tbl, c("common_horizon_vector"), default = paste(common_horizons_year, collapse = "|")),
      common_threshold_vector = first_existing_character(tbl, c("common_threshold_vector"), default = NA_character_),
      stage6_final_class = normalize_stage6_flag_text(first_existing_character(tbl, c("stage6_final_class", "cure_model_eligibility_flag", "final_decision_flag"), default = NA_character_)),
      carry_forward_stage7 = dplyr::coalesce(col_or_lgl(tbl, "carry_forward_stage7"), final_decision_flag != "unsupportive"),
      carry_forward_stage8 = dplyr::coalesce(col_or_lgl(tbl, "carry_forward_stage8"), final_decision_flag != "unsupportive"),
      screening_repeated_in_stage8 = dplyr::coalesce(col_or_lgl(tbl, "screening_repeated_in_stage8"), FALSE)
    )
  
  out %>%
    select(names(make_empty_stage6_carry_forward_tbl()))
}

read_stage6_carry_forward <- function(path) {
  assert_file_exists(path, "Stage 6 screening file")
  ext <- tolower(tools::file_ext(path))
  
  if (ext == "csv") {
    tbl <- safe_read_csv(path)
  } else if (ext == "rds") {
    obj <- readRDS(path)
    if (is.data.frame(obj)) {
      tbl <- tibble::as_tibble(obj)
    } else if (is.list(obj) && !is.null(obj$outputs) && !is.null(obj$outputs$carry_forward_flag_table)) {
      tbl <- tibble::as_tibble(obj$outputs$carry_forward_flag_table)
    } else {
      stop("Unsupported Stage 6 RDS structure. Expected a data frame or screening bundle with outputs$carry_forward_flag_table.", call. = FALSE)
    }
  } else {
    stop("Stage 6 screening file must be either CSV or RDS.", call. = FALSE)
  }
  
  normalize_stage6_carry_forward(tbl)
}

load_stage6_carry_forward_optional <- function(path = NA_character_) {
  if (is.null(path) || length(path) == 0L || is.na(path) || !nzchar(path)) {
    return(make_empty_stage6_carry_forward_tbl())
  }
  
  if (!file.exists(path)) {
    warning(
      sprintf("Stage 6 screening file not found: %s. Proceeding without Stage 6 carry-forward.", path),
      call. = FALSE
    )
    return(make_empty_stage6_carry_forward_tbl())
  }
  
  read_stage6_carry_forward(path)
}

make_stage6_join_key <- function(dataset, site_branch = NA_character_) {
  dataset_chr <- as.character(dataset)
  site_branch_chr <- as.character(site_branch)
  
  dplyr::case_when(
    dataset_chr %in% c("PNU", "SNU") ~ dataset_chr,
    dataset_chr == "merged" & site_branch_chr == "site_adjusted" ~ "merged__site_adjusted",
    dataset_chr == "merged" ~ "merged__site_free",
    TRUE ~ NA_character_
  )
}

make_empty_stage6_registry <- function() {
  out <- make_empty_stage6_carry_forward_tbl() %>%
    mutate(
      stage6_join_key = character(),
      stage6_flag_joined = logical()
    ) %>%
    rename_with(
      .fn = ~ ifelse(.x %in% c("stage6_join_key", "stage6_flag_joined"), .x, paste0("stage6__", .x))
    )
  
  out
}

standardize_stage6_registry_for_join <- function(stage6_registry) {
  template <- make_empty_stage6_registry()
  
  if (is.null(stage6_registry)) {
    return(template)
  }
  
  out <- tibble::as_tibble(stage6_registry)
  if (!"stage6_join_key" %in% names(out)) {
    out <- out %>%
      normalize_stage6_carry_forward() %>%
      build_stage6_registry()
  }
  
  for (nm in names(template)) {
    if (!nm %in% names(out)) {
      out[[nm]] <- empty_like_n(template[[nm]], nrow(out))
    }
  }
  
  out %>%
    mutate(
      stage6_join_key = as.character(stage6_join_key),
      stage6_flag_joined = dplyr::coalesce(as.logical(stage6_flag_joined), FALSE)
    ) %>%
    select(all_of(names(template)))
}

build_stage6_registry <- function(stage6_tbl) {
  stage6_tbl %>%
    normalize_stage6_carry_forward() %>%
    mutate(
      stage6_join_key = as.character(dataset_key),
      stage6_flag_joined = TRUE
    ) %>%
    rename_with(
      .fn = ~ ifelse(.x %in% c("stage6_join_key", "stage6_flag_joined"), .x, paste0("stage6__", .x))
    ) %>%
    standardize_stage6_registry_for_join()
}

join_stage6_flags <- function(df, stage6_registry) {
  df <- tibble::as_tibble(df)
  df <- ensure_columns_exist(df, list(dataset = NA_character_, site_branch = NA_character_))
  stage6_registry <- standardize_stage6_registry_for_join(stage6_registry)
  
  existing_stage6_cols <- grep("^stage6", names(df), value = TRUE)
  join_key <- make_stage6_join_key(dataset = col_or_chr(df, "dataset"), site_branch = col_or_chr(df, "site_branch"))
  
  out <- df %>%
    select(-any_of(existing_stage6_cols)) %>%
    mutate(stage6_join_key = join_key) %>%
    left_join(stage6_registry, by = "stage6_join_key")
  
  if (!"stage6_flag_joined" %in% names(out)) {
    out$stage6_flag_joined <- rep(FALSE, nrow(out))
    joined_candidates <- intersect(c("stage6_flag_joined.x", "stage6_flag_joined.y"), names(out))
    for (nm in joined_candidates) {
      out$stage6_flag_joined <- dplyr::coalesce(as.logical(out[[nm]]), as.logical(out$stage6_flag_joined))
    }
  }
  
  if (!"stage6__screening_note" %in% names(out)) {
    out$stage6__screening_note <- rep(NA_character_, nrow(out))
    note_candidates <- intersect(c("stage6__screening_note.x", "stage6__screening_note.y"), names(out))
    for (nm in note_candidates) {
      out$stage6__screening_note <- dplyr::coalesce(as.character(out[[nm]]), as.character(out$stage6__screening_note))
    }
  }
  
  out %>%
    mutate(
      stage6_flag_joined = dplyr::coalesce(as.logical(stage6_flag_joined), FALSE),
      stage6__screening_note = dplyr::coalesce(
        as.character(stage6__screening_note),
        if_else(stage6_flag_joined, NA_character_, "Stage 6 carry-forward was not joined in this refresh run.")
      )
    ) %>%
    select(-any_of(c(
      "stage6_flag_joined.x", "stage6_flag_joined.y",
      "stage6__screening_note.x", "stage6__screening_note.y"
    )))
}

# 🔴 Define: Stage 7 artifact readers ===============================
load_stage7_outputs <- function(root_dir) {
  resolved_root_dir <- resolve_stage7_root_dir(root_dir)
  
  required_files <- c(
    "stage7_fit_registry.csv",
    "stage7_coefficients.csv",
    "stage7_subject_predictions.csv",
    "stage7_risk_summary.csv",
    "stage7_delta_risk.csv",
    "stage7_threshold_metrics.csv"
  )
  
  missing_required <- required_files[!file.exists(file.path(resolved_root_dir, required_files))]
  if (length(missing_required) > 0L) {
    stop(
      paste0(
        "The resolved Stage 7 folder is missing required output file(s): ",
        paste(missing_required, collapse = ", "),
        "\nResolved folder: ", resolved_root_dir
      ),
      call. = FALSE
    )
  }
  
  run_manifest_tbl <- safe_read_csv_if_exists(file.path(resolved_root_dir, "stage7_run_manifest.csv"))
  if (is.null(run_manifest_tbl)) {
    message("ℹ️  stage7_run_manifest.csv not found; proceeding without it.")
    run_manifest_tbl <- tibble()
  }
  
  bootstrap_model_qc_tbl <- safe_read_csv_if_exists(file.path(resolved_root_dir, "stage7_bootstrap_model_qc.csv"))
  if (is.null(bootstrap_model_qc_tbl)) {
    bootstrap_model_qc_tbl <- tibble()
  }
  
  shard_registry_tbl <- safe_read_csv_if_exists(file.path(resolved_root_dir, "stage7_shard_registry.csv"))
  if (is.null(shard_registry_tbl)) {
    shard_registry_tbl <- tibble()
  }
  
  list(
    resolved_root_dir = resolved_root_dir,
    run_manifest = run_manifest_tbl,
    fit_registry = safe_read_csv(file.path(resolved_root_dir, "stage7_fit_registry.csv")),
    coefficients = safe_read_csv(file.path(resolved_root_dir, "stage7_coefficients.csv")),
    subject_predictions = safe_read_csv(file.path(resolved_root_dir, "stage7_subject_predictions.csv")),
    risk_summary = safe_read_csv(file.path(resolved_root_dir, "stage7_risk_summary.csv")),
    delta_risk = safe_read_csv(file.path(resolved_root_dir, "stage7_delta_risk.csv")),
    threshold_metrics = safe_read_csv(file.path(resolved_root_dir, "stage7_threshold_metrics.csv")),
    bootstrap_model_qc = bootstrap_model_qc_tbl,
    shard_registry = shard_registry_tbl,
    merged_bootstrap_rds_file = file.path(resolved_root_dir, "stage7_bootstrap_merged_results.rds")
  )
}

# 🔴 Define: shared metadata enrichers ===============================
derive_site_placement_label <- function(dataset, site_branch, existing_label = NA_character_, incidence_uses_site = NA, latency_uses_site = NA) {
  existing_label <- as.character(existing_label)
  dataset <- as.character(dataset)
  site_branch <- as.character(site_branch)
  
  out <- rep(NA_character_, length(dataset))
  keep_existing <- !is.na(existing_label) & nzchar(existing_label)
  out[keep_existing] <- existing_label[keep_existing]
  
  need <- !keep_existing
  if (any(need)) {
    inc <- as.logical(incidence_uses_site)
    lat <- as.logical(latency_uses_site)
    
    out[need & dataset != "merged"] <- "not_applicable_single_site"
    out[need & dataset == "merged" & !is.na(inc) & !is.na(lat) & inc & lat] <- "site_in_both"
    out[need & dataset == "merged" & !is.na(inc) & !is.na(lat) & inc & !lat] <- "site_in_incidence_only"
    out[need & dataset == "merged" & !is.na(inc) & !is.na(lat) & !inc & lat] <- "site_in_latency_only"
    out[need & dataset == "merged" & !is.na(inc) & !is.na(lat) & !inc & !lat] <- "site_in_neither"
    out[need & dataset == "merged" & is.na(out) & site_branch == "site_adjusted"] <- "site_in_both"
    out[need & dataset == "merged" & is.na(out) & site_branch != "site_adjusted"] <- "site_in_neither"
  }
  
  out
}

add_common_stage7_metadata <- function(df, formula_lookup, risk_scale_default = main_risk_scale) {
  df <- tibble::as_tibble(df)
  
  df <- ensure_columns_exist(
    df,
    list(
      dataset = NA_character_,
      formula_variant = NA_character_,
      formula_label = NA_character_,
      site_branch = NA_character_,
      interaction_branch = NA_character_,
      risk_scale = NA_character_,
      stage = NA_character_,
      branch = NA_character_,
      formula_scope = NA_character_,
      site_term_interpretation = NA_character_,
      uses_site = NA,
      uses_age_sex_interaction = NA,
      site_placement_label = NA_character_,
      incidence_uses_site = NA,
      latency_uses_site = NA
    )
  )
  
  if (nrow(df) == 0L) {
    return(df)
  }
  
  if ("formula_name" %in% names(df)) {
    df$formula_variant <- dplyr::coalesce(as.character(df$formula_variant), as.character(df$formula_name))
  }
  
  formula_join <- formula_lookup %>%
    rename(
      formula_label_stage1 = formula_label,
      formula_rhs_stage1 = formula_rhs,
      formula_full_stage1 = formula_full,
      uses_site_stage1 = uses_site,
      uses_age_sex_interaction_stage1 = uses_age_sex_interaction,
      site_branch_stage1 = site_branch,
      interaction_branch_stage1 = interaction_branch,
      risk_scale_stage1 = risk_scale,
      formula_scope_stage1 = formula_scope,
      site_term_interpretation_stage1 = site_term_interpretation
    )
  
  out <- df %>%
    mutate(
      dataset = as.character(dataset),
      formula_variant = as.character(formula_variant)
    ) %>%
    left_join(formula_join, by = c("dataset", "formula_variant")) %>%
    mutate(
      formula_label = dplyr::coalesce(as.character(formula_label), formula_label_stage1),
      site_branch = dplyr::coalesce(as.character(site_branch), site_branch_stage1),
      interaction_branch = dplyr::coalesce(as.character(interaction_branch), interaction_branch_stage1),
      risk_scale = dplyr::coalesce(as.character(risk_scale), risk_scale_stage1, risk_scale_default),
      formula_scope = dplyr::coalesce(as.character(formula_scope), formula_scope_stage1),
      site_term_interpretation = dplyr::coalesce(as.character(site_term_interpretation), site_term_interpretation_stage1),
      uses_site = dplyr::coalesce(as.logical(uses_site), as.logical(uses_site_stage1)),
      uses_age_sex_interaction = dplyr::coalesce(as.logical(uses_age_sex_interaction), as.logical(uses_age_sex_interaction_stage1)),
      incidence_uses_site = dplyr::coalesce(as.logical(incidence_uses_site), as.logical(uses_site_stage1)),
      latency_uses_site = dplyr::coalesce(as.logical(latency_uses_site), as.logical(uses_site_stage1)),
      stage = dplyr::coalesce(as.character(stage), "Stage 7"),
      branch = dplyr::coalesce(as.character(branch), "Stage7"),
      dataset_key = make_stage6_join_key(dataset = dataset, site_branch = site_branch),
      site_placement_label = derive_site_placement_label(
        dataset = dataset,
        site_branch = site_branch,
        existing_label = site_placement_label,
        incidence_uses_site = incidence_uses_site,
        latency_uses_site = latency_uses_site
      )
    ) %>%
    select(-any_of(c(
      "formula_label_stage1",
      "formula_rhs_stage1",
      "formula_full_stage1",
      "uses_site_stage1",
      "uses_age_sex_interaction_stage1",
      "site_branch_stage1",
      "interaction_branch_stage1",
      "risk_scale_stage1",
      "formula_scope_stage1",
      "site_term_interpretation_stage1"
    )))
  
  out
}

join_horizon_support <- function(df, horizon_lookup) {
  if (nrow(df) == 0L || !"horizon_year" %in% names(df)) {
    return(df)
  }
  
  existing_cols <- intersect(
    c("support_tier", "interpretation_tier", "horizon_evidence_class", "claim_restriction_flag", "interpretation_note"),
    names(df)
  )
  
  df %>%
    mutate(
      dataset = as.character(dataset),
      horizon_year = as.integer(horizon_year)
    ) %>%
    select(-any_of(existing_cols)) %>%
    left_join(horizon_lookup, by = c("dataset", "horizon_year"))
}

# 🔴 Define: IPCW rebuilding logic ===============================
compute_ipcw_frame <- function(df, horizons_year) {
  censor_fit <- survival::survfit(survival::Surv(time_year, censor_main) ~ 1, data = df)
  
  step_value <- function(x, y, at, yleft = 1, yright = NULL, lower = 0, upper = 1) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    at <- as.numeric(at)
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
    if (length(x) == 0L || length(y) == 0L) {
      out <- rep(yleft, length(at))
      return(pmin(pmax(out, lower), upper))
    }
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    if (is.null(yright)) {
      yright <- tail(y, 1)
    }
    out <- approx(
      x = x,
      y = y,
      xout = at,
      method = "constant",
      f = 0,
      yleft = yleft,
      yright = yright,
      ties = "ordered"
    )$y
    pmin(pmax(out, lower), upper)
  }
  
  G_at_t <- step_value(censor_fit$time, censor_fit$surv, horizons_year, yleft = 1, lower = 0, upper = 1)
  G_at_obs <- step_value(censor_fit$time, censor_fit$surv, df$time_year, yleft = 1, lower = 0, upper = 1)
  
  bind_rows(lapply(seq_along(horizons_year), function(idx) {
    hh <- as.numeric(horizons_year[idx])
    G_h <- pmax(G_at_t[idx], 1e-10)
    is_case <- df$event_main == 1L & df$time_year <= hh
    is_control <- df$time_year > hh
    
    tibble(
      dataset = unique(as.character(df$dataset)),
      unique_person_id = df$unique_person_id,
      horizon_year = as.integer(hh),
      event_by_horizon = as.integer(is_case),
      known_case_or_control = as.integer(is_case | is_control),
      weight_case = if_else(is_case, 1 / pmax(G_at_obs, 1e-10), 0),
      weight_control = if_else(is_control, 1 / G_h, 0)
    )
  }))
}

rebuild_ipcw_registry <- function(stage1_analysis_datasets, horizons_year) {
  required_names <- c("PNU", "SNU", "merged")
  if (!is.list(stage1_analysis_datasets) || !all(required_names %in% names(stage1_analysis_datasets))) {
    stop("stage1_analysis_datasets must contain PNU, SNU, and merged.", call. = FALSE)
  }
  
  bind_rows(lapply(required_names, function(nm) {
    df <- tibble::as_tibble(stage1_analysis_datasets[[nm]])
    if (!"dataset" %in% names(df)) {
      df$dataset <- nm
    }
    if (!"unique_person_id" %in% names(df)) {
      df$unique_person_id <- paste(df$site, df$id, sep = "::")
    }
    if (!"time_year" %in% names(df)) {
      df$time_year <- suppressWarnings(as.numeric(df$days_followup)) / 365.25
    }
    if (!"event_main" %in% names(df)) {
      df$event_main <- as.integer(df$status_num == 1L)
    }
    if (!"censor_main" %in% names(df)) {
      df$censor_main <- as.integer(df$status_num %in% c(0L, 2L))
    }
    compute_ipcw_frame(df, horizons_year = horizons_year)
  }))
}

# 🔴 Define: threshold suppression rules ===============================
apply_threshold_reporting_rules <- function(threshold_metrics_tbl) {
  metric_cols_to_blank <- c(
    "positive_classification_rate",
    "weighted_tp",
    "weighted_fn",
    "weighted_fp",
    "weighted_tn",
    "sensitivity",
    "specificity",
    "ppv",
    "false_positive_weighted_per_n",
    "false_positive_burden_non_event",
    "false_positive_burden_primary",
    "false_positive_burden_all",
    "unnecessary_high_risk_per_100_population",
    "unnecessary_high_risk_per_100_non_event",
    "net_benefit",
    "false_positive_burden_primary_ci_lower",
    "false_positive_burden_primary_ci_upper",
    "false_positive_weighted_per_n_ci_lower",
    "false_positive_weighted_per_n_ci_upper",
    "net_benefit_ci_lower",
    "net_benefit_ci_upper",
    "false_positive_burden_ci_lower",
    "false_positive_burden_ci_upper"
  )
  
  existing_metric_cols <- intersect(metric_cols_to_blank, names(threshold_metrics_tbl))
  
  suppression_flag <- threshold_metrics_tbl$dataset == "PNU" & (
    (!is.na(threshold_metrics_tbl$support_tier) & threshold_metrics_tbl$support_tier == "projection") |
      as.integer(threshold_metrics_tbl$horizon_year) >= as.integer(pnu_threshold_suppress_from_year)
  )
  
  out <- threshold_metrics_tbl %>%
    mutate(
      threshold_projection_suppressed_flag = suppression_flag,
      threshold_reporting_status = if_else(threshold_projection_suppressed_flag, "suppressed_projection", "reported"),
      threshold_suppression_reason = if_else(
        threshold_projection_suppressed_flag,
        "PNU threshold-based metrics suppressed from year 3 onward because they are projection-only and not primary-supported.",
        NA_character_
      )
    )
  
  if (length(existing_metric_cols) > 0L) {
    out[existing_metric_cols] <- lapply(out[existing_metric_cols], function(col) {
      col[suppression_flag] <- NA
      col
    })
  }
  
  out
}

# 🔴 Define: supportive decomposition exports ===============================
standardize_fit_metric_columns <- function(df) {
  if (nrow(df) == 0L) {
    return(df)
  }
  
  df %>%
    mutate(
      std_logLik = first_existing_value(df, c("logLik", "fit_logLik", "model_logLik", "maximized_logLik", "max_logLik", "log_likelihood", "ll")),
      std_AIC = first_existing_value(df, c("AIC", "fit_AIC", "model_AIC", "aic")),
      std_BIC = first_existing_value(df, c("BIC", "fit_BIC", "model_BIC", "bic")),
      std_n_parameters = first_existing_value(df, c("n_parameters", "parameter_count", "df_model", "effective_parameter_n")),
      std_convergence_status = first_existing_character(df, c("overall_status", "fit_status", "status", "convergence_status"), default = NA_character_),
      std_model_role = dplyr::case_when(
        stringr::str_detect(tolower(first_existing_character(df, c("model_block", "model_class", "model_id"), default = "")), "no[-_ ]?cure") ~ "no_cure",
        stringr::str_detect(tolower(first_existing_character(df, c("model_block", "model_class", "model_id"), default = "")), "cure") ~ "cure",
        TRUE ~ NA_character_
      )
    )
}

load_stage5_fit_metrics_optional <- function(stage5_root_dir) {
  if (is.null(stage5_root_dir) || length(stage5_root_dir) == 0L || is.na(stage5_root_dir) || !nzchar(stage5_root_dir) || !dir.exists(stage5_root_dir)) {
    return(tibble())
  }
  
  candidate_files <- c(
    "stage5_fit_registry.csv",
    "stage5_model_registry.csv",
    "stage5_family_fit_summary.csv",
    "stage5_model_performance.csv"
  )
  
  existing <- candidate_files[file.exists(file.path(stage5_root_dir, candidate_files))]
  if (length(existing) == 0L) {
    return(tibble())
  }
  
  safe_read_csv(file.path(stage5_root_dir, existing[[1]]))
}

build_family_matched_fit_contrast <- function(fit_registry_tbl, stage5_fit_tbl = tibble()) {
  output_schema <- tibble(
    dataset_key = character(),
    dataset = character(),
    formula_variant = character(),
    site_branch = character(),
    family_pair = character(),
    cure_model_id = character(),
    matched_noncure_model_id = character(),
    logLik_noncure = double(),
    logLik_cure = double(),
    delta_logLik_cure_minus_noncure = double(),
    LR_2delta_logLik = double(),
    delta_AIC_cure_minus_noncure = double(),
    delta_BIC_cure_minus_noncure = double(),
    lrt_calibration_status = character(),
    lrt_pvalue_bootstrap = double(),
    same_family_fit_gain_signal = character(),
    convergence_pair_flag = logical(),
    risk_scale = character()
  )
  
  if (nrow(fit_registry_tbl) == 0L) {
    return(output_schema)
  }
  
  cure_df <- standardize_fit_metric_columns(fit_registry_tbl) %>%
    filter(std_model_role != "no_cure") %>%
    mutate(
      family_pair = dplyr::coalesce(as.character(model_family), as.character(latency_type), as.character(model_id)),
      matched_noncure_model_id = as.character(matched_nocure_model_id),
      dataset_key = dplyr::coalesce(as.character(dataset_key), make_stage6_join_key(dataset, site_branch))
    )
  
  noncure_stage7 <- standardize_fit_metric_columns(fit_registry_tbl) %>%
    filter(std_model_role == "no_cure") %>%
    mutate(
      family_pair = dplyr::coalesce(as.character(model_family), as.character(latency_type), as.character(model_id)),
      dataset_key = dplyr::coalesce(as.character(dataset_key), make_stage6_join_key(dataset, site_branch))
    )
  
  stage5_fit_tbl <- ensure_columns_exist(
    stage5_fit_tbl,
    list(
      dataset = NA_character_,
      source_dataset = NA_character_,
      formula_variant = NA_character_,
      formula_name = NA_character_,
      analysis_variant = NA_character_,
      site_branch = NA_character_,
      model_family = NA_character_,
      latency_type = NA_character_,
      model_id = NA_character_,
      dataset_key = NA_character_
    )
  )
  
  stage5_formula_variant_fallback <- first_existing_character(
    stage5_fit_tbl,
    c("formula_variant", "formula_name", "analysis_variant"),
    default = NA_character_
  )
  
  noncure_stage5 <- standardize_fit_metric_columns(stage5_fit_tbl) %>%
    mutate(
      dataset = dplyr::coalesce(
        as.character(dataset),
        first_existing_character(stage5_fit_tbl, c("source_dataset"), default = NA_character_)
      ),
      formula_variant = dplyr::coalesce(
        as.character(formula_variant),
        stage5_formula_variant_fallback
      ),
      site_branch = dplyr::coalesce(
        as.character(site_branch),
        dplyr::if_else(
          as.character(dataset) == "merged" &
            stringr::str_detect(dplyr::coalesce(stage5_formula_variant_fallback, ""), "site"),
          "site_adjusted",
          "site_free",
          missing = NA_character_
        )
      ),
      family_pair = dplyr::coalesce(as.character(model_family), as.character(latency_type), as.character(model_id)),
      dataset_key = dplyr::coalesce(as.character(dataset_key), make_stage6_join_key(dataset, site_branch))
    )
  
  noncure_df <- bind_rows(noncure_stage7, noncure_stage5) %>%
    filter(!is.na(dataset), !is.na(formula_variant)) %>%
    arrange(dataset_key, formula_variant, family_pair, desc(!is.na(std_logLik))) %>%
    group_by(dataset_key, formula_variant, family_pair, model_id) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  if (nrow(cure_df) == 0L) {
    return(output_schema)
  }
  
  by_model_id <- cure_df %>%
    filter(!is.na(matched_noncure_model_id), matched_noncure_model_id != "") %>%
    left_join(
      noncure_df %>%
        select(
          dataset_key,
          formula_variant,
          model_id,
          std_logLik,
          std_AIC,
          std_BIC,
          std_n_parameters,
          std_convergence_status
        ) %>%
        rename(
          noncure_model_id = model_id,
          logLik_noncure = std_logLik,
          AIC_noncure = std_AIC,
          BIC_noncure = std_BIC,
          n_parameters_noncure = std_n_parameters,
          convergence_status_noncure = std_convergence_status
        ),
      by = c("dataset_key", "formula_variant", "matched_noncure_model_id" = "noncure_model_id")
    )
  
  remaining <- cure_df %>%
    filter(is.na(matched_noncure_model_id) | matched_noncure_model_id == "") %>%
    left_join(
      noncure_df %>%
        select(
          dataset_key,
          formula_variant,
          family_pair,
          model_id,
          std_logLik,
          std_AIC,
          std_BIC,
          std_n_parameters,
          std_convergence_status
        ) %>%
        rename(
          matched_noncure_model_id = model_id,
          logLik_noncure = std_logLik,
          AIC_noncure = std_AIC,
          BIC_noncure = std_BIC,
          n_parameters_noncure = std_n_parameters,
          convergence_status_noncure = std_convergence_status
        ),
      by = c("dataset_key", "formula_variant", "family_pair")
    )
  
  out <- bind_rows(by_model_id, remaining) %>%
    transmute(
      dataset_key = dataset_key,
      dataset = as.character(dataset),
      formula_variant = as.character(formula_variant),
      site_branch = as.character(site_branch),
      family_pair = as.character(family_pair),
      cure_model_id = as.character(model_id),
      matched_noncure_model_id = as.character(matched_noncure_model_id),
      logLik_noncure = as.numeric(logLik_noncure),
      logLik_cure = as.numeric(std_logLik),
      delta_logLik_cure_minus_noncure = logLik_cure - logLik_noncure,
      LR_2delta_logLik = 2 * (logLik_cure - logLik_noncure),
      delta_AIC_cure_minus_noncure = as.numeric(std_AIC) - as.numeric(AIC_noncure),
      delta_BIC_cure_minus_noncure = as.numeric(std_BIC) - as.numeric(BIC_noncure),
      lrt_calibration_status = dplyr::case_when(
        is.na(logLik_noncure) | is.na(logLik_cure) ~ "missing_likelihood_metrics",
        TRUE ~ "not_reported_nonregular_problem"
      ),
      lrt_pvalue_bootstrap = NA_real_,
      same_family_fit_gain_signal = dplyr::case_when(
        is.na(delta_logLik_cure_minus_noncure) ~ "cannot_assess",
        delta_logLik_cure_minus_noncure > 0 & !is.na(delta_AIC_cure_minus_noncure) & delta_AIC_cure_minus_noncure < 0 & !is.na(delta_BIC_cure_minus_noncure) & delta_BIC_cure_minus_noncure < 0 ~ "coherent_cure_fit_gain",
        delta_logLik_cure_minus_noncure > 0 ~ "partial_cure_fit_gain",
        TRUE ~ "no_cure_fit_gain"
      ),
      convergence_pair_flag = is_success_status(std_convergence_status) & is_success_status(convergence_status_noncure),
      risk_scale = dplyr::coalesce(as.character(risk_scale), main_risk_scale)
    ) %>%
    arrange(dataset_key, formula_variant, family_pair, cure_model_id)
  
  bind_rows(output_schema[0, ], out)
}

classify_hazard_shape <- function(hazard_values, cumulative_risk_values) {
  h <- as.numeric(hazard_values)
  F <- as.numeric(cumulative_risk_values)
  
  if (length(h) == 0L || all(is.na(h))) {
    return("insufficient_grid")
  }
  if (any(diff(F) < -1e-6, na.rm = TRUE)) {
    return("nonmonotone_risk")
  }
  
  h_nonmiss <- h[is.finite(h)]
  if (length(h_nonmiss) <= 1L) {
    return("insufficient_grid")
  }
  
  d <- diff(h_nonmiss)
  tol <- max(0.01, stats::median(abs(h_nonmiss), na.rm = TRUE) * 0.10)
  
  if (all(abs(d) <= tol, na.rm = TRUE)) {
    return("flat")
  }
  if (all(d >= -tol, na.rm = TRUE)) {
    return("increasing")
  }
  if (all(d <= tol, na.rm = TRUE)) {
    return("decreasing")
  }
  
  peak_idx <- which.max(h_nonmiss)
  if (peak_idx > 1L && peak_idx < length(h_nonmiss)) {
    left_ok <- all(diff(h_nonmiss[seq_len(peak_idx)]) >= -tol, na.rm = TRUE)
    right_ok <- all(diff(h_nonmiss[peak_idx:length(h_nonmiss)]) <= tol, na.rm = TRUE)
    if (left_ok && right_ok) {
      return("unimodal")
    }
  }
  
  "irregular"
}

build_hazard_shape_plausibility <- function(risk_summary_tbl) {
  output_schema <- tibble(
    dataset = character(),
    dataset_key = character(),
    formula_variant = character(),
    formula_label = character(),
    site_branch = character(),
    interaction_branch = character(),
    model_id = character(),
    model_class = character(),
    model_block = character(),
    model_family = character(),
    latency_type = character(),
    risk_scale = character(),
    hazard_target = double(),
    hazard_1y = double(),
    hazard_2y = double(),
    hazard_3y = double(),
    hazard_4y = double(),
    hazard_5y = double(),
    hazard_6y = double(),
    hazard_7y = double(),
    hazard_8y = double(),
    hazard_9y = double(),
    hazard_10y = double(),
    hazard_ratio_10y_vs_1y = double(),
    shape_class = character(),
    hazard_estimation_method = character(),
    monotonicity_violation_flag = logical()
  )
  
  if (nrow(risk_summary_tbl) == 0L || !"horizon_year" %in% names(risk_summary_tbl)) {
    return(output_schema)
  }
  
  risk_col <- first_existing_column_name(risk_summary_tbl, c("mean_risk_overall", "risk_overall", "overall_risk", "predicted_risk", "risk"))
  if (is.na(risk_col)) {
    return(output_schema)
  }
  
  group_cols <- intersect(
    c("dataset", "dataset_key", "formula_variant", "formula_label", "site_branch", "interaction_branch", "model_id", "model_class", "model_block", "model_family", "latency_type", "risk_scale"),
    names(risk_summary_tbl)
  )
  
  out <- risk_summary_tbl %>%
    mutate(
      horizon_year = as.integer(horizon_year),
      cumulative_risk = pmin(pmax(suppressWarnings(as.numeric(.data[[risk_col]])), 0), 0.999999)
    ) %>%
    filter(!is.na(horizon_year)) %>%
    group_by(across(all_of(group_cols))) %>%
    arrange(horizon_year, .by_group = TRUE) %>%
    summarise(
      hazard_tbl = list({
        df <- tibble(
          horizon_year = horizon_year,
          cumulative_risk = cumulative_risk
        ) %>%
          group_by(horizon_year) %>%
          summarise(cumulative_risk = dplyr::last(cumulative_risk), .groups = "drop") %>%
          arrange(horizon_year) %>%
          mutate(
            prev_risk = dplyr::lag(cumulative_risk, default = 0),
            discrete_hazard = ifelse(cumulative_risk >= prev_risk, (cumulative_risk - prev_risk) / pmax(1 - prev_risk, 1e-8), NA_real_)
          )
        df
      }),
      .groups = "drop"
    ) %>%
    mutate(
      monotonicity_violation_flag = purrr::map_lgl(hazard_tbl, ~ any(diff(.x$cumulative_risk) < -1e-6, na.rm = TRUE)),
      hazard_target = purrr::map_dbl(hazard_tbl, ~ {
        vals <- .x$discrete_hazard[is.finite(.x$discrete_hazard)]
        if (length(vals) == 0L) NA_real_ else tail(vals, 1)
      }),
      hazard_1y = purrr::map_dbl(hazard_tbl, ~ dplyr::coalesce(.x$discrete_hazard[match(1L, .x$horizon_year)], NA_real_)),
      hazard_2y = purrr::map_dbl(hazard_tbl, ~ dplyr::coalesce(.x$discrete_hazard[match(2L, .x$horizon_year)], NA_real_)),
      hazard_3y = purrr::map_dbl(hazard_tbl, ~ dplyr::coalesce(.x$discrete_hazard[match(3L, .x$horizon_year)], NA_real_)),
      hazard_4y = purrr::map_dbl(hazard_tbl, ~ dplyr::coalesce(.x$discrete_hazard[match(4L, .x$horizon_year)], NA_real_)),
      hazard_5y = purrr::map_dbl(hazard_tbl, ~ dplyr::coalesce(.x$discrete_hazard[match(5L, .x$horizon_year)], NA_real_)),
      hazard_6y = purrr::map_dbl(hazard_tbl, ~ dplyr::coalesce(.x$discrete_hazard[match(6L, .x$horizon_year)], NA_real_)),
      hazard_7y = purrr::map_dbl(hazard_tbl, ~ dplyr::coalesce(.x$discrete_hazard[match(7L, .x$horizon_year)], NA_real_)),
      hazard_8y = purrr::map_dbl(hazard_tbl, ~ dplyr::coalesce(.x$discrete_hazard[match(8L, .x$horizon_year)], NA_real_)),
      hazard_9y = purrr::map_dbl(hazard_tbl, ~ dplyr::coalesce(.x$discrete_hazard[match(9L, .x$horizon_year)], NA_real_)),
      hazard_10y = purrr::map_dbl(hazard_tbl, ~ dplyr::coalesce(.x$discrete_hazard[match(10L, .x$horizon_year)], NA_real_)),
      hazard_ratio_10y_vs_1y = ifelse(is.finite(hazard_10y) & is.finite(hazard_1y) & hazard_1y > 0, hazard_10y / hazard_1y, NA_real_),
      shape_class = purrr::map_chr(hazard_tbl, ~ classify_hazard_shape(.x$discrete_hazard, .x$cumulative_risk)),
      hazard_estimation_method = "discrete_increment_from_cumulative_mean_risk"
    ) %>%
    select(-hazard_tbl)
  
  bind_rows(output_schema[0, ], out)
}

build_cure_supporting_decomposition <- function(fit_registry_tbl, risk_summary_tbl) {
  output_schema <- tibble(
    dataset = character(),
    dataset_key = character(),
    formula_variant = character(),
    formula_label = character(),
    site_branch = character(),
    interaction_branch = character(),
    site_placement_label = character(),
    model_id = character(),
    model_class = character(),
    model_family = character(),
    latency_type = character(),
    risk_scale = character(),
    horizon_year = integer(),
    cure_fraction = double(),
    susceptible_fraction = double(),
    uncured_survival = double(),
    uncured_risk = double(),
    MSTu = double(),
    uncured_mean_support_flag = character(),
    source_uncured_survival_col = character(),
    source_uncured_risk_col = character(),
    stage6__cure_model_eligibility_flag = character(),
    stage6__followup_not_contradicted_flag = logical()
  )
  
  if (nrow(fit_registry_tbl) == 0L && nrow(risk_summary_tbl) == 0L) {
    return(output_schema)
  }
  
  model_cols <- intersect(
    c("dataset", "dataset_key", "formula_variant", "formula_label", "site_branch", "interaction_branch", "site_placement_label", "model_id", "model_class", "model_family", "latency_type", "risk_scale", "stage6__cure_model_eligibility_flag", "stage6__followup_not_contradicted_flag"),
    names(fit_registry_tbl)
  )
  
  model_level <- fit_registry_tbl %>%
    mutate(
      cure_fraction = first_existing_value(fit_registry_tbl, c("cure_fraction", "estimated_cure_fraction", "mean_cure_fraction", "pi_hat", "cure_rate", "cure_prob", "prob_cure")),
      susceptible_fraction = first_existing_value(fit_registry_tbl, c("susceptible_fraction", "estimated_susceptible_fraction", "mean_susceptible_fraction")),
      MSTu = first_existing_value(fit_registry_tbl, c("MSTu", "mstu", "uncured_mean_survival", "mean_survival_uncured", "mst_uncured"))
    ) %>%
    mutate(
      susceptible_fraction = ifelse(is.na(susceptible_fraction) & !is.na(cure_fraction), 1 - cure_fraction, susceptible_fraction),
      cure_fraction = ifelse(is.na(cure_fraction) & !is.na(susceptible_fraction), 1 - susceptible_fraction, cure_fraction)
    ) %>%
    select(all_of(c(model_cols, "cure_fraction", "susceptible_fraction", "MSTu"))) %>%
    distinct()
  
  uncured_survival_col <- first_existing_column_name(risk_summary_tbl, c("mean_survival_uncured", "uncured_survival", "mean_survival_susceptible", "susceptible_survival", "survival_uncured", "survival_susceptible"))
  uncured_risk_col <- first_existing_column_name(risk_summary_tbl, c("mean_risk_uncured", "uncured_risk", "mean_risk_susceptible", "susceptible_risk", "risk_uncured", "risk_susceptible"))
  
  annual_tbl <- tibble()
  if (!is.na(uncured_survival_col) || !is.na(uncured_risk_col)) {
    annual_tbl <- risk_summary_tbl %>%
      transmute(
        dataset = if ("dataset" %in% names(risk_summary_tbl)) as.character(dataset) else NA_character_,
        dataset_key = if ("dataset_key" %in% names(risk_summary_tbl)) as.character(dataset_key) else make_stage6_join_key(as.character(dataset), as.character(site_branch)),
        formula_variant = if ("formula_variant" %in% names(risk_summary_tbl)) as.character(formula_variant) else NA_character_,
        formula_label = if ("formula_label" %in% names(risk_summary_tbl)) as.character(formula_label) else NA_character_,
        site_branch = if ("site_branch" %in% names(risk_summary_tbl)) as.character(site_branch) else NA_character_,
        interaction_branch = if ("interaction_branch" %in% names(risk_summary_tbl)) as.character(interaction_branch) else NA_character_,
        site_placement_label = if ("site_placement_label" %in% names(risk_summary_tbl)) as.character(site_placement_label) else NA_character_,
        model_id = if ("model_id" %in% names(risk_summary_tbl)) as.character(model_id) else NA_character_,
        model_class = if ("model_class" %in% names(risk_summary_tbl)) as.character(model_class) else NA_character_,
        model_family = if ("model_family" %in% names(risk_summary_tbl)) as.character(model_family) else NA_character_,
        latency_type = if ("latency_type" %in% names(risk_summary_tbl)) as.character(latency_type) else NA_character_,
        risk_scale = if ("risk_scale" %in% names(risk_summary_tbl)) as.character(risk_scale) else main_risk_scale,
        horizon_year = if ("horizon_year" %in% names(risk_summary_tbl)) as.integer(horizon_year) else NA_integer_,
        uncured_survival = if (!is.na(uncured_survival_col)) suppressWarnings(as.numeric(.data[[uncured_survival_col]])) else NA_real_,
        uncured_risk = if (!is.na(uncured_risk_col)) suppressWarnings(as.numeric(.data[[uncured_risk_col]])) else NA_real_
      )
  }
  
  combined <- if (nrow(annual_tbl) > 0L) {
    annual_tbl %>%
      left_join(model_level, by = c("dataset", "dataset_key", "formula_variant", "formula_label", "site_branch", "interaction_branch", "site_placement_label", "model_id", "model_class", "model_family", "latency_type", "risk_scale"))
  } else {
    model_level %>%
      mutate(horizon_year = NA_integer_, uncured_survival = NA_real_, uncured_risk = NA_real_)
  }
  
  out <- combined %>%
    mutate(
      source_uncured_survival_col = uncured_survival_col,
      source_uncured_risk_col = uncured_risk_col,
      uncured_mean_support_flag = dplyr::case_when(
        !is.na(MSTu) & dplyr::coalesce(stage6__followup_not_contradicted_flag, FALSE) ~ "reported",
        !is.na(MSTu) ~ "available_but_not_supported",
        (!is.na(uncured_survival) | !is.na(uncured_risk)) ~ "annual_grid_only",
        TRUE ~ "not_available_from_saved_stage7_outputs"
      )
    ) %>%
    select(names(output_schema))
  
  bind_rows(output_schema[0, ], out)
}

# 🔴 Define: plot builders and writers ===============================
make_plot_group_label <- function(dataset, formula_variant, formula_label = NA_character_) {
  lbl <- ifelse(!is.na(formula_label) & nzchar(formula_label), formula_label, formula_variant)
  paste(dataset, lbl, sep = " | ")
}

build_stage7_visual_plot_registry <- function(risk_summary_tbl, delta_risk_tbl, threshold_metrics_tbl, common_horizons_year, threshold_plot_horizons_year) {
  plots <- list()
  
  if (nrow(risk_summary_tbl) > 0L && all(c("mean_risk_overall", "mean_risk_ci_lower", "mean_risk_ci_upper", "horizon_year", "model_id", "dataset", "formula_variant") %in% names(risk_summary_tbl))) {
    plot_risk_tbl <- risk_summary_tbl %>%
      mutate(plot_group = make_plot_group_label(dataset, formula_variant, formula_label))
    
    if (nrow(plot_risk_tbl) > 0L) {
      plots[["stage7_visual_overall_risk"]] <- list(
        source_table = "stage7_risk_summary.csv",
        description = "Overall predicted transition risk by horizon.",
        plot = ggplot(plot_risk_tbl, aes(x = horizon_year, y = mean_risk_overall, color = model_id, group = model_id)) +
          geom_line(linewidth = 0.8, na.rm = TRUE) +
          geom_point(size = 1.2, na.rm = TRUE) +
          geom_ribbon(
            aes(ymin = mean_risk_ci_lower, ymax = mean_risk_ci_upper, fill = model_id),
            alpha = 0.12,
            linewidth = 0,
            inherit.aes = TRUE,
            na.rm = TRUE,
            show.legend = FALSE
          ) +
          facet_wrap(~ plot_group, scales = "free_y") +
          scale_x_continuous(breaks = common_horizons_year) +
          labs(
            title = "Stage 7: overall predicted transition risk by horizon",
            x = "Years after cohort entry",
            y = "Mean predicted overall risk",
            color = "Model"
          ) +
          theme_bw()
      )
    }
  }
  
  if (nrow(delta_risk_tbl) > 0L && all(c("mean_delta_risk_nc_minus_cure", "mean_delta_risk_ci_lower", "mean_delta_risk_ci_upper", "horizon_year", "model_id", "dataset", "formula_variant") %in% names(delta_risk_tbl))) {
    plot_delta_tbl <- delta_risk_tbl %>%
      mutate(plot_group = make_plot_group_label(dataset, formula_variant, formula_label))
    
    if (nrow(plot_delta_tbl) > 0L) {
      plots[["stage7_visual_delta_risk"]] <- list(
        source_table = "stage7_delta_risk.csv",
        description = "DeltaRisk_NC_C(t) across horizons.",
        plot = ggplot(plot_delta_tbl, aes(x = horizon_year, y = mean_delta_risk_nc_minus_cure, color = model_id, group = model_id)) +
          geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
          geom_line(linewidth = 0.8, na.rm = TRUE) +
          geom_point(size = 1.2, na.rm = TRUE) +
          geom_ribbon(
            aes(ymin = mean_delta_risk_ci_lower, ymax = mean_delta_risk_ci_upper, fill = model_id),
            alpha = 0.12,
            linewidth = 0,
            inherit.aes = TRUE,
            na.rm = TRUE,
            show.legend = FALSE
          ) +
          facet_wrap(~ plot_group, scales = "free_y") +
          scale_x_continuous(breaks = common_horizons_year) +
          labs(
            title = "Stage 7: DeltaRisk_NC_C(t) = risk(no-cure) - risk(cure)",
            x = "Years after cohort entry",
            y = "Mean delta risk (no-cure minus cure)",
            color = "Cure model"
          ) +
          theme_bw()
      )
    }
  }
  
  if (nrow(threshold_metrics_tbl) > 0L && all(c("false_positive_burden_primary", "false_positive_burden_primary_ci_lower", "false_positive_burden_primary_ci_upper", "net_benefit", "net_benefit_ci_lower", "net_benefit_ci_upper", "threshold", "horizon_year", "model_id", "dataset", "formula_variant") %in% names(threshold_metrics_tbl))) {
    plot_threshold_tbl <- threshold_metrics_tbl %>%
      filter(horizon_year %in% threshold_plot_horizons_year) %>%
      filter(!threshold_projection_suppressed_flag) %>%
      mutate(plot_group = make_plot_group_label(dataset, formula_variant, formula_label))
    
    if (nrow(plot_threshold_tbl) > 0L) {
      plots[["stage7_visual_false_positive_burden"]] <- list(
        source_table = "stage7_threshold_metrics.csv",
        description = "False-positive burden among non-events across thresholds.",
        plot = ggplot(plot_threshold_tbl, aes(x = threshold, y = false_positive_burden_primary, color = model_id, group = model_id)) +
          geom_line(linewidth = 0.8, na.rm = TRUE) +
          geom_point(size = 1.2, na.rm = TRUE) +
          geom_ribbon(
            aes(ymin = false_positive_burden_primary_ci_lower, ymax = false_positive_burden_primary_ci_upper, fill = model_id),
            alpha = 0.12,
            linewidth = 0,
            inherit.aes = TRUE,
            na.rm = TRUE,
            show.legend = FALSE
          ) +
          facet_grid(horizon_year ~ plot_group, scales = "free_y") +
          labs(
            title = "Stage 7: false-positive burden among non-events",
            x = "Risk threshold",
            y = "False-positive burden among non-events (IPCW-weighted)",
            color = "Model"
          ) +
          theme_bw()
      )
      
      plots[["stage7_visual_net_benefit"]] <- list(
        source_table = "stage7_threshold_metrics.csv",
        description = "Net benefit across thresholds.",
        plot = ggplot(plot_threshold_tbl, aes(x = threshold, y = net_benefit, color = model_id, group = model_id)) +
          geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
          geom_line(linewidth = 0.8, na.rm = TRUE) +
          geom_point(size = 1.2, na.rm = TRUE) +
          geom_ribbon(
            aes(ymin = net_benefit_ci_lower, ymax = net_benefit_ci_upper, fill = model_id),
            alpha = 0.12,
            linewidth = 0,
            inherit.aes = TRUE,
            na.rm = TRUE,
            show.legend = FALSE
          ) +
          facet_grid(horizon_year ~ plot_group, scales = "free_y") +
          labs(
            title = "Stage 7: net benefit across thresholds",
            x = "Risk threshold",
            y = "Net benefit (IPCW-weighted)",
            color = "Model"
          ) +
          theme_bw()
      )
    }
  }
  
  plots
}

save_plot_registry_outputs <- function(plot_registry, combined_pdf_file, force_rebuild_pdf = FALSE, write_pngs = TRUE, force_rebuild_pngs = FALSE, width = plot_width_in, height = plot_height_in, dpi = plot_dpi) {
  output_records <- tibble(
    plot_name = character(),
    source_table = character(),
    file_type = character(),
    file_name = character(),
    file_path = character(),
    description = character(),
    exists = logical()
  )
  
  if (length(plot_registry) == 0L) {
    return(output_records)
  }
  
  if (should_rebuild_output(combined_pdf_file, force_rebuild = force_rebuild_pdf)) {
    grDevices::pdf(combined_pdf_file, width = width, height = height)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    for (plot_name in names(plot_registry)) {
      print(plot_registry[[plot_name]]$plot)
    }
    grDevices::dev.off()
  }
  
  output_records <- bind_rows(
    output_records,
    tibble(
      plot_name = "stage7_visual_summary",
      source_table = paste(unique(vapply(plot_registry, `[[`, character(1), "source_table")), collapse = "|"),
      file_type = "pdf",
      file_name = basename(combined_pdf_file),
      file_path = normalize_existing_path(combined_pdf_file),
      description = "Combined Stage 7 visual summary PDF.",
      exists = file.exists(combined_pdf_file)
    )
  )
  
  if (isTRUE(write_pngs)) {
    for (plot_name in names(plot_registry)) {
      png_file <- make_output_path(paste0(plot_name, ".png"))
      if (should_rebuild_output(png_file, force_rebuild = force_rebuild_pngs)) {
        ggplot2::ggsave(filename = png_file, plot = plot_registry[[plot_name]]$plot, width = width, height = height, dpi = dpi, units = "in")
      }
      output_records <- bind_rows(
        output_records,
        tibble(
          plot_name = plot_name,
          source_table = plot_registry[[plot_name]]$source_table,
          file_type = "png",
          file_name = basename(png_file),
          file_path = normalize_existing_path(png_file),
          description = plot_registry[[plot_name]]$description,
          exists = file.exists(png_file)
        )
      )
    }
  }
  
  output_records
}

# 🔴 Define: bootstrap extraction readers ===============================
required_bootstrap_tables <- c(
  "fit_registry",
  "coefficients",
  "risk_summary",
  "delta_risk",
  "threshold_metrics"
)

standardize_bootstrap_results <- function(x) {
  if (!is.list(x)) {
    stop("`stage7_bootstrap_merged_results.rds` must contain a list.", call. = FALSE)
  }
  for (nm in required_bootstrap_tables) {
    if (!nm %in% names(x) || is.null(x[[nm]])) {
      x[[nm]] <- tibble()
    }
    x[[nm]] <- tibble::as_tibble(x[[nm]])
  }
  x[required_bootstrap_tables]
}

normalize_fit_registry <- function(fit_registry_raw) {
  if (nrow(fit_registry_raw) == 0L) {
    return(
      tibble(
        rep_id = integer(),
        dataset = character(),
        formula_variant = character(),
        site_branch = character(),
        interaction_branch = character(),
        model_id = character(),
        model_class = character(),
        model_block = character(),
        model_family = character(),
        latency_type = character(),
        matched_nocure_model_id = character(),
        fit_status = character(),
        prediction_status = character(),
        overall_status = character(),
        fit_success = logical(),
        prediction_success = logical(),
        overall_success = logical(),
        failure_type = character(),
        convergence_code = double(),
        warning_count_total = double(),
        primary_error_message = character(),
        primary_warning_message = character(),
        model_key = character(),
        failure_key = character()
      )
    )
  }
  
  fit_success_fallback <- is_success_status(col_or_chr(fit_registry_raw, "fit_status"))
  prediction_success_fallback <- is_success_status(col_or_chr(fit_registry_raw, "prediction_status"))
  
  fit_success_vec <- dplyr::coalesce(
    col_or_lgl(fit_registry_raw, "fit_component_success"),
    fit_success_fallback
  )
  
  prediction_success_vec <- dplyr::coalesce(
    col_or_lgl(fit_registry_raw, "prediction_component_success"),
    prediction_success_fallback
  )
  
  warning_count_total_vec <- dplyr::coalesce(
    col_or_num(fit_registry_raw, "warning_count_total"),
    col_or_num(fit_registry_raw, "fit_warning_count") + col_or_num(fit_registry_raw, "prediction_warning_count"),
    col_or_num(fit_registry_raw, "fit_warning_count"),
    col_or_num(fit_registry_raw, "prediction_warning_count"),
    rep(0, nrow(fit_registry_raw))
  )
  
  tibble(
    rep_id = as.integer(col_or_num(fit_registry_raw, "rep_id")),
    dataset = col_or_chr(fit_registry_raw, "dataset", "<missing_dataset>"),
    formula_variant = col_or_chr(fit_registry_raw, "formula_variant", "<missing_formula>"),
    site_branch = col_or_chr(fit_registry_raw, "site_branch"),
    interaction_branch = col_or_chr(fit_registry_raw, "interaction_branch"),
    model_id = col_or_chr(fit_registry_raw, "model_id", "<missing_model>"),
    model_class = col_or_chr(fit_registry_raw, "model_class"),
    model_block = col_or_chr(fit_registry_raw, "model_block"),
    model_family = col_or_chr(fit_registry_raw, "model_family"),
    latency_type = col_or_chr(fit_registry_raw, "latency_type"),
    matched_nocure_model_id = col_or_chr(fit_registry_raw, "matched_nocure_model_id"),
    fit_status = col_or_chr(fit_registry_raw, "fit_status"),
    prediction_status = col_or_chr(fit_registry_raw, "prediction_status"),
    overall_status = col_or_chr(fit_registry_raw, "overall_status"),
    fit_success = as.logical(fit_success_vec),
    prediction_success = as.logical(prediction_success_vec),
    convergence_code = col_or_num(fit_registry_raw, "convergence_code"),
    warning_count_total = as.numeric(warning_count_total_vec),
    primary_error_message = coalesce_character(
      col_or_chr(fit_registry_raw, "bootstrap_error"),
      col_or_chr(fit_registry_raw, "error_message"),
      col_or_chr(fit_registry_raw, "fit_error_message"),
      col_or_chr(fit_registry_raw, "prediction_error_message"),
      col_or_chr(fit_registry_raw, "fit_message"),
      col_or_chr(fit_registry_raw, "prediction_message"),
      col_or_chr(fit_registry_raw, "message")
    ),
    primary_warning_message = coalesce_character(
      col_or_chr(fit_registry_raw, "warning_message"),
      col_or_chr(fit_registry_raw, "fit_warning_message"),
      col_or_chr(fit_registry_raw, "prediction_warning_message")
    )
  ) %>%
    mutate(
      overall_success = fit_success & prediction_success,
      failure_type = dplyr::case_when(
        !fit_success & !prediction_success ~ "fit_and_prediction_failed",
        !fit_success ~ "fit_failed",
        fit_success & !prediction_success ~ "prediction_failed",
        TRUE ~ "success"
      ),
      model_key = paste(dataset, formula_variant, model_id, sep = "__"),
      failure_key = ifelse(
        failure_type == "success",
        NA_character_,
        paste(dataset, formula_variant, model_id, failure_type, sep = "::")
      )
    )
}

make_failure_events <- function(fit_registry_norm) {
  fit_registry_norm %>%
    filter(!overall_success) %>%
    transmute(
      rep_id,
      dataset,
      formula_variant,
      site_branch,
      interaction_branch,
      model_id,
      model_class,
      model_block,
      model_family,
      latency_type,
      matched_nocure_model_id,
      failure_type,
      fit_status,
      prediction_status,
      overall_status,
      convergence_code,
      warning_count_total,
      primary_error_message,
      primary_warning_message
    ) %>%
    arrange(rep_id, dataset, formula_variant, model_id)
}

make_replicate_overview <- function(fit_registry_norm) {
  fit_registry_norm %>%
    group_by(rep_id) %>%
    summarise(
      n_model_rows = n(),
      n_fail_total = sum(!overall_success, na.rm = TRUE),
      n_fit_fail = sum(!fit_success, na.rm = TRUE),
      n_prediction_fail = sum(fit_success & !prediction_success, na.rm = TRUE),
      n_warning_rows = sum(warning_count_total > 0, na.rm = TRUE),
      fit_success_rate = mean(fit_success, na.rm = TRUE),
      prediction_success_rate = mean(prediction_success, na.rm = TRUE),
      overall_success_rate = mean(overall_success, na.rm = TRUE),
      failed_dataset_set = collapse_unique_sorted(dataset[!overall_success], empty_value = ""),
      failed_model_set = collapse_unique_sorted(model_id[!overall_success], empty_value = ""),
      failure_signature = ifelse(
        sum(!overall_success, na.rm = TRUE) == 0,
        "NO_FAILURE",
        paste(sort(unique(failure_key[!is.na(failure_key)])), collapse = " | ")
      ),
      pnu_fail_n = sum(dataset == "PNU" & !overall_success, na.rm = TRUE),
      snu_fail_n = sum(dataset == "SNU" & !overall_success, na.rm = TRUE),
      merged_fail_n = sum(dataset == "merged" & !overall_success, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(rep_id)
}

make_failure_pattern_summary <- function(replicate_overview) {
  replicate_overview %>%
    group_by(failure_signature) %>%
    summarise(
      n_replicates = n(),
      n_replicates_with_failures = sum(n_fail_total > 0, na.rm = TRUE),
      mean_failures_per_replicate = mean(n_fail_total, na.rm = TRUE),
      min_rep_id = min(rep_id, na.rm = TRUE),
      max_rep_id = max(rep_id, na.rm = TRUE),
      example_rep_ids = paste(head(rep_id, 20), collapse = "|"),
      .groups = "drop"
    ) %>%
    arrange(desc(n_replicates), desc(mean_failures_per_replicate), failure_signature)
}

make_failure_by_model <- function(failure_events) {
  if (nrow(failure_events) == 0L) {
    return(
      tibble(
        dataset = character(),
        formula_variant = character(),
        model_id = character(),
        failure_type = character(),
        n_failure_rows = integer(),
        n_replicates = integer(),
        example_rep_ids = character()
      )
    )
  }
  
  failure_events %>%
    group_by(dataset, formula_variant, model_id, failure_type) %>%
    summarise(
      n_failure_rows = n(),
      n_replicates = n_distinct(rep_id),
      example_rep_ids = paste(head(sort(unique(rep_id)), 20), collapse = "|"),
      .groups = "drop"
    ) %>%
    arrange(desc(n_replicates), desc(n_failure_rows), dataset, formula_variant, model_id)
}

metric_id_cols <- c(
  "dataset", "formula_variant", "site_branch", "interaction_branch",
  "model_id", "model_class", "model_block", "model_family", "latency_type",
  "matched_nocure_model_id", "component", "term", "horizon_year", "threshold"
)

numeric_metric_cols <- function(df, structural_cols = character(), extra_exclude = character()) {
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  out <- setdiff(numeric_cols, structural_cols)
  out <- setdiff(out, extra_exclude)
  out
}

complete_metric_id_columns <- function(df) {
  default_values <- list(
    dataset = NA_character_,
    formula_variant = NA_character_,
    site_branch = NA_character_,
    interaction_branch = NA_character_,
    model_id = NA_character_,
    model_class = NA_character_,
    model_block = NA_character_,
    model_family = NA_character_,
    latency_type = NA_character_,
    matched_nocure_model_id = NA_character_,
    component = NA_character_,
    term = NA_character_,
    horizon_year = NA_integer_,
    threshold = NA_real_
  )
  
  for (nm in names(default_values)) {
    if (!nm %in% names(df)) {
      df[[nm]] <- default_values[[nm]]
    }
  }
  
  df
}

get_metric_specs <- function(bootstrap_results) {
  list(
    coefficients = list(
      df = bootstrap_results$coefficients,
      structural_cols = c("rep_id"),
      metric_cols = intersect("estimate", names(bootstrap_results$coefficients))
    ),
    risk_summary = list(
      df = bootstrap_results$risk_summary,
      structural_cols = c("rep_id", "horizon_year"),
      metric_cols = numeric_metric_cols(
        bootstrap_results$risk_summary,
        structural_cols = c("rep_id", "horizon_year"),
        extra_exclude = c("n_subject")
      )
    ),
    delta_risk = list(
      df = bootstrap_results$delta_risk,
      structural_cols = c("rep_id", "horizon_year"),
      metric_cols = numeric_metric_cols(
        bootstrap_results$delta_risk,
        structural_cols = c("rep_id", "horizon_year"),
        extra_exclude = c("n_subject")
      )
    ),
    threshold_metrics = list(
      df = bootstrap_results$threshold_metrics,
      structural_cols = c("rep_id", "horizon_year", "threshold"),
      metric_cols = numeric_metric_cols(
        bootstrap_results$threshold_metrics,
        structural_cols = c("rep_id", "horizon_year", "threshold"),
        extra_exclude = c("n_subject")
      )
    )
  )
}

summarize_one_metric_column <- function(df, table_name, metric_col) {
  if (nrow(df) == 0L || !metric_col %in% names(df) || !"rep_id" %in% names(df)) {
    return(tibble())
  }
  
  group_cols <- intersect(metric_id_cols, names(df))
  
  work_df <- df %>%
    mutate(
      metric_value = suppressWarnings(as.numeric(.data[[metric_col]])),
      rep_id = as.integer(rep_id)
    ) %>%
    select(all_of(c("rep_id", group_cols, "metric_value")))
  
  if (all(is.na(work_df$metric_value))) {
    return(tibble())
  }
  
  work_df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      table_name = table_name,
      metric_name = metric_col,
      n_replicates_nonmissing = n_distinct(rep_id[!is.na(metric_value)]),
      mean_value = mean(metric_value, na.rm = TRUE),
      sd_value = stats::sd(metric_value, na.rm = TRUE),
      min_value = suppressWarnings(min(metric_value, na.rm = TRUE)),
      q25_value = safe_quantile(metric_value, 0.25),
      median_value = safe_quantile(metric_value, 0.50),
      q75_value = safe_quantile(metric_value, 0.75),
      max_value = suppressWarnings(max(metric_value, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      min_value = ifelse(is.infinite(min_value), NA_real_, min_value),
      max_value = ifelse(is.infinite(max_value), NA_real_, max_value)
    )
}

summarize_metric_table <- function(df, table_name, metric_cols) {
  if (nrow(df) == 0L || length(metric_cols) == 0L) {
    return(tibble())
  }
  
  bind_rows(lapply(metric_cols, function(metric_col) {
    summarize_one_metric_column(df, table_name, metric_col)
  }))
}

extract_top_metric_values <- function(top_cells, bootstrap_results) {
  if (nrow(top_cells) == 0L) {
    return(
      tibble(
        metric_cell_id = integer(),
        table_name = character(),
        metric_name = character(),
        rep_id = integer(),
        dataset = character(),
        formula_variant = character(),
        site_branch = character(),
        interaction_branch = character(),
        model_id = character(),
        model_class = character(),
        model_block = character(),
        model_family = character(),
        latency_type = character(),
        matched_nocure_model_id = character(),
        component = character(),
        term = character(),
        horizon_year = integer(),
        threshold = double(),
        metric_value = double(),
        metric_key = character()
      )
    )
  }
  
  out <- list()
  
  for (table_name in unique(top_cells$table_name)) {
    df <- bootstrap_results[[table_name]]
    if (is.null(df) || nrow(df) == 0L) next
    
    table_cells <- top_cells %>% filter(table_name == !!table_name)
    
    for (metric_name in unique(table_cells$metric_name)) {
      if (!metric_name %in% names(df)) next
      
      cells_sub <- table_cells %>% filter(metric_name == !!metric_name)
      join_cols <- intersect(metric_id_cols, names(df))
      
      src <- df %>%
        mutate(
          rep_id = as.integer(rep_id),
          metric_value = suppressWarnings(as.numeric(.data[[metric_name]]))
        ) %>%
        select(all_of(c("rep_id", join_cols, "metric_value")))
      
      cell_map <- cells_sub %>%
        select(metric_cell_id, metric_key, all_of(join_cols)) %>%
        distinct()
      
      joined <- dplyr::inner_join(
        src,
        cell_map,
        by = join_cols,
        na_matches = "na"
      ) %>%
        mutate(
          table_name = table_name,
          metric_name = metric_name
        ) %>%
        complete_metric_id_columns() %>%
        select(
          metric_cell_id, table_name, metric_name, rep_id,
          dataset, formula_variant, site_branch, interaction_branch,
          model_id, model_class, model_block, model_family, latency_type,
          matched_nocure_model_id, component, term, horizon_year, threshold,
          metric_value, metric_key
        )
      
      out[[paste(table_name, metric_name, sep = "__")]] <- joined
    }
  }
  
  bind_rows(out) %>% arrange(metric_cell_id, rep_id)
}

write_metric_values_long_streaming <- function(bootstrap_results, out_file) {
  if (file.exists(out_file)) file.remove(out_file)
  metric_specs <- get_metric_specs(bootstrap_results)
  
  for (table_name in names(metric_specs)) {
    spec <- metric_specs[[table_name]]
    df <- spec$df
    metric_cols <- spec$metric_cols
    
    if (nrow(df) == 0L || length(metric_cols) == 0L || !"rep_id" %in% names(df)) next
    
    for (metric_name in metric_cols) {
      if (!metric_name %in% names(df)) next
      
      chunk_df <- df %>%
        transmute(
          table_name = table_name,
          metric_name = metric_name,
          rep_id = as.integer(rep_id),
          dataset = if ("dataset" %in% names(df)) as.character(dataset) else NA_character_,
          formula_variant = if ("formula_variant" %in% names(df)) as.character(formula_variant) else NA_character_,
          site_branch = if ("site_branch" %in% names(df)) as.character(site_branch) else NA_character_,
          interaction_branch = if ("interaction_branch" %in% names(df)) as.character(interaction_branch) else NA_character_,
          model_id = if ("model_id" %in% names(df)) as.character(model_id) else NA_character_,
          model_class = if ("model_class" %in% names(df)) as.character(model_class) else NA_character_,
          model_block = if ("model_block" %in% names(df)) as.character(model_block) else NA_character_,
          model_family = if ("model_family" %in% names(df)) as.character(model_family) else NA_character_,
          latency_type = if ("latency_type" %in% names(df)) as.character(latency_type) else NA_character_,
          matched_nocure_model_id = if ("matched_nocure_model_id" %in% names(df)) as.character(matched_nocure_model_id) else NA_character_,
          component = if ("component" %in% names(df)) as.character(component) else NA_character_,
          term = if ("term" %in% names(df)) as.character(term) else NA_character_,
          horizon_year = if ("horizon_year" %in% names(df)) as.integer(horizon_year) else NA_integer_,
          threshold = if ("threshold" %in% names(df)) as.numeric(threshold) else NA_real_,
          metric_value = suppressWarnings(as.numeric(.data[[metric_name]]))
        )
      
      readr::write_csv(chunk_df, out_file, append = file.exists(out_file))
    }
  }
  
  normalize_existing_path(out_file)
}

build_bootstrap_plot_registry <- function(replicate_overview, failure_pattern_summary, failure_by_model, top_metric_values_long) {
  plots <- list()
  
  if (nrow(replicate_overview) > 0L) {
    plots[["stage7_bootstrap_failure_histogram"]] <- list(
      source_table = "stage7_replicate_overview.csv",
      description = "Bootstrap replicate failure count distribution.",
      plot = ggplot(replicate_overview, aes(x = n_fail_total)) +
        geom_histogram(binwidth = 1, boundary = -0.5, closed = "right") +
        labs(
          title = "Bootstrap replicate failure count distribution",
          x = "Failed model rows per replicate",
          y = "Replicate count"
        ) +
        theme_minimal()
    )
    
    plots[["stage7_bootstrap_failure_by_replicate"]] <- list(
      source_table = "stage7_replicate_overview.csv",
      description = "Failure count by bootstrap replicate.",
      plot = ggplot(replicate_overview, aes(x = rep_id, y = n_fail_total)) +
        geom_point(alpha = 0.7) +
        geom_line(alpha = 0.35) +
        labs(
          title = "Failure count by bootstrap replicate",
          x = "Replicate ID",
          y = "Failed model rows"
        ) +
        theme_minimal()
    )
  }
  
  if (nrow(failure_pattern_summary) > 0L) {
    plot_df <- failure_pattern_summary %>%
      slice_head(n = top_n_failure_patterns) %>%
      mutate(pattern_label = stringr::str_wrap(failure_signature, width = 80))
    
    plots[["stage7_bootstrap_failure_patterns"]] <- list(
      source_table = "stage7_replicate_failure_pattern_summary.csv",
      description = sprintf("Top %s bootstrap failure signatures.", top_n_failure_patterns),
      plot = ggplot(plot_df, aes(x = reorder(pattern_label, n_replicates), y = n_replicates)) +
        geom_col() +
        coord_flip() +
        labs(
          title = sprintf("Top %s failure signatures", top_n_failure_patterns),
          x = "Failure signature",
          y = "Replicate count"
        ) +
        theme_minimal()
    )
  }
  
  if (nrow(failure_by_model) > 0L) {
    plot_df <- failure_by_model %>%
      slice_head(n = top_n_failure_models) %>%
      mutate(model_label = stringr::str_wrap(paste(dataset, formula_variant, model_id, failure_type, sep = " | "), width = 70))
    
    plots[["stage7_bootstrap_failure_models"]] <- list(
      source_table = "stage7_replicate_failure_by_model.csv",
      description = sprintf("Top %s failing model patterns.", top_n_failure_models),
      plot = ggplot(plot_df, aes(x = reorder(model_label, n_replicates), y = n_replicates)) +
        geom_col() +
        coord_flip() +
        labs(
          title = sprintf("Top %s failing model patterns", top_n_failure_models),
          x = "Dataset | formula | model | failure type",
          y = "Replicate count"
        ) +
        theme_minimal()
    )
  }
  
  if (nrow(top_metric_values_long) > 0L) {
    plot_df <- top_metric_values_long %>%
      mutate(metric_key_wrapped = stringr::str_wrap(metric_key, width = 70))
    
    plots[["stage7_bootstrap_top_variable_metrics"]] <- list(
      source_table = "stage7_top_variable_metric_values_long.csv",
      description = sprintf("Top %s most variable metric cells.", top_n_variable_metric_cells),
      plot = ggplot(plot_df, aes(x = metric_key_wrapped, y = metric_value)) +
        geom_boxplot(outlier.alpha = 0.25) +
        coord_flip() +
        labs(
          title = sprintf("Top %s most variable metric cells", top_n_variable_metric_cells),
          x = "Metric cell",
          y = "Metric value across replicates"
        ) +
        theme_minimal()
    )
  }
  
  plots
}

# 🔴 Define: full Stage 7 core fitting pipeline ===============================
stage7_core_required_files <- c(
  "stage7_fit_registry.csv",
  "stage7_coefficients.csv",
  "stage7_subject_predictions.csv",
  "stage7_risk_summary.csv",
  "stage7_delta_risk.csv",
  "stage7_threshold_metrics.csv"
)

stage7_core_outputs_exist <- function(root_dir) {
  root_dir <- normalize_existing_path(root_dir)
  all(file.exists(file.path(root_dir, stage7_core_required_files)))
}

resolve_stage7_run_root_or_create <- function(root_dir) {
  root_dir <- normalizePath(root_dir, winslash = "/", mustWork = FALSE)
  
  if (dir.exists(root_dir) && stage7_core_outputs_exist(root_dir)) {
    return(root_dir)
  }
  
  if (dir.exists(root_dir)) {
    subdirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
    candidate_dirs <- subdirs[vapply(subdirs, stage7_core_outputs_exist, logical(1))]
    if (length(candidate_dirs) == 1L) {
      message(sprintf("ℹ️  Stage 7 run root auto-resolved to existing nested run directory: %s", candidate_dirs[[1]]))
      return(normalize_existing_path(candidate_dirs[[1]]))
    }
    if (length(candidate_dirs) > 1L) {
      stop(
        paste0(
          "Multiple existing Stage 7 run directories were found under the supplied run_root_dir. ",
          "Please point run_root_dir to the exact folder you want to reuse or overwrite.\nCandidates:\n- ",
          paste(candidate_dirs, collapse = "\n- ")
        ),
        call. = FALSE
      )
    }
  }
  
  dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)
  normalize_existing_path(root_dir)
}

stage7_capture_with_warnings <- function(expr) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  list(
    value = if (inherits(value, "error")) NULL else value,
    warnings = warnings,
    error_message = if (inherits(value, "error")) conditionMessage(value) else NA_character_
  )
}

stage7_trim_rhs <- function(rhs) {
  rhs <- stringr::str_squish(as.character(rhs))
  rhs[is.na(rhs)] <- NA_character_
  rhs
}

stage7_remove_site_term_rhs <- function(rhs) {
  rhs <- stage7_trim_rhs(rhs)
  rhs <- stringr::str_replace_all(rhs, "\\s*\\+\\s*site\\b", "")
  rhs <- stringr::str_replace_all(rhs, "\\bsite\\b\\s*\\+\\s*", "")
  rhs <- stringr::str_replace_all(rhs, "\\s*\\+\\s*\\+\\s*", " + ")
  rhs <- stringr::str_replace_all(rhs, "^\\s*\\+\\s*", "")
  rhs <- stringr::str_replace_all(rhs, "\\s*\\+\\s*$", "")
  stringr::str_squish(rhs)
}

stage7_has_site_term <- function(rhs) {
  grepl("\\bsite\\b", stage7_trim_rhs(rhs))
}

stage7_has_age_sex_interaction <- function(rhs) {
  grepl("age_s:sex_num", stage7_trim_rhs(rhs), fixed = TRUE)
}

stage7_formula_label_from_variant <- function(formula_variant) {
  dplyr::case_when(
    formula_variant %in% c("base", "base_sitefree") ~ "Base",
    formula_variant %in% c("interaction", "interaction_sitefree") ~ "Interaction",
    formula_variant == "base_siteadjusted" ~ "Site-adjusted",
    formula_variant == "interaction_siteadjusted" ~ "Site + interaction",
    formula_variant == "benchmark" ~ "Benchmark KM",
    TRUE ~ as.character(formula_variant)
  )
}

stage7_prepare_formula_plan <- function(stage1_inputs) {
  formula_tbl <- stage1_inputs$formula_registry %>%
    mutate(
      formula_variant = dplyr::case_when(
        dataset != "merged" & formula_name == "base" ~ "base",
        dataset != "merged" & formula_name == "interaction" ~ "interaction",
        dataset == "merged" & formula_name == "base" ~ "base_sitefree",
        dataset == "merged" & formula_name == "interaction" ~ "interaction_sitefree",
        dataset == "merged" & formula_name == "site_added" ~ "base_siteadjusted",
        dataset == "merged" & formula_name == "site_interaction" ~ "interaction_siteadjusted",
        TRUE ~ as.character(formula_name)
      ),
      incidence_rhs = dplyr::case_when(
        dataset == "merged" & formula_name %in% c("site_added", "site_interaction") ~ stage7_remove_site_term_rhs(formula_rhs),
        TRUE ~ stage7_trim_rhs(formula_rhs)
      ),
      latency_rhs = stage7_trim_rhs(formula_rhs),
      incidence_uses_site = stage7_has_site_term(incidence_rhs),
      latency_uses_site = stage7_has_site_term(latency_rhs),
      uses_site = incidence_uses_site | latency_uses_site,
      uses_age_sex_interaction = stage7_has_age_sex_interaction(incidence_rhs) | stage7_has_age_sex_interaction(latency_rhs),
      site_branch = if_else(latency_uses_site, "site_adjusted", "site_free"),
      interaction_branch = if_else(uses_age_sex_interaction, "age_sex_interaction", "no_age_sex_interaction"),
      formula_label = stage7_formula_label_from_variant(formula_variant),
      formula_scope = "main_transition_only_scale",
      site_term_interpretation = if_else(dataset == "merged" & uses_site, "structural_context_proxy_not_causal_treatment_effect", "not_applicable"),
      site_placement_label = derive_site_placement_label(
        dataset = dataset,
        site_branch = site_branch,
        existing_label = NA_character_,
        incidence_uses_site = incidence_uses_site,
        latency_uses_site = latency_uses_site
      ),
      stage = "Stage 7",
      branch = "Stage7",
      incidence_link = "logit"
    ) %>%
    filter(
      (dataset != "merged" & formula_variant %in% c("base", "interaction")) |
        (dataset == "merged" & formula_variant %in% c("base_sitefree", "base_siteadjusted", "interaction_sitefree", "interaction_siteadjusted"))
    ) %>%
    select(
      dataset, formula_variant, formula_label,
      site_branch, interaction_branch,
      incidence_rhs, latency_rhs,
      formula_scope, site_term_interpretation,
      uses_site, uses_age_sex_interaction,
      incidence_uses_site, latency_uses_site,
      site_placement_label, stage, branch, incidence_link
    )
  
  benchmark_tbl <- tibble(
    dataset = c("PNU", "SNU", "merged"),
    formula_variant = "benchmark",
    formula_label = "Benchmark KM",
    site_branch = "benchmark",
    interaction_branch = "benchmark",
    incidence_rhs = NA_character_,
    latency_rhs = NA_character_,
    formula_scope = "benchmark_only",
    site_term_interpretation = "not_applicable",
    uses_site = FALSE,
    uses_age_sex_interaction = FALSE,
    incidence_uses_site = FALSE,
    latency_uses_site = FALSE,
    site_placement_label = c("not_applicable_single_site", "not_applicable_single_site", "site_in_neither"),
    stage = "Stage 7",
    branch = "Stage7",
    incidence_link = NA_character_
  )
  
  list(model_plan = formula_tbl, benchmark_plan = benchmark_tbl)
}

stage7_model_catalog <- function() {
  tibble::tribble(
    ~model_id, ~fit_method, ~model_class, ~model_block, ~model_family, ~latency_type, ~fit_engine, ~is_cure_model, ~matched_nocure_model_id,
    "benchmark_km", "benchmark_km", "km_benchmark", "benchmark", "km", "benchmark", "survival::survfit", FALSE, NA_character_,
    "nocure_exp", "nocure_exp", "no_cure", "no_cure", "exp", "parametric_aft", "survival::survreg", FALSE, NA_character_,
    "nocure_weibull", "nocure_weibull", "no_cure", "no_cure", "weibull", "parametric_aft", "survival::survreg", FALSE, NA_character_,
    "nocure_lnorm", "nocure_lnorm", "no_cure", "no_cure", "lnorm", "parametric_aft", "survival::survreg", FALSE, NA_character_,
    "nocure_llogis", "nocure_llogis", "no_cure", "no_cure", "llogis", "parametric_aft", "survival::survreg", FALSE, NA_character_,
    "nocure_coxph", "nocure_coxph", "no_cure", "no_cure", "coxph", "cox_ph", "survival::coxph", FALSE, NA_character_,
    "cure_exp", "cure_exp", "parametric_mixture_cure", "frequentist_cure_main", "exp", "parametric_aft", "custom_mle", TRUE, "nocure_exp",
    "cure_weibull", "cure_weibull", "parametric_mixture_cure", "frequentist_cure_main", "weibull", "parametric_aft", "custom_mle", TRUE, "nocure_weibull",
    "cure_lnorm", "cure_lnorm", "parametric_mixture_cure", "frequentist_cure_main", "lnorm", "parametric_aft", "custom_mle", TRUE, "nocure_lnorm",
    "cure_llogis", "cure_llogis", "parametric_mixture_cure", "frequentist_cure_main", "llogis", "parametric_aft", "custom_mle", TRUE, "nocure_llogis",
    "cure_coxlatency", "cure_coxlatency", "semiparametric_cox_latency_cure", "frequentist_cure_main", "coxph", "cox_ph", "smcure", TRUE, "nocure_coxph",
    "cure_aft_sensitivity", "cure_aft_sensitivity", "aft_latency_cure_sensitivity", "frequentist_cure_sensitivity", "aft", "aft_semiparametric", "smcure", TRUE, NA_character_
  )
}

stage7_make_surv_formula <- function(rhs) {
  rhs <- stage7_trim_rhs(rhs)
  if (is.na(rhs) || rhs == "") {
    as.formula("survival::Surv(time_year, event_main) ~ 1")
  } else {
    as.formula(paste("survival::Surv(time_year, event_main) ~", rhs))
  }
}

stage7_make_rhs_formula <- function(rhs) {
  rhs <- stage7_trim_rhs(rhs)
  if (is.na(rhs) || rhs == "") {
    as.formula("~ 1")
  } else {
    as.formula(paste("~", rhs))
  }
}

stage7_model_matrix <- function(df, rhs, drop_intercept = FALSE) {
  mm <- stats::model.matrix(stage7_make_rhs_formula(rhs), data = df)
  if (isTRUE(drop_intercept) && "(Intercept)" %in% colnames(mm)) {
    mm <- mm[, setdiff(colnames(mm), "(Intercept)"), drop = FALSE]
  }
  mm
}

stage7_align_mm_to_coef <- function(mm, coef_names) {
  mm <- as.data.frame(as.matrix(mm))
  coef_names <- as.character(coef_names)
  
  out <- sapply(coef_names, function(nm) {
    if (identical(nm, "(Intercept)")) {
      rep(1, nrow(mm))
    } else if (nm %in% names(mm)) {
      mm[[nm]]
    } else {
      rep(0, nrow(mm))
    }
  })
  
  out <- as.matrix(out)
  if (is.null(dim(out))) out <- matrix(out, nrow = nrow(mm))
  colnames(out) <- coef_names
  out
}

stage7_standardize_term_name <- function(x) {
  x <- as.character(x)
  out <- vapply(strsplit(x, ":", fixed = TRUE), function(parts) {
    parts <- vapply(parts, function(part) {
      if (identical(part, "(Intercept)")) return(part)
      if (grepl("^site", part)) {
        suffix <- sub("^site", "", part)
        suffix <- gsub("[^A-Za-z0-9]+", "_", suffix)
        suffix <- tolower(gsub("^_+|_+$", "", suffix))
        return(paste0("site_", suffix))
      }
      part
    }, character(1))
    paste(parts, collapse = ":")
  }, character(1))
  out
}

stage7_step_eval <- function(time, value, at, yleft = 1, yright = NULL) {
  time <- as.numeric(time)
  value <- as.numeric(value)
  at <- as.numeric(at)
  keep <- is.finite(time) & is.finite(value)
  time <- time[keep]
  value <- value[keep]
  if (length(time) == 0L) return(rep(yleft, length(at)))
  ord <- order(time)
  time <- time[ord]
  value <- value[ord]
  if (is.null(yright)) yright <- tail(value, 1)
  approx(x = time, y = value, xout = at, method = "constant", f = 0, yleft = yleft, yright = yright, ties = "ordered")$y
}

stage7_clamp_prob <- function(x, eps = 1e-12) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

stage7_parametric_surv_density <- function(time, lp, family, aux_log = NA_real_) {
  time <- pmax(as.numeric(time), 1e-10)
  lp <- as.numeric(lp)
  if (family == "exp") {
    scale_time <- exp(lp)
    surv <- exp(-time / scale_time)
    dens <- (1 / scale_time) * surv
  } else if (family == "weibull") {
    sigma <- exp(aux_log)
    shape <- 1 / sigma
    scale_time <- exp(lp)
    z <- (time / scale_time)^shape
    surv <- exp(-z)
    dens <- (shape / scale_time) * (time / scale_time)^(shape - 1) * surv
  } else if (family == "lnorm") {
    sigma <- exp(aux_log)
    z <- (log(time) - lp) / sigma
    surv <- 1 - stats::pnorm(z)
    dens <- stats::dnorm(z) / (sigma * time)
  } else if (family == "llogis") {
    sigma <- exp(aux_log)
    z <- (log(time) - lp) / sigma
    surv <- 1 - stats::plogis(z)
    dens <- (stats::plogis(z) * (1 - stats::plogis(z))) / (sigma * time)
  } else {
    stop(sprintf("Unsupported parametric family: %s", family), call. = FALSE)
  }
  list(surv = pmin(pmax(surv, 0), 1), dens = pmax(dens, 1e-300))
}

stage7_build_prediction_rows <- function(df, horizons, overall_surv_mat, susceptible_surv_mat, uncured_prob_vec, meta_row) {
  sex_label <- if ("sex_label" %in% names(df)) {
    as.character(df$sex_label)
  } else {
    dplyr::case_when(
      suppressWarnings(as.numeric(df$sex_num)) == 0 ~ "Male",
      suppressWarnings(as.numeric(df$sex_num)) == 1 ~ "Female",
      TRUE ~ NA_character_
    )
  }
  
  bind_rows(lapply(seq_along(horizons), function(j) {
    tibble(
      unique_person_id = as.character(df$unique_person_id),
      site = as.character(df$site),
      sex_num = suppressWarnings(as.numeric(df$sex_num)),
      sex_label = sex_label,
      horizon_year = as.integer(horizons[[j]]),
      predicted_survival_overall = as.numeric(overall_surv_mat[, j]),
      predicted_risk_overall = 1 - as.numeric(overall_surv_mat[, j]),
      predicted_survival_susceptible = as.numeric(susceptible_surv_mat[, j]),
      predicted_risk_susceptible = 1 - as.numeric(susceptible_surv_mat[, j]),
      predicted_uncured_fraction = as.numeric(uncured_prob_vec),
      predicted_cure_fraction = 1 - as.numeric(uncured_prob_vec),
      dataset = meta_row$dataset,
      formula_variant = meta_row$formula_variant,
      site_branch = meta_row$site_branch,
      interaction_branch = meta_row$interaction_branch,
      incidence_rhs = meta_row$incidence_rhs,
      latency_rhs = meta_row$latency_rhs,
      formula_label = meta_row$formula_label,
      formula_scope = meta_row$formula_scope,
      site_term_interpretation = meta_row$site_term_interpretation,
      uses_site = meta_row$uses_site,
      uses_age_sex_interaction = meta_row$uses_age_sex_interaction,
      site_placement_label = meta_row$site_placement_label,
      incidence_uses_site = meta_row$incidence_uses_site,
      latency_uses_site = meta_row$latency_uses_site,
      stage = meta_row$stage,
      branch = meta_row$branch,
      incidence_link = meta_row$incidence_link,
      model_id = meta_row$model_id,
      model_class = meta_row$model_class,
      model_block = meta_row$model_block,
      model_family = meta_row$model_family,
      latency_type = meta_row$latency_type,
      is_cure_model = meta_row$is_cure_model,
      matched_nocure_model_id = meta_row$matched_nocure_model_id
    )
  }))
}

stage7_fit_benchmark_km <- function(df, horizons) {
  fit_capture <- stage7_capture_with_warnings(survival::survfit(survival::Surv(time_year, event_main) ~ 1, data = df))
  if (is.null(fit_capture$value)) {
    return(list(success = FALSE, error_message = fit_capture$error_message, warnings = fit_capture$warnings))
  }
  fit <- fit_capture$value
  surv_at <- stage7_step_eval(fit$time, fit$surv, horizons, yleft = 1, yright = tail(fit$surv, 1))
  overall_surv_mat <- matrix(rep(surv_at, each = nrow(df)), nrow = nrow(df), ncol = length(horizons))
  susceptible_surv_mat <- overall_surv_mat
  list(
    success = TRUE,
    fit_object = fit,
    warnings = fit_capture$warnings,
    coefficients = tibble(),
    predictions = list(
      overall_survival = overall_surv_mat,
      susceptible_survival = susceptible_surv_mat,
      uncured_prob = rep(1, nrow(df))
    ),
    loglik = NA_real_,
    n_parameters = 0,
    convergence_code = 0,
    has_converged_solution = TRUE,
    optimizer_method = "Kaplan-Meier"
  )
}

stage7_fit_nocure_survreg <- function(df, latency_rhs, family) {
  dist_map <- c(exp = "exponential", weibull = "weibull", lnorm = "lognormal", llogis = "loglogistic")
  surv_formula <- stage7_make_surv_formula(latency_rhs)
  fit_capture <- stage7_capture_with_warnings(
    survival::survreg(surv_formula, data = df, dist = unname(dist_map[[family]]))
  )
  if (is.null(fit_capture$value)) {
    return(list(success = FALSE, error_message = fit_capture$error_message, warnings = fit_capture$warnings))
  }
  fit <- fit_capture$value
  mm <- stats::model.matrix(delete.response(stats::terms(fit)), data = df)
  lp <- drop(mm %*% stats::coef(fit))
  aux_log <- if (family == "exp") NA_real_ else log(fit$scale)
  list(
    success = TRUE,
    fit_object = fit,
    warnings = fit_capture$warnings,
    coefficients = bind_rows(
      tibble(component = "latency", term = stage7_standardize_term_name(names(stats::coef(fit))), estimate = as.numeric(stats::coef(fit))),
      if (family == "exp") tibble() else tibble(component = "latency_auxiliary", term = "log_shape_or_sigma", estimate = aux_log)
    ),
    predictions = list(
      lp = lp,
      aux_log = aux_log,
      family = family
    ),
    loglik = if (!is.null(fit$loglik) && length(fit$loglik) >= 2L) as.numeric(fit$loglik[[2]]) else NA_real_,
    n_parameters = length(stats::coef(fit)) + ifelse(family == "exp", 0, 1),
    convergence_code = 0,
    has_converged_solution = TRUE,
    optimizer_method = "survival::survreg"
  )
}

stage7_predict_nocure_survreg <- function(fit_result, horizons, n_subject) {
  lp <- fit_result$predictions$lp
  family <- fit_result$predictions$family
  aux_log <- fit_result$predictions$aux_log
  surv_mat <- sapply(horizons, function(h) stage7_parametric_surv_density(rep(h, n_subject), lp, family, aux_log)$surv)
  if (is.null(dim(surv_mat))) surv_mat <- matrix(surv_mat, nrow = n_subject, ncol = length(horizons))
  list(
    overall_survival = surv_mat,
    susceptible_survival = surv_mat,
    uncured_prob = rep(1, n_subject)
  )
}

stage7_fit_nocure_coxph <- function(df, latency_rhs) {
  surv_formula <- stage7_make_surv_formula(latency_rhs)
  fit_capture <- stage7_capture_with_warnings(
    survival::coxph(surv_formula, data = df, x = TRUE, model = TRUE)
  )
  if (is.null(fit_capture$value)) {
    return(list(success = FALSE, error_message = fit_capture$error_message, warnings = fit_capture$warnings))
  }
  fit <- fit_capture$value
  basehaz_df <- survival::basehaz(fit, centered = FALSE)
  lp <- stats::predict(fit, newdata = df, type = "lp")
  list(
    success = TRUE,
    fit_object = fit,
    warnings = fit_capture$warnings,
    coefficients = tibble(component = "latency", term = stage7_standardize_term_name(names(stats::coef(fit))), estimate = as.numeric(stats::coef(fit))),
    predictions = list(basehaz = basehaz_df, lp = lp),
    loglik = if (!is.null(fit$loglik) && length(fit$loglik) >= 2L) as.numeric(fit$loglik[[2]]) else NA_real_,
    n_parameters = length(stats::coef(fit)),
    convergence_code = 0,
    has_converged_solution = TRUE,
    optimizer_method = "survival::coxph"
  )
}

stage7_predict_nocure_coxph <- function(fit_result, horizons, n_subject) {
  H0 <- stage7_step_eval(fit_result$predictions$basehaz$time, fit_result$predictions$basehaz$hazard, horizons, yleft = 0, yright = tail(fit_result$predictions$basehaz$hazard, 1))
  surv_mat <- sapply(H0, function(h0) exp(-h0 * exp(fit_result$predictions$lp)))
  if (is.null(dim(surv_mat))) surv_mat <- matrix(surv_mat, nrow = n_subject, ncol = length(horizons))
  list(
    overall_survival = surv_mat,
    susceptible_survival = surv_mat,
    uncured_prob = rep(1, n_subject)
  )
}

stage7_fit_parametric_cure <- function(df, incidence_rhs, latency_rhs, family) {
  X_inc <- stage7_model_matrix(df, incidence_rhs, drop_intercept = FALSE)
  X_lat <- stage7_model_matrix(df, latency_rhs, drop_intercept = FALSE)
  y_time <- pmax(as.numeric(df$time_year), 1e-10)
  y_event <- as.integer(df$event_main)
  
  start_gamma <- tryCatch({
    stats::coef(stats::glm(stats::as.formula(paste("event_main ~", stage7_trim_rhs(incidence_rhs))), data = df, family = stats::binomial(link = "logit")))
  }, error = function(e) rep(0, ncol(X_inc)))
  start_gamma <- start_gamma[colnames(X_inc)]
  start_gamma[is.na(start_gamma)] <- 0
  
  start_survreg <- tryCatch({
    dist_map <- c(exp = "exponential", weibull = "weibull", lnorm = "lognormal", llogis = "loglogistic")
    survival::survreg(stage7_make_surv_formula(latency_rhs), data = df, dist = unname(dist_map[[family]]))
  }, error = function(e) NULL)
  
  if (!is.null(start_survreg)) {
    start_beta <- stats::coef(start_survreg)[colnames(X_lat)]
    start_beta[is.na(start_beta)] <- 0
    start_aux <- if (family == "exp") numeric(0) else log(start_survreg$scale)
  } else {
    start_beta <- rep(0, ncol(X_lat))
    names(start_beta) <- colnames(X_lat)
    start_aux <- if (family == "exp") numeric(0) else 0
  }
  
  par_start <- c(start_gamma, start_beta, start_aux)
  
  neg_loglik <- function(par) {
    gamma <- par[seq_len(ncol(X_inc))]
    beta <- par[ncol(X_inc) + seq_len(ncol(X_lat))]
    aux_log <- if (family == "exp") NA_real_ else par[[length(par)]]
    uncured_prob <- stage7_clamp_prob(plogis(drop(X_inc %*% gamma)))
    lp_lat <- drop(X_lat %*% beta)
    dens_surv <- stage7_parametric_surv_density(y_time, lp_lat, family, aux_log)
    surv_uncured <- pmax(dens_surv$surv, 1e-300)
    dens_uncured <- pmax(dens_surv$dens, 1e-300)
    ll <- sum(y_event * (log(uncured_prob) + log(dens_uncured)) + (1 - y_event) * log((1 - uncured_prob) + uncured_prob * surv_uncured))
    if (!is.finite(ll)) return(1e12)
    -ll
  }
  
  fit_attempts <- list()
  fit_attempts[[1]] <- stage7_capture_with_warnings(stats::optim(par_start, neg_loglik, method = "BFGS", control = list(maxit = 2000, reltol = 1e-10)))
  fit_attempts[[2]] <- stage7_capture_with_warnings(stats::optim(par_start, neg_loglik, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-10)))
  
  valid_attempts <- keep(fit_attempts, ~ !is.null(.x$value) && is.list(.x$value) && is.finite(.x$value$value))
  if (length(valid_attempts) == 0L) {
    all_warnings <- unlist(map(fit_attempts, "warnings"), use.names = FALSE)
    all_errors <- unlist(map(fit_attempts, "error_message"), use.names = FALSE)
    return(list(success = FALSE, error_message = paste(na.omit(all_errors), collapse = " | "), warnings = all_warnings))
  }
  
  best_idx <- which.min(vapply(valid_attempts, function(x) x$value$value, numeric(1)))
  fit_capture <- valid_attempts[[best_idx]]
  opt <- fit_capture$value
  par_hat <- opt$par
  gamma_hat <- par_hat[seq_len(ncol(X_inc))]
  beta_hat <- par_hat[ncol(X_inc) + seq_len(ncol(X_lat))]
  aux_log_hat <- if (family == "exp") NA_real_ else par_hat[[length(par_hat)]]
  uncured_prob <- stage7_clamp_prob(plogis(drop(X_inc %*% gamma_hat)))
  lp_lat <- drop(X_lat %*% beta_hat)
  
  coef_tbl <- bind_rows(
    tibble(component = "incidence_uncured_logit", term = stage7_standardize_term_name(colnames(X_inc)), estimate = as.numeric(gamma_hat)),
    tibble(component = "latency", term = stage7_standardize_term_name(colnames(X_lat)), estimate = as.numeric(beta_hat)),
    if (family == "exp") tibble() else tibble(component = "latency_auxiliary", term = "log_shape_or_sigma", estimate = aux_log_hat)
  )
  
  list(
    success = TRUE,
    fit_object = list(par = par_hat, family = family, X_inc = X_inc, X_lat = X_lat, opt = opt),
    warnings = fit_capture$warnings,
    coefficients = coef_tbl,
    predictions = list(gamma = gamma_hat, beta = beta_hat, aux_log = aux_log_hat, family = family, uncured_prob = uncured_prob, lp_lat = lp_lat),
    loglik = -opt$value,
    n_parameters = length(par_hat),
    convergence_code = as.numeric(opt$convergence),
    has_converged_solution = identical(as.numeric(opt$convergence), 0),
    optimizer_method = paste0("optim::", opt$method %||% "optim")
  )
}

stage7_predict_parametric_cure <- function(fit_result, horizons, n_subject) {
  pred <- fit_result$predictions
  surv_mat <- sapply(horizons, function(h) stage7_parametric_surv_density(rep(h, n_subject), pred$lp_lat, pred$family, pred$aux_log)$surv)
  if (is.null(dim(surv_mat))) surv_mat <- matrix(surv_mat, nrow = n_subject, ncol = length(horizons))
  overall_surv_mat <- (1 - pred$uncured_prob) + pred$uncured_prob * surv_mat
  list(
    overall_survival = overall_surv_mat,
    susceptible_survival = surv_mat,
    uncured_prob = pred$uncured_prob
  )
}

stage7_fit_smcure_model <- function(df, incidence_rhs, latency_rhs, model_type) {
  if (!requireNamespace("smcure", quietly = TRUE)) {
    return(list(success = FALSE, error_message = "Package `smcure` is required for semiparametric cure models but is not installed.", warnings = character()))
  }
  
  surv_formula <- stage7_make_surv_formula(latency_rhs)
  cure_formula <- stage7_make_rhs_formula(incidence_rhs)
  fit_capture <- stage7_capture_with_warnings(
    smcure::smcure(
      formula = surv_formula,
      cureform = cure_formula,
      data = df,
      model = ifelse(model_type == "ph", "ph", "aft"),
      link = "logit",
      Var = FALSE,
      emmax = 200,
      eps = 1e-07
    )
  )
  if (is.null(fit_capture$value)) {
    return(list(success = FALSE, error_message = fit_capture$error_message, warnings = fit_capture$warnings))
  }
  fit <- fit_capture$value
  b_names <- fit$bnm %||% names(fit$b) %||% rep("incidence", length(fit$b))
  beta_names <- fit$betanm %||% names(fit$beta) %||% rep("latency", length(fit$beta))
  coef_tbl <- bind_rows(
    tibble(component = "incidence_uncured_logit", term = stage7_standardize_term_name(b_names), estimate = as.numeric(fit$b)),
    tibble(component = "latency", term = stage7_standardize_term_name(beta_names), estimate = as.numeric(fit$beta))
  )
  list(
    success = TRUE,
    fit_object = fit,
    warnings = fit_capture$warnings,
    coefficients = coef_tbl,
    predictions = list(model_type = model_type),
    loglik = NA_real_,
    n_parameters = length(fit$b) + length(fit$beta),
    convergence_code = 0,
    has_converged_solution = TRUE,
    optimizer_method = paste0("smcure::smcure_", model_type)
  )
}

stage7_predict_smcure <- function(fit_result, df, incidence_rhs, latency_rhs, horizons) {
  fit <- fit_result$fit_object
  model_type <- fit_result$predictions$model_type
  
  inc_coef_names <- fit$bnm %||% names(fit$b)
  mm_inc_raw <- stage7_model_matrix(df, incidence_rhs, drop_intercept = FALSE)
  Z <- stage7_align_mm_to_coef(mm_inc_raw, inc_coef_names)
  b_vec <- as.numeric(fit$b)
  uncured_prob <- stage7_clamp_prob(plogis(drop(Z %*% b_vec)))
  
  if (model_type == "ph") {
    lat_coef_names <- fit$betanm %||% names(fit$beta)
    mm_lat_raw <- stage7_model_matrix(df, latency_rhs, drop_intercept = TRUE)
    X <- stage7_align_mm_to_coef(mm_lat_raw, lat_coef_names)
    beta_vec <- as.numeric(fit$beta)
    
    base_tbl <- tibble(time = as.numeric(fit$Time), surv = as.numeric(fit$s)) %>%
      filter(is.finite(time), is.finite(surv)) %>%
      group_by(time) %>%
      summarise(surv = min(surv), .groups = "drop") %>%
      arrange(time)
    
    lp <- drop(X %*% beta_vec)
    base_surv_at <- stage7_step_eval(base_tbl$time, base_tbl$surv, horizons, yleft = 1, yright = tail(base_tbl$surv, 1))
    susceptible_surv_mat <- sapply(base_surv_at, function(s0) pmin(pmax(s0, 0), 1) ^ exp(lp))
    if (is.null(dim(susceptible_surv_mat))) {
      susceptible_surv_mat <- matrix(susceptible_surv_mat, nrow = nrow(df), ncol = length(horizons))
    }
  } else {
    lat_coef_names <- fit$betanm %||% names(fit$beta)
    mm_lat_raw <- stage7_model_matrix(df, latency_rhs, drop_intercept = FALSE)
    X <- stage7_align_mm_to_coef(mm_lat_raw, lat_coef_names)
    beta_vec <- as.numeric(fit$beta)
    
    s0 <- as.numeric(fit$s)
    err <- as.numeric(fit$error)
    ebetaX <- drop(exp(X %*% beta_vec))
    
    susceptible_surv_mat <- sapply(seq_along(ebetaX), function(j) {
      subj_time <- ebetaX[[j]] * exp(err)
      ord <- order(subj_time)
      stage7_step_eval(subj_time[ord], s0[ord], horizons, yleft = 1, yright = tail(s0[ord], 1))
    })
    susceptible_surv_mat <- t(susceptible_surv_mat)
  }
  
  overall_surv_mat <- (1 - uncured_prob) + uncured_prob * susceptible_surv_mat
  list(
    overall_survival = overall_surv_mat,
    susceptible_survival = susceptible_surv_mat,
    uncured_prob = uncured_prob
  )
}

stage7_summarise_risk_summary <- function(subject_predictions_tbl) {
  if (nrow(subject_predictions_tbl) == 0L) return(tibble())
  subject_predictions_tbl %>%
    group_by(dataset, formula_variant, site_branch, interaction_branch, model_id, model_class, model_block, model_family, latency_type, horizon_year, matched_nocure_model_id) %>%
    summarise(
      n_subject = dplyr::n(),
      mean_survival_overall = mean(predicted_survival_overall, na.rm = TRUE),
      mean_risk_overall = mean(predicted_risk_overall, na.rm = TRUE),
      median_risk_overall = safe_quantile(predicted_risk_overall, 0.50),
      q25_risk_overall = safe_quantile(predicted_risk_overall, 0.25),
      q75_risk_overall = safe_quantile(predicted_risk_overall, 0.75),
      mean_survival_susceptible = mean(predicted_survival_susceptible, na.rm = TRUE),
      mean_risk_susceptible = mean(predicted_risk_susceptible, na.rm = TRUE),
      mean_cure_fraction = mean(predicted_cure_fraction, na.rm = TRUE),
      median_cure_fraction = safe_quantile(predicted_cure_fraction, 0.50),
      q25_cure_fraction = safe_quantile(predicted_cure_fraction, 0.25),
      q75_cure_fraction = safe_quantile(predicted_cure_fraction, 0.75),
      boot_n_success_mean_risk = NA_real_,
      boot_success_rate_mean_risk = NA_real_,
      mean_risk_ci_lower = NA_real_,
      mean_risk_ci_upper = NA_real_,
      boot_n_success_mean_cure = NA_real_,
      boot_success_rate_mean_cure = NA_real_,
      mean_cure_fraction_ci_lower = NA_real_,
      mean_cure_fraction_ci_upper = NA_real_,
      .groups = "drop"
    )
}

stage7_summarise_delta_risk <- function(subject_predictions_tbl) {
  if (nrow(subject_predictions_tbl) == 0L) return(tibble())
  
  cure_preds <- subject_predictions_tbl %>%
    filter(as.logical(is_cure_model), !is.na(matched_nocure_model_id), matched_nocure_model_id != "") %>%
    select(
      unique_person_id, dataset, formula_variant, site_branch, interaction_branch, horizon_year,
      model_id, matched_nocure_model_id, predicted_risk_overall
    ) %>%
    rename(cure_model_id = model_id, cure_predicted_risk = predicted_risk_overall)
  
  noncure_preds <- subject_predictions_tbl %>%
    filter(!as.logical(is_cure_model)) %>%
    select(unique_person_id, dataset, formula_variant, site_branch, interaction_branch, horizon_year, model_id, predicted_risk_overall) %>%
    rename(matched_nocure_model_id = model_id, nocure_predicted_risk = predicted_risk_overall)
  
  joined <- cure_preds %>%
    inner_join(
      noncure_preds,
      by = c("unique_person_id", "dataset", "formula_variant", "site_branch", "interaction_branch", "horizon_year", "matched_nocure_model_id")
    ) %>%
    mutate(delta_risk_nc_minus_cure = nocure_predicted_risk - cure_predicted_risk)
  
  if (nrow(joined) == 0L) {
    warning("stage7_summarise_delta_risk(): cure/non-cure join returned 0 rows.", call. = FALSE)
    return(tibble())
  }
  
  subject_predictions_tbl %>%
    filter(model_id %in% unique(joined$cure_model_id)) %>%
    distinct(dataset, formula_variant, site_branch, interaction_branch, model_id, matched_nocure_model_id) %>%
    inner_join(
      joined %>%
        group_by(dataset, formula_variant, site_branch, interaction_branch, cure_model_id, matched_nocure_model_id, horizon_year) %>%
        summarise(
          n_subject = dplyr::n(),
          mean_delta_risk_nc_minus_cure = mean(delta_risk_nc_minus_cure, na.rm = TRUE),
          median_delta_risk_nc_minus_cure = safe_quantile(delta_risk_nc_minus_cure, 0.50),
          q25_delta_risk_nc_minus_cure = safe_quantile(delta_risk_nc_minus_cure, 0.25),
          q75_delta_risk_nc_minus_cure = safe_quantile(delta_risk_nc_minus_cure, 0.75),
          boot_n_success_delta = NA_real_,
          boot_success_rate_delta = NA_real_,
          mean_delta_risk_ci_lower = NA_real_,
          mean_delta_risk_ci_upper = NA_real_,
          .groups = "drop"
        ) %>%
        rename(model_id = cure_model_id),
      by = c("dataset", "formula_variant", "site_branch", "interaction_branch", "model_id", "matched_nocure_model_id")
    )
}

stage7_compute_threshold_metrics_raw <- function(subject_predictions_tbl, ipcw_tbl, thresholds) {
  if (nrow(subject_predictions_tbl) == 0L || nrow(ipcw_tbl) == 0L) return(tibble())
  pred_tbl <- subject_predictions_tbl %>%
    select(unique_person_id, dataset, formula_variant, site_branch, interaction_branch, model_id, model_class, model_block, model_family, latency_type, matched_nocure_model_id, horizon_year, predicted_risk_overall) %>%
    inner_join(ipcw_tbl, by = c("dataset", "unique_person_id", "horizon_year"))
  
  if (nrow(pred_tbl) == 0L) return(tibble())
  
  bind_rows(lapply(as.numeric(thresholds), function(thr) {
    pred_tbl %>%
      mutate(pred_positive = predicted_risk_overall >= thr) %>%
      group_by(dataset, formula_variant, site_branch, interaction_branch, model_id, model_class, model_block, model_family, latency_type, matched_nocure_model_id, horizon_year) %>%
      summarise(
        threshold = thr,
        n_subject = dplyr::n(),
        positive_classification_rate = mean(pred_positive, na.rm = TRUE),
        weighted_tp = sum(ifelse(pred_positive, weight_case, 0), na.rm = TRUE),
        weighted_fn = sum(ifelse(!pred_positive, weight_case, 0), na.rm = TRUE),
        weighted_fp = sum(ifelse(pred_positive, weight_control, 0), na.rm = TRUE),
        weighted_tn = sum(ifelse(!pred_positive, weight_control, 0), na.rm = TRUE),
        sensitivity = safe_divide(weighted_tp, weighted_tp + weighted_fn),
        specificity = safe_divide(weighted_tn, weighted_tn + weighted_fp),
        ppv = safe_divide(weighted_tp, weighted_tp + weighted_fp),
        false_positive_weighted_per_n = safe_divide(weighted_fp, n_subject),
        false_positive_burden_non_event = safe_divide(weighted_fp, weighted_fp + weighted_tn),
        false_positive_burden_primary = safe_divide(weighted_fp, weighted_fp + weighted_tn),
        false_positive_burden_all = safe_divide(weighted_fp, n_subject),
        unnecessary_high_risk_per_100_population = 100 * safe_divide(weighted_fp, n_subject),
        unnecessary_high_risk_per_100_non_event = 100 * safe_divide(weighted_fp, weighted_fp + weighted_tn),
        net_benefit = safe_divide(weighted_tp, n_subject) - safe_divide(weighted_fp, n_subject) * (thr / pmax(1 - thr, 1e-8)),
        boot_n_success_fp_primary = NA_real_,
        boot_success_rate_fp_primary = NA_real_,
        false_positive_burden_primary_ci_lower = NA_real_,
        false_positive_burden_primary_ci_upper = NA_real_,
        boot_n_success_fp_weighted = NA_real_,
        boot_success_rate_fp_weighted = NA_real_,
        false_positive_weighted_per_n_ci_lower = NA_real_,
        false_positive_weighted_per_n_ci_upper = NA_real_,
        boot_n_success_nb = NA_real_,
        boot_success_rate_nb = NA_real_,
        net_benefit_ci_lower = NA_real_,
        net_benefit_ci_upper = NA_real_,
        false_positive_burden_ci_lower = NA_real_,
        false_positive_burden_ci_upper = NA_real_,
        .groups = "drop"
      )
  }))
}

stage7_build_bootstrap_qc_placeholder <- function(fit_registry_tbl) {
  if (nrow(fit_registry_tbl) == 0L) return(tibble())
  fit_registry_tbl %>%
    transmute(
      dataset, formula_variant, site_branch, interaction_branch, model_id, model_class, model_block, model_family, latency_type, matched_nocure_model_id,
      expected_reps = as.integer(stage7_expected_bootstrap_reps),
      n_rep_observed = 0L,
      n_rep_fit_success = 0L,
      n_rep_prediction_success = 0L,
      n_rep_overall_success = 0L,
      n_rep_fit_failed = 0L,
      n_rep_prediction_failed = 0L,
      n_rep_warning = 0L,
      n_rep_nonzero_convergence = 0L,
      fit_success_rate = NA_real_,
      prediction_success_rate = NA_real_,
      overall_success_rate = NA_real_,
      ci_min_success_required = NA_real_,
      ci_eligible = FALSE,
      qc_flag = ifelse(stage7_expected_bootstrap_reps > 0L, "bootstrap_not_run", "bootstrap_not_requested")
    )
}

stage7_fit_core_outputs <- function(run_root_dir, stage1_inputs) {
  run_root_dir <- normalize_existing_path(run_root_dir)
  dir.create(run_root_dir, recursive = TRUE, showWarnings = FALSE)
  
  formula_plan_obj <- stage7_prepare_formula_plan(stage1_inputs)
  model_plan <- formula_plan_obj$model_plan
  benchmark_plan <- formula_plan_obj$benchmark_plan
  model_catalog <- stage7_model_catalog()
  horizons <- sort(unique(as.integer(stage1_inputs$horizon_registry$horizon_year)))
  thresholds <- sort(unique(as.numeric(stage1_inputs$threshold_registry$threshold)))
  
  fit_rows <- list()
  coef_rows <- list()
  pred_rows <- list()
  fit_counter <- 0L
  
  stage7_prepare_dataset <- function(df, dataset_name) {
    df <- tibble::as_tibble(df)
    if (!"dataset" %in% names(df)) df$dataset <- dataset_name
    if (!"unique_person_id" %in% names(df)) {
      if (all(c("site", "id") %in% names(df))) {
        df$unique_person_id <- paste0(df$site, "__", df$id)
      } else {
        stop(sprintf("Stage 1 dataset `%s` does not contain `unique_person_id` or the fallback `site` + `id` columns.", dataset_name), call. = FALSE)
      }
    }
    if (!"site" %in% names(df)) df$site <- dataset_name
    df$site <- factor(as.character(df$site))
    df$sex_num <- suppressWarnings(as.numeric(df$sex_num))
    df$time_year <- suppressWarnings(as.numeric(df$time_year))
    df$event_main <- as.integer(df$event_main)
    if (!"censor_main" %in% names(df)) df$censor_main <- as.integer(df$event_main == 0)
    df
  }
  
  run_one_model <- function(df, meta_row, model_row) {
    t0 <- proc.time()[[3]]
    fit_result <- switch(
      model_row$fit_method,
      benchmark_km = stage7_fit_benchmark_km(df, horizons),
      nocure_exp = stage7_fit_nocure_survreg(df, meta_row$latency_rhs, "exp"),
      nocure_weibull = stage7_fit_nocure_survreg(df, meta_row$latency_rhs, "weibull"),
      nocure_lnorm = stage7_fit_nocure_survreg(df, meta_row$latency_rhs, "lnorm"),
      nocure_llogis = stage7_fit_nocure_survreg(df, meta_row$latency_rhs, "llogis"),
      nocure_coxph = stage7_fit_nocure_coxph(df, meta_row$latency_rhs),
      cure_exp = stage7_fit_parametric_cure(df, meta_row$incidence_rhs, meta_row$latency_rhs, "exp"),
      cure_weibull = stage7_fit_parametric_cure(df, meta_row$incidence_rhs, meta_row$latency_rhs, "weibull"),
      cure_lnorm = stage7_fit_parametric_cure(df, meta_row$incidence_rhs, meta_row$latency_rhs, "lnorm"),
      cure_llogis = stage7_fit_parametric_cure(df, meta_row$incidence_rhs, meta_row$latency_rhs, "llogis"),
      cure_coxlatency = stage7_fit_smcure_model(df, meta_row$incidence_rhs, meta_row$latency_rhs, model_type = "ph"),
      cure_aft_sensitivity = stage7_fit_smcure_model(df, meta_row$incidence_rhs, meta_row$latency_rhs, model_type = "aft"),
      stop(sprintf("Unsupported fit method: %s", model_row$fit_method), call. = FALSE)
    )
    fit_elapsed <- proc.time()[[3]] - t0
    
    prediction_result <- list(value = NULL, warnings = character(), error_message = NA_character_)
    if (isTRUE(fit_result$success)) {
      t1 <- proc.time()[[3]]
      prediction_result <- switch(
        model_row$fit_method,
        benchmark_km = list(value = fit_result$predictions, warnings = character(), error_message = NA_character_),
        nocure_exp = list(value = stage7_predict_nocure_survreg(fit_result, horizons, nrow(df)), warnings = character(), error_message = NA_character_),
        nocure_weibull = list(value = stage7_predict_nocure_survreg(fit_result, horizons, nrow(df)), warnings = character(), error_message = NA_character_),
        nocure_lnorm = list(value = stage7_predict_nocure_survreg(fit_result, horizons, nrow(df)), warnings = character(), error_message = NA_character_),
        nocure_llogis = list(value = stage7_predict_nocure_survreg(fit_result, horizons, nrow(df)), warnings = character(), error_message = NA_character_),
        nocure_coxph = list(value = stage7_predict_nocure_coxph(fit_result, horizons, nrow(df)), warnings = character(), error_message = NA_character_),
        cure_exp = list(value = stage7_predict_parametric_cure(fit_result, horizons, nrow(df)), warnings = character(), error_message = NA_character_),
        cure_weibull = list(value = stage7_predict_parametric_cure(fit_result, horizons, nrow(df)), warnings = character(), error_message = NA_character_),
        cure_lnorm = list(value = stage7_predict_parametric_cure(fit_result, horizons, nrow(df)), warnings = character(), error_message = NA_character_),
        cure_llogis = list(value = stage7_predict_parametric_cure(fit_result, horizons, nrow(df)), warnings = character(), error_message = NA_character_),
        cure_coxlatency = stage7_capture_with_warnings(stage7_predict_smcure(fit_result, df, meta_row$incidence_rhs, meta_row$latency_rhs, horizons)),
        cure_aft_sensitivity = stage7_capture_with_warnings(stage7_predict_smcure(fit_result, df, meta_row$incidence_rhs, meta_row$latency_rhs, horizons)),
        stop(sprintf("Unsupported prediction method: %s", model_row$fit_method), call. = FALSE)
      )
      prediction_elapsed <- proc.time()[[3]] - t1
    } else {
      prediction_elapsed <- NA_real_
    }
    
    model_meta <- meta_row %>%
      mutate(
        model_id = model_row$model_id,
        model_class = model_row$model_class,
        model_block = model_row$model_block,
        model_family = model_row$model_family,
        latency_type = model_row$latency_type,
        fit_engine = model_row$fit_engine,
        is_cure_model = model_row$is_cure_model,
        matched_nocure_model_id = model_row$matched_nocure_model_id
      )
    
    save_file <- NA_character_
    if (isTRUE(save_stage7_fit_rds) && isTRUE(fit_result$success) && !is.null(fit_result$fit_object)) {
      save_file <- file.path(
        run_root_dir,
        paste0(stage7_model_rds_prefix, "__", model_meta$dataset, "__", model_meta$formula_variant, "__", model_meta$model_id, ".rds")
      )
      saveRDS(fit_result$fit_object, save_file)
    }
    
    fit_ok <- isTRUE(fit_result$success)
    pred_ok <- fit_ok && !is.null(prediction_result$value)
    fit_warning_text <- collapse_unique_sorted(fit_result$warnings, empty_value = "")
    pred_warning_text <- collapse_unique_sorted(prediction_result$warnings, empty_value = "")
    
    fit_row <- tibble(
      dataset = model_meta$dataset,
      formula_variant = model_meta$formula_variant,
      site_branch = model_meta$site_branch,
      interaction_branch = model_meta$interaction_branch,
      model_id = model_meta$model_id,
      model_class = model_meta$model_class,
      model_block = model_meta$model_block,
      model_family = model_meta$model_family,
      latency_type = model_meta$latency_type,
      fit_engine = model_meta$fit_engine,
      is_cure_model = model_meta$is_cure_model,
      matched_nocure_model_id = model_meta$matched_nocure_model_id,
      fit_status = ifelse(fit_ok, "ok_fit", "fit_error"),
      prediction_status = ifelse(pred_ok, "ok_prediction", ifelse(fit_ok, "prediction_error", "prediction_not_attempted")),
      overall_status = ifelse(fit_ok && pred_ok, "ok", ifelse(fit_ok, "prediction_error", "fit_error")),
      fit_component_success = fit_ok,
      prediction_component_success = pred_ok,
      convergence_code = ifelse(fit_ok, fit_result$convergence_code %||% NA_real_, NA_real_),
      has_converged_solution = ifelse(fit_ok, fit_result$has_converged_solution %||% NA, NA),
      optimizer_method = ifelse(fit_ok, fit_result$optimizer_method %||% NA_character_, NA_character_),
      fit_source = ifelse(fit_ok, "fresh_core_fit", "fresh_core_fit_failed"),
      fit_cache_validation = ifelse(fit_ok, "freshly_fitted", "failed_during_fresh_fit"),
      loglik = ifelse(fit_ok, fit_result$loglik %||% NA_real_, NA_real_),
      n_parameters = ifelse(fit_ok, fit_result$n_parameters %||% NA_real_, NA_real_),
      AIC = ifelse(fit_ok && is.finite(fit_result$loglik %||% NA_real_) && is.finite(fit_result$n_parameters %||% NA_real_), -2 * fit_result$loglik + 2 * fit_result$n_parameters, NA_real_),
      BIC = ifelse(fit_ok && is.finite(fit_result$loglik %||% NA_real_) && is.finite(fit_result$n_parameters %||% NA_real_), -2 * fit_result$loglik + log(nrow(df)) * fit_result$n_parameters, NA_real_),
      fit_message = dplyr::coalesce(fit_result$error_message, ""),
      prediction_message = dplyr::coalesce(prediction_result$error_message, ""),
      message = paste(c(na.omit(c(fit_result$error_message, prediction_result$error_message))), collapse = " | "),
      error_message = paste(c(na.omit(c(fit_result$error_message, prediction_result$error_message))), collapse = " | "),
      fit_warning_count = length(fit_result$warnings),
      fit_warning_message = fit_warning_text,
      prediction_warning_count = length(prediction_result$warnings),
      prediction_warning_message = pred_warning_text,
      warning_count = length(fit_result$warnings) + length(prediction_result$warnings),
      warning_message = collapse_unique_sorted(c(fit_result$warnings, prediction_result$warnings), empty_value = ""),
      warning_count_total = length(fit_result$warnings) + length(prediction_result$warnings),
      fit_elapsed_sec = fit_elapsed,
      prediction_elapsed_sec = prediction_elapsed,
      elapsed_sec = fit_elapsed + dplyr::coalesce(prediction_elapsed, 0),
      prediction_rows = ifelse(pred_ok, nrow(df) * length(horizons), 0),
      expected_prediction_rows = nrow(df) * length(horizons),
      saved_rds_file = save_file,
      n = nrow(df),
      n_transition = sum(df$event_main == 1, na.rm = TRUE),
      n_right_censor = if ("status_num" %in% names(df)) sum(df$status_num == 0, na.rm = TRUE) else sum(df$event_main == 0 & df$censor_main == 1, na.rm = TRUE),
      n_remission = if ("status_num" %in% names(df)) sum(df$status_num == 2, na.rm = TRUE) else NA_real_,
      incidence_rhs = model_meta$incidence_rhs,
      latency_rhs = model_meta$latency_rhs,
      formula_label = model_meta$formula_label,
      formula_scope = model_meta$formula_scope,
      site_term_interpretation = model_meta$site_term_interpretation,
      uses_site = model_meta$uses_site,
      uses_age_sex_interaction = model_meta$uses_age_sex_interaction,
      site_placement_label = model_meta$site_placement_label,
      incidence_uses_site = model_meta$incidence_uses_site,
      latency_uses_site = model_meta$latency_uses_site,
      stage = model_meta$stage,
      branch = model_meta$branch,
      incidence_link = model_meta$incidence_link,
      dataset_key = make_stage6_join_key(model_meta$dataset, site_branch = model_meta$site_branch)
    )
    
    coef_tbl <- if (fit_ok && nrow(fit_result$coefficients) > 0L) {
      fit_result$coefficients %>%
        mutate(
          dataset = model_meta$dataset,
          formula_variant = model_meta$formula_variant,
          site_branch = model_meta$site_branch,
          interaction_branch = model_meta$interaction_branch,
          model_id = model_meta$model_id,
          model_class = model_meta$model_class,
          model_block = model_meta$model_block,
          model_family = model_meta$model_family,
          latency_type = model_meta$latency_type,
          matched_nocure_model_id = model_meta$matched_nocure_model_id,
          formula_label = model_meta$formula_label,
          formula_scope = model_meta$formula_scope,
          site_term_interpretation = model_meta$site_term_interpretation,
          uses_site = model_meta$uses_site,
          uses_age_sex_interaction = model_meta$uses_age_sex_interaction,
          site_placement_label = model_meta$site_placement_label,
          incidence_uses_site = model_meta$incidence_uses_site,
          latency_uses_site = model_meta$latency_uses_site,
          stage = model_meta$stage,
          branch = model_meta$branch,
          incidence_link = model_meta$incidence_link,
          dataset_key = make_stage6_join_key(model_meta$dataset, site_branch = model_meta$site_branch),
          boot_n_success_estimate = NA_real_,
          boot_success_rate_estimate = NA_real_,
          estimate_ci_lower = NA_real_,
          estimate_ci_upper = NA_real_
        ) %>%
        select(dataset, formula_variant, site_branch, interaction_branch, model_id, model_class, model_block, model_family, latency_type, component, term, estimate, boot_n_success_estimate, boot_success_rate_estimate, estimate_ci_lower, estimate_ci_upper, formula_label, formula_scope, site_term_interpretation, uses_site, uses_age_sex_interaction, site_placement_label, incidence_uses_site, latency_uses_site, stage, branch, incidence_link, matched_nocure_model_id, dataset_key)
    } else {
      tibble()
    }
    
    pred_tbl <- if (pred_ok) {
      stage7_build_prediction_rows(
        df = df,
        horizons = horizons,
        overall_surv_mat = prediction_result$value$overall_survival,
        susceptible_surv_mat = prediction_result$value$susceptible_survival,
        uncured_prob_vec = prediction_result$value$uncured_prob,
        meta_row = model_meta
      )
    } else {
      tibble()
    }
    
    list(fit_row = fit_row, coefficients = coef_tbl, predictions = pred_tbl)
  }
  
  for (dataset_name in names(stage1_inputs$analysis_datasets)) {
    df <- stage7_prepare_dataset(stage1_inputs$analysis_datasets[[dataset_name]], dataset_name)
    
    benchmark_meta <- benchmark_plan %>% filter(dataset == dataset_name)
    benchmark_spec <- model_catalog %>% filter(model_id == "benchmark_km")
    fit_counter <- fit_counter + 1L
    message(sprintf("[%03d] Fitting %s | %s | %s", fit_counter, dataset_name, benchmark_meta$formula_variant[[1]], benchmark_spec$model_id[[1]]))
    result <- run_one_model(df, benchmark_meta, benchmark_spec)
    fit_rows[[length(fit_rows) + 1L]] <- result$fit_row
    if (nrow(result$coefficients) > 0L) coef_rows[[length(coef_rows) + 1L]] <- result$coefficients
    if (nrow(result$predictions) > 0L) pred_rows[[length(pred_rows) + 1L]] <- result$predictions
    
    dataset_plan <- model_plan %>% filter(dataset == dataset_name)
    nonbenchmark_catalog <- model_catalog %>% filter(model_id != "benchmark_km")
    for (ii in seq_len(nrow(dataset_plan))) {
      meta_row <- dataset_plan[ii, , drop = FALSE]
      for (jj in seq_len(nrow(nonbenchmark_catalog))) {
        model_row <- nonbenchmark_catalog[jj, , drop = FALSE]
        fit_counter <- fit_counter + 1L
        message(sprintf("[%03d] Fitting %s | %s | %s", fit_counter, dataset_name, meta_row$formula_variant[[1]], model_row$model_id[[1]]))
        result <- run_one_model(df, meta_row, model_row)
        fit_rows[[length(fit_rows) + 1L]] <- result$fit_row
        if (nrow(result$coefficients) > 0L) coef_rows[[length(coef_rows) + 1L]] <- result$coefficients
        if (nrow(result$predictions) > 0L) pred_rows[[length(pred_rows) + 1L]] <- result$predictions
      }
    }
  }
  
  fit_registry_tbl <- bind_rows(fit_rows)
  coefficients_tbl <- bind_rows(coef_rows)
  subject_predictions_tbl <- bind_rows(pred_rows)
  
  risk_summary_tbl <- stage7_summarise_risk_summary(subject_predictions_tbl)
  delta_risk_tbl <- stage7_summarise_delta_risk(subject_predictions_tbl)
  ipcw_tbl <- bind_rows(lapply(names(stage1_inputs$analysis_datasets), function(nm) {
    df <- stage1_inputs$analysis_datasets[[nm]]
    df <- tibble::as_tibble(df)
    if (!"dataset" %in% names(df)) df$dataset <- nm
    compute_ipcw_frame(df, horizons_year = horizons)
  }))
  threshold_metrics_tbl <- stage7_compute_threshold_metrics_raw(subject_predictions_tbl, ipcw_tbl, thresholds)
  bootstrap_model_qc_tbl <- stage7_build_bootstrap_qc_placeholder(fit_registry_tbl)
  
  safe_write_csv(fit_registry_tbl, file.path(run_root_dir, "stage7_fit_registry.csv"))
  safe_write_csv(coefficients_tbl, file.path(run_root_dir, "stage7_coefficients.csv"))
  safe_write_csv(subject_predictions_tbl, file.path(run_root_dir, "stage7_subject_predictions.csv"))
  safe_write_csv(risk_summary_tbl, file.path(run_root_dir, "stage7_risk_summary.csv"))
  safe_write_csv(delta_risk_tbl, file.path(run_root_dir, "stage7_delta_risk.csv"))
  safe_write_csv(threshold_metrics_tbl, file.path(run_root_dir, "stage7_threshold_metrics.csv"))
  safe_write_csv(bootstrap_model_qc_tbl, file.path(run_root_dir, "stage7_bootstrap_model_qc.csv"))
  
  run_manifest_tbl <- tibble(
    file_name = c(
      "stage7_fit_registry.csv",
      "stage7_coefficients.csv",
      "stage7_subject_predictions.csv",
      "stage7_risk_summary.csv",
      "stage7_delta_risk.csv",
      "stage7_threshold_metrics.csv",
      "stage7_bootstrap_model_qc.csv"
    ),
    description = c(
      "Raw Stage 7 fit registry generated from the full fitting pipeline.",
      "Raw Stage 7 coefficient table generated from the full fitting pipeline.",
      "Raw Stage 7 subject-level prediction table generated from the full fitting pipeline.",
      "Raw Stage 7 risk summary generated from the full fitting pipeline.",
      "Raw Stage 7 DeltaRisk_NC_C summary generated from the full fitting pipeline.",
      "Raw Stage 7 threshold metrics generated from the full fitting pipeline.",
      "Bootstrap QC placeholder generated when no Stage 7 bootstrap bundle is available."
    ),
    file_path = file.path(run_root_dir, c(
      "stage7_fit_registry.csv",
      "stage7_coefficients.csv",
      "stage7_subject_predictions.csv",
      "stage7_risk_summary.csv",
      "stage7_delta_risk.csv",
      "stage7_threshold_metrics.csv",
      "stage7_bootstrap_model_qc.csv"
    )),
    generation_mode = "fresh_core_fit"
  )
  safe_write_csv(run_manifest_tbl, file.path(run_root_dir, "stage7_run_manifest.csv"))
  
  invisible(list(
    fit_registry = fit_registry_tbl,
    coefficients = coefficients_tbl,
    subject_predictions = subject_predictions_tbl,
    risk_summary = risk_summary_tbl,
    delta_risk = delta_risk_tbl,
    threshold_metrics = threshold_metrics_tbl,
    bootstrap_model_qc = bootstrap_model_qc_tbl
  ))
}

ensure_stage7_core_outputs <- function(run_root_dir, stage1_inputs) {
  resolved_run_root_dir <- resolve_stage7_run_root_or_create(run_root_dir)
  if (stage7_core_outputs_exist(resolved_run_root_dir) && isTRUE(skip_stage7_core_refit_if_outputs_exist) && !isTRUE(force_stage7_core_refit)) {
    message(sprintf("ℹ️  Reusing existing Stage 7 core outputs in: %s", resolved_run_root_dir))
    return(list(resolved_run_root_dir = resolved_run_root_dir, reused_existing = TRUE))
  }
  message(sprintf("ℹ️  Stage 7 core outputs were not fully available; fitting core models into: %s", resolved_run_root_dir))
  stage7_fit_core_outputs(resolved_run_root_dir, stage1_inputs = stage1_inputs)
  list(resolved_run_root_dir = resolved_run_root_dir, reused_existing = FALSE)
}

# 🔴 Load: current artifacts and supporting inputs ===============================
stage1_inputs <- load_stage1_backbone(stage1_data_path)
formula_lookup <- prepare_formula_lookup(stage1_inputs$formula_registry)
horizon_lookup <- prepare_horizon_lookup(stage1_inputs$horizon_registry)
stage6_carry_forward_tbl <- load_stage6_carry_forward_optional(stage6_screening_file)
stage6_registry <- build_stage6_registry(stage6_carry_forward_tbl)
stage7_core_status <- ensure_stage7_core_outputs(run_root_dir, stage1_inputs = stage1_inputs)
stage7_objects <- load_stage7_outputs(stage7_core_status$resolved_run_root_dir)
run_root_dir <- stage7_objects$resolved_root_dir
export_path <- resolve_stage7_export_path(
  export_path = export_path,
  resolved_run_root_dir = run_root_dir,
  suffix = stage7_refresh_export_suffix
)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)
stage5_fit_registry_tbl <- load_stage5_fit_metrics_optional(stage5_root_dir)

# 🔴 Transform: patch Stage 7 tables with Stage 1 and Stage 6 metadata ===============================
fit_registry_file_target <- make_output_path("stage7_fit_registry.csv")
fit_registry_tbl <- load_or_rebuild_csv_output(
  path = fit_registry_file_target,
  force_rebuild = refresh_patched_exports_even_when_outputs_exist,
  required_cols = c("dataset", "formula_variant", "risk_scale", "site_placement_label"),
  builder_fn = function() {
    stage7_objects$fit_registry %>%
      add_common_stage7_metadata(formula_lookup = formula_lookup, risk_scale_default = stage1_inputs$main_risk_scale) %>%
      join_stage6_flags(stage6_registry = stage6_registry)
  }
)

coefficients_file_target <- make_output_path("stage7_coefficients.csv")
coefficients_tbl <- load_or_rebuild_csv_output(
  path = coefficients_file_target,
  force_rebuild = refresh_patched_exports_even_when_outputs_exist,
  required_cols = c("dataset", "formula_variant", "risk_scale", "site_placement_label"),
  builder_fn = function() {
    stage7_objects$coefficients %>%
      add_common_stage7_metadata(formula_lookup = formula_lookup, risk_scale_default = stage1_inputs$main_risk_scale) %>%
      join_stage6_flags(stage6_registry = stage6_registry)
  }
)

subject_predictions_file_target <- make_output_path("stage7_subject_predictions.csv")
subject_predictions_tbl <- load_or_rebuild_csv_output(
  path = subject_predictions_file_target,
  force_rebuild = refresh_patched_exports_even_when_outputs_exist,
  required_cols = c("dataset", "formula_variant", "risk_scale", "support_tier", "horizon_evidence_class", "claim_restriction_flag"),
  builder_fn = function() {
    stage7_objects$subject_predictions %>%
      add_common_stage7_metadata(formula_lookup = formula_lookup, risk_scale_default = stage1_inputs$main_risk_scale) %>%
      join_stage6_flags(stage6_registry = stage6_registry) %>%
      join_horizon_support(horizon_lookup = horizon_lookup)
  }
)

risk_summary_file_target <- make_output_path("stage7_risk_summary.csv")
risk_summary_tbl <- load_or_rebuild_csv_output(
  path = risk_summary_file_target,
  force_rebuild = refresh_patched_exports_even_when_outputs_exist,
  required_cols = c("dataset", "formula_variant", "risk_scale", "support_tier", "horizon_evidence_class", "claim_restriction_flag"),
  builder_fn = function() {
    stage7_objects$risk_summary %>%
      add_common_stage7_metadata(formula_lookup = formula_lookup, risk_scale_default = stage1_inputs$main_risk_scale) %>%
      join_stage6_flags(stage6_registry = stage6_registry) %>%
      join_horizon_support(horizon_lookup = horizon_lookup)
  }
)

delta_risk_file_target <- make_output_path("stage7_delta_risk.csv")
delta_risk_tbl <- load_or_rebuild_csv_output(
  path = delta_risk_file_target,
  force_rebuild = refresh_patched_exports_even_when_outputs_exist,
  required_cols = c("dataset", "formula_variant", "risk_scale", "support_tier", "horizon_evidence_class", "claim_restriction_flag"),
  builder_fn = function() {
    stage7_objects$delta_risk %>%
      add_common_stage7_metadata(formula_lookup = formula_lookup, risk_scale_default = stage1_inputs$main_risk_scale) %>%
      join_stage6_flags(stage6_registry = stage6_registry) %>%
      join_horizon_support(horizon_lookup = horizon_lookup)
  }
)

threshold_metrics_file_target <- make_output_path("stage7_threshold_metrics.csv")
threshold_metrics_tbl <- load_or_rebuild_csv_output(
  path = threshold_metrics_file_target,
  force_rebuild = refresh_patched_exports_even_when_outputs_exist,
  required_cols = c("dataset", "formula_variant", "risk_scale", "support_tier", "horizon_evidence_class", "claim_restriction_flag", "threshold_projection_suppressed_flag"),
  builder_fn = function() {
    stage7_objects$threshold_metrics %>%
      add_common_stage7_metadata(formula_lookup = formula_lookup, risk_scale_default = stage1_inputs$main_risk_scale) %>%
      join_stage6_flags(stage6_registry = stage6_registry) %>%
      join_horizon_support(horizon_lookup = horizon_lookup) %>%
      apply_threshold_reporting_rules()
  }
)

bootstrap_model_qc_file_target <- make_output_path("stage7_bootstrap_model_qc.csv")
bootstrap_model_qc_tbl <- load_or_rebuild_csv_output(
  path = bootstrap_model_qc_file_target,
  force_rebuild = refresh_patched_exports_even_when_outputs_exist,
  required_cols = c("dataset", "formula_variant", "risk_scale", "site_placement_label"),
  builder_fn = function() {
    stage7_objects$bootstrap_model_qc %>%
      add_common_stage7_metadata(formula_lookup = formula_lookup, risk_scale_default = stage1_inputs$main_risk_scale) %>%
      join_stage6_flags(stage6_registry = stage6_registry)
  }
)

# 🔴 Rebuild: IPCW registry and new Stage 7 support tables ===============================
ipcw_registry_file_target <- make_output_path("stage7_ipcw_registry.csv")
if (should_rebuild_output(ipcw_registry_file_target, force_rebuild = refresh_ipcw_registry_even_when_output_exists, required_cols = c("dataset", "unique_person_id", "horizon_year", "weight_case", "weight_control"))) {
  ipcw_registry_tbl <- rebuild_ipcw_registry(stage1_inputs$analysis_datasets, horizons_year = stage1_inputs$common_horizons_year)
} else {
  ipcw_registry_tbl <- safe_read_csv(ipcw_registry_file_target)
}

hazard_shape_file_target <- make_output_path("stage7_hazard_shape_plausibility.csv")
if (should_rebuild_output(hazard_shape_file_target, force_rebuild = refresh_new_spec_tables_even_when_outputs_exist, required_cols = c("dataset", "formula_variant", "model_id", "hazard_1y", "hazard_10y", "shape_class", "risk_scale"))) {
  hazard_shape_tbl <- build_hazard_shape_plausibility(risk_summary_tbl)
} else {
  hazard_shape_tbl <- safe_read_csv(hazard_shape_file_target)
}

fit_contrast_file_target <- make_output_path("stage7_family_matched_fit_contrast.csv")
if (should_rebuild_output(fit_contrast_file_target, force_rebuild = refresh_new_spec_tables_even_when_outputs_exist, required_cols = c("dataset_key", "family_pair", "LR_2delta_logLik", "lrt_calibration_status", "risk_scale"))) {
  family_matched_fit_contrast_tbl <- build_family_matched_fit_contrast(fit_registry_tbl, stage5_fit_registry_tbl)
} else {
  family_matched_fit_contrast_tbl <- safe_read_csv(fit_contrast_file_target)
}

cure_decomp_file_target <- make_output_path("stage7_cure_supporting_decomposition.csv")
if (should_rebuild_output(cure_decomp_file_target, force_rebuild = refresh_new_spec_tables_even_when_outputs_exist, required_cols = c("dataset", "formula_variant", "model_id", "susceptible_fraction", "risk_scale", "uncured_mean_support_flag"))) {
  cure_supporting_decomposition_tbl <- build_cure_supporting_decomposition(fit_registry_tbl, risk_summary_tbl)
} else {
  cure_supporting_decomposition_tbl <- safe_read_csv(cure_decomp_file_target)
}

refresh_metadata_tbl <- tibble(
  setting = c(
    "requested_run_root_dir",
    "resolved_run_root_dir",
    "requested_export_path",
    "export_path",
    "stage1_data_path",
    "stage6_screening_file",
    "stage5_root_dir",
    "main_risk_scale",
    "common_horizons_year",
    "threshold_plot_horizons_year",
    "refresh_patched_exports_even_when_outputs_exist",
    "refresh_ipcw_registry_even_when_output_exists",
    "refresh_new_spec_tables_even_when_outputs_exist",
    "regenerate_visuals_even_when_outputs_exist",
    "run_bootstrap_extraction",
    "refresh_bootstrap_extraction_even_when_outputs_exist",
    "write_full_metric_values_long_csv",
    "write_individual_plot_pngs"
  ),
  value = c(
    normalize_existing_path(requested_run_root_dir),
    normalize_existing_path(run_root_dir),
    ifelse(is.na(requested_export_path), NA_character_, normalize_existing_path(requested_export_path)),
    normalize_existing_path(export_path),
    normalize_existing_path(stage1_data_path),
    ifelse(is.na(stage6_screening_file), NA_character_, stage6_screening_file),
    ifelse(is.na(stage5_root_dir), NA_character_, stage5_root_dir),
    stage1_inputs$main_risk_scale,
    paste(stage1_inputs$common_horizons_year, collapse = ","),
    paste(threshold_plot_horizons_year, collapse = ","),
    as.character(refresh_patched_exports_even_when_outputs_exist),
    as.character(refresh_ipcw_registry_even_when_output_exists),
    as.character(refresh_new_spec_tables_even_when_outputs_exist),
    as.character(regenerate_visuals_even_when_outputs_exist),
    as.character(run_bootstrap_extraction),
    as.character(refresh_bootstrap_extraction_even_when_outputs_exist),
    as.character(write_full_metric_values_long_csv),
    as.character(write_individual_plot_pngs)
  )
)

# 🔴 Regenerate: visual summary PDF and PNG exports ===============================
visual_plot_registry <- build_stage7_visual_plot_registry(
  risk_summary_tbl = risk_summary_tbl,
  delta_risk_tbl = delta_risk_tbl,
  threshold_metrics_tbl = threshold_metrics_tbl,
  common_horizons_year = stage1_inputs$common_horizons_year,
  threshold_plot_horizons_year = threshold_plot_horizons_year
)

visual_summary_pdf <- make_output_path("stage7_visual_summary.pdf")
visual_plot_manifest <- save_plot_registry_outputs(
  plot_registry = visual_plot_registry,
  combined_pdf_file = visual_summary_pdf,
  force_rebuild_pdf = regenerate_visuals_even_when_outputs_exist,
  write_pngs = write_individual_plot_pngs,
  force_rebuild_pngs = regenerate_visuals_even_when_outputs_exist,
  width = plot_width_in,
  height = plot_height_in,
  dpi = plot_dpi
)

# 🔴 Export: patched Stage 7 tables ===============================
fit_registry_file <- normalize_existing_path(fit_registry_file_target)
coefficients_file <- normalize_existing_path(coefficients_file_target)
subject_predictions_file <- normalize_existing_path(subject_predictions_file_target)
risk_summary_file <- normalize_existing_path(risk_summary_file_target)
delta_risk_file <- normalize_existing_path(delta_risk_file_target)
threshold_metrics_file <- normalize_existing_path(threshold_metrics_file_target)
ipcw_registry_file <- safe_write_csv(ipcw_registry_tbl, ipcw_registry_file_target)
bootstrap_model_qc_file <- normalize_existing_path(bootstrap_model_qc_file_target)
stage6_registry_file <- safe_write_csv(stage6_carry_forward_tbl, make_output_path("stage7_stage6_carry_forward_registry.csv"))
hazard_shape_file <- safe_write_csv(hazard_shape_tbl, hazard_shape_file_target)
family_matched_fit_contrast_file <- safe_write_csv(family_matched_fit_contrast_tbl, fit_contrast_file_target)
cure_supporting_decomposition_file <- safe_write_csv(cure_supporting_decomposition_tbl, cure_decomp_file_target)
refresh_metadata_file <- safe_write_csv(refresh_metadata_tbl, make_output_path("stage7_refresh_metadata.csv"))

# 🔴 Extract: memory-safe bootstrap summaries ===============================
replicate_overview_file <- NA_character_
failure_events_file <- NA_character_
failure_pattern_file <- NA_character_
failure_by_model_file <- NA_character_
metric_distribution_file <- NA_character_
top_metric_cells_file <- NA_character_
top_metric_values_file <- NA_character_
diagnostic_pdf_file <- NA_character_
bootstrap_plot_manifest <- tibble(
  plot_name = character(),
  source_table = character(),
  file_type = character(),
  file_name = character(),
  file_path = character(),
  description = character(),
  exists = logical()
)

bootstrap_expected_csvs <- c(
  make_output_path("stage7_replicate_overview.csv"),
  make_output_path("stage7_replicate_failure_events.csv"),
  make_output_path("stage7_replicate_failure_pattern_summary.csv"),
  make_output_path("stage7_replicate_failure_by_model.csv"),
  make_output_path("stage7_metric_distribution_summary.csv"),
  make_output_path("stage7_top_variable_metric_cells.csv"),
  make_output_path("stage7_top_variable_metric_values_long.csv")
)

if (isTRUE(run_bootstrap_extraction) && file.exists(stage7_objects$merged_bootstrap_rds_file)) {
  existing_bootstrap_csvs_ready <- all(file.exists(bootstrap_expected_csvs))
  
  if (!refresh_bootstrap_extraction_even_when_outputs_exist && existing_bootstrap_csvs_ready) {
    replicate_overview <- safe_read_csv(make_output_path("stage7_replicate_overview.csv"))
    failure_events <- safe_read_csv(make_output_path("stage7_replicate_failure_events.csv"))
    failure_pattern_summary <- safe_read_csv(make_output_path("stage7_replicate_failure_pattern_summary.csv"))
    failure_by_model <- safe_read_csv(make_output_path("stage7_replicate_failure_by_model.csv"))
    metric_distribution_summary <- safe_read_csv(make_output_path("stage7_metric_distribution_summary.csv"))
    top_variable_metric_cells <- safe_read_csv(make_output_path("stage7_top_variable_metric_cells.csv"))
    top_metric_values_long <- safe_read_csv(make_output_path("stage7_top_variable_metric_values_long.csv"))
  } else {
    bootstrap_results <- standardize_bootstrap_results(readRDS(stage7_objects$merged_bootstrap_rds_file))
    
    fit_registry_norm <- normalize_fit_registry(bootstrap_results$fit_registry)
    failure_events <- make_failure_events(fit_registry_norm)
    replicate_overview <- make_replicate_overview(fit_registry_norm)
    failure_pattern_summary <- make_failure_pattern_summary(replicate_overview)
    failure_by_model <- make_failure_by_model(failure_events)
    
    metric_specs <- get_metric_specs(bootstrap_results)
    metric_distribution_summary <- bind_rows(
      lapply(names(metric_specs), function(table_name) {
        spec <- metric_specs[[table_name]]
        summarize_metric_table(spec$df, table_name, spec$metric_cols)
      })
    ) %>%
      mutate(metric_key = make_metric_key(.)) %>%
      arrange(table_name, dataset, formula_variant, model_id, component, term, horizon_year, threshold, metric_name)
    
    top_variable_metric_cells <- metric_distribution_summary %>%
      filter(n_replicates_nonmissing >= 5, is.finite(sd_value), !is.na(sd_value)) %>%
      arrange(desc(sd_value), desc(n_replicates_nonmissing), metric_key) %>%
      slice_head(n = top_n_variable_metric_cells) %>%
      mutate(metric_cell_id = row_number())
    
    top_metric_values_long <- extract_top_metric_values(top_variable_metric_cells, bootstrap_results)
    
    if (isTRUE(write_full_metric_values_long_csv)) {
      write_metric_values_long_streaming(
        bootstrap_results = bootstrap_results,
        out_file = make_output_path("stage7_metric_values_long_full.csv")
      )
    }
    
    replicate_overview_file <- safe_write_csv(replicate_overview, make_output_path("stage7_replicate_overview.csv"))
    failure_events_file <- safe_write_csv(failure_events, make_output_path("stage7_replicate_failure_events.csv"))
    failure_pattern_file <- safe_write_csv(failure_pattern_summary, make_output_path("stage7_replicate_failure_pattern_summary.csv"))
    failure_by_model_file <- safe_write_csv(failure_by_model, make_output_path("stage7_replicate_failure_by_model.csv"))
    metric_distribution_file <- safe_write_csv(metric_distribution_summary, make_output_path("stage7_metric_distribution_summary.csv"))
    top_metric_cells_file <- safe_write_csv(top_variable_metric_cells, make_output_path("stage7_top_variable_metric_cells.csv"))
    top_metric_values_file <- safe_write_csv(top_metric_values_long, make_output_path("stage7_top_variable_metric_values_long.csv"))
  }
  
  if (is.na(replicate_overview_file)) replicate_overview_file <- normalize_existing_path(make_output_path("stage7_replicate_overview.csv"))
  if (is.na(failure_events_file)) failure_events_file <- normalize_existing_path(make_output_path("stage7_replicate_failure_events.csv"))
  if (is.na(failure_pattern_file)) failure_pattern_file <- normalize_existing_path(make_output_path("stage7_replicate_failure_pattern_summary.csv"))
  if (is.na(failure_by_model_file)) failure_by_model_file <- normalize_existing_path(make_output_path("stage7_replicate_failure_by_model.csv"))
  if (is.na(metric_distribution_file)) metric_distribution_file <- normalize_existing_path(make_output_path("stage7_metric_distribution_summary.csv"))
  if (is.na(top_metric_cells_file)) top_metric_cells_file <- normalize_existing_path(make_output_path("stage7_top_variable_metric_cells.csv"))
  if (is.na(top_metric_values_file)) top_metric_values_file <- normalize_existing_path(make_output_path("stage7_top_variable_metric_values_long.csv"))
  
  bootstrap_plot_registry <- build_bootstrap_plot_registry(
    replicate_overview = replicate_overview,
    failure_pattern_summary = failure_pattern_summary,
    failure_by_model = failure_by_model,
    top_metric_values_long = top_metric_values_long
  )
  
  diagnostic_pdf_file <- make_output_path("stage7_bootstrap_extraction_plots.pdf")
  bootstrap_plot_manifest <- save_plot_registry_outputs(
    plot_registry = bootstrap_plot_registry,
    combined_pdf_file = diagnostic_pdf_file,
    force_rebuild_pdf = refresh_bootstrap_extraction_even_when_outputs_exist,
    write_pngs = write_individual_plot_pngs,
    force_rebuild_pngs = refresh_bootstrap_extraction_even_when_outputs_exist,
    width = plot_width_in,
    height = plot_height_in,
    dpi = plot_dpi
  )
}

# 🔴 Assemble: refreshed Stage 7 output manifest ===============================
main_manifest <- tibble(
  file_name = c(
    basename(fit_registry_file),
    basename(coefficients_file),
    basename(subject_predictions_file),
    basename(risk_summary_file),
    basename(delta_risk_file),
    basename(threshold_metrics_file),
    basename(ipcw_registry_file),
    basename(bootstrap_model_qc_file),
    basename(stage6_registry_file),
    basename(hazard_shape_file),
    basename(family_matched_fit_contrast_file),
    basename(cure_supporting_decomposition_file),
    basename(refresh_metadata_file),
    if (file.exists(visual_summary_pdf)) basename(visual_summary_pdf) else NA_character_,
    if (!is.na(replicate_overview_file)) basename(replicate_overview_file) else NA_character_,
    if (!is.na(failure_events_file)) basename(failure_events_file) else NA_character_,
    if (!is.na(failure_pattern_file)) basename(failure_pattern_file) else NA_character_,
    if (!is.na(failure_by_model_file)) basename(failure_by_model_file) else NA_character_,
    if (!is.na(metric_distribution_file)) basename(metric_distribution_file) else NA_character_,
    if (!is.na(top_metric_cells_file)) basename(top_metric_cells_file) else NA_character_,
    if (!is.na(top_metric_values_file)) basename(top_metric_values_file) else NA_character_,
    if (!is.na(diagnostic_pdf_file)) basename(diagnostic_pdf_file) else NA_character_,
    basename(stage7_objects$merged_bootstrap_rds_file)
  ),
  file_path = c(
    fit_registry_file,
    coefficients_file,
    subject_predictions_file,
    risk_summary_file,
    delta_risk_file,
    threshold_metrics_file,
    ipcw_registry_file,
    bootstrap_model_qc_file,
    stage6_registry_file,
    hazard_shape_file,
    family_matched_fit_contrast_file,
    cure_supporting_decomposition_file,
    refresh_metadata_file,
    if (file.exists(visual_summary_pdf)) normalize_existing_path(visual_summary_pdf) else NA_character_,
    replicate_overview_file,
    failure_events_file,
    failure_pattern_file,
    failure_by_model_file,
    metric_distribution_file,
    top_metric_cells_file,
    top_metric_values_file,
    if (!is.na(diagnostic_pdf_file)) normalize_existing_path(diagnostic_pdf_file) else NA_character_,
    normalize_existing_path(stage7_objects$merged_bootstrap_rds_file)
  ),
  description = c(
    "Patched fit registry with Stage 1 metadata and Stage 6 carry-forward flags.",
    "Patched coefficient table with Stage 1 metadata and Stage 6 carry-forward flags.",
    "Patched subject-by-horizon predictions with Stage 1 metadata and Stage 6 carry-forward flags.",
    "Patched risk summary with risk-scale and horizon-support metadata.",
    "Patched delta-risk table with risk-scale and horizon-support metadata.",
    "Patched threshold metrics with horizon-support metadata and PNU projection suppression.",
    "Rebuilt or reused IPCW registry from the Stage 1 analysis datasets.",
    "Patched bootstrap model QC with Stage 1 metadata and Stage 6 carry-forward flags.",
    "Canonical Stage 6 carry-forward table used for joins.",
    "Hazard-shape plausibility table on the annual 1-10 year grid.",
    "Family-matched cure-versus-non-cure fit-contrast table.",
    "Cure-model-only supporting decomposition table with susceptible-fraction and uncured-only summaries when available.",
    "Refresh metadata and rerun-safe settings used for this Stage 7 patch run.",
    "Combined Stage 7 visual summary PDF.",
    "Bootstrap replicate overview.",
    "Bootstrap failure events by replicate and model row.",
    "Bootstrap failure-pattern frequency table.",
    "Bootstrap failure counts aggregated by model.",
    "Across-replicate summary statistics for each metric cell.",
    "Top variable metric cells selected from the summary table.",
    "Replicate-level raw values for the top variable metric cells.",
    "Bootstrap extraction diagnostic PDF.",
    "Saved merged bootstrap results RDS reused without rerunning bootstrap fitting."
  )
) %>%
  filter(!is.na(file_name))

png_manifest <- bind_rows(visual_plot_manifest, bootstrap_plot_manifest) %>%
  filter(file_type == "png" | file_type == "pdf") %>%
  transmute(
    file_name = file_name,
    file_path = file_path,
    description = description
  )

output_manifest <- bind_rows(
  main_manifest,
  if (isTRUE(include_pngs_in_manifest)) png_manifest else png_manifest %>% filter(file_name %in% basename(c(visual_summary_pdf, diagnostic_pdf_file)))
) %>%
  distinct(file_name, .keep_all = TRUE) %>%
  mutate(
    exists = file.exists(file_path),
    size_bytes = ifelse(exists, file.info(file_path)$size, NA_real_)
  )

output_manifest_file <- safe_write_csv(output_manifest, make_output_path("stage7_output_manifest.csv"))

message("✅ Stage 7 revised refresh completed using existing Stage 7 results whenever possible, with PDF and per-plot PNG exports preserved.")
print(output_manifest)
