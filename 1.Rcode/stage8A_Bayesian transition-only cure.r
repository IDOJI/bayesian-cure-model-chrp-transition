# 🔴 Configure: Stage-8A paths, controls, and thresholds ===============================
# 🟧 Configure: OS-aware root path ===============================
sys_name <- Sys.info()[["sysname"]]

project_root <- switch(
  sys_name,
  "Darwin"  = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  "Windows" = "C:/Users/clair/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  stop("Unsupported OS: ", sys_name)
)

# 🟧 Configure: Stage-8A paths, controls, and thresholds ===============================
stage1_export_path <- file.path(project_root, "stage1_Backbone lock")
merged_data_path   <- file.path(project_root, "MERGED_dataset3_pnu_snu.csv")
export_path        <- file.path(project_root, "stage8A_Bayesian transition-only cure")
stage5_export_path <- file.path(project_root, "stage5_Individualized no-cure comparator")
stage6_export_path <- file.path(project_root, "stage6_Cure-appropriateness screening")
stage8b_export_path <- NULL

stage1_bundle_file <- file.path(stage1_export_path, "stage1_backbone_bundle.rds")
stage1_datasets_file <- file.path(stage1_export_path, "stage1_analysis_datasets.rds")
stage1_formula_registry_file <- file.path(stage1_export_path, "stage1_formula_registry.csv")
stage1_horizon_registry_file <- file.path(stage1_export_path, "stage1_horizon_registry.csv")
stage1_threshold_registry_file <- file.path(stage1_export_path, "stage1_threshold_registry.csv")
stage1_modeling_registry_file <- file.path(stage1_export_path, "stage1_modeling_registry.csv")
stage1_metadata_registry_file <- file.path(stage1_export_path, "stage1_metadata_registry.csv")

stage6_screening_flag_csv <- file.path(stage6_export_path, "stage6_carry_forward_flag_table.csv")
stage5_nocure_cohort_csv <- file.path(stage5_export_path, "stage5_model_performance_long.csv")
stage5_nocure_classification_csv <- file.path(stage5_export_path, "stage5_model_performance_long.csv")

pnu_site_label <- "PNU"
snu_site_label <- "SNU"

common_horizons_year_fallback <- 1:10
risk_thresholds_fallback <- c(0.05, 0.10, 0.15, 0.20)

include_merged_incidence_site_supplementary <- TRUE
run_model_ids <- NULL

reuse_existing_stage8_outputs <- TRUE
require_existing_rds_to_skip_fit <- FALSE
preserve_existing_diagnostic_pdf <- FALSE
always_regenerate_plot_bundle <- TRUE

stan_chains <- 4L
stan_iter <- 2000L
stan_warmup <- 1000L
stan_thin <- 1L
stan_seed <- 20260328L
stan_adapt_delta <- 0.95
stan_max_treedepth <- 12L
stan_refresh <- 0L

prior_predictive_draws <- 500L
posterior_prediction_draws <- 400L
save_full_stanfit_rds <- FALSE

ppc_tolerance_abs <- 0.15
ess_min_threshold <- 400
degenerate_draw_fraction_threshold <- 0.90
degenerate_subject_fraction_threshold <- 0.95
tiny_susceptible_prob <- 0.01
huge_susceptible_prob <- 0.99
tiny_median_years <- 0.05
huge_median_years <- 50

plot_width <- 11
plot_height <- 8.5
plot_dpi <- 300

prior_materiality_thresholds <- list(
  risk = 0.03,
  cure_fraction = 0.03,
  false_positive_burden = 0.03,
  FP100 = 3,
  NB = 0.01,
  PPV = 0.05,
  TPR = 0.05
)

# 🔴 Initialize: packages, options, and runtime defaults ===============================
core_packages <- c(
  "readr", "dplyr", "tibble", "tidyr", "ggplot2", "survival", "matrixStats"
)

missing_core_packages <- core_packages[
  !vapply(core_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_core_packages) > 0L) {
  stop(
    "Install required packages before running this script: ",
    paste(missing_core_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(survival)
  library(matrixStats)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

# 🔴 Define: low-level helpers, I/O, and coercion ===============================
## 🟠 Define: file readers, path normalizers, and safe writers ===============================
normalize_dataset_label <- function(x) {
  x_chr <- trimws(as.character(x))
  x_up <- toupper(x_chr)
  out <- ifelse(
    x_up == "MERGED",
    "merged",
    ifelse(x_up == "PNU", "PNU", ifelse(x_up == "SNU", "SNU", x_chr))
  )
  out[is.na(x)] <- NA_character_
  unname(out)
}

normalize_existing_path <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(NA_character_)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

safe_integer <- function(x) {
  suppressWarnings(as.integer(as.character(x)))
}

first_existing_name <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) {
    return(NULL)
  }
  hit[[1L]]
}

read_delimited_or_rds <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(NULL)
  }
  if (!file.exists(path)) {
    stop(sprintf("Input file does not exist: %s", path), call. = FALSE)
  }
  ext <- tolower(tools::file_ext(path))
  path_low <- tolower(path)

  if (grepl("\\.(csv|txt)(\\.gz)?$", path_low)) {
    return(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  }
  if (ext == "rds") {
    return(readRDS(path))
  }
  stop(sprintf("Unsupported input extension for `%s`.", path), call. = FALSE)
}

read_existing_table_if_present <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  read_delimited_or_rds(path)
}

make_temp_output_path <- function(path, tag = "tmp") {
  dir <- dirname(path)
  base <- basename(path)
  stamp <- paste0(
    format(Sys.time(), "%Y%m%d%H%M%S"),
    "_",
    sprintf("%08d", sample.int(99999999, 1L))
  )

  if (grepl("\\.csv\\.gz$", base, ignore.case = TRUE)) {
    base <- sub("\\.csv\\.gz$", paste0("_", tag, "_", stamp, ".csv.gz"), base, ignore.case = TRUE)
  } else if (grepl("\\.csv$", base, ignore.case = TRUE)) {
    base <- sub("\\.csv$", paste0("_", tag, "_", stamp, ".csv"), base, ignore.case = TRUE)
  } else if (grepl("\\.pdf$", base, ignore.case = TRUE)) {
    base <- sub("\\.pdf$", paste0("_", tag, "_", stamp, ".pdf"), base, ignore.case = TRUE)
  } else if (grepl("\\.rds$", base, ignore.case = TRUE)) {
    base <- sub("\\.rds$", paste0("_", tag, "_", stamp, ".rds"), base, ignore.case = TRUE)
  } else if (grepl("\\.png$", base, ignore.case = TRUE)) {
    base <- sub("\\.png$", paste0("_", tag, "_", stamp, ".png"), base, ignore.case = TRUE)
  } else {
    base <- paste0(base, "_", tag, "_", stamp)
  }

  file.path(dir, base)
}

safe_promote_file <- function(tmp_path, final_path) {
  if (!file.exists(tmp_path)) {
    stop(sprintf("Temporary file does not exist: %s", tmp_path), call. = FALSE)
  }
  dir.create(dirname(final_path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(final_path)) {
    unlink(final_path)
  }
  ok <- file.rename(tmp_path, final_path)
  if (!isTRUE(ok)) {
    ok_copy <- file.copy(tmp_path, final_path, overwrite = TRUE)
    unlink(tmp_path)
    if (!isTRUE(ok_copy)) {
      stop(sprintf("Failed to replace file: %s", final_path), call. = FALSE)
    }
  }
  invisible(TRUE)
}

write_csv_preserve_schema <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- make_temp_output_path(path, tag = "tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  readr::write_csv(df, tmp)
  safe_promote_file(tmp, path)
}

save_rds_preserve_schema <- function(obj, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- make_temp_output_path(path, tag = "tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(obj, tmp)
  safe_promote_file(tmp, path)
}

pdf_file_is_usable <- function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  info <- file.info(path)
  isTRUE(!is.na(info$size[[1]]) && info$size[[1]] > 0)
}

## 🟠 Define: table utilities and defensive join helpers ===============================
flatten_scalar_list_col <- function(x) {
  if (!is.list(x)) {
    return(x)
  }

  values <- lapply(x, function(one) {
    if (is.null(one) || length(one) == 0L) {
      return(NA)
    }
    if (inherits(one, "data.frame")) {
      return(paste(unlist(one, recursive = TRUE, use.names = FALSE), collapse = "|"))
    }
    if (length(one) == 1L) {
      return(one[[1L]])
    }
    paste(as.character(unlist(one, recursive = TRUE, use.names = FALSE)), collapse = "|")
  })

  non_missing <- values[!vapply(values, function(one) length(one) == 1L && is.na(one), logical(1))]
  if (length(non_missing) == 0L) {
    return(rep(NA_character_, length(values)))
  }

  if (all(vapply(non_missing, is.logical, logical(1)))) {
    return(vapply(values, function(one) if (is.na(one)) NA else as.logical(one), logical(1)))
  }

  if (all(vapply(non_missing, function(one) is.numeric(one) || is.integer(one), logical(1)))) {
    return(vapply(values, function(one) if (is.na(one)) NA_real_ else as.numeric(one), numeric(1)))
  }

  vapply(values, function(one) if (is.na(one)) NA_character_ else as.character(one), character(1))
}

simplify_scalar_list_cols <- function(df, cols = names(df)) {
  if (is.null(df) || !inherits(df, "data.frame") || length(cols) == 0L) {
    return(df)
  }
  out <- df
  target_cols <- intersect(cols, names(out))
  for (nm in target_cols) {
    out[[nm]] <- flatten_scalar_list_col(out[[nm]])
  }
  out
}

bind_rows_safe <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) == 0L) {
    return(tibble())
  }
  bind_rows(x)
}

nrow_or_zero <- function(df) {
  if (is.null(df) || !inherits(df, "data.frame")) {
    return(0L)
  }
  nrow(df)
}

drop_existing_columns <- function(df, cols) {
  if (is.null(df) || !inherits(df, "data.frame") || length(cols) == 0L) {
    return(df)
  }
  df %>% select(-any_of(cols))
}

left_join_replacing_columns <- function(x, y, by) {
  if (is.null(x) || !inherits(x, "data.frame")) {
    return(x)
  }
  if (is.null(y) || !inherits(y, "data.frame") || nrow(y) == 0L) {
    return(x)
  }
  by_cols <- unname(by)
  if (length(by_cols) == 0L) {
    return(x)
  }
  if (any(!by_cols %in% names(x)) || any(!by_cols %in% names(y))) {
    return(x)
  }
  join_cols <- setdiff(names(y), by_cols)
  x %>% drop_existing_columns(join_cols) %>% left_join(y, by = by)
}

subset_model_table <- function(df, model_id) {
  if (is.null(df) || !inherits(df, "data.frame") || !("model_id" %in% names(df))) {
    return(tibble())
  }
  df %>% filter(.data$model_id == .env$model_id)
}

subset_model_registry_row <- function(df, model_id) {
  out <- subset_model_table(df, model_id)
  if (nrow(out) == 0L) {
    return(NULL)
  }
  out
}

pick_reuse_table <- function(model_id, existing_df = NULL, bundle_df = NULL) {
  existing_sub <- subset_model_table(existing_df, model_id)
  if (nrow(existing_sub) > 0L) {
    return(existing_sub)
  }
  bundle_sub <- subset_model_table(bundle_df, model_id)
  if (nrow(bundle_sub) > 0L) {
    return(bundle_sub)
  }
  if (inherits(bundle_df, "data.frame") && nrow(bundle_df) > 0L && !("model_id" %in% names(bundle_df))) {
    return(bundle_df)
  }
  tibble()
}

make_threshold_key <- function(x) {
  out <- rep("__NA__", length(x))
  idx <- !is.na(x)
  out[idx] <- sprintf("%.10f", as.numeric(x[idx]))
  out
}

format_stage8_progress <- function(done, total) {
  if (is.na(total) || total <= 0L) {
    return("[0/0 | 100.0%]")
  }
  sprintf("[%d/%d | %5.1f%%]", as.integer(done), as.integer(total), 100 * as.numeric(done) / as.numeric(total))
}

format_stage8_number <- function(x, digits = 3L) {
  if (length(x) == 0L || is.na(x) || !is.finite(x)) {
    return("NA")
  }
  formatC(as.numeric(x), format = "f", digits = digits)
}

emit_stage8_progress <- function(done, total, model_id, detail) {
  message(format_stage8_progress(done, total), " ", as.character(model_id), " ", as.character(detail))
}

elapsed_stage8_seconds <- function(start_time) {
  as.numeric(difftime(Sys.time(), start_time, units = "secs"))
}

is_localhost_connection_warning <- function(w) {
  grepl("closing unused connection .*<-localhost:", conditionMessage(w))
}
# 🔴 Define: Stage-1 backbone alignment and fallback dataset preparation ===============================
## 🟠 Define: horizon support governance and Stage-1 coercion helpers ===============================
derive_support_tier <- function(dataset_name, horizon_year) {
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
  stop(sprintf("Unknown dataset `%s` in support-tier derivation.", dataset_name), call. = FALSE)
}

derive_horizon_evidence_class <- function(dataset_name, horizon_year, support_tier = NULL) {
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
  stop(sprintf("Unknown dataset `%s` in horizon-evidence derivation.", dataset_name), call. = FALSE)
}

derive_claim_restriction_flag <- function(horizon_evidence_class) {
  if (is.na(horizon_evidence_class)) {
    return(NA_character_)
  }
  if (identical(horizon_evidence_class, "directly_observed_data_supported")) {
    return("primary_claim_allowed")
  }
  if (identical(horizon_evidence_class, "partly_model_dependent")) {
    return("secondary_or_sensitivity_only")
  }
  "projection_only"
}

derive_interpretation_note <- function(dataset_name, horizon_year, support_tier, horizon_evidence_class, claim_restriction_flag) {
  if (support_tier == "primary_supported" && horizon_evidence_class == "directly_observed_data_supported") {
    return("Primary supported horizon with comparatively direct follow-up support.")
  }
  if (dataset_name == "PNU" && as.integer(horizon_year) == 2L) {
    return("Sensitivity horizon for PNU; partly model-dependent and not for primary claims.")
  }
  if (support_tier == "secondary" && horizon_evidence_class == "partly_model_dependent") {
    return("Secondary horizon with growing tail uncertainty; interpret with model-dependence explicitly acknowledged.")
  }
  if (claim_restriction_flag == "projection_only") {
    return("Projection-dominant horizon; mostly extrapolated and not eligible for primary claims.")
  }
  "Common comparison horizon retained for cross-model comparability."
}

augment_horizon_registry_fields <- function(horizon_registry) {
  if (is.null(horizon_registry) || !inherits(horizon_registry, "data.frame") || nrow(horizon_registry) == 0L) {
    return(tibble())
  }

  out <- horizon_registry %>%
    mutate(
      dataset = normalize_dataset_label(.data$dataset),
      horizon_year = as.integer(safe_numeric(.data$horizon_year))
    ) %>%
    filter(dataset %in% c("PNU", "SNU", "merged"), !is.na(horizon_year))

  if (!"support_tier" %in% names(out)) {
    out$support_tier <- NA_character_
  }
  if (!"horizon_evidence_class" %in% names(out)) {
    out$horizon_evidence_class <- NA_character_
  }
  if (!"claim_restriction_flag" %in% names(out)) {
    out$claim_restriction_flag <- NA_character_
  }
  if (!"interpretation_note" %in% names(out)) {
    out$interpretation_note <- NA_character_
  }
  if (!"interpretation_tier" %in% names(out)) {
    out$interpretation_tier <- NA_character_
  }
  if (!"horizon_support_label" %in% names(out)) {
    out$horizon_support_label <- NA_character_
  }
  if (!"primary_supported_flag" %in% names(out)) {
    out$primary_supported_flag <- NA
  }
  if (!"horizon_days" %in% names(out)) {
    out$horizon_days <- out$horizon_year * 365.25
  }

  out <- out %>%
    rowwise() %>%
    mutate(
      support_tier = {
        current <- as.character(support_tier %||% NA_character_)
        current <- ifelse(is.na(current) || !nzchar(trimws(current)), NA_character_, current)
        if (!is.na(current)) {
          current
        } else {
          legacy_tier <- as.character(interpretation_tier %||% NA_character_)
          legacy_label <- as.character(horizon_support_label %||% NA_character_)
          if (!is.na(legacy_tier) && grepl("primary", legacy_tier, ignore.case = TRUE)) {
            "primary_supported"
          } else if (!is.na(legacy_label) && grepl("primary", legacy_label, ignore.case = TRUE)) {
            "primary_supported"
          } else if (!is.na(legacy_tier) && grepl("secondary", legacy_tier, ignore.case = TRUE)) {
            "secondary"
          } else if (!is.na(legacy_tier) && grepl("sensitivity", legacy_tier, ignore.case = TRUE)) {
            "sensitivity"
          } else if (!is.na(legacy_tier) && grepl("projection", legacy_tier, ignore.case = TRUE)) {
            "projection"
          } else {
            derive_support_tier(dataset, horizon_year)
          }
        }
      },
      horizon_evidence_class = {
        current <- as.character(horizon_evidence_class %||% NA_character_)
        current <- ifelse(is.na(current) || !nzchar(trimws(current)), NA_character_, current)
        if (!is.na(current)) current else derive_horizon_evidence_class(dataset, horizon_year, support_tier)
      },
      claim_restriction_flag = {
        current <- as.character(claim_restriction_flag %||% NA_character_)
        current <- ifelse(is.na(current) || !nzchar(trimws(current)), NA_character_, current)
        if (!is.na(current)) current else derive_claim_restriction_flag(horizon_evidence_class)
      },
      interpretation_tier = {
        current <- as.character(interpretation_tier %||% NA_character_)
        current <- ifelse(is.na(current) || !nzchar(trimws(current)), NA_character_, current)
        if (!is.na(current)) current else ifelse(support_tier == "primary_supported", "primary-supported", support_tier)
      },
      horizon_support_label = {
        current <- as.character(horizon_support_label %||% NA_character_)
        current <- ifelse(is.na(current) || !nzchar(trimws(current)), NA_character_, current)
        if (!is.na(current)) current else support_tier
      },
      primary_supported_flag = ifelse(is.na(primary_supported_flag), support_tier == "primary_supported", as.logical(primary_supported_flag)),
      interpretation_note = {
        current <- as.character(interpretation_note %||% NA_character_)
        current <- ifelse(is.na(current) || !nzchar(trimws(current)), NA_character_, current)
        if (!is.na(current)) current else derive_interpretation_note(dataset, horizon_year, support_tier, horizon_evidence_class, claim_restriction_flag)
      },
      horizon_days = as.numeric(safe_numeric(horizon_days))
    ) %>%
    ungroup() %>%
    distinct(dataset, horizon_year, .keep_all = TRUE) %>%
    arrange(factor(dataset, levels = c("PNU", "SNU", "merged")), horizon_year)

  out
}

build_fallback_formula_registry <- function() {
  bind_rows(
    tibble(
      dataset = c("PNU", "PNU", "SNU", "SNU"),
      dataset_key = c("PNU", "PNU", "SNU", "SNU"),
      formula_name = c("base", "interaction", "base", "interaction"),
      formula_label = c("Base", "Interaction", "Base", "Interaction"),
      formula_rhs = c(
        "age_s + sex_num",
        "age_s + sex_num + age_s:sex_num",
        "age_s + sex_num",
        "age_s + sex_num + age_s:sex_num"
      )
    ),
    tibble(
      dataset = rep("merged", 4L),
      dataset_key = rep("merged", 4L),
      formula_name = c("base", "interaction", "site_added", "site_interaction"),
      formula_label = c("Base", "Interaction", "Site-added", "Site + interaction"),
      formula_rhs = c(
        "age_s + sex_num",
        "age_s + sex_num + age_s:sex_num",
        "age_s + sex_num + site",
        "age_s + sex_num + age_s:sex_num + site"
      )
    )
  ) %>%
    mutate(
      formula_id = paste(dataset, formula_name, sep = "__"),
      formula_full = paste("~", formula_rhs),
      uses_site = grepl("\\bsite\\b", formula_rhs),
      uses_age_sex_interaction = grepl("age_s:sex_num", formula_rhs, fixed = TRUE),
      site_branch = if_else(uses_site, "site_adjusted", "site_free"),
      interaction_branch = if_else(uses_age_sex_interaction, "age_sex_interaction", "no_age_sex_interaction"),
      risk_scale = "transition_only_main",
      formula_scope = "main_transition_only_scale",
      site_term_interpretation = if_else(uses_site, "structural_context_proxy_not_causal_treatment_effect", "not_applicable")
    )
}

build_fallback_threshold_registry <- function(thresholds) {
  tibble(
    threshold = sort(unique(as.numeric(thresholds))),
    threshold_label = paste0(sprintf("%.0f", 100 * sort(unique(as.numeric(thresholds)))), "%")
  )
}

build_fallback_horizon_registry <- function(horizons_year) {
  bind_rows(lapply(c("PNU", "SNU", "merged"), function(ds) {
    tibble(dataset = ds, horizon_year = as.integer(horizons_year))
  })) %>%
    mutate(horizon_days = horizon_year * 365.25) %>%
    augment_horizon_registry_fields()
}

build_fallback_metadata_registry <- function() {
  tibble(
    metadata_group = c("stage", "stage", "risk_scale", "risk_scale"),
    metadata_name = c("stage_name", "stage_role", "main_risk_scale", "supplementary_risk_scale"),
    metadata_value = c("Stage 1 fallback", "Fallback backbone reconstructed directly from merged data.", "transition_only_main", "transition_cif_competing")
  )
}

build_dataset_registry <- function(analysis_datasets) {
  bind_rows(lapply(names(analysis_datasets), function(dataset_name) {
    df <- analysis_datasets[[dataset_name]]
    tibble(
      dataset = dataset_name,
      dataset_key = dataset_name,
      n_rows = nrow(df),
      n_unique_person_id = dplyr::n_distinct(df$unique_person_id),
      n_site_levels = dplyr::n_distinct(df$site),
      site_values = paste(sort(unique(df$site)), collapse = "|"),
      unique_person_id_variable = "unique_person_id",
      unique_person_id_rule = "site + id",
      time_variable = "days_followup",
      reporting_time_variable = "time_year",
      event_definition = "status_num == 1",
      censoring_definition = "status_num %in% c(0, 2)",
      main_risk_scale = "transition_only_main",
      supplementary_risk_scale = "transition_cif_competing"
    )
  }))
}

coerce_stage1_formula_registry <- function(df) {
  out <- if (is.null(df) || !inherits(df, "data.frame") || nrow(df) == 0L) {
    build_fallback_formula_registry()
  } else {
    tibble::as_tibble(df)
  }

  if (!"dataset" %in% names(out)) {
    stop("Stage 1 formula registry must contain `dataset`.", call. = FALSE)
  }
  if (!"formula_name" %in% names(out)) {
    stop("Stage 1 formula registry must contain `formula_name`.", call. = FALSE)
  }

  out <- out %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      dataset_key = if ("dataset_key" %in% names(out)) normalize_dataset_label(dataset_key) else dataset,
      formula_label = if ("formula_label" %in% names(out)) as.character(formula_label) else tools::toTitleCase(gsub("_", " ", formula_name)),
      formula_rhs = if ("formula_rhs" %in% names(out)) as.character(formula_rhs) else NA_character_,
      formula_id = if ("formula_id" %in% names(out)) as.character(formula_id) else paste(dataset, formula_name, sep = "__"),
      formula_full = if ("formula_full" %in% names(out)) as.character(formula_full) else ifelse(!is.na(formula_rhs), paste("~", formula_rhs), NA_character_),
      uses_site = if ("uses_site" %in% names(out)) as.logical(uses_site) else grepl("\\bsite\\b", formula_rhs),
      uses_age_sex_interaction = if ("uses_age_sex_interaction" %in% names(out)) as.logical(uses_age_sex_interaction) else grepl("age_s:sex_num", formula_rhs, fixed = TRUE),
      site_branch = if ("site_branch" %in% names(out)) as.character(site_branch) else ifelse(uses_site, "site_adjusted", "site_free"),
      interaction_branch = if ("interaction_branch" %in% names(out)) as.character(interaction_branch) else ifelse(uses_age_sex_interaction, "age_sex_interaction", "no_age_sex_interaction"),
      risk_scale = if ("risk_scale" %in% names(out)) as.character(risk_scale) else "transition_only_main",
      formula_scope = if ("formula_scope" %in% names(out)) as.character(formula_scope) else "main_transition_only_scale",
      site_term_interpretation = if ("site_term_interpretation" %in% names(out)) as.character(site_term_interpretation) else ifelse(uses_site, "structural_context_proxy_not_causal_treatment_effect", "not_applicable")
    ) %>%
    distinct(dataset, formula_name, .keep_all = TRUE) %>%
    arrange(factor(dataset, levels = c("PNU", "SNU", "merged")), formula_name)

  out
}

coerce_stage1_threshold_registry <- function(df, thresholds_fallback) {
  out <- if (is.null(df) || !inherits(df, "data.frame") || nrow(df) == 0L) {
    build_fallback_threshold_registry(thresholds_fallback)
  } else {
    tibble::as_tibble(df)
  }

  threshold_col <- first_existing_name(out, c("threshold", "risk_threshold"))
  if (is.null(threshold_col)) {
    out <- build_fallback_threshold_registry(thresholds_fallback)
  } else {
    out <- out %>% mutate(threshold = as.numeric(safe_numeric(.data[[threshold_col]])))
    if (!"threshold_label" %in% names(out)) {
      out$threshold_label <- paste0(sprintf("%.0f", 100 * out$threshold), "%")
    }
    out <- out %>% filter(!is.na(threshold)) %>% distinct(threshold, .keep_all = TRUE) %>% arrange(threshold)
  }
  out
}

coerce_stage1_horizon_registry <- function(df, horizons_fallback) {
  out <- if (is.null(df) || !inherits(df, "data.frame") || nrow(df) == 0L) {
    build_fallback_horizon_registry(horizons_fallback)
  } else {
    tibble::as_tibble(df)
  }

  dataset_col <- first_existing_name(out, c("dataset", "dataset_key"))
  horizon_col <- first_existing_name(out, c("horizon_year", "horizon", "year"))
  if (is.null(dataset_col) || is.null(horizon_col)) {
    out <- build_fallback_horizon_registry(horizons_fallback)
  } else {
    out <- out %>% rename(dataset = !!dataset_col, horizon_year = !!horizon_col)
    out <- augment_horizon_registry_fields(out)
  }
  out
}

## 🟠 Define: fallback raw-data preparation when Stage-1 artifacts are absent ===============================
standardize_known_site_labels <- function(df, pnu_label, snu_label) {
  if (!"site" %in% names(df)) {
    stop("Merged input must contain column `site`.", call. = FALSE)
  }
  df %>%
    mutate(
      site = trimws(as.character(site)),
      site = case_when(
        toupper(site) == toupper(pnu_label) ~ pnu_label,
        toupper(site) == toupper(snu_label) ~ snu_label,
        TRUE ~ site
      )
    )
}

prepare_analysis_dataset <- function(df, dataset_name, pnu_label, snu_label) {
  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(sprintf("[%s] Missing required columns: %s", dataset_name, paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  out <- df %>%
    mutate(
      id = trimws(as.character(id)),
      site = trimws(as.character(site)),
      site = case_when(
        toupper(site) == toupper(pnu_label) ~ pnu_label,
        toupper(site) == toupper(snu_label) ~ snu_label,
        TRUE ~ site
      ),
      sex_num = as.integer(safe_numeric(sex_num)),
      age_exact_entry = safe_numeric(age_exact_entry),
      days_followup = safe_numeric(days_followup),
      status_num = as.integer(safe_numeric(status_num))
    )

  if (nrow(out) == 0L) {
    stop(sprintf("[%s] Dataset has zero rows.", dataset_name), call. = FALSE)
  }
  if (anyNA(out[required_cols])) {
    stop(sprintf("[%s] Missing values detected in required backbone columns.", dataset_name), call. = FALSE)
  }
  if (any(out$id == "")) {
    stop(sprintf("[%s] Blank IDs detected.", dataset_name), call. = FALSE)
  }
  if (any(!out$sex_num %in% c(0L, 1L))) {
    stop(sprintf("[%s] `sex_num` must be coded 0/1.", dataset_name), call. = FALSE)
  }
  if (any(!out$status_num %in% c(0L, 1L, 2L))) {
    stop(sprintf("[%s] `status_num` must be coded 0/1/2.", dataset_name), call. = FALSE)
  }
  if (any(out$days_followup < 0)) {
    stop(sprintf("[%s] Negative `days_followup` values detected.", dataset_name), call. = FALSE)
  }

  if (dataset_name != "merged") {
    out <- out %>% mutate(site = dataset_name)
  }

  age_sd <- stats::sd(out$age_exact_entry)
  if (is.na(age_sd) || age_sd <= 0) {
    stop(sprintf("[%s] `age_exact_entry` must have positive SD.", dataset_name), call. = FALSE)
  }
  age_mean <- mean(out$age_exact_entry)

  out <- out %>%
    mutate(
      unique_person_id = paste(site, id, sep = "_"),
      time_year = pmax(days_followup / 365.25, 1e-8),
      event_main = as.integer(status_num == 1L),
      right_censor_flag = as.integer(status_num == 0L),
      censor_main = as.integer(status_num %in% c(0L, 2L)),
      remission_flag = as.integer(status_num == 2L),
      sex_label = factor(ifelse(sex_num == 0L, "Female", "Male"), levels = c("Female", "Male")),
      status_label = factor(
        case_when(
          status_num == 0L ~ "right_censoring",
          status_num == 2L ~ "remission",
          TRUE ~ "transition"
        ),
        levels = c("right_censoring", "remission", "transition")
      ),
      age_s = (age_exact_entry - age_mean) / (2 * age_sd)
    )

  if (nrow(out) != dplyr::n_distinct(out$unique_person_id)) {
    stop(sprintf("[%s] `site + id` is not unique.", dataset_name), call. = FALSE)
  }
  if (dataset_name == "merged" && dplyr::n_distinct(out$site) < 2L) {
    stop("[merged] Merged dataset must contain at least two site levels.", call. = FALSE)
  }

  scaling_row <- tibble(
    dataset = dataset_name,
    dataset_key = dataset_name,
    variable = "age_exact_entry",
    scaled_variable = "age_s",
    center_mean = age_mean,
    scale_sd = age_sd,
    scale_two_sd = 2 * age_sd,
    scaling_rule = "(age_exact_entry - mean(age_exact_entry)) / (2 * sd(age_exact_entry))"
  )

  list(data = out, scaling = scaling_row)
}

load_stage8_datasets_from_merged <- function(merged_path, pnu_label, snu_label) {
  merged_raw <- read_delimited_or_rds(merged_path)
  if (is.null(merged_raw)) {
    stop("Provide `merged_data_path` when Stage 1 artifacts are unavailable.", call. = FALSE)
  }
  merged_raw <- standardize_known_site_labels(merged_raw, pnu_label, snu_label)
  pnu_raw <- merged_raw %>% filter(site == pnu_label)
  snu_raw <- merged_raw %>% filter(site == snu_label)
  if (nrow(pnu_raw) == 0L) {
    stop(sprintf("No rows found for PNU site label `%s` in merged data.", pnu_label), call. = FALSE)
  }
  if (nrow(snu_raw) == 0L) {
    stop(sprintf("No rows found for SNU site label `%s` in merged data.", snu_label), call. = FALSE)
  }
  prepared_pnu <- prepare_analysis_dataset(pnu_raw, "PNU", pnu_label, snu_label)
  prepared_snu <- prepare_analysis_dataset(snu_raw, "SNU", pnu_label, snu_label)
  prepared_merged <- prepare_analysis_dataset(merged_raw, "merged", pnu_label, snu_label)
  list(
    datasets = list(PNU = prepared_pnu$data, SNU = prepared_snu$data, merged = prepared_merged$data),
    scaling_registry = bind_rows(prepared_pnu$scaling, prepared_snu$scaling, prepared_merged$scaling)
  )
}

load_stage1_backbone_or_fallback <- function() {
  stage1_bundle <- NULL
  stage1_loaded_flag <- FALSE

  if (file.exists(stage1_bundle_file) || file.exists(stage1_datasets_file)) {
    if (file.exists(stage1_bundle_file)) {
      stage1_bundle <- tryCatch(readRDS(stage1_bundle_file), error = function(e) NULL)
    }

    analysis_datasets <- if (file.exists(stage1_datasets_file)) {
      tryCatch(readRDS(stage1_datasets_file), error = function(e) NULL)
    } else {
      NULL
    }

    if (is.null(analysis_datasets) && !is.null(stage1_bundle) && is.list(stage1_bundle$datasets)) {
      analysis_datasets <- stage1_bundle$datasets
    }

    if (!is.null(analysis_datasets) && all(c("PNU", "SNU", "merged") %in% names(analysis_datasets))) {
      analysis_datasets <- lapply(analysis_datasets[c("PNU", "SNU", "merged")], tibble::as_tibble)
      stage1_loaded_flag <- TRUE

      formula_registry <- read_existing_table_if_present(stage1_formula_registry_file)
      if (is.null(formula_registry) && !is.null(stage1_bundle)) {
        formula_registry <- stage1_bundle$registries$formula_registry %||% NULL
      }
      horizon_registry <- read_existing_table_if_present(stage1_horizon_registry_file)
      if (is.null(horizon_registry) && !is.null(stage1_bundle)) {
        horizon_registry <- stage1_bundle$registries$horizon_registry %||% NULL
      }
      threshold_registry <- read_existing_table_if_present(stage1_threshold_registry_file)
      if (is.null(threshold_registry) && !is.null(stage1_bundle)) {
        threshold_registry <- stage1_bundle$registries$threshold_registry %||% NULL
      }
      modeling_registry <- read_existing_table_if_present(stage1_modeling_registry_file)
      if (is.null(modeling_registry) && !is.null(stage1_bundle)) {
        modeling_registry <- stage1_bundle$registries$modeling_registry %||% NULL
      }
      metadata_registry <- read_existing_table_if_present(stage1_metadata_registry_file)
      if (is.null(metadata_registry) && !is.null(stage1_bundle)) {
        metadata_registry <- stage1_bundle$registries$metadata_registry %||% NULL
      }
      scaling_registry <- if (!is.null(stage1_bundle)) stage1_bundle$registries$scaling_registry %||% NULL else NULL
      if (is.null(scaling_registry)) {
        scaling_registry <- bind_rows(lapply(names(analysis_datasets), function(ds) {
          tibble(
            dataset = ds,
            dataset_key = ds,
            variable = "age_exact_entry",
            scaled_variable = "age_s",
            center_mean = mean(analysis_datasets[[ds]]$age_exact_entry),
            scale_sd = sd(analysis_datasets[[ds]]$age_exact_entry),
            scale_two_sd = 2 * sd(analysis_datasets[[ds]]$age_exact_entry),
            scaling_rule = "(age_exact_entry - mean(age_exact_entry)) / (2 * sd(age_exact_entry))"
          )
        }))
      }

      formula_registry <- coerce_stage1_formula_registry(formula_registry)
      threshold_registry <- coerce_stage1_threshold_registry(threshold_registry, risk_thresholds_fallback)
      horizon_registry <- coerce_stage1_horizon_registry(horizon_registry, common_horizons_year_fallback)
      metadata_registry <- if (is.null(metadata_registry) || !inherits(metadata_registry, "data.frame") || nrow(metadata_registry) == 0L) build_fallback_metadata_registry() else tibble::as_tibble(metadata_registry)
      modeling_registry <- if (is.null(modeling_registry) || !inherits(modeling_registry, "data.frame")) tibble() else tibble::as_tibble(modeling_registry)
      dataset_registry <- build_dataset_registry(analysis_datasets)

      return(list(
        stage1_loaded_flag = TRUE,
        stage1_bundle = stage1_bundle,
        analysis_datasets = analysis_datasets,
        scaling_registry = scaling_registry,
        formula_registry = formula_registry,
        horizon_registry = horizon_registry,
        threshold_registry = threshold_registry,
        modeling_registry = modeling_registry,
        metadata_registry = metadata_registry,
        dataset_registry = dataset_registry
      ))
    }
  }

  fallback_loaded <- load_stage8_datasets_from_merged(merged_data_path, pnu_site_label, snu_site_label)
  analysis_datasets <- fallback_loaded$datasets
  scaling_registry <- fallback_loaded$scaling_registry
  formula_registry <- build_fallback_formula_registry()
  threshold_registry <- build_fallback_threshold_registry(risk_thresholds_fallback)
  horizon_registry <- build_fallback_horizon_registry(common_horizons_year_fallback)
  metadata_registry <- build_fallback_metadata_registry()
  dataset_registry <- build_dataset_registry(analysis_datasets)

  list(
    stage1_loaded_flag = FALSE,
    stage1_bundle = NULL,
    analysis_datasets = analysis_datasets,
    scaling_registry = scaling_registry,
    formula_registry = formula_registry,
    horizon_registry = horizon_registry,
    threshold_registry = threshold_registry,
    modeling_registry = tibble(),
    metadata_registry = metadata_registry,
    dataset_registry = dataset_registry
  )
}

# 🔴 Define: Stage-5 / Stage-6 linkage readers and canonical carry-forward fields ===============================
## 🟠 Define: Stage-5 and Stage-6 normalization helpers ===============================
extract_stage5_performance_table <- function(obj) {
  if (is.null(obj)) {
    return(NULL)
  }
  if (inherits(obj, "data.frame")) {
    return(obj)
  }
  if (is.list(obj)) {
    if (inherits(obj$model_performance_long, "data.frame")) {
      return(obj$model_performance_long)
    }
    if (is.list(obj$outputs) && inherits(obj$outputs$model_performance_long, "data.frame")) {
      return(obj$outputs$model_performance_long)
    }
  }
  NULL
}

extract_stage6_screening_table <- function(obj) {
  if (is.null(obj)) {
    return(NULL)
  }
  if (inherits(obj, "data.frame")) {
    return(obj)
  }
  if (is.list(obj)) {
    if (inherits(obj$carry_forward_flag_table, "data.frame")) {
      return(obj$carry_forward_flag_table)
    }
    if (inherits(obj$screening_summary, "data.frame")) {
      return(obj$screening_summary)
    }
    if (is.list(obj$outputs) && inherits(obj$outputs$carry_forward_flag_table, "data.frame")) {
      return(obj$outputs$carry_forward_flag_table)
    }
    if (is.list(obj$outputs) && inherits(obj$outputs$screening_summary, "data.frame")) {
      return(obj$outputs$screening_summary)
    }
  }
  NULL
}

map_stage6_variant_to_formula_anchor <- function(dataset_key = NA_character_, source_dataset = NA_character_, analysis_variant = NA_character_, hsu_formula_branch = NA_character_) {
  dataset_key <- as.character(dataset_key %||% NA_character_)
  source_dataset <- normalize_dataset_label(source_dataset %||% NA_character_)
  analysis_variant <- as.character(analysis_variant %||% NA_character_)
  hsu_formula_branch <- as.character(hsu_formula_branch %||% NA_character_)

  if (!is.na(dataset_key) && nzchar(dataset_key)) {
    if (identical(dataset_key, "merged__site_free")) {
      return(c("base", "interaction"))
    }
    if (identical(dataset_key, "merged__site_adjusted")) {
      return(c("site_added", "site_interaction"))
    }
    if (dataset_key %in% c("PNU", "SNU", "merged")) {
      return("ALL")
    }
  }

  if (identical(source_dataset, "merged")) {
    branch_text <- paste(analysis_variant, hsu_formula_branch)
    if (grepl("site_adjusted", branch_text, fixed = TRUE)) {
      return(c("site_added", "site_interaction"))
    }
    if (grepl("site_free", branch_text, fixed = TRUE)) {
      return(c("base", "interaction"))
    }
  }
  "ALL"
}

read_screening_flags <- function(path) {
  empty_out <- tibble(
    dataset = character(),
    formula_anchor = character(),
    cure_model_eligibility_flag = character(),
    primary_gate_method = character(),
    primary_gate_flag = character(),
    receus_primary_class = character(),
    presence_modifier_flag = character(),
    cure_presence_support_flag = character(),
    followup_contradiction_flag = character(),
    followup_not_contradicted_flag = character(),
    screening_note = character()
  )

  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(empty_out)
  }

  raw_obj <- read_delimited_or_rds(path)
  df <- extract_stage6_screening_table(raw_obj)
  if (is.null(df) || !inherits(df, "data.frame")) {
    return(empty_out)
  }
  df <- tibble::as_tibble(df)

  if (("dataset_key" %in% names(df)) && ("cure_model_eligibility_flag" %in% names(df) || "final_decision_flag" %in% names(df))) {
    dataset_source_col <- first_existing_name(df, c("source_dataset", "dataset", "cohort"))
    dataset_key_col <- first_existing_name(df, c("dataset_key"))
    analysis_variant_col <- first_existing_name(df, c("analysis_variant"))
    branch_col <- first_existing_name(df, c("hsu_formula_branch", "analysis_variant"))

    out_stage6 <- bind_rows(lapply(seq_len(nrow(df)), function(i) {
      one <- df[i, , drop = FALSE]
      dataset_value <- if (!is.null(dataset_source_col)) one[[dataset_source_col]][[1L]] else NA_character_
      dataset_std <- normalize_dataset_label(dataset_value)
      formula_anchors <- map_stage6_variant_to_formula_anchor(
        dataset_key = if (!is.null(dataset_key_col)) one[[dataset_key_col]][[1L]] else NA_character_,
        source_dataset = dataset_std,
        analysis_variant = if (!is.null(analysis_variant_col)) one[[analysis_variant_col]][[1L]] else NA_character_,
        hsu_formula_branch = if (!is.null(branch_col)) one[[branch_col]][[1L]] else NA_character_
      )
      tibble(
        dataset = dataset_std,
        formula_anchor = formula_anchors,
        cure_model_eligibility_flag = as.character((one$cure_model_eligibility_flag %||% one$final_decision_flag)[[1L]]),
        primary_gate_method = as.character((one$primary_gate_method %||% NA_character_)[[1L]]),
        primary_gate_flag = as.character((one$primary_gate_flag %||% NA_character_)[[1L]]),
        receus_primary_class = as.character((one$receus_primary_class %||% NA_character_)[[1L]]),
        presence_modifier_flag = as.character((one$presence_modifier_flag %||% NA_character_)[[1L]]),
        cure_presence_support_flag = as.character((one$cure_presence_support_flag %||% NA_character_)[[1L]]),
        followup_contradiction_flag = as.character((one$followup_contradiction_flag %||% NA_character_)[[1L]]),
        followup_not_contradicted_flag = as.character((one$followup_not_contradicted_flag %||% NA_character_)[[1L]]),
        screening_note = as.character((one$screening_note %||% one$screening_detail %||% NA_character_)[[1L]])
      )
    })) %>%
      filter(dataset %in% c("PNU", "SNU", "merged"))

    if (nrow(out_stage6) > 0L) {
      return(out_stage6)
    }
  }

  dataset_col <- first_existing_name(df, c("dataset", "cohort", "source_dataset"))
  flag_col <- first_existing_name(df, c("cure_model_eligibility_flag", "final_decision_flag", "decision_flag", "screening_flag"))
  note_col <- first_existing_name(df, c("screening_note", "screening_detail", "detail"))
  formula_col <- first_existing_name(df, c("formula_anchor", "formula_name", "formula_type"))

  if (is.null(dataset_col) || is.null(flag_col)) {
    return(empty_out)
  }

  out <- df %>%
    transmute(
      dataset = normalize_dataset_label(.data[[dataset_col]]),
      formula_anchor = if (!is.null(formula_col)) as.character(.data[[formula_col]]) else "ALL",
      cure_model_eligibility_flag = as.character(.data[[flag_col]]),
      primary_gate_method = if ("primary_gate_method" %in% names(df)) as.character(.data$primary_gate_method) else NA_character_,
      primary_gate_flag = if ("primary_gate_flag" %in% names(df)) as.character(.data$primary_gate_flag) else NA_character_,
      receus_primary_class = if ("receus_primary_class" %in% names(df)) as.character(.data$receus_primary_class) else NA_character_,
      presence_modifier_flag = if ("presence_modifier_flag" %in% names(df)) as.character(.data$presence_modifier_flag) else NA_character_,
      cure_presence_support_flag = if ("cure_presence_support_flag" %in% names(df)) as.character(.data$cure_presence_support_flag) else NA_character_,
      followup_contradiction_flag = if ("followup_contradiction_flag" %in% names(df)) as.character(.data$followup_contradiction_flag) else NA_character_,
      followup_not_contradicted_flag = if ("followup_not_contradicted_flag" %in% names(df)) as.character(.data$followup_not_contradicted_flag) else NA_character_,
      screening_note = if (!is.null(note_col)) as.character(.data[[note_col]]) else NA_character_
    ) %>%
    filter(dataset %in% c("PNU", "SNU", "merged"))

  if (nrow(out) == 0L) empty_out else out
}

build_screening_model_lookup <- function(screening_flags, model_grid) {
  model_base <- model_grid %>% distinct(model_id, dataset, formula_anchor)
  if (is.null(screening_flags) || nrow(screening_flags) == 0L) {
    return(model_base %>% mutate(
      cure_model_eligibility_flag = NA_character_,
      primary_gate_method = NA_character_,
      primary_gate_flag = NA_character_,
      receus_primary_class = NA_character_,
      presence_modifier_flag = NA_character_,
      cure_presence_support_flag = NA_character_,
      followup_contradiction_flag = NA_character_,
      followup_not_contradicted_flag = NA_character_,
      screening_note = NA_character_
    ))
  }

  exact_lookup <- model_base %>%
    left_join(
      screening_flags %>%
        filter(formula_anchor != "ALL") %>%
        distinct(dataset, formula_anchor, .keep_all = TRUE),
      by = c("dataset", "formula_anchor")
    )

  all_lookup <- model_base %>%
    left_join(
      screening_flags %>%
        filter(formula_anchor == "ALL") %>%
        distinct(dataset, .keep_all = TRUE) %>%
        rename_with(~paste0(.x, "_all"), -dataset),
      by = "dataset"
    )

  exact_lookup %>%
    left_join(select(all_lookup, -dataset, -formula_anchor), by = "model_id") %>%
    transmute(
      model_id = model_id,
      cure_model_eligibility_flag = dplyr::coalesce(cure_model_eligibility_flag, cure_model_eligibility_flag_all),
      primary_gate_method = dplyr::coalesce(primary_gate_method, primary_gate_method_all),
      primary_gate_flag = dplyr::coalesce(primary_gate_flag, primary_gate_flag_all),
      receus_primary_class = dplyr::coalesce(receus_primary_class, receus_primary_class_all),
      presence_modifier_flag = dplyr::coalesce(presence_modifier_flag, presence_modifier_flag_all),
      cure_presence_support_flag = dplyr::coalesce(cure_presence_support_flag, cure_presence_support_flag_all),
      followup_contradiction_flag = dplyr::coalesce(followup_contradiction_flag, followup_contradiction_flag_all),
      followup_not_contradicted_flag = dplyr::coalesce(followup_not_contradicted_flag, followup_not_contradicted_flag_all),
      screening_note = dplyr::coalesce(screening_note, screening_note_all)
    ) %>%
    mutate(
      screening_flag = cure_model_eligibility_flag,
      screening_detail = screening_note,
      carry_forward_stage8 = ifelse(is.na(cure_model_eligibility_flag), NA, cure_model_eligibility_flag != "unsupportive")
    )
}

normalize_nocure_cohort <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  raw_obj <- read_delimited_or_rds(path)
  df <- extract_stage5_performance_table(raw_obj)
  if (is.null(df) || !inherits(df, "data.frame")) {
    return(NULL)
  }
  df <- tibble::as_tibble(df)

  if (all(c("metric_domain", "metric_name", "metric_value") %in% names(df))) {
    out <- df %>%
      filter(model_class == "non_cure", metric_domain == "horizon_summary", metric_name == "mean_predicted_risk") %>%
      transmute(
        dataset = normalize_dataset_label(dataset),
        formula_anchor = if ("formula_name" %in% names(df)) as.character(formula_name) else "ALL",
        no_cure_model_id = as.character(model_id),
        horizon_year = as.integer(safe_numeric(horizon_year)),
        metric = "meanRisk",
        value = safe_numeric(metric_value)
      ) %>%
      filter(!is.na(horizon_year), !is.na(value))
    if (nrow(out) > 0L) {
      return(out)
    }
  }

  dataset_col <- first_existing_name(df, c("dataset", "cohort"))
  horizon_col <- first_existing_name(df, c("horizon_year", "horizon", "year"))
  formula_col <- first_existing_name(df, c("formula_anchor", "formula_name", "formula_type"))
  model_col <- first_existing_name(df, c("no_cure_model_id", "model_id", "family", "model"))
  risk_col <- first_existing_name(df, c("meanRisk_no_cure", "meanRisk", "mean_risk", "risk_mean"))

  if (is.null(dataset_col) || is.null(horizon_col) || is.null(risk_col)) {
    return(NULL)
  }

  tibble(
    dataset = normalize_dataset_label(df[[dataset_col]]),
    formula_anchor = if (!is.null(formula_col)) as.character(df[[formula_col]]) else "ALL",
    no_cure_model_id = if (!is.null(model_col)) as.character(df[[model_col]]) else "NO_CURE_REFERENCE",
    horizon_year = as.integer(safe_numeric(df[[horizon_col]])),
    metric = "meanRisk",
    value = safe_numeric(df[[risk_col]])
  ) %>%
    filter(!is.na(horizon_year), !is.na(value))
}

normalize_nocure_classification <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  raw_obj <- read_delimited_or_rds(path)
  df <- extract_stage5_performance_table(raw_obj)
  if (is.null(df) || !inherits(df, "data.frame")) {
    return(NULL)
  }
  df <- tibble::as_tibble(df)

  if (all(c("metric_domain", "metric_name", "metric_value") %in% names(df))) {
    metric_name_map <- c(
      false_positive_burden_nonevents = "false_positive_burden",
      false_positive_per_100 = "FP100",
      net_benefit = "NB",
      PPV = "PPV",
      TPR = "TPR",
      FPR = "FPR"
    )
    out <- df %>%
      filter(model_class == "non_cure", metric_domain == "threshold_summary", metric_name %in% names(metric_name_map)) %>%
      transmute(
        dataset = normalize_dataset_label(dataset),
        formula_anchor = if ("formula_name" %in% names(df)) as.character(formula_name) else "ALL",
        no_cure_model_id = as.character(model_id),
        horizon_year = as.integer(safe_numeric(horizon_year)),
        threshold = as.numeric(safe_numeric(threshold)),
        metric = unname(metric_name_map[metric_name]),
        value = safe_numeric(metric_value)
      ) %>%
      filter(!is.na(horizon_year), !is.na(threshold), !is.na(value))
    if (nrow(out) > 0L) {
      return(out)
    }
  }

  dataset_col <- first_existing_name(df, c("dataset", "cohort"))
  horizon_col <- first_existing_name(df, c("horizon_year", "horizon", "year"))
  threshold_col <- first_existing_name(df, c("threshold", "risk_threshold"))
  formula_col <- first_existing_name(df, c("formula_anchor", "formula_name", "formula_type"))
  model_col <- first_existing_name(df, c("no_cure_model_id", "model_id", "family", "model"))

  metric_candidates <- c(
    false_positive_burden = first_existing_name(df, c("false_positive_burden", "fp_burden")),
    FP100 = first_existing_name(df, c("FP100", "fp100", "false_positives_per_100")),
    NB = first_existing_name(df, c("NB", "net_benefit", "nb_mean")),
    PPV = first_existing_name(df, c("PPV", "ppv")),
    TPR = first_existing_name(df, c("TPR", "tpr")),
    FPR = first_existing_name(df, c("FPR", "fpr"))
  )
  metric_candidates <- metric_candidates[!vapply(metric_candidates, is.null, logical(1))]

  if (is.null(dataset_col) || is.null(horizon_col) || is.null(threshold_col) || length(metric_candidates) == 0L) {
    return(NULL)
  }

  out_base <- tibble(
    dataset = normalize_dataset_label(df[[dataset_col]]),
    formula_anchor = if (!is.null(formula_col)) as.character(df[[formula_col]]) else "ALL",
    no_cure_model_id = if (!is.null(model_col)) as.character(df[[model_col]]) else "NO_CURE_REFERENCE",
    horizon_year = as.integer(safe_numeric(df[[horizon_col]])),
    threshold = as.numeric(safe_numeric(df[[threshold_col]]))
  )

  bind_rows(lapply(names(metric_candidates), function(metric_name) {
    tibble(
      dataset = out_base$dataset,
      formula_anchor = out_base$formula_anchor,
      no_cure_model_id = out_base$no_cure_model_id,
      horizon_year = out_base$horizon_year,
      threshold = out_base$threshold,
      metric = metric_name,
      value = safe_numeric(df[[metric_candidates[[metric_name]]]])
    )
  })) %>%
    filter(!is.na(horizon_year), !is.na(threshold), !is.na(value))
}
# 🔴 Define: model-grid construction, reuse bundles, and prior encoding ===============================
## 🟠 Define: RDS reuse extraction and compatibility checks ===============================
empty_reuse_bundle <- function() {
  list(
    coefficient_summary = NULL,
    diagnostics_parameter_level = NULL,
    ppc_summary = NULL,
    posterior_subject_profile = NULL,
    posterior_subject_yearly = NULL,
    posterior_cohort_yearly = NULL,
    posterior_classification = NULL,
    prior_predictive_summary = NULL,
    selected_parameter_draws = NULL,
    fit_object = NULL,
    rds_path = NA_character_
  )
}

get_model_rds_candidates <- function(model_id, export_dir, registry_row = NULL) {
  cands <- c(
    file.path(export_dir, paste0(model_id, "__bayes_stage8_fit.rds")),
    if (!is.null(registry_row) && "rds_path" %in% names(registry_row)) as.character(registry_row$rds_path[[1L]]) else NA_character_
  )
  cands <- unique(cands)
  cands <- cands[!is.na(cands) & nzchar(cands)]
  unname(cands)
}

extract_stage8_bundle_tables <- function(obj) {
  out <- empty_reuse_bundle()
  if (is.null(obj)) {
    return(out)
  }
  if (inherits(obj, "stanfit")) {
    out$fit_object <- obj
    return(out)
  }
  if (is.list(obj)) {
    if (inherits(obj$coefficient_summary, "data.frame")) out$coefficient_summary <- obj$coefficient_summary
    if (inherits(obj$diagnostics_parameter_level, "data.frame")) out$diagnostics_parameter_level <- obj$diagnostics_parameter_level
    if (inherits(obj$ppc_summary, "data.frame")) out$ppc_summary <- obj$ppc_summary
    if (inherits(obj$posterior_subject_profile, "data.frame")) out$posterior_subject_profile <- obj$posterior_subject_profile
    if (inherits(obj$posterior_subject_yearly, "data.frame")) out$posterior_subject_yearly <- obj$posterior_subject_yearly
    if (inherits(obj$posterior_cohort_yearly, "data.frame")) out$posterior_cohort_yearly <- obj$posterior_cohort_yearly
    if (inherits(obj$posterior_classification, "data.frame")) out$posterior_classification <- obj$posterior_classification
    if (inherits(obj$prior_predictive_summary, "data.frame")) out$prior_predictive_summary <- obj$prior_predictive_summary
    if (inherits(obj$selected_parameter_draws, "data.frame")) out$selected_parameter_draws <- obj$selected_parameter_draws
    if (inherits(obj$fit, "stanfit")) out$fit_object <- obj$fit
    if (inherits(obj$stanfit, "stanfit")) out$fit_object <- obj$stanfit
  }
  out
}

reuse_bundle_cache <- new.env(parent = emptyenv())

get_reuse_bundle <- function(model_id, export_dir, registry_row = NULL) {
  cache_key <- paste0("bundle__", model_id)
  if (exists(cache_key, envir = reuse_bundle_cache, inherits = FALSE)) {
    return(get(cache_key, envir = reuse_bundle_cache, inherits = FALSE))
  }
  cands <- get_model_rds_candidates(model_id, export_dir, registry_row)
  out <- empty_reuse_bundle()
  if (length(cands) > 0L) {
    for (cand in cands) {
      if (!file.exists(cand)) next
      obj <- tryCatch(readRDS(cand), error = function(e) NULL)
      if (is.null(obj)) next
      out <- extract_stage8_bundle_tables(obj)
      out$rds_path <- cand
      break
    }
  }
  assign(cache_key, out, envir = reuse_bundle_cache)
  out
}

normalize_prior_branch_value <- function(x) {
  x <- as.character(x %||% NA_character_)
  if (is.na(x) || !nzchar(x)) {
    return(NA_character_)
  }
  if (identical(x, "neutral_weakly_informative")) {
    return("neutral_no_external_info")
  }
  x
}

registry_field_value <- function(reg, field, default = NA) {
  if (is.null(reg) || !inherits(reg, "data.frame") || nrow(reg) == 0L || !(field %in% names(reg))) {
    return(default)
  }
  value <- reg[[field]][[1L]]
  if (is.null(value) || length(value) == 0L || (length(value) == 1L && is.na(value))) {
    return(default)
  }
  value
}

load_existing_stage8_exports <- function(export_dir) {
  list(
    model_registry = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_model_registry.csv")),
    coefficient_summary = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_coefficient_summary.csv")),
    posterior_subject_profile = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_posterior_subject_profile.csv.gz")),
    posterior_subject_yearly = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_posterior_subject_yearly.csv.gz")),
    posterior_cohort_yearly = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_posterior_cohort_yearly.csv")),
    posterior_classification = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_posterior_classification.csv")),
    posterior_delta_vs_nocure = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_posterior_delta_vs_nocure.csv")),
    diagnostics_parameter_level = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_diagnostics_parameter_level.csv")),
    ppc_summary = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_ppc_summary.csv")),
    prior_predictive_summary = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_prior_predictive_summary.csv")),
    prior_branch_delta = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_prior_branch_delta.csv")),
    incidence_anchor_update = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_incidence_anchor_update.csv")),
    hazard_plausibility = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_hazard_plausibility.csv")),
    uncured_decomposition = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_uncured_decomposition.csv")),
    stage8A_vs_stage8B_delta = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_stage8A_vs_stage8B_delta.csv")),
    diagnostic_pdf_exists = pdf_file_is_usable(file.path(export_dir, "bayes_stage8_diagnostic_plots.pdf"))
  )
}

is_model_reusable <- function(model_row, dataset_df, existing_exports, export_dir, require_existing_rds = FALSE) {
  reg <- subset_model_registry_row(existing_exports$model_registry, model_row$model_id[[1L]])
  if (is.null(reg) || nrow(reg) != 1L) {
    return(FALSE)
  }

  reg_dataset <- as.character(registry_field_value(reg, "dataset", registry_field_value(reg, "dataset_key", NA_character_)))
  reg_family <- as.character(registry_field_value(reg, "family_code", NA_character_))
  reg_formula <- as.character(registry_field_value(reg, "formula_anchor", NA_character_))
  reg_branch <- as.character(registry_field_value(reg, "branch", "Stage8A"))
  reg_risk_scale <- as.character(registry_field_value(reg, "risk_scale", "transition_only_main"))
  reg_prior_branch <- normalize_prior_branch_value(registry_field_value(reg, "prior_branch", "anchor_informed"))
  reg_site_prior <- as.character(registry_field_value(reg, "site_prior_family", ifelse(model_row$has_site_term[[1L]], "normal_0_1_main", "not_applicable")))

  if (!identical(reg_dataset, as.character(model_row$dataset[[1L]]))) return(FALSE)
  if (!identical(reg_family, as.character(model_row$family_code[[1L]]))) return(FALSE)
  if (!identical(reg_formula, as.character(model_row$formula_anchor[[1L]]))) return(FALSE)
  if (!identical(reg_branch, "Stage8A")) return(FALSE)
  if (!identical(reg_risk_scale, "transition_only_main")) return(FALSE)
  if (!identical(reg_prior_branch, as.character(model_row$prior_branch[[1L]]))) return(FALSE)
  if (!identical(reg_site_prior, as.character(model_row$site_prior_family[[1L]]))) return(FALSE)
  if (!identical(as.character(registry_field_value(reg, "fit_status", NA_character_)), "ok")) return(FALSE)

  n_event_now <- sum(dataset_df$event_main)
  n_censor_now <- sum(dataset_df$censor_main)
  n_remission_now <- sum(dataset_df$remission_flag)
  if (!isTRUE(nrow(dataset_df) == as.integer(safe_numeric(registry_field_value(reg, "n", NA_integer_))))) return(FALSE)
  if (!isTRUE(n_event_now == as.integer(safe_numeric(registry_field_value(reg, "n_event", NA_integer_))))) return(FALSE)
  if (!isTRUE(n_censor_now == as.integer(safe_numeric(registry_field_value(reg, "n_censor_main", NA_integer_))))) return(FALSE)
  if (!isTRUE(n_remission_now == as.integer(safe_numeric(registry_field_value(reg, "n_remission", NA_integer_))))) return(FALSE)

  rds_candidates <- get_model_rds_candidates(model_row$model_id[[1L]], export_dir, reg)
  has_rds <- any(file.exists(rds_candidates))
  if (isTRUE(require_existing_rds) && !has_rds) {
    return(FALSE)
  }

  reuse_bundle <- if (has_rds) get_reuse_bundle(model_row$model_id[[1L]], export_dir, reg) else empty_reuse_bundle()
  required_nonadmissible <- c("coefficient_summary", "diagnostics_parameter_level", "ppc_summary")
  for (nm in required_nonadmissible) {
    has_existing <- nrow(subset_model_table(existing_exports[[nm]], model_row$model_id[[1L]])) > 0L
    has_bundle <- nrow(subset_model_table(reuse_bundle[[nm]], model_row$model_id[[1L]])) > 0L
    if (!has_existing && !has_bundle) {
      return(FALSE)
    }
  }

  reg_admissible <- isTRUE(as.logical(registry_field_value(reg, "admissible_flag", FALSE)))
  if (reg_admissible) {
    required_admissible <- c("posterior_subject_profile", "posterior_subject_yearly", "posterior_cohort_yearly", "posterior_classification")
    for (nm in required_admissible) {
      has_existing <- nrow(subset_model_table(existing_exports[[nm]], model_row$model_id[[1L]])) > 0L
      has_bundle <- nrow(subset_model_table(reuse_bundle[[nm]], model_row$model_id[[1L]])) > 0L
      if (!has_existing && !has_bundle) {
        return(FALSE)
      }
    }
  }
  TRUE
}

## 🟠 Define: Stage-8A structural grid, prior branches, and design bundles ===============================
build_model_grid <- function(include_supplementary = TRUE) {
  base_rows <- tibble::tribble(
    ~dataset, ~structural_model_id, ~latency_branch, ~formula_anchor, ~incidence_site_indicator, ~latency_site_indicator, ~latency_interaction_indicator, ~is_supplementary_branch,
    "PNU", "PNU-L0", "L0", "base", FALSE, FALSE, FALSE, FALSE,
    "PNU", "PNU-L1", "L1", "interaction", FALSE, FALSE, TRUE, FALSE,
    "SNU", "SNU-L0", "L0", "base", FALSE, FALSE, FALSE, FALSE,
    "SNU", "SNU-L1", "L1", "interaction", FALSE, FALSE, TRUE, FALSE,
    "merged", "MERGED-L0S0", "L0S0", "base", FALSE, FALSE, FALSE, FALSE,
    "merged", "MERGED-L1S0", "L1S0", "interaction", FALSE, FALSE, TRUE, FALSE,
    "merged", "MERGED-L0S1", "L0S1", "site_added", FALSE, TRUE, FALSE, FALSE,
    "merged", "MERGED-L1S1", "L1S1", "site_interaction", FALSE, TRUE, TRUE, FALSE
  )

  if (isTRUE(include_supplementary)) {
    base_rows <- bind_rows(
      base_rows,
      tibble::tribble(
        ~dataset, ~structural_model_id, ~latency_branch, ~formula_anchor, ~incidence_site_indicator, ~latency_site_indicator, ~latency_interaction_indicator, ~is_supplementary_branch,
        "merged", "MERGED-Isite-L0S0", "L0S0", "base", TRUE, FALSE, FALSE, TRUE,
        "merged", "MERGED-Isite-L1S0", "L1S0", "interaction", TRUE, FALSE, TRUE, TRUE,
        "merged", "MERGED-Isite-L0S1", "L0S1", "site_added", TRUE, TRUE, FALSE, TRUE,
        "merged", "MERGED-Isite-L1S1", "L1S1", "site_interaction", TRUE, TRUE, TRUE, TRUE
      )
    )
  }

  family_rows <- tibble::tribble(
    ~family_code, ~latency_family, ~family_id,
    "E", "exponential", 1L,
    "W", "weibull", 2L,
    "LN", "lognormal", 3L,
    "LL", "loglogistic", 4L
  )

  prior_branch_rows <- tibble::tribble(
    ~prior_branch,
    "anchor_informed",
    "neutral_no_external_info"
  )

  out <- tidyr::crossing(base_rows, family_rows, prior_branch_rows) %>%
    rowwise() %>%
    mutate(
      has_site_term = isTRUE(incidence_site_indicator) || isTRUE(latency_site_indicator),
      site_prior_family = list(if (has_site_term) c("normal_0_1_main", "student_t3_0_1_sensitivity") else "not_applicable")
    ) %>%
    tidyr::unnest(site_prior_family) %>%
    ungroup() %>%
    mutate(
      branch = "Stage8A",
      risk_scale = "transition_only_main",
      dataset_key = dataset,
      site_placement_label = case_when(
        incidence_site_indicator & latency_site_indicator ~ "site_in_both",
        incidence_site_indicator & !latency_site_indicator ~ "site_in_incidence_only",
        !incidence_site_indicator & latency_site_indicator ~ "site_in_latency_only",
        TRUE ~ "site_in_neither"
      ),
      legacy_model_id = paste(structural_model_id, family_code, sep = "-"),
      model_id = dplyr::case_when(
        prior_branch == "anchor_informed" & site_prior_family %in% c("normal_0_1_main", "not_applicable") ~ legacy_model_id,
        prior_branch == "neutral_no_external_info" & site_prior_family == "not_applicable" ~ paste(legacy_model_id, "PB-neutral", sep = "__"),
        prior_branch == "neutral_no_external_info" & site_prior_family == "normal_0_1_main" ~ paste(legacy_model_id, "PB-neutral", sep = "__"),
        prior_branch == "anchor_informed" & site_prior_family == "student_t3_0_1_sensitivity" ~ paste(legacy_model_id, "SP-t3", sep = "__"),
        TRUE ~ paste(legacy_model_id, "PB-neutral", "SP-t3", sep = "__")
      )
    ) %>%
    arrange(factor(dataset, levels = c("PNU", "SNU", "merged")), structural_model_id, family_code, prior_branch, site_prior_family)

  out
}

build_prior_specs <- function() {
  list(
    anchor_informed = list(
      prior_branch = "anchor_informed",
      alpha_gp = -9.581369553169,
      mu_beta_inc_base = c(0.419871845822, 0.907608052926, 0.586202561451, 0.466865123863, 0.037997248763),
      sd_beta_inc_base = c(0.132789397422, 0.173731076538, 0.191221553945, 0.270393197518, 0.302838606651),
      sd_beta_inc_site = 1.0,
      sd_delta = 3.5,
      sd_gamma0 = 2.5,
      sd_gamma = 1.0,
      site_prior_df = 3.0,
      site_prior_scale = 1.0,
      weibull_shape_sd = 0.35,
      lognormal_shape_sd = 0.50,
      loglogistic_shape_sd = 0.50
    ),
    neutral_no_external_info = list(
      prior_branch = "neutral_no_external_info",
      alpha_gp = 0.0,
      mu_beta_inc_base = rep(0, 5L),
      sd_beta_inc_base = rep(2.0, 5L),
      sd_beta_inc_site = 1.0,
      sd_delta = 3.5,
      sd_gamma0 = 2.5,
      sd_gamma = 1.0,
      site_prior_df = 3.0,
      site_prior_scale = 1.0,
      weibull_shape_sd = 0.35,
      lognormal_shape_sd = 0.50,
      loglogistic_shape_sd = 0.50
    )
  )
}

make_design_bundle <- function(df, model_row, prior_spec, snu_label) {
  z_i <- as.integer(df$sex_num)
  x20_i <- as.integer(df$age_exact_entry >= 20 & df$age_exact_entry < 30)
  x30_i <- as.integer(df$age_exact_entry >= 30)
  s_i <- as.integer(df$site == snu_label)

  X_inc <- cbind(
    sex_num = z_i,
    age20_29 = x20_i,
    age30plus = x30_i,
    sex_x_age20_29 = z_i * x20_i,
    sex_x_age30plus = z_i * x30_i
  )
  mu_beta_inc <- prior_spec$mu_beta_inc_base
  sd_beta_inc <- prior_spec$sd_beta_inc_base
  incidence_site_position <- 0L

  if (isTRUE(model_row$incidence_site_indicator[[1L]])) {
    X_inc <- cbind(X_inc, site_SNU = s_i)
    mu_beta_inc <- c(mu_beta_inc, 0)
    sd_beta_inc <- c(sd_beta_inc, prior_spec$sd_beta_inc_site)
    incidence_site_position <- ncol(X_inc)
  }

  a_i <- as.numeric(df$age_s)
  az_i <- a_i * z_i
  latency_site_position <- 0L
  X_lat <- switch(
    model_row$latency_branch[[1L]],
    L0 = cbind(age_s = a_i, sex_num = z_i),
    L1 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i),
    L0S0 = cbind(age_s = a_i, sex_num = z_i),
    L1S0 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i),
    L0S1 = {
      out <- cbind(age_s = a_i, sex_num = z_i, site_SNU = s_i)
      latency_site_position <- 3L
      out
    },
    L1S1 = {
      out <- cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i, site_SNU = s_i)
      latency_site_position <- 4L
      out
    },
    stop(sprintf("Unknown latency branch `%s`.", model_row$latency_branch[[1L]]), call. = FALSE)
  )

  shape_prior_sd <- switch(
    model_row$family_code[[1L]],
    E = 1.0,
    W = prior_spec$weibull_shape_sd,
    LN = prior_spec$lognormal_shape_sd,
    LL = prior_spec$loglogistic_shape_sd,
    1.0
  )

  list(
    X_inc = unclass(as.matrix(X_inc)),
    X_lat = unclass(as.matrix(X_lat)),
    mu_beta_inc = as.numeric(mu_beta_inc),
    sd_beta_inc = as.numeric(sd_beta_inc),
    sd_gamma_lat = rep(prior_spec$sd_gamma, ncol(X_lat)),
    alpha_gp = prior_spec$alpha_gp,
    time = as.numeric(df$time_year),
    event = as.integer(df$event_main),
    id_df = df %>% select(unique_person_id, id, site, sex_num, age_exact_entry, age_s),
    incidence_site_position = as.integer(incidence_site_position),
    latency_site_position = as.integer(latency_site_position),
    use_student_t_site_prior = as.integer(model_row$site_prior_family[[1L]] == "student_t3_0_1_sensitivity"),
    site_prior_df = as.numeric(prior_spec$site_prior_df),
    site_prior_scale = as.numeric(prior_spec$site_prior_scale),
    shape_prior_sd = as.numeric(shape_prior_sd),
    x20_i = x20_i,
    x30_i = x30_i,
    s_i = s_i,
    z_i = z_i
  )
}
# 🔴 Define: survival math, IPCW metrics, and Bayesian diagnostics ===============================
## 🟠 Define: posterior summaries, prior predictive simulation, and threshold metrics ===============================
km_eval <- function(survfit_obj, times) {
  base_times <- survfit_obj$time
  base_surv <- survfit_obj$surv
  vapply(times, FUN.VALUE = numeric(1), FUN = function(tt) {
    idx <- max(c(0L, which(base_times <= tt)))
    if (idx == 0L) 1 else base_surv[[idx]]
  })
}

build_ipcw_reference <- function(df, horizons) {
  event_fit <- survival::survfit(survival::Surv(time_year, event_main) ~ 1, data = df)
  censor_fit <- survival::survfit(survival::Surv(time_year, censor_main) ~ 1, data = df)

  horizon_rows <- lapply(horizons, function(h) {
    G_t <- pmax(km_eval(censor_fit, h), 1e-8)
    G_tminus <- pmax(km_eval(censor_fit, pmax(df$time_year - 1e-10, 0)), 1e-8)
    w_case <- ifelse(df$event_main == 1L & df$time_year <= h, 1 / G_tminus, 0)
    w_control <- ifelse(df$time_year > h, 1 / G_t, 0)
    prevalence <- 1 - km_eval(event_fit, h)

    tibble(
      horizon_year = h,
      observed_km_risk = prevalence,
      denom_case = sum(w_case),
      denom_control = sum(w_control),
      G_t = G_t,
      n_subject = nrow(df),
      w_case = list(w_case),
      w_control = list(w_control)
    )
  })

  bind_rows(horizon_rows)
}

summary_scalar <- function(x) {
  tibble(
    mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    q025 = stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    q50 = stats::quantile(x, 0.500, na.rm = TRUE, names = FALSE),
    q975 = stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE)
  )
}

summarize_cols_matrix <- function(mat) {
  q <- matrixStats::colQuantiles(mat, probs = c(0.025, 0.500, 0.975), na.rm = TRUE)
  tibble(
    mean = matrixStats::colMeans2(mat, na.rm = TRUE),
    q025 = q[, 1],
    q50 = q[, 2],
    q975 = q[, 3]
  )
}

compute_linear_terms <- function(draws, X_inc, X_lat, alpha_gp) {
  eta_inc <- draws$beta_inc %*% t(X_inc)
  eta_inc <- sweep(eta_inc, 1, alpha_gp + draws$delta0, FUN = "+")
  pi_mat <- plogis(eta_inc)

  mu_lat <- draws$gamma_lat %*% t(X_lat)
  mu_lat <- sweep(mu_lat, 1, draws$gamma0, FUN = "+")
  median_mat <- exp(mu_lat)

  list(
    pi_mat = pi_mat,
    cure_prob_mat = 1 - pi_mat,
    mu_lat_mat = mu_lat,
    median_mat = median_mat
  )
}

family_survival_hazard <- function(horizon_year, family_code, mu_lat_mat, median_mat, draws) {
  if (family_code == "E") {
    lambda <- median_mat / log(2)
    Su <- exp(-horizon_year / lambda)
    haz <- 1 / lambda
  } else if (family_code == "W") {
    k <- exp(draws$rho_W)
    lambda <- median_mat / ((log(2))^(1 / k))
    ratio <- horizon_year / lambda
    Su <- exp(-(ratio^k))
    haz <- (k / lambda) * (ratio^(k - 1))
  } else if (family_code == "LN") {
    sigma <- exp(draws$log_sigma_LN)
    z <- (log(horizon_year) - mu_lat_mat) / sigma
    Su <- 1 - pnorm(z)
    log_pdf <- dnorm(z, log = TRUE) - log(horizon_year) - log(sigma)
    haz <- exp(log_pdf) / pmax(Su, 1e-12)
  } else if (family_code == "LL") {
    k <- exp(-draws$psi_LL)
    lambda <- median_mat
    ratio <- horizon_year / lambda
    Su <- 1 / (1 + ratio^k)
    haz <- (k / lambda) * (ratio^(k - 1)) / (1 + ratio^k)
  } else {
    stop(sprintf("Unknown family code `%s`.", family_code), call. = FALSE)
  }

  Su <- pmin(pmax(Su, 1e-12), 1 - 1e-12)
  haz <- pmax(haz, 1e-12)
  list(Su = Su, haz = haz)
}

extract_draws_compact <- function(fit, K_inc, K_lat) {
  ext <- rstan::extract(
    fit,
    pars = c("delta0", "beta_inc", "gamma0", "gamma_lat", "rho_W", "log_sigma_LN", "psi_LL", "log_lik"),
    permuted = TRUE,
    inc_warmup = FALSE
  )

  beta_inc <- ext$beta_inc
  gamma_lat <- ext$gamma_lat
  log_lik <- ext$log_lik

  if (is.null(dim(beta_inc))) beta_inc <- matrix(beta_inc, ncol = K_inc)
  if (is.null(dim(gamma_lat))) gamma_lat <- matrix(gamma_lat, ncol = K_lat)
  if (is.null(dim(log_lik))) log_lik <- matrix(log_lik, nrow = length(ext$delta0))

  list(
    delta0 = as.numeric(ext$delta0),
    beta_inc = beta_inc,
    gamma0 = as.numeric(ext$gamma0),
    gamma_lat = gamma_lat,
    rho_W = as.numeric(ext$rho_W),
    log_sigma_LN = as.numeric(ext$log_sigma_LN),
    psi_LL = as.numeric(ext$psi_LL),
    log_lik = log_lik
  )
}

compute_degeneracy <- function(pi_mat, median_mat, supported_risk_list) {
  near_zero_pi <- rowMeans(pi_mat < tiny_susceptible_prob) > degenerate_subject_fraction_threshold
  near_one_pi <- rowMeans(pi_mat > huge_susceptible_prob) > degenerate_subject_fraction_threshold
  near_zero_median <- rowMeans(median_mat < tiny_median_years) > degenerate_subject_fraction_threshold

  if (length(supported_risk_list) > 0L) {
    risk_mean_supported <- Reduce("+", supported_risk_list) / length(supported_risk_list)
    huge_median_and_flat <- (rowMeans(median_mat > huge_median_years) > degenerate_subject_fraction_threshold) &
      (rowMeans(risk_mean_supported < 0.01) > degenerate_subject_fraction_threshold)
  } else {
    huge_median_and_flat <- rep(FALSE, nrow(pi_mat))
  }

  any_problem <- near_zero_pi | near_one_pi | near_zero_median | huge_median_and_flat

  tibble(
    near_zero_pi_rate = mean(near_zero_pi),
    near_one_pi_rate = mean(near_one_pi),
    near_zero_median_rate = mean(near_zero_median),
    huge_median_flat_risk_rate = mean(huge_median_and_flat),
    degenerate_flag = mean(any_problem) > degenerate_draw_fraction_threshold
  )
}

simulate_prior_predictive <- function(df, design_bundle, model_row, prior_spec, n_draws, horizons_eval) {
  K_inc <- ncol(design_bundle$X_inc)
  K_lat <- ncol(design_bundle$X_lat)

  sim_draws <- list(
    delta0 = rnorm(n_draws, 0, prior_spec$sd_delta),
    beta_inc = matrix(
      rnorm(
        n_draws * K_inc,
        mean = rep(design_bundle$mu_beta_inc, each = n_draws),
        sd = rep(design_bundle$sd_beta_inc, each = n_draws)
      ),
      nrow = n_draws,
      byrow = FALSE
    ),
    gamma0 = rnorm(n_draws, 0, prior_spec$sd_gamma0),
    gamma_lat = matrix(
      rnorm(n_draws * K_lat, mean = 0, sd = rep(design_bundle$sd_gamma_lat, each = n_draws)),
      nrow = n_draws,
      byrow = FALSE
    ),
    rho_W = rnorm(n_draws, 0, prior_spec$weibull_shape_sd),
    log_sigma_LN = rnorm(n_draws, 0, prior_spec$lognormal_shape_sd),
    psi_LL = rnorm(n_draws, 0, prior_spec$loglogistic_shape_sd)
  )

  lin <- compute_linear_terms(sim_draws, design_bundle$X_inc, design_bundle$X_lat, design_bundle$alpha_gp)
  supported_horizons <- if (model_row$dataset[[1L]] == "PNU") c(1, 2) else c(1, 2, 5)
  supported_risk_list <- lapply(supported_horizons, function(h) {
    fh <- family_survival_hazard(h, model_row$family_code[[1L]], lin$mu_lat_mat, lin$median_mat, sim_draws)
    1 - ((1 - lin$pi_mat) + lin$pi_mat * fh$Su)
  })
  degeneracy <- compute_degeneracy(lin$pi_mat, lin$median_mat, supported_risk_list)

  prior_rows <- list(
    tibble(dataset = model_row$dataset[[1L]], model_id = model_row$model_id[[1L]], branch = "Stage8A", risk_scale = "transition_only_main", prior_branch = model_row$prior_branch[[1L]], metric = "cohort_mean_susceptible_probability", horizon_year = NA_integer_, summary_scalar(rowMeans(lin$pi_mat)), prior_degenerate_flag = degeneracy$degenerate_flag),
    tibble(dataset = model_row$dataset[[1L]], model_id = model_row$model_id[[1L]], branch = "Stage8A", risk_scale = "transition_only_main", prior_branch = model_row$prior_branch[[1L]], metric = "cohort_mean_cure_probability", horizon_year = NA_integer_, summary_scalar(rowMeans(1 - lin$pi_mat)), prior_degenerate_flag = degeneracy$degenerate_flag),
    tibble(dataset = model_row$dataset[[1L]], model_id = model_row$model_id[[1L]], branch = "Stage8A", risk_scale = "transition_only_main", prior_branch = model_row$prior_branch[[1L]], metric = "cohort_median_susceptible_time", horizon_year = NA_integer_, summary_scalar(apply(lin$median_mat, 1, stats::median)), prior_degenerate_flag = degeneracy$degenerate_flag)
  )

  risk_rows <- lapply(horizons_eval, function(h) {
    fh <- family_survival_hazard(h, model_row$family_code[[1L]], lin$mu_lat_mat, lin$median_mat, sim_draws)
    risk_mat <- 1 - ((1 - lin$pi_mat) + lin$pi_mat * fh$Su)
    tibble(dataset = model_row$dataset[[1L]], model_id = model_row$model_id[[1L]], branch = "Stage8A", risk_scale = "transition_only_main", prior_branch = model_row$prior_branch[[1L]], metric = "cohort_mean_risk", horizon_year = h, summary_scalar(rowMeans(risk_mat)), prior_degenerate_flag = degeneracy$degenerate_flag)
  })

  list(summary = bind_rows(prior_rows, risk_rows), degeneracy = degeneracy)
}

weighted_auc_horizon <- function(score, horizon_row) {
  w_case <- unlist(horizon_row$w_case)
  w_control <- unlist(horizon_row$w_control)
  case_idx <- which(w_case > 0)
  control_idx <- which(w_control > 0)
  if (length(case_idx) == 0L || length(control_idx) == 0L) {
    return(NA_real_)
  }

  controls <- tibble(score = score[control_idx], w = w_control[control_idx]) %>%
    group_by(score) %>%
    summarise(w = sum(w), .groups = "drop") %>%
    arrange(score)
  controls$cum_less <- c(0, head(cumsum(controls$w), -1L))
  match_idx <- match(score[case_idx], controls$score)
  less_w <- controls$cum_less[match_idx]
  tie_w <- controls$w[match_idx]
  num <- sum(w_case[case_idx] * (less_w + 0.5 * tie_w))
  den <- sum(w_case[case_idx]) * sum(controls$w)
  if (!is.finite(den) || den <= 0) {
    return(NA_real_)
  }
  num / den
}

ipcw_brier_horizon <- function(score, horizon_row) {
  w_case <- unlist(horizon_row$w_case)
  w_control <- unlist(horizon_row$w_control)
  n_subject <- length(score)
  if (n_subject == 0L) return(NA_real_)
  sum(w_case * (1 - score)^2 + w_control * score^2) / n_subject
}

compute_classification_summary <- function(risk_draws, horizon_row, thresholds) {
  prevalence <- as.numeric(horizon_row$observed_km_risk)
  w_case <- unlist(horizon_row$w_case)
  w_control <- unlist(horizon_row$w_control)
  denom_case <- as.numeric(horizon_row$denom_case)
  denom_control <- as.numeric(horizon_row$denom_control)
  n_subject <- ncol(risk_draws)

  out_abs <- vector("list", length(thresholds))

  for (j in seq_along(thresholds)) {
    thr <- thresholds[[j]]
    H_mat <- (risk_draws >= thr) * 1

    if (denom_case > 0) {
      tpr_draw <- as.vector(H_mat %*% w_case) / denom_case
    } else {
      tpr_draw <- rep(NA_real_, nrow(risk_draws))
    }

    if (denom_control > 0) {
      fpr_draw <- as.vector(H_mat %*% w_control) / denom_control
    } else {
      fpr_draw <- rep(NA_real_, nrow(risk_draws))
    }

    pos_rate_draw <- prevalence * tpr_draw + (1 - prevalence) * fpr_draw
    pos_count_draw <- pos_rate_draw * n_subject
    ppv_draw <- ifelse(pos_rate_draw > 0, prevalence * tpr_draw / pos_rate_draw, NA_real_)
    fdp_draw <- 1 - ppv_draw
    fp_burden_draw <- (1 - prevalence) * fpr_draw
    fp_count_draw <- fp_burden_draw * n_subject
    fp100_draw <- 100 * fp_burden_draw
    nb_draw <- prevalence * tpr_draw - (1 - prevalence) * fpr_draw * (thr / (1 - thr))

    out_abs[[j]] <- tibble(
      threshold = thr,
      positive_classification_rate_mean = mean(pos_rate_draw, na.rm = TRUE),
      positive_classification_rate_q025 = stats::quantile(pos_rate_draw, 0.025, na.rm = TRUE, names = FALSE),
      positive_classification_rate_q50 = stats::quantile(pos_rate_draw, 0.500, na.rm = TRUE, names = FALSE),
      positive_classification_rate_q975 = stats::quantile(pos_rate_draw, 0.975, na.rm = TRUE, names = FALSE),
      positive_classification_count_mean = mean(pos_count_draw, na.rm = TRUE),
      positive_classification_count_q025 = stats::quantile(pos_count_draw, 0.025, na.rm = TRUE, names = FALSE),
      positive_classification_count_q50 = stats::quantile(pos_count_draw, 0.500, na.rm = TRUE, names = FALSE),
      positive_classification_count_q975 = stats::quantile(pos_count_draw, 0.975, na.rm = TRUE, names = FALSE),
      FPR_mean = mean(fpr_draw, na.rm = TRUE),
      FPR_q025 = stats::quantile(fpr_draw, 0.025, na.rm = TRUE, names = FALSE),
      FPR_q50 = stats::quantile(fpr_draw, 0.500, na.rm = TRUE, names = FALSE),
      FPR_q975 = stats::quantile(fpr_draw, 0.975, na.rm = TRUE, names = FALSE),
      false_positive_burden_mean = mean(fp_burden_draw, na.rm = TRUE),
      false_positive_burden_q025 = stats::quantile(fp_burden_draw, 0.025, na.rm = TRUE, names = FALSE),
      false_positive_burden_q50 = stats::quantile(fp_burden_draw, 0.500, na.rm = TRUE, names = FALSE),
      false_positive_burden_q975 = stats::quantile(fp_burden_draw, 0.975, na.rm = TRUE, names = FALSE),
      false_positive_count_mean = mean(fp_count_draw, na.rm = TRUE),
      false_positive_count_q025 = stats::quantile(fp_count_draw, 0.025, na.rm = TRUE, names = FALSE),
      false_positive_count_q50 = stats::quantile(fp_count_draw, 0.500, na.rm = TRUE, names = FALSE),
      false_positive_count_q975 = stats::quantile(fp_count_draw, 0.975, na.rm = TRUE, names = FALSE),
      FP100_mean = mean(fp100_draw, na.rm = TRUE),
      FP100_q025 = stats::quantile(fp100_draw, 0.025, na.rm = TRUE, names = FALSE),
      FP100_q50 = stats::quantile(fp100_draw, 0.500, na.rm = TRUE, names = FALSE),
      FP100_q975 = stats::quantile(fp100_draw, 0.975, na.rm = TRUE, names = FALSE),
      PPV_mean = mean(ppv_draw, na.rm = TRUE),
      PPV_q025 = stats::quantile(ppv_draw, 0.025, na.rm = TRUE, names = FALSE),
      PPV_q50 = stats::quantile(ppv_draw, 0.500, na.rm = TRUE, names = FALSE),
      PPV_q975 = stats::quantile(ppv_draw, 0.975, na.rm = TRUE, names = FALSE),
      FDP_mean = mean(fdp_draw, na.rm = TRUE),
      FDP_q025 = stats::quantile(fdp_draw, 0.025, na.rm = TRUE, names = FALSE),
      FDP_q50 = stats::quantile(fdp_draw, 0.500, na.rm = TRUE, names = FALSE),
      FDP_q975 = stats::quantile(fdp_draw, 0.975, na.rm = TRUE, names = FALSE),
      TPR_mean = mean(tpr_draw, na.rm = TRUE),
      TPR_q025 = stats::quantile(tpr_draw, 0.025, na.rm = TRUE, names = FALSE),
      TPR_q50 = stats::quantile(tpr_draw, 0.500, na.rm = TRUE, names = FALSE),
      TPR_q975 = stats::quantile(tpr_draw, 0.975, na.rm = TRUE, names = FALSE),
      NB_mean = mean(nb_draw, na.rm = TRUE),
      NB_q025 = stats::quantile(nb_draw, 0.025, na.rm = TRUE, names = FALSE),
      NB_q50 = stats::quantile(nb_draw, 0.500, na.rm = TRUE, names = FALSE),
      NB_q975 = stats::quantile(nb_draw, 0.975, na.rm = TRUE, names = FALSE)
    )
  }
  bind_rows(out_abs)
}

compute_information_criteria <- function(log_lik) {
  out <- list(
    waic = NA_real_,
    p_waic = NA_real_,
    p_waic_high_n = NA_integer_,
    p_waic_high_pct = NA_real_,
    looic = NA_real_,
    p_loo = NA_real_,
    pareto_k_max = NA_real_,
    pareto_k_bad_n = NA_integer_,
    pareto_k_bad_pct = NA_real_,
    pareto_k_very_bad_n = NA_integer_,
    waic_warning_flag = FALSE,
    loo_warning_flag = FALSE,
    info_criteria_warning_detail = NA_character_
  )
  if (is.null(log_lik) || !requireNamespace("loo", quietly = TRUE)) {
    return(out)
  }

  warning_texts <- character()
  collect_warning <- function(w) {
    warning_texts <<- c(warning_texts, conditionMessage(w))
    tryInvokeRestart("muffleWarning")
  }

  waic_obj <- withCallingHandlers(tryCatch(loo::waic(log_lik), error = function(e) e), warning = collect_warning)
  if (!inherits(waic_obj, "error")) {
    if ("waic" %in% rownames(waic_obj$estimates)) out$waic <- as.numeric(waic_obj$estimates["waic", "Estimate"])
    if ("p_waic" %in% rownames(waic_obj$estimates)) out$p_waic <- as.numeric(waic_obj$estimates["p_waic", "Estimate"])
    if (!is.null(waic_obj$pointwise) && "p_waic" %in% colnames(waic_obj$pointwise)) {
      p_waic_pointwise <- as.numeric(waic_obj$pointwise[, "p_waic"])
      out$p_waic_high_n <- as.integer(sum(p_waic_pointwise > 0.4, na.rm = TRUE))
      out$p_waic_high_pct <- if (length(p_waic_pointwise) > 0L) 100 * out$p_waic_high_n / length(p_waic_pointwise) else NA_real_
    }
  } else {
    warning_texts <- c(warning_texts, paste0("waic_error: ", conditionMessage(waic_obj)))
  }

  loo_obj <- withCallingHandlers(tryCatch(loo::loo(log_lik), error = function(e) e), warning = collect_warning)
  if (!inherits(loo_obj, "error")) {
    if ("looic" %in% rownames(loo_obj$estimates)) out$looic <- as.numeric(loo_obj$estimates["looic", "Estimate"])
    if ("p_loo" %in% rownames(loo_obj$estimates)) out$p_loo <- as.numeric(loo_obj$estimates["p_loo", "Estimate"])
    if (!is.null(loo_obj$diagnostics$pareto_k)) {
      pareto_k <- as.numeric(loo_obj$diagnostics$pareto_k)
      out$pareto_k_max <- if (length(pareto_k) > 0L && any(is.finite(pareto_k))) max(pareto_k, na.rm = TRUE) else NA_real_
      out$pareto_k_bad_n <- as.integer(sum(pareto_k > 0.7, na.rm = TRUE))
      out$pareto_k_very_bad_n <- as.integer(sum(pareto_k > 1.0, na.rm = TRUE))
      out$pareto_k_bad_pct <- if (length(pareto_k) > 0L) 100 * out$pareto_k_bad_n / length(pareto_k) else NA_real_
    }
  } else {
    warning_texts <- c(warning_texts, paste0("loo_error: ", conditionMessage(loo_obj)))
  }

  out$waic_warning_flag <- any(grepl("p_waic", warning_texts, fixed = TRUE))
  out$loo_warning_flag <- any(grepl("Pareto k", warning_texts, fixed = TRUE))
  if (length(warning_texts) > 0L) out$info_criteria_warning_detail <- paste(unique(warning_texts), collapse = " | ")
  out
}

normalize_prior_predictive_summary <- function(df) {
  required_cols <- c("dataset", "model_id", "branch", "risk_scale", "prior_branch", "metric", "horizon_year", "mean", "q025", "q50", "q975", "prior_degenerate_flag")
  if (is.null(df) || !inherits(df, "data.frame") || nrow(df) == 0L) {
    out <- as_tibble(setNames(rep(list(logical()), length(required_cols)), required_cols))[0, ]
    out$dataset <- character(); out$model_id <- character(); out$branch <- character(); out$risk_scale <- character(); out$prior_branch <- character(); out$metric <- character(); out$horizon_year <- integer(); out$mean <- double(); out$q025 <- double(); out$q50 <- double(); out$q975 <- double(); out$prior_degenerate_flag <- logical()
    return(out)
  }
  out <- tibble::as_tibble(df)
  for (nm in setdiff(required_cols, names(out))) {
    out[[nm]] <- if (nm == "horizon_year") NA_integer_ else if (nm %in% c("mean", "q025", "q50", "q975")) NA_real_ else if (nm == "prior_degenerate_flag") NA else NA_character_
  }
  out %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      model_id = as.character(model_id),
      branch = ifelse(is.na(branch) | !nzchar(branch), "Stage8A", as.character(branch)),
      risk_scale = ifelse(is.na(risk_scale) | !nzchar(risk_scale), "transition_only_main", as.character(risk_scale)),
      prior_branch = normalize_prior_branch_value(ifelse(is.na(prior_branch) | !nzchar(prior_branch), "anchor_informed", prior_branch)),
      metric = as.character(metric),
      horizon_year = as.integer(safe_numeric(horizon_year)),
      mean = safe_numeric(mean),
      q025 = safe_numeric(q025),
      q50 = safe_numeric(q50),
      q975 = safe_numeric(q975),
      prior_degenerate_flag = as.logical(prior_degenerate_flag)
    ) %>%
    distinct(dataset, model_id, metric, horizon_year, .keep_all = TRUE)
}

# 🔴 Define: Stan model compilation with site-prior sensitivity support ===============================
## 🟠 Define: compiled Stage-8A Stan program and trace helpers ===============================
compile_stage8_stan_model <- function() {
  stan_code <- r"(
functions {
  real loglogistic_lpdf_custom(real t, real lambda, real k) {
    real ratio;
    ratio = t / lambda;
    return log(k) - log(lambda) + (k - 1) * log(ratio) - 2 * log1p(pow(ratio, k));
  }
  real loglogistic_lccdf_custom(real t, real lambda, real k) {
    return -log1p(pow(t / lambda, k));
  }
}
data {
  int<lower=1> N;
  vector<lower=1e-8>[N] time;
  int<lower=0, upper=1> event[N];
  int<lower=1> K_inc;
  matrix[N, K_inc] X_inc;
  int<lower=1> K_lat;
  matrix[N, K_lat] X_lat;
  int<lower=1, upper=4> family_id;
  real alpha_gp;
  vector[K_inc] mu_beta_inc;
  vector<lower=0>[K_inc] sd_beta_inc;
  real<lower=0> sd_delta;
  real<lower=0> sd_gamma0;
  vector<lower=0>[K_lat] sd_gamma_lat;
  real<lower=0> sd_shape;
  int<lower=0> incidence_site_position;
  int<lower=0> latency_site_position;
  int<lower=0, upper=1> use_student_t_site_prior;
  real<lower=0> site_prior_df;
  real<lower=0> site_prior_scale;
}
parameters {
  real delta0;
  vector[K_inc] beta_inc;
  real gamma0;
  vector[K_lat] gamma_lat;
  real rho_W;
  real log_sigma_LN;
  real psi_LL;
}
model {
  delta0 ~ normal(0, sd_delta);
  gamma0 ~ normal(0, sd_gamma0);
  rho_W ~ normal(0, sd_shape);
  log_sigma_LN ~ normal(0, sd_shape);
  psi_LL ~ normal(0, sd_shape);

  for (j in 1:K_inc) {
    if (use_student_t_site_prior == 1 && j == incidence_site_position) {
      beta_inc[j] ~ student_t(site_prior_df, 0, site_prior_scale);
    } else {
      beta_inc[j] ~ normal(mu_beta_inc[j], sd_beta_inc[j]);
    }
  }
  for (j in 1:K_lat) {
    if (use_student_t_site_prior == 1 && j == latency_site_position) {
      gamma_lat[j] ~ student_t(site_prior_df, 0, site_prior_scale);
    } else {
      gamma_lat[j] ~ normal(0, sd_gamma_lat[j]);
    }
  }

  for (i in 1:N) {
    real eta_inc;
    real pi_i;
    real mu_lat;
    real m_i;
    real logS;
    real logf;

    eta_inc = alpha_gp + delta0 + dot_product(row(X_inc, i), beta_inc);
    pi_i = inv_logit(eta_inc);
    mu_lat = gamma0 + dot_product(row(X_lat, i), gamma_lat);
    m_i = exp(mu_lat);

    if (family_id == 1) {
      real lambda;
      lambda = m_i / log(2);
      logS = -time[i] / lambda;
      logf = exponential_lpdf(time[i] | 1 / lambda);
    } else if (family_id == 2) {
      real kW;
      real lambda;
      kW = exp(rho_W);
      lambda = m_i / pow(log(2), 1 / kW);
      logS = -pow(time[i] / lambda, kW);
      logf = weibull_lpdf(time[i] | kW, lambda);
    } else if (family_id == 3) {
      real sigmaLN;
      sigmaLN = exp(log_sigma_LN);
      logS = normal_lccdf(log(time[i]) | mu_lat, sigmaLN);
      logf = lognormal_lpdf(time[i] | mu_lat, sigmaLN);
    } else {
      real kLL;
      real lambda;
      kLL = exp(-psi_LL);
      lambda = m_i;
      logS = loglogistic_lccdf_custom(time[i], lambda, kLL);
      logf = loglogistic_lpdf_custom(time[i], lambda, kLL);
    }

    if (event[i] == 1) {
      target += log(pi_i) + logf;
    } else {
      target += log_sum_exp(log1m(pi_i), log(pi_i) + logS);
    }
  }
}
generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    real eta_inc;
    real pi_i;
    real mu_lat;
    real m_i;
    real logS;
    real logf;

    eta_inc = alpha_gp + delta0 + dot_product(row(X_inc, i), beta_inc);
    pi_i = inv_logit(eta_inc);
    mu_lat = gamma0 + dot_product(row(X_lat, i), gamma_lat);
    m_i = exp(mu_lat);

    if (family_id == 1) {
      real lambda;
      lambda = m_i / log(2);
      logS = -time[i] / lambda;
      logf = exponential_lpdf(time[i] | 1 / lambda);
    } else if (family_id == 2) {
      real kW;
      real lambda;
      kW = exp(rho_W);
      lambda = m_i / pow(log(2), 1 / kW);
      logS = -pow(time[i] / lambda, kW);
      logf = weibull_lpdf(time[i] | kW, lambda);
    } else if (family_id == 3) {
      real sigmaLN;
      sigmaLN = exp(log_sigma_LN);
      logS = normal_lccdf(log(time[i]) | mu_lat, sigmaLN);
      logf = lognormal_lpdf(time[i] | mu_lat, sigmaLN);
    } else {
      real kLL;
      real lambda;
      kLL = exp(-psi_LL);
      lambda = m_i;
      logS = loglogistic_lccdf_custom(time[i], lambda, kLL);
      logf = loglogistic_lpdf_custom(time[i], lambda, kLL);
    }

    if (event[i] == 1) {
      log_lik[i] = log(pi_i) + logf;
    } else {
      log_lik[i] = log_sum_exp(log1m(pi_i), log(pi_i) + logS);
    }
  }
}
)"
  rstan::stan_model(model_code = stan_code, model_name = "stage8A_bayesian_cure_block")
}

select_trace_parameters <- function(family_code, K_inc, K_lat) {
  selected_pars <- c("delta0", "gamma0")
  if (K_inc >= 1L) selected_pars <- c(selected_pars, "beta_inc[1]")
  if (family_code == "E") {
    if (K_lat >= 1L) selected_pars <- c(selected_pars, "gamma_lat[1]")
  } else if (family_code == "W") {
    selected_pars <- c(selected_pars, "rho_W")
  } else if (family_code == "LN") {
    selected_pars <- c(selected_pars, "log_sigma_LN")
  } else if (family_code == "LL") {
    selected_pars <- c(selected_pars, "psi_LL")
  }
  unique(selected_pars)
}

make_trace_record <- function(fit, model_id, family_code, K_inc, K_lat) {
  selected_pars <- select_trace_parameters(family_code = family_code, K_inc = K_inc, K_lat = K_lat)
  list(model_id = model_id, selected_pars = selected_pars, arr = as.array(fit, pars = selected_pars))
}

plot_trace_record <- function(trace_record) {
  selected_pars <- trace_record$selected_pars
  arr <- trace_record$arr
  model_id <- trace_record$model_id
  par_old <- par(no.readonly = TRUE)
  on.exit(par(par_old), add = TRUE)
  par(mfrow = c(2, 2), mar = c(3, 3, 3, 1))
  for (j in seq_len(min(length(selected_pars), 4L))) {
    matplot(arr[, , j], type = "l", lty = 1, col = seq_len(dim(arr)[2]), main = paste(model_id, selected_pars[[j]]), xlab = "Iteration", ylab = "")
  }
  if (length(selected_pars) < 4L) {
    for (j in seq_len(4L - length(selected_pars))) plot.new()
  }
}
# 🔴 Define: post-fit derivations, deltas, and figure-source builders ===============================
## 🟠 Define: no-cure deltas, anchor updates, hazard plausibility, and plots ===============================
empty_delta_template <- function() {
  tibble(
    dataset = character(),
    branch = character(),
    risk_scale = character(),
    prior_branch = character(),
    site_prior_family = character(),
    model_id = character(),
    no_cure_model_id = character(),
    latency_family = character(),
    formula_anchor = character(),
    horizon_year = integer(),
    threshold = double(),
    metric = character(),
    no_cure_value = double(),
    delta_mean = double(),
    delta_q025 = double(),
    delta_q50 = double(),
    delta_q975 = double(),
    support_tier = character(),
    horizon_evidence_class = character(),
    claim_restriction_flag = character(),
    prior_tail_sensitive = logical(),
    classification_estimable_flag = logical(),
    denom_case = double(),
    denom_control = double(),
    admissible_flag = logical(),
    admissibility_reasons = character(),
    screening_flag = character(),
    screening_detail = character()
  )
}

rebuild_delta_vs_nocure <- function(posterior_cohort_yearly, posterior_classification, nocure_cohort_long = NULL, nocure_class_long = NULL) {
  out_rows <- list()

  if (!is.null(nocure_cohort_long) && nrow(nocure_cohort_long) > 0L && nrow(posterior_cohort_yearly) > 0L) {
    cohort_join <- posterior_cohort_yearly %>%
      select(dataset, branch, risk_scale, prior_branch, site_prior_family, model_id, latency_family, formula_anchor, horizon_year, meanRisk_Bayes_mean, meanRisk_Bayes_q025, meanRisk_Bayes_q50, meanRisk_Bayes_q975) %>%
      bind_rows(
        posterior_cohort_yearly %>%
          select(dataset, branch, risk_scale, prior_branch, site_prior_family, model_id, latency_family, formula_anchor, horizon_year, meanRisk_Bayes_mean, meanRisk_Bayes_q025, meanRisk_Bayes_q50, meanRisk_Bayes_q975) %>%
          mutate(formula_anchor = "ALL")
      ) %>%
      inner_join(
        nocure_cohort_long %>% filter(metric == "meanRisk") %>% select(dataset, no_cure_model_id, formula_anchor, horizon_year, metric, value),
        by = c("dataset", "formula_anchor", "horizon_year")
      ) %>%
      transmute(
        dataset = dataset,
        branch = branch,
        risk_scale = risk_scale,
        prior_branch = prior_branch,
        site_prior_family = site_prior_family,
        model_id = model_id,
        no_cure_model_id = no_cure_model_id,
        latency_family = latency_family,
        formula_anchor = formula_anchor,
        horizon_year = horizon_year,
        threshold = NA_real_,
        metric = metric,
        no_cure_value = value,
        delta_mean = value - meanRisk_Bayes_mean,
        delta_q025 = value - meanRisk_Bayes_q975,
        delta_q50 = value - meanRisk_Bayes_q50,
        delta_q975 = value - meanRisk_Bayes_q025
      ) %>%
      distinct()
    out_rows[[length(out_rows) + 1L]] <- cohort_join
  }

  if (!is.null(nocure_class_long) && nrow(nocure_class_long) > 0L && nrow(posterior_classification) > 0L) {
    metric_map <- tibble(
      metric = c("false_positive_burden", "FP100", "NB", "PPV", "TPR", "FPR"),
      mean_col = c("false_positive_burden_mean", "FP100_mean", "NB_mean", "PPV_mean", "TPR_mean", "FPR_mean"),
      q025_col = c("false_positive_burden_q025", "FP100_q025", "NB_q025", "PPV_q025", "TPR_q025", "FPR_q025"),
      q50_col = c("false_positive_burden_q50", "FP100_q50", "NB_q50", "PPV_q50", "TPR_q50", "FPR_q50"),
      q975_col = c("false_positive_burden_q975", "FP100_q975", "NB_q975", "PPV_q975", "TPR_q975", "FPR_q975")
    )

    class_long <- bind_rows(lapply(seq_len(nrow(metric_map)), function(i) {
      one <- metric_map[i, , drop = FALSE]
      posterior_classification %>%
        transmute(
          dataset = dataset,
          branch = branch,
          risk_scale = risk_scale,
          prior_branch = prior_branch,
          site_prior_family = site_prior_family,
          model_id = model_id,
          latency_family = latency_family,
          formula_anchor = formula_anchor,
          horizon_year = horizon_year,
          threshold = threshold,
          metric = one$metric[[1L]],
          bayes_mean = .data[[one$mean_col[[1L]]]],
          bayes_q025 = .data[[one$q025_col[[1L]]]],
          bayes_q50 = .data[[one$q50_col[[1L]]]],
          bayes_q975 = .data[[one$q975_col[[1L]]]]
        )
    }))

    class_join <- class_long %>%
      bind_rows(class_long %>% mutate(formula_anchor = "ALL")) %>%
      inner_join(nocure_class_long %>% select(dataset, no_cure_model_id, formula_anchor, horizon_year, threshold, metric, value), by = c("dataset", "formula_anchor", "horizon_year", "threshold", "metric")) %>%
      transmute(
        dataset = dataset,
        branch = branch,
        risk_scale = risk_scale,
        prior_branch = prior_branch,
        site_prior_family = site_prior_family,
        model_id = model_id,
        no_cure_model_id = no_cure_model_id,
        latency_family = latency_family,
        formula_anchor = formula_anchor,
        horizon_year = horizon_year,
        threshold = threshold,
        metric = metric,
        no_cure_value = value,
        delta_mean = value - bayes_mean,
        delta_q025 = value - bayes_q975,
        delta_q50 = value - bayes_q50,
        delta_q975 = value - bayes_q025
      ) %>%
      distinct()

    out_rows[[length(out_rows) + 1L]] <- class_join
  }

  out <- bind_rows_safe(out_rows)
  if (nrow(out) == 0L) {
    return(empty_delta_template())
  }
  out
}

build_cohort_point_metrics <- function(posterior_subject_yearly, ipcw_registry) {
  if (is.null(posterior_subject_yearly) || nrow(posterior_subject_yearly) == 0L) {
    return(tibble())
  }
  bind_rows(lapply(split(posterior_subject_yearly, list(posterior_subject_yearly$model_id, posterior_subject_yearly$horizon_year), drop = TRUE), function(df) {
    if (nrow(df) == 0L) return(NULL)
    dataset_name <- as.character(df$dataset[[1L]])
    horizon_now <- as.integer(df$horizon_year[[1L]])
    horizon_row <- ipcw_registry[[dataset_name]] %>% filter(horizon_year == horizon_now)
    if (nrow(horizon_row) == 0L) return(NULL)
    score <- safe_numeric(df$risk_mean)
    mean_score <- mean(score, na.rm = TRUE)
    tibble(
      dataset = dataset_name,
      model_id = as.character(df$model_id[[1L]]),
      horizon_year = horizon_now,
      auc_ipcw = weighted_auc_horizon(score, horizon_row),
      brier_ipcw = ipcw_brier_horizon(score, horizon_row),
      calibration_in_the_large = as.numeric(horizon_row$observed_km_risk[[1L]]) - mean_score,
      calibration_ratio = ifelse(as.numeric(horizon_row$observed_km_risk[[1L]]) > 0, mean_score / as.numeric(horizon_row$observed_km_risk[[1L]]), NA_real_),
      calibration_abs_error = abs(as.numeric(horizon_row$observed_km_risk[[1L]]) - mean_score),
      observed_km_risk = as.numeric(horizon_row$observed_km_risk[[1L]])
    )
  }))
}

build_prior_branch_delta <- function(model_registry, posterior_cohort_yearly, posterior_classification) {
  model_meta <- model_registry %>%
    select(model_id, dataset, branch, risk_scale, structural_model_id, family_code, formula_anchor, site_prior_family, prior_branch)

  cohort_base <- posterior_cohort_yearly %>%
    left_join(model_meta, by = c("model_id", "dataset", "branch", "risk_scale", "formula_anchor", "site_prior_family"))

  cohort_long <- bind_rows(
    cohort_base %>%
      transmute(
        dataset = dataset,
        branch = branch,
        risk_scale = risk_scale,
        structural_model_id = structural_model_id,
        family_code = family_code,
        formula_anchor = formula_anchor,
        site_prior_family = site_prior_family,
        model_id = model_id,
        prior_branch = prior_branch,
        horizon_year = horizon_year,
        threshold = NA_real_,
        metric = "risk",
        mean = meanRisk_Bayes_mean,
        q025 = meanRisk_Bayes_q025,
        q50 = meanRisk_Bayes_q50,
        q975 = meanRisk_Bayes_q975
      ),
    cohort_base %>%
      transmute(
        dataset = dataset,
        branch = branch,
        risk_scale = risk_scale,
        structural_model_id = structural_model_id,
        family_code = family_code,
        formula_anchor = formula_anchor,
        site_prior_family = site_prior_family,
        model_id = model_id,
        prior_branch = prior_branch,
        horizon_year = horizon_year,
        threshold = NA_real_,
        metric = "cure_fraction",
        mean = cohort_mean_cure_fraction_mean,
        q025 = cohort_mean_cure_fraction_q025,
        q50 = cohort_mean_cure_fraction_q50,
        q975 = cohort_mean_cure_fraction_q975
      )
  )

  metric_map <- tibble(
    metric = c("false_positive_burden", "FP100", "NB", "PPV", "TPR"),
    mean_col = c("false_positive_burden_mean", "FP100_mean", "NB_mean", "PPV_mean", "TPR_mean"),
    q025_col = c("false_positive_burden_q025", "FP100_q025", "NB_q025", "PPV_q025", "TPR_q025"),
    q50_col = c("false_positive_burden_q50", "FP100_q50", "NB_q50", "PPV_q50", "TPR_q50"),
    q975_col = c("false_positive_burden_q975", "FP100_q975", "NB_q975", "PPV_q975", "TPR_q975")
  )
  class_long <- bind_rows(lapply(seq_len(nrow(metric_map)), function(i) {
    one <- metric_map[i, , drop = FALSE]
    posterior_classification %>%
      transmute(
        dataset = dataset,
        branch = branch,
        risk_scale = risk_scale,
        structural_model_id = structural_model_id,
        family_code = family_code,
        formula_anchor = formula_anchor,
        site_prior_family = site_prior_family,
        model_id = model_id,
        prior_branch = prior_branch,
        horizon_year = horizon_year,
        threshold = threshold,
        metric = one$metric[[1L]],
        mean = .data[[one$mean_col[[1L]]]],
        q025 = .data[[one$q025_col[[1L]]]],
        q50 = .data[[one$q50_col[[1L]]]],
        q975 = .data[[one$q975_col[[1L]]]]
      )
  }))

  source_long <- bind_rows(cohort_long, class_long)
  anchor_df <- source_long %>% filter(prior_branch == "anchor_informed") %>% rename_with(~paste0(.x, "_anchor"), all_of(c("model_id", "mean", "q025", "q50", "q975")))
  neutral_df <- source_long %>% filter(prior_branch == "neutral_no_external_info") %>% rename_with(~paste0(.x, "_neutral"), all_of(c("model_id", "mean", "q025", "q50", "q975")))

  out <- anchor_df %>%
    inner_join(
      neutral_df,
      by = c("dataset", "branch", "risk_scale", "structural_model_id", "family_code", "formula_anchor", "site_prior_family", "horizon_year", "threshold", "metric")
    ) %>%
    mutate(
      delta_mean_anchor_minus_neutral = mean_anchor - mean_neutral,
      delta_q025_anchor_minus_neutral = q025_anchor - q975_neutral,
      delta_q50_anchor_minus_neutral = q50_anchor - q50_neutral,
      delta_q975_anchor_minus_neutral = q975_anchor - q025_neutral,
      materiality_threshold = dplyr::case_when(
        metric == "risk" ~ prior_materiality_thresholds$risk,
        metric == "cure_fraction" ~ prior_materiality_thresholds$cure_fraction,
        metric == "false_positive_burden" ~ prior_materiality_thresholds$false_positive_burden,
        metric == "FP100" ~ prior_materiality_thresholds$FP100,
        metric == "NB" ~ prior_materiality_thresholds$NB,
        metric == "PPV" ~ prior_materiality_thresholds$PPV,
        metric == "TPR" ~ prior_materiality_thresholds$TPR,
        TRUE ~ NA_real_
      ),
      conclusion_sign_change_flag = metric == "NB" & is.finite(mean_anchor) & is.finite(mean_neutral) & sign(mean_anchor) != sign(mean_neutral),
      materiality_exceeded_flag = ifelse(is.na(materiality_threshold), FALSE, abs(delta_mean_anchor_minus_neutral) >= materiality_threshold),
      prior_tail_sensitive = materiality_exceeded_flag | conclusion_sign_change_flag
    ) %>%
    transmute(
      dataset = dataset,
      branch = branch,
      risk_scale = risk_scale,
      structural_model_id = structural_model_id,
      family_code = family_code,
      formula_anchor = formula_anchor,
      site_prior_family = site_prior_family,
      horizon_year = horizon_year,
      threshold = threshold,
      metric = metric,
      model_id_anchor = model_id_anchor,
      model_id_neutral = model_id_neutral,
      anchor_mean = mean_anchor,
      neutral_mean = mean_neutral,
      delta_mean_anchor_minus_neutral = delta_mean_anchor_minus_neutral,
      delta_q025_anchor_minus_neutral = delta_q025_anchor_minus_neutral,
      delta_q50_anchor_minus_neutral = delta_q50_anchor_minus_neutral,
      delta_q975_anchor_minus_neutral = delta_q975_anchor_minus_neutral,
      materiality_threshold = materiality_threshold,
      materiality_exceeded_flag = materiality_exceeded_flag,
      conclusion_sign_change_flag = conclusion_sign_change_flag,
      prior_tail_sensitive = prior_tail_sensitive
    )

  if (nrow(out) == 0L) {
    return(empty_prior_branch_delta())
  }
  out
}

extract_parameter_draws_for_anchor_update <- function(model_id, export_dir, coefficient_summary) {
  bundle <- get_reuse_bundle(model_id, export_dir)
  if (inherits(bundle$selected_parameter_draws, "data.frame") && nrow(bundle$selected_parameter_draws) > 0L) {
    return(bundle$selected_parameter_draws)
  }

  coef_tbl <- subset_model_table(coefficient_summary, model_id)
  if (nrow(coef_tbl) == 0L) {
    return(NULL)
  }
  target_params <- coef_tbl %>% filter(parameter %in% c("delta0", paste0("beta_inc[", 1:5, "]")))
  if (nrow(target_params) == 0L) {
    return(NULL)
  }
  n_draws <- 400L
  sim_list <- lapply(seq_len(nrow(target_params)), function(i) {
    tibble::tibble(!!target_params$parameter[[i]] := rnorm(n_draws, mean = safe_numeric(target_params$mean[[i]]), sd = pmax(safe_numeric(target_params$sd[[i]]), 1e-6)))
  })
  bind_cols(sim_list)
}

build_incidence_anchor_update <- function(model_registry, coefficient_summary, export_dir) {
  if (nrow(model_registry) == 0L) {
    return(empty_incidence_anchor_update())
  }

  age_sex_cells <- tibble(
    age_sex_anchor_cell = c("Female_<20", "Female_20_29", "Female_30plus", "Male_<20", "Male_20_29", "Male_30plus"),
    sex_num = c(0, 0, 0, 1, 1, 1),
    age20_29 = c(0, 1, 0, 0, 1, 0),
    age30plus = c(0, 0, 1, 0, 0, 1),
    sex_x_age20_29 = c(0, 0, 0, 0, 1, 0),
    sex_x_age30plus = c(0, 0, 0, 0, 0, 1)
  )
  alpha_anchor <- -9.581369553169
  mu_anchor <- c(0.419871845822, 0.907608052926, 0.586202561451, 0.466865123863, 0.037997248763)
  external_logit <- alpha_anchor + as.matrix(age_sex_cells[, c("sex_num", "age20_29", "age30plus", "sex_x_age20_29", "sex_x_age30plus")]) %*% mu_anchor
  external_one_year_risk <- plogis(external_logit)
  external_incidence_rate_per10k <- -10000 * log(1 - external_one_year_risk)

  bind_rows(lapply(seq_len(nrow(model_registry)), function(i) {
    reg <- model_registry[i, , drop = FALSE]
    param_draws <- extract_parameter_draws_for_anchor_update(reg$model_id[[1L]], export_dir, coefficient_summary)
    if (is.null(param_draws) || nrow(param_draws) == 0L) return(NULL)
    req_cols <- c("delta0", paste0("beta_inc[", 1:5, "]"))
    if (!all(req_cols %in% names(param_draws))) return(NULL)

    alpha_gp <- if (normalize_prior_branch_value(reg$prior_branch[[1L]]) == "anchor_informed") -9.581369553169 else 0
    X_cell <- as.matrix(age_sex_cells[, c("sex_num", "age20_29", "age30plus", "sex_x_age20_29", "sex_x_age30plus")])
    beta_mat <- as.matrix(param_draws[, paste0("beta_inc[", 1:5, "]"), drop = FALSE])
    eta <- sweep(beta_mat %*% t(X_cell), 1, alpha_gp + safe_numeric(param_draws$delta0), FUN = "+")
    prior_center_logit <- as.numeric(alpha_gp + X_cell %*% if (normalize_prior_branch_value(reg$prior_branch[[1L]]) == "anchor_informed") mu_anchor else rep(0, 5L))
    prior_center_risk <- plogis(prior_center_logit)
    q <- matrixStats::colQuantiles(eta, probs = c(0.025, 0.975), na.rm = TRUE)
    risk_q <- matrixStats::colQuantiles(plogis(eta), probs = c(0.025, 0.975), na.rm = TRUE)

    tibble(
      dataset = reg$dataset[[1L]],
      dataset_key = reg$dataset[[1L]],
      branch = reg$branch[[1L]],
      risk_scale = reg$risk_scale[[1L]],
      prior_branch = reg$prior_branch[[1L]],
      site_prior_family = reg$site_prior_family[[1L]],
      model_id = reg$model_id[[1L]],
      retained_fit_id = reg$model_id[[1L]],
      latency_family = reg$latency_family[[1L]],
      formula_anchor = reg$formula_anchor[[1L]],
      age_sex_anchor_cell = age_sex_cells$age_sex_anchor_cell,
      external_incidence_rate_per10k = as.numeric(external_incidence_rate_per10k),
      external_one_year_risk = as.numeric(external_one_year_risk),
      prior_center_logit = prior_center_logit,
      posterior_mean_logit = matrixStats::colMeans2(eta, na.rm = TRUE),
      posterior_lower_logit = q[, 1],
      posterior_upper_logit = q[, 2],
      posterior_mean_one_year_risk = matrixStats::colMeans2(plogis(eta), na.rm = TRUE),
      posterior_lower_one_year_risk = risk_q[, 1],
      posterior_upper_one_year_risk = risk_q[, 2],
      posterior_minus_prior_logit = matrixStats::colMeans2(eta, na.rm = TRUE) - prior_center_logit,
      posterior_minus_prior_risk = matrixStats::colMeans2(plogis(eta), na.rm = TRUE) - prior_center_risk
    )
  }))
}

classify_hazard_shape <- function(hazard_values) {
  hazard_values <- as.numeric(hazard_values)
  if (length(hazard_values) <= 1L || anyNA(hazard_values)) return(NA_character_)
  rng <- max(hazard_values) - min(hazard_values)
  tol <- max(1e-8, 0.02 * max(hazard_values))
  diffs <- diff(hazard_values)
  if (rng < tol) return("approximately_flat")
  if (all(diffs >= -tol)) return("monotone_increasing")
  if (all(diffs <= tol)) return("monotone_decreasing")
  peak_idx <- which.max(hazard_values)
  if (peak_idx > 1L && peak_idx < length(hazard_values)) {
    left_ok <- all(diff(hazard_values[seq_len(peak_idx)]) >= -tol)
    right_ok <- all(diff(hazard_values[peak_idx:length(hazard_values)]) <= tol)
    if (left_ok && right_ok) return("unimodal")
  }
  "complex_nonmonotone"
}

build_hazard_plausibility <- function(model_registry, posterior_cohort_yearly) {
  if (nrow(posterior_cohort_yearly) == 0L) return(empty_hazard_plausibility())
  posterior_cohort_yearly %>%
    group_by(model_id, dataset, branch, risk_scale, prior_branch, site_prior_family, latency_family, formula_anchor) %>%
    summarise(
      hazard_target = "posterior_mean_population_hazard",
      hazard_1y = meanHazard_mean[horizon_year == 1][1],
      hazard_2y = meanHazard_mean[horizon_year == 2][1],
      hazard_3y = meanHazard_mean[horizon_year == 3][1],
      hazard_4y = meanHazard_mean[horizon_year == 4][1],
      hazard_5y = meanHazard_mean[horizon_year == 5][1],
      hazard_6y = meanHazard_mean[horizon_year == 6][1],
      hazard_7y = meanHazard_mean[horizon_year == 7][1],
      hazard_8y = meanHazard_mean[horizon_year == 8][1],
      hazard_9y = meanHazard_mean[horizon_year == 9][1],
      hazard_10y = meanHazard_mean[horizon_year == 10][1],
      hazard_ratio_10y_vs_1y = ifelse(is.finite(hazard_1y) & hazard_1y > 0, hazard_10y / hazard_1y, NA_real_),
      shape_class = classify_hazard_shape(c(hazard_1y, hazard_2y, hazard_3y, hazard_4y, hazard_5y, hazard_6y, hazard_7y, hazard_8y, hazard_9y, hazard_10y)),
      .groups = "drop"
    ) %>%
    left_join(model_registry %>% select(model_id, structural_model_id, family_code, site_placement_label, admissible_flag, admissibility_reasons), by = "model_id")
}

build_uncured_decomposition <- function(posterior_cohort_yearly) {
  if (nrow(posterior_cohort_yearly) == 0L) return(empty_uncured_decomposition())
  posterior_cohort_yearly %>%
    transmute(
      dataset = dataset,
      dataset_key = dataset,
      branch = branch,
      risk_scale = risk_scale,
      prior_branch = prior_branch,
      site_prior_family = site_prior_family,
      model_id = model_id,
      latency_family = latency_family,
      formula_anchor = formula_anchor,
      horizon_year = horizon_year,
      support_tier = support_tier,
      horizon_evidence_class = horizon_evidence_class,
      claim_restriction_flag = claim_restriction_flag,
      cure_fraction_mean = cohort_mean_cure_fraction_mean,
      cure_fraction_q025 = cohort_mean_cure_fraction_q025,
      cure_fraction_q50 = cohort_mean_cure_fraction_q50,
      cure_fraction_q975 = cohort_mean_cure_fraction_q975,
      susceptible_fraction_mean = 1 - cohort_mean_cure_fraction_mean,
      susceptible_fraction_q025 = 1 - cohort_mean_cure_fraction_q975,
      susceptible_fraction_q50 = 1 - cohort_mean_cure_fraction_q50,
      susceptible_fraction_q975 = 1 - cohort_mean_cure_fraction_q025,
      uncured_survival_mean = meanSusceptibleSurvival_mean,
      uncured_survival_q025 = meanSusceptibleSurvival_q025,
      uncured_survival_q50 = meanSusceptibleSurvival_q50,
      uncured_survival_q975 = meanSusceptibleSurvival_q975,
      uncured_risk_mean = 1 - meanSusceptibleSurvival_mean,
      uncured_risk_q025 = 1 - meanSusceptibleSurvival_q975,
      uncured_risk_q50 = 1 - meanSusceptibleSurvival_q50,
      uncured_risk_q975 = 1 - meanSusceptibleSurvival_q025,
      MSTu = NA_real_,
      uncured_mean_support_flag = FALSE,
      uncured_mean_support_note = "Single-number uncured-only mean was not exported because adequate tail support was not guaranteed under the revised support-governance rule."
    )
}

empty_prior_branch_delta <- function() {
  tibble(
    dataset = character(),
    branch = character(),
    risk_scale = character(),
    structural_model_id = character(),
    family_code = character(),
    formula_anchor = character(),
    site_prior_family = character(),
    horizon_year = integer(),
    threshold = double(),
    metric = character(),
    model_id_anchor = character(),
    model_id_neutral = character(),
    anchor_mean = double(),
    neutral_mean = double(),
    delta_mean_anchor_minus_neutral = double(),
    delta_q025_anchor_minus_neutral = double(),
    delta_q50_anchor_minus_neutral = double(),
    delta_q975_anchor_minus_neutral = double(),
    materiality_threshold = double(),
    materiality_exceeded_flag = logical(),
    conclusion_sign_change_flag = logical(),
    prior_tail_sensitive = logical()
  )
}

empty_incidence_anchor_update <- function() {
  tibble(
    dataset = character(), dataset_key = character(), branch = character(), risk_scale = character(), prior_branch = character(), site_prior_family = character(), model_id = character(), retained_fit_id = character(), latency_family = character(), formula_anchor = character(), age_sex_anchor_cell = character(), external_incidence_rate_per10k = double(), external_one_year_risk = double(), prior_center_logit = double(), posterior_mean_logit = double(), posterior_lower_logit = double(), posterior_upper_logit = double(), posterior_mean_one_year_risk = double(), posterior_lower_one_year_risk = double(), posterior_upper_one_year_risk = double(), posterior_minus_prior_logit = double(), posterior_minus_prior_risk = double()
  )
}

empty_hazard_plausibility <- function() {
  tibble(
    model_id = character(), dataset = character(), branch = character(), risk_scale = character(), prior_branch = character(), site_prior_family = character(), latency_family = character(), formula_anchor = character(), hazard_target = character(), hazard_1y = double(), hazard_2y = double(), hazard_3y = double(), hazard_4y = double(), hazard_5y = double(), hazard_6y = double(), hazard_7y = double(), hazard_8y = double(), hazard_9y = double(), hazard_10y = double(), hazard_ratio_10y_vs_1y = double(), shape_class = character(), structural_model_id = character(), family_code = character(), site_placement_label = character(), admissible_flag = logical(), admissibility_reasons = character()
  )
}

empty_uncured_decomposition <- function() {
  tibble(
    dataset = character(), dataset_key = character(), branch = character(), risk_scale = character(), prior_branch = character(), site_prior_family = character(), model_id = character(), latency_family = character(), formula_anchor = character(), horizon_year = integer(), support_tier = character(), horizon_evidence_class = character(), claim_restriction_flag = character(), cure_fraction_mean = double(), cure_fraction_q025 = double(), cure_fraction_q50 = double(), cure_fraction_q975 = double(), susceptible_fraction_mean = double(), susceptible_fraction_q025 = double(), susceptible_fraction_q50 = double(), susceptible_fraction_q975 = double(), uncured_survival_mean = double(), uncured_survival_q025 = double(), uncured_survival_q50 = double(), uncured_survival_q975 = double(), uncured_risk_mean = double(), uncured_risk_q025 = double(), uncured_risk_q50 = double(), uncured_risk_q975 = double(), MSTu = double(), uncured_mean_support_flag = logical(), uncured_mean_support_note = character()
  )
}

empty_stage8A_vs_stage8B_delta <- function() {
  tibble(
    dataset = character(),
    branch_8A = character(),
    branch_8B = character(),
    risk_scale_8A = character(),
    risk_scale_8B = character(),
    model_id_8A = character(),
    model_id_8B = character(),
    horizon_year = integer(),
    threshold = double(),
    metric = character(),
    delta_mean_8B_minus_8A = double(),
    delta_q025_8B_minus_8A = double(),
    delta_q50_8B_minus_8A = double(),
    delta_q975_8B_minus_8A = double(),
    support_tier = character(),
    horizon_evidence_class = character(),
    claim_restriction_flag = character()
  )
}

make_stage8_output_audit <- function(model_grid, model_registry, posterior_cohort_yearly, posterior_classification, ppc_summary, posterior_delta_vs_nocure, diagnostic_pdf_path) {
  model_registry <- simplify_scalar_list_cols(model_registry)
  posterior_cohort_yearly <- simplify_scalar_list_cols(posterior_cohort_yearly)
  posterior_classification <- simplify_scalar_list_cols(posterior_classification)
  ppc_summary <- simplify_scalar_list_cols(ppc_summary)
  posterior_delta_vs_nocure <- simplify_scalar_list_cols(posterior_delta_vs_nocure)

  admissible_n <- sum(as.logical(model_registry$admissible_flag), na.rm = TRUE)
  registry_dup_n <- if ("model_id" %in% names(model_registry) && nrow(model_registry) > 0L) sum(duplicated(model_registry[c("model_id")])) else 0L
  cohort_dup_n <- if (all(c("model_id", "horizon_year") %in% names(posterior_cohort_yearly)) && nrow(posterior_cohort_yearly) > 0L) sum(duplicated(posterior_cohort_yearly[c("model_id", "horizon_year")])) else 0L
  class_dup_n <- if (all(c("model_id", "horizon_year", "threshold") %in% names(posterior_classification)) && nrow(posterior_classification) > 0L) sum(duplicated(posterior_classification[c("model_id", "horizon_year", "threshold")])) else 0L
  ppc_dup_n <- if (all(c("model_id", "horizon_year") %in% names(ppc_summary)) && nrow(ppc_summary) > 0L) sum(duplicated(ppc_summary[c("model_id", "horizon_year")])) else 0L

  bind_rows(
    tibble(check_name = "model_registry_row_count", status = ifelse(nrow(model_registry) == nrow(model_grid), "pass", "fail"), observed_value = as.character(nrow(model_registry)), expected_value = as.character(nrow(model_grid)), detail = "One registry row should exist per model in the Stage 8A grid."),
    tibble(check_name = "model_registry_duplicate_model_id", status = ifelse(registry_dup_n == 0L, "pass", "fail"), observed_value = as.character(registry_dup_n), expected_value = "0", detail = "Model IDs must be unique in the model registry."),
    tibble(check_name = "posterior_cohort_yearly_row_count", status = ifelse(nrow(posterior_cohort_yearly) == admissible_n * 10L, "pass", "fail"), observed_value = as.character(nrow(posterior_cohort_yearly)), expected_value = as.character(admissible_n * 10L), detail = "Admissible models should contribute one cohort-level row per horizon."),
    tibble(check_name = "posterior_classification_duplicate_keys", status = ifelse(class_dup_n == 0L, "pass", "fail"), observed_value = as.character(class_dup_n), expected_value = "0", detail = "Each model-horizon-threshold key should appear at most once in the classification table."),
    tibble(check_name = "posterior_cohort_yearly_duplicate_keys", status = ifelse(cohort_dup_n == 0L, "pass", "fail"), observed_value = as.character(cohort_dup_n), expected_value = "0", detail = "Each model-horizon key should appear at most once in the cohort-level table."),
    tibble(check_name = "ppc_summary_duplicate_keys", status = ifelse(ppc_dup_n == 0L, "pass", "fail"), observed_value = as.character(ppc_dup_n), expected_value = "0", detail = "Each model-horizon key should appear at most once in the PPC table."),
    tibble(check_name = "posterior_delta_vs_nocure_schema_present", status = ifelse(is.data.frame(posterior_delta_vs_nocure), "pass", "fail"), observed_value = ifelse(is.data.frame(posterior_delta_vs_nocure), "TRUE", "FALSE"), expected_value = "TRUE", detail = "Delta-vs-no-cure table should exist, even if empty."),
    tibble(check_name = "diagnostic_pdf_exists", status = ifelse(pdf_file_is_usable(diagnostic_pdf_path), "pass", "fail"), observed_value = ifelse(pdf_file_is_usable(diagnostic_pdf_path), "TRUE", "FALSE"), expected_value = "TRUE", detail = "Diagnostic PDF should exist after the run unless preserved intentionally after a plotting error.")
  )
}

save_plot_object <- function(plot_obj, png_path) {
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  tmp <- make_temp_output_path(png_path, tag = "tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  ggplot2::ggsave(filename = tmp, plot = plot_obj, width = plot_width, height = plot_height, dpi = plot_dpi)
  safe_promote_file(tmp, png_path)
}

safe_generate_diagnostic_pdf <- function(trace_records, plot_objects, final_path) {
  dir.create(dirname(final_path), recursive = TRUE, showWarnings = FALSE)
  tmp_pdf <- make_temp_output_path(final_path, tag = "tmp")
  pdf_open <- FALSE
  close_pdf <- function() {
    if (isTRUE(pdf_open)) {
      try(grDevices::dev.off(), silent = TRUE)
      pdf_open <<- FALSE
    }
  }
  on.exit({ close_pdf(); if (file.exists(tmp_pdf)) unlink(tmp_pdf) }, add = TRUE)
  grDevices::pdf(tmp_pdf, width = plot_width, height = plot_height, onefile = TRUE)
  pdf_open <- TRUE
  if (length(trace_records) > 0L) {
    for (trace_record in trace_records) plot_trace_record(trace_record)
  }
  if (length(plot_objects) > 0L) {
    for (one in plot_objects) print(one)
  }
  if (length(trace_records) == 0L && length(plot_objects) == 0L) {
    plot.new(); text(0.5, 0.5, "No diagnostic pages were available for this Stage 8A run.")
  }
  close_pdf()
  if (!pdf_file_is_usable(tmp_pdf)) stop("Temporary diagnostic PDF was not created correctly.", call. = FALSE)
  safe_promote_file(tmp_pdf, final_path)
  invisible(TRUE)
}
# 🔴 Load: Stage-1 backbone, optional comparators, and frozen horizon metadata ===============================
## 🟠 Load: analysis-ready datasets and cross-stage linkage tables ===============================
backbone_objects <- load_stage1_backbone_or_fallback()
analysis_datasets <- backbone_objects$analysis_datasets
scaling_registry <- tibble::as_tibble(backbone_objects$scaling_registry)
formula_registry <- tibble::as_tibble(backbone_objects$formula_registry)
horizon_registry <- tibble::as_tibble(backbone_objects$horizon_registry)
threshold_registry <- tibble::as_tibble(backbone_objects$threshold_registry)
metadata_registry_stage1 <- tibble::as_tibble(backbone_objects$metadata_registry)
dataset_registry <- tibble::as_tibble(backbone_objects$dataset_registry)

horizons_year <- sort(unique(as.integer(horizon_registry$horizon_year)))
risk_thresholds <- sort(unique(as.numeric(threshold_registry$threshold)))
if (!identical(horizons_year, 1:10)) {
  stop("Stage 8A requires the common horizon vector 1:10 years.", call. = FALSE)
}
if (length(risk_thresholds) == 0L || anyNA(risk_thresholds) || any(risk_thresholds <= 0 | risk_thresholds >= 1)) {
  stop("Threshold registry must contain probabilities strictly between 0 and 1.", call. = FALSE)
}

screening_flags <- read_screening_flags(stage6_screening_flag_csv)
nocure_cohort_long <- normalize_nocure_cohort(stage5_nocure_cohort_csv)
nocure_class_long <- normalize_nocure_classification(stage5_nocure_classification_csv)

if (nrow(screening_flags) == 0L) {
  warning("Stage 6 screening linkage was not recovered. Carry-forward screening fields will remain missing.", call. = FALSE)
}
if (is.null(nocure_cohort_long) || nrow(nocure_cohort_long) == 0L) {
  warning("Stage 5 cohort comparison linkage was not recovered. `bayes_stage8_posterior_delta_vs_nocure.csv` may be partially or fully empty.", call. = FALSE)
}
if (is.null(nocure_class_long) || nrow(nocure_class_long) == 0L) {
  warning("Stage 5 threshold comparison linkage was not recovered. `bayes_stage8_posterior_delta_vs_nocure.csv` may be partially or fully empty.", call. = FALSE)
}

ipcw_registry <- lapply(names(analysis_datasets), function(ds) build_ipcw_reference(analysis_datasets[[ds]], horizons_year))
names(ipcw_registry) <- names(analysis_datasets)

horizon_metadata_registry <- bind_rows(lapply(names(ipcw_registry), function(ds) {
  ipcw_registry[[ds]] %>% mutate(dataset = ds)
})) %>%
  left_join(horizon_registry, by = c("dataset", "horizon_year")) %>%
  mutate(
    dataset_key = dataset,
    support_tier = as.character(support_tier),
    horizon_evidence_class = as.character(horizon_evidence_class),
    claim_restriction_flag = as.character(claim_restriction_flag),
    classification_estimable_flag = (denom_case > 0) & (denom_control > 0),
    support_priority = dplyr::case_when(
      support_tier == "primary_supported" ~ 1L,
      support_tier == "secondary" ~ 2L,
      support_tier == "sensitivity" ~ 3L,
      TRUE ~ 4L
    )
  ) %>%
  arrange(factor(dataset, levels = c("PNU", "SNU", "merged")), horizon_year)

# 🔴 Prepare: Stage-8A model grid, reuse plan, and Stan availability ===============================
## 🟠 Prepare: grid expansion for prior branches and site-prior families ===============================
prior_specs <- build_prior_specs()
model_grid <- build_model_grid(include_merged_incidence_site_supplementary)
if (!is.null(run_model_ids)) {
  model_grid <- model_grid %>% filter(model_id %in% run_model_ids | legacy_model_id %in% run_model_ids)
  if (nrow(model_grid) == 0L) {
    stop("`run_model_ids` filtered out all Stage 8A models.", call. = FALSE)
  }
}
model_order <- model_grid$model_id
screening_lookup <- build_screening_model_lookup(screening_flags, model_grid)

existing_exports <- if (isTRUE(reuse_existing_stage8_outputs)) load_existing_stage8_exports(export_path) else list(
  model_registry = NULL,
  coefficient_summary = NULL,
  posterior_subject_profile = NULL,
  posterior_subject_yearly = NULL,
  posterior_cohort_yearly = NULL,
  posterior_classification = NULL,
  posterior_delta_vs_nocure = NULL,
  diagnostics_parameter_level = NULL,
  ppc_summary = NULL,
  prior_predictive_summary = NULL,
  prior_branch_delta = NULL,
  incidence_anchor_update = NULL,
  hazard_plausibility = NULL,
  uncured_decomposition = NULL,
  stage8A_vs_stage8B_delta = NULL,
  diagnostic_pdf_exists = FALSE
)

reuse_flags <- vapply(seq_len(nrow(model_grid)), function(ii) {
  model_row <- model_grid[ii, , drop = FALSE]
  dataset_df <- analysis_datasets[[model_row$dataset[[1L]]]]
  is_model_reusable(model_row, dataset_df, existing_exports, export_path, require_existing_rds_to_skip_fit)
}, logical(1))
model_grid$reuse_existing_fit <- reuse_flags
n_models_reused <- sum(model_grid$reuse_existing_fit)
n_models_to_fit <- sum(!model_grid$reuse_existing_fit)
message("Stage 8A reuse plan: ", n_models_reused, " model(s) reused; ", n_models_to_fit, " model(s) require fitting.")

need_fit_packages <- n_models_to_fit > 0L
if (need_fit_packages) {
  fit_packages <- c("rstan", "posterior", "loo")
  missing_fit_packages <- fit_packages[!vapply(fit_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_fit_packages) > 0L) {
    stop("Install required fitting packages before running this script: ", paste(missing_fit_packages, collapse = ", "), call. = FALSE)
  }
  suppressPackageStartupMessages({
    library(rstan)
    library(posterior)
  })
  rstan_options(auto_write = TRUE)
  options(mc.cores = max(1L, min(stan_chains, parallel::detectCores(logical = TRUE))))
  stan_model_compiled <- compile_stage8_stan_model()
} else {
  stan_model_compiled <- NULL
}

diagnostic_pdf_path <- file.path(export_path, "bayes_stage8_diagnostic_plots.pdf")
trace_records <- list()

# 🔴 Run: Stage-8A model loop with reuse-first execution ===============================
## 🟠 Run: fit or reuse each expanded Stage-8A specification ===============================
main_registry_rows <- list()
coef_rows <- list()
subject_profile_rows <- list()
subject_yearly_rows <- list()
cohort_yearly_rows <- list()
class_rows <- list()
diagnostic_rows <- list()
ppc_rows_all <- list()
prior_predictive_rows <- list()
total_models <- nrow(model_grid)

withCallingHandlers({
for (ii in seq_len(nrow(model_grid))) {
  model_row <- model_grid[ii, , drop = FALSE]
  model_id_now <- model_row$model_id[[1L]]
  dataset_name <- model_row$dataset[[1L]]
  dataset_df <- analysis_datasets[[dataset_name]]
  model_started_at <- Sys.time()

  emit_stage8_progress(ii - 1L, total_models, model_id_now, paste0("starting Stage 8A model (dataset=", dataset_name, ", family=", model_row$family_code[[1L]], ", prior=", model_row$prior_branch[[1L]], ", site_prior=", model_row$site_prior_family[[1L]], ", reuse=", isTRUE(model_row$reuse_existing_fit[[1L]]), ")"))

  screening_row <- screening_lookup %>% filter(model_id == model_id_now)
  if (nrow(screening_row) == 0L) {
    screening_row <- tibble(
      model_id = model_id_now,
      cure_model_eligibility_flag = NA_character_,
      primary_gate_method = NA_character_,
      primary_gate_flag = NA_character_,
      receus_primary_class = NA_character_,
      presence_modifier_flag = NA_character_,
      cure_presence_support_flag = NA_character_,
      followup_contradiction_flag = NA_character_,
      followup_not_contradicted_flag = NA_character_,
      screening_note = NA_character_,
      screening_flag = NA_character_,
      screening_detail = NA_character_,
      carry_forward_stage8 = NA
    )
  }

  existing_reg <- subset_model_registry_row(existing_exports$model_registry, model_id_now)
  reuse_bundle <- if (isTRUE(model_row$reuse_existing_fit[[1L]])) get_reuse_bundle(model_id_now, export_path, existing_reg) else empty_reuse_bundle()

  prior_spec <- prior_specs[[model_row$prior_branch[[1L]]]]
  existing_prior_model <- normalize_prior_predictive_summary(bind_rows_safe(list(subset_model_table(existing_exports$prior_predictive_summary, model_id_now), subset_model_table(reuse_bundle$prior_predictive_summary, model_id_now))))
  need_prior_regeneration <- nrow(existing_prior_model) == 0L
  design_main <- make_design_bundle(dataset_df, model_row, prior_spec, snu_site_label)

  if (need_prior_regeneration) {
    set.seed(stan_seed + ii)
    prior_pred_main <- simulate_prior_predictive(df = dataset_df, design_bundle = design_main, model_row = model_row, prior_spec = prior_spec, n_draws = prior_predictive_draws, horizons_eval = c(1, 2, 5, 10))
    prior_predictive_rows[[length(prior_predictive_rows) + 1L]] <- prior_pred_main$summary
  } else {
    prior_pred_main <- list(degeneracy = tibble(degenerate_flag = any(existing_prior_model$prior_degenerate_flag, na.rm = TRUE)))
    prior_predictive_rows[[length(prior_predictive_rows) + 1L]] <- existing_prior_model
  }

  if (isTRUE(model_row$reuse_existing_fit[[1L]])) {
    if (!is.null(existing_reg)) {
      main_registry_rows[[length(main_registry_rows) + 1L]] <- existing_reg
    }
    coef_rows[[length(coef_rows) + 1L]] <- pick_reuse_table(model_id_now, existing_exports$coefficient_summary, reuse_bundle$coefficient_summary)
    diagnostic_rows[[length(diagnostic_rows) + 1L]] <- pick_reuse_table(model_id_now, existing_exports$diagnostics_parameter_level, reuse_bundle$diagnostics_parameter_level)
    ppc_rows_all[[length(ppc_rows_all) + 1L]] <- pick_reuse_table(model_id_now, existing_exports$ppc_summary, reuse_bundle$ppc_summary)

    existing_reg_admissible <- !is.null(existing_reg) && isTRUE(as.logical(registry_field_value(existing_reg, "admissible_flag", FALSE)))
    if (isTRUE(existing_reg_admissible)) {
      subject_profile_rows[[length(subject_profile_rows) + 1L]] <- pick_reuse_table(model_id_now, existing_exports$posterior_subject_profile, reuse_bundle$posterior_subject_profile)
      subject_yearly_rows[[length(subject_yearly_rows) + 1L]] <- pick_reuse_table(model_id_now, existing_exports$posterior_subject_yearly, reuse_bundle$posterior_subject_yearly)
      cohort_yearly_rows[[length(cohort_yearly_rows) + 1L]] <- pick_reuse_table(model_id_now, existing_exports$posterior_cohort_yearly, reuse_bundle$posterior_cohort_yearly)
      class_rows[[length(class_rows) + 1L]] <- pick_reuse_table(model_id_now, existing_exports$posterior_classification, reuse_bundle$posterior_classification)
    }

    emit_stage8_progress(ii, total_models, model_id_now, paste0("reused existing outputs; elapsed=", format_stage8_number(elapsed_stage8_seconds(model_started_at), 1L), "s"))
    next
  }

  stan_data <- list(
    N = nrow(dataset_df),
    time = as.numeric(design_main$time),
    event = as.integer(design_main$event),
    K_inc = ncol(design_main$X_inc),
    X_inc = design_main$X_inc,
    K_lat = ncol(design_main$X_lat),
    X_lat = design_main$X_lat,
    family_id = as.integer(model_row$family_id[[1L]]),
    alpha_gp = as.numeric(design_main$alpha_gp),
    mu_beta_inc = as.numeric(design_main$mu_beta_inc),
    sd_beta_inc = as.numeric(design_main$sd_beta_inc),
    sd_delta = as.numeric(prior_spec$sd_delta),
    sd_gamma0 = as.numeric(prior_spec$sd_gamma0),
    sd_gamma_lat = as.numeric(design_main$sd_gamma_lat),
    sd_shape = as.numeric(design_main$shape_prior_sd),
    incidence_site_position = as.integer(design_main$incidence_site_position),
    latency_site_position = as.integer(design_main$latency_site_position),
    use_student_t_site_prior = as.integer(design_main$use_student_t_site_prior),
    site_prior_df = as.numeric(design_main$site_prior_df),
    site_prior_scale = as.numeric(design_main$site_prior_scale)
  )

  fit_status <- "ok"
  fit_error_message <- NA_character_
  fit <- tryCatch(
    rstan::sampling(
      object = stan_model_compiled,
      data = stan_data,
      chains = stan_chains,
      iter = stan_iter,
      warmup = stan_warmup,
      thin = stan_thin,
      seed = stan_seed + ii,
      refresh = stan_refresh,
      control = list(adapt_delta = stan_adapt_delta, max_treedepth = stan_max_treedepth)
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    fit_status <- "sampling_error"
    fit_error_message <- conditionMessage(fit)
    main_registry_rows[[length(main_registry_rows) + 1L]] <- tibble(
      dataset = dataset_name,
      dataset_key = dataset_name,
      model_id = model_id_now,
      legacy_model_id = model_row$legacy_model_id[[1L]],
      branch = "Stage8A",
      risk_scale = "transition_only_main",
      prior_branch = model_row$prior_branch[[1L]],
      site_prior_family = model_row$site_prior_family[[1L]],
      site_placement_label = model_row$site_placement_label[[1L]],
      structural_model_id = model_row$structural_model_id[[1L]],
      latency_family = model_row$latency_family[[1L]],
      family_code = model_row$family_code[[1L]],
      formula_anchor = model_row$formula_anchor[[1L]],
      incidence_site_indicator = model_row$incidence_site_indicator[[1L]],
      latency_site_indicator = model_row$latency_site_indicator[[1L]],
      latency_interaction_indicator = model_row$latency_interaction_indicator[[1L]],
      is_supplementary_branch = model_row$is_supplementary_branch[[1L]],
      fit_status = fit_status,
      fit_error_message = fit_error_message,
      admissible_flag = FALSE,
      admissibility_reasons = "sampling_error",
      prior_degenerate_flag = prior_pred_main$degeneracy$degenerate_flag[[1L]],
      posterior_degenerate_flag = NA,
      ppc_gross_contradiction_flag = NA,
      divergences = NA_integer_,
      max_rhat = NA_real_,
      min_bulk_ess = NA_real_,
      min_tail_ess = NA_real_,
      treedepth_exceeded = NA_integer_,
      waic = NA_real_,
      looic = NA_real_,
      p_waic = NA_real_,
      p_waic_high_n = NA_integer_,
      p_waic_high_pct = NA_real_,
      p_loo = NA_real_,
      pareto_k_max = NA_real_,
      pareto_k_bad_n = NA_integer_,
      pareto_k_bad_pct = NA_real_,
      pareto_k_very_bad_n = NA_integer_,
      waic_warning_flag = NA,
      loo_warning_flag = NA,
      info_criteria_warning_detail = NA_character_,
      n = nrow(dataset_df),
      n_event = sum(dataset_df$event_main),
      n_censor_main = sum(dataset_df$censor_main),
      n_remission = sum(dataset_df$remission_flag),
      rds_path = NA_character_
    )
    emit_stage8_progress(ii, total_models, model_id_now, paste0("sampling error after ", format_stage8_number(elapsed_stage8_seconds(model_started_at), 1L), "s: ", fit_error_message))
    next
  }

  trace_records[[length(trace_records) + 1L]] <- make_trace_record(fit, model_id_now, model_row$family_code[[1L]], ncol(design_main$X_inc), ncol(design_main$X_lat))

  param_names <- c("delta0", "gamma0", paste0("beta_inc[", seq_len(ncol(design_main$X_inc)), "]"), paste0("gamma_lat[", seq_len(ncol(design_main$X_lat)), "]"))
  if (model_row$family_code[[1L]] == "W") param_names <- c(param_names, "rho_W")
  if (model_row$family_code[[1L]] == "LN") param_names <- c(param_names, "log_sigma_LN")
  if (model_row$family_code[[1L]] == "LL") param_names <- c(param_names, "psi_LL")

  param_array <- posterior::as_draws_array(as.array(fit, pars = param_names))
  param_diag_tbl <- posterior::summarise_draws(param_array, mean = base::mean, sd = stats::sd, rhat = posterior::rhat, ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail)
  param_draws_mat <- posterior::as_draws_matrix(param_array)

  coef_tbl <- tibble(
    dataset = dataset_name,
    dataset_key = dataset_name,
    model_id = model_id_now,
    parameter = colnames(param_draws_mat),
    mean = apply(param_draws_mat, 2, mean),
    sd = apply(param_draws_mat, 2, stats::sd),
    q025 = apply(param_draws_mat, 2, stats::quantile, probs = 0.025, names = FALSE),
    q50 = apply(param_draws_mat, 2, stats::quantile, probs = 0.500, names = FALSE),
    q975 = apply(param_draws_mat, 2, stats::quantile, probs = 0.975, names = FALSE)
  )
  coef_rows[[length(coef_rows) + 1L]] <- coef_tbl

  diagnostic_param_tbl_model <- tibble(
    dataset = dataset_name,
    dataset_key = dataset_name,
    model_id = model_id_now,
    parameter = param_diag_tbl$variable,
    mean = param_diag_tbl$mean,
    sd = param_diag_tbl$sd,
    rhat = param_diag_tbl$rhat,
    ess_bulk = param_diag_tbl$ess_bulk,
    ess_tail = param_diag_tbl$ess_tail
  )

  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- sum(vapply(sampler_params, function(x) sum(x[, "divergent__"]), numeric(1)))
  treedepth_exceeded <- sum(vapply(sampler_params, function(x) sum(x[, "treedepth__"] >= stan_max_treedepth), numeric(1)))

  draws_compact <- extract_draws_compact(fit, K_inc = ncol(design_main$X_inc), K_lat = ncol(design_main$X_lat))
  total_draws <- length(draws_compact$delta0)
  set.seed(stan_seed + 200000L + ii)
  keep_draw_idx <- if (total_draws <= posterior_prediction_draws) seq_len(total_draws) else sort(sample(seq_len(total_draws), size = posterior_prediction_draws, replace = FALSE))

  draws_pred <- list(
    delta0 = draws_compact$delta0[keep_draw_idx],
    beta_inc = draws_compact$beta_inc[keep_draw_idx, , drop = FALSE],
    gamma0 = draws_compact$gamma0[keep_draw_idx],
    gamma_lat = draws_compact$gamma_lat[keep_draw_idx, , drop = FALSE],
    rho_W = draws_compact$rho_W[keep_draw_idx],
    log_sigma_LN = draws_compact$log_sigma_LN[keep_draw_idx],
    psi_LL = draws_compact$psi_LL[keep_draw_idx]
  )

  linear_terms <- compute_linear_terms(draws_pred, design_main$X_inc, design_main$X_lat, design_main$alpha_gp)
  supported_horizons <- if (dataset_name == "PNU") c(1, 2) else c(1, 2, 5)
  supported_risk_list <- lapply(supported_horizons, function(h) {
    fh_sup <- family_survival_hazard(h, model_row$family_code[[1L]], linear_terms$mu_lat_mat, linear_terms$median_mat, draws_pred)
    1 - ((1 - linear_terms$pi_mat) + linear_terms$pi_mat * fh_sup$Su)
  })
  posterior_degeneracy <- compute_degeneracy(linear_terms$pi_mat, linear_terms$median_mat, supported_risk_list)
  info_criteria <- compute_information_criteria(draws_compact$log_lik)

  subject_profile_summary <- summarize_cols_matrix(linear_terms$cure_prob_mat)
  susceptible_prob_summary <- summarize_cols_matrix(linear_terms$pi_mat)
  median_susc_summary <- summarize_cols_matrix(linear_terms$median_mat)
  subject_profile_tbl <- bind_cols(
    tibble(dataset = dataset_name, dataset_key = dataset_name, model_id = model_id_now),
    design_main$id_df,
    tibble(
      cure_prob_mean = subject_profile_summary$mean,
      cure_prob_q025 = subject_profile_summary$q025,
      cure_prob_q50 = subject_profile_summary$q50,
      cure_prob_q975 = subject_profile_summary$q975,
      susceptible_prob_mean = susceptible_prob_summary$mean,
      susceptible_prob_q025 = susceptible_prob_summary$q025,
      susceptible_prob_q50 = susceptible_prob_summary$q50,
      susceptible_prob_q975 = susceptible_prob_summary$q975,
      median_susc_time_mean = median_susc_summary$mean,
      median_susc_time_q025 = median_susc_summary$q025,
      median_susc_time_q50 = median_susc_summary$q50,
      median_susc_time_q975 = median_susc_summary$q975
    )
  )

  ppc_model_rows <- list(); cohort_model_rows <- list(); subject_year_model_rows <- list(); class_model_rows <- list()

  for (h in horizons_year) {
    fh <- family_survival_hazard(h, model_row$family_code[[1L]], linear_terms$mu_lat_mat, linear_terms$median_mat, draws_pred)
    pop_surv <- (1 - linear_terms$pi_mat) + linear_terms$pi_mat * fh$Su
    risk_mat <- 1 - pop_surv
    meta_row <- horizon_metadata_registry %>% filter(dataset == dataset_name, horizon_year == h)
    if (nrow(meta_row) == 0L) {
      meta_row <- tibble(dataset = dataset_name, horizon_year = h, support_tier = derive_support_tier(dataset_name, h), horizon_evidence_class = derive_horizon_evidence_class(dataset_name, h), claim_restriction_flag = derive_claim_restriction_flag(derive_horizon_evidence_class(dataset_name, h)), interpretation_note = derive_interpretation_note(dataset_name, h, derive_support_tier(dataset_name, h), derive_horizon_evidence_class(dataset_name, h), derive_claim_restriction_flag(derive_horizon_evidence_class(dataset_name, h))), observed_km_risk = NA_real_, denom_case = NA_real_, denom_control = NA_real_, classification_estimable_flag = NA)
    }
    horizon_ref <- ipcw_registry[[dataset_name]] %>% filter(horizon_year == h)

    subj_surv_summary <- summarize_cols_matrix(pop_surv)
    subj_risk_summary <- summarize_cols_matrix(risk_mat)
    subj_su_summary <- summarize_cols_matrix(fh$Su)
    point_risk_subject <- matrixStats::colMeans2(risk_mat, na.rm = TRUE)
    mean_risk_draw <- rowMeans(risk_mat)
    mean_surv_draw <- rowMeans(pop_surv)
    mean_sus_surv_draw <- rowMeans(fh$Su)
    mean_hazard_draw <- rowMeans(fh$haz)
    mean_cure_draw <- rowMeans(linear_terms$cure_prob_mat)

    subject_year_model_rows[[length(subject_year_model_rows) + 1L]] <- bind_cols(
      tibble(dataset = dataset_name, dataset_key = dataset_name, branch = "Stage8A", risk_scale = "transition_only_main", prior_branch = model_row$prior_branch[[1L]], site_prior_family = model_row$site_prior_family[[1L]], model_id = model_id_now, latency_family = model_row$latency_family[[1L]], formula_anchor = model_row$formula_anchor[[1L]], horizon_year = h, support_tier = meta_row$support_tier[[1L]], horizon_evidence_class = meta_row$horizon_evidence_class[[1L]], claim_restriction_flag = meta_row$claim_restriction_flag[[1L]]),
      design_main$id_df,
      tibble(
        S_pop_mean = subj_surv_summary$mean,
        S_pop_q025 = subj_surv_summary$q025,
        S_pop_q50 = subj_surv_summary$q50,
        S_pop_q975 = subj_surv_summary$q975,
        risk_mean = subj_risk_summary$mean,
        risk_q025 = subj_risk_summary$q025,
        risk_q50 = subj_risk_summary$q50,
        risk_q975 = subj_risk_summary$q975,
        S_u_mean = subj_su_summary$mean,
        S_u_q025 = subj_su_summary$q025,
        S_u_q50 = subj_su_summary$q50,
        S_u_q975 = subj_su_summary$q975
      )
    )

    cohort_model_rows[[length(cohort_model_rows) + 1L]] <- tibble(
      dataset = dataset_name,
      dataset_key = dataset_name,
      branch = "Stage8A",
      risk_scale = "transition_only_main",
      prior_branch = model_row$prior_branch[[1L]],
      site_prior_family = model_row$site_prior_family[[1L]],
      model_id = model_id_now,
      latency_family = model_row$latency_family[[1L]],
      formula_anchor = model_row$formula_anchor[[1L]],
      horizon_year = h,
      support_tier = meta_row$support_tier[[1L]],
      horizon_evidence_class = meta_row$horizon_evidence_class[[1L]],
      claim_restriction_flag = meta_row$claim_restriction_flag[[1L]],
      interpretation_note = meta_row$interpretation_note[[1L]],
      observed_km_risk = horizon_ref$observed_km_risk[[1L]],
      meanRisk_Bayes_mean = mean(mean_risk_draw),
      meanRisk_Bayes_q025 = stats::quantile(mean_risk_draw, 0.025, names = FALSE),
      meanRisk_Bayes_q50 = stats::quantile(mean_risk_draw, 0.500, names = FALSE),
      meanRisk_Bayes_q975 = stats::quantile(mean_risk_draw, 0.975, names = FALSE),
      meanSurvival_mean = mean(mean_surv_draw),
      meanSurvival_q025 = stats::quantile(mean_surv_draw, 0.025, names = FALSE),
      meanSurvival_q50 = stats::quantile(mean_surv_draw, 0.500, names = FALSE),
      meanSurvival_q975 = stats::quantile(mean_surv_draw, 0.975, names = FALSE),
      meanSusceptibleSurvival_mean = mean(mean_sus_surv_draw),
      meanSusceptibleSurvival_q025 = stats::quantile(mean_sus_surv_draw, 0.025, names = FALSE),
      meanSusceptibleSurvival_q50 = stats::quantile(mean_sus_surv_draw, 0.500, names = FALSE),
      meanSusceptibleSurvival_q975 = stats::quantile(mean_sus_surv_draw, 0.975, names = FALSE),
      meanHazard_mean = mean(mean_hazard_draw),
      meanHazard_q025 = stats::quantile(mean_hazard_draw, 0.025, names = FALSE),
      meanHazard_q50 = stats::quantile(mean_hazard_draw, 0.500, names = FALSE),
      meanHazard_q975 = stats::quantile(mean_hazard_draw, 0.975, names = FALSE),
      cohort_mean_cure_fraction_mean = mean(mean_cure_draw),
      cohort_mean_cure_fraction_q025 = stats::quantile(mean_cure_draw, 0.025, names = FALSE),
      cohort_mean_cure_fraction_q50 = stats::quantile(mean_cure_draw, 0.500, names = FALSE),
      cohort_mean_cure_fraction_q975 = stats::quantile(mean_cure_draw, 0.975, names = FALSE),
      auc_ipcw = weighted_auc_horizon(point_risk_subject, horizon_ref),
      brier_ipcw = ipcw_brier_horizon(point_risk_subject, horizon_ref),
      calibration_in_the_large = horizon_ref$observed_km_risk[[1L]] - mean(point_risk_subject, na.rm = TRUE),
      calibration_ratio = ifelse(horizon_ref$observed_km_risk[[1L]] > 0, mean(point_risk_subject, na.rm = TRUE) / horizon_ref$observed_km_risk[[1L]], NA_real_),
      calibration_abs_error = abs(horizon_ref$observed_km_risk[[1L]] - mean(point_risk_subject, na.rm = TRUE))
    )

    ppc_model_rows[[length(ppc_model_rows) + 1L]] <- tibble(
      dataset = dataset_name,
      dataset_key = dataset_name,
      branch = "Stage8A",
      risk_scale = "transition_only_main",
      prior_branch = model_row$prior_branch[[1L]],
      site_prior_family = model_row$site_prior_family[[1L]],
      model_id = model_id_now,
      horizon_year = h,
      support_tier = meta_row$support_tier[[1L]],
      horizon_evidence_class = meta_row$horizon_evidence_class[[1L]],
      claim_restriction_flag = meta_row$claim_restriction_flag[[1L]],
      observed_km_risk = horizon_ref$observed_km_risk[[1L]],
      posterior_mean_risk = mean(mean_risk_draw),
      posterior_q025_risk = stats::quantile(mean_risk_draw, 0.025, names = FALSE),
      posterior_q975_risk = stats::quantile(mean_risk_draw, 0.975, names = FALSE),
      absolute_difference = abs(mean(mean_risk_draw) - horizon_ref$observed_km_risk[[1L]]),
      gross_contradiction_flag = ((h %in% supported_horizons) && (horizon_ref$observed_km_risk[[1L]] < stats::quantile(mean_risk_draw, 0.025, names = FALSE) || horizon_ref$observed_km_risk[[1L]] > stats::quantile(mean_risk_draw, 0.975, names = FALSE)) && abs(mean(mean_risk_draw) - horizon_ref$observed_km_risk[[1L]]) > ppc_tolerance_abs)
    )

    class_out <- compute_classification_summary(risk_draws = risk_mat, horizon_row = horizon_ref, thresholds = risk_thresholds)
    if (nrow(class_out) > 0L) {
      class_model_rows[[length(class_model_rows) + 1L]] <- class_out %>%
        mutate(
          dataset = dataset_name,
          dataset_key = dataset_name,
          branch = "Stage8A",
          risk_scale = "transition_only_main",
          prior_branch = model_row$prior_branch[[1L]],
          site_prior_family = model_row$site_prior_family[[1L]],
          model_id = model_id_now,
          latency_family = model_row$latency_family[[1L]],
          formula_anchor = model_row$formula_anchor[[1L]],
          horizon_year = h,
          support_tier = meta_row$support_tier[[1L]],
          horizon_evidence_class = meta_row$horizon_evidence_class[[1L]],
          claim_restriction_flag = meta_row$claim_restriction_flag[[1L]],
          interpretation_note = meta_row$interpretation_note[[1L]],
          observed_km_risk = horizon_ref$observed_km_risk[[1L]],
          denom_case = horizon_ref$denom_case[[1L]],
          denom_control = horizon_ref$denom_control[[1L]],
          classification_estimable_flag = meta_row$classification_estimable_flag[[1L]]
        ) %>%
        relocate(dataset, dataset_key, branch, risk_scale, prior_branch, site_prior_family, model_id, latency_family, formula_anchor, horizon_year, threshold)
    }
  }

  ppc_model_tbl <- bind_rows(ppc_model_rows)
  cohort_model_tbl <- bind_rows(cohort_model_rows)
  subject_year_model_tbl <- bind_rows(subject_year_model_rows)
  class_model_tbl <- bind_rows(class_model_rows)
  ppc_gross_contradiction_flag <- any(ppc_model_tbl$gross_contradiction_flag, na.rm = TRUE)

  max_rhat <- max(param_diag_tbl$rhat, na.rm = TRUE)
  min_bulk_ess <- min(param_diag_tbl$ess_bulk, na.rm = TRUE)
  min_tail_ess <- min(param_diag_tbl$ess_tail, na.rm = TRUE)

  admissibility_reasons <- character()
  if (prior_pred_main$degeneracy$degenerate_flag[[1L]]) admissibility_reasons <- c(admissibility_reasons, "prior_degenerate")
  if (posterior_degeneracy$degenerate_flag[[1L]]) admissibility_reasons <- c(admissibility_reasons, "posterior_degenerate")
  if (!is.finite(max_rhat) || max_rhat >= 1.01) admissibility_reasons <- c(admissibility_reasons, "rhat")
  if (!is.finite(min_bulk_ess) || min_bulk_ess < ess_min_threshold) admissibility_reasons <- c(admissibility_reasons, "bulk_ess")
  if (!is.finite(min_tail_ess) || min_tail_ess < ess_min_threshold) admissibility_reasons <- c(admissibility_reasons, "tail_ess")
  if (divergences > 0) admissibility_reasons <- c(admissibility_reasons, "divergences")
  if (treedepth_exceeded > 0) admissibility_reasons <- c(admissibility_reasons, "treedepth")
  if (isTRUE(ppc_gross_contradiction_flag)) admissibility_reasons <- c(admissibility_reasons, "ppc")
  admissible_flag <- length(admissibility_reasons) == 0L

  set.seed(stan_seed + 300000L + ii)
  selected_parameter_draws <- as.data.frame(param_draws_mat[sample(seq_len(nrow(param_draws_mat)), size = min(400L, nrow(param_draws_mat))), , drop = FALSE])

  compact_rds <- list(
    dataset = dataset_name,
    model_id = model_id_now,
    branch = "Stage8A",
    risk_scale = "transition_only_main",
    prior_branch = model_row$prior_branch[[1L]],
    site_prior_family = model_row$site_prior_family[[1L]],
    site_placement_label = model_row$site_placement_label[[1L]],
    latency_family = model_row$latency_family[[1L]],
    formula_anchor = model_row$formula_anchor[[1L]],
    fit_status = fit_status,
    admissible_flag = admissible_flag,
    admissibility_reasons = admissibility_reasons,
    coefficient_summary = coef_tbl,
    diagnostics_parameter_level = diagnostic_param_tbl_model,
    ppc_summary = ppc_model_tbl,
    posterior_subject_profile = if (isTRUE(admissible_flag)) subject_profile_tbl else tibble(),
    posterior_subject_yearly = if (isTRUE(admissible_flag)) subject_year_model_tbl else tibble(),
    posterior_cohort_yearly = if (isTRUE(admissible_flag)) cohort_model_tbl else tibble(),
    posterior_classification = if (isTRUE(admissible_flag)) class_model_tbl else tibble(),
    prior_predictive_summary = prior_pred_main$summary,
    selected_parameter_draws = selected_parameter_draws
  )
  if (isTRUE(save_full_stanfit_rds)) compact_rds$fit <- fit
  rds_path <- file.path(export_path, paste0(model_id_now, "__bayes_stage8_fit.rds"))
  save_rds_preserve_schema(compact_rds, rds_path)

  main_registry_rows[[length(main_registry_rows) + 1L]] <- tibble(
    dataset = dataset_name,
    dataset_key = dataset_name,
    model_id = model_id_now,
    legacy_model_id = model_row$legacy_model_id[[1L]],
    branch = "Stage8A",
    risk_scale = "transition_only_main",
    prior_branch = model_row$prior_branch[[1L]],
    site_prior_family = model_row$site_prior_family[[1L]],
    site_placement_label = model_row$site_placement_label[[1L]],
    structural_model_id = model_row$structural_model_id[[1L]],
    latency_family = model_row$latency_family[[1L]],
    family_code = model_row$family_code[[1L]],
    formula_anchor = model_row$formula_anchor[[1L]],
    incidence_site_indicator = model_row$incidence_site_indicator[[1L]],
    latency_site_indicator = model_row$latency_site_indicator[[1L]],
    latency_interaction_indicator = model_row$latency_interaction_indicator[[1L]],
    is_supplementary_branch = model_row$is_supplementary_branch[[1L]],
    fit_status = fit_status,
    fit_error_message = fit_error_message,
    admissible_flag = admissible_flag,
    admissibility_reasons = if (length(admissibility_reasons) == 0L) "" else paste(admissibility_reasons, collapse = "|"),
    prior_degenerate_flag = prior_pred_main$degeneracy$degenerate_flag[[1L]],
    posterior_degenerate_flag = posterior_degeneracy$degenerate_flag[[1L]],
    ppc_gross_contradiction_flag = ppc_gross_contradiction_flag,
    divergences = divergences,
    max_rhat = max_rhat,
    min_bulk_ess = min_bulk_ess,
    min_tail_ess = min_tail_ess,
    treedepth_exceeded = treedepth_exceeded,
    waic = info_criteria$waic,
    looic = info_criteria$looic,
    p_waic = info_criteria$p_waic,
    p_waic_high_n = info_criteria$p_waic_high_n,
    p_waic_high_pct = info_criteria$p_waic_high_pct,
    p_loo = info_criteria$p_loo,
    pareto_k_max = info_criteria$pareto_k_max,
    pareto_k_bad_n = info_criteria$pareto_k_bad_n,
    pareto_k_bad_pct = info_criteria$pareto_k_bad_pct,
    pareto_k_very_bad_n = info_criteria$pareto_k_very_bad_n,
    waic_warning_flag = info_criteria$waic_warning_flag,
    loo_warning_flag = info_criteria$loo_warning_flag,
    info_criteria_warning_detail = info_criteria$info_criteria_warning_detail,
    n = nrow(dataset_df),
    n_event = sum(dataset_df$event_main),
    n_censor_main = sum(dataset_df$censor_main),
    n_remission = sum(dataset_df$remission_flag),
    cohort_mean_cure_fraction_mean = mean(rowMeans(linear_terms$cure_prob_mat)),
    cohort_mean_cure_fraction_q025 = stats::quantile(rowMeans(linear_terms$cure_prob_mat), 0.025, names = FALSE),
    cohort_mean_cure_fraction_q50 = stats::quantile(rowMeans(linear_terms$cure_prob_mat), 0.500, names = FALSE),
    cohort_mean_cure_fraction_q975 = stats::quantile(rowMeans(linear_terms$cure_prob_mat), 0.975, names = FALSE),
    rds_path = rds_path
  )

  diagnostic_rows[[length(diagnostic_rows) + 1L]] <- diagnostic_param_tbl_model
  ppc_rows_all[[length(ppc_rows_all) + 1L]] <- ppc_model_tbl
  if (isTRUE(admissible_flag)) {
    subject_profile_rows[[length(subject_profile_rows) + 1L]] <- subject_profile_tbl
    subject_yearly_rows[[length(subject_yearly_rows) + 1L]] <- subject_year_model_tbl
    cohort_yearly_rows[[length(cohort_yearly_rows) + 1L]] <- cohort_model_tbl
    class_rows[[length(class_rows) + 1L]] <- class_model_tbl
  }

  rm(fit, draws_compact, draws_pred, linear_terms, param_array, param_draws_mat)
  gc(verbose = FALSE)
  emit_stage8_progress(ii, total_models, model_id_now, paste0("completed; elapsed=", format_stage8_number(elapsed_stage8_seconds(model_started_at), 1L), "s; WAIC=", format_stage8_number(info_criteria$waic, 2L), "; LOOIC=", format_stage8_number(info_criteria$looic, 2L), "; Pareto k max=", format_stage8_number(info_criteria$pareto_k_max, 3L)))
}
}, warning = function(w) {
  if (is_localhost_connection_warning(w)) tryInvokeRestart("muffleWarning")
})
# 🔴 Assemble: final Stage-8A tables, prior-sensitivity flags, and plot sources ===============================
## 🟠 Assemble: bind, annotate, and enrich all exported Stage-8A objects ===============================
model_annotation <- model_grid %>%
  select(model_id, legacy_model_id, branch, risk_scale, prior_branch, site_prior_family, site_placement_label, structural_model_id, latency_family, family_code, family_id, formula_anchor, incidence_site_indicator, latency_site_indicator, latency_interaction_indicator, is_supplementary_branch, dataset, dataset_key, reuse_existing_fit)

model_registry <- bind_rows_safe(main_registry_rows) %>%
  left_join_replacing_columns(model_annotation, by = "model_id") %>%
  left_join_replacing_columns(screening_lookup, by = "model_id") %>%
  mutate(
    prior_branch = normalize_prior_branch_value(prior_branch),
    screening_flag = dplyr::coalesce(screening_flag, cure_model_eligibility_flag),
    screening_detail = dplyr::coalesce(screening_detail, screening_note),
    fit_reused_flag = as.logical(reuse_existing_fit),
    branch = ifelse(is.na(branch) | !nzchar(branch), "Stage8A", as.character(branch)),
    risk_scale = ifelse(is.na(risk_scale) | !nzchar(risk_scale), "transition_only_main", as.character(risk_scale)),
    site_prior_family = ifelse(is.na(site_prior_family) | !nzchar(site_prior_family), ifelse(incidence_site_indicator | latency_site_indicator, "normal_0_1_main", "not_applicable"), as.character(site_prior_family))
  ) %>%
  arrange(factor(model_id, levels = model_order))

coefficient_summary <- bind_rows_safe(coef_rows) %>%
  left_join_replacing_columns(model_annotation %>% select(model_id, dataset, dataset_key, branch, risk_scale, prior_branch, site_prior_family, latency_family, formula_anchor), by = "model_id") %>%
  arrange(factor(model_id, levels = model_order), parameter)

posterior_subject_profile <- bind_rows_safe(subject_profile_rows) %>%
  left_join_replacing_columns(model_annotation %>% select(model_id, dataset, dataset_key, branch, risk_scale, prior_branch, site_prior_family, structural_model_id, family_code, latency_family, formula_anchor, site_placement_label), by = "model_id") %>%
  arrange(factor(model_id, levels = model_order), unique_person_id)

posterior_subject_yearly <- bind_rows_safe(subject_yearly_rows) %>%
  left_join_replacing_columns(model_annotation %>% select(model_id, dataset, dataset_key, branch, risk_scale, prior_branch, site_prior_family, structural_model_id, family_code, latency_family, formula_anchor, site_placement_label), by = "model_id") %>%
  left_join_replacing_columns(horizon_metadata_registry %>% select(dataset, horizon_year, support_tier, horizon_evidence_class, claim_restriction_flag, interpretation_note, observed_km_risk), by = c("dataset", "horizon_year")) %>%
  arrange(factor(model_id, levels = model_order), horizon_year, unique_person_id)

posterior_cohort_yearly <- bind_rows_safe(cohort_yearly_rows) %>%
  left_join_replacing_columns(model_annotation %>% select(model_id, dataset, dataset_key, branch, risk_scale, prior_branch, site_prior_family, structural_model_id, family_code, latency_family, formula_anchor, site_placement_label), by = "model_id") %>%
  left_join_replacing_columns(horizon_metadata_registry %>% select(dataset, horizon_year, support_tier, horizon_evidence_class, claim_restriction_flag, interpretation_note, observed_km_risk), by = c("dataset", "horizon_year")) %>%
  arrange(factor(model_id, levels = model_order), horizon_year)

posterior_classification <- bind_rows_safe(class_rows) %>%
  left_join_replacing_columns(model_annotation %>% select(model_id, dataset, dataset_key, branch, risk_scale, prior_branch, site_prior_family, structural_model_id, family_code, latency_family, formula_anchor, site_placement_label), by = "model_id") %>%
  left_join_replacing_columns(horizon_metadata_registry %>% select(dataset, horizon_year, support_tier, horizon_evidence_class, claim_restriction_flag, interpretation_note, observed_km_risk, denom_case, denom_control, classification_estimable_flag), by = c("dataset", "horizon_year")) %>%
  arrange(factor(model_id, levels = model_order), horizon_year, threshold)

diagnostics_parameter_level <- bind_rows_safe(diagnostic_rows) %>%
  left_join_replacing_columns(model_annotation %>% select(model_id, dataset, dataset_key, branch, risk_scale, prior_branch, site_prior_family, latency_family, formula_anchor), by = "model_id") %>%
  arrange(factor(model_id, levels = model_order), parameter)

ppc_summary <- bind_rows_safe(ppc_rows_all) %>%
  left_join_replacing_columns(model_annotation %>% select(model_id, dataset, dataset_key, branch, risk_scale, prior_branch, site_prior_family, structural_model_id, family_code, latency_family, formula_anchor, site_placement_label), by = "model_id") %>%
  left_join_replacing_columns(horizon_metadata_registry %>% select(dataset, horizon_year, support_tier, horizon_evidence_class, claim_restriction_flag, interpretation_note, observed_km_risk), by = c("dataset", "horizon_year")) %>%
  arrange(factor(model_id, levels = model_order), horizon_year)

prior_predictive_summary <- normalize_prior_predictive_summary(bind_rows_safe(prior_predictive_rows)) %>%
  left_join_replacing_columns(model_annotation %>% select(model_id, dataset, dataset_key, branch, risk_scale, prior_branch, latency_family, formula_anchor, site_prior_family), by = c("model_id", "dataset", "branch", "risk_scale", "prior_branch")) %>%
  arrange(factor(model_id, levels = model_order), metric, horizon_year)

cohort_point_metrics <- build_cohort_point_metrics(posterior_subject_yearly, ipcw_registry)
posterior_cohort_yearly <- posterior_cohort_yearly %>%
  left_join_replacing_columns(cohort_point_metrics, by = c("dataset", "model_id", "horizon_year")) %>%
  group_by(model_id) %>%
  arrange(horizon_year, .by_group = TRUE) %>%
  mutate(ibs_from_1y_to_h = cumsum(replace_na(brier_ipcw, 0)) / seq_along(horizon_year)) %>%
  ungroup()

prior_branch_delta <- build_prior_branch_delta(model_registry, posterior_cohort_yearly, posterior_classification) %>%
  left_join(horizon_metadata_registry %>% select(dataset, horizon_year, support_tier, horizon_evidence_class, claim_restriction_flag), by = c("dataset", "horizon_year")) %>%
  arrange(dataset, structural_model_id, family_code, horizon_year, threshold, metric)

cohort_prior_flag <- prior_branch_delta %>%
  filter(is.na(threshold)) %>%
  group_by(dataset, branch, risk_scale, structural_model_id, family_code, formula_anchor, site_prior_family, horizon_year) %>%
  summarise(prior_tail_sensitive = any(prior_tail_sensitive, na.rm = TRUE), .groups = "drop")

class_prior_flag <- prior_branch_delta %>%
  filter(!is.na(threshold)) %>%
  group_by(dataset, branch, risk_scale, structural_model_id, family_code, formula_anchor, site_prior_family, horizon_year, threshold) %>%
  summarise(prior_tail_sensitive_threshold = any(prior_tail_sensitive, na.rm = TRUE), .groups = "drop")

posterior_cohort_yearly <- posterior_cohort_yearly %>%
  left_join(cohort_prior_flag, by = c("dataset", "branch", "risk_scale", "structural_model_id", "family_code", "formula_anchor", "site_prior_family", "horizon_year")) %>%
  mutate(
    prior_tail_sensitive = dplyr::coalesce(prior_tail_sensitive, FALSE),
    claim_restriction_flag = ifelse(prior_tail_sensitive & horizon_evidence_class == "mostly_extrapolated", "projection_plus_prior_sensitive", claim_restriction_flag)
  ) %>%
  arrange(factor(model_id, levels = model_order), horizon_year)

posterior_classification <- posterior_classification %>%
  left_join(cohort_prior_flag, by = c("dataset", "branch", "risk_scale", "structural_model_id", "family_code", "formula_anchor", "site_prior_family", "horizon_year")) %>%
  left_join(class_prior_flag, by = c("dataset", "branch", "risk_scale", "structural_model_id", "family_code", "formula_anchor", "site_prior_family", "horizon_year", "threshold")) %>%
  mutate(
    prior_tail_sensitive = dplyr::coalesce(prior_tail_sensitive_threshold, FALSE) | dplyr::coalesce(prior_tail_sensitive, FALSE),
    claim_restriction_flag = ifelse(prior_tail_sensitive & horizon_evidence_class == "mostly_extrapolated", "projection_plus_prior_sensitive", claim_restriction_flag)
  ) %>%
  select(-any_of(c("prior_tail_sensitive_threshold"))) %>%
  arrange(factor(model_id, levels = model_order), horizon_year, threshold)

ppc_summary <- ppc_summary %>%
  left_join(cohort_prior_flag, by = c("dataset", "branch", "risk_scale", "structural_model_id", "family_code", "formula_anchor", "site_prior_family", "horizon_year")) %>%
  mutate(
    prior_tail_sensitive = dplyr::coalesce(prior_tail_sensitive, FALSE),
    claim_restriction_flag = ifelse(prior_tail_sensitive & horizon_evidence_class == "mostly_extrapolated", "projection_plus_prior_sensitive", claim_restriction_flag)
  ) %>%
  arrange(factor(model_id, levels = model_order), horizon_year)

posterior_subject_yearly <- posterior_subject_yearly %>%
  left_join(cohort_prior_flag, by = c("dataset", "branch", "risk_scale", "structural_model_id", "family_code", "formula_anchor", "site_prior_family", "horizon_year")) %>%
  mutate(
    prior_tail_sensitive = dplyr::coalesce(prior_tail_sensitive, FALSE),
    claim_restriction_flag = ifelse(prior_tail_sensitive & horizon_evidence_class == "mostly_extrapolated", "projection_plus_prior_sensitive", claim_restriction_flag)
  ) %>%
  arrange(factor(model_id, levels = model_order), horizon_year, unique_person_id)

model_prior_flag <- bind_rows(
  posterior_cohort_yearly %>% select(model_id, prior_tail_sensitive),
  posterior_classification %>% select(model_id, prior_tail_sensitive)
) %>% group_by(model_id) %>% summarise(prior_tail_sensitive_any = any(prior_tail_sensitive, na.rm = TRUE), .groups = "drop")

model_registry <- model_registry %>%
  left_join(model_prior_flag, by = "model_id") %>%
  mutate(prior_tail_sensitive_any = dplyr::coalesce(prior_tail_sensitive_any, FALSE)) %>%
  arrange(factor(model_id, levels = model_order))

posterior_delta_vs_nocure <- tryCatch(
  rebuild_delta_vs_nocure(posterior_cohort_yearly, posterior_classification, nocure_cohort_long, nocure_class_long),
  error = function(e) {
    warning(paste0("Stage 8A could not rebuild posterior_delta_vs_nocure; an empty schema-preserving table will be written instead. Reason: ", conditionMessage(e)), call. = FALSE)
    empty_delta_template()
  }
)

posterior_delta_vs_nocure <- posterior_delta_vs_nocure %>%
  left_join_replacing_columns(model_annotation %>% select(model_id, dataset, dataset_key, branch, risk_scale, prior_branch, site_prior_family, structural_model_id, family_code, latency_family, formula_anchor, site_placement_label), by = c("model_id", "dataset", "branch", "risk_scale", "prior_branch", "site_prior_family", "latency_family", "formula_anchor")) %>%
  left_join_replacing_columns(model_registry %>% select(model_id, admissible_flag, admissibility_reasons, screening_flag, screening_detail), by = "model_id") %>%
  mutate(threshold_key = make_threshold_key(threshold)) %>%
  left_join_replacing_columns(
    bind_rows(
      posterior_cohort_yearly %>% transmute(dataset, model_id, formula_anchor, horizon_year, threshold_key = make_threshold_key(NA_real_), support_tier, horizon_evidence_class, claim_restriction_flag, prior_tail_sensitive, classification_estimable_flag = NA, denom_case = NA_real_, denom_control = NA_real_),
      posterior_classification %>% transmute(dataset, model_id, formula_anchor, horizon_year, threshold_key = make_threshold_key(threshold), support_tier, horizon_evidence_class, claim_restriction_flag, prior_tail_sensitive, classification_estimable_flag, denom_case, denom_control)
    ) %>% distinct(),
    by = c("dataset", "model_id", "formula_anchor", "horizon_year", "threshold_key")
  ) %>%
  select(-threshold_key) %>%
  arrange(dataset, factor(model_id, levels = model_order), horizon_year, threshold, metric, no_cure_model_id)

incidence_anchor_update <- build_incidence_anchor_update(model_registry, coefficient_summary, export_path)
if (ncol(incidence_anchor_update) == 0L) incidence_anchor_update <- empty_incidence_anchor_update()
incidence_anchor_update <- incidence_anchor_update %>% arrange(factor(model_id, levels = model_order), age_sex_anchor_cell)

hazard_plausibility <- build_hazard_plausibility(model_registry, posterior_cohort_yearly)
if (ncol(hazard_plausibility) == 0L) hazard_plausibility <- empty_hazard_plausibility()
hazard_plausibility <- hazard_plausibility %>% arrange(factor(model_id, levels = model_order))

uncured_decomposition <- build_uncured_decomposition(posterior_cohort_yearly)
if (ncol(uncured_decomposition) == 0L) uncured_decomposition <- empty_uncured_decomposition()
uncured_decomposition <- uncured_decomposition %>% arrange(factor(model_id, levels = model_order), horizon_year)

stage8A_vs_stage8B_delta <- empty_stage8A_vs_stage8B_delta()

carryforward_metadata <- model_registry %>%
  select(
    dataset, dataset_key, model_id, branch, risk_scale, prior_branch, site_prior_family, site_placement_label,
    latency_family, formula_anchor, incidence_site_indicator, latency_site_indicator, latency_interaction_indicator,
    admissible_flag, admissibility_reasons, fit_reused_flag, prior_tail_sensitive_any,
    cure_model_eligibility_flag, primary_gate_method, primary_gate_flag, receus_primary_class,
    presence_modifier_flag, cure_presence_support_flag, followup_contradiction_flag, followup_not_contradicted_flag,
    screening_note, screening_flag, screening_detail, carry_forward_stage8
  )

stage8_metadata_registry <- bind_rows(
  tibble(metadata_group = "stage", metadata_name = c("stage_name", "stage_role", "branch", "risk_scale"), metadata_value = c("Stage 8A", "Bayesian transition-only cure branch aligned to revised v5 specification.", "Stage8A", "transition_only_main")),
  tibble(metadata_group = "paths", metadata_name = c("stage1_export_path", "stage1_bundle_file", "stage1_datasets_file", "stage5_export_path", "stage6_export_path", "stage8b_export_path", "export_path"), metadata_value = c(normalize_existing_path(stage1_export_path), normalize_existing_path(stage1_bundle_file), normalize_existing_path(stage1_datasets_file), normalize_existing_path(stage5_export_path), normalize_existing_path(stage6_export_path), normalize_existing_path(stage8b_export_path), normalize_existing_path(export_path))),
  tibble(metadata_group = "reuse", metadata_name = c("reuse_existing_stage8_outputs", "require_existing_rds_to_skip_fit", "preserve_existing_diagnostic_pdf"), metadata_value = c(as.character(reuse_existing_stage8_outputs), as.character(require_existing_rds_to_skip_fit), as.character(preserve_existing_diagnostic_pdf))),
  tibble(metadata_group = "stan", metadata_name = c("stan_chains", "stan_iter", "stan_warmup", "stan_thin", "stan_seed", "stan_adapt_delta", "stan_max_treedepth", "posterior_prediction_draws"), metadata_value = c(stan_chains, stan_iter, stan_warmup, stan_thin, stan_seed, stan_adapt_delta, stan_max_treedepth, posterior_prediction_draws)),
  tibble(metadata_group = "prior", metadata_name = c("anchor_prior_branch_label", "neutral_prior_branch_label", "site_prior_main_label", "site_prior_sensitivity_label", "anchor_alpha_gp_female_under20", "anchor_materiality_risk", "anchor_materiality_cure_fraction", "anchor_materiality_false_positive_burden", "anchor_materiality_FP100", "anchor_materiality_NB", "anchor_materiality_PPV", "anchor_materiality_TPR"), metadata_value = c("anchor_informed", "neutral_no_external_info", "normal_0_1_main", "student_t3_0_1_sensitivity", -9.581369553169, prior_materiality_thresholds$risk, prior_materiality_thresholds$cure_fraction, prior_materiality_thresholds$false_positive_burden, prior_materiality_thresholds$FP100, prior_materiality_thresholds$NB, prior_materiality_thresholds$PPV, prior_materiality_thresholds$TPR)),
  tibble(metadata_group = "grid", metadata_name = c("include_merged_incidence_site_supplementary", "common_horizons_year", "risk_thresholds", "n_models_total", "n_models_reused", "n_models_to_fit"), metadata_value = c(as.character(include_merged_incidence_site_supplementary), paste(horizons_year, collapse = ","), paste(format(risk_thresholds, trim = TRUE, scientific = FALSE), collapse = ","), nrow(model_grid), n_models_reused, n_models_to_fit))
)

output_audit <- make_stage8_output_audit(model_grid, model_registry, posterior_cohort_yearly, posterior_classification, ppc_summary, posterior_delta_vs_nocure, diagnostic_pdf_path)

# 🔴 Create: figure objects from exported data frames ===============================
## 🟠 Create: combined PDF pages and per-plot PNG graphics ===============================
plot_objects <- list()

if (nrow(posterior_cohort_yearly) > 0L) {
  g_risk <- posterior_cohort_yearly %>%
    ggplot(aes(x = horizon_year, y = meanRisk_Bayes_mean, color = model_id, fill = model_id)) +
    geom_ribbon(aes(ymin = meanRisk_Bayes_q025, ymax = meanRisk_Bayes_q975), alpha = 0.15, linewidth = 0) +
    geom_line(linewidth = 0.7) +
    facet_grid(dataset ~ prior_branch, scales = "free_y") +
    labs(title = "Stage 8A posterior cohort mean risk trajectories", x = "Horizon (years)", y = "Posterior mean risk") +
    theme_bw() + theme(legend.position = "none")
  plot_objects[["bayes_stage8_plot01_risk_trajectories"]] <- g_risk

  g_hazard <- posterior_cohort_yearly %>%
    ggplot(aes(x = horizon_year, y = meanHazard_mean, color = model_id, fill = model_id)) +
    geom_ribbon(aes(ymin = meanHazard_q025, ymax = meanHazard_q975), alpha = 0.15, linewidth = 0) +
    geom_line(linewidth = 0.7) +
    facet_grid(dataset ~ prior_branch, scales = "free_y") +
    labs(title = "Stage 8A posterior cohort mean hazard trajectories", x = "Horizon (years)", y = "Posterior mean hazard") +
    theme_bw() + theme(legend.position = "none")
  plot_objects[["bayes_stage8_plot02_hazard_trajectories"]] <- g_hazard
}

if (nrow(posterior_classification) > 0L) {
  nb_plot_df <- posterior_classification %>% filter(horizon_year %in% c(1, 2, 5))
  if (nrow(nb_plot_df) > 0L) {
    g_nb <- nb_plot_df %>%
      ggplot(aes(x = threshold, y = NB_mean, color = model_id, fill = model_id)) +
      geom_ribbon(aes(ymin = NB_q025, ymax = NB_q975), alpha = 0.15, linewidth = 0) +
      geom_line(linewidth = 0.7) +
      facet_grid(dataset + prior_branch ~ horizon_year, scales = "free_y") +
      labs(title = "Stage 8A net benefit by threshold", x = "Risk threshold", y = "Net benefit") +
      theme_bw() + theme(legend.position = "none")
    plot_objects[["bayes_stage8_plot03_net_benefit"]] <- g_nb
  }
}

if (nrow(ppc_summary) > 0L) {
  g_ppc <- ppc_summary %>%
    ggplot(aes(x = horizon_year, y = posterior_mean_risk, color = model_id)) +
    geom_errorbar(aes(ymin = posterior_q025_risk, ymax = posterior_q975_risk), width = 0.12, alpha = 0.6) +
    geom_line(linewidth = 0.6) +
    geom_point(linewidth = 0.6) +
    geom_point(aes(y = observed_km_risk), shape = 4, size = 2.0, stroke = 0.9, color = "black") +
    facet_grid(dataset ~ prior_branch, scales = "free_y") +
    labs(title = "Stage 8A posterior predictive checks against observed KM risk", x = "Horizon (years)", y = "Risk") +
    theme_bw() + theme(legend.position = "none")
  plot_objects[["bayes_stage8_plot04_ppc"]] <- g_ppc
}

if (nrow(prior_branch_delta) > 0L) {
  g_anchor_delta <- prior_branch_delta %>%
    filter(metric %in% c("risk", "cure_fraction"), is.na(threshold)) %>%
    ggplot(aes(x = horizon_year, y = delta_mean_anchor_minus_neutral, color = metric, fill = metric)) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_ribbon(aes(ymin = delta_q025_anchor_minus_neutral, ymax = delta_q975_anchor_minus_neutral), alpha = 0.15, linewidth = 0) +
    geom_line(linewidth = 0.7) +
    facet_grid(dataset ~ metric, scales = "free_y") +
    labs(title = "Anchor-informed minus neutral prior branch deltas", x = "Horizon (years)", y = "Delta") +
    theme_bw() + theme(legend.position = "none")
  plot_objects[["bayes_stage8_plot05_anchor_vs_neutral_delta"]] <- g_anchor_delta
}

if (nrow(incidence_anchor_update) > 0L) {
  g_anchor_update <- incidence_anchor_update %>%
    ggplot(aes(x = age_sex_anchor_cell, y = posterior_mean_one_year_risk, ymin = posterior_lower_one_year_risk, ymax = posterior_upper_one_year_risk)) +
    geom_pointrange() +
    facet_grid(dataset ~ prior_branch, scales = "free_y") +
    labs(title = "Stage 8A incidence anchor update", x = "Age-sex anchor cell", y = "Posterior one-year susceptibility risk") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot_objects[["bayes_stage8_plot06_incidence_anchor_update"]] <- g_anchor_update
}

if (nrow(uncured_decomposition) > 0L) {
  g_uncured <- uncured_decomposition %>%
    ggplot(aes(x = horizon_year, y = uncured_survival_mean, color = model_id, fill = model_id)) +
    geom_ribbon(aes(ymin = uncured_survival_q025, ymax = uncured_survival_q975), alpha = 0.15, linewidth = 0) +
    geom_line(linewidth = 0.7) +
    facet_grid(dataset ~ prior_branch, scales = "free_y") +
    labs(title = "Stage 8A uncured-only survival trajectories", x = "Horizon (years)", y = "Uncured-only survival") +
    theme_bw() + theme(legend.position = "none")
  plot_objects[["bayes_stage8_plot07_uncured_decomposition"]] <- g_uncured
}
# 🔴 Export: CSV tables, plot files, and final manifest ===============================
## 🟠 Export: write schema-preserving Stage-8A artifacts ===============================
model_registry <- simplify_scalar_list_cols(model_registry)
coefficient_summary <- simplify_scalar_list_cols(coefficient_summary)
posterior_subject_profile <- simplify_scalar_list_cols(posterior_subject_profile)
posterior_subject_yearly <- simplify_scalar_list_cols(posterior_subject_yearly)
posterior_cohort_yearly <- simplify_scalar_list_cols(posterior_cohort_yearly)
posterior_classification <- simplify_scalar_list_cols(posterior_classification)
posterior_delta_vs_nocure <- simplify_scalar_list_cols(posterior_delta_vs_nocure)
diagnostics_parameter_level <- simplify_scalar_list_cols(diagnostics_parameter_level)
ppc_summary <- simplify_scalar_list_cols(ppc_summary)
prior_predictive_summary <- simplify_scalar_list_cols(prior_predictive_summary)
prior_branch_delta <- simplify_scalar_list_cols(prior_branch_delta)
incidence_anchor_update <- simplify_scalar_list_cols(incidence_anchor_update)
hazard_plausibility <- simplify_scalar_list_cols(hazard_plausibility)
uncured_decomposition <- simplify_scalar_list_cols(uncured_decomposition)
stage8A_vs_stage8B_delta <- simplify_scalar_list_cols(stage8A_vs_stage8B_delta)
stage8_metadata_registry <- simplify_scalar_list_cols(stage8_metadata_registry)
output_audit <- simplify_scalar_list_cols(output_audit)
carryforward_metadata <- simplify_scalar_list_cols(carryforward_metadata)

write_csv_preserve_schema(model_registry, file.path(export_path, "bayes_stage8_model_registry.csv"))
write_csv_preserve_schema(coefficient_summary, file.path(export_path, "bayes_stage8_coefficient_summary.csv"))
write_csv_preserve_schema(posterior_subject_profile, file.path(export_path, "bayes_stage8_posterior_subject_profile.csv.gz"))
write_csv_preserve_schema(posterior_subject_yearly, file.path(export_path, "bayes_stage8_posterior_subject_yearly.csv.gz"))
write_csv_preserve_schema(posterior_cohort_yearly, file.path(export_path, "bayes_stage8_posterior_cohort_yearly.csv"))
write_csv_preserve_schema(posterior_classification, file.path(export_path, "bayes_stage8_posterior_classification.csv"))
write_csv_preserve_schema(posterior_delta_vs_nocure, file.path(export_path, "bayes_stage8_posterior_delta_vs_nocure.csv"))
write_csv_preserve_schema(diagnostics_parameter_level, file.path(export_path, "bayes_stage8_diagnostics_parameter_level.csv"))
write_csv_preserve_schema(ppc_summary, file.path(export_path, "bayes_stage8_ppc_summary.csv"))
write_csv_preserve_schema(prior_predictive_summary, file.path(export_path, "bayes_stage8_prior_predictive_summary.csv"))
write_csv_preserve_schema(carryforward_metadata, file.path(export_path, "bayes_stage8_carryforward_metadata.csv"))
write_csv_preserve_schema(output_audit, file.path(export_path, "bayes_stage8_output_audit.csv"))
write_csv_preserve_schema(prior_branch_delta, file.path(export_path, "bayes_stage8_prior_branch_delta.csv"))
write_csv_preserve_schema(incidence_anchor_update, file.path(export_path, "bayes_stage8_incidence_anchor_update.csv"))
write_csv_preserve_schema(hazard_plausibility, file.path(export_path, "bayes_stage8_hazard_plausibility.csv"))
write_csv_preserve_schema(uncured_decomposition, file.path(export_path, "bayes_stage8_uncured_decomposition.csv"))
write_csv_preserve_schema(stage8_metadata_registry, file.path(export_path, "bayes_stage8_metadata_registry.csv"))
write_csv_preserve_schema(stage8A_vs_stage8B_delta, file.path(export_path, "bayes_stage8_stage8A_vs_stage8B_delta.csv"))

if (isTRUE(always_regenerate_plot_bundle) || n_models_to_fit > 0L || !pdf_file_is_usable(diagnostic_pdf_path) || !isTRUE(preserve_existing_diagnostic_pdf)) {
  pdf_ok <- tryCatch({
    safe_generate_diagnostic_pdf(trace_records = trace_records, plot_objects = plot_objects, final_path = diagnostic_pdf_path)
    TRUE
  }, error = function(e) {
    warning(paste0("Diagnostic PDF was not regenerated. Existing PDF was preserved if present. Reason: ", conditionMessage(e)), call. = FALSE)
    FALSE
  })
  if (isTRUE(pdf_ok)) {
    message("Stage 8A diagnostic PDF regenerated safely: ", diagnostic_pdf_path)
  }
}

if (length(plot_objects) > 0L) {
  for (nm in names(plot_objects)) {
    png_path <- file.path(export_path, paste0(nm, ".png"))
    try(save_plot_object(plot_objects[[nm]], png_path), silent = FALSE)
  }
}

export_file_names <- c(
  "bayes_stage8_model_registry.csv",
  "bayes_stage8_coefficient_summary.csv",
  "bayes_stage8_posterior_subject_profile.csv.gz",
  "bayes_stage8_posterior_subject_yearly.csv.gz",
  "bayes_stage8_posterior_cohort_yearly.csv",
  "bayes_stage8_posterior_classification.csv",
  "bayes_stage8_posterior_delta_vs_nocure.csv",
  "bayes_stage8_diagnostics_parameter_level.csv",
  "bayes_stage8_ppc_summary.csv",
  "bayes_stage8_prior_predictive_summary.csv",
  "bayes_stage8_carryforward_metadata.csv",
  "bayes_stage8_output_audit.csv",
  "bayes_stage8_prior_branch_delta.csv",
  "bayes_stage8_incidence_anchor_update.csv",
  "bayes_stage8_hazard_plausibility.csv",
  "bayes_stage8_uncured_decomposition.csv",
  "bayes_stage8_metadata_registry.csv",
  "bayes_stage8_stage8A_vs_stage8B_delta.csv",
  "bayes_stage8_export_manifest.csv",
  "bayes_stage8_diagnostic_plots.pdf"
)

export_object_names <- c(
  "model_registry",
  "coefficient_summary",
  "posterior_subject_profile",
  "posterior_subject_yearly",
  "posterior_cohort_yearly",
  "posterior_classification",
  "posterior_delta_vs_nocure",
  "diagnostics_parameter_level",
  "ppc_summary",
  "prior_predictive_summary",
  "carryforward_metadata",
  "output_audit",
  "prior_branch_delta",
  "incidence_anchor_update",
  "hazard_plausibility",
  "uncured_decomposition",
  "stage8_metadata_registry",
  "stage8A_vs_stage8B_delta",
  "export_manifest",
  "diagnostic_plot_pdf"
)

export_descriptions <- c(
  "Stage 8A model-level diagnostics, admissibility, prior-branch, and site-prior registry",
  "Stage 8A posterior coefficient summaries",
  "Subject-level posterior cure/susceptible profile summaries",
  "Subject-year posterior prediction table on the common 1-10 year grid",
  "Horizon-level Stage 8A cohort risk/survival/hazard summaries with calibration/Brier/AUC fields",
  "Threshold-level Stage 8A classification and clinical-usefulness summaries",
  "Bayesian cure versus Stage 5 non-cure comparison table",
  "Parameter-level convergence diagnostics",
  "Posterior predictive check summary against observed KM risk",
  "Prior-predictive summary table for the actual fitted prior branch",
  "Canonical Stage 6 carry-forward fields attached to each retained Stage 8A fit",
  "Stage 8A export audit and structural QC summary",
  "Anchor-informed versus neutral prior-branch delta table",
  "Prior-to-posterior incidence anchor update source table",
  "Hazard-shape plausibility table on the 1-10 year annual grid",
  "Uncured-only / non-cure-fraction supporting decomposition table",
  "Implementation metadata, thresholds, paths, and prior labels used in Stage 8A",
  "Schema-preserving Stage 8A-versus-Stage 8B delta placeholder for later integration",
  "Manifest of Stage 8A exported non-PNG files",
  "Combined diagnostic and figure PDF; each plot is also saved as a standalone PNG"
)

export_manifest <- tibble(
  file_name = export_file_names,
  object_name = export_object_names,
  description = export_descriptions,
  file_path = file.path(export_path, export_file_names)
)
write_csv_preserve_schema(export_manifest, file.path(export_path, "bayes_stage8_export_manifest.csv"))

message("Stage 8A completed. Reused models: ", n_models_reused, "; newly fitted models: ", n_models_to_fit, ".")
