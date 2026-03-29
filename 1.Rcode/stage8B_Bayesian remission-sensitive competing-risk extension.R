# 🔴 Configure: Stage8B paths, controls, and fixed grids ===============================
# 🟧 Configure: OS-aware root path ===============================
sys_name <- Sys.info()[["sysname"]]

project_root <- switch(
  sys_name,
  "Darwin"  = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  "Windows" = "C:/Users/clair/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  stop("Unsupported OS: ", sys_name)
)

# 🟧 Configure: Stage8B paths ===============================
stage1_export_path  <- file.path(project_root, "stage1_Backbone lock")
export_path         <- file.path(project_root, "stage8B_Bayesian remission-sensitive competing-risk extension")
stage6_export_path  <- file.path(project_root, "stage6_Cure-appropriateness screening")
stage8a_export_path <- file.path(project_root, "stage8A_Bayesian transition-only cure")




stage1_bundle_file <- file.path(stage1_export_path, "stage1_backbone_bundle.rds")
stage1_analysis_datasets_file <- file.path(stage1_export_path, "stage1_analysis_datasets.rds")
stage1_horizon_registry_file <- file.path(stage1_export_path, "stage1_horizon_registry.csv")
stage1_threshold_registry_file <- file.path(stage1_export_path, "stage1_threshold_registry.csv")
stage1_scaling_registry_file <- file.path(stage1_export_path, "stage1_scaling_registry.csv")
stage1_dataset_registry_file <- file.path(stage1_export_path, "stage1_dataset_registry.csv")
stage1_export_manifest_file <- file.path(stage1_export_path, "stage1_export_manifest.csv")

stage6_screening_flag_csv <- file.path(stage6_export_path, "stage6_carry_forward_flag_table.csv")

stage8a_model_registry_csv <- file.path(stage8a_export_path, "bayes_stage8_model_registry.csv")
stage8a_cohort_yearly_csv <- file.path(stage8a_export_path, "bayes_stage8_posterior_cohort_yearly.csv")
stage8a_classification_csv <- file.path(stage8a_export_path, "bayes_stage8_posterior_classification.csv")

pnu_site_label <- "PNU"
snu_site_label <- "SNU"

run_model_ids <- NULL
include_merged_incidence_site_supplementary <- TRUE
fit_prior_branches <- c("anchor_informed", "neutral_no_external_info")
fit_site_prior_families <- c("normal_0_1_main", "student_t3_0_1_sensitivity")

reuse_existing_stage8b_rds <- TRUE
save_full_stanfit_rds <- FALSE

stan_chains <- 4L
stan_iter <- 2000L
stan_warmup <- 1000L
stan_thin <- 1L
stan_seed <- 20260328L
stan_adapt_delta <- 0.95
stan_max_treedepth <- 12L
stan_refresh <- 0L

prior_predictive_draws <- 250L
posterior_prediction_draws <- 300L
remission_cut_years <- c(0, 1, 2, 5, 10, 100)
integration_step_year <- 0.05

ppc_tolerance_abs <- 0.12
ess_min_threshold <- 400
rhat_max_threshold <- 1.01
degenerate_draw_fraction_threshold <- 0.90
degenerate_subject_fraction_threshold <- 0.95
tiny_susceptible_prob <- 0.01
huge_susceptible_prob <- 0.99
tiny_median_years <- 0.05
huge_median_years <- 50

prior_predictive_horizons <- c(1L, 2L, 5L, 10L)
prior_tail_warning_mean_years <- 100
prior_tail_warning_q975_years <- 200

prior_materiality_risk <- 0.03
prior_materiality_false_positive_burden <- 0.03
prior_materiality_nb <- 0.01

horizons_year <- NULL
risk_thresholds <- NULL

# 🔴 Initialize: packages, options, and safe defaults ===============================
core_packages <- c(
  "readr", "dplyr", "tibble", "tidyr", "ggplot2", "survival", "matrixStats", "purrr"
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
  library(purrr)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

# 🔴 Define: backbone-safe I/O and scalar helpers ===============================
## 🟠 Define: path, coercion, and join utilities ===============================
assert_existing_file <- function(path, label = "File") {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
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
  stop(sprintf("Unsupported file extension for `%s`.", path), call. = FALSE)
}

make_temp_output_path <- function(path, tag = "tmp") {
  dir <- dirname(path)
  base <- basename(path)
  stamp <- paste0(
    format(Sys.time(), "%Y%m%d%H%M%S"),
    "_",
    sprintf("%08d", sample.int(99999999, 1))
  )
  if (grepl("\\.csv\\.gz$", base, ignore.case = TRUE)) {
    base <- sub("\\.csv\\.gz$", paste0("_", tag, "_", stamp, ".csv.gz"), base, ignore.case = TRUE)
  } else if (grepl("\\.csv$", base, ignore.case = TRUE)) {
    base <- sub("\\.csv$", paste0("_", tag, "_", stamp, ".csv"), base, ignore.case = TRUE)
  } else if (grepl("\\.pdf$", base, ignore.case = TRUE)) {
    base <- sub("\\.pdf$", paste0("_", tag, "_", stamp, ".pdf"), base, ignore.case = TRUE)
  } else if (grepl("\\.rds$", base, ignore.case = TRUE)) {
    base <- sub("\\.rds$", paste0("_", tag, "_", stamp, ".rds"), base, ignore.case = TRUE)
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

pdf_file_is_usable <- function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  info <- file.info(path)
  isTRUE(!is.na(info$size[[1]]) && info$size[[1]] > 0L)
}

nrow_or_zero <- function(df) {
  if (is.null(df) || !inherits(df, "data.frame")) {
    return(0L)
  }
  nrow(df)
}

bind_rows_safe <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) == 0L) {
    return(tibble())
  }
  bind_rows(x)
}

drop_existing_columns <- function(df, cols) {
  if (is.null(df) || !inherits(df, "data.frame") || length(cols) == 0L) {
    return(df)
  }
  dplyr::select(df, -any_of(cols))
}

left_join_replacing_columns <- function(x, y, by) {
  if (is.null(x) || !inherits(x, "data.frame")) {
    return(x)
  }
  if (is.null(y) || !inherits(y, "data.frame") || nrow(y) == 0L) {
    return(x)
  }
  by_cols <- unname(by)
  if (any(!by_cols %in% names(x)) || any(!by_cols %in% names(y))) {
    return(x)
  }
  join_cols <- setdiff(names(y), by_cols)
  x %>%
    drop_existing_columns(join_cols) %>%
    left_join(y, by = by)
}

first_existing_name <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) {
    return(NULL)
  }
  hit[[1L]]
}

format_number <- function(x, digits = 3L) {
  if (length(x) == 0L || is.na(x) || !is.finite(x)) {
    return("NA")
  }
  formatC(as.numeric(x), format = "f", digits = digits)
}

elapsed_seconds <- function(start_time) {
  as.numeric(difftime(Sys.time(), start_time, units = "secs"))
}

emit_progress <- function(done, total, model_id, detail) {
  pct <- if (total <= 0L) 100 else 100 * done / total
  message(sprintf("[%d/%d | %5.1f%%] %s %s", as.integer(done), as.integer(total), pct, model_id, detail))
}

is_localhost_connection_warning <- function(w) {
  grepl("closing unused connection .*<-localhost:", conditionMessage(w))
}

# 🔴 Define: Stage1 import and governance helpers ===============================
## 🟠 Define: Stage1 bundle readers and horizon governance ===============================
empty_stage1_bundle <- function() {
  list(
    datasets = NULL,
    registries = list(
      dataset_registry = tibble(),
      horizon_registry = tibble(),
      threshold_registry = tibble(),
      scaling_registry = tibble(),
      metadata_registry = tibble()
    ),
    config = list()
  )
}

load_stage1_backbone <- function(
    bundle_file,
    datasets_file,
    dataset_registry_file,
    horizon_registry_file,
    threshold_registry_file,
    scaling_registry_file
) {
  out <- empty_stage1_bundle()
  
  if (file.exists(bundle_file)) {
    obj <- readRDS(bundle_file)
    if (is.list(obj)) {
      out$config <- obj$config %||% list()
      out$datasets <- obj$datasets %||% NULL
      regs <- obj$registries %||% list()
      out$registries$dataset_registry <- tibble::as_tibble(regs$dataset_registry %||% tibble())
      out$registries$horizon_registry <- tibble::as_tibble(regs$horizon_registry %||% tibble())
      out$registries$threshold_registry <- tibble::as_tibble(regs$threshold_registry %||% tibble())
      out$registries$scaling_registry <- tibble::as_tibble(regs$scaling_registry %||% tibble())
      out$registries$metadata_registry <- tibble::as_tibble(regs$metadata_registry %||% tibble())
    }
  }
  
  if (is.null(out$datasets) && file.exists(datasets_file)) {
    out$datasets <- readRDS(datasets_file)
  }
  
  if (nrow_or_zero(out$registries$dataset_registry) == 0L && file.exists(dataset_registry_file)) {
    out$registries$dataset_registry <- read_delimited_or_rds(dataset_registry_file)
  }
  if (nrow_or_zero(out$registries$horizon_registry) == 0L && file.exists(horizon_registry_file)) {
    out$registries$horizon_registry <- read_delimited_or_rds(horizon_registry_file)
  }
  if (nrow_or_zero(out$registries$threshold_registry) == 0L && file.exists(threshold_registry_file)) {
    out$registries$threshold_registry <- read_delimited_or_rds(threshold_registry_file)
  }
  if (nrow_or_zero(out$registries$scaling_registry) == 0L && file.exists(scaling_registry_file)) {
    out$registries$scaling_registry <- read_delimited_or_rds(scaling_registry_file)
  }
  
  if (is.null(out$datasets) || !is.list(out$datasets) || length(out$datasets) == 0L) {
    stop("Stage 1 datasets could not be loaded from `stage1_backbone_bundle.rds` or `stage1_analysis_datasets.rds`.", call. = FALSE)
  }
  
  req_datasets <- c("PNU", "SNU", "merged")
  if (any(!req_datasets %in% names(out$datasets))) {
    stop("Stage 1 datasets must contain `PNU`, `SNU`, and `merged`.", call. = FALSE)
  }
  
  out$registries$horizon_registry <- tibble::as_tibble(out$registries$horizon_registry)
  out$registries$threshold_registry <- tibble::as_tibble(out$registries$threshold_registry)
  out$registries$dataset_registry <- tibble::as_tibble(out$registries$dataset_registry)
  out$registries$scaling_registry <- tibble::as_tibble(out$registries$scaling_registry)
  
  out
}

build_horizon_annotation_from_stage1 <- function(horizon_registry, datasets, horizons) {
  out <- tibble::as_tibble(horizon_registry)
  
  if (nrow(out) == 0L) {
    out <- tidyr::crossing(
      dataset = names(datasets),
      horizon_year = as.integer(horizons)
    ) %>%
      mutate(
        dataset_key = dataset,
        support_tier = case_when(
          dataset == "PNU" & horizon_year == 1L ~ "primary_supported",
          dataset == "PNU" & horizon_year == 2L ~ "sensitivity",
          dataset == "PNU" & horizon_year >= 3L ~ "projection",
          dataset %in% c("SNU", "merged") & horizon_year %in% c(1L, 2L) ~ "primary_supported",
          dataset %in% c("SNU", "merged") & horizon_year <= 5L ~ "secondary",
          TRUE ~ "projection"
        ),
        horizon_evidence_class = case_when(
          dataset == "PNU" & horizon_year == 1L ~ "directly_observed_data_supported",
          dataset == "PNU" & horizon_year == 2L ~ "partly_model_dependent",
          dataset == "PNU" & horizon_year >= 3L ~ "mostly_extrapolated",
          dataset %in% c("SNU", "merged") & horizon_year %in% c(1L, 2L) ~ "directly_observed_data_supported",
          dataset %in% c("SNU", "merged") & horizon_year <= 5L ~ "partly_model_dependent",
          TRUE ~ "mostly_extrapolated"
        ),
        claim_restriction_flag = case_when(
          horizon_evidence_class == "directly_observed_data_supported" ~ "primary_claim_allowed",
          horizon_evidence_class == "partly_model_dependent" ~ "secondary_or_sensitivity_only",
          TRUE ~ "projection_only"
        ),
        interpretation_note = case_when(
          claim_restriction_flag == "primary_claim_allowed" ~ "Primary supported horizon with comparatively direct follow-up support.",
          claim_restriction_flag == "secondary_or_sensitivity_only" ~ "Secondary or sensitivity horizon with explicit model dependence.",
          TRUE ~ "Projection-dominant horizon."
        )
      )
  }
  
  out %>%
    transmute(
      dataset = as.character(dataset %||% dataset_key),
      dataset_key = as.character(dataset_key %||% dataset),
      horizon_year = as.integer(safe_numeric(horizon_year)),
      support_tier = as.character(support_tier),
      horizon_evidence_class = as.character(horizon_evidence_class),
      claim_restriction_flag = as.character(claim_restriction_flag),
      interpretation_note = as.character(
        coalesce(
          .data$interpretation_note,
          if ("support_tier_standard" %in% names(out)) as.character(out$support_tier_standard) else NA_character_
        )
      )
    ) %>%
    filter(dataset %in% c("PNU", "SNU", "merged"), horizon_year %in% horizons) %>%
    distinct(dataset, horizon_year, .keep_all = TRUE) %>%
    arrange(factor(dataset, levels = c("PNU", "SNU", "merged")), horizon_year)
}

build_threshold_vector_from_stage1 <- function(threshold_registry, fallback = c(0.05, 0.10, 0.15, 0.20)) {
  out <- tibble::as_tibble(threshold_registry)
  if (nrow(out) == 0L || !("threshold" %in% names(out))) {
    return(as.numeric(fallback))
  }
  sort(unique(safe_numeric(out$threshold)))
}

# 🔴 Define: Stage6 carry-forward readers ===============================
## 🟠 Define: Stage6 screening import and model lookup ===============================
empty_stage6_lookup <- function() {
  tibble(
    dataset = character(),
    formula_anchor = character(),
    cure_model_eligibility_flag = character(),
    receus_primary_class = character(),
    cure_presence_support_flag = character(),
    followup_not_contradicted_flag = character(),
    screening_note = character()
  )
}

normalize_stage6_formula_anchor <- function(x) {
  x <- trimws(as.character(x))
  out <- case_when(
    x %in% c("base", "Base") ~ "base",
    x %in% c("interaction", "Interaction") ~ "interaction",
    x %in% c("site_added", "Site-added", "site-added") ~ "site_added",
    x %in% c("site_interaction", "Site + interaction", "site+interaction") ~ "site_interaction",
    TRUE ~ x
  )
  out
}

expand_stage6_dataset_key <- function(dataset_value, formula_anchor_value = NA_character_) {
  dataset_value <- trimws(as.character(dataset_value %||% NA_character_))
  formula_anchor_value <- trimws(as.character(formula_anchor_value %||% NA_character_))
  
  if (identical(dataset_value, "merged__site_free")) {
    return(tibble(dataset = "merged", formula_anchor = c("base", "interaction")))
  }
  if (identical(dataset_value, "merged__site_adjusted")) {
    return(tibble(dataset = "merged", formula_anchor = c("site_added", "site_interaction")))
  }
  if (dataset_value %in% c("PNU", "SNU", "merged")) {
    return(tibble(
      dataset = dataset_value,
      formula_anchor = if (!is.na(formula_anchor_value) && nzchar(formula_anchor_value)) formula_anchor_value else "ALL"
    ))
  }
  
  tibble(
    dataset = dataset_value,
    formula_anchor = if (!is.na(formula_anchor_value) && nzchar(formula_anchor_value)) formula_anchor_value else "ALL"
  )
}

read_stage6_screening_lookup <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(empty_stage6_lookup())
  }
  df <- read_delimited_or_rds(path)
  df <- tibble::as_tibble(df)
  
  dataset_col <- first_existing_name(df, c("dataset_key", "dataset", "source_dataset", "cohort", "screening_dataset_key"))
  formula_col <- first_existing_name(df, c("formula_anchor", "formula_name", "formula_type"))
  eligibility_col <- first_existing_name(df, c("cure_model_eligibility_flag", "stage6_final_class", "final_decision_flag", "screening_flag"))
  receus_col <- first_existing_name(df, c("receus_primary_class", "primary_gate_flag"))
  presence_col <- first_existing_name(df, c("cure_presence_support_flag", "presence_modifier_flag"))
  followup_col <- first_existing_name(df, c("followup_not_contradicted_flag", "followup_contradiction_flag"))
  note_col <- first_existing_name(df, c("screening_note", "screening_detail", "primary_gate_method", "screening_context"))
  
  if (is.null(dataset_col) || is.null(eligibility_col)) {
    return(empty_stage6_lookup())
  }
  
  raw_out <- tibble(
    dataset_value = as.character(df[[dataset_col]]),
    formula_anchor_value = if (!is.null(formula_col)) normalize_stage6_formula_anchor(df[[formula_col]]) else NA_character_,
    cure_model_eligibility_flag = as.character(df[[eligibility_col]]),
    receus_primary_class = if (!is.null(receus_col)) as.character(df[[receus_col]]) else NA_character_,
    cure_presence_support_flag = if (!is.null(presence_col)) as.character(df[[presence_col]]) else NA_character_,
    followup_not_contradicted_flag = if (!is.null(followup_col)) as.character(df[[followup_col]]) else NA_character_,
    screening_note = if (!is.null(note_col)) as.character(df[[note_col]]) else NA_character_
  )
  
  out <- bind_rows(lapply(seq_len(nrow(raw_out)), function(i) {
    one <- raw_out[i, , drop = FALSE]
    expanded <- expand_stage6_dataset_key(one$dataset_value[[1]], one$formula_anchor_value[[1]])
    bind_cols(
      expanded,
      one %>% select(-dataset_value, -formula_anchor_value)[rep(1L, nrow(expanded)), , drop = FALSE]
    )
  })) %>%
    filter(dataset %in% c("PNU", "SNU", "merged")) %>%
    mutate(formula_anchor = ifelse(is.na(formula_anchor) | formula_anchor == "", "ALL", formula_anchor)) %>%
    group_by(dataset, formula_anchor) %>%
    summarise(
      cure_model_eligibility_flag = paste(unique(na.omit(cure_model_eligibility_flag)), collapse = "|"),
      receus_primary_class = paste(unique(na.omit(receus_primary_class)), collapse = "|"),
      cure_presence_support_flag = paste(unique(na.omit(cure_presence_support_flag)), collapse = "|"),
      followup_not_contradicted_flag = paste(unique(na.omit(followup_not_contradicted_flag)), collapse = "|"),
      screening_note = paste(unique(na.omit(screening_note)), collapse = "; "),
      .groups = "drop"
    )
  
  if (nrow(out) == 0L) {
    return(empty_stage6_lookup())
  }
  out
}

build_stage6_model_lookup <- function(screening_lookup, model_grid) {
  if (nrow_or_zero(screening_lookup) == 0L) {
    return(model_grid %>%
             transmute(
               model_id = model_id,
               cure_model_eligibility_flag = NA_character_,
               receus_primary_class = NA_character_,
               cure_presence_support_flag = NA_character_,
               followup_not_contradicted_flag = NA_character_,
               screening_note = NA_character_
             ))
  }
  
  model_base <- model_grid %>%
    distinct(model_id, dataset, formula_anchor)
  
  exact_lookup <- model_base %>%
    left_join(
      screening_lookup %>% filter(formula_anchor != "ALL"),
      by = c("dataset", "formula_anchor")
    )
  
  all_lookup <- model_base %>%
    left_join(
      screening_lookup %>%
        filter(formula_anchor == "ALL") %>%
        select(dataset, ends_with("flag"), screening_note, receus_primary_class) %>%
        rename_with(~ paste0(.x, "_all"), -dataset),
      by = "dataset"
    )
  
  exact_lookup %>%
    left_join(all_lookup, by = c("model_id", "dataset", "formula_anchor")) %>%
    transmute(
      model_id = model_id,
      cure_model_eligibility_flag = coalesce(cure_model_eligibility_flag, cure_model_eligibility_flag_all),
      receus_primary_class = coalesce(receus_primary_class, receus_primary_class_all),
      cure_presence_support_flag = coalesce(cure_presence_support_flag, cure_presence_support_flag_all),
      followup_not_contradicted_flag = coalesce(followup_not_contradicted_flag, followup_not_contradicted_flag_all),
      screening_note = coalesce(screening_note, screening_note_all)
    )
}

# 🔴 Define: Stage8A import helpers for 8B-versus-8A deltas ===============================
## 🟠 Define: Stage8A readers and harmonizers ===============================
load_stage8a_outputs <- function(model_registry_csv, cohort_yearly_csv, classification_csv) {
  out <- list(
    model_registry = tibble(),
    posterior_cohort_yearly = tibble(),
    posterior_classification = tibble()
  )
  if (file.exists(model_registry_csv)) {
    out$model_registry <- tibble::as_tibble(read_delimited_or_rds(model_registry_csv))
  }
  if (file.exists(cohort_yearly_csv)) {
    out$posterior_cohort_yearly <- tibble::as_tibble(read_delimited_or_rds(cohort_yearly_csv))
  }
  if (file.exists(classification_csv)) {
    out$posterior_classification <- tibble::as_tibble(read_delimited_or_rds(classification_csv))
  }
  out
}

augment_stage8a_with_model_keys <- function(stage8a_outputs) {
  reg <- tibble::as_tibble(stage8a_outputs$model_registry)
  cohort <- tibble::as_tibble(stage8a_outputs$posterior_cohort_yearly)
  class_tbl <- tibble::as_tibble(stage8a_outputs$posterior_classification)
  
  if (nrow_or_zero(reg) == 0L) {
    return(list(
      posterior_cohort_yearly = tibble(),
      posterior_classification = tibble()
    ))
  }
  
  key_tbl <- reg %>%
    transmute(
      dataset = as.character(dataset),
      model_id = as.character(model_id),
      structural_model_id = as.character(structural_model_id),
      family_code = as.character(family_code),
      formula_anchor = as.character(formula_anchor)
    )
  
  cohort_out <- cohort %>%
    left_join(key_tbl, by = c("dataset", "model_id", "formula_anchor")) %>%
    mutate(horizon_year = as.integer(safe_numeric(horizon_year)))
  
  class_out <- class_tbl %>%
    left_join(key_tbl, by = c("dataset", "model_id", "formula_anchor")) %>%
    mutate(
      horizon_year = as.integer(safe_numeric(horizon_year)),
      threshold = as.numeric(safe_numeric(threshold))
    )
  
  list(
    posterior_cohort_yearly = cohort_out,
    posterior_classification = class_out
  )
}

# 🔴 Define: Stage8B structural grid and priors ===============================
## 🟠 Define: model grid, prior specs, and design bundles ===============================
build_stage8b_model_grid <- function(
    include_merged_incidence_site_supplementary = TRUE,
    fit_prior_branches = c("anchor_informed", "neutral_no_external_info"),
    fit_site_prior_families = c("normal_0_1_main", "student_t3_0_1_sensitivity")
) {
  single_rows <- tibble::tribble(
    ~dataset, ~structural_model_id, ~formula_anchor, ~transition_latency_branch, ~remission_branch, ~incidence_site_indicator, ~latency_site_indicator, ~remission_site_indicator, ~interaction_indicator, ~site_placement_label, ~is_supplementary_branch,
    "PNU", "PNU-L0R0", "base", "L0", "R0", FALSE, FALSE, FALSE, FALSE, "site_in_neither", FALSE,
    "PNU", "PNU-L1R1", "interaction", "L1", "R1", FALSE, FALSE, FALSE, TRUE, "site_in_neither", FALSE,
    "SNU", "SNU-L0R0", "base", "L0", "R0", FALSE, FALSE, FALSE, FALSE, "site_in_neither", FALSE,
    "SNU", "SNU-L1R1", "interaction", "L1", "R1", FALSE, FALSE, FALSE, TRUE, "site_in_neither", FALSE
  )
  
  merged_rows <- tibble::tribble(
    ~dataset, ~structural_model_id, ~formula_anchor, ~transition_latency_branch, ~remission_branch, ~incidence_site_indicator, ~latency_site_indicator, ~remission_site_indicator, ~interaction_indicator, ~site_placement_label, ~is_supplementary_branch,
    "merged", "MERGED-L0S0-R0S0", "base", "L0S0", "R0S0", FALSE, FALSE, FALSE, FALSE, "site_in_neither", FALSE,
    "merged", "MERGED-L1S0-R1S0", "interaction", "L1S0", "R1S0", FALSE, FALSE, FALSE, TRUE, "site_in_neither", FALSE,
    "merged", "MERGED-L0S1-R0S1", "site_added", "L0S1", "R0S1", FALSE, TRUE, TRUE, FALSE, "site_in_latency_only", FALSE,
    "merged", "MERGED-L1S1-R1S1", "site_interaction", "L1S1", "R1S1", FALSE, TRUE, TRUE, TRUE, "site_in_latency_only", FALSE
  )
  
  if (isTRUE(include_merged_incidence_site_supplementary)) {
    merged_rows <- bind_rows(
      merged_rows,
      tibble::tribble(
        ~dataset, ~structural_model_id, ~formula_anchor, ~transition_latency_branch, ~remission_branch, ~incidence_site_indicator, ~latency_site_indicator, ~remission_site_indicator, ~interaction_indicator, ~site_placement_label, ~is_supplementary_branch,
        "merged", "MERGED-I1-L0S0-R0S0", "base", "L0S0", "R0S0", TRUE, FALSE, FALSE, FALSE, "site_in_incidence_only", TRUE,
        "merged", "MERGED-I1-L1S0-R1S0", "interaction", "L1S0", "R1S0", TRUE, FALSE, FALSE, TRUE, "site_in_incidence_only", TRUE,
        "merged", "MERGED-I1-L0S1-R0S1", "site_added", "L0S1", "R0S1", TRUE, TRUE, TRUE, FALSE, "site_in_both", TRUE,
        "merged", "MERGED-I1-L1S1-R1S1", "site_interaction", "L1S1", "R1S1", TRUE, TRUE, TRUE, TRUE, "site_in_both", TRUE
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
  
  base_grid <- tidyr::crossing(bind_rows(single_rows, merged_rows), family_rows) %>%
    mutate(
      has_any_site_term = incidence_site_indicator | latency_site_indicator | remission_site_indicator
    )
  
  rows_out <- list()
  
  for (ii in seq_len(nrow(base_grid))) {
    row_i <- base_grid[ii, , drop = FALSE]
    site_prior_candidates <- if (isTRUE(row_i$has_any_site_term[[1]])) fit_site_prior_families else "not_applicable"
    cross_i <- tidyr::crossing(
      row_i,
      tibble(
        prior_branch = fit_prior_branches,
        site_prior_family = site_prior_candidates
      )
    )
    rows_out[[ii]] <- cross_i
  }
  
  bind_rows(rows_out) %>%
    mutate(
      dataset_key = dataset,
      branch = "Stage8B",
      risk_scale = "transition_cif_competing",
      retained_fit_id = paste(structural_model_id, family_code, site_prior_family, sep = "__"),
      model_id = paste(structural_model_id, family_code, prior_branch, site_prior_family, sep = "__")
    ) %>%
    arrange(factor(dataset, levels = c("PNU", "SNU", "merged")), structural_model_id, family_code, prior_branch, site_prior_family)
}

build_stage8b_prior_specs <- function() {
  list(
    anchor_informed = list(
      prior_branch = "anchor_informed",
      alpha_prior_center = -9.581369553169,
      mu_beta_inc_base = c(0.419871845822, 0.907608052926, 0.586202561451, 0.466865123863, 0.037997248763),
      sd_beta_inc_base = c(0.132789397422, 0.173731076538, 0.191221553945, 0.270393197518, 0.302838606651),
      sd_delta = 3.5,
      sd_gamma0 = 2.5,
      sd_gamma = 1.0,
      sd_rho_W = 0.35,
      sd_log_sigma_LN = 0.50,
      sd_psi_LL = 0.50,
      sd_xi0 = 2.5,
      sd_rem = 1.5,
      sd_log_rho_rem = 1.5,
      sd_site = 1.0
    ),
    neutral_no_external_info = list(
      prior_branch = "neutral_no_external_info",
      alpha_prior_center = 0.0,
      mu_beta_inc_base = rep(0, 5),
      sd_beta_inc_base = rep(2.0, 5),
      sd_delta = 3.5,
      sd_gamma0 = 2.5,
      sd_gamma = 1.0,
      sd_rho_W = 0.35,
      sd_log_sigma_LN = 0.50,
      sd_psi_LL = 0.50,
      sd_xi0 = 2.5,
      sd_rem = 1.5,
      sd_log_rho_rem = 1.5,
      sd_site = 1.0
    )
  )
}

make_stage8b_design_bundle <- function(df, model_row, prior_spec, snu_label, remission_cut_years) {
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
  
  if (isTRUE(model_row$incidence_site_indicator)) {
    X_inc <- cbind(X_inc, site_SNU = s_i)
    mu_beta_inc <- c(mu_beta_inc, 0)
    sd_beta_inc <- c(sd_beta_inc, prior_spec$sd_site)
  }
  
  a_i <- as.numeric(df$age_s)
  az_i <- a_i * z_i
  
  X_lat <- switch(
    model_row$transition_latency_branch,
    L0 = cbind(age_s = a_i, sex_num = z_i),
    L1 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i),
    L0S0 = cbind(age_s = a_i, sex_num = z_i),
    L1S0 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i),
    L0S1 = cbind(age_s = a_i, sex_num = z_i, site_SNU = s_i),
    L1S1 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i, site_SNU = s_i),
    stop(sprintf("Unknown transition latency branch `%s`.", model_row$transition_latency_branch), call. = FALSE)
  )
  
  X_rem <- switch(
    model_row$remission_branch,
    R0 = cbind(age_s = a_i, sex_num = z_i),
    R1 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i),
    R0S0 = cbind(age_s = a_i, sex_num = z_i),
    R1S0 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i),
    R0S1 = cbind(age_s = a_i, sex_num = z_i, site_SNU = s_i),
    R1S1 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i, site_SNU = s_i),
    stop(sprintf("Unknown remission branch `%s`.", model_row$remission_branch), call. = FALSE)
  )
  
  sd_gamma_lat <- rep(prior_spec$sd_gamma, ncol(X_lat))
  sd_beta_rem <- rep(prior_spec$sd_rem, ncol(X_rem))
  if (grepl("S1$", model_row$transition_latency_branch) && isTRUE(model_row$latency_site_indicator)) {
    sd_gamma_lat[ncol(X_lat)] <- prior_spec$sd_site
  }
  if (grepl("S1$", model_row$remission_branch) && isTRUE(model_row$remission_site_indicator)) {
    sd_beta_rem[ncol(X_rem)] <- prior_spec$sd_site
  }
  
  rem_cuts <- as.numeric(remission_cut_years)
  if (length(rem_cuts) < 2L || any(diff(rem_cuts) <= 0)) {
    stop("`remission_cut_years` must be a strictly increasing numeric vector.", call. = FALSE)
  }
  
  list(
    X_inc = unclass(as.matrix(X_inc)),
    X_lat = unclass(as.matrix(X_lat)),
    X_rem = unclass(as.matrix(X_rem)),
    mu_beta_inc = as.numeric(mu_beta_inc),
    sd_beta_inc = as.numeric(sd_beta_inc),
    alpha_prior_center = as.numeric(prior_spec$alpha_prior_center),
    sd_gamma_lat = as.numeric(sd_gamma_lat),
    sd_beta_rem = as.numeric(sd_beta_rem),
    rem_cuts = rem_cuts,
    time = as.numeric(df$time_year),
    status = as.integer(df$status_num),
    id_df = df %>% select(unique_person_id, id, site, sex_num, age_exact_entry, age_s)
  )
}


# 🔴 Define: observed competing-risk references and support metadata ===============================
## 🟠 Define: Aalen-Johansen, censoring IPCW, and support joins ===============================
km_step_eval <- function(survfit_obj, times) {
  base_times <- survfit_obj$time
  base_surv <- survfit_obj$surv
  vapply(
    times,
    FUN.VALUE = numeric(1),
    FUN = function(tt) {
      idx <- max(c(0L, which(base_times <= tt)))
      if (idx == 0L) {
        1
      } else {
        base_surv[[idx]]
      }
    }
  )
}

build_aalen_johansen_reference <- function(df, horizons) {
  event_times <- sort(unique(df$time_year[df$status_num %in% c(1L, 2L)]))
  if (length(event_times) == 0L) {
    return(
      tibble(
        horizon_year = as.integer(horizons),
        observed_transition_cif = 0,
        observed_remission_cif = 0,
        observed_all_event_free = 1
      )
    )
  }
  
  n_total <- nrow(df)
  S_all <- 1
  F1 <- 0
  F2 <- 0
  
  aj_tbl <- tibble(
    time_year = event_times,
    observed_transition_cif = NA_real_,
    observed_remission_cif = NA_real_,
    observed_all_event_free = NA_real_
  )
  
  n_risk <- n_total
  for (ii in seq_along(event_times)) {
    tt <- event_times[[ii]]
    d1 <- sum(df$status_num == 1L & df$time_year == tt)
    d2 <- sum(df$status_num == 2L & df$time_year == tt)
    c0 <- sum(df$status_num == 0L & df$time_year == tt)
    
    if (n_risk <= 0L) {
      aj_tbl$observed_transition_cif[[ii]] <- F1
      aj_tbl$observed_remission_cif[[ii]] <- F2
      aj_tbl$observed_all_event_free[[ii]] <- S_all
      next
    }
    
    F1 <- F1 + S_all * d1 / n_risk
    F2 <- F2 + S_all * d2 / n_risk
    S_all <- S_all * (1 - (d1 + d2) / n_risk)
    
    aj_tbl$observed_transition_cif[[ii]] <- F1
    aj_tbl$observed_remission_cif[[ii]] <- F2
    aj_tbl$observed_all_event_free[[ii]] <- S_all
    
    n_risk <- n_risk - d1 - d2 - c0
  }
  
  eval_step <- function(values, times, horizons_eval) {
    vapply(horizons_eval, function(hh) {
      idx <- max(c(0L, which(times <= hh)))
      if (idx == 0L) {
        return(0)
      }
      values[[idx]]
    }, numeric(1))
  }
  
  tibble(
    horizon_year = as.integer(horizons),
    observed_transition_cif = eval_step(aj_tbl$observed_transition_cif, aj_tbl$time_year, horizons),
    observed_remission_cif = eval_step(aj_tbl$observed_remission_cif, aj_tbl$time_year, horizons),
    observed_all_event_free = {
      vals <- vapply(horizons, function(hh) {
        idx <- max(c(0L, which(aj_tbl$time_year <= hh)))
        if (idx == 0L) {
          return(1)
        }
        aj_tbl$observed_all_event_free[[idx]]
      }, numeric(1))
      vals
    }
  )
}

build_stage8b_ipcw_reference <- function(df, horizons) {
  censor_fit <- survival::survfit(survival::Surv(time_year, status_num == 0L) ~ 1, data = df)
  aj_ref <- build_aalen_johansen_reference(df, horizons)
  
  out_rows <- lapply(horizons, function(hh) {
    G_h <- pmax(km_step_eval(censor_fit, hh), 1e-8)
    
    G_tminus <- pmax(km_step_eval(censor_fit, pmax(df$time_year - 1e-10, 0)), 1e-8)
    G_t <- pmax(km_step_eval(censor_fit, df$time_year), 1e-8)
    
    w_case <- ifelse(df$status_num == 1L & df$time_year <= hh, 1 / G_tminus, 0)
    w_control <- ifelse(
      df$status_num == 2L & df$time_year <= hh,
      1 / G_t,
      ifelse(df$time_year > hh, 1 / G_h, 0)
    )
    
    aj_row <- aj_ref %>% filter(horizon_year == hh)
    
    tibble(
      horizon_year = as.integer(hh),
      observed_transition_cif = aj_row$observed_transition_cif[[1]],
      observed_remission_cif = aj_row$observed_remission_cif[[1]],
      observed_all_event_free = aj_row$observed_all_event_free[[1]],
      denom_case = sum(w_case),
      denom_control = sum(w_control),
      G_h = G_h,
      w_case = list(w_case),
      w_control = list(w_control)
    )
  })
  
  bind_rows(out_rows)
}

make_stage8b_support_registry <- function(horizon_registry_stage1, datasets = c("PNU", "SNU", "merged"), horizons = 1:10) {
  out <- tibble::as_tibble(horizon_registry_stage1)
  if (nrow(out) == 0L) {
    return(
      tidyr::expand_grid(dataset = datasets, horizon_year = as.integer(horizons)) %>%
        mutate(
          support_tier = NA_character_,
          support_tier_standard = NA_character_,
          horizon_evidence_class = NA_character_,
          claim_restriction_flag = NA_character_,
          interpretation_note = NA_character_
        )
    )
  }
  
  out %>%
    transmute(
      dataset = as.character(dataset),
      dataset_key = as.character(dataset_key %||% dataset),
      horizon_year = as.integer(safe_numeric(horizon_year)),
      support_tier = as.character(support_tier),
      support_tier_standard = as.character(support_tier_standard %||% support_tier),
      horizon_evidence_class = as.character(horizon_evidence_class),
      claim_restriction_flag = as.character(claim_restriction_flag),
      interpretation_note = as.character(interpretation_note)
    ) %>%
    filter(dataset %in% datasets, horizon_year %in% horizons) %>%
    distinct(dataset, horizon_year, .keep_all = TRUE) %>%
    arrange(factor(dataset, levels = datasets), horizon_year)
}

# 🔴 Define: Stage8B posterior math and numerical integration ===============================
## 🟠 Define: transition-family, remission, and CIF integration helpers ===============================
interval_exposure_vector <- function(t, cuts) {
  lower <- cuts[-length(cuts)]
  upper <- cuts[-1]
  pmax(pmin(t, upper) - lower, 0)
}

interval_index_for_time <- function(t, cuts) {
  lower <- cuts[-length(cuts)]
  upper <- cuts[-1]
  idx <- which(t >= lower & t < upper)
  if (length(idx) == 0L) {
    return(length(lower))
  }
  idx[[1L]]
}

compute_stage8b_linear_terms <- function(draws, X_inc, X_lat, X_rem, alpha_prior_center) {
  eta_inc <- draws$beta_inc %*% t(X_inc)
  eta_inc <- sweep(eta_inc, 1, alpha_prior_center + draws$delta0, FUN = "+")
  pi_mat <- plogis(eta_inc)
  
  mu_lat <- draws$gamma_lat %*% t(X_lat)
  mu_lat <- sweep(mu_lat, 1, draws$gamma0, FUN = "+")
  median_mat <- exp(mu_lat)
  
  eta_rem <- draws$beta_rem %*% t(X_rem)
  eta_rem <- sweep(eta_rem, 1, draws$xi0, FUN = "+")
  
  list(
    pi_mat = pmin(pmax(pi_mat, 1e-12), 1 - 1e-12),
    cure_prob_mat = pmin(pmax(1 - pi_mat, 1e-12), 1 - 1e-12),
    mu_lat_mat = mu_lat,
    median_mat = pmin(pmax(median_mat, 1e-8), 1e8),
    eta_rem_mat = eta_rem
  )
}

transition_family_components <- function(t, family_code, mu_lat_mat, median_mat, draws) {
  if (family_code == "E") {
    lambda <- median_mat / log(2)
    Su <- exp(-t / lambda)
    haz <- 1 / lambda
    fu <- haz * Su
  } else if (family_code == "W") {
    k <- exp(draws$rho_W)
    lambda <- median_mat / ((log(2))^(1 / k))
    ratio <- t / lambda
    Su <- exp(-(ratio^k))
    haz <- (k / lambda) * (ratio^(k - 1))
    fu <- haz * Su
  } else if (family_code == "LN") {
    sigma <- exp(draws$log_sigma_LN)
    z <- (log(t) - mu_lat_mat) / sigma
    Su <- 1 - pnorm(z)
    log_pdf <- dnorm(z, log = TRUE) - log(t) - log(sigma)
    fu <- exp(log_pdf)
    haz <- fu / pmax(Su, 1e-12)
  } else if (family_code == "LL") {
    k <- exp(-draws$psi_LL)
    lambda <- median_mat
    ratio <- t / lambda
    Su <- 1 / (1 + ratio^k)
    haz <- (k / lambda) * (ratio^(k - 1)) / (1 + ratio^k)
    fu <- haz * Su
  } else {
    stop(sprintf("Unknown family code `%s`.", family_code), call. = FALSE)
  }
  
  list(
    Su = pmin(pmax(Su, 1e-12), 1 - 1e-12),
    haz = pmax(haz, 1e-12),
    fu = pmax(fu, 1e-12)
  )
}

remission_components_piecewise <- function(t, eta_rem_mat, log_rho_rem_mat, rem_cuts) {
  expo <- interval_exposure_vector(t, rem_cuts)
  base_draw <- as.vector(exp(log_rho_rem_mat) %*% expo)
  H2 <- sweep(exp(eta_rem_mat), 1, base_draw, FUN = "*")
  
  interval_idx <- interval_index_for_time(t, rem_cuts)
  rho_t_draw <- exp(log_rho_rem_mat[, interval_idx])
  h2 <- sweep(exp(eta_rem_mat), 1, rho_t_draw, FUN = "*")
  
  list(
    H2 = pmax(H2, 0),
    G2 = pmin(pmax(exp(-H2), 1e-12), 1),
    h2 = pmax(h2, 1e-12)
  )
}

build_stage8b_prediction_trajectories <- function(state, draws, model_row, horizons, integration_step_year, rem_cuts) {
  max_h <- max(horizons)
  bounds <- seq(0, max_h, by = integration_step_year)
  if (tail(bounds, 1) < max_h) {
    bounds <- c(bounds, max_h)
  }
  mids <- (bounds[-1] + bounds[-length(bounds)]) / 2
  widths <- diff(bounds)
  
  n_draw <- nrow(state$pi_mat)
  n_subj <- ncol(state$pi_mat)
  
  F1_cum <- matrix(0, nrow = n_draw, ncol = n_subj)
  F2_cum <- matrix(0, nrow = n_draw, ncol = n_subj)
  
  horizon_results <- vector("list", length(horizons))
  names(horizon_results) <- as.character(horizons)
  
  step_to_horizon <- findInterval(horizons, vec = bounds[-1], left.open = FALSE)
  
  horizon_counter <- 1L
  for (jj in seq_along(mids)) {
    tt <- max(mids[[jj]], 1e-8)
    comp1 <- transition_family_components(tt, model_row$family_code[[1]], state$mu_lat_mat, state$median_mat, draws)
    comp2 <- remission_components_piecewise(tt, state$eta_rem_mat, draws$log_rho_rem, rem_cuts)
    
    M1 <- (1 - state$pi_mat) + state$pi_mat * comp1$Su
    trans_inc <- state$pi_mat * comp1$fu * comp2$G2 * widths[[jj]]
    rem_inc <- comp2$h2 * comp2$G2 * M1 * widths[[jj]]
    
    F1_cum <- F1_cum + trans_inc
    F2_cum <- F2_cum + rem_inc
    
    while (horizon_counter <= length(horizons) && step_to_horizon[[horizon_counter]] == jj) {
      hh <- horizons[[horizon_counter]]
      comp1_h <- transition_family_components(max(hh, 1e-8), model_row$family_code[[1]], state$mu_lat_mat, state$median_mat, draws)
      comp2_h <- remission_components_piecewise(max(hh, 1e-8), state$eta_rem_mat, draws$log_rho_rem, rem_cuts)
      
      M1_h <- (1 - state$pi_mat) + state$pi_mat * comp1_h$Su
      S_all_h <- comp2_h$G2 * M1_h
      
      hazard_trans_pop <- (state$pi_mat * comp1_h$fu) / pmax(M1_h, 1e-12)
      
      horizon_results[[as.character(hh)]] <- list(
        transition_cif = pmin(pmax(F1_cum, 0), 1),
        remission_cif = pmin(pmax(F2_cum, 0), 1),
        all_event_free = pmin(pmax(S_all_h, 0), 1),
        uncured_survival = comp1_h$Su,
        transition_population_hazard = pmax(hazard_trans_pop, 0),
        remission_hazard = comp2_h$h2
      )
      horizon_counter <- horizon_counter + 1L
    }
  }
  
  horizon_results
}

compute_stage8b_degeneracy <- function(pi_mat, median_mat, transition_cif_supported) {
  near_zero_pi <- rowMeans(pi_mat < tiny_susceptible_prob) > degenerate_subject_fraction_threshold
  near_one_pi <- rowMeans(pi_mat > huge_susceptible_prob) > degenerate_subject_fraction_threshold
  near_zero_median <- rowMeans(median_mat < tiny_median_years) > degenerate_subject_fraction_threshold
  huge_median_flat <- if (length(transition_cif_supported) > 0L) {
    mean_risk_supported <- Reduce("+", transition_cif_supported) / length(transition_cif_supported)
    (rowMeans(median_mat > huge_median_years) > degenerate_subject_fraction_threshold) &
      (rowMeans(mean_risk_supported < 0.01) > degenerate_subject_fraction_threshold)
  } else {
    rep(FALSE, nrow(pi_mat))
  }
  any_problem <- near_zero_pi | near_one_pi | near_zero_median | huge_median_flat
  
  tibble(
    near_zero_pi_rate = mean(near_zero_pi),
    near_one_pi_rate = mean(near_one_pi),
    near_zero_median_rate = mean(near_zero_median),
    huge_median_flat_risk_rate = mean(huge_median_flat),
    degenerate_flag = mean(any_problem) > degenerate_draw_fraction_threshold
  )
}

compute_stage8b_classification_summary <- function(risk_draws, horizon_row, thresholds) {
  prevalence <- as.numeric(horizon_row$observed_transition_cif)
  w_case <- unlist(horizon_row$w_case)
  w_control <- unlist(horizon_row$w_control)
  denom_case <- as.numeric(horizon_row$denom_case)
  denom_control <- as.numeric(horizon_row$denom_control)
  n_subject <- ncol(risk_draws)
  
  out_list <- vector("list", length(thresholds))
  for (jj in seq_along(thresholds)) {
    thr <- thresholds[[jj]]
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
    ppv_draw <- ifelse(pos_rate_draw > 0, prevalence * tpr_draw / pos_rate_draw, NA_real_)
    fp_burden_draw <- (1 - prevalence) * fpr_draw
    fp100_draw <- 100 * fp_burden_draw
    fp_count_draw <- n_subject * fp_burden_draw
    nb_draw <- prevalence * tpr_draw - (1 - prevalence) * fpr_draw * (thr / (1 - thr))
    
    out_list[[jj]] <- tibble(
      threshold = thr,
      positive_rate_mean = mean(pos_rate_draw, na.rm = TRUE),
      positive_rate_q025 = stats::quantile(pos_rate_draw, 0.025, na.rm = TRUE, names = FALSE),
      positive_rate_q50 = stats::quantile(pos_rate_draw, 0.500, na.rm = TRUE, names = FALSE),
      positive_rate_q975 = stats::quantile(pos_rate_draw, 0.975, na.rm = TRUE, names = FALSE),
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
  
  bind_rows(out_list)
}

classify_hazard_shape <- function(hazard_values, tol = 1e-6) {
  hv <- as.numeric(hazard_values)
  if (length(hv) < 2L || all(is.na(hv))) {
    return("undetermined")
  }
  diffs <- diff(hv)
  diffs[abs(diffs) < tol] <- 0
  
  if (all(diffs >= 0)) {
    if (all(diffs == 0)) return("flat")
    return("monotone_increasing")
  }
  if (all(diffs <= 0)) {
    if (all(diffs == 0)) return("flat")
    return("monotone_decreasing")
  }
  peak_idx <- which.max(hv)
  if (peak_idx > 1L && peak_idx < length(hv)) {
    left_ok <- all(diff(hv[1:peak_idx]) >= -tol)
    right_ok <- all(diff(hv[peak_idx:length(hv)]) <= tol)
    if (left_ok && right_ok) {
      return("unimodal")
    }
  }
  "irregular"
}


# 🔴 Define: prior predictive summaries, anchor updates, and diagnostics ===============================
## 🟠 Define: summary helpers, prior predictive simulation, and anchor-update tables ===============================
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

ppc_horizons_for_dataset <- function(dataset_name) {
  if (dataset_name == "PNU") {
    return(c(1L, 2L))
  }
  c(1L, 2L, 5L)
}

get_stage8b_site_prior_indices <- function(model_row, design_bundle) {
  list(
    inc_idx = if (isTRUE(model_row$incidence_site_indicator[[1]])) ncol(design_bundle$X_inc) else NA_integer_,
    lat_idx = if (isTRUE(model_row$latency_site_indicator[[1]])) ncol(design_bundle$X_lat) else NA_integer_,
    rem_idx = if (isTRUE(model_row$remission_site_indicator[[1]])) ncol(design_bundle$X_rem) else NA_integer_
  )
}

apply_stage8b_site_prior_draws <- function(mat, idx, site_prior_family, scale) {
  if (is.na(idx) || idx <= 0L || idx > ncol(mat)) {
    return(mat)
  }
  if (identical(site_prior_family, "student_t3_0_1_sensitivity")) {
    mat[, idx] <- stats::rt(nrow(mat), df = 3) * scale
  }
  mat
}

empty_stage8b_prior_predictive_summary <- function() {
  tibble(
    dataset_key = character(),
    branch = character(),
    model_id = character(),
    retained_fit_id = character(),
    risk_scale = character(),
    prior_branch = character(),
    site_prior_family = character(),
    metric = character(),
    horizon = integer(),
    mean = double(),
    sd = double(),
    q025 = double(),
    q50 = double(),
    q975 = double(),
    prior_degenerate_flag = logical()
  )
}

simulate_stage8b_prior_predictive <- function(dataset_df, model_row, design_bundle, prior_spec, n_draws, horizons_eval) {
  idx_info <- get_stage8b_site_prior_indices(model_row, design_bundle)
  K_inc <- ncol(design_bundle$X_inc)
  K_lat <- ncol(design_bundle$X_lat)
  K_rem <- ncol(design_bundle$X_rem)
  J_rem <- length(design_bundle$rem_cuts) - 1L
  
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
    xi0 = rnorm(n_draws, 0, prior_spec$sd_xi0),
    beta_rem = matrix(
      rnorm(n_draws * K_rem, mean = 0, sd = rep(design_bundle$sd_beta_rem, each = n_draws)),
      nrow = n_draws,
      byrow = FALSE
    ),
    log_rho_rem = matrix(
      rnorm(n_draws * J_rem, mean = 0, sd = prior_spec$sd_log_rho_rem),
      nrow = n_draws,
      byrow = FALSE
    ),
    rho_W = rnorm(n_draws, 0, prior_spec$sd_rho_W),
    log_sigma_LN = rnorm(n_draws, 0, prior_spec$sd_log_sigma_LN),
    psi_LL = rnorm(n_draws, 0, prior_spec$sd_psi_LL)
  )
  
  sim_draws$beta_inc <- apply_stage8b_site_prior_draws(
    sim_draws$beta_inc,
    idx_info$inc_idx,
    model_row$site_prior_family[[1]],
    prior_spec$sd_site
  )
  sim_draws$gamma_lat <- apply_stage8b_site_prior_draws(
    sim_draws$gamma_lat,
    idx_info$lat_idx,
    model_row$site_prior_family[[1]],
    prior_spec$sd_site
  )
  sim_draws$beta_rem <- apply_stage8b_site_prior_draws(
    sim_draws$beta_rem,
    idx_info$rem_idx,
    model_row$site_prior_family[[1]],
    prior_spec$sd_site
  )
  
  state <- compute_stage8b_linear_terms(
    draws = sim_draws,
    X_inc = design_bundle$X_inc,
    X_lat = design_bundle$X_lat,
    X_rem = design_bundle$X_rem,
    alpha_prior_center = design_bundle$alpha_prior_center
  )
  
  pred <- build_stage8b_prediction_trajectories(
    state = state,
    draws = sim_draws,
    model_row = model_row,
    horizons = sort(unique(as.integer(horizons_eval))),
    integration_step_year = integration_step_year,
    rem_cuts = design_bundle$rem_cuts
  )
  
  supported_horizons <- intersect(as.character(ppc_horizons_for_dataset(model_row$dataset[[1]])), names(pred))
  transition_supported <- lapply(supported_horizons, function(hh) pred[[hh]]$transition_cif)
  degeneracy <- compute_stage8b_degeneracy(state$pi_mat, state$median_mat, transition_supported)
  
  out_rows <- list(
    tibble(
      dataset_key = model_row$dataset[[1]],
      branch = "Stage8B",
      model_id = model_row$model_id[[1]],
      retained_fit_id = model_row$retained_fit_id[[1]],
      risk_scale = model_row$risk_scale[[1]],
      prior_branch = model_row$prior_branch[[1]],
      site_prior_family = model_row$site_prior_family[[1]],
      metric = "susceptible_fraction",
      horizon = NA_integer_,
      summary_scalar(rowMeans(state$pi_mat)),
      prior_degenerate_flag = degeneracy$degenerate_flag[[1]]
    ),
    tibble(
      dataset_key = model_row$dataset[[1]],
      branch = "Stage8B",
      model_id = model_row$model_id[[1]],
      retained_fit_id = model_row$retained_fit_id[[1]],
      risk_scale = model_row$risk_scale[[1]],
      prior_branch = model_row$prior_branch[[1]],
      site_prior_family = model_row$site_prior_family[[1]],
      metric = "cure_fraction",
      horizon = NA_integer_,
      summary_scalar(rowMeans(1 - state$pi_mat)),
      prior_degenerate_flag = degeneracy$degenerate_flag[[1]]
    ),
    tibble(
      dataset_key = model_row$dataset[[1]],
      branch = "Stage8B",
      model_id = model_row$model_id[[1]],
      retained_fit_id = model_row$retained_fit_id[[1]],
      risk_scale = model_row$risk_scale[[1]],
      prior_branch = model_row$prior_branch[[1]],
      site_prior_family = model_row$site_prior_family[[1]],
      metric = "median_susceptible_time",
      horizon = NA_integer_,
      summary_scalar(apply(state$median_mat, 1, stats::median)),
      prior_degenerate_flag = degeneracy$degenerate_flag[[1]]
    )
  )
  
  for (hh in sort(unique(as.integer(horizons_eval)))) {
    pred_h <- pred[[as.character(hh)]]
    out_rows[[length(out_rows) + 1L]] <- tibble(
      dataset_key = model_row$dataset[[1]],
      branch = "Stage8B",
      model_id = model_row$model_id[[1]],
      retained_fit_id = model_row$retained_fit_id[[1]],
      risk_scale = model_row$risk_scale[[1]],
      prior_branch = model_row$prior_branch[[1]],
      site_prior_family = model_row$site_prior_family[[1]],
      metric = "transition_cif",
      horizon = hh,
      summary_scalar(rowMeans(pred_h$transition_cif)),
      prior_degenerate_flag = degeneracy$degenerate_flag[[1]]
    )
    out_rows[[length(out_rows) + 1L]] <- tibble(
      dataset_key = model_row$dataset[[1]],
      branch = "Stage8B",
      model_id = model_row$model_id[[1]],
      retained_fit_id = model_row$retained_fit_id[[1]],
      risk_scale = model_row$risk_scale[[1]],
      prior_branch = model_row$prior_branch[[1]],
      site_prior_family = model_row$site_prior_family[[1]],
      metric = "remission_cif",
      horizon = hh,
      summary_scalar(rowMeans(pred_h$remission_cif)),
      prior_degenerate_flag = degeneracy$degenerate_flag[[1]]
    )
    out_rows[[length(out_rows) + 1L]] <- tibble(
      dataset_key = model_row$dataset[[1]],
      branch = "Stage8B",
      model_id = model_row$model_id[[1]],
      retained_fit_id = model_row$retained_fit_id[[1]],
      risk_scale = model_row$risk_scale[[1]],
      prior_branch = model_row$prior_branch[[1]],
      site_prior_family = model_row$site_prior_family[[1]],
      metric = "all_event_free",
      horizon = hh,
      summary_scalar(rowMeans(pred_h$all_event_free)),
      prior_degenerate_flag = degeneracy$degenerate_flag[[1]]
    )
  }
  
  bind_rows(out_rows)
}

annotate_stage8b_prior_predictive <- function(df) {
  if (nrow_or_zero(df) == 0L) {
    return(empty_stage8b_prior_predictive_summary() %>%
             mutate(
               prior_tail_warning_flag = logical(),
               prior_tail_warning_detail = character()
             ))
  }
  
  tibble::as_tibble(df) %>%
    mutate(
      horizon = as.integer(safe_numeric(horizon)),
      mean = safe_numeric(mean),
      q025 = safe_numeric(q025),
      q50 = safe_numeric(q50),
      q975 = safe_numeric(q975),
      prior_degenerate_flag = as.logical(prior_degenerate_flag),
      prior_tail_warning_flag = metric == "median_susceptible_time" &
        ((mean > prior_tail_warning_mean_years) | (q975 > prior_tail_warning_q975_years)),
      prior_tail_warning_detail = case_when(
        metric == "median_susceptible_time" & mean > prior_tail_warning_mean_years & q975 > prior_tail_warning_q975_years ~
          paste0("Prior median susceptible time mean > ", prior_tail_warning_mean_years, " years and q975 > ", prior_tail_warning_q975_years, " years."),
        metric == "median_susceptible_time" & mean > prior_tail_warning_mean_years ~
          paste0("Prior median susceptible time mean > ", prior_tail_warning_mean_years, " years."),
        metric == "median_susceptible_time" & q975 > prior_tail_warning_q975_years ~
          paste0("Prior median susceptible time q975 > ", prior_tail_warning_q975_years, " years."),
        TRUE ~ NA_character_
      )
    )
}

make_stage8b_anchor_cells <- function(model_row) {
  age_band_df <- tibble::tribble(
    ~sex_num, ~age20_29, ~age30plus, ~cell_label,
    0L, 0L, 0L, "Female_lt20",
    1L, 0L, 0L, "Male_lt20",
    0L, 1L, 0L, "Female_20_29",
    1L, 1L, 0L, "Male_20_29",
    0L, 0L, 1L, "Female_30plus",
    1L, 0L, 1L, "Male_30plus"
  )
  
  site_levels <- if (model_row$dataset[[1]] == "merged" && isTRUE(model_row$incidence_site_indicator[[1]])) {
    c("PNU", "SNU")
  } else if (model_row$dataset[[1]] == "merged") {
    "pooled"
  } else {
    model_row$dataset[[1]]
  }
  
  tidyr::crossing(age_band_df, tibble(site_level = site_levels)) %>%
    mutate(
      sex_x_age20_29 = sex_num * age20_29,
      sex_x_age30plus = sex_num * age30plus,
      site_SNU = as.integer(site_level == "SNU"),
      age_sex_anchor_cell = ifelse(site_level %in% c("pooled", "PNU", "SNU"),
                                   paste0(cell_label, "__", site_level),
                                   cell_label)
    )
}

make_stage8b_incidence_anchor_update <- function(model_row, design_bundle, prior_spec, draws_compact) {
  if (!identical(model_row$prior_branch[[1]], "anchor_informed")) {
    return(tibble())
  }
  
  cells <- make_stage8b_anchor_cells(model_row)
  out_rows <- vector("list", nrow(cells))
  
  for (ii in seq_len(nrow(cells))) {
    one <- cells[ii, , drop = FALSE]
    x_vec <- c(
      sex_num = one$sex_num[[1]],
      age20_29 = one$age20_29[[1]],
      age30plus = one$age30plus[[1]],
      sex_x_age20_29 = one$sex_x_age20_29[[1]],
      sex_x_age30plus = one$sex_x_age30plus[[1]]
    )
    if (isTRUE(model_row$incidence_site_indicator[[1]])) {
      x_vec <- c(x_vec, site_SNU = one$site_SNU[[1]])
    }
    
    prior_center_logit <- design_bundle$alpha_prior_center + sum(design_bundle$mu_beta_inc * x_vec)
    external_one_year_risk <- plogis(prior_center_logit)
    external_incidence_rate_per10k <- -10000 * log(pmax(1 - external_one_year_risk, 1e-12))
    
    post_logit <- design_bundle$alpha_prior_center + draws_compact$delta0 + as.vector(draws_compact$beta_inc %*% x_vec)
    post_risk <- plogis(post_logit)
    
    out_rows[[ii]] <- tibble(
      dataset_key = model_row$dataset[[1]],
      branch = "Stage8B",
      risk_scale = model_row$risk_scale[[1]],
      prior_branch = model_row$prior_branch[[1]],
      retained_fit_id = model_row$retained_fit_id[[1]],
      model_id = model_row$model_id[[1]],
      site_prior_family = model_row$site_prior_family[[1]],
      age_sex_anchor_cell = one$age_sex_anchor_cell[[1]],
      site_level = one$site_level[[1]],
      external_incidence_rate_per10k = external_incidence_rate_per10k,
      external_one_year_risk = external_one_year_risk,
      prior_center_logit = prior_center_logit,
      posterior_mean_logit = mean(post_logit),
      posterior_lower_logit = stats::quantile(post_logit, 0.025, names = FALSE),
      posterior_upper_logit = stats::quantile(post_logit, 0.975, names = FALSE),
      posterior_mean_one_year_risk = mean(post_risk),
      posterior_lower_one_year_risk = stats::quantile(post_risk, 0.025, names = FALSE),
      posterior_upper_one_year_risk = stats::quantile(post_risk, 0.975, names = FALSE),
      posterior_minus_prior_logit = mean(post_logit) - prior_center_logit,
      posterior_minus_prior_risk = mean(post_risk) - external_one_year_risk
    )
  }
  
  bind_rows(out_rows)
}

# 🔴 Define: Stan compilation, draw extraction, and diagnostics ===============================
## 🟠 Define: Stan model, compact draws, and information criteria ===============================
compile_stage8b_stan_model <- function() {
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
  real rem_base_hazard(real t, vector cuts, vector log_rho_rem) {
    int J;
    real out;
    J = num_elements(log_rho_rem);
    out = 0;
    for (j in 1:J) {
      real lower;
      real upper;
      real exposure;
      lower = cuts[j];
      upper = cuts[j + 1];
      exposure = fmax(fmin(t, upper) - lower, 0);
      out += exp(log_rho_rem[j]) * exposure;
    }
    return out;
  }
  int rem_interval_index(real t, vector cuts) {
    int J;
    J = num_elements(cuts) - 1;
    for (j in 1:J) {
      if (t >= cuts[j] && t < cuts[j + 1]) {
        return j;
      }
    }
    return J;
  }
}
data {
  int<lower=1> N;
  vector<lower=1e-8>[N] time;
  int<lower=0, upper=2> status[N];
  int<lower=1> K_inc;
  matrix[N, K_inc] X_inc;
  int<lower=1> K_lat;
  matrix[N, K_lat] X_lat;
  int<lower=1> K_rem;
  matrix[N, K_rem] X_rem;
  int<lower=1, upper=4> family_id;
  real alpha_prior_center;
  vector[K_inc] mu_beta_inc;
  vector<lower=0>[K_inc] sd_beta_inc;
  real<lower=0> sd_delta;
  real<lower=0> sd_gamma0;
  vector<lower=0>[K_lat] sd_gamma_lat;
  real<lower=0> sd_xi0;
  vector<lower=0>[K_rem] sd_beta_rem;
  real<lower=0> sd_log_rho_rem;
  real<lower=0> sd_rho_W;
  real<lower=0> sd_log_sigma_LN;
  real<lower=0> sd_psi_LL;
  int<lower=1> J_rem;
  vector[J_rem + 1] rem_cuts;
  real<lower=0> sd_site_prior;
  int<lower=0> inc_site_index;
  int<lower=0> lat_site_index;
  int<lower=0> rem_site_index;
  int<lower=0, upper=1> use_t_prior_inc_site;
  int<lower=0, upper=1> use_t_prior_lat_site;
  int<lower=0, upper=1> use_t_prior_rem_site;
}
parameters {
  real delta0;
  vector[K_inc] beta_inc;
  real gamma0;
  vector[K_lat] gamma_lat;
  real xi0;
  vector[K_rem] beta_rem;
  vector[J_rem] log_rho_rem;
  real rho_W;
  real log_sigma_LN;
  real psi_LL;
}
model {
  delta0 ~ normal(0, sd_delta);
  gamma0 ~ normal(0, sd_gamma0);
  xi0 ~ normal(0, sd_xi0);
  log_rho_rem ~ normal(0, sd_log_rho_rem);
  rho_W ~ normal(0, sd_rho_W);
  log_sigma_LN ~ normal(0, sd_log_sigma_LN);
  psi_LL ~ normal(0, sd_psi_LL);

  for (j in 1:K_inc) {
    if (use_t_prior_inc_site == 1 && inc_site_index > 0 && j == inc_site_index) {
      target += student_t_lpdf(beta_inc[j] | 3, 0, sd_site_prior);
    } else {
      target += normal_lpdf(beta_inc[j] | mu_beta_inc[j], sd_beta_inc[j]);
    }
  }

  for (j in 1:K_lat) {
    if (use_t_prior_lat_site == 1 && lat_site_index > 0 && j == lat_site_index) {
      target += student_t_lpdf(gamma_lat[j] | 3, 0, sd_site_prior);
    } else {
      target += normal_lpdf(gamma_lat[j] | 0, sd_gamma_lat[j]);
    }
  }

  for (j in 1:K_rem) {
    if (use_t_prior_rem_site == 1 && rem_site_index > 0 && j == rem_site_index) {
      target += student_t_lpdf(beta_rem[j] | 3, 0, sd_site_prior);
    } else {
      target += normal_lpdf(beta_rem[j] | 0, sd_beta_rem[j]);
    }
  }

  for (i in 1:N) {
    real eta_inc;
    real pi_i;
    real mu_lat;
    real m_i;
    real logS1;
    real logf1;
    real eta_rem;
    real H2_base;
    real logG2;
    real h2;
    real logM1;

    eta_inc = alpha_prior_center + delta0 + dot_product(row(X_inc, i), beta_inc);
    pi_i = inv_logit(eta_inc);
    mu_lat = gamma0 + dot_product(row(X_lat, i), gamma_lat);
    m_i = exp(mu_lat);
    eta_rem = xi0 + dot_product(row(X_rem, i), beta_rem);
    H2_base = rem_base_hazard(time[i], rem_cuts, log_rho_rem);
    logG2 = -exp(eta_rem) * H2_base;
    h2 = exp(eta_rem + log_rho_rem[rem_interval_index(time[i], rem_cuts)]);

    if (family_id == 1) {
      real lambda;
      lambda = m_i / log(2);
      logS1 = -time[i] / lambda;
      logf1 = exponential_lpdf(time[i] | 1 / lambda);
    } else if (family_id == 2) {
      real kW;
      real lambda;
      kW = exp(rho_W);
      lambda = m_i / pow(log(2), 1 / kW);
      logS1 = -pow(time[i] / lambda, kW);
      logf1 = weibull_lpdf(time[i] | kW, lambda);
    } else if (family_id == 3) {
      real sigmaLN;
      sigmaLN = exp(log_sigma_LN);
      logS1 = normal_lccdf(log(time[i]) | mu_lat, sigmaLN);
      logf1 = lognormal_lpdf(time[i] | mu_lat, sigmaLN);
    } else {
      real kLL;
      real lambda;
      kLL = exp(-psi_LL);
      lambda = m_i;
      logS1 = loglogistic_lccdf_custom(time[i], lambda, kLL);
      logf1 = loglogistic_lpdf_custom(time[i], lambda, kLL);
    }

    logM1 = log_sum_exp(log1m(pi_i), log(pi_i) + logS1);

    if (status[i] == 1) {
      target += log(pi_i) + logf1 + logG2;
    } else if (status[i] == 2) {
      target += log(h2) + logG2 + logM1;
    } else {
      target += logG2 + logM1;
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
    real logS1;
    real logf1;
    real eta_rem;
    real H2_base;
    real logG2;
    real h2;
    real logM1;

    eta_inc = alpha_prior_center + delta0 + dot_product(row(X_inc, i), beta_inc);
    pi_i = inv_logit(eta_inc);
    mu_lat = gamma0 + dot_product(row(X_lat, i), gamma_lat);
    m_i = exp(mu_lat);
    eta_rem = xi0 + dot_product(row(X_rem, i), beta_rem);
    H2_base = rem_base_hazard(time[i], rem_cuts, log_rho_rem);
    logG2 = -exp(eta_rem) * H2_base;
    h2 = exp(eta_rem + log_rho_rem[rem_interval_index(time[i], rem_cuts)]);

    if (family_id == 1) {
      real lambda;
      lambda = m_i / log(2);
      logS1 = -time[i] / lambda;
      logf1 = exponential_lpdf(time[i] | 1 / lambda);
    } else if (family_id == 2) {
      real kW;
      real lambda;
      kW = exp(rho_W);
      lambda = m_i / pow(log(2), 1 / kW);
      logS1 = -pow(time[i] / lambda, kW);
      logf1 = weibull_lpdf(time[i] | kW, lambda);
    } else if (family_id == 3) {
      real sigmaLN;
      sigmaLN = exp(log_sigma_LN);
      logS1 = normal_lccdf(log(time[i]) | mu_lat, sigmaLN);
      logf1 = lognormal_lpdf(time[i] | mu_lat, sigmaLN);
    } else {
      real kLL;
      real lambda;
      kLL = exp(-psi_LL);
      lambda = m_i;
      logS1 = loglogistic_lccdf_custom(time[i], lambda, kLL);
      logf1 = loglogistic_lpdf_custom(time[i], lambda, kLL);
    }

    logM1 = log_sum_exp(log1m(pi_i), log(pi_i) + logS1);

    if (status[i] == 1) {
      log_lik[i] = log(pi_i) + logf1 + logG2;
    } else if (status[i] == 2) {
      log_lik[i] = log(h2) + logG2 + logM1;
    } else {
      log_lik[i] = logG2 + logM1;
    }
  }
}
)"
rstan::stan_model(model_code = stan_code, model_name = "stage8b_bayesian_competing_risk_cure")
}

extract_stage8b_draws_compact <- function(fit, K_inc, K_lat, K_rem, J_rem) {
  ext <- rstan::extract(
    fit,
    pars = c("delta0", "beta_inc", "gamma0", "gamma_lat", "xi0", "beta_rem", "log_rho_rem", "rho_W", "log_sigma_LN", "psi_LL", "log_lik"),
    permuted = TRUE,
    inc_warmup = FALSE
  )
  
  beta_inc <- ext$beta_inc
  gamma_lat <- ext$gamma_lat
  beta_rem <- ext$beta_rem
  log_rho_rem <- ext$log_rho_rem
  log_lik <- ext$log_lik
  
  if (is.null(dim(beta_inc))) beta_inc <- matrix(beta_inc, ncol = K_inc)
  if (is.null(dim(gamma_lat))) gamma_lat <- matrix(gamma_lat, ncol = K_lat)
  if (is.null(dim(beta_rem))) beta_rem <- matrix(beta_rem, ncol = K_rem)
  if (is.null(dim(log_rho_rem))) log_rho_rem <- matrix(log_rho_rem, ncol = J_rem)
  if (is.null(dim(log_lik))) log_lik <- matrix(log_lik, nrow = length(ext$delta0))
  
  list(
    delta0 = as.numeric(ext$delta0),
    beta_inc = beta_inc,
    gamma0 = as.numeric(ext$gamma0),
    gamma_lat = gamma_lat,
    xi0 = as.numeric(ext$xi0),
    beta_rem = beta_rem,
    log_rho_rem = log_rho_rem,
    rho_W = as.numeric(ext$rho_W),
    log_sigma_LN = as.numeric(ext$log_sigma_LN),
    psi_LL = as.numeric(ext$psi_LL),
    log_lik = log_lik
  )
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
  
  waic_obj <- withCallingHandlers(
    tryCatch(loo::waic(log_lik), error = function(e) e),
    warning = collect_warning
  )
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
  
  loo_obj <- withCallingHandlers(
    tryCatch(loo::loo(log_lik), error = function(e) e),
    warning = collect_warning
  )
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
  if (length(warning_texts) > 0L) {
    out$info_criteria_warning_detail <- paste(unique(warning_texts), collapse = " | ")
  }
  out
}

select_stage8b_trace_parameters <- function(family_code, K_inc, K_lat, K_rem, J_rem) {
  out <- c("delta0", "gamma0", "xi0")
  if (K_inc >= 1) out <- c(out, "beta_inc[1]")
  if (K_lat >= 1) out <- c(out, "gamma_lat[1]")
  if (K_rem >= 1) out <- c(out, "beta_rem[1]")
  if (J_rem >= 1) out <- c(out, "log_rho_rem[1]")
  if (family_code == "W") out <- c(out, "rho_W")
  if (family_code == "LN") out <- c(out, "log_sigma_LN")
  if (family_code == "LL") out <- c(out, "psi_LL")
  unique(out[seq_len(min(length(out), 6L))])
}

make_stage8b_trace_record <- function(fit, model_id, family_code, K_inc, K_lat, K_rem, J_rem) {
  pars <- select_stage8b_trace_parameters(family_code, K_inc, K_lat, K_rem, J_rem)
  list(
    model_id = model_id,
    selected_pars = pars,
    arr = as.array(fit, pars = pars)
  )
}

plot_stage8b_trace_record <- function(trace_record) {
  selected_pars <- trace_record$selected_pars
  arr <- trace_record$arr
  par_old <- par(no.readonly = TRUE)
  on.exit(par(par_old), add = TRUE)
  par(mfrow = c(2, 2), mar = c(3, 3, 3, 1))
  for (jj in seq_len(min(length(selected_pars), 4L))) {
    matplot(
      arr[, , jj],
      type = "l",
      lty = 1,
      col = seq_len(dim(arr)[2]),
      main = paste(trace_record$model_id, selected_pars[[jj]]),
      xlab = "Iteration",
      ylab = ""
    )
  }
  if (length(selected_pars) < 4L) {
    for (jj in seq_len(4L - length(selected_pars))) plot.new()
  }
}

safe_generate_stage8b_diagnostic_pdf <- function(
    trace_records,
    posterior_cohort_yearly,
    posterior_classification,
    ppc_summary,
    final_path
) {
  dir.create(dirname(final_path), recursive = TRUE, showWarnings = FALSE)
  tmp_pdf <- make_temp_output_path(final_path, tag = "tmp")
  pdf_open <- FALSE
  
  close_pdf <- function() {
    if (isTRUE(pdf_open)) {
      try(grDevices::dev.off(), silent = TRUE)
      pdf_open <<- FALSE
    }
  }
  
  on.exit({
    close_pdf()
    if (file.exists(tmp_pdf)) unlink(tmp_pdf)
  }, add = TRUE)
  
  grDevices::pdf(tmp_pdf, width = 11, height = 8.5, onefile = TRUE)
  pdf_open <- TRUE
  
  if (length(trace_records) > 0L) {
    for (tr in trace_records) {
      plot_stage8b_trace_record(tr)
    }
  }
  
  if (nrow_or_zero(posterior_cohort_yearly) > 0L) {
    g1 <- posterior_cohort_yearly %>%
      ggplot(aes(x = horizon, y = transition_cif_mean, color = model_id, fill = model_id)) +
      geom_ribbon(aes(ymin = transition_cif_q025, ymax = transition_cif_q975), alpha = 0.15, linewidth = 0) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~ dataset_key, scales = "free_y") +
      labs(title = "Stage 8B posterior transition CIF trajectories", x = "Horizon (years)", y = "Transition CIF") +
      theme_bw() +
      theme(legend.position = "none")
    print(g1)
    
    g2 <- posterior_cohort_yearly %>%
      ggplot(aes(x = horizon, y = remission_cif_mean, color = model_id, fill = model_id)) +
      geom_ribbon(aes(ymin = remission_cif_q025, ymax = remission_cif_q975), alpha = 0.15, linewidth = 0) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~ dataset_key, scales = "free_y") +
      labs(title = "Stage 8B posterior remission CIF trajectories", x = "Horizon (years)", y = "Remission CIF") +
      theme_bw() +
      theme(legend.position = "none")
    print(g2)
  }
  
  if (nrow_or_zero(posterior_classification) > 0L) {
    nb_df <- posterior_classification %>% filter(horizon %in% c(1L, 2L, 5L))
    if (nrow(nb_df) > 0L) {
      g_nb <- nb_df %>%
        ggplot(aes(x = threshold, y = NB_mean, color = model_id, fill = model_id)) +
        geom_ribbon(aes(ymin = NB_q025, ymax = NB_q975), alpha = 0.15, linewidth = 0) +
        geom_line(linewidth = 0.7) +
        facet_grid(dataset_key ~ horizon, scales = "free_y") +
        labs(title = "Stage 8B net benefit by threshold", x = "Threshold", y = "Net benefit") +
        theme_bw() +
        theme(legend.position = "none")
      print(g_nb)
    }
  }
  
  if (nrow_or_zero(ppc_summary) > 0L) {
    g_ppc <- ppc_summary %>%
      ggplot(aes(x = horizon, y = posterior_mean_transition_cif, color = model_id)) +
      geom_errorbar(aes(ymin = posterior_q025_transition_cif, ymax = posterior_q975_transition_cif), width = 0.12, alpha = 0.6) +
      geom_line(linewidth = 0.6) +
      geom_point(linewidth = 0.6) +
      geom_point(aes(y = observed_transition_cif), shape = 4, size = 2.0, stroke = 0.9, color = "black") +
      facet_wrap(~ dataset_key, scales = "free_y") +
      labs(title = "Stage 8B posterior predictive checks vs observed transition CIF", x = "Horizon (years)", y = "Transition CIF") +
      theme_bw() +
      theme(legend.position = "none")
    print(g_ppc)
  }
  
  if (length(trace_records) == 0L && nrow_or_zero(posterior_cohort_yearly) == 0L && nrow_or_zero(posterior_classification) == 0L && nrow_or_zero(ppc_summary) == 0L) {
    plot.new()
    text(0.5, 0.5, "No diagnostic pages were available for this Stage 8B run.")
  }
  
  close_pdf()
  if (!pdf_file_is_usable(tmp_pdf)) {
    stop("Temporary Stage8B diagnostic PDF was not created correctly.", call. = FALSE)
  }
  safe_promote_file(tmp_pdf, final_path)
  invisible(TRUE)
}

load_stage8b_reuse_bundle <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  tryCatch(readRDS(path), error = function(e) NULL)
}

is_stage8b_bundle_reusable <- function(bundle, model_row, dataset_df) {
  if (is.null(bundle) || !is.list(bundle) || is.null(bundle$model_registry_row)) {
    return(FALSE)
  }
  reg <- tibble::as_tibble(bundle$model_registry_row)
  if (nrow(reg) != 1L) {
    return(FALSE)
  }
  identical(as.character(reg$model_id[[1]]), as.character(model_row$model_id[[1]])) &&
    identical(as.character(reg$dataset_key[[1]]), as.character(model_row$dataset[[1]])) &&
    identical(as.character(reg$fit_status[[1]]), "ok") &&
    identical(as.integer(reg$n[[1]]), nrow(dataset_df))
}


# 🔴 Load: Stage1 backbone, Stage6 carry-forward, and Stage8A comparators ===============================
## 🟠 Define: backbone objects, common grids, and comparator registries ===============================
stage1_backbone <- load_stage1_backbone(
  bundle_file = stage1_bundle_file,
  datasets_file = stage1_analysis_datasets_file,
  dataset_registry_file = stage1_dataset_registry_file,
  horizon_registry_file = stage1_horizon_registry_file,
  threshold_registry_file = stage1_threshold_registry_file,
  scaling_registry_file = stage1_scaling_registry_file
)

analysis_datasets <- stage1_backbone$datasets
dataset_registry_stage1 <- tibble::as_tibble(stage1_backbone$registries$dataset_registry)
scaling_registry <- tibble::as_tibble(stage1_backbone$registries$scaling_registry)
horizon_registry_stage1 <- build_horizon_annotation_from_stage1(
  horizon_registry = stage1_backbone$registries$horizon_registry,
  datasets = analysis_datasets,
  horizons = 1:10
)
thresholds_from_stage1 <- build_threshold_vector_from_stage1(stage1_backbone$registries$threshold_registry)

horizons_year <- sort(unique(as.integer(horizon_registry_stage1$horizon_year)))
risk_thresholds <- sort(unique(as.numeric(thresholds_from_stage1)))

if (!identical(horizons_year, 1:10)) {
  stop("Stage8B requires Stage1 horizon grid `1:10`.", call. = FALSE)
}
if (length(risk_thresholds) == 0L || anyNA(risk_thresholds) || any(risk_thresholds <= 0 | risk_thresholds >= 1)) {
  stop("Stage8B risk thresholds must be probabilities strictly between 0 and 1.", call. = FALSE)
}

screening_lookup_stage6 <- read_stage6_screening_lookup(stage6_screening_flag_csv)

stage8a_outputs_raw <- load_stage8a_outputs(
  model_registry_csv = stage8a_model_registry_csv,
  cohort_yearly_csv = stage8a_cohort_yearly_csv,
  classification_csv = stage8a_classification_csv
)
stage8a_outputs_aug <- augment_stage8a_with_model_keys(stage8a_outputs_raw)

dataset_registry_stage8b <- bind_rows(lapply(names(analysis_datasets), function(ds) {
  dat <- analysis_datasets[[ds]]
  tibble(
    dataset_key = ds,
    n = nrow(dat),
    n_transition = sum(dat$status_num == 1L),
    n_remission = sum(dat$status_num == 2L),
    n_right_censoring = sum(dat$status_num == 0L),
    person_time_years = sum(dat$time_year),
    median_followup_years = stats::median(dat$time_year),
    max_followup_years = max(dat$time_year)
  )
}))

aj_registry <- lapply(names(analysis_datasets), function(ds) {
  build_aalen_johansen_reference(analysis_datasets[[ds]], horizons_year) %>%
    mutate(dataset_key = ds)
})
names(aj_registry) <- names(analysis_datasets)

ipcw_registry <- lapply(names(analysis_datasets), function(ds) {
  build_stage8b_ipcw_reference(analysis_datasets[[ds]], horizons_year) %>%
    mutate(dataset_key = ds)
})
names(ipcw_registry) <- names(analysis_datasets)

support_registry <- make_stage8b_support_registry(horizon_registry_stage1, datasets = names(analysis_datasets), horizons = horizons_year)

# 🔴 Prepare: model grid, prior objects, and optional reuse plan ===============================
## 🟠 Define: grid expansion, screening joins, and Stan setup ===============================
prior_specs <- build_stage8b_prior_specs()

model_grid <- build_stage8b_model_grid(
  include_merged_incidence_site_supplementary = include_merged_incidence_site_supplementary,
  fit_prior_branches = fit_prior_branches,
  fit_site_prior_families = fit_site_prior_families
)

if (!is.null(run_model_ids)) {
  model_grid <- model_grid %>% filter(model_id %in% run_model_ids)
  if (nrow(model_grid) == 0L) {
    stop("`run_model_ids` filtered out all Stage8B models.", call. = FALSE)
  }
}

stage6_model_lookup <- build_stage6_model_lookup(screening_lookup_stage6, model_grid)
model_grid <- model_grid %>%
  left_join(stage6_model_lookup, by = "model_id")

model_grid$stage8b_rds_path <- file.path(export_path, paste0(model_grid$model_id, "__bayes_stage8b_fit.rds"))
model_grid$reuse_existing_fit <- FALSE

if (isTRUE(reuse_existing_stage8b_rds)) {
  for (ii in seq_len(nrow(model_grid))) {
    bundle_i <- load_stage8b_reuse_bundle(model_grid$stage8b_rds_path[[ii]])
    model_grid$reuse_existing_fit[[ii]] <- is_stage8b_bundle_reusable(bundle_i, model_grid[ii, , drop = FALSE], analysis_datasets[[model_grid$dataset[[ii]]]])
  }
}

n_models_reused <- sum(model_grid$reuse_existing_fit)
n_models_to_fit <- sum(!model_grid$reuse_existing_fit)

message(
  "Stage8B reuse plan: ",
  n_models_reused,
  " model(s) reused; ",
  n_models_to_fit,
  " model(s) require fitting."
)

if (n_models_to_fit > 0L) {
  fit_packages <- c("rstan", "posterior", "loo")
  missing_fit_packages <- fit_packages[!vapply(fit_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_fit_packages) > 0L) {
    stop(
      "Install required fitting packages before running this script: ",
      paste(missing_fit_packages, collapse = ", "),
      call. = FALSE
    )
  }
  
  suppressPackageStartupMessages({
    library(rstan)
    library(posterior)
  })
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = max(1L, min(stan_chains, parallel::detectCores(logical = TRUE))))
  stan_model_compiled <- compile_stage8b_stan_model()
} else {
  stan_model_compiled <- NULL
}

diagnostic_pdf_path <- file.path(export_path, "bayes_stage8b_diagnostic_plots.pdf")
trace_records <- list()

# 🔴 Run: Stage8B model loop with reuse-first execution ===============================
## 🟠 Define: model fitting, posterior summaries, and reusable RDS bundles ===============================
registry_rows <- list()
coef_rows <- list()
diag_param_rows <- list()
ppc_rows <- list()
subject_profile_rows <- list()
subject_yearly_rows <- list()
cohort_rows <- list()
classification_rows <- list()
prior_predictive_rows <- list()
hazard_shape_rows <- list()
uncured_support_rows <- list()
anchor_update_rows <- list()

withCallingHandlers({
  for (ii in seq_len(nrow(model_grid))) {
    model_row <- model_grid[ii, , drop = FALSE]
    model_id_now <- model_row$model_id[[1]]
    dataset_name <- model_row$dataset[[1]]
    dataset_df <- analysis_datasets[[dataset_name]]
    model_started_at <- Sys.time()
    
    emit_progress(
      ii - 1L,
      nrow(model_grid),
      model_id_now,
      paste0("starting Stage8B (dataset=", dataset_name, ", family=", model_row$family_code[[1]], ", prior=", model_row$prior_branch[[1]], ", site-prior=", model_row$site_prior_family[[1]], ", reuse=", isTRUE(model_row$reuse_existing_fit[[1]]), ")")
    )
    
    prior_spec <- prior_specs[[model_row$prior_branch[[1]]]]
    design_bundle <- make_stage8b_design_bundle(dataset_df, model_row, prior_spec, snu_site_label, remission_cut_years)
    
    set.seed(stan_seed + ii)
    prior_predictive_tbl <- simulate_stage8b_prior_predictive(
      dataset_df = dataset_df,
      model_row = model_row,
      design_bundle = design_bundle,
      prior_spec = prior_spec,
      n_draws = prior_predictive_draws,
      horizons_eval = prior_predictive_horizons
    )
    prior_predictive_rows[[length(prior_predictive_rows) + 1L]] <- prior_predictive_tbl
    
    if (isTRUE(model_row$reuse_existing_fit[[1]])) {
      bundle_reuse <- load_stage8b_reuse_bundle(model_row$stage8b_rds_path[[1]])
      registry_rows[[length(registry_rows) + 1L]] <- tibble::as_tibble(bundle_reuse$model_registry_row)
      coef_rows[[length(coef_rows) + 1L]] <- tibble::as_tibble(bundle_reuse$coefficient_summary)
      diag_param_rows[[length(diag_param_rows) + 1L]] <- tibble::as_tibble(bundle_reuse$diagnostics_parameter_level)
      ppc_rows[[length(ppc_rows) + 1L]] <- tibble::as_tibble(bundle_reuse$ppc_summary)
      subject_profile_rows[[length(subject_profile_rows) + 1L]] <- tibble::as_tibble(bundle_reuse$posterior_subject_profile)
      subject_yearly_rows[[length(subject_yearly_rows) + 1L]] <- tibble::as_tibble(bundle_reuse$posterior_subject_yearly)
      cohort_rows[[length(cohort_rows) + 1L]] <- tibble::as_tibble(bundle_reuse$posterior_cohort_yearly)
      classification_rows[[length(classification_rows) + 1L]] <- tibble::as_tibble(bundle_reuse$posterior_classification)
      hazard_shape_rows[[length(hazard_shape_rows) + 1L]] <- tibble::as_tibble(bundle_reuse$hazard_shape_plausibility)
      uncured_support_rows[[length(uncured_support_rows) + 1L]] <- tibble::as_tibble(bundle_reuse$uncured_supporting_decomposition)
      anchor_update_rows[[length(anchor_update_rows) + 1L]] <- tibble::as_tibble(bundle_reuse$incidence_anchor_update)
      
      emit_progress(
        ii,
        nrow(model_grid),
        model_id_now,
        paste0("reused existing Stage8B fit; elapsed=", format_number(elapsed_seconds(model_started_at), digits = 1L), "s")
      )
      next
    }
    
    site_idx <- get_stage8b_site_prior_indices(model_row, design_bundle)
    
    stan_data <- list(
      N = nrow(dataset_df),
      time = as.numeric(design_bundle$time),
      status = as.integer(design_bundle$status),
      K_inc = ncol(design_bundle$X_inc),
      X_inc = design_bundle$X_inc,
      K_lat = ncol(design_bundle$X_lat),
      X_lat = design_bundle$X_lat,
      K_rem = ncol(design_bundle$X_rem),
      X_rem = design_bundle$X_rem,
      family_id = as.integer(model_row$family_id[[1]]),
      alpha_prior_center = as.numeric(design_bundle$alpha_prior_center),
      mu_beta_inc = as.numeric(design_bundle$mu_beta_inc),
      sd_beta_inc = as.numeric(design_bundle$sd_beta_inc),
      sd_delta = as.numeric(prior_spec$sd_delta),
      sd_gamma0 = as.numeric(prior_spec$sd_gamma0),
      sd_gamma_lat = as.numeric(design_bundle$sd_gamma_lat),
      sd_xi0 = as.numeric(prior_spec$sd_xi0),
      sd_beta_rem = as.numeric(design_bundle$sd_beta_rem),
      sd_log_rho_rem = as.numeric(prior_spec$sd_log_rho_rem),
      sd_rho_W = as.numeric(prior_spec$sd_rho_W),
      sd_log_sigma_LN = as.numeric(prior_spec$sd_log_sigma_LN),
      sd_psi_LL = as.numeric(prior_spec$sd_psi_LL),
      J_rem = length(design_bundle$rem_cuts) - 1L,
      rem_cuts = as.numeric(design_bundle$rem_cuts),
      sd_site_prior = as.numeric(prior_spec$sd_site),
      inc_site_index = ifelse(is.na(site_idx$inc_idx), 0L, as.integer(site_idx$inc_idx)),
      lat_site_index = ifelse(is.na(site_idx$lat_idx), 0L, as.integer(site_idx$lat_idx)),
      rem_site_index = ifelse(is.na(site_idx$rem_idx), 0L, as.integer(site_idx$rem_idx)),
      use_t_prior_inc_site = as.integer(!is.na(site_idx$inc_idx) && identical(model_row$site_prior_family[[1]], "student_t3_0_1_sensitivity")),
      use_t_prior_lat_site = as.integer(!is.na(site_idx$lat_idx) && identical(model_row$site_prior_family[[1]], "student_t3_0_1_sensitivity")),
      use_t_prior_rem_site = as.integer(!is.na(site_idx$rem_idx) && identical(model_row$site_prior_family[[1]], "student_t3_0_1_sensitivity"))
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
      
      registry_rows[[length(registry_rows) + 1L]] <- tibble(
        dataset_key = dataset_name,
        model_id = model_id_now,
        retained_fit_id = model_row$retained_fit_id[[1]],
        structural_model_id = model_row$structural_model_id[[1]],
        formula_anchor = model_row$formula_anchor[[1]],
        transition_latency_branch = model_row$transition_latency_branch[[1]],
        remission_branch = model_row$remission_branch[[1]],
        site_placement_label = model_row$site_placement_label[[1]],
        branch = "Stage8B",
        risk_scale = model_row$risk_scale[[1]],
        prior_branch = model_row$prior_branch[[1]],
        site_prior_family = model_row$site_prior_family[[1]],
        latency_family = model_row$latency_family[[1]],
        family_code = model_row$family_code[[1]],
        fit_status = fit_status,
        fit_error_message = fit_error_message,
        admissibility_flag = FALSE,
        admissibility_reasons = "sampling_error",
        prior_degenerate_flag = any(prior_predictive_tbl$prior_degenerate_flag %in% TRUE),
        posterior_degenerate_flag = NA,
        ppc_gross_contradiction_flag = NA,
        coherence_violation_flag = NA,
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
        n_transition = sum(dataset_df$status_num == 1L),
        n_remission = sum(dataset_df$status_num == 2L),
        n_right_censoring = sum(dataset_df$status_num == 0L),
        cure_model_eligibility_flag = model_row$cure_model_eligibility_flag[[1]],
        receus_primary_class = model_row$receus_primary_class[[1]],
        cure_presence_support_flag = model_row$cure_presence_support_flag[[1]],
        followup_not_contradicted_flag = model_row$followup_not_contradicted_flag[[1]],
        screening_note = model_row$screening_note[[1]],
        rds_path = NA_character_,
        fit_reused_flag = FALSE
      )
      
      emit_progress(
        ii,
        nrow(model_grid),
        model_id_now,
        paste0("sampling error after ", format_number(elapsed_seconds(model_started_at), digits = 1L), "s: ", fit_error_message)
      )
      next
    }
    
    trace_records[[length(trace_records) + 1L]] <- make_stage8b_trace_record(
      fit = fit,
      model_id = model_id_now,
      family_code = model_row$family_code[[1]],
      K_inc = ncol(design_bundle$X_inc),
      K_lat = ncol(design_bundle$X_lat),
      K_rem = ncol(design_bundle$X_rem),
      J_rem = length(design_bundle$rem_cuts) - 1L
    )
    
    param_names <- c(
      "delta0", "gamma0", "xi0",
      paste0("beta_inc[", seq_len(ncol(design_bundle$X_inc)), "]"),
      paste0("gamma_lat[", seq_len(ncol(design_bundle$X_lat)), "]"),
      paste0("beta_rem[", seq_len(ncol(design_bundle$X_rem)), "]"),
      paste0("log_rho_rem[", seq_len(length(design_bundle$rem_cuts) - 1L), "]")
    )
    if (model_row$family_code[[1]] == "W") param_names <- c(param_names, "rho_W")
    if (model_row$family_code[[1]] == "LN") param_names <- c(param_names, "log_sigma_LN")
    if (model_row$family_code[[1]] == "LL") param_names <- c(param_names, "psi_LL")
    
    param_array <- posterior::as_draws_array(as.array(fit, pars = param_names))
    param_diag_tbl <- posterior::summarise_draws(
      param_array,
      mean = base::mean,
      sd = stats::sd,
      rhat = posterior::rhat,
      ess_bulk = posterior::ess_bulk,
      ess_tail = posterior::ess_tail
    )
    param_draws_mat <- posterior::as_draws_matrix(param_array)
    
    coef_tbl <- tibble(
      dataset_key = dataset_name,
      model_id = model_id_now,
      retained_fit_id = model_row$retained_fit_id[[1]],
      parameter = colnames(param_draws_mat),
      mean = apply(param_draws_mat, 2, mean),
      sd = apply(param_draws_mat, 2, stats::sd),
      q025 = apply(param_draws_mat, 2, stats::quantile, probs = 0.025, names = FALSE),
      q50 = apply(param_draws_mat, 2, stats::quantile, probs = 0.500, names = FALSE),
      q975 = apply(param_draws_mat, 2, stats::quantile, probs = 0.975, names = FALSE)
    )
    
    diag_param_tbl_model <- tibble(
      dataset_key = dataset_name,
      model_id = model_id_now,
      retained_fit_id = model_row$retained_fit_id[[1]],
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
    
    draws_compact <- extract_stage8b_draws_compact(
      fit = fit,
      K_inc = ncol(design_bundle$X_inc),
      K_lat = ncol(design_bundle$X_lat),
      K_rem = ncol(design_bundle$X_rem),
      J_rem = length(design_bundle$rem_cuts) - 1L
    )
    
    info_criteria <- compute_information_criteria(draws_compact$log_lik)
    
    total_draws <- length(draws_compact$delta0)
    set.seed(stan_seed + 200000L + ii)
    keep_draw_idx <- if (total_draws <= posterior_prediction_draws) {
      seq_len(total_draws)
    } else {
      sort(sample(seq_len(total_draws), size = posterior_prediction_draws, replace = FALSE))
    }
    
    draws_pred <- list(
      delta0 = draws_compact$delta0[keep_draw_idx],
      beta_inc = draws_compact$beta_inc[keep_draw_idx, , drop = FALSE],
      gamma0 = draws_compact$gamma0[keep_draw_idx],
      gamma_lat = draws_compact$gamma_lat[keep_draw_idx, , drop = FALSE],
      xi0 = draws_compact$xi0[keep_draw_idx],
      beta_rem = draws_compact$beta_rem[keep_draw_idx, , drop = FALSE],
      log_rho_rem = draws_compact$log_rho_rem[keep_draw_idx, , drop = FALSE],
      rho_W = draws_compact$rho_W[keep_draw_idx],
      log_sigma_LN = draws_compact$log_sigma_LN[keep_draw_idx],
      psi_LL = draws_compact$psi_LL[keep_draw_idx]
    )
    
    state <- compute_stage8b_linear_terms(
      draws = draws_pred,
      X_inc = design_bundle$X_inc,
      X_lat = design_bundle$X_lat,
      X_rem = design_bundle$X_rem,
      alpha_prior_center = design_bundle$alpha_prior_center
    )
    
    pred_trajectories <- build_stage8b_prediction_trajectories(
      state = state,
      draws = draws_pred,
      model_row = model_row,
      horizons = horizons_year,
      integration_step_year = integration_step_year,
      rem_cuts = design_bundle$rem_cuts
    )
    
    supported_horizons <- ppc_horizons_for_dataset(dataset_name)
    transition_supported <- lapply(as.character(supported_horizons), function(hh) pred_trajectories[[hh]]$transition_cif)
    posterior_degeneracy <- compute_stage8b_degeneracy(state$pi_mat, state$median_mat, transition_supported)
    
    subject_profile_tbl <- bind_cols(
      tibble(
        dataset_key = dataset_name,
        model_id = model_id_now,
        retained_fit_id = model_row$retained_fit_id[[1]],
        branch = "Stage8B",
        risk_scale = model_row$risk_scale[[1]],
        prior_branch = model_row$prior_branch[[1]],
        site_prior_family = model_row$site_prior_family[[1]]
      ),
      design_bundle$id_df,
      tibble(
        cure_fraction_mean = summarize_cols_matrix(state$cure_prob_mat)$mean,
        cure_fraction_q025 = summarize_cols_matrix(state$cure_prob_mat)$q025,
        cure_fraction_q50 = summarize_cols_matrix(state$cure_prob_mat)$q50,
        cure_fraction_q975 = summarize_cols_matrix(state$cure_prob_mat)$q975,
        susceptible_fraction_mean = summarize_cols_matrix(state$pi_mat)$mean,
        susceptible_fraction_q025 = summarize_cols_matrix(state$pi_mat)$q025,
        susceptible_fraction_q50 = summarize_cols_matrix(state$pi_mat)$q50,
        susceptible_fraction_q975 = summarize_cols_matrix(state$pi_mat)$q975,
        median_susceptible_time_mean = summarize_cols_matrix(state$median_mat)$mean,
        median_susceptible_time_q025 = summarize_cols_matrix(state$median_mat)$q025,
        median_susceptible_time_q50 = summarize_cols_matrix(state$median_mat)$q50,
        median_susceptible_time_q975 = summarize_cols_matrix(state$median_mat)$q975
      )
    )
    
    horizon_refs <- ipcw_registry[[dataset_name]]
    support_refs <- support_registry %>% filter(dataset == dataset_name)
    
    subject_year_rows_model <- list()
    cohort_rows_model <- list()
    class_rows_model <- list()
    ppc_rows_model <- list()
    
    coherence_error_max <- 0
    hazard_means <- rep(NA_real_, length(horizons_year))
    names(hazard_means) <- paste0("hazard_", horizons_year, "y")
    
    for (hh in horizons_year) {
      pred_h <- pred_trajectories[[as.character(hh)]]
      horizon_ref <- horizon_refs %>% filter(horizon_year == hh)
      support_ref <- support_refs %>% filter(horizon_year == hh)
      
      trans_sum <- summarize_cols_matrix(pred_h$transition_cif)
      rem_sum <- summarize_cols_matrix(pred_h$remission_cif)
      free_sum <- summarize_cols_matrix(pred_h$all_event_free)
      uncured_surv_sum <- summarize_cols_matrix(pred_h$uncured_survival)
      uncured_risk_sum <- summarize_cols_matrix(1 - pred_h$uncured_survival)
      
      coherence_error_max <- max(
        coherence_error_max,
        max(abs(pred_h$transition_cif + pred_h$remission_cif + pred_h$all_event_free - 1), na.rm = TRUE)
      )
      
      subject_year_rows_model[[length(subject_year_rows_model) + 1L]] <- bind_cols(
        tibble(
          dataset_key = dataset_name,
          model_id = model_id_now,
          retained_fit_id = model_row$retained_fit_id[[1]],
          branch = "Stage8B",
          risk_scale = model_row$risk_scale[[1]],
          prior_branch = model_row$prior_branch[[1]],
          site_prior_family = model_row$site_prior_family[[1]],
          horizon = hh
        ),
        design_bundle$id_df,
        tibble(
          transition_cif_mean = trans_sum$mean,
          transition_cif_q025 = trans_sum$q025,
          transition_cif_q50 = trans_sum$q50,
          transition_cif_q975 = trans_sum$q975,
          remission_cif_mean = rem_sum$mean,
          remission_cif_q025 = rem_sum$q025,
          remission_cif_q50 = rem_sum$q50,
          remission_cif_q975 = rem_sum$q975,
          all_event_free_mean = free_sum$mean,
          all_event_free_q025 = free_sum$q025,
          all_event_free_q50 = free_sum$q50,
          all_event_free_q975 = free_sum$q975,
          uncured_survival_mean = uncured_surv_sum$mean,
          uncured_survival_q025 = uncured_surv_sum$q025,
          uncured_survival_q50 = uncured_surv_sum$q50,
          uncured_survival_q975 = uncured_surv_sum$q975,
          uncured_risk_mean = uncured_risk_sum$mean,
          uncured_risk_q025 = uncured_risk_sum$q025,
          uncured_risk_q50 = uncured_risk_sum$q50,
          uncured_risk_q975 = uncured_risk_sum$q975
        )
      )
      
      transition_cif_draw <- rowMeans(pred_h$transition_cif)
      remission_cif_draw <- rowMeans(pred_h$remission_cif)
      free_draw <- rowMeans(pred_h$all_event_free)
      hazard_draw <- rowMeans(pred_h$transition_population_hazard)
      cure_draw <- rowMeans(state$cure_prob_mat)
      susc_draw <- rowMeans(state$pi_mat)
      uncured_surv_draw <- rowMeans(pred_h$uncured_survival)
      uncured_risk_draw <- rowMeans(1 - pred_h$uncured_survival)
      
      hazard_means[[paste0("hazard_", hh, "y")]] <- mean(hazard_draw, na.rm = TRUE)
      
      cohort_rows_model[[length(cohort_rows_model) + 1L]] <- tibble(
        dataset_key = dataset_name,
        model_id = model_id_now,
        retained_fit_id = model_row$retained_fit_id[[1]],
        structural_model_id = model_row$structural_model_id[[1]],
        formula_anchor = model_row$formula_anchor[[1]],
        family_code = model_row$family_code[[1]],
        latency_family = model_row$latency_family[[1]],
        site_placement_label = model_row$site_placement_label[[1]],
        branch = "Stage8B",
        risk_scale = model_row$risk_scale[[1]],
        prior_branch = model_row$prior_branch[[1]],
        site_prior_family = model_row$site_prior_family[[1]],
        horizon = hh,
        transition_cif_mean = mean(transition_cif_draw),
        transition_cif_q025 = stats::quantile(transition_cif_draw, 0.025, names = FALSE),
        transition_cif_q50 = stats::quantile(transition_cif_draw, 0.500, names = FALSE),
        transition_cif_q975 = stats::quantile(transition_cif_draw, 0.975, names = FALSE),
        remission_cif_mean = mean(remission_cif_draw),
        remission_cif_q025 = stats::quantile(remission_cif_draw, 0.025, names = FALSE),
        remission_cif_q50 = stats::quantile(remission_cif_draw, 0.500, names = FALSE),
        remission_cif_q975 = stats::quantile(remission_cif_draw, 0.975, names = FALSE),
        all_event_free_mean = mean(free_draw),
        all_event_free_q025 = stats::quantile(free_draw, 0.025, names = FALSE),
        all_event_free_q50 = stats::quantile(free_draw, 0.500, names = FALSE),
        all_event_free_q975 = stats::quantile(free_draw, 0.975, names = FALSE),
        cohort_mean_cure_fraction_mean = mean(cure_draw),
        cohort_mean_cure_fraction_q025 = stats::quantile(cure_draw, 0.025, names = FALSE),
        cohort_mean_cure_fraction_q50 = stats::quantile(cure_draw, 0.500, names = FALSE),
        cohort_mean_cure_fraction_q975 = stats::quantile(cure_draw, 0.975, names = FALSE),
        cohort_mean_susceptible_fraction_mean = mean(susc_draw),
        cohort_mean_susceptible_fraction_q025 = stats::quantile(susc_draw, 0.025, names = FALSE),
        cohort_mean_susceptible_fraction_q50 = stats::quantile(susc_draw, 0.500, names = FALSE),
        cohort_mean_susceptible_fraction_q975 = stats::quantile(susc_draw, 0.975, names = FALSE),
        mean_uncured_survival_mean = mean(uncured_surv_draw),
        mean_uncured_survival_q025 = stats::quantile(uncured_surv_draw, 0.025, names = FALSE),
        mean_uncured_survival_q50 = stats::quantile(uncured_surv_draw, 0.500, names = FALSE),
        mean_uncured_survival_q975 = stats::quantile(uncured_surv_draw, 0.975, names = FALSE),
        mean_uncured_risk_mean = mean(uncured_risk_draw),
        mean_uncured_risk_q025 = stats::quantile(uncured_risk_draw, 0.025, names = FALSE),
        mean_uncured_risk_q50 = stats::quantile(uncured_risk_draw, 0.500, names = FALSE),
        mean_uncured_risk_q975 = stats::quantile(uncured_risk_draw, 0.975, names = FALSE),
        mean_transition_hazard_mean = mean(hazard_draw),
        mean_transition_hazard_q025 = stats::quantile(hazard_draw, 0.025, names = FALSE),
        mean_transition_hazard_q50 = stats::quantile(hazard_draw, 0.500, names = FALSE),
        mean_transition_hazard_q975 = stats::quantile(hazard_draw, 0.975, names = FALSE)
      )
      
      ppc_rows_model[[length(ppc_rows_model) + 1L]] <- tibble(
        dataset_key = dataset_name,
        model_id = model_id_now,
        retained_fit_id = model_row$retained_fit_id[[1]],
        structural_model_id = model_row$structural_model_id[[1]],
        formula_anchor = model_row$formula_anchor[[1]],
        family_code = model_row$family_code[[1]],
        branch = "Stage8B",
        risk_scale = model_row$risk_scale[[1]],
        prior_branch = model_row$prior_branch[[1]],
        site_prior_family = model_row$site_prior_family[[1]],
        horizon = hh,
        observed_transition_cif = horizon_ref$observed_transition_cif[[1]],
        observed_remission_cif = horizon_ref$observed_remission_cif[[1]],
        observed_all_event_free = horizon_ref$observed_all_event_free[[1]],
        posterior_mean_transition_cif = mean(transition_cif_draw),
        posterior_q025_transition_cif = stats::quantile(transition_cif_draw, 0.025, names = FALSE),
        posterior_q975_transition_cif = stats::quantile(transition_cif_draw, 0.975, names = FALSE),
        posterior_mean_remission_cif = mean(remission_cif_draw),
        posterior_q025_remission_cif = stats::quantile(remission_cif_draw, 0.025, names = FALSE),
        posterior_q975_remission_cif = stats::quantile(remission_cif_draw, 0.975, names = FALSE),
        absolute_difference_transition_cif = abs(mean(transition_cif_draw) - horizon_ref$observed_transition_cif[[1]]),
        gross_contradiction_flag = (
          (hh %in% ppc_horizons_for_dataset(dataset_name)) &&
            (
              horizon_ref$observed_transition_cif[[1]] < stats::quantile(transition_cif_draw, 0.025, names = FALSE) ||
                horizon_ref$observed_transition_cif[[1]] > stats::quantile(transition_cif_draw, 0.975, names = FALSE)
            ) &&
            abs(mean(transition_cif_draw) - horizon_ref$observed_transition_cif[[1]]) > ppc_tolerance_abs
        )
      )
      
      class_tbl_h <- compute_stage8b_classification_summary(
        risk_draws = pred_h$transition_cif,
        horizon_row = horizon_ref,
        thresholds = risk_thresholds
      ) %>%
        mutate(
          dataset_key = dataset_name,
          model_id = model_id_now,
          retained_fit_id = model_row$retained_fit_id[[1]],
          structural_model_id = model_row$structural_model_id[[1]],
          formula_anchor = model_row$formula_anchor[[1]],
          family_code = model_row$family_code[[1]],
          latency_family = model_row$latency_family[[1]],
          branch = "Stage8B",
          risk_scale = model_row$risk_scale[[1]],
          prior_branch = model_row$prior_branch[[1]],
          site_prior_family = model_row$site_prior_family[[1]],
          horizon = hh,
          observed_transition_cif = horizon_ref$observed_transition_cif[[1]]
        ) %>%
        relocate(dataset_key, model_id, retained_fit_id, structural_model_id, formula_anchor, family_code, branch, risk_scale, prior_branch, site_prior_family, horizon, threshold)
      
      class_rows_model[[length(class_rows_model) + 1L]] <- class_tbl_h
    }
    
    ppc_model_tbl <- bind_rows(ppc_rows_model)
    subject_year_tbl <- bind_rows(subject_year_rows_model)
    cohort_model_tbl <- bind_rows(cohort_rows_model)
    class_model_tbl <- bind_rows(class_rows_model)
    
    ppc_gross_contradiction_flag <- any(ppc_model_tbl$gross_contradiction_flag, na.rm = TRUE)
    coherence_violation_flag <- is.finite(coherence_error_max) && coherence_error_max > 0.02
    
    max_rhat <- max(param_diag_tbl$rhat, na.rm = TRUE)
    min_bulk_ess <- min(param_diag_tbl$ess_bulk, na.rm = TRUE)
    min_tail_ess <- min(param_diag_tbl$ess_tail, na.rm = TRUE)
    
    admissibility_reasons <- character()
    if (any(prior_predictive_tbl$prior_degenerate_flag %in% TRUE)) admissibility_reasons <- c(admissibility_reasons, "prior_degenerate")
    if (posterior_degeneracy$degenerate_flag[[1]]) admissibility_reasons <- c(admissibility_reasons, "posterior_degenerate")
    if (!is.finite(max_rhat) || max_rhat >= rhat_max_threshold) admissibility_reasons <- c(admissibility_reasons, "rhat")
    if (!is.finite(min_bulk_ess) || min_bulk_ess < ess_min_threshold) admissibility_reasons <- c(admissibility_reasons, "bulk_ess")
    if (!is.finite(min_tail_ess) || min_tail_ess < ess_min_threshold) admissibility_reasons <- c(admissibility_reasons, "tail_ess")
    if (divergences > 0) admissibility_reasons <- c(admissibility_reasons, "divergences")
    if (treedepth_exceeded > 0) admissibility_reasons <- c(admissibility_reasons, "treedepth")
    if (isTRUE(ppc_gross_contradiction_flag)) admissibility_reasons <- c(admissibility_reasons, "ppc")
    if (isTRUE(coherence_violation_flag)) admissibility_reasons <- c(admissibility_reasons, "coherence")
    
    admissibility_flag <- length(admissibility_reasons) == 0L
    
    hazard_shape_tbl <- tibble(
      dataset_key = dataset_name,
      model_id = model_id_now,
      retained_fit_id = model_row$retained_fit_id[[1]],
      structural_model_id = model_row$structural_model_id[[1]],
      formula_anchor = model_row$formula_anchor[[1]],
      family_code = model_row$family_code[[1]],
      latency_family = model_row$latency_family[[1]],
      branch = "Stage8B",
      risk_scale = model_row$risk_scale[[1]],
      prior_branch = model_row$prior_branch[[1]],
      site_prior_family = model_row$site_prior_family[[1]],
      hazard_target = "transition_population_cause_specific_hazard",
      hazard_1y = hazard_means[["hazard_1y"]],
      hazard_2y = hazard_means[["hazard_2y"]],
      hazard_3y = hazard_means[["hazard_3y"]],
      hazard_4y = hazard_means[["hazard_4y"]],
      hazard_5y = hazard_means[["hazard_5y"]],
      hazard_6y = hazard_means[["hazard_6y"]],
      hazard_7y = hazard_means[["hazard_7y"]],
      hazard_8y = hazard_means[["hazard_8y"]],
      hazard_9y = hazard_means[["hazard_9y"]],
      hazard_10y = hazard_means[["hazard_10y"]],
      hazard_ratio_10y_vs_1y = safe_numeric(hazard_means[["hazard_10y"]]) / pmax(safe_numeric(hazard_means[["hazard_1y"]]), 1e-12),
      shape_class = classify_hazard_shape(unname(hazard_means))
    )
    
    cure_draw <- rowMeans(state$cure_prob_mat)
    susc_draw <- rowMeans(state$pi_mat)
    uncured_support_tbl <- bind_rows(lapply(horizons_year, function(hh) {
      pred_h <- pred_trajectories[[as.character(hh)]]
      uncured_surv_draw <- rowMeans(pred_h$uncured_survival)
      uncured_risk_draw <- rowMeans(1 - pred_h$uncured_survival)
      tibble(
        dataset_key = dataset_name,
        model_id = model_id_now,
        retained_fit_id = model_row$retained_fit_id[[1]],
        structural_model_id = model_row$structural_model_id[[1]],
        formula_anchor = model_row$formula_anchor[[1]],
        family_code = model_row$family_code[[1]],
        branch = "Stage8B",
        risk_scale = model_row$risk_scale[[1]],
        prior_branch = model_row$prior_branch[[1]],
        site_prior_family = model_row$site_prior_family[[1]],
        horizon = hh,
        cure_fraction_mean = mean(cure_draw),
        cure_fraction_q025 = stats::quantile(cure_draw, 0.025, names = FALSE),
        cure_fraction_q50 = stats::quantile(cure_draw, 0.500, names = FALSE),
        cure_fraction_q975 = stats::quantile(cure_draw, 0.975, names = FALSE),
        susceptible_fraction_mean = mean(susc_draw),
        susceptible_fraction_q025 = stats::quantile(susc_draw, 0.025, names = FALSE),
        susceptible_fraction_q50 = stats::quantile(susc_draw, 0.500, names = FALSE),
        susceptible_fraction_q975 = stats::quantile(susc_draw, 0.975, names = FALSE),
        uncured_survival_mean = mean(uncured_surv_draw),
        uncured_survival_q025 = stats::quantile(uncured_surv_draw, 0.025, names = FALSE),
        uncured_survival_q50 = stats::quantile(uncured_surv_draw, 0.500, names = FALSE),
        uncured_survival_q975 = stats::quantile(uncured_surv_draw, 0.975, names = FALSE),
        uncured_risk_mean = mean(uncured_risk_draw),
        uncured_risk_q025 = stats::quantile(uncured_risk_draw, 0.025, names = FALSE),
        uncured_risk_q50 = stats::quantile(uncured_risk_draw, 0.500, names = FALSE),
        uncured_risk_q975 = stats::quantile(uncured_risk_draw, 0.975, names = FALSE),
        MSTu_mean = NA_real_,
        MSTu_q025 = NA_real_,
        MSTu_q50 = NA_real_,
        MSTu_q975 = NA_real_,
        uncured_mean_support_flag = FALSE
      )
    }))
    
    anchor_update_tbl <- make_stage8b_incidence_anchor_update(model_row, design_bundle, prior_spec, draws_compact)
    
    rds_path_now <- model_row$stage8b_rds_path[[1]]
    registry_row <- tibble(
      dataset_key = dataset_name,
      model_id = model_id_now,
      retained_fit_id = model_row$retained_fit_id[[1]],
      structural_model_id = model_row$structural_model_id[[1]],
      formula_anchor = model_row$formula_anchor[[1]],
      transition_latency_branch = model_row$transition_latency_branch[[1]],
      remission_branch = model_row$remission_branch[[1]],
      site_placement_label = model_row$site_placement_label[[1]],
      branch = "Stage8B",
      risk_scale = model_row$risk_scale[[1]],
      prior_branch = model_row$prior_branch[[1]],
      site_prior_family = model_row$site_prior_family[[1]],
      latency_family = model_row$latency_family[[1]],
      family_code = model_row$family_code[[1]],
      fit_status = fit_status,
      fit_error_message = fit_error_message,
      admissibility_flag = admissibility_flag,
      admissibility_reasons = if (length(admissibility_reasons) == 0L) "" else paste(admissibility_reasons, collapse = "|"),
      prior_degenerate_flag = any(prior_predictive_tbl$prior_degenerate_flag %in% TRUE),
      posterior_degenerate_flag = posterior_degeneracy$degenerate_flag[[1]],
      ppc_gross_contradiction_flag = ppc_gross_contradiction_flag,
      coherence_violation_flag = coherence_violation_flag,
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
      n_transition = sum(dataset_df$status_num == 1L),
      n_remission = sum(dataset_df$status_num == 2L),
      n_right_censoring = sum(dataset_df$status_num == 0L),
      cure_model_eligibility_flag = model_row$cure_model_eligibility_flag[[1]],
      receus_primary_class = model_row$receus_primary_class[[1]],
      cure_presence_support_flag = model_row$cure_presence_support_flag[[1]],
      followup_not_contradicted_flag = model_row$followup_not_contradicted_flag[[1]],
      screening_note = model_row$screening_note[[1]],
      cohort_mean_cure_fraction_mean = mean(cure_draw),
      cohort_mean_cure_fraction_q025 = stats::quantile(cure_draw, 0.025, names = FALSE),
      cohort_mean_cure_fraction_q50 = stats::quantile(cure_draw, 0.500, names = FALSE),
      cohort_mean_cure_fraction_q975 = stats::quantile(cure_draw, 0.975, names = FALSE),
      rds_path = rds_path_now,
      fit_reused_flag = FALSE
    )
    
    save_obj <- list(
      version = "stage8b_v1",
      model_registry_row = registry_row,
      coefficient_summary = coef_tbl,
      diagnostics_parameter_level = diag_param_tbl_model,
      ppc_summary = ppc_model_tbl,
      posterior_subject_profile = if (isTRUE(admissibility_flag)) subject_profile_tbl else tibble(),
      posterior_subject_yearly = if (isTRUE(admissibility_flag)) subject_year_tbl else tibble(),
      posterior_cohort_yearly = if (isTRUE(admissibility_flag)) cohort_model_tbl else tibble(),
      posterior_classification = if (isTRUE(admissibility_flag)) class_model_tbl else tibble(),
      prior_predictive_summary = prior_predictive_tbl,
      hazard_shape_plausibility = if (isTRUE(admissibility_flag)) hazard_shape_tbl else tibble(),
      uncured_supporting_decomposition = if (isTRUE(admissibility_flag)) uncured_support_tbl else tibble(),
      incidence_anchor_update = if (isTRUE(admissibility_flag)) anchor_update_tbl else tibble()
    )
    if (isTRUE(save_full_stanfit_rds)) {
      save_obj$fit <- fit
    }
    saveRDS(save_obj, rds_path_now)
    
    registry_rows[[length(registry_rows) + 1L]] <- registry_row
    coef_rows[[length(coef_rows) + 1L]] <- coef_tbl
    diag_param_rows[[length(diag_param_rows) + 1L]] <- diag_param_tbl_model
    ppc_rows[[length(ppc_rows) + 1L]] <- ppc_model_tbl
    if (isTRUE(admissibility_flag)) {
      subject_profile_rows[[length(subject_profile_rows) + 1L]] <- subject_profile_tbl
      subject_yearly_rows[[length(subject_yearly_rows) + 1L]] <- subject_year_tbl
      cohort_rows[[length(cohort_rows) + 1L]] <- cohort_model_tbl
      classification_rows[[length(classification_rows) + 1L]] <- class_model_tbl
      hazard_shape_rows[[length(hazard_shape_rows) + 1L]] <- hazard_shape_tbl
      uncured_support_rows[[length(uncured_support_rows) + 1L]] <- uncured_support_tbl
      anchor_update_rows[[length(anchor_update_rows) + 1L]] <- anchor_update_tbl
    }
    
    rm(fit, param_array, param_draws_mat, draws_compact, draws_pred, state, pred_trajectories)
    gc(verbose = FALSE)
    
    emit_progress(
      ii,
      nrow(model_grid),
      model_id_now,
      paste0(
        "completed; elapsed=", format_number(elapsed_seconds(model_started_at), digits = 1L),
        "s; WAIC=", format_number(info_criteria$waic, digits = 2L),
        "; LOOIC=", format_number(info_criteria$looic, digits = 2L),
        "; Pareto k max=", format_number(info_criteria$pareto_k_max, digits = 3L)
      )
    )
  }
}, warning = function(w) {
  if (is_localhost_connection_warning(w)) {
    tryInvokeRestart("muffleWarning")
  }
})


# 🔴 Assemble: final Stage8B tables and governance metadata ===============================
## 🟠 Define: combined tables, support joins, and delta reconstructions ===============================
model_order <- model_grid$model_id

model_registry <- bind_rows_safe(registry_rows) %>%
  arrange(factor(model_id, levels = model_order))

coefficient_summary <- bind_rows_safe(coef_rows) %>%
  arrange(factor(model_id, levels = model_order), parameter)

diagnostics_parameter_level <- bind_rows_safe(diag_param_rows) %>%
  arrange(factor(model_id, levels = model_order), parameter)

ppc_summary <- bind_rows_safe(ppc_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon)

posterior_subject_profile <- bind_rows_safe(subject_profile_rows) %>%
  arrange(factor(model_id, levels = model_order), unique_person_id)

posterior_subject_yearly <- bind_rows_safe(subject_yearly_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon, unique_person_id)

posterior_cohort_yearly <- bind_rows_safe(cohort_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon)

posterior_classification <- bind_rows_safe(classification_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon, threshold)

hazard_shape_plausibility <- bind_rows_safe(hazard_shape_rows) %>%
  arrange(factor(model_id, levels = model_order))

uncured_supporting_decomposition <- bind_rows_safe(uncured_support_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon)

incidence_anchor_update <- bind_rows_safe(anchor_update_rows) %>%
  arrange(factor(model_id, levels = model_order), age_sex_anchor_cell)

prior_predictive_summary <- annotate_stage8b_prior_predictive(bind_rows_safe(prior_predictive_rows)) %>%
  arrange(factor(model_id, levels = model_order), metric, prior_branch, horizon)

model_annotation <- model_registry %>%
  select(
    model_id,
    retained_fit_id,
    admissibility_flag,
    admissibility_reasons,
    cure_model_eligibility_flag,
    receus_primary_class,
    cure_presence_support_flag,
    followup_not_contradicted_flag,
    screening_note,
    fit_reused_flag
  )

horizon_annotation <- support_registry %>%
  transmute(
    dataset_key = dataset,
    horizon = horizon_year,
    support_tier = support_tier,
    horizon_evidence_class = horizon_evidence_class,
    claim_restriction_flag = claim_restriction_flag,
    interpretation_note = interpretation_note
  )

posterior_subject_yearly <- posterior_subject_yearly %>%
  left_join_replacing_columns(horizon_annotation, by = c("dataset_key", "horizon")) %>%
  left_join_replacing_columns(model_annotation, by = c("model_id", "retained_fit_id"))

posterior_cohort_yearly <- posterior_cohort_yearly %>%
  left_join_replacing_columns(horizon_annotation, by = c("dataset_key", "horizon")) %>%
  left_join_replacing_columns(model_annotation, by = c("model_id", "retained_fit_id"))

posterior_classification <- posterior_classification %>%
  left_join_replacing_columns(horizon_annotation, by = c("dataset_key", "horizon")) %>%
  left_join_replacing_columns(model_annotation, by = c("model_id", "retained_fit_id"))

ppc_summary <- ppc_summary %>%
  left_join_replacing_columns(horizon_annotation, by = c("dataset_key", "horizon")) %>%
  left_join_replacing_columns(model_annotation, by = c("model_id", "retained_fit_id"))

uncured_supporting_decomposition <- uncured_supporting_decomposition %>%
  left_join_replacing_columns(horizon_annotation, by = c("dataset_key", "horizon")) %>%
  left_join_replacing_columns(model_annotation, by = c("model_id", "retained_fit_id"))

prior_warning_by_model <- prior_predictive_summary %>%
  group_by(model_id, retained_fit_id) %>%
  summarise(
    prior_tail_warning_flag = any(prior_tail_warning_flag %in% TRUE),
    prior_tail_warning_detail = paste(unique(na.omit(prior_tail_warning_detail[prior_tail_warning_flag %in% TRUE])), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(prior_tail_warning_detail = dplyr::na_if(prior_tail_warning_detail, ""))

model_registry <- model_registry %>%
  left_join(prior_warning_by_model, by = c("model_id", "retained_fit_id")) %>%
  mutate(prior_tail_warning_flag = coalesce(prior_tail_warning_flag, FALSE))

cohort_delta_long <- bind_rows(
  posterior_cohort_yearly %>%
    transmute(
      dataset_key,
      retained_fit_id,
      structural_model_id,
      formula_anchor,
      family_code,
      site_prior_family,
      model_id,
      prior_branch,
      horizon,
      threshold = NA_real_,
      metric = "transition_cif",
      estimate = transition_cif_mean
    ),
  posterior_cohort_yearly %>%
    transmute(
      dataset_key,
      retained_fit_id,
      structural_model_id,
      formula_anchor,
      family_code,
      site_prior_family,
      model_id,
      prior_branch,
      horizon,
      threshold = NA_real_,
      metric = "cure_fraction",
      estimate = cohort_mean_cure_fraction_mean
    )
)

class_delta_long <- bind_rows(
  posterior_classification %>%
    transmute(dataset_key, retained_fit_id, structural_model_id, formula_anchor, family_code, site_prior_family, model_id, prior_branch, horizon, threshold, metric = "false_positive_burden", estimate = false_positive_burden_mean),
  posterior_classification %>%
    transmute(dataset_key, retained_fit_id, structural_model_id, formula_anchor, family_code, site_prior_family, model_id, prior_branch, horizon, threshold, metric = "FP100", estimate = FP100_mean),
  posterior_classification %>%
    transmute(dataset_key, retained_fit_id, structural_model_id, formula_anchor, family_code, site_prior_family, model_id, prior_branch, horizon, threshold, metric = "NB", estimate = NB_mean),
  posterior_classification %>%
    transmute(dataset_key, retained_fit_id, structural_model_id, formula_anchor, family_code, site_prior_family, model_id, prior_branch, horizon, threshold, metric = "PPV", estimate = PPV_mean),
  posterior_classification %>%
    transmute(dataset_key, retained_fit_id, structural_model_id, formula_anchor, family_code, site_prior_family, model_id, prior_branch, horizon, threshold, metric = "TPR", estimate = TPR_mean)
)

anchor_vs_neutral_delta <- bind_rows(cohort_delta_long, class_delta_long) %>%
  mutate(
    threshold_key = ifelse(is.na(threshold), "__NA__", sprintf("%.10f", threshold))
  ) %>%
  select(-threshold) %>%
  tidyr::pivot_wider(
    names_from = prior_branch,
    values_from = c(model_id, estimate),
    names_sep = "__"
  ) %>%
  transmute(
    dataset_key = dataset_key,
    branch = "Stage8B",
    risk_scale = "transition_cif_competing",
    retained_fit_id = retained_fit_id,
    structural_model_id = structural_model_id,
    formula_anchor = formula_anchor,
    family_code = family_code,
    site_prior_family = site_prior_family,
    horizon = horizon,
    threshold = ifelse(threshold_key == "__NA__", NA_real_, safe_numeric(threshold_key)),
    metric = metric,
    model_id_anchor = model_id__anchor_informed,
    model_id_neutral = model_id__neutral_no_external_info,
    anchor_estimate = estimate__anchor_informed,
    neutral_estimate = estimate__neutral_no_external_info,
    delta_anchor_minus_neutral = estimate__anchor_informed - estimate__neutral_no_external_info
  ) %>%
  left_join_replacing_columns(horizon_annotation, by = c("dataset_key", "horizon")) %>%
  arrange(retained_fit_id, horizon, threshold, metric)

prior_tail_sensitive_lookup <- anchor_vs_neutral_delta %>%
  mutate(
    threshold_key = ifelse(is.na(threshold), "__NA__", sprintf("%.10f", threshold)),
    material_flag = case_when(
      metric %in% c("transition_cif", "cure_fraction", "PPV", "TPR") ~ abs(delta_anchor_minus_neutral) > prior_materiality_risk,
      metric == "false_positive_burden" ~ abs(delta_anchor_minus_neutral) > prior_materiality_false_positive_burden,
      metric == "FP100" ~ abs(delta_anchor_minus_neutral) > 100 * prior_materiality_false_positive_burden,
      metric == "NB" ~ abs(delta_anchor_minus_neutral) > prior_materiality_nb,
      TRUE ~ FALSE
    )
  ) %>%
  group_by(dataset_key, retained_fit_id, horizon, threshold_key) %>%
  summarise(prior_tail_sensitive = any(material_flag %in% TRUE), .groups = "drop")

apply_prior_tail_sensitive <- function(df, has_threshold = FALSE) {
  if (nrow_or_zero(df) == 0L) {
    return(df)
  }
  out <- df %>%
    mutate(threshold_key = if (has_threshold) ifelse(is.na(threshold), "__NA__", sprintf("%.10f", threshold)) else "__NA__") %>%
    left_join(prior_tail_sensitive_lookup, by = c("dataset_key", "retained_fit_id", "horizon", "threshold_key")) %>%
    left_join(prior_warning_by_model %>% select(model_id, retained_fit_id, prior_tail_warning_flag), by = c("model_id", "retained_fit_id")) %>%
    mutate(
      prior_tail_sensitive = coalesce(prior_tail_sensitive, FALSE) | coalesce(prior_tail_warning_flag, FALSE),
      claim_restriction_flag = ifelse(
        claim_restriction_flag == "projection_only" & prior_tail_sensitive,
        "projection_plus_prior_sensitive",
        claim_restriction_flag
      )
    ) %>%
    select(-threshold_key, -prior_tail_warning_flag)
  out
}

posterior_subject_yearly <- apply_prior_tail_sensitive(posterior_subject_yearly, has_threshold = FALSE)
posterior_cohort_yearly <- apply_prior_tail_sensitive(posterior_cohort_yearly, has_threshold = FALSE)
posterior_classification <- apply_prior_tail_sensitive(posterior_classification, has_threshold = TRUE)
ppc_summary <- apply_prior_tail_sensitive(ppc_summary, has_threshold = FALSE)
uncured_supporting_decomposition <- apply_prior_tail_sensitive(uncured_supporting_decomposition, has_threshold = FALSE)

stage8a_cohort_key <- tibble::as_tibble(stage8a_outputs_aug$posterior_cohort_yearly) %>%
  transmute(
    dataset_key = dataset,
    structural_model_id = structural_model_id,
    formula_anchor = formula_anchor,
    family_code = family_code,
    horizon = as.integer(horizon_year),
    stage8a_transition_risk = meanRisk_Bayes_mean,
    stage8a_cure_fraction = cohort_mean_cure_fraction_mean
  )

stage8a_class_key <- tibble::as_tibble(stage8a_outputs_aug$posterior_classification) %>%
  transmute(
    dataset_key = dataset,
    structural_model_id = structural_model_id,
    formula_anchor = formula_anchor,
    family_code = family_code,
    horizon = as.integer(horizon_year),
    threshold = as.numeric(safe_numeric(threshold)),
    stage8a_false_positive_burden = false_positive_burden_mean,
    stage8a_FP100 = FP100_mean,
    stage8a_NB = NB_mean,
    stage8a_PPV = PPV_mean,
    stage8a_TPR = TPR_mean
  )

delta_vs_stage8a <- bind_rows(
  posterior_cohort_yearly %>%
    left_join(stage8a_cohort_key, by = c("dataset_key", "structural_model_id", "formula_anchor", "family_code", "horizon")) %>%
    transmute(
      dataset_key,
      branch = "Stage8B",
      risk_scale = "transition_cif_competing",
      model_id,
      retained_fit_id,
      structural_model_id,
      formula_anchor,
      family_code,
      prior_branch,
      site_prior_family,
      horizon,
      threshold = NA_real_,
      metric = "transition_cif",
      stage8b_estimate = transition_cif_mean,
      stage8a_estimate = stage8a_transition_risk,
      delta_8B_minus_8A = transition_cif_mean - stage8a_transition_risk
    ),
  posterior_cohort_yearly %>%
    left_join(stage8a_cohort_key, by = c("dataset_key", "structural_model_id", "formula_anchor", "family_code", "horizon")) %>%
    transmute(
      dataset_key,
      branch = "Stage8B",
      risk_scale = "transition_cif_competing",
      model_id,
      retained_fit_id,
      structural_model_id,
      formula_anchor,
      family_code,
      prior_branch,
      site_prior_family,
      horizon,
      threshold = NA_real_,
      metric = "cure_fraction",
      stage8b_estimate = cohort_mean_cure_fraction_mean,
      stage8a_estimate = stage8a_cure_fraction,
      delta_8B_minus_8A = cohort_mean_cure_fraction_mean - stage8a_cure_fraction
    ),
  posterior_classification %>%
    left_join(stage8a_class_key, by = c("dataset_key", "structural_model_id", "formula_anchor", "family_code", "horizon", "threshold")) %>%
    transmute(
      dataset_key,
      branch = "Stage8B",
      risk_scale = "transition_cif_competing",
      model_id,
      retained_fit_id,
      structural_model_id,
      formula_anchor,
      family_code,
      prior_branch,
      site_prior_family,
      horizon,
      threshold,
      metric = "false_positive_burden",
      stage8b_estimate = false_positive_burden_mean,
      stage8a_estimate = stage8a_false_positive_burden,
      delta_8B_minus_8A = false_positive_burden_mean - stage8a_false_positive_burden
    ),
  posterior_classification %>%
    left_join(stage8a_class_key, by = c("dataset_key", "structural_model_id", "formula_anchor", "family_code", "horizon", "threshold")) %>%
    transmute(
      dataset_key,
      branch = "Stage8B",
      risk_scale = "transition_cif_competing",
      model_id,
      retained_fit_id,
      structural_model_id,
      formula_anchor,
      family_code,
      prior_branch,
      site_prior_family,
      horizon,
      threshold,
      metric = "FP100",
      stage8b_estimate = FP100_mean,
      stage8a_estimate = stage8a_FP100,
      delta_8B_minus_8A = FP100_mean - stage8a_FP100
    ),
  posterior_classification %>%
    left_join(stage8a_class_key, by = c("dataset_key", "structural_model_id", "formula_anchor", "family_code", "horizon", "threshold")) %>%
    transmute(
      dataset_key,
      branch = "Stage8B",
      risk_scale = "transition_cif_competing",
      model_id,
      retained_fit_id,
      structural_model_id,
      formula_anchor,
      family_code,
      prior_branch,
      site_prior_family,
      horizon,
      threshold,
      metric = "NB",
      stage8b_estimate = NB_mean,
      stage8a_estimate = stage8a_NB,
      delta_8B_minus_8A = NB_mean - stage8a_NB
    ),
  posterior_classification %>%
    left_join(stage8a_class_key, by = c("dataset_key", "structural_model_id", "formula_anchor", "family_code", "horizon", "threshold")) %>%
    transmute(
      dataset_key,
      branch = "Stage8B",
      risk_scale = "transition_cif_competing",
      model_id,
      retained_fit_id,
      structural_model_id,
      formula_anchor,
      family_code,
      prior_branch,
      site_prior_family,
      horizon,
      threshold,
      metric = "PPV",
      stage8b_estimate = PPV_mean,
      stage8a_estimate = stage8a_PPV,
      delta_8B_minus_8A = PPV_mean - stage8a_PPV
    ),
  posterior_classification %>%
    left_join(stage8a_class_key, by = c("dataset_key", "structural_model_id", "formula_anchor", "family_code", "horizon", "threshold")) %>%
    transmute(
      dataset_key,
      branch = "Stage8B",
      risk_scale = "transition_cif_competing",
      model_id,
      retained_fit_id,
      structural_model_id,
      formula_anchor,
      family_code,
      prior_branch,
      site_prior_family,
      horizon,
      threshold,
      metric = "TPR",
      stage8b_estimate = TPR_mean,
      stage8a_estimate = stage8a_TPR,
      delta_8B_minus_8A = TPR_mean - stage8a_TPR
    )
) %>%
  left_join_replacing_columns(horizon_annotation, by = c("dataset_key", "horizon")) %>%
  arrange(factor(model_id, levels = model_order), horizon, threshold, metric)

horizon_support_panel <- posterior_cohort_yearly %>%
  select(
    dataset_key,
    model_id,
    retained_fit_id,
    structural_model_id,
    formula_anchor,
    family_code,
    prior_branch,
    site_prior_family,
    horizon,
    transition_cif_mean,
    transition_cif_q025,
    transition_cif_q975,
    support_tier,
    horizon_evidence_class,
    claim_restriction_flag,
    prior_tail_sensitive,
    admissibility_flag
  ) %>%
  arrange(factor(model_id, levels = model_order), horizon)

# 🔴 Check: output consistency and export manifest ===============================
## 🟠 Define: audit summary, PDF generation, and CSV/RDS manifest ===============================
output_audit <- bind_rows(
  tibble(
    check_name = "model_registry_row_count",
    status = ifelse(nrow(model_registry) == nrow(model_grid), "pass", "fail"),
    observed_value = as.character(nrow(model_registry)),
    expected_value = as.character(nrow(model_grid)),
    detail = "One model-registry row should exist per Stage8B model."
  ),
  tibble(
    check_name = "posterior_cohort_yearly_row_count",
    status = ifelse(nrow(posterior_cohort_yearly) == sum(model_registry$admissibility_flag %in% TRUE) * length(horizons_year), "pass", "fail"),
    observed_value = as.character(nrow(posterior_cohort_yearly)),
    expected_value = as.character(sum(model_registry$admissibility_flag %in% TRUE) * length(horizons_year)),
    detail = "Admissible fits should contribute one cohort row per horizon."
  ),
  tibble(
    check_name = "posterior_classification_row_count",
    status = ifelse(nrow(posterior_classification) == sum(model_registry$admissibility_flag %in% TRUE) * length(horizons_year) * length(risk_thresholds), "pass", "fail"),
    observed_value = as.character(nrow(posterior_classification)),
    expected_value = as.character(sum(model_registry$admissibility_flag %in% TRUE) * length(horizons_year) * length(risk_thresholds)),
    detail = "Admissible fits should contribute one classification row per horizon-threshold pair."
  ),
  tibble(
    check_name = "anchor_vs_neutral_delta_nonempty",
    status = ifelse(nrow(anchor_vs_neutral_delta) > 0L, "pass", "warn"),
    observed_value = as.character(nrow(anchor_vs_neutral_delta)),
    expected_value = ">0",
    detail = "Anchor-versus-neutral delta table should be populated when both prior branches were fitted."
  ),
  tibble(
    check_name = "delta_vs_stage8a_nonempty",
    status = ifelse(nrow(delta_vs_stage8a) > 0L, "pass", "warn"),
    observed_value = as.character(nrow(delta_vs_stage8a)),
    expected_value = ">0",
    detail = "Stage8B-versus-Stage8A delta table is populated only when matching Stage8A outputs were found."
  )
)

pdf_ok <- tryCatch(
  {
    safe_generate_stage8b_diagnostic_pdf(
      trace_records = trace_records,
      posterior_cohort_yearly = posterior_cohort_yearly,
      posterior_classification = posterior_classification,
      ppc_summary = ppc_summary,
      final_path = diagnostic_pdf_path
    )
    TRUE
  },
  error = function(e) {
    warning(
      paste0("Stage8B diagnostic PDF generation failed: ", conditionMessage(e)),
      call. = FALSE
    )
    FALSE
  }
)

output_audit <- bind_rows(
  output_audit,
  tibble(
    check_name = "diagnostic_pdf_exists",
    status = ifelse(pdf_file_is_usable(diagnostic_pdf_path), "pass", "fail"),
    observed_value = ifelse(pdf_file_is_usable(diagnostic_pdf_path), "TRUE", "FALSE"),
    expected_value = "TRUE",
    detail = "Diagnostic PDF should exist after Stage8B export."
  )
)

metadata_registry <- tibble::tribble(
  ~metadata_group, ~metadata_name, ~metadata_value,
  "stage", "stage_name", "Stage 8B remission-sensitive Bayesian competing-risk cure extension",
  "stage", "branch", "Stage8B",
  "stage", "risk_scale", "transition_cif_competing",
  "stage", "main_stage8_branch_reference", "Stage8A transition-only cure branch remains primary",
  "inputs", "stage1_export_path", stage1_export_path,
  "inputs", "stage6_screening_flag_csv", stage6_screening_flag_csv,
  "inputs", "stage8a_export_path", stage8a_export_path,
  "thresholds", "common_horizons_year", paste(horizons_year, collapse = ","),
  "thresholds", "risk_thresholds", paste(format(risk_thresholds, trim = TRUE, scientific = FALSE), collapse = ","),
  "model", "transition_families", paste(c("E", "W", "LN", "LL"), collapse = "|"),
  "model", "prior_branches", paste(unique(model_grid$prior_branch), collapse = "|"),
  "model", "site_prior_families", paste(sort(unique(model_grid$site_prior_family)), collapse = "|"),
  "model", "remission_cut_years", paste(remission_cut_years, collapse = "|"),
  "model", "remission_model", "piecewise_exponential",
  "stan", "stan_chains", as.character(stan_chains),
  "stan", "stan_iter", as.character(stan_iter),
  "stan", "stan_warmup", as.character(stan_warmup),
  "stan", "posterior_prediction_draws", as.character(posterior_prediction_draws),
  "interpretation", "cure_fraction_meaning", "non_susceptibility_to_transition",
  "interpretation", "site_effect_rule", "site_terms_are_structural_context_proxies_not_causal_treatment_effects",
  "interpretation", "cross_scale_rule", "Stage8B does not replace Stage8A; cross-scale deltas are remission-aware change summaries only"
)

export_manifest <- tibble(
  file_name = c(
    "bayes_stage8b_model_registry.csv",
    "bayes_stage8b_coefficient_summary.csv",
    "bayes_stage8b_diagnostics_parameter_level.csv",
    "bayes_stage8b_ppc_summary.csv",
    "bayes_stage8b_posterior_subject_profile.csv.gz",
    "bayes_stage8b_posterior_subject_yearly.csv.gz",
    "bayes_stage8b_posterior_cohort_yearly.csv",
    "bayes_stage8b_posterior_classification.csv",
    "bayes_stage8b_prior_predictive_summary.csv",
    "bayes_stage8b_hazard_plausibility.csv",
    "bayes_stage8b_uncured_decomposition.csv",
    "bayes_stage8b_anchor_vs_neutral_delta.csv",
    "bayes_stage8b_incidence_anchor_update.csv",
    "bayes_stage8b_stage8a_vs_stage8b_delta.csv",
    "bayes_stage8b_horizon_support_panel.csv",
    "bayes_stage8b_output_audit.csv",
    "bayes_stage8b_metadata_registry.csv",
    "bayes_stage8b_diagnostic_plots.pdf",
    "bayes_stage8b_bundle.rds",
    "bayes_stage8b_export_manifest.csv"
  ),
  object_name = c(
    "model_registry",
    "coefficient_summary",
    "diagnostics_parameter_level",
    "ppc_summary",
    "posterior_subject_profile",
    "posterior_subject_yearly",
    "posterior_cohort_yearly",
    "posterior_classification",
    "prior_predictive_summary",
    "hazard_shape_plausibility",
    "uncured_supporting_decomposition",
    "anchor_vs_neutral_delta",
    "incidence_anchor_update",
    "delta_vs_stage8a",
    "horizon_support_panel",
    "output_audit",
    "metadata_registry",
    "diagnostic_pdf",
    "stage8b_bundle",
    "export_manifest"
  ),
  description = c(
    "Stage8B model-level registry with admissibility, diagnostics, and carry-forward screening fields",
    "Parameter posterior summaries for all fitted Stage8B models",
    "Parameter-level convergence diagnostics for Stage8B fits",
    "Posterior predictive checks against observed transition/remission CIF quantities",
    "Subject-level posterior cure/susceptibility profile table",
    "Subject-by-horizon posterior transition/remission/all-event-free table",
    "Cohort-by-horizon posterior Stage8B summaries on the transition CIF scale",
    "Threshold-based Stage8B classification and decision-curve summaries using transition CIF",
    "Prior predictive summaries for Stage8B fits",
    "Hazard-shape plausibility table on the 1-10 year annual grid",
    "Cure-model-only supporting decomposition table with uncured-only annual summaries",
    "Mandatory anchor-informed versus neutral prior delta table",
    "Mandatory prior-to-posterior incidence-shape update source table",
    "Mandatory remission-aware Stage8B versus Stage8A delta table",
    "Figure-ready source table for the horizon-support panel",
    "Automated Stage8B output audit",
    "Stage8B metadata registry aligned to the revised integrated master specification",
    "Diagnostic PDF generated from exported Stage8B summary tables",
    "Reusable Stage8B bundle with key outputs and configuration",
    "Manifest of Stage8B exported files"
  ),
  file_path = file.path(export_path, c(
    "bayes_stage8b_model_registry.csv",
    "bayes_stage8b_coefficient_summary.csv",
    "bayes_stage8b_diagnostics_parameter_level.csv",
    "bayes_stage8b_ppc_summary.csv",
    "bayes_stage8b_posterior_subject_profile.csv.gz",
    "bayes_stage8b_posterior_subject_yearly.csv.gz",
    "bayes_stage8b_posterior_cohort_yearly.csv",
    "bayes_stage8b_posterior_classification.csv",
    "bayes_stage8b_prior_predictive_summary.csv",
    "bayes_stage8b_hazard_plausibility.csv",
    "bayes_stage8b_uncured_decomposition.csv",
    "bayes_stage8b_anchor_vs_neutral_delta.csv",
    "bayes_stage8b_incidence_anchor_update.csv",
    "bayes_stage8b_stage8a_vs_stage8b_delta.csv",
    "bayes_stage8b_horizon_support_panel.csv",
    "bayes_stage8b_output_audit.csv",
    "bayes_stage8b_metadata_registry.csv",
    "bayes_stage8b_diagnostic_plots.pdf",
    "bayes_stage8b_bundle.rds",
    "bayes_stage8b_export_manifest.csv"
  ))
)

# 🔴 Export: Stage8B CSV outputs and final manifest ===============================
## 🟠 Define: write tables, save manifest, and finish run ===============================
write_csv_preserve_schema(model_registry, file.path(export_path, "bayes_stage8b_model_registry.csv"))
write_csv_preserve_schema(coefficient_summary, file.path(export_path, "bayes_stage8b_coefficient_summary.csv"))
write_csv_preserve_schema(diagnostics_parameter_level, file.path(export_path, "bayes_stage8b_diagnostics_parameter_level.csv"))
write_csv_preserve_schema(ppc_summary, file.path(export_path, "bayes_stage8b_ppc_summary.csv"))
write_csv_preserve_schema(posterior_subject_profile, file.path(export_path, "bayes_stage8b_posterior_subject_profile.csv.gz"))
write_csv_preserve_schema(posterior_subject_yearly, file.path(export_path, "bayes_stage8b_posterior_subject_yearly.csv.gz"))
write_csv_preserve_schema(posterior_cohort_yearly, file.path(export_path, "bayes_stage8b_posterior_cohort_yearly.csv"))
write_csv_preserve_schema(posterior_classification, file.path(export_path, "bayes_stage8b_posterior_classification.csv"))
write_csv_preserve_schema(prior_predictive_summary, file.path(export_path, "bayes_stage8b_prior_predictive_summary.csv"))
write_csv_preserve_schema(hazard_shape_plausibility, file.path(export_path, "bayes_stage8b_hazard_plausibility.csv"))
write_csv_preserve_schema(uncured_supporting_decomposition, file.path(export_path, "bayes_stage8b_uncured_decomposition.csv"))
write_csv_preserve_schema(anchor_vs_neutral_delta, file.path(export_path, "bayes_stage8b_anchor_vs_neutral_delta.csv"))
write_csv_preserve_schema(incidence_anchor_update, file.path(export_path, "bayes_stage8b_incidence_anchor_update.csv"))
write_csv_preserve_schema(delta_vs_stage8a, file.path(export_path, "bayes_stage8b_stage8a_vs_stage8b_delta.csv"))
write_csv_preserve_schema(horizon_support_panel, file.path(export_path, "bayes_stage8b_horizon_support_panel.csv"))
write_csv_preserve_schema(output_audit, file.path(export_path, "bayes_stage8b_output_audit.csv"))
write_csv_preserve_schema(metadata_registry, file.path(export_path, "bayes_stage8b_metadata_registry.csv"))
saveRDS(
  list(
    stage = "Stage8B",
    created_at = as.character(Sys.time()),
    session_info = utils::sessionInfo(),
    config = list(
      stage1_export_path = stage1_export_path,
      export_path = export_path,
      stage6_export_path = stage6_export_path,
      stage8a_export_path = stage8a_export_path,
      pnu_site_label = pnu_site_label,
      snu_site_label = snu_site_label,
      horizons_year = horizons_year,
      risk_thresholds = risk_thresholds,
      remission_cut_years = remission_cut_years,
      fit_prior_branches = fit_prior_branches,
      fit_site_prior_families = fit_site_prior_families
    ),
    outputs = list(
      model_registry = model_registry,
      posterior_cohort_yearly = posterior_cohort_yearly,
      posterior_classification = posterior_classification,
      prior_predictive_summary = prior_predictive_summary,
      hazard_plausibility = hazard_shape_plausibility,
      uncured_decomposition = uncured_supporting_decomposition,
      anchor_vs_neutral_delta = anchor_vs_neutral_delta,
      incidence_anchor_update = incidence_anchor_update,
      stage8a_vs_stage8b_delta = delta_vs_stage8a,
      horizon_support_panel = horizon_support_panel,
      ppc_summary = ppc_summary
    )
  ),
  file.path(export_path, "bayes_stage8b_bundle.rds")
)
write_csv_preserve_schema(export_manifest, file.path(export_path, "bayes_stage8b_export_manifest.csv"))

message(
  "Stage8B completed. Reused models: ",
  n_models_reused,
  "; newly fitted models: ",
  n_models_to_fit,
  "; diagnostic PDF ok: ",
  pdf_ok,
  "."
)
