
# 🔴 Configure: Stage-8A paths, controls, and scientific constants ===============================
sys_name <- Sys.info()[["sysname"]]

project_root <- switch(
  sys_name,
  "Darwin"  = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  "Windows" = "C:/Users/clair/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  stop("Unsupported OS: ", sys_name)
)

merged_data_path   <- file.path(project_root, "data", "MERGED_dataset3_pnu_snu.csv")
export_path        <- file.path(project_root, "stage8A_Bayesian transition-only cure")
stage5_export_path <- file.path(project_root, "stage5_Individualized no-cure comparator")
stage6_export_path <- file.path(project_root, "stage6_Cure-appropriateness screening")

pnu_site_label <- "PNU"
snu_site_label <- "SNU"

horizons_year <- 1:10
risk_thresholds <- c(0.05, 0.10, 0.15, 0.20)

run_model_ids <- NULL
save_model_rds <- TRUE
reuse_existing_rds <- TRUE

stan_chains <- 4L
stan_iter <- 2000L
stan_warmup <- 1000L
stan_thin <- 1L
stan_seed <- 20260329L
stan_adapt_delta <- 0.95
stan_max_treedepth <- 12L
stan_refresh <- 0L

prior_predictive_draws <- 500L
posterior_prediction_draws <- 400L

ppc_tolerance_abs <- 0.15
ess_min_threshold <- 400
rhat_max_threshold <- 1.01

degenerate_draw_fraction_threshold <- 0.90
degenerate_subject_fraction_threshold <- 0.95
tiny_susceptible_prob <- 0.01
huge_susceptible_prob <- 0.99
tiny_median_years <- 0.05
huge_median_years <- 50

prior_tail_warning_mean_years <- 100
prior_tail_warning_q975_years <- 200

materiality_delta_risk <- 0.03
materiality_delta_cure_fraction <- 0.03
materiality_delta_false_positive_burden <- 0.03
materiality_delta_FP100 <- 3.0
materiality_delta_NB <- 0.01
materiality_delta_PPV <- 0.03
materiality_delta_TPR <- 0.03

stage8_branch_label <- "Stage8A"
stage8_risk_scale <- "transition_only_main"
stage8b_risk_scale <- "transition_cif_competing"
neutral_prior_branch_label <- "neutral_no_external_info"
anchor_prior_branch_label <- "anchor_informed"
site_prior_family_main <- "normal_0_1_main"
site_prior_family_sensitivity <- "student_t3_0_1_sensitivity"

stage6_screening_flag_csv <- file.path(stage6_export_path, "stage6_carry_forward_flag_table.csv")
stage5_model_performance_csv <- file.path(stage5_export_path, "stage5_model_performance_long.csv")

# 🔴 Initialize: packages, options, and output folder ===============================
core_packages <- c(
  "readr", "dplyr", "tibble", "tidyr", "ggplot2", "survival",
  "matrixStats", "rstan", "posterior", "loo"
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
  library(rstan)
  library(posterior)
  library(loo)
})

options(stringsAsFactors = FALSE, scipen = 999)
rstan_options(auto_write = TRUE)
options(mc.cores = max(1L, min(stan_chains, parallel::detectCores(logical = TRUE))))

dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

if (!identical(as.integer(horizons_year), 1:10)) {
  stop("`horizons_year` must be exactly 1:10.", call. = FALSE)
}
if (length(risk_thresholds) == 0L || anyNA(risk_thresholds) || any(risk_thresholds <= 0 | risk_thresholds >= 1)) {
  stop("`risk_thresholds` must be probabilities strictly between 0 and 1.", call. = FALSE)
}
if (toupper(pnu_site_label) == toupper(snu_site_label)) {
  stop("`pnu_site_label` and `snu_site_label` must differ.", call. = FALSE)
}

# 🔴 Define: low-level helpers for robust I/O and coercion ===============================
`%||%` <- function(x, y) if (is.null(x)) y else x

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

read_delimited_or_rds <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(NULL)
  }
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))
  if (grepl("\\.csv\\.gz$", path, ignore.case = TRUE)) {
    return(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  }
  if (ext == "csv") {
    return(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  }
  if (ext == "rds") {
    return(readRDS(path))
  }

  stop("Unsupported extension for file: ", path, call. = FALSE)
}

make_temp_output_path <- function(path, tag = "tmp") {
  stamp <- paste0(format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%08d", sample.int(99999999, 1)))
  dir <- dirname(path)
  base <- basename(path)

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
  dir.create(dirname(final_path), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(final_path)) {
    unlink(final_path)
  }

  ok <- file.rename(tmp_path, final_path)
  if (!isTRUE(ok)) {
    ok_copy <- file.copy(tmp_path, final_path, overwrite = TRUE)
    unlink(tmp_path)
    if (!isTRUE(ok_copy)) {
      stop("Failed to replace file: ", final_path, call. = FALSE)
    }
  }

  invisible(TRUE)
}

write_csv_preserve_schema <- function(df, path) {
  tmp <- make_temp_output_path(path, tag = "tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  readr::write_csv(df, tmp)
  safe_promote_file(tmp, path)
}

bind_rows_safe <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) == 0L) {
    return(tibble())
  }
  bind_rows(x)
}

first_existing_name <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) {
    return(NULL)
  }
  hit[[1L]]
}

nrow_or_zero <- function(df) {
  if (is.null(df) || !inherits(df, "data.frame")) {
    return(0L)
  }
  nrow(df)
}

pdf_file_is_usable <- function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  finfo <- file.info(path)
  isTRUE(!is.na(finfo$size[[1L]]) && finfo$size[[1L]] > 0)
}

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

simplify_scalar_list_cols <- function(df) {
  if (is.null(df) || !inherits(df, "data.frame")) {
    return(df)
  }
  out <- df
  for (nm in names(out)) {
    out[[nm]] <- flatten_scalar_list_col(out[[nm]])
  }
  out
}

stage8a_model_rds_path <- function(model_id, export_dir = export_path) {
  file.path(export_dir, paste0(model_id, "__bayes_stage8a_fit.rds"))
}

bundle_get <- function(bundle, name) {
  if (is.null(bundle) || !is.list(bundle)) {
    return(NULL)
  }
  bundle[[name]] %||% NULL
}

bundle_table_or_empty <- function(bundle, element, model_id = NULL) {
  obj <- bundle[[element]]
  if (is.null(obj) || !inherits(obj, "data.frame") || nrow(obj) == 0L) {
    return(tibble())
  }
  out <- tibble::as_tibble(obj)
  if (!is.null(model_id) && !("model_id" %in% names(out))) {
    out$model_id <- model_id
  }
  out
}

load_existing_stage8a_bundle <- function(path, model_row = NULL, dataset_df = NULL) {
  missing_bundle <- function(found = FALSE, reason = "missing_rds", bundle = NULL, registry = NULL, admissible = FALSE, reusable = FALSE) {
    list(found = found, reason = reason, bundle = bundle, registry = registry, admissible = admissible, reusable = reusable)
  }

  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(missing_bundle())
  }

  bundle <- tryCatch(readRDS(path), error = identity)
  if (inherits(bundle, "error")) {
    return(
      missing_bundle(
        found = TRUE,
        reason = paste0("rds_read_error: ", conditionMessage(bundle))
      )
    )
  }
  if (!is.list(bundle)) {
    return(missing_bundle(found = TRUE, reason = "invalid_bundle_object"))
  }

  registry <- bundle_get(bundle, "model_registry_row") %||% bundle_get(bundle, "model_registry")
  if (is.null(registry) || !inherits(registry, "data.frame") || nrow(registry) == 0L) {
    return(missing_bundle(found = TRUE, reason = "missing_model_registry_row", bundle = bundle))
  }
  registry <- tibble::as_tibble(registry[1L, , drop = FALSE])

  admissible_flag <- isTRUE(as.logical(registry$admissible_flag[[1L]] %||% FALSE))
  registry_model_id <- as.character(registry$model_id[[1L]] %||% NA_character_)
  registry_dataset_key <- as.character(registry$dataset_key[[1L]] %||% NA_character_)
  registry_fit_status <- as.character(registry$fit_status[[1L]] %||% NA_character_)
  registry_n <- suppressWarnings(as.integer(registry$n[[1L]] %||% NA_integer_))

  required_names <- c(
    "coefficient_summary",
    "diagnostics_parameter_level",
    "ppc_summary",
    "prior_predictive_summary",
    "prediction_long",
    "cohort_yearly",
    "classification",
    "uncured_decomposition",
    "hazard_plausibility"
  )
  has_required_outputs <- all(vapply(required_names, function(nm) {
    !is.null(bundle_get(bundle, nm))
  }, logical(1)))

  inc_draws <- bundle_get(bundle, "posterior_incidence_draws")
  has_required_draws <- !isTRUE(admissible_flag) || (
    is.list(inc_draws) &&
      all(c("alpha_inc", "beta_inc") %in% names(inc_draws))
  )

  expected_model_id <- if (!is.null(model_row) && "model_id" %in% names(model_row)) {
    as.character(model_row$model_id[[1L]])
  } else {
    NA_character_
  }
  expected_dataset_key <- if (!is.null(model_row) && "dataset_key" %in% names(model_row)) {
    as.character(model_row$dataset_key[[1L]])
  } else {
    NA_character_
  }
  expected_n <- if (!is.null(dataset_df)) nrow(dataset_df) else NA_integer_

  reusable_reason <- dplyr::case_when(
    !is.na(expected_model_id) && is.na(registry_model_id) ~ "missing_model_id",
    !is.na(expected_model_id) && !is.na(registry_model_id) && !identical(registry_model_id, expected_model_id) ~ "model_id_mismatch",
    !is.na(expected_dataset_key) && is.na(registry_dataset_key) ~ "missing_dataset_key",
    !is.na(expected_dataset_key) && !is.na(registry_dataset_key) && !identical(registry_dataset_key, expected_dataset_key) ~ "dataset_key_mismatch",
    !is.na(expected_n) && !is.na(registry_n) && !identical(registry_n, expected_n) ~ "dataset_size_mismatch",
    !identical(registry_fit_status, "ok") ~ "fit_status_not_ok",
    !has_required_outputs ~ "missing_required_outputs",
    !has_required_draws ~ "missing_posterior_incidence_draws",
    TRUE ~ "existing_rds"
  )

  reusable_flag <- identical(reusable_reason, "existing_rds")

  missing_bundle(
    found = TRUE,
    reason = reusable_reason,
    bundle = bundle,
    registry = registry,
    admissible = admissible_flag,
    reusable = reusable_flag
  )
}

is_stage8a_bundle_reusable <- function(bundle_check) {
  isTRUE(bundle_check$reusable)
}

default_stage8a_bundle_check <- function(reason = "reuse_not_evaluated") {
  list(
    found = FALSE,
    reason = reason,
    bundle = NULL,
    registry = NULL,
    admissible = FALSE,
    reusable = FALSE
  )
}

# 🔴 Define: Stage-8A semantic helpers and support governance ===============================
normalize_dataset_label <- function(x) {
  x_chr <- trimws(as.character(x))
  x_up <- toupper(x_chr)
  out <- dplyr::case_when(
    x_up %in% c("PNU", "P") ~ "PNU",
    x_up %in% c("SNU", "S") ~ "SNU",
    x_up %in% c("MERGED", "M", "BOTH", "ALL_SITES") ~ "merged",
    TRUE ~ x_chr
  )
  out[is.na(x)] <- NA_character_
  out
}

normalize_stage6_formula_anchor <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr[is.na(x_chr) | x_chr == ""] <- "ALL"
  x_key <- toupper(gsub("[^A-Z0-9]+", "", x_chr))
  out <- dplyr::case_when(
    x_key %in% c("ALL", "ANY", "GLOBAL") ~ "ALL",
    x_key %in% c("BASE", "MAIN", "BASELINE") ~ "base",
    x_key %in% c("INTERACTION") ~ "interaction",
    x_key %in% c("SITEADDED", "SITE") ~ "site_added",
    x_key %in% c("SITEINTERACTION") ~ "site_interaction",
    TRUE ~ x_chr
  )
  out[is.na(out) | out == ""] <- "ALL"
  out
}

stage8_collapse_unique_text <- function(x, sep = "|") {
  x_chr <- trimws(as.character(x))
  x_chr <- unique(x_chr[!is.na(x_chr) & nzchar(x_chr) & !(toupper(x_chr) %in% c("NA", "NULL"))])
  if (length(x_chr) == 0L) {
    return(NA_character_)
  }
  paste(x_chr, collapse = sep)
}

stage8_parse_logicalish <- function(x) {
  x_chr <- trimws(toupper(as.character(x)))
  out <- rep(NA, length(x_chr))
  out[x_chr %in% c("TRUE", "T", "1", "YES", "Y")] <- TRUE
  out[x_chr %in% c("FALSE", "F", "0", "NO", "N")] <- FALSE
  as.logical(out)
}

screening_value_or_na <- function(screening_row, field) {
  if (is.null(screening_row) || !inherits(screening_row, "data.frame") || nrow(screening_row) == 0L) {
    return(NA_character_)
  }
  if (!(field %in% names(screening_row))) {
    return(NA_character_)
  }
  value <- screening_row[[field]][[1L]]
  if (is.null(value) || length(value) == 0L) {
    return(NA_character_)
  }
  as.character(value)
}

screening_logical_or_na <- function(screening_row, field) {
  if (is.null(screening_row) || !inherits(screening_row, "data.frame") || nrow(screening_row) == 0L) {
    return(NA)
  }
  if (!(field %in% names(screening_row))) {
    return(NA)
  }
  value <- screening_row[[field]][[1L]]
  if (is.null(value) || length(value) == 0L) {
    return(NA)
  }
  if (is.logical(value)) {
    return(as.logical(value))
  }
  stage8_parse_logicalish(value)[[1L]]
}

derive_support_tier <- function(dataset_name, horizon_year) {
  h <- as.integer(horizon_year)
  if (dataset_name == "PNU") {
    if (h == 1L) return("primary_supported")
    if (h == 2L) return("sensitivity")
    return("projection")
  }
  if (dataset_name %in% c("SNU", "merged")) {
    if (h %in% c(1L, 2L)) return("primary_supported")
    if (h %in% c(3L, 4L, 5L)) return("secondary")
    return("projection")
  }
  "projection"
}

derive_horizon_evidence_class <- function(dataset_name, horizon_year) {
  h <- as.integer(horizon_year)
  if (dataset_name == "PNU") {
    if (h == 1L) return("directly_observed_data_supported")
    if (h == 2L) return("partly_model_dependent")
    return("mostly_extrapolated")
  }
  if (dataset_name %in% c("SNU", "merged")) {
    if (h %in% c(1L, 2L)) return("directly_observed_data_supported")
    if (h %in% c(3L, 4L, 5L)) return("partly_model_dependent")
    return("mostly_extrapolated")
  }
  "mostly_extrapolated"
}

derive_claim_restriction_flag <- function(horizon_evidence_class) {
  if (is.na(horizon_evidence_class)) {
    return(NA_character_)
  }
  if (horizon_evidence_class == "directly_observed_data_supported") {
    return("primary_claim_allowed")
  }
  if (horizon_evidence_class == "partly_model_dependent") {
    return("secondary_or_sensitivity_only")
  }
  "projection_only"
}

derive_interpretation_note <- function(dataset_name, horizon_year, support_tier, horizon_evidence_class, claim_restriction_flag) {
  h <- as.integer(horizon_year)

  if (support_tier == "primary_supported" && horizon_evidence_class == "directly_observed_data_supported") {
    return("Primary supported horizon with comparatively direct follow-up support.")
  }
  if (dataset_name == "PNU" && h == 2L) {
    return("Sensitivity horizon for PNU; partly model-dependent and not for primary claims.")
  }
  if (support_tier == "secondary" && horizon_evidence_class == "partly_model_dependent") {
    return("Secondary horizon with growing tail uncertainty; interpret with explicit model-dependence.")
  }
  if (claim_restriction_flag == "projection_only") {
    return("Projection-dominant horizon; mostly extrapolated and not eligible for primary claims.")
  }

  "Common comparison horizon retained for cross-model comparability."
}

build_horizon_metadata <- function(analysis_datasets, horizons, ipcw_registry) {
  out_rows <- vector("list", length(analysis_datasets))
  dataset_names <- names(analysis_datasets)

  for (ii in seq_along(dataset_names)) {
    dataset_name <- dataset_names[[ii]]
    ipcw_df <- ipcw_registry[[dataset_name]]

    out_rows[[ii]] <- tibble(
      dataset_key = dataset_name,
      horizon = as.integer(horizons),
      support_tier = vapply(horizons, function(h) derive_support_tier(dataset_name, h), character(1)),
      horizon_evidence_class = vapply(horizons, function(h) derive_horizon_evidence_class(dataset_name, h), character(1)),
      claim_restriction_flag = vapply(horizons, function(h) derive_claim_restriction_flag(derive_horizon_evidence_class(dataset_name, h)), character(1)),
      interpretation_note = mapply(
        derive_interpretation_note,
        dataset_name,
        horizons,
        vapply(horizons, function(h) derive_support_tier(dataset_name, h), character(1)),
        vapply(horizons, function(h) derive_horizon_evidence_class(dataset_name, h), character(1)),
        vapply(horizons, function(h) derive_claim_restriction_flag(derive_horizon_evidence_class(dataset_name, h)), character(1)),
        USE.NAMES = FALSE
      )
    ) %>%
      left_join(
        ipcw_df %>%
          select(horizon, observed_km_risk, denom_case, denom_control),
        by = "horizon"
      ) %>%
      mutate(
        risk_set_support_flag = denom_case > 0 & denom_control > 0
      )
  }

  bind_rows(out_rows) %>%
    arrange(match(dataset_key, names(analysis_datasets)), horizon)
}


# 🔴 Define: upstream Stage-5 and Stage-6 normalization helpers ===============================
read_screening_flags <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(tibble())
  }

  raw_obj <- read_delimited_or_rds(path)
  df <- if (inherits(raw_obj, "data.frame")) tibble::as_tibble(raw_obj) else tibble()

  if (nrow(df) == 0L) {
    return(tibble())
  }

  dataset_col <- first_existing_name(df, c("source_dataset", "dataset", "cohort"))
  dataset_key_col <- first_existing_name(df, c("dataset_key", "screening_dataset_key"))
  formula_col <- first_existing_name(df, c("formula_anchor", "formula_name", "formula_type"))
  eligibility_col <- first_existing_name(df, c("cure_model_eligibility_flag", "stage6_final_class", "final_decision_flag", "decision_flag", "screening_flag"))
  primary_gate_method_col <- first_existing_name(df, c("primary_gate_method", "method", "screening_method"))
  primary_gate_flag_col <- first_existing_name(df, c("primary_gate_flag"))
  receus_col <- first_existing_name(df, c("receus_primary_class"))
  presence_modifier_col <- first_existing_name(df, c("presence_modifier_flag"))
  cure_presence_col <- first_existing_name(df, c("cure_presence_support_flag"))
  followup_contradiction_col <- first_existing_name(df, c("followup_contradiction_flag"))
  followup_not_col <- first_existing_name(df, c("followup_not_contradicted_flag"))
  note_col <- first_existing_name(df, c("screening_note", "screening_detail", "screening_context", "detail"))
  carry_forward_col <- first_existing_name(df, c("carry_forward_stage8"))

  if (is.null(eligibility_col)) {
    return(tibble())
  }

  dataset_value <- if (!is.null(dataset_col)) {
    normalize_dataset_label(df[[dataset_col]])
  } else if (!is.null(dataset_key_col)) {
    normalize_dataset_label(sub("__.*$", "", as.character(df[[dataset_key_col]])))
  } else {
    rep(NA_character_, nrow(df))
  }

  formula_value <- if (!is.null(formula_col)) {
    normalize_stage6_formula_anchor(df[[formula_col]])
  } else {
    rep("ALL", nrow(df))
  }

  out <- tibble(
    dataset_key = dataset_value,
    formula_anchor = formula_value,
    cure_model_eligibility_flag = as.character(df[[eligibility_col]]),
    primary_gate_method = if (!is.null(primary_gate_method_col)) as.character(df[[primary_gate_method_col]]) else NA_character_,
    primary_gate_flag = if (!is.null(primary_gate_flag_col)) as.character(df[[primary_gate_flag_col]]) else NA_character_,
    receus_primary_class = if (!is.null(receus_col)) as.character(df[[receus_col]]) else NA_character_,
    presence_modifier_flag = if (!is.null(presence_modifier_col)) as.character(df[[presence_modifier_col]]) else NA_character_,
    cure_presence_support_flag = if (!is.null(cure_presence_col)) as.character(df[[cure_presence_col]]) else NA_character_,
    followup_contradiction_flag = if (!is.null(followup_contradiction_col)) as.character(df[[followup_contradiction_col]]) else NA_character_,
    followup_not_contradicted_flag = if (!is.null(followup_not_col)) as.character(df[[followup_not_col]]) else NA_character_,
    screening_note = if (!is.null(note_col)) as.character(df[[note_col]]) else NA_character_,
    carry_forward_stage8 = if (!is.null(carry_forward_col)) stage8_parse_logicalish(df[[carry_forward_col]]) else NA
  ) %>%
    filter(dataset_key %in% c("PNU", "SNU", "merged")) %>%
    group_by(dataset_key, formula_anchor) %>%
    summarise(
      cure_model_eligibility_flag = stage8_collapse_unique_text(cure_model_eligibility_flag),
      primary_gate_method = stage8_collapse_unique_text(primary_gate_method),
      primary_gate_flag = stage8_collapse_unique_text(primary_gate_flag),
      receus_primary_class = stage8_collapse_unique_text(receus_primary_class),
      presence_modifier_flag = stage8_collapse_unique_text(presence_modifier_flag),
      cure_presence_support_flag = stage8_collapse_unique_text(cure_presence_support_flag),
      followup_contradiction_flag = stage8_collapse_unique_text(followup_contradiction_flag),
      followup_not_contradicted_flag = stage8_collapse_unique_text(followup_not_contradicted_flag),
      screening_note = stage8_collapse_unique_text(screening_note, sep = "; "),
      carry_forward_stage8 = if (all(is.na(carry_forward_stage8))) NA else any(carry_forward_stage8 %in% TRUE),
      .groups = "drop"
    )

  out
}

build_screening_model_lookup <- function(screening_flags, model_grid) {
  base_lookup <- model_grid %>%
    distinct(model_id, dataset_key, formula_anchor)

  if (nrow(screening_flags) == 0L) {
    return(
      base_lookup %>%
        mutate(
          cure_model_eligibility_flag = NA_character_,
          primary_gate_method = NA_character_,
          primary_gate_flag = NA_character_,
          receus_primary_class = NA_character_,
          presence_modifier_flag = NA_character_,
          cure_presence_support_flag = NA_character_,
          followup_contradiction_flag = NA_character_,
          followup_not_contradicted_flag = NA_character_,
          screening_note = NA_character_,
          carry_forward_stage8 = NA
        )
    )
  }

  exact_lookup <- base_lookup %>%
    left_join(
      screening_flags %>%
        filter(formula_anchor != "ALL") %>%
        distinct(dataset_key, formula_anchor, .keep_all = TRUE),
      by = c("dataset_key", "formula_anchor")
    )

  all_lookup <- base_lookup %>%
    left_join(
      screening_flags %>%
        filter(formula_anchor == "ALL") %>%
        distinct(dataset_key, .keep_all = TRUE) %>%
        rename_with(~paste0(.x, "_all"), -dataset_key),
      by = "dataset_key"
    )

  out <- exact_lookup %>%
    select(model_id)

  merge_fields <- c(
    "cure_model_eligibility_flag", "primary_gate_method", "primary_gate_flag",
    "receus_primary_class", "presence_modifier_flag", "cure_presence_support_flag",
    "followup_contradiction_flag", "followup_not_contradicted_flag",
    "screening_note", "carry_forward_stage8"
  )

  for (nm in merge_fields) {
    out[[nm]] <- dplyr::coalesce(exact_lookup[[nm]], all_lookup[[paste0(nm, "_all")]])
  }

  out
}

normalize_nocure_cohort <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(tibble())
  }

  df <- read_delimited_or_rds(path)
  if (!inherits(df, "data.frame")) {
    return(tibble())
  }
  df <- tibble::as_tibble(df)

  if (all(c("metric_domain", "metric_name", "metric_value") %in% names(df))) {
    out <- df %>%
      filter(
        model_class == "non_cure",
        metric_domain == "horizon_summary",
        metric_name == "mean_predicted_risk"
      ) %>%
      transmute(
        dataset_key = normalize_dataset_label(dataset),
        formula_anchor = normalize_stage6_formula_anchor(.data$formula_name %||% "ALL"),
        no_cure_model_id = as.character(model_id),
        horizon = as.integer(safe_numeric(horizon_year)),
        metric = "risk_mean",
        value = safe_numeric(metric_value)
      ) %>%
      filter(!is.na(horizon), !is.na(value))

    return(out)
  }

  tibble()
}

normalize_nocure_classification <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(tibble())
  }

  df <- read_delimited_or_rds(path)
  if (!inherits(df, "data.frame")) {
    return(tibble())
  }
  df <- tibble::as_tibble(df)

  if (all(c("metric_domain", "metric_name", "metric_value") %in% names(df))) {
    metric_map <- c(
      false_positive_burden_nonevents = "false_positive_burden_mean",
      false_positive_per_100 = "FP100_mean",
      net_benefit = "NB_mean",
      ppv = "PPV_mean",
      tpr = "TPR_mean"
    )

    out <- df %>%
      filter(
        model_class == "non_cure",
        metric_domain == "threshold_summary",
        metric_name %in% names(metric_map)
      ) %>%
      transmute(
        dataset_key = normalize_dataset_label(dataset),
        formula_anchor = normalize_stage6_formula_anchor(.data$formula_name %||% "ALL"),
        no_cure_model_id = as.character(model_id),
        horizon = as.integer(safe_numeric(horizon_year)),
        threshold = safe_numeric(threshold),
        metric = unname(metric_map[metric_name]),
        value = safe_numeric(metric_value)
      ) %>%
      filter(!is.na(horizon), !is.na(threshold), !is.na(value))

    return(out)
  }

  tibble()
}

# 🔴 Define: backbone loaders, prior artifacts, and model grid ===============================
prepare_analysis_dataset <- function(df, dataset_name, pnu_label, snu_label) {
  if (is.null(df) || !inherits(df, "data.frame")) {
    stop("Input dataset is NULL or not a data frame for ", dataset_name, call. = FALSE)
  }

  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop("[", dataset_name, "] Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
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
    stop("[", dataset_name, "] Dataset has zero rows.", call. = FALSE)
  }
  if (anyNA(out[required_cols])) {
    stop("[", dataset_name, "] Missing values detected in required backbone columns.", call. = FALSE)
  }
  if (any(out$id == "")) {
    stop("[", dataset_name, "] Blank `id` values detected.", call. = FALSE)
  }
  if (any(!out$sex_num %in% c(0L, 1L))) {
    stop("[", dataset_name, "] `sex_num` must be coded 0/1.", call. = FALSE)
  }
  if (any(!out$status_num %in% c(0L, 1L, 2L))) {
    stop("[", dataset_name, "] `status_num` must be coded 0/1/2.", call. = FALSE)
  }
  if (any(out$days_followup < 0)) {
    stop("[", dataset_name, "] Negative `days_followup` values detected.", call. = FALSE)
  }

  if (dataset_name != "merged") {
    out <- out %>% mutate(site = dataset_name)
  }

  age_mean <- mean(out$age_exact_entry)
  age_sd <- stats::sd(out$age_exact_entry)
  if (is.na(age_sd) || age_sd <= 0) {
    stop("[", dataset_name, "] `age_exact_entry` must have positive SD.", call. = FALSE)
  }

  out <- out %>%
    mutate(
      unique_person_id = paste(site, id, sep = "_"),
      time_year = pmax(days_followup / 365.25, 1e-8),
      event_main = as.integer(status_num == 1L),
      right_censor_flag = as.integer(status_num == 0L),
      remission_flag = as.integer(status_num == 2L),
      censor_main = as.integer(status_num %in% c(0L, 2L)),
      age_s = (age_exact_entry - age_mean) / (2 * age_sd),
      sex_label = factor(if_else(sex_num == 0L, "Female", "Male"), levels = c("Female", "Male"))
    )

  if (dplyr::n_distinct(out$unique_person_id) != nrow(out)) {
    stop("[", dataset_name, "] `site + id` is not unique.", call. = FALSE)
  }
  if (dataset_name == "merged" && dplyr::n_distinct(out$site) < 2L) {
    stop("[merged] Merged dataset must contain at least two site levels.", call. = FALSE)
  }

  list(
    data = out,
    summary = tibble(
      dataset_key = dataset_name,
      n = nrow(out),
      n_transition = sum(out$event_main),
      n_remission = sum(out$remission_flag),
      n_right_censoring = sum(out$right_censor_flag),
      n_main_censoring = sum(out$censor_main),
      age_mean = age_mean,
      age_sd = age_sd,
      max_followup_years = max(out$time_year),
      median_followup_years = stats::median(out$time_year),
      site_values = paste(sort(unique(out$site)), collapse = "|")
    )
  )
}

load_stage8a_datasets <- function(merged_path, pnu_label, snu_label) {
  merged_raw <- read_delimited_or_rds(merged_path)
  if (!inherits(merged_raw, "data.frame")) {
    stop("Merged dataset could not be loaded.", call. = FALSE)
  }

  merged_raw <- tibble::as_tibble(merged_raw)
  if (!("site" %in% names(merged_raw))) {
    stop("Merged dataset must contain a `site` column.", call. = FALSE)
  }

  merged_raw <- merged_raw %>%
    mutate(
      site = trimws(as.character(site)),
      site = case_when(
        toupper(site) == toupper(pnu_label) ~ pnu_label,
        toupper(site) == toupper(snu_label) ~ snu_label,
        TRUE ~ site
      )
    )

  pnu_raw <- merged_raw %>% filter(site == pnu_label)
  snu_raw <- merged_raw %>% filter(site == snu_label)

  if (nrow(pnu_raw) == 0L) {
    stop("No rows found for PNU site label `", pnu_label, "`.", call. = FALSE)
  }
  if (nrow(snu_raw) == 0L) {
    stop("No rows found for SNU site label `", snu_label, "`.", call. = FALSE)
  }

  pnu_obj <- prepare_analysis_dataset(pnu_raw, "PNU", pnu_label, snu_label)
  snu_obj <- prepare_analysis_dataset(snu_raw, "SNU", pnu_label, snu_label)
  merged_obj <- prepare_analysis_dataset(merged_raw, "merged", pnu_label, snu_label)

  list(
    datasets = list(
      PNU = pnu_obj$data,
      SNU = snu_obj$data,
      merged = merged_obj$data
    ),
    dataset_registry = bind_rows(pnu_obj$summary, snu_obj$summary, merged_obj$summary)
  )
}

build_prior_artifacts <- function() {
  alpha_gp_vdw <- -9.581369553169
  mu_beta_inc_anchor <- c(
    0.419871845822, 0.907608052926, 0.586202561451, 0.466865123863, 0.037997248763
  )
  sd_beta_inc_anchor <- c(
    0.132789397422, 0.173731076538, 0.191221553945, 0.270393197518, 0.302838606651
  )

  anchor_cell_df <- tidyr::crossing(
    sex_num = c(0L, 1L),
    age_band = c("<20", "20_29", "30plus")
  ) %>%
    mutate(
      age20_29 = as.integer(age_band == "20_29"),
      age30plus = as.integer(age_band == "30plus"),
      sex_x_age20_29 = sex_num * age20_29,
      sex_x_age30plus = sex_num * age30plus,
      age_sex_anchor_cell = paste0(if_else(sex_num == 0L, "Female", "Male"), "_", age_band),
      prior_center_logit_anchor = alpha_gp_vdw +
        sex_num * mu_beta_inc_anchor[[1L]] +
        age20_29 * mu_beta_inc_anchor[[2L]] +
        age30plus * mu_beta_inc_anchor[[3L]] +
        sex_x_age20_29 * mu_beta_inc_anchor[[4L]] +
        sex_x_age30plus * mu_beta_inc_anchor[[5L]],
      external_one_year_risk = plogis(prior_center_logit_anchor),
      external_incidence_rate_per10k = -log(pmax(1 - external_one_year_risk, 1e-12)) * 10000
    ) %>%
    select(
      age_sex_anchor_cell,
      sex_num,
      age_band,
      external_incidence_rate_per10k,
      external_one_year_risk,
      prior_center_logit_anchor
    )

  list(
    alpha_gp_vdw = alpha_gp_vdw,
    mu_beta_inc_anchor = mu_beta_inc_anchor,
    sd_beta_inc_anchor = sd_beta_inc_anchor,
    alpha_inc_neutral_mean = 0,
    alpha_inc_sd = 3.5,
    sd_beta_inc_neutral = 2.0,
    sd_gamma0 = 2.5,
    sd_gamma_lat = 1.0,
    sd_shape_W = 0.35,
    sd_log_sigma_LN = 0.50,
    sd_psi_LL = 0.50,
    anchor_cells = anchor_cell_df
  )
}

site_placement_label <- function(incidence_site_indicator, latency_site_indicator) {
  dplyr::case_when(
    isTRUE(incidence_site_indicator) & isTRUE(latency_site_indicator) ~ "site_in_both",
    isTRUE(incidence_site_indicator) & !isTRUE(latency_site_indicator) ~ "site_in_incidence_only",
    !isTRUE(incidence_site_indicator) & isTRUE(latency_site_indicator) ~ "site_in_latency_only",
    TRUE ~ "site_in_neither"
  )
}

build_model_grid <- function() {
  base_rows <- tibble::tribble(
    ~dataset_key, ~structural_model_id, ~incidence_site_indicator, ~latency_site_indicator, ~latency_interaction_indicator, ~formula_anchor,
    "PNU",        "PNU-L0",              FALSE,                      FALSE,                     FALSE,                         "base",
    "PNU",        "PNU-L1",              FALSE,                      FALSE,                     TRUE,                          "interaction",
    "SNU",        "SNU-L0",              FALSE,                      FALSE,                     FALSE,                         "base",
    "SNU",        "SNU-L1",              FALSE,                      FALSE,                     TRUE,                          "interaction",
    "merged",     "MERGED-I0-L0S0",      FALSE,                      FALSE,                     FALSE,                         "base",
    "merged",     "MERGED-I0-L1S0",      FALSE,                      FALSE,                     TRUE,                          "interaction",
    "merged",     "MERGED-I0-L0S1",      FALSE,                      TRUE,                      FALSE,                         "site_added",
    "merged",     "MERGED-I0-L1S1",      FALSE,                      TRUE,                      TRUE,                          "site_interaction",
    "merged",     "MERGED-I1-L0S0",      TRUE,                       FALSE,                     FALSE,                         "base",
    "merged",     "MERGED-I1-L1S0",      TRUE,                       FALSE,                     TRUE,                          "interaction",
    "merged",     "MERGED-I1-L0S1",      TRUE,                       TRUE,                      FALSE,                         "site_added",
    "merged",     "MERGED-I1-L1S1",      TRUE,                       TRUE,                      TRUE,                          "site_interaction"
  )

  family_rows <- tibble::tribble(
    ~family_code, ~latency_family, ~family_id,
    "E",          "exponential",   1L,
    "W",          "weibull",       2L,
    "LN",         "lognormal",     3L,
    "LL",         "loglogistic",   4L
  )

  seed_grid <- tidyr::crossing(
    base_rows,
    family_rows,
    prior_branch = c(anchor_prior_branch_label, neutral_prior_branch_label)
  ) %>%
    mutate(
      has_site_term = dataset_key == "merged" & (incidence_site_indicator | latency_site_indicator),
      site_placement = mapply(site_placement_label, incidence_site_indicator, latency_site_indicator, USE.NAMES = FALSE)
    )

  out_rows <- vector("list", nrow(seed_grid))
  for (ii in seq_len(nrow(seed_grid))) {
    one <- seed_grid[ii, , drop = FALSE]
    site_prior_values <- if (isTRUE(one$has_site_term[[1L]])) {
      c(site_prior_family_main, site_prior_family_sensitivity)
    } else {
      NA_character_
    }

    out_rows[[ii]] <- tibble(
      dataset_key = one$dataset_key[[1L]],
      structural_model_id = one$structural_model_id[[1L]],
      incidence_site_indicator = one$incidence_site_indicator[[1L]],
      latency_site_indicator = one$latency_site_indicator[[1L]],
      latency_interaction_indicator = one$latency_interaction_indicator[[1L]],
      formula_anchor = one$formula_anchor[[1L]],
      family_code = one$family_code[[1L]],
      latency_family = one$latency_family[[1L]],
      family_id = one$family_id[[1L]],
      prior_branch = one$prior_branch[[1L]],
      site_prior_family = site_prior_values,
      has_site_term = one$has_site_term[[1L]],
      site_placement = one$site_placement[[1L]]
    )
  }

  bind_rows(out_rows) %>%
    mutate(
      branch = stage8_branch_label,
      risk_scale = stage8_risk_scale,
      site_prior_family_export = if_else(has_site_term, site_prior_family, NA_character_),
      prior_branch_code = if_else(prior_branch == anchor_prior_branch_label, "ANCH", "NEUT"),
      site_prior_code = case_when(
        !has_site_term ~ "NOSITE",
        site_prior_family == site_prior_family_main ~ "N01",
        site_prior_family == site_prior_family_sensitivity ~ "T31",
        TRUE ~ "UNK"
      ),
      model_id = paste(structural_model_id, family_code, prior_branch_code, site_prior_code, sep = "-")
    ) %>%
    arrange(dataset_key, structural_model_id, family_code, prior_branch, site_prior_code)
}

make_design_bundle <- function(df, model_row, prior_artifacts, snu_label) {
  z_i <- as.integer(df$sex_num)
  x20_i <- as.integer(df$age_exact_entry >= 20 & df$age_exact_entry < 30)
  x30_i <- as.integer(df$age_exact_entry >= 30)
  s_i <- as.integer(df$site == snu_label)

  X_inc_base <- cbind(
    sex_num = z_i,
    age20_29 = x20_i,
    age30plus = x30_i,
    sex_x_age20_29 = z_i * x20_i,
    sex_x_age30plus = z_i * x30_i
  )

  if (identical(model_row$prior_branch[[1L]], anchor_prior_branch_label)) {
    mu_beta_inc <- prior_artifacts$mu_beta_inc_anchor
    sd_beta_inc <- prior_artifacts$sd_beta_inc_anchor
    alpha_inc_prior_mean <- prior_artifacts$alpha_gp_vdw
  } else {
    mu_beta_inc <- rep(0, 5L)
    sd_beta_inc <- rep(prior_artifacts$sd_beta_inc_neutral, 5L)
    alpha_inc_prior_mean <- prior_artifacts$alpha_inc_neutral_mean
  }

  X_inc <- X_inc_base
  inc_is_site <- rep(0L, ncol(X_inc))
  if (isTRUE(model_row$incidence_site_indicator[[1L]])) {
    X_inc <- cbind(X_inc, site_SNU = s_i)
    mu_beta_inc <- c(mu_beta_inc, 0)
    sd_beta_inc <- c(sd_beta_inc, 1)
    inc_is_site <- c(inc_is_site, 1L)
  }

  a_i <- as.numeric(df$age_s)
  az_i <- a_i * z_i

  X_lat <- switch(
    model_row$structural_model_id[[1L]],
    "PNU-L0" = cbind(age_s = a_i, sex_num = z_i),
    "PNU-L1" = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i),
    "SNU-L0" = cbind(age_s = a_i, sex_num = z_i),
    "SNU-L1" = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i),
    "MERGED-I0-L0S0" = cbind(age_s = a_i, sex_num = z_i),
    "MERGED-I0-L1S0" = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i),
    "MERGED-I0-L0S1" = cbind(age_s = a_i, sex_num = z_i, site_SNU = s_i),
    "MERGED-I0-L1S1" = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i, site_SNU = s_i),
    "MERGED-I1-L0S0" = cbind(age_s = a_i, sex_num = z_i),
    "MERGED-I1-L1S0" = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i),
    "MERGED-I1-L0S1" = cbind(age_s = a_i, sex_num = z_i, site_SNU = s_i),
    "MERGED-I1-L1S1" = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i, site_SNU = s_i),
    stop("Unknown structural model id: ", model_row$structural_model_id[[1L]], call. = FALSE)
  )

  lat_is_site <- as.integer(colnames(X_lat) == "site_SNU")
  sd_gamma_lat <- rep(prior_artifacts$sd_gamma_lat, ncol(X_lat))

  list(
    X_inc = unclass(as.matrix(X_inc)),
    inc_is_site = as.integer(inc_is_site),
    mu_beta_inc = as.numeric(mu_beta_inc),
    sd_beta_inc = as.numeric(sd_beta_inc),
    alpha_inc_prior_mean = as.numeric(alpha_inc_prior_mean),
    alpha_inc_prior_sd = as.numeric(prior_artifacts$alpha_inc_sd),
    X_lat = unclass(as.matrix(X_lat)),
    lat_is_site = as.integer(lat_is_site),
    sd_gamma_lat = as.numeric(sd_gamma_lat),
    sd_gamma0 = as.numeric(prior_artifacts$sd_gamma0),
    sd_shape_W = as.numeric(prior_artifacts$sd_shape_W),
    sd_log_sigma_LN = as.numeric(prior_artifacts$sd_log_sigma_LN),
    sd_psi_LL = as.numeric(prior_artifacts$sd_psi_LL),
    site_prior_family_id = if (!isTRUE(model_row$has_site_term[[1L]])) {
      0L
    } else if (identical(model_row$site_prior_family[[1L]], site_prior_family_main)) {
      1L
    } else {
      2L
    },
    time = as.numeric(df$time_year),
    event = as.integer(df$event_main),
    id_df = df %>%
      select(unique_person_id, id, site, sex_num, age_exact_entry, age_s)
  )
}


# 🔴 Define: censoring-aware benchmarks, survival math, and Bayesian kernels ===============================
km_eval <- function(survfit_obj, times) {
  base_times <- survfit_obj$time
  base_surv <- survfit_obj$surv

  vapply(times, function(tt) {
    idx <- max(c(0L, which(base_times <= tt)))
    if (idx == 0L) {
      1
    } else {
      base_surv[[idx]]
    }
  }, numeric(1))
}

build_ipcw_reference <- function(df, horizons) {
  event_fit <- survival::survfit(survival::Surv(time_year, event_main) ~ 1, data = df)
  censor_fit <- survival::survfit(survival::Surv(time_year, censor_main) ~ 1, data = df)

  out_rows <- vector("list", length(horizons))
  for (ii in seq_along(horizons)) {
    h <- horizons[[ii]]
    G_t <- pmax(km_eval(censor_fit, h), 1e-8)
    G_tminus <- pmax(km_eval(censor_fit, pmax(df$time_year - 1e-10, 0)), 1e-8)
    w_case <- ifelse(df$event_main == 1L & df$time_year <= h, 1 / G_tminus, 0)
    w_control <- ifelse(df$time_year > h, 1 / G_t, 0)
    prevalence <- 1 - km_eval(event_fit, h)

    out_rows[[ii]] <- tibble(
      horizon = as.integer(h),
      observed_km_risk = as.numeric(prevalence),
      denom_case = sum(w_case),
      denom_control = sum(w_control),
      w_case = list(w_case),
      w_control = list(w_control)
    )
  }

  bind_rows(out_rows)
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
    q025 = q[, 1L],
    q50 = q[, 2L],
    q975 = q[, 3L]
  )
}

compute_linear_terms <- function(draws, X_inc, X_lat) {
  eta_inc <- draws$alpha_inc + draws$beta_inc %*% t(X_inc)
  pi_mat <- plogis(eta_inc)

  mu_lat <- draws$gamma0 + draws$gamma_lat %*% t(X_lat)
  median_mat <- exp(mu_lat)

  list(
    eta_inc_mat = eta_inc,
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
    lambda <- sweep(median_mat, 1L, (log(2))^(1 / k), FUN = "/")
    ratio <- sweep(lambda, 1L, horizon_year, FUN = function(lam, tt) tt / lam)
    k_mat <- matrix(k, nrow = length(k), ncol = ncol(lambda))
    Su <- exp(-(ratio^k_mat))
    haz <- (k_mat / lambda) * (ratio^(k_mat - 1))
  } else if (family_code == "LN") {
    sigma <- exp(draws$log_sigma_LN)
    sigma_mat <- matrix(sigma, nrow = length(sigma), ncol = ncol(mu_lat_mat))
    z <- (log(horizon_year) - mu_lat_mat) / sigma_mat
    Su <- 1 - pnorm(z)
    log_pdf <- dnorm(z, log = TRUE) - log(horizon_year) - log(sigma_mat)
    haz <- exp(log_pdf) / pmax(Su, 1e-12)
  } else if (family_code == "LL") {
    k <- exp(-draws$psi_LL)
    lambda <- median_mat
    ratio <- sweep(lambda, 1L, horizon_year, FUN = function(lam, tt) tt / lam)
    k_mat <- matrix(k, nrow = length(k), ncol = ncol(lambda))
    Su <- 1 / (1 + ratio^k_mat)
    haz <- (k_mat / lambda) * (ratio^(k_mat - 1)) / (1 + ratio^k_mat)
  } else {
    stop("Unknown family code: ", family_code, call. = FALSE)
  }

  Su <- pmin(pmax(Su, 1e-12), 1 - 1e-12)
  haz <- pmax(haz, 1e-12)

  list(Su = Su, haz = haz)
}

compute_mstu_draws <- function(family_code, mu_lat_mat, median_mat, draws) {
  if (family_code == "E") {
    mstu_mat <- median_mat / log(2)
  } else if (family_code == "W") {
    k <- exp(draws$rho_W)
    lambda <- sweep(median_mat, 1L, (log(2))^(1 / k), FUN = "/")
    k_mat <- matrix(k, nrow = length(k), ncol = ncol(lambda))
    mstu_mat <- lambda * gamma(1 + 1 / k_mat)
  } else if (family_code == "LN") {
    sigma <- exp(draws$log_sigma_LN)
    sigma_mat <- matrix(sigma, nrow = length(sigma), ncol = ncol(mu_lat_mat))
    mstu_mat <- exp(mu_lat_mat + 0.5 * sigma_mat^2)
  } else if (family_code == "LL") {
    k <- exp(-draws$psi_LL)
    k_mat <- matrix(k, nrow = length(k), ncol = ncol(median_mat))
    lambda <- median_mat
    safe_mean <- ifelse(
      k_mat > 1,
      lambda * (pi / k_mat) / sin(pi / k_mat),
      Inf
    )
    mstu_mat <- safe_mean
  } else {
    stop("Unknown family code: ", family_code, call. = FALSE)
  }

  rowMeans(mstu_mat, na.rm = TRUE)
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

simulate_prior_predictive <- function(design_bundle, model_row, prior_artifacts, n_draws, horizons_eval) {
  K_inc <- ncol(design_bundle$X_inc)
  K_lat <- ncol(design_bundle$X_lat)

  beta_inc <- matrix(NA_real_, nrow = n_draws, ncol = K_inc)
  for (j in seq_len(K_inc)) {
    if (design_bundle$inc_is_site[[j]] == 1L && design_bundle$site_prior_family_id %in% c(1L, 2L)) {
      if (design_bundle$site_prior_family_id == 1L) {
        beta_inc[, j] <- rnorm(n_draws, 0, 1)
      } else {
        beta_inc[, j] <- rt(n_draws, df = 3)
      }
    } else {
      beta_inc[, j] <- rnorm(
        n_draws,
        mean = design_bundle$mu_beta_inc[[j]],
        sd = design_bundle$sd_beta_inc[[j]]
      )
    }
  }

  gamma_lat <- matrix(NA_real_, nrow = n_draws, ncol = K_lat)
  for (j in seq_len(K_lat)) {
    if (design_bundle$lat_is_site[[j]] == 1L && design_bundle$site_prior_family_id %in% c(1L, 2L)) {
      if (design_bundle$site_prior_family_id == 1L) {
        gamma_lat[, j] <- rnorm(n_draws, 0, 1)
      } else {
        gamma_lat[, j] <- rt(n_draws, df = 3)
      }
    } else {
      gamma_lat[, j] <- rnorm(n_draws, 0, design_bundle$sd_gamma_lat[[j]])
    }
  }

  sim_draws <- list(
    alpha_inc = rnorm(n_draws, design_bundle$alpha_inc_prior_mean, design_bundle$alpha_inc_prior_sd),
    beta_inc = beta_inc,
    gamma0 = rnorm(n_draws, 0, design_bundle$sd_gamma0),
    gamma_lat = gamma_lat,
    rho_W = rnorm(n_draws, 0, prior_artifacts$sd_shape_W),
    log_sigma_LN = rnorm(n_draws, 0, prior_artifacts$sd_log_sigma_LN),
    psi_LL = rnorm(n_draws, 0, prior_artifacts$sd_psi_LL)
  )

  linear_terms <- compute_linear_terms(sim_draws, design_bundle$X_inc, design_bundle$X_lat)

  supported_horizons <- if (model_row$dataset_key[[1L]] == "PNU") c(1L, 2L) else c(1L, 2L, 5L)
  supported_risk_list <- lapply(supported_horizons, function(h) {
    fh <- family_survival_hazard(h, model_row$family_code[[1L]], linear_terms$mu_lat_mat, linear_terms$median_mat, sim_draws)
    1 - ((1 - linear_terms$pi_mat) + linear_terms$pi_mat * fh$Su)
  })

  degeneracy <- compute_degeneracy(linear_terms$pi_mat, linear_terms$median_mat, supported_risk_list)

  prior_rows <- list(
    tibble(
      dataset_key = model_row$dataset_key[[1L]],
      model_id = model_row$model_id[[1L]],
      branch = stage8_branch_label,
      risk_scale = stage8_risk_scale,
      prior_branch = model_row$prior_branch[[1L]],
      site_prior_family = model_row$site_prior_family_export[[1L]],
      metric = "cohort_mean_susceptible_probability",
      horizon = NA_integer_,
      summary_scalar(rowMeans(linear_terms$pi_mat)),
      prior_degenerate_flag = degeneracy$degenerate_flag[[1L]]
    ),
    tibble(
      dataset_key = model_row$dataset_key[[1L]],
      model_id = model_row$model_id[[1L]],
      branch = stage8_branch_label,
      risk_scale = stage8_risk_scale,
      prior_branch = model_row$prior_branch[[1L]],
      site_prior_family = model_row$site_prior_family_export[[1L]],
      metric = "cohort_mean_cure_probability",
      horizon = NA_integer_,
      summary_scalar(rowMeans(linear_terms$cure_prob_mat)),
      prior_degenerate_flag = degeneracy$degenerate_flag[[1L]]
    ),
    tibble(
      dataset_key = model_row$dataset_key[[1L]],
      model_id = model_row$model_id[[1L]],
      branch = stage8_branch_label,
      risk_scale = stage8_risk_scale,
      prior_branch = model_row$prior_branch[[1L]],
      site_prior_family = model_row$site_prior_family_export[[1L]],
      metric = "cohort_median_susceptible_time",
      horizon = NA_integer_,
      summary_scalar(apply(linear_terms$median_mat, 1L, stats::median)),
      prior_degenerate_flag = degeneracy$degenerate_flag[[1L]]
    )
  )

  risk_rows <- lapply(horizons_eval, function(h) {
    fh <- family_survival_hazard(h, model_row$family_code[[1L]], linear_terms$mu_lat_mat, linear_terms$median_mat, sim_draws)
    risk_mat <- 1 - ((1 - linear_terms$pi_mat) + linear_terms$pi_mat * fh$Su)
    tibble(
      dataset_key = model_row$dataset_key[[1L]],
      model_id = model_row$model_id[[1L]],
      branch = stage8_branch_label,
      risk_scale = stage8_risk_scale,
      prior_branch = model_row$prior_branch[[1L]],
      site_prior_family = model_row$site_prior_family_export[[1L]],
      metric = "cohort_mean_risk",
      horizon = as.integer(h),
      summary_scalar(rowMeans(risk_mat)),
      prior_degenerate_flag = degeneracy$degenerate_flag[[1L]]
    )
  })

  bind_rows(prior_rows, risk_rows)
}

annotate_prior_predictive_summary <- function(df) {
  if (nrow(df) == 0L) {
    return(df)
  }

  df %>%
    mutate(
      prior_tail_warning_flag = metric == "cohort_median_susceptible_time" &
        ((safe_numeric(mean) > prior_tail_warning_mean_years) | (safe_numeric(q975) > prior_tail_warning_q975_years)),
      prior_tail_warning_detail = case_when(
        metric == "cohort_median_susceptible_time" &
          safe_numeric(mean) > prior_tail_warning_mean_years &
          safe_numeric(q975) > prior_tail_warning_q975_years ~
            paste0("Prior-predictive median susceptible time exceeds both mean>", prior_tail_warning_mean_years, " and q975>", prior_tail_warning_q975_years, " years."),
        metric == "cohort_median_susceptible_time" &
          safe_numeric(mean) > prior_tail_warning_mean_years ~
            paste0("Prior-predictive median susceptible time mean exceeds ", prior_tail_warning_mean_years, " years."),
        metric == "cohort_median_susceptible_time" &
          safe_numeric(q975) > prior_tail_warning_q975_years ~
            paste0("Prior-predictive median susceptible time q975 exceeds ", prior_tail_warning_q975_years, " years."),
        TRUE ~ NA_character_
      )
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

  if (is.null(log_lik)) {
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
    if ("waic" %in% rownames(waic_obj$estimates)) {
      out$waic <- as.numeric(waic_obj$estimates["waic", "Estimate"])
    }
    if ("p_waic" %in% rownames(waic_obj$estimates)) {
      out$p_waic <- as.numeric(waic_obj$estimates["p_waic", "Estimate"])
    }
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
    if ("looic" %in% rownames(loo_obj$estimates)) {
      out$looic <- as.numeric(loo_obj$estimates["looic", "Estimate"])
    }
    if ("p_loo" %in% rownames(loo_obj$estimates)) {
      out$p_loo <- as.numeric(loo_obj$estimates["p_loo", "Estimate"])
    }
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

compute_classification_summary <- function(risk_draws, horizon_row, thresholds, cohort_n) {
  prevalence <- as.numeric(horizon_row$observed_km_risk)
  w_case <- unlist(horizon_row$w_case)
  w_control <- unlist(horizon_row$w_control)
  denom_case <- as.numeric(horizon_row$denom_case)
  denom_control <- as.numeric(horizon_row$denom_control)

  out_rows <- vector("list", length(thresholds))
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
    fp_count_draw <- cohort_n * fp_burden_draw
    nb_draw <- prevalence * tpr_draw - (1 - prevalence) * fpr_draw * (thr / (1 - thr))

    out_rows[[jj]] <- tibble(
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
      NB_q975 = stats::quantile(nb_draw, 0.975, na.rm = TRUE, names = FALSE),
      classification_estimable_flag = denom_case > 0 & denom_control > 0
    )
  }

  bind_rows(out_rows)
}

classify_hazard_shape <- function(hazard_values) {
  h <- as.numeric(hazard_values)
  if (length(h) < 2L || any(!is.finite(h))) {
    return("unstable_or_nonfinite")
  }

  d <- diff(h)
  if (all(d <= 0)) {
    return("monotone_decreasing")
  }
  if (all(d >= 0)) {
    return("monotone_increasing")
  }

  peak <- which.max(h)
  trough <- which.min(h)

  if (peak > 1L && peak < length(h) && all(diff(h[1:peak]) >= 0) && all(diff(h[peak:length(h)]) <= 0)) {
    return("unimodal")
  }
  if (trough > 1L && trough < length(h) && all(diff(h[1:trough]) <= 0) && all(diff(h[trough:length(h)]) >= 0)) {
    return("u_shaped")
  }

  "complex_nonmonotone"
}

compile_stage8a_stan_model <- function() {
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
  array[K_inc] int<lower=0, upper=1> inc_is_site;
  int<lower=1> K_lat;
  matrix[N, K_lat] X_lat;
  array[K_lat] int<lower=0, upper=1> lat_is_site;
  int<lower=1, upper=4> family_id;
  int<lower=0, upper=2> site_prior_family_id;
  real alpha_inc_prior_mean;
  real<lower=0> alpha_inc_prior_sd;
  vector[K_inc] mu_beta_inc;
  vector<lower=0>[K_inc] sd_beta_inc;
  real<lower=0> sd_gamma0;
  vector<lower=0>[K_lat] sd_gamma_lat;
  real<lower=0> sd_shape_W;
  real<lower=0> sd_log_sigma_LN;
  real<lower=0> sd_psi_LL;
}
parameters {
  real alpha_inc;
  vector[K_inc] beta_inc;
  real gamma0;
  vector[K_lat] gamma_lat;
  real rho_W;
  real log_sigma_LN;
  real psi_LL;
}
model {
  alpha_inc ~ normal(alpha_inc_prior_mean, alpha_inc_prior_sd);

  for (j in 1:K_inc) {
    if (inc_is_site[j] == 1) {
      if (site_prior_family_id == 2) {
        beta_inc[j] ~ student_t(3, 0, 1);
      } else {
        beta_inc[j] ~ normal(0, 1);
      }
    } else {
      beta_inc[j] ~ normal(mu_beta_inc[j], sd_beta_inc[j]);
    }
  }

  gamma0 ~ normal(0, sd_gamma0);
  for (j in 1:K_lat) {
    if (lat_is_site[j] == 1) {
      if (site_prior_family_id == 2) {
        gamma_lat[j] ~ student_t(3, 0, 1);
      } else {
        gamma_lat[j] ~ normal(0, 1);
      }
    } else {
      gamma_lat[j] ~ normal(0, sd_gamma_lat[j]);
    }
  }

  rho_W ~ normal(0, sd_shape_W);
  log_sigma_LN ~ normal(0, sd_log_sigma_LN);
  psi_LL ~ normal(0, sd_psi_LL);

  for (i in 1:N) {
    real eta_inc;
    real pi_i;
    real mu_lat;
    real m_i;
    real logS;
    real logf;

    eta_inc = alpha_inc + dot_product(row(X_inc, i), beta_inc);
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

    eta_inc = alpha_inc + dot_product(row(X_inc, i), beta_inc);
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
  rstan::stan_model(model_code = stan_code, model_name = "stage8a_bayesian_transition_only_cure")
}

extract_draws_compact <- function(fit, K_inc, K_lat) {
  ext <- rstan::extract(
    fit,
    pars = c("alpha_inc", "beta_inc", "gamma0", "gamma_lat", "rho_W", "log_sigma_LN", "psi_LL", "log_lik"),
    permuted = TRUE,
    inc_warmup = FALSE
  )

  beta_inc <- ext$beta_inc
  gamma_lat <- ext$gamma_lat
  log_lik <- ext$log_lik

  if (is.null(dim(beta_inc))) {
    beta_inc <- matrix(beta_inc, ncol = K_inc)
  }
  if (is.null(dim(gamma_lat))) {
    gamma_lat <- matrix(gamma_lat, ncol = K_lat)
  }
  if (is.null(dim(log_lik))) {
    log_lik <- matrix(log_lik, nrow = length(ext$alpha_inc))
  }

  list(
    alpha_inc = as.numeric(ext$alpha_inc),
    beta_inc = beta_inc,
    gamma0 = as.numeric(ext$gamma0),
    gamma_lat = gamma_lat,
    rho_W = as.numeric(ext$rho_W),
    log_sigma_LN = as.numeric(ext$log_sigma_LN),
    psi_LL = as.numeric(ext$psi_LL),
    log_lik = log_lik
  )
}


# 🔴 Define: diagnostic plotting, prior-tail logic, and delta builders ===============================
format_stage8_progress <- function(done, total) {
  sprintf("[%d/%d | %5.1f%%]", as.integer(done), as.integer(total), 100 * as.numeric(done) / as.numeric(total))
}

emit_stage8_progress <- function(done, total, model_id, detail) {
  message(format_stage8_progress(done, total), " ", model_id, " ", detail)
}

select_trace_parameters <- function(family_code, K_inc, K_lat) {
  selected_pars <- c("alpha_inc", "gamma0")
  if (K_inc >= 1L) {
    selected_pars <- c(selected_pars, "beta_inc[1]")
  }
  if (K_lat >= 1L) {
    selected_pars <- c(selected_pars, "gamma_lat[1]")
  }
  if (family_code == "W") selected_pars <- c(selected_pars, "rho_W")
  if (family_code == "LN") selected_pars <- c(selected_pars, "log_sigma_LN")
  if (family_code == "LL") selected_pars <- c(selected_pars, "psi_LL")
  unique(selected_pars)
}

make_trace_record <- function(fit, model_id, family_code, K_inc, K_lat) {
  selected_pars <- select_trace_parameters(family_code, K_inc, K_lat)
  list(
    model_id = model_id,
    selected_pars = selected_pars,
    arr = as.array(fit, pars = selected_pars)
  )
}

plot_trace_record <- function(trace_record) {
  selected_pars <- trace_record$selected_pars
  arr <- trace_record$arr
  model_id <- trace_record$model_id

  par_old <- par(no.readonly = TRUE)
  on.exit(par(par_old), add = TRUE)
  par(mfrow = c(2, 2), mar = c(3, 3, 3, 1))

  for (jj in seq_len(min(length(selected_pars), 4L))) {
    matplot(
      arr[, , jj],
      type = "l",
      lty = 1,
      col = seq_len(dim(arr)[2]),
      main = paste(model_id, selected_pars[[jj]]),
      xlab = "Iteration",
      ylab = ""
    )
  }
  if (length(selected_pars) < 4L) {
    for (jj in seq_len(4L - length(selected_pars))) {
      plot.new()
    }
  }
}

safe_generate_diagnostic_pdf <- function(trace_records, cohort_df, class_df, ppc_df, final_path) {
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
    for (one in trace_records) {
      plot_trace_record(one)
    }
  }

  if (nrow(cohort_df) > 0L) {
    risk_plot_df <- cohort_df %>%
      filter(
        is.finite(horizon),
        is.finite(risk_mean),
        is.finite(risk_q025),
        is.finite(risk_q975)
      )

    if (nrow(risk_plot_df) > 0L) {
      g_risk <- risk_plot_df %>%
      ggplot(aes(x = horizon, y = risk_mean, color = model_id, fill = model_id)) +
      geom_ribbon(aes(ymin = risk_q025, ymax = risk_q975), alpha = 0.15, linewidth = 0) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~ dataset_key, scales = "free_y") +
      labs(
        title = "Stage 8A posterior cohort mean risk trajectories",
        x = "Horizon (years)",
        y = "Posterior mean risk"
      ) +
        theme_bw() +
        theme(legend.position = "none")
      print(g_risk)
    }

    hazard_plot_df <- cohort_df %>%
      filter(
        is.finite(horizon),
        is.finite(hazard_mean),
        is.finite(hazard_q025),
        is.finite(hazard_q975)
      )

    if (nrow(hazard_plot_df) > 0L) {
      g_haz <- hazard_plot_df %>%
      ggplot(aes(x = horizon, y = hazard_mean, color = model_id, fill = model_id)) +
      geom_ribbon(aes(ymin = hazard_q025, ymax = hazard_q975), alpha = 0.15, linewidth = 0) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~ dataset_key, scales = "free_y") +
      labs(
        title = "Stage 8A posterior cohort mean hazard trajectories",
        x = "Horizon (years)",
        y = "Posterior mean hazard"
      ) +
        theme_bw() +
        theme(legend.position = "none")
      print(g_haz)
    }
  }

  if (nrow(class_df) > 0L) {
    nb_plot_df <- class_df %>%
      filter(horizon %in% c(1L, 2L, 5L)) %>%
      filter(
        is.finite(threshold),
        is.finite(NB_mean),
        is.finite(NB_q025),
        is.finite(NB_q975)
      )

    if (nrow(nb_plot_df) > 0L) {
      g_nb <- nb_plot_df %>%
      ggplot(aes(x = threshold, y = NB_mean, color = model_id, fill = model_id)) +
      geom_ribbon(aes(ymin = NB_q025, ymax = NB_q975), alpha = 0.15, linewidth = 0) +
      geom_line(linewidth = 0.7) +
      facet_grid(dataset_key ~ horizon, scales = "free_y") +
      labs(
        title = "Stage 8A net benefit by threshold",
        x = "Risk threshold",
        y = "Net benefit"
      ) +
        theme_bw() +
        theme(legend.position = "none")
      print(g_nb)
    }
  }

  if (nrow(ppc_df) > 0L) {
    ppc_plot_df <- ppc_df %>%
      filter(
        is.finite(horizon),
        is.finite(posterior_mean_risk),
        is.finite(posterior_q025_risk),
        is.finite(posterior_q975_risk)
      )

    if (nrow(ppc_plot_df) > 0L) {
      g_ppc <- ppc_plot_df %>%
      ggplot(aes(x = horizon, y = posterior_mean_risk, color = model_id)) +
      geom_errorbar(aes(ymin = posterior_q025_risk, ymax = posterior_q975_risk), width = 0.12, alpha = 0.6) +
      geom_line(linewidth = 0.6) +
      geom_point(size = 1.2) +
      geom_point(
        data = ppc_plot_df %>% filter(is.finite(observed_km_risk)),
        aes(x = horizon, y = observed_km_risk),
        shape = 4,
        size = 2,
        stroke = 0.9,
        color = "black",
        inherit.aes = FALSE
      ) +
      facet_wrap(~ dataset_key, scales = "free_y") +
      labs(
        title = "Stage 8A posterior predictive checks against observed KM risk",
        x = "Horizon (years)",
        y = "Risk"
      ) +
        theme_bw() +
        theme(legend.position = "none")
      print(g_ppc)
    }
  }

  if (length(trace_records) == 0L && nrow(cohort_df) == 0L && nrow(class_df) == 0L && nrow(ppc_df) == 0L) {
    plot.new()
    text(0.5, 0.5, "No diagnostic pages were available for this Stage 8A run.")
  }

  close_pdf()
  safe_promote_file(tmp_pdf, final_path)
  invisible(TRUE)
}

make_threshold_key <- function(x) {
  out <- rep("__NA__", length(x))
  idx <- !is.na(x)
  out[idx] <- sprintf("%.10f", as.numeric(x[idx]))
  out
}

compute_prior_tail_flags <- function(anchor_vs_neutral_delta_panel) {
  if (nrow(anchor_vs_neutral_delta_panel) == 0L) {
    return(empty_prior_tail_by_horizon())
  }

  safe_max_abs <- function(x) {
    x <- abs(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0L) {
      return(0)
    }
    max(x)
  }

  anchor_vs_neutral_delta_panel %>%
    group_by(dataset_key, structural_model_id, family_code, site_prior_family, horizon) %>%
    summarise(
      max_abs_delta_risk = safe_max_abs(delta_risk_anchor_minus_neutral),
      max_abs_delta_cure_fraction = safe_max_abs(delta_cure_fraction_anchor_minus_neutral),
      max_abs_delta_false_positive_burden = safe_max_abs(delta_false_positive_burden_anchor_minus_neutral),
      max_abs_delta_FP100 = safe_max_abs(delta_FP100_anchor_minus_neutral),
      max_abs_delta_NB = safe_max_abs(delta_NB_anchor_minus_neutral),
      max_abs_delta_PPV = safe_max_abs(delta_PPV_anchor_minus_neutral),
      max_abs_delta_TPR = safe_max_abs(delta_TPR_anchor_minus_neutral),
      nb_sign_changed = any(anchor_NB_mean * neutral_NB_mean < 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      prior_tail_sensitive = (
        max_abs_delta_risk >= materiality_delta_risk |
        max_abs_delta_cure_fraction >= materiality_delta_cure_fraction |
        max_abs_delta_false_positive_burden >= materiality_delta_false_positive_burden |
        max_abs_delta_FP100 >= materiality_delta_FP100 |
        max_abs_delta_NB >= materiality_delta_NB |
        max_abs_delta_PPV >= materiality_delta_PPV |
        max_abs_delta_TPR >= materiality_delta_TPR |
        nb_sign_changed
      )
    ) %>%
    select(dataset_key, structural_model_id, family_code, site_prior_family, horizon, prior_tail_sensitive)
}

empty_anchor_vs_neutral_delta_panel <- function() {
  tibble(
    dataset_key = character(),
    structural_model_id = character(),
    family_code = character(),
    site_prior_family = character(),
    site_placement = character(),
    branch = character(),
    risk_scale = character(),
    horizon = integer(),
    support_tier = character(),
    horizon_evidence_class = character(),
    claim_restriction_flag = character(),
    anchor_risk_mean = double(),
    anchor_cure_fraction_mean = double(),
    neutral_risk_mean = double(),
    neutral_cure_fraction_mean = double(),
    delta_risk_anchor_minus_neutral = double(),
    delta_cure_fraction_anchor_minus_neutral = double(),
    threshold = double(),
    anchor_false_positive_burden_mean = double(),
    neutral_false_positive_burden_mean = double(),
    delta_false_positive_burden_anchor_minus_neutral = double(),
    anchor_FP100_mean = double(),
    neutral_FP100_mean = double(),
    delta_FP100_anchor_minus_neutral = double(),
    anchor_NB_mean = double(),
    neutral_NB_mean = double(),
    delta_NB_anchor_minus_neutral = double(),
    anchor_PPV_mean = double(),
    neutral_PPV_mean = double(),
    delta_PPV_anchor_minus_neutral = double(),
    anchor_TPR_mean = double(),
    neutral_TPR_mean = double(),
    delta_TPR_anchor_minus_neutral = double(),
    prior_tail_sensitive = logical()
  )
}

empty_prior_tail_by_horizon <- function() {
  tibble(
    dataset_key = character(),
    structural_model_id = character(),
    family_code = character(),
    site_prior_family = character(),
    horizon = integer(),
    prior_tail_sensitive = logical()
  )
}

empty_incidence_anchor_update_panel <- function() {
  tibble(
    dataset_key = character(),
    branch = character(),
    risk_scale = character(),
    prior_branch = character(),
    retained_fit_id = character(),
    age_sex_anchor_cell = character(),
    external_incidence_rate_per10k = double(),
    external_one_year_risk = double(),
    prior_center_logit = double(),
    posterior_mean_logit = double(),
    posterior_lower_logit = double(),
    posterior_upper_logit = double(),
    posterior_mean_one_year_risk = double(),
    posterior_lower_one_year_risk = double(),
    posterior_upper_one_year_risk = double(),
    posterior_minus_prior_logit = double(),
    posterior_minus_prior_risk = double()
  )
}

empty_delta_vs_nocure <- function() {
  tibble(
    dataset_key = character(),
    model_id = character(),
    no_cure_model_id = character(),
    formula_anchor = character(),
    horizon = integer(),
    threshold = double(),
    metric = character(),
    value = double()
  )
}

build_anchor_vs_neutral_delta_panel <- function(cohort_df, class_df) {
  if (nrow(cohort_df) == 0L) {
    return(empty_anchor_vs_neutral_delta_panel())
  }

  cohort_anchor <- cohort_df %>%
    filter(prior_branch == anchor_prior_branch_label) %>%
    select(
      dataset_key, structural_model_id, family_code, site_prior_family, site_placement,
      branch, risk_scale, horizon, support_tier, horizon_evidence_class, claim_restriction_flag,
      anchor_risk_mean = risk_mean,
      anchor_cure_fraction_mean = cure_fraction_mean
    )

  cohort_neutral <- cohort_df %>%
    filter(prior_branch == neutral_prior_branch_label) %>%
    select(
      dataset_key, structural_model_id, family_code, site_prior_family, horizon,
      neutral_risk_mean = risk_mean,
      neutral_cure_fraction_mean = cure_fraction_mean
    )

  cohort_pair <- cohort_anchor %>%
    inner_join(
      cohort_neutral,
      by = c("dataset_key", "structural_model_id", "family_code", "site_prior_family", "horizon")
    ) %>%
    mutate(
      delta_risk_anchor_minus_neutral = anchor_risk_mean - neutral_risk_mean,
      delta_cure_fraction_anchor_minus_neutral = anchor_cure_fraction_mean - neutral_cure_fraction_mean
    )

  if (nrow(class_df) == 0L) {
    return(
      cohort_pair %>%
        mutate(
          threshold = NA_real_,
          anchor_false_positive_burden_mean = NA_real_,
          neutral_false_positive_burden_mean = NA_real_,
          delta_false_positive_burden_anchor_minus_neutral = NA_real_,
          anchor_FP100_mean = NA_real_,
          neutral_FP100_mean = NA_real_,
          delta_FP100_anchor_minus_neutral = NA_real_,
          anchor_NB_mean = NA_real_,
          neutral_NB_mean = NA_real_,
          delta_NB_anchor_minus_neutral = NA_real_,
          anchor_PPV_mean = NA_real_,
          neutral_PPV_mean = NA_real_,
          delta_PPV_anchor_minus_neutral = NA_real_,
          anchor_TPR_mean = NA_real_,
          neutral_TPR_mean = NA_real_,
          delta_TPR_anchor_minus_neutral = NA_real_
        )
    )
  }

  class_anchor <- class_df %>%
    filter(prior_branch == anchor_prior_branch_label) %>%
    select(
      dataset_key, structural_model_id, family_code, site_prior_family, horizon, threshold,
      anchor_false_positive_burden_mean = false_positive_burden_mean,
      anchor_FP100_mean = FP100_mean,
      anchor_NB_mean = NB_mean,
      anchor_PPV_mean = PPV_mean,
      anchor_TPR_mean = TPR_mean
    )

  class_neutral <- class_df %>%
    filter(prior_branch == neutral_prior_branch_label) %>%
    select(
      dataset_key, structural_model_id, family_code, site_prior_family, horizon, threshold,
      neutral_false_positive_burden_mean = false_positive_burden_mean,
      neutral_FP100_mean = FP100_mean,
      neutral_NB_mean = NB_mean,
      neutral_PPV_mean = PPV_mean,
      neutral_TPR_mean = TPR_mean
    )

  class_pair <- class_anchor %>%
    inner_join(
      class_neutral,
      by = c("dataset_key", "structural_model_id", "family_code", "site_prior_family", "horizon", "threshold")
    ) %>%
    mutate(
      delta_false_positive_burden_anchor_minus_neutral = anchor_false_positive_burden_mean - neutral_false_positive_burden_mean,
      delta_FP100_anchor_minus_neutral = anchor_FP100_mean - neutral_FP100_mean,
      delta_NB_anchor_minus_neutral = anchor_NB_mean - neutral_NB_mean,
      delta_PPV_anchor_minus_neutral = anchor_PPV_mean - neutral_PPV_mean,
      delta_TPR_anchor_minus_neutral = anchor_TPR_mean - neutral_TPR_mean
    )

  cohort_pair %>%
    left_join(
      class_pair,
      by = c("dataset_key", "structural_model_id", "family_code", "site_prior_family", "horizon")
    ) %>%
    arrange(dataset_key, structural_model_id, family_code, site_prior_family, horizon, threshold)
}

build_incidence_anchor_update_panel <- function(model_grid, retained_model_ids, fit_draw_map, prior_artifacts) {
  if (length(retained_model_ids) == 0L) {
    return(empty_incidence_anchor_update_panel())
  }

  out_rows <- list()
  retained_grid <- model_grid %>% filter(model_id %in% retained_model_ids)

  for (ii in seq_len(nrow(retained_grid))) {
    model_row <- retained_grid[ii, , drop = FALSE]
    fit_draws <- fit_draw_map[[model_row$model_id[[1L]]]]
    if (is.null(fit_draws)) {
      next
    }

    base_beta <- fit_draws$beta_inc[, seq_len(5L), drop = FALSE]
    alpha_draw <- fit_draws$alpha_inc

    for (jj in seq_len(nrow(prior_artifacts$anchor_cells))) {
      cell_row <- prior_artifacts$anchor_cells[jj, , drop = FALSE]
      x_cell <- c(
        cell_row$sex_num[[1L]],
        as.integer(cell_row$age_band[[1L]] == "20_29"),
        as.integer(cell_row$age_band[[1L]] == "30plus"),
        cell_row$sex_num[[1L]] * as.integer(cell_row$age_band[[1L]] == "20_29"),
        cell_row$sex_num[[1L]] * as.integer(cell_row$age_band[[1L]] == "30plus")
      )

      posterior_logit <- alpha_draw + as.vector(base_beta %*% x_cell)
      posterior_risk <- plogis(posterior_logit)

      prior_center_logit <- if (identical(model_row$prior_branch[[1L]], anchor_prior_branch_label)) {
        cell_row$prior_center_logit_anchor[[1L]]
      } else {
        0
      }

      out_rows[[length(out_rows) + 1L]] <- tibble(
        dataset_key = model_row$dataset_key[[1L]],
        branch = stage8_branch_label,
        risk_scale = stage8_risk_scale,
        prior_branch = model_row$prior_branch[[1L]],
        retained_fit_id = model_row$model_id[[1L]],
        age_sex_anchor_cell = cell_row$age_sex_anchor_cell[[1L]],
        external_incidence_rate_per10k = cell_row$external_incidence_rate_per10k[[1L]],
        external_one_year_risk = cell_row$external_one_year_risk[[1L]],
        prior_center_logit = prior_center_logit,
        posterior_mean_logit = mean(posterior_logit),
        posterior_lower_logit = stats::quantile(posterior_logit, 0.025, names = FALSE),
        posterior_upper_logit = stats::quantile(posterior_logit, 0.975, names = FALSE),
        posterior_mean_one_year_risk = mean(posterior_risk),
        posterior_lower_one_year_risk = stats::quantile(posterior_risk, 0.025, names = FALSE),
        posterior_upper_one_year_risk = stats::quantile(posterior_risk, 0.975, names = FALSE),
        posterior_minus_prior_logit = mean(posterior_logit) - prior_center_logit,
        posterior_minus_prior_risk = mean(posterior_risk) - plogis(prior_center_logit)
      )
    }
  }

  bind_rows_safe(out_rows)
}


build_delta_vs_nocure <- function(performance_long, nocure_cohort_long, nocure_class_long) {
  if (nrow(performance_long) == 0L) {
    return(empty_delta_vs_nocure())
  }

  cohort_rows <- performance_long %>%
    filter(is.na(threshold)) %>%
    select(dataset_key, model_id, formula_anchor, horizon, risk_mean) %>%
    distinct()

  nocure_cohort_ref <- nocure_cohort_long %>%
    filter(metric == "risk_mean") %>%
    select(dataset_key, formula_anchor, no_cure_model_id, horizon, metric, value) %>%
    distinct()

  nocure_class_ref <- nocure_class_long %>%
    select(dataset_key, formula_anchor, no_cure_model_id, horizon, threshold, metric, value) %>%
    distinct()

  out_rows <- list()

  if (nrow(nocure_cohort_ref) > 0L) {
    cohort_join <- cohort_rows %>%
      inner_join(
        nocure_cohort_ref,
        by = c("dataset_key", "formula_anchor", "horizon"),
        relationship = "many-to-many"
      ) %>%
      transmute(
        dataset_key = dataset_key,
        model_id = model_id,
        no_cure_model_id = no_cure_model_id,
        formula_anchor = formula_anchor,
        horizon = horizon,
        threshold = NA_real_,
        metric = "delta_risk_nocure_minus_stage8a",
        value = value - risk_mean
      )
    out_rows[[length(out_rows) + 1L]] <- cohort_join
  }

  if (nrow(nocure_class_ref) > 0L) {
    metric_map <- tribble(
      ~metric,                        ~stage8_col,
      "false_positive_burden_mean",   "false_positive_burden_mean",
      "FP100_mean",                   "FP100_mean",
      "NB_mean",                      "NB_mean",
      "PPV_mean",                     "PPV_mean",
      "TPR_mean",                     "TPR_mean"
    )

    for (ii in seq_len(nrow(metric_map))) {
      one <- metric_map[ii, , drop = FALSE]
      stage_col <- one$stage8_col[[1L]]

      stage8_sub <- performance_long %>%
        filter(!is.na(threshold)) %>%
        transmute(
          dataset_key = dataset_key,
          model_id = model_id,
          formula_anchor = formula_anchor,
          horizon = horizon,
          threshold = threshold,
          stage8_value = .data[[stage_col]]
        ) %>%
        distinct()

      nocure_sub <- nocure_class_ref %>%
        filter(metric == one$metric[[1L]]) %>%
        distinct()

      join_df <- stage8_sub %>%
        inner_join(
          nocure_sub,
          by = c("dataset_key", "formula_anchor", "horizon", "threshold"),
          relationship = "many-to-many"
        ) %>%
        transmute(
          dataset_key = dataset_key,
          model_id = model_id,
          no_cure_model_id = no_cure_model_id,
          formula_anchor = formula_anchor,
          horizon = horizon,
          threshold = threshold,
          metric = paste0("delta_", one$metric[[1L]], "_nocure_minus_stage8a"),
          value = value - stage8_value
        )

      out_rows[[length(out_rows) + 1L]] <- join_df
    }
  }

  bind_rows_safe(out_rows)
}

# 🔴 Define: metadata registries and schema-preserving placeholders ===============================
build_metadata_registry <- function() {
  tibble::tribble(
    ~metadata_group, ~metadata_name, ~metadata_value,
    "stage", "stage_name", "Stage 8A Bayesian transition-only cure",
    "stage", "branch", stage8_branch_label,
    "stage", "risk_scale", stage8_risk_scale,
    "stage", "main_event_definition", "status_num == 1",
    "stage", "main_censoring_definition", "status_num %in% c(0, 2)",
    "stage", "supplementary_competing_scale", stage8b_risk_scale,
    "stage", "pnu_site_label", pnu_site_label,
    "stage", "snu_site_label", snu_site_label,
    "stage", "horizon_vector", paste(horizons_year, collapse = ","),
    "stage", "threshold_vector", paste(risk_thresholds, collapse = ","),
    "stage", "prior_branch_anchor", anchor_prior_branch_label,
    "stage", "prior_branch_neutral", neutral_prior_branch_label,
    "stage", "site_prior_family_main", site_prior_family_main,
    "stage", "site_prior_family_sensitivity", site_prior_family_sensitivity,
    "priors", "alpha_gp_vdw", as.character(-9.581369553169),
    "priors", "mu_beta_inc_anchor", paste(c(0.419871845822, 0.907608052926, 0.586202561451, 0.466865123863, 0.037997248763), collapse = ","),
    "priors", "sd_beta_inc_anchor", paste(c(0.132789397422, 0.173731076538, 0.191221553945, 0.270393197518, 0.302838606651), collapse = ","),
    "priors", "neutral_beta_sd", "2",
    "priors", "latency_gamma0_sd", "2.5",
    "priors", "latency_gamma_sd", "1",
    "priors", "weibull_shape_sd", "0.35",
    "priors", "lognormal_scale_sd", "0.50",
    "priors", "loglogistic_shape_sd", "0.50",
    "materiality", "delta_risk", as.character(materiality_delta_risk),
    "materiality", "delta_cure_fraction", as.character(materiality_delta_cure_fraction),
    "materiality", "delta_false_positive_burden", as.character(materiality_delta_false_positive_burden),
    "materiality", "delta_FP100", as.character(materiality_delta_FP100),
    "materiality", "delta_NB", as.character(materiality_delta_NB),
    "materiality", "delta_PPV", as.character(materiality_delta_PPV),
    "materiality", "delta_TPR", as.character(materiality_delta_TPR),
    "mcmc", "stan_chains", as.character(stan_chains),
    "mcmc", "stan_iter", as.character(stan_iter),
    "mcmc", "stan_warmup", as.character(stan_warmup),
    "mcmc", "stan_seed", as.character(stan_seed),
    "mcmc", "stan_adapt_delta", as.character(stan_adapt_delta),
    "mcmc", "stan_max_treedepth", as.character(stan_max_treedepth),
    "diagnostics", "ess_min_threshold", as.character(ess_min_threshold),
    "diagnostics", "rhat_max_threshold", as.character(rhat_max_threshold),
    "diagnostics", "ppc_tolerance_abs", as.character(ppc_tolerance_abs),
    "inputs", "merged_data_path", merged_data_path,
    "inputs", "stage5_model_performance_csv", stage5_model_performance_csv,
    "inputs", "stage6_screening_flag_csv", stage6_screening_flag_csv,
    "outputs", "export_path", export_path,
    "documents", "governing_shared_master", "Integrated_Modeling_Master_Specification_English_REVISED_v5.md",
    "documents", "governing_stage8_companion", "Bayesian_Modeling_Specification_Stage8_REVISED_v5.md",
    "documents", "governing_code_rules", "Rules Before Generating R Code_🇬🇧ENG.md",
    "notes", "site_effect_interpretation", "Merged site terms are structural context proxies, not causal treatment effects.",
    "notes", "retention_rule", "Retain all admissible Stage 8A fits; do not force a single Bayesian winner.",
    "notes", "stage8b_delta_policy", "Stage 8A does not compute 8A-vs-8B deltas; Stage 8B owns that comparison."
  )
}

build_output_audit <- function(model_grid, model_registry, prediction_long, performance_long, reporting_metadata, anchor_delta, anchor_update, uncured_decomp, hazard_plausibility, diagnostic_pdf_path) {
  expected_models <- nrow(model_grid)
  actual_models <- nrow(model_registry)
  admissible_n <- sum(model_registry$admissible_flag %in% TRUE, na.rm = TRUE)

  checks <- list(
    tibble(
      check_name = "model_registry_row_count",
      status = ifelse(actual_models == expected_models, "pass", "fail"),
      observed_value = as.character(actual_models),
      expected_value = as.character(expected_models),
      detail = "One model registry row should exist for every Stage 8A grid row."
    ),
    tibble(
      check_name = "prediction_long_nonempty_when_admissible",
      status = ifelse(admissible_n == 0L || nrow(prediction_long) > 0L, "pass", "fail"),
      observed_value = as.character(nrow(prediction_long)),
      expected_value = ">0 if admissible fits exist",
      detail = "Admissible fits should populate the long-format prediction source table."
    ),
    tibble(
      check_name = "performance_long_nonempty_when_admissible",
      status = ifelse(admissible_n == 0L || nrow(performance_long) > 0L, "pass", "fail"),
      observed_value = as.character(nrow(performance_long)),
      expected_value = ">0 if admissible fits exist",
      detail = "Admissible fits should populate the long-format performance/classification table."
    ),
    tibble(
      check_name = "reporting_metadata_nonempty_when_admissible",
      status = ifelse(admissible_n == 0L || nrow(reporting_metadata) > 0L, "pass", "fail"),
      observed_value = as.character(nrow(reporting_metadata)),
      expected_value = ">0 if admissible fits exist",
      detail = "Reporting metadata should be directly joinable to performance rows."
    ),
    tibble(
      check_name = "anchor_update_panel_nonempty_when_admissible",
      status = ifelse(admissible_n == 0L || nrow(anchor_update) > 0L, "pass", "fail"),
      observed_value = as.character(nrow(anchor_update)),
      expected_value = ">0 if admissible fits exist",
      detail = "Incidence anchor update panel is mandatory for admissible Stage 8A fits."
    ),
    tibble(
      check_name = "uncured_decomposition_nonempty_when_admissible",
      status = ifelse(admissible_n == 0L || nrow(uncured_decomp) > 0L, "pass", "fail"),
      observed_value = as.character(nrow(uncured_decomp)),
      expected_value = ">0 if admissible fits exist",
      detail = "Uncured-only supporting decomposition should be exported for admissible fits."
    ),
    tibble(
      check_name = "hazard_plausibility_nonempty_when_admissible",
      status = ifelse(admissible_n == 0L || nrow(hazard_plausibility) > 0L, "pass", "fail"),
      observed_value = as.character(nrow(hazard_plausibility)),
      expected_value = ">0 if admissible fits exist",
      detail = "Hazard-shape plausibility export should be available for admissible fits."
    ),
    tibble(
      check_name = "diagnostic_pdf_exists",
      status = ifelse(pdf_file_is_usable(diagnostic_pdf_path), "pass", "fail"),
      observed_value = ifelse(pdf_file_is_usable(diagnostic_pdf_path), "TRUE", "FALSE"),
      expected_value = "TRUE",
      detail = "The aggregate Stage 8A diagnostic PDF should exist."
    ),
    tibble(
      check_name = "anchor_delta_panel_present",
      status = ifelse(admissible_n == 0L || nrow(anchor_delta) > 0L, "pass", "fail"),
      observed_value = as.character(nrow(anchor_delta)),
      expected_value = ">0 if paired anchor/neutral admissible fits exist",
      detail = "Anchor-versus-neutral delta table is mandatory when both prior branches are retained."
    )
  )

  bind_rows(checks)
}


# 🔴 Load: backbone data, carry-forward screening, and no-cure comparators ===============================
prior_artifacts <- build_prior_artifacts()

loaded_objects <- load_stage8a_datasets(
  merged_path = merged_data_path,
  pnu_label = pnu_site_label,
  snu_label = snu_site_label
)

analysis_datasets <- loaded_objects$datasets
dataset_registry <- loaded_objects$dataset_registry

screening_flags <- read_screening_flags(stage6_screening_flag_csv)
nocure_cohort_long <- normalize_nocure_cohort(stage5_model_performance_csv)
nocure_class_long <- normalize_nocure_classification(stage5_model_performance_csv)

ipcw_registry <- lapply(names(analysis_datasets), function(ds) {
  build_ipcw_reference(analysis_datasets[[ds]], horizons_year)
})
names(ipcw_registry) <- names(analysis_datasets)

horizon_metadata <- build_horizon_metadata(
  analysis_datasets = analysis_datasets,
  horizons = horizons_year,
  ipcw_registry = ipcw_registry
)

metadata_registry <- build_metadata_registry()

# 🔴 Prepare: model grid, Stage-6 lookup, and Stan compiler ===============================
model_grid <- build_model_grid() %>%
  mutate(
    site_prior_family = if_else(has_site_term, site_prior_family, "not_applicable"),
    site_prior_family_export = site_prior_family
  )

if (!is.null(run_model_ids)) {
  model_grid <- model_grid %>% filter(model_id %in% run_model_ids)
  if (nrow(model_grid) == 0L) {
    stop("`run_model_ids` filtered out all Stage 8A models.", call. = FALSE)
  }
}

screening_lookup <- build_screening_model_lookup(screening_flags, model_grid)
model_grid$stage8a_rds_path <- if (isTRUE(save_model_rds)) {
  stage8a_model_rds_path(model_grid$model_id, export_path)
} else {
  NA_character_
}
model_grid$reuse_existing_fit <- FALSE
stage8a_reuse_checks <- vector("list", nrow(model_grid))

if (isTRUE(reuse_existing_rds) && isTRUE(save_model_rds)) {
  for (ii in seq_len(nrow(model_grid))) {
    stage8a_reuse_checks[[ii]] <- load_existing_stage8a_bundle(
      path = model_grid$stage8a_rds_path[[ii]],
      model_row = model_grid[ii, , drop = FALSE],
      dataset_df = analysis_datasets[[model_grid$dataset_key[[ii]]]]
    )
    model_grid$reuse_existing_fit[[ii]] <- is_stage8a_bundle_reusable(stage8a_reuse_checks[[ii]])
  }
}

n_models_reused <- sum(model_grid$reuse_existing_fit %in% TRUE)
n_models_to_fit <- sum(!(model_grid$reuse_existing_fit %in% TRUE))
message("Stage8A reuse plan: ", n_models_reused, " model(s) reused; ", n_models_to_fit, " model(s) require fitting.")

stan_model_compiled <- if (n_models_to_fit > 0L) compile_stage8a_stan_model() else NULL

trace_records <- list()
fit_draw_map <- list()

model_registry_rows <- list()
coefficient_rows <- list()
diagnostic_parameter_rows <- list()
ppc_rows <- list()
prior_predictive_rows <- list()
prediction_rows <- list()
cohort_rows <- list()
classification_rows <- list()
uncured_rows <- list()
hazard_rows <- list()

total_models <- nrow(model_grid)
model_order <- model_grid$model_id

# 🔴 Run: Stage-8A Bayesian model grid with anchor-vs-neutral priors ===============================
for (ii in seq_len(total_models)) {
  model_row <- model_grid[ii, , drop = FALSE]
  model_id_now <- model_row$model_id[[1L]]
  dataset_name <- model_row$dataset_key[[1L]]
  dataset_df <- analysis_datasets[[dataset_name]]
  model_started_at <- Sys.time()

  emit_stage8_progress(
    ii - 1L,
    total_models,
    model_id_now,
    paste0(
      "starting fit (dataset=", dataset_name,
      ", family=", model_row$family_code[[1L]],
      ", prior=", model_row$prior_branch[[1L]],
      ", site_prior=", model_row$site_prior_family_export[[1L]], ")"
    )
  )

  screening_row <- screening_lookup %>%
    filter(model_id == model_id_now)

  existing_rds_path <- model_row$stage8a_rds_path[[1L]] %||% NA_character_
  existing_bundle_check <- stage8a_reuse_checks[[ii]] %||%
    default_stage8a_bundle_check(if (isTRUE(save_model_rds)) "reuse_not_precomputed" else "rds_saving_disabled")

  if (is_stage8a_bundle_reusable(existing_bundle_check)) {
    reused_registry <- tibble::as_tibble(existing_bundle_check$registry) %>%
      mutate(
        dataset_key = dataset_name,
        model_id = model_id_now,
        branch = stage8_branch_label,
        risk_scale = stage8_risk_scale,
        cure_model_eligibility_flag = screening_value_or_na(screening_row, "cure_model_eligibility_flag"),
        primary_gate_method = screening_value_or_na(screening_row, "primary_gate_method"),
        primary_gate_flag = screening_value_or_na(screening_row, "primary_gate_flag"),
        receus_primary_class = screening_value_or_na(screening_row, "receus_primary_class"),
        presence_modifier_flag = screening_value_or_na(screening_row, "presence_modifier_flag"),
        cure_presence_support_flag = screening_value_or_na(screening_row, "cure_presence_support_flag"),
        followup_contradiction_flag = screening_value_or_na(screening_row, "followup_contradiction_flag"),
        followup_not_contradicted_flag = screening_value_or_na(screening_row, "followup_not_contradicted_flag"),
        screening_note = screening_value_or_na(screening_row, "screening_note"),
        carry_forward_stage8 = screening_logical_or_na(screening_row, "carry_forward_stage8"),
        reused_existing_rds_flag = TRUE,
        rds_validation_status = "existing_rds_skipped",
        rds_validation_reason = existing_bundle_check$reason,
        rds_path = existing_rds_path
      )

    model_registry_rows[[length(model_registry_rows) + 1L]] <- reused_registry
    coefficient_rows[[length(coefficient_rows) + 1L]] <- bundle_table_or_empty(existing_bundle_check$bundle, "coefficient_summary", model_id_now)
    diagnostic_parameter_rows[[length(diagnostic_parameter_rows) + 1L]] <- bundle_table_or_empty(existing_bundle_check$bundle, "diagnostics_parameter_level", model_id_now)
    ppc_rows[[length(ppc_rows) + 1L]] <- bundle_table_or_empty(existing_bundle_check$bundle, "ppc_summary", model_id_now)
    prior_predictive_rows[[length(prior_predictive_rows) + 1L]] <- bundle_table_or_empty(existing_bundle_check$bundle, "prior_predictive_summary", model_id_now)

    if (isTRUE(existing_bundle_check$admissible)) {
      prediction_rows[[length(prediction_rows) + 1L]] <- bundle_table_or_empty(existing_bundle_check$bundle, "prediction_long", model_id_now)
      cohort_rows[[length(cohort_rows) + 1L]] <- bundle_table_or_empty(existing_bundle_check$bundle, "cohort_yearly", model_id_now)
      classification_rows[[length(classification_rows) + 1L]] <- bundle_table_or_empty(existing_bundle_check$bundle, "classification", model_id_now)
      uncured_rows[[length(uncured_rows) + 1L]] <- bundle_table_or_empty(existing_bundle_check$bundle, "uncured_decomposition", model_id_now)
      hazard_rows[[length(hazard_rows) + 1L]] <- bundle_table_or_empty(existing_bundle_check$bundle, "hazard_plausibility", model_id_now)

      inc_draws <- existing_bundle_check$bundle[["posterior_incidence_draws"]]
      if (is.list(inc_draws)) {
        fit_draw_map[[model_id_now]] <- inc_draws
      }
    }

    emit_stage8_progress(ii, total_models, model_id_now, paste0("existing RDS detected; skipped refit: ", basename(existing_rds_path)))
    next
  }

  if (isTRUE(existing_bundle_check$found)) {
    emit_stage8_progress(
      ii - 1L,
      total_models,
      model_id_now,
      paste0("existing RDS found but refit required (reason=", existing_bundle_check$reason, ")")
    )
  }

  if (is.null(stan_model_compiled)) {
    stop("Internal error: Stage 8A Stan model was not compiled for a model that requires refitting.", call. = FALSE)
  }

  design_bundle <- make_design_bundle(
    df = dataset_df,
    model_row = model_row,
    prior_artifacts = prior_artifacts,
    snu_label = snu_site_label
  )

  set.seed(stan_seed + ii)
  prior_predictive_rows[[length(prior_predictive_rows) + 1L]] <- simulate_prior_predictive(
    design_bundle = design_bundle,
    model_row = model_row,
    prior_artifacts = prior_artifacts,
    n_draws = prior_predictive_draws,
    horizons_eval = c(1L, 2L, 5L, 10L)
  )

  stan_data <- list(
    N = nrow(dataset_df),
    time = as.numeric(design_bundle$time),
    event = as.integer(design_bundle$event),
    K_inc = ncol(design_bundle$X_inc),
    X_inc = design_bundle$X_inc,
    inc_is_site = as.integer(design_bundle$inc_is_site),
    K_lat = ncol(design_bundle$X_lat),
    X_lat = design_bundle$X_lat,
    lat_is_site = as.integer(design_bundle$lat_is_site),
    family_id = as.integer(model_row$family_id[[1L]]),
    site_prior_family_id = as.integer(design_bundle$site_prior_family_id),
    alpha_inc_prior_mean = as.numeric(design_bundle$alpha_inc_prior_mean),
    alpha_inc_prior_sd = as.numeric(design_bundle$alpha_inc_prior_sd),
    mu_beta_inc = as.numeric(design_bundle$mu_beta_inc),
    sd_beta_inc = as.numeric(design_bundle$sd_beta_inc),
    sd_gamma0 = as.numeric(design_bundle$sd_gamma0),
    sd_gamma_lat = as.numeric(design_bundle$sd_gamma_lat),
    sd_shape_W = as.numeric(design_bundle$sd_shape_W),
    sd_log_sigma_LN = as.numeric(design_bundle$sd_log_sigma_LN),
    sd_psi_LL = as.numeric(design_bundle$sd_psi_LL)
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

    model_registry_rows[[length(model_registry_rows) + 1L]] <- tibble(
      dataset_key = dataset_name,
      model_id = model_id_now,
      retained_fit_id = model_id_now,
      branch = stage8_branch_label,
      risk_scale = stage8_risk_scale,
      structural_model_id = model_row$structural_model_id[[1L]],
      formula_anchor = model_row$formula_anchor[[1L]],
      family_code = model_row$family_code[[1L]],
      latency_family = model_row$latency_family[[1L]],
      prior_branch = model_row$prior_branch[[1L]],
      site_prior_family = model_row$site_prior_family_export[[1L]],
      site_placement = model_row$site_placement[[1L]],
      fit_status = fit_status,
      fit_error_message = fit_error_message,
      admissible_flag = FALSE,
      admissibility_reasons = "sampling_error",
      divergences = NA_integer_,
      treedepth_exceeded = NA_integer_,
      max_rhat = NA_real_,
      min_bulk_ess = NA_real_,
      min_tail_ess = NA_real_,
      waic = NA_real_,
      p_waic = NA_real_,
      looic = NA_real_,
      p_loo = NA_real_,
      p_waic_high_n = NA_integer_,
      p_waic_high_pct = NA_real_,
      pareto_k_max = NA_real_,
      pareto_k_bad_n = NA_integer_,
      pareto_k_bad_pct = NA_real_,
      pareto_k_very_bad_n = NA_integer_,
      waic_warning_flag = NA,
      loo_warning_flag = NA,
      info_criteria_warning_detail = NA_character_,
      prior_degenerate_flag = NA,
      posterior_degenerate_flag = NA,
      ppc_gross_contradiction_flag = NA,
      n = nrow(dataset_df),
      n_transition = sum(dataset_df$event_main),
      n_remission = sum(dataset_df$remission_flag),
      n_main_censoring = sum(dataset_df$censor_main),
      cure_model_eligibility_flag = screening_value_or_na(screening_row, "cure_model_eligibility_flag"),
      primary_gate_method = screening_value_or_na(screening_row, "primary_gate_method"),
      primary_gate_flag = screening_value_or_na(screening_row, "primary_gate_flag"),
      receus_primary_class = screening_value_or_na(screening_row, "receus_primary_class"),
      presence_modifier_flag = screening_value_or_na(screening_row, "presence_modifier_flag"),
      cure_presence_support_flag = screening_value_or_na(screening_row, "cure_presence_support_flag"),
      followup_contradiction_flag = screening_value_or_na(screening_row, "followup_contradiction_flag"),
      followup_not_contradicted_flag = screening_value_or_na(screening_row, "followup_not_contradicted_flag"),
      screening_note = screening_value_or_na(screening_row, "screening_note"),
      carry_forward_stage8 = screening_logical_or_na(screening_row, "carry_forward_stage8"),
      reused_existing_rds_flag = FALSE,
      rds_validation_status = "new_fit",
      rds_validation_reason = NA_character_,
      prior_tail_warning_flag = NA,
      prior_tail_warning_detail = NA_character_,
      prior_tail_sensitive_any = NA,
      rds_path = if (save_model_rds) existing_rds_path else NA_character_
    )

    emit_stage8_progress(ii, total_models, model_id_now, paste0("sampling error: ", fit_error_message))
    next
  }

  trace_records[[length(trace_records) + 1L]] <- make_trace_record(
    fit = fit,
    model_id = model_id_now,
    family_code = model_row$family_code[[1L]],
    K_inc = ncol(design_bundle$X_inc),
    K_lat = ncol(design_bundle$X_lat)
  )

  parameter_names <- c("alpha_inc", "gamma0")
  parameter_names <- c(parameter_names, paste0("beta_inc[", seq_len(ncol(design_bundle$X_inc)), "]"))
  parameter_names <- c(parameter_names, paste0("gamma_lat[", seq_len(ncol(design_bundle$X_lat)), "]"))
  if (model_row$family_code[[1L]] == "W") parameter_names <- c(parameter_names, "rho_W")
  if (model_row$family_code[[1L]] == "LN") parameter_names <- c(parameter_names, "log_sigma_LN")
  if (model_row$family_code[[1L]] == "LL") parameter_names <- c(parameter_names, "psi_LL")

  param_array <- posterior::as_draws_array(as.array(fit, pars = parameter_names))
  param_diag_tbl <- posterior::summarise_draws(
    param_array,
    mean = base::mean,
    sd = stats::sd,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )

  param_draws_mat <- posterior::as_draws_matrix(param_array)

  coefficient_rows[[length(coefficient_rows) + 1L]] <- tibble(
    dataset_key = dataset_name,
    model_id = model_id_now,
    retained_fit_id = model_id_now,
    branch = stage8_branch_label,
    risk_scale = stage8_risk_scale,
    prior_branch = model_row$prior_branch[[1L]],
    site_prior_family = model_row$site_prior_family_export[[1L]],
    structural_model_id = model_row$structural_model_id[[1L]],
    formula_anchor = model_row$formula_anchor[[1L]],
    family_code = model_row$family_code[[1L]],
    latency_family = model_row$latency_family[[1L]],
    parameter = colnames(param_draws_mat),
    mean = apply(param_draws_mat, 2L, mean),
    sd = apply(param_draws_mat, 2L, stats::sd),
    q025 = apply(param_draws_mat, 2L, stats::quantile, probs = 0.025, names = FALSE),
    q50 = apply(param_draws_mat, 2L, stats::quantile, probs = 0.500, names = FALSE),
    q975 = apply(param_draws_mat, 2L, stats::quantile, probs = 0.975, names = FALSE)
  )

  diagnostic_parameter_rows[[length(diagnostic_parameter_rows) + 1L]] <- tibble(
    dataset_key = dataset_name,
    model_id = model_id_now,
    retained_fit_id = model_id_now,
    branch = stage8_branch_label,
    risk_scale = stage8_risk_scale,
    prior_branch = model_row$prior_branch[[1L]],
    site_prior_family = model_row$site_prior_family_export[[1L]],
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

  draws_compact <- extract_draws_compact(
    fit = fit,
    K_inc = ncol(design_bundle$X_inc),
    K_lat = ncol(design_bundle$X_lat)
  )

  total_draws <- length(draws_compact$alpha_inc)
  set.seed(stan_seed + 100000L + ii)
  keep_draw_idx <- if (total_draws <= posterior_prediction_draws) {
    seq_len(total_draws)
  } else {
    sort(sample(seq_len(total_draws), size = posterior_prediction_draws, replace = FALSE))
  }

  draws_pred <- list(
    alpha_inc = draws_compact$alpha_inc[keep_draw_idx],
    beta_inc = draws_compact$beta_inc[keep_draw_idx, , drop = FALSE],
    gamma0 = draws_compact$gamma0[keep_draw_idx],
    gamma_lat = draws_compact$gamma_lat[keep_draw_idx, , drop = FALSE],
    rho_W = draws_compact$rho_W[keep_draw_idx],
    log_sigma_LN = draws_compact$log_sigma_LN[keep_draw_idx],
    psi_LL = draws_compact$psi_LL[keep_draw_idx]
  )

  linear_terms <- compute_linear_terms(
    draws = draws_pred,
    X_inc = design_bundle$X_inc,
    X_lat = design_bundle$X_lat
  )

  supported_horizons <- if (dataset_name == "PNU") c(1L, 2L) else c(1L, 2L, 5L)
  supported_risk_list <- lapply(supported_horizons, function(h) {
    fh_sup <- family_survival_hazard(h, model_row$family_code[[1L]], linear_terms$mu_lat_mat, linear_terms$median_mat, draws_pred)
    1 - ((1 - linear_terms$pi_mat) + linear_terms$pi_mat * fh_sup$Su)
  })

  posterior_degeneracy <- compute_degeneracy(
    pi_mat = linear_terms$pi_mat,
    median_mat = linear_terms$median_mat,
    supported_risk_list = supported_risk_list
  )

  info_criteria <- compute_information_criteria(draws_compact$log_lik)

  cure_prob_summary <- summarize_cols_matrix(linear_terms$cure_prob_mat)
  susceptible_prob_summary <- summarize_cols_matrix(linear_terms$pi_mat)

  mean_cure_draw <- rowMeans(linear_terms$cure_prob_mat)
  mean_susceptible_draw <- rowMeans(linear_terms$pi_mat)
  mstu_draws <- compute_mstu_draws(
    family_code = model_row$family_code[[1L]],
    mu_lat_mat = linear_terms$mu_lat_mat,
    median_mat = linear_terms$median_mat,
    draws = draws_pred
  )

  followup_not_contradicted_flag_now <- screening_value_or_na(screening_row, "followup_not_contradicted_flag")
  followup_ok <- stage8_parse_logicalish(followup_not_contradicted_flag_now)[[1L]]
  uncured_mean_support_flag_now <- isTRUE(followup_ok) && isTRUE(max(dataset_df$time_year, na.rm = TRUE) >= 5)

  prediction_model_rows <- list()
  cohort_model_rows <- list()
  class_model_rows <- list()
  ppc_model_rows <- list()
  uncured_model_rows <- list()

  for (h in horizons_year) {
    fh <- family_survival_hazard(
      horizon_year = h,
      family_code = model_row$family_code[[1L]],
      mu_lat_mat = linear_terms$mu_lat_mat,
      median_mat = linear_terms$median_mat,
      draws = draws_pred
    )

    pop_surv <- (1 - linear_terms$pi_mat) + linear_terms$pi_mat * fh$Su
    risk_mat <- 1 - pop_surv

    subj_surv_summary <- summarize_cols_matrix(pop_surv)
    subj_risk_summary <- summarize_cols_matrix(risk_mat)
    subj_su_summary <- summarize_cols_matrix(fh$Su)

    prediction_model_rows[[length(prediction_model_rows) + 1L]] <- bind_cols(
      tibble(
        dataset_key = dataset_name,
        model_id = model_id_now,
        retained_fit_id = model_id_now,
        branch = stage8_branch_label,
        risk_scale = stage8_risk_scale,
        prior_branch = model_row$prior_branch[[1L]],
        site_prior_family = model_row$site_prior_family_export[[1L]],
        structural_model_id = model_row$structural_model_id[[1L]],
        formula_anchor = model_row$formula_anchor[[1L]],
        family_code = model_row$family_code[[1L]],
        latency_family = model_row$latency_family[[1L]],
        site_placement = model_row$site_placement[[1L]],
        horizon = h
      ),
      horizon_metadata %>%
        filter(dataset_key == dataset_name, horizon == h) %>%
        select(support_tier, horizon_evidence_class, claim_restriction_flag, interpretation_note),
      design_bundle$id_df,
      tibble(
        cure_prob_mean = cure_prob_summary$mean,
        cure_prob_q025 = cure_prob_summary$q025,
        cure_prob_q50 = cure_prob_summary$q50,
        cure_prob_q975 = cure_prob_summary$q975,
        susceptible_prob_mean = susceptible_prob_summary$mean,
        susceptible_prob_q025 = susceptible_prob_summary$q025,
        susceptible_prob_q50 = susceptible_prob_summary$q50,
        susceptible_prob_q975 = susceptible_prob_summary$q975,
        S_pop_mean = subj_surv_summary$mean,
        S_pop_q025 = subj_surv_summary$q025,
        S_pop_q50 = subj_surv_summary$q50,
        S_pop_q975 = subj_surv_summary$q975,
        risk_mean = subj_risk_summary$mean,
        risk_q025 = subj_risk_summary$q025,
        risk_q50 = subj_risk_summary$q50,
        risk_q975 = subj_risk_summary$q975,
        uncured_survival_mean = subj_su_summary$mean,
        uncured_survival_q025 = subj_su_summary$q025,
        uncured_survival_q50 = subj_su_summary$q50,
        uncured_survival_q975 = subj_su_summary$q975
      )
    )

    mean_risk_draw <- rowMeans(risk_mat)
    mean_surv_draw <- rowMeans(pop_surv)
    mean_uncured_surv_draw <- rowMeans(fh$Su)
    mean_hazard_draw <- rowMeans(fh$haz)

    cohort_model_rows[[length(cohort_model_rows) + 1L]] <- tibble(
      dataset_key = dataset_name,
      model_id = model_id_now,
      retained_fit_id = model_id_now,
      branch = stage8_branch_label,
      risk_scale = stage8_risk_scale,
      prior_branch = model_row$prior_branch[[1L]],
      site_prior_family = model_row$site_prior_family_export[[1L]],
      structural_model_id = model_row$structural_model_id[[1L]],
      formula_anchor = model_row$formula_anchor[[1L]],
      family_code = model_row$family_code[[1L]],
      latency_family = model_row$latency_family[[1L]],
      site_placement = model_row$site_placement[[1L]],
      horizon = h,
      risk_mean = mean(mean_risk_draw),
      risk_q025 = stats::quantile(mean_risk_draw, 0.025, names = FALSE),
      risk_q50 = stats::quantile(mean_risk_draw, 0.500, names = FALSE),
      risk_q975 = stats::quantile(mean_risk_draw, 0.975, names = FALSE),
      survival_mean = mean(mean_surv_draw),
      survival_q025 = stats::quantile(mean_surv_draw, 0.025, names = FALSE),
      survival_q50 = stats::quantile(mean_surv_draw, 0.500, names = FALSE),
      survival_q975 = stats::quantile(mean_surv_draw, 0.975, names = FALSE),
      uncured_survival_mean = mean(mean_uncured_surv_draw),
      uncured_survival_q025 = stats::quantile(mean_uncured_surv_draw, 0.025, names = FALSE),
      uncured_survival_q50 = stats::quantile(mean_uncured_surv_draw, 0.500, names = FALSE),
      uncured_survival_q975 = stats::quantile(mean_uncured_surv_draw, 0.975, names = FALSE),
      uncured_risk_mean = mean(1 - mean_uncured_surv_draw),
      uncured_risk_q025 = stats::quantile(1 - mean_uncured_surv_draw, 0.025, names = FALSE),
      uncured_risk_q50 = stats::quantile(1 - mean_uncured_surv_draw, 0.500, names = FALSE),
      uncured_risk_q975 = stats::quantile(1 - mean_uncured_surv_draw, 0.975, names = FALSE),
      hazard_mean = mean(mean_hazard_draw),
      hazard_q025 = stats::quantile(mean_hazard_draw, 0.025, names = FALSE),
      hazard_q50 = stats::quantile(mean_hazard_draw, 0.500, names = FALSE),
      hazard_q975 = stats::quantile(mean_hazard_draw, 0.975, names = FALSE),
      cure_fraction_mean = mean(mean_cure_draw),
      cure_fraction_q025 = stats::quantile(mean_cure_draw, 0.025, names = FALSE),
      cure_fraction_q50 = stats::quantile(mean_cure_draw, 0.500, names = FALSE),
      cure_fraction_q975 = stats::quantile(mean_cure_draw, 0.975, names = FALSE),
      susceptible_fraction_mean = mean(mean_susceptible_draw),
      susceptible_fraction_q025 = stats::quantile(mean_susceptible_draw, 0.025, names = FALSE),
      susceptible_fraction_q50 = stats::quantile(mean_susceptible_draw, 0.500, names = FALSE),
      susceptible_fraction_q975 = stats::quantile(mean_susceptible_draw, 0.975, names = FALSE),
      discrimination_AUC = NA_real_,
      calibration_intercept = NA_real_,
      calibration_slope = NA_real_,
      brier_score = NA_real_,
      IBS = NA_real_,
      MSTu_mean = if (uncured_mean_support_flag_now) mean(mstu_draws[is.finite(mstu_draws)], na.rm = TRUE) else NA_real_,
      MSTu_q025 = if (uncured_mean_support_flag_now) stats::quantile(mstu_draws[is.finite(mstu_draws)], 0.025, names = FALSE, na.rm = TRUE) else NA_real_,
      MSTu_q50 = if (uncured_mean_support_flag_now) stats::quantile(mstu_draws[is.finite(mstu_draws)], 0.500, names = FALSE, na.rm = TRUE) else NA_real_,
      MSTu_q975 = if (uncured_mean_support_flag_now) stats::quantile(mstu_draws[is.finite(mstu_draws)], 0.975, names = FALSE, na.rm = TRUE) else NA_real_,
      uncured_mean_support_flag = uncured_mean_support_flag_now
    ) %>%
      left_join(
        horizon_metadata %>% filter(dataset_key == dataset_name, horizon == h),
        by = c("dataset_key", "horizon")
      )

    horizon_ref <- ipcw_registry[[dataset_name]] %>% filter(horizon == h)

    ppc_model_rows[[length(ppc_model_rows) + 1L]] <- tibble(
      dataset_key = dataset_name,
      model_id = model_id_now,
      retained_fit_id = model_id_now,
      branch = stage8_branch_label,
      risk_scale = stage8_risk_scale,
      prior_branch = model_row$prior_branch[[1L]],
      site_prior_family = model_row$site_prior_family_export[[1L]],
      structural_model_id = model_row$structural_model_id[[1L]],
      family_code = model_row$family_code[[1L]],
      horizon = h,
      observed_km_risk = horizon_ref$observed_km_risk[[1L]],
      posterior_mean_risk = mean(mean_risk_draw),
      posterior_q025_risk = stats::quantile(mean_risk_draw, 0.025, names = FALSE),
      posterior_q975_risk = stats::quantile(mean_risk_draw, 0.975, names = FALSE),
      absolute_difference = abs(mean(mean_risk_draw) - horizon_ref$observed_km_risk[[1L]]),
      gross_contradiction_flag = (
        h %in% supported_horizons &&
          (
            horizon_ref$observed_km_risk[[1L]] < stats::quantile(mean_risk_draw, 0.025, names = FALSE) ||
              horizon_ref$observed_km_risk[[1L]] > stats::quantile(mean_risk_draw, 0.975, names = FALSE)
          ) &&
          abs(mean(mean_risk_draw) - horizon_ref$observed_km_risk[[1L]]) > ppc_tolerance_abs
      )
    )

    class_out <- compute_classification_summary(
      risk_draws = risk_mat,
      horizon_row = horizon_ref,
      thresholds = risk_thresholds,
      cohort_n = nrow(dataset_df)
    )

    class_model_rows[[length(class_model_rows) + 1L]] <- class_out %>%
      mutate(
        dataset_key = dataset_name,
        model_id = model_id_now,
        retained_fit_id = model_id_now,
        branch = stage8_branch_label,
        risk_scale = stage8_risk_scale,
        prior_branch = model_row$prior_branch[[1L]],
        site_prior_family = model_row$site_prior_family_export[[1L]],
        structural_model_id = model_row$structural_model_id[[1L]],
        formula_anchor = model_row$formula_anchor[[1L]],
        family_code = model_row$family_code[[1L]],
        latency_family = model_row$latency_family[[1L]],
        site_placement = model_row$site_placement[[1L]],
        horizon = h,
        observed_km_risk = horizon_ref$observed_km_risk[[1L]]
      ) %>%
      left_join(
        horizon_metadata %>% filter(dataset_key == dataset_name, horizon == h),
        by = c("dataset_key", "horizon")
      )

    uncured_model_rows[[length(uncured_model_rows) + 1L]] <- tibble(
      dataset_key = dataset_name,
      model_id = model_id_now,
      retained_fit_id = model_id_now,
      branch = stage8_branch_label,
      risk_scale = stage8_risk_scale,
      prior_branch = model_row$prior_branch[[1L]],
      site_prior_family = model_row$site_prior_family_export[[1L]],
      structural_model_id = model_row$structural_model_id[[1L]],
      formula_anchor = model_row$formula_anchor[[1L]],
      family_code = model_row$family_code[[1L]],
      latency_family = model_row$latency_family[[1L]],
      site_placement = model_row$site_placement[[1L]],
      horizon = h,
      cure_fraction = mean(mean_cure_draw),
      susceptible_fraction = mean(mean_susceptible_draw),
      uncured_survival = mean(mean_uncured_surv_draw),
      uncured_survival_q025 = stats::quantile(mean_uncured_surv_draw, 0.025, names = FALSE),
      uncured_survival_q50 = stats::quantile(mean_uncured_surv_draw, 0.500, names = FALSE),
      uncured_survival_q975 = stats::quantile(mean_uncured_surv_draw, 0.975, names = FALSE),
      uncured_risk = mean(1 - mean_uncured_surv_draw),
      uncured_risk_q025 = stats::quantile(1 - mean_uncured_surv_draw, 0.025, names = FALSE),
      uncured_risk_q50 = stats::quantile(1 - mean_uncured_surv_draw, 0.500, names = FALSE),
      uncured_risk_q975 = stats::quantile(1 - mean_uncured_surv_draw, 0.975, names = FALSE),
      MSTu_mean = if (uncured_mean_support_flag_now) mean(mstu_draws[is.finite(mstu_draws)], na.rm = TRUE) else NA_real_,
      MSTu_q025 = if (uncured_mean_support_flag_now) stats::quantile(mstu_draws[is.finite(mstu_draws)], 0.025, names = FALSE, na.rm = TRUE) else NA_real_,
      MSTu_q50 = if (uncured_mean_support_flag_now) stats::quantile(mstu_draws[is.finite(mstu_draws)], 0.500, names = FALSE, na.rm = TRUE) else NA_real_,
      MSTu_q975 = if (uncured_mean_support_flag_now) stats::quantile(mstu_draws[is.finite(mstu_draws)], 0.975, names = FALSE, na.rm = TRUE) else NA_real_,
      uncured_mean_support_flag = uncured_mean_support_flag_now
    )
  }

  ppc_model_tbl <- bind_rows_safe(ppc_model_rows)
  cohort_model_tbl <- bind_rows_safe(cohort_model_rows)
  class_model_tbl <- bind_rows_safe(class_model_rows)
  prediction_model_tbl <- bind_rows_safe(prediction_model_rows)
  uncured_model_tbl <- bind_rows_safe(uncured_model_rows)

  ppc_gross_contradiction_flag <- any(ppc_model_tbl$gross_contradiction_flag %in% TRUE, na.rm = TRUE)
  max_rhat <- max(param_diag_tbl$rhat, na.rm = TRUE)
  min_bulk_ess <- min(param_diag_tbl$ess_bulk, na.rm = TRUE)
  min_tail_ess <- min(param_diag_tbl$ess_tail, na.rm = TRUE)

  admissibility_reasons <- character()
  prior_deg_flag <- any(prior_predictive_rows[[length(prior_predictive_rows)]]$prior_degenerate_flag %in% TRUE, na.rm = TRUE)
  if (isTRUE(prior_deg_flag)) admissibility_reasons <- c(admissibility_reasons, "prior_degenerate")
  if (posterior_degeneracy$degenerate_flag[[1L]]) admissibility_reasons <- c(admissibility_reasons, "posterior_degenerate")
  if (!is.finite(max_rhat) || max_rhat >= rhat_max_threshold) admissibility_reasons <- c(admissibility_reasons, "rhat")
  if (!is.finite(min_bulk_ess) || min_bulk_ess < ess_min_threshold) admissibility_reasons <- c(admissibility_reasons, "bulk_ess")
  if (!is.finite(min_tail_ess) || min_tail_ess < ess_min_threshold) admissibility_reasons <- c(admissibility_reasons, "tail_ess")
  if (divergences > 0L) admissibility_reasons <- c(admissibility_reasons, "divergences")
  if (treedepth_exceeded > 0L) admissibility_reasons <- c(admissibility_reasons, "treedepth")
  if (isTRUE(ppc_gross_contradiction_flag)) admissibility_reasons <- c(admissibility_reasons, "ppc")

  admissible_flag <- length(admissibility_reasons) == 0L

  hazard_values <- cohort_model_tbl$hazard_mean[match(horizons_year, cohort_model_tbl$horizon)]
  hazard_rows[[length(hazard_rows) + 1L]] <- if (isTRUE(admissible_flag)) {
    tibble(
      dataset_key = dataset_name,
      model_id = model_id_now,
      retained_fit_id = model_id_now,
      branch = stage8_branch_label,
      risk_scale = stage8_risk_scale,
      prior_branch = model_row$prior_branch[[1L]],
      site_prior_family = model_row$site_prior_family_export[[1L]],
      structural_model_id = model_row$structural_model_id[[1L]],
      formula_anchor = model_row$formula_anchor[[1L]],
      family_code = model_row$family_code[[1L]],
      latency_family = model_row$latency_family[[1L]],
      site_placement = model_row$site_placement[[1L]],
      hazard_target = "cohort_mean_population_hazard",
      hazard_1y = hazard_values[[1L]],
      hazard_2y = hazard_values[[2L]],
      hazard_3y = hazard_values[[3L]],
      hazard_4y = hazard_values[[4L]],
      hazard_5y = hazard_values[[5L]],
      hazard_6y = hazard_values[[6L]],
      hazard_7y = hazard_values[[7L]],
      hazard_8y = hazard_values[[8L]],
      hazard_9y = hazard_values[[9L]],
      hazard_10y = hazard_values[[10L]],
      hazard_ratio_10y_vs_1y = hazard_values[[10L]] / hazard_values[[1L]],
      shape_class = classify_hazard_shape(hazard_values)
    )
  } else {
    NULL
  }

  model_registry_rows[[length(model_registry_rows) + 1L]] <- tibble(
    dataset_key = dataset_name,
    model_id = model_id_now,
    retained_fit_id = model_id_now,
    branch = stage8_branch_label,
    risk_scale = stage8_risk_scale,
    structural_model_id = model_row$structural_model_id[[1L]],
    formula_anchor = model_row$formula_anchor[[1L]],
    family_code = model_row$family_code[[1L]],
    latency_family = model_row$latency_family[[1L]],
    prior_branch = model_row$prior_branch[[1L]],
    site_prior_family = model_row$site_prior_family_export[[1L]],
    site_placement = model_row$site_placement[[1L]],
    fit_status = fit_status,
    fit_error_message = fit_error_message,
    admissible_flag = admissible_flag,
    admissibility_reasons = if (length(admissibility_reasons) == 0L) "" else paste(admissibility_reasons, collapse = "|"),
    divergences = divergences,
    treedepth_exceeded = treedepth_exceeded,
    max_rhat = max_rhat,
    min_bulk_ess = min_bulk_ess,
    min_tail_ess = min_tail_ess,
    waic = info_criteria$waic,
    p_waic = info_criteria$p_waic,
    looic = info_criteria$looic,
    p_loo = info_criteria$p_loo,
    p_waic_high_n = info_criteria$p_waic_high_n,
    p_waic_high_pct = info_criteria$p_waic_high_pct,
    pareto_k_max = info_criteria$pareto_k_max,
    pareto_k_bad_n = info_criteria$pareto_k_bad_n,
    pareto_k_bad_pct = info_criteria$pareto_k_bad_pct,
    pareto_k_very_bad_n = info_criteria$pareto_k_very_bad_n,
    waic_warning_flag = info_criteria$waic_warning_flag,
    loo_warning_flag = info_criteria$loo_warning_flag,
    info_criteria_warning_detail = info_criteria$info_criteria_warning_detail,
    prior_degenerate_flag = prior_deg_flag,
    posterior_degenerate_flag = posterior_degeneracy$degenerate_flag[[1L]],
    ppc_gross_contradiction_flag = ppc_gross_contradiction_flag,
    n = nrow(dataset_df),
    n_transition = sum(dataset_df$event_main),
    n_remission = sum(dataset_df$remission_flag),
    n_main_censoring = sum(dataset_df$censor_main),
    cure_model_eligibility_flag = screening_value_or_na(screening_row, "cure_model_eligibility_flag"),
    primary_gate_method = screening_value_or_na(screening_row, "primary_gate_method"),
    primary_gate_flag = screening_value_or_na(screening_row, "primary_gate_flag"),
    receus_primary_class = screening_value_or_na(screening_row, "receus_primary_class"),
    presence_modifier_flag = screening_value_or_na(screening_row, "presence_modifier_flag"),
    cure_presence_support_flag = screening_value_or_na(screening_row, "cure_presence_support_flag"),
    followup_contradiction_flag = screening_value_or_na(screening_row, "followup_contradiction_flag"),
    followup_not_contradicted_flag = screening_value_or_na(screening_row, "followup_not_contradicted_flag"),
    screening_note = screening_value_or_na(screening_row, "screening_note"),
    carry_forward_stage8 = screening_logical_or_na(screening_row, "carry_forward_stage8"),
    reused_existing_rds_flag = FALSE,
    rds_validation_status = "new_fit",
    rds_validation_reason = NA_character_,
    prior_tail_warning_flag = NA,
    prior_tail_warning_detail = NA_character_,
    prior_tail_sensitive_any = NA,
    rds_path = if (save_model_rds) existing_rds_path else NA_character_
  )

  ppc_rows[[length(ppc_rows) + 1L]] <- ppc_model_tbl

  if (isTRUE(admissible_flag)) {
    prediction_rows[[length(prediction_rows) + 1L]] <- prediction_model_tbl
    cohort_rows[[length(cohort_rows) + 1L]] <- cohort_model_tbl
    classification_rows[[length(classification_rows) + 1L]] <- class_model_tbl
    uncured_rows[[length(uncured_rows) + 1L]] <- uncured_model_tbl
    fit_draw_map[[model_id_now]] <- list(
      alpha_inc = draws_pred$alpha_inc,
      beta_inc = draws_pred$beta_inc
    )
  }

  if (isTRUE(save_model_rds)) {
    bundle <- list(
      bundle_version = "stage8a_v1",
      model_registry_row = model_registry_rows[[length(model_registry_rows)]],
      coefficient_summary = coefficient_rows[[length(coefficient_rows)]],
      diagnostics_parameter_level = diagnostic_parameter_rows[[length(diagnostic_parameter_rows)]],
      ppc_summary = ppc_model_tbl,
      prior_predictive_summary = prior_predictive_rows[[length(prior_predictive_rows)]],
      posterior_incidence_draws = if (isTRUE(admissible_flag)) list(
        alpha_inc = draws_pred$alpha_inc,
        beta_inc = draws_pred$beta_inc
      ) else NULL,
      prediction_long = if (isTRUE(admissible_flag)) prediction_model_tbl else tibble(),
      cohort_yearly = if (isTRUE(admissible_flag)) cohort_model_tbl else tibble(),
      classification = if (isTRUE(admissible_flag)) class_model_tbl else tibble(),
      uncured_decomposition = if (isTRUE(admissible_flag)) uncured_model_tbl else tibble(),
      hazard_plausibility = if (isTRUE(admissible_flag)) hazard_rows[[length(hazard_rows)]] else tibble()
    )
    saveRDS(bundle, existing_rds_path)
  }

  rm(fit, draws_compact, draws_pred, linear_terms, param_array, param_draws_mat)
  gc(verbose = FALSE)

  emit_stage8_progress(
    ii,
    total_models,
    model_id_now,
    paste0(
      "completed; admissible=", admissible_flag,
      "; WAIC=", formatC(info_criteria$waic, digits = 2L, format = "f"),
      "; LOOIC=", formatC(info_criteria$looic, digits = 2L, format = "f")
    )
  )
}


# 🔴 Finalize: assemble retained tables, prior-tail flags, and reporting layers ===============================
model_registry <- bind_rows_safe(model_registry_rows) %>%
  arrange(factor(model_id, levels = model_order)) %>%
  mutate(
    reused_existing_rds_flag = dplyr::coalesce(as.logical(reused_existing_rds_flag), FALSE),
    rds_validation_status = dplyr::coalesce(as.character(rds_validation_status), ifelse(reused_existing_rds_flag, "existing_rds_skipped", "new_fit")),
    rds_validation_reason = as.character(rds_validation_reason)
  )

coefficient_summary <- bind_rows_safe(coefficient_rows) %>%
  arrange(factor(model_id, levels = model_order), parameter)

diagnostics_parameter_level <- bind_rows_safe(diagnostic_parameter_rows) %>%
  arrange(factor(model_id, levels = model_order), parameter)

ppc_summary <- bind_rows_safe(ppc_rows) %>%
  left_join(
    horizon_metadata %>%
      select(dataset_key, horizon, support_tier, horizon_evidence_class, claim_restriction_flag),
    by = c("dataset_key", "horizon")
  ) %>%
  arrange(factor(model_id, levels = model_order), horizon)

prior_predictive_summary <- bind_rows_safe(prior_predictive_rows) %>%
  annotate_prior_predictive_summary() %>%
  arrange(factor(model_id, levels = model_order), metric, horizon)

prediction_long <- bind_rows_safe(prediction_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon, unique_person_id)

cohort_yearly <- bind_rows_safe(cohort_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon)

classification_table <- bind_rows_safe(classification_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon, threshold)

uncured_only_decomposition_panel <- bind_rows_safe(uncured_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon)

hazard_plausibility <- bind_rows_safe(hazard_rows) %>%
  arrange(factor(model_id, levels = model_order))

prior_warning_by_model <- prior_predictive_summary %>%
  group_by(dataset_key, model_id) %>%
  summarise(
    prior_tail_warning_flag = any(prior_tail_warning_flag %in% TRUE, na.rm = TRUE),
    prior_tail_warning_detail = paste(unique(na.omit(prior_tail_warning_detail[prior_tail_warning_flag %in% TRUE])), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(prior_tail_warning_detail = dplyr::na_if(prior_tail_warning_detail, ""))

model_registry <- model_registry %>%
  left_join(prior_warning_by_model, by = c("dataset_key", "model_id")) %>%
  mutate(
    prior_tail_warning_flag = dplyr::coalesce(prior_tail_warning_flag.y, prior_tail_warning_flag.x, FALSE),
    prior_tail_warning_detail = dplyr::coalesce(prior_tail_warning_detail.y, prior_tail_warning_detail.x)
  ) %>%
  select(-any_of(c("prior_tail_warning_flag.x", "prior_tail_warning_flag.y", "prior_tail_warning_detail.x", "prior_tail_warning_detail.y")))

anchor_vs_neutral_delta_panel <- build_anchor_vs_neutral_delta_panel(
  cohort_df = cohort_yearly,
  class_df = classification_table
)

prior_tail_by_horizon <- compute_prior_tail_flags(anchor_vs_neutral_delta_panel)

if (nrow(cohort_yearly) > 0L) {
  cohort_yearly <- cohort_yearly %>%
    left_join(
      prior_tail_by_horizon,
      by = c("dataset_key", "structural_model_id", "family_code", "site_prior_family", "horizon")
    ) %>%
    mutate(
      prior_tail_sensitive = dplyr::coalesce(prior_tail_sensitive, FALSE),
      claim_restriction_flag = if_else(
        prior_tail_sensitive & horizon_evidence_class == "mostly_extrapolated",
        "projection_plus_prior_sensitive",
        claim_restriction_flag
      ),
      admissibility_flag = TRUE
    )
}

if (nrow(classification_table) > 0L) {
  classification_table <- classification_table %>%
    left_join(
      prior_tail_by_horizon,
      by = c("dataset_key", "structural_model_id", "family_code", "site_prior_family", "horizon")
    ) %>%
    mutate(
      prior_tail_sensitive = dplyr::coalesce(prior_tail_sensitive, FALSE),
      claim_restriction_flag = if_else(
        prior_tail_sensitive & horizon_evidence_class == "mostly_extrapolated",
        "projection_plus_prior_sensitive",
        claim_restriction_flag
      ),
      admissibility_flag = TRUE
    )
}

if (nrow(prediction_long) > 0L) {
  prediction_long <- prediction_long %>%
    left_join(
      prior_tail_by_horizon,
      by = c("dataset_key", "structural_model_id", "family_code", "site_prior_family", "horizon")
    ) %>%
    mutate(
      prior_tail_sensitive = dplyr::coalesce(prior_tail_sensitive, FALSE),
      claim_restriction_flag = if_else(
        prior_tail_sensitive & horizon_evidence_class == "mostly_extrapolated",
        "projection_plus_prior_sensitive",
        claim_restriction_flag
      ),
      admissibility_flag = TRUE
    )
}

if (nrow(uncured_only_decomposition_panel) > 0L && nrow(cohort_yearly) > 0L) {
  uncured_only_decomposition_panel <- uncured_only_decomposition_panel %>%
    left_join(
      cohort_yearly %>%
        select(
          model_id, horizon, support_tier, horizon_evidence_class,
          claim_restriction_flag, prior_tail_sensitive
        ),
      by = c("model_id", "horizon")
    )
}

prior_tail_by_model <- prior_tail_by_horizon %>%
  group_by(dataset_key, structural_model_id, family_code, site_prior_family) %>%
  summarise(prior_tail_sensitive_any = any(prior_tail_sensitive %in% TRUE), .groups = "drop")

if (nrow(prior_tail_by_model) > 0L) {
  model_registry <- model_registry %>%
    left_join(
      prior_tail_by_model,
      by = c("dataset_key", "structural_model_id", "family_code", "site_prior_family")
    ) %>%
    mutate(
      prior_tail_sensitive_any = dplyr::coalesce(
        prior_tail_sensitive_any.y,
        prior_tail_sensitive_any.x,
        FALSE
      )
    ) %>%
    select(-any_of(c("prior_tail_sensitive_any.x", "prior_tail_sensitive_any.y")))
} else {
  model_registry <- model_registry %>%
    mutate(prior_tail_sensitive_any = dplyr::coalesce(prior_tail_sensitive_any, FALSE))
}

if (nrow(anchor_vs_neutral_delta_panel) > 0L) {
  anchor_vs_neutral_delta_panel <- anchor_vs_neutral_delta_panel %>%
    left_join(
      prior_tail_by_horizon,
      by = c("dataset_key", "structural_model_id", "family_code", "site_prior_family", "horizon")
    ) %>%
    mutate(
      prior_tail_sensitive = dplyr::coalesce(prior_tail_sensitive, FALSE),
      claim_restriction_flag = if_else(
        prior_tail_sensitive & horizon_evidence_class == "mostly_extrapolated",
        "projection_plus_prior_sensitive",
        claim_restriction_flag
      )
    )
} else {
  anchor_vs_neutral_delta_panel <- empty_anchor_vs_neutral_delta_panel()
}

incidence_anchor_update_panel <- build_incidence_anchor_update_panel(
  model_grid = model_grid,
  retained_model_ids = cohort_yearly$model_id,
  fit_draw_map = fit_draw_map,
  prior_artifacts = prior_artifacts
)

if (nrow(incidence_anchor_update_panel) == 0L) {
  incidence_anchor_update_panel <- empty_incidence_anchor_update_panel()
}

horizon_support_panel <- if (nrow(cohort_yearly) > 0L) {
  horizon_metadata %>%
    left_join(
      cohort_yearly %>%
        group_by(dataset_key, horizon) %>%
        summarise(prior_tail_sensitive = any(prior_tail_sensitive %in% TRUE), .groups = "drop"),
      by = c("dataset_key", "horizon")
    ) %>%
    mutate(
      # Horizon support is dataset-level, so retain a single flag when any retained fit is prior-tail-sensitive.
      prior_tail_sensitive = dplyr::coalesce(prior_tail_sensitive, FALSE),
      prior_tail_sensitive_any_fit = prior_tail_sensitive
    ) %>%
    arrange(dataset_key, horizon)
} else {
  horizon_metadata %>%
    mutate(
      prior_tail_sensitive = FALSE,
      prior_tail_sensitive_any_fit = prior_tail_sensitive
    ) %>%
    arrange(dataset_key, horizon)
}

cohort_part_long <- if (nrow(cohort_yearly) > 0L) {
  cohort_yearly %>%
    mutate(
      threshold = NA_real_,
      positive_rate_mean = NA_real_,
      positive_rate_q025 = NA_real_,
      positive_rate_q50 = NA_real_,
      positive_rate_q975 = NA_real_,
      FPR_mean = NA_real_,
      FPR_q025 = NA_real_,
      FPR_q50 = NA_real_,
      FPR_q975 = NA_real_,
      false_positive_burden_mean = NA_real_,
      false_positive_burden_q025 = NA_real_,
      false_positive_burden_q50 = NA_real_,
      false_positive_burden_q975 = NA_real_,
      false_positive_count_mean = NA_real_,
      false_positive_count_q025 = NA_real_,
      false_positive_count_q50 = NA_real_,
      false_positive_count_q975 = NA_real_,
      FP100_mean = NA_real_,
      FP100_q025 = NA_real_,
      FP100_q50 = NA_real_,
      FP100_q975 = NA_real_,
      PPV_mean = NA_real_,
      PPV_q025 = NA_real_,
      PPV_q50 = NA_real_,
      PPV_q975 = NA_real_,
      TPR_mean = NA_real_,
      TPR_q025 = NA_real_,
      TPR_q50 = NA_real_,
      TPR_q975 = NA_real_,
      NB_mean = NA_real_,
      NB_q025 = NA_real_,
      NB_q50 = NA_real_,
      NB_q975 = NA_real_,
      classification_estimable_flag = NA
    )
} else {
  tibble()
}

class_part_long <- if (nrow(classification_table) > 0L && nrow(cohort_yearly) > 0L) {
  classification_table %>%
    left_join(
      cohort_yearly %>%
        select(
          dataset_key, model_id, retained_fit_id, horizon,
          risk_mean, risk_q025, risk_q50, risk_q975,
          survival_mean, survival_q025, survival_q50, survival_q975,
          uncured_survival_mean, uncured_survival_q025, uncured_survival_q50, uncured_survival_q975,
          uncured_risk_mean, uncured_risk_q025, uncured_risk_q50, uncured_risk_q975,
          hazard_mean, hazard_q025, hazard_q50, hazard_q975,
          cure_fraction_mean, cure_fraction_q025, cure_fraction_q50, cure_fraction_q975,
          susceptible_fraction_mean, susceptible_fraction_q025, susceptible_fraction_q50, susceptible_fraction_q975,
          discrimination_AUC, calibration_intercept, calibration_slope, brier_score, IBS,
          MSTu_mean, MSTu_q025, MSTu_q50, MSTu_q975, uncured_mean_support_flag
        ),
      by = c("dataset_key", "model_id", "retained_fit_id", "horizon")
    )
} else {
  tibble()
}

performance_classification_long <- bind_rows_safe(list(cohort_part_long, class_part_long))
if (nrow(performance_classification_long) > 0L) {
  performance_classification_long <- performance_classification_long %>%
    mutate(
      instability_marker = claim_restriction_flag %in% c("projection_only", "projection_plus_prior_sensitive")
    ) %>%
    arrange(factor(model_id, levels = model_order), horizon, threshold)
} else {
  performance_classification_long <- tibble()
}

reporting_metadata <- if (nrow(performance_classification_long) > 0L) {
  performance_classification_long %>%
    transmute(
      dataset_key = dataset_key,
      model_id = model_id,
      retained_fit_id = retained_fit_id,
      branch = branch,
      risk_scale = risk_scale,
      prior_branch = prior_branch,
      site_prior_family = site_prior_family,
      horizon = horizon,
      threshold = threshold,
      support_tier = support_tier,
      horizon_evidence_class = horizon_evidence_class,
      claim_restriction_flag = claim_restriction_flag,
      prior_tail_sensitive = prior_tail_sensitive,
      admissibility_flag = admissibility_flag
    ) %>%
    distinct()
} else {
  tibble()
}

posterior_delta_vs_nocure <- if (nrow(performance_classification_long) > 0L) {
  build_delta_vs_nocure(
    performance_long = performance_classification_long,
    nocure_cohort_long = nocure_cohort_long,
    nocure_class_long = nocure_class_long
  )
} else {
  tibble()
}

diagnostic_pdf_path <- file.path(export_path, "bayes_stage8a_diagnostic_plots.pdf")
safe_generate_diagnostic_pdf(
  trace_records = trace_records,
  cohort_df = cohort_yearly,
  class_df = classification_table,
  ppc_df = ppc_summary,
  final_path = diagnostic_pdf_path
)

output_audit <- build_output_audit(
  model_grid = model_grid,
  model_registry = model_registry,
  prediction_long = prediction_long,
  performance_long = performance_classification_long,
  reporting_metadata = reporting_metadata,
  anchor_delta = anchor_vs_neutral_delta_panel,
  anchor_update = incidence_anchor_update_panel,
  uncured_decomp = uncured_only_decomposition_panel,
  hazard_plausibility = hazard_plausibility,
  diagnostic_pdf_path = diagnostic_pdf_path
)

# 🔴 Export: Stage-8A source-of-truth tables and figure-ready panels ===============================
write_csv_preserve_schema(simplify_scalar_list_cols(metadata_registry), file.path(export_path, "bayes_stage8a_metadata_registry.csv"))
write_csv_preserve_schema(simplify_scalar_list_cols(dataset_registry), file.path(export_path, "bayes_stage8a_dataset_registry.csv"))
write_csv_preserve_schema(simplify_scalar_list_cols(horizon_metadata), file.path(export_path, "bayes_stage8a_horizon_metadata.csv"))
write_csv_preserve_schema(simplify_scalar_list_cols(model_registry), file.path(export_path, "bayes_stage8a_model_registry.csv"))
write_csv_preserve_schema(simplify_scalar_list_cols(coefficient_summary), file.path(export_path, "bayes_stage8a_coefficient_summary.csv"))
write_csv_preserve_schema(simplify_scalar_list_cols(diagnostics_parameter_level), file.path(export_path, "bayes_stage8a_diagnostics_parameter_level.csv"))
write_csv_preserve_schema(simplify_scalar_list_cols(ppc_summary), file.path(export_path, "bayes_stage8a_ppc_summary.csv"))
write_csv_preserve_schema(simplify_scalar_list_cols(prior_predictive_summary), file.path(export_path, "bayes_stage8a_prior_predictive_summary.csv"))
write_csv_preserve_schema(simplify_scalar_list_cols(prediction_long), file.path(export_path, "bayes_stage8a_prediction_long.csv.gz"))
write_csv_preserve_schema(simplify_scalar_list_cols(performance_classification_long), file.path(export_path, "bayes_stage8a_performance_classification_long.csv"))
write_csv_preserve_schema(simplify_scalar_list_cols(reporting_metadata), file.path(export_path, "bayes_stage8a_reporting_metadata.csv"))
write_csv_preserve_schema(simplify_scalar_list_cols(anchor_vs_neutral_delta_panel), file.path(export_path, "bayes_stage8a_anchor_vs_neutral_delta_panel.csv"))
write_csv_preserve_schema(simplify_scalar_list_cols(incidence_anchor_update_panel), file.path(export_path, "bayes_stage8a_incidence_anchor_update_panel.csv"))
write_csv_preserve_schema(simplify_scalar_list_cols(uncured_only_decomposition_panel), file.path(export_path, "bayes_stage8a_uncured_only_decomposition_panel.csv"))
write_csv_preserve_schema(simplify_scalar_list_cols(hazard_plausibility), file.path(export_path, "bayes_stage8a_hazard_plausibility.csv"))
write_csv_preserve_schema(simplify_scalar_list_cols(posterior_delta_vs_nocure), file.path(export_path, "bayes_stage8a_posterior_delta_vs_nocure.csv"))
write_csv_preserve_schema(simplify_scalar_list_cols(horizon_support_panel), file.path(export_path, "bayes_stage8a_horizon_support_panel.csv"))
write_csv_preserve_schema(simplify_scalar_list_cols(output_audit), file.path(export_path, "bayes_stage8a_output_audit.csv"))

message(
  "Stage 8A completed. Total models attempted: ",
  total_models,
  "; admissible fits retained: ",
  sum(model_registry$admissible_flag %in% TRUE, na.rm = TRUE),
  "."
)
