# 🔴 Configure: Stage-8 paths, controls, and reporting thresholds ===============================
# 🟧 Configure: OS-aware project root ===============================
sys_name <- Sys.info()[["sysname"]]

project_root <- switch(
  sys_name,
  "Darwin"  = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  "Windows" = "C:/Users/clair/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  stop("Unsupported OS: ", sys_name)
)

# 🟧 Configure: Stage folder paths ===============================
merged_data_path   <- file.path(project_root, "MERGED_dataset3_pnu_snu.csv")
export_path        <- file.path(project_root, "stage8B_Bayesian remission-sensitive competing-risk extension")
stage5_export_path <- file.path(project_root, "stage5_Individualized no-cure comparator")
stage6_export_path <- file.path(project_root, "stage6_Cure-appropriateness screening")
stage8a_export_path <- file.path(project_root, "stage8A_Bayesian transition-only cure")

# 🟧 Helpers: find files inside stage folders ===============================
find_stage_file <- function(stage_path, pattern, recursive = TRUE) {
  if (!dir.exists(stage_path)) {
    stop("Directory does not exist: ", stage_path)
  }
  
  hits <- list.files(
    path = stage_path,
    pattern = pattern,
    full.names = TRUE,
    recursive = recursive,
    ignore.case = TRUE
  )
  
  if (length(hits) == 0) {
    stop("No matching file found in: ", stage_path, " | pattern: ", pattern)
  }
  
  if (length(hits) > 1) {
    warning("Multiple matching files found. Using the first one: ", hits[1])
  }
  
  hits[1]
}

# 🟧 Locate required CSV files from each stage folder ===============================
stage6_screening_flag_csv <- find_stage_file(
  stage6_export_path,
  "^stage6_carry_forward_flag_table\\.csv$"
)

stage5_nocure_cohort_csv <- find_stage_file(
  stage5_export_path,
  "^stage5_model_performance_long\\.csv$"
)

stage5_nocure_classification_csv <- find_stage_file(
  stage5_export_path,
  "^stage5_model_performance_long\\.csv$"
)



pnu_site_label <- "PNU"
snu_site_label <- "SNU"

horizons_year <- 1:10
risk_thresholds <- c(0.05, 0.10, 0.15, 0.20)

include_merged_incidence_site_supplementary <- TRUE
run_model_ids <- NULL

reuse_existing_stage8_outputs <- TRUE
require_existing_rds_to_skip_fit <- TRUE
preserve_existing_diagnostic_pdf <- TRUE

stan_chains <- 4L
stan_iter <- 2000L
stan_warmup <- 1000L
stan_thin <- 1L
stan_seed <- 20260322L
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

prior_predictive_horizons <- c(1, 2, 5, 10)
prior_tail_warning_mean_years <- 100
prior_tail_warning_q975_years <- 200

# 🔴 Initialize: packages and runtime defaults ===============================
core_packages <- c(
  "readr", "dplyr", "tibble", "tidyr", "ggplot2", "survival", "matrixStats"
)

missing_core_packages <- core_packages[
  !vapply(core_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_core_packages) > 0) {
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

if (!identical(as.integer(horizons_year), 1:10)) {
  stop("`horizons_year` must be exactly 1:10 for Stage 8.", call. = FALSE)
}
if (anyNA(risk_thresholds) || any(risk_thresholds <= 0 | risk_thresholds >= 1)) {
  stop("`risk_thresholds` must be probabilities strictly between 0 and 1.", call. = FALSE)
}
if (toupper(pnu_site_label) == toupper(snu_site_label)) {
  stop("`pnu_site_label` and `snu_site_label` must differ.", call. = FALSE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# 🔴 Define: reusable helpers and Bayesian engine ===============================
## 🟠 Define: I/O, coercion, reuse, and join helpers ===============================
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

read_delimited_or_rds <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(NULL)
  }
  if (!file.exists(path)) {
    stop(sprintf("Input file does not exist: %s", path), call. = FALSE)
  }
  path_low <- tolower(path)
  ext <- tolower(tools::file_ext(path))
  if (grepl("\\.(csv|txt)(\\.gz)?$", path_low)) {
    return(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  }
  if (ext == "rds") {
    return(readRDS(path))
  }
  stop(sprintf("Unsupported file extension for `%s`.", path), call. = FALSE)
}

read_existing_table_if_present <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  read_delimited_or_rds(path)
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
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
      return(one[[1]])
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

  vapply(
    values,
    function(one) {
      if (is.na(one)) {
        return(NA_character_)
      }
      as.character(one)
    },
    character(1)
  )
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

first_existing_name <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    return(NULL)
  }
  hit[[1]]
}

bind_rows_safe <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) == 0) {
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
  x %>%
    drop_existing_columns(join_cols) %>%
    left_join(y, by = by)
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
  isTRUE(!is.na(info$size[[1]]) && info$size[[1]] > 0)
}

subset_model_table <- function(df, model_id) {
  if (is.null(df) || !("model_id" %in% names(df))) {
    return(tibble())
  }
  df %>% filter(.data$model_id == .env$model_id)
}

subset_model_registry_row <- function(df, model_id) {
  out <- subset_model_table(df, model_id)
  if (nrow(out) == 0) {
    return(NULL)
  }
  out
}

make_threshold_key <- function(x) {
  out <- rep("__NA__", length(x))
  idx <- !is.na(x)
  out[idx] <- sprintf("%.10f", as.numeric(x[idx]))
  out
}

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
    fit_object = NULL,
    rds_path = NA_character_
  )
}

get_model_rds_candidates <- function(model_id, export_dir, registry_row = NULL) {
  cands <- c(
    file.path(export_dir, paste0(model_id, "__bayes_stage8_fit.rds")),
    if (!is.null(registry_row) && "rds_path" %in% names(registry_row)) as.character(registry_row$rds_path[[1]]) else NA_character_
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
      if (!file.exists(cand)) {
        next
      }
      obj <- tryCatch(readRDS(cand), error = function(e) NULL)
      if (is.null(obj)) {
        next
      }
      out <- extract_stage8_bundle_tables(obj)
      out$rds_path <- cand
      break
    }
  }

  assign(cache_key, out, envir = reuse_bundle_cache)
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

empty_delta_template <- function() {
  tibble(
    dataset = character(),
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
    horizon_support_label = character(),
    support_priority = integer(),
    reporting_status = character(),
    reporting_note = character(),
    classification_estimable_flag = logical(),
    denom_case = double(),
    denom_control = double(),
    admissible_flag = logical(),
    admissibility_reasons = character(),
    screening_flag = character(),
    screening_detail = character()
  )
}

empty_screening_template <- function() {
  tibble(
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
    screening_note = character(),
    screening_flag = character(),
    screening_detail = character(),
    carry_forward_stage8 = logical()
  )
}

empty_screening_lookup_row <- function(model_id = NA_character_) {
  tibble(
    model_id = as.character(model_id),
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

screening_value_or_na <- function(screening_row, field) {
  if (is.null(screening_row) || !inherits(screening_row, "data.frame") || nrow(screening_row) == 0L) {
    return(NA_character_)
  }
  if (!(field %in% names(screening_row))) {
    return(NA_character_)
  }
  value <- screening_row[[field]][[1]]
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
  value <- screening_row[[field]][[1]]
  if (is.null(value) || length(value) == 0L) {
    return(NA)
  }
  as.logical(value)
}

stage8_normalize_dataset_token <- function(x) {
  x_chr <- trimws(as.character(x))
  x_up <- toupper(x_chr)
  out <- dplyr::case_when(
    x_up %in% c("PNU", "P") ~ "PNU",
    x_up %in% c("SNU", "S") ~ "SNU",
    x_up %in% c("MERGED", "M", "BOTH", "ALL_SITES") ~ "merged",
    TRUE ~ normalize_dataset_label(x_chr)
  )
  out[is.na(x)] <- NA_character_
  out
}

stage8_dataset_from_stage6_key <- function(x) {
  x_chr <- trimws(as.character(x))
  x_base <- sub("__.*$", "", x_chr)
  stage8_normalize_dataset_token(x_base)
}

normalize_stage6_formula_anchor <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr[is.na(x_chr) | x_chr == ""] <- "ALL"
  x_key <- toupper(gsub("[^A-Z0-9]+", "", x_chr))
  out <- dplyr::case_when(
    x_key %in% c("ALL", "ANY", "GLOBAL") ~ "ALL",
    x_key %in% c("BASE", "MAIN", "BASELINE") ~ "base",
    TRUE ~ x_chr
  )
  out[is.na(out) | out == ""] <- "ALL"
  out
}

stage8_screening_lookup_fields <- function() {
  c(
    "cure_model_eligibility_flag",
    "primary_gate_method",
    "primary_gate_flag",
    "receus_primary_class",
    "presence_modifier_flag",
    "cure_presence_support_flag",
    "followup_contradiction_flag",
    "followup_not_contradicted_flag",
    "screening_note",
    "screening_flag",
    "screening_detail",
    "carry_forward_stage8"
  )
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

stage8_flag_chr <- function(x, preserve_unparsed = TRUE) {
  raw_chr <- trimws(as.character(x))
  parsed <- stage8_parse_logicalish(raw_chr)
  out <- ifelse(is.na(parsed), NA_character_, ifelse(parsed, "TRUE", "FALSE"))

  if (isTRUE(preserve_unparsed)) {
    keep_raw <- is.na(parsed) & !is.na(raw_chr) & nzchar(raw_chr) & !(toupper(raw_chr) %in% c("NA", "NULL"))
    out[keep_raw] <- raw_chr[keep_raw]
  }

  out
}

stage8_invert_flag_chr <- function(x) {
  parsed <- stage8_parse_logicalish(x)
  ifelse(is.na(parsed), NA_character_, ifelse(!parsed, "TRUE", "FALSE"))
}

stage8_extract_followup_flags <- function(df, followup_not_col = NULL, followup_contradiction_col = NULL) {
  followup_not_contradicted_flag <- if (!is.null(followup_not_col)) {
    stage8_flag_chr(df[[followup_not_col]], preserve_unparsed = TRUE)
  } else if (!is.null(followup_contradiction_col)) {
    stage8_invert_flag_chr(df[[followup_contradiction_col]])
  } else {
    rep(NA_character_, nrow(df))
  }

  followup_contradiction_flag <- if (!is.null(followup_contradiction_col)) {
    stage8_flag_chr(df[[followup_contradiction_col]], preserve_unparsed = TRUE)
  } else if (!is.null(followup_not_col)) {
    stage8_invert_flag_chr(df[[followup_not_col]])
  } else {
    rep(NA_character_, nrow(df))
  }

  tibble(
    followup_contradiction_flag = followup_contradiction_flag,
    followup_not_contradicted_flag = followup_not_contradicted_flag
  )
}

empty_prior_predictive_summary <- function() {
  tibble(
    dataset = character(),
    model_id = character(),
    metric = character(),
    horizon_year = integer(),
    prior_set = character(),
    mean = double(),
    q025 = double(),
    q50 = double(),
    q975 = double(),
    prior_degenerate_flag = logical()
  )
}

normalize_prior_predictive_summary <- function(df) {
  out <- if (is.null(df) || !inherits(df, "data.frame")) {
    empty_prior_predictive_summary()
  } else {
    tibble::as_tibble(df)
  }

  required_cols <- c(
    "dataset",
    "model_id",
    "metric",
    "horizon_year",
    "prior_set",
    "mean",
    "q025",
    "q50",
    "q975",
    "prior_degenerate_flag"
  )

  if (nrow(out) == 0L && length(setdiff(required_cols, names(out))) == length(required_cols)) {
    return(empty_prior_predictive_summary())
  }

  for (nm in setdiff(required_cols, names(out))) {
    if (nm == "horizon_year") {
      out[[nm]] <- NA_integer_
    } else if (nm %in% c("mean", "q025", "q50", "q975")) {
      out[[nm]] <- NA_real_
    } else if (nm == "prior_degenerate_flag") {
      out[[nm]] <- NA
    } else {
      out[[nm]] <- NA_character_
    }
  }

  out %>%
    select(any_of(required_cols), everything()) %>%
    mutate(
      dataset = as.character(.data$dataset),
      model_id = as.character(.data$model_id),
      metric = as.character(.data$metric),
      horizon_year = as.integer(safe_numeric(.data$horizon_year)),
      prior_set = as.character(.data$prior_set),
      mean = safe_numeric(.data$mean),
      q025 = safe_numeric(.data$q025),
      q50 = safe_numeric(.data$q50),
      q975 = safe_numeric(.data$q975),
      prior_degenerate_flag = as.logical(.data$prior_degenerate_flag)
    ) %>%
    distinct(dataset, model_id, metric, horizon_year, prior_set, .keep_all = TRUE)
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
  message(
    format_stage8_progress(done, total),
    " ",
    as.character(model_id),
    " ",
    as.character(detail)
  )
}

elapsed_stage8_seconds <- function(start_time) {
  as.numeric(difftime(Sys.time(), start_time, units = "secs"))
}

is_localhost_connection_warning <- function(w) {
  grepl("closing unused connection .*<-localhost:", conditionMessage(w))
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
    if ("waic" %in% rownames(waic_obj$estimates)) {
      out$waic <- as.numeric(waic_obj$estimates["waic", "Estimate"])
    }
    if ("p_waic" %in% rownames(waic_obj$estimates)) {
      out$p_waic <- as.numeric(waic_obj$estimates["p_waic", "Estimate"])
    }
    if (!is.null(waic_obj$pointwise) && "p_waic" %in% colnames(waic_obj$pointwise)) {
      p_waic_pointwise <- as.numeric(waic_obj$pointwise[, "p_waic"])
      out$p_waic_high_n <- as.integer(sum(p_waic_pointwise > 0.4, na.rm = TRUE))
      out$p_waic_high_pct <- if (length(p_waic_pointwise) > 0L) {
        100 * out$p_waic_high_n / length(p_waic_pointwise)
      } else {
        NA_real_
      }
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
      out$pareto_k_max <- if (length(pareto_k) > 0L && any(is.finite(pareto_k))) {
        max(pareto_k, na.rm = TRUE)
      } else {
        NA_real_
      }
      out$pareto_k_bad_n <- as.integer(sum(pareto_k > 0.7, na.rm = TRUE))
      out$pareto_k_very_bad_n <- as.integer(sum(pareto_k > 1.0, na.rm = TRUE))
      out$pareto_k_bad_pct <- if (length(pareto_k) > 0L) {
        100 * out$pareto_k_bad_n / length(pareto_k)
      } else {
        NA_real_
      }
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

load_existing_stage8_exports <- function(export_dir) {
  list(
    model_registry = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_model_registry.csv")),
    coefficient_summary = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_coefficient_summary.csv")),
    posterior_subject_profile = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_posterior_subject_profile.csv.gz")),
    posterior_subject_yearly = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_posterior_subject_yearly.csv.gz")),
    posterior_cohort_yearly = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_posterior_cohort_yearly.csv")),
    posterior_classification = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_posterior_classification.csv")),
    diagnostics_parameter_level = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_diagnostics_parameter_level.csv")),
    ppc_summary = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_ppc_summary.csv")),
    prior_predictive_summary = read_existing_table_if_present(file.path(export_dir, "bayes_stage8_prior_predictive_summary.csv")),
    diagnostic_pdf_exists = pdf_file_is_usable(file.path(export_dir, "bayes_stage8_diagnostic_plots.pdf"))
  )
}

is_model_reusable <- function(model_row, dataset_df, existing_exports, export_dir, require_existing_rds = TRUE) {
  reg <- subset_model_registry_row(existing_exports$model_registry, model_row$model_id[[1]])
  if (is.null(reg) || nrow(reg) != 1L) {
    return(FALSE)
  }

  if (!identical(as.character(reg$dataset[[1]]), as.character(model_row$dataset[[1]]))) {
    return(FALSE)
  }
  if (!identical(as.character(reg$family_code[[1]]), as.character(model_row$family_code[[1]]))) {
    return(FALSE)
  }
  if (!identical(as.character(reg$formula_anchor[[1]]), as.character(model_row$formula_anchor[[1]]))) {
    return(FALSE)
  }
  if (!identical(as.character(reg$fit_status[[1]]), "ok")) {
    return(FALSE)
  }

  n_event_now <- sum(dataset_df$event_main)
  n_censor_now <- sum(dataset_df$censor_main)
  n_remission_now <- sum(dataset_df$remission_flag)

  if (!isTRUE(nrow(dataset_df) == reg$n[[1]])) {
    return(FALSE)
  }
  if (!isTRUE(n_event_now == reg$n_event[[1]])) {
    return(FALSE)
  }
  if (!isTRUE(n_censor_now == reg$n_censor_main[[1]])) {
    return(FALSE)
  }
  if (!isTRUE(n_remission_now == reg$n_remission[[1]])) {
    return(FALSE)
  }

  rds_candidates <- get_model_rds_candidates(model_row$model_id[[1]], export_dir, reg)
  has_rds <- any(file.exists(rds_candidates))
  if (isTRUE(require_existing_rds) && !has_rds) {
    return(FALSE)
  }

  reuse_bundle <- if (has_rds) {
    get_reuse_bundle(model_row$model_id[[1]], export_dir, reg)
  } else {
    empty_reuse_bundle()
  }

  required_nonadmissible <- c("coefficient_summary", "diagnostics_parameter_level", "ppc_summary")
  for (nm in required_nonadmissible) {
    has_existing <- nrow(subset_model_table(existing_exports[[nm]], model_row$model_id[[1]])) > 0L
    has_bundle <- nrow(subset_model_table(reuse_bundle[[nm]], model_row$model_id[[1]])) > 0L
    if (!has_existing && !has_bundle) {
      return(FALSE)
    }
  }

  if (isTRUE(as.logical(reg$admissible_flag[[1]]))) {
    required_admissible <- c(
      "posterior_subject_profile",
      "posterior_subject_yearly",
      "posterior_cohort_yearly",
      "posterior_classification"
    )
    for (nm in required_admissible) {
      has_existing <- nrow(subset_model_table(existing_exports[[nm]], model_row$model_id[[1]])) > 0L
      has_bundle <- nrow(subset_model_table(reuse_bundle[[nm]], model_row$model_id[[1]])) > 0L
      if (!has_existing && !has_bundle) {
        return(FALSE)
      }
    }
  }

  TRUE
}

## 🟠 Define: backbone data preparation ===============================
prepare_analysis_dataset <- function(df, dataset_name, pnu_label, snu_label) {
  if (is.null(df)) {
    stop(sprintf("[%s] Input data frame is NULL.", dataset_name), call. = FALSE)
  }

  if (!("site" %in% names(df))) {
    df$site <- dataset_name
  }

  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      sprintf("[%s] Missing required columns: %s", dataset_name, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
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

  if (nrow(out) == 0) {
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
      censor_main = as.integer(status_num %in% c(0L, 2L)),
      remission_flag = as.integer(status_num == 2L),
      age_s = (age_exact_entry - age_mean) / (2 * age_sd)
    )

  if (nrow(out) != dplyr::n_distinct(out$unique_person_id)) {
    stop(sprintf("[%s] `site + id` is not unique.", dataset_name), call. = FALSE)
  }

  if (dataset_name == "merged" && dplyr::n_distinct(out$site) < 2L) {
    stop("[merged] Merged dataset must contain at least two site levels.", call. = FALSE)
  }

  list(
    data = out,
    scaling = tibble(
      dataset = dataset_name,
      center_mean_age_exact_entry = age_mean,
      sd_age_exact_entry = age_sd,
      scale_rule = "(age_exact_entry - mean(age_exact_entry)) / (2 * sd(age_exact_entry))"
    )
  )
}

load_stage8_datasets <- function(merged_path, pnu_label, snu_label) {
  merged_raw <- read_delimited_or_rds(merged_path)

  if (is.null(merged_raw)) {
    stop("Provide `merged_data_path`.", call. = FALSE)
  }
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

  if (nrow(pnu_raw) == 0) {
    stop(sprintf("No rows found for PNU site label `%s` in merged data.", pnu_label), call. = FALSE)
  }
  if (nrow(snu_raw) == 0) {
    stop(sprintf("No rows found for SNU site label `%s` in merged data.", snu_label), call. = FALSE)
  }

  prepared_pnu <- prepare_analysis_dataset(pnu_raw, "PNU", pnu_label, snu_label)
  prepared_snu <- prepare_analysis_dataset(snu_raw, "SNU", pnu_label, snu_label)
  prepared_merged <- prepare_analysis_dataset(merged_raw, "merged", pnu_label, snu_label)

  list(
    datasets = list(
      PNU = prepared_pnu$data,
      SNU = prepared_snu$data,
      merged = prepared_merged$data
    ),
    scaling = bind_rows(prepared_pnu$scaling, prepared_snu$scaling, prepared_merged$scaling)
  )
}

## 🟠 Define: stage-link normalization ===============================
read_screening_flags <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(empty_screening_template())
  }

  raw_obj <- read_delimited_or_rds(path)
  df <- extract_stage6_screening_table(raw_obj)
  if (is.null(df) && inherits(raw_obj, "data.frame")) {
    df <- raw_obj
  }
  if (is.null(df) || !inherits(df, "data.frame")) {
    return(empty_screening_template())
  }
  df <- tibble::as_tibble(df)

  dataset_source_col <- first_existing_name(df, c("source_dataset", "dataset", "cohort"))
  dataset_key_col <- first_existing_name(df, c("dataset_key", "screening_dataset_key"))
  analysis_variant_col <- first_existing_name(df, c("analysis_variant"))
  branch_col <- first_existing_name(df, c("hsu_formula_branch", "analysis_variant"))
  formula_col <- first_existing_name(df, c("formula_anchor", "formula_name", "formula_type"))
  eligibility_col <- first_existing_name(df, c("cure_model_eligibility_flag", "stage6_final_class", "final_decision_flag", "decision_flag", "screening_flag"))
  primary_gate_method_col <- first_existing_name(df, c("primary_gate_method", "method", "screening_method", "analysis_variant"))
  primary_gate_flag_col <- first_existing_name(df, c("primary_gate_flag"))
  receus_col <- first_existing_name(df, c("receus_primary_class"))
  presence_modifier_col <- first_existing_name(df, c("presence_modifier_flag"))
  cure_presence_col <- first_existing_name(df, c("cure_presence_support_flag"))
  followup_not_col <- first_existing_name(df, c("followup_not_contradicted_flag"))
  followup_contradiction_col <- first_existing_name(df, c("followup_contradiction_flag"))
  note_col <- first_existing_name(df, c("screening_note", "screening_detail", "screening_context", "detail"))

  if (is.null(eligibility_col) || (is.null(dataset_source_col) && is.null(dataset_key_col))) {
    return(empty_screening_template())
  }

  followup_flags <- stage8_extract_followup_flags(
    df = df,
    followup_not_col = followup_not_col,
    followup_contradiction_col = followup_contradiction_col
  )

  has_stage6_variant_map <- !is.null(dataset_key_col) || !is.null(dataset_source_col)

  out <- if (isTRUE(has_stage6_variant_map)) {
    bind_rows(lapply(seq_len(nrow(df)), function(i) {
      one <- df[i, , drop = FALSE]
      dataset_std <- if (!is.null(dataset_source_col)) {
        stage8_normalize_dataset_token(one[[dataset_source_col]][[1]])
      } else {
        stage8_dataset_from_stage6_key(one[[dataset_key_col]][[1]])
      }
      formula_anchors <- if (!is.null(formula_col)) {
        normalize_stage6_formula_anchor(one[[formula_col]][[1]])
      } else {
        map_stage6_variant_to_formula_anchor(
          dataset_key = if (!is.null(dataset_key_col)) one[[dataset_key_col]][[1]] else NA_character_,
          source_dataset = dataset_std,
          analysis_variant = if (!is.null(analysis_variant_col)) one[[analysis_variant_col]][[1]] else NA_character_,
          hsu_formula_branch = if (!is.null(branch_col)) one[[branch_col]][[1]] else NA_character_
        )
      }

      note_parts <- c(
        if (!is.null(primary_gate_method_col)) paste0("method:", as.character(one[[primary_gate_method_col]][[1]])) else NA_character_,
        if (!is.null(note_col)) as.character(one[[note_col]][[1]]) else NA_character_
      )

      tibble(
        dataset = dataset_std,
        formula_anchor = formula_anchors,
        cure_model_eligibility_flag = as.character(one[[eligibility_col]][[1]]),
        primary_gate_method = if (!is.null(primary_gate_method_col)) as.character(one[[primary_gate_method_col]][[1]]) else NA_character_,
        primary_gate_flag = if (!is.null(primary_gate_flag_col)) as.character(one[[primary_gate_flag_col]][[1]]) else NA_character_,
        receus_primary_class = if (!is.null(receus_col)) as.character(one[[receus_col]][[1]]) else NA_character_,
        presence_modifier_flag = if (!is.null(presence_modifier_col)) as.character(one[[presence_modifier_col]][[1]]) else NA_character_,
        cure_presence_support_flag = if (!is.null(cure_presence_col)) as.character(one[[cure_presence_col]][[1]]) else NA_character_,
        followup_contradiction_flag = followup_flags$followup_contradiction_flag[[i]],
        followup_not_contradicted_flag = followup_flags$followup_not_contradicted_flag[[i]],
        screening_note = if (!is.null(note_col)) as.character(one[[note_col]][[1]]) else NA_character_,
        screening_detail_component = stage8_collapse_unique_text(note_parts, sep = "; ")
      )
    }))
  } else {
    tibble(
      dataset = stage8_normalize_dataset_token(df[[dataset_source_col]]),
      formula_anchor = if (!is.null(formula_col)) normalize_stage6_formula_anchor(df[[formula_col]]) else "ALL",
      cure_model_eligibility_flag = as.character(df[[eligibility_col]]),
      primary_gate_method = if (!is.null(primary_gate_method_col)) as.character(df[[primary_gate_method_col]]) else NA_character_,
      primary_gate_flag = if (!is.null(primary_gate_flag_col)) as.character(df[[primary_gate_flag_col]]) else NA_character_,
      receus_primary_class = if (!is.null(receus_col)) as.character(df[[receus_col]]) else NA_character_,
      presence_modifier_flag = if (!is.null(presence_modifier_col)) as.character(df[[presence_modifier_col]]) else NA_character_,
      cure_presence_support_flag = if (!is.null(cure_presence_col)) as.character(df[[cure_presence_col]]) else NA_character_,
      followup_contradiction_flag = followup_flags$followup_contradiction_flag,
      followup_not_contradicted_flag = followup_flags$followup_not_contradicted_flag,
      screening_note = if (!is.null(note_col)) as.character(df[[note_col]]) else NA_character_,
      screening_detail_component = stage8_collapse_unique_text(
        c(
          if (!is.null(primary_gate_method_col)) paste0("method:", as.character(df[[primary_gate_method_col]])) else NA_character_,
          if (!is.null(note_col)) as.character(df[[note_col]]) else NA_character_
        ),
        sep = "; "
      )
    )
  }

  out <- out %>%
    mutate(
      dataset = stage8_normalize_dataset_token(dataset),
      formula_anchor = normalize_stage6_formula_anchor(formula_anchor)
    ) %>%
    filter(dataset %in% c("PNU", "SNU", "merged")) %>%
    group_by(dataset, formula_anchor) %>%
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
      screening_detail = stage8_collapse_unique_text(screening_detail_component, sep = "; "),
      .groups = "drop"
    ) %>%
    mutate(
      screening_flag = cure_model_eligibility_flag,
      carry_forward_stage8 = dplyr::case_when(
        is.na(cure_model_eligibility_flag) ~ NA,
        grepl("\\|", cure_model_eligibility_flag, fixed = TRUE) ~ NA,
        cure_model_eligibility_flag == "unsupportive" ~ FALSE,
        TRUE ~ TRUE
      )
    ) %>%
    select(
      dataset,
      formula_anchor,
      all_of(stage8_screening_lookup_fields())
    )

  if (nrow(out) == 0L) {
    return(empty_screening_template())
  }

  out
}

build_screening_model_lookup <- function(screening_flags, model_grid) {
  model_base <- model_grid %>%
    distinct(model_id, dataset, formula_anchor)

  lookup_fields <- stage8_screening_lookup_fields()

  if (is.null(screening_flags) || nrow(screening_flags) == 0) {
    out <- model_base %>% select(model_id)
    template <- empty_screening_lookup_row()
    for (nm in setdiff(names(template), "model_id")) {
      out[[nm]] <- template[[nm]][[1]]
    }
    return(out)
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
    ) %>%
    select(-dataset, -any_of(c("formula_anchor", "formula_anchor_all")))

  out <- exact_lookup %>%
    left_join(all_lookup, by = "model_id") %>%
    transmute(model_id = model_id)

  for (nm in lookup_fields) {
    out[[nm]] <- dplyr::coalesce(exact_lookup[[nm]], all_lookup[[paste0(nm, "_all")]])
  }

  out
}

read_stage6_screening_lookup <- function(path) {
  read_screening_flags(path)
}

stage8b_support_tier_from_grid <- function(dataset, horizon_year) {
  dataset <- stage8_normalize_dataset_token(dataset)
  horizon_year <- as.integer(safe_numeric(horizon_year))

  dplyr::case_when(
    dataset == "PNU" & horizon_year == 1L ~ "primary_supported",
    dataset == "PNU" & horizon_year == 2L ~ "sensitivity",
    dataset == "PNU" & horizon_year >= 3L ~ "projection",
    dataset %in% c("SNU", "merged") & horizon_year %in% c(1L, 2L) ~ "primary_supported",
    dataset %in% c("SNU", "merged") & horizon_year <= 5L ~ "secondary",
    TRUE ~ "projection"
  )
}

stage8b_horizon_evidence_from_grid <- function(dataset, horizon_year) {
  dataset <- stage8_normalize_dataset_token(dataset)
  horizon_year <- as.integer(safe_numeric(horizon_year))

  dplyr::case_when(
    dataset == "PNU" & horizon_year == 1L ~ "directly_observed_data_supported",
    dataset == "PNU" & horizon_year == 2L ~ "partly_model_dependent",
    dataset == "PNU" & horizon_year >= 3L ~ "mostly_extrapolated",
    dataset %in% c("SNU", "merged") & horizon_year %in% c(1L, 2L) ~ "directly_observed_data_supported",
    dataset %in% c("SNU", "merged") & horizon_year <= 5L ~ "partly_model_dependent",
    TRUE ~ "mostly_extrapolated"
  )
}

stage8b_claim_restriction_from_evidence <- function(horizon_evidence_class) {
  dplyr::case_when(
    horizon_evidence_class == "directly_observed_data_supported" ~ "primary_claim_allowed",
    horizon_evidence_class == "partly_model_dependent" ~ "secondary_or_sensitivity_only",
    TRUE ~ "projection_only"
  )
}

stage8b_interpretation_note_from_grid <- function(dataset, horizon_year, support_tier, horizon_evidence_class, claim_restriction_flag) {
  dataset <- stage8_normalize_dataset_token(dataset)
  horizon_year <- as.integer(safe_numeric(horizon_year))

  dplyr::case_when(
    claim_restriction_flag == "primary_claim_allowed" ~ "Primary supported horizon with comparatively direct follow-up support.",
    dataset == "PNU" & horizon_year == 2L ~ "Sensitivity horizon for PNU; partly model-dependent and not for primary claims.",
    claim_restriction_flag == "secondary_or_sensitivity_only" ~ "Secondary or sensitivity horizon with explicit model dependence.",
    TRUE ~ "Projection-dominant horizon."
  )
}

build_horizon_annotation_from_stage1 <- function(horizon_registry, datasets, horizons) {
  dataset_levels <- if (is.list(datasets) && !is.null(names(datasets)) && any(nzchar(names(datasets)))) {
    names(datasets)
  } else {
    as.character(datasets)
  }
  dataset_levels <- unique(stage8_normalize_dataset_token(dataset_levels))

  out <- if (is.null(horizon_registry)) tibble() else tibble::as_tibble(horizon_registry)

  if (nrow(out) == 0L) {
    out <- tidyr::crossing(
      dataset = dataset_levels,
      horizon_year = as.integer(horizons)
    )
  }

  if (!("dataset" %in% names(out))) {
    if ("dataset_key" %in% names(out)) {
      out$dataset <- out$dataset_key
    } else {
      stop("Stage1 horizon registry must contain `dataset` or `dataset_key`.", call. = FALSE)
    }
  }

  if (!("dataset_key" %in% names(out))) {
    out$dataset_key <- out$dataset
  }

  if (!("horizon_year" %in% names(out))) {
    horizon_col <- first_existing_name(out, c("horizon", "year"))
    if (is.null(horizon_col)) {
      stop("Stage1 horizon registry must contain `horizon_year` (or `horizon` / `year`).", call. = FALSE)
    }
    out$horizon_year <- out[[horizon_col]]
  }

  if (!("support_tier" %in% names(out))) out$support_tier <- NA_character_
  if (!("support_tier_standard" %in% names(out))) out$support_tier_standard <- NA_character_
  if (!("horizon_evidence_class" %in% names(out))) out$horizon_evidence_class <- NA_character_
  if (!("claim_restriction_flag" %in% names(out))) out$claim_restriction_flag <- NA_character_
  if (!("interpretation_note" %in% names(out))) out$interpretation_note <- NA_character_

  out %>%
    mutate(
      dataset = stage8_normalize_dataset_token(dataset),
      dataset_key = dplyr::coalesce(stage8_normalize_dataset_token(dataset_key), stage8_normalize_dataset_token(dataset)),
      horizon_year = as.integer(safe_numeric(horizon_year)),
      support_tier = dplyr::coalesce(
        dplyr::na_if(as.character(support_tier), ""),
        dplyr::na_if(as.character(support_tier_standard), ""),
        stage8b_support_tier_from_grid(dataset, horizon_year)
      ),
      horizon_evidence_class = dplyr::coalesce(
        dplyr::na_if(as.character(horizon_evidence_class), ""),
        stage8b_horizon_evidence_from_grid(dataset, horizon_year)
      ),
      claim_restriction_flag = dplyr::coalesce(
        dplyr::na_if(as.character(claim_restriction_flag), ""),
        stage8b_claim_restriction_from_evidence(horizon_evidence_class)
      ),
      interpretation_note = dplyr::coalesce(
        dplyr::na_if(as.character(interpretation_note), ""),
        stage8b_interpretation_note_from_grid(dataset, horizon_year, support_tier, horizon_evidence_class, claim_restriction_flag)
      )
    ) %>%
    transmute(
      dataset = as.character(dataset),
      dataset_key = as.character(dataset_key),
      horizon_year = as.integer(horizon_year),
      support_tier = as.character(support_tier),
      horizon_evidence_class = as.character(horizon_evidence_class),
      claim_restriction_flag = as.character(claim_restriction_flag),
      interpretation_note = as.character(interpretation_note)
    ) %>%
    filter(dataset %in% dataset_levels, horizon_year %in% as.integer(horizons)) %>%
    distinct(dataset, horizon_year, .keep_all = TRUE) %>%
    arrange(factor(dataset, levels = dataset_levels), horizon_year)
}

make_stage8b_support_registry <- function(horizon_registry_stage1, datasets = c("PNU", "SNU", "merged"), horizons = 1:10) {
  dataset_levels <- unique(stage8_normalize_dataset_token(datasets))

  out <- build_horizon_annotation_from_stage1(
    horizon_registry = horizon_registry_stage1,
    datasets = dataset_levels,
    horizons = horizons
  )

  if (nrow(out) == 0L) {
    return(
      tidyr::expand_grid(dataset = dataset_levels, horizon_year = as.integer(horizons)) %>%
        mutate(
          dataset_key = dataset,
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
      dataset_key = as.character(dataset_key),
      horizon_year = as.integer(safe_numeric(horizon_year)),
      support_tier = as.character(support_tier),
      support_tier_standard = as.character(support_tier),
      horizon_evidence_class = as.character(horizon_evidence_class),
      claim_restriction_flag = as.character(claim_restriction_flag),
      interpretation_note = as.character(interpretation_note)
    ) %>%
    filter(dataset %in% dataset_levels, horizon_year %in% as.integer(horizons)) %>%
    distinct(dataset, horizon_year, .keep_all = TRUE) %>%
    arrange(factor(dataset, levels = dataset_levels), horizon_year)
}

normalize_nocure_cohort <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  raw_obj <- read_delimited_or_rds(path)
  df <- extract_stage5_performance_table(raw_obj)
  if (is.null(df) && inherits(raw_obj, "data.frame")) {
    df <- raw_obj
  }
  if (is.null(df) || !inherits(df, "data.frame")) {
    return(NULL)
  }

  if (all(c("metric_domain", "metric_name", "metric_value") %in% names(df))) {
    out <- df %>%
      filter(
        model_class == "non_cure",
        metric_domain == "horizon_summary",
        metric_name == "mean_predicted_risk"
      ) %>%
      transmute(
        dataset = normalize_dataset_label(dataset),
        formula_anchor = if ("formula_name" %in% names(df)) as.character(formula_name) else "ALL",
        no_cure_model_id = as.character(model_id),
        horizon_year = as.integer(safe_numeric(horizon_year)),
        metric = "meanRisk",
        value = safe_numeric(metric_value)
      ) %>%
      filter(!is.na(horizon_year), !is.na(value))

    if (nrow(out) > 0) {
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
  if (is.null(df) && inherits(raw_obj, "data.frame")) {
    df <- raw_obj
  }
  if (is.null(df) || !inherits(df, "data.frame")) {
    return(NULL)
  }

  if (all(c("metric_domain", "metric_name", "metric_value") %in% names(df))) {
    metric_name_map <- c(
      false_positive_burden_nonevents = "FPR",
      false_positive_per_100 = "FP100",
      net_benefit = "NB",
      false_positive_burden_all = "false_positive_burden"
    )

    out <- df %>%
      filter(
        model_class == "non_cure",
        metric_domain == "threshold_summary",
        metric_name %in% names(metric_name_map)
      ) %>%
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

    if (nrow(out) > 0) {
      return(out)
    }
  }

  dataset_col <- first_existing_name(df, c("dataset", "cohort"))
  horizon_col <- first_existing_name(df, c("horizon_year", "horizon", "year"))
  threshold_col <- first_existing_name(df, c("threshold", "risk_threshold"))
  formula_col <- first_existing_name(df, c("formula_anchor", "formula_name", "formula_type"))
  model_col <- first_existing_name(df, c("no_cure_model_id", "model_id", "family", "model"))

  metric_candidates <- c(
    FPR = first_existing_name(df, c("FPR", "fpr", "fpr_mean")),
    FP100 = first_existing_name(df, c("FP100", "fp100", "false_positives_per_100")),
    NB = first_existing_name(df, c("NB", "net_benefit", "nb_mean")),
    false_positive_burden = first_existing_name(df, c("false_positive_burden", "fp_burden"))
  )
  metric_candidates <- metric_candidates[!vapply(metric_candidates, is.null, logical(1))]

  if (is.null(dataset_col) || is.null(horizon_col) || is.null(threshold_col) || length(metric_candidates) == 0) {
    return(NULL)
  }

  out <- tibble(
    dataset = normalize_dataset_label(df[[dataset_col]]),
    formula_anchor = if (!is.null(formula_col)) as.character(df[[formula_col]]) else "ALL",
    no_cure_model_id = if (!is.null(model_col)) as.character(df[[model_col]]) else "NO_CURE_REFERENCE",
    horizon_year = as.integer(safe_numeric(df[[horizon_col]])),
    threshold = as.numeric(safe_numeric(df[[threshold_col]]))
  )

  metric_long <- lapply(names(metric_candidates), function(metric_name) {
    tibble(
      dataset = out$dataset,
      formula_anchor = out$formula_anchor,
      no_cure_model_id = out$no_cure_model_id,
      horizon_year = out$horizon_year,
      threshold = out$threshold,
      metric = metric_name,
      value = safe_numeric(df[[metric_candidates[[metric_name]]]])
    )
  })

  bind_rows(metric_long) %>%
    filter(!is.na(horizon_year), !is.na(threshold), !is.na(value))
}

## 🟠 Define: structural grid and prior objects ===============================
build_model_grid <- function(include_supplementary = TRUE) {
  base_rows <- tibble::tribble(
    ~dataset,  ~structural_model_id, ~latency_branch, ~formula_anchor,    ~incidence_site_indicator, ~latency_site_indicator, ~latency_interaction_indicator, ~is_supplementary_branch,
    "PNU",     "PNU-L0",             "L0",            "base",             FALSE,                      FALSE,                     FALSE,                         FALSE,
    "PNU",     "PNU-L1",             "L1",            "interaction",      FALSE,                      FALSE,                     TRUE,                          FALSE,
    "SNU",     "SNU-L0",             "L0",            "base",             FALSE,                      FALSE,                     FALSE,                         FALSE,
    "SNU",     "SNU-L1",             "L1",            "interaction",      FALSE,                      FALSE,                     TRUE,                          FALSE,
    "merged",  "MERGED-L0S0",        "L0S0",          "base",             FALSE,                      FALSE,                     FALSE,                         FALSE,
    "merged",  "MERGED-L1S0",        "L1S0",          "interaction",      FALSE,                      FALSE,                     TRUE,                          FALSE,
    "merged",  "MERGED-L0S1",        "L0S1",          "site_added",       FALSE,                      TRUE,                      FALSE,                         FALSE,
    "merged",  "MERGED-L1S1",        "L1S1",          "site_interaction", FALSE,                      TRUE,                      TRUE,                          FALSE
  )

  if (isTRUE(include_supplementary)) {
    base_rows <- bind_rows(
      base_rows,
      tibble::tribble(
        ~dataset,  ~structural_model_id, ~latency_branch, ~formula_anchor,    ~incidence_site_indicator, ~latency_site_indicator, ~latency_interaction_indicator, ~is_supplementary_branch,
        "merged",  "MERGED-Isite-L0S0",  "L0S0",          "base",             TRUE,                       FALSE,                     FALSE,                         TRUE,
        "merged",  "MERGED-Isite-L1S0",  "L1S0",          "interaction",      TRUE,                       FALSE,                     TRUE,                          TRUE,
        "merged",  "MERGED-Isite-L0S1",  "L0S1",          "site_added",       TRUE,                       TRUE,                      FALSE,                         TRUE,
        "merged",  "MERGED-Isite-L1S1",  "L1S1",          "site_interaction", TRUE,                       TRUE,                      TRUE,                          TRUE
      )
    )
  }

  family_rows <- tibble::tribble(
    ~family_code, ~latency_family, ~family_id,
    "E",          "exponential",   1L,
    "W",          "weibull",       2L,
    "LN",         "lognormal",     3L,
    "LL",         "loglogistic",   4L
  )

  tidyr::crossing(base_rows, family_rows) %>%
    mutate(
      model_id = paste(structural_model_id, family_code, sep = "-"),
      prior_set = "main"
    )
}

build_prior_specs <- function() {
  list(
    main = list(
      prior_set = "main",
      alpha_gp = -9.581369553169,
      mu_beta_inc_base = c(0.419871845822, 0.907608052926, 0.586202561451, 0.466865123863, 0.037997248763),
      sd_beta_inc_base = c(0.132789397422, 0.173731076538, 0.191221553945, 0.270393197518, 0.302838606651),
      sd_beta_inc_site = 0.5,
      sd_delta = 3.5,
      sd_gamma0 = 2.5,
      sd_gamma = 1.0,
      sd_shape = 1.0
    ),
    wider_latency = list(
      prior_set = "wider_latency",
      alpha_gp = -9.581369553169,
      mu_beta_inc_base = c(0.419871845822, 0.907608052926, 0.586202561451, 0.466865123863, 0.037997248763),
      sd_beta_inc_base = c(0.132789397422, 0.173731076538, 0.191221553945, 0.270393197518, 0.302838606651),
      sd_beta_inc_site = 0.5,
      sd_delta = 3.5,
      sd_gamma0 = 3.5,
      sd_gamma = 1.5,
      sd_shape = 1.5
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

  if (isTRUE(model_row$incidence_site_indicator)) {
    X_inc <- cbind(X_inc, site_SNU = s_i)
    mu_beta_inc <- c(mu_beta_inc, 0)
    sd_beta_inc <- c(sd_beta_inc, prior_spec$sd_beta_inc_site)
  }

  a_i <- as.numeric(df$age_s)
  az_i <- a_i * z_i

  X_lat <- switch(
    model_row$latency_branch,
    L0 = cbind(age_s = a_i, sex_num = z_i),
    L1 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i),
    L0S0 = cbind(age_s = a_i, sex_num = z_i),
    L1S0 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i),
    L0S1 = cbind(age_s = a_i, sex_num = z_i, site_SNU = s_i),
    L1S1 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i, site_SNU = s_i),
    stop(sprintf("Unknown latency branch `%s`.", model_row$latency_branch), call. = FALSE)
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
    id_df = df %>% select(unique_person_id, id, site, sex_num, age_exact_entry, age_s)
  )
}

## 🟠 Define: survival math, support labels, and metadata registries ===============================
support_label <- function(dataset_name, horizon_year) {
  if (dataset_name == "PNU") {
    if (horizon_year == 1) return("primary-supported")
    if (horizon_year == 2) return("secondary")
    return("projection")
  }
  if (dataset_name %in% c("SNU", "merged")) {
    if (horizon_year %in% c(1, 2)) return("primary-supported")
    if (horizon_year <= 5) return("secondary")
    return("projection")
  }
  "projection"
}

support_priority_from_label <- function(x) {
  x <- as.character(x)
  out <- rep.int(9L, length(x))
  out[x == "primary-supported"] <- 1L
  out[x == "secondary"] <- 2L
  out[x == "projection"] <- 3L
  out
}

ppc_horizons_for_dataset <- function(dataset_name) {
  if (dataset_name == "PNU") {
    return(c(1, 2))
  }
  c(1, 2, 5)
}

km_eval <- function(survfit_obj, times) {
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
      G_t = G_t
    ) %>%
      mutate(
        w_case = list(w_case),
        w_control = list(w_control)
      )
  })

  bind_rows(horizon_rows)
}

make_horizon_metadata_registry <- function(ipcw_registry) {
  bind_rows(lapply(names(ipcw_registry), function(ds) {
    ref <- ipcw_registry[[ds]]
    labels <- vapply(ref$horizon_year, function(h) support_label(ds, h), character(1))
    estimable <- (ref$denom_case > 0) & (ref$denom_control > 0)

    classification_note <- case_when(
      !estimable & ref$denom_case <= 0 & ref$denom_control <= 0 & labels == "projection" ~
        "IPCW case and control denominators are zero at this horizon. Threshold-based classification is not estimable; this horizon is also projection-only.",
      !estimable & ref$denom_control <= 0 & labels == "projection" ~
        "IPCW control denominator is zero at this horizon. Threshold-based classification is not estimable; this horizon is also projection-only.",
      !estimable & ref$denom_case <= 0 & labels == "projection" ~
        "IPCW case denominator is zero at this horizon. Threshold-based classification is not estimable; this horizon is also projection-only.",
      !estimable & ref$denom_case <= 0 & ref$denom_control <= 0 ~
        "IPCW case and control denominators are zero at this horizon. Threshold-based classification is not estimable.",
      !estimable & ref$denom_control <= 0 ~
        "IPCW control denominator is zero at this horizon. Threshold-based classification is not estimable.",
      !estimable & ref$denom_case <= 0 ~
        "IPCW case denominator is zero at this horizon. Threshold-based classification is not estimable.",
      labels == "projection" ~
        "Projection horizon. Threshold-based classification is reported for completeness and should not be prioritized in main interpretation.",
      labels == "secondary" ~
        "Secondary-support horizon.",
      TRUE ~
        "Primary-supported horizon."
    )

    cohort_note <- case_when(
      labels == "projection" ~
        "Projection horizon. Cohort-level posterior trajectories are reported, but primary interpretation should remain anchored to supported horizons.",
      labels == "secondary" ~
        "Secondary-support horizon.",
      TRUE ~
        "Primary-supported horizon."
    )

    tibble(
      dataset = ds,
      horizon_year = ref$horizon_year,
      horizon_support_label = labels,
      support_priority = support_priority_from_label(labels),
      observed_km_risk = ref$observed_km_risk,
      denom_case = ref$denom_case,
      denom_control = ref$denom_control,
      classification_estimable_flag = estimable,
      classification_reporting_status = case_when(
        !estimable ~ "not_estimable",
        labels == "projection" ~ "reported_with_projection_caution",
        TRUE ~ "reported"
      ),
      classification_reporting_note = classification_note,
      cohort_reporting_status = case_when(
        labels == "projection" ~ "reported_with_projection_caution",
        TRUE ~ "reported"
      ),
      cohort_reporting_note = cohort_note
    )
  })) %>%
    arrange(dataset, horizon_year)
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

  if (is.null(dim(beta_inc))) {
    beta_inc <- matrix(beta_inc, ncol = K_inc)
  }
  if (is.null(dim(gamma_lat))) {
    gamma_lat <- matrix(gamma_lat, ncol = K_lat)
  }
  if (is.null(dim(log_lik))) {
    log_lik <- matrix(log_lik, nrow = length(ext$delta0))
  }

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

  if (length(supported_risk_list) > 0) {
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
    rho_W = rnorm(n_draws, 0, prior_spec$sd_shape),
    log_sigma_LN = rnorm(n_draws, 0, prior_spec$sd_shape),
    psi_LL = rnorm(n_draws, 0, prior_spec$sd_shape)
  )

  lin <- compute_linear_terms(sim_draws, design_bundle$X_inc, design_bundle$X_lat, design_bundle$alpha_gp)

  supported_horizons <- ppc_horizons_for_dataset(model_row$dataset)
  supported_risk_list <- lapply(supported_horizons, function(h) {
    fh <- family_survival_hazard(h, model_row$family_code, lin$mu_lat_mat, lin$median_mat, sim_draws)
    1 - ((1 - lin$pi_mat) + lin$pi_mat * fh$Su)
  })

  degeneracy <- compute_degeneracy(lin$pi_mat, lin$median_mat, supported_risk_list)

  prior_rows <- list(
    tibble(
      dataset = model_row$dataset,
      model_id = model_row$model_id,
      prior_set = prior_spec$prior_set,
      metric = "cohort_mean_susceptible_probability",
      horizon_year = NA_integer_,
      summary_scalar(rowMeans(lin$pi_mat)),
      prior_degenerate_flag = degeneracy$degenerate_flag
    ),
    tibble(
      dataset = model_row$dataset,
      model_id = model_row$model_id,
      prior_set = prior_spec$prior_set,
      metric = "cohort_mean_cure_probability",
      horizon_year = NA_integer_,
      summary_scalar(rowMeans(1 - lin$pi_mat)),
      prior_degenerate_flag = degeneracy$degenerate_flag
    ),
    tibble(
      dataset = model_row$dataset,
      model_id = model_row$model_id,
      prior_set = prior_spec$prior_set,
      metric = "LNTF",
      horizon_year = NA_integer_,
      summary_scalar(rowMeans(1 - lin$pi_mat)),
      prior_degenerate_flag = degeneracy$degenerate_flag
    ),
    tibble(
      dataset = model_row$dataset,
      model_id = model_row$model_id,
      prior_set = prior_spec$prior_set,
      metric = "cohort_median_susceptible_time",
      horizon_year = NA_integer_,
      summary_scalar(apply(lin$median_mat, 1, stats::median)),
      prior_degenerate_flag = degeneracy$degenerate_flag
    )
  )

  risk_rows <- lapply(horizons_eval, function(h) {
    fh <- family_survival_hazard(h, model_row$family_code, lin$mu_lat_mat, lin$median_mat, sim_draws)
    risk_mat <- 1 - ((1 - lin$pi_mat) + lin$pi_mat * fh$Su)
    tibble(
      dataset = model_row$dataset,
      model_id = model_row$model_id,
      prior_set = prior_spec$prior_set,
      metric = "cohort_mean_risk",
      horizon_year = h,
      summary_scalar(rowMeans(risk_mat)),
      prior_degenerate_flag = degeneracy$degenerate_flag
    )
  })

  list(
    summary = bind_rows(prior_rows, risk_rows),
    degeneracy = degeneracy
  )
}

annotate_prior_predictive_summary <- function(df, mean_threshold_years, q975_threshold_years) {
  if (is.null(df) || !inherits(df, "data.frame")) {
    return(
      empty_prior_predictive_summary() %>%
        mutate(
          prior_tail_warning_flag = logical(),
          prior_tail_warning_detail = character()
        )
    )
  }

  normalize_prior_predictive_summary(df) %>%
    drop_existing_columns(c("prior_tail_warning_flag", "prior_tail_warning_detail")) %>%
    mutate(
      prior_tail_warning_flag = .data$metric == "cohort_median_susceptible_time" &
        ((safe_numeric(.data$mean) > mean_threshold_years) | (safe_numeric(.data$q975) > q975_threshold_years)),
      prior_tail_warning_detail = case_when(
        .data$metric == "cohort_median_susceptible_time" &
          safe_numeric(.data$mean) > mean_threshold_years &
          safe_numeric(.data$q975) > q975_threshold_years ~
          paste0(
            "Prior-predictive cohort median susceptible time is very long (mean > ",
            mean_threshold_years,
            " years and q975 > ",
            q975_threshold_years,
            " years)."
          ),
        .data$metric == "cohort_median_susceptible_time" &
          safe_numeric(.data$mean) > mean_threshold_years ~
          paste0(
            "Prior-predictive cohort median susceptible time mean exceeds ",
            mean_threshold_years,
            " years."
          ),
        .data$metric == "cohort_median_susceptible_time" &
          safe_numeric(.data$q975) > q975_threshold_years ~
          paste0(
            "Prior-predictive cohort median susceptible time q975 exceeds ",
            q975_threshold_years,
            " years."
          ),
        TRUE ~ NA_character_
      )
    )
}

## 🟠 Define: censoring-aware classification summaries ===============================
compute_classification_summary <- function(risk_draws, horizon_row, thresholds) {
  prevalence <- as.numeric(horizon_row$observed_km_risk)
  w_case <- unlist(horizon_row$w_case)
  w_control <- unlist(horizon_row$w_control)
  denom_case <- as.numeric(horizon_row$denom_case)
  denom_control <- as.numeric(horizon_row$denom_control)

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
    ppv_draw <- ifelse(pos_rate_draw > 0, prevalence * tpr_draw / pos_rate_draw, NA_real_)
    fdp_draw <- 1 - ppv_draw
    fp_burden_draw <- (1 - prevalence) * fpr_draw
    fp100_draw <- 100 * fp_burden_draw
    nb_draw <- prevalence * tpr_draw - (1 - prevalence) * fpr_draw * (thr / (1 - thr))

    out_abs[[j]] <- tibble(
      threshold = thr,
      pos_rate_mean = mean(pos_rate_draw, na.rm = TRUE),
      pos_rate_q025 = stats::quantile(pos_rate_draw, 0.025, na.rm = TRUE, names = FALSE),
      pos_rate_q50 = stats::quantile(pos_rate_draw, 0.500, na.rm = TRUE, names = FALSE),
      pos_rate_q975 = stats::quantile(pos_rate_draw, 0.975, na.rm = TRUE, names = FALSE),
      FPR_mean = mean(fpr_draw, na.rm = TRUE),
      FPR_q025 = stats::quantile(fpr_draw, 0.025, na.rm = TRUE, names = FALSE),
      FPR_q50 = stats::quantile(fpr_draw, 0.500, na.rm = TRUE, names = FALSE),
      FPR_q975 = stats::quantile(fpr_draw, 0.975, na.rm = TRUE, names = FALSE),
      false_positive_burden_mean = mean(fp_burden_draw, na.rm = TRUE),
      false_positive_burden_q025 = stats::quantile(fp_burden_draw, 0.025, na.rm = TRUE, names = FALSE),
      false_positive_burden_q50 = stats::quantile(fp_burden_draw, 0.500, na.rm = TRUE, names = FALSE),
      false_positive_burden_q975 = stats::quantile(fp_burden_draw, 0.975, na.rm = TRUE, names = FALSE),
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

## 🟠 Define: Stan compilation and PDF rendering ===============================
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
  beta_inc ~ normal(mu_beta_inc, sd_beta_inc);
  gamma0 ~ normal(0, sd_gamma0);
  gamma_lat ~ normal(rep_vector(0, K_lat), sd_gamma_lat);
  rho_W ~ normal(0, sd_shape);
  log_sigma_LN ~ normal(0, sd_shape);
  psi_LL ~ normal(0, sd_shape);

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
  rstan::stan_model(model_code = stan_code, model_name = "stage8_bayesian_cure_block")
}

select_trace_parameters <- function(family_code, K_inc, K_lat) {
  selected_pars <- c("delta0", "gamma0")
  if (K_inc >= 1) {
    selected_pars <- c(selected_pars, "beta_inc[1]")
  }
  if (family_code == "E") {
    if (K_lat >= 1) {
      selected_pars <- c(selected_pars, "gamma_lat[1]")
    }
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

  for (j in seq_len(min(length(selected_pars), 4L))) {
    matplot(
      arr[, , j],
      type = "l",
      lty = 1,
      col = seq_len(dim(arr)[2]),
      main = paste(model_id, selected_pars[[j]]),
      xlab = "Iteration",
      ylab = ""
    )
  }
  if (length(selected_pars) < 4L) {
    for (j in seq_len(4L - length(selected_pars))) {
      plot.new()
    }
  }
}

safe_generate_diagnostic_pdf <- function(trace_records, n_models_reused, n_models_to_fit, posterior_cohort_yearly, posterior_classification, ppc_summary, final_path) {
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

  if (n_models_reused > 0L && n_models_to_fit > 0L) {
    plot.new()
    text(0.5, 0.65, "Stage 8 diagnostic PDF note", cex = 1.25, font = 2)
    text(0.5, 0.48, paste0("Reused models from prior outputs: ", n_models_reused))
    text(0.5, 0.40, paste0("Newly fitted models in this run: ", n_models_to_fit))
    text(0.5, 0.28, "Trace pages below are shown only for models fitted in this run.")
    text(0.5, 0.20, "Aggregate pages at the end use the full final result tables.")
  }

  if (n_models_reused > 0L && n_models_to_fit == 0L) {
    plot.new()
    text(0.5, 0.62, "Stage 8 diagnostic PDF note", cex = 1.25, font = 2)
    text(0.5, 0.45, "All model-level fits were reused from existing Stage 8 outputs.")
    text(0.5, 0.32, "This PDF therefore contains aggregate pages only.")
  }

  if (length(trace_records) > 0L) {
    for (trace_record in trace_records) {
      plot_trace_record(trace_record)
    }
  }

  if (nrow(posterior_cohort_yearly) > 0) {
    g_risk <- posterior_cohort_yearly %>%
      ggplot(aes(x = horizon_year, y = meanRisk_Bayes_mean, color = model_id, fill = model_id)) +
      geom_ribbon(aes(ymin = meanRisk_Bayes_q025, ymax = meanRisk_Bayes_q975), alpha = 0.15, linewidth = 0) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~ dataset, scales = "free_y") +
      labs(
        title = "Stage 8 Bayesian cure: posterior cohort mean risk trajectories",
        x = "Horizon (years)",
        y = "Posterior mean risk"
      ) +
      theme_bw() +
      theme(legend.position = "none")
    print(g_risk)

    g_hazard <- posterior_cohort_yearly %>%
      ggplot(aes(x = horizon_year, y = meanHazard_mean, color = model_id, fill = model_id)) +
      geom_ribbon(aes(ymin = meanHazard_q025, ymax = meanHazard_q975), alpha = 0.15, linewidth = 0) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~ dataset, scales = "free_y") +
      labs(
        title = "Stage 8 Bayesian cure: posterior cohort mean hazard trajectories",
        x = "Horizon (years)",
        y = "Posterior mean hazard"
      ) +
      theme_bw() +
      theme(legend.position = "none")
    print(g_hazard)
  }

  if (nrow(posterior_classification) > 0) {
    nb_plot_df <- posterior_classification %>%
      filter(horizon_year %in% c(1, 2, 5))
    if (nrow(nb_plot_df) > 0) {
      g_nb <- nb_plot_df %>%
        ggplot(aes(x = threshold, y = NB_mean, color = model_id, fill = model_id)) +
        geom_ribbon(aes(ymin = NB_q025, ymax = NB_q975), alpha = 0.15, linewidth = 0) +
        geom_line(linewidth = 0.7) +
        facet_grid(dataset ~ horizon_year, scales = "free_y") +
        labs(
          title = "Stage 8 Bayesian cure: net benefit by threshold",
          x = "Risk threshold",
          y = "Net benefit"
        ) +
        theme_bw() +
        theme(legend.position = "none")
      print(g_nb)
    }
  }

  if (nrow(ppc_summary) > 0) {
    g_ppc <- ppc_summary %>%
      ggplot(aes(x = horizon_year, y = posterior_mean_risk, color = model_id)) +
      geom_errorbar(aes(ymin = posterior_q025_risk, ymax = posterior_q975_risk), width = 0.12, alpha = 0.6) +
      geom_line(linewidth = 0.6) +
      geom_point(linewidth = 0.6) +
      geom_point(aes(y = observed_km_risk), shape = 4, size = 2.0, stroke = 0.9, color = "black") +
      facet_wrap(~ dataset, scales = "free_y") +
      labs(
        title = "Stage 8 Bayesian cure: posterior predictive checks against observed KM risk",
        x = "Horizon (years)",
        y = "Risk"
      ) +
      theme_bw() +
      theme(legend.position = "none")
    print(g_ppc)
  }

  if (length(trace_records) == 0L && nrow(posterior_cohort_yearly) == 0L && nrow(posterior_classification) == 0L && nrow(ppc_summary) == 0L) {
    plot.new()
    text(0.5, 0.5, "No diagnostic pages were available for this Stage 8 run.")
  }

  close_pdf()
  if (!pdf_file_is_usable(tmp_pdf)) {
    stop("Temporary diagnostic PDF was not created correctly.", call. = FALSE)
  }
  safe_promote_file(tmp_pdf, final_path)
  invisible(TRUE)
}

## 🟠 Define: post-hoc Stage-10 comparison reconstruction ===============================
rebuild_delta_vs_nocure <- function(posterior_cohort_yearly, posterior_classification, nocure_cohort_long = NULL, nocure_class_long = NULL) {
  out_rows <- list()

  if (!is.null(nocure_cohort_long) && nrow(nocure_cohort_long) > 0 && nrow(posterior_cohort_yearly) > 0) {
    posterior_cohort_aug <- posterior_cohort_yearly %>%
      select(dataset, model_id, formula_anchor, horizon_year, meanRisk_Bayes_mean, meanRisk_Bayes_q025, meanRisk_Bayes_q50, meanRisk_Bayes_q975) %>%
      mutate(formula_anchor_original = formula_anchor) %>%
      bind_rows(
        posterior_cohort_yearly %>%
          select(dataset, model_id, formula_anchor, horizon_year, meanRisk_Bayes_mean, meanRisk_Bayes_q025, meanRisk_Bayes_q50, meanRisk_Bayes_q975) %>%
          mutate(formula_anchor_original = formula_anchor, formula_anchor = "ALL")
      )

    cohort_join <- posterior_cohort_aug %>%
      inner_join(
        nocure_cohort_long %>%
          filter(metric == "meanRisk") %>%
          select(dataset, no_cure_model_id, formula_anchor, horizon_year, metric, value),
        by = c("dataset", "formula_anchor", "horizon_year")
      ) %>%
      transmute(
        dataset = dataset,
        model_id = model_id,
        no_cure_model_id = no_cure_model_id,
        formula_anchor = formula_anchor_original,
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

  if (!is.null(nocure_class_long) && nrow(nocure_class_long) > 0 && nrow(posterior_classification) > 0) {
    metric_map <- tribble(
      ~metric, ~mean_col, ~q025_col, ~q50_col, ~q975_col,
      "FPR", "FPR_mean", "FPR_q025", "FPR_q50", "FPR_q975",
      "FP100", "FP100_mean", "FP100_q025", "FP100_q50", "FP100_q975",
      "NB", "NB_mean", "NB_q025", "NB_q50", "NB_q975",
      "false_positive_burden", "false_positive_burden_mean", "false_positive_burden_q025", "false_positive_burden_q50", "false_positive_burden_q975"
    )

    class_long <- bind_rows(lapply(seq_len(nrow(metric_map)), function(i) {
      one <- metric_map[i, , drop = FALSE]
      posterior_classification %>%
        transmute(
          dataset = dataset,
          model_id = model_id,
          formula_anchor = formula_anchor,
          horizon_year = horizon_year,
          threshold = threshold,
          metric = one$metric[[1]],
          bayes_mean = .data[[one$mean_col[[1]]]],
          bayes_q025 = .data[[one$q025_col[[1]]]],
          bayes_q50 = .data[[one$q50_col[[1]]]],
          bayes_q975 = .data[[one$q975_col[[1]]]]
        )
    }))

    class_long_aug <- class_long %>%
      mutate(formula_anchor_original = formula_anchor) %>%
      bind_rows(
        class_long %>%
          mutate(formula_anchor_original = formula_anchor, formula_anchor = "ALL")
      )

    class_join <- class_long_aug %>%
      inner_join(
        nocure_class_long %>%
          select(dataset, no_cure_model_id, formula_anchor, horizon_year, threshold, metric, value),
        by = c("dataset", "formula_anchor", "horizon_year", "threshold", "metric")
      ) %>%
      transmute(
        dataset = dataset,
        model_id = model_id,
        no_cure_model_id = no_cure_model_id,
        formula_anchor = formula_anchor_original,
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
  if (nrow(out) == 0) {
    return(empty_delta_template())
  }

  out %>%
    arrange(dataset, model_id, horizon_year, threshold, metric, no_cure_model_id)
}

stage8_search_scalar_in_object <- function(obj, candidate_names, max_depth = 3L) {
  if (max_depth < 0L || is.null(obj)) {
    return(NULL)
  }

  if (is.data.frame(obj)) {
    hit <- candidate_names[candidate_names %in% names(obj)]
    if (length(hit) > 0L && nrow(obj) >= 1L) {
      return(obj[[hit[[1L]]]][[1L]])
    }
    return(NULL)
  }

  if (is.list(obj)) {
    nms <- names(obj)
    if (!is.null(nms)) {
      hit <- candidate_names[candidate_names %in% nms]
      if (length(hit) > 0L) {
        val <- obj[[hit[[1L]]]]
        if (is.data.frame(val) && nrow(val) >= 1L) {
          return(val[[1L]][[1L]])
        }
        if (is.atomic(val) && length(val) >= 1L && is.null(dim(val))) {
          return(val[[1L]])
        }
      }

      priority_names <- c(
        "fit_qc", "qc", "diagnostics", "diag", "fit_summary", "summary",
        "information_criteria", "info_criteria", "waic", "loo", "waic_result",
        "loo_result", "metrics", "fit_metrics", "diagnostic_summary",
        "export_cache", "export_payload", "registry_row", "model_registry_row",
        "metadata"
      )

      for (nm in priority_names[priority_names %in% nms]) {
        val <- stage8_search_scalar_in_object(obj[[nm]], candidate_names, max_depth = max_depth - 1L)
        if (!is.null(val)) {
          return(val)
        }
      }
    }
  }

  NULL
}

stage8_scalar_numeric <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(as.character(x[[1L]])))
}

stage8_scalar_character <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(NA_character_)
  }
  as.character(x[[1L]])
}

stage8_scalar_logical <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(NA)
  }
  if (is.logical(x[[1L]])) {
    return(as.logical(x[[1L]]))
  }
  stage8_parse_logicalish(x[[1L]])[[1L]]
}

stage8_read_fit_qc_from_rds <- function(rds_path) {
  blank <- tibble(
    divergences = NA_real_,
    max_rhat = NA_real_,
    min_bulk_ess = NA_real_,
    min_tail_ess = NA_real_,
    treedepth_exceeded = NA_real_,
    waic = NA_real_,
    looic = NA_real_,
    p_waic = NA_real_,
    p_waic_high_n = NA_real_,
    p_waic_high_pct = NA_real_,
    p_loo = NA_real_,
    pareto_k_max = NA_real_,
    pareto_k_bad_n = NA_real_,
    pareto_k_bad_pct = NA_real_,
    pareto_k_very_bad_n = NA_real_,
    waic_warning_flag = NA,
    loo_warning_flag = NA,
    info_criteria_warning_detail = NA_character_
  )

  if (is.na(rds_path) || !nzchar(as.character(rds_path)) || !file.exists(rds_path)) {
    return(blank)
  }

  obj <- tryCatch(readRDS(rds_path), error = function(e) e)
  if (inherits(obj, "error")) {
    return(blank)
  }

  get_num <- function(candidates) {
    stage8_scalar_numeric(stage8_search_scalar_in_object(obj, candidates, max_depth = 4L))
  }
  get_chr <- function(candidates) {
    stage8_scalar_character(stage8_search_scalar_in_object(obj, candidates, max_depth = 4L))
  }
  get_lgl <- function(candidates) {
    stage8_scalar_logical(stage8_search_scalar_in_object(obj, candidates, max_depth = 4L))
  }

  waic <- get_num(c("waic", "waic_ic", "waic_value"))
  if (is.na(waic)) {
    elpd_waic <- get_num(c("elpd_waic"))
    if (!is.na(elpd_waic)) {
      waic <- -2 * elpd_waic
    }
  }

  looic <- get_num(c("looic", "loo_ic", "loo_value"))
  if (is.na(looic)) {
    elpd_loo <- get_num(c("elpd_loo"))
    if (!is.na(elpd_loo)) {
      looic <- -2 * elpd_loo
    }
  }

  out <- tibble(
    divergences = get_num(c("divergences", "n_divergences", "n_divergent", "num_divergent", "n_divergent_transitions")),
    max_rhat = get_num(c("max_rhat", "rhat_max")),
    min_bulk_ess = get_num(c("min_bulk_ess", "bulk_ess_min", "ess_bulk_min")),
    min_tail_ess = get_num(c("min_tail_ess", "tail_ess_min", "ess_tail_min")),
    treedepth_exceeded = get_num(c("treedepth_exceeded", "treedepth_exceeded_n", "n_treedepth_exceeded")),
    waic = waic,
    looic = looic,
    p_waic = get_num(c("p_waic")),
    p_waic_high_n = get_num(c("p_waic_high_n")),
    p_waic_high_pct = get_num(c("p_waic_high_pct")),
    p_loo = get_num(c("p_loo")),
    pareto_k_max = get_num(c("pareto_k_max", "k_max")),
    pareto_k_bad_n = get_num(c("pareto_k_bad_n", "k_bad_n", "pareto_bad_n")),
    pareto_k_bad_pct = get_num(c("pareto_k_bad_pct", "k_bad_pct", "pareto_bad_pct")),
    pareto_k_very_bad_n = get_num(c("pareto_k_very_bad_n", "k_very_bad_n", "pareto_very_bad_n")),
    waic_warning_flag = get_lgl(c("waic_warning_flag", "waic_warning")),
    loo_warning_flag = get_lgl(c("loo_warning_flag", "loo_warning")),
    info_criteria_warning_detail = get_chr(c("info_criteria_warning_detail", "ic_warning_detail", "warning_detail"))
  )

  if (is.na(out$waic_warning_flag) && !is.na(out$p_waic_high_n)) {
    out$waic_warning_flag <- out$p_waic_high_n > 0
  }
  if (is.na(out$loo_warning_flag) && (!is.na(out$pareto_k_bad_n) || !is.na(out$pareto_k_very_bad_n))) {
    out$loo_warning_flag <- dplyr::coalesce(out$pareto_k_bad_n, 0) > 0 | dplyr::coalesce(out$pareto_k_very_bad_n, 0) > 0
  }

  out
}

stage8_ensure_model_registry_columns <- function(model_registry) {
  out <- tibble::as_tibble(model_registry)

  required_cols <- c(
    "fit_status", "rds_path", "reuse_existing_fit", "fit_reused_flag",
    "divergences", "max_rhat", "min_bulk_ess", "min_tail_ess", "treedepth_exceeded",
    "waic", "looic", "p_waic", "p_waic_high_n", "p_waic_high_pct", "p_loo",
    "pareto_k_max", "pareto_k_bad_n", "pareto_k_bad_pct", "pareto_k_very_bad_n",
    "waic_warning_flag", "loo_warning_flag", "info_criteria_warning_detail"
  )

  for (nm in setdiff(required_cols, names(out))) {
    if (nm %in% c("fit_status", "rds_path", "info_criteria_warning_detail")) {
      out[[nm]] <- NA_character_
    } else if (nm %in% c("waic_warning_flag", "loo_warning_flag", "reuse_existing_fit", "fit_reused_flag")) {
      out[[nm]] <- NA
    } else {
      out[[nm]] <- NA_real_
    }
  }

  out
}

fill_stage8_model_registry_reuse_qc <- function(model_registry) {
  out <- stage8_ensure_model_registry_columns(model_registry)

  out$fit_reused_flag <- dplyr::coalesce(
    as.logical(out$fit_reused_flag),
    as.logical(out$reuse_existing_fit),
    FALSE
  )

  core_missing_flag <- is.na(out$divergences) | is.na(out$waic) | is.na(out$looic)
  needs_harvest <- out$fit_status == "ok" &
    out$fit_reused_flag &
    core_missing_flag &
    !is.na(out$rds_path) &
    nzchar(as.character(out$rds_path))

  if (!any(needs_harvest)) {
    return(out)
  }

  harvested <- bind_rows(lapply(out$rds_path[needs_harvest], stage8_read_fit_qc_from_rds))
  target_idx <- which(needs_harvest)
  overlapping_cols <- intersect(names(harvested), names(out))

  for (nm in overlapping_cols) {
    out[[nm]][target_idx] <- dplyr::coalesce(out[[nm]][target_idx], harvested[[nm]])
  }

  if ("waic_warning_flag" %in% names(out)) {
    out$waic_warning_flag <- dplyr::coalesce(
      as.logical(out$waic_warning_flag),
      (!is.na(out$p_waic_high_n) & out$p_waic_high_n > 0)
    )
  }

  if ("loo_warning_flag" %in% names(out)) {
    out$loo_warning_flag <- dplyr::coalesce(
      as.logical(out$loo_warning_flag),
      (!is.na(out$pareto_k_bad_n) & out$pareto_k_bad_n > 0) |
        (!is.na(out$pareto_k_very_bad_n) & out$pareto_k_very_bad_n > 0)
    )
  }

  out
}

## 🟠 Define: self-audit checks for downstream safety ===============================
make_stage8_output_audit <- function(model_grid, model_registry, posterior_cohort_yearly, posterior_classification, ppc_summary, posterior_delta_vs_nocure, diagnostic_pdf_path, horizons_year, risk_thresholds) {
  model_registry <- fill_stage8_model_registry_reuse_qc(simplify_scalar_list_cols(model_registry))
  posterior_cohort_yearly <- simplify_scalar_list_cols(posterior_cohort_yearly)
  posterior_classification <- simplify_scalar_list_cols(posterior_classification)
  ppc_summary <- simplify_scalar_list_cols(ppc_summary)
  posterior_delta_vs_nocure <- simplify_scalar_list_cols(posterior_delta_vs_nocure)

  if ("model_id" %in% names(model_registry)) {
    model_registry$model_id <- as.character(model_registry$model_id)
  }
  if (all(c("model_id", "horizon_year") %in% names(posterior_cohort_yearly))) {
    posterior_cohort_yearly$model_id <- as.character(posterior_cohort_yearly$model_id)
    posterior_cohort_yearly$horizon_year <- safe_numeric(posterior_cohort_yearly$horizon_year)
  }
  if (all(c("model_id", "horizon_year", "threshold") %in% names(posterior_classification))) {
    posterior_classification$model_id <- as.character(posterior_classification$model_id)
    posterior_classification$horizon_year <- safe_numeric(posterior_classification$horizon_year)
    posterior_classification$threshold <- safe_numeric(posterior_classification$threshold)
  }
  if (all(c("model_id", "horizon_year") %in% names(ppc_summary))) {
    ppc_summary$model_id <- as.character(ppc_summary$model_id)
    ppc_summary$horizon_year <- safe_numeric(ppc_summary$horizon_year)
  }

  ok_fit_flag <- model_registry$fit_status == "ok"
  fit_reused_flag <- dplyr::coalesce(as.logical(model_registry$fit_reused_flag), as.logical(model_registry$reuse_existing_fit), FALSE)
  core_qc_missing_flag <- is.na(model_registry$divergences) | is.na(model_registry$waic) | is.na(model_registry$looic)
  core_qc_required_flag <- ok_fit_flag & !fit_reused_flag
  missing_required_core_qc_n <- sum(core_qc_required_flag & core_qc_missing_flag, na.rm = TRUE)
  missing_reused_core_qc_n <- sum(ok_fit_flag & fit_reused_flag & core_qc_missing_flag, na.rm = TRUE)

  admissible_n <- sum(as.logical(model_registry$admissible_flag), na.rm = TRUE)

  registry_dup_n <- if ("model_id" %in% names(model_registry) && nrow(model_registry) > 0L) {
    sum(duplicated(model_registry[c("model_id")]))
  } else {
    0L
  }

  cohort_dup_n <- if (all(c("model_id", "horizon_year") %in% names(posterior_cohort_yearly)) && nrow(posterior_cohort_yearly) > 0L) {
    sum(duplicated(posterior_cohort_yearly[c("model_id", "horizon_year")]))
  } else {
    0L
  }

  class_dup_n <- if (all(c("model_id", "horizon_year", "threshold") %in% names(posterior_classification)) && nrow(posterior_classification) > 0L) {
    sum(duplicated(posterior_classification[c("model_id", "horizon_year", "threshold")]))
  } else {
    0L
  }

  ppc_dup_n <- if (all(c("model_id", "horizon_year") %in% names(ppc_summary)) && nrow(ppc_summary) > 0L) {
    sum(duplicated(ppc_summary[c("model_id", "horizon_year")]))
  } else {
    0L
  }

  cohort_missing_meta_n <- if (nrow(posterior_cohort_yearly) == 0L) {
    0L
  } else {
    sum(
      is.na(posterior_cohort_yearly$screening_flag) |
        is.na(posterior_cohort_yearly$admissible_flag) |
        is.na(posterior_cohort_yearly$support_priority)
    )
  }

  class_missing_meta_n <- if (nrow(posterior_classification) == 0L) {
    0L
  } else {
    sum(
      is.na(posterior_classification$classification_reporting_status) |
        is.na(posterior_classification$support_priority)
    )
  }

  delta_missing_meta_n <- if (nrow(posterior_delta_vs_nocure) == 0L) {
    0L
  } else {
    sum(
      is.na(posterior_delta_vs_nocure$horizon_support_label) |
        is.na(posterior_delta_vs_nocure$support_priority)
    )
  }

  bind_rows(
    tibble(
      check_name = "model_registry_row_count",
      status = ifelse(nrow(model_registry) == nrow(model_grid), "pass", "fail"),
      observed_value = as.character(nrow(model_registry)),
      expected_value = as.character(nrow(model_grid)),
      detail = "One registry row should exist per model in the Stage 8 grid."
    ),
    tibble(
      check_name = "model_registry_duplicate_model_id",
      status = ifelse(registry_dup_n == 0L, "pass", "fail"),
      observed_value = as.character(registry_dup_n),
      expected_value = "0",
      detail = "Model IDs must be unique in the model registry."
    ),
    tibble(
      check_name = "model_registry_core_qc_fields_nonmissing",
      status = ifelse(missing_required_core_qc_n == 0L, "pass", "fail"),
      observed_value = as.character(missing_required_core_qc_n),
      expected_value = "0",
      detail = paste0(
        "Newly fit successful rows must have populated core QC summary fields. Reused successful rows still missing core QC after harvest: ",
        missing_reused_core_qc_n,
        " (allowed for legacy RDS without stored IC summaries)."
      )
    ),
    tibble(
      check_name = "posterior_cohort_yearly_row_count",
      status = ifelse(nrow(posterior_cohort_yearly) == admissible_n * length(horizons_year), "pass", "fail"),
      observed_value = as.character(nrow(posterior_cohort_yearly)),
      expected_value = as.character(admissible_n * length(horizons_year)),
      detail = "Admissible models should contribute one cohort-level row per horizon."
    ),
    tibble(
      check_name = "posterior_classification_row_count",
      status = ifelse(nrow(posterior_classification) == admissible_n * length(horizons_year) * length(risk_thresholds), "pass", "fail"),
      observed_value = as.character(nrow(posterior_classification)),
      expected_value = as.character(admissible_n * length(horizons_year) * length(risk_thresholds)),
      detail = "Admissible models should contribute one threshold row per horizon and threshold."
    ),
    tibble(
      check_name = "ppc_summary_row_count",
      status = ifelse(nrow(ppc_summary) == nrow(model_grid) * length(horizons_year), "pass", "fail"),
      observed_value = as.character(nrow(ppc_summary)),
      expected_value = as.character(nrow(model_grid) * length(horizons_year)),
      detail = "Every model should contribute one PPC row per horizon."
    ),
    tibble(
      check_name = "posterior_cohort_yearly_duplicate_keys",
      status = ifelse(cohort_dup_n == 0L, "pass", "fail"),
      observed_value = as.character(cohort_dup_n),
      expected_value = "0",
      detail = "Each model-horizon key should appear at most once in the cohort-level table."
    ),
    tibble(
      check_name = "posterior_classification_duplicate_keys",
      status = ifelse(class_dup_n == 0L, "pass", "fail"),
      observed_value = as.character(class_dup_n),
      expected_value = "0",
      detail = "Each model-horizon-threshold key should appear at most once in the classification table."
    ),
    tibble(
      check_name = "ppc_summary_duplicate_keys",
      status = ifelse(ppc_dup_n == 0L, "pass", "fail"),
      observed_value = as.character(ppc_dup_n),
      expected_value = "0",
      detail = "Each model-horizon key should appear at most once in the PPC table."
    ),
    tibble(
      check_name = "posterior_cohort_yearly_missing_carryforward_metadata",
      status = ifelse(cohort_missing_meta_n == 0L, "pass", "fail"),
      observed_value = as.character(cohort_missing_meta_n),
      expected_value = "0",
      detail = "Cohort-level output should retain screening, admissibility, and support metadata."
    ),
    tibble(
      check_name = "posterior_classification_missing_reporting_metadata",
      status = ifelse(class_missing_meta_n == 0L, "pass", "fail"),
      observed_value = as.character(class_missing_meta_n),
      expected_value = "0",
      detail = "Threshold-based output should retain estimability and reporting metadata."
    ),
    tibble(
      check_name = "posterior_delta_vs_nocure_missing_support_metadata",
      status = ifelse(delta_missing_meta_n == 0L, "pass", "fail"),
      observed_value = as.character(delta_missing_meta_n),
      expected_value = "0",
      detail = "Bayesian-vs-no-cure comparison rows should retain horizon support metadata."
    ),
    tibble(
      check_name = "diagnostic_pdf_exists",
      status = ifelse(pdf_file_is_usable(diagnostic_pdf_path), "pass", "fail"),
      observed_value = ifelse(pdf_file_is_usable(diagnostic_pdf_path), "TRUE", "FALSE"),
      expected_value = "TRUE",
      detail = "Diagnostic PDF should exist after the run, unless PDF generation was intentionally skipped."
    )
  )
}

stage8_finalize_registry_and_audit <- function(model_grid,
                                               model_registry,
                                               posterior_cohort_yearly,
                                               posterior_classification,
                                               ppc_summary,
                                               posterior_delta_vs_nocure,
                                               diagnostic_pdf_path,
                                               horizons_year,
                                               risk_thresholds) {
  model_registry_fixed <- fill_stage8_model_registry_reuse_qc(model_registry)

  output_audit_fixed <- make_stage8_output_audit(
    model_grid = model_grid,
    model_registry = model_registry_fixed,
    posterior_cohort_yearly = posterior_cohort_yearly,
    posterior_classification = posterior_classification,
    ppc_summary = ppc_summary,
    posterior_delta_vs_nocure = posterior_delta_vs_nocure,
    diagnostic_pdf_path = diagnostic_pdf_path,
    horizons_year = horizons_year,
    risk_thresholds = risk_thresholds
  )

  list(
    model_registry = model_registry_fixed,
    output_audit = output_audit_fixed
  )
}

stage8_repair_existing_stage8_exports <- function(export_path,
                                                  diagnostic_plot_pdf_name = "bayes_stage8_diagnostic_plots.pdf") {
  model_registry_path <- file.path(export_path, "bayes_stage8_model_registry.csv")
  posterior_cohort_yearly_path <- file.path(export_path, "bayes_stage8_posterior_cohort_yearly.csv")
  posterior_classification_path <- file.path(export_path, "bayes_stage8_posterior_classification.csv")
  ppc_summary_path <- file.path(export_path, "bayes_stage8_ppc_summary.csv")
  posterior_delta_vs_nocure_path <- file.path(export_path, "bayes_stage8_posterior_delta_vs_nocure.csv")
  output_audit_path <- file.path(export_path, "bayes_stage8_output_audit.csv")
  diagnostic_plot_pdf_path <- file.path(export_path, diagnostic_plot_pdf_name)

  model_registry <- read_delimited_or_rds(model_registry_path)
  posterior_cohort_yearly <- read_delimited_or_rds(posterior_cohort_yearly_path)
  posterior_classification <- read_delimited_or_rds(posterior_classification_path)
  ppc_summary <- read_delimited_or_rds(ppc_summary_path)
  posterior_delta_vs_nocure <- if (file.exists(posterior_delta_vs_nocure_path)) {
    read_delimited_or_rds(posterior_delta_vs_nocure_path)
  } else {
    empty_delta_template()
  }

  horizons_year_now <- if (nrow_or_zero(posterior_cohort_yearly) > 0L) {
    sort(unique(as.integer(safe_numeric(posterior_cohort_yearly$horizon_year))))
  } else {
    integer()
  }
  risk_thresholds_now <- if (nrow_or_zero(posterior_classification) > 0L) {
    sort(unique(safe_numeric(posterior_classification$threshold)))
  } else {
    numeric()
  }

  fixed <- stage8_finalize_registry_and_audit(
    model_grid = tibble(model_id = unique(model_registry$model_id)),
    model_registry = model_registry,
    posterior_cohort_yearly = posterior_cohort_yearly,
    posterior_classification = posterior_classification,
    ppc_summary = ppc_summary,
    posterior_delta_vs_nocure = posterior_delta_vs_nocure,
    diagnostic_pdf_path = diagnostic_plot_pdf_path,
    horizons_year = horizons_year_now,
    risk_thresholds = risk_thresholds_now
  )

  write_csv_preserve_schema(fixed$model_registry, model_registry_path)
  write_csv_preserve_schema(fixed$output_audit, output_audit_path)
  invisible(fixed)
}

# 🔴 Load: analysis-ready datasets and optional comparators ===============================
## 🟠 Load: backbone data objects and Stage-5/6 linkage ===============================
loaded_objects <- load_stage8_datasets(
  merged_path = merged_data_path,
  pnu_label = pnu_site_label,
  snu_label = snu_site_label
)

analysis_datasets <- loaded_objects$datasets
scaling_registry <- loaded_objects$scaling

screening_flags <- read_screening_flags(stage6_screening_flag_csv)
nocure_cohort_long <- normalize_nocure_cohort(stage5_nocure_cohort_csv)
nocure_class_long <- normalize_nocure_classification(stage5_nocure_classification_csv)

if (nrow(screening_flags) == 0L) {
  warning("Stage 6 screening linkage was not recovered. `screening_flag` and `screening_detail` may remain missing.", call. = FALSE)
}
if (is.null(nocure_cohort_long) || nrow(nocure_cohort_long) == 0L) {
  warning("Stage 5 cohort comparison linkage was not recovered. `bayes_stage8_posterior_delta_vs_nocure.csv` may be partially or fully empty.", call. = FALSE)
}
if (is.null(nocure_class_long) || nrow(nocure_class_long) == 0L) {
  warning("Stage 5 threshold comparison linkage was not recovered. `bayes_stage8_posterior_delta_vs_nocure.csv` may be partially or fully empty.", call. = FALSE)
}

dataset_registry <- bind_rows(
  lapply(names(analysis_datasets), function(ds) {
    dat <- analysis_datasets[[ds]]
    tibble(
      dataset = ds,
      n = nrow(dat),
      n_event = sum(dat$event_main),
      n_censor_main = sum(dat$censor_main),
      n_remission = sum(dat$remission_flag),
      mean_age_exact_entry = mean(dat$age_exact_entry),
      prop_male = mean(dat$sex_num == 1L)
    )
  })
)

## 🟠 Load: IPCW benchmarks and horizon metadata ===============================
ipcw_registry <- lapply(names(analysis_datasets), function(ds) {
  dat <- analysis_datasets[[ds]]
  ref <- build_ipcw_reference(dat, horizons_year)
  ref$dataset <- ds
  ref
})
names(ipcw_registry) <- names(analysis_datasets)

horizon_metadata_registry <- make_horizon_metadata_registry(ipcw_registry)

observed_km_risk <- horizon_metadata_registry %>%
  select(
    dataset,
    horizon_year,
    observed_km_risk,
    horizon_support_label,
    support_priority,
    cohort_reporting_status,
    cohort_reporting_note
  )

# 🔴 Prepare: model grid, export reuse plan, and Stan availability ===============================
prior_specs <- build_prior_specs()
model_grid <- build_model_grid(include_merged_incidence_site_supplementary)

if (!is.null(run_model_ids)) {
  model_grid <- model_grid %>% filter(model_id %in% run_model_ids)
  if (nrow(model_grid) == 0) {
    stop("`run_model_ids` filtered out all models.", call. = FALSE)
  }
}

screening_lookup <- build_screening_model_lookup(screening_flags, model_grid)

existing_exports <- if (isTRUE(reuse_existing_stage8_outputs)) {
  load_existing_stage8_exports(export_path)
} else {
  list(
    model_registry = NULL,
    coefficient_summary = NULL,
    posterior_subject_profile = NULL,
    posterior_subject_yearly = NULL,
    posterior_cohort_yearly = NULL,
    posterior_classification = NULL,
    diagnostics_parameter_level = NULL,
    ppc_summary = NULL,
    prior_predictive_summary = NULL,
    diagnostic_pdf_exists = FALSE
  )
}

reuse_flags <- vapply(seq_len(nrow(model_grid)), function(ii) {
  model_row <- model_grid[ii, , drop = FALSE]
  dataset_df <- analysis_datasets[[model_row$dataset[[1]]]]
  is_model_reusable(
    model_row = model_row,
    dataset_df = dataset_df,
    existing_exports = existing_exports,
    export_dir = export_path,
    require_existing_rds = require_existing_rds_to_skip_fit
  )
}, logical(1))
model_grid$reuse_existing_fit <- reuse_flags

n_models_reused <- sum(model_grid$reuse_existing_fit)
n_models_to_fit <- sum(!model_grid$reuse_existing_fit)

message(
  "Stage 8 reuse plan: ",
  n_models_reused,
  " model(s) reused from existing outputs/RDS; ",
  n_models_to_fit,
  " model(s) require fitting."
)

need_fit_packages <- n_models_to_fit > 0L
if (need_fit_packages) {
  fit_packages <- c("rstan", "posterior", "loo")
  missing_fit_packages <- fit_packages[
    !vapply(fit_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_fit_packages) > 0) {
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
  stan_model_compiled <- compile_stage8_stan_model()
} else {
  stan_model_compiled <- NULL
}

diagnostic_pdf_path <- file.path(export_path, "bayes_stage8_diagnostic_plots.pdf")
create_diagnostic_pdf <- (!isTRUE(reuse_existing_stage8_outputs)) ||
  (n_models_to_fit > 0L) ||
  (!isTRUE(preserve_existing_diagnostic_pdf)) ||
  (!pdf_file_is_usable(diagnostic_pdf_path))

trace_records <- list()

# 🔴 Run: model grid with reuse-first execution ===============================
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
  model_id_now <- model_row$model_id[[1]]
  dataset_name <- model_row$dataset[[1]]
  dataset_df <- analysis_datasets[[dataset_name]]
  model_started_at <- Sys.time()

  emit_stage8_progress(
    ii - 1L,
    total_models,
    model_id_now,
    paste0(
      "starting model fit (dataset=",
      dataset_name,
      ", family=",
      model_row$family_code[[1]],
      ", reuse=",
      isTRUE(model_row$reuse_existing_fit[[1]]),
      ")"
    )
  )

  screening_row <- screening_lookup %>%
    filter(model_id == model_id_now)
  if (nrow(screening_row) == 0L) {
    screening_row <- empty_screening_lookup_row(model_id_now)
  }
  screening_flag_now <- screening_value_or_na(screening_row, "screening_flag")
  screening_detail_now <- screening_value_or_na(screening_row, "screening_detail")
  cure_model_eligibility_flag_now <- screening_value_or_na(screening_row, "cure_model_eligibility_flag")
  primary_gate_method_now <- screening_value_or_na(screening_row, "primary_gate_method")
  primary_gate_flag_now <- screening_value_or_na(screening_row, "primary_gate_flag")
  receus_primary_class_now <- screening_value_or_na(screening_row, "receus_primary_class")
  presence_modifier_flag_now <- screening_value_or_na(screening_row, "presence_modifier_flag")
  cure_presence_support_flag_now <- screening_value_or_na(screening_row, "cure_presence_support_flag")
  followup_contradiction_flag_now <- screening_value_or_na(screening_row, "followup_contradiction_flag")
  followup_not_contradicted_flag_now <- screening_value_or_na(screening_row, "followup_not_contradicted_flag")
  screening_note_now <- screening_value_or_na(screening_row, "screening_note")
  carry_forward_stage8_now <- screening_logical_or_na(screening_row, "carry_forward_stage8")

  existing_reg <- subset_model_registry_row(existing_exports$model_registry, model_id_now)
  reuse_bundle <- if (isTRUE(model_row$reuse_existing_fit[[1]])) {
    get_reuse_bundle(model_id_now, export_path, existing_reg)
  } else {
    empty_reuse_bundle()
  }

  prior_main <- prior_specs$main
  prior_wider <- prior_specs$wider_latency

  existing_prior_model <- bind_rows_safe(list(
    subset_model_table(existing_exports$prior_predictive_summary, model_id_now),
    subset_model_table(reuse_bundle$prior_predictive_summary, model_id_now)
  )) %>%
    select(-any_of(c("prior_tail_warning_flag", "prior_tail_warning_detail"))) %>%
    normalize_prior_predictive_summary()

  need_prior_regeneration <- !(isTRUE(model_row$reuse_existing_fit[[1]]) && nrow(existing_prior_model) > 0L)

  design_main <- NULL
  design_wider <- NULL
  prior_pred_main <- NULL
  prior_pred_wider <- NULL

  if (need_prior_regeneration || !isTRUE(model_row$reuse_existing_fit[[1]])) {
    design_main <- make_design_bundle(dataset_df, model_row, prior_main, snu_site_label)
    design_wider <- make_design_bundle(dataset_df, model_row, prior_wider, snu_site_label)
  }

  if (need_prior_regeneration) {
    set.seed(stan_seed + ii)
    prior_pred_main <- simulate_prior_predictive(
      df = dataset_df,
      design_bundle = design_main,
      model_row = model_row,
      prior_spec = prior_main,
      n_draws = prior_predictive_draws,
      horizons_eval = prior_predictive_horizons
    )

    set.seed(stan_seed + 100000L + ii)
    prior_pred_wider <- simulate_prior_predictive(
      df = dataset_df,
      design_bundle = design_wider,
      model_row = model_row,
      prior_spec = prior_wider,
      n_draws = prior_predictive_draws,
      horizons_eval = prior_predictive_horizons
    )

    prior_predictive_rows[[length(prior_predictive_rows) + 1L]] <- prior_pred_main$summary
    prior_predictive_rows[[length(prior_predictive_rows) + 1L]] <- prior_pred_wider$summary
  } else {
    prior_predictive_rows[[length(prior_predictive_rows) + 1L]] <- existing_prior_model
  }

  if (isTRUE(model_row$reuse_existing_fit[[1]])) {
    if (!is.null(existing_reg)) {
      main_registry_rows[[length(main_registry_rows) + 1L]] <- existing_reg
    }

    coef_rows[[length(coef_rows) + 1L]] <- pick_reuse_table(
      model_id = model_id_now,
      existing_df = existing_exports$coefficient_summary,
      bundle_df = reuse_bundle$coefficient_summary
    )

    diagnostic_rows[[length(diagnostic_rows) + 1L]] <- pick_reuse_table(
      model_id = model_id_now,
      existing_df = existing_exports$diagnostics_parameter_level,
      bundle_df = reuse_bundle$diagnostics_parameter_level
    )

    ppc_rows_all[[length(ppc_rows_all) + 1L]] <- pick_reuse_table(
      model_id = model_id_now,
      existing_df = existing_exports$ppc_summary,
      bundle_df = reuse_bundle$ppc_summary
    )

    existing_reg_admissible <- !is.null(existing_reg) && isTRUE(as.logical(existing_reg$admissible_flag[[1]]))

    if (isTRUE(existing_reg_admissible)) {
      subject_profile_rows[[length(subject_profile_rows) + 1L]] <- pick_reuse_table(
        model_id = model_id_now,
        existing_df = existing_exports$posterior_subject_profile,
        bundle_df = reuse_bundle$posterior_subject_profile
      )

      subject_yearly_rows[[length(subject_yearly_rows) + 1L]] <- pick_reuse_table(
        model_id = model_id_now,
        existing_df = existing_exports$posterior_subject_yearly,
        bundle_df = reuse_bundle$posterior_subject_yearly
      )

      cohort_yearly_rows[[length(cohort_yearly_rows) + 1L]] <- pick_reuse_table(
        model_id = model_id_now,
        existing_df = existing_exports$posterior_cohort_yearly,
        bundle_df = reuse_bundle$posterior_cohort_yearly
      )

      class_rows[[length(class_rows) + 1L]] <- pick_reuse_table(
        model_id = model_id_now,
        existing_df = existing_exports$posterior_classification,
        bundle_df = reuse_bundle$posterior_classification
      )
    }

    emit_stage8_progress(
      ii,
      total_models,
      model_id_now,
      paste0("reused existing fit; elapsed=", format_stage8_number(elapsed_stage8_seconds(model_started_at), digits = 1L), "s")
    )

    next
  }

  if (is.null(design_main) || is.null(design_wider)) {
    design_main <- make_design_bundle(dataset_df, model_row, prior_main, snu_site_label)
    design_wider <- make_design_bundle(dataset_df, model_row, prior_wider, snu_site_label)
  }

  stan_data <- list(
    N = nrow(dataset_df),
    time = as.numeric(design_main$time),
    event = as.integer(design_main$event),
    K_inc = ncol(design_main$X_inc),
    X_inc = design_main$X_inc,
    K_lat = ncol(design_main$X_lat),
    X_lat = design_main$X_lat,
    family_id = as.integer(model_row$family_id[[1]]),
    alpha_gp = as.numeric(design_main$alpha_gp),
    mu_beta_inc = as.numeric(design_main$mu_beta_inc),
    sd_beta_inc = as.numeric(design_main$sd_beta_inc),
    sd_delta = as.numeric(prior_main$sd_delta),
    sd_gamma0 = as.numeric(prior_main$sd_gamma0),
    sd_gamma_lat = as.numeric(design_main$sd_gamma_lat),
    sd_shape = as.numeric(prior_main$sd_shape)
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
      model_id = model_id_now,
      structural_model_id = model_row$structural_model_id[[1]],
      latency_family = model_row$latency_family[[1]],
      family_code = model_row$family_code[[1]],
      formula_anchor = model_row$formula_anchor[[1]],
      incidence_site_indicator = model_row$incidence_site_indicator[[1]],
      latency_site_indicator = model_row$latency_site_indicator[[1]],
      latency_interaction_indicator = model_row$latency_interaction_indicator[[1]],
      is_supplementary_branch = model_row$is_supplementary_branch[[1]],
      fit_status = fit_status,
      fit_error_message = fit_error_message,
      admissible_flag = FALSE,
      admissibility_reasons = "sampling_error",
      prior_degenerate_flag = prior_pred_main$degeneracy$degenerate_flag[[1]],
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
      cure_model_eligibility_flag = cure_model_eligibility_flag_now,
      primary_gate_method = primary_gate_method_now,
      primary_gate_flag = primary_gate_flag_now,
      receus_primary_class = receus_primary_class_now,
      presence_modifier_flag = presence_modifier_flag_now,
      cure_presence_support_flag = cure_presence_support_flag_now,
      followup_contradiction_flag = followup_contradiction_flag_now,
      followup_not_contradicted_flag = followup_not_contradicted_flag_now,
      screening_note = screening_note_now,
      screening_flag = screening_flag_now,
      screening_detail = screening_detail_now,
      carry_forward_stage8 = carry_forward_stage8_now,
      rds_path = NA_character_
    )
    emit_stage8_progress(
      ii,
      total_models,
      model_id_now,
      paste0("sampling error after ", format_stage8_number(elapsed_stage8_seconds(model_started_at), digits = 1L), "s: ", fit_error_message)
    )
    next
  }

  if (create_diagnostic_pdf) {
    trace_records[[length(trace_records) + 1L]] <- make_trace_record(
      fit = fit,
      model_id = model_id_now,
      family_code = model_row$family_code[[1]],
      K_inc = ncol(design_main$X_inc),
      K_lat = ncol(design_main$X_lat)
    )
  }

  param_names <- c("delta0", "gamma0")
  param_names <- c(param_names, paste0("beta_inc[", seq_len(ncol(design_main$X_inc)), "]"))
  param_names <- c(param_names, paste0("gamma_lat[", seq_len(ncol(design_main$X_lat)), "]"))
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
    dataset = dataset_name,
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
    rho_W = draws_compact$rho_W[keep_draw_idx],
    log_sigma_LN = draws_compact$log_sigma_LN[keep_draw_idx],
    psi_LL = draws_compact$psi_LL[keep_draw_idx]
  )

  linear_terms <- compute_linear_terms(draws_pred, design_main$X_inc, design_main$X_lat, design_main$alpha_gp)

  supported_horizons <- ppc_horizons_for_dataset(dataset_name)
  supported_risk_list <- lapply(supported_horizons, function(h) {
    fh_sup <- family_survival_hazard(h, model_row$family_code[[1]], linear_terms$mu_lat_mat, linear_terms$median_mat, draws_pred)
    1 - ((1 - linear_terms$pi_mat) + linear_terms$pi_mat * fh_sup$Su)
  })
  posterior_degeneracy <- compute_degeneracy(linear_terms$pi_mat, linear_terms$median_mat, supported_risk_list)

  info_criteria <- compute_information_criteria(draws_compact$log_lik)
  waic_val <- info_criteria$waic
  looic_val <- info_criteria$looic

  subject_profile_summary <- summarize_cols_matrix(linear_terms$cure_prob_mat)
  susceptible_prob_summary <- summarize_cols_matrix(linear_terms$pi_mat)
  median_susc_summary <- summarize_cols_matrix(linear_terms$median_mat)

  subject_profile_tbl <- bind_cols(
    tibble(
      dataset = dataset_name,
      model_id = model_id_now
    ),
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

  ppc_model_rows <- list()
  cohort_model_rows <- list()
  subject_year_model_rows <- list()
  class_model_rows <- list()

  for (h in horizons_year) {
    fh <- family_survival_hazard(h, model_row$family_code[[1]], linear_terms$mu_lat_mat, linear_terms$median_mat, draws_pred)
    pop_surv <- (1 - linear_terms$pi_mat) + linear_terms$pi_mat * fh$Su
    risk_mat <- 1 - pop_surv

    subj_surv_summary <- summarize_cols_matrix(pop_surv)
    subj_risk_summary <- summarize_cols_matrix(risk_mat)
    subj_su_summary <- summarize_cols_matrix(fh$Su)

    subject_year_model_rows[[length(subject_year_model_rows) + 1L]] <- bind_cols(
      tibble(
        dataset = dataset_name,
        model_id = model_id_now,
        horizon_year = h,
        horizon_support_label = support_label(dataset_name, h)
      ),
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

    mean_risk_draw <- rowMeans(risk_mat)
    mean_surv_draw <- rowMeans(pop_surv)
    mean_sus_surv_draw <- rowMeans(fh$Su)
    mean_hazard_draw <- rowMeans(fh$haz)
    mean_cure_draw <- rowMeans(linear_terms$cure_prob_mat)

    cohort_model_rows[[length(cohort_model_rows) + 1L]] <- tibble(
      dataset = dataset_name,
      model_id = model_id_now,
      latency_family = model_row$latency_family[[1]],
      formula_anchor = model_row$formula_anchor[[1]],
      horizon_year = h,
      horizon_support_label = support_label(dataset_name, h),
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
      cohort_mean_cure_fraction_q975 = stats::quantile(mean_cure_draw, 0.975, names = FALSE)
    )

    horizon_ref <- ipcw_registry[[dataset_name]] %>% filter(horizon_year == h)

    ppc_model_rows[[length(ppc_model_rows) + 1L]] <- tibble(
      dataset = dataset_name,
      model_id = model_id_now,
      horizon_year = h,
      horizon_support_label = support_label(dataset_name, h),
      observed_km_risk = horizon_ref$observed_km_risk[[1]],
      posterior_mean_risk = mean(mean_risk_draw),
      posterior_q025_risk = stats::quantile(mean_risk_draw, 0.025, names = FALSE),
      posterior_q975_risk = stats::quantile(mean_risk_draw, 0.975, names = FALSE),
      absolute_difference = abs(mean(mean_risk_draw) - horizon_ref$observed_km_risk[[1]]),
      gross_contradiction_flag = (
        (h %in% ppc_horizons_for_dataset(dataset_name)) &&
          (
            horizon_ref$observed_km_risk[[1]] < stats::quantile(mean_risk_draw, 0.025, names = FALSE) ||
              horizon_ref$observed_km_risk[[1]] > stats::quantile(mean_risk_draw, 0.975, names = FALSE)
          ) &&
          abs(mean(mean_risk_draw) - horizon_ref$observed_km_risk[[1]]) > ppc_tolerance_abs
      )
    )

    class_out <- compute_classification_summary(
      risk_draws = risk_mat,
      horizon_row = horizon_ref,
      thresholds = risk_thresholds
    )

    if (nrow(class_out) > 0) {
      class_model_rows[[length(class_model_rows) + 1L]] <- class_out %>%
        mutate(
          dataset = dataset_name,
          model_id = model_id_now,
          latency_family = model_row$latency_family[[1]],
          formula_anchor = model_row$formula_anchor[[1]],
          horizon_year = h,
          horizon_support_label = support_label(dataset_name, h),
          observed_km_risk = horizon_ref$observed_km_risk[[1]]
        ) %>%
        relocate(dataset, model_id, latency_family, formula_anchor, horizon_year, horizon_support_label, threshold)
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

  admissibility_reasons <- c()
  if (prior_pred_main$degeneracy$degenerate_flag[[1]]) admissibility_reasons <- c(admissibility_reasons, "prior_degenerate")
  if (posterior_degeneracy$degenerate_flag[[1]]) admissibility_reasons <- c(admissibility_reasons, "posterior_degenerate")
  if (!is.finite(max_rhat) || max_rhat >= 1.01) admissibility_reasons <- c(admissibility_reasons, "rhat")
  if (!is.finite(min_bulk_ess) || min_bulk_ess < ess_min_threshold) admissibility_reasons <- c(admissibility_reasons, "bulk_ess")
  if (!is.finite(min_tail_ess) || min_tail_ess < ess_min_threshold) admissibility_reasons <- c(admissibility_reasons, "tail_ess")
  if (divergences > 0) admissibility_reasons <- c(admissibility_reasons, "divergences")
  if (treedepth_exceeded > 0) admissibility_reasons <- c(admissibility_reasons, "treedepth")
  if (isTRUE(ppc_gross_contradiction_flag)) admissibility_reasons <- c(admissibility_reasons, "ppc")
  admissible_flag <- length(admissibility_reasons) == 0

  prior_predictive_bundle <- bind_rows(
    prior_pred_main$summary,
    prior_pred_wider$summary
  )

  compact_rds <- list(
    dataset = dataset_name,
    model_id = model_id_now,
    latency_family = model_row$latency_family[[1]],
    formula_anchor = model_row$formula_anchor[[1]],
    fit_status = fit_status,
    admissible_flag = admissible_flag,
    admissibility_reasons = admissibility_reasons,
    diagnostics = tibble(
      divergences = divergences,
      max_rhat = max_rhat,
      min_bulk_ess = min_bulk_ess,
      min_tail_ess = min_tail_ess,
      treedepth_exceeded = treedepth_exceeded,
      waic = waic_val,
      looic = looic_val,
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
      ppc_gross_contradiction_flag = ppc_gross_contradiction_flag,
      prior_degenerate_flag = prior_pred_main$degeneracy$degenerate_flag[[1]],
      posterior_degenerate_flag = posterior_degeneracy$degenerate_flag[[1]]
    ),
    coefficient_summary = coef_tbl,
    diagnostics_parameter_level = diagnostic_param_tbl_model,
    ppc_summary = ppc_model_tbl,
    posterior_subject_profile = if (isTRUE(admissible_flag)) subject_profile_tbl else tibble(),
    posterior_subject_yearly = if (isTRUE(admissible_flag)) subject_year_model_tbl else tibble(),
    posterior_cohort_yearly = if (isTRUE(admissible_flag)) cohort_model_tbl else tibble(),
    posterior_classification = if (isTRUE(admissible_flag)) class_model_tbl else tibble(),
    prior_predictive_summary = prior_predictive_bundle,
    selected_parameter_draws = as.data.frame(
      {
        set.seed(stan_seed + 300000L + ii)
        param_draws_mat[
          sample(seq_len(nrow(param_draws_mat)), size = min(100L, nrow(param_draws_mat))),
          ,
          drop = FALSE
        ]
      }
    )
  )

  rds_path <- file.path(export_path, paste0(model_id_now, "__bayes_stage8_fit.rds"))
  save_obj <- compact_rds
  if (isTRUE(save_full_stanfit_rds)) {
    save_obj$fit <- fit
  }
  saveRDS(save_obj, rds_path)

  main_registry_rows[[length(main_registry_rows) + 1L]] <- tibble(
    dataset = dataset_name,
    model_id = model_id_now,
    structural_model_id = model_row$structural_model_id[[1]],
    latency_family = model_row$latency_family[[1]],
    family_code = model_row$family_code[[1]],
    formula_anchor = model_row$formula_anchor[[1]],
    incidence_site_indicator = model_row$incidence_site_indicator[[1]],
    latency_site_indicator = model_row$latency_site_indicator[[1]],
    latency_interaction_indicator = model_row$latency_interaction_indicator[[1]],
    is_supplementary_branch = model_row$is_supplementary_branch[[1]],
    fit_status = fit_status,
    fit_error_message = fit_error_message,
    admissible_flag = admissible_flag,
    admissibility_reasons = if (length(admissibility_reasons) == 0) "" else paste(admissibility_reasons, collapse = "|"),
    prior_degenerate_flag = prior_pred_main$degeneracy$degenerate_flag[[1]],
    posterior_degenerate_flag = posterior_degeneracy$degenerate_flag[[1]],
    ppc_gross_contradiction_flag = ppc_gross_contradiction_flag,
    divergences = divergences,
    max_rhat = max_rhat,
    min_bulk_ess = min_bulk_ess,
    min_tail_ess = min_tail_ess,
    treedepth_exceeded = treedepth_exceeded,
    waic = waic_val,
    looic = looic_val,
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
    cure_model_eligibility_flag = cure_model_eligibility_flag_now,
    primary_gate_method = primary_gate_method_now,
    primary_gate_flag = primary_gate_flag_now,
    receus_primary_class = receus_primary_class_now,
    presence_modifier_flag = presence_modifier_flag_now,
    cure_presence_support_flag = cure_presence_support_flag_now,
    followup_contradiction_flag = followup_contradiction_flag_now,
    followup_not_contradicted_flag = followup_not_contradicted_flag_now,
    screening_note = screening_note_now,
    screening_flag = screening_flag_now,
    screening_detail = screening_detail_now,
    carry_forward_stage8 = carry_forward_stage8_now,
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

  emit_stage8_progress(
    ii,
    total_models,
    model_id_now,
    paste0(
      "completed; elapsed=",
      format_stage8_number(elapsed_stage8_seconds(model_started_at), digits = 1L),
      "s; WAIC=",
      format_stage8_number(waic_val, digits = 2L),
      "; pWAIC=",
      format_stage8_number(info_criteria$p_waic, digits = 3L),
      "; LOOIC=",
      format_stage8_number(looic_val, digits = 2L),
      "; pLOO=",
      format_stage8_number(info_criteria$p_loo, digits = 3L),
      "; Pareto k max=",
      format_stage8_number(info_criteria$pareto_k_max, digits = 3L)
    )
  )
}
}, warning = function(w) {
  if (is_localhost_connection_warning(w)) {
    tryInvokeRestart("muffleWarning")
  }
})

# 🔴 Export: final tables, metadata enrichment, and diagnostics ===============================
## 🟠 Export: bind primary result tables before annotation ===============================
model_order <- model_grid$model_id

model_registry <- bind_rows_safe(main_registry_rows) %>%
  select(-any_of(c(
    "screening_flag_carry",
    "screening_detail_carry",
    "fit_reused_flag",
    "prior_tail_warning_flag",
    "prior_tail_warning_detail"
  ))) %>%
  mutate(
    rds_path = ifelse(
      file.exists(file.path(export_path, paste0(model_id, "__bayes_stage8_fit.rds"))),
      file.path(export_path, paste0(model_id, "__bayes_stage8_fit.rds")),
      rds_path
    )
  ) %>%
  {
    out <- .
    for (nm in stage8_screening_lookup_fields()) {
      if (!(nm %in% names(out))) {
        out[[nm]] <- if (nm == "carry_forward_stage8") NA else NA_character_
      }
    }
    out
  } %>%
  left_join(
    model_grid %>%
      select(model_id, reuse_existing_fit),
    by = "model_id"
  ) %>%
  rename(fit_reused_flag = reuse_existing_fit) %>%
  left_join(
    screening_lookup %>%
      rename_with(~paste0(.x, "_import"), -model_id),
    by = "model_id"
  ) %>%
  mutate(
    cure_model_eligibility_flag = dplyr::coalesce(cure_model_eligibility_flag, cure_model_eligibility_flag_import),
    primary_gate_method = dplyr::coalesce(primary_gate_method, primary_gate_method_import),
    primary_gate_flag = dplyr::coalesce(primary_gate_flag, primary_gate_flag_import),
    receus_primary_class = dplyr::coalesce(receus_primary_class, receus_primary_class_import),
    presence_modifier_flag = dplyr::coalesce(presence_modifier_flag, presence_modifier_flag_import),
    cure_presence_support_flag = dplyr::coalesce(cure_presence_support_flag, cure_presence_support_flag_import),
    followup_contradiction_flag = dplyr::coalesce(followup_contradiction_flag, followup_contradiction_flag_import),
    followup_not_contradicted_flag = dplyr::coalesce(followup_not_contradicted_flag, followup_not_contradicted_flag_import),
    screening_note = dplyr::coalesce(screening_note, screening_note_import),
    screening_flag = dplyr::coalesce(screening_flag, screening_flag_import),
    screening_detail = dplyr::coalesce(screening_detail, screening_detail_import),
    carry_forward_stage8 = dplyr::coalesce(as.logical(carry_forward_stage8), as.logical(carry_forward_stage8_import))
  ) %>%
  select(-any_of(paste0(stage8_screening_lookup_fields(), "_import"))) %>%
  arrange(factor(model_id, levels = model_order))

coefficient_summary <- bind_rows_safe(coef_rows) %>%
  arrange(factor(model_id, levels = model_order), parameter)

posterior_subject_profile <- bind_rows_safe(subject_profile_rows) %>%
  arrange(factor(model_id, levels = model_order), unique_person_id)

posterior_subject_yearly <- bind_rows_safe(subject_yearly_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon_year, unique_person_id)

posterior_cohort_yearly <- bind_rows_safe(cohort_yearly_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon_year)

posterior_classification <- bind_rows_safe(class_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon_year, threshold)

diagnostics_parameter_level <- bind_rows_safe(diagnostic_rows) %>%
  arrange(factor(model_id, levels = model_order), parameter)

ppc_summary <- bind_rows_safe(ppc_rows_all) %>%
  arrange(factor(model_id, levels = model_order), horizon_year)

prior_predictive_summary <- normalize_prior_predictive_summary(bind_rows_safe(prior_predictive_rows)) %>%
  arrange(factor(model_id, levels = model_order), metric, prior_set, horizon_year)

## 🟠 Export: add prior-warning columns without forcing refits ===============================
prior_predictive_summary <- annotate_prior_predictive_summary(
  df = prior_predictive_summary,
  mean_threshold_years = prior_tail_warning_mean_years,
  q975_threshold_years = prior_tail_warning_q975_years
) %>%
  mutate(
    mean = safe_numeric(.data$mean),
    q025 = safe_numeric(.data$q025),
    q50 = safe_numeric(.data$q50),
    q975 = safe_numeric(.data$q975),
    prior_degenerate_flag = as.logical(.data$prior_degenerate_flag),
    prior_tail_warning_flag = as.logical(.data$prior_tail_warning_flag)
  ) %>%
  distinct(dataset, model_id, metric, horizon_year, prior_set, .keep_all = TRUE) %>%
  arrange(factor(model_id, levels = model_order), metric, prior_set, horizon_year)

prior_sensitivity_summary <- prior_predictive_summary %>%
  select(
    dataset,
    model_id,
    metric,
    horizon_year,
    prior_set,
    mean,
    q025,
    q975,
    prior_degenerate_flag,
    prior_tail_warning_flag
  ) %>%
  tidyr::pivot_wider(
    names_from = prior_set,
    values_from = c(mean, q025, q975, prior_degenerate_flag, prior_tail_warning_flag),
    names_sep = "__"
  ) %>%
  mutate(
    delta_mean_wider_minus_main = mean__wider_latency - mean__main,
    prior_tail_warning_flag_any = dplyr::coalesce(prior_tail_warning_flag__main, FALSE) |
      dplyr::coalesce(prior_tail_warning_flag__wider_latency, FALSE)
  ) %>%
  arrange(factor(model_id, levels = model_order), metric, horizon_year)

prior_warning_by_model <- prior_predictive_summary %>%
  group_by(dataset, model_id) %>%
  summarise(
    prior_tail_warning_flag = any(dplyr::coalesce(prior_tail_warning_flag, FALSE)),
    prior_tail_warning_detail = paste(
      unique(na.omit(prior_tail_warning_detail[prior_tail_warning_flag %in% TRUE])),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  mutate(
    prior_tail_warning_detail = dplyr::na_if(prior_tail_warning_detail, "")
  )

model_registry <- model_registry %>%
  left_join(prior_warning_by_model, by = c("dataset", "model_id")) %>%
  mutate(
    prior_tail_warning_flag = dplyr::coalesce(prior_tail_warning_flag, FALSE)
  ) %>%
  arrange(factor(model_id, levels = model_order))

## 🟠 Export: enrich outputs with carry-forward and reporting metadata ===============================
model_annotation <- model_registry %>%
  select(
    model_id,
    admissible_flag,
    admissibility_reasons,
    screening_flag,
    screening_detail,
    fit_reused_flag
  )

cohort_horizon_annotation <- horizon_metadata_registry %>%
  select(
    dataset,
    horizon_year,
    horizon_support_label,
    support_priority,
    cohort_reporting_status,
    cohort_reporting_note
  )

classification_horizon_annotation <- horizon_metadata_registry %>%
  select(
    dataset,
    horizon_year,
    horizon_support_label,
    support_priority,
    denom_case,
    denom_control,
    classification_estimable_flag,
    classification_reporting_status,
    classification_reporting_note
  )

posterior_subject_yearly <- posterior_subject_yearly %>%
  left_join_replacing_columns(
    cohort_horizon_annotation %>%
      select(dataset, horizon_year, horizon_support_label, support_priority),
    by = c("dataset", "horizon_year", "horizon_support_label")
  ) %>%
  arrange(factor(model_id, levels = model_order), horizon_year, unique_person_id)

posterior_cohort_yearly <- posterior_cohort_yearly %>%
  left_join_replacing_columns(
    cohort_horizon_annotation,
    by = c("dataset", "horizon_year", "horizon_support_label")
  ) %>%
  left_join_replacing_columns(
    model_annotation,
    by = "model_id"
  ) %>%
  arrange(factor(model_id, levels = model_order), horizon_year)

posterior_classification <- posterior_classification %>%
  left_join_replacing_columns(
    classification_horizon_annotation,
    by = c("dataset", "horizon_year", "horizon_support_label")
  ) %>%
  left_join_replacing_columns(
    model_annotation,
    by = "model_id"
  ) %>%
  arrange(factor(model_id, levels = model_order), horizon_year, threshold)

ppc_summary <- ppc_summary %>%
  left_join_replacing_columns(
    cohort_horizon_annotation,
    by = c("dataset", "horizon_year", "horizon_support_label")
  ) %>%
  left_join_replacing_columns(
    model_annotation,
    by = "model_id"
  ) %>%
  arrange(factor(model_id, levels = model_order), horizon_year)

## 🟠 Export: rebuild Stage-10 comparison file with support metadata ===============================
posterior_delta_vs_nocure <- tryCatch(
  rebuild_delta_vs_nocure(
    posterior_cohort_yearly = posterior_cohort_yearly,
    posterior_classification = posterior_classification,
    nocure_cohort_long = nocure_cohort_long,
    nocure_class_long = nocure_class_long
  ),
  error = function(e) {
    warning(
      paste0(
        "Stage 8 could not rebuild `posterior_delta_vs_nocure`; an empty schema-preserving table will be written instead. Reason: ",
        conditionMessage(e)
      ),
      call. = FALSE
    )
    empty_delta_template()
  }
)

delta_horizon_meta <- bind_rows_safe(list(
  if (nrow_or_zero(posterior_cohort_yearly) > 0L) {
    posterior_cohort_yearly %>%
      transmute(
        dataset = dataset,
        model_id = model_id,
        formula_anchor = formula_anchor,
        horizon_year = horizon_year,
        threshold_key = make_threshold_key(NA_real_),
        latency_family = latency_family,
        horizon_support_label = horizon_support_label,
        support_priority = support_priority,
        reporting_status = cohort_reporting_status,
        reporting_note = cohort_reporting_note,
        classification_estimable_flag = NA,
        denom_case = NA_real_,
        denom_control = NA_real_
      )
  } else {
    NULL
  },
  if (nrow_or_zero(posterior_classification) > 0L) {
    posterior_classification %>%
      transmute(
        dataset = dataset,
        model_id = model_id,
        formula_anchor = formula_anchor,
        horizon_year = horizon_year,
        threshold_key = make_threshold_key(threshold),
        latency_family = latency_family,
        horizon_support_label = horizon_support_label,
        support_priority = support_priority,
        reporting_status = classification_reporting_status,
        reporting_note = classification_reporting_note,
        classification_estimable_flag = classification_estimable_flag,
        denom_case = denom_case,
        denom_control = denom_control
      )
  } else {
    NULL
  }
)) %>%
  distinct()

posterior_delta_vs_nocure <- posterior_delta_vs_nocure %>%
  mutate(threshold_key = make_threshold_key(threshold)) %>%
  left_join_replacing_columns(
    delta_horizon_meta,
    by = c("dataset", "model_id", "formula_anchor", "horizon_year", "threshold_key")
  ) %>%
  select(-any_of("threshold_key")) %>%
  left_join_replacing_columns(
    model_annotation %>%
      select(model_id, admissible_flag, admissibility_reasons, screening_flag, screening_detail),
    by = "model_id"
  ) %>%
  arrange(
    dataset,
    factor(model_id, levels = model_order),
    horizon_year,
    threshold,
    metric,
    no_cure_model_id
  )

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
prior_sensitivity_summary <- simplify_scalar_list_cols(prior_sensitivity_summary)

## 🟠 Export: regenerate diagnostic PDF safely without refitting ===============================
if (create_diagnostic_pdf) {
  pdf_ok <- tryCatch(
    {
      safe_generate_diagnostic_pdf(
        trace_records = trace_records,
        n_models_reused = n_models_reused,
        n_models_to_fit = n_models_to_fit,
        posterior_cohort_yearly = posterior_cohort_yearly,
        posterior_classification = posterior_classification,
        ppc_summary = ppc_summary,
        final_path = diagnostic_pdf_path
      )
      TRUE
    },
    error = function(e) {
      warning(
        paste0(
          "Diagnostic PDF was not regenerated. Existing PDF was preserved. Reason: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
      FALSE
    }
  )

  if (isTRUE(pdf_ok)) {
    message("Stage 8 diagnostic PDF regenerated safely: ", diagnostic_pdf_path)
  }
}

## 🟠 Export: finalize reused-fit QC harvest, carry-forward metadata, and audit ===============================
final_stage8 <- stage8_finalize_registry_and_audit(
  model_grid = model_grid,
  model_registry = model_registry,
  posterior_cohort_yearly = posterior_cohort_yearly,
  posterior_classification = posterior_classification,
  ppc_summary = ppc_summary,
  posterior_delta_vs_nocure = posterior_delta_vs_nocure,
  diagnostic_pdf_path = diagnostic_pdf_path,
  horizons_year = horizons_year,
  risk_thresholds = risk_thresholds
)

model_registry <- simplify_scalar_list_cols(final_stage8$model_registry)
output_audit <- simplify_scalar_list_cols(final_stage8$output_audit)

carryforward_metadata <- model_registry %>%
  select(
    dataset,
    model_id,
    latency_family,
    formula_anchor,
    incidence_site_indicator,
    latency_site_indicator,
    latency_interaction_indicator,
    admissible_flag,
    admissibility_reasons,
    cure_model_eligibility_flag,
    primary_gate_method,
    primary_gate_flag,
    receus_primary_class,
    presence_modifier_flag,
    cure_presence_support_flag,
    followup_contradiction_flag,
    followup_not_contradicted_flag,
    screening_note,
    screening_flag,
    screening_detail,
    carry_forward_stage8,
    fit_reused_flag,
    prior_tail_warning_flag,
    prior_tail_warning_detail
  )

## 🟠 Export: write CSV outputs ===============================
write_csv_preserve_schema(dataset_registry, file.path(export_path, "bayes_stage8_dataset_registry.csv"))
write_csv_preserve_schema(scaling_registry, file.path(export_path, "bayes_stage8_scaling_registry.csv"))
write_csv_preserve_schema(horizon_metadata_registry, file.path(export_path, "bayes_stage8_horizon_metadata.csv"))
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
write_csv_preserve_schema(prior_sensitivity_summary, file.path(export_path, "bayes_stage8_prior_sensitivity_summary.csv"))
write_csv_preserve_schema(carryforward_metadata, file.path(export_path, "bayes_stage8_carryforward_metadata.csv"))
write_csv_preserve_schema(output_audit, file.path(export_path, "bayes_stage8_output_audit.csv"))

message(
  "Stage 8 completed. Reused models: ",
  n_models_reused,
  "; newly fitted models: ",
  n_models_to_fit,
  "."
)
