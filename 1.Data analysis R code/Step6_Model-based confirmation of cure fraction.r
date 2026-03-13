# 🔴 Configure: Step6 runtime and bundle paths ===============================

## 🟠 Declare: user-editable bundle directories and boundary-test controls ===============================
options(stringsAsFactors = FALSE, scipen = 999)

auto_install_packages <- FALSE
verbose <- TRUE

step4_bundle_dir <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦4.Step4_Non-cure benchmark models/attachments'

step5_bundle_dir <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦5.Step5_MLE MCM/attachments'

export_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦6.Step6_formal cure fraction test/attachments'

step3_gate_status <- c(
  merged = "pass",
  pnu = "pass",
  snu = "pass"
)

run_all_pairs <- TRUE
strict_bundle_alignment <- TRUE

boundary_alpha <- 0.05

strict_bootstrap_success_target_primary <- 999L
strict_bootstrap_success_target_sensitivity <- 999L

relaxed_bootstrap_success_target_primary <- 999L
relaxed_bootstrap_success_target_sensitivity <- 999L

bootstrap_max_attempt_multiplier <- 4L
censoring_bootstrap_method <- "reverse_km"
failure_log_max_rows <- 5000L

zero_time_offset_fallback_days <- 1e-03
cure_optim_method <- "BFGS"
optim_maxit <- 1000L
survreg_maxiter <- 100L
alignment_tolerance <- 1e-10
random_seed <- 20260312L
save_rds_compression <- "xz"

dataset_labels <- c("merged", "pnu", "snu")
expected_step4_dataset_label <- c(merged = "merged", pnu = "PNU", snu = "SNU")
expected_step5_dataset_label <- c(merged = "merged", pnu = "pnu", snu = "snu")

file_dataset_manifest <- file.path(export_path, "step6_dataset_manifest.csv")
file_primary_results <- file.path(export_path, "step6_cure_fraction_primary.csv")
file_allpairs_results <- file.path(export_path, "step6_cure_fraction_all_pairs.csv")

dataset_rds_paths <- c(
  merged = file.path(export_path, "step6_test_merged.rds"),
  pnu = file.path(export_path, "step6_test_pnu.rds"),
  snu = file.path(export_path, "step6_test_snu.rds")
)

## 🟠 Build: automatic dataset-level bundle paths from directory inputs ===============================
build_bundle_paths_from_dir <- function(bundle_dir, file_prefix) {
  if (!dir.exists(bundle_dir)) {
    stop("Bundle directory not found: ", bundle_dir, call. = FALSE)
  }
  
  out <- c(
    merged = file.path(bundle_dir, paste0(file_prefix, "_merged.rds")),
    pnu = file.path(bundle_dir, paste0(file_prefix, "_pnu.rds")),
    snu = file.path(bundle_dir, paste0(file_prefix, "_snu.rds"))
  )
  
  missing_files <- out[!file.exists(out)]
  
  if (length(missing_files) > 0L) {
    stop(
      paste0(
        "Required bundle file(s) not found under ",
        bundle_dir,
        ": ",
        paste(basename(missing_files), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  out
}

step4_bundle_paths <- build_bundle_paths_from_dir(
  bundle_dir = step4_bundle_dir,
  file_prefix = "step4_nc"
)

step5_bundle_paths <- build_bundle_paths_from_dir(
  bundle_dir = step5_bundle_dir,
  file_prefix = "step5_models"
)

if (!dir.exists(export_path)) {
  dir.create(export_path, recursive = TRUE, showWarnings = FALSE)
}

## 🟠 Validate: package namespaces and control inputs ===============================
required_packages <- c(
  "survival",
  "flexsurvcure"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_packages) > 0L) {
  if (isTRUE(auto_install_packages)) {
    install.packages(missing_packages)
  } else {
    stop(
      paste0(
        "Missing required package(s): ",
        paste(missing_packages, collapse = ", "),
        ". Install them first or set auto_install_packages <- TRUE."
      ),
      call. = FALSE
    )
  }
}

library(survival)
library(flexsurvcure)

if (!setequal(names(step4_bundle_paths), dataset_labels)) {
  stop("step4_bundle_paths must be a named character vector with names: merged, pnu, snu", call. = FALSE)
}

if (!setequal(names(step5_bundle_paths), dataset_labels)) {
  stop("step5_bundle_paths must be a named character vector with names: merged, pnu, snu", call. = FALSE)
}

if (!setequal(names(step3_gate_status), dataset_labels)) {
  stop("step3_gate_status must be a named character vector with names: merged, pnu, snu", call. = FALSE)
}

if (!isTRUE(is.logical(run_all_pairs)) || length(run_all_pairs) != 1L) {
  stop("run_all_pairs must be a single logical value.", call. = FALSE)
}

if (!isTRUE(is.logical(strict_bundle_alignment)) || length(strict_bundle_alignment) != 1L) {
  stop("strict_bundle_alignment must be a single logical value.", call. = FALSE)
}

if (!is.finite(boundary_alpha) || boundary_alpha <= 0 || boundary_alpha >= 1) {
  stop("boundary_alpha must be a single number between 0 and 1.", call. = FALSE)
}

for (nm in c(
  "strict_bootstrap_success_target_primary",
  "strict_bootstrap_success_target_sensitivity",
  "relaxed_bootstrap_success_target_primary",
  "relaxed_bootstrap_success_target_sensitivity"
)) {
  val <- get(nm, inherits = FALSE)
  if (!is.finite(val) || val < 1L) {
    stop(nm, " must be a positive integer.", call. = FALSE)
  }
}

if (!is.finite(bootstrap_max_attempt_multiplier) || bootstrap_max_attempt_multiplier < 1L) {
  stop("bootstrap_max_attempt_multiplier must be at least 1.", call. = FALSE)
}

if (!(tolower(censoring_bootstrap_method) %in% c("reverse_km", "empirical_censoring_resample"))) {
  stop("censoring_bootstrap_method must be one of: reverse_km, empirical_censoring_resample", call. = FALSE)
}

if (!is.finite(zero_time_offset_fallback_days) || zero_time_offset_fallback_days <= 0) {
  stop("zero_time_offset_fallback_days must be a positive number.", call. = FALSE)
}

if (!is.finite(optim_maxit) || optim_maxit < 1L) {
  stop("optim_maxit must be a positive integer.", call. = FALSE)
}

if (!is.finite(survreg_maxiter) || survreg_maxiter < 1L) {
  stop("survreg_maxiter must be a positive integer.", call. = FALSE)
}

if (!is.finite(alignment_tolerance) || alignment_tolerance < 0) {
  stop("alignment_tolerance must be a non-negative number.", call. = FALSE)
}

if (!is.finite(random_seed)) {
  stop("random_seed must be finite.", call. = FALSE)
}

if (!(save_rds_compression %in% c("gzip", "bzip2", "xz", FALSE, TRUE))) {
  stop("save_rds_compression must be one of gzip, bzip2, xz, TRUE, FALSE.", call. = FALSE)
}

for (ds in dataset_labels) {
  if (!file.exists(step4_bundle_paths[[ds]])) {
    stop("Step4 bundle not found for dataset ", ds, ": ", step4_bundle_paths[[ds]], call. = FALSE)
  }
  if (!file.exists(step5_bundle_paths[[ds]])) {
    stop("Step5 bundle not found for dataset ", ds, ": ", step5_bundle_paths[[ds]], call. = FALSE)
  }
}

set.seed(random_seed)

# 🔴 Define: shared utilities for bundle validation and result rows ===============================

## 🟠 Build: coercion helpers, binders, and typed row templates ===============================
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) {
    return(y)
  }
  x
}

vmessage <- function(...) {
  if (isTRUE(verbose)) {
    message(...)
  }
}

blank_to_na <- function(x) {
  x <- as.character(x[!is.na(x)])
  x <- trimws(x)
  x <- x[nzchar(x)]
  x <- x[!toupper(x) %in% c("NA", "NULL")]
  x <- unique(x)
  if (length(x) == 0L) {
    return(NA_character_)
  }
  paste(x, collapse = " | ")
}

truncate_text <- function(x, width = 160L) {
  x <- as.character(x)
  x <- ifelse(
    is.na(x),
    NA_character_,
    ifelse(nchar(x) > width, paste0(substr(x, 1L, width - 3L), "..."), x)
  )
  x
}

note_if_false <- function(flag, text) {
  if (isFALSE(flag)) {
    return(text)
  }
  NA_character_
}

bind_rows_fill <- function(df_list) {
  df_list <- Filter(Negate(is.null), df_list)
  if (length(df_list) == 0L) {
    return(data.frame())
  }
  
  all_names <- unique(unlist(lapply(df_list, names), use.names = FALSE))
  
  df_list_filled <- lapply(df_list, function(d) {
    d <- as.data.frame(d, stringsAsFactors = FALSE)
    missing_names <- setdiff(all_names, names(d))
    if (length(missing_names) > 0L) {
      for (nm in missing_names) {
        d[[nm]] <- NA
      }
    }
    d <- d[, all_names, drop = FALSE]
    rownames(d) <- NULL
    d
  })
  
  out <- do.call(rbind, df_list_filled)
  rownames(out) <- NULL
  out
}

normalize_gate_status <- function(x) {
  y <- tolower(trimws(as.character(x)[1L]))
  if (!y %in% c("pass", "borderline", "fail")) {
    stop("Invalid gate status: ", x, ". Allowed values are pass, borderline, fail.", call. = FALSE)
  }
  y
}

make_subject_key <- function(data) {
  paste(
    toupper(trimws(as.character(data$site))),
    trimws(as.character(data$id)),
    sep = "::"
  )
}

detect_instability_message <- function(text) {
  if (length(text) == 0L || all(is.na(text))) {
    return(FALSE)
  }
  
  grepl(
    pattern = paste(
      c(
        "Hessian",
        "not positive definite",
        "non-positive definite",
        "false convergence",
        "failed",
        "singular",
        "NaNs produced",
        "non-finite",
        "cannot be evaluated",
        "did not converge",
        "ran out of iterations",
        "optimization failed"
      ),
      collapse = "|"
    ),
    x = text,
    ignore.case = TRUE
  )
}

new_empty_failure_log <- function() {
  data.frame(
    attempt = integer(0),
    stage = character(0),
    reason = character(0),
    stringsAsFactors = FALSE
  )
}

new_empty_failure_counts <- function() {
  data.frame(
    stage = character(0),
    reason = character(0),
    n = integer(0),
    stringsAsFactors = FALSE
  )
}

new_bootstrap_result <- function(target_B) {
  list(
    d_boot = numeric(0),
    attempts = 0L,
    success_B = 0L,
    failed_B = 0L,
    target_B = as.integer(target_B),
    target_met = FALSE,
    failure_log = new_empty_failure_log(),
    failure_counts = new_empty_failure_counts(),
    failure_summary = NA_character_,
    accepted_unstable_null_B = 0L,
    accepted_unstable_alt_B = 0L,
    accepted_any_unstable_B = 0L,
    accepted_all_stable_B = 0L
  )
}

pair_result_template <- function() {
  list(
    dataset = NA_character_,
    gate_status_step3 = NA_character_,
    step5_gate_pass = NA,
    step5_selected_best_cure_id = NA_character_,
    is_primary_pair = NA,
    alt_model_id = NA_character_,
    null_model_id = NA_character_,
    latency_dist = NA_character_,
    null_survreg_dist = NA_character_,
    step4_family = NA_character_,
    covariate_structure = NA_character_,
    analysis_n = as.integer(NA),
    n_event = as.integer(NA),
    n_censor = as.integer(NA),
    zero_time_offset_days = NA_real_,
    bootstrap_method = NA_character_,
    null_fit_ok = NA,
    null_converged = NA,
    null_stable = NA,
    alt_fit_ok = NA,
    alt_converged = NA,
    alt_stable = NA,
    step4_matched_ok = NA,
    step4_matched_converged = NA,
    step4_step5_subject_set_match = NA,
    step4_step5_endpoint_match = NA,
    step4_step5_transformed_covariate_match = NA,
    step4_step5_scaling_match = NA,
    logLik_null = NA_real_,
    logLik_alt = NA_real_,
    df_model_null = NA_real_,
    df_model_alt = NA_real_,
    AIC_null = NA_real_,
    AIC_alt = NA_real_,
    BIC_null = NA_real_,
    BIC_alt = NA_real_,
    d_obs = NA_real_,
    strict_observed_evaluable = NA,
    strict_bootstrap_target_B = as.integer(NA),
    strict_bootstrap_success_B = as.integer(NA),
    strict_bootstrap_failed_B = as.integer(NA),
    strict_bootstrap_attempts = as.integer(NA),
    strict_bootstrap_target_met = NA,
    strict_bootstrap_failure_rate = NA_real_,
    strict_bootstrap_null_mean = NA_real_,
    strict_bootstrap_null_median = NA_real_,
    strict_bootstrap_null_q95 = NA_real_,
    strict_bootstrap_null_q975 = NA_real_,
    strict_p_value_boundary = NA_real_,
    strict_alpha = NA_real_,
    strict_decision = NA_character_,
    strict_reading_role = NA_character_,
    strict_bootstrap_failure_summary = NA_character_,
    relaxed_observed_evaluable = NA,
    relaxed_bootstrap_target_B = as.integer(NA),
    relaxed_bootstrap_success_B = as.integer(NA),
    relaxed_bootstrap_failed_B = as.integer(NA),
    relaxed_bootstrap_attempts = as.integer(NA),
    relaxed_bootstrap_target_met = NA,
    relaxed_bootstrap_failure_rate = NA_real_,
    relaxed_bootstrap_null_mean = NA_real_,
    relaxed_bootstrap_null_median = NA_real_,
    relaxed_bootstrap_null_q95 = NA_real_,
    relaxed_bootstrap_null_q975 = NA_real_,
    relaxed_p_value_boundary = NA_real_,
    relaxed_alpha = NA_real_,
    relaxed_decision = NA_character_,
    relaxed_reading_role = NA_character_,
    relaxed_bootstrap_failure_summary = NA_character_,
    relaxed_bootstrap_accepted_unstable_null_B = as.integer(NA),
    relaxed_bootstrap_accepted_unstable_alt_B = as.integer(NA),
    relaxed_bootstrap_accepted_any_unstable_B = as.integer(NA),
    relaxed_bootstrap_accepted_all_stable_B = as.integer(NA),
    null_warning_message = NA_character_,
    alt_warning_message = NA_character_,
    null_error_message = NA_character_,
    alt_error_message = NA_character_,
    note = NA_character_
  )
}

build_pair_result_row <- function(values = list()) {
  template <- pair_result_template()
  row_list <- utils::modifyList(template, values)
  
  row_list <- Map(
    f = function(val, tmpl) {
      if (is.null(val) || length(val) == 0L) {
        return(tmpl)
      }
      val[1L]
    },
    val = row_list,
    tmpl = template
  )
  
  as.data.frame(row_list, stringsAsFactors = FALSE, check.names = FALSE)
}

new_empty_pair_result <- function() {
  out <- build_pair_result_row()
  out[FALSE, , drop = FALSE]
}

summarize_bootstrap_result <- function(bootstrap_res) {
  d <- bootstrap_res$d_boot
  
  list(
    failure_rate = if (bootstrap_res$attempts > 0L) bootstrap_res$failed_B / bootstrap_res$attempts else NA_real_,
    null_mean = if (length(d) > 0L) mean(d) else NA_real_,
    null_median = if (length(d) > 0L) stats::median(d) else NA_real_,
    null_q95 = if (length(d) > 0L) as.numeric(stats::quantile(d, probs = 0.95, names = FALSE, type = 7)) else NA_real_,
    null_q975 = if (length(d) > 0L) as.numeric(stats::quantile(d, probs = 0.975, names = FALSE, type = 7)) else NA_real_
  )
}

# 🔴 Define: Step6 model pair catalog and formula generators ===============================

## 🟠 Build: matched boundary-test catalog and formula makers ===============================
make_pair_catalog <- function() {
  data.frame(
    alt_model_id = c(
      "cure_weibull_main",
      "cure_weibull_interaction",
      "cure_lnorm_main",
      "cure_lnorm_interaction",
      "cure_llogis_main",
      "cure_llogis_interaction"
    ),
    null_model_id = c(
      "nc_weibull_main",
      "nc_weibull_interaction",
      "nc_lognormal_main",
      "nc_lognormal_interaction",
      "nc_loglogistic_main",
      "nc_loglogistic_interaction"
    ),
    latency_dist = c(
      "weibull",
      "weibull",
      "lnorm",
      "lnorm",
      "llogis",
      "llogis"
    ),
    null_survreg_dist = c(
      "weibull",
      "weibull",
      "lognormal",
      "lognormal",
      "loglogistic",
      "loglogistic"
    ),
    step4_family = c(
      "weibull_aft",
      "weibull_aft",
      "lognormal_aft",
      "lognormal_aft",
      "loglogistic_aft",
      "loglogistic_aft"
    ),
    step4_specification = c(
      "main_effects",
      "sex_age_interaction",
      "main_effects",
      "sex_age_interaction",
      "main_effects",
      "sex_age_interaction"
    ),
    covariate_structure = c(
      "main",
      "interaction",
      "main",
      "interaction",
      "main",
      "interaction"
    ),
    latency_param = c(
      "scale",
      "scale",
      "meanlog",
      "meanlog",
      "scale",
      "scale"
    ),
    stringsAsFactors = FALSE
  )
}

make_rhs_text <- function(covariate_structure) {
  if (identical(covariate_structure, "main")) {
    return("c_sex + z_age")
  }
  if (identical(covariate_structure, "interaction")) {
    return("c_sex + z_age + int_sex_age")
  }
  stop("Unknown covariate_structure: ", covariate_structure, call. = FALSE)
}

make_surv_formula <- function(covariate_structure) {
  stats::as.formula(
    paste("survival::Surv(time_model, event_primary) ~", make_rhs_text(covariate_structure))
  )
}

make_anc_list <- function(latency_param, covariate_structure) {
  anc_formula <- stats::as.formula(paste("~", make_rhs_text(covariate_structure)))
  setNames(list(anc_formula), latency_param)
}

select_pair_queue <- function(pair_catalog, primary_alt_id, run_all_pairs) {
  if (!isTRUE(run_all_pairs)) {
    if (is.na(primary_alt_id) || !primary_alt_id %in% pair_catalog$alt_model_id) {
      return(pair_catalog[FALSE, , drop = FALSE])
    }
    out <- pair_catalog[pair_catalog$alt_model_id == primary_alt_id, , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  
  if (!is.na(primary_alt_id) && primary_alt_id %in% pair_catalog$alt_model_id) {
    first_row <- pair_catalog[pair_catalog$alt_model_id == primary_alt_id, , drop = FALSE]
    other_rows <- pair_catalog[pair_catalog$alt_model_id != primary_alt_id, , drop = FALSE]
    out <- rbind(first_row, other_rows)
    rownames(out) <- NULL
    return(out)
  }
  
  pair_catalog
}

make_pair_seed <- function(base_seed, dataset_name, alt_model_id, mode_label) {
  seed_raw <- as.integer(base_seed) + sum(utf8ToInt(paste(dataset_name, alt_model_id, mode_label, sep = "::")))
  seed_mod <- seed_raw %% .Machine$integer.max
  if (!is.finite(seed_mod) || seed_mod <= 0L) {
    seed_mod <- 1L
  }
  as.integer(seed_mod)
}

# 🔴 Define: bundle readers and alignment checks ===============================

## 🟠 Build: RDS readers, structure checks, and aligned analysis frames ===============================
assert_step4_bundle_structure <- function(bundle, dataset_name) {
  required_names <- c(
    "dataset",
    "data_complete_case",
    "standardization_constants",
    "fit_nc_list"
  )
  missing_names <- setdiff(required_names, names(bundle))
  if (length(missing_names) > 0L) {
    stop(
      "Step4 bundle for dataset ",
      dataset_name,
      " is missing required element(s): ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }
  
  expected_label <- expected_step4_dataset_label[[dataset_name]]
  actual_label <- as.character(bundle$dataset[1L])
  
  if (!identical(actual_label, expected_label)) {
    stop(
      "Step4 bundle dataset label mismatch for ",
      dataset_name,
      ". Expected ",
      expected_label,
      " but found ",
      actual_label,
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

assert_step5_bundle_structure <- function(bundle, dataset_name) {
  required_names <- c(
    "dataset_label",
    "analysis_frame",
    "standardization",
    "selected_models",
    "gate"
  )
  missing_names <- setdiff(required_names, names(bundle))
  if (length(missing_names) > 0L) {
    stop(
      "Step5 bundle for dataset ",
      dataset_name,
      " is missing required element(s): ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }
  
  expected_label <- expected_step5_dataset_label[[dataset_name]]
  actual_label <- as.character(bundle$dataset_label[1L])
  
  if (!identical(actual_label, expected_label)) {
    stop(
      "Step5 bundle dataset label mismatch for ",
      dataset_name,
      ". Expected ",
      expected_label,
      " but found ",
      actual_label,
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

extract_step5_selected_cure_id <- function(step5_bundle) {
  selected_id <- blank_to_na(step5_bundle$selected_models$best_cure_id %||% NA_character_)
  if (!is.na(selected_id)) {
    return(selected_id)
  }
  
  if (!is.null(step5_bundle$model_summary)) {
    model_summary <- as.data.frame(step5_bundle$model_summary, stringsAsFactors = FALSE)
    hit <- model_summary[model_summary$is_best_cure %in% TRUE, , drop = FALSE]
    if (nrow(hit) > 0L) {
      return(as.character(hit$model_id[1L]))
    }
  }
  
  NA_character_
}

extract_step5_selection_rule <- function(step5_bundle) {
  blank_to_na(step5_bundle$selected_models$best_cure_selection_rule %||% NA_character_)
}

prepare_step5_analysis_frame <- function(step5_bundle, zero_time_offset_fallback_days) {
  dat <- as.data.frame(step5_bundle$analysis_frame, stringsAsFactors = FALSE)
  
  required_cols <- c(
    "id",
    "site",
    "event_primary",
    "c_sex",
    "z_age"
  )
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0L) {
    stop(
      "Step5 analysis_frame is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  zero_time_offset_days <- suppressWarnings(as.numeric(
    step5_bundle$charter$zero_time_offset_days %||%
      step5_bundle$standardization$zero_time_offset_days[1L] %||%
      zero_time_offset_fallback_days
  ))
  
  if (!is.finite(zero_time_offset_days) || zero_time_offset_days <= 0) {
    zero_time_offset_days <- zero_time_offset_fallback_days
  }
  
  dat$id <- as.character(dat$id)
  dat$site <- toupper(trimws(as.character(dat$site)))
  
  if (!("time_primary" %in% names(dat)) && ("time_model" %in% names(dat))) {
    dat$time_primary <- as.numeric(dat$time_model)
  }
  
  if (!("time_model" %in% names(dat)) && ("time_primary" %in% names(dat))) {
    dat$time_model <- pmax(as.numeric(dat$time_primary), zero_time_offset_days)
  }
  
  if (!(("time_primary" %in% names(dat)) && ("time_model" %in% names(dat)))) {
    stop("Step5 analysis_frame must contain time_primary or time_model.", call. = FALSE)
  }
  
  dat$time_primary <- as.numeric(dat$time_primary)
  dat$time_model <- pmax(as.numeric(dat$time_model), zero_time_offset_days)
  dat$event_primary <- as.integer(as.numeric(dat$event_primary))
  dat$c_sex <- as.numeric(dat$c_sex)
  dat$z_age <- as.numeric(dat$z_age)
  
  if (!("int_sex_age" %in% names(dat))) {
    dat$int_sex_age <- dat$c_sex * dat$z_age
  } else {
    dat$int_sex_age <- as.numeric(dat$int_sex_age)
  }
  
  if (!("status_primary" %in% names(dat))) {
    dat$status_primary <- ifelse(dat$event_primary == 1L, "transition", "censored")
  } else {
    dat$status_primary <- as.character(dat$status_primary)
  }
  
  dat$subject_key <- make_subject_key(dat)
  
  if (anyDuplicated(dat$subject_key) > 0L) {
    duplicated_example <- dat$subject_key[duplicated(dat$subject_key)][1L]
    stop("Duplicated subject key detected in Step5 analysis_frame: ", duplicated_example, call. = FALSE)
  }
  
  if (any(!is.finite(dat$time_model)) || any(dat$time_model <= 0)) {
    stop("Step5 analysis_frame contains non-finite or non-positive time_model values.", call. = FALSE)
  }
  
  if (any(!dat$event_primary %in% c(0L, 1L))) {
    stop("Step5 analysis_frame contains invalid event_primary values.", call. = FALSE)
  }
  
  if (any(!is.finite(dat$c_sex)) || any(!is.finite(dat$z_age)) || any(!is.finite(dat$int_sex_age))) {
    stop("Step5 analysis_frame contains non-finite transformed covariates.", call. = FALSE)
  }
  
  dat <- dat[order(dat$subject_key), , drop = FALSE]
  rownames(dat) <- NULL
  
  list(
    data = dat,
    zero_time_offset_days = zero_time_offset_days
  )
}

extract_step4_scaling_vector <- function(step4_bundle) {
  df <- as.data.frame(step4_bundle$standardization_constants, stringsAsFactors = FALSE)
  c(
    age_mean_merged = as.numeric(df$age_mean_merged[1L]),
    age_sd_merged = as.numeric(df$age_sd_merged[1L]),
    sex_mean_merged = as.numeric(df$sex_mean_merged[1L])
  )
}

extract_step5_scaling_vector <- function(step5_bundle) {
  df <- as.data.frame(step5_bundle$standardization, stringsAsFactors = FALSE)
  c(
    age_mean_merged = as.numeric(df$age_mean_merged[1L]),
    age_sd_merged = as.numeric(df$age_sd_merged[1L]),
    sex_mean_merged = as.numeric(df$sex_mean_merged[1L])
  )
}

validate_step4_step5_alignment <- function(step4_bundle, step5_bundle, step5_data, tolerance = 1e-10) {
  step4_data <- as.data.frame(step4_bundle$data_complete_case, stringsAsFactors = FALSE)
  
  required_step4_cols <- c(
    "id",
    "site",
    "time_primary",
    "event_primary",
    "c_sex",
    "z_age"
  )
  missing_step4_cols <- setdiff(required_step4_cols, names(step4_data))
  if (length(missing_step4_cols) > 0L) {
    stop(
      "Step4 data_complete_case is missing required column(s): ",
      paste(missing_step4_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  step4_data$id <- as.character(step4_data$id)
  step4_data$site <- toupper(trimws(as.character(step4_data$site)))
  step4_data$time_primary <- as.numeric(step4_data$time_primary)
  step4_data$event_primary <- as.integer(as.numeric(step4_data$event_primary))
  step4_data$c_sex <- as.numeric(step4_data$c_sex)
  step4_data$z_age <- as.numeric(step4_data$z_age)
  
  if (!("sex_age_interaction" %in% names(step4_data))) {
    step4_data$sex_age_interaction <- step4_data$c_sex * step4_data$z_age
  } else {
    step4_data$sex_age_interaction <- as.numeric(step4_data$sex_age_interaction)
  }
  
  step4_data$subject_key <- make_subject_key(step4_data)
  
  if (anyDuplicated(step4_data$subject_key) > 0L) {
    duplicated_example <- step4_data$subject_key[duplicated(step4_data$subject_key)][1L]
    stop("Duplicated subject key detected in Step4 data_complete_case: ", duplicated_example, call. = FALSE)
  }
  
  subject_set_match <- setequal(step4_data$subject_key, step5_data$subject_key)
  
  endpoint_match <- FALSE
  transformed_covariate_match <- FALSE
  
  if (isTRUE(subject_set_match)) {
    step4_small <- step4_data[, c("subject_key", "time_primary", "event_primary", "c_sex", "z_age", "sex_age_interaction"), drop = FALSE]
    step5_small <- step5_data[, c("subject_key", "time_primary", "event_primary", "c_sex", "z_age", "int_sex_age"), drop = FALSE]
    
    merged_check <- merge(
      step4_small,
      step5_small,
      by = "subject_key",
      all = FALSE,
      sort = FALSE,
      suffixes = c("_step4", "_step5")
    )
    
    endpoint_match <- nrow(merged_check) == nrow(step5_small) &&
      all(abs(merged_check$time_primary_step4 - merged_check$time_primary_step5) <= tolerance) &&
      all(merged_check$event_primary_step4 == merged_check$event_primary_step5)
    
    transformed_covariate_match <- nrow(merged_check) == nrow(step5_small) &&
      all(abs(merged_check$c_sex_step4 - merged_check$c_sex_step5) <= tolerance) &&
      all(abs(merged_check$z_age_step4 - merged_check$z_age_step5) <= tolerance) &&
      all(abs(merged_check$sex_age_interaction - merged_check$int_sex_age) <= tolerance)
  }
  
  scaling_step4 <- extract_step4_scaling_vector(step4_bundle)
  scaling_step5 <- extract_step5_scaling_vector(step5_bundle)
  
  scaling_match <- all(is.finite(scaling_step4)) &&
    all(is.finite(scaling_step5)) &&
    all(abs(scaling_step4 - scaling_step5) <= tolerance)
  
  note <- blank_to_na(c(
    note_if_false(subject_set_match, "Step4 and Step5 subject sets differ."),
    note_if_false(endpoint_match, "Step4 and Step5 endpoint variables differ."),
    note_if_false(transformed_covariate_match, "Step4 and Step5 transformed covariates differ."),
    note_if_false(scaling_match, "Step4 and Step5 scaling constants differ.")
  ))
  
  list(
    subject_set_match = subject_set_match,
    endpoint_match = endpoint_match,
    transformed_covariate_match = transformed_covariate_match,
    scaling_match = scaling_match,
    note = note
  )
}

extract_step4_matched_record <- function(step4_bundle, pair_spec) {
  fit_nc_list <- step4_bundle$fit_nc_list %||% NULL
  if (is.null(fit_nc_list)) {
    return(NULL)
  }
  
  specification_key <- as.character(pair_spec$step4_specification[1L])
  family_key <- as.character(pair_spec$step4_family[1L])
  
  if (!specification_key %in% names(fit_nc_list)) {
    return(NULL)
  }
  
  if (!family_key %in% names(fit_nc_list[[specification_key]])) {
    return(NULL)
  }
  
  fit_nc_list[[specification_key]][[family_key]]
}

extract_step4_matched_status <- function(step4_record) {
  list(
    ok = if (!is.null(step4_record)) isTRUE(step4_record$ok) else NA,
    converged = if (!is.null(step4_record)) isTRUE(step4_record$converged) else NA
  )
}

# 🔴 Define: fit wrappers and likelihood summaries ===============================

## 🟠 Build: guarded fit runners, convergence checks, and evaluation rules ===============================
capture_fit <- function(expr_fun) {
  warnings_caught <- character(0)
  t0 <- proc.time()[["elapsed"]]
  
  fit_obj <- tryCatch(
    withCallingHandlers(
      expr_fun(),
      warning = function(w) {
        warnings_caught <<- c(warnings_caught, conditionMessage(w))
        tryInvokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  
  elapsed_sec <- proc.time()[["elapsed"]] - t0
  
  if (inherits(fit_obj, "error")) {
    return(list(
      fit = NULL,
      ok = FALSE,
      warnings = unique(warnings_caught),
      error = conditionMessage(fit_obj),
      elapsed_sec = elapsed_sec
    ))
  }
  
  list(
    fit = fit_obj,
    ok = TRUE,
    warnings = unique(warnings_caught),
    error = NA_character_,
    elapsed_sec = elapsed_sec
  )
}

get_survreg_converged <- function(fit) {
  if (is.null(fit)) {
    return(FALSE)
  }
  
  coef_ok <- tryCatch({
    cf <- stats::coef(fit)
    all(is.finite(as.numeric(cf)))
  }, error = function(e) FALSE)
  
  scale_ok <- isTRUE(is.finite(fit$scale)) && fit$scale > 0
  
  fail_ok <- is.null(fit$fail) ||
    identical(fit$fail, FALSE) ||
    identical(fit$fail, 0L) ||
    identical(fit$fail, 0)
  
  isTRUE(coef_ok) && isTRUE(scale_ok) && isTRUE(fail_ok)
}

get_flexsurvcure_converged <- function(fit) {
  if (is.null(fit)) {
    return(FALSE)
  }
  
  coef_ok <- tryCatch({
    cf <- stats::coef(fit)
    all(is.finite(as.numeric(cf)))
  }, error = function(e) FALSE)
  
  if (!coef_ok) {
    return(FALSE)
  }
  
  if (!is.null(fit$converged)) {
    return(isTRUE(as.logical(fit$converged)) && coef_ok)
  }
  
  if (!is.null(fit$optim$convergence)) {
    return((identical(fit$optim$convergence, 0L) || identical(fit$optim$convergence, 0)) && coef_ok)
  }
  
  TRUE
}

finalize_fit_info <- function(fit_info, model_family) {
  converged_flag <- FALSE
  
  if (identical(model_family, "null_survreg")) {
    converged_flag <- get_survreg_converged(fit_info$fit)
  } else if (identical(model_family, "alt_cure")) {
    converged_flag <- get_flexsurvcure_converged(fit_info$fit)
  } else {
    stop("Unknown model_family in finalize_fit_info: ", model_family, call. = FALSE)
  }
  
  warning_message <- blank_to_na(fit_info$warnings)
  error_message <- blank_to_na(fit_info$error)
  
  stable_flag <- isTRUE(fit_info$ok) &&
    isTRUE(converged_flag) &&
    !detect_instability_message(paste(c(warning_message, error_message), collapse = " | "))
  
  fit_info$converged <- converged_flag
  fit_info$stable <- stable_flag
  fit_info$warning_message <- warning_message
  fit_info$error_message <- error_message
  fit_info
}

safe_fit_metrics <- function(fit, n_obs) {
  loglik_obj <- tryCatch(stats::logLik(fit), error = function(e) NULL)
  
  logLik_value <- if (!is.null(loglik_obj)) as.numeric(loglik_obj) else NA_real_
  df_model <- if (!is.null(loglik_obj)) {
    as.numeric(attr(loglik_obj, "df"))
  } else {
    tryCatch(length(stats::coef(fit)), error = function(e) NA_real_)
  }
  
  aic_value <- tryCatch(as.numeric(stats::AIC(fit)), error = function(e) NA_real_)
  bic_value <- tryCatch(as.numeric(stats::BIC(fit)), error = function(e) NA_real_)
  
  if (!is.finite(aic_value) && is.finite(logLik_value) && is.finite(df_model)) {
    aic_value <- -2 * logLik_value + 2 * df_model
  }
  
  if (!is.finite(bic_value) && is.finite(logLik_value) && is.finite(df_model) && is.finite(n_obs)) {
    bic_value <- -2 * logLik_value + log(n_obs) * df_model
  }
  
  list(
    logLik = logLik_value,
    df_model = df_model,
    AIC = aic_value,
    BIC = bic_value
  )
}

fit_null_no_cure_model <- function(data, pair_spec, survreg_maxiter) {
  formula_obj <- make_surv_formula(pair_spec$covariate_structure[1L])
  
  fit_info <- capture_fit(function() {
    survival::survreg(
      formula = formula_obj,
      data = data,
      dist = pair_spec$null_survreg_dist[1L],
      control = survival::survreg.control(maxiter = survreg_maxiter),
      x = TRUE,
      y = TRUE,
      model = TRUE
    )
  })
  
  finalize_fit_info(fit_info, model_family = "null_survreg")
}

fit_alt_cure_model <- function(data, pair_spec, cure_optim_method, optim_maxit, survreg_maxiter) {
  formula_obj <- make_surv_formula(pair_spec$covariate_structure[1L])
  anc_list <- make_anc_list(
    latency_param = pair_spec$latency_param[1L],
    covariate_structure = pair_spec$covariate_structure[1L]
  )
  
  fit_info <- capture_fit(function() {
    flexsurvcure::flexsurvcure(
      formula = formula_obj,
      data = data,
      dist = pair_spec$latency_dist[1L],
      anc = anc_list,
      link = "logistic",
      mixture = TRUE,
      na.action = stats::na.omit,
      sr.control = survival::survreg.control(maxiter = survreg_maxiter),
      method = cure_optim_method,
      control = list(maxit = optim_maxit)
    )
  })
  
  finalize_fit_info(fit_info, model_family = "alt_cure")
}

is_pair_evaluable_mode <- function(null_fit_info, alt_fit_info, null_metrics, alt_metrics, require_stable = TRUE) {
  fit_part <- if (isTRUE(require_stable)) {
    isTRUE(null_fit_info$stable) && isTRUE(alt_fit_info$stable)
  } else {
    isTRUE(null_fit_info$ok) &&
      isTRUE(null_fit_info$converged) &&
      isTRUE(alt_fit_info$ok) &&
      isTRUE(alt_fit_info$converged)
  }
  
  fit_part &&
    is.finite(null_metrics$logLik) &&
    is.finite(alt_metrics$logLik)
}

compose_fit_failure_reason <- function(fit_info, require_stable = TRUE) {
  if (isFALSE(fit_info$ok)) {
    return(blank_to_na(c("fit_error", truncate_text(fit_info$error_message, 160L))))
  }
  
  if (!isTRUE(fit_info$converged)) {
    return(blank_to_na(c("nonconverged", truncate_text(fit_info$warning_message, 160L), truncate_text(fit_info$error_message, 160L))))
  }
  
  if (isTRUE(require_stable) && !isTRUE(fit_info$stable)) {
    return(blank_to_na(c("numerical_instability", truncate_text(fit_info$warning_message, 160L), truncate_text(fit_info$error_message, 160L))))
  }
  
  "unknown_fit_failure"
}

# 🔴 Define: censoring samplers and bootstrap engines ===============================

## 🟠 Build: reverse-KM censoring samplers and mode-specific bootstrap runners ===============================
build_censoring_sampler <- function(data, method) {
  method <- tolower(method[1L])
  censor_times_observed <- as.numeric(data$time_model[data$event_primary == 0L])
  
  if (identical(method, "empirical_censoring_resample")) {
    if (length(censor_times_observed) == 0L) {
      return(list(
        sample = function(n) rep(Inf, n),
        summary = data.frame(
          method_requested = method,
          method_used = "all_infinite_no_observed_censoring",
          n_observed_censor = 0L,
          unique_censor_times = 0L,
          total_finite_mass = 0,
          infinite_mass = 1,
          stringsAsFactors = FALSE
        )
      ))
    }
    
    return(list(
      sample = function(n) sample(censor_times_observed, size = n, replace = TRUE),
      summary = data.frame(
        method_requested = method,
        method_used = "empirical_censoring_resample",
        n_observed_censor = length(censor_times_observed),
        unique_censor_times = length(unique(censor_times_observed)),
        total_finite_mass = 1,
        infinite_mass = 0,
        stringsAsFactors = FALSE
      )
    ))
  }
  
  status_censor <- 1L - as.integer(data$event_primary)
  
  if (sum(status_censor, na.rm = TRUE) == 0L) {
    return(list(
      sample = function(n) rep(Inf, n),
      summary = data.frame(
        method_requested = method,
        method_used = "reverse_km_no_observed_censoring",
        n_observed_censor = 0L,
        unique_censor_times = 0L,
        total_finite_mass = 0,
        infinite_mass = 1,
        stringsAsFactors = FALSE
      )
    ))
  }
  
  sf <- survival::survfit(
    survival::Surv(time_model, status_censor) ~ 1,
    data = data
  )
  
  times <- as.numeric(sf$time)
  surv <- as.numeric(sf$surv)
  
  if (length(times) == 0L || length(surv) == 0L) {
    if (length(censor_times_observed) == 0L) {
      return(list(
        sample = function(n) rep(Inf, n),
        summary = data.frame(
          method_requested = method,
          method_used = "reverse_km_empty_all_infinite",
          n_observed_censor = 0L,
          unique_censor_times = 0L,
          total_finite_mass = 0,
          infinite_mass = 1,
          stringsAsFactors = FALSE
        )
      ))
    }
    
    return(list(
      sample = function(n) sample(censor_times_observed, size = n, replace = TRUE),
      summary = data.frame(
        method_requested = method,
        method_used = "reverse_km_fallback_empirical",
        n_observed_censor = length(censor_times_observed),
        unique_censor_times = length(unique(censor_times_observed)),
        total_finite_mass = 1,
        infinite_mass = 0,
        stringsAsFactors = FALSE
      )
    ))
  }
  
  s_prev <- c(1, head(surv, -1L))
  jump_prob <- pmax(0, s_prev - surv)
  
  if (sum(jump_prob) > 1 + 1e-10) {
    jump_prob <- jump_prob / sum(jump_prob)
  }
  
  total_finite_mass <- min(1, sum(jump_prob))
  infinite_mass <- max(0, 1 - total_finite_mass)
  cum_prob <- cumsum(jump_prob)
  
  sample_fun <- function(n) {
    u <- stats::runif(n)
    out <- rep(Inf, n)
    
    if (length(times) > 0L && total_finite_mass > 0) {
      finite_idx <- which(u <= total_finite_mass)
      if (length(finite_idx) > 0L) {
        j <- findInterval(u[finite_idx], cum_prob) + 1L
        j[j < 1L] <- 1L
        j[j > length(times)] <- length(times)
        out[finite_idx] <- times[j]
      }
    }
    
    out
  }
  
  list(
    sample = sample_fun,
    summary = data.frame(
      method_requested = method,
      method_used = "reverse_km",
      n_observed_censor = length(censor_times_observed),
      unique_censor_times = length(times),
      total_finite_mass = total_finite_mass,
      infinite_mass = infinite_mass,
      stringsAsFactors = FALSE
    )
  )
}

simulate_event_times_from_null <- function(null_fit, newdata, survreg_dist, zero_time_offset_days) {
  lp <- as.numeric(stats::predict(null_fit, newdata = newdata, type = "lp"))
  u <- stats::runif(length(lp))
  
  sim_time <- as.numeric(
    survival::qsurvreg(
      p = u,
      mean = lp,
      scale = null_fit$scale,
      distribution = survreg_dist
    )
  )
  
  if (length(sim_time) != nrow(newdata)) {
    stop("Simulated event-time vector length does not match nrow(newdata).", call. = FALSE)
  }
  
  if (any(!is.finite(sim_time))) {
    stop("Non-finite event times were generated from the null model.", call. = FALSE)
  }
  
  pmax(sim_time, zero_time_offset_days)
}

make_bootstrap_dataset <- function(data_template, null_fit, pair_spec, censor_sampler, zero_time_offset_days) {
  event_time <- simulate_event_times_from_null(
    null_fit = null_fit,
    newdata = data_template,
    survreg_dist = pair_spec$null_survreg_dist[1L],
    zero_time_offset_days = zero_time_offset_days
  )
  
  censor_time <- as.numeric(censor_sampler$sample(nrow(data_template)))
  
  if (length(censor_time) != nrow(data_template)) {
    stop("Censoring-time vector length does not match nrow(data_template).", call. = FALSE)
  }
  
  observed_time <- pmin(event_time, censor_time)
  event_observed <- as.integer(event_time <= censor_time | !is.finite(censor_time))
  
  observed_time <- pmax(as.numeric(observed_time), zero_time_offset_days)
  
  dat_b <- data_template
  dat_b$time_model <- observed_time
  dat_b$time_primary <- observed_time
  dat_b$event_primary <- event_observed
  dat_b$status_primary <- ifelse(dat_b$event_primary == 1L, "transition", "censored")
  dat_b$time_model_offset_applied <- as.integer(dat_b$time_model <= zero_time_offset_days + sqrt(.Machine$double.eps))
  
  dat_b
}

run_boundary_bootstrap <- function(
    data_template,
    pair_spec,
    null_fit_observed,
    censor_sampler,
    zero_time_offset_days,
    bootstrap_target_B,
    bootstrap_max_attempt_multiplier,
    cure_optim_method,
    optim_maxit,
    survreg_maxiter,
    failure_log_max_rows,
    require_stable = TRUE
) {
  d_boot <- numeric(0)
  failure_log <- new_empty_failure_log()
  failure_counter <- integer(0)
  
  accepted_unstable_null_B <- 0L
  accepted_unstable_alt_B <- 0L
  accepted_any_unstable_B <- 0L
  accepted_all_stable_B <- 0L
  
  max_attempts <- max(1L, as.integer(bootstrap_target_B) * as.integer(bootstrap_max_attempt_multiplier))
  attempt <- 0L
  
  add_failure <- function(stage, reason) {
    key <- paste(stage, reason, sep = "||")
    
    if (key %in% names(failure_counter)) {
      failure_counter[key] <<- failure_counter[key] + 1L
    } else {
      failure_counter[key] <<- 1L
    }
    
    if (nrow(failure_log) < failure_log_max_rows) {
      failure_log <<- rbind(
        failure_log,
        data.frame(
          attempt = as.integer(attempt),
          stage = as.character(stage),
          reason = as.character(reason),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  while (length(d_boot) < bootstrap_target_B && attempt < max_attempts) {
    attempt <- attempt + 1L
    
    dat_b <- tryCatch(
      make_bootstrap_dataset(
        data_template = data_template,
        null_fit = null_fit_observed,
        pair_spec = pair_spec,
        censor_sampler = censor_sampler,
        zero_time_offset_days = zero_time_offset_days
      ),
      error = function(e) e
    )
    
    if (inherits(dat_b, "error")) {
      add_failure("data_generation", truncate_text(conditionMessage(dat_b), 160L))
      next
    }
    
    if (sum(dat_b$event_primary == 1L, na.rm = TRUE) == 0L) {
      add_failure("data_generation", "no_events_in_bootstrap_sample")
      next
    }
    
    null_fit_b <- fit_null_no_cure_model(
      data = dat_b,
      pair_spec = pair_spec,
      survreg_maxiter = survreg_maxiter
    )
    
    null_metrics_b <- if (!is.null(null_fit_b$fit)) {
      safe_fit_metrics(null_fit_b$fit, n_obs = nrow(dat_b))
    } else {
      list(logLik = NA_real_, df_model = NA_real_, AIC = NA_real_, BIC = NA_real_)
    }
    
    alt_fit_b <- fit_alt_cure_model(
      data = dat_b,
      pair_spec = pair_spec,
      cure_optim_method = cure_optim_method,
      optim_maxit = optim_maxit,
      survreg_maxiter = survreg_maxiter
    )
    
    alt_metrics_b <- if (!is.null(alt_fit_b$fit)) {
      safe_fit_metrics(alt_fit_b$fit, n_obs = nrow(dat_b))
    } else {
      list(logLik = NA_real_, df_model = NA_real_, AIC = NA_real_, BIC = NA_real_)
    }
    
    accepted_flag <- is_pair_evaluable_mode(
      null_fit_info = null_fit_b,
      alt_fit_info = alt_fit_b,
      null_metrics = null_metrics_b,
      alt_metrics = alt_metrics_b,
      require_stable = require_stable
    )
    
    if (!isTRUE(accepted_flag)) {
      if (!(isTRUE(null_fit_b$ok) && isTRUE(null_fit_b$converged) && (!isTRUE(require_stable) || isTRUE(null_fit_b$stable)))) {
        add_failure("null_refit", compose_fit_failure_reason(null_fit_b, require_stable = require_stable))
        next
      }
      
      if (!(isTRUE(alt_fit_b$ok) && isTRUE(alt_fit_b$converged) && (!isTRUE(require_stable) || isTRUE(alt_fit_b$stable)))) {
        add_failure("alt_refit", compose_fit_failure_reason(alt_fit_b, require_stable = require_stable))
        next
      }
      
      if (!is.finite(null_metrics_b$logLik) || !is.finite(alt_metrics_b$logLik)) {
        add_failure("logLik", "non_finite_loglik")
        next
      }
      
      add_failure("other", "not_accepted_for_mode")
      next
    }
    
    unstable_null_flag <- !isTRUE(null_fit_b$stable)
    unstable_alt_flag <- !isTRUE(alt_fit_b$stable)
    
    accepted_unstable_null_B <- accepted_unstable_null_B + as.integer(unstable_null_flag)
    accepted_unstable_alt_B <- accepted_unstable_alt_B + as.integer(unstable_alt_flag)
    accepted_any_unstable_B <- accepted_any_unstable_B + as.integer(unstable_null_flag || unstable_alt_flag)
    accepted_all_stable_B <- accepted_all_stable_B + as.integer(!(unstable_null_flag || unstable_alt_flag))
    
    d_boot <- c(d_boot, max(0, 2 * (alt_metrics_b$logLik - null_metrics_b$logLik)))
  }
  
  failure_counts <- if (length(failure_counter) == 0L) {
    new_empty_failure_counts()
  } else {
    split_names <- strsplit(names(failure_counter), "||", fixed = TRUE)
    out <- data.frame(
      stage = vapply(split_names, `[`, character(1), 1L),
      reason = vapply(split_names, `[`, character(1), 2L),
      n = as.integer(failure_counter),
      stringsAsFactors = FALSE
    )
    out[order(-out$n, out$stage, out$reason), , drop = FALSE]
  }
  
  out <- new_bootstrap_result(target_B = bootstrap_target_B)
  out$d_boot <- d_boot
  out$attempts <- attempt
  out$success_B <- length(d_boot)
  out$failed_B <- if (length(failure_counter) == 0L) 0L else sum(failure_counter)
  out$target_met <- length(d_boot) >= bootstrap_target_B
  out$failure_log <- failure_log
  out$failure_counts <- failure_counts
  out$failure_summary <- blank_to_na(
    if (nrow(failure_counts) > 0L) {
      paste0(failure_counts$stage, ":", failure_counts$reason, "=", failure_counts$n)
    } else {
      NA_character_
    }
  )
  out$accepted_unstable_null_B <- accepted_unstable_null_B
  out$accepted_unstable_alt_B <- accepted_unstable_alt_B
  out$accepted_any_unstable_B <- accepted_any_unstable_B
  out$accepted_all_stable_B <- accepted_all_stable_B
  out
}

# 🔴 Define: pair evaluators and decision rules ===============================

## 🟠 Build: strict and relaxed decision rules, pair evaluators, and placeholders ===============================
classify_strict_reading_role <- function(gate_status, evaluable_flag) {
  if (!isTRUE(evaluable_flag)) {
    return("not_evaluable")
  }
  
  if (identical(gate_status, "pass")) {
    return("confirmatory")
  }
  
  if (identical(gate_status, "borderline")) {
    return("supportive")
  }
  
  if (identical(gate_status, "fail")) {
    return("exploratory")
  }
  
  "not_evaluable"
}

classify_relaxed_reading_role <- function(evaluable_flag) {
  if (!isTRUE(evaluable_flag)) {
    return("not_evaluable")
  }
  
  "exploratory_relaxed"
}

classify_decision <- function(p_value, alpha, evaluable_flag) {
  if (!isTRUE(evaluable_flag) || !is.finite(p_value)) {
    return("not_evaluable")
  }
  
  if (p_value < alpha) {
    return("reject_H0_support_nonzero_cure")
  }
  
  "do_not_reject_H0"
}

evaluate_boundary_pair <- function(
    data_template,
    dataset_name,
    gate_status,
    step5_gate_pass,
    step5_selected_cure_id,
    pair_spec,
    is_primary_pair,
    dataset_validation,
    step4_record,
    alpha,
    strict_bootstrap_target_B,
    relaxed_bootstrap_target_B,
    bootstrap_max_attempt_multiplier,
    censoring_bootstrap_method,
    cure_optim_method,
    optim_maxit,
    survreg_maxiter,
    zero_time_offset_days,
    random_seed_base,
    failure_log_max_rows
) {
  analysis_n <- nrow(data_template)
  n_event <- sum(data_template$event_primary == 1L, na.rm = TRUE)
  n_censor <- sum(data_template$event_primary == 0L, na.rm = TRUE)
  
  step4_status <- extract_step4_matched_status(step4_record)
  
  strict_bootstrap_res <- new_bootstrap_result(strict_bootstrap_target_B)
  relaxed_bootstrap_res <- new_bootstrap_result(relaxed_bootstrap_target_B)
  
  strict_pair_seed <- NA_integer_
  relaxed_pair_seed <- NA_integer_
  censor_sampler_summary <- data.frame(
    method_requested = censoring_bootstrap_method,
    method_used = NA_character_,
    n_observed_censor = n_censor,
    unique_censor_times = NA_integer_,
    total_finite_mass = NA_real_,
    infinite_mass = NA_real_,
    stringsAsFactors = FALSE
  )
  
  null_fit_info <- list(
    fit = NULL,
    ok = FALSE,
    converged = FALSE,
    stable = FALSE,
    warning_message = NA_character_,
    error_message = NA_character_
  )
  
  alt_fit_info <- list(
    fit = NULL,
    ok = FALSE,
    converged = FALSE,
    stable = FALSE,
    warning_message = NA_character_,
    error_message = NA_character_
  )
  
  null_metrics <- list(logLik = NA_real_, df_model = NA_real_, AIC = NA_real_, BIC = NA_real_)
  alt_metrics <- list(logLik = NA_real_, df_model = NA_real_, AIC = NA_real_, BIC = NA_real_)
  
  strict_observed_evaluable <- FALSE
  relaxed_observed_evaluable <- FALSE
  d_obs <- NA_real_
  
  if (analysis_n > 0L && n_event > 0L) {
    null_fit_info <- fit_null_no_cure_model(
      data = data_template,
      pair_spec = pair_spec,
      survreg_maxiter = survreg_maxiter
    )
    
    alt_fit_info <- fit_alt_cure_model(
      data = data_template,
      pair_spec = pair_spec,
      cure_optim_method = cure_optim_method,
      optim_maxit = optim_maxit,
      survreg_maxiter = survreg_maxiter
    )
    
    null_metrics <- if (!is.null(null_fit_info$fit)) {
      safe_fit_metrics(null_fit_info$fit, n_obs = analysis_n)
    } else {
      list(logLik = NA_real_, df_model = NA_real_, AIC = NA_real_, BIC = NA_real_)
    }
    
    alt_metrics <- if (!is.null(alt_fit_info$fit)) {
      safe_fit_metrics(alt_fit_info$fit, n_obs = analysis_n)
    } else {
      list(logLik = NA_real_, df_model = NA_real_, AIC = NA_real_, BIC = NA_real_)
    }
    
    strict_observed_evaluable <- is_pair_evaluable_mode(
      null_fit_info = null_fit_info,
      alt_fit_info = alt_fit_info,
      null_metrics = null_metrics,
      alt_metrics = alt_metrics,
      require_stable = TRUE
    )
    
    relaxed_observed_evaluable <- is_pair_evaluable_mode(
      null_fit_info = null_fit_info,
      alt_fit_info = alt_fit_info,
      null_metrics = null_metrics,
      alt_metrics = alt_metrics,
      require_stable = FALSE
    )
    
    if (isTRUE(relaxed_observed_evaluable)) {
      d_obs <- max(0, 2 * (alt_metrics$logLik - null_metrics$logLik))
    }
    
    if (isTRUE(strict_observed_evaluable) || isTRUE(relaxed_observed_evaluable)) {
      censor_sampler <- build_censoring_sampler(
        data = data_template,
        method = censoring_bootstrap_method
      )
      censor_sampler_summary <- censor_sampler$summary
      
      if (isTRUE(strict_observed_evaluable)) {
        strict_pair_seed <- make_pair_seed(
          base_seed = random_seed_base,
          dataset_name = dataset_name,
          alt_model_id = pair_spec$alt_model_id[1L],
          mode_label = "strict"
        )
        
        set.seed(strict_pair_seed)
        
        strict_bootstrap_res <- run_boundary_bootstrap(
          data_template = data_template,
          pair_spec = pair_spec,
          null_fit_observed = null_fit_info$fit,
          censor_sampler = censor_sampler,
          zero_time_offset_days = zero_time_offset_days,
          bootstrap_target_B = strict_bootstrap_target_B,
          bootstrap_max_attempt_multiplier = bootstrap_max_attempt_multiplier,
          cure_optim_method = cure_optim_method,
          optim_maxit = optim_maxit,
          survreg_maxiter = survreg_maxiter,
          failure_log_max_rows = failure_log_max_rows,
          require_stable = TRUE
        )
      }
      
      if (isTRUE(relaxed_observed_evaluable)) {
        relaxed_pair_seed <- make_pair_seed(
          base_seed = random_seed_base,
          dataset_name = dataset_name,
          alt_model_id = pair_spec$alt_model_id[1L],
          mode_label = "relaxed"
        )
        
        set.seed(relaxed_pair_seed)
        
        relaxed_bootstrap_res <- run_boundary_bootstrap(
          data_template = data_template,
          pair_spec = pair_spec,
          null_fit_observed = null_fit_info$fit,
          censor_sampler = censor_sampler,
          zero_time_offset_days = zero_time_offset_days,
          bootstrap_target_B = relaxed_bootstrap_target_B,
          bootstrap_max_attempt_multiplier = bootstrap_max_attempt_multiplier,
          cure_optim_method = cure_optim_method,
          optim_maxit = optim_maxit,
          survreg_maxiter = survreg_maxiter,
          failure_log_max_rows = failure_log_max_rows,
          require_stable = FALSE
        )
      }
    }
  }
  
  strict_p_value_boundary <- if (isTRUE(strict_observed_evaluable) && strict_bootstrap_res$success_B > 0L) {
    (1 + sum(strict_bootstrap_res$d_boot >= d_obs)) / (strict_bootstrap_res$success_B + 1)
  } else {
    NA_real_
  }
  
  relaxed_p_value_boundary <- if (isTRUE(relaxed_observed_evaluable) && relaxed_bootstrap_res$success_B > 0L) {
    (1 + sum(relaxed_bootstrap_res$d_boot >= d_obs)) / (relaxed_bootstrap_res$success_B + 1)
  } else {
    NA_real_
  }
  
  strict_summary <- summarize_bootstrap_result(strict_bootstrap_res)
  relaxed_summary <- summarize_bootstrap_result(relaxed_bootstrap_res)
  
  strict_evaluable_final <- isTRUE(strict_observed_evaluable) && isTRUE(strict_bootstrap_res$target_met)
  relaxed_evaluable_final <- isTRUE(relaxed_observed_evaluable) && isTRUE(relaxed_bootstrap_res$target_met)
  
  strict_reading_role <- classify_strict_reading_role(
    gate_status = gate_status,
    evaluable_flag = strict_evaluable_final
  )
  
  relaxed_reading_role <- classify_relaxed_reading_role(
    evaluable_flag = relaxed_evaluable_final
  )
  
  strict_decision <- classify_decision(
    p_value = strict_p_value_boundary,
    alpha = alpha,
    evaluable_flag = strict_evaluable_final
  )
  
  relaxed_decision <- classify_decision(
    p_value = relaxed_p_value_boundary,
    alpha = alpha,
    evaluable_flag = relaxed_evaluable_final
  )
  
  note <- blank_to_na(c(
    note_if_false(dataset_validation$subject_set_match, "Step4 and Step5 subject sets differ."),
    note_if_false(dataset_validation$endpoint_match, "Step4 and Step5 endpoint variables differ."),
    note_if_false(dataset_validation$transformed_covariate_match, "Step4 and Step5 transformed covariates differ."),
    note_if_false(dataset_validation$scaling_match, "Step4 and Step5 scaling constants differ."),
    if (analysis_n == 0L) "No rows available in the common Step6 analysis frame." else NA_character_,
    if (analysis_n > 0L && n_event == 0L) "No transition events are available in the common Step6 analysis frame." else NA_character_,
    if (!isTRUE(relaxed_observed_evaluable) && analysis_n > 0L && n_event > 0L) {
      "Observed matched-pair refit was not evaluable even under the relaxed rule."
    } else {
      NA_character_
    },
    if (!isTRUE(strict_observed_evaluable) && isTRUE(relaxed_observed_evaluable)) {
      "Observed matched-pair refit failed the strict stability rule but passed the relaxed exploratory rule."
    } else {
      NA_character_
    },
    if (isTRUE(relaxed_observed_evaluable) && !isTRUE(alt_fit_info$stable)) {
      "Observed alternative cure fit had a stability warning; relaxed results must be interpreted as exploratory."
    } else {
      NA_character_
    },
    if (isTRUE(strict_observed_evaluable) && !isTRUE(strict_bootstrap_res$target_met)) {
      paste0(
        "Strict bootstrap success target was not reached (",
        strict_bootstrap_res$success_B,
        "/",
        strict_bootstrap_target_B,
        ")."
      )
    } else {
      NA_character_
    },
    if (isTRUE(relaxed_observed_evaluable) && !isTRUE(relaxed_bootstrap_res$target_met)) {
      paste0(
        "Relaxed bootstrap success target was not reached (",
        relaxed_bootstrap_res$success_B,
        "/",
        relaxed_bootstrap_target_B,
        ")."
      )
    } else {
      NA_character_
    },
    if (isTRUE(relaxed_evaluable_final) && relaxed_bootstrap_res$accepted_any_unstable_B > 0L) {
      paste0(
        "Relaxed accepted bootstrap replicates included unstable fits (any unstable accepted = ",
        relaxed_bootstrap_res$accepted_any_unstable_B,
        ")."
      )
    } else {
      NA_character_
    }
  ))
  
  row <- build_pair_result_row(
    values = list(
      dataset = dataset_name,
      gate_status_step3 = gate_status,
      step5_gate_pass = step5_gate_pass,
      step5_selected_best_cure_id = step5_selected_cure_id,
      is_primary_pair = is_primary_pair,
      alt_model_id = as.character(pair_spec$alt_model_id[1L]),
      null_model_id = as.character(pair_spec$null_model_id[1L]),
      latency_dist = as.character(pair_spec$latency_dist[1L]),
      null_survreg_dist = as.character(pair_spec$null_survreg_dist[1L]),
      step4_family = as.character(pair_spec$step4_family[1L]),
      covariate_structure = as.character(pair_spec$covariate_structure[1L]),
      analysis_n = as.integer(analysis_n),
      n_event = as.integer(n_event),
      n_censor = as.integer(n_censor),
      zero_time_offset_days = as.numeric(zero_time_offset_days),
      bootstrap_method = as.character(censoring_bootstrap_method),
      null_fit_ok = null_fit_info$ok,
      null_converged = null_fit_info$converged,
      null_stable = null_fit_info$stable,
      alt_fit_ok = alt_fit_info$ok,
      alt_converged = alt_fit_info$converged,
      alt_stable = alt_fit_info$stable,
      step4_matched_ok = step4_status$ok,
      step4_matched_converged = step4_status$converged,
      step4_step5_subject_set_match = dataset_validation$subject_set_match,
      step4_step5_endpoint_match = dataset_validation$endpoint_match,
      step4_step5_transformed_covariate_match = dataset_validation$transformed_covariate_match,
      step4_step5_scaling_match = dataset_validation$scaling_match,
      logLik_null = null_metrics$logLik,
      logLik_alt = alt_metrics$logLik,
      df_model_null = null_metrics$df_model,
      df_model_alt = alt_metrics$df_model,
      AIC_null = null_metrics$AIC,
      AIC_alt = alt_metrics$AIC,
      BIC_null = null_metrics$BIC,
      BIC_alt = alt_metrics$BIC,
      d_obs = d_obs,
      strict_observed_evaluable = strict_observed_evaluable,
      strict_bootstrap_target_B = strict_bootstrap_res$target_B,
      strict_bootstrap_success_B = strict_bootstrap_res$success_B,
      strict_bootstrap_failed_B = strict_bootstrap_res$failed_B,
      strict_bootstrap_attempts = strict_bootstrap_res$attempts,
      strict_bootstrap_target_met = strict_bootstrap_res$target_met,
      strict_bootstrap_failure_rate = strict_summary$failure_rate,
      strict_bootstrap_null_mean = strict_summary$null_mean,
      strict_bootstrap_null_median = strict_summary$null_median,
      strict_bootstrap_null_q95 = strict_summary$null_q95,
      strict_bootstrap_null_q975 = strict_summary$null_q975,
      strict_p_value_boundary = strict_p_value_boundary,
      strict_alpha = alpha,
      strict_decision = strict_decision,
      strict_reading_role = strict_reading_role,
      strict_bootstrap_failure_summary = strict_bootstrap_res$failure_summary,
      relaxed_observed_evaluable = relaxed_observed_evaluable,
      relaxed_bootstrap_target_B = relaxed_bootstrap_res$target_B,
      relaxed_bootstrap_success_B = relaxed_bootstrap_res$success_B,
      relaxed_bootstrap_failed_B = relaxed_bootstrap_res$failed_B,
      relaxed_bootstrap_attempts = relaxed_bootstrap_res$attempts,
      relaxed_bootstrap_target_met = relaxed_bootstrap_res$target_met,
      relaxed_bootstrap_failure_rate = relaxed_summary$failure_rate,
      relaxed_bootstrap_null_mean = relaxed_summary$null_mean,
      relaxed_bootstrap_null_median = relaxed_summary$null_median,
      relaxed_bootstrap_null_q95 = relaxed_summary$null_q95,
      relaxed_bootstrap_null_q975 = relaxed_summary$null_q975,
      relaxed_p_value_boundary = relaxed_p_value_boundary,
      relaxed_alpha = alpha,
      relaxed_decision = relaxed_decision,
      relaxed_reading_role = relaxed_reading_role,
      relaxed_bootstrap_failure_summary = relaxed_bootstrap_res$failure_summary,
      relaxed_bootstrap_accepted_unstable_null_B = relaxed_bootstrap_res$accepted_unstable_null_B,
      relaxed_bootstrap_accepted_unstable_alt_B = relaxed_bootstrap_res$accepted_unstable_alt_B,
      relaxed_bootstrap_accepted_any_unstable_B = relaxed_bootstrap_res$accepted_any_unstable_B,
      relaxed_bootstrap_accepted_all_stable_B = relaxed_bootstrap_res$accepted_all_stable_B,
      null_warning_message = null_fit_info$warning_message,
      alt_warning_message = alt_fit_info$warning_message,
      null_error_message = null_fit_info$error_message,
      alt_error_message = alt_fit_info$error_message,
      note = note
    )
  )
  
  list(
    row = row,
    store = list(
      pair_spec = pair_spec,
      step4_matched_record = step4_record,
      observed = list(
        null_fit = null_fit_info$fit,
        alt_fit = alt_fit_info$fit,
        null_fit_info = null_fit_info,
        alt_fit_info = alt_fit_info,
        null_metrics = null_metrics,
        alt_metrics = alt_metrics,
        d_obs = d_obs,
        strict_observed_evaluable = strict_observed_evaluable,
        relaxed_observed_evaluable = relaxed_observed_evaluable
      ),
      censoring_sampler_summary = censor_sampler_summary,
      strict = list(
        pair_seed = strict_pair_seed,
        bootstrap = strict_bootstrap_res
      ),
      relaxed = list(
        pair_seed = relaxed_pair_seed,
        bootstrap = relaxed_bootstrap_res
      )
    )
  )
}

make_primary_placeholder_row <- function(
    dataset_name,
    gate_status,
    step5_gate_pass,
    step5_selected_cure_id,
    analysis_n,
    n_event,
    n_censor,
    zero_time_offset_days,
    bootstrap_method,
    alpha,
    note
) {
  pair_catalog_local <- make_pair_catalog()
  pair_spec <- if (!is.na(step5_selected_cure_id) && step5_selected_cure_id %in% pair_catalog_local$alt_model_id) {
    pair_catalog_local[pair_catalog_local$alt_model_id == step5_selected_cure_id, , drop = FALSE]
  } else {
    NULL
  }
  
  build_pair_result_row(
    values = list(
      dataset = dataset_name,
      gate_status_step3 = gate_status,
      step5_gate_pass = step5_gate_pass,
      step5_selected_best_cure_id = step5_selected_cure_id,
      is_primary_pair = TRUE,
      alt_model_id = if (!is.null(pair_spec)) as.character(pair_spec$alt_model_id[1L]) else NA_character_,
      null_model_id = if (!is.null(pair_spec)) as.character(pair_spec$null_model_id[1L]) else NA_character_,
      latency_dist = if (!is.null(pair_spec)) as.character(pair_spec$latency_dist[1L]) else NA_character_,
      null_survreg_dist = if (!is.null(pair_spec)) as.character(pair_spec$null_survreg_dist[1L]) else NA_character_,
      step4_family = if (!is.null(pair_spec)) as.character(pair_spec$step4_family[1L]) else NA_character_,
      covariate_structure = if (!is.null(pair_spec)) as.character(pair_spec$covariate_structure[1L]) else NA_character_,
      analysis_n = as.integer(analysis_n),
      n_event = as.integer(n_event),
      n_censor = as.integer(n_censor),
      zero_time_offset_days = as.numeric(zero_time_offset_days),
      bootstrap_method = as.character(bootstrap_method),
      strict_observed_evaluable = FALSE,
      strict_bootstrap_target_B = strict_bootstrap_success_target_primary,
      strict_bootstrap_success_B = 0L,
      strict_bootstrap_failed_B = 0L,
      strict_bootstrap_attempts = 0L,
      strict_bootstrap_target_met = FALSE,
      strict_alpha = alpha,
      strict_decision = "not_evaluable",
      strict_reading_role = "not_evaluable",
      relaxed_observed_evaluable = FALSE,
      relaxed_bootstrap_target_B = relaxed_bootstrap_success_target_primary,
      relaxed_bootstrap_success_B = 0L,
      relaxed_bootstrap_failed_B = 0L,
      relaxed_bootstrap_attempts = 0L,
      relaxed_bootstrap_target_met = FALSE,
      relaxed_alpha = alpha,
      relaxed_decision = "not_evaluable",
      relaxed_reading_role = "not_evaluable",
      note = note
    )
  )
}

# 🔴 Read: dataset-level Step4 and Step5 bundles ===============================

## 🟠 Import: aligned bundles and derive dataset-level prerequisites ===============================
pair_catalog <- make_pair_catalog()

step4_bundles <- list()
step5_bundles <- list()
prepared_analysis <- list()
dataset_validation_list <- list()

for (ds in dataset_labels) {
  step4_bundles[[ds]] <- readRDS(step4_bundle_paths[[ds]])
  step5_bundles[[ds]] <- readRDS(step5_bundle_paths[[ds]])
  
  assert_step4_bundle_structure(step4_bundles[[ds]], dataset_name = ds)
  assert_step5_bundle_structure(step5_bundles[[ds]], dataset_name = ds)
  
  prepared_analysis[[ds]] <- prepare_step5_analysis_frame(
    step5_bundle = step5_bundles[[ds]],
    zero_time_offset_fallback_days = zero_time_offset_fallback_days
  )
  
  dataset_validation_list[[ds]] <- validate_step4_step5_alignment(
    step4_bundle = step4_bundles[[ds]],
    step5_bundle = step5_bundles[[ds]],
    step5_data = prepared_analysis[[ds]]$data,
    tolerance = alignment_tolerance
  )
  
  if (isTRUE(strict_bundle_alignment)) {
    failed_checks <- c(
      if (isFALSE(dataset_validation_list[[ds]]$subject_set_match)) "subject_set_match" else NA_character_,
      if (isFALSE(dataset_validation_list[[ds]]$endpoint_match)) "endpoint_match" else NA_character_,
      if (isFALSE(dataset_validation_list[[ds]]$transformed_covariate_match)) "transformed_covariate_match" else NA_character_,
      if (isFALSE(dataset_validation_list[[ds]]$scaling_match)) "scaling_match" else NA_character_
    )
    failed_checks <- failed_checks[!is.na(failed_checks)]
    
    if (length(failed_checks) > 0L) {
      stop(
        "Strict Step4-Step5 alignment failed for dataset ",
        ds,
        ". Failed checks: ",
        paste(failed_checks, collapse = ", "),
        call. = FALSE
      )
    }
  }
}

# 🔴 Execute: Step6 boundary tests by dataset ===============================

## 🟠 Loop: strict confirmatory and relaxed exploratory bootstrap calibration ===============================
all_pairs_results_list <- list()
primary_results_list <- list()
manifest_rows <- list()
dataset_bundle_objects <- list()

for (ds in dataset_labels) {
  vmessage("Running Step6 boundary test for dataset: ", ds)
  
  step4_bundle <- step4_bundles[[ds]]
  step5_bundle <- step5_bundles[[ds]]
  data_i <- prepared_analysis[[ds]]$data
  zero_time_offset_days_i <- prepared_analysis[[ds]]$zero_time_offset_days
  
  gate_status_i <- normalize_gate_status(step3_gate_status[[ds]])
  step5_gate_pass_i <- isTRUE(step5_bundle$gate$pass)
  step5_gate_reason_i <- blank_to_na(step5_bundle$gate$reason %||% NA_character_)
  step5_selected_cure_id_i <- extract_step5_selected_cure_id(step5_bundle)
  step5_selection_rule_i <- extract_step5_selection_rule(step5_bundle)
  
  analysis_n_i <- nrow(data_i)
  n_event_i <- sum(data_i$event_primary == 1L, na.rm = TRUE)
  n_censor_i <- sum(data_i$event_primary == 0L, na.rm = TRUE)
  
  pair_queue_i <- select_pair_queue(
    pair_catalog = pair_catalog,
    primary_alt_id = step5_selected_cure_id_i,
    run_all_pairs = run_all_pairs
  )
  
  pair_rows_i <- list()
  pair_store_i <- list()
  
  if (nrow(pair_queue_i) > 0L) {
    for (i in seq_len(nrow(pair_queue_i))) {
      pair_spec_i <- pair_queue_i[i, , drop = FALSE]
      alt_model_id_i <- pair_spec_i$alt_model_id[1L]
      is_primary_pair_i <- !is.na(step5_selected_cure_id_i) && identical(alt_model_id_i, step5_selected_cure_id_i)
      
      strict_bootstrap_target_B_i <- if (isTRUE(is_primary_pair_i)) {
        strict_bootstrap_success_target_primary
      } else {
        strict_bootstrap_success_target_sensitivity
      }
      
      relaxed_bootstrap_target_B_i <- if (isTRUE(is_primary_pair_i)) {
        relaxed_bootstrap_success_target_primary
      } else {
        relaxed_bootstrap_success_target_sensitivity
      }
      
      vmessage(
        "  - Pair ",
        alt_model_id_i,
        " (primary = ",
        is_primary_pair_i,
        ", strict target B = ",
        strict_bootstrap_target_B_i,
        ", relaxed target B = ",
        relaxed_bootstrap_target_B_i,
        ")"
      )
      
      step4_record_i <- extract_step4_matched_record(
        step4_bundle = step4_bundle,
        pair_spec = pair_spec_i
      )
      
      eval_out_i <- evaluate_boundary_pair(
        data_template = data_i,
        dataset_name = ds,
        gate_status = gate_status_i,
        step5_gate_pass = step5_gate_pass_i,
        step5_selected_cure_id = step5_selected_cure_id_i,
        pair_spec = pair_spec_i,
        is_primary_pair = is_primary_pair_i,
        dataset_validation = dataset_validation_list[[ds]],
        step4_record = step4_record_i,
        alpha = boundary_alpha,
        strict_bootstrap_target_B = strict_bootstrap_target_B_i,
        relaxed_bootstrap_target_B = relaxed_bootstrap_target_B_i,
        bootstrap_max_attempt_multiplier = bootstrap_max_attempt_multiplier,
        censoring_bootstrap_method = censoring_bootstrap_method,
        cure_optim_method = cure_optim_method,
        optim_maxit = optim_maxit,
        survreg_maxiter = survreg_maxiter,
        zero_time_offset_days = zero_time_offset_days_i,
        random_seed_base = random_seed,
        failure_log_max_rows = failure_log_max_rows
      )
      
      pair_rows_i[[alt_model_id_i]] <- eval_out_i$row
      pair_store_i[[alt_model_id_i]] <- eval_out_i$store
      
      vmessage(
        "    d_obs = ",
        ifelse(is.finite(eval_out_i$row$d_obs[1L]), format(eval_out_i$row$d_obs[1L], digits = 6), "NA"),
        " | strict success = ",
        eval_out_i$row$strict_bootstrap_success_B[1L],
        "/",
        eval_out_i$row$strict_bootstrap_target_B[1L],
        ", p = ",
        ifelse(is.finite(eval_out_i$row$strict_p_value_boundary[1L]), format(eval_out_i$row$strict_p_value_boundary[1L], digits = 6), "NA"),
        " | relaxed success = ",
        eval_out_i$row$relaxed_bootstrap_success_B[1L],
        "/",
        eval_out_i$row$relaxed_bootstrap_target_B[1L],
        ", p = ",
        ifelse(is.finite(eval_out_i$row$relaxed_p_value_boundary[1L]), format(eval_out_i$row$relaxed_p_value_boundary[1L], digits = 6), "NA")
      )
    }
  }
  
  all_pairs_df_i <- if (length(pair_rows_i) > 0L) {
    bind_rows_fill(pair_rows_i)
  } else {
    new_empty_pair_result()
  }
  
  primary_df_i <- all_pairs_df_i[all_pairs_df_i$is_primary_pair %in% TRUE, , drop = FALSE]
  
  if (nrow(primary_df_i) == 0L) {
    primary_note_i <- if (!is.na(step5_selected_cure_id_i)) {
      paste0(
        "Primary Step5 best_cure_id was not evaluable in Step6: ",
        step5_selected_cure_id_i,
        "."
      )
    } else {
      "No Step5 best_cure_id was available for primary Step6 evaluation."
    }
    
    primary_df_i <- make_primary_placeholder_row(
      dataset_name = ds,
      gate_status = gate_status_i,
      step5_gate_pass = step5_gate_pass_i,
      step5_selected_cure_id = step5_selected_cure_id_i,
      analysis_n = analysis_n_i,
      n_event = n_event_i,
      n_censor = n_censor_i,
      zero_time_offset_days = zero_time_offset_days_i,
      bootstrap_method = censoring_bootstrap_method,
      alpha = boundary_alpha,
      note = primary_note_i
    )
  }
  
  primary_row_i <- primary_df_i[1L, , drop = FALSE]
  primary_pair_path_i <- dataset_rds_paths[[ds]]
  
  manifest_rows[[ds]] <- data.frame(
    dataset = ds,
    step3_gate_status = gate_status_i,
    step5_gate_pass = step5_gate_pass_i,
    step5_gate_reason = step5_gate_reason_i,
    step5_selected_best_cure_id = step5_selected_cure_id_i,
    step5_best_cure_selection_rule = step5_selection_rule_i,
    analysis_n = analysis_n_i,
    n_event = n_event_i,
    n_censor = n_censor_i,
    zero_time_offset_days = zero_time_offset_days_i,
    strict_bundle_alignment = strict_bundle_alignment,
    step4_step5_subject_set_match = dataset_validation_list[[ds]]$subject_set_match,
    step4_step5_endpoint_match = dataset_validation_list[[ds]]$endpoint_match,
    step4_step5_transformed_covariate_match = dataset_validation_list[[ds]]$transformed_covariate_match,
    step4_step5_scaling_match = dataset_validation_list[[ds]]$scaling_match,
    run_all_pairs = run_all_pairs,
    pairs_run_n = nrow(all_pairs_df_i),
    primary_alt_model_id = as.character(primary_row_i$alt_model_id[1L]),
    primary_null_model_id = as.character(primary_row_i$null_model_id[1L]),
    primary_latency_dist = as.character(primary_row_i$latency_dist[1L]),
    primary_covariate_structure = as.character(primary_row_i$covariate_structure[1L]),
    primary_d_obs = as.numeric(primary_row_i$d_obs[1L]),
    primary_strict_observed_evaluable = as.logical(primary_row_i$strict_observed_evaluable[1L]),
    primary_strict_bootstrap_target_B = as.integer(primary_row_i$strict_bootstrap_target_B[1L]),
    primary_strict_bootstrap_success_B = as.integer(primary_row_i$strict_bootstrap_success_B[1L]),
    primary_strict_bootstrap_target_met = as.logical(primary_row_i$strict_bootstrap_target_met[1L]),
    primary_strict_p_value_boundary = as.numeric(primary_row_i$strict_p_value_boundary[1L]),
    primary_strict_decision = as.character(primary_row_i$strict_decision[1L]),
    primary_strict_reading_role = as.character(primary_row_i$strict_reading_role[1L]),
    primary_relaxed_observed_evaluable = as.logical(primary_row_i$relaxed_observed_evaluable[1L]),
    primary_relaxed_bootstrap_target_B = as.integer(primary_row_i$relaxed_bootstrap_target_B[1L]),
    primary_relaxed_bootstrap_success_B = as.integer(primary_row_i$relaxed_bootstrap_success_B[1L]),
    primary_relaxed_bootstrap_target_met = as.logical(primary_row_i$relaxed_bootstrap_target_met[1L]),
    primary_relaxed_p_value_boundary = as.numeric(primary_row_i$relaxed_p_value_boundary[1L]),
    primary_relaxed_decision = as.character(primary_row_i$relaxed_decision[1L]),
    primary_relaxed_reading_role = as.character(primary_row_i$relaxed_reading_role[1L]),
    primary_note = as.character(primary_row_i$note[1L]),
    step4_bundle_path = normalizePath(step4_bundle_paths[[ds]], winslash = "/", mustWork = FALSE),
    step5_bundle_path = normalizePath(step5_bundle_paths[[ds]], winslash = "/", mustWork = FALSE),
    step6_bundle_path = normalizePath(primary_pair_path_i, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE
  )
  
  dataset_bundle_objects[[ds]] <- list(
    step = "Step6",
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    dataset_label = ds,
    step3_gate_status = gate_status_i,
    step4_bundle_path = normalizePath(step4_bundle_paths[[ds]], winslash = "/", mustWork = FALSE),
    step5_bundle_path = normalizePath(step5_bundle_paths[[ds]], winslash = "/", mustWork = FALSE),
    strict_bundle_alignment = strict_bundle_alignment,
    boundary_test_specification = list(
      primary_pair_rule = "Step5 best cure model versus matched no-cure counterpart with the same latency family and the same covariate structure.",
      null_model_engine = "survival::survreg",
      alternative_model_engine = "flexsurvcure::flexsurvcure",
      cox_used_in_step6 = FALSE,
      boundary_statistic = "d_obs = max(0, 2 * (logLik_alt - logLik_null))",
      strict_rule = "Observed and bootstrap refits must be stable: fit_ok + converged + no instability warning + finite logLik.",
      relaxed_rule = "Observed and bootstrap refits may be unstable, but must have fit_ok + converged + finite logLik. Relaxed results are exploratory only.",
      alpha = boundary_alpha,
      censoring_bootstrap_method = censoring_bootstrap_method,
      strict_bootstrap_success_target_primary = strict_bootstrap_success_target_primary,
      strict_bootstrap_success_target_sensitivity = strict_bootstrap_success_target_sensitivity,
      relaxed_bootstrap_success_target_primary = relaxed_bootstrap_success_target_primary,
      relaxed_bootstrap_success_target_sensitivity = relaxed_bootstrap_success_target_sensitivity,
      bootstrap_max_attempt_multiplier = bootstrap_max_attempt_multiplier
    ),
    step5_primary_selection = list(
      gate_pass = step5_gate_pass_i,
      gate_reason = step5_gate_reason_i,
      selected_best_cure_id = step5_selected_cure_id_i,
      selected_best_cure_rule = step5_selection_rule_i
    ),
    alignment_validation = dataset_validation_list[[ds]],
    zero_time_offset_days = zero_time_offset_days_i,
    analysis_frame = data_i,
    pair_catalog = pair_catalog,
    all_pairs_results = all_pairs_df_i,
    primary_result = primary_df_i,
    pair_store = pair_store_i,
    session_info = utils::sessionInfo()
  )
  
  all_pairs_results_list[[ds]] <- all_pairs_df_i
  primary_results_list[[ds]] <- primary_df_i
}

test_cure_fraction_all <- bind_rows_fill(all_pairs_results_list)
test_cure_fraction <- bind_rows_fill(primary_results_list)
step6_dataset_manifest <- bind_rows_fill(manifest_rows)

if (nrow(test_cure_fraction_all) == 0L) {
  test_cure_fraction_all <- new_empty_pair_result()
}

if (nrow(test_cure_fraction) == 0L) {
  test_cure_fraction <- new_empty_pair_result()
}

# 🔴 Export: CSV tables and self-contained RDS bundles ===============================

## 🟠 Write: source-of-truth Step6 summary tables ===============================
utils::write.csv(step6_dataset_manifest, file_dataset_manifest, row.names = FALSE, na = "")
utils::write.csv(test_cure_fraction, file_primary_results, row.names = FALSE, na = "")
utils::write.csv(test_cure_fraction_all, file_allpairs_results, row.names = FALSE, na = "")

## 🟠 Save: dataset-level boundary-test bundles ===============================
saveRDS(dataset_bundle_objects$merged, file = dataset_rds_paths[["merged"]], compress = save_rds_compression)
saveRDS(dataset_bundle_objects$pnu, file = dataset_rds_paths[["pnu"]], compress = save_rds_compression)
saveRDS(dataset_bundle_objects$snu, file = dataset_rds_paths[["snu"]], compress = save_rds_compression)

## 🟠 Assign: Step6 workflow objects into current session ===============================
assign("test_cure_fraction", test_cure_fraction, envir = .GlobalEnv)
assign("test_cure_fraction_all", test_cure_fraction_all, envir = .GlobalEnv)
assign("step6_dataset_manifest", step6_dataset_manifest, envir = .GlobalEnv)
assign("step6_pair_catalog", pair_catalog, envir = .GlobalEnv)
assign("step6_dataset_bundles", dataset_bundle_objects, envir = .GlobalEnv)

invisible(list(
  export_path = normalizePath(export_path, winslash = "/", mustWork = FALSE),
  dataset_manifest = step6_dataset_manifest,
  test_cure_fraction = test_cure_fraction,
  test_cure_fraction_all = test_cure_fraction_all
))