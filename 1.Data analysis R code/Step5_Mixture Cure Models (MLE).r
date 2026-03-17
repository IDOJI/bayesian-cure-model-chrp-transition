# 🔴 Rebuild: Step5 MLE cure script after anomaly detection ===============================
# Prior outputs showed clear anomalies:
# 1) fit_registry$success was FALSE for all 80 models despite 78 converged fits
# 2) cure_only.csv and latency_summary.csv were empty
# 3) mean_predictions and subject_predictions contained only Cox-latency models
# 4) flexsurvcure prediction export failed with summary parsing errors
# This rewritten script fixes those issues and adds strict post-fit validation.

# 🔴 Configure: paths and run switches ===============================

DATA_PATH <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
EXPORT_DIR <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦5.Step5_MLE MCM/attachments'

REFIT_MODELS <- TRUE
AUTO_INSTALL_PACKAGES <- FALSE
EXPORT_SUBJECT_PREDICTIONS <- TRUE
STOP_IF_ANY_PREDICTION_FAILURE <- FALSE

SEED <- 20260317L
SCRIPT_VERSION <- "step5_mle_cure_rebuild_v2"

PNU_SITE_VALUE <- "PNU"
SNU_SITE_VALUE <- "SNU"
DATASETS_TO_RUN <- c("PNU", "SNU", "MERGED")
SITE_REFERENCE <- PNU_SITE_VALUE

ID_VAR <- "id"
SITE_VAR <- "site"
SEX_VAR <- "sex_num"
AGE_VAR <- "age_exact_entry"
AGE_FALLBACK_VAR <- "age_int"
TIME_VAR <- "days_followup"
STATUS_VAR <- "status_num"

EVENT_CODE <- 1L
RIGHT_CENSOR_CODES <- c(0L, 2L)

TIME_UNIT_DIVISOR <- 365.25
TIME_EPS_YEARS <- 1e-08
PREDICTION_YEARS <- 1:10
BOUNDARY_EPS <- 1e-04

PARAMETRIC_FAMILIES <- c("exp", "weibull", "llogis", "lnorm")
FLEXSURVCURE_LINK <- "logistic"
FLEXSURV_MAXIT <- 5000L

SMCURE_LINK <- "logit"
SMCURE_VAR <- FALSE
SMCURE_NBOOT <- 0L
SMCURE_EMMAX <- 200L
SMCURE_EPS <- 1e-07

MASTER_RDS_NAME <- "step5_mle_cure_master.rds"
MODEL_GRID_CSV_NAME <- "step5_mle_cure_model_grid.csv"
PREPROC_CSV_NAME <- "step5_mle_cure_preprocessing.csv"
ANALYSIS_SPEC_CSV_NAME <- "step5_mle_cure_analysis_spec.csv"
FIT_REGISTRY_CSV_NAME <- "step5_mle_cure_fit_registry.csv"
COEFFICIENTS_CSV_NAME <- "step5_mle_cure_coefficients.csv"
CURE_ONLY_CSV_NAME <- "step5_mle_cure_cure_only.csv"
MEAN_PRED_CSV_NAME <- "step5_mle_cure_mean_predictions.csv"
LATENCY_SUMMARY_CSV_NAME <- "step5_mle_cure_latency_summary.csv"
SUBJECT_PRED_CSV_NAME <- "step5_mle_cure_subject_predictions.csv"
MANIFEST_CSV_NAME <- "step5_mle_cure_manifest.csv"

# 🔴 Initialize: session and file map ===============================

options(stringsAsFactors = FALSE)
set.seed(SEED)

DATA_PATH <- path.expand(DATA_PATH)
EXPORT_DIR <- path.expand(EXPORT_DIR)

if (!dir.exists(EXPORT_DIR)) {
  dir.create(EXPORT_DIR, recursive = TRUE, showWarnings = FALSE)
}

MASTER_RDS_PATH <- file.path(EXPORT_DIR, MASTER_RDS_NAME)
MODEL_GRID_CSV_PATH <- file.path(EXPORT_DIR, MODEL_GRID_CSV_NAME)
PREPROC_CSV_PATH <- file.path(EXPORT_DIR, PREPROC_CSV_NAME)
ANALYSIS_SPEC_CSV_PATH <- file.path(EXPORT_DIR, ANALYSIS_SPEC_CSV_NAME)
FIT_REGISTRY_CSV_PATH <- file.path(EXPORT_DIR, FIT_REGISTRY_CSV_NAME)
COEFFICIENTS_CSV_PATH <- file.path(EXPORT_DIR, COEFFICIENTS_CSV_NAME)
CURE_ONLY_CSV_PATH <- file.path(EXPORT_DIR, CURE_ONLY_CSV_NAME)
MEAN_PRED_CSV_PATH <- file.path(EXPORT_DIR, MEAN_PRED_CSV_NAME)
LATENCY_SUMMARY_CSV_PATH <- file.path(EXPORT_DIR, LATENCY_SUMMARY_CSV_NAME)
SUBJECT_PRED_CSV_PATH <- file.path(EXPORT_DIR, SUBJECT_PRED_CSV_NAME)
MANIFEST_CSV_PATH <- file.path(EXPORT_DIR, MANIFEST_CSV_NAME)

# 🔴 Check: packages and utility operators ===============================

## 🟠 Ensure: package availability ===============================

ensure_packages <- function(pkgs, auto_install = FALSE) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0L) {
    if (isTRUE(auto_install)) {
      install.packages(missing_pkgs)
    } else {
      stop(
        "Missing required packages: ",
        paste(missing_pkgs, collapse = ", "),
        ". Set AUTO_INSTALL_PACKAGES <- TRUE to install automatically."
      )
    }
  }
  invisible(TRUE)
}

required_pkgs <- c("survival", "flexsurv", "flexsurvcure", "smcure")
ensure_packages(required_pkgs, auto_install = AUTO_INSTALL_PACKAGES)

suppressPackageStartupMessages({
  library(survival)
  library(flexsurv)
  library(flexsurvcure)
  library(smcure)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

# 🔴 Define: empty templates and helper utilities ===============================

## 🟠 Build: empty export tables ===============================

empty_coeff_df <- function() {
  data.frame(
    dataset = character(),
    fit_id = character(),
    component = character(),
    term = character(),
    estimate = numeric(),
    std_error = numeric(),
    zvalue = numeric(),
    pvalue = numeric(),
    lower95 = numeric(),
    upper95 = numeric(),
    source_scale = character(),
    stringsAsFactors = FALSE
  )
}

empty_cure_only_df <- function() {
  data.frame(
    dataset = character(),
    fit_id = character(),
    lane = character(),
    family = character(),
    incidence_target = character(),
    cure_fraction_mean = numeric(),
    cure_fraction_median = numeric(),
    cure_fraction_min = numeric(),
    cure_fraction_max = numeric(),
    uncured_probability_mean = numeric(),
    boundary_flag = logical(),
    stringsAsFactors = FALSE
  )
}

empty_mean_pred_df <- function() {
  data.frame(
    dataset = character(),
    fit_id = character(),
    year = integer(),
    meanS = numeric(),
    meanRisk = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_latency_summary_df <- function() {
  data.frame(
    dataset = character(),
    fit_id = character(),
    year = integer(),
    meanSu = numeric(),
    medianSu = numeric(),
    minSu = numeric(),
    maxSu = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_subject_pred_df <- function() {
  data.frame(
    dataset = character(),
    fit_id = character(),
    subject_uid = character(),
    subject_rowid_dataset = integer(),
    site_id = character(),
    site = character(),
    id = character(),
    year = integer(),
    Shat_pop = numeric(),
    Riskhat_pop = numeric(),
    Shat_uncured = numeric(),
    Riskhat_uncured = numeric(),
    cure_fraction = numeric(),
    uncured_probability = numeric(),
    stringsAsFactors = FALSE
  )
}

## 🟠 Create: generic helpers ===============================

collapse_messages <- function(x) {
  x <- unique(stats::na.omit(as.character(x)))
  if (length(x) == 0L) "" else paste(x, collapse = " | ")
}

combine_rows <- function(x, empty_df_fun = NULL) {
  if (length(x) == 0L) {
    if (is.null(empty_df_fun)) {
      return(data.frame(stringsAsFactors = FALSE))
    }
    return(empty_df_fun())
  }
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

stop_if_missing_columns <- function(df, cols, object_name = "data") {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(object_name, " is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}

safe_pkg_version <- function(pkg) {
  as.character(utils::packageVersion(pkg))
}

make_rhs <- function(vars) {
  vars <- vars[!duplicated(vars)]
  if (length(vars) == 0L) "1" else paste(vars, collapse = " + ")
}

split_var_string <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) {
    character(0)
  } else {
    strsplit(x, "|", fixed = TRUE)[[1L]]
  }
}

safe_loglik <- function(fit) {
  out <- tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_)
  if (!is.finite(out) && !is.null(fit$loglik)) {
    out <- suppressWarnings(as.numeric(fit$loglik))
  }
  out
}

safe_aic_from_loglik <- function(loglik, k) {
  if (is.finite(loglik) && is.finite(k)) -2 * loglik + 2 * k else NA_real_
}

safe_bic_from_loglik <- function(loglik, k, n) {
  if (is.finite(loglik) && is.finite(k) && is.finite(n)) -2 * loglik + log(n) * k else NA_real_
}

clamp01 <- function(x) {
  pmin(pmax(as.numeric(x), 0), 1)
}

find_first_matching_name <- function(x, patterns) {
  x_lower <- tolower(x)
  for (pat in patterns) {
    idx <- which(grepl(pat, x_lower, perl = TRUE))
    if (length(idx) > 0L) {
      return(x[idx[1L]])
    }
  }
  NA_character_
}

coerce_matrix <- function(x, nrow_expected = NULL, ncol_expected = NULL) {
  out <- as.matrix(x)
  storage.mode(out) <- "double"
  if (!is.null(nrow_expected) && nrow(out) != nrow_expected) {
    stop("Unexpected matrix nrow: ", nrow(out), " vs ", nrow_expected)
  }
  if (!is.null(ncol_expected) && ncol(out) != ncol_expected) {
    stop("Unexpected matrix ncol: ", ncol(out), " vs ", ncol_expected)
  }
  out
}

# 🔴 Define: data ingest and preprocessing ===============================

## 🟠 Read: merged source csv ===============================

read_and_clean_input <- function(data_path) {
  if (!file.exists(data_path)) {
    stop("DATA_PATH does not exist: ", data_path)
  }
  
  raw_df <- utils::read.csv(
    file = data_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required_cols <- c(ID_VAR, SITE_VAR, SEX_VAR, TIME_VAR, STATUS_VAR)
  stop_if_missing_columns(raw_df, required_cols, object_name = "Input CSV")
  
  age_var_used <- if (AGE_VAR %in% names(raw_df)) AGE_VAR else AGE_FALLBACK_VAR
  if (!age_var_used %in% names(raw_df)) {
    stop("Neither AGE_VAR nor AGE_FALLBACK_VAR is present in the input data.")
  }
  
  raw_df[[SITE_VAR]] <- trimws(as.character(raw_df[[SITE_VAR]]))
  raw_df[[ID_VAR]] <- as.character(raw_df[[ID_VAR]])
  raw_df[[SEX_VAR]] <- suppressWarnings(as.numeric(as.character(raw_df[[SEX_VAR]])))
  raw_df[[TIME_VAR]] <- suppressWarnings(as.numeric(as.character(raw_df[[TIME_VAR]])))
  raw_df[[STATUS_VAR]] <- suppressWarnings(as.integer(as.character(raw_df[[STATUS_VAR]])))
  raw_df[[age_var_used]] <- suppressWarnings(as.numeric(as.character(raw_df[[age_var_used]])))
  
  essential_cols <- c(SITE_VAR, ID_VAR, SEX_VAR, TIME_VAR, STATUS_VAR, age_var_used)
  keep <- stats::complete.cases(raw_df[, essential_cols, drop = FALSE]) &
    !is.na(raw_df[[TIME_VAR]]) &
    raw_df[[TIME_VAR]] >= 0
  
  dropped_n <- sum(!keep)
  cleaned_df <- raw_df[keep, , drop = FALSE]
  
  cleaned_df$age_raw <- cleaned_df[[age_var_used]]
  cleaned_df$event <- as.integer(cleaned_df[[STATUS_VAR]] == EVENT_CODE)
  cleaned_df$time_years <- pmax(cleaned_df[[TIME_VAR]] / TIME_UNIT_DIVISOR, TIME_EPS_YEARS)
  cleaned_df$site_id <- paste(cleaned_df[[SITE_VAR]], cleaned_df[[ID_VAR]], sep = "::")
  cleaned_df$subject_uid <- sprintf("UID_%06d", seq_len(nrow(cleaned_df)))
  
  attr(cleaned_df, "age_var_used") <- age_var_used
  attr(cleaned_df, "dropped_n") <- dropped_n
  attr(cleaned_df, "raw_n") <- nrow(raw_df)
  cleaned_df
}

## 🟠 Prepare: dataset-specific scaling and site coding ===============================

prepare_analysis_dataset <- function(core_df, dataset_key) {
  dataset_key <- toupper(dataset_key)
  
  if (dataset_key == "PNU") {
    df <- core_df[core_df[[SITE_VAR]] == PNU_SITE_VALUE, , drop = FALSE]
  } else if (dataset_key == "SNU") {
    df <- core_df[core_df[[SITE_VAR]] == SNU_SITE_VALUE, , drop = FALSE]
  } else if (dataset_key == "MERGED") {
    df <- core_df
  } else {
    stop("Unknown dataset_key: ", dataset_key)
  }
  
  if (nrow(df) == 0L) {
    stop("No rows available for dataset: ", dataset_key)
  }
  
  if (sum(df$event, na.rm = TRUE) < 1L) {
    stop("Dataset ", dataset_key, " has no transition events after preprocessing.")
  }
  
  age_center <- mean(df$age_raw, na.rm = TRUE)
  age_sd <- stats::sd(df$age_raw, na.rm = TRUE)
  if (!is.finite(age_sd) || age_sd <= 0) {
    stop("Dataset ", dataset_key, " has non-positive age SD; cannot compute age_s.")
  }
  
  age_scale_2sd <- 2 * age_sd
  df$age_s <- (df$age_raw - age_center) / age_scale_2sd
  df$sex_num <- df[[SEX_VAR]]
  df$age_sex_int <- df$age_s * df$sex_num
  
  site_reference <- ""
  site_nonreference <- ""
  
  if (dataset_key == "MERGED") {
    site_levels <- sort(unique(df[[SITE_VAR]]))
    if (length(site_levels) != 2L) {
      stop(
        "Merged dataset currently supports exactly 2 site levels for Step5. Found: ",
        paste(site_levels, collapse = ", ")
      )
    }
    if (!SITE_REFERENCE %in% site_levels) {
      stop("SITE_REFERENCE (", SITE_REFERENCE, ") is not present in merged data.")
    }
    site_reference <- SITE_REFERENCE
    site_nonreference <- setdiff(site_levels, site_reference)
    if (length(site_nonreference) != 1L) {
      stop("Merged dataset could not determine the non-reference site uniquely.")
    }
    site_nonreference <- site_nonreference[1L]
    df$site_dummy <- as.integer(df[[SITE_VAR]] == site_nonreference)
  } else {
    df$site_dummy <- 0L
  }
  
  df$dataset_key <- dataset_key
  df$subject_rowid_dataset <- seq_len(nrow(df))
  
  preproc_row <- data.frame(
    dataset = dataset_key,
    n = nrow(df),
    n_event = sum(df$event),
    n_censor = sum(df$event == 0L),
    age_var_used = attr(core_df, "age_var_used"),
    age_center = age_center,
    age_scale_2sd = age_scale_2sd,
    site_reference = site_reference,
    site_nonreference = site_nonreference,
    min_followup_years = min(df$time_years, na.rm = TRUE),
    max_followup_years = max(df$time_years, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  list(data = df, preproc = preproc_row)
}

# 🔴 Define: model grid and coefficient extraction ===============================

## 🟠 Build: Step5 cure-model grid ===============================

build_component_vars <- function(dataset_key, interaction_flag, site_flag) {
  vars <- c("age_s", "sex_num")
  if (as.integer(interaction_flag) == 1L) {
    vars <- c(vars, "age_sex_int")
  }
  if (toupper(dataset_key) == "MERGED" && as.integer(site_flag) == 1L) {
    vars <- c(vars, "site_dummy")
  }
  unique(vars)
}

build_step5_grid <- function(dataset_key) {
  dataset_key <- toupper(dataset_key)
  family_values <- c(PARAMETRIC_FAMILIES, "cox_latency")
  incidence_flags <- c(0L, 1L)
  latency_flags <- c(0L, 1L)
  site_flags <- if (dataset_key == "MERGED") c(0L, 1L) else 0L
  
  out <- vector("list", length = length(family_values) * length(incidence_flags) * length(latency_flags) * length(site_flags))
  idx <- 1L
  
  for (family_i in family_values) {
    for (inc_i in incidence_flags) {
      for (lat_i in latency_flags) {
        for (site_i in site_flags) {
          inc_vars <- build_component_vars(dataset_key, inc_i, site_i)
          lat_vars <- build_component_vars(dataset_key, lat_i, site_i)
          
          cov_spec <- if (dataset_key == "MERGED") {
            sprintf("C%d%dS%d", inc_i, lat_i, site_i)
          } else {
            sprintf("C%d%d", inc_i, lat_i)
          }
          
          lane <- if (identical(family_i, "cox_latency")) "PH_semiparametric" else "AFT_parametric"
          engine <- if (identical(family_i, "cox_latency")) "smcure" else "flexsurvcure"
          anc_param <- if (identical(family_i, "cox_latency")) {
            ""
          } else {
            switch(
              family_i,
              exp = "rate",
              weibull = "scale",
              llogis = "scale",
              lnorm = "meanlog",
              stop("Unknown parametric family: ", family_i)
            )
          }
          
          fit_id <- paste(dataset_key, lane, family_i, cov_spec, sep = "__")
          
          out[[idx]] <- data.frame(
            dataset = dataset_key,
            model_class = "cure_MLE",
            lane = lane,
            engine = engine,
            family = family_i,
            anc_param = anc_param,
            incidence_int_flag = inc_i,
            latency_int_flag = lat_i,
            site_flag = site_i,
            cov_spec = cov_spec,
            incidence_terms = make_rhs(inc_vars),
            latency_terms = make_rhs(lat_vars),
            incidence_vars_str = paste(inc_vars, collapse = "|"),
            latency_vars_str = paste(lat_vars, collapse = "|"),
            fit_id = fit_id,
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
    }
  }
  
  grid_df <- combine_rows(out)
  grid_df$fit_index <- seq_len(nrow(grid_df))
  grid_df
}

## 🟠 Extract: flexsurvcure coefficient table ===============================

extract_flexsurvcure_res_table <- function(fit) {
  if (is.null(fit$res.t)) {
    return(NULL)
  }
  
  df <- as.data.frame(fit$res.t, stringsAsFactors = FALSE)
  df$term <- rownames(fit$res.t)
  rownames(df) <- NULL
  
  original_names <- names(df)
  
  est_col <- find_first_matching_name(original_names, c("^est$", "^estimate$"))
  se_col <- find_first_matching_name(original_names, c("^se$", "std", "stderr", "std.error"))
  lower_col <- find_first_matching_name(original_names, c("^l95", "^lower"))
  upper_col <- find_first_matching_name(original_names, c("^u95", "^upper"))
  z_col <- find_first_matching_name(original_names, c("^z$", "zvalue"))
  p_col <- find_first_matching_name(original_names, c("^p$", "pvalue"))
  
  out <- data.frame(
    term = df$term,
    estimate = if (!is.na(est_col)) suppressWarnings(as.numeric(df[[est_col]])) else NA_real_,
    std_error = if (!is.na(se_col)) suppressWarnings(as.numeric(df[[se_col]])) else NA_real_,
    zvalue = if (!is.na(z_col)) suppressWarnings(as.numeric(df[[z_col]])) else NA_real_,
    pvalue = if (!is.na(p_col)) suppressWarnings(as.numeric(df[[p_col]])) else NA_real_,
    lower95 = if (!is.na(lower_col)) suppressWarnings(as.numeric(df[[lower_col]])) else NA_real_,
    upper95 = if (!is.na(upper_col)) suppressWarnings(as.numeric(df[[upper_col]])) else NA_real_,
    stringsAsFactors = FALSE
  )
  
  if (all(is.na(out$lower95)) && any(is.finite(out$std_error))) {
    out$lower95 <- out$estimate - 1.96 * out$std_error
  }
  if (all(is.na(out$upper95)) && any(is.finite(out$std_error))) {
    out$upper95 <- out$estimate + 1.96 * out$std_error
  }
  if (all(is.na(out$zvalue)) && any(is.finite(out$std_error))) {
    out$zvalue <- out$estimate / out$std_error
  }
  if (all(is.na(out$pvalue)) && any(is.finite(out$zvalue))) {
    out$pvalue <- 2 * stats::pnorm(-abs(out$zvalue))
  }
  
  out
}

extract_flexsurvcure_coef_table <- function(fit, dataset_key, fit_id, incidence_vars) {
  res_df <- extract_flexsurvcure_res_table(fit)
  if (is.null(res_df)) {
    return(empty_coeff_df())
  }
  
  latency_parameter_prefix <- c("rate", "scale", "shape", "meanlog", "sdlog")
  
  res_df$component <- vapply(
    res_df$term,
    FUN.VALUE = character(1L),
    FUN = function(term_i) {
      if (term_i == "theta" || term_i %in% incidence_vars) {
        return("incidence")
      }
      if (term_i %in% latency_parameter_prefix) {
        return("latency")
      }
      if (any(startsWith(term_i, paste0(latency_parameter_prefix, "(")))) {
        return("latency")
      }
      "latency"
    }
  )
  
  res_df$source_scale <- ifelse(
    res_df$component == "incidence",
    "logit_cure_probability",
    "package_native_real_line"
  )
  res_df$dataset <- dataset_key
  res_df$fit_id <- fit_id
  
  res_df[, c("dataset", "fit_id", "component", "term", "estimate", "std_error", "zvalue", "pvalue", "lower95", "upper95", "source_scale"), drop = FALSE]
}

## 🟠 Extract: smcure coefficient table ===============================

extract_smcure_coef_table <- function(fit, dataset_key, fit_id) {
  incidence_df <- data.frame(
    dataset = dataset_key,
    fit_id = fit_id,
    component = "incidence",
    term = fit$bnm,
    estimate = as.numeric(fit$b),
    std_error = NA_real_,
    zvalue = NA_real_,
    pvalue = NA_real_,
    lower95 = NA_real_,
    upper95 = NA_real_,
    source_scale = "logit_uncured_probability",
    stringsAsFactors = FALSE
  )
  
  latency_df <- data.frame(
    dataset = dataset_key,
    fit_id = fit_id,
    component = "latency",
    term = fit$betanm,
    estimate = as.numeric(fit$beta),
    std_error = NA_real_,
    zvalue = NA_real_,
    pvalue = NA_real_,
    lower95 = NA_real_,
    upper95 = NA_real_,
    source_scale = "cox_ph_log_hazard_ratio",
    stringsAsFactors = FALSE
  )
  
  if (!is.null(fit$b_sd)) {
    incidence_df$std_error <- as.numeric(fit$b_sd)
    incidence_df$zvalue <- as.numeric(fit$b_zvalue)
    incidence_df$pvalue <- as.numeric(fit$b_pvalue)
    incidence_df$lower95 <- incidence_df$estimate - 1.96 * incidence_df$std_error
    incidence_df$upper95 <- incidence_df$estimate + 1.96 * incidence_df$std_error
  }
  
  if (!is.null(fit$beta_sd)) {
    latency_df$std_error <- as.numeric(fit$beta_sd)
    latency_df$zvalue <- as.numeric(fit$beta_zvalue)
    latency_df$pvalue <- as.numeric(fit$beta_pvalue)
    latency_df$lower95 <- latency_df$estimate - 1.96 * latency_df$std_error
    latency_df$upper95 <- latency_df$estimate + 1.96 * latency_df$std_error
  }
  
  out <- rbind(incidence_df, latency_df)
  rownames(out) <- NULL
  out
}

# 🔴 Define: prediction helpers ===============================

## 🟠 Parse: flexible summary output from flexsurvcure ===============================

extract_single_flexsurv_survival <- function(summary_obj, eval_times) {
  candidate <- NULL
  
  if (is.matrix(summary_obj) || is.data.frame(summary_obj)) {
    candidate <- summary_obj
  } else if (is.list(summary_obj)) {
    obj_names <- names(summary_obj)
    keep_idx <- seq_along(summary_obj)
    if (!is.null(obj_names)) {
      keep_idx <- which(obj_names != "X")
      if (length(keep_idx) == 0L) {
        keep_idx <- seq_along(summary_obj)
      }
    }
    for (ii in keep_idx) {
      if (is.matrix(summary_obj[[ii]]) || is.data.frame(summary_obj[[ii]])) {
        candidate <- summary_obj[[ii]]
        break
      }
    }
  }
  
  if (is.null(candidate)) {
    stop("Could not parse summary.flexsurvreg output for a single-row prediction.")
  }
  
  df <- as.data.frame(candidate, stringsAsFactors = FALSE)
  nms <- names(df)
  
  est_col <- if ("est" %in% nms) {
    "est"
  } else {
    num_cols <- nms[vapply(df, is.numeric, logical(1))]
    if (length(num_cols) == 0L) {
      stop("No numeric estimate column found in parsed summary output.")
    }
    num_cols[1L]
  }
  
  if ("time" %in% nms) {
    time_vec <- suppressWarnings(as.numeric(df$time))
    row_idx <- vapply(
      eval_times,
      FUN.VALUE = integer(1L),
      FUN = function(tt) {
        exact_idx <- which(abs(time_vec - tt) < 1e-10)
        if (length(exact_idx) > 0L) exact_idx[1L] else which.min(abs(time_vec - tt))
      }
    )
    est <- suppressWarnings(as.numeric(df[[est_col]][row_idx]))
  } else {
    if (nrow(df) < length(eval_times)) {
      stop("Parsed summary output has fewer rows than requested prediction times.")
    }
    est <- suppressWarnings(as.numeric(df[[est_col]][seq_len(length(eval_times))]))
  }
  
  clamp01(est)
}

## 🟠 Predict: parametric overall survival matrix ===============================

predict_parametric_survival_matrix <- function(fit, data_df, eval_times) {
  n_subject <- nrow(data_df)
  n_times <- length(eval_times)
  out <- matrix(NA_real_, nrow = n_subject, ncol = n_times)
  fail_messages <- character(0)
  
  for (i in seq_len(n_subject)) {
    out[i, ] <- tryCatch(
      {
        summ_i <- suppressWarnings(
          summary(
            fit,
            newdata = data_df[i, , drop = FALSE],
            type = "survival",
            t = eval_times,
            start = 0,
            ci = FALSE,
            se = FALSE,
            B = 0,
            tidy = FALSE
          )
        )
        extract_single_flexsurv_survival(summ_i, eval_times)
      },
      error = function(e) {
        fail_messages <<- c(fail_messages, paste0("row ", i, ": ", conditionMessage(e)))
        rep(NA_real_, n_times)
      }
    )
  }
  
  list(
    matrix = clamp01(out),
    n_failed_rows = sum(apply(is.na(out), 1L, all)),
    failure_messages = unique(fail_messages)
  )
}

## 🟠 Compute: parametric cure fraction from fitted coefficients ===============================

compute_parametric_cure_fraction <- function(fit, data_df, incidence_vars) {
  res_df <- extract_flexsurvcure_res_table(fit)
  if (is.null(res_df)) {
    stop("Could not retrieve flexsurvcure coefficient table for cure-fraction computation.")
  }
  
  intercept <- res_df$estimate[res_df$term == "theta"]
  if (length(intercept) == 0L || !is.finite(intercept[1L])) {
    intercept <- 0
  } else {
    intercept <- intercept[1L]
  }
  
  eta <- rep(intercept, nrow(data_df))
  
  if (length(incidence_vars) > 0L) {
    x_mat <- as.matrix(data_df[, incidence_vars, drop = FALSE])
    storage.mode(x_mat) <- "double"
    
    beta_vec <- vapply(
      incidence_vars,
      FUN.VALUE = numeric(1L),
      FUN = function(v) {
        tmp <- res_df$estimate[res_df$term == v]
        if (length(tmp) == 0L || !is.finite(tmp[1L])) 0 else tmp[1L]
      }
    )
    
    eta <- as.numeric(eta + x_mat %*% beta_vec)
  }
  
  clamp01(stats::plogis(eta))
}

## 🟠 Evaluate: semiparametric step-function predictions ===============================

step_eval_single <- function(time_vec, surv_vec, eval_times) {
  tt <- c(0, as.numeric(time_vec))
  ss <- c(1, as.numeric(surv_vec))
  
  ord <- order(tt, seq_along(tt))
  tt <- tt[ord]
  ss <- ss[ord]
  
  keep <- !duplicated(tt, fromLast = TRUE)
  tt <- tt[keep]
  ss <- ss[keep]
  
  idx <- findInterval(eval_times, tt, rightmost.closed = TRUE, all.inside = TRUE)
  idx[idx < 1L] <- 1L
  idx[idx > length(ss)] <- length(ss)
  
  clamp01(ss[idx])
}

step_eval_matrix <- function(time_vec, surv_mat, eval_times) {
  n_subject <- ncol(surv_mat)
  out <- matrix(NA_real_, nrow = n_subject, ncol = length(eval_times))
  for (j in seq_len(n_subject)) {
    out[j, ] <- step_eval_single(time_vec = time_vec, surv_vec = surv_mat[, j], eval_times = eval_times)
  }
  out
}

compute_uncured_survival <- function(pop_surv_mat, cure_fraction) {
  pop_surv_mat <- coerce_matrix(pop_surv_mat)
  cure_fraction <- as.numeric(cure_fraction)
  uncured_prob <- pmax(1 - cure_fraction, .Machine$double.eps)
  out <- sweep(pop_surv_mat, 1L, cure_fraction, FUN = "-")
  out <- sweep(out, 1L, uncured_prob, FUN = "/")
  clamp01(out)
}

## 🟠 Predict: semiparametric Cox-latency matrices ===============================

predict_smcure_matrices <- function(fit, data_df, incidence_vars, latency_vars, eval_times) {
  newX <- stats::model.matrix(
    object = stats::as.formula(paste0("~ ", make_rhs(latency_vars))),
    data = data_df
  )
  if (ncol(newX) <= 1L) {
    stop("Latency design matrix has no covariate columns after removing the intercept.")
  }
  newX <- newX[, -1L, drop = FALSE]
  
  newZ <- as.matrix(data_df[, incidence_vars, drop = FALSE])
  storage.mode(newX) <- "double"
  storage.mode(newZ) <- "double"
  
  pred_obj <- smcure::predictsmcure(
    object = fit,
    newX = newX,
    newZ = newZ,
    model = "ph"
  )
  
  pred_mat <- as.matrix(pred_obj$prediction)
  storage.mode(pred_mat) <- "double"
  
  time_vec <- pred_mat[, ncol(pred_mat)]
  surv_curves <- pred_mat[, -ncol(pred_mat), drop = FALSE]
  
  ord <- order(time_vec, seq_along(time_vec))
  time_vec <- time_vec[ord]
  surv_curves <- surv_curves[ord, , drop = FALSE]
  
  pop_surv_mat <- step_eval_matrix(
    time_vec = time_vec,
    surv_mat = surv_curves,
    eval_times = eval_times
  )
  
  uncure_prob <- clamp01(as.numeric(pred_obj$newuncureprob))
  cure_fraction <- clamp01(1 - uncure_prob)
  uncured_surv_mat <- compute_uncured_survival(pop_surv_mat = pop_surv_mat, cure_fraction = cure_fraction)
  
  list(
    pop_surv_mat = coerce_matrix(pop_surv_mat, nrow_expected = nrow(data_df), ncol_expected = length(eval_times)),
    uncured_surv_mat = coerce_matrix(uncured_surv_mat, nrow_expected = nrow(data_df), ncol_expected = length(eval_times)),
    cure_fraction = cure_fraction,
    n_failed_rows = 0L,
    failure_messages = character(0)
  )
}

# 🔴 Define: output table builders ===============================

## 🟠 Create: subject-level long predictions ===============================

make_subject_prediction_long <- function(data_df, dataset_key, fit_id, pop_surv_mat, uncured_surv_mat, cure_fraction) {
  pop_surv_mat <- coerce_matrix(pop_surv_mat, nrow_expected = nrow(data_df), ncol_expected = length(PREDICTION_YEARS))
  uncured_surv_mat <- coerce_matrix(uncured_surv_mat, nrow_expected = nrow(data_df), ncol_expected = length(PREDICTION_YEARS))
  cure_fraction <- as.numeric(cure_fraction)
  
  n_subject <- nrow(data_df)
  n_years <- ncol(pop_surv_mat)
  
  data.frame(
    dataset = dataset_key,
    fit_id = fit_id,
    subject_uid = rep(data_df$subject_uid, each = n_years),
    subject_rowid_dataset = rep(data_df$subject_rowid_dataset, each = n_years),
    site_id = rep(data_df$site_id, each = n_years),
    site = rep(data_df[[SITE_VAR]], each = n_years),
    id = rep(data_df[[ID_VAR]], each = n_years),
    year = rep(PREDICTION_YEARS, times = n_subject),
    Shat_pop = as.vector(t(pop_surv_mat)),
    Riskhat_pop = 1 - as.vector(t(pop_surv_mat)),
    Shat_uncured = as.vector(t(uncured_surv_mat)),
    Riskhat_uncured = 1 - as.vector(t(uncured_surv_mat)),
    cure_fraction = rep(cure_fraction, each = n_years),
    uncured_probability = rep(1 - cure_fraction, each = n_years),
    stringsAsFactors = FALSE
  )
}

## 🟠 Create: cohort-average and latency summaries ===============================

make_mean_prediction_df <- function(dataset_key, fit_id, pop_surv_mat) {
  pop_surv_mat <- coerce_matrix(pop_surv_mat, ncol_expected = length(PREDICTION_YEARS))
  mean_s <- colMeans(pop_surv_mat, na.rm = TRUE)
  data.frame(
    dataset = dataset_key,
    fit_id = fit_id,
    year = PREDICTION_YEARS,
    meanS = mean_s,
    meanRisk = 1 - mean_s,
    stringsAsFactors = FALSE
  )
}

make_latency_summary_df <- function(dataset_key, fit_id, uncured_surv_mat) {
  uncured_surv_mat <- coerce_matrix(uncured_surv_mat, ncol_expected = length(PREDICTION_YEARS))
  data.frame(
    dataset = dataset_key,
    fit_id = fit_id,
    year = PREDICTION_YEARS,
    meanSu = colMeans(uncured_surv_mat, na.rm = TRUE),
    medianSu = apply(uncured_surv_mat, 2L, stats::median, na.rm = TRUE),
    minSu = apply(uncured_surv_mat, 2L, min, na.rm = TRUE),
    maxSu = apply(uncured_surv_mat, 2L, max, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

make_cure_only_row <- function(dataset_key, fit_id, lane, family, incidence_target, cure_fraction, boundary_flag) {
  cure_fraction <- as.numeric(cure_fraction)
  data.frame(
    dataset = dataset_key,
    fit_id = fit_id,
    lane = lane,
    family = family,
    incidence_target = incidence_target,
    cure_fraction_mean = mean(cure_fraction, na.rm = TRUE),
    cure_fraction_median = stats::median(cure_fraction, na.rm = TRUE),
    cure_fraction_min = min(cure_fraction, na.rm = TRUE),
    cure_fraction_max = max(cure_fraction, na.rm = TRUE),
    uncured_probability_mean = mean(1 - cure_fraction, na.rm = TRUE),
    boundary_flag = isTRUE(boundary_flag),
    stringsAsFactors = FALSE
  )
}

# 🔴 Define: model-fitting wrappers ===============================

## 🟠 Fit: one parametric AFT cure model ===============================

fit_one_parametric_cure_model <- function(spec_row, data_df) {
  fit_started_at <- Sys.time()
  fit_warning_messages <- character(0)
  fit_error_messages <- character(0)
  pred_warning_messages <- character(0)
  pred_error_messages <- character(0)
  
  incidence_vars <- split_var_string(spec_row$incidence_vars_str)
  latency_vars <- split_var_string(spec_row$latency_vars_str)
  
  surv_formula <- stats::as.formula(
    paste0("survival::Surv(time_years, event) ~ ", make_rhs(incidence_vars))
  )
  
  anc_list <- stats::setNames(
    object = list(stats::as.formula(paste0("~ ", make_rhs(latency_vars)))),
    nm = spec_row$anc_param
  )
  
  fit_obj <- tryCatch(
    withCallingHandlers(
      flexsurvcure::flexsurvcure(
        formula = surv_formula,
        data = data_df,
        dist = spec_row$family,
        link = FLEXSURVCURE_LINK,
        mixture = TRUE,
        anc = anc_list,
        na.action = stats::na.omit,
        control = list(maxit = FLEXSURV_MAXIT)
      ),
      warning = function(w) {
        fit_warning_messages <<- c(fit_warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      structure(list(message = conditionMessage(e)), class = "step5_fit_error")
    }
  )
  
  fit_time_sec <- as.numeric(difftime(Sys.time(), fit_started_at, units = "secs"))
  
  if (inherits(fit_obj, "step5_fit_error")) {
    registry_row <- data.frame(
      dataset = spec_row$dataset,
      fit_id = spec_row$fit_id,
      model_class = spec_row$model_class,
      lane = spec_row$lane,
      engine = spec_row$engine,
      family = spec_row$family,
      cov_spec = spec_row$cov_spec,
      incidence_int_flag = spec_row$incidence_int_flag,
      latency_int_flag = spec_row$latency_int_flag,
      site_flag = spec_row$site_flag,
      incidence_terms = spec_row$incidence_terms,
      latency_terms = spec_row$latency_terms,
      n = nrow(data_df),
      n_event = sum(data_df$event),
      n_censor = sum(data_df$event == 0L),
      n_parameters = NA_real_,
      logLik = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      converged_flag = FALSE,
      incidence_target = "cure_probability",
      mean_cure_fraction = NA_real_,
      boundary_flag = NA,
      fit_success = FALSE,
      prediction_success = FALSE,
      success = FALSE,
      n_prediction_failed_rows = NA_integer_,
      fit_time_sec = fit_time_sec,
      prediction_time_sec = NA_real_,
      error_stage = "fit",
      error_message = fit_obj$message,
      warning_messages = collapse_messages(fit_warning_messages),
      stringsAsFactors = FALSE
    )
    
    return(list(
      fit_bundle = list(
        spec = spec_row,
        fit = NULL,
        engine = spec_row$engine,
        success = FALSE,
        fit_success = FALSE,
        prediction_success = FALSE,
        warnings = fit_warning_messages,
        error_message = fit_obj$message
      ),
      registry = registry_row,
      coefficient_table = empty_coeff_df(),
      subject_predictions = empty_subject_pred_df(),
      mean_predictions = empty_mean_pred_df(),
      latency_summary = empty_latency_summary_df(),
      cure_only = empty_cure_only_df()
    ))
  }
  
  coefficient_table <- extract_flexsurvcure_coef_table(
    fit = fit_obj,
    dataset_key = spec_row$dataset,
    fit_id = spec_row$fit_id,
    incidence_vars = incidence_vars
  )
  
  pred_started_at <- Sys.time()
  
  pop_surv_mat <- NULL
  uncured_surv_mat <- NULL
  cure_fraction <- NULL
  n_failed_rows <- NA_integer_
  
  pred_success <- tryCatch(
    withCallingHandlers(
      {
        pop_pred <- predict_parametric_survival_matrix(
          fit = fit_obj,
          data_df = data_df,
          eval_times = PREDICTION_YEARS
        )
        pop_surv_mat <- pop_pred$matrix
        n_failed_rows <- pop_pred$n_failed_rows
        if (length(pop_pred$failure_messages) > 0L) {
          pred_error_messages <<- c(pred_error_messages, pop_pred$failure_messages)
        }
        
        cure_fraction <- compute_parametric_cure_fraction(
          fit = fit_obj,
          data_df = data_df,
          incidence_vars = incidence_vars
        )
        uncured_surv_mat <- compute_uncured_survival(
          pop_surv_mat = pop_surv_mat,
          cure_fraction = cure_fraction
        )
        TRUE
      },
      warning = function(w) {
        pred_warning_messages <<- c(pred_warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      pred_error_messages <<- c(pred_error_messages, conditionMessage(e))
      FALSE
    }
  )
  
  prediction_time_sec <- as.numeric(difftime(Sys.time(), pred_started_at, units = "secs"))
  
  if (!isTRUE(pred_success)) {
    pop_surv_mat <- NULL
    uncured_surv_mat <- NULL
    cure_fraction <- NULL
  }
  
  subject_predictions <- if (isTRUE(pred_success)) {
    make_subject_prediction_long(
      data_df = data_df,
      dataset_key = spec_row$dataset,
      fit_id = spec_row$fit_id,
      pop_surv_mat = pop_surv_mat,
      uncured_surv_mat = uncured_surv_mat,
      cure_fraction = cure_fraction
    )
  } else {
    empty_subject_pred_df()
  }
  
  mean_predictions <- if (isTRUE(pred_success)) {
    make_mean_prediction_df(
      dataset_key = spec_row$dataset,
      fit_id = spec_row$fit_id,
      pop_surv_mat = pop_surv_mat
    )
  } else {
    empty_mean_pred_df()
  }
  
  latency_summary <- if (isTRUE(pred_success)) {
    make_latency_summary_df(
      dataset_key = spec_row$dataset,
      fit_id = spec_row$fit_id,
      uncured_surv_mat = uncured_surv_mat
    )
  } else {
    empty_latency_summary_df()
  }
  
  cure_only <- if (isTRUE(pred_success)) {
    make_cure_only_row(
      dataset_key = spec_row$dataset,
      fit_id = spec_row$fit_id,
      lane = spec_row$lane,
      family = spec_row$family,
      incidence_target = "cure_probability",
      cure_fraction = cure_fraction,
      boundary_flag = any(cure_fraction < BOUNDARY_EPS | cure_fraction > (1 - BOUNDARY_EPS))
    )
  } else {
    empty_cure_only_df()
  }
  
  loglik_val <- safe_loglik(fit_obj)
  n_params <- if (nrow(coefficient_table) > 0L) nrow(coefficient_table) else suppressWarnings(length(stats::coef(fit_obj)))
  aic_val <- safe_aic_from_loglik(loglik_val, n_params)
  bic_val <- safe_bic_from_loglik(loglik_val, n_params, nrow(data_df))
  
  converged_flag <- TRUE
  if (!is.null(fit_obj$converged)) {
    converged_flag <- isTRUE(fit_obj$converged)
  }
  if (!is.null(fit_obj$opt$convergence)) {
    converged_flag <- converged_flag && identical(as.integer(fit_obj$opt$convergence), 0L)
  }
  
  prediction_complete <- isTRUE(pred_success) && identical(n_failed_rows, 0L)
  
  registry_row <- data.frame(
    dataset = spec_row$dataset,
    fit_id = spec_row$fit_id,
    model_class = spec_row$model_class,
    lane = spec_row$lane,
    engine = spec_row$engine,
    family = spec_row$family,
    cov_spec = spec_row$cov_spec,
    incidence_int_flag = spec_row$incidence_int_flag,
    latency_int_flag = spec_row$latency_int_flag,
    site_flag = spec_row$site_flag,
    incidence_terms = spec_row$incidence_terms,
    latency_terms = spec_row$latency_terms,
    n = nrow(data_df),
    n_event = sum(data_df$event),
    n_censor = sum(data_df$event == 0L),
    n_parameters = n_params,
    logLik = loglik_val,
    AIC = aic_val,
    BIC = bic_val,
    converged_flag = converged_flag,
    incidence_target = "cure_probability",
    mean_cure_fraction = if (isTRUE(pred_success)) mean(cure_fraction, na.rm = TRUE) else NA_real_,
    boundary_flag = if (isTRUE(pred_success)) any(cure_fraction < BOUNDARY_EPS | cure_fraction > (1 - BOUNDARY_EPS)) else NA,
    fit_success = TRUE,
    prediction_success = prediction_complete,
    success = TRUE && prediction_complete,
    n_prediction_failed_rows = n_failed_rows,
    fit_time_sec = fit_time_sec,
    prediction_time_sec = prediction_time_sec,
    error_stage = if (isTRUE(pred_success)) "" else "prediction",
    error_message = collapse_messages(pred_error_messages),
    warning_messages = collapse_messages(c(fit_warning_messages, pred_warning_messages)),
    stringsAsFactors = FALSE
  )
  
  list(
    fit_bundle = list(
      spec = spec_row,
      fit = fit_obj,
      engine = spec_row$engine,
      success = TRUE && prediction_complete,
      fit_success = TRUE,
      prediction_success = prediction_complete,
      warnings = c(fit_warning_messages, pred_warning_messages),
      error_message = collapse_messages(pred_error_messages)
    ),
    registry = registry_row,
    coefficient_table = coefficient_table,
    subject_predictions = subject_predictions,
    mean_predictions = mean_predictions,
    latency_summary = latency_summary,
    cure_only = cure_only
  )
}

## 🟠 Fit: one semiparametric Cox-latency cure model ===============================

fit_one_smcure_ph_model <- function(spec_row, data_df) {
  fit_started_at <- Sys.time()
  fit_warning_messages <- character(0)
  fit_error_messages <- character(0)
  pred_warning_messages <- character(0)
  pred_error_messages <- character(0)
  
  incidence_vars <- split_var_string(spec_row$incidence_vars_str)
  latency_vars <- split_var_string(spec_row$latency_vars_str)
  
  surv_formula <- stats::as.formula(
    paste0("survival::Surv(time_years, event) ~ ", make_rhs(latency_vars))
  )
  cure_formula <- stats::as.formula(
    paste0("~ ", make_rhs(incidence_vars))
  )
  
  fit_obj <- tryCatch(
    withCallingHandlers(
      {
        invisible(capture.output(
          tmp_fit <- smcure::smcure(
            formula = surv_formula,
            cureform = cure_formula,
            data = data_df,
            na.action = stats::na.omit,
            model = "ph",
            link = SMCURE_LINK,
            Var = SMCURE_VAR,
            emmax = SMCURE_EMMAX,
            eps = SMCURE_EPS,
            nboot = if (SMCURE_VAR) SMCURE_NBOOT else 0L
          )
        ))
        tmp_fit
      },
      warning = function(w) {
        fit_warning_messages <<- c(fit_warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      structure(list(message = conditionMessage(e)), class = "step5_fit_error")
    }
  )
  
  fit_time_sec <- as.numeric(difftime(Sys.time(), fit_started_at, units = "secs"))
  
  if (inherits(fit_obj, "step5_fit_error")) {
    registry_row <- data.frame(
      dataset = spec_row$dataset,
      fit_id = spec_row$fit_id,
      model_class = spec_row$model_class,
      lane = spec_row$lane,
      engine = spec_row$engine,
      family = spec_row$family,
      cov_spec = spec_row$cov_spec,
      incidence_int_flag = spec_row$incidence_int_flag,
      latency_int_flag = spec_row$latency_int_flag,
      site_flag = spec_row$site_flag,
      incidence_terms = spec_row$incidence_terms,
      latency_terms = spec_row$latency_terms,
      n = nrow(data_df),
      n_event = sum(data_df$event),
      n_censor = sum(data_df$event == 0L),
      n_parameters = NA_real_,
      logLik = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      converged_flag = FALSE,
      incidence_target = "uncured_probability",
      mean_cure_fraction = NA_real_,
      boundary_flag = NA,
      fit_success = FALSE,
      prediction_success = FALSE,
      success = FALSE,
      n_prediction_failed_rows = NA_integer_,
      fit_time_sec = fit_time_sec,
      prediction_time_sec = NA_real_,
      error_stage = "fit",
      error_message = fit_obj$message,
      warning_messages = collapse_messages(fit_warning_messages),
      stringsAsFactors = FALSE
    )
    
    return(list(
      fit_bundle = list(
        spec = spec_row,
        fit = NULL,
        engine = spec_row$engine,
        success = FALSE,
        fit_success = FALSE,
        prediction_success = FALSE,
        warnings = fit_warning_messages,
        error_message = fit_obj$message
      ),
      registry = registry_row,
      coefficient_table = empty_coeff_df(),
      subject_predictions = empty_subject_pred_df(),
      mean_predictions = empty_mean_pred_df(),
      latency_summary = empty_latency_summary_df(),
      cure_only = empty_cure_only_df()
    ))
  }
  
  coefficient_table <- extract_smcure_coef_table(
    fit = fit_obj,
    dataset_key = spec_row$dataset,
    fit_id = spec_row$fit_id
  )
  
  pred_started_at <- Sys.time()
  
  pop_surv_mat <- NULL
  uncured_surv_mat <- NULL
  cure_fraction <- NULL
  n_failed_rows <- 0L
  
  pred_success <- tryCatch(
    withCallingHandlers(
      {
        pred_obj <- predict_smcure_matrices(
          fit = fit_obj,
          data_df = data_df,
          incidence_vars = incidence_vars,
          latency_vars = latency_vars,
          eval_times = PREDICTION_YEARS
        )
        pop_surv_mat <- pred_obj$pop_surv_mat
        uncured_surv_mat <- pred_obj$uncured_surv_mat
        cure_fraction <- pred_obj$cure_fraction
        n_failed_rows <- pred_obj$n_failed_rows
        if (length(pred_obj$failure_messages) > 0L) {
          pred_error_messages <<- c(pred_error_messages, pred_obj$failure_messages)
        }
        TRUE
      },
      warning = function(w) {
        pred_warning_messages <<- c(pred_warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      pred_error_messages <<- c(pred_error_messages, conditionMessage(e))
      FALSE
    }
  )
  
  prediction_time_sec <- as.numeric(difftime(Sys.time(), pred_started_at, units = "secs"))
  prediction_complete <- isTRUE(pred_success) && identical(n_failed_rows, 0L)
  
  subject_predictions <- if (isTRUE(pred_success)) {
    make_subject_prediction_long(
      data_df = data_df,
      dataset_key = spec_row$dataset,
      fit_id = spec_row$fit_id,
      pop_surv_mat = pop_surv_mat,
      uncured_surv_mat = uncured_surv_mat,
      cure_fraction = cure_fraction
    )
  } else {
    empty_subject_pred_df()
  }
  
  mean_predictions <- if (isTRUE(pred_success)) {
    make_mean_prediction_df(
      dataset_key = spec_row$dataset,
      fit_id = spec_row$fit_id,
      pop_surv_mat = pop_surv_mat
    )
  } else {
    empty_mean_pred_df()
  }
  
  latency_summary <- if (isTRUE(pred_success)) {
    make_latency_summary_df(
      dataset_key = spec_row$dataset,
      fit_id = spec_row$fit_id,
      uncured_surv_mat = uncured_surv_mat
    )
  } else {
    empty_latency_summary_df()
  }
  
  cure_only <- if (isTRUE(pred_success)) {
    make_cure_only_row(
      dataset_key = spec_row$dataset,
      fit_id = spec_row$fit_id,
      lane = spec_row$lane,
      family = spec_row$family,
      incidence_target = "uncured_probability",
      cure_fraction = cure_fraction,
      boundary_flag = any(cure_fraction < BOUNDARY_EPS | cure_fraction > (1 - BOUNDARY_EPS))
    )
  } else {
    empty_cure_only_df()
  }
  
  n_params <- length(fit_obj$b %||% numeric()) + length(fit_obj$beta %||% numeric())
  converged_flag <- all(is.finite(c(fit_obj$b %||% numeric(), fit_obj$beta %||% numeric())))
  
  registry_row <- data.frame(
    dataset = spec_row$dataset,
    fit_id = spec_row$fit_id,
    model_class = spec_row$model_class,
    lane = spec_row$lane,
    engine = spec_row$engine,
    family = spec_row$family,
    cov_spec = spec_row$cov_spec,
    incidence_int_flag = spec_row$incidence_int_flag,
    latency_int_flag = spec_row$latency_int_flag,
    site_flag = spec_row$site_flag,
    incidence_terms = spec_row$incidence_terms,
    latency_terms = spec_row$latency_terms,
    n = nrow(data_df),
    n_event = sum(data_df$event),
    n_censor = sum(data_df$event == 0L),
    n_parameters = n_params,
    logLik = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    converged_flag = converged_flag,
    incidence_target = "uncured_probability",
    mean_cure_fraction = if (isTRUE(pred_success)) mean(cure_fraction, na.rm = TRUE) else NA_real_,
    boundary_flag = if (isTRUE(pred_success)) any(cure_fraction < BOUNDARY_EPS | cure_fraction > (1 - BOUNDARY_EPS)) else NA,
    fit_success = TRUE,
    prediction_success = prediction_complete,
    success = TRUE && prediction_complete,
    n_prediction_failed_rows = n_failed_rows,
    fit_time_sec = fit_time_sec,
    prediction_time_sec = prediction_time_sec,
    error_stage = if (isTRUE(pred_success)) "" else "prediction",
    error_message = collapse_messages(pred_error_messages),
    warning_messages = collapse_messages(c(fit_warning_messages, pred_warning_messages)),
    stringsAsFactors = FALSE
  )
  
  list(
    fit_bundle = list(
      spec = spec_row,
      fit = fit_obj,
      engine = spec_row$engine,
      success = TRUE && prediction_complete,
      fit_success = TRUE,
      prediction_success = prediction_complete,
      warnings = c(fit_warning_messages, pred_warning_messages),
      error_message = collapse_messages(pred_error_messages)
    ),
    registry = registry_row,
    coefficient_table = coefficient_table,
    subject_predictions = subject_predictions,
    mean_predictions = mean_predictions,
    latency_summary = latency_summary,
    cure_only = cure_only
  )
}

# 🔴 Define: validation helpers for non-broken exports ===============================

## 🟠 Validate: internal consistency of exported tables ===============================

validate_step5_outputs <- function(model_grid, fit_registry, cure_only, mean_predictions, latency_summary, subject_predictions, analysis_data_list) {
  if (nrow(model_grid) == 0L) {
    stop("Model grid is empty.")
  }
  
  if (nrow(fit_registry) != nrow(model_grid)) {
    stop("fit_registry row count does not match model_grid row count.")
  }
  
  if (!all(model_grid$fit_id %in% fit_registry$fit_id)) {
    stop("Some fit_id values in model_grid are missing from fit_registry.")
  }
  
  required_fit_cols <- c("fit_success", "prediction_success", "success", "n_prediction_failed_rows")
  missing_fit_cols <- setdiff(required_fit_cols, names(fit_registry))
  if (length(missing_fit_cols) > 0L) {
    stop("fit_registry is missing required QC columns: ", paste(missing_fit_cols, collapse = ", "))
  }
  
  pred_ok_ids <- fit_registry$fit_id[isTRUE(fit_registry$prediction_success) | fit_registry$prediction_success == TRUE]
  pred_ok_ids <- unique(pred_ok_ids)
  
  if (length(pred_ok_ids) == 0L) {
    stop("No model achieved prediction_success = TRUE. Exports would be broken.")
  }
  
  if (nrow(cure_only) != length(pred_ok_ids)) {
    stop("cure_only row count does not equal the number of prediction_success models.")
  }
  
  expected_mean_rows <- length(pred_ok_ids) * length(PREDICTION_YEARS)
  if (nrow(mean_predictions) != expected_mean_rows) {
    stop("mean_predictions row count mismatch: ", nrow(mean_predictions), " vs ", expected_mean_rows)
  }
  
  if (nrow(latency_summary) != expected_mean_rows) {
    stop("latency_summary row count mismatch: ", nrow(latency_summary), " vs ", expected_mean_rows)
  }
  
  if (EXPORT_SUBJECT_PREDICTIONS) {
    fit_to_dataset <- fit_registry[, c("fit_id", "dataset"), drop = FALSE]
    dataset_sizes <- vapply(analysis_data_list, nrow, integer(1))
    expected_subject_rows <- sum(dataset_sizes[fit_to_dataset$dataset] * length(PREDICTION_YEARS) * as.integer(fit_to_dataset$fit_id %in% pred_ok_ids))
    if (nrow(subject_predictions) != expected_subject_rows) {
      stop("subject_predictions row count mismatch: ", nrow(subject_predictions), " vs ", expected_subject_rows)
    }
  }
  
  param_ok_n <- sum(fit_registry$engine == "flexsurvcure" & fit_registry$prediction_success)
  cox_ok_n <- sum(fit_registry$engine == "smcure" & fit_registry$prediction_success)
  
  if (param_ok_n == 0L) {
    stop("No parametric flexsurvcure model achieved prediction_success = TRUE.")
  }
  
  if (cox_ok_n == 0L) {
    stop("No semiparametric smcure model achieved prediction_success = TRUE.")
  }
  
  invisible(
    data.frame(
      metric = c("prediction_success_total", "prediction_success_parametric", "prediction_success_cox"),
      value = c(length(pred_ok_ids), param_ok_n, cox_ok_n),
      stringsAsFactors = FALSE
    )
  )
}

# 🔴 Execute: data ingest, grid building, fitting, and QC ===============================

## 🟠 Ingest: source data and dataset branches ===============================

core_data <- read_and_clean_input(DATA_PATH)
raw_n <- attr(core_data, "raw_n")
dropped_n <- attr(core_data, "dropped_n")
age_var_used_global <- attr(core_data, "age_var_used")
duplicate_site_id_n <- sum(duplicated(core_data$site_id))

analysis_data_list <- list()
preproc_rows <- list()

for (ds in DATASETS_TO_RUN) {
  prepared <- prepare_analysis_dataset(core_df = core_data, dataset_key = ds)
  analysis_data_list[[toupper(ds)]] <- prepared$data
  preproc_rows[[toupper(ds)]] <- prepared$preproc
}

preproc_spec <- combine_rows(preproc_rows)
model_grid <- combine_rows(lapply(DATASETS_TO_RUN, build_step5_grid))

## 🟠 Run: fit or reload master object ===============================

if (!REFIT_MODELS && !file.exists(MASTER_RDS_PATH)) {
  stop("REFIT_MODELS is FALSE, but no existing master RDS was found: ", MASTER_RDS_PATH)
}

if (!REFIT_MODELS) {
  master_obj <- readRDS(MASTER_RDS_PATH)
  
  analysis_spec <- master_obj$analysis_spec
  preproc_spec <- master_obj$preprocessing
  model_grid <- master_obj$model_grid
  fit_registry <- master_obj$fit_registry
  coefficient_table <- master_obj$coefficient_table
  cure_only <- master_obj$cure_only
  mean_predictions <- master_obj$mean_predictions
  latency_summary <- master_obj$latency_summary
  subject_predictions <- master_obj$subject_predictions
  fit_objects <- master_obj$fit_objects
} else {
  analysis_spec <- data.frame(
    script_version = SCRIPT_VERSION,
    analysis_run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    data_path = DATA_PATH,
    export_dir = EXPORT_DIR,
    refit_models = REFIT_MODELS,
    datasets_run = paste(DATASETS_TO_RUN, collapse = "|"),
    age_var_requested = AGE_VAR,
    age_var_fallback = AGE_FALLBACK_VAR,
    age_var_used_global = age_var_used_global,
    time_var = TIME_VAR,
    status_var = STATUS_VAR,
    event_code = EVENT_CODE,
    right_censor_codes = paste(RIGHT_CENSOR_CODES, collapse = "|"),
    time_unit_divisor = TIME_UNIT_DIVISOR,
    time_eps_years = TIME_EPS_YEARS,
    prediction_years = paste(PREDICTION_YEARS, collapse = "|"),
    pnu_site_value = PNU_SITE_VALUE,
    snu_site_value = SNU_SITE_VALUE,
    site_reference = SITE_REFERENCE,
    boundary_eps = BOUNDARY_EPS,
    flexsurvcure_link = FLEXSURVCURE_LINK,
    flexsurv_maxit = FLEXSURV_MAXIT,
    smcure_link = SMCURE_LINK,
    smcure_var = SMCURE_VAR,
    smcure_nboot = if (SMCURE_VAR) SMCURE_NBOOT else 0L,
    smcure_emmax = SMCURE_EMMAX,
    smcure_eps = SMCURE_EPS,
    raw_rows = raw_n,
    cleaned_rows = nrow(core_data),
    dropped_rows_global_cleaning = dropped_n,
    duplicated_site_id_n = duplicate_site_id_n,
    stop_if_any_prediction_failure = STOP_IF_ANY_PREDICTION_FAILURE,
    package_survival = safe_pkg_version("survival"),
    package_flexsurv = safe_pkg_version("flexsurv"),
    package_flexsurvcure = safe_pkg_version("flexsurvcure"),
    package_smcure = safe_pkg_version("smcure"),
    R_version = R.version.string,
    stringsAsFactors = FALSE
  )
  
  fit_objects <- vector("list", length = nrow(model_grid))
  names(fit_objects) <- model_grid$fit_id
  
  registry_rows <- vector("list", length = nrow(model_grid))
  coefficient_tables <- list()
  cure_only_rows <- list()
  mean_prediction_rows <- list()
  latency_summary_rows <- list()
  subject_prediction_rows <- list()
  
  for (i in seq_len(nrow(model_grid))) {
    spec_row <- as.list(model_grid[i, , drop = FALSE])
    data_df <- analysis_data_list[[spec_row$dataset]]
    
    if (identical(spec_row$engine, "flexsurvcure")) {
      res_i <- fit_one_parametric_cure_model(spec_row = spec_row, data_df = data_df)
    } else if (identical(spec_row$engine, "smcure")) {
      res_i <- fit_one_smcure_ph_model(spec_row = spec_row, data_df = data_df)
    } else {
      stop("Unknown engine in model grid: ", spec_row$engine)
    }
    
    fit_objects[[spec_row$fit_id]] <- res_i$fit_bundle
    registry_rows[[i]] <- res_i$registry
    
    if (nrow(res_i$coefficient_table) > 0L) {
      coefficient_tables[[length(coefficient_tables) + 1L]] <- res_i$coefficient_table
    }
    if (nrow(res_i$cure_only) > 0L) {
      cure_only_rows[[length(cure_only_rows) + 1L]] <- res_i$cure_only
    }
    if (nrow(res_i$mean_predictions) > 0L) {
      mean_prediction_rows[[length(mean_prediction_rows) + 1L]] <- res_i$mean_predictions
    }
    if (nrow(res_i$latency_summary) > 0L) {
      latency_summary_rows[[length(latency_summary_rows) + 1L]] <- res_i$latency_summary
    }
    if (nrow(res_i$subject_predictions) > 0L) {
      subject_prediction_rows[[length(subject_prediction_rows) + 1L]] <- res_i$subject_predictions
    }
  }
  
  fit_registry <- combine_rows(registry_rows)
  coefficient_table <- combine_rows(coefficient_tables, empty_df_fun = empty_coeff_df)
  cure_only <- combine_rows(cure_only_rows, empty_df_fun = empty_cure_only_df)
  mean_predictions <- combine_rows(mean_prediction_rows, empty_df_fun = empty_mean_pred_df)
  latency_summary <- combine_rows(latency_summary_rows, empty_df_fun = empty_latency_summary_df)
  subject_predictions <- combine_rows(subject_prediction_rows, empty_df_fun = empty_subject_pred_df)
  
  qc_summary <- validate_step5_outputs(
    model_grid = model_grid,
    fit_registry = fit_registry,
    cure_only = cure_only,
    mean_predictions = mean_predictions,
    latency_summary = latency_summary,
    subject_predictions = subject_predictions,
    analysis_data_list = analysis_data_list
  )
  
  if (STOP_IF_ANY_PREDICTION_FAILURE && any(!fit_registry$prediction_success)) {
    bad_ids <- fit_registry$fit_id[!fit_registry$prediction_success]
    stop(
      "At least one model did not achieve prediction_success = TRUE: ",
      paste(bad_ids, collapse = ", ")
    )
  }
  
  master_obj <- list(
    analysis_spec = analysis_spec,
    preprocessing = preproc_spec,
    model_grid = model_grid,
    fit_registry = fit_registry,
    coefficient_table = coefficient_table,
    cure_only = cure_only,
    mean_predictions = mean_predictions,
    latency_summary = latency_summary,
    subject_predictions = subject_predictions,
    analysis_data = analysis_data_list,
    fit_objects = fit_objects,
    qc_summary = qc_summary
  )
  
  saveRDS(master_obj, MASTER_RDS_PATH, compress = "xz")
}

# 🔴 Export: csv outputs and manifest ===============================

## 🟠 Order: final tables consistently ===============================

if (nrow(model_grid) > 0L) {
  model_grid <- model_grid[order(model_grid$dataset, model_grid$lane, model_grid$family, model_grid$cov_spec), , drop = FALSE]
}

if (nrow(preproc_spec) > 0L) {
  preproc_spec <- preproc_spec[order(preproc_spec$dataset), , drop = FALSE]
}

if (nrow(fit_registry) > 0L) {
  fit_registry <- fit_registry[order(fit_registry$dataset, fit_registry$lane, fit_registry$family, fit_registry$cov_spec), , drop = FALSE]
}

if (nrow(coefficient_table) > 0L) {
  coefficient_table <- coefficient_table[order(coefficient_table$dataset, coefficient_table$fit_id, coefficient_table$component, coefficient_table$term), , drop = FALSE]
}

if (nrow(cure_only) > 0L) {
  cure_only <- cure_only[order(cure_only$dataset, cure_only$fit_id), , drop = FALSE]
}

if (nrow(mean_predictions) > 0L) {
  mean_predictions <- mean_predictions[order(mean_predictions$dataset, mean_predictions$fit_id, mean_predictions$year), , drop = FALSE]
}

if (nrow(latency_summary) > 0L) {
  latency_summary <- latency_summary[order(latency_summary$dataset, latency_summary$fit_id, latency_summary$year), , drop = FALSE]
}

if (nrow(subject_predictions) > 0L) {
  subject_predictions <- subject_predictions[order(subject_predictions$dataset, subject_predictions$fit_id, subject_predictions$subject_uid, subject_predictions$year), , drop = FALSE]
}

## 🟠 Write: source-of-truth files ===============================

utils::write.csv(model_grid, MODEL_GRID_CSV_PATH, row.names = FALSE, na = "")
utils::write.csv(preproc_spec, PREPROC_CSV_PATH, row.names = FALSE, na = "")
utils::write.csv(analysis_spec, ANALYSIS_SPEC_CSV_PATH, row.names = FALSE, na = "")
utils::write.csv(fit_registry, FIT_REGISTRY_CSV_PATH, row.names = FALSE, na = "")
utils::write.csv(coefficient_table, COEFFICIENTS_CSV_PATH, row.names = FALSE, na = "")
utils::write.csv(cure_only, CURE_ONLY_CSV_PATH, row.names = FALSE, na = "")
utils::write.csv(mean_predictions, MEAN_PRED_CSV_PATH, row.names = FALSE, na = "")
utils::write.csv(latency_summary, LATENCY_SUMMARY_CSV_PATH, row.names = FALSE, na = "")

if (EXPORT_SUBJECT_PREDICTIONS) {
  utils::write.csv(subject_predictions, SUBJECT_PRED_CSV_PATH, row.names = FALSE, na = "")
}

manifest_rows <- list(
  data.frame(
    file_name = basename(MASTER_RDS_PATH),
    purpose = "Master RDS with fit objects, prepared data, exported tables, and QC summary",
    rows = NA_integer_,
    stringsAsFactors = FALSE
  ),
  data.frame(
    file_name = basename(MODEL_GRID_CSV_PATH),
    purpose = "Step5 cure-MLE model grid",
    rows = nrow(model_grid),
    stringsAsFactors = FALSE
  ),
  data.frame(
    file_name = basename(PREPROC_CSV_PATH),
    purpose = "Dataset-specific preprocessing and scaling summary",
    rows = nrow(preproc_spec),
    stringsAsFactors = FALSE
  ),
  data.frame(
    file_name = basename(ANALYSIS_SPEC_CSV_PATH),
    purpose = "Global analysis specification and package versions",
    rows = nrow(analysis_spec),
    stringsAsFactors = FALSE
  ),
  data.frame(
    file_name = basename(FIT_REGISTRY_CSV_PATH),
    purpose = "Fit registry with fit_success, prediction_success, and QC fields",
    rows = nrow(fit_registry),
    stringsAsFactors = FALSE
  ),
  data.frame(
    file_name = basename(COEFFICIENTS_CSV_PATH),
    purpose = "Incidence and latency coefficient table",
    rows = nrow(coefficient_table),
    stringsAsFactors = FALSE
  ),
  data.frame(
    file_name = basename(CURE_ONLY_CSV_PATH),
    purpose = "Cure-only summary by fitted model",
    rows = nrow(cure_only),
    stringsAsFactors = FALSE
  ),
  data.frame(
    file_name = basename(MEAN_PRED_CSV_PATH),
    purpose = "Cohort-average population survival and risk by year",
    rows = nrow(mean_predictions),
    stringsAsFactors = FALSE
  ),
  data.frame(
    file_name = basename(LATENCY_SUMMARY_CSV_PATH),
    purpose = "Latency survival summary by year",
    rows = nrow(latency_summary),
    stringsAsFactors = FALSE
  )
)

if (EXPORT_SUBJECT_PREDICTIONS) {
  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
    file_name = basename(SUBJECT_PRED_CSV_PATH),
    purpose = "Subject-level yearly predictions",
    rows = nrow(subject_predictions),
    stringsAsFactors = FALSE
  )
}

manifest <- combine_rows(manifest_rows)
utils::write.csv(manifest, MANIFEST_CSV_PATH, row.names = FALSE, na = "")