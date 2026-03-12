# 🔴 Configure: Step4 runtime and paths ===============================

## 🟠 Define: user-editable options and file names ===============================

options(stringsAsFactors = FALSE)

auto_install_packages <- FALSE


data_path <-"/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"

export_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦4.Step4_Non-cure benchmark models/attachments'


time_unit_days_per_year <- 365.25

required_packages <- c(
  "survival",
  "flexsurv"
)

allowed_sites <- c("PNU", "SNU")
allowed_status_num <- c(0, 1, 2)

dataset_labels <- c("merged", "PNU", "SNU")
specification_labels <- c("main_effects", "sex_age_interaction")
family_labels <- c("cox_ph", "weibull_aft", "lognormal_aft", "loglogistic_aft")
parametric_family_labels <- c("weibull_aft", "lognormal_aft", "loglogistic_aft")

file_standardization_constants <- file.path(export_path, "step4_standardization_constants.csv")
file_dataset_summary <- file.path(export_path, "step4_dataset_summary.csv")
file_model_summary <- file.path(export_path, "step4_model_summary.csv")
file_best_summary <- file.path(export_path, "step4_best_parametric_summary.csv")
file_coefficient_summary <- file.path(export_path, "step4_coefficients.csv")

file_rds_merged <- file.path(export_path, "step4_nc_merged.rds")
file_rds_pnu <- file.path(export_path, "step4_nc_pnu.rds")
file_rds_snu <- file.path(export_path, "step4_nc_snu.rds")

if (!dir.exists(export_path)) {
  dir.create(export_path, recursive = TRUE, showWarnings = FALSE)
}

## 🟠 Load: required libraries safely ===============================

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  if (isTRUE(auto_install_packages)) {
    install.packages(missing_packages)
  } else {
    stop(
      sprintf(
        "Missing required packages: %s. Install them first or set auto_install_packages <- TRUE.",
        paste(missing_packages, collapse = ", ")
      )
    )
  }
}

library(survival)
library(flexsurv)

# 🔴 Define: helper functions for validation and transformation ===============================

## 🟠 Build: input-check helpers and coercion utilities ===============================

assert_required_columns <- function(data, required_cols) {
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
}

trim_upper_char <- function(x) {
  out <- trimws(as.character(x))
  out[out == ""] <- NA_character_
  toupper(out)
}

trim_char <- function(x) {
  out <- trimws(as.character(x))
  out[out == ""] <- NA_character_
  out
}

coerce_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

collapse_messages <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x) == 0L) {
    return(NA_character_)
  }
  paste(x, collapse = " | ")
}

safe_deparse_formula <- function(fml) {
  paste(deparse(fml, width.cutoff = 500L), collapse = " ")
}

is_complete_nonmissing <- function(...) {
  mats <- lapply(list(...), function(x) !is.na(x))
  Reduce(`&`, mats)
}

new_empty_model_summary <- function() {
  data.frame(
    dataset = character(0),
    specification = character(0),
    family = character(0),
    model_class = character(0),
    model_role = character(0),
    formula_text = character(0),
    n = integer(0),
    n_event = integer(0),
    n_censored = integer(0),
    n_unique_event_times = integer(0),
    ok = logical(0),
    converged = logical(0),
    df_model = numeric(0),
    logLik = numeric(0),
    AIC = numeric(0),
    BIC = numeric(0),
    warnings = character(0),
    error = character(0),
    stringsAsFactors = FALSE
  )
}

new_empty_best_summary <- function() {
  data.frame(
    dataset = character(0),
    specification = character(0),
    selected_family = character(0),
    model_class = character(0),
    model_role = character(0),
    formula_text = character(0),
    n = integer(0),
    n_event = integer(0),
    n_censored = integer(0),
    n_unique_event_times = integer(0),
    logLik = numeric(0),
    AIC = numeric(0),
    BIC = numeric(0),
    aic_rank_within_spec = integer(0),
    delta_AIC = numeric(0),
    delta_BIC = numeric(0),
    selection_pool = character(0),
    ok = logical(0),
    converged = logical(0),
    warnings = character(0),
    error = character(0),
    stringsAsFactors = FALSE
  )
}

new_empty_coefficient_summary <- function() {
  data.frame(
    dataset = character(0),
    specification = character(0),
    family = character(0),
    model_class = character(0),
    model_role = character(0),
    term = character(0),
    term_role = character(0),
    estimate = numeric(0),
    std_error = numeric(0),
    statistic = numeric(0),
    p_value = numeric(0),
    exp_estimate = numeric(0),
    effect_scale = character(0),
    stringsAsFactors = FALSE
  )
}

## 🟠 Build: standardization helpers for sex and age ===============================

compute_standardization_constants <- function(data_merged_complete) {
  age_mean <- mean(data_merged_complete$age_exact_entry, na.rm = TRUE)
  age_sd <- stats::sd(data_merged_complete$age_exact_entry, na.rm = TRUE)
  sex_mean <- mean(data_merged_complete$sex_num, na.rm = TRUE)
  
  if (!is.finite(age_mean) || !is.finite(age_sd) || age_sd <= 0) {
    stop("Merged complete-case age_exact_entry has non-finite mean or non-positive SD.")
  }
  
  if (!is.finite(sex_mean)) {
    stop("Merged complete-case sex_num mean is non-finite.")
  }
  
  data.frame(
    age_mean_merged = age_mean,
    age_sd_merged = age_sd,
    age_scale_denominator = 2 * age_sd,
    sex_mean_merged = sex_mean,
    n_reference = nrow(data_merged_complete),
    stringsAsFactors = FALSE
  )
}

apply_standardization_constants <- function(data, constants_df) {
  age_mean <- constants_df$age_mean_merged[1]
  age_sd <- constants_df$age_sd_merged[1]
  sex_mean <- constants_df$sex_mean_merged[1]
  
  data$z_age <- (data$age_exact_entry - age_mean) / (2 * age_sd)
  data$c_sex <- data$sex_num - sex_mean
  data$sex_age_interaction <- data$c_sex * data$z_age
  
  data
}

# 🔴 Define: helper functions for model fitting and summaries ===============================

## 🟠 Build: formulas, distributions, and guarded fit wrapper ===============================

get_model_formula <- function(specification_label) {
  if (identical(specification_label, "main_effects")) {
    return(Surv(time_primary, event_primary) ~ c_sex + z_age)
  }
  if (identical(specification_label, "sex_age_interaction")) {
    return(Surv(time_primary, event_primary) ~ c_sex * z_age)
  }
  stop(sprintf("Unknown specification label: %s", specification_label))
}

get_flexsurv_dist <- function(family_label) {
  if (identical(family_label, "weibull_aft")) return("weibull")
  if (identical(family_label, "lognormal_aft")) return("lnorm")
  if (identical(family_label, "loglogistic_aft")) return("llogis")
  stop(sprintf("Unknown parametric family label: %s", family_label))
}

get_model_class <- function(family_label) {
  if (identical(family_label, "cox_ph")) return("coxph")
  if (family_label %in% parametric_family_labels) return("flexsurvreg")
  stop(sprintf("Unknown family label: %s", family_label))
}

get_model_role <- function(family_label) {
  if (identical(family_label, "cox_ph")) return("semiparametric_PH_benchmark")
  if (family_label %in% parametric_family_labels) return("parametric_AFT_benchmark")
  stop(sprintf("Unknown family label: %s", family_label))
}

capture_fit <- function(expr) {
  warnings_collected <- character(0)
  
  fit_obj <- withCallingHandlers(
    tryCatch(
      expr,
      error = function(e) e
    ),
    warning = function(w) {
      warnings_collected <<- c(warnings_collected, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  if (inherits(fit_obj, "error")) {
    return(list(
      ok = FALSE,
      fit = NULL,
      warnings = unique(warnings_collected),
      error = conditionMessage(fit_obj)
    ))
  }
  
  list(
    ok = TRUE,
    fit = fit_obj,
    warnings = unique(warnings_collected),
    error = NA_character_
  )
}

fit_one_step4_model <- function(data_complete, dataset_label, specification_label, family_label) {
  formula_obj <- get_model_formula(specification_label)
  formula_text <- safe_deparse_formula(formula_obj)
  model_class <- get_model_class(family_label)
  model_role <- get_model_role(family_label)
  covariate_names <- setdiff(colnames(stats::model.matrix(formula_obj, data = data_complete)), "(Intercept)")
  
  n_total <- nrow(data_complete)
  n_event <- sum(data_complete$event_primary == 1, na.rm = TRUE)
  n_censored <- sum(data_complete$event_primary == 0, na.rm = TRUE)
  n_unique_event_times <- length(unique(data_complete$time_primary[data_complete$event_primary == 1]))
  
  if (n_total == 0L) {
    return(list(
      dataset = dataset_label,
      specification = specification_label,
      family = family_label,
      model_class = model_class,
      model_role = model_role,
      formula_text = formula_text,
      covariate_names = covariate_names,
      n = n_total,
      n_event = n_event,
      n_censored = n_censored,
      n_unique_event_times = n_unique_event_times,
      ok = FALSE,
      converged = FALSE,
      warnings = character(0),
      error = "No complete-case rows available for this dataset.",
      fit = NULL
    ))
  }
  
  if (n_event == 0L) {
    return(list(
      dataset = dataset_label,
      specification = specification_label,
      family = family_label,
      model_class = model_class,
      model_role = model_role,
      formula_text = formula_text,
      covariate_names = covariate_names,
      n = n_total,
      n_event = n_event,
      n_censored = n_censored,
      n_unique_event_times = n_unique_event_times,
      ok = FALSE,
      converged = FALSE,
      warnings = character(0),
      error = "No transition events available for this dataset.",
      fit = NULL
    ))
  }
  
  fit_result <- if (identical(family_label, "cox_ph")) {
    capture_fit(
      survival::coxph(
        formula = formula_obj,
        data = data_complete,
        ties = "efron",
        x = TRUE,
        y = TRUE,
        model = TRUE,
        singular.ok = FALSE
      )
    )
  } else {
    capture_fit(
      flexsurv::flexsurvreg(
        formula = formula_obj,
        data = data_complete,
        dist = get_flexsurv_dist(family_label)
      )
    )
  }
  
  converged_flag <- FALSE
  
  if (isTRUE(fit_result$ok) && !is.null(fit_result$fit)) {
    if (inherits(fit_result$fit, "coxph")) {
      converged_flag <- is.null(fit_result$fit$fail)
    } else if (inherits(fit_result$fit, "flexsurvreg")) {
      if (!is.null(fit_result$fit$optim$convergence)) {
        converged_flag <- identical(fit_result$fit$optim$convergence, 0L) || identical(fit_result$fit$optim$convergence, 0)
      } else if (!is.null(fit_result$fit$converged)) {
        converged_flag <- isTRUE(fit_result$fit$converged)
      } else {
        converged_flag <- TRUE
      }
    }
  }
  
  list(
    dataset = dataset_label,
    specification = specification_label,
    family = family_label,
    model_class = model_class,
    model_role = model_role,
    formula_text = formula_text,
    covariate_names = covariate_names,
    n = n_total,
    n_event = n_event,
    n_censored = n_censored,
    n_unique_event_times = n_unique_event_times,
    ok = isTRUE(fit_result$ok),
    converged = converged_flag,
    warnings = fit_result$warnings,
    error = fit_result$error,
    fit = fit_result$fit
  )
}

## 🟠 Build: record flatteners and tidy summary extractors ===============================

flatten_fit_records <- function(fit_list_nested) {
  out <- list()
  
  for (dataset_label in names(fit_list_nested)) {
    for (spec_label in names(fit_list_nested[[dataset_label]])) {
      for (family_label in names(fit_list_nested[[dataset_label]][[spec_label]])) {
        out[[length(out) + 1L]] <- fit_list_nested[[dataset_label]][[spec_label]][[family_label]]
      }
    }
  }
  
  out
}

extract_model_summary_row <- function(fit_record) {
  loglik_value <- NA_real_
  aic_value <- NA_real_
  bic_value <- NA_real_
  df_model <- NA_real_
  
  if (isTRUE(fit_record$ok) && !is.null(fit_record$fit)) {
    loglik_obj <- tryCatch(stats::logLik(fit_record$fit), error = function(e) NULL)
    
    if (!is.null(loglik_obj)) {
      loglik_value <- as.numeric(loglik_obj)
      df_model <- attr(loglik_obj, "df")
    }
    
    aic_value <- tryCatch(stats::AIC(fit_record$fit), error = function(e) NA_real_)
    bic_value <- tryCatch(stats::BIC(fit_record$fit), error = function(e) NA_real_)
  }
  
  data.frame(
    dataset = fit_record$dataset,
    specification = fit_record$specification,
    family = fit_record$family,
    model_class = fit_record$model_class,
    model_role = fit_record$model_role,
    formula_text = fit_record$formula_text,
    n = fit_record$n,
    n_event = fit_record$n_event,
    n_censored = fit_record$n_censored,
    n_unique_event_times = fit_record$n_unique_event_times,
    ok = fit_record$ok,
    converged = fit_record$converged,
    df_model = df_model,
    logLik = loglik_value,
    AIC = aic_value,
    BIC = bic_value,
    warnings = collapse_messages(fit_record$warnings),
    error = fit_record$error,
    stringsAsFactors = FALSE
  )
}

extract_coefficient_table <- function(fit_record) {
  if (!isTRUE(fit_record$ok) || is.null(fit_record$fit)) {
    return(NULL)
  }
  
  if (inherits(fit_record$fit, "coxph")) {
    sx <- summary(fit_record$fit)
    coef_mat <- sx$coefficients
    
    if (is.null(coef_mat) || length(coef_mat) == 0L) {
      return(NULL)
    }
    
    if (is.null(dim(coef_mat))) {
      coef_mat <- matrix(coef_mat, nrow = 1L)
      rownames(coef_mat) <- fit_record$covariate_names[1]
      colnames(coef_mat) <- names(summary(fit_record$fit)$coefficients)
    }
    
    out <- data.frame(
      dataset = fit_record$dataset,
      specification = fit_record$specification,
      family = fit_record$family,
      model_class = fit_record$model_class,
      model_role = fit_record$model_role,
      term = rownames(coef_mat),
      term_role = "covariate",
      estimate = unname(coef_mat[, "coef"]),
      std_error = unname(coef_mat[, "se(coef)"]),
      statistic = unname(coef_mat[, "z"]),
      p_value = unname(coef_mat[, "Pr(>|z|)"]),
      exp_estimate = unname(coef_mat[, "exp(coef)"]),
      effect_scale = "hazard_ratio",
      stringsAsFactors = FALSE
    )
    
    rownames(out) <- NULL
    return(out)
  }
  
  if (inherits(fit_record$fit, "flexsurvreg")) {
    est <- tryCatch(stats::coef(fit_record$fit), error = function(e) NULL)
    vc <- tryCatch(stats::vcov(fit_record$fit), error = function(e) NULL)
    
    if (is.null(est) || length(est) == 0L) {
      return(NULL)
    }
    
    se <- rep(NA_real_, length(est))
    if (!is.null(vc)) {
      se_try <- tryCatch(sqrt(diag(vc)), error = function(e) rep(NA_real_, length(est)))
      if (length(se_try) == length(est)) {
        se <- se_try
      }
    }
    
    statistic <- est / se
    p_value <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
    term_role <- ifelse(names(est) %in% fit_record$covariate_names, "covariate", "distribution_parameter")
    exp_estimate <- ifelse(term_role == "covariate", exp(est), NA_real_)
    effect_scale <- ifelse(term_role == "covariate", "time_ratio", "not_applicable")
    
    out <- data.frame(
      dataset = fit_record$dataset,
      specification = fit_record$specification,
      family = fit_record$family,
      model_class = fit_record$model_class,
      model_role = fit_record$model_role,
      term = names(est),
      term_role = term_role,
      estimate = unname(est),
      std_error = unname(se),
      statistic = unname(statistic),
      p_value = unname(p_value),
      exp_estimate = unname(exp_estimate),
      effect_scale = effect_scale,
      stringsAsFactors = FALSE
    )
    
    rownames(out) <- NULL
    return(out)
  }
  
  NULL
}

build_best_parametric_summary <- function(model_summary_df) {
  eligible <- model_summary_df[
    model_summary_df$family %in% parametric_family_labels &
      model_summary_df$ok &
      model_summary_df$converged &
      is.finite(model_summary_df$AIC),
    ,
    drop = FALSE
  ]
  
  if (nrow(eligible) == 0L) {
    return(new_empty_best_summary())
  }
  
  split_key <- paste(eligible$dataset, eligible$specification, sep = "||")
  split_list <- split(eligible, split_key, drop = TRUE)
  
  best_rows <- lapply(split_list, function(df_one) {
    family_order <- match(df_one$family, parametric_family_labels)
    ord <- order(df_one$AIC, df_one$BIC, family_order, na.last = TRUE)
    df_ord <- df_one[ord, , drop = FALSE]
    
    best <- df_ord[1, , drop = FALSE]
    
    data.frame(
      dataset = best$dataset,
      specification = best$specification,
      selected_family = best$family,
      model_class = best$model_class,
      model_role = best$model_role,
      formula_text = best$formula_text,
      n = best$n,
      n_event = best$n_event,
      n_censored = best$n_censored,
      n_unique_event_times = best$n_unique_event_times,
      logLik = best$logLik,
      AIC = best$AIC,
      BIC = best$BIC,
      aic_rank_within_spec = 1L,
      delta_AIC = 0,
      delta_BIC = 0,
      selection_pool = "parametric_AFT_only",
      ok = best$ok,
      converged = best$converged,
      warnings = best$warnings,
      error = best$error,
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, best_rows)
  rownames(out) <- NULL
  out
}

build_fit_nc_best <- function(fit_nc_list, best_parametric_summary) {
  out <- list()
  
  for (dataset_label in dataset_labels) {
    out[[dataset_label]] <- list()
    
    for (specification_label in specification_labels) {
      hit <- best_parametric_summary[
        best_parametric_summary$dataset == dataset_label &
          best_parametric_summary$specification == specification_label,
        ,
        drop = FALSE
      ]
      
      if (nrow(hit) == 0L) {
        out[[dataset_label]][[specification_label]] <- NULL
      } else {
        selected_family <- hit$selected_family[1]
        out[[dataset_label]][[specification_label]] <- fit_nc_list[[dataset_label]][[specification_label]][[selected_family]]
      }
    }
  }
  
  out
}

# 🔴 Read: merged input and derive Step4 analysis data ===============================

## 🟠 Import: merged source file and enforce core variables ===============================

if (!file.exists(data_path)) {
  stop(sprintf("Input file not found: %s", data_path))
}

dat_raw <- utils::read.csv(
  file = data_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "site",
  "id",
  "date_entry",
  "days_followup",
  "status_num",
  "sex_num",
  "age_exact_entry"
)

assert_required_columns(dat_raw, required_columns)

## 🟠 Check: endpoint coding, site labels, and key integrity ===============================

dat_analysis <- dat_raw

dat_analysis$site <- trim_upper_char(dat_analysis$site)
dat_analysis$id <- trim_char(dat_analysis$id)
dat_analysis$date_entry <- trim_char(dat_analysis$date_entry)

dat_analysis$days_followup <- coerce_numeric(dat_analysis$days_followup)
dat_analysis$status_num <- coerce_numeric(dat_analysis$status_num)
dat_analysis$sex_num <- coerce_numeric(dat_analysis$sex_num)
dat_analysis$age_exact_entry <- coerce_numeric(dat_analysis$age_exact_entry)

if (any(is.na(dat_analysis$site) | !(dat_analysis$site %in% allowed_sites))) {
  stop(sprintf("Found invalid site values. Allowed values are: %s", paste(allowed_sites, collapse = ", ")))
}

if (any(!is.na(dat_analysis$status_num) & !(dat_analysis$status_num %in% allowed_status_num))) {
  stop(sprintf("Found invalid status_num values. Allowed values are: %s", paste(allowed_status_num, collapse = ", ")))
}

if (any(!is.na(dat_analysis$sex_num) & !(dat_analysis$sex_num %in% c(0, 1)))) {
  stop("Found invalid sex_num values. Allowed values are 0 and 1.")
}

if (any(dat_analysis$days_followup < 0, na.rm = TRUE)) {
  stop("Found negative values in days_followup.")
}

dat_analysis$site_id <- paste(dat_analysis$site, dat_analysis$id, sep = "__")

if (anyDuplicated(dat_analysis$site_id) > 0L) {
  duplicated_example <- dat_analysis$site_id[duplicated(dat_analysis$site_id)][1]
  stop(sprintf("Found duplicated site + id key. Example: %s", duplicated_example))
}

## 🟠 Transform: primary endpoint and merged-reference standardization ===============================

dat_analysis$time_primary <- dat_analysis$days_followup
dat_analysis$event_primary <- ifelse(
  dat_analysis$status_num == 1,
  1,
  ifelse(dat_analysis$status_num %in% c(0, 2), 0, NA_real_)
)
dat_analysis$status_primary <- ifelse(
  dat_analysis$event_primary == 1,
  "transition",
  ifelse(dat_analysis$event_primary == 0, "censored", NA_character_)
)

dat_analysis$step4_complete_case <- is_complete_nonmissing(
  dat_analysis$site,
  dat_analysis$id,
  dat_analysis$date_entry,
  dat_analysis$time_primary,
  dat_analysis$event_primary,
  dat_analysis$sex_num,
  dat_analysis$age_exact_entry
)

dat_analysis$analysis_row_id <- seq_len(nrow(dat_analysis))

merged_reference <- dat_analysis[dat_analysis$step4_complete_case, , drop = FALSE]

if (nrow(merged_reference) == 0L) {
  stop("No complete-case rows available in merged Step4 analysis set.")
}

standardization_constants <- compute_standardization_constants(merged_reference)
dat_analysis <- apply_standardization_constants(dat_analysis, standardization_constants)

## 🟠 Split: merged, PNU, and SNU complete-case analysis sets ===============================

dataset_data_full <- list(
  merged = dat_analysis,
  PNU = dat_analysis[dat_analysis$site == "PNU", , drop = FALSE],
  SNU = dat_analysis[dat_analysis$site == "SNU", , drop = FALSE]
)

dataset_data_complete <- list(
  merged = dat_analysis[dat_analysis$step4_complete_case, , drop = FALSE],
  PNU = dat_analysis[dat_analysis$step4_complete_case & dat_analysis$site == "PNU", , drop = FALSE],
  SNU = dat_analysis[dat_analysis$step4_complete_case & dat_analysis$site == "SNU", , drop = FALSE]
)

dataset_summary_list <- lapply(dataset_labels, function(dataset_label) {
  df_full <- dataset_data_full[[dataset_label]]
  df_cc <- dataset_data_complete[[dataset_label]]
  
  data.frame(
    dataset = dataset_label,
    n_total_subset = nrow(df_full),
    n_complete_case = nrow(df_cc),
    n_excluded_missing = nrow(df_full) - nrow(df_cc),
    n_event_complete = sum(df_cc$event_primary == 1, na.rm = TRUE),
    n_censored_complete = sum(df_cc$event_primary == 0, na.rm = TRUE),
    followup_median_days = if (nrow(df_cc) > 0L) stats::median(df_cc$time_primary, na.rm = TRUE) else NA_real_,
    followup_max_days = if (nrow(df_cc) > 0L) max(df_cc$time_primary, na.rm = TRUE) else NA_real_,
    age_mean_complete = if (nrow(df_cc) > 0L) mean(df_cc$age_exact_entry, na.rm = TRUE) else NA_real_,
    age_sd_complete = if (nrow(df_cc) > 0L) stats::sd(df_cc$age_exact_entry, na.rm = TRUE) else NA_real_,
    female_pct_complete = if (nrow(df_cc) > 0L) mean(df_cc$sex_num == 0, na.rm = TRUE) * 100 else NA_real_,
    male_pct_complete = if (nrow(df_cc) > 0L) mean(df_cc$sex_num == 1, na.rm = TRUE) * 100 else NA_real_,
    stringsAsFactors = FALSE
  )
})

dataset_summary <- do.call(rbind, dataset_summary_list)
rownames(dataset_summary) <- NULL

# 🔴 Fit: Step4 non-cure benchmark models ===============================

## 🟠 Fit: all dataset-by-specification-by-family combinations ===============================

fit_nc_list <- list()

for (dataset_label in dataset_labels) {
  fit_nc_list[[dataset_label]] <- list()
  
  for (specification_label in specification_labels) {
    fit_nc_list[[dataset_label]][[specification_label]] <- list()
    
    for (family_label in family_labels) {
      fit_nc_list[[dataset_label]][[specification_label]][[family_label]] <- fit_one_step4_model(
        data_complete = dataset_data_complete[[dataset_label]],
        dataset_label = dataset_label,
        specification_label = specification_label,
        family_label = family_label
      )
    }
  }
}

## 🟠 Summarize: tidy tables for models, coefficients, and best parametric fits ===============================

all_fit_records <- flatten_fit_records(fit_nc_list)

model_summary_rows <- lapply(all_fit_records, extract_model_summary_row)
model_summary <- if (length(model_summary_rows) > 0L) {
  do.call(rbind, model_summary_rows)
} else {
  new_empty_model_summary()
}
rownames(model_summary) <- NULL

coefficient_rows <- lapply(all_fit_records, extract_coefficient_table)
coefficient_rows <- Filter(Negate(is.null), coefficient_rows)
coefficient_summary <- if (length(coefficient_rows) > 0L) {
  do.call(rbind, coefficient_rows)
} else {
  new_empty_coefficient_summary()
}
rownames(coefficient_summary) <- NULL

best_parametric_summary <- build_best_parametric_summary(model_summary)
fit_nc_best <- build_fit_nc_best(fit_nc_list, best_parametric_summary)

# 🔴 Export: CSV tables and dataset-level RDS bundles ===============================

## 🟠 Write: source-of-truth CSV outputs ===============================

utils::write.csv(standardization_constants, file_standardization_constants, row.names = FALSE, na = "")
utils::write.csv(dataset_summary, file_dataset_summary, row.names = FALSE, na = "")
utils::write.csv(model_summary, file_model_summary, row.names = FALSE, na = "")
utils::write.csv(best_parametric_summary, file_best_summary, row.names = FALSE, na = "")
utils::write.csv(coefficient_summary, file_coefficient_summary, row.names = FALSE, na = "")

## 🟠 Save: self-contained RDS bundles for merged, PNU, and SNU ===============================

make_dataset_bundle <- function(dataset_label, rds_path) {
  dataset_bundle <- list(
    step = "Step4_non_cure_benchmark",
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    source_data_path = normalizePath(data_path, winslash = "/", mustWork = FALSE),
    dataset = dataset_label,
    endpoint_definition = list(
      time_origin = "date_entry",
      time_variable = "time_primary",
      event_variable = "event_primary",
      event_definition = "status_num == 1",
      censoring_definition = "status_num %in% c(0, 2)",
      time_unit = "days",
      years_to_days_multiplier = time_unit_days_per_year
    ),
    covariate_definition = list(
      sex_variable = "sex_num",
      age_variable = "age_exact_entry",
      sex_coding = "0 = Female, 1 = Male",
      site_in_formula = FALSE,
      age_standardization = "(age_exact_entry - mean_merged) / (2 * sd_merged)",
      sex_standardization = "sex_num - mean_merged",
      interaction_term = "c_sex * z_age",
      standardization_reference = "merged complete-case Step4 analysis set"
    ),
    standardization_constants = standardization_constants,
    data_full_subset = dataset_data_full[[dataset_label]],
    data_complete_case = dataset_data_complete[[dataset_label]],
    fit_nc_list = fit_nc_list[[dataset_label]],
    fit_nc_best = fit_nc_best[[dataset_label]],
    dataset_summary = dataset_summary[dataset_summary$dataset == dataset_label, , drop = FALSE],
    model_summary = model_summary[model_summary$dataset == dataset_label, , drop = FALSE],
    best_parametric_summary = best_parametric_summary[best_parametric_summary$dataset == dataset_label, , drop = FALSE],
    coefficient_summary = coefficient_summary[coefficient_summary$dataset == dataset_label, , drop = FALSE],
    session_info = utils::sessionInfo()
  )
  
  saveRDS(dataset_bundle, file = rds_path)
}

make_dataset_bundle("merged", file_rds_merged)
make_dataset_bundle("PNU", file_rds_pnu)
make_dataset_bundle("SNU", file_rds_snu)

## 🟠 Assign: workflow objects into current session ===============================

assign("dat_analysis", dat_analysis, envir = .GlobalEnv)
assign("fit_nc_list", fit_nc_list, envir = .GlobalEnv)
assign("fit_nc_best", fit_nc_best, envir = .GlobalEnv)
assign("step4_standardization_constants", standardization_constants, envir = .GlobalEnv)
assign("step4_dataset_summary", dataset_summary, envir = .GlobalEnv)
assign("step4_model_summary", model_summary, envir = .GlobalEnv)
assign("step4_best_parametric_summary", best_parametric_summary, envir = .GlobalEnv)
assign("step4_coefficient_summary", coefficient_summary, envir = .GlobalEnv)