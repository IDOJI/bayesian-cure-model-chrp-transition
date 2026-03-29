# 🔴 Configure: paths and controls ===============================
data_path <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/stage1_Backbone lock'
export_path <-'/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/stage5_Individualized no-cure comparator'

stage1_bundle_file <- file.path(data_path, "stage1_backbone_bundle.rds")
stage1_datasets_file <- file.path(data_path, "stage1_analysis_datasets.rds")
stage1_dataset_registry_file <- file.path(data_path, "stage1_dataset_registry.csv")
stage1_formula_registry_file <- file.path(data_path, "stage1_formula_registry.csv")
stage1_horizon_registry_file <- file.path(data_path, "stage1_horizon_registry.csv")
stage1_threshold_registry_file <- file.path(data_path, "stage1_threshold_registry.csv")
stage1_scaling_registry_file <- file.path(data_path, "stage1_scaling_registry.csv")
stage1_metadata_registry_file <- file.path(data_path, "stage1_metadata_registry.csv")

calibration_group_count <- 10L
calibration_min_supported_bins_for_regression <- 2L
prediction_horizons_for_plots <- c(1, 2, 5, 10)
probability_clip_epsilon <- 1e-06
time_origin_epsilon_year <- 1e-08

options(stringsAsFactors = FALSE, scipen = 999)

# 🔴 Build: helper functions ===============================
## 🟠 Declare: package guards ===============================
required_packages <- c(
  "dplyr",
  "tidyr",
  "purrr",
  "tibble",
  "readr",
  "ggplot2",
  "survival",
  "flexsurv"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop(
    sprintf(
      "Install required packages before running this script: %s",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(survival)
  library(flexsurv)
})

has_survAUC <- requireNamespace("survAUC", quietly = TRUE)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

clip_prob <- function(x, eps = probability_clip_epsilon) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

safe_char <- function(x) {
  as.character(x %||% NA_character_)
}

safe_num <- function(x) {
  as.numeric(x %||% NA_real_)
}

bind_rows_from_df_lists <- function(..., object_name = "input_lists") {
  out <- c(...)
  out <- purrr::compact(out)
  
  if (length(out) == 0) {
    return(tibble())
  }
  
  is_df <- vapply(out, inherits, logical(1), what = "data.frame")
  
  if (!all(is_df)) {
    bad_positions <- which(!is_df)
    stop(
      sprintf(
        "`%s` contains at least one non-data.frame object at positions: %s",
        object_name,
        paste(bad_positions, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  dplyr::bind_rows(out)
}

## 🟠 Create: bundle readers ===============================
read_stage1_bundle <- function(bundle_path, datasets_path = NULL) {
  if (!file.exists(bundle_path)) {
    stop(sprintf("Stage 1 backbone bundle not found: %s", bundle_path), call. = FALSE)
  }
  
  bundle <- readRDS(bundle_path)
  
  if (!is.list(bundle)) {
    stop("`stage1_backbone_bundle.rds` must be an R list.", call. = FALSE)
  }
  
  if (is.null(bundle$datasets) && !is.null(datasets_path) && file.exists(datasets_path)) {
    bundle$datasets <- readRDS(datasets_path)
  }
  
  bundle
}

extract_stage1_piece <- function(bundle, piece_name) {
  find_first_data_frame <- function(x) {
    if (inherits(x, "data.frame")) {
      return(tibble::as_tibble(x))
    }
    
    if (!is.list(x) || length(x) == 0L) {
      return(NULL)
    }
    
    preferred_names <- c("data", "df", "table", "registry", "value")
    
    for (nm in preferred_names) {
      if (!is.null(x[[nm]])) {
        out <- find_first_data_frame(x[[nm]])
        if (!is.null(out)) {
          return(out)
        }
      }
    }
    
    for (idx in seq_along(x)) {
      out <- find_first_data_frame(x[[idx]])
      if (!is.null(out)) {
        return(out)
      }
    }
    
    NULL
  }
  
  csv_lookup <- list(
    dataset_registry = stage1_dataset_registry_file,
    formula_registry = stage1_formula_registry_file,
    horizon_registry = stage1_horizon_registry_file,
    threshold_registry = stage1_threshold_registry_file,
    scaling_registry = stage1_scaling_registry_file,
    metadata_registry = stage1_metadata_registry_file
  )
  
  csv_path <- csv_lookup[[piece_name]] %||% NULL
  
  if (!is.null(csv_path) && file.exists(csv_path)) {
    piece <- readr::read_csv(csv_path, show_col_types = FALSE)
  } else {
    piece <- bundle$registries[[piece_name]] %||% NULL
    
    if (is.null(piece)) {
      stop(sprintf("Stage 1 bundle is missing registry `%s`.", piece_name), call. = FALSE)
    }
    
    piece <- find_first_data_frame(piece)
  }
  
  if (!inherits(piece, "data.frame")) {
    stop(
      sprintf(
        "Stage 1 registry `%s` could not be recovered as a rectangular data frame.",
        piece_name
      ),
      call. = FALSE
    )
  }
  
  list_cols <- names(piece)[vapply(piece, is.list, logical(1))]
  
  if (length(list_cols) > 0L) {
    stop(
      sprintf(
        "Stage 1 registry `%s` contains list-columns that cannot be grouped safely: %s",
        piece_name,
        paste(list_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  tibble::as_tibble(piece)
}

extract_stage1_datasets <- function(bundle) {
  datasets <- bundle$datasets %||% NULL
  if (is.null(datasets)) {
    stop("Stage 1 bundle does not contain `datasets`.", call. = FALSE)
  }
  expected_names <- c("PNU", "SNU", "merged")
  if (!all(expected_names %in% names(datasets))) {
    stop("Stage 1 datasets must contain `PNU`, `SNU`, and `merged`.", call. = FALSE)
  }
  datasets
}

extract_stage1_thresholds <- function(bundle, threshold_registry) {
  out <- bundle$config$risk_thresholds %||% NULL
  if (is.null(out)) {
    out <- threshold_registry$threshold
  }
  out <- sort(unique(as.numeric(out)))
  if (length(out) == 0 || anyNA(out) || any(out <= 0 | out >= 1)) {
    stop("Stage 1 threshold vector is invalid.", call. = FALSE)
  }
  out
}

extract_stage1_horizons <- function(bundle, horizon_registry) {
  out <- bundle$config$common_horizons_year %||% NULL
  if (is.null(out)) {
    out <- sort(unique(as.integer(horizon_registry$horizon_year)))
  }
  out <- as.integer(out)
  if (!identical(out, 1:10)) {
    stop("Stage 5 requires the Stage 1 common horizon grid to be exactly 1:10 years.", call. = FALSE)
  }
  out
}

## 🟠 Create: dataset normalizers ===============================
normalize_dataset_for_stage5 <- function(df, dataset_name) {
  required_cols <- c(
    "unique_person_id",
    "site",
    "id",
    "sex_num",
    "age_exact_entry",
    "age_s",
    "days_followup",
    "time_year",
    "event_main",
    "censor_main",
    "status_num"
  )
  
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      sprintf("[%s] Stage 1 dataset is missing required columns: %s", dataset_name, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  
  out <- df %>%
    mutate(
      unique_person_id = as.character(unique_person_id),
      site = factor(as.character(site)),
      id = as.character(id),
      sex_num = as.integer(sex_num),
      age_exact_entry = as.numeric(age_exact_entry),
      age_s = as.numeric(age_s),
      days_followup = as.numeric(days_followup),
      time_year = as.numeric(time_year),
      time_year_model = pmax(as.numeric(time_year), time_origin_epsilon_year),
      event_main = as.integer(event_main),
      censor_main = as.integer(censor_main),
      status_num = as.integer(status_num)
    )
  
  if (anyNA(out$unique_person_id) || any(out$unique_person_id == "")) {
    stop(sprintf("[%s] Invalid `unique_person_id` values detected.", dataset_name), call. = FALSE)
  }
  
  if (nrow(out) != dplyr::n_distinct(out$unique_person_id)) {
    stop(sprintf("[%s] `unique_person_id` must be unique.", dataset_name), call. = FALSE)
  }
  
  if (anyNA(out$sex_num) || any(!out$sex_num %in% c(0L, 1L))) {
    stop(sprintf("[%s] `sex_num` must be 0/1 only.", dataset_name), call. = FALSE)
  }
  
  if (anyNA(out$event_main) || any(!out$event_main %in% c(0L, 1L))) {
    stop(sprintf("[%s] `event_main` must be 0/1 only.", dataset_name), call. = FALSE)
  }
  
  if (anyNA(out$censor_main) || any(!out$censor_main %in% c(0L, 1L))) {
    stop(sprintf("[%s] `censor_main` must be 0/1 only.", dataset_name), call. = FALSE)
  }
  
  if (any(out$days_followup < 0, na.rm = TRUE) || any(out$time_year < 0, na.rm = TRUE)) {
    stop(sprintf("[%s] Negative follow-up time detected.", dataset_name), call. = FALSE)
  }
  
  out
}

assign_anchor_horizon <- function(dataset_name) {
  if (dataset_name == "PNU") {
    return(1L)
  }
  if (dataset_name %in% c("SNU", "merged")) {
    return(2L)
  }
  stop(sprintf("Unknown dataset for anchor horizon: %s", dataset_name), call. = FALSE)
}

## 🟠 Create: survival estimators ===============================
evaluate_survfit <- function(fit, times) {
  times <- as.numeric(times)
  
  if (length(times) == 0) {
    return(numeric(0))
  }
  
  if (is.null(fit) || length(fit$time) == 0) {
    return(rep(1, length(times)))
  }
  
  summary_obj <- summary(fit, times = times, extend = TRUE)
  surv_vec <- as.numeric(summary_obj$surv)
  
  if (length(surv_vec) != length(times)) {
    stop("Unexpected length mismatch while evaluating a survfit object.", call. = FALSE)
  }
  
  surv_vec
}

fit_event_survfit <- function(time_year, event_main) {
  survival::survfit(survival::Surv(time_year, event_main) ~ 1)
}

fit_censor_survfit <- function(time_year, censor_main) {
  survival::survfit(survival::Surv(time_year, censor_main) ~ 1)
}

build_ipcw_cache <- function(data, horizons_year) {
  censor_fit <- fit_censor_survfit(data$time_year, data$censor_main)
  event_fit <- fit_event_survfit(data$time_year, data$event_main)
  
  cache <- vector("list", length(horizons_year))
  names(cache) <- paste0("year_", horizons_year)
  
  for (idx in seq_along(horizons_year)) {
    horizon <- as.numeric(horizons_year[idx])
    
    g_horizon <- evaluate_survfit(censor_fit, horizon)
    g_left <- evaluate_survfit(censor_fit, pmax(data$time_year - 1e-08, 0))
    observed_km_risk <- 1 - evaluate_survfit(event_fit, horizon)
    
    case_ind <- data$event_main == 1L & data$time_year <= horizon
    control_ind <- data$time_year > horizon
    
    w_case <- ifelse(case_ind, 1 / pmax(g_left, time_origin_epsilon_year), 0)
    w_control <- ifelse(control_ind, 1 / pmax(g_horizon, time_origin_epsilon_year), 0)
    
    y_horizon <- ifelse(case_ind, 1, ifelse(control_ind, 0, NA_real_))
    w_ipcw <- ifelse(case_ind, w_case, ifelse(control_ind, w_control, NA_real_))
    
    case_count_ipcw <- sum(w_case, na.rm = TRUE)
    nonevent_count_ipcw <- sum(w_control, na.rm = TRUE)
    case_support_flag <- is.finite(case_count_ipcw) && case_count_ipcw > 0
    nonevent_support_flag <- is.finite(nonevent_count_ipcw) && nonevent_count_ipcw > 0
    binary_outcome_support_flag <- case_support_flag && nonevent_support_flag
    observed_ipcw_risk <- if ((case_count_ipcw + nonevent_count_ipcw) > 0) {
      case_count_ipcw / (case_count_ipcw + nonevent_count_ipcw)
    } else {
      NA_real_
    }
    
    cache[[idx]] <- list(
      horizon_year = horizon,
      censor_fit = censor_fit,
      event_fit = event_fit,
      censor_survival_at_horizon = as.numeric(g_horizon),
      observed_km_risk = as.numeric(observed_km_risk),
      observed_ipcw_risk = as.numeric(observed_ipcw_risk),
      case_ind = case_ind,
      control_ind = control_ind,
      y_horizon = y_horizon,
      w_case = w_case,
      w_control = w_control,
      w_ipcw = w_ipcw,
      case_count_ipcw = case_count_ipcw,
      nonevent_count_ipcw = nonevent_count_ipcw,
      known_count = sum(case_ind | control_ind),
      case_support_flag = as.logical(case_support_flag),
      nonevent_support_flag = as.logical(nonevent_support_flag),
      binary_outcome_support_flag = as.logical(binary_outcome_support_flag)
    )
  }
  
  list(
    censor_fit = censor_fit,
    event_fit = event_fit,
    horizon_cache = cache
  )
}

extract_flexsurv_matrix <- function(summary_list, target_times) {
  if (!is.list(summary_list)) {
    summary_list <- list(summary_list)
  }
  
  pred_mat <- t(
    vapply(
      summary_list,
      FUN.VALUE = numeric(length(target_times)),
      FUN = function(x) {
        xx <- as.data.frame(x)
        if ("time" %in% names(xx)) {
          xx <- xx[match(target_times, xx$time), , drop = FALSE]
        }
        if ("est" %in% names(xx)) {
          return(as.numeric(xx$est))
        }
        as.numeric(xx[[ncol(xx)]])
      }
    )
  )
  
  pred_mat
}

extract_flexsurv_vector <- function(summary_list) {
  if (!is.list(summary_list)) {
    summary_list <- list(summary_list)
  }
  
  vapply(
    summary_list,
    FUN.VALUE = numeric(1),
    FUN = function(x) {
      xx <- as.data.frame(x)
      if ("est" %in% names(xx)) {
        return(as.numeric(xx$est[1]))
      }
      as.numeric(xx[[ncol(xx)]][1])
    }
  )
}

predict_parametric_survival <- function(fit, newdata, horizons_year) {
  surv_list <- summary(
    fit,
    newdata = newdata,
    t = horizons_year,
    type = "survival",
    ci = FALSE,
    se = FALSE
  )
  
  median_list <- tryCatch(
    summary(
      fit,
      newdata = newdata,
      type = "quantile",
      quantiles = 0.5,
      ci = FALSE,
      se = FALSE
    ),
    error = function(e) NULL
  )
  
  surv_mat <- extract_flexsurv_matrix(surv_list, horizons_year)
  median_vec <- if (is.null(median_list)) {
    rep(NA_real_, nrow(newdata))
  } else {
    extract_flexsurv_vector(median_list)
  }
  
  list(
    survival = surv_mat,
    risk = 1 - surv_mat,
    median_survival = median_vec,
    score = NULL
  )
}

predict_cox_survival <- function(fit, newdata, horizons_year) {
  linear_predictor <- as.numeric(stats::predict(fit, newdata = newdata, type = "lp"))
  baseline_hazard <- survival::basehaz(fit, centered = FALSE)
  
  if (nrow(baseline_hazard) == 0) {
    cumulative_baseline <- rep(0, length(horizons_year))
  } else {
    baseline_hazard <- baseline_hazard[order(baseline_hazard$time), , drop = FALSE]
    cumulative_fun <- stats::stepfun(baseline_hazard$time, c(0, baseline_hazard$hazard), right = TRUE)
    cumulative_baseline <- pmax(cumulative_fun(horizons_year), 0)
  }
  
  survival_mat <- exp(-outer(exp(linear_predictor), cumulative_baseline, `*`))
  
  list(
    survival = survival_mat,
    risk = 1 - survival_mat,
    median_survival = rep(NA_real_, nrow(newdata)),
    score = linear_predictor
  )
}

## 🟠 Create: discrimination scorers ===============================
compute_weighted_auc_cd <- function(score, w_case, w_control) {
  valid <- !is.na(score) & ((w_case > 0) | (w_control > 0))
  
  if (sum(w_case[valid] > 0) == 0 || sum(w_control[valid] > 0) == 0) {
    return(NA_real_)
  }
  
  auc_df <- tibble(
    score = as.numeric(score[valid]),
    w_case = as.numeric(w_case[valid]),
    w_control = as.numeric(w_control[valid])
  ) %>%
    group_by(score) %>%
    summarise(
      w_case = sum(w_case, na.rm = TRUE),
      w_control = sum(w_control, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(score)
  
  cum_control_before <- c(0, head(cumsum(auc_df$w_control), -1))
  concordant_weight <- sum(auc_df$w_case * cum_control_before, na.rm = TRUE) +
    0.5 * sum(auc_df$w_case * auc_df$w_control, na.rm = TRUE)
  
  denominator <- sum(auc_df$w_case, na.rm = TRUE) * sum(auc_df$w_control, na.rm = TRUE)
  
  if (denominator <= 0) {
    return(NA_real_)
  }
  
  concordant_weight / denominator
}

compute_harrell_c <- function(time_year, event_main, risk_score) {
  concordance_df <- tibble(
    time_year = as.numeric(time_year),
    event_main = as.integer(event_main),
    risk_score = as.numeric(risk_score)
  )
  
  result <- tryCatch(
    survival::concordance(
      survival::Surv(time_year, event_main) ~ risk_score,
      data = concordance_df,
      reverse = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(result)) {
    return(list(estimate = NA_real_, se = NA_real_))
  }
  
  estimate_val <- tryCatch(
    {
      if (is.list(result) && !is.null(result[["concordance"]])) {
        as.numeric(result[["concordance"]])[1]
      } else if (!is.null(names(result)) && "concordance" %in% names(result)) {
        as.numeric(result[["concordance"]])[1]
      } else {
        as.numeric(result)[1]
      }
    },
    error = function(e) NA_real_
  )
  
  var_val <- tryCatch(
    {
      if (is.list(result) && !is.null(result[["var"]])) {
        as.numeric(result[["var"]])[1]
      } else if (!is.null(names(result)) && "var" %in% names(result)) {
        as.numeric(result[["var"]])[1]
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_
  )
  
  se_val <- if (is.finite(var_val) && var_val >= 0) sqrt(var_val) else NA_real_
  
  list(
    estimate = estimate_val,
    se = as.numeric(se_val)
  )
}

normalize_harrell_c <- function(x) {
  if (is.list(x)) {
    estimate_val <- suppressWarnings(
      as.numeric((x$estimate %||% x$concordance %||% NA_real_))[1]
    )
    se_val <- suppressWarnings(
      as.numeric((x$se %||% if (!is.null(x$var)) sqrt(as.numeric(x$var)[1]) else NA_real_))[1]
    )
  } else {
    estimate_val <- suppressWarnings(as.numeric(x)[1])
    se_val <- NA_real_
  }
  
  list(
    estimate = estimate_val,
    se = se_val
  )
}

compute_uno_c <- function(time_year, event_main, risk_score) {
  if (!has_survAUC) {
    return(list(estimate = NA_real_))
  }
  
  surv_object <- survival::Surv(time_year, event_main)
  
  result <- tryCatch(
    survAUC::UnoC(
      Surv.rsp = surv_object,
      Surv.rsp.new = surv_object,
      lpnew = as.numeric(risk_score)
    ),
    error = function(e) NA_real_
  )
  
  list(estimate = as.numeric(result))
}

## 🟠 Create: calibration scorers ===============================
safe_weighted_glm <- function(formula_obj, data, weights_vec, offset_vec = NULL) {
  warning_store <- character(0)
  
  fit <- withCallingHandlers(
    tryCatch(
      stats::glm(
        formula = formula_obj,
        family = stats::binomial(),
        data = data,
        weights = weights_vec,
        offset = offset_vec
      ),
      error = function(e) e
    ),
    warning = function(w) {
      warning_store <<- c(warning_store, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  list(
    fit = fit,
    warning_text = paste(unique(warning_store), collapse = " | ")
  )
}

compute_calibration_statistics <- function(pred_risk, ipcw_info, observed_reference_risk) {
  pred_risk <- clip_prob(pred_risk)
  
  known_df <- tibble(
    y = ipcw_info$y_horizon,
    lp = qlogis(pred_risk),
    w = ipcw_info$w_ipcw
  ) %>%
    filter(!is.na(y), !is.na(lp), is.finite(lp), !is.na(w), is.finite(w), w > 0)
  
  calibration_case_weight <- sum(known_df$w * known_df$y, na.rm = TRUE)
  calibration_nonevent_weight <- sum(known_df$w * (1 - known_df$y), na.rm = TRUE)
  
  binary_outcome_support_flag <- as.integer(
    is.finite(calibration_case_weight) &&
      is.finite(calibration_nonevent_weight) &&
      calibration_case_weight > 0 &&
      calibration_nonevent_weight > 0
  )
  
  calibration_minimum_n_flag <- as.integer(nrow(known_df) >= 10)
  
  grouped_bins_for_support <- compute_calibration_bins(
    pred_risk = pred_risk,
    ipcw_info = ipcw_info,
    group_count = calibration_group_count
  )
  
  calibration_supported_bin_count <- sum(
    grouped_bins_for_support$bin_binary_support_flag == 1L,
    na.rm = TRUE
  )
  
  calibration_supported_bin_min_flag <- as.integer(
    calibration_supported_bin_count >= calibration_min_supported_bins_for_regression
  )
  
  calibration_support_flag <- as.integer(
    binary_outcome_support_flag == 1L &&
      calibration_minimum_n_flag == 1L &&
      calibration_supported_bin_min_flag == 1L
  )
  
  mean_predicted_risk_value <- mean(pred_risk, na.rm = TRUE)
  mean_calibration_difference_value <- if (
    binary_outcome_support_flag == 1L &&
    is.finite(observed_reference_risk)
  ) {
    observed_reference_risk - mean_predicted_risk_value
  } else {
    NA_real_
  }
  
  out <- list(
    mean_predicted_risk = mean_predicted_risk_value,
    mean_calibration_difference = mean_calibration_difference_value,
    calibration_intercept_offset = NA_real_,
    calibration_intercept_free = NA_real_,
    calibration_slope = NA_real_,
    calibration_known_n = nrow(known_df),
    calibration_case_weight = calibration_case_weight,
    calibration_nonevent_weight = calibration_nonevent_weight,
    binary_outcome_support_flag = binary_outcome_support_flag,
    calibration_minimum_n_flag = calibration_minimum_n_flag,
    calibration_supported_bin_count = as.integer(calibration_supported_bin_count),
    calibration_supported_bin_min_flag = calibration_supported_bin_min_flag,
    calibration_support_flag = calibration_support_flag,
    calibration_fit_success_flag = 0L,
    calibration_regression_suppressed_flag = 0L,
    calibration_warning = NA_character_
  )
  
  if (calibration_support_flag != 1L) {
    return(out)
  }
  
  intercept_fit <- safe_weighted_glm(
    formula_obj = y ~ 1,
    data = known_df,
    weights_vec = known_df$w,
    offset_vec = known_df$lp
  )
  
  slope_fit <- safe_weighted_glm(
    formula_obj = y ~ lp,
    data = known_df,
    weights_vec = known_df$w,
    offset_vec = NULL
  )
  
  intercept_offset <- if (inherits(intercept_fit$fit, "error")) {
    NA_real_
  } else {
    as.numeric(stats::coef(intercept_fit$fit)[1])
  }
  
  intercept_free <- if (inherits(slope_fit$fit, "error")) {
    NA_real_
  } else {
    as.numeric(stats::coef(slope_fit$fit)[1])
  }
  
  slope_estimate <- if (inherits(slope_fit$fit, "error")) {
    NA_real_
  } else {
    as.numeric(stats::coef(slope_fit$fit)[2])
  }
  
  fit_warning <- paste(
    unique(c(intercept_fit$warning_text, slope_fit$warning_text)),
    collapse = " | "
  )
  if (!nzchar(trimws(fit_warning))) {
    fit_warning <- NA_character_
  }
  
  calibration_fit_success_flag <- as.integer(
    is.finite(intercept_offset) &&
      is.finite(intercept_free) &&
      is.finite(slope_estimate)
  )
  
  out$calibration_intercept_offset <- intercept_offset
  out$calibration_intercept_free <- intercept_free
  out$calibration_slope <- slope_estimate
  out$calibration_fit_success_flag <- calibration_fit_success_flag
  out$calibration_regression_suppressed_flag <- as.integer(calibration_fit_success_flag != 1L)
  out$calibration_warning <- fit_warning
  
  if (calibration_fit_success_flag != 1L) {
    out$calibration_intercept_offset <- NA_real_
    out$calibration_intercept_free <- NA_real_
    out$calibration_slope <- NA_real_
  }
  
  out
}

derive_calibration_failure_reason <- function(calibration_stats) {
  support_flag <- isTRUE((calibration_stats$calibration_support_flag %||% 0L) == 1L)
  fit_success_flag <- isTRUE((calibration_stats$calibration_fit_success_flag %||% 0L) == 1L)
  
  if (support_flag && fit_success_flag) {
    return("none")
  }
  
  reasons <- character(0)
  
  if (!isTRUE((calibration_stats$binary_outcome_support_flag %||% 0L) == 1L)) {
    reasons <- c(reasons, "binary_outcome_not_supported")
  }
  
  if (!isTRUE((calibration_stats$calibration_minimum_n_flag %||% 1L) == 1L)) {
    reasons <- c(reasons, "known_n_below_minimum")
  }
  
  if (!isTRUE((calibration_stats$calibration_supported_bin_min_flag %||% 1L) == 1L)) {
    reasons <- c(reasons, "supported_bins_below_minimum")
  }
  
  if (support_flag && !fit_success_flag) {
    reasons <- c(reasons, "calibration_regression_fit_failed")
  }
  
  if (length(reasons) == 0L) {
    reasons <- c(reasons, "unspecified_failure")
  }
  
  paste(unique(reasons), collapse = "; ")
}

compute_calibration_bins <- function(pred_risk, ipcw_info, group_count = calibration_group_count) {
  pred_risk <- clip_prob(pred_risk)
  unique_pred_n <- dplyr::n_distinct(pred_risk)
  overall_binary_support_flag <- if (!is.null(ipcw_info$binary_outcome_support_flag)) {
    isTRUE(ipcw_info$binary_outcome_support_flag)
  } else {
    sum(ipcw_info$w_case, na.rm = TRUE) > 0 && sum(ipcw_info$w_control, na.rm = TRUE) > 0
  }
  
  if (unique_pred_n <= 1) {
    bin_factor <- factor(rep("All", length(pred_risk)))
  } else {
    probs <- seq(0, 1, length.out = min(group_count, unique_pred_n) + 1)
    breaks_vec <- unique(stats::quantile(pred_risk, probs = probs, na.rm = TRUE, type = 2))
    if (length(breaks_vec) <= 2) {
      bin_factor <- factor(rep("All", length(pred_risk)))
    } else {
      bin_factor <- cut(pred_risk, breaks = breaks_vec, include.lowest = TRUE, dig.lab = 12)
    }
  }
  
  tibble(
    bin_label = bin_factor,
    predicted_risk = pred_risk,
    w_case = ipcw_info$w_case,
    w_control = ipcw_info$w_control
  ) %>%
    group_by(bin_label) %>%
    summarise(
      n_subject = n(),
      mean_predicted_risk = mean(predicted_risk, na.rm = TRUE),
      observed_risk_ipcw = {
        denom <- sum(w_case + w_control, na.rm = TRUE)
        if (overall_binary_support_flag && denom > 0) sum(w_case, na.rm = TRUE) / denom else NA_real_
      },
      case_count_ipcw = sum(w_case, na.rm = TRUE),
      nonevent_count_ipcw = sum(w_control, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      bin_index = row_number(),
      case_support_flag = as.integer(case_count_ipcw > 0),
      nonevent_support_flag = as.integer(nonevent_count_ipcw > 0),
      bin_binary_support_flag = as.integer(case_count_ipcw > 0 & nonevent_count_ipcw > 0),
      binary_outcome_support_flag = as.integer(overall_binary_support_flag)
    )
}

## 🟠 Create: threshold scorers ===============================
compute_threshold_vector <- function(pred_positive, ipcw_info, threshold, n_total) {
  pred_positive <- as.logical(pred_positive)
  
  tp <- sum(ipcw_info$w_case[pred_positive], na.rm = TRUE)
  fp <- sum(ipcw_info$w_control[pred_positive], na.rm = TRUE)
  fn <- sum(ipcw_info$w_case[!pred_positive], na.rm = TRUE)
  tn <- sum(ipcw_info$w_control[!pred_positive], na.rm = TRUE)
  
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  ppv <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  npv <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
  fdp <- if ((tp + fp) > 0) fp / (tp + fp) else NA_real_
  prevalence <- (tp + fn) / n_total
  net_benefit <- (tp / n_total) - (fp / n_total) * (threshold / (1 - threshold))
  net_benefit_all <- prevalence - (1 - prevalence) * (threshold / (1 - threshold))
  net_reduction_unnecessary_per_100 <- 100 * (net_benefit - net_benefit_all) / (threshold / (1 - threshold))
  
  c(
    positive_classification_rate = mean(pred_positive),
    true_positive_count_ipcw = tp,
    false_positive_count_ipcw = fp,
    false_negative_count_ipcw = fn,
    true_negative_count_ipcw = tn,
    false_positive_burden_all = fp / n_total,
    false_positive_burden_nonevents = if ((tn + fp) > 0) fp / (tn + fp) else NA_real_,
    false_positive_per_100 = 100 * fp / n_total,
    sensitivity = sensitivity,
    specificity = specificity,
    ppv = ppv,
    npv = npv,
    false_discovery_proportion = fdp,
    event_prevalence_ipcw = prevalence,
    net_benefit = net_benefit,
    net_benefit_treat_all = net_benefit_all,
    net_reduction_unnecessary_per_100 = net_reduction_unnecessary_per_100
  )
}

## 🟠 Create: export helpers ===============================
horizon_label_lookup <- function(dataset_name, horizon_year, horizon_registry) {
  horizon_registry %>%
    filter(dataset == dataset_name, horizon_year == !!horizon_year) %>%
    slice(1) %>%
    select(horizon_id, interpretation_tier, primary_supported_flag, interpretation_note) %>%
    mutate(
      horizon_year = as.integer(horizon_year),
      horizon_label = paste0("Year ", horizon_year)
    )
}

make_metric_rows <- function(meta_row, horizon_info, threshold_value, metric_domain, metrics_named) {
  tibble(
    stage = "Stage 5",
    dataset = meta_row$dataset %||% NA_character_,
    model_class = meta_row$model_class %||% NA_character_,
    model_family = meta_row$model_family %||% NA_character_,
    model_id = meta_row$model_id %||% NA_character_,
    strategy_family = meta_row$strategy_family %||% NA_character_,
    formula_id = meta_row$formula_id %||% NA_character_,
    formula_name = meta_row$formula_name %||% NA_character_,
    formula_label = meta_row$formula_label %||% NA_character_,
    site_branch = meta_row$site_branch %||% NA_character_,
    interaction_branch = meta_row$interaction_branch %||% NA_character_,
    anchor_horizon_year = meta_row$anchor_horizon_year %||% NA_integer_,
    horizon_id = horizon_info$horizon_id %||% NA_character_,
    horizon_year = horizon_info$horizon_year %||% NA_integer_,
    horizon_label = horizon_info$horizon_label %||% NA_character_,
    interpretation_tier = horizon_info$interpretation_tier %||% NA_character_,
    primary_supported_flag = horizon_info$primary_supported_flag %||% NA,
    interpretation_note = horizon_info$interpretation_note %||% NA_character_,
    threshold = threshold_value,
    metric_domain = metric_domain,
    metric_name = names(metrics_named),
    metric_value = as.numeric(unname(metrics_named))
  )
}

## 🟠 Create: model utilities ===============================
fit_single_noncure_model <- function(data, formula_rhs, model_family, flexsurv_dist = NA_character_) {
  surv_formula <- stats::as.formula(
    paste0("survival::Surv(time_year_model, event_main) ~ ", formula_rhs)
  )
  
  if (model_family == "coxph") {
    fit_object <- tryCatch(
      survival::coxph(
        formula = surv_formula,
        data = data,
        ties = "breslow",
        x = TRUE,
        y = TRUE,
        model = TRUE,
        singular.ok = TRUE
      ),
      error = function(e) e
    )
    
    if (inherits(fit_object, "error")) {
      return(list(
        fit = NULL,
        converged = FALSE,
        convergence_code = NA_integer_,
        error_message = conditionMessage(fit_object)
      ))
    }
    
    return(list(
      fit = fit_object,
      converged = is.null(fit_object$fail),
      convergence_code = if (is.null(fit_object$fail)) 0L else 1L,
      error_message = safe_char(fit_object$fail)
    ))
  }
  
  fit_object <- tryCatch(
    flexsurv::flexsurvreg(
      formula = surv_formula,
      data = data,
      dist = flexsurv_dist
    ),
    error = function(e) e
  )
  
  if (inherits(fit_object, "error")) {
    return(list(
      fit = NULL,
      converged = FALSE,
      convergence_code = NA_integer_,
      error_message = conditionMessage(fit_object)
    ))
  }
  
  conv_code <- fit_object$opt$convergence %||% fit_object$optim$convergence %||% NA_integer_
  conv_flag <- if (is.na(conv_code)) TRUE else conv_code == 0L
  
  list(
    fit = fit_object,
    converged = conv_flag,
    convergence_code = conv_code,
    error_message = NA_character_
  )
}

extract_fit_statistics <- function(fit_object, data, model_family, convergence_code, converged, error_message) {
  if (is.null(fit_object)) {
    return(tibble(
      converged = FALSE,
      convergence_code = convergence_code,
      error_message = error_message,
      n_obs = nrow(data),
      n_event = sum(data$event_main, na.rm = TRUE),
      n_parameter = NA_integer_,
      logLik = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      extrapolation_rule = ifelse(model_family == "coxph", "baseline_carried_forward_after_last_event", "parametric_extrapolation")
    ))
  }
  
  n_parameter <- length(stats::coef(fit_object))
  loglik_value <- tryCatch(as.numeric(stats::logLik(fit_object))[1], error = function(e) NA_real_)
  aic_value <- tryCatch(stats::AIC(fit_object), error = function(e) NA_real_)
  bic_value <- if (!is.na(loglik_value)) {
    -2 * loglik_value + log(nrow(data)) * n_parameter
  } else {
    NA_real_
  }
  
  tibble(
    converged = as.logical(converged),
    convergence_code = as.integer(convergence_code),
    error_message = safe_char(error_message),
    n_obs = nrow(data),
    n_event = sum(data$event_main, na.rm = TRUE),
    n_parameter = as.integer(n_parameter),
    logLik = loglik_value,
    AIC = as.numeric(aic_value),
    BIC = as.numeric(bic_value),
    extrapolation_rule = ifelse(model_family == "coxph", "baseline_carried_forward_after_last_event", "parametric_extrapolation")
  )
}

# 🔴 Load: stage-one artifacts ===============================
## 🟠 Read: backbone bundle ===============================
stage1_bundle <- read_stage1_bundle(
  bundle_path = stage1_bundle_file,
  datasets_path = stage1_datasets_file
)

## 🟠 Extract: analysis objects ===============================
formula_registry <- extract_stage1_piece(stage1_bundle, "formula_registry")
dataset_registry <- extract_stage1_piece(stage1_bundle, "dataset_registry")
horizon_registry <- extract_stage1_piece(stage1_bundle, "horizon_registry")
threshold_registry <- extract_stage1_piece(stage1_bundle, "threshold_registry")
scaling_registry <- extract_stage1_piece(stage1_bundle, "scaling_registry")
metadata_registry <- extract_stage1_piece(stage1_bundle, "metadata_registry")
analysis_datasets <- extract_stage1_datasets(stage1_bundle)

analysis_datasets <- lapply(
  names(analysis_datasets),
  function(dataset_name) normalize_dataset_for_stage5(analysis_datasets[[dataset_name]], dataset_name)
)
names(analysis_datasets) <- c("PNU", "SNU", "merged")

common_horizons_year <- extract_stage1_horizons(stage1_bundle, horizon_registry)
risk_thresholds <- extract_stage1_thresholds(stage1_bundle, threshold_registry)
prediction_horizons_for_plots <- intersect(common_horizons_year, prediction_horizons_for_plots)

dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# 🔴 Check results ===============================
cat("stage1_bundle_file:\n", normalizePath(stage1_bundle_file, winslash = "/", mustWork = FALSE), "\n\n")
cat("stage1_datasets_file:\n", normalizePath(stage1_datasets_file, winslash = "/", mustWork = FALSE), "\n\n")

cat("formula_registry dataset raw:\n")
print(unique(formula_registry$dataset))

cat("\nformula_registry dataset normalized:\n")
print(sort(unique(trimws(as.character(formula_registry$dataset)))))

cat("\nformula count by dataset:\n")
print(as.data.frame(table(trimws(as.character(formula_registry$dataset)), useNA = "ifany")))

cat("\nstage1_bundle datasets names:\n")
print(names(stage1_bundle$datasets))

# 🔴 Validate: frozen inputs ===============================
## 🟠 Check: stage-one consistency ===============================
if (!all(c("PNU", "SNU", "merged") %in% dataset_registry$dataset)) {
  stop("Stage 1 dataset registry must contain PNU, SNU, and merged.", call. = FALSE)
}

normalize_dataset_label <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    toupper(x) == "PNU" ~ "PNU",
    toupper(x) == "SNU" ~ "SNU",
    tolower(x) == "merged" ~ "merged",
    TRUE ~ x
  )
}

required_datasets <- c("PNU", "SNU", "merged")

formula_registry <- formula_registry %>%
  dplyr::mutate(dataset = normalize_dataset_label(dataset))

if (!setequal(unique(formula_registry$dataset), required_datasets)) {
  stop(
    sprintf(
      "Stage 1 formula registry datasets mismatch. observed = [%s]",
      paste(sort(unique(formula_registry$dataset)), collapse = ", ")
    ),
    call. = FALSE
  )
}

observed_formula_counts <- table(formula_registry$dataset)
expected_formula_counts <- c(PNU = 2L, SNU = 2L, merged = 4L)

if (!all(names(expected_formula_counts) %in% names(observed_formula_counts))) {
  stop("Stage 1 formula registry is missing one or more required dataset blocks.", call. = FALSE)
}

if (any(as.integer(observed_formula_counts[names(expected_formula_counts)]) < expected_formula_counts)) {
  stop(
    sprintf(
      "Stage 1 formula registry has incomplete formula rows. observed counts = [%s]",
      paste(sprintf("%s:%s", names(observed_formula_counts), as.integer(observed_formula_counts)), collapse = ", ")
    ),
    call. = FALSE
  )
}

horizon_year_grid <- sort(unique(as.integer(horizon_registry$horizon_year)))

if (anyNA(horizon_year_grid) || !identical(horizon_year_grid, 1:10)) {
  stop(
    sprintf(
      "Stage 1 horizon registry is not fixed to 1:10 years. observed = [%s]",
      paste(horizon_year_grid, collapse = ", ")
    ),
    call. = FALSE
  )
}

if (!all(risk_thresholds %in% threshold_registry$threshold)) {
  stop("Stage 1 threshold registry and extracted threshold vector are inconsistent.", call. = FALSE)
}

scaling_check <- scaling_registry %>%
  filter(variable == "age_exact_entry", scaled_variable == "age_s")

if (nrow(scaling_check) < 3) {
  stop("Stage 1 scaling registry must contain dataset-specific age_s scaling rows.", call. = FALSE)
}

## 🟠 Check: modeling plan ===============================
family_registry <- tibble(
  model_family = c("exponential", "weibull", "lognormal", "loglogistic", "coxph"),
  fit_engine = c("flexsurv", "flexsurv", "flexsurv", "flexsurv", "survival"),
  flexsurv_dist = c("exp", "weibull", "lnorm", "llogis", NA_character_),
  model_class = "non_cure"
)

model_plan <- tidyr::crossing(
  formula_registry %>%
    select(
      dataset,
      formula_id,
      formula_name,
      formula_label,
      formula_rhs,
      site_branch,
      interaction_branch
    ),
  family_registry
) %>%
  mutate(
    model_id = paste(dataset, formula_name, model_family, sep = "__"),
    anchor_horizon_year = vapply(dataset, assign_anchor_horizon, integer(1)),
    stage = "Stage 5"
  ) %>%
  arrange(match(dataset, c("PNU", "SNU", "merged")), formula_name, model_family)

if (nrow(model_plan) == 0) {
  stop("No Stage 5 non-cure models were created in the model plan.", call. = FALSE)
}

# 🔴 Fit: non-cure models ===============================
## 🟠 Prepare: dataset references ===============================
dataset_reference_rows <- list()
dataset_threshold_reference_rows <- list()
dataset_ipcw_registry <- list()

for (dataset_name in c("PNU", "SNU", "merged")) {
  dataset_data <- analysis_datasets[[dataset_name]]
  dataset_ipcw <- build_ipcw_cache(dataset_data, common_horizons_year)
  dataset_ipcw_registry[[dataset_name]] <- dataset_ipcw
  
  observed_horizon_rows <- list()
  observed_threshold_rows <- list()
  
  for (horizon_year in common_horizons_year) {
    horizon_key <- paste0("year_", horizon_year)
    ipcw_info <- dataset_ipcw$horizon_cache[[horizon_key]]
    horizon_info <- horizon_label_lookup(dataset_name, horizon_year, horizon_registry)
    
    observed_meta <- list(
      dataset = dataset_name,
      model_class = "benchmark",
      model_family = "observed_km",
      model_id = paste(dataset_name, "observed_km", sep = "__"),
      strategy_family = "observed_km",
      formula_id = NA_character_,
      formula_name = NA_character_,
      formula_label = NA_character_,
      site_branch = NA_character_,
      interaction_branch = NA_character_,
      anchor_horizon_year = assign_anchor_horizon(dataset_name)
    )
    
    observed_metrics <- c(
      observed_km_risk = ipcw_info$observed_km_risk,
      observed_ipcw_risk = ipcw_info$observed_ipcw_risk,
      case_count_ipcw = ipcw_info$case_count_ipcw,
      nonevent_count_ipcw = ipcw_info$nonevent_count_ipcw,
      known_count = ipcw_info$known_count,
      cohort_n = nrow(dataset_data)
    )
    
    observed_horizon_rows[[horizon_key]] <- make_metric_rows(
      meta_row = observed_meta,
      horizon_info = horizon_info,
      threshold_value = NA_real_,
      metric_domain = "horizon_reference",
      metrics_named = observed_metrics
    )
    
    for (threshold_value in risk_thresholds) {
      treat_all_meta <- observed_meta
      treat_all_meta$model_class <- "reference"
      treat_all_meta$model_family <- "treat_all"
      treat_all_meta$model_id <- paste(dataset_name, "treat_all", sep = "__")
      treat_all_meta$strategy_family <- "treat_all"
      
      treat_none_meta <- observed_meta
      treat_none_meta$model_class <- "reference"
      treat_none_meta$model_family <- "treat_none"
      treat_none_meta$model_id <- paste(dataset_name, "treat_none", sep = "__")
      treat_none_meta$strategy_family <- "treat_none"
      
      treat_all_metrics <- compute_threshold_vector(
        pred_positive = rep(TRUE, nrow(dataset_data)),
        ipcw_info = ipcw_info,
        threshold = threshold_value,
        n_total = nrow(dataset_data)
      )
      
      treat_none_metrics <- compute_threshold_vector(
        pred_positive = rep(FALSE, nrow(dataset_data)),
        ipcw_info = ipcw_info,
        threshold = threshold_value,
        n_total = nrow(dataset_data)
      )
      
      observed_threshold_rows[[paste0(horizon_key, "_all_", threshold_value)]] <- make_metric_rows(
        meta_row = treat_all_meta,
        horizon_info = horizon_info,
        threshold_value = threshold_value,
        metric_domain = "threshold_reference",
        metrics_named = treat_all_metrics
      )
      
      observed_threshold_rows[[paste0(horizon_key, "_none_", threshold_value)]] <- make_metric_rows(
        meta_row = treat_none_meta,
        horizon_info = horizon_info,
        threshold_value = threshold_value,
        metric_domain = "threshold_reference",
        metrics_named = treat_none_metrics
      )
    }
  }
  
  dataset_reference_rows[[dataset_name]] <- bind_rows(observed_horizon_rows)
  dataset_threshold_reference_rows[[dataset_name]] <- bind_rows(observed_threshold_rows)
}

## 🟠 Iterate: model fitting ===============================
fitted_models <- list()
model_registry_rows <- list()
risk_long_rows <- list()
performance_rows <- list()
calibration_bin_rows <- list()
calibration_diagnostic_rows <- list()

for (plan_idx in seq_len(nrow(model_plan))) {
  plan_row <- model_plan[plan_idx, , drop = FALSE]
  dataset_name <- plan_row$dataset[[1]]
  formula_rhs <- plan_row$formula_rhs[[1]]
  model_family <- plan_row$model_family[[1]]
  flexsurv_dist <- plan_row$flexsurv_dist[[1]]
  
  dataset_data <- analysis_datasets[[dataset_name]]
  dataset_ipcw <- dataset_ipcw_registry[[dataset_name]]
  anchor_horizon_year <- plan_row$anchor_horizon_year[[1]]
  anchor_horizon_index <- match(anchor_horizon_year, common_horizons_year)
  
  fit_result <- fit_single_noncure_model(
    data = dataset_data,
    formula_rhs = formula_rhs,
    model_family = model_family,
    flexsurv_dist = flexsurv_dist
  )
  
  fit_stats <- extract_fit_statistics(
    fit_object = fit_result$fit,
    data = dataset_data,
    model_family = model_family,
    convergence_code = fit_result$convergence_code,
    converged = fit_result$converged,
    error_message = fit_result$error_message
  )
  
  registry_row <- plan_row %>%
    mutate(
      stage = "Stage 5",
      analysis_time_variable = "days_followup",
      reporting_time_variable = "time_year",
      event_definition = "status_num == 1",
      censoring_definition = "status_num %in% c(0, 2)"
    ) %>%
    bind_cols(fit_stats)
  
  model_registry_rows[[plan_row$model_id[[1]]]] <- registry_row
  
  if (!isTRUE(fit_result$converged) || is.null(fit_result$fit)) {
    next
  }
  
  if (model_family == "coxph") {
    pred_result <- predict_cox_survival(
      fit = fit_result$fit,
      newdata = dataset_data,
      horizons_year = common_horizons_year
    )
  } else {
    pred_result <- tryCatch(
      predict_parametric_survival(
        fit = fit_result$fit,
        newdata = dataset_data,
        horizons_year = common_horizons_year
      ),
      error = function(e) {
        NULL
      }
    )
    
    if (is.null(pred_result)) {
      model_registry_rows[[plan_row$model_id[[1]]]] <- model_registry_rows[[plan_row$model_id[[1]]]] %>%
        mutate(
          converged = FALSE,
          error_message = "Prediction extraction failed after successful model fit."
        )
      next
    }
  }
  
  fitted_models[[plan_row$model_id[[1]]]] <- fit_result$fit
  
  prediction_long <- tibble(
    stage = "Stage 5",
    dataset = dataset_name,
    model_class = "non_cure",
    model_family = model_family,
    model_id = plan_row$model_id[[1]],
    formula_id = plan_row$formula_id[[1]],
    formula_name = plan_row$formula_name[[1]],
    formula_label = plan_row$formula_label[[1]],
    site_branch = plan_row$site_branch[[1]],
    interaction_branch = plan_row$interaction_branch[[1]],
    unique_person_id = rep(dataset_data$unique_person_id, each = length(common_horizons_year)),
    site = rep(as.character(dataset_data$site), each = length(common_horizons_year)),
    id = rep(dataset_data$id, each = length(common_horizons_year)),
    sex_num = rep(dataset_data$sex_num, each = length(common_horizons_year)),
    age_exact_entry = rep(dataset_data$age_exact_entry, each = length(common_horizons_year)),
    age_s = rep(dataset_data$age_s, each = length(common_horizons_year)),
    horizon_year = rep(common_horizons_year, times = nrow(dataset_data)),
    survival_pred = as.vector(t(pred_result$survival)),
    risk_pred = as.vector(t(pred_result$risk)),
    median_survival_pred = rep(pred_result$median_survival, each = length(common_horizons_year)),
    anchor_horizon_year = anchor_horizon_year,
    c_index_score = rep(pred_result$risk[, anchor_horizon_index], each = length(common_horizons_year))
  ) %>%
    left_join(
      horizon_registry %>%
        select(dataset, horizon_year, horizon_id, interpretation_tier, primary_supported_flag, interpretation_note),
      by = c("dataset", "horizon_year")
    )
  
  risk_long_rows[[plan_row$model_id[[1]]]] <- prediction_long
  
  harrell_c <- compute_harrell_c(
    time_year = dataset_data$time_year,
    event_main = dataset_data$event_main,
    risk_score = pred_result$risk[, anchor_horizon_index]
  )
  
  harrell_c <- normalize_harrell_c(harrell_c)
  harrell_c_estimate <- harrell_c[["estimate"]]
  harrell_c_se <- harrell_c[["se"]]
  
  uno_c <- compute_uno_c(
    time_year = dataset_data$time_year,
    event_main = dataset_data$event_main,
    risk_score = pred_result$risk[, anchor_horizon_index]
  )
  
  cindex_meta <- list(
    dataset = dataset_name,
    model_class = "non_cure",
    model_family = model_family,
    model_id = plan_row$model_id[[1]],
    strategy_family = model_family,
    formula_id = plan_row$formula_id[[1]],
    formula_name = plan_row$formula_name[[1]],
    formula_label = plan_row$formula_label[[1]],
    site_branch = plan_row$site_branch[[1]],
    interaction_branch = plan_row$interaction_branch[[1]],
    anchor_horizon_year = anchor_horizon_year
  )
  
  cindex_horizon_info <- horizon_label_lookup(dataset_name, anchor_horizon_year, horizon_registry)
  
  performance_rows[[paste0(plan_row$model_id[[1]], "_cindex")]] <- bind_rows(
    make_metric_rows(
      meta_row = cindex_meta,
      horizon_info = cindex_horizon_info,
      threshold_value = NA_real_,
      metric_domain = "time_range_discrimination",
      metrics_named = c(
        harrell_c = harrell_c_estimate,
        harrell_c_se = harrell_c_se,
        uno_c = uno_c$estimate
      )
    )
  )
  
  brier_by_horizon <- numeric(length(common_horizons_year))
  brier_by_horizon[] <- NA_real_
  
  for (h_idx in seq_along(common_horizons_year)) {
    horizon_year <- common_horizons_year[h_idx]
    horizon_key <- paste0("year_", horizon_year)
    horizon_info <- horizon_label_lookup(dataset_name, horizon_year, horizon_registry)
    ipcw_info <- dataset_ipcw$horizon_cache[[horizon_key]]
    pred_risk_h <- pred_result$risk[, h_idx]
    pred_risk_clip <- clip_prob(pred_risk_h)
    
    auc_est <- compute_weighted_auc_cd(
      score = pred_risk_h,
      w_case = ipcw_info$w_case,
      w_control = ipcw_info$w_control
    )
    
    y_known <- ipcw_info$y_horizon
    w_known <- ipcw_info$w_ipcw
    known_index <- !is.na(y_known) & !is.na(w_known) & w_known > 0
    
    brier_score <- if (any(known_index)) {
      sum(w_known[known_index] * (y_known[known_index] - pred_risk_h[known_index])^2, na.rm = TRUE) / nrow(dataset_data)
    } else {
      NA_real_
    }
    
    brier_null <- if (any(known_index)) {
      sum(w_known[known_index] * (y_known[known_index] - ipcw_info$observed_km_risk)^2, na.rm = TRUE) / nrow(dataset_data)
    } else {
      NA_real_
    }
    
    scaled_brier <- if (!is.na(brier_null) && brier_null > 0) {
      1 - (brier_score / brier_null)
    } else {
      NA_real_
    }
    
    brier_by_horizon[h_idx] <- brier_score
    
    calibration_stats <- compute_calibration_statistics(
      pred_risk = pred_risk_h,
      ipcw_info = ipcw_info,
      observed_reference_risk = ipcw_info$observed_km_risk
    )
    
    horizon_metric_values <- c(
      mean_predicted_risk = calibration_stats$mean_predicted_risk,
      observed_km_risk = ipcw_info$observed_km_risk,
      observed_ipcw_risk = ipcw_info$observed_ipcw_risk,
      mean_calibration_difference = calibration_stats$mean_calibration_difference,
      calibration_intercept_offset = calibration_stats$calibration_intercept_offset,
      calibration_intercept_free = calibration_stats$calibration_intercept_free,
      calibration_slope = calibration_stats$calibration_slope,
      binary_outcome_support_flag = calibration_stats$binary_outcome_support_flag,
      calibration_known_n = calibration_stats$calibration_known_n,
      calibration_case_weight = calibration_stats$calibration_case_weight,
      calibration_nonevent_weight = calibration_stats$calibration_nonevent_weight,
      calibration_minimum_n_flag = calibration_stats$calibration_minimum_n_flag,
      calibration_supported_bin_count = calibration_stats$calibration_supported_bin_count,
      calibration_supported_bin_min_flag = calibration_stats$calibration_supported_bin_min_flag,
      calibration_support_flag = calibration_stats$calibration_support_flag,
      calibration_fit_success_flag = calibration_stats$calibration_fit_success_flag,
      calibration_regression_suppressed_flag = calibration_stats$calibration_regression_suppressed_flag,
      time_dependent_auc = auc_est,
      brier_score = brier_score,
      scaled_brier_score = scaled_brier,
      case_count_ipcw = ipcw_info$case_count_ipcw,
      nonevent_count_ipcw = ipcw_info$nonevent_count_ipcw,
      known_count = ipcw_info$known_count
    )
    
    performance_rows[[paste0(plan_row$model_id[[1]], "_horizon_", horizon_year)]] <- bind_rows(
      performance_rows[[paste0(plan_row$model_id[[1]], "_horizon_", horizon_year)]],
      make_metric_rows(
        meta_row = cindex_meta,
        horizon_info = horizon_info,
        threshold_value = NA_real_,
        metric_domain = "horizon_summary",
        metrics_named = horizon_metric_values
      )
    )
    
    calibration_diagnostic_rows[[paste0(plan_row$model_id[[1]], "_caldiag_", horizon_year)]] <- tibble(
      stage = "Stage 5",
      dataset = dataset_name,
      model_class = "non_cure",
      model_family = model_family,
      model_id = plan_row$model_id[[1]],
      strategy_family = model_family,
      formula_id = plan_row$formula_id[[1]],
      formula_name = plan_row$formula_name[[1]],
      formula_label = plan_row$formula_label[[1]],
      site_branch = plan_row$site_branch[[1]],
      interaction_branch = plan_row$interaction_branch[[1]],
      anchor_horizon_year = anchor_horizon_year,
      horizon_id = horizon_info$horizon_id,
      horizon_year = horizon_year,
      horizon_label = horizon_info$horizon_label,
      interpretation_tier = horizon_info$interpretation_tier,
      primary_supported_flag = horizon_info$primary_supported_flag,
      interpretation_note = horizon_info$interpretation_note,
      mean_predicted_risk = calibration_stats$mean_predicted_risk,
      mean_calibration_difference = calibration_stats$mean_calibration_difference,
      calibration_intercept_offset = calibration_stats$calibration_intercept_offset,
      calibration_intercept_free = calibration_stats$calibration_intercept_free,
      calibration_slope = calibration_stats$calibration_slope,
      calibration_known_n = calibration_stats$calibration_known_n,
      calibration_case_weight = calibration_stats$calibration_case_weight,
      calibration_nonevent_weight = calibration_stats$calibration_nonevent_weight,
      binary_outcome_support_flag = calibration_stats$binary_outcome_support_flag,
      calibration_minimum_n_flag = calibration_stats$calibration_minimum_n_flag,
      calibration_supported_bin_count = calibration_stats$calibration_supported_bin_count,
      calibration_supported_bin_min_flag = calibration_stats$calibration_supported_bin_min_flag,
      calibration_support_flag = calibration_stats$calibration_support_flag,
      calibration_fit_success_flag = calibration_stats$calibration_fit_success_flag,
      calibration_regression_suppressed_flag = calibration_stats$calibration_regression_suppressed_flag,
      calibration_failure_reason = derive_calibration_failure_reason(calibration_stats),
      calibration_warning = calibration_stats$calibration_warning
    )
    
    calibration_bins <- compute_calibration_bins(
      pred_risk = pred_risk_h,
      ipcw_info = ipcw_info,
      group_count = calibration_group_count
    ) %>%
      mutate(
        stage = "Stage 5",
        dataset = dataset_name,
        model_class = "non_cure",
        model_family = model_family,
        model_id = plan_row$model_id[[1]],
        formula_id = plan_row$formula_id[[1]],
        formula_name = plan_row$formula_name[[1]],
        formula_label = plan_row$formula_label[[1]],
        site_branch = plan_row$site_branch[[1]],
        interaction_branch = plan_row$interaction_branch[[1]],
        anchor_horizon_year = anchor_horizon_year,
        horizon_year = horizon_year,
        horizon_id = horizon_info$horizon_id,
        interpretation_tier = horizon_info$interpretation_tier,
        primary_supported_flag = horizon_info$primary_supported_flag,
        interpretation_note = horizon_info$interpretation_note
      ) %>%
      relocate(
        stage, dataset, model_class, model_family, model_id,
        formula_id, formula_name, formula_label,
        site_branch, interaction_branch, anchor_horizon_year,
        horizon_id, horizon_year, interpretation_tier,
        primary_supported_flag, interpretation_note
      )
    
    calibration_bin_rows[[paste0(plan_row$model_id[[1]], "_calbin_", horizon_year)]] <- calibration_bins
    
    for (threshold_value in risk_thresholds) {
      threshold_metrics <- compute_threshold_vector(
        pred_positive = pred_risk_h >= threshold_value,
        ipcw_info = ipcw_info,
        threshold = threshold_value,
        n_total = nrow(dataset_data)
      )
      
      performance_rows[[paste0(plan_row$model_id[[1]], "_threshold_", horizon_year, "_", threshold_value)]] <- bind_rows(
        performance_rows[[paste0(plan_row$model_id[[1]], "_threshold_", horizon_year, "_", threshold_value)]],
        make_metric_rows(
          meta_row = cindex_meta,
          horizon_info = horizon_info,
          threshold_value = threshold_value,
          metric_domain = "threshold_summary",
          metrics_named = threshold_metrics
        )
      )
    }
  }
  
  ibs_valid <- !is.na(brier_by_horizon)
  ibs_value <- if (sum(ibs_valid) >= 2) {
    valid_horizons <- common_horizons_year[ibs_valid]
    valid_brier <- brier_by_horizon[ibs_valid]
    sum(diff(valid_horizons) * (head(valid_brier, -1) + tail(valid_brier, -1)) / 2) /
      (max(valid_horizons) - min(valid_horizons))
  } else {
    NA_real_
  }
  
  ibs_horizon_info <- horizon_label_lookup(dataset_name, anchor_horizon_year, horizon_registry)
  
  performance_rows[[paste0(plan_row$model_id[[1]], "_ibs")]] <- bind_rows(
    performance_rows[[paste0(plan_row$model_id[[1]], "_ibs")]],
    make_metric_rows(
      meta_row = cindex_meta,
      horizon_info = ibs_horizon_info,
      threshold_value = NA_real_,
      metric_domain = "integrated_summary",
      metrics_named = c(integrated_brier_score = ibs_value)
    )
  )
}

# 🔴 Assemble: downstream tables ===============================
## 🟠 Combine: prediction matrix ===============================
subject_horizon_risk_long <- bind_rows_from_df_lists(
  risk_long_rows,
  object_name = "risk_long_rows"
) %>%
  arrange(
    match(dataset, c("PNU", "SNU", "merged")),
    formula_name,
    model_family,
    unique_person_id,
    horizon_year
  )

## 🟠 Combine: performance stack ===============================
model_registry <- bind_rows_from_df_lists(
  model_registry_rows,
  object_name = "model_registry_rows"
) %>%
  arrange(
    match(dataset, c("PNU", "SNU", "merged")),
    formula_name,
    model_family
  )

all_performance_rows <- c(
  dataset_reference_rows,
  dataset_threshold_reference_rows,
  performance_rows
)

all_performance_rows <- purrr::compact(all_performance_rows)

if (!all(vapply(all_performance_rows, inherits, logical(1), what = "data.frame"))) {
  stop("`all_performance_rows` contains at least one non-data.frame object.", call. = FALSE)
}

model_performance_long <- bind_rows(all_performance_rows) %>%
  arrange(
    match(dataset, c("PNU", "SNU", "merged")),
    formula_name,
    model_family,
    horizon_year,
    threshold,
    metric_domain,
    metric_name
  )

calibration_bins_long <- bind_rows_from_df_lists(
  calibration_bin_rows,
  object_name = "calibration_bin_rows"
) %>%
  arrange(
    match(dataset, c("PNU", "SNU", "merged")),
    formula_name,
    model_family,
    horizon_year,
    bin_index
  )

calibration_diagnostics_long <- bind_rows_from_df_lists(
  calibration_diagnostic_rows,
  object_name = "calibration_diagnostic_rows"
) %>%
  arrange(
    match(dataset, c("PNU", "SNU", "merged")),
    formula_name,
    model_family,
    horizon_year
  )

## 🟠 Draw: pdf figures ===============================
pdf_file <- file.path(export_path, "stage5_summary_plots.pdf")
grDevices::pdf(pdf_file, width = 12, height = 8, onefile = TRUE)

risk_plot_source <- subject_horizon_risk_long %>%
  group_by(dataset, formula_label, model_family, horizon_year, interpretation_tier) %>%
  summarise(
    median_risk = median(risk_pred, na.rm = TRUE),
    risk_p25 = stats::quantile(risk_pred, probs = 0.25, na.rm = TRUE, names = FALSE),
    risk_p75 = stats::quantile(risk_pred, probs = 0.75, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  )

if (nrow(risk_plot_source) > 0) {
  print(
    ggplot(risk_plot_source, aes(x = horizon_year, y = median_risk, color = model_family, fill = model_family)) +
      geom_line(linewidth = 0.8) +
      geom_ribbon(aes(ymin = risk_p25, ymax = risk_p75), alpha = 0.15, linewidth = 0) +
      facet_grid(dataset ~ formula_label, scales = "free_y") +
      scale_x_continuous(breaks = common_horizons_year) +
      labs(
        title = "Stage 5 non-cure predicted risk summary",
        x = "Prediction horizon (years)",
        y = "Median predicted risk"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
  )
}

auc_plot_source <- model_performance_long %>%
  filter(
    metric_name == "time_dependent_auc",
    model_class == "non_cure",
    !is.na(metric_value)
  ) %>%
  arrange(dataset, formula_label, model_family, horizon_year) %>%
  group_by(dataset, formula_label, model_family) %>%
  mutate(
    plot_segment_id = cumsum(c(TRUE, diff(horizon_year) != 1L))
  ) %>%
  ungroup()

auc_line_source <- auc_plot_source %>%
  add_count(dataset, formula_label, model_family, plot_segment_id, name = "segment_n") %>%
  filter(segment_n >= 2L)

if (nrow(auc_plot_source) > 0) {
  print(
    ggplot(auc_plot_source, aes(x = horizon_year, y = metric_value, color = model_family)) +
      geom_line(
        data = auc_line_source,
        aes(group = interaction(model_family, plot_segment_id)),
        linewidth = 0.8
      ) +
      geom_point(size = 1.5) +
      facet_grid(dataset ~ formula_label) +
      scale_x_continuous(breaks = common_horizons_year) +
      labs(
        title = "Stage 5 cumulative/dynamic time-dependent AUC",
        x = "Prediction horizon (years)",
        y = "AUC"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
  )
}

brier_plot_source <- model_performance_long %>%
  filter(
    metric_name == "brier_score",
    model_class == "non_cure",
    !is.na(metric_value)
  ) %>%
  arrange(dataset, formula_label, model_family, horizon_year) %>%
  group_by(dataset, formula_label, model_family) %>%
  mutate(
    plot_segment_id = cumsum(c(TRUE, diff(horizon_year) != 1L))
  ) %>%
  ungroup()

brier_line_source <- brier_plot_source %>%
  add_count(dataset, formula_label, model_family, plot_segment_id, name = "segment_n") %>%
  filter(segment_n >= 2L)

if (nrow(brier_plot_source) > 0) {
  print(
    ggplot(brier_plot_source, aes(x = horizon_year, y = metric_value, color = model_family)) +
      geom_line(
        data = brier_line_source,
        aes(group = interaction(model_family, plot_segment_id)),
        linewidth = 0.8
      ) +
      geom_point(size = 1.5) +
      facet_grid(dataset ~ formula_label) +
      scale_x_continuous(breaks = common_horizons_year) +
      labs(
        title = "Stage 5 horizon-specific Brier score",
        x = "Prediction horizon (years)",
        y = "Brier score"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
  )
}

calibration_plot_source <- calibration_bins_long %>%
  filter(
    horizon_year %in% prediction_horizons_for_plots,
    binary_outcome_support_flag == 1,
    bin_binary_support_flag == 1,
    !is.na(mean_predicted_risk),
    !is.na(observed_risk_ipcw)
  ) %>%
  add_count(
    dataset,
    formula_label,
    model_family,
    horizon_year,
    name = "supported_bin_count_for_plot"
  ) %>%
  filter(
    supported_bin_count_for_plot >= calibration_min_supported_bins_for_regression
  ) %>%
  arrange(dataset, formula_label, model_family, horizon_year, bin_index) %>%
  group_by(dataset, formula_label, model_family, horizon_year) %>%
  mutate(
    calibration_line_segment_id = cumsum(c(TRUE, diff(bin_index) != 1L))
  ) %>%
  ungroup()

calibration_line_source <- calibration_plot_source %>%
  add_count(
    dataset, formula_label, model_family, horizon_year, calibration_line_segment_id,
    name = "segment_n"
  ) %>%
  filter(segment_n >= 2L)

if (nrow(calibration_plot_source) > 0) {
  print(
    ggplot(calibration_plot_source, aes(x = mean_predicted_risk, y = observed_risk_ipcw, color = model_family)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      geom_point(size = 1.5) +
      geom_line(
        data = calibration_line_source,
        aes(group = interaction(model_family, calibration_line_segment_id)),
        linewidth = 0.8
      ) +
      facet_grid(dataset + formula_label ~ horizon_year, scales = "free") +
      labs(
        title = "Stage 5 grouped calibration by prediction horizon (supported bins only)",
        subtitle = "Grouped calibration is shown only when at least two bin-level case/nonevent-supported groups exist.",
        x = "Mean predicted risk",
        y = "Observed risk (IPCW)"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
  )
}

net_benefit_plot_source <- model_performance_long %>%
  filter(
    metric_name == "net_benefit",
    horizon_year %in% prediction_horizons_for_plots,
    model_family %in% c("exponential", "weibull", "lognormal", "loglogistic", "coxph", "treat_all", "treat_none")
  ) %>%
  mutate(
    formula_label = if_else(is.na(formula_label), "Reference", formula_label)
  )

if (nrow(net_benefit_plot_source) > 0) {
  print(
    ggplot(net_benefit_plot_source, aes(x = threshold, y = metric_value, color = model_family)) +
      geom_line(linewidth = 0.8) +
      facet_grid(dataset + formula_label ~ horizon_year, scales = "free_y") +
      labs(
        title = "Stage 5 decision-curve net benefit",
        x = "Risk threshold",
        y = "Net benefit"
      ) +
      theme_bw() +
      theme(legend.position = "bottom")
  )
}

grDevices::dev.off()

# 🔴 Save: stage-five outputs ===============================
## 🟠 Write: csv sources ===============================
model_registry_file <- file.path(export_path, "stage5_model_registry.csv")
risk_long_file <- file.path(export_path, "stage5_subject_horizon_risk_long.csv")
performance_long_file <- file.path(export_path, "stage5_model_performance_long.csv")
calibration_bins_file <- file.path(export_path, "stage5_calibration_bins_long.csv")
calibration_diagnostics_file <- file.path(export_path, "stage5_calibration_diagnostics_long.csv")

readr::write_csv(model_registry, model_registry_file)
readr::write_csv(subject_horizon_risk_long, risk_long_file)
readr::write_csv(model_performance_long, performance_long_file)
readr::write_csv(calibration_bins_long, calibration_bins_file)
readr::write_csv(calibration_diagnostics_long, calibration_diagnostics_file)

## 🟠 Persist: rds bundles ===============================
fitted_models_file <- file.path(export_path, "stage5_fitted_models.rds")
results_bundle_file <- file.path(export_path, "stage5_results_bundle.rds")

saveRDS(fitted_models, fitted_models_file)

stage5_results_bundle <- list(
  stage = "Stage 5",
  created_at = as.character(Sys.time()),
  data_path = normalizePath(data_path, winslash = "/", mustWork = FALSE),
  export_path = normalizePath(export_path, winslash = "/", mustWork = FALSE),
  stage1_bundle_file = normalizePath(stage1_bundle_file, winslash = "/", mustWork = FALSE),
  stage1_datasets_file = normalizePath(stage1_datasets_file, winslash = "/", mustWork = FALSE),
  common_horizons_year = common_horizons_year,
  risk_thresholds = risk_thresholds,
  calibration_min_supported_bins_for_regression = calibration_min_supported_bins_for_regression,
  package_versions = vapply(required_packages, function(pkg) as.character(utils::packageVersion(pkg)), character(1)),
  model_registry = model_registry,
  subject_horizon_risk_long = subject_horizon_risk_long,
  model_performance_long = model_performance_long,
  calibration_bins_long = calibration_bins_long,
  calibration_diagnostics_long = calibration_diagnostics_long
)

saveRDS(stage5_results_bundle, results_bundle_file)

## 🟠 Record: file manifest ===============================
export_manifest <- tibble(
  file_name = c(
    "stage5_model_registry.csv",
    "stage5_subject_horizon_risk_long.csv",
    "stage5_model_performance_long.csv",
    "stage5_calibration_bins_long.csv",
    "stage5_calibration_diagnostics_long.csv",
    "stage5_summary_plots.pdf",
    "stage5_fitted_models.rds",
    "stage5_results_bundle.rds",
    "stage5_export_manifest.csv"
  ),
  object_name = c(
    "model_registry",
    "subject_horizon_risk_long",
    "model_performance_long",
    "calibration_bins_long",
    "calibration_diagnostics_long",
    "summary_plots_pdf",
    "fitted_models",
    "stage5_results_bundle",
    "export_manifest"
  ),
  description = c(
    "Model-level fit registry for all Stage 5 non-cure models",
    "Subject-by-horizon predicted survival and risk long table",
    "Long-format performance table containing horizon summaries, threshold summaries, and benchmark references",
    "Grouped calibration plot source-of-truth table",
    "Model-by-horizon calibration diagnostics with sparse-bin gating and failure reasons",
    "Stage 5 summary figure PDF generated from exported data frames",
    "Fitted non-cure model objects as an RDS list",
    "Reusable Stage 5 results bundle with exported tables",
    "Manifest of all Stage 5 exported files"
  ),
  file_path = c(
    model_registry_file,
    risk_long_file,
    performance_long_file,
    calibration_bins_file,
    calibration_diagnostics_file,
    pdf_file,
    fitted_models_file,
    results_bundle_file,
    file.path(export_path, "stage5_export_manifest.csv")
  )
)

readr::write_csv(export_manifest, file.path(export_path, "stage5_export_manifest.csv"))
