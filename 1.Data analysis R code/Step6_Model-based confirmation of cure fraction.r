# 🔴 Configure: Step6 paths and switches ===============================
DATA_PATH <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
STEP4_EXPORT_DIR <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦4.Step4_Non-cure benchmark models/attachments"
STEP4_FILE_PREFIX <- "step4_nocure"
STEP5_EXPORT_DIR <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦5.Step5_MLE MCM/attachments"



STEP6_EXPORT_DIR <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦6.Step6_formal cure fraction test/attachments'

STEP6_FILE_PREFIX <- "step6_evidence"
PREDICTION_YEARS <- 1:10
SITE_REFERENCE <- "PNU"
RECOMPUTE_STEP5_PERFORMANCE <- TRUE
STEP6_SEED <- 20260316L

STEP5_MASTER_RDS_NAME <- "step5_mle_cure_master.rds"
STEP5_MODEL_GRID_CSV_NAME <- "step5_mle_cure_model_grid.csv"
STEP5_PREPROC_CSV_NAME <- "step5_mle_cure_preprocessing.csv"
STEP5_ANALYSIS_SPEC_CSV_NAME <- "step5_mle_cure_analysis_spec.csv"
STEP5_FIT_REGISTRY_CSV_NAME <- "step5_mle_cure_fit_registry.csv"
STEP5_COEFFICIENTS_CSV_NAME <- "step5_mle_cure_coefficients.csv"
STEP5_CURE_ONLY_CSV_NAME <- "step5_mle_cure_cure_only.csv"
STEP5_MEAN_PRED_CSV_NAME <- "step5_mle_cure_mean_predictions.csv"
STEP5_LATENCY_SUMMARY_CSV_NAME <- "step5_mle_cure_latency_summary.csv"
STEP5_SUBJECT_PRED_CSV_NAME <- "step5_mle_cure_subject_predictions.csv"
STEP5_MANIFEST_CSV_NAME <- "step5_mle_cure_manifest.csv"

# 🔴 Initialize: runtime options and file targets ===============================
options(stringsAsFactors = FALSE, scipen = 999)
set.seed(STEP6_SEED)

DATA_PATH <- path.expand(DATA_PATH)
STEP4_EXPORT_DIR <- path.expand(STEP4_EXPORT_DIR)
STEP5_EXPORT_DIR <- path.expand(STEP5_EXPORT_DIR)
STEP6_EXPORT_DIR <- path.expand(STEP6_EXPORT_DIR)

if (!dir.exists(STEP6_EXPORT_DIR)) {
  dir.create(STEP6_EXPORT_DIR, recursive = TRUE, showWarnings = FALSE)
}

STEP4_ANALYSIS_DATA_PATH <- file.path(STEP4_EXPORT_DIR, paste0(STEP4_FILE_PREFIX, "_analysis_data.csv"))
STEP4_MODEL_REGISTRY_PATH <- file.path(STEP4_EXPORT_DIR, paste0(STEP4_FILE_PREFIX, "_model_registry.csv"))
STEP4_SUBJECT_PRED_PATH <- file.path(STEP4_EXPORT_DIR, paste0(STEP4_FILE_PREFIX, "_subject_yearly_predictions.csv"))
STEP4_MEAN_PRED_PATH <- file.path(STEP4_EXPORT_DIR, paste0(STEP4_FILE_PREFIX, "_mean_yearly_predictions.csv"))
STEP4_PERF_YEARLY_PATH <- file.path(STEP4_EXPORT_DIR, paste0(STEP4_FILE_PREFIX, "_performance_yearly.csv"))
STEP4_PERF_SUMMARY_PATH <- file.path(STEP4_EXPORT_DIR, paste0(STEP4_FILE_PREFIX, "_performance_summary.csv"))
STEP4_FIT_BUNDLE_PATH <- file.path(STEP4_EXPORT_DIR, paste0(STEP4_FILE_PREFIX, "_fit_bundle.rds"))
STEP4_MANIFEST_PATH <- file.path(STEP4_EXPORT_DIR, paste0(STEP4_FILE_PREFIX, "_manifest.csv"))

STEP5_MASTER_RDS_PATH <- file.path(STEP5_EXPORT_DIR, STEP5_MASTER_RDS_NAME)
STEP5_MODEL_GRID_PATH <- file.path(STEP5_EXPORT_DIR, STEP5_MODEL_GRID_CSV_NAME)
STEP5_PREPROC_PATH <- file.path(STEP5_EXPORT_DIR, STEP5_PREPROC_CSV_NAME)
STEP5_ANALYSIS_SPEC_PATH <- file.path(STEP5_EXPORT_DIR, STEP5_ANALYSIS_SPEC_CSV_NAME)
STEP5_FIT_REGISTRY_PATH <- file.path(STEP5_EXPORT_DIR, STEP5_FIT_REGISTRY_CSV_NAME)
STEP5_COEFFICIENTS_PATH <- file.path(STEP5_EXPORT_DIR, STEP5_COEFFICIENTS_CSV_NAME)
STEP5_CURE_ONLY_PATH <- file.path(STEP5_EXPORT_DIR, STEP5_CURE_ONLY_CSV_NAME)
STEP5_MEAN_PRED_PATH <- file.path(STEP5_EXPORT_DIR, STEP5_MEAN_PRED_CSV_NAME)
STEP5_LATENCY_SUMMARY_PATH <- file.path(STEP5_EXPORT_DIR, STEP5_LATENCY_SUMMARY_CSV_NAME)
STEP5_SUBJECT_PRED_PATH <- file.path(STEP5_EXPORT_DIR, STEP5_SUBJECT_PRED_CSV_NAME)
STEP5_MANIFEST_PATH <- file.path(STEP5_EXPORT_DIR, STEP5_MANIFEST_CSV_NAME)

ANALYSIS_SPEC_OUT <- file.path(STEP6_EXPORT_DIR, paste0(STEP6_FILE_PREFIX, "_analysis_spec.csv"))
IMPORT_AUDIT_OUT <- file.path(STEP6_EXPORT_DIR, paste0(STEP6_FILE_PREFIX, "_import_audit.csv"))
KM_BENCHMARK_OUT <- file.path(STEP6_EXPORT_DIR, paste0(STEP6_FILE_PREFIX, "_km_benchmark.csv"))
MODEL_YEARLY_OUT <- file.path(STEP6_EXPORT_DIR, paste0(STEP6_FILE_PREFIX, "_model_yearly_metrics.csv"))
MODEL_SUMMARY_OUT <- file.path(STEP6_EXPORT_DIR, paste0(STEP6_FILE_PREFIX, "_model_summary_metrics.csv"))
PAIR_REGISTRY_OUT <- file.path(STEP6_EXPORT_DIR, paste0(STEP6_FILE_PREFIX, "_pair_registry.csv"))
PAIR_YEARLY_OUT <- file.path(STEP6_EXPORT_DIR, paste0(STEP6_FILE_PREFIX, "_pair_yearly.csv"))
PAIR_SUMMARY_OUT <- file.path(STEP6_EXPORT_DIR, paste0(STEP6_FILE_PREFIX, "_pair_summary.csv"))
COEFFICIENTS_OUT <- file.path(STEP6_EXPORT_DIR, paste0(STEP6_FILE_PREFIX, "_coefficients.csv"))
LANE_YEARLY_OUT <- file.path(STEP6_EXPORT_DIR, paste0(STEP6_FILE_PREFIX, "_lane_yearly.csv"))
LANE_SUMMARY_OUT <- file.path(STEP6_EXPORT_DIR, paste0(STEP6_FILE_PREFIX, "_lane_summary.csv"))
MASTER_RDS_OUT <- file.path(STEP6_EXPORT_DIR, paste0(STEP6_FILE_PREFIX, "_master.rds"))
MANIFEST_OUT <- file.path(STEP6_EXPORT_DIR, paste0(STEP6_FILE_PREFIX, "_manifest.csv"))

# 🔴 Validate: required packages and imported artifacts ===============================
## 🟠 Check: required packages ===============================
required_pkgs <- c("survival")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("다음 패키지가 필요합니다: ", paste(missing_pkgs, collapse = ", "))
}
suppressPackageStartupMessages({
  library(survival)
})

## 🟠 Check: Step4 and Step5 exported files ===============================
required_input_paths <- c(
  DATA_PATH,
  STEP4_ANALYSIS_DATA_PATH,
  STEP4_MODEL_REGISTRY_PATH,
  STEP4_SUBJECT_PRED_PATH,
  STEP4_MEAN_PRED_PATH,
  STEP4_PERF_YEARLY_PATH,
  STEP4_PERF_SUMMARY_PATH,
  STEP4_FIT_BUNDLE_PATH,
  STEP4_MANIFEST_PATH,
  STEP5_MASTER_RDS_PATH,
  STEP5_MODEL_GRID_PATH,
  STEP5_PREPROC_PATH,
  STEP5_ANALYSIS_SPEC_PATH,
  STEP5_FIT_REGISTRY_PATH,
  STEP5_COEFFICIENTS_PATH,
  STEP5_CURE_ONLY_PATH,
  STEP5_MEAN_PRED_PATH,
  STEP5_LATENCY_SUMMARY_PATH,
  STEP5_SUBJECT_PRED_PATH,
  STEP5_MANIFEST_PATH
)

missing_input_paths <- required_input_paths[!file.exists(required_input_paths)]
if (length(missing_input_paths) > 0L) {
  stop("다음 입력 파일이 존재하지 않습니다: ", paste(missing_input_paths, collapse = " | "))
}

# 🔴 Define: helper functions for harmonization and metrics ===============================
## 🟠 Helpers: row binding and audit metadata ===============================
bind_rows_safe <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) == 0L) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

bind_rows_fill <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) == 0L) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  all_names <- unique(unlist(lapply(x, names)))
  filled <- lapply(x, function(df) {
    missing_cols <- setdiff(all_names, names(df))
    if (length(missing_cols) > 0L) {
      for (nm in missing_cols) {
        df[[nm]] <- NA
      }
    }
    df <- df[, all_names, drop = FALSE]
    rownames(df) <- NULL
    df
  })
  out <- do.call(rbind, filled)
  rownames(out) <- NULL
  out
}

read_csv_safely <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

object_nrow <- function(obj) {
  if (is.data.frame(obj) || is.matrix(obj)) {
    return(nrow(obj))
  }
  if (is.list(obj)) {
    return(length(obj))
  }
  NA_integer_
}

object_ncol <- function(obj) {
  if (is.data.frame(obj) || is.matrix(obj)) {
    return(ncol(obj))
  }
  NA_integer_
}

collapse_unique <- function(x) {
  x <- unique(stats::na.omit(as.character(x)))
  if (length(x) == 0L) "" else paste(sort(x), collapse = "|")
}

make_audit_row <- function(source_step, object_name, file_path, object_ref, used_directly = TRUE) {
  data.frame(
    source_step = source_step,
    object_name = object_name,
    file_name = basename(file_path),
    file_path = file_path,
    exists = file.exists(file_path),
    object_class = paste(class(object_ref), collapse = "|"),
    n_rows = object_nrow(object_ref),
    n_cols = object_ncol(object_ref),
    file_size_bytes = if (file.exists(file_path)) unname(file.info(file_path)$size) else NA_real_,
    used_directly = used_directly,
    stringsAsFactors = FALSE
  )
}

to_numeric_if_present <- function(df, cols) {
  for (nm in intersect(cols, names(df))) {
    df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
  }
  df
}

to_logical_if_present <- function(df, cols) {
  for (nm in intersect(cols, names(df))) {
    if (is.logical(df[[nm]])) next
    v <- trimws(as.character(df[[nm]]))
    df[[nm]] <- ifelse(
      v %in% c("TRUE", "T", "1"),
      TRUE,
      ifelse(v %in% c("FALSE", "F", "0"), FALSE, NA)
    )
  }
  df
}

safe_min <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}

safe_max <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

gap_to_range <- function(value, lower, upper) {
  if (!is.finite(value) || !is.finite(lower) || !is.finite(upper)) {
    return(NA_real_)
  }
  if (value < lower) {
    return(value - lower)
  }
  if (value > upper) {
    return(value - upper)
  }
  0
}

safe_value_at_year <- function(df, year, colname) {
  hit <- df[df$year == year, colname, drop = TRUE]
  hit <- hit[is.finite(hit)]
  if (length(hit) == 0L) NA_real_ else as.numeric(hit[1L])
}

## 🟠 Helpers: standardize identifiers and metadata ===============================
standardize_dataset <- function(x) {
  y <- toupper(trimws(as.character(x)))
  out <- ifelse(
    y == "MERGED",
    "merged",
    ifelse(y %in% c("PNU", "SNU"), y, y)
  )
  out
}

standardize_lane <- function(x) {
  y <- trimws(as.character(x))
  out <- ifelse(
    y %in% c("AFT_param", "AFT_parametric"),
    "AFT_param",
    ifelse(y %in% c("PH_semiparam", "PH_semiparametric"), "PH_semiparam", y)
  )
  out
}

standardize_family <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- ifelse(
    y %in% c("exp", "exponential"),
    "exponential",
    ifelse(
      y == "weibull",
      "weibull",
      ifelse(
        y %in% c("llogis", "loglogistic"),
        "loglogistic",
        ifelse(
          y %in% c("lnorm", "lognormal"),
          "lognormal",
          ifelse(
            y %in% c("coxph", "cox_latency", "cox"),
            "cox",
            y
          )
        )
      )
    )
  )
  out
}

standardize_step4_registry <- function(df) {
  df$dataset_id <- standardize_dataset(df$dataset_id)
  df$model_id_std <- as.character(df$model_id)
  df$lane_std <- standardize_lane(df$lane)
  df$family_std <- standardize_family(df$family)
  df$model_class <- "nocure"
  df$source_step <- "step4"
  df <- to_numeric_if_present(
    df,
    c("n_subject", "n_event", "n_censored", "n_params", "logLik_full", "AIC", "BIC",
      "partial_logLik", "AIC_partial", "max_followup_years", "max_event_years")
  )
  df <- to_logical_if_present(df, c("fit_success", "convergence_flag", "singular_hessian_flag", "boundary_flag"))
  df
}

standardize_step4_mean <- function(df) {
  df$dataset_id <- standardize_dataset(df$dataset_id)
  df$model_id_std <- as.character(df$model_id)
  df$lane_std <- standardize_lane(df$lane)
  df$family_std <- standardize_family(df$family)
  df$model_class <- "nocure"
  df$source_step <- "step4"
  df$subgroup_id <- "overall"
  df <- to_numeric_if_present(df, c("year", "mean_surv", "mean_risk"))
  df
}

standardize_step4_perf_yearly <- function(df) {
  df$dataset_id <- standardize_dataset(df$dataset_id)
  df$model_id_std <- as.character(df$model_id)
  df$lane_std <- standardize_lane(df$lane)
  df$family_std <- standardize_family(df$family)
  df$model_class <- "nocure"
  df$source_step <- "step4"
  df <- to_numeric_if_present(df, c("year", "auc_cd", "brier", "n_cases", "n_controls", "n_censored_before_t", "censor_surv_at_t"))
  df
}

standardize_step4_perf_summary <- function(df) {
  df$dataset_id <- standardize_dataset(df$dataset_id)
  df$model_id_std <- as.character(df$model_id)
  df$lane_std <- standardize_lane(df$lane)
  df$family_std <- standardize_family(df$family)
  df$model_class <- "nocure"
  df$source_step <- "step4"
  df <- to_numeric_if_present(df, c("ibs_1to10", "ibs_available", "ibs_available_from", "ibs_available_to", "auc_mean_available", "brier_mean_available", "n_evaluable_years"))
  df
}

standardize_step4_subject_predictions <- function(df) {
  df$dataset_id <- standardize_dataset(df$dataset_id)
  df$model_id_std <- as.character(df$model_id)
  df$subject_key <- as.character(df$subject_uid)
  df <- to_numeric_if_present(df, c("year", "predicted_surv", "predicted_risk"))
  df
}

standardize_step5_registry <- function(df) {
  df$dataset_id <- standardize_dataset(df$dataset)
  df$model_id_std <- as.character(df$fit_id)
  df$lane_std <- standardize_lane(df$lane)
  df$family_std <- standardize_family(df$family)
  df$source_step <- "step5"
  df <- to_numeric_if_present(
    df,
    c("n", "n_event", "n_censor", "n_parameters", "logLik", "AIC", "BIC", "mean_cure_fraction", "fit_time_sec")
  )
  df <- to_logical_if_present(df, c("success", "converged_flag", "boundary_flag"))
  df$fit_success <- df$success
  df$convergence_flag <- df$converged_flag
  df$n_subject <- df$n
  df$n_censored <- df$n_censor
  df$logLik_full <- df$logLik
  df$likelihood_type <- ifelse(df$family_std == "cox", "partial", "full")
  df$partial_logLik <- NA_real_
  df$AIC_partial <- NA_real_
  df
}

standardize_step5_mean <- function(df, registry_std) {
  df$dataset_id <- standardize_dataset(df$dataset)
  df$model_id_std <- as.character(df$fit_id)
  df$model_class <- "cure_MLE"
  df$source_step <- "step5"
  df$subgroup_id <- "overall"
  df$mean_surv <- suppressWarnings(as.numeric(df$meanS))
  df$mean_risk <- suppressWarnings(as.numeric(df$meanRisk))
  df$year <- suppressWarnings(as.numeric(df$year))
  keep_cols <- c("model_id_std", "lane_std", "family_std", "cov_spec")
  df <- merge(df, registry_std[, keep_cols, drop = FALSE], by = "model_id_std", all.x = TRUE, sort = FALSE)
  df
}

standardize_step5_subject_predictions <- function(df, registry_std) {
  df$dataset_id <- standardize_dataset(df$dataset)
  df$model_id_std <- as.character(df$fit_id)
  df$model_class <- "cure_MLE"
  df$source_step <- "step5"
  df$subject_key <- as.character(df$site_id)
  df$predicted_surv <- suppressWarnings(as.numeric(df$Shat_pop))
  df$predicted_risk <- suppressWarnings(as.numeric(df$Riskhat_pop))
  df$year <- suppressWarnings(as.numeric(df$year))
  keep_cols <- c("model_id_std", "lane_std", "family_std", "cov_spec")
  df <- merge(df, registry_std[, keep_cols, drop = FALSE], by = "model_id_std", all.x = TRUE, sort = FALSE)
  df
}

standardize_step5_cure_only <- function(df) {
  df$dataset_id <- standardize_dataset(df$dataset)
  df$model_id_std <- as.character(df$fit_id)
  df$family_std <- standardize_family(df$family)
  df$lane_std <- standardize_lane(df$lane)
  df <- to_numeric_if_present(
    df,
    c("cure_fraction_mean", "cure_fraction_median", "cure_fraction_min", "cure_fraction_max", "uncured_probability_mean")
  )
  df <- to_logical_if_present(df, c("boundary_flag"))
  df
}

standardize_step5_latency_summary <- function(df) {
  df$dataset_id <- standardize_dataset(df$dataset)
  df$model_id_std <- as.character(df$fit_id)
  df <- to_numeric_if_present(df, c("year", "meanSu", "medianSu", "minSu", "maxSu"))
  df
}

standardize_step5_coefficients <- function(df, registry_std) {
  df$dataset_id <- standardize_dataset(df$dataset)
  df$model_id_std <- as.character(df$fit_id)
  df <- to_numeric_if_present(df, c("estimate", "std_error", "zvalue", "pvalue", "lower95", "upper95"))
  keep_cols <- c("model_id_std", "lane_std", "family_std", "cov_spec", "model_class")
  df <- merge(df, registry_std[, keep_cols, drop = FALSE], by = "model_id_std", all.x = TRUE, sort = FALSE)
  df$source_step <- "step5"
  df
}

## 🟠 Helpers: KM and time-dependent performance ===============================
fit_censoring_km <- function(time, event) {
  fit <- survival::survfit(survival::Surv(time, 1 - event) ~ 1)
  ctime <- fit$time
  csurv <- fit$surv
  
  eval_fun <- function(tvals, left = FALSE) {
    tt <- as.numeric(tvals)
    if (left) {
      tt <- pmax(tt - 1e-10, 0)
    }
    idx <- findInterval(tt, ctime)
    ifelse(idx <= 0, 1, csurv[idx])
  }
  
  list(time = ctime, surv = csurv, eval = eval_fun)
}

weighted_auc_cd <- function(marker, case_w, control_w) {
  keep_case <- is.finite(marker) & is.finite(case_w) & (case_w > 0)
  keep_ctrl <- is.finite(marker) & is.finite(control_w) & (control_w > 0)
  
  if (sum(keep_case) == 0 || sum(keep_ctrl) == 0) {
    return(NA_real_)
  }
  
  m_case <- marker[keep_case]
  w_case <- case_w[keep_case]
  m_ctrl <- marker[keep_ctrl]
  w_ctrl <- control_w[keep_ctrl]
  
  ctrl_key <- signif(m_ctrl, 12)
  ctrl_agg <- stats::aggregate(w_ctrl, by = list(marker = ctrl_key), FUN = sum)
  ctrl_agg <- ctrl_agg[order(ctrl_agg$marker), , drop = FALSE]
  rownames(ctrl_agg) <- NULL
  
  cum_le <- cumsum(ctrl_agg$x)
  case_key <- signif(m_case, 12)
  
  idx_le <- findInterval(case_key, ctrl_agg$marker)
  idx_eq <- match(case_key, ctrl_agg$marker)
  
  weight_le <- ifelse(idx_le <= 0, 0, cum_le[idx_le])
  weight_eq <- ifelse(is.na(idx_eq), 0, ctrl_agg$x[idx_eq])
  weight_lt <- weight_le - weight_eq
  
  auc_num <- sum(w_case * (weight_lt + 0.5 * weight_eq))
  auc_den <- sum(w_case) * sum(w_ctrl)
  
  if (!is.finite(auc_den) || auc_den <= 0) {
    return(NA_real_)
  }
  
  auc_num / auc_den
}

brier_ipcw <- function(surv_prob, time, event, t_eval, Gfit) {
  g_t <- Gfit$eval(t_eval, left = FALSE)
  g_y_left <- Gfit$eval(time, left = TRUE)
  
  invalid_case <- (event == 1 & time <= t_eval) & (!is.finite(g_y_left) | g_y_left <= 0)
  if (!is.finite(g_t) || g_t <= 0 || any(invalid_case)) {
    return(NA_real_)
  }
  
  term_event <- ifelse(event == 1 & time <= t_eval, (surv_prob^2) / g_y_left, 0)
  term_survive <- ifelse(time > t_eval, ((1 - surv_prob)^2) / g_t, 0)
  
  mean(term_event + term_survive)
}

trapz_mean <- function(x, y) {
  if (length(x) < 2L || length(y) < 2L) {
    return(NA_real_)
  }
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  if (max(x) <= min(x)) {
    return(NA_real_)
  }
  sum(diff(x) * (y[-length(y)] + y[-1]) / 2) / (max(x) - min(x))
}

evaluate_model_performance <- function(time, event, surv_mat, risk_mat, times, Gfit) {
  yearly_rows <- vector("list", length(times))
  
  for (j in seq_along(times)) {
    t_eval <- times[j]
    g_t <- Gfit$eval(t_eval, left = FALSE)
    g_y_left <- Gfit$eval(time, left = TRUE)
    
    case_w <- ifelse(event == 1 & time <= t_eval, 1 / g_y_left, 0)
    control_w <- ifelse(time > t_eval, 1 / g_t, 0)
    
    if (!is.finite(g_t) || g_t <= 0 || any((event == 1 & time <= t_eval) & (!is.finite(g_y_left) | g_y_left <= 0))) {
      auc_t <- NA_real_
      brier_t <- NA_real_
      evaluable <- FALSE
    } else {
      auc_t <- weighted_auc_cd(risk_mat[, j], case_w, control_w)
      brier_t <- brier_ipcw(surv_mat[, j], time, event, t_eval, Gfit)
      evaluable <- is.finite(auc_t) && is.finite(brier_t)
    }
    
    yearly_rows[[j]] <- data.frame(
      year = t_eval,
      auc_cd = auc_t,
      brier = brier_t,
      n_cases = sum(event == 1 & time <= t_eval),
      n_controls = sum(time > t_eval),
      n_censored_before_t = sum(event == 0 & time <= t_eval),
      censor_surv_at_t = ifelse(is.finite(g_t), g_t, NA_real_),
      evaluable = evaluable,
      stringsAsFactors = FALSE
    )
  }
  
  yearly_df <- bind_rows_safe(yearly_rows)
  brier_ok <- is.finite(yearly_df$brier)
  auc_ok <- is.finite(yearly_df$auc_cd)
  
  ibs_1to10 <- if (all(brier_ok)) trapz_mean(yearly_df$year, yearly_df$brier) else NA_real_
  ibs_available <- if (sum(brier_ok) >= 2L) trapz_mean(yearly_df$year[brier_ok], yearly_df$brier[brier_ok]) else NA_real_
  
  summary_df <- data.frame(
    ibs_1to10 = ibs_1to10,
    ibs_available = ibs_available,
    ibs_available_from = ifelse(any(brier_ok), min(yearly_df$year[brier_ok]), NA_real_),
    ibs_available_to = ifelse(any(brier_ok), max(yearly_df$year[brier_ok]), NA_real_),
    auc_mean_available = ifelse(any(auc_ok), mean(yearly_df$auc_cd[auc_ok]), NA_real_),
    brier_mean_available = ifelse(any(brier_ok), mean(yearly_df$brier[brier_ok]), NA_real_),
    n_evaluable_years = sum(brier_ok),
    stringsAsFactors = FALSE
  )
  
  list(yearly = yearly_df, summary = summary_df)
}

km_yearly_summary <- function(branch_df, years) {
  fit <- survival::survfit(
    survival::Surv(time_years_model, event_transition) ~ 1,
    data = branch_df
  )
  sm <- summary(fit, times = years, extend = TRUE)
  
  surv_vals <- rep(NA_real_, length(years))
  lower_vals <- rep(NA_real_, length(years))
  upper_vals <- rep(NA_real_, length(years))
  nrisk_vals <- rep(NA_real_, length(years))
  
  if (length(sm$surv) > 0L) {
    surv_vals[seq_along(sm$surv)] <- sm$surv
  }
  if (!is.null(sm$lower) && length(sm$lower) > 0L) {
    lower_vals[seq_along(sm$lower)] <- sm$lower
  }
  if (!is.null(sm$upper) && length(sm$upper) > 0L) {
    upper_vals[seq_along(sm$upper)] <- sm$upper
  }
  if (!is.null(sm$n.risk) && length(sm$n.risk) > 0L) {
    nrisk_vals[seq_along(sm$n.risk)] <- sm$n.risk
  }
  
  max_event_year <- if (sum(branch_df$event_transition == 1L) > 0L) max(branch_df$time_years_model[branch_df$event_transition == 1L]) else NA_real_
  max_followup_year <- max(branch_df$time_years_model, na.rm = TRUE)
  tail_reference_year <- if (any(years <= max_event_year, na.rm = TRUE)) max(years[years <= max_event_year]) else NA_real_
  
  data.frame(
    dataset_id = unique(branch_df$dataset_id)[1],
    year = years,
    km_surv = surv_vals,
    km_risk = 1 - surv_vals,
    km_lcl = lower_vals,
    km_ucl = upper_vals,
    n_risk = nrisk_vals,
    n_subject = nrow(branch_df),
    n_event = sum(branch_df$event_transition == 1L),
    n_censored = sum(branch_df$event_transition == 0L),
    last_event_year = max_event_year,
    max_followup_year = max_followup_year,
    tail_reference_year = tail_reference_year,
    support_flag = ifelse(is.finite(max_event_year) & years <= max_event_year, "within_last_event", "beyond_last_event"),
    beyond_max_followup_flag = as.integer(years > max_followup_year),
    stringsAsFactors = FALSE
  )
}

attach_km_metrics <- function(model_yearly_df, km_df) {
  out <- merge(model_yearly_df, km_df, by = c("dataset_id", "year"), all.x = TRUE, sort = FALSE)
  out$shortfall_model_minus_km <- out$mean_risk - out$km_risk
  out$ratio_km_to_model <- ifelse(is.finite(out$mean_risk) & out$mean_risk > 0, out$km_risk / out$mean_risk, NA_real_)
  out$missedpos100 <- 100 * out$shortfall_model_minus_km
  out
}

## 🟠 Helpers: coefficient extraction and comparison builders ===============================
long_to_matrix <- function(pred_long, analysis_df, years, value_col) {
  keys <- as.character(analysis_df$subject_key)
  out <- matrix(
    NA_real_,
    nrow = length(keys),
    ncol = length(years),
    dimnames = list(keys, paste0("year_", years))
  )
  
  for (yr in years) {
    sub <- pred_long[pred_long$year == yr, c("subject_key", value_col), drop = FALSE]
    sub <- sub[!duplicated(sub$subject_key), , drop = FALSE]
    vals <- as.numeric(sub[[value_col]])
    names(vals) <- as.character(sub$subject_key)
    out[, paste0("year_", yr)] <- vals[keys]
  }
  
  out
}

compute_step5_performance <- function(subject_pred_std, fit_registry_std, analysis_branches, years) {
  yearly_rows <- list()
  summary_rows <- list()
  
  for (model_id in unique(subject_pred_std$model_id_std)) {
    pred_sub <- subject_pred_std[subject_pred_std$model_id_std == model_id, , drop = FALSE]
    meta <- fit_registry_std[fit_registry_std$model_id_std == model_id, , drop = FALSE]
    
    if (nrow(pred_sub) == 0L || nrow(meta) == 0L) next
    
    dataset_id <- meta$dataset_id[1]
    analysis_df <- analysis_branches[[dataset_id]]
    if (is.null(analysis_df) || nrow(analysis_df) == 0L) next
    
    surv_mat <- long_to_matrix(pred_sub, analysis_df, years, value_col = "predicted_surv")
    risk_mat <- 1 - surv_mat
    Gfit <- fit_censoring_km(analysis_df$time_years_model, analysis_df$event_transition)
    
    perf <- evaluate_model_performance(
      time = analysis_df$time_years_model,
      event = analysis_df$event_transition,
      surv_mat = surv_mat,
      risk_mat = risk_mat,
      times = years,
      Gfit = Gfit
    )
    
    meta_yearly <- data.frame(
      dataset_id = dataset_id,
      subgroup_id = "overall",
      model_id_std = model_id,
      model_class = meta$model_class[1],
      lane_std = meta$lane_std[1],
      family_std = meta$family_std[1],
      cov_spec = meta$cov_spec[1],
      source_step = "step5",
      stringsAsFactors = FALSE
    )
    
    yearly_rows[[length(yearly_rows) + 1L]] <- cbind(
      meta_yearly[rep(1, nrow(perf$yearly)), , drop = FALSE],
      perf$yearly
    )
    
    summary_rows[[length(summary_rows) + 1L]] <- cbind(
      meta_yearly,
      perf$summary
    )
  }
  
  list(
    yearly = bind_rows_fill(yearly_rows),
    summary = bind_rows_fill(summary_rows)
  )
}

summarize_model_yearly <- function(model_yearly_df) {
  rows <- list()
  ids <- unique(model_yearly_df$model_id_std)
  
  for (model_id in ids) {
    sub <- model_yearly_df[model_yearly_df$model_id_std == model_id, , drop = FALSE]
    if (nrow(sub) == 0L) next
    
    tail_year <- unique(stats::na.omit(sub$tail_reference_year))
    tail_year <- if (length(tail_year) == 0L) NA_real_ else as.numeric(tail_year[1])
    
    tail_gap <- NA_real_
    tail_abs_gap <- NA_real_
    if (is.finite(tail_year)) {
      tail_hit <- sub[sub$year == tail_year, , drop = FALSE]
      if (nrow(tail_hit) > 0L) {
        tail_gap <- as.numeric(tail_hit$shortfall_model_minus_km[1])
        tail_abs_gap <- abs(tail_gap)
      }
    }
    
    finite_shortfall <- which(is.finite(sub$shortfall_model_minus_km))
    max_abs_shortfall <- if (length(finite_shortfall) > 0L) max(abs(sub$shortfall_model_minus_km[finite_shortfall])) else NA_real_
    year_max_abs_shortfall <- if (length(finite_shortfall) > 0L) sub$year[finite_shortfall][which.max(abs(sub$shortfall_model_minus_km[finite_shortfall]))] else NA_real_
    
    rows[[length(rows) + 1L]] <- data.frame(
      dataset_id = sub$dataset_id[1],
      model_id_std = model_id,
      risk_1y = safe_value_at_year(sub, 1, "mean_risk"),
      risk_5y = safe_value_at_year(sub, 5, "mean_risk"),
      risk_10y = safe_value_at_year(sub, 10, "mean_risk"),
      shortfall_1y = safe_value_at_year(sub, 1, "shortfall_model_minus_km"),
      shortfall_5y = safe_value_at_year(sub, 5, "shortfall_model_minus_km"),
      shortfall_10y = safe_value_at_year(sub, 10, "shortfall_model_minus_km"),
      auc_1y = safe_value_at_year(sub, 1, "auc_cd"),
      auc_5y = safe_value_at_year(sub, 5, "auc_cd"),
      auc_10y = safe_value_at_year(sub, 10, "auc_cd"),
      brier_1y = safe_value_at_year(sub, 1, "brier"),
      brier_5y = safe_value_at_year(sub, 5, "brier"),
      brier_10y = safe_value_at_year(sub, 10, "brier"),
      tail_reference_year = tail_year,
      tail_gap = tail_gap,
      tail_abs_gap = tail_abs_gap,
      max_abs_shortfall = max_abs_shortfall,
      year_max_abs_shortfall = year_max_abs_shortfall,
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows_fill(rows)
}

summarize_latency_yearly <- function(latency_df) {
  rows <- list()
  ids <- unique(latency_df$model_id_std)
  
  for (model_id in ids) {
    sub <- latency_df[latency_df$model_id_std == model_id, , drop = FALSE]
    if (nrow(sub) == 0L) next
    
    rows[[length(rows) + 1L]] <- data.frame(
      dataset_id = sub$dataset_id[1],
      model_id_std = model_id,
      latency_meanSu_1y = safe_value_at_year(sub, 1, "meanSu"),
      latency_meanSu_5y = safe_value_at_year(sub, 5, "meanSu"),
      latency_meanSu_10y = safe_value_at_year(sub, 10, "meanSu"),
      latency_medianSu_10y = safe_value_at_year(sub, 10, "medianSu"),
      latency_minSu_10y = safe_value_at_year(sub, 10, "minSu"),
      latency_maxSu_10y = safe_value_at_year(sub, 10, "maxSu"),
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows_fill(rows)
}

validate_mean_predictions <- function(subject_long_std, mean_df_std) {
  if (nrow(subject_long_std) == 0L || nrow(mean_df_std) == 0L) {
    return(NA_real_)
  }
  
  agg_surv <- stats::aggregate(predicted_surv ~ dataset_id + model_id_std + year, data = subject_long_std, FUN = mean)
  agg_risk <- stats::aggregate(predicted_risk ~ dataset_id + model_id_std + year, data = subject_long_std, FUN = mean)
  names(agg_surv)[4] <- "agg_mean_surv"
  names(agg_risk)[4] <- "agg_mean_risk"
  
  chk <- merge(agg_surv, agg_risk, by = c("dataset_id", "model_id_std", "year"), all = TRUE)
  chk <- merge(
    chk,
    mean_df_std[, c("dataset_id", "model_id_std", "year", "mean_surv", "mean_risk"), drop = FALSE],
    by = c("dataset_id", "model_id_std", "year"),
    all = TRUE
  )
  
  diffs <- c(abs(chk$agg_mean_surv - chk$mean_surv), abs(chk$agg_mean_risk - chk$mean_risk))
  diffs <- diffs[is.finite(diffs)]
  if (length(diffs) == 0L) NA_real_ else max(diffs)
}

extract_step4_coefficients <- function(fit_bundle, registry_std) {
  fit_list <- fit_bundle$fit_objects
  if (is.null(fit_list) || length(fit_list) == 0L) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  
  rows <- list()
  
  for (model_id in names(fit_list)) {
    entry <- fit_list[[model_id]]
    fit_obj <- entry$fit
    if (is.null(fit_obj)) next
    
    coef_vec <- tryCatch(stats::coef(fit_obj), error = function(e) NULL)
    if (is.null(coef_vec) || length(coef_vec) == 0L) next
    
    vc <- tryCatch(stats::vcov(fit_obj), error = function(e) NULL)
    se <- rep(NA_real_, length(coef_vec))
    names(se) <- names(coef_vec)
    
    if (!is.null(vc)) {
      se_all <- suppressWarnings(sqrt(diag(vc)))
      se_all <- as.numeric(se_all)
      names(se_all) <- rownames(vc)
      if (all(is.na(names(se_all)))) {
        names(se_all) <- names(coef_vec)[seq_along(se_all)]
      }
      se[names(coef_vec)] <- se_all[names(coef_vec)]
    }
    
    zvalue <- coef_vec / se
    pvalue <- 2 * stats::pnorm(-abs(zvalue))
    lower95 <- coef_vec - 1.96 * se
    upper95 <- coef_vec + 1.96 * se
    
    meta <- registry_std[registry_std$model_id_std == model_id, , drop = FALSE]
    
    rows[[length(rows) + 1L]] <- data.frame(
      dataset_id = if (nrow(meta) > 0L) meta$dataset_id[1] else NA_character_,
      model_id_std = model_id,
      model_class = "nocure",
      lane_std = if (nrow(meta) > 0L) meta$lane_std[1] else NA_character_,
      family_std = if (nrow(meta) > 0L) meta$family_std[1] else NA_character_,
      cov_spec = if (nrow(meta) > 0L) meta$cov_spec[1] else NA_character_,
      component = "nocure",
      term = names(coef_vec),
      estimate = as.numeric(coef_vec),
      std_error = as.numeric(se),
      zvalue = as.numeric(zvalue),
      pvalue = as.numeric(pvalue),
      lower95 = as.numeric(lower95),
      upper95 = as.numeric(upper95),
      source_scale = if (inherits(fit_obj, "coxph")) "cox_ph_log_hazard_ratio" else "survreg_location",
      source_step = "step4",
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows_fill(rows)
}

lookup_model_id <- function(model_summary_df, dataset_id, model_class, family_std, cov_spec, lane_std = NULL) {
  sub <- model_summary_df[
    model_summary_df$dataset_id == dataset_id &
      model_summary_df$model_class == model_class &
      model_summary_df$family_std == family_std &
      model_summary_df$cov_spec == cov_spec,
    ,
    drop = FALSE
  ]
  
  if (!is.null(lane_std)) {
    sub <- sub[sub$lane_std == lane_std, , drop = FALSE]
  }
  
  if ("fit_success" %in% names(sub)) {
    sub <- sub[is.na(sub$fit_success) | sub$fit_success, , drop = FALSE]
  }
  
  if (nrow(sub) == 0L) {
    return(NA_character_)
  }
  
  as.character(sub$model_id_std[1])
}

build_pairwise_yearly <- function(pair_registry_df, model_yearly_df) {
  rows <- list()
  
  for (i in seq_len(nrow(pair_registry_df))) {
    pair_i <- pair_registry_df[i, , drop = FALSE]
    
    left_df <- model_yearly_df[model_yearly_df$model_id_std == pair_i$left_model_id, , drop = FALSE]
    right_df <- model_yearly_df[model_yearly_df$model_id_std == pair_i$right_model_id, , drop = FALSE]
    
    if (nrow(left_df) == 0L || nrow(right_df) == 0L) next
    
    left_keep <- left_df[
      ,
      c(
        "dataset_id", "year", "km_surv", "km_risk", "km_lcl", "km_ucl", "n_risk", "tail_reference_year",
        "mean_surv", "mean_risk", "auc_cd", "brier", "shortfall_model_minus_km",
        "ratio_km_to_model", "missedpos100"
      ),
      drop = FALSE
    ]
    
    right_keep <- right_df[
      ,
      c(
        "dataset_id", "year", "mean_surv", "mean_risk", "auc_cd", "brier", "shortfall_model_minus_km",
        "ratio_km_to_model", "missedpos100"
      ),
      drop = FALSE
    ]
    
    names(left_keep)[9:ncol(left_keep)] <- paste0(names(left_keep)[9:ncol(left_keep)], "_left")
    names(right_keep)[3:ncol(right_keep)] <- paste0(names(right_keep)[3:ncol(right_keep)], "_right")
    
    cmp_df <- merge(left_keep, right_keep, by = c("dataset_id", "year"), all = FALSE, sort = FALSE)
    
    cmp_df$delta_risk_right_minus_left <- cmp_df$mean_risk_right - cmp_df$mean_risk_left
    cmp_df$delta_risk_left_minus_right <- -cmp_df$delta_risk_right_minus_left
    
    cmp_df$delta_auc_right_minus_left <- cmp_df$auc_cd_right - cmp_df$auc_cd_left
    cmp_df$delta_auc_left_minus_right <- -cmp_df$delta_auc_right_minus_left
    
    cmp_df$delta_brier_right_minus_left <- cmp_df$brier_right - cmp_df$brier_left
    cmp_df$delta_brier_left_minus_right <- -cmp_df$delta_brier_right_minus_left
    
    cmp_df$delta_shortfall_right_minus_left <- cmp_df$shortfall_model_minus_km_right - cmp_df$shortfall_model_minus_km_left
    cmp_df$delta_shortfall_left_minus_right <- -cmp_df$delta_shortfall_right_minus_left
    
    cmp_df$delta_ratio_km_to_model_right_minus_left <- cmp_df$ratio_km_to_model_right - cmp_df$ratio_km_to_model_left
    cmp_df$delta_ratio_km_to_model_left_minus_right <- -cmp_df$delta_ratio_km_to_model_right_minus_left
    
    cmp_df$delta_missedpos100_right_minus_left <- cmp_df$missedpos100_right - cmp_df$missedpos100_left
    cmp_df$delta_missedpos100_left_minus_right <- -cmp_df$delta_missedpos100_right_minus_left
    
    meta_rep <- pair_i[rep(1L, nrow(cmp_df)), , drop = FALSE]
    rows[[length(rows) + 1L]] <- cbind(meta_rep, cmp_df)
  }
  
  bind_rows_fill(rows)
}

summarize_pairwise <- function(pair_registry_df, pair_yearly_df, model_summary_df) {
  rows <- list()
  
  for (i in seq_len(nrow(pair_registry_df))) {
    pair_i <- pair_registry_df[i, , drop = FALSE]
    sub <- pair_yearly_df[pair_yearly_df$pair_id == pair_i$pair_id, , drop = FALSE]
    if (nrow(sub) == 0L) next
    
    left_sum <- model_summary_df[model_summary_df$model_id_std == pair_i$left_model_id, , drop = FALSE]
    right_sum <- model_summary_df[model_summary_df$model_id_std == pair_i$right_model_id, , drop = FALSE]
    
    primary_risk_col <- if (pair_i$primary_sign == "right_minus_left") "delta_risk_right_minus_left" else "delta_risk_left_minus_right"
    primary_shortfall_col <- if (pair_i$primary_sign == "right_minus_left") "delta_shortfall_right_minus_left" else "delta_shortfall_left_minus_right"
    
    finite_primary <- which(is.finite(sub[[primary_risk_col]]))
    max_abs_delta_risk <- if (length(finite_primary) > 0L) max(abs(sub[[primary_risk_col]][finite_primary])) else NA_real_
    year_max_abs_delta_risk <- if (length(finite_primary) > 0L) sub$year[finite_primary][which.max(abs(sub[[primary_risk_col]][finite_primary]))] else NA_real_
    
    left_tail_gap <- if (nrow(left_sum) > 0L) left_sum$tail_gap[1] else NA_real_
    right_tail_gap <- if (nrow(right_sum) > 0L) right_sum$tail_gap[1] else NA_real_
    left_tail_abs_gap <- if (nrow(left_sum) > 0L) left_sum$tail_abs_gap[1] else NA_real_
    right_tail_abs_gap <- if (nrow(right_sum) > 0L) right_sum$tail_abs_gap[1] else NA_real_
    
    delta_ibs_right_minus_left <- if (nrow(left_sum) > 0L && nrow(right_sum) > 0L) right_sum$ibs_available[1] - left_sum$ibs_available[1] else NA_real_
    delta_auc_mean_right_minus_left <- if (nrow(left_sum) > 0L && nrow(right_sum) > 0L) right_sum$auc_mean_available[1] - left_sum$auc_mean_available[1] else NA_real_
    delta_brier_mean_right_minus_left <- if (nrow(left_sum) > 0L && nrow(right_sum) > 0L) right_sum$brier_mean_available[1] - left_sum$brier_mean_available[1] else NA_real_
    
    delta_ibs_left_minus_right <- -delta_ibs_right_minus_left
    delta_auc_mean_left_minus_right <- -delta_auc_mean_right_minus_left
    delta_brier_mean_left_minus_right <- -delta_brier_mean_right_minus_left
    
    primary_delta_ibs <- if (pair_i$primary_sign == "right_minus_left") delta_ibs_right_minus_left else delta_ibs_left_minus_right
    primary_delta_auc_mean <- if (pair_i$primary_sign == "right_minus_left") delta_auc_mean_right_minus_left else delta_auc_mean_left_minus_right
    primary_delta_brier_mean <- if (pair_i$primary_sign == "right_minus_left") delta_brier_mean_right_minus_left else delta_brier_mean_left_minus_right
    
    tail_gap_right_minus_left <- right_tail_gap - left_tail_gap
    tail_gap_left_minus_right <- -tail_gap_right_minus_left
    tail_abs_gap_right_minus_left <- right_tail_abs_gap - left_tail_abs_gap
    tail_abs_gap_left_minus_right <- -tail_abs_gap_right_minus_left
    
    rows[[length(rows) + 1L]] <- data.frame(
      pair_id = pair_i$pair_id,
      section_code = pair_i$section_code,
      contrast_label = pair_i$contrast_label,
      dataset_id = pair_i$dataset_id,
      family_std = pair_i$family_std,
      lane_std = pair_i$lane_std,
      left_model_id = pair_i$left_model_id,
      right_model_id = pair_i$right_model_id,
      left_model_class = pair_i$left_model_class,
      right_model_class = pair_i$right_model_class,
      left_cov_spec = pair_i$left_cov_spec,
      right_cov_spec = pair_i$right_cov_spec,
      primary_sign = pair_i$primary_sign,
      comparison_note = pair_i$comparison_note,
      years_compared_n = nrow(sub),
      primary_delta_risk_1y = safe_value_at_year(sub, 1, primary_risk_col),
      primary_delta_risk_5y = safe_value_at_year(sub, 5, primary_risk_col),
      primary_delta_risk_10y = safe_value_at_year(sub, 10, primary_risk_col),
      primary_delta_shortfall_1y = safe_value_at_year(sub, 1, primary_shortfall_col),
      primary_delta_shortfall_5y = safe_value_at_year(sub, 5, primary_shortfall_col),
      primary_delta_shortfall_10y = safe_value_at_year(sub, 10, primary_shortfall_col),
      max_abs_primary_delta_risk = max_abs_delta_risk,
      year_max_abs_primary_delta_risk = year_max_abs_delta_risk,
      primary_delta_ibs_available = primary_delta_ibs,
      primary_delta_auc_mean_available = primary_delta_auc_mean,
      primary_delta_brier_mean_available = primary_delta_brier_mean,
      tail_gap_right_minus_left = tail_gap_right_minus_left,
      tail_gap_left_minus_right = tail_gap_left_minus_right,
      primary_tail_gap_delta = if (pair_i$primary_sign == "right_minus_left") tail_gap_right_minus_left else tail_gap_left_minus_right,
      tail_abs_gap_right_minus_left = tail_abs_gap_right_minus_left,
      tail_abs_gap_left_minus_right = tail_abs_gap_left_minus_right,
      primary_tail_abs_gap_delta = if (pair_i$primary_sign == "right_minus_left") tail_abs_gap_right_minus_left else tail_abs_gap_left_minus_right,
      AIC_right_minus_left = if (nrow(left_sum) > 0L && nrow(right_sum) > 0L) right_sum$AIC[1] - left_sum$AIC[1] else NA_real_,
      BIC_right_minus_left = if (nrow(left_sum) > 0L && nrow(right_sum) > 0L) right_sum$BIC[1] - left_sum$BIC[1] else NA_real_,
      AIC_left_minus_right = if (nrow(left_sum) > 0L && nrow(right_sum) > 0L) left_sum$AIC[1] - right_sum$AIC[1] else NA_real_,
      BIC_left_minus_right = if (nrow(left_sum) > 0L && nrow(right_sum) > 0L) left_sum$BIC[1] - right_sum$BIC[1] else NA_real_,
      left_boundary_flag = if (nrow(left_sum) > 0L) left_sum$boundary_flag[1] else NA,
      right_boundary_flag = if (nrow(right_sum) > 0L) right_sum$boundary_flag[1] else NA,
      left_convergence_flag = if (nrow(left_sum) > 0L) left_sum$convergence_flag[1] else NA,
      right_convergence_flag = if (nrow(right_sum) > 0L) right_sum$convergence_flag[1] else NA,
      left_cure_fraction_mean = if (nrow(left_sum) > 0L) left_sum$cure_fraction_mean[1] else NA_real_,
      right_cure_fraction_mean = if (nrow(right_sum) > 0L) right_sum$cure_fraction_mean[1] else NA_real_,
      boundary_test_p = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows_fill(rows)
}

build_lane_sensitivity <- function(model_summary_df, model_yearly_df) {
  aft_families <- c("exponential", "weibull", "loglogistic", "lognormal")
  covspec_map <- list(
    PNU = c("C00", "C10", "C01", "C11"),
    SNU = c("C00", "C10", "C01", "C11"),
    merged = c("C00S0", "C10S0", "C01S0", "C11S0", "C00S1", "C10S1", "C01S1", "C11S1")
  )
  
  yearly_rows <- list()
  summary_rows <- list()
  
  for (dataset_id in names(covspec_map)) {
    for (cov_spec_i in covspec_map[[dataset_id]]) {
      cox_id <- lookup_model_id(
        model_summary_df = model_summary_df,
        dataset_id = dataset_id,
        model_class = "cure_MLE",
        family_std = "cox",
        cov_spec = cov_spec_i,
        lane_std = "PH_semiparam"
      )
      
      aft_ids <- unique(stats::na.omit(unname(vapply(
        aft_families,
        function(fam) {
          lookup_model_id(
            model_summary_df = model_summary_df,
            dataset_id = dataset_id,
            model_class = "cure_MLE",
            family_std = fam,
            cov_spec = cov_spec_i,
            lane_std = "AFT_param"
          )
        },
        FUN.VALUE = character(1)
      ))))
      
      aft_ids <- aft_ids[nzchar(aft_ids)]
      
      if (is.na(cox_id) || length(aft_ids) == 0L) next
      
      cox_y <- model_yearly_df[model_yearly_df$model_id_std == cox_id, , drop = FALSE]
      aft_y <- model_yearly_df[model_yearly_df$model_id_std %in% aft_ids, , drop = FALSE]
      if (nrow(cox_y) == 0L || nrow(aft_y) == 0L) next
      
      years <- sort(unique(cox_y$year))
      cov_year_rows <- list()
      
      for (yr in years) {
        cox_row <- cox_y[cox_y$year == yr, , drop = FALSE]
        aft_row <- aft_y[aft_y$year == yr, , drop = FALSE]
        if (nrow(cox_row) == 0L || nrow(aft_row) == 0L) next
        
        aft_min_risk <- safe_min(aft_row$mean_risk)
        aft_max_risk <- safe_max(aft_row$mean_risk)
        aft_mean_risk <- safe_mean(aft_row$mean_risk)
        
        aft_min_auc <- safe_min(aft_row$auc_cd)
        aft_max_auc <- safe_max(aft_row$auc_cd)
        aft_mean_auc <- safe_mean(aft_row$auc_cd)
        
        aft_min_brier <- safe_min(aft_row$brier)
        aft_max_brier <- safe_max(aft_row$brier)
        aft_mean_brier <- safe_mean(aft_row$brier)
        
        aft_min_shortfall <- safe_min(aft_row$shortfall_model_minus_km)
        aft_max_shortfall <- safe_max(aft_row$shortfall_model_minus_km)
        aft_mean_shortfall <- safe_mean(aft_row$shortfall_model_minus_km)
        
        cov_year_rows[[length(cov_year_rows) + 1L]] <- data.frame(
          section_code = "6F",
          dataset_id = dataset_id,
          cov_spec = cov_spec_i,
          cox_model_id = cox_id,
          aft_model_ids = collapse_unique(aft_ids),
          aft_model_n = length(aft_ids),
          year = yr,
          cox_mean_risk = cox_row$mean_risk[1],
          aft_min_risk = aft_min_risk,
          aft_max_risk = aft_max_risk,
          aft_mean_risk = aft_mean_risk,
          aft_range_width_risk = aft_max_risk - aft_min_risk,
          cox_gap_from_aft_risk_range = gap_to_range(cox_row$mean_risk[1], aft_min_risk, aft_max_risk),
          delta_risk_cox_minus_aft_mean = cox_row$mean_risk[1] - aft_mean_risk,
          cox_auc_cd = cox_row$auc_cd[1],
          aft_min_auc = aft_min_auc,
          aft_max_auc = aft_max_auc,
          aft_mean_auc = aft_mean_auc,
          cox_gap_from_aft_auc_range = gap_to_range(cox_row$auc_cd[1], aft_min_auc, aft_max_auc),
          cox_brier = cox_row$brier[1],
          aft_min_brier = aft_min_brier,
          aft_max_brier = aft_max_brier,
          aft_mean_brier = aft_mean_brier,
          cox_gap_from_aft_brier_range = gap_to_range(cox_row$brier[1], aft_min_brier, aft_max_brier),
          cox_shortfall = cox_row$shortfall_model_minus_km[1],
          aft_min_shortfall = aft_min_shortfall,
          aft_max_shortfall = aft_max_shortfall,
          aft_mean_shortfall = aft_mean_shortfall,
          cox_gap_from_aft_shortfall_range = gap_to_range(cox_row$shortfall_model_minus_km[1], aft_min_shortfall, aft_max_shortfall),
          stringsAsFactors = FALSE
        )
      }
      
      cov_year_df <- bind_rows_fill(cov_year_rows)
      yearly_rows[[length(yearly_rows) + 1L]] <- cov_year_df
      
      cox_sum <- model_summary_df[model_summary_df$model_id_std == cox_id, , drop = FALSE]
      aft_sum <- model_summary_df[model_summary_df$model_id_std %in% aft_ids, , drop = FALSE]
      
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        section_code = "6F",
        dataset_id = dataset_id,
        cov_spec = cov_spec_i,
        cox_model_id = cox_id,
        aft_model_ids = collapse_unique(aft_ids),
        aft_model_n = length(aft_ids),
        cox_gap_risk_1y = safe_value_at_year(cov_year_df, 1, "cox_gap_from_aft_risk_range"),
        cox_gap_risk_5y = safe_value_at_year(cov_year_df, 5, "cox_gap_from_aft_risk_range"),
        cox_gap_risk_10y = safe_value_at_year(cov_year_df, 10, "cox_gap_from_aft_risk_range"),
        years_outside_aft_risk_range_n = sum(is.finite(cov_year_df$cox_gap_from_aft_risk_range) & abs(cov_year_df$cox_gap_from_aft_risk_range) > 0),
        pct_years_outside_aft_risk_range = mean(is.finite(cov_year_df$cox_gap_from_aft_risk_range) & abs(cov_year_df$cox_gap_from_aft_risk_range) > 0),
        max_abs_cox_gap_from_aft_risk_range = if (any(is.finite(cov_year_df$cox_gap_from_aft_risk_range))) max(abs(cov_year_df$cox_gap_from_aft_risk_range[is.finite(cov_year_df$cox_gap_from_aft_risk_range)])) else NA_real_,
        year_of_max_abs_cox_gap_risk = if (any(is.finite(cov_year_df$cox_gap_from_aft_risk_range))) cov_year_df$year[which.max(abs(cov_year_df$cox_gap_from_aft_risk_range))] else NA_real_,
        mean_aft_range_width_risk = safe_mean(cov_year_df$aft_range_width_risk),
        cox_ibs_available = if (nrow(cox_sum) > 0L) cox_sum$ibs_available[1] else NA_real_,
        aft_min_ibs_available = safe_min(aft_sum$ibs_available),
        aft_max_ibs_available = safe_max(aft_sum$ibs_available),
        aft_mean_ibs_available = safe_mean(aft_sum$ibs_available),
        cox_gap_from_aft_ibs_range = if (nrow(cox_sum) > 0L) gap_to_range(cox_sum$ibs_available[1], safe_min(aft_sum$ibs_available), safe_max(aft_sum$ibs_available)) else NA_real_,
        cox_auc_mean_available = if (nrow(cox_sum) > 0L) cox_sum$auc_mean_available[1] else NA_real_,
        aft_min_auc_mean_available = safe_min(aft_sum$auc_mean_available),
        aft_max_auc_mean_available = safe_max(aft_sum$auc_mean_available),
        aft_mean_auc_mean_available = safe_mean(aft_sum$auc_mean_available),
        cox_gap_from_aft_auc_mean_range = if (nrow(cox_sum) > 0L) gap_to_range(cox_sum$auc_mean_available[1], safe_min(aft_sum$auc_mean_available), safe_max(aft_sum$auc_mean_available)) else NA_real_,
        cox_brier_mean_available = if (nrow(cox_sum) > 0L) cox_sum$brier_mean_available[1] else NA_real_,
        aft_min_brier_mean_available = safe_min(aft_sum$brier_mean_available),
        aft_max_brier_mean_available = safe_max(aft_sum$brier_mean_available),
        aft_mean_brier_mean_available = safe_mean(aft_sum$brier_mean_available),
        cox_gap_from_aft_brier_mean_range = if (nrow(cox_sum) > 0L) gap_to_range(cox_sum$brier_mean_available[1], safe_min(aft_sum$brier_mean_available), safe_max(aft_sum$brier_mean_available)) else NA_real_,
        comparison_note = "AFT cure family envelope versus Cox-latency cure model within same dataset and covariate structure",
        stringsAsFactors = FALSE
      )
    }
  }
  
  list(
    yearly = bind_rows_fill(yearly_rows),
    summary = bind_rows_fill(summary_rows)
  )
}

# 🔴 Import: Step4 and Step5 exported artifacts ===============================
## 🟠 Read: raw input and export bundles ===============================
raw_input <- read_csv_safely(DATA_PATH)

step4_analysis_data <- read_csv_safely(STEP4_ANALYSIS_DATA_PATH)
step4_model_registry <- read_csv_safely(STEP4_MODEL_REGISTRY_PATH)
step4_subject_predictions <- read_csv_safely(STEP4_SUBJECT_PRED_PATH)
step4_mean_predictions <- read_csv_safely(STEP4_MEAN_PRED_PATH)
step4_performance_yearly <- read_csv_safely(STEP4_PERF_YEARLY_PATH)
step4_performance_summary <- read_csv_safely(STEP4_PERF_SUMMARY_PATH)
step4_fit_bundle <- readRDS(STEP4_FIT_BUNDLE_PATH)
step4_manifest <- read_csv_safely(STEP4_MANIFEST_PATH)

step5_model_grid <- read_csv_safely(STEP5_MODEL_GRID_PATH)
step5_preproc <- read_csv_safely(STEP5_PREPROC_PATH)
step5_analysis_spec <- read_csv_safely(STEP5_ANALYSIS_SPEC_PATH)
step5_fit_registry <- read_csv_safely(STEP5_FIT_REGISTRY_PATH)
step5_coefficients <- read_csv_safely(STEP5_COEFFICIENTS_PATH)
step5_cure_only <- read_csv_safely(STEP5_CURE_ONLY_PATH)
step5_mean_predictions <- read_csv_safely(STEP5_MEAN_PRED_PATH)
step5_latency_summary <- read_csv_safely(STEP5_LATENCY_SUMMARY_PATH)
step5_subject_predictions <- read_csv_safely(STEP5_SUBJECT_PRED_PATH)
step5_master <- readRDS(STEP5_MASTER_RDS_PATH)
step5_manifest <- read_csv_safely(STEP5_MANIFEST_PATH)

import_audit <- bind_rows_fill(list(
  make_audit_row("step4", "analysis_data", STEP4_ANALYSIS_DATA_PATH, step4_analysis_data, TRUE),
  make_audit_row("step4", "model_registry", STEP4_MODEL_REGISTRY_PATH, step4_model_registry, TRUE),
  make_audit_row("step4", "subject_predictions", STEP4_SUBJECT_PRED_PATH, step4_subject_predictions, TRUE),
  make_audit_row("step4", "mean_predictions", STEP4_MEAN_PRED_PATH, step4_mean_predictions, TRUE),
  make_audit_row("step4", "performance_yearly", STEP4_PERF_YEARLY_PATH, step4_performance_yearly, TRUE),
  make_audit_row("step4", "performance_summary", STEP4_PERF_SUMMARY_PATH, step4_performance_summary, TRUE),
  make_audit_row("step4", "fit_bundle", STEP4_FIT_BUNDLE_PATH, step4_fit_bundle, TRUE),
  make_audit_row("step4", "manifest", STEP4_MANIFEST_PATH, step4_manifest, TRUE),
  make_audit_row("step5", "model_grid", STEP5_MODEL_GRID_PATH, step5_model_grid, TRUE),
  make_audit_row("step5", "preprocessing", STEP5_PREPROC_PATH, step5_preproc, TRUE),
  make_audit_row("step5", "analysis_spec", STEP5_ANALYSIS_SPEC_PATH, step5_analysis_spec, TRUE),
  make_audit_row("step5", "fit_registry", STEP5_FIT_REGISTRY_PATH, step5_fit_registry, TRUE),
  make_audit_row("step5", "coefficients", STEP5_COEFFICIENTS_PATH, step5_coefficients, TRUE),
  make_audit_row("step5", "cure_only", STEP5_CURE_ONLY_PATH, step5_cure_only, TRUE),
  make_audit_row("step5", "mean_predictions", STEP5_MEAN_PRED_PATH, step5_mean_predictions, TRUE),
  make_audit_row("step5", "latency_summary", STEP5_LATENCY_SUMMARY_PATH, step5_latency_summary, TRUE),
  make_audit_row("step5", "subject_predictions", STEP5_SUBJECT_PRED_PATH, step5_subject_predictions, TRUE),
  make_audit_row("step5", "master_rds", STEP5_MASTER_RDS_PATH, step5_master, TRUE),
  make_audit_row("step5", "manifest", STEP5_MANIFEST_PATH, step5_manifest, TRUE)
))

## 🟠 Transform: canonical model catalogs and analysis branches ===============================
step4_analysis_data$dataset_id <- standardize_dataset(step4_analysis_data$dataset_id)
step4_analysis_data$subject_key <- as.character(step4_analysis_data$subject_uid)
step4_analysis_data <- to_numeric_if_present(step4_analysis_data, c("time_years_model", "event_transition"))

analysis_core <- unique(
  step4_analysis_data[, c("dataset_id", "subject_key", "time_years_model", "event_transition"), drop = FALSE]
)
analysis_core <- analysis_core[order(analysis_core$dataset_id, analysis_core$subject_key), , drop = FALSE]
analysis_branches <- split(analysis_core, analysis_core$dataset_id)

step4_registry_std <- standardize_step4_registry(step4_model_registry)
step4_mean_std <- standardize_step4_mean(step4_mean_predictions)
step4_perf_yearly_std <- standardize_step4_perf_yearly(step4_performance_yearly)
step4_perf_summary_std <- standardize_step4_perf_summary(step4_performance_summary)
step4_subject_pred_std <- standardize_step4_subject_predictions(step4_subject_predictions)

step5_registry_std <- standardize_step5_registry(step5_fit_registry)
step5_mean_std <- standardize_step5_mean(step5_mean_predictions, step5_registry_std)
step5_subject_pred_std <- standardize_step5_subject_predictions(step5_subject_predictions, step5_registry_std)
step5_cure_only_std <- standardize_step5_cure_only(step5_cure_only)
step5_latency_std <- standardize_step5_latency_summary(step5_latency_summary)
step5_coef_std <- standardize_step5_coefficients(step5_coefficients, step5_registry_std)

step4_mean_validation_max_abs_diff <- validate_mean_predictions(step4_subject_pred_std, step4_mean_std)
step5_mean_validation_max_abs_diff <- validate_mean_predictions(step5_subject_pred_std, step5_mean_std)

# 🔴 Compute: model-level metrics and harmonized tables ===============================
## 🟠 Build: Kaplan-Meier benchmarks by dataset ===============================
km_benchmark <- bind_rows_fill(lapply(analysis_branches, km_yearly_summary, years = PREDICTION_YEARS))
km_benchmark <- km_benchmark[order(km_benchmark$dataset_id, km_benchmark$year), , drop = FALSE]

## 🟠 Build: Step4 harmonized model metrics ===============================
step4_model_yearly <- merge(
  step4_mean_std[, c("dataset_id", "subgroup_id", "model_id_std", "model_class", "lane_std", "family_std", "cov_spec", "source_step", "year", "mean_surv", "mean_risk"), drop = FALSE],
  step4_perf_yearly_std[, c("dataset_id", "model_id_std", "year", "auc_cd", "brier"), drop = FALSE],
  by = c("dataset_id", "model_id_std", "year"),
  all = TRUE,
  sort = FALSE
)
step4_model_yearly <- attach_km_metrics(step4_model_yearly, km_benchmark)
step4_yearly_bits <- summarize_model_yearly(step4_model_yearly)

step4_model_summary <- merge(
  step4_registry_std[, c(
    "dataset_id", "model_id_std", "model_class", "lane_std", "family_std", "cov_spec", "source_step",
    "fit_success", "convergence_flag", "singular_hessian_flag", "boundary_flag", "n_subject",
    "n_event", "n_censored", "n_params", "logLik_full", "AIC", "BIC", "partial_logLik",
    "AIC_partial", "likelihood_type", "interaction_flag", "site_flag", "warning_messages", "error_message"
  ), drop = FALSE],
  step4_perf_summary_std[, c("dataset_id", "model_id_std", "ibs_1to10", "ibs_available", "ibs_available_from", "ibs_available_to", "auc_mean_available", "brier_mean_available", "n_evaluable_years"), drop = FALSE],
  by = c("dataset_id", "model_id_std"),
  all = TRUE,
  sort = FALSE
)
step4_model_summary$n_parameters <- step4_model_summary$n_params
step4_model_summary$incidence_int_flag <- NA_integer_
step4_model_summary$latency_int_flag <- NA_integer_
step4_model_summary$cure_fraction_mean <- NA_real_
step4_model_summary$cure_fraction_median <- NA_real_
step4_model_summary$cure_fraction_min <- NA_real_
step4_model_summary$cure_fraction_max <- NA_real_
step4_model_summary$uncured_probability_mean <- NA_real_
step4_model_summary$incidence_target <- NA_character_
step4_model_summary <- merge(step4_model_summary, step4_yearly_bits, by = c("dataset_id", "model_id_std"), all = TRUE, sort = FALSE)

## 🟠 Build: Step5 apparent performance and harmonized model metrics ===============================
if (isTRUE(RECOMPUTE_STEP5_PERFORMANCE)) {
  step5_perf <- compute_step5_performance(step5_subject_pred_std, step5_registry_std, analysis_branches, PREDICTION_YEARS)
} else {
  step5_perf <- list(yearly = data.frame(stringsAsFactors = FALSE), summary = data.frame(stringsAsFactors = FALSE))
}

step5_model_yearly <- merge(
  step5_mean_std[, c("dataset_id", "subgroup_id", "model_id_std", "model_class", "lane_std", "family_std", "cov_spec", "source_step", "year", "mean_surv", "mean_risk"), drop = FALSE],
  step5_perf$yearly[, c("dataset_id", "model_id_std", "year", "auc_cd", "brier"), drop = FALSE],
  by = c("dataset_id", "model_id_std", "year"),
  all = TRUE,
  sort = FALSE
)
step5_model_yearly <- attach_km_metrics(step5_model_yearly, km_benchmark)
step5_yearly_bits <- summarize_model_yearly(step5_model_yearly)
step5_latency_bits <- summarize_latency_yearly(step5_latency_std)

step5_model_summary <- merge(
  step5_registry_std[, c(
    "dataset_id", "model_id_std", "model_class", "lane_std", "family_std", "cov_spec", "source_step",
    "fit_success", "convergence_flag", "boundary_flag", "n_subject", "n_event", "n_censored",
    "n_parameters", "logLik_full", "AIC", "BIC", "likelihood_type", "partial_logLik",
    "AIC_partial", "incidence_int_flag", "latency_int_flag", "site_flag", "warning_messages",
    "error_message", "incidence_target", "mean_cure_fraction"
  ), drop = FALSE],
  step5_perf$summary[, c("dataset_id", "model_id_std", "ibs_1to10", "ibs_available", "ibs_available_from", "ibs_available_to", "auc_mean_available", "brier_mean_available", "n_evaluable_years"), drop = FALSE],
  by = c("dataset_id", "model_id_std"),
  all = TRUE,
  sort = FALSE
)

step5_model_summary <- merge(
  step5_model_summary,
  step5_cure_only_std[, c("dataset_id", "model_id_std", "cure_fraction_mean", "cure_fraction_median", "cure_fraction_min", "cure_fraction_max", "uncured_probability_mean"), drop = FALSE],
  by = c("dataset_id", "model_id_std"),
  all = TRUE,
  sort = FALSE
)
step5_model_summary <- merge(step5_model_summary, step5_yearly_bits, by = c("dataset_id", "model_id_std"), all = TRUE, sort = FALSE)
step5_model_summary <- merge(step5_model_summary, step5_latency_bits, by = c("dataset_id", "model_id_std"), all = TRUE, sort = FALSE)

## 🟠 Build: unified coefficient catalog ===============================
step4_coef_std <- extract_step4_coefficients(step4_fit_bundle, step4_registry_std)
coefficients_combined <- bind_rows_fill(list(step4_coef_std, step5_coef_std))
coefficients_combined$term_role <- ifelse(
  coefficients_combined$term %in% c("age_s:sex_num", "age_sex_int"),
  "age_sex_interaction",
  ifelse(
    grepl("^site", coefficients_combined$term),
    "site_effect",
    ifelse(coefficients_combined$component %in% c("incidence", "latency"), paste0(coefficients_combined$component, "_term"), "other")
  )
)
coefficients_combined <- coefficients_combined[order(coefficients_combined$dataset_id, coefficients_combined$model_class, coefficients_combined$model_id_std, coefficients_combined$component, coefficients_combined$term), , drop = FALSE]

## 🟠 Build: unified model dashboards ===============================
model_yearly_metrics <- bind_rows_fill(list(step4_model_yearly, step5_model_yearly))
model_yearly_metrics$model_label <- paste(model_yearly_metrics$dataset_id, model_yearly_metrics$model_class, model_yearly_metrics$family_std, model_yearly_metrics$cov_spec, sep = "__")
model_yearly_metrics <- model_yearly_metrics[order(model_yearly_metrics$dataset_id, model_yearly_metrics$model_class, model_yearly_metrics$family_std, model_yearly_metrics$cov_spec, model_yearly_metrics$year), , drop = FALSE]

model_summary_metrics <- bind_rows_fill(list(step4_model_summary, step5_model_summary))
model_summary_metrics$model_label <- paste(model_summary_metrics$dataset_id, model_summary_metrics$model_class, model_summary_metrics$family_std, model_summary_metrics$cov_spec, sep = "__")
model_summary_metrics <- model_summary_metrics[order(model_summary_metrics$dataset_id, model_summary_metrics$model_class, model_summary_metrics$family_std, model_summary_metrics$cov_spec), , drop = FALSE]

# 🔴 Compare: Step6 evidence matrices ===============================
## 🟠 Assemble: no-cure internal contrasts and cure matched contrasts ===============================
pair_rows <- list()

add_pair_row <- function(section_code,
                         contrast_label,
                         dataset_id,
                         family_std,
                         lane_std,
                         left_model_id,
                         right_model_id,
                         left_model_class,
                         right_model_class,
                         left_cov_spec,
                         right_cov_spec,
                         primary_sign,
                         comparison_note) {
  if (any(is.na(c(left_model_id, right_model_id)))) return(invisible(NULL))
  if (!nzchar(left_model_id) || !nzchar(right_model_id)) return(invisible(NULL))
  
  pair_rows[[length(pair_rows) + 1L]] <<- data.frame(
    pair_id = sprintf("PAIR_%04d", length(pair_rows) + 1L),
    section_code = section_code,
    contrast_label = contrast_label,
    dataset_id = dataset_id,
    family_std = family_std,
    lane_std = lane_std,
    left_model_id = left_model_id,
    right_model_id = right_model_id,
    left_model_class = left_model_class,
    right_model_class = right_model_class,
    left_cov_spec = left_cov_spec,
    right_cov_spec = right_cov_spec,
    primary_sign = primary_sign,
    comparison_note = comparison_note,
    stringsAsFactors = FALSE
  )
}

all_nocure_families <- c("exponential", "weibull", "loglogistic", "lognormal", "cox")
all_cure_families <- c("exponential", "weibull", "loglogistic", "lognormal", "cox")

for (dataset_id in c("PNU", "SNU")) {
  for (family_std in all_nocure_families) {
    lane_std <- if (family_std == "cox") "PH_semiparam" else "AFT_param"
    
    add_pair_row(
      section_code = "6A",
      contrast_label = "nc_internal_interaction",
      dataset_id = dataset_id,
      family_std = family_std,
      lane_std = lane_std,
      left_model_id = lookup_model_id(model_summary_metrics, dataset_id, "nocure", family_std, "N0", lane_std),
      right_model_id = lookup_model_id(model_summary_metrics, dataset_id, "nocure", family_std, "N1", lane_std),
      left_model_class = "nocure",
      right_model_class = "nocure",
      left_cov_spec = "N0",
      right_cov_spec = "N1",
      primary_sign = "right_minus_left",
      comparison_note = "No-cure family 내부 비교: interaction 추가 효과 (N1 - N0)"
    )
  }
}

for (family_std in all_nocure_families) {
  lane_std <- if (family_std == "cox") "PH_semiparam" else "AFT_param"
  
  add_pair_row(
    section_code = "6A|6E",
    contrast_label = "nc_site_effect",
    dataset_id = "merged",
    family_std = family_std,
    lane_std = lane_std,
    left_model_id = lookup_model_id(model_summary_metrics, "merged", "nocure", family_std, "N0S0", lane_std),
    right_model_id = lookup_model_id(model_summary_metrics, "merged", "nocure", family_std, "N0S1", lane_std),
    left_model_class = "nocure",
    right_model_class = "nocure",
    left_cov_spec = "N0S0",
    right_cov_spec = "N0S1",
    primary_sign = "right_minus_left",
    comparison_note = "No-cure merged site 효과: no interaction 고정 (N0S1 - N0S0)"
  )
  
  add_pair_row(
    section_code = "6A|6E",
    contrast_label = "nc_site_effect",
    dataset_id = "merged",
    family_std = family_std,
    lane_std = lane_std,
    left_model_id = lookup_model_id(model_summary_metrics, "merged", "nocure", family_std, "N1S0", lane_std),
    right_model_id = lookup_model_id(model_summary_metrics, "merged", "nocure", family_std, "N1S1", lane_std),
    left_model_class = "nocure",
    right_model_class = "nocure",
    left_cov_spec = "N1S0",
    right_cov_spec = "N1S1",
    primary_sign = "right_minus_left",
    comparison_note = "No-cure merged site 효과: interaction 유지 (N1S1 - N1S0)"
  )
  
  add_pair_row(
    section_code = "6A",
    contrast_label = "nc_internal_interaction",
    dataset_id = "merged",
    family_std = family_std,
    lane_std = lane_std,
    left_model_id = lookup_model_id(model_summary_metrics, "merged", "nocure", family_std, "N0S0", lane_std),
    right_model_id = lookup_model_id(model_summary_metrics, "merged", "nocure", family_std, "N1S0", lane_std),
    left_model_class = "nocure",
    right_model_class = "nocure",
    left_cov_spec = "N0S0",
    right_cov_spec = "N1S0",
    primary_sign = "right_minus_left",
    comparison_note = "No-cure merged interaction 효과: site 없음 고정 (N1S0 - N0S0)"
  )
  
  add_pair_row(
    section_code = "6A",
    contrast_label = "nc_internal_interaction",
    dataset_id = "merged",
    family_std = family_std,
    lane_std = lane_std,
    left_model_id = lookup_model_id(model_summary_metrics, "merged", "nocure", family_std, "N0S1", lane_std),
    right_model_id = lookup_model_id(model_summary_metrics, "merged", "nocure", family_std, "N1S1", lane_std),
    left_model_class = "nocure",
    right_model_class = "nocure",
    left_cov_spec = "N0S1",
    right_cov_spec = "N1S1",
    primary_sign = "right_minus_left",
    comparison_note = "No-cure merged interaction 효과: site 포함 고정 (N1S1 - N0S1)"
  )
}

for (dataset_id in c("PNU", "SNU")) {
  for (family_std in all_cure_families) {
    lane_std <- if (family_std == "cox") "PH_semiparam" else "AFT_param"
    
    add_pair_row(
      section_code = "6B",
      contrast_label = "nc_vs_cure_latency_matched",
      dataset_id = dataset_id,
      family_std = family_std,
      lane_std = lane_std,
      left_model_id = lookup_model_id(model_summary_metrics, dataset_id, "nocure", family_std, "N0", lane_std),
      right_model_id = lookup_model_id(model_summary_metrics, dataset_id, "cure_MLE", family_std, "C00", lane_std),
      left_model_class = "nocure",
      right_model_class = "cure_MLE",
      left_cov_spec = "N0",
      right_cov_spec = "C00",
      primary_sign = "left_minus_right",
      comparison_note = "Latency L0 매칭: no-cure N0 versus cure C00"
    )
    
    add_pair_row(
      section_code = "6B",
      contrast_label = "nc_vs_cure_latency_matched",
      dataset_id = dataset_id,
      family_std = family_std,
      lane_std = lane_std,
      left_model_id = lookup_model_id(model_summary_metrics, dataset_id, "nocure", family_std, "N0", lane_std),
      right_model_id = lookup_model_id(model_summary_metrics, dataset_id, "cure_MLE", family_std, "C10", lane_std),
      left_model_class = "nocure",
      right_model_class = "cure_MLE",
      left_cov_spec = "N0",
      right_cov_spec = "C10",
      primary_sign = "left_minus_right",
      comparison_note = "Latency L0 매칭: no-cure N0 versus cure C10"
    )
    
    add_pair_row(
      section_code = "6B",
      contrast_label = "nc_vs_cure_latency_matched",
      dataset_id = dataset_id,
      family_std = family_std,
      lane_std = lane_std,
      left_model_id = lookup_model_id(model_summary_metrics, dataset_id, "nocure", family_std, "N1", lane_std),
      right_model_id = lookup_model_id(model_summary_metrics, dataset_id, "cure_MLE", family_std, "C01", lane_std),
      left_model_class = "nocure",
      right_model_class = "cure_MLE",
      left_cov_spec = "N1",
      right_cov_spec = "C01",
      primary_sign = "left_minus_right",
      comparison_note = "Latency L1 매칭: no-cure N1 versus cure C01"
    )
    
    add_pair_row(
      section_code = "6B",
      contrast_label = "nc_vs_cure_latency_matched",
      dataset_id = dataset_id,
      family_std = family_std,
      lane_std = lane_std,
      left_model_id = lookup_model_id(model_summary_metrics, dataset_id, "nocure", family_std, "N1", lane_std),
      right_model_id = lookup_model_id(model_summary_metrics, dataset_id, "cure_MLE", family_std, "C11", lane_std),
      left_model_class = "nocure",
      right_model_class = "cure_MLE",
      left_cov_spec = "N1",
      right_cov_spec = "C11",
      primary_sign = "left_minus_right",
      comparison_note = "Latency L1 매칭: no-cure N1 versus cure C11"
    )
  }
}

for (family_std in all_cure_families) {
  lane_std <- if (family_std == "cox") "PH_semiparam" else "AFT_param"
  
  add_pair_row("6B", "nc_vs_cure_latency_matched", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "nocure", family_std, "N0S0", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C00S0", lane_std),
               "nocure", "cure_MLE", "N0S0", "C00S0", "left_minus_right",
               "Merged latency L0/site S0 매칭: N0S0 versus C00S0")
  
  add_pair_row("6B", "nc_vs_cure_latency_matched", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "nocure", family_std, "N0S0", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C10S0", lane_std),
               "nocure", "cure_MLE", "N0S0", "C10S0", "left_minus_right",
               "Merged latency L0/site S0 매칭: N0S0 versus C10S0")
  
  add_pair_row("6B", "nc_vs_cure_latency_matched", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "nocure", family_std, "N1S0", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C01S0", lane_std),
               "nocure", "cure_MLE", "N1S0", "C01S0", "left_minus_right",
               "Merged latency L1/site S0 매칭: N1S0 versus C01S0")
  
  add_pair_row("6B", "nc_vs_cure_latency_matched", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "nocure", family_std, "N1S0", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C11S0", lane_std),
               "nocure", "cure_MLE", "N1S0", "C11S0", "left_minus_right",
               "Merged latency L1/site S0 매칭: N1S0 versus C11S0")
  
  add_pair_row("6B", "nc_vs_cure_latency_matched", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "nocure", family_std, "N0S1", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C00S1", lane_std),
               "nocure", "cure_MLE", "N0S1", "C00S1", "left_minus_right",
               "Merged latency L0/site S1 매칭: N0S1 versus C00S1")
  
  add_pair_row("6B", "nc_vs_cure_latency_matched", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "nocure", family_std, "N0S1", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C10S1", lane_std),
               "nocure", "cure_MLE", "N0S1", "C10S1", "left_minus_right",
               "Merged latency L0/site S1 매칭: N0S1 versus C10S1")
  
  add_pair_row("6B", "nc_vs_cure_latency_matched", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "nocure", family_std, "N1S1", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C01S1", lane_std),
               "nocure", "cure_MLE", "N1S1", "C01S1", "left_minus_right",
               "Merged latency L1/site S1 매칭: N1S1 versus C01S1")
  
  add_pair_row("6B", "nc_vs_cure_latency_matched", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "nocure", family_std, "N1S1", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C11S1", lane_std),
               "nocure", "cure_MLE", "N1S1", "C11S1", "left_minus_right",
               "Merged latency L1/site S1 매칭: N1S1 versus C11S1")
}

## 🟠 Assemble: incidence latency and site contrasts ===============================
for (dataset_id in c("PNU", "SNU")) {
  for (family_std in all_cure_families) {
    lane_std <- if (family_std == "cox") "PH_semiparam" else "AFT_param"
    
    add_pair_row("6C", "cure_incidence_interaction", dataset_id, family_std, lane_std,
                 lookup_model_id(model_summary_metrics, dataset_id, "cure_MLE", family_std, "C00", lane_std),
                 lookup_model_id(model_summary_metrics, dataset_id, "cure_MLE", family_std, "C10", lane_std),
                 "cure_MLE", "cure_MLE", "C00", "C10", "right_minus_left",
                 "Cure incidence interaction 효과: latency L0 고정 (C10 - C00)")
    
    add_pair_row("6C", "cure_incidence_interaction", dataset_id, family_std, lane_std,
                 lookup_model_id(model_summary_metrics, dataset_id, "cure_MLE", family_std, "C01", lane_std),
                 lookup_model_id(model_summary_metrics, dataset_id, "cure_MLE", family_std, "C11", lane_std),
                 "cure_MLE", "cure_MLE", "C01", "C11", "right_minus_left",
                 "Cure incidence interaction 효과: latency L1 고정 (C11 - C01)")
    
    add_pair_row("6D", "cure_latency_interaction", dataset_id, family_std, lane_std,
                 lookup_model_id(model_summary_metrics, dataset_id, "cure_MLE", family_std, "C00", lane_std),
                 lookup_model_id(model_summary_metrics, dataset_id, "cure_MLE", family_std, "C01", lane_std),
                 "cure_MLE", "cure_MLE", "C00", "C01", "right_minus_left",
                 "Cure latency interaction 효과: incidence I0 고정 (C01 - C00)")
    
    add_pair_row("6D", "cure_latency_interaction", dataset_id, family_std, lane_std,
                 lookup_model_id(model_summary_metrics, dataset_id, "cure_MLE", family_std, "C10", lane_std),
                 lookup_model_id(model_summary_metrics, dataset_id, "cure_MLE", family_std, "C11", lane_std),
                 "cure_MLE", "cure_MLE", "C10", "C11", "right_minus_left",
                 "Cure latency interaction 효과: incidence I1 고정 (C11 - C10)")
  }
}

for (family_std in all_cure_families) {
  lane_std <- if (family_std == "cox") "PH_semiparam" else "AFT_param"
  
  add_pair_row("6C", "cure_incidence_interaction", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C00S0", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C10S0", lane_std),
               "cure_MLE", "cure_MLE", "C00S0", "C10S0", "right_minus_left",
               "Cure incidence interaction 효과: site S0, latency L0 고정 (C10S0 - C00S0)")
  
  add_pair_row("6C", "cure_incidence_interaction", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C01S0", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C11S0", lane_std),
               "cure_MLE", "cure_MLE", "C01S0", "C11S0", "right_minus_left",
               "Cure incidence interaction 효과: site S0, latency L1 고정 (C11S0 - C01S0)")
  
  add_pair_row("6C", "cure_incidence_interaction", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C00S1", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C10S1", lane_std),
               "cure_MLE", "cure_MLE", "C00S1", "C10S1", "right_minus_left",
               "Cure incidence interaction 효과: site S1, latency L0 고정 (C10S1 - C00S1)")
  
  add_pair_row("6C", "cure_incidence_interaction", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C01S1", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C11S1", lane_std),
               "cure_MLE", "cure_MLE", "C01S1", "C11S1", "right_minus_left",
               "Cure incidence interaction 효과: site S1, latency L1 고정 (C11S1 - C01S1)")
  
  add_pair_row("6D", "cure_latency_interaction", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C00S0", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C01S0", lane_std),
               "cure_MLE", "cure_MLE", "C00S0", "C01S0", "right_minus_left",
               "Cure latency interaction 효과: site S0, incidence I0 고정 (C01S0 - C00S0)")
  
  add_pair_row("6D", "cure_latency_interaction", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C10S0", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C11S0", lane_std),
               "cure_MLE", "cure_MLE", "C10S0", "C11S0", "right_minus_left",
               "Cure latency interaction 효과: site S0, incidence I1 고정 (C11S0 - C10S0)")
  
  add_pair_row("6D", "cure_latency_interaction", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C00S1", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C01S1", lane_std),
               "cure_MLE", "cure_MLE", "C00S1", "C01S1", "right_minus_left",
               "Cure latency interaction 효과: site S1, incidence I0 고정 (C01S1 - C00S1)")
  
  add_pair_row("6D", "cure_latency_interaction", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C10S1", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C11S1", lane_std),
               "cure_MLE", "cure_MLE", "C10S1", "C11S1", "right_minus_left",
               "Cure latency interaction 효과: site S1, incidence I1 고정 (C11S1 - C10S1)")
  
  add_pair_row("6E", "cure_site_effect", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C00S0", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C00S1", lane_std),
               "cure_MLE", "cure_MLE", "C00S0", "C00S1", "right_minus_left",
               "Cure merged site 효과: incidence I0 latency L0 고정 (C00S1 - C00S0)")
  
  add_pair_row("6E", "cure_site_effect", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C10S0", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C10S1", lane_std),
               "cure_MLE", "cure_MLE", "C10S0", "C10S1", "right_minus_left",
               "Cure merged site 효과: incidence I1 latency L0 고정 (C10S1 - C10S0)")
  
  add_pair_row("6E", "cure_site_effect", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C01S0", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C01S1", lane_std),
               "cure_MLE", "cure_MLE", "C01S0", "C01S1", "right_minus_left",
               "Cure merged site 효과: incidence I0 latency L1 고정 (C01S1 - C01S0)")
  
  add_pair_row("6E", "cure_site_effect", "merged", family_std, lane_std,
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C11S0", lane_std),
               lookup_model_id(model_summary_metrics, "merged", "cure_MLE", family_std, "C11S1", lane_std),
               "cure_MLE", "cure_MLE", "C11S0", "C11S1", "right_minus_left",
               "Cure merged site 효과: incidence I1 latency L1 고정 (C11S1 - C11S0)")
}

pair_registry <- bind_rows_fill(pair_rows)
pair_registry <- pair_registry[order(pair_registry$section_code, pair_registry$dataset_id, pair_registry$family_std, pair_registry$left_cov_spec, pair_registry$right_cov_spec), , drop = FALSE]

pair_yearly <- build_pairwise_yearly(pair_registry, model_yearly_metrics)
pair_yearly <- pair_yearly[order(pair_yearly$section_code, pair_yearly$dataset_id, pair_yearly$family_std, pair_yearly$pair_id, pair_yearly$year), , drop = FALSE]

pair_summary <- summarize_pairwise(pair_registry, pair_yearly, model_summary_metrics)
pair_summary <- pair_summary[order(pair_summary$section_code, pair_summary$dataset_id, pair_summary$family_std, pair_summary$pair_id), , drop = FALSE]

## 🟠 Assemble: AFT envelope versus Cox-latency contrasts ===============================
lane_sensitivity <- build_lane_sensitivity(model_summary_metrics, model_yearly_metrics)
lane_yearly <- lane_sensitivity$yearly
lane_summary <- lane_sensitivity$summary

if (nrow(lane_yearly) > 0L) {
  lane_yearly <- lane_yearly[order(lane_yearly$dataset_id, lane_yearly$cov_spec, lane_yearly$year), , drop = FALSE]
}
if (nrow(lane_summary) > 0L) {
  lane_summary <- lane_summary[order(lane_summary$dataset_id, lane_summary$cov_spec), , drop = FALSE]
}

# 🔴 Export: Step6 source-of-truth files and master bundle ===============================
## 🟠 Write: CSV artifacts and RDS bundle ===============================
analysis_spec <- data.frame(
  analysis_run_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  data_path = DATA_PATH,
  step4_export_dir = STEP4_EXPORT_DIR,
  step4_file_prefix = STEP4_FILE_PREFIX,
  step5_export_dir = STEP5_EXPORT_DIR,
  step6_export_dir = STEP6_EXPORT_DIR,
  prediction_years = paste(PREDICTION_YEARS, collapse = "|"),
  site_reference = SITE_REFERENCE,
  recompute_step5_performance = RECOMPUTE_STEP5_PERFORMANCE,
  raw_rows = nrow(raw_input),
  raw_site_levels = if ("site" %in% names(raw_input)) collapse_unique(raw_input$site) else "",
  raw_unique_site_id_n = if (all(c("site", "id") %in% names(raw_input))) length(unique(paste(raw_input$site, raw_input$id, sep = "::"))) else NA_integer_,
  step4_analysis_rows = nrow(step4_analysis_data),
  step4_models = nrow(step4_registry_std),
  step5_models = nrow(step5_registry_std),
  combined_model_yearly_rows = nrow(model_yearly_metrics),
  combined_model_summary_rows = nrow(model_summary_metrics),
  pair_registry_rows = nrow(pair_registry),
  pair_yearly_rows = nrow(pair_yearly),
  pair_summary_rows = nrow(pair_summary),
  lane_yearly_rows = nrow(lane_yearly),
  lane_summary_rows = nrow(lane_summary),
  step4_mean_prediction_max_abs_diff = step4_mean_validation_max_abs_diff,
  step5_mean_prediction_max_abs_diff = step5_mean_validation_max_abs_diff,
  km_tail_reference_definition = "last whole year <= last observed event time within each dataset branch",
  pair_delta_note = "pair_yearly stores both right-left and left-right deltas; pair_summary marks primary_sign",
  package_survival = as.character(utils::packageVersion("survival")),
  R_version = R.version.string,
  stringsAsFactors = FALSE
)

master_bundle <- list(
  script_stage = "Step6_evidence_matrix",
  analysis_spec = analysis_spec,
  imported = list(
    step4 = list(
      analysis_data = step4_analysis_data,
      model_registry = step4_model_registry,
      subject_predictions = step4_subject_predictions,
      mean_predictions = step4_mean_predictions,
      performance_yearly = step4_performance_yearly,
      performance_summary = step4_performance_summary,
      fit_bundle = step4_fit_bundle,
      manifest = step4_manifest
    ),
    step5 = list(
      model_grid = step5_model_grid,
      preprocessing = step5_preproc,
      analysis_spec = step5_analysis_spec,
      fit_registry = step5_fit_registry,
      coefficients = step5_coefficients,
      cure_only = step5_cure_only,
      mean_predictions = step5_mean_predictions,
      latency_summary = step5_latency_summary,
      subject_predictions = step5_subject_predictions,
      master_rds = step5_master,
      manifest = step5_manifest
    )
  ),
  import_audit = import_audit,
  km_benchmark = km_benchmark,
  model_yearly_metrics = model_yearly_metrics,
  model_summary_metrics = model_summary_metrics,
  coefficients = coefficients_combined,
  pair_registry = pair_registry,
  pair_yearly = pair_yearly,
  pair_summary = pair_summary,
  lane_yearly = lane_yearly,
  lane_summary = lane_summary,
  session_info = paste(capture.output(sessionInfo()), collapse = "\n"),
  created_at = as.character(Sys.time())
)

utils::write.csv(analysis_spec, ANALYSIS_SPEC_OUT, row.names = FALSE, na = "")
utils::write.csv(import_audit, IMPORT_AUDIT_OUT, row.names = FALSE, na = "")
utils::write.csv(km_benchmark, KM_BENCHMARK_OUT, row.names = FALSE, na = "")
utils::write.csv(model_yearly_metrics, MODEL_YEARLY_OUT, row.names = FALSE, na = "")
utils::write.csv(model_summary_metrics, MODEL_SUMMARY_OUT, row.names = FALSE, na = "")
utils::write.csv(pair_registry, PAIR_REGISTRY_OUT, row.names = FALSE, na = "")
utils::write.csv(pair_yearly, PAIR_YEARLY_OUT, row.names = FALSE, na = "")
utils::write.csv(pair_summary, PAIR_SUMMARY_OUT, row.names = FALSE, na = "")
utils::write.csv(coefficients_combined, COEFFICIENTS_OUT, row.names = FALSE, na = "")
utils::write.csv(lane_yearly, LANE_YEARLY_OUT, row.names = FALSE, na = "")
utils::write.csv(lane_summary, LANE_SUMMARY_OUT, row.names = FALSE, na = "")
saveRDS(master_bundle, MASTER_RDS_OUT, compress = "xz")

manifest <- data.frame(
  file_name = c(
    basename(ANALYSIS_SPEC_OUT),
    basename(IMPORT_AUDIT_OUT),
    basename(KM_BENCHMARK_OUT),
    basename(MODEL_YEARLY_OUT),
    basename(MODEL_SUMMARY_OUT),
    basename(PAIR_REGISTRY_OUT),
    basename(PAIR_YEARLY_OUT),
    basename(PAIR_SUMMARY_OUT),
    basename(COEFFICIENTS_OUT),
    basename(LANE_YEARLY_OUT),
    basename(LANE_SUMMARY_OUT),
    basename(MASTER_RDS_OUT),
    basename(MANIFEST_OUT)
  ),
  file_path = c(
    ANALYSIS_SPEC_OUT,
    IMPORT_AUDIT_OUT,
    KM_BENCHMARK_OUT,
    MODEL_YEARLY_OUT,
    MODEL_SUMMARY_OUT,
    PAIR_REGISTRY_OUT,
    PAIR_YEARLY_OUT,
    PAIR_SUMMARY_OUT,
    COEFFICIENTS_OUT,
    LANE_YEARLY_OUT,
    LANE_SUMMARY_OUT,
    MASTER_RDS_OUT,
    MANIFEST_OUT
  ),
  purpose = c(
    "Step6 실행 사양, 입력 경로, 검증 메타데이터",
    "Step4/Step5 export 파일 import 감사표",
    "데이터셋별 KM benchmark와 연도별 위험도",
    "모든 Step4/Step5 모델의 연도별 harmonized 성능 및 shortfall",
    "모든 Step4/Step5 모델의 요약 성능, fit 지표, tail gap",
    "Step6 pairwise contrast registry",
    "Step6 pairwise contrast yearly source-of-truth",
    "Step6 pairwise contrast summary",
    "Step4+Step5 통합 coefficient catalog",
    "AFT family envelope versus Cox-latency yearly contrasts",
    "AFT family envelope versus Cox-latency summary",
    "Step6 master RDS with imported fits and derived outputs",
    "Step6 output manifest"
  ),
  rows = c(
    nrow(analysis_spec),
    nrow(import_audit),
    nrow(km_benchmark),
    nrow(model_yearly_metrics),
    nrow(model_summary_metrics),
    nrow(pair_registry),
    nrow(pair_yearly),
    nrow(pair_summary),
    nrow(coefficients_combined),
    nrow(lane_yearly),
    nrow(lane_summary),
    NA_integer_,
    NA_integer_
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(manifest, MANIFEST_OUT, row.names = FALSE, na = "")
manifest$file_size_bytes <- unname(file.info(manifest$file_path)$size)
utils::write.csv(manifest, MANIFEST_OUT, row.names = FALSE, na = "")