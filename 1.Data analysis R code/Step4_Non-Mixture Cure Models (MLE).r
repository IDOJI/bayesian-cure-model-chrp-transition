# 🔴 Configure: paths and analysis options ===============================
DATA_PATH <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
EXPORT_PATH <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦4.Step4_Non-cure benchmark models/attachments'

INSTALL_MISSING_PACKAGES <- FALSE
RANDOM_SEED <- 20260301L

ANALYSIS_BRANCHES <- c("merged", "PNU", "SNU")
PREDICTION_HORIZONS_YEARS <- 1:10

CANDIDATE_MODELS <- c(
  "weibull",
  "lnorm",
  "llogis",
  "gengamma",
  "flexsurvspline_k1"
)

FLEXSURVSPLINE_K <- 1L
FLEXSURVSPLINE_SCALE <- "hazard"

BASELINE_COVARIATES <- c("age_model_c", "sex_fact", "site")
MIN_ROWS_PER_ANALYSIS <- 30L
MIN_EVENTS_PER_ANALYSIS <- 5L
TIME_SCALE_DENOMINATOR <- 365.25
EPSILON_TIME <- 1e-08
MIN_GHAT <- 1e-06

# 🟠 Attach: packages and startup checks ===============================
required_pkgs <- c("survival", "flexsurv", "dplyr", "tidyr", "purrr", "timeROC")

load_or_install_packages <- function(pkgs, install_missing = FALSE) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0 && install_missing) {
    install.packages(missing_pkgs)
  }
  still_missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing) > 0) {
    stop(
      "다음 패키지가 설치되어 있지 않습니다: ",
      paste(still_missing, collapse = ", "),
      "\nINSTALL_MISSING_PACKAGES <- TRUE 로 바꾸거나 직접 설치한 뒤 다시 실행하세요."
    )
  }
  invisible(lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

set.seed(RANDOM_SEED)
options(stringsAsFactors = FALSE)
if (!dir.exists(EXPORT_PATH)) dir.create(EXPORT_PATH, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(DATA_PATH)) stop("DATA_PATH에 지정한 파일이 존재하지 않습니다: ", DATA_PATH)

load_or_install_packages(required_pkgs, install_missing = INSTALL_MISSING_PACKAGES)

# 🔴 Define: helper utilities ===============================
## 🟠 Build: reusable functions for step4 ===============================
choose_age_source <- function(df) {
  candidates <- c("age_exact_entry", "age_int")
  picked <- candidates[candidates %in% names(df)]
  if (length(picked) == 0) return(NA_character_)
  picked[1]
}

standardize_status_to_transition <- function(df) {
  if ("status_num" %in% names(df)) {
    status_num_chr <- as.character(df$status_num)
    suppressWarnings(status_num_num <- as.numeric(status_num_chr))
    event_transition <- ifelse(!is.na(status_num_num), as.integer(status_num_num == 1), NA_integer_)
    return(event_transition)
  }
  if ("status" %in% names(df)) {
    status_chr <- tolower(trimws(as.character(df$status)))
    event_transition <- ifelse(status_chr == "transition", 1L,
                               ifelse(status_chr %in% c("right_censoring", "remission"), 0L, NA_integer_))
    return(event_transition)
  }
  stop("status_num 또는 status 컬럼이 필요합니다.")
}

capture_warnings <- function(expr) {
  warnings <- character(0)
  result <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(result = result, warnings = unique(warnings))
}

safe_fit_noncure_model <- function(model_key, surv_formula, dat) {
  fit_attempt <- capture_warnings(
    tryCatch(
      {
        if (identical(model_key, "flexsurvspline_k1")) {
          flexsurv::flexsurvspline(
            formula = surv_formula,
            data = dat,
            k = FLEXSURVSPLINE_K,
            scale = FLEXSURVSPLINE_SCALE
          )
        } else {
          flexsurv::flexsurvreg(
            formula = surv_formula,
            data = dat,
            dist = model_key
          )
        }
      },
      error = function(e) e
    )
  )
  
  if (inherits(fit_attempt$result, "error")) {
    return(list(
      success = FALSE,
      fit = NULL,
      error_message = conditionMessage(fit_attempt$result),
      warning_message = paste(fit_attempt$warnings, collapse = " | ")
    ))
  }
  
  list(
    success = TRUE,
    fit = fit_attempt$result,
    error_message = NA_character_,
    warning_message = if (length(fit_attempt$warnings) == 0) NA_character_ else paste(fit_attempt$warnings, collapse = " | ")
  )
}

extract_fit_metrics <- function(fit_obj) {
  out <- list(
    logLik = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    npar = NA_real_
  )
  ll <- tryCatch(logLik(fit_obj), error = function(e) NULL)
  if (!is.null(ll)) {
    out$logLik <- as.numeric(ll)
    out$npar <- attr(ll, "df")
  }
  out$AIC <- tryCatch(AIC(fit_obj), error = function(e) NA_real_)
  out$BIC <- tryCatch(BIC(fit_obj), error = function(e) NA_real_)
  out
}

extract_one_subject_survival <- function(fit_obj, one_row_df, times) {
  pred_try <- tryCatch(
    predict(
      fit_obj,
      newdata = one_row_df,
      type = "survival",
      times = times,
      ci = FALSE,
      se.fit = FALSE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(pred_try)) {
    pred_vec <- as.numeric(pred_try)
    if (length(pred_vec) == length(times)) {
      return(pmin(pmax(pred_vec, 0), 1))
    }
  }
  
  sum_try <- tryCatch(
    summary(
      fit_obj,
      newdata = one_row_df,
      type = "survival",
      t = times,
      ci = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.data.frame(sum_try) && "est" %in% names(sum_try)) {
    est <- as.numeric(sum_try$est)
    if (length(est) == length(times)) return(pmin(pmax(est, 0), 1))
  }
  
  if (is.list(sum_try) && length(sum_try) >= 1) {
    first_obj <- sum_try[[1]]
    if (is.data.frame(first_obj) && "est" %in% names(first_obj)) {
      est <- as.numeric(first_obj$est)
      if (length(est) == length(times)) return(pmin(pmax(est, 0), 1))
    }
    if (is.list(first_obj) && !is.null(first_obj$est)) {
      est <- as.numeric(first_obj$est)
      if (length(est) == length(times)) return(pmin(pmax(est, 0), 1))
    }
  }
  
  stop("flexsurv 예측값 파싱에 실패했습니다.")
}

predict_survival_matrix <- function(fit_obj, newdata_df, times) {
  pred_try <- tryCatch(
    predict(
      fit_obj,
      newdata = newdata_df,
      type = "survival",
      times = times,
      ci = FALSE,
      se.fit = FALSE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(pred_try)) {
    pred_mat <- as.matrix(pred_try)
    if (nrow(pred_mat) == nrow(newdata_df) && ncol(pred_mat) == length(times)) {
      pred_mat <- pmin(pmax(pred_mat, 0), 1)
      colnames(pred_mat) <- paste0("year_", times)
      rownames(pred_mat) <- newdata_df$site_id
      return(pred_mat)
    }
    if (nrow(pred_mat) == length(times) && ncol(pred_mat) == nrow(newdata_df)) {
      pred_mat <- t(pred_mat)
      pred_mat <- pmin(pmax(pred_mat, 0), 1)
      colnames(pred_mat) <- paste0("year_", times)
      rownames(pred_mat) <- newdata_df$site_id
      return(pred_mat)
    }
  }
  
  out <- matrix(NA_real_, nrow = nrow(newdata_df), ncol = length(times))
  for (i in seq_len(nrow(newdata_df))) {
    out[i, ] <- extract_one_subject_survival(
      fit_obj = fit_obj,
      one_row_df = newdata_df[i, , drop = FALSE],
      times = times
    )
  }
  out <- pmin(pmax(out, 0), 1)
  colnames(out) <- paste0("year_", times)
  rownames(out) <- newdata_df$site_id
  out
}

make_censoring_survival_evaluator <- function(time, event) {
  censor_event <- 1L - as.integer(event)
  if (sum(censor_event, na.rm = TRUE) == 0L) {
    return(function(t, left_limit = FALSE) rep(1, length(t)))
  }
  
  sf_cens <- survival::survfit(survival::Surv(time, censor_event) ~ 1)
  
  function(t, left_limit = FALSE) {
    tt <- as.numeric(t)
    if (left_limit) tt <- pmax(tt - EPSILON_TIME, 0)
    surv_vals <- summary(sf_cens, times = tt, extend = TRUE)$surv
    if (length(surv_vals) == 0L) surv_vals <- rep(1, length(tt))
    surv_vals[is.na(surv_vals)] <- 1
    pmax(as.numeric(surv_vals), MIN_GHAT)
  }
}

compute_ipcw_brier <- function(time, event, surv_pred_t, horizon_t, Gfun) {
  time <- as.numeric(time)
  event <- as.integer(event)
  surv_pred_t <- as.numeric(surv_pred_t)
  
  G_t <- Gfun(horizon_t, left_limit = FALSE)
  G_y_minus <- Gfun(time, left_limit = TRUE)
  
  obs_surv_status <- ifelse(time > horizon_t, 1,
                            ifelse(time <= horizon_t & event == 1L, 0, NA_real_))
  
  weights <- ifelse(time <= horizon_t & event == 1L, 1 / G_y_minus,
                    ifelse(time > horizon_t, 1 / G_t, NA_real_))
  
  brier_val <- mean(weights * (obs_surv_status - surv_pred_t)^2, na.rm = TRUE)
  as.numeric(brier_val)
}

compute_auc_at_time <- function(time, event, risk_pred_t, horizon_t) {
  time <- as.numeric(time)
  event <- as.integer(event)
  risk_pred_t <- as.numeric(risk_pred_t)
  
  n_cases <- sum(time <= horizon_t & event == 1L, na.rm = TRUE)
  n_controls <- sum(time > horizon_t, na.rm = TRUE)
  
  if (n_cases < 1L || n_controls < 1L) return(NA_real_)
  
  if (sd(risk_pred_t, na.rm = TRUE) == 0) return(0.5)
  
  roc_obj <- tryCatch(
    timeROC::timeROC(
      T = time,
      delta = event,
      marker = risk_pred_t,
      cause = 1,
      weighting = "marginal",
      times = horizon_t,
      iid = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(roc_obj)) return(NA_real_)
  auc_val <- roc_obj$AUC
  if (length(auc_val) == 0L) return(NA_real_)
  as.numeric(auc_val[1])
}

trapz_scalar <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) == 0L) return(NA_real_)
  if (length(x) == 1L) return(y[1])
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

build_branch_data <- function(df, branch_name) {
  if (identical(branch_name, "merged")) return(df)
  df[df$site == branch_name, , drop = FALSE]
}

build_model_frame <- function(df_branch) {
  usable_covars <- BASELINE_COVARIATES[BASELINE_COVARIATES %in% names(df_branch)]
  
  if ("site" %in% usable_covars && length(unique(df_branch$site[!is.na(df_branch$site)])) <= 1L) {
    usable_covars <- setdiff(usable_covars, "site")
  }
  
  non_constant <- usable_covars[
    vapply(
      usable_covars,
      function(v) length(unique(df_branch[[v]][!is.na(df_branch[[v]])])) > 1L,
      logical(1)
    )
  ]
  
  keep_vars <- unique(c("site_id", "site", "id", "time_years", "event_transition", non_constant))
  dat_cc <- df_branch[, keep_vars, drop = FALSE]
  cc_flag <- complete.cases(dat_cc)
  dat_cc <- dat_cc[cc_flag, , drop = FALSE]
  
  non_constant_after_cc <- non_constant[
    vapply(
      non_constant,
      function(v) length(unique(dat_cc[[v]][!is.na(dat_cc[[v]])])) > 1L,
      logical(1)
    )
  ]
  
  rhs <- if (length(non_constant_after_cc) == 0L) "1" else paste(non_constant_after_cc, collapse = " + ")
  surv_formula <- as.formula(paste0("survival::Surv(time_years, event_transition) ~ ", rhs))
  
  list(
    data = dat_cc,
    formula = surv_formula,
    used_covariates = non_constant_after_cc,
    excluded_due_to_missing = sum(!cc_flag)
  )
}

# 🔴 Import: merged cohort and prepare transition endpoint ===============================
## 🟠 Read: raw csv and validate schema ===============================
raw_df <- read.csv(
  file = DATA_PATH,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

required_any_of <- list(
  id_site = c("id", "site"),
  followup = c("days_followup"),
  status = c("status_num", "status")
)

missing_core <- c(
  setdiff(required_any_of$id_site, names(raw_df)),
  setdiff(required_any_of$followup, names(raw_df))
)
if (length(intersect(required_any_of$status, names(raw_df))) == 0L) {
  missing_core <- c(missing_core, "status_num/status")
}
if (length(missing_core) > 0L) {
  stop("필수 컬럼이 없습니다: ", paste(unique(missing_core), collapse = ", "))
}

## 🟠 Transform: transition-only analysis dataset ===============================
analysis_df <- raw_df

analysis_df$site <- trimws(as.character(analysis_df$site))
analysis_df$id <- as.character(analysis_df$id)
analysis_df$site_id <- paste0(analysis_df$site, "__", analysis_df$id)

analysis_df$days_followup <- suppressWarnings(as.numeric(as.character(analysis_df$days_followup)))
analysis_df$event_transition <- standardize_status_to_transition(analysis_df)

age_source <- choose_age_source(analysis_df)
if (!is.na(age_source)) {
  analysis_df$age_model <- suppressWarnings(as.numeric(as.character(analysis_df[[age_source]])))
} else {
  analysis_df$age_model <- NA_real_
}
analysis_df$age_model_c <- analysis_df$age_model - mean(analysis_df$age_model, na.rm = TRUE)

if (!("sex_fact" %in% names(analysis_df))) {
  if ("sex_num" %in% names(analysis_df)) {
    sex_num_tmp <- suppressWarnings(as.numeric(as.character(analysis_df$sex_num)))
    analysis_df$sex_fact <- ifelse(sex_num_tmp == 0, "Female",
                                   ifelse(sex_num_tmp == 1, "Male", NA_character_))
  } else {
    analysis_df$sex_fact <- NA_character_
  }
}
analysis_df$sex_fact <- factor(analysis_df$sex_fact, levels = c("Female", "Male"))
analysis_df$site <- factor(analysis_df$site)

analysis_df$time_years <- analysis_df$days_followup / TIME_SCALE_DENOMINATOR

initial_exclusion_df <- data.frame(
  row_index = seq_len(nrow(analysis_df)),
  site_id = analysis_df$site_id,
  exclude_reason = NA_character_,
  stringsAsFactors = FALSE
)

initial_exclusion_df$exclude_reason[is.na(analysis_df$id) | trimws(analysis_df$id) == ""] <- "missing_id"
initial_exclusion_df$exclude_reason[is.na(initial_exclusion_df$exclude_reason) & (is.na(analysis_df$site) | trimws(as.character(analysis_df$site)) == "")] <- "missing_site"
initial_exclusion_df$exclude_reason[is.na(initial_exclusion_df$exclude_reason) & is.na(analysis_df$days_followup)] <- "missing_days_followup"
initial_exclusion_df$exclude_reason[is.na(initial_exclusion_df$exclude_reason) & analysis_df$days_followup < 0] <- "negative_days_followup"
initial_exclusion_df$exclude_reason[is.na(initial_exclusion_df$exclude_reason) & is.na(analysis_df$event_transition)] <- "missing_or_invalid_status"

analysis_df_valid <- analysis_df[is.na(initial_exclusion_df$exclude_reason), , drop = FALSE]
initial_exclusion_df <- initial_exclusion_df[!is.na(initial_exclusion_df$exclude_reason), , drop = FALSE]

# 🔴 Fit: Step4 non-cure benchmark models ===============================
## 🟠 Partition: dataset branches for modeling ===============================
analysis_registry_rows <- list()
fit_message_rows <- list()
pred_subject_rows <- list()
pred_mean_rows <- list()
perf_time_rows <- list()
perf_overall_rows <- list()
fit_nc_list <- list()
fit_nc_best <- list()

for (branch_name in ANALYSIS_BRANCHES) {
  branch_df_raw <- build_branch_data(analysis_df_valid, branch_name)
  analysis_id <- paste(branch_name, "overall", sep = "__")
  
  if (nrow(branch_df_raw) == 0L) {
    fit_message_rows[[length(fit_message_rows) + 1L]] <- data.frame(
      analysis_id = analysis_id,
      dataset_branch = branch_name,
      subgroup = "overall",
      model_name = NA_character_,
      message_type = "skip",
      message_text = "분석 대상 행이 0개라서 건너뜀",
      stringsAsFactors = FALSE
    )
    next
  }
  
  model_frame_info <- build_model_frame(branch_df_raw)
  branch_df <- model_frame_info$data
  
  if (nrow(branch_df) < MIN_ROWS_PER_ANALYSIS) {
    fit_message_rows[[length(fit_message_rows) + 1L]] <- data.frame(
      analysis_id = analysis_id,
      dataset_branch = branch_name,
      subgroup = "overall",
      model_name = NA_character_,
      message_type = "skip",
      message_text = paste0("완전케이스 기준 행 수가 ", nrow(branch_df), "개로 최소 기준 미만이라 건너뜀"),
      stringsAsFactors = FALSE
    )
    next
  }
  
  if (sum(branch_df$event_transition, na.rm = TRUE) < MIN_EVENTS_PER_ANALYSIS) {
    fit_message_rows[[length(fit_message_rows) + 1L]] <- data.frame(
      analysis_id = analysis_id,
      dataset_branch = branch_name,
      subgroup = "overall",
      model_name = NA_character_,
      message_type = "skip",
      message_text = paste0("event 수가 ", sum(branch_df$event_transition, na.rm = TRUE), "개로 최소 기준 미만이라 건너뜀"),
      stringsAsFactors = FALSE
    )
    next
  }
  
  fit_nc_list[[analysis_id]] <- list()
  surv_formula <- model_frame_info$formula
  formula_chr <- paste(deparse(surv_formula), collapse = " ")
  
  Gfun <- make_censoring_survival_evaluator(
    time = branch_df$time_years,
    event = branch_df$event_transition
  )
  
  branch_max_followup <- max(branch_df$time_years, na.rm = TRUE)
  branch_last_event <- if (sum(branch_df$event_transition, na.rm = TRUE) > 0L) {
    max(branch_df$time_years[branch_df$event_transition == 1L], na.rm = TRUE)
  } else {
    NA_real_
  }
  
  for (model_key in CANDIDATE_MODELS) {
    fit_res <- safe_fit_noncure_model(
      model_key = model_key,
      surv_formula = surv_formula,
      dat = branch_df
    )
    
    if (!fit_res$success) {
      analysis_registry_rows[[length(analysis_registry_rows) + 1L]] <- data.frame(
        analysis_id = analysis_id,
        dataset_branch = branch_name,
        subgroup = "overall",
        model_name = model_key,
        formula = formula_chr,
        used_covariates = paste(model_frame_info$used_covariates, collapse = " + "),
        n = nrow(branch_df),
        n_event = sum(branch_df$event_transition, na.rm = TRUE),
        n_censor = sum(branch_df$event_transition == 0L, na.rm = TRUE),
        excluded_due_to_missing_covariates = model_frame_info$excluded_due_to_missing,
        max_followup_years = branch_max_followup,
        last_event_years = branch_last_event,
        logLik = NA_real_,
        AIC = NA_real_,
        BIC = NA_real_,
        npar = NA_real_,
        fit_success = FALSE,
        error_message = fit_res$error_message,
        warning_message = fit_res$warning_message,
        stringsAsFactors = FALSE
      )
      
      fit_message_rows[[length(fit_message_rows) + 1L]] <- data.frame(
        analysis_id = analysis_id,
        dataset_branch = branch_name,
        subgroup = "overall",
        model_name = model_key,
        message_type = "fit_error",
        message_text = fit_res$error_message,
        stringsAsFactors = FALSE
      )
      next
    }
    
    fit_obj <- fit_res$fit
    fit_nc_list[[analysis_id]][[model_key]] <- fit_obj
    
    fit_metrics <- extract_fit_metrics(fit_obj)
    
    surv_mat <- predict_survival_matrix(
      fit_obj = fit_obj,
      newdata_df = branch_df,
      times = PREDICTION_HORIZONS_YEARS
    )
    risk_mat <- 1 - surv_mat
    
    subject_pred_df <- data.frame(
      site_id = rep(branch_df$site_id, each = length(PREDICTION_HORIZONS_YEARS)),
      site = rep(as.character(branch_df$site), each = length(PREDICTION_HORIZONS_YEARS)),
      id = rep(branch_df$id, each = length(PREDICTION_HORIZONS_YEARS)),
      dataset_branch = branch_name,
      subgroup = "overall",
      analysis_id = analysis_id,
      model_name = model_key,
      year = rep(PREDICTION_HORIZONS_YEARS, times = nrow(branch_df)),
      surv_prob = as.numeric(t(surv_mat)),
      risk_prob = as.numeric(t(risk_mat)),
      stringsAsFactors = FALSE
    )
    pred_subject_rows[[length(pred_subject_rows) + 1L]] <- subject_pred_df
    
    mean_pred_df <- data.frame(
      dataset_branch = branch_name,
      subgroup = "overall",
      analysis_id = analysis_id,
      model_name = model_key,
      year = PREDICTION_HORIZONS_YEARS,
      mean_surv_prob = colMeans(surv_mat, na.rm = TRUE),
      mean_risk_prob = colMeans(risk_mat, na.rm = TRUE),
      prediction_source = "mean_subject_prediction",
      stringsAsFactors = FALSE
    )
    pred_mean_rows[[length(pred_mean_rows) + 1L]] <- mean_pred_df
    
    perf_by_time_one_model <- lapply(seq_along(PREDICTION_HORIZONS_YEARS), function(j) {
      horizon_t <- PREDICTION_HORIZONS_YEARS[j]
      risk_t <- risk_mat[, j]
      surv_t <- surv_mat[, j]
      
      n_cases_t <- sum(branch_df$time_years <= horizon_t & branch_df$event_transition == 1L, na.rm = TRUE)
      n_controls_t <- sum(branch_df$time_years > horizon_t, na.rm = TRUE)
      
      auc_t <- if (horizon_t <= branch_max_followup) {
        compute_auc_at_time(
          time = branch_df$time_years,
          event = branch_df$event_transition,
          risk_pred_t = risk_t,
          horizon_t = horizon_t
        )
      } else {
        NA_real_
      }
      
      brier_t <- if (horizon_t <= branch_max_followup) {
        compute_ipcw_brier(
          time = branch_df$time_years,
          event = branch_df$event_transition,
          surv_pred_t = surv_t,
          horizon_t = horizon_t,
          Gfun = Gfun
        )
      } else {
        NA_real_
      }
      
      data.frame(
        dataset_branch = branch_name,
        subgroup = "overall",
        analysis_id = analysis_id,
        model_name = model_key,
        year = horizon_t,
        n_cases_t = n_cases_t,
        n_controls_t = n_controls_t,
        max_followup_years = branch_max_followup,
        auc_t = auc_t,
        brier_t = brier_t,
        metric_in_observed_support = horizon_t <= branch_max_followup,
        stringsAsFactors = FALSE
      )
    })
    perf_by_time_one_model <- dplyr::bind_rows(perf_by_time_one_model)
    perf_time_rows[[length(perf_time_rows) + 1L]] <- perf_by_time_one_model
    
    valid_brier_df <- perf_by_time_one_model[is.finite(perf_by_time_one_model$brier_t), , drop = FALSE]
    ibs_grid <- if (nrow(valid_brier_df) == 0L) {
      NA_real_
    } else if (nrow(valid_brier_df) == 1L) {
      valid_brier_df$brier_t[1]
    } else {
      trapz_scalar(valid_brier_df$year, valid_brier_df$brier_t) /
        (max(valid_brier_df$year) - min(valid_brier_df$year))
    }
    
    perf_overall_rows[[length(perf_overall_rows) + 1L]] <- data.frame(
      analysis_id = analysis_id,
      dataset_branch = branch_name,
      subgroup = "overall",
      model_name = model_key,
      formula = formula_chr,
      used_covariates = paste(model_frame_info$used_covariates, collapse = " + "),
      n = nrow(branch_df),
      n_event = sum(branch_df$event_transition, na.rm = TRUE),
      n_censor = sum(branch_df$event_transition == 0L, na.rm = TRUE),
      excluded_due_to_missing_covariates = model_frame_info$excluded_due_to_missing,
      max_followup_years = branch_max_followup,
      last_event_years = branch_last_event,
      logLik = fit_metrics$logLik,
      AIC = fit_metrics$AIC,
      BIC = fit_metrics$BIC,
      npar = fit_metrics$npar,
      IBS_grid_1to10y = ibs_grid,
      fit_success = TRUE,
      warning_message = fit_res$warning_message,
      stringsAsFactors = FALSE
    )
    
    analysis_registry_rows[[length(analysis_registry_rows) + 1L]] <- data.frame(
      analysis_id = analysis_id,
      dataset_branch = branch_name,
      subgroup = "overall",
      model_name = model_key,
      formula = formula_chr,
      used_covariates = paste(model_frame_info$used_covariates, collapse = " + "),
      n = nrow(branch_df),
      n_event = sum(branch_df$event_transition, na.rm = TRUE),
      n_censor = sum(branch_df$event_transition == 0L, na.rm = TRUE),
      excluded_due_to_missing_covariates = model_frame_info$excluded_due_to_missing,
      max_followup_years = branch_max_followup,
      last_event_years = branch_last_event,
      logLik = fit_metrics$logLik,
      AIC = fit_metrics$AIC,
      BIC = fit_metrics$BIC,
      npar = fit_metrics$npar,
      fit_success = TRUE,
      error_message = NA_character_,
      warning_message = fit_res$warning_message,
      stringsAsFactors = FALSE
    )
  }
}

## 🟠 Select: best models by information criteria ===============================
analysis_registry_df <- dplyr::bind_rows(analysis_registry_rows)
fit_messages_df <- dplyr::bind_rows(fit_message_rows)
pred_subject_nc <- dplyr::bind_rows(pred_subject_rows)
pred_mean_nc <- dplyr::bind_rows(pred_mean_rows)
perf_nc_time <- dplyr::bind_rows(perf_time_rows)
perf_nc_overall <- dplyr::bind_rows(perf_overall_rows)

if (nrow(analysis_registry_df) > 0L) {
  analysis_registry_df <- analysis_registry_df %>%
    dplyr::group_by(analysis_id) %>%
    dplyr::mutate(
      best_by_AIC = fit_success & is.finite(AIC) & (AIC == min(AIC[fit_success & is.finite(AIC)], na.rm = TRUE)),
      best_by_BIC = fit_success & is.finite(BIC) & (BIC == min(BIC[fit_success & is.finite(BIC)], na.rm = TRUE))
    ) %>%
    dplyr::ungroup()
  
  perf_nc_overall <- perf_nc_overall %>%
    dplyr::left_join(
      analysis_registry_df %>%
        dplyr::select(analysis_id, model_name, best_by_AIC, best_by_BIC),
      by = c("analysis_id", "model_name")
    )
  
  successful_best_aic <- analysis_registry_df %>%
    dplyr::filter(fit_success, best_by_AIC) %>%
    dplyr::select(analysis_id, model_name)
  
  if (nrow(successful_best_aic) > 0L) {
    for (i in seq_len(nrow(successful_best_aic))) {
      aid <- successful_best_aic$analysis_id[i]
      mname <- successful_best_aic$model_name[i]
      fit_nc_best[[aid]] <- fit_nc_list[[aid]][[mname]]
    }
  }
}

# 🔴 Export: Step4 artifacts and reusable objects ===============================
## 🟠 Write: csv summaries and rds bundle ===============================
step4_bundle <- list(
  meta = list(
    step = "Step4_non_cure_benchmark_MLE",
    data_path = DATA_PATH,
    export_path = EXPORT_PATH,
    analysis_branches = ANALYSIS_BRANCHES,
    prediction_horizons_years = PREDICTION_HORIZONS_YEARS,
    candidate_models = CANDIDATE_MODELS,
    remission_handling = "remission_as_censoring",
    time_unit = "years_from_days_followup_div_365.25",
    age_source = age_source,
    random_seed = RANDOM_SEED
  ),
  analysis_data_prepared = analysis_df_valid,
  initial_exclusions = initial_exclusion_df,
  fit_nc_list = fit_nc_list,
  fit_nc_best = fit_nc_best,
  analysis_registry = analysis_registry_df,
  pred_subject_nc = pred_subject_nc,
  pred_mean_nc = pred_mean_nc,
  perf_nc_time = perf_nc_time,
  perf_nc_overall = perf_nc_overall,
  fit_messages = fit_messages_df
)

saveRDS(
  object = step4_bundle,
  file = file.path(EXPORT_PATH, "step4_noncure_bundle.rds")
)

utils::write.csv(
  analysis_registry_df,
  file = file.path(EXPORT_PATH, "step4_noncure_fit_registry.csv"),
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  perf_nc_time,
  file = file.path(EXPORT_PATH, "step4_noncure_performance_by_time.csv"),
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  perf_nc_overall,
  file = file.path(EXPORT_PATH, "step4_noncure_performance_overall.csv"),
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  pred_mean_nc,
  file = file.path(EXPORT_PATH, "step4_noncure_mean_predictions_by_time.csv"),
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  pred_subject_nc,
  file = file.path(EXPORT_PATH, "step4_noncure_subject_predictions_long.csv"),
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  initial_exclusion_df,
  file = file.path(EXPORT_PATH, "step4_noncure_initial_exclusions.csv"),
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  fit_messages_df,
  file = file.path(EXPORT_PATH, "step4_noncure_fit_messages.csv"),
  row.names = FALSE,
  na = ""
)