# 🔴 설정: 경로와 Step4 옵션 ===============================
data_path <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
export_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦4.Step4_Non-cure benchmark models/attachments'

file_prefix <- "step4_nocure"
age_var_preferred <- "age_exact_entry"   # 없으면 age_int로 자동 대체
site_reference <- "PNU"
prediction_years <- 1:10
days_per_year <- 365.25
time_zero_correction_days <- 0.5         # survreg용 0일 추적 보정
survreg_maxiter <- 100
cox_ties <- "efron"

options(stringsAsFactors = FALSE, scipen = 999)

if (!dir.exists(export_path)) {
  dir.create(export_path, recursive = TRUE, showWarnings = FALSE)
}

analysis_data_file <- file.path(export_path, paste0(file_prefix, "_analysis_data.csv"))
model_registry_file <- file.path(export_path, paste0(file_prefix, "_model_registry.csv"))
subject_predictions_file <- file.path(export_path, paste0(file_prefix, "_subject_yearly_predictions.csv"))
mean_predictions_file <- file.path(export_path, paste0(file_prefix, "_mean_yearly_predictions.csv"))
performance_yearly_file <- file.path(export_path, paste0(file_prefix, "_performance_yearly.csv"))
performance_summary_file <- file.path(export_path, paste0(file_prefix, "_performance_summary.csv"))
fit_bundle_file <- file.path(export_path, paste0(file_prefix, "_fit_bundle.rds"))
manifest_file <- file.path(export_path, paste0(file_prefix, "_manifest.csv"))

# 🟠 확인: 패키지와 기본 환경 ===============================
required_packages <- c("survival")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("다음 패키지가 필요합니다: ", paste(missing_packages, collapse = ", "))
}
library(survival)

# 🟠 정의: 공통 헬퍼 함수 ===============================
bind_rows_safe <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) == 0) {
    return(data.frame())
  }
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

# 🟡 정규화: 사이트와 기본 변수 ===============================
normalize_site <- function(x) {
  y <- as.character(x)
  y <- trimws(y)
  y <- toupper(y)
  y <- gsub("[[:space:]_\\-]", "", y)
  
  out <- ifelse(y == "PNU", "PNU",
                ifelse(y %in% c("SNU", "CHRP"), "SNU", NA_character_))
  
  uniq <- unique(y[!is.na(y)])
  if (any(is.na(out)) && length(uniq) == 2 && "PNU" %in% uniq) {
    other <- setdiff(uniq, "PNU")
    out[y == other] <- "SNU"
  }
  out
}

choose_age_var <- function(colnames_raw, preferred = "age_exact_entry") {
  candidates <- c(preferred, "age_exact_entry", "age_int")
  candidates <- unique(candidates)
  hit <- candidates[candidates %in% colnames_raw]
  if (length(hit) == 0) {
    stop("age 변수를 찾을 수 없습니다. age_exact_entry 또는 age_int가 필요합니다.")
  }
  hit[1]
}

validate_binary <- function(x, var_name) {
  ok <- unique(x[!is.na(x)])
  if (!all(ok %in% c(0, 1))) {
    stop(var_name, " 는 0/1 binary coding 이어야 합니다.")
  }
  invisible(TRUE)
}

derive_status_num <- function(raw) {
  if ("status_num" %in% names(raw)) {
    return(as.integer(raw$status_num))
  }
  if (!"status" %in% names(raw)) {
    stop("status_num 또는 status 컬럼이 필요합니다.")
  }
  status_chr <- tolower(trimws(as.character(raw$status)))
  out <- ifelse(status_chr == "transition", 1L,
                ifelse(status_chr == "remission", 2L,
                       ifelse(status_chr == "right_censoring", 0L, NA_integer_)))
  as.integer(out)
}

derive_sex_num <- function(raw) {
  if ("sex_num" %in% names(raw)) {
    return(as.integer(raw$sex_num))
  }
  if (!"sex_fact" %in% names(raw)) {
    stop("sex_num 또는 sex_fact 컬럼이 필요합니다.")
  }
  sex_chr <- tolower(trimws(as.character(raw$sex_fact)))
  out <- ifelse(sex_chr == "female", 0L,
                ifelse(sex_chr == "male", 1L, NA_integer_))
  as.integer(out)
}

prepare_base_data <- function(raw,
                              age_var_preferred,
                              days_per_year,
                              time_zero_correction_days,
                              site_reference) {
  age_var_used <- choose_age_var(names(raw), preferred = age_var_preferred)
  
  id_chr <- as.character(raw$id)
  site_std <- normalize_site(raw$site)
  age_raw <- suppressWarnings(as.numeric(raw[[age_var_used]]))
  sex_num <- suppressWarnings(as.integer(derive_sex_num(raw)))
  status_num <- suppressWarnings(as.integer(derive_status_num(raw)))
  days_followup <- suppressWarnings(as.numeric(raw$days_followup))
  
  status_label <- if ("status" %in% names(raw)) {
    as.character(raw$status)
  } else {
    ifelse(status_num == 1L, "transition",
           ifelse(status_num == 2L, "remission", "right_censoring"))
  }
  
  keep <- !is.na(id_chr) &
    nzchar(id_chr) &
    !is.na(site_std) &
    is.finite(age_raw) &
    !is.na(sex_num) &
    !is.na(status_num) &
    is.finite(days_followup)
  
  dat <- data.frame(
    id = id_chr[keep],
    site_raw = as.character(raw$site[keep]),
    site = site_std[keep],
    age_raw = age_raw[keep],
    sex_num = sex_num[keep],
    status_num = status_num[keep],
    status = status_label[keep],
    days_followup = days_followup[keep],
    age_var_used = age_var_used,
    stringsAsFactors = FALSE
  )
  
  if (nrow(dat) == 0) {
    stop("분석 가능한 행이 없습니다.")
  }
  
  if (any(dat$days_followup < 0, na.rm = TRUE)) {
    stop("days_followup 에 음수가 존재합니다.")
  }
  
  validate_binary(dat$sex_num, "sex_num")
  
  if (any(!dat$status_num %in% c(0L, 1L, 2L))) {
    stop("status_num 은 0/1/2 coding 이어야 합니다.")
  }
  
  dat$event_transition <- as.integer(dat$status_num == 1L)
  dat$time_years_raw <- dat$days_followup / days_per_year
  dat$time_was_zero_adjusted <- as.integer(dat$time_years_raw <= 0)
  dat$time_years_model <- ifelse(
    dat$time_years_raw <= 0,
    time_zero_correction_days / days_per_year,
    dat$time_years_raw
  )
  
  all_site_levels <- unique(c(site_reference, sort(unique(dat$site))))
  dat$site <- factor(dat$site, levels = all_site_levels)
  dat$subject_uid <- paste0(as.character(dat$site), "::", dat$id)
  
  if (anyDuplicated(dat$subject_uid) > 0) {
    dup_ids <- unique(dat$subject_uid[duplicated(dat$subject_uid)])
    stop("site + id 기준 중복 행이 존재합니다: ", paste(head(dup_ids, 10), collapse = ", "))
  }
  
  dat
}

# 🟡 준비: 브랜치 분석 데이터 ===============================
prepare_branch_data <- function(base_dat, dataset_id, site_reference) {
  if (dataset_id == "merged") {
    dat <- base_dat
  } else {
    dat <- base_dat[as.character(base_dat$site) == dataset_id, , drop = FALSE]
  }
  
  if (nrow(dat) == 0) {
    stop(dataset_id, " 브랜치 데이터가 비어 있습니다.")
  }
  
  dat <- dat[order(dat$subject_uid), , drop = FALSE]
  
  age_center <- mean(dat$age_raw, na.rm = TRUE)
  age_sd <- stats::sd(dat$age_raw, na.rm = TRUE)
  
  if (!is.finite(age_sd) || age_sd <= 0) {
    stop(dataset_id, " 브랜치에서 age SD가 0 또는 비정상입니다.")
  }
  
  dat$dataset_id <- dataset_id
  dat$subgroup_id <- "overall"
  dat$age_center_branch <- age_center
  dat$age_scale_2sd_branch <- 2 * age_sd
  dat$age_s <- (dat$age_raw - age_center) / (2 * age_sd)
  
  if (dataset_id == "merged") {
    dat$site <- factor(as.character(dat$site))
    if (site_reference %in% levels(dat$site)) {
      dat$site <- stats::relevel(dat$site, ref = site_reference)
    }
  } else {
    dat$site <- droplevels(factor(as.character(dat$site)))
  }
  
  dat
}

build_branch_list <- function(base_dat, site_reference) {
  out <- list(
    PNU = prepare_branch_data(base_dat, "PNU", site_reference = site_reference),
    SNU = prepare_branch_data(base_dat, "SNU", site_reference = site_reference),
    merged = prepare_branch_data(base_dat, "merged", site_reference = site_reference)
  )
  out
}

# 🟡 생성: no-cure 모델 그리드 ===============================
build_nocure_specs <- function(dataset_id) {
  if (dataset_id %in% c("PNU", "SNU")) {
    return(list(
      N0 = list(rhs = "age_s + sex_num", interaction_flag = 0L, site_flag = 0L),
      N1 = list(rhs = "age_s * sex_num", interaction_flag = 1L, site_flag = 0L)
    ))
  }
  
  if (dataset_id == "merged") {
    return(list(
      N0S0 = list(rhs = "age_s + sex_num", interaction_flag = 0L, site_flag = 0L),
      N1S0 = list(rhs = "age_s * sex_num", interaction_flag = 1L, site_flag = 0L),
      N0S1 = list(rhs = "age_s + sex_num + site", interaction_flag = 0L, site_flag = 1L),
      N1S1 = list(rhs = "age_s * sex_num + site", interaction_flag = 1L, site_flag = 1L)
    ))
  }
  
  stop("알 수 없는 dataset_id: ", dataset_id)
}

build_nocure_grid <- function(dataset_id) {
  specs <- build_nocure_specs(dataset_id)
  aft_families <- c("exponential", "weibull", "loglogistic", "lognormal")
  
  rows <- list()
  idx <- 1L
  
  for (spec_name in names(specs)) {
    spec <- specs[[spec_name]]
    
    for (fam in aft_families) {
      model_id <- paste(dataset_id, "AFT_param", fam, spec_name, sep = "__")
      rows[[idx]] <- data.frame(
        dataset_id = dataset_id,
        lane = "AFT_param",
        family = fam,
        cov_spec = spec_name,
        formula_rhs = spec$rhs,
        interaction_flag = spec$interaction_flag,
        site_flag = spec$site_flag,
        model_id = model_id,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
    
    model_id <- paste(dataset_id, "PH_semiparam", "coxph", spec_name, sep = "__")
    rows[[idx]] <- data.frame(
      dataset_id = dataset_id,
      lane = "PH_semiparam",
      family = "coxph",
      cov_spec = spec_name,
      formula_rhs = spec$rhs,
      interaction_flag = spec$interaction_flag,
      site_flag = spec$site_flag,
      model_id = model_id,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
  
  out <- bind_rows_safe(rows)
  rownames(out) <- NULL
  out
}

# 🟡 적합: AFT 및 Cox 모형 ===============================
collect_warnings <- function(expr) {
  warnings <- character(0)
  result <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(result = result, warnings = unique(warnings))
}

fit_survreg_model <- function(formula_obj, data, dist_name, maxiter) {
  fit_res <- collect_warnings(
    survival::survreg(
      formula = formula_obj,
      data = data,
      dist = dist_name,
      control = survival::survreg.control(maxiter = maxiter),
      model = TRUE,
      x = TRUE,
      y = TRUE
    )
  )
  
  if (inherits(fit_res$result, "error")) {
    return(list(fit = NULL, error_message = conditionMessage(fit_res$result), warnings = fit_res$warnings))
  }
  
  list(fit = fit_res$result, error_message = NA_character_, warnings = fit_res$warnings)
}

fit_coxph_model <- function(formula_obj, data, ties) {
  fit_res <- collect_warnings(
    survival::coxph(
      formula = formula_obj,
      data = data,
      ties = ties,
      x = TRUE,
      y = TRUE,
      model = TRUE,
      singular.ok = TRUE
    )
  )
  
  if (inherits(fit_res$result, "error")) {
    return(list(fit = NULL, error_message = conditionMessage(fit_res$result), warnings = fit_res$warnings))
  }
  
  list(fit = fit_res$result, error_message = NA_character_, warnings = fit_res$warnings)
}

safe_vcov_diagnostics <- function(fit, warnings, dist_name = NA_character_) {
  vc <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  singular_flag <- FALSE
  
  if (is.null(vc) || any(!is.finite(vc))) {
    singular_flag <- TRUE
  } else if (length(vc) > 0) {
    eigvals <- tryCatch(
      eigen((vc + t(vc)) / 2, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) NA_real_
    )
    singular_flag <- any(!is.finite(eigvals)) || min(abs(eigvals)) < 1e-12
  }
  
  boundary_flag <- any(grepl("boundary|infinite|monotone likelihood", warnings, ignore.case = TRUE))
  if (!is.na(dist_name) && inherits(fit, "survreg")) {
    if (dist_name != "exponential" && is.finite(fit$scale)) {
      boundary_flag <- boundary_flag || (fit$scale < 1e-06) || (fit$scale > 1e+06)
    }
  }
  
  list(
    singular_hessian_flag = singular_flag,
    boundary_flag = boundary_flag
  )
}

# 🟡 예측: 생존확률과 위험도 ===============================
build_linear_predictor <- function(fit, newdata) {
  term_obj <- if (!is.null(fit$terms)) fit$terms else stats::terms(fit)
  mm <- stats::model.matrix(
    object = stats::delete.response(term_obj),
    data = newdata,
    xlev = fit$xlevels
  )
  
  coef_vec <- stats::coef(fit)
  
  if ("(Intercept)" %in% colnames(mm) && !("(Intercept)" %in% names(coef_vec))) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }
  
  missing_cols <- setdiff(names(coef_vec), colnames(mm))
  if (length(missing_cols) > 0) {
    zero_block <- matrix(
      0,
      nrow = nrow(mm),
      ncol = length(missing_cols),
      dimnames = list(NULL, missing_cols)
    )
    mm <- cbind(mm, zero_block)
  }
  
  mm <- mm[, names(coef_vec), drop = FALSE]
  as.numeric(mm %*% coef_vec)
}

predict_survreg_survival <- function(fit, newdata, times, dist_name) {
  lp <- build_linear_predictor(fit, newdata)
  scale_param <- if (dist_name == "exponential") 1 else fit$scale
  
  out <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(times))
  colnames(out) <- paste0("year_", times)
  
  for (j in seq_along(times)) {
    tt <- pmax(times[j], .Machine$double.eps)
    
    if (dist_name %in% c("weibull", "exponential")) {
      out[, j] <- exp(- (tt / exp(lp))^(1 / scale_param))
    } else if (dist_name == "lognormal") {
      out[, j] <- 1 - stats::pnorm((log(tt) - lp) / scale_param)
    } else if (dist_name == "loglogistic") {
      out[, j] <- 1 / (1 + exp((log(tt) - lp) / scale_param))
    } else {
      stop("지원하지 않는 survreg 분포입니다: ", dist_name)
    }
  }
  
  pmin(pmax(out, 0), 1)
}

prepare_basehaz_fun <- function(cox_fit) {
  bh <- survival::basehaz(cox_fit, centered = FALSE)
  if (nrow(bh) == 0) {
    return(function(t) rep(0, length(t)))
  }
  
  bh <- bh[order(bh$time), c("time", "hazard")]
  bh <- bh[!duplicated(bh$time, fromLast = TRUE), , drop = FALSE]
  
  function(t) {
    idx <- findInterval(t, bh$time)
    ifelse(idx <= 0, 0, bh$hazard[idx])
  }
}

predict_cox_survival <- function(fit, newdata, times) {
  lp <- build_linear_predictor(fit, newdata)
  H0_fun <- prepare_basehaz_fun(fit)
  H0_vals <- H0_fun(times)
  out <- exp(-outer(exp(lp), H0_vals, `*`))
  colnames(out) <- paste0("year_", times)
  pmin(pmax(out, 0), 1)
}

# 🟡 계산: IPCW 성능 지표 ===============================
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
  if (length(x) < 2 || length(y) < 2) {
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
  ibs_available <- if (sum(brier_ok) >= 2) trapz_mean(yearly_df$year[brier_ok], yearly_df$brier[brier_ok]) else NA_real_
  
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

matrix_to_long <- function(surv_mat, data, times, model_meta) {
  risk_mat <- 1 - surv_mat
  n_subject <- nrow(data)
  n_time <- length(times)
  
  out <- data.frame(
    dataset_id = model_meta$dataset_id,
    subgroup_id = data$subgroup_id,
    model_id = model_meta$model_id,
    lane = model_meta$lane,
    family = model_meta$family,
    cov_spec = model_meta$cov_spec,
    subject_uid = rep(data$subject_uid, each = n_time),
    site = rep(as.character(data$site), each = n_time),
    sex_num = rep(data$sex_num, each = n_time),
    age_raw = rep(data$age_raw, each = n_time),
    age_s = rep(data$age_s, each = n_time),
    year = rep(times, times = n_subject),
    predicted_surv = as.vector(t(surv_mat)),
    predicted_risk = as.vector(t(risk_mat)),
    is_beyond_max_followup = rep(times > max(data$time_years_model), times = n_subject),
    is_beyond_last_event = rep(times > max(data$time_years_model[data$event_transition == 1]), times = n_subject),
    stringsAsFactors = FALSE
  )
  
  rownames(out) <- NULL
  out
}

cohort_mean_from_long <- function(pred_long) {
  if (nrow(pred_long) == 0) {
    return(data.frame())
  }
  
  grp_cols <- c("dataset_id", "subgroup_id", "model_id", "lane", "family", "cov_spec", "year")
  
  surv_mean <- stats::aggregate(
    predicted_surv ~ dataset_id + subgroup_id + model_id + lane + family + cov_spec + year,
    data = pred_long,
    FUN = mean
  )
  
  risk_mean <- stats::aggregate(
    predicted_risk ~ dataset_id + subgroup_id + model_id + lane + family + cov_spec + year,
    data = pred_long,
    FUN = mean
  )
  
  out <- merge(surv_mean, risk_mean, by = grp_cols, all = TRUE)
  names(out)[names(out) == "predicted_surv"] <- "mean_surv"
  names(out)[names(out) == "predicted_risk"] <- "mean_risk"
  out
}

# 🟠 읽기: 원시 데이터와 브랜치 객체 ===============================
raw_dat <- utils::read.csv(data_path, check.names = FALSE)
base_dat <- prepare_base_data(
  raw = raw_dat,
  age_var_preferred = age_var_preferred,
  days_per_year = days_per_year,
  time_zero_correction_days = time_zero_correction_days,
  site_reference = site_reference
)

if (!all(c("PNU", "SNU") %in% unique(as.character(base_dat$site)))) {
  stop("site 표준화 후 PNU와 SNU 두 브랜치를 모두 찾지 못했습니다.")
}

branch_list <- build_branch_list(base_dat, site_reference = site_reference)

analysis_data_export <- bind_rows_safe(lapply(names(branch_list), function(dataset_id) {
  dat <- branch_list[[dataset_id]]
  data.frame(
    dataset_id = dataset_id,
    subgroup_id = dat$subgroup_id,
    subject_uid = dat$subject_uid,
    id = dat$id,
    site = as.character(dat$site),
    site_raw = dat$site_raw,
    age_var_used = dat$age_var_used,
    age_raw = dat$age_raw,
    age_center_branch = dat$age_center_branch,
    age_scale_2sd_branch = dat$age_scale_2sd_branch,
    age_s = dat$age_s,
    sex_num = dat$sex_num,
    status_num = dat$status_num,
    status = dat$status,
    event_transition = dat$event_transition,
    days_followup = dat$days_followup,
    time_years_raw = dat$time_years_raw,
    time_years_model = dat$time_years_model,
    time_was_zero_adjusted = dat$time_was_zero_adjusted,
    stringsAsFactors = FALSE
  )
}))

# 🟠 실행: 전체 no-cure 그리드 적합 ===============================
all_fit_objects <- list()
registry_rows <- list()
prediction_long_rows <- list()
performance_yearly_rows <- list()
performance_summary_rows <- list()
model_grid_rows <- list()

for (dataset_id in names(branch_list)) {
  dat <- branch_list[[dataset_id]]
  censoring_fit <- fit_censoring_km(dat$time_years_model, dat$event_transition)
  grid_df <- build_nocure_grid(dataset_id)
  model_grid_rows[[dataset_id]] <- grid_df
  
  branch_n <- nrow(dat)
  branch_event_n <- sum(dat$event_transition == 1L)
  branch_censored_n <- sum(dat$event_transition == 0L)
  branch_max_followup <- max(dat$time_years_model)
  branch_max_event <- if (branch_event_n > 0) max(dat$time_years_model[dat$event_transition == 1L]) else NA_real_
  branch_site_levels <- paste(levels(dat$site), collapse = "|")
  
  for (i in seq_len(nrow(grid_df))) {
    spec_row <- grid_df[i, , drop = FALSE]
    model_id <- spec_row$model_id
    formula_text <- paste0("survival::Surv(time_years_model, event_transition) ~ ", spec_row$formula_rhs)
    formula_obj <- stats::as.formula(formula_text)
    
    if (spec_row$family != "coxph") {
      fit_res <- fit_survreg_model(
        formula_obj = formula_obj,
        data = dat,
        dist_name = spec_row$family,
        maxiter = survreg_maxiter
      )
    } else {
      fit_res <- fit_coxph_model(
        formula_obj = formula_obj,
        data = dat,
        ties = cox_ties
      )
    }
    
    fit_obj <- fit_res$fit
    error_message <- fit_res$error_message
    warning_messages <- fit_res$warnings
    
    convergence_warning_flag <- any(grepl(
      "did not converge|didn't converge|ran out of iterations|loglik converged before variable|inner loop failed",
      warning_messages,
      ignore.case = TRUE
    ))
    
    fit_success <- !is.null(fit_obj) && is.na(error_message)
    convergence_flag <- fit_success && !convergence_warning_flag
    
    if (fit_success) {
      diag_info <- safe_vcov_diagnostics(
        fit = fit_obj,
        warnings = warning_messages,
        dist_name = if (spec_row$family == "coxph") NA_character_ else spec_row$family
      )
    } else {
      diag_info <- list(
        singular_hessian_flag = NA,
        boundary_flag = NA
      )
    }
    
    if (fit_success && spec_row$family != "coxph") {
      loglik_val <- as.numeric(stats::logLik(fit_obj))
      n_params <- attr(stats::logLik(fit_obj), "df")
      aic_val <- -2 * loglik_val + 2 * n_params
      bic_val <- -2 * loglik_val + log(branch_n) * n_params
      partial_loglik_val <- NA_real_
      aic_partial_val <- NA_real_
      likelihood_type <- "full"
    } else if (fit_success && spec_row$family == "coxph") {
      n_params <- sum(is.finite(stats::coef(fit_obj)))
      partial_loglik_val <- fit_obj$loglik[2]
      aic_partial_val <- -2 * partial_loglik_val + 2 * n_params
      loglik_val <- NA_real_
      aic_val <- NA_real_
      bic_val <- NA_real_
      likelihood_type <- "partial"
    } else {
      n_params <- NA_real_
      loglik_val <- NA_real_
      aic_val <- NA_real_
      bic_val <- NA_real_
      partial_loglik_val <- NA_real_
      aic_partial_val <- NA_real_
      likelihood_type <- if (spec_row$family == "coxph") "partial" else "full"
    }
    
    registry_row <- data.frame(
      dataset_id = dataset_id,
      subgroup_id = "overall",
      model_id = model_id,
      model_class = "nocure",
      lane = spec_row$lane,
      family = spec_row$family,
      cov_spec = spec_row$cov_spec,
      formula_full = formula_text,
      formula_rhs = spec_row$formula_rhs,
      interaction_flag = spec_row$interaction_flag,
      site_flag = spec_row$site_flag,
      likelihood_type = likelihood_type,
      fit_success = fit_success,
      convergence_flag = convergence_flag,
      singular_hessian_flag = diag_info$singular_hessian_flag,
      boundary_flag = diag_info$boundary_flag,
      n_subject = branch_n,
      n_event = branch_event_n,
      n_censored = branch_censored_n,
      max_followup_years = branch_max_followup,
      max_event_years = branch_max_event,
      zero_time_adjusted_n = sum(dat$time_was_zero_adjusted),
      age_var_used = unique(dat$age_var_used),
      age_center_branch = unique(dat$age_center_branch),
      age_scale_2sd_branch = unique(dat$age_scale_2sd_branch),
      site_levels = branch_site_levels,
      n_params = n_params,
      logLik_full = loglik_val,
      AIC = aic_val,
      BIC = bic_val,
      partial_logLik = partial_loglik_val,
      AIC_partial = aic_partial_val,
      warning_messages = if (length(warning_messages) > 0) paste(unique(warning_messages), collapse = " | ") else NA_character_,
      error_message = error_message,
      stringsAsFactors = FALSE
    )
    
    registry_rows[[length(registry_rows) + 1L]] <- registry_row
    
    all_fit_objects[[model_id]] <- list(
      meta = registry_row,
      fit = fit_obj,
      warnings = warning_messages,
      error_message = error_message
    )
    
    if (!fit_success) {
      next
    }
    
    if (spec_row$family != "coxph") {
      surv_mat <- predict_survreg_survival(
        fit = fit_obj,
        newdata = dat,
        times = prediction_years,
        dist_name = spec_row$family
      )
    } else {
      surv_mat <- predict_cox_survival(
        fit = fit_obj,
        newdata = dat,
        times = prediction_years
      )
    }
    
    risk_mat <- 1 - surv_mat
    
    model_meta <- list(
      dataset_id = dataset_id,
      model_id = model_id,
      lane = spec_row$lane,
      family = spec_row$family,
      cov_spec = spec_row$cov_spec
    )
    
    pred_long <- matrix_to_long(
      surv_mat = surv_mat,
      data = dat,
      times = prediction_years,
      model_meta = model_meta
    )
    
    prediction_long_rows[[length(prediction_long_rows) + 1L]] <- pred_long
    
    perf_list <- evaluate_model_performance(
      time = dat$time_years_model,
      event = dat$event_transition,
      surv_mat = surv_mat,
      risk_mat = risk_mat,
      times = prediction_years,
      Gfit = censoring_fit
    )
    
    perf_yearly <- cbind(
      data.frame(
        dataset_id = dataset_id,
        subgroup_id = "overall",
        model_id = model_id,
        lane = spec_row$lane,
        family = spec_row$family,
        cov_spec = spec_row$cov_spec,
        performance_set = "apparent",
        stringsAsFactors = FALSE
      )[rep(1, nrow(perf_list$yearly)), , drop = FALSE],
      perf_list$yearly
    )
    
    perf_summary <- cbind(
      data.frame(
        dataset_id = dataset_id,
        subgroup_id = "overall",
        model_id = model_id,
        lane = spec_row$lane,
        family = spec_row$family,
        cov_spec = spec_row$cov_spec,
        performance_set = "apparent",
        stringsAsFactors = FALSE
      ),
      perf_list$summary
    )
    
    performance_yearly_rows[[length(performance_yearly_rows) + 1L]] <- perf_yearly
    performance_summary_rows[[length(performance_summary_rows) + 1L]] <- perf_summary
  }
}

model_registry_df <- bind_rows_safe(registry_rows)
subject_predictions_df <- bind_rows_safe(prediction_long_rows)
mean_predictions_df <- cohort_mean_from_long(subject_predictions_df)
performance_yearly_df <- bind_rows_safe(performance_yearly_rows)
performance_summary_df <- bind_rows_safe(performance_summary_rows)
model_grid_df <- bind_rows_safe(model_grid_rows)

# 🟠 정리: 결과 객체와 요약 테이블 ===============================
fit_bundle <- list(
  script_stage = "Step4_no_cure_full_grid",
  config = list(
    data_path = data_path,
    export_path = export_path,
    file_prefix = file_prefix,
    age_var_preferred = age_var_preferred,
    age_var_resolved = unique(base_dat$age_var_used),
    site_reference = site_reference,
    prediction_years = prediction_years,
    days_per_year = days_per_year,
    time_zero_correction_days = time_zero_correction_days,
    survreg_maxiter = survreg_maxiter,
    cox_ties = cox_ties
  ),
  analysis_data = analysis_data_export,
  model_grid = model_grid_df,
  model_registry = model_registry_df,
  fit_objects = all_fit_objects,
  subject_yearly_predictions = subject_predictions_df,
  mean_yearly_predictions = mean_predictions_df,
  performance_yearly = performance_yearly_df,
  performance_summary = performance_summary_df,
  session_info = paste(capture.output(sessionInfo()), collapse = "\n"),
  created_at = as.character(Sys.time())
)

# 🟠 저장: CSV와 RDS 산출물 ===============================
utils::write.csv(analysis_data_export, analysis_data_file, row.names = FALSE, na = "")
utils::write.csv(model_registry_df, model_registry_file, row.names = FALSE, na = "")
utils::write.csv(subject_predictions_df, subject_predictions_file, row.names = FALSE, na = "")
utils::write.csv(mean_predictions_df, mean_predictions_file, row.names = FALSE, na = "")
utils::write.csv(performance_yearly_df, performance_yearly_file, row.names = FALSE, na = "")
utils::write.csv(performance_summary_df, performance_summary_file, row.names = FALSE, na = "")
saveRDS(fit_bundle, fit_bundle_file)

manifest_df <- data.frame(
  file_name = c(
    basename(analysis_data_file),
    basename(model_registry_file),
    basename(subject_predictions_file),
    basename(mean_predictions_file),
    basename(performance_yearly_file),
    basename(performance_summary_file),
    basename(fit_bundle_file),
    basename(manifest_file)
  ),
  file_path = c(
    analysis_data_file,
    model_registry_file,
    subject_predictions_file,
    mean_predictions_file,
    performance_yearly_file,
    performance_summary_file,
    fit_bundle_file,
    manifest_file
  ),
  description = c(
    "브랜치별(PNU/SNU/merged) Step4 분석 데이터와 전처리 결과",
    "모형별 적합 상태, formula, fit 통계량, 경고/오류 요약",
    "모형별/개체별/연도별 예측 생존확률 및 위험도",
    "모형별/연도별 cohort-average mean survival/risk",
    "모형별/연도별 cumulative-dynamic AUC 및 Brier score",
    "모형별 IBS(가능 범위), 평균 AUC/Brier 요약",
    "적합 객체와 모든 핵심 산출물을 담은 RDS 번들",
    "export된 파일 목록"
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(manifest_df, manifest_file, row.names = FALSE, na = "")
manifest_df$file_size_bytes <- unname(file.info(manifest_df$file_path)$size)
utils::write.csv(manifest_df, manifest_file, row.names = FALSE, na = "")