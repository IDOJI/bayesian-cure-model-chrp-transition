# rm(list=ls())
# 🔴 Configure: runtime options and package checks ===============================

## 🟠 Declare: user-editable paths and control flags ===============================
options(stringsAsFactors = FALSE, scipen = 999)

data_file_full <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
export_path <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦5.Step5_MLE MCM/attachments"


gate_pass <- c(
  merged = TRUE,
  pnu = TRUE,
  snu = TRUE
)

k_years <- c(1, 3, 5, 10)
days_per_year <- 365.25

zero_time_offset_days <- 1e-03
fitplot_grid_n <- 250L
extrapolation_grid_n <- 350L

cox_ties <- "efron"
cure_optim_method <- "BFGS"
optim_maxit <- 1000L
survreg_maxiter <- 100L

write_fitplot_png <- TRUE
write_extrapolation_png <- TRUE
random_seed <- 20250312L

set.seed(random_seed)

if (!dir.exists(export_path)) {
  dir.create(export_path, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(export_path)) {
  stop("export_path could not be created: ", export_path, call. = FALSE)
}

if (!setequal(names(gate_pass), c("merged", "pnu", "snu"))) {
  stop("gate_pass must be a named logical vector with names: merged, pnu, snu", call. = FALSE)
}


if (!file.exists(data_file_full)) {
  stop("Input CSV not found: ", data_file_full, call. = FALSE)
}

if (length(k_years) == 0L || all(!is.finite(k_years))) {
  stop("k_years must contain at least one finite value.", call. = FALSE)
}

extrapolation_horizon_years <- max(k_years[is.finite(k_years)], na.rm = TRUE)
extrapolation_horizon_days <- extrapolation_horizon_years * days_per_year

## 🟠 Validate: required package namespaces ===============================
required_pkgs <- c("survival", "flexsurvcure")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0L) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(missing_pkgs, collapse = ", "),
      ". Install them before running this script."
    ),
    call. = FALSE
  )
}

# 🔴 Define: reusable helper utilities ===============================

## 🟠 Build: coercion and QC helpers ===============================
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

blank_to_na <- function(x) {
  if (length(x) == 0L || all(is.na(x))) {
    return(NA_character_)
  }
  y <- paste(x, collapse = " | ")
  y <- trimws(y)
  if (!nzchar(y)) NA_character_ else y
}

first_or_na <- function(x) {
  if (length(x) == 0L) NA_character_ else as.character(x[1L])
}

is_missing_string <- function(x) {
  is.na(x) | trimws(as.character(x)) == ""
}

collapse_head <- function(x, n = 10L) {
  x <- utils::head(x, n)
  paste(x, collapse = ", ")
}

coerce_numeric_clean <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", trimws(as.character(x)))))
}

normalize_site <- function(x) {
  toupper(trimws(as.character(x)))
}

coerce_status_num <- function(x) {
  y_raw <- trimws(tolower(as.character(x)))
  out <- rep.int(NA_integer_, length(y_raw))
  
  out[y_raw %in% c("0", "0.0", "right_censoring", "right censoring", "censoring", "censored")] <- 0L
  out[y_raw %in% c("1", "1.0", "transition", "event")] <- 1L
  out[y_raw %in% c("2", "2.0", "remission")] <- 2L
  
  suppressWarnings(y_num <- as.integer(as.numeric(y_raw)))
  idx <- is.na(out) & !is.na(y_num) & y_num %in% c(0L, 1L, 2L)
  out[idx] <- y_num[idx]
  
  out
}

coerce_binary_sex <- function(x) {
  y_raw <- trimws(tolower(as.character(x)))
  out <- rep.int(NA_real_, length(y_raw))
  
  out[y_raw %in% c("0", "0.0", "female", "f")] <- 0
  out[y_raw %in% c("1", "1.0", "male", "m")] <- 1
  
  suppressWarnings(y_num <- as.numeric(y_raw))
  idx <- is.na(out) & !is.na(y_num) & y_num %in% c(0, 1)
  out[idx] <- y_num[idx]
  
  out
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

step_interp <- function(x, y, xout) {
  xout <- as.numeric(xout)
  if (length(xout) == 0L) {
    return(numeric(0))
  }
  
  ok <- is.finite(x) & is.finite(y)
  x <- as.numeric(x[ok])
  y <- as.numeric(y[ok])
  
  if (length(x) == 0L) {
    return(rep(NA_real_, length(xout)))
  }
  
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  
  keep <- !duplicated(x, fromLast = TRUE)
  x <- x[keep]
  y <- y[keep]
  
  if (length(x) == 1L) {
    return(rep(y[1L], length(xout)))
  }
  
  stats::approx(
    x = x,
    y = y,
    xout = xout,
    method = "constant",
    f = 0,
    rule = 2
  )$y
}

make_regular_times <- function(max_time, n_points) {
  if (!is.finite(max_time) || max_time <= 0) {
    return(c(0))
  }
  unique(as.numeric(c(0, seq(0, max_time, length.out = max(2L, n_points)), max_time)))
}

compute_scaling_reference <- function(dat_merged_complete) {
  age_mean <- mean(dat_merged_complete$age_exact_entry, na.rm = TRUE)
  age_sd <- stats::sd(dat_merged_complete$age_exact_entry, na.rm = TRUE)
  sex_mean <- mean(dat_merged_complete$sex_num, na.rm = TRUE)
  
  if (!is.finite(age_mean) || !is.finite(age_sd) || age_sd <= 0) {
    stop("Merged analysis set produced a non-finite or non-positive age SD.", call. = FALSE)
  }
  if (!is.finite(sex_mean)) {
    stop("Merged analysis set produced a non-finite sex mean.", call. = FALSE)
  }
  
  list(
    age_mean = age_mean,
    age_sd = age_sd,
    age_two_sd = 2 * age_sd,
    sex_mean = sex_mean
  )
}

apply_step5_scaling <- function(dat, scaling_ref, zero_time_offset_days) {
  dat$z_age <- (dat$age_exact_entry - scaling_ref$age_mean) / scaling_ref$age_two_sd
  dat$c_sex <- dat$sex_num - scaling_ref$sex_mean
  dat$int_sex_age <- dat$c_sex * dat$z_age
  dat$time_model <- ifelse(dat$time_primary <= 0, zero_time_offset_days, dat$time_primary)
  dat$time_model_offset_applied <- as.integer(dat$time_primary <= 0)
  dat
}

make_rhs_text <- function(covariate_structure) {
  if (identical(covariate_structure, "main")) {
    "c_sex + z_age"
  } else if (identical(covariate_structure, "interaction")) {
    "c_sex + z_age + int_sex_age"
  } else {
    stop("Unknown covariate_structure: ", covariate_structure, call. = FALSE)
  }
}

make_formula_texts <- function(spec) {
  rhs_text <- make_rhs_text(spec$covariate_structure[1L])
  list(
    rhs = rhs_text,
    surv_formula = paste0("survival::Surv(time_model, event_primary) ~ ", rhs_text),
    incidence_formula = if (identical(spec$model_family[1L], "mixture_cure")) paste0("theta(cure_fraction) ~ ", rhs_text) else NA_character_,
    latency_formula = if (identical(spec$model_family[1L], "mixture_cure")) paste0(spec$latency_param[1L], " ~ ", rhs_text) else NA_character_
  )
}

validate_core_inputs <- function(dat) {
  issue_list <- list(
    missing_id = which(is_missing_string(dat$id)),
    missing_site = which(is_missing_string(dat$site)),
    missing_date_entry = which(is_missing_string(dat$date_entry)),
    missing_days_followup = which(!is.finite(dat$days_followup)),
    missing_status_num = which(is.na(dat$status_num)),
    invalid_site = which(!is_missing_string(dat$site) & !(dat$site %in% c("PNU", "SNU"))),
    invalid_days_followup = which(is.finite(dat$days_followup) & dat$days_followup < 0),
    invalid_status_num = which(!is.na(dat$status_num) & !(dat$status_num %in% c(0L, 1L, 2L))),
    invalid_age_exact_entry = which(!dat$age_exact_entry_missing & !is.finite(dat$age_exact_entry)),
    negative_age_exact_entry = which(is.finite(dat$age_exact_entry) & dat$age_exact_entry < 0),
    invalid_sex_num = which(!dat$sex_num_missing & is.na(dat$sex_num)),
    duplicate_subject_key = which(duplicated(dat$subject_key))
  )
  
  bad_names <- names(issue_list)[vapply(issue_list, length, integer(1)) > 0L]
  
  if (length(bad_names) > 0L) {
    msg_parts <- vapply(
      bad_names,
      function(nm) {
        paste0(nm, " [n=", length(issue_list[[nm]]), "] rows=", collapse_head(issue_list[[nm]]))
      },
      character(1)
    )
    
    stop(
      paste(c("Core input validation failed.", msg_parts), collapse = "\n"),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

empty_curve_df <- function() {
  data.frame(
    dataset = character(),
    plot_type = character(),
    curve_role = character(),
    model_id = character(),
    model_family = character(),
    selection_pool = character(),
    latency_dist = character(),
    covariate_structure = character(),
    time_days = numeric(),
    time_years = numeric(),
    survival = numeric(),
    max_event_time_days = numeric(),
    max_followup_days = numeric(),
    support_limit_type = character(),
    is_beyond_max_event = logical(),
    is_beyond_max_followup = logical(),
    is_extrapolated = logical(),
    prediction_regime = character(),
    stringsAsFactors = FALSE
  )
}

empty_coeff_df <- function() {
  data.frame(
    dataset = character(),
    model_id = character(),
    model_family = character(),
    component = character(),
    term = character(),
    estimate = numeric(),
    lower95 = numeric(),
    upper95 = numeric(),
    exp_estimate = numeric(),
    exp_lower95 = numeric(),
    exp_upper95 = numeric(),
    se = numeric(),
    statistic = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_kyear_df <- function() {
  data.frame(
    dataset = character(),
    model_id = character(),
    model_family = character(),
    selection_pool = character(),
    latency_dist = character(),
    covariate_structure = character(),
    k_year = numeric(),
    time_days = numeric(),
    time_years = numeric(),
    survival = numeric(),
    max_event_time_days = numeric(),
    max_followup_days = numeric(),
    support_limit_type = character(),
    is_beyond_max_event = logical(),
    is_beyond_max_followup = logical(),
    is_extrapolated = logical(),
    prediction_regime = character(),
    prediction_note = character(),
    stringsAsFactors = FALSE
  )
}

sanitize_name_key <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", x))
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
        "cannot be evaluated"
      ),
      collapse = "|"
    ),
    x = text,
    ignore.case = TRUE
  )
}

# 🔴 Define: model fitting and extraction helpers ===============================

## 🟠 Build: specification catalog and fitting wrappers ===============================
make_model_specs <- function() {
  data.frame(
    model_id = c(
      "cox_main",
      "cox_interaction",
      "cure_weibull_main",
      "cure_weibull_interaction",
      "cure_lnorm_main",
      "cure_lnorm_interaction",
      "cure_llogis_main",
      "cure_llogis_interaction"
    ),
    model_family = c("cox", "cox", rep("mixture_cure", 6L)),
    selection_pool = c("benchmark", "benchmark", rep("cure_candidate", 6L)),
    latency_dist = c(NA_character_, NA_character_, "weibull", "weibull", "lnorm", "lnorm", "llogis", "llogis"),
    latency_param = c(NA_character_, NA_character_, "scale", "scale", "meanlog", "meanlog", "scale", "scale"),
    covariate_structure = c("main", "interaction", "main", "interaction", "main", "interaction", "main", "interaction"),
    link = c(NA_character_, NA_character_, rep("logistic", 6L)),
    incidence_target = c(NA_character_, NA_character_, rep("cure_fraction_theta", 6L)),
    incidence_note = c(
      NA_character_,
      NA_character_,
      rep("flexsurvcure parameterizes theta as cure fraction; susceptible proportion = 1 - theta", 6L)
    ),
    stringsAsFactors = FALSE
  )
}

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
      converged = FALSE,
      error = conditionMessage(fit_obj),
      warnings = blank_to_na(unique(warnings_caught)),
      elapsed_sec = elapsed_sec
    ))
  }
  
  list(
    fit = fit_obj,
    ok = TRUE,
    converged = NA,
    error = NA_character_,
    warnings = blank_to_na(unique(warnings_caught)),
    elapsed_sec = elapsed_sec
  )
}

get_fit_converged <- function(fit, model_family) {
  if (is.null(fit)) {
    return(FALSE)
  }
  
  coef_ok <- tryCatch({
    cf <- stats::coef(fit)
    all(is.finite(as.numeric(cf)))
  }, error = function(e) TRUE)
  
  if (identical(model_family, "cox")) {
    return(is.null(fit$fail) && coef_ok)
  }
  
  conv_flag <- if (!is.null(fit$converged)) isTRUE(as.logical(fit$converged)) else TRUE
  conv_flag && coef_ok
}

fit_model_from_spec <- function(spec, dat, cox_ties, cure_optim_method, optim_maxit, survreg_maxiter) {
  rhs_text <- make_rhs_text(spec$covariate_structure[1L])
  surv_formula <- stats::as.formula(paste("survival::Surv(time_model, event_primary) ~", rhs_text))
  
  if (identical(spec$model_family[1L], "cox")) {
    fit_info <- capture_fit(function() {
      survival::coxph(
        formula = surv_formula,
        data = dat,
        ties = cox_ties,
        x = TRUE,
        y = TRUE,
        model = TRUE
      )
    })
  } else {
    anc_formula <- stats::as.formula(paste("~", rhs_text))
    anc_list <- setNames(list(anc_formula), spec$latency_param[1L])
    
    fit_info <- capture_fit(function() {
      flexsurvcure::flexsurvcure(
        formula = surv_formula,
        data = dat,
        dist = spec$latency_dist[1L],
        anc = anc_list,
        link = spec$link[1L],
        mixture = TRUE,
        na.action = stats::na.omit,
        sr.control = survival::survreg.control(maxiter = survreg_maxiter),
        method = cure_optim_method,
        control = list(maxit = optim_maxit)
      )
    })
  }
  
  fit_info$converged <- if (!is.null(fit_info$fit)) get_fit_converged(fit_info$fit, spec$model_family[1L]) else FALSE
  fit_info
}

extract_fit_summary <- function(spec, fit_info, dataset_name, dataset_stats) {
  formula_texts <- make_formula_texts(spec)
  
  fit <- fit_info$fit
  fit_ok <- isTRUE(fit_info$ok) && !is.null(fit)
  converged <- isTRUE(fit_info$converged)
  
  loglik_obj <- if (fit_ok) tryCatch(stats::logLik(fit), error = function(e) NULL) else NULL
  loglik_val <- if (!is.null(loglik_obj)) as.numeric(loglik_obj) else NA_real_
  df_model <- if (!is.null(loglik_obj)) as.numeric(attr(loglik_obj, "df")) else NA_real_
  
  if (!is.finite(df_model) && fit_ok) {
    df_model <- tryCatch(length(stats::coef(fit)), error = function(e) NA_real_)
  }
  
  aic_val <- if (fit_ok) tryCatch(as.numeric(stats::AIC(fit)), error = function(e) NA_real_) else NA_real_
  bic_val <- if (fit_ok) tryCatch(as.numeric(stats::BIC(fit)), error = function(e) NA_real_) else NA_real_
  
  if (!is.finite(aic_val) && is.finite(loglik_val) && is.finite(df_model)) {
    aic_val <- -2 * loglik_val + 2 * df_model
  }
  
  if (!is.finite(bic_val) && is.finite(loglik_val) && is.finite(df_model) && is.finite(dataset_stats$n_fit)) {
    bic_val <- -2 * loglik_val + log(dataset_stats$n_fit) * df_model
  }
  
  warning_message <- blank_to_na(fit_info$warnings)
  error_message <- blank_to_na(fit_info$error)
  is_numerically_stable <- fit_ok && converged && !detect_instability_message(paste(warning_message, error_message, collapse = " | "))
  
  data.frame(
    dataset = dataset_name,
    model_id = spec$model_id[1L],
    model_family = spec$model_family[1L],
    selection_pool = spec$selection_pool[1L],
    latency_dist = spec$latency_dist[1L],
    latency_param = spec$latency_param[1L],
    covariate_structure = spec$covariate_structure[1L],
    link = spec$link[1L],
    incidence_target = spec$incidence_target[1L],
    incidence_note = spec$incidence_note[1L],
    surv_formula = formula_texts$surv_formula,
    incidence_formula = formula_texts$incidence_formula,
    latency_formula = formula_texts$latency_formula,
    fit_status = if (!fit_ok) "error" else if (!converged) "nonconverged" else "fitted",
    fit_ok = fit_ok,
    converged = converged,
    is_numerically_stable = is_numerically_stable,
    warning_message = warning_message,
    error_message = error_message,
    elapsed_sec = fit_info$elapsed_sec,
    n_fit = dataset_stats$n_fit,
    n_event = dataset_stats$n_event,
    n_censor = dataset_stats$n_censor,
    n_zero_time = dataset_stats$zero_time_count,
    max_followup_days = dataset_stats$max_followup_days,
    max_event_time_days = dataset_stats$max_event_time_days,
    logLik = loglik_val,
    df_model = df_model,
    AIC = aic_val,
    BIC = bic_val,
    eligible_for_best_cure = FALSE,
    eligible_for_best_cox = FALSE,
    is_best_cure = FALSE,
    is_best_cox = FALSE,
    stringsAsFactors = FALSE
  )
}

classify_cure_term <- function(term, spec) {
  cov_terms <- all.vars(stats::as.formula(paste("~", make_rhs_text(spec$covariate_structure[1L]))))
  latency_param <- spec$latency_param[1L]
  
  if (!is.na(latency_param) && grepl(paste0("^", latency_param, "(\\b|\\(|\\.)"), term)) {
    return("latency_location")
  }
  
  if (term %in% c("(Intercept)", cov_terms)) {
    return("incidence_cure_fraction")
  }
  
  if (grepl("^theta(\\b|\\(|\\.)", term)) {
    return("incidence_cure_fraction")
  }
  
  if (grepl("shape|sigma", term, ignore.case = TRUE)) {
    return("auxiliary_distribution")
  }
  
  "auxiliary_or_other"
}

extract_coefficients <- function(spec, fit, dataset_name) {
  if (is.null(fit)) {
    return(NULL)
  }
  
  model_id <- spec$model_id[1L]
  model_family <- spec$model_family[1L]
  
  if (identical(model_family, "cox")) {
    cf <- tryCatch(stats::coef(fit), error = function(e) NULL)
    if (is.null(cf) || length(cf) == 0L) {
      return(NULL)
    }
    
    vc <- tryCatch(stats::vcov(fit), error = function(e) NULL)
    se <- if (!is.null(vc)) sqrt(diag(vc)) else rep(NA_real_, length(cf))
    z_stat <- cf / se
    p_val <- 2 * stats::pnorm(-abs(z_stat))
    ci <- tryCatch(stats::confint(fit), error = function(e) matrix(NA_real_, nrow = length(cf), ncol = 2L))
    if (!is.matrix(ci) || nrow(ci) != length(cf)) {
      ci <- matrix(NA_real_, nrow = length(cf), ncol = 2L)
    }
    
    out <- data.frame(
      dataset = dataset_name,
      model_id = model_id,
      model_family = model_family,
      component = "hazard_ph",
      term = names(cf),
      estimate = as.numeric(cf),
      lower95 = as.numeric(ci[, 1L]),
      upper95 = as.numeric(ci[, 2L]),
      exp_estimate = as.numeric(exp(cf)),
      exp_lower95 = as.numeric(exp(ci[, 1L])),
      exp_upper95 = as.numeric(exp(ci[, 2L])),
      se = as.numeric(se),
      statistic = as.numeric(z_stat),
      p_value = as.numeric(p_val),
      stringsAsFactors = FALSE
    )
    
    rownames(out) <- NULL
    return(out)
  }
  
  res <- fit$res
  if (is.null(res)) {
    cf <- tryCatch(stats::coef(fit), error = function(e) NULL)
    if (is.null(cf)) {
      return(NULL)
    }
    
    out <- data.frame(
      dataset = dataset_name,
      model_id = model_id,
      model_family = model_family,
      component = vapply(names(cf), function(x) classify_cure_term(x, spec), character(1)),
      term = names(cf),
      estimate = as.numeric(cf),
      lower95 = NA_real_,
      upper95 = NA_real_,
      exp_estimate = NA_real_,
      exp_lower95 = NA_real_,
      exp_upper95 = NA_real_,
      se = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      stringsAsFactors = FALSE
    )
    
    rownames(out) <- NULL
    return(out)
  }
  
  res_df <- as.data.frame(res, stringsAsFactors = FALSE)
  term <- rownames(res_df)
  rownames(res_df) <- NULL
  
  key <- sanitize_name_key(names(res_df))
  est_col <- names(res_df)[match(TRUE, key %in% c("est", "estimate"))]
  lower_col <- names(res_df)[match(TRUE, grepl("^l95|^lower95|^lcl", key))]
  upper_col <- names(res_df)[match(TRUE, grepl("^u95|^upper95|^ucl", key))]
  se_col <- names(res_df)[match(TRUE, key %in% c("se", "stderr", "stderror", "stderror"))]
  stat_col <- names(res_df)[match(TRUE, key %in% c("z", "statistic", "wald", "teststatistic"))]
  p_col <- names(res_df)[match(TRUE, key %in% c("p", "pvalue", "pr", "prz"))]
  
  estimate <- if (!is.na(est_col)) as.numeric(res_df[[est_col]]) else NA_real_
  lower95 <- if (!is.na(lower_col)) as.numeric(res_df[[lower_col]]) else NA_real_
  upper95 <- if (!is.na(upper_col)) as.numeric(res_df[[upper_col]]) else NA_real_
  se <- if (!is.na(se_col)) as.numeric(res_df[[se_col]]) else NA_real_
  statistic <- if (!is.na(stat_col)) as.numeric(res_df[[stat_col]]) else NA_real_
  p_value <- if (!is.na(p_col)) as.numeric(res_df[[p_col]]) else NA_real_
  
  out <- data.frame(
    dataset = dataset_name,
    model_id = model_id,
    model_family = model_family,
    component = vapply(term, function(x) classify_cure_term(x, spec), character(1)),
    term = term,
    estimate = estimate,
    lower95 = lower95,
    upper95 = upper95,
    exp_estimate = NA_real_,
    exp_lower95 = NA_real_,
    exp_upper95 = NA_real_,
    se = se,
    statistic = statistic,
    p_value = p_value,
    stringsAsFactors = FALSE
  )
  
  rownames(out) <- NULL
  out
}

make_skipped_summary <- function(model_specs, dataset_name, reason, dataset_stats) {
  rows <- lapply(seq_len(nrow(model_specs)), function(i) {
    spec <- model_specs[i, , drop = FALSE]
    formula_texts <- make_formula_texts(spec)
    
    data.frame(
      dataset = dataset_name,
      model_id = spec$model_id[1L],
      model_family = spec$model_family[1L],
      selection_pool = spec$selection_pool[1L],
      latency_dist = spec$latency_dist[1L],
      latency_param = spec$latency_param[1L],
      covariate_structure = spec$covariate_structure[1L],
      link = spec$link[1L],
      incidence_target = spec$incidence_target[1L],
      incidence_note = spec$incidence_note[1L],
      surv_formula = formula_texts$surv_formula,
      incidence_formula = formula_texts$incidence_formula,
      latency_formula = formula_texts$latency_formula,
      fit_status = reason,
      fit_ok = FALSE,
      converged = FALSE,
      is_numerically_stable = FALSE,
      warning_message = NA_character_,
      error_message = reason,
      elapsed_sec = NA_real_,
      n_fit = dataset_stats$n_fit,
      n_event = dataset_stats$n_event,
      n_censor = dataset_stats$n_censor,
      n_zero_time = dataset_stats$zero_time_count,
      max_followup_days = dataset_stats$max_followup_days,
      max_event_time_days = dataset_stats$max_event_time_days,
      logLik = NA_real_,
      df_model = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      eligible_for_best_cure = FALSE,
      eligible_for_best_cox = FALSE,
      is_best_cure = FALSE,
      is_best_cox = FALSE,
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows_fill(rows)
}

select_best_model <- function(summary_df, selection_pool) {
  base_mask <- summary_df$selection_pool == selection_pool &
    summary_df$fit_ok &
    summary_df$converged &
    is.finite(summary_df$AIC)
  
  stable_mask <- base_mask & summary_df$is_numerically_stable
  
  if (any(stable_mask)) {
    idx <- which(stable_mask)[which.min(summary_df$AIC[stable_mask])]
    return(list(model_id = summary_df$model_id[idx], selection_rule = "stable_aic"))
  }
  
  if (any(base_mask)) {
    idx <- which(base_mask)[which.min(summary_df$AIC[base_mask])]
    return(list(model_id = summary_df$model_id[idx], selection_rule = "fallback_aic_nonstable"))
  }
  
  list(model_id = NA_character_, selection_rule = "no_eligible_model")
}

# 🔴 Define: prediction and plotting helpers ===============================

## 🟠 Build: marginal survival prediction utilities ===============================
get_time_est_df <- function(x, default_times = NULL) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  
  nms <- names(x)
  key <- sanitize_name_key(nms)
  
  time_col <- if (any(key == "time")) nms[which(key == "time")[1L]] else NA_character_
  est_col <- if (any(key %in% c("est", "estimate"))) nms[which(key %in% c("est", "estimate"))[1L]] else NA_character_
  
  if (is.na(time_col)) {
    if (is.null(default_times)) {
      stop("No time column found and default_times is NULL.", call. = FALSE)
    }
    x$time <- default_times
    time_col <- "time"
  }
  
  if (is.na(est_col)) {
    numeric_cols <- nms[vapply(x, is.numeric, FUN.VALUE = logical(1))]
    numeric_cols <- setdiff(numeric_cols, time_col)
    if (length(numeric_cols) == 0L) {
      stop("No estimate column found in prediction output.", call. = FALSE)
    }
    est_col <- numeric_cols[1L]
  }
  
  out <- data.frame(
    time_days = as.numeric(x[[time_col]]),
    survival = as.numeric(x[[est_col]]),
    stringsAsFactors = FALSE
  )
  
  out <- out[order(out$time_days), , drop = FALSE]
  rownames(out) <- NULL
  out
}

predict_marginal_survival_cox <- function(fit, newdata, times) {
  times <- sort(unique(as.numeric(times)))
  max_event_time <- if (!is.null(fit$y) && nrow(fit$y) > 0L) {
    event_rows <- fit$y[, 2L] == 1
    if (any(event_rows)) max(fit$y[event_rows, 1L], na.rm = TRUE) else NA_real_
  } else {
    NA_real_
  }
  
  bh <- survival::basehaz(fit, centered = FALSE)
  if (nrow(bh) == 0L || !is.finite(max_event_time)) {
    return(data.frame(time_days = times, survival = rep(NA_real_, length(times)), stringsAsFactors = FALSE))
  }
  
  bh <- stats::aggregate(hazard ~ time, data = bh, FUN = max)
  bh <- bh[order(bh$time), , drop = FALSE]
  
  base_time <- c(0, bh$time)
  base_hazard <- c(0, bh$hazard)
  
  H0 <- step_interp(base_time, base_hazard, pmin(times, max_event_time))
  lp <- stats::predict(fit, newdata = newdata, type = "lp")
  surv_mat <- exp(-outer(H0, exp(lp), `*`))
  marginal_surv <- rowMeans(surv_mat, na.rm = FALSE)
  marginal_surv[times > max_event_time] <- NA_real_
  
  data.frame(
    time_days = times,
    survival = as.numeric(marginal_surv),
    stringsAsFactors = FALSE
  )
}

predict_marginal_survival_cure <- function(fit, newdata, times) {
  times <- sort(unique(as.numeric(times)))
  
  pred_tidy <- tryCatch(
    summary(
      fit,
      newdata = newdata,
      type = "survival",
      t = times,
      ci = FALSE,
      B = 0,
      tidy = TRUE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(pred_tidy)) {
    pred_tidy <- as.data.frame(pred_tidy, stringsAsFactors = FALSE)
    pred_df <- get_time_est_df(pred_tidy, default_times = times)
    pred_df <- stats::aggregate(survival ~ time_days, data = pred_df, FUN = mean)
    pred_df <- pred_df[order(pred_df$time_days), , drop = FALSE]
    rownames(pred_df) <- NULL
    return(pred_df)
  }
  
  pred_list <- lapply(seq_len(nrow(newdata)), function(i) {
    pred_i <- summary(
      fit,
      newdata = newdata[i, , drop = FALSE],
      type = "survival",
      t = times,
      ci = FALSE,
      B = 0,
      tidy = FALSE
    )
    
    if (is.list(pred_i)) {
      pred_i <- pred_i[[1L]]
    }
    
    get_time_est_df(pred_i, default_times = times)
  })
  
  pred_mat <- do.call(cbind, lapply(pred_list, function(d) d$survival))
  data.frame(
    time_days = times,
    survival = rowMeans(pred_mat, na.rm = FALSE),
    stringsAsFactors = FALSE
  )
}

predict_marginal_survival_from_spec <- function(spec, fit, newdata, times) {
  if (identical(spec$model_family[1L], "cox")) {
    predict_marginal_survival_cox(
      fit = fit,
      newdata = newdata,
      times = times
    )
  } else {
    predict_marginal_survival_cure(
      fit = fit,
      newdata = newdata,
      times = times
    )
  }
}

build_curve_export_df <- function(pred_df, dataset_name, plot_type, curve_role, spec, dataset_stats) {
  if (is.null(pred_df) || nrow(pred_df) == 0L) {
    return(empty_curve_df())
  }
  
  time_days <- as.numeric(pred_df$time_days)
  survival <- as.numeric(pred_df$survival)
  
  max_event_time_days <- dataset_stats$max_event_time_days
  max_followup_days <- dataset_stats$max_followup_days
  
  is_beyond_max_event <- if (is.finite(max_event_time_days)) time_days > max_event_time_days else rep(NA, length(time_days))
  is_beyond_max_followup <- if (is.finite(max_followup_days)) time_days > max_followup_days else rep(NA, length(time_days))
  is_extrapolated <- if (identical(spec$model_family[1L], "mixture_cure")) is_beyond_max_followup else rep(FALSE, length(time_days))
  
  support_limit_type <- if (identical(spec$model_family[1L], "cox")) {
    rep("max_event_time", length(time_days))
  } else {
    rep("parametric_extrapolation_allowed", length(time_days))
  }
  
  prediction_regime <- if (identical(spec$model_family[1L], "cox")) {
    ifelse(is_beyond_max_event, "not_available_beyond_max_event", "within_event_supported_range")
  } else {
    ifelse(is_beyond_max_followup, "parametric_extrapolation_beyond_max_followup", "within_observed_followup")
  }
  
  data.frame(
    dataset = dataset_name,
    plot_type = plot_type,
    curve_role = curve_role,
    model_id = spec$model_id[1L],
    model_family = spec$model_family[1L],
    selection_pool = spec$selection_pool[1L],
    latency_dist = spec$latency_dist[1L],
    covariate_structure = spec$covariate_structure[1L],
    time_days = time_days,
    time_years = time_days / days_per_year,
    survival = survival,
    max_event_time_days = max_event_time_days,
    max_followup_days = max_followup_days,
    support_limit_type = support_limit_type,
    is_beyond_max_event = is_beyond_max_event,
    is_beyond_max_followup = is_beyond_max_followup,
    is_extrapolated = is_extrapolated,
    prediction_regime = prediction_regime,
    stringsAsFactors = FALSE
  )
}

build_kyear_export_df <- function(pred_df, dataset_name, spec, dataset_stats, k_years) {
  time_days <- as.numeric(k_years * days_per_year)
  surv_vals <- step_interp(pred_df$time_days, pred_df$survival, time_days)
  
  max_event_time_days <- dataset_stats$max_event_time_days
  max_followup_days <- dataset_stats$max_followup_days
  
  is_beyond_max_event <- if (is.finite(max_event_time_days)) time_days > max_event_time_days else rep(NA, length(time_days))
  is_beyond_max_followup <- if (is.finite(max_followup_days)) time_days > max_followup_days else rep(NA, length(time_days))
  
  if (identical(spec$model_family[1L], "cox")) {
    surv_vals[is_beyond_max_event] <- NA_real_
    support_limit_type <- rep("max_event_time", length(time_days))
    prediction_regime <- ifelse(is_beyond_max_event, "not_available_beyond_max_event", "within_event_supported_range")
    prediction_note <- ifelse(is_beyond_max_event, "NA beyond last observed event time for Cox benchmark.", NA_character_)
    is_extrapolated <- rep(FALSE, length(time_days))
  } else {
    support_limit_type <- rep("parametric_extrapolation_allowed", length(time_days))
    prediction_regime <- ifelse(is_beyond_max_followup, "parametric_extrapolation_beyond_max_followup", "within_observed_followup")
    prediction_note <- ifelse(is_beyond_max_followup, "Model-based extrapolation beyond maximum observed follow-up.", NA_character_)
    is_extrapolated <- is_beyond_max_followup
  }
  
  data.frame(
    dataset = dataset_name,
    model_id = spec$model_id[1L],
    model_family = spec$model_family[1L],
    selection_pool = spec$selection_pool[1L],
    latency_dist = spec$latency_dist[1L],
    covariate_structure = spec$covariate_structure[1L],
    k_year = k_years,
    time_days = time_days,
    time_years = k_years,
    survival = surv_vals,
    max_event_time_days = max_event_time_days,
    max_followup_days = max_followup_days,
    support_limit_type = support_limit_type,
    is_beyond_max_event = is_beyond_max_event,
    is_beyond_max_followup = is_beyond_max_followup,
    is_extrapolated = is_extrapolated,
    prediction_regime = prediction_regime,
    prediction_note = prediction_note,
    stringsAsFactors = FALSE
  )
}

plot_fit_df <- function(curve_df, file_path, dataset_label, max_event_time_days, max_followup_days) {
  if (nrow(curve_df) == 0L) {
    return(invisible(NULL))
  }
  
  tryCatch({
    grDevices::png(filename = file_path, width = 2000, height = 1600, res = 220)
    on.exit(grDevices::dev.off(), add = TRUE)
    
    x_all <- curve_df$time_years[is.finite(curve_df$time_years)]
    y_all <- curve_df$survival[is.finite(curve_df$survival)]
    
    xlim <- if (length(x_all) > 0L) range(x_all, na.rm = TRUE) else c(0, 1)
    ylim <- if (length(y_all) > 0L) range(c(0, 1, y_all), na.rm = TRUE) else c(0, 1)
    
    graphics::plot(
      NA,
      xlim = xlim,
      ylim = ylim,
      xlab = "Years since cohort entry",
      ylab = "Survival probability",
      main = paste0("Step5 fit-range: ", toupper(dataset_label))
    )
    
    if (is.finite(max_event_time_days)) {
      graphics::abline(v = max_event_time_days / days_per_year, col = "gray50", lty = 3, lwd = 1.5)
    }
    if (is.finite(max_followup_days)) {
      graphics::abline(v = max_followup_days / days_per_year, col = "gray70", lty = 2, lwd = 1.5)
    }
    
    cox_df <- curve_df[curve_df$curve_role == "best_cox" & is.finite(curve_df$survival), , drop = FALSE]
    if (nrow(cox_df) > 0L) {
      graphics::lines(cox_df$time_years, cox_df$survival, col = "royalblue3", lty = 2, lwd = 2.5)
    }
    
    cure_df <- curve_df[curve_df$curve_role == "best_cure" & is.finite(curve_df$survival), , drop = FALSE]
    if (nrow(cure_df) > 0L) {
      graphics::lines(cure_df$time_years, cure_df$survival, col = "firebrick3", lty = 1, lwd = 2.8)
    }
    
    legend_labels <- c()
    legend_cols <- c()
    legend_lty <- c()
    legend_lwd <- c()
    
    if (nrow(cox_df) > 0L) {
      legend_labels <- c(legend_labels, paste0("Best Cox: ", unique(cox_df$model_id)[1L]))
      legend_cols <- c(legend_cols, "royalblue3")
      legend_lty <- c(legend_lty, 2)
      legend_lwd <- c(legend_lwd, 2.5)
    }
    if (nrow(cure_df) > 0L) {
      legend_labels <- c(legend_labels, paste0("Best Cure: ", unique(cure_df$model_id)[1L]))
      legend_cols <- c(legend_cols, "firebrick3")
      legend_lty <- c(legend_lty, 1)
      legend_lwd <- c(legend_lwd, 2.8)
    }
    if (is.finite(max_event_time_days)) {
      legend_labels <- c(legend_labels, "Max event time")
      legend_cols <- c(legend_cols, "gray50")
      legend_lty <- c(legend_lty, 3)
      legend_lwd <- c(legend_lwd, 1.5)
    }
    if (is.finite(max_followup_days)) {
      legend_labels <- c(legend_labels, "Max follow-up")
      legend_cols <- c(legend_cols, "gray70")
      legend_lty <- c(legend_lty, 2)
      legend_lwd <- c(legend_lwd, 1.5)
    }
    
    if (length(legend_labels) > 0L) {
      graphics::legend(
        "topright",
        legend = legend_labels,
        col = legend_cols,
        lty = legend_lty,
        lwd = legend_lwd,
        bty = "n"
      )
    }
  }, error = function(e) {
    message("Fit-range PNG could not be written for ", dataset_label, ": ", conditionMessage(e))
    invisible(NULL)
  })
}

plot_extrapolation_df <- function(curve_df, file_path, dataset_label, max_event_time_days, max_followup_days) {
  if (nrow(curve_df) == 0L) {
    return(invisible(NULL))
  }
  
  tryCatch({
    grDevices::png(filename = file_path, width = 2000, height = 1600, res = 220)
    on.exit(grDevices::dev.off(), add = TRUE)
    
    x_all <- curve_df$time_years[is.finite(curve_df$time_years)]
    y_all <- curve_df$survival[is.finite(curve_df$survival)]
    
    xlim <- if (length(x_all) > 0L) range(x_all, na.rm = TRUE) else c(0, 1)
    ylim <- if (length(y_all) > 0L) range(c(0, 1, y_all), na.rm = TRUE) else c(0, 1)
    
    graphics::plot(
      NA,
      xlim = xlim,
      ylim = ylim,
      xlab = "Years since cohort entry",
      ylab = "Survival probability",
      main = paste0("Step5 parametric extrapolation: ", toupper(dataset_label))
    )
    
    if (is.finite(max_event_time_days)) {
      graphics::abline(v = max_event_time_days / days_per_year, col = "gray50", lty = 3, lwd = 1.5)
    }
    if (is.finite(max_followup_days)) {
      graphics::abline(v = max_followup_days / days_per_year, col = "gray70", lty = 2, lwd = 1.5)
    }
    
    within_df <- curve_df[!curve_df$is_extrapolated & is.finite(curve_df$survival), , drop = FALSE]
    extra_df <- curve_df[curve_df$is_extrapolated & is.finite(curve_df$survival), , drop = FALSE]
    
    if (nrow(within_df) > 0L) {
      graphics::lines(within_df$time_years, within_df$survival, col = "firebrick3", lty = 1, lwd = 2.8)
    }
    if (nrow(extra_df) > 0L) {
      graphics::lines(extra_df$time_years, extra_df$survival, col = "firebrick3", lty = 2, lwd = 2.8)
    }
    
    legend_labels <- c(
      paste0("Best Cure: ", unique(curve_df$model_id)[1L], " (within follow-up)"),
      paste0("Best Cure: ", unique(curve_df$model_id)[1L], " (extrapolated)")
    )
    legend_cols <- c("firebrick3", "firebrick3")
    legend_lty <- c(1, 2)
    legend_lwd <- c(2.8, 2.8)
    
    if (is.finite(max_event_time_days)) {
      legend_labels <- c(legend_labels, "Max event time")
      legend_cols <- c(legend_cols, "gray50")
      legend_lty <- c(legend_lty, 3)
      legend_lwd <- c(legend_lwd, 1.5)
    }
    if (is.finite(max_followup_days)) {
      legend_labels <- c(legend_labels, "Max follow-up")
      legend_cols <- c(legend_cols, "gray70")
      legend_lty <- c(legend_lty, 2)
      legend_lwd <- c(legend_lwd, 1.5)
    }
    
    graphics::legend(
      "topright",
      legend = legend_labels,
      col = legend_cols,
      lty = legend_lty,
      lwd = legend_lwd,
      bty = "n"
    )
  }, error = function(e) {
    message("Extrapolation PNG could not be written for ", dataset_label, ": ", conditionMessage(e))
    invisible(NULL)
  })
}

# 🔴 Read: merged source and create primary endpoint ===============================

## 🟠 Import: merged CSV and normalize raw fields ===============================
read_col_classes <- c(
  id = "character",
  site = "character",
  date_entry = "character",
  days_followup = "character",
  status_num = "character",
  sex_num = "character",
  age_exact_entry = "character"
)

dat_raw <- utils::read.csv(
  file = data_file_full,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN", ".", "NULL"),
  colClasses = read_col_classes
)

required_cols <- c(
  "id",
  "site",
  "date_entry",
  "days_followup",
  "status_num",
  "sex_num",
  "age_exact_entry"
)

missing_cols <- setdiff(required_cols, names(dat_raw))
if (length(missing_cols) > 0L) {
  stop("Required column(s) missing from input CSV: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

dat_raw$site_original <- dat_raw$site
dat_raw$site <- normalize_site(dat_raw$site)
dat_raw$id <- as.character(dat_raw$id)
dat_raw$date_entry <- as.character(dat_raw$date_entry)

dat_raw$days_followup <- coerce_numeric_clean(dat_raw$days_followup)
dat_raw$status_num <- coerce_status_num(dat_raw$status_num)
dat_raw$sex_num_missing <- is_missing_string(dat_raw$sex_num)
dat_raw$sex_num <- coerce_binary_sex(dat_raw$sex_num)
dat_raw$age_exact_entry_missing <- is_missing_string(dat_raw$age_exact_entry)
dat_raw$age_exact_entry <- coerce_numeric_clean(dat_raw$age_exact_entry)
dat_raw$subject_key <- paste(dat_raw$site, dat_raw$id, sep = "::")

## 🟠 Check: core analysis integrity ===============================
validate_core_inputs(dat_raw)

## 🟠 Derive: transition-only endpoint variables ===============================
dat_core <- dat_raw
dat_core$time_primary <- dat_core$days_followup
dat_core$event_primary <- ifelse(dat_core$status_num == 1L, 1L, 0L)
dat_core$status_primary <- ifelse(dat_core$event_primary == 1L, "transition", "censored")
dat_core$step5_complete_case <- !is.na(dat_core$sex_num) & !is.na(dat_core$age_exact_entry)

# 🔴 Split: frozen scaling reference and dataset branches ===============================

## 🟠 Compute: merged-reference standardization constants ===============================
analysis_frame_merged <- dat_core[dat_core$step5_complete_case, , drop = FALSE]
analysis_frame_merged <- analysis_frame_merged[order(analysis_frame_merged$site, analysis_frame_merged$id), , drop = FALSE]
rownames(analysis_frame_merged) <- NULL

if (nrow(analysis_frame_merged) == 0L) {
  stop("No complete cases remained for Step5 after applying sex_num and age_exact_entry requirements.", call. = FALSE)
}

scaling_ref <- compute_scaling_reference(analysis_frame_merged)

scaling_reference_df <- data.frame(
  reference_dataset = "merged_complete_case_analysis_set",
  reference_n = nrow(analysis_frame_merged),
  age_variable = "age_exact_entry",
  sex_variable = "sex_num",
  age_mean_merged = scaling_ref$age_mean,
  age_sd_merged = scaling_ref$age_sd,
  age_two_sd_merged = scaling_ref$age_two_sd,
  sex_mean_merged = scaling_ref$sex_mean,
  zero_time_offset_days = zero_time_offset_days,
  scaling_note = "Continuous age centered and divided by 2 SD; binary sex centered only; interaction built from transformed inputs; merged reference applied to merged/PNU/SNU.",
  stringsAsFactors = FALSE
)

## 🟠 Assemble: merged PNU SNU analysis frames ===============================
analysis_frame_merged <- apply_step5_scaling(analysis_frame_merged, scaling_ref, zero_time_offset_days)

analysis_list <- list(
  merged = analysis_frame_merged,
  pnu = analysis_frame_merged[analysis_frame_merged$site == "PNU", , drop = FALSE],
  snu = analysis_frame_merged[analysis_frame_merged$site == "SNU", , drop = FALSE]
)

core_list <- list(
  merged = dat_core,
  pnu = dat_core[dat_core$site == "PNU", , drop = FALSE],
  snu = dat_core[dat_core$site == "SNU", , drop = FALSE]
)

dataset_manifest <- bind_rows_fill(lapply(names(core_list), function(ds) {
  core_ds <- core_list[[ds]]
  fit_ds <- analysis_list[[ds]]
  
  data.frame(
    dataset = ds,
    gate_pass_input = isTRUE(gate_pass[[ds]]),
    n_core_valid_source = nrow(core_ds),
    n_missing_sex_step5 = sum(is.na(core_ds$sex_num)),
    n_missing_age_step5 = sum(is.na(core_ds$age_exact_entry)),
    n_excluded_missing_step5_covariates = nrow(core_ds) - nrow(fit_ds),
    n_fit_analysis = nrow(fit_ds),
    n_event = sum(fit_ds$event_primary == 1L),
    n_censor = sum(fit_ds$event_primary == 0L),
    n_status0_source = sum(core_ds$status_num == 0L),
    n_status1_source = sum(core_ds$status_num == 1L),
    n_status2_source = sum(core_ds$status_num == 2L),
    zero_time_count = sum(fit_ds$time_primary <= 0),
    max_followup_days = if (nrow(fit_ds) > 0L) max(fit_ds$time_primary, na.rm = TRUE) else NA_real_,
    max_event_time_days = if (sum(fit_ds$event_primary == 1L) > 0L) max(fit_ds$time_primary[fit_ds$event_primary == 1L], na.rm = TRUE) else NA_real_,
    fitplot_csv = file.path(export_path, paste0("step5_fitplot_", ds, ".csv")),
    fitplot_png = file.path(export_path, paste0("step5_fitplot_", ds, ".png")),
    extrapolation_csv = file.path(export_path, paste0("step5_extrapolation_", ds, ".csv")),
    extrapolation_png = file.path(export_path, paste0("step5_extrapolation_", ds, ".png")),
    bundle_rds = file.path(export_path, paste0("step5_models_", ds, ".rds")),
    best_cox_id = NA_character_,
    best_cox_selection_rule = NA_character_,
    best_cure_id = NA_character_,
    best_cure_selection_rule = NA_character_,
    stringsAsFactors = FALSE
  )
}))

# 🔴 Register: Step5 model grids and storage objects ===============================

## 🟠 Define: benchmark and cure candidate catalog ===============================
model_specs <- make_model_specs()

model_summary_parts <- list()
coef_parts <- list()
kyear_parts <- list()

# 🔴 Execute: dataset-wise fitting and assessment ===============================

## 🟠 Loop: fit models within each cohort branch ===============================
for (dataset_name in names(analysis_list)) {
  data_i <- analysis_list[[dataset_name]]
  prediction_newdata <- data_i[, c("c_sex", "z_age", "int_sex_age"), drop = FALSE]
  
  dataset_stats <- list(
    n_fit = nrow(data_i),
    n_event = sum(data_i$event_primary == 1L),
    n_censor = sum(data_i$event_primary == 0L),
    zero_time_count = sum(data_i$time_primary <= 0),
    max_followup_days = if (nrow(data_i) > 0L) max(data_i$time_primary, na.rm = TRUE) else NA_real_,
    max_event_time_days = if (sum(data_i$event_primary == 1L) > 0L) max(data_i$time_primary[data_i$event_primary == 1L], na.rm = TRUE) else NA_real_
  )
  
  fitplot_csv_path <- dataset_manifest$fitplot_csv[dataset_manifest$dataset == dataset_name]
  fitplot_png_path <- dataset_manifest$fitplot_png[dataset_manifest$dataset == dataset_name]
  extrap_csv_path <- dataset_manifest$extrapolation_csv[dataset_manifest$dataset == dataset_name]
  extrap_png_path <- dataset_manifest$extrapolation_png[dataset_manifest$dataset == dataset_name]
  bundle_rds_path <- dataset_manifest$bundle_rds[dataset_manifest$dataset == dataset_name]
  
  fits_this <- setNames(vector("list", length = nrow(model_specs)), model_specs$model_id)
  summary_rows <- list()
  coef_rows <- list()
  kyear_rows <- list()
  fitplot_df <- empty_curve_df()
  extrapolation_df <- empty_curve_df()
  
  skip_reason <- NULL
  if (!isTRUE(gate_pass[[dataset_name]])) {
    skip_reason <- "skipped_gate"
  } else if (nrow(data_i) == 0L) {
    skip_reason <- "skipped_no_complete_cases"
  } else if (dataset_stats$n_event == 0L) {
    skip_reason <- "skipped_no_events"
  }
  
  if (!is.null(skip_reason)) {
    summary_this <- make_skipped_summary(
      model_specs = model_specs,
      dataset_name = dataset_name,
      reason = skip_reason,
      dataset_stats = dataset_stats
    )
    
    model_summary_parts[[dataset_name]] <- summary_this
    coef_parts[[dataset_name]] <- empty_coeff_df()
    kyear_parts[[dataset_name]] <- empty_kyear_df()
    
    utils::write.csv(fitplot_df, fitplot_csv_path, row.names = FALSE)
    utils::write.csv(extrapolation_df, extrap_csv_path, row.names = FALSE)
    
    bundle_obj <- list(
      step = "Step5",
      dataset_label = dataset_name,
      created_at = as.character(Sys.time()),
      input_file = normalizePath(data_file_full, winslash = "/", mustWork = FALSE),
      gate = list(pass = isTRUE(gate_pass[[dataset_name]]), reason = skip_reason),
      charter = list(
        time_origin = "date_entry",
        primary_time_variable = "days_followup",
        primary_model_time_variable = "time_model",
        primary_event = "transition",
        censoring = c("right_censoring", "remission"),
        endpoint_note = "status_num == 1 treated as event; status_num in c(0,2) treated as censoring.",
        zero_time_offset_days = zero_time_offset_days
      ),
      standardization = scaling_reference_df,
      dataset_stats = dataset_stats,
      model_specs = model_specs,
      analysis_frame = data_i,
      fits = fits_this,
      model_summary = summary_this,
      coefficient_table = empty_coeff_df(),
      kyear_predictions = empty_kyear_df(),
      fitplot_data = fitplot_df,
      extrapolation_data = extrapolation_df,
      selected_models = list(best_cure_id = NA_character_, best_cox_id = NA_character_),
      session_info = utils::capture.output(sessionInfo())
    )
    
    saveRDS(bundle_obj, bundle_rds_path, compress = "xz")
    next
  }
  
  for (i in seq_len(nrow(model_specs))) {
    spec <- model_specs[i, , drop = FALSE]
    model_id <- spec$model_id[1L]
    
    fit_info <- fit_model_from_spec(
      spec = spec,
      dat = data_i,
      cox_ties = cox_ties,
      cure_optim_method = cure_optim_method,
      optim_maxit = optim_maxit,
      survreg_maxiter = survreg_maxiter
    )
    
    fits_this[[model_id]] <- fit_info$fit
    
    summary_rows[[model_id]] <- extract_fit_summary(
      spec = spec,
      fit_info = fit_info,
      dataset_name = dataset_name,
      dataset_stats = dataset_stats
    )
    
    coef_rows[[model_id]] <- extract_coefficients(
      spec = spec,
      fit = fit_info$fit,
      dataset_name = dataset_name
    )
  }
  
  summary_this <- bind_rows_fill(summary_rows)
  coef_this <- bind_rows_fill(coef_rows)
  
  best_cure_sel <- select_best_model(summary_this, "cure_candidate")
  best_cox_sel <- select_best_model(summary_this, "benchmark")
  
  if (nrow(summary_this) > 0L) {
    summary_this$eligible_for_best_cure <- with(summary_this, selection_pool == "cure_candidate" & fit_ok & converged & is.finite(AIC) & is_numerically_stable)
    summary_this$eligible_for_best_cox <- with(summary_this, selection_pool == "benchmark" & fit_ok & converged & is.finite(AIC) & is_numerically_stable)
    
    if (!is.na(best_cure_sel$model_id)) {
      summary_this$is_best_cure[summary_this$model_id == best_cure_sel$model_id] <- TRUE
    }
    if (!is.na(best_cox_sel$model_id)) {
      summary_this$is_best_cox[summary_this$model_id == best_cox_sel$model_id] <- TRUE
    }
  }
  
  best_cure_id <- best_cure_sel$model_id
  best_cox_id <- best_cox_sel$model_id
  
  dataset_manifest$best_cure_id[dataset_manifest$dataset == dataset_name] <- best_cure_id
  dataset_manifest$best_cure_selection_rule[dataset_manifest$dataset == dataset_name] <- best_cure_sel$selection_rule
  dataset_manifest$best_cox_id[dataset_manifest$dataset == dataset_name] <- best_cox_id
  dataset_manifest$best_cox_selection_rule[dataset_manifest$dataset == dataset_name] <- best_cox_sel$selection_rule
  
  prediction_times_days <- as.numeric(k_years * days_per_year)
  
  for (i in seq_len(nrow(model_specs))) {
    spec <- model_specs[i, , drop = FALSE]
    model_id <- spec$model_id[1L]
    fit_obj <- fits_this[[model_id]]
    
    model_row <- summary_this[summary_this$model_id == model_id, , drop = FALSE]
    if (nrow(model_row) == 0L || !isTRUE(model_row$fit_ok[1L]) || !isTRUE(model_row$converged[1L])) {
      next
    }
    
    pred_k <- tryCatch(
      predict_marginal_survival_from_spec(
        spec = spec,
        fit = fit_obj,
        newdata = prediction_newdata,
        times = prediction_times_days
      ),
      error = function(e) e
    )
    
    if (inherits(pred_k, "error")) {
      kyear_rows[[model_id]] <- data.frame(
        dataset = dataset_name,
        model_id = model_id,
        model_family = spec$model_family[1L],
        selection_pool = spec$selection_pool[1L],
        latency_dist = spec$latency_dist[1L],
        covariate_structure = spec$covariate_structure[1L],
        k_year = k_years,
        time_days = prediction_times_days,
        time_years = k_years,
        survival = NA_real_,
        max_event_time_days = dataset_stats$max_event_time_days,
        max_followup_days = dataset_stats$max_followup_days,
        support_limit_type = if (identical(spec$model_family[1L], "cox")) "max_event_time" else "parametric_extrapolation_allowed",
        is_beyond_max_event = prediction_times_days > dataset_stats$max_event_time_days,
        is_beyond_max_followup = prediction_times_days > dataset_stats$max_followup_days,
        is_extrapolated = if (identical(spec$model_family[1L], "mixture_cure")) prediction_times_days > dataset_stats$max_followup_days else FALSE,
        prediction_regime = "prediction_error",
        prediction_note = conditionMessage(pred_k),
        stringsAsFactors = FALSE
      )
      next
    }
    
    kyear_rows[[model_id]] <- build_kyear_export_df(
      pred_df = pred_k,
      dataset_name = dataset_name,
      spec = spec,
      dataset_stats = dataset_stats,
      k_years = k_years
    )
  }
  
  kyear_this <- bind_rows_fill(kyear_rows)
  
  fit_times <- make_regular_times(dataset_stats$max_followup_days, fitplot_grid_n)
  fitplot_parts <- list()
  
  if (!is.na(best_cox_id) && !is.null(fits_this[[best_cox_id]])) {
    spec_best_cox <- model_specs[model_specs$model_id == best_cox_id, , drop = FALSE]
    pred_best_cox_fit <- tryCatch(
      predict_marginal_survival_from_spec(
        spec = spec_best_cox,
        fit = fits_this[[best_cox_id]],
        newdata = prediction_newdata,
        times = fit_times
      ),
      error = function(e) NULL
    )
    
    if (!is.null(pred_best_cox_fit)) {
      fitplot_parts[["best_cox"]] <- build_curve_export_df(
        pred_df = pred_best_cox_fit,
        dataset_name = dataset_name,
        plot_type = "fitplot",
        curve_role = "best_cox",
        spec = spec_best_cox,
        dataset_stats = dataset_stats
      )
    }
  }
  
  if (!is.na(best_cure_id) && !is.null(fits_this[[best_cure_id]])) {
    spec_best_cure <- model_specs[model_specs$model_id == best_cure_id, , drop = FALSE]
    pred_best_cure_fit <- tryCatch(
      predict_marginal_survival_from_spec(
        spec = spec_best_cure,
        fit = fits_this[[best_cure_id]],
        newdata = prediction_newdata,
        times = fit_times
      ),
      error = function(e) NULL
    )
    
    if (!is.null(pred_best_cure_fit)) {
      fitplot_parts[["best_cure"]] <- build_curve_export_df(
        pred_df = pred_best_cure_fit,
        dataset_name = dataset_name,
        plot_type = "fitplot",
        curve_role = "best_cure",
        spec = spec_best_cure,
        dataset_stats = dataset_stats
      )
    }
  }
  
  fitplot_df <- bind_rows_fill(fitplot_parts)
  utils::write.csv(fitplot_df, fitplot_csv_path, row.names = FALSE)
  
  if (isTRUE(write_fitplot_png) && nrow(fitplot_df) > 0L) {
    plot_fit_df(
      curve_df = fitplot_df,
      file_path = fitplot_png_path,
      dataset_label = dataset_name,
      max_event_time_days = dataset_stats$max_event_time_days,
      max_followup_days = dataset_stats$max_followup_days
    )
  }
  
  extrap_horizon_days <- max(dataset_stats$max_followup_days, extrapolation_horizon_days, na.rm = TRUE)
  extrap_times <- make_regular_times(extrap_horizon_days, extrapolation_grid_n)
  extrapolation_parts <- list()
  
  if (!is.na(best_cure_id) && !is.null(fits_this[[best_cure_id]])) {
    spec_best_cure <- model_specs[model_specs$model_id == best_cure_id, , drop = FALSE]
    pred_best_cure_extra <- tryCatch(
      predict_marginal_survival_from_spec(
        spec = spec_best_cure,
        fit = fits_this[[best_cure_id]],
        newdata = prediction_newdata,
        times = extrap_times
      ),
      error = function(e) NULL
    )
    
    if (!is.null(pred_best_cure_extra)) {
      extrapolation_parts[["best_cure"]] <- build_curve_export_df(
        pred_df = pred_best_cure_extra,
        dataset_name = dataset_name,
        plot_type = "extrapolation",
        curve_role = "best_cure",
        spec = spec_best_cure,
        dataset_stats = dataset_stats
      )
    }
  }
  
  extrapolation_df <- bind_rows_fill(extrapolation_parts)
  utils::write.csv(extrapolation_df, extrap_csv_path, row.names = FALSE)
  
  if (isTRUE(write_extrapolation_png) && nrow(extrapolation_df) > 0L) {
    plot_extrapolation_df(
      curve_df = extrapolation_df,
      file_path = extrap_png_path,
      dataset_label = dataset_name,
      max_event_time_days = dataset_stats$max_event_time_days,
      max_followup_days = dataset_stats$max_followup_days
    )
  }
  
  model_summary_parts[[dataset_name]] <- summary_this
  coef_parts[[dataset_name]] <- if (nrow(coef_this) > 0L) coef_this else empty_coeff_df()
  kyear_parts[[dataset_name]] <- if (nrow(kyear_this) > 0L) kyear_this else empty_kyear_df()
  
  bundle_obj <- list(
    step = "Step5",
    dataset_label = dataset_name,
    created_at = as.character(Sys.time()),
    input_file = normalizePath(data_file_full, winslash = "/", mustWork = FALSE),
    gate = list(pass = isTRUE(gate_pass[[dataset_name]]), reason = NA_character_),
    charter = list(
      time_origin = "date_entry",
      primary_time_variable = "days_followup",
      primary_model_time_variable = "time_model",
      primary_event = "transition",
      censoring = c("right_censoring", "remission"),
      endpoint_note = "status_num == 1 treated as event; status_num in c(0,2) treated as censoring.",
      zero_time_offset_days = zero_time_offset_days,
      zero_time_offset_note = "time_model = pmax(time_primary, zero_time_offset_days) was used for fitting to guard positive-support parametric latency models."
    ),
    standardization = scaling_reference_df,
    dataset_stats = dataset_stats,
    model_specs = model_specs,
    analysis_frame = data_i,
    fits = fits_this,
    model_summary = summary_this,
    coefficient_table = if (nrow(coef_this) > 0L) coef_this else empty_coeff_df(),
    kyear_predictions = if (nrow(kyear_this) > 0L) kyear_this else empty_kyear_df(),
    fitplot_data = fitplot_df,
    extrapolation_data = extrapolation_df,
    selected_models = list(
      best_cure_id = best_cure_id,
      best_cure_selection_rule = best_cure_sel$selection_rule,
      best_cox_id = best_cox_id,
      best_cox_selection_rule = best_cox_sel$selection_rule
    ),
    prediction_notes = list(
      cure_parameterization = "flexsurvcure theta = cure fraction",
      susceptible_probability = "susceptible proportion = 1 - theta",
      cox_note = "Cox benchmark predictions are exported as NA beyond the last observed event time.",
      cure_note = "Parametric cure predictions are retained beyond maximum follow-up and flagged as extrapolated."
    ),
    session_info = utils::capture.output(sessionInfo())
  )
  
  saveRDS(bundle_obj, bundle_rds_path, compress = "xz")
}

# 🔴 Export: combined summaries and self-contained bundles ===============================

## 🟠 Write: CSV outputs across all datasets ===============================
model_summary_all <- bind_rows_fill(model_summary_parts)
coef_all <- bind_rows_fill(coef_parts)
kyear_all <- bind_rows_fill(kyear_parts)

if (nrow(coef_all) == 0L) {
  coef_all <- empty_coeff_df()
}

if (nrow(kyear_all) == 0L) {
  kyear_all <- empty_kyear_df()
}

utils::write.csv(
  scaling_reference_df,
  file.path(export_path, "step5_scaling_reference.csv"),
  row.names = FALSE
)

utils::write.csv(
  dataset_manifest,
  file.path(export_path, "step5_dataset_manifest.csv"),
  row.names = FALSE
)

utils::write.csv(
  model_summary_all,
  file.path(export_path, "step5_model_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  coef_all,
  file.path(export_path, "step5_model_coefficients.csv"),
  row.names = FALSE
)

utils::write.csv(
  kyear_all,
  file.path(export_path, "step5_model_kyear_survival.csv"),
  row.names = FALSE
)

invisible(list(
  input_file = normalizePath(data_file_full, winslash = "/", mustWork = FALSE),
  export_path = normalizePath(export_path, winslash = "/", mustWork = FALSE),
  scaling_reference = scaling_reference_df,
  dataset_manifest = dataset_manifest,
  model_summary = model_summary_all,
  model_coefficients = coef_all,
  model_kyear_survival = kyear_all
))