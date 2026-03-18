# 🔴 Configure: paths and execution switches ===============================

## 🟠 Define: data path and export location ===============================
data_path <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
export_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦5.Step5_MLE MCM/attachments'

## 🟠 Set: user-facing analysis parameters ===============================
script_version <- "step5_mle_cure_grid_v2_0"
random_seed <- 20260317L

age_var_preferred <- "age_exact_entry"
id_var <- "id"
site_var <- "site"
sex_var <- "sex_num"
time_var <- "days_followup"
status_var <- "status_num"

pnu_site_label <- "PNU"
snu_site_label <- "SNU"
site_reference_label <- "PNU"

days_per_year <- 365.25
time_zero_epsilon_days <- 0.5
prediction_years <- 1:10
prob_floor <- 1e-12
boundary_prob_threshold <- 1e-4
big_nll <- 1e20

reuse_existing_model_rds <- TRUE
overwrite_model_rds <- FALSE
force_refit_if_script_version_mismatch <- TRUE

parametric_maxit <- 2000L
parametric_reltol <- 1e-10
parametric_gamma_intercept_shifts <- c(0, 1.5, -1.5)
parametric_gamma_mean_uncured_starts <- c(0.50, 0.80, 0.95)

em_max_iter <- 200L
em_tol <- 1e-6
coxph_iter_max <- 50L

write_subject_level_predictions <- TRUE
write_mean_level_predictions <- TRUE
write_param_estimates <- TRUE
write_cure_summary <- TRUE
write_vcov_long <- TRUE
write_cox_baseline_curve <- TRUE
write_fit_diagnostics <- TRUE
write_perf_metrics <- TRUE

# 🔴 Initialize: packages and global options ===============================

## 🟠 Resolve: dependencies and options ===============================
required_pkgs <- c("survival")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("다음 패키지를 먼저 설치해 주세요: ", paste(missing_pkgs, collapse = ", "))
}

set.seed(random_seed)
options(stringsAsFactors = FALSE, scipen = 999)

if (!dir.exists(export_path)) {
  dir.create(export_path, recursive = TRUE, showWarnings = FALSE)
}

# 🔴 Build: helper utilities for fitting and export ===============================

## 🟠 Create: small utilities and safe wrappers ===============================
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

collapse_terms <- function(x) {
  if (length(x) == 0L) "(none)" else paste(x, collapse = " + ")
}

write_csv_utf8 <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(path)
}

bind_rows_or_template <- function(lst, template_df) {
  if (length(lst) == 0L) template_df else do.call(rbind, lst)
}

capture_warnings <- function(expr) {
  warns <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warns <<- unique(c(warns, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warns)
}

record_manifest_row <- function(manifest_list, file_name, file_type, n_rows = NA_integer_, note = NA_character_, model_id = NA_character_) {
  manifest_list[[length(manifest_list) + 1L]] <- data.frame(
    file_name = file_name,
    file_type = file_type,
    n_rows = n_rows,
    model_id = model_id,
    note = note,
    stringsAsFactors = FALSE
  )
  manifest_list
}

trapz_base <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

safe_inverse_hessian <- function(h, tol = 1e-8) {
  out <- list(
    vcov = if (is.null(h)) matrix(NA_real_, 0L, 0L) else matrix(NA_real_, nrow(h), ncol(h)),
    singular = TRUE,
    condition_number = NA_real_
  )
  
  if (is.null(h) || length(h) == 0L || any(!is.finite(h))) {
    return(out)
  }
  
  hs <- (h + t(h)) / 2
  out$condition_number <- tryCatch(kappa(hs), error = function(e) NA_real_)
  eigvals <- tryCatch(eigen(hs, symmetric = TRUE, only.values = TRUE)$values, error = function(e) rep(NA_real_, nrow(hs)))
  
  if (all(is.finite(eigvals)) && all(eigvals > tol)) {
    vc <- tryCatch(solve(hs), error = function(e) NULL)
    if (!is.null(vc) && all(is.finite(diag(vc)))) {
      out$vcov <- vc
      out$singular <- FALSE
    }
  }
  
  out
}

norm_p_value <- function(z) {
  ifelse(is.finite(z), 2 * stats::pnorm(abs(z), lower.tail = FALSE), NA_real_)
}

wald_ci <- function(est, se, level = 0.95) {
  alpha <- 1 - level
  z <- stats::qnorm(1 - alpha / 2)
  lcl <- est - z * se
  ucl <- est + z * se
  cbind(lcl = lcl, ucl = ucl)
}

## 🟠 Build: model matrix creation with rank pruning ===============================
make_model_matrix_with_pruning <- function(data, terms, include_intercept = TRUE, prefix = "x") {
  rhs <- if (length(terms) == 0L) "1" else paste(terms, collapse = " + ")
  fml <- stats::as.formula(paste("~", rhs))
  mm0 <- stats::model.matrix(fml, data = data)
  
  if (!include_intercept && "(Intercept)" %in% colnames(mm0)) {
    mm0 <- mm0[, colnames(mm0) != "(Intercept)", drop = FALSE]
  }
  
  original_names <- colnames(mm0)
  dropped_terms <- character()
  
  if (ncol(mm0) > 0L) {
    qr_mm <- qr(mm0)
    keep_idx <- sort(qr_mm$pivot[seq_len(qr_mm$rank)])
    if (length(keep_idx) < ncol(mm0)) {
      dropped_terms <- original_names[-keep_idx]
    }
    mm0 <- mm0[, keep_idx, drop = FALSE]
    original_names <- colnames(mm0)
  }
  
  new_names <- if (ncol(mm0) > 0L) sprintf("%s_%02d", prefix, seq_len(ncol(mm0))) else character(0)
  colnames(mm0) <- new_names
  
  map <- data.frame(
    coef_name = new_names,
    term_label = original_names,
    stringsAsFactors = FALSE
  )
  
  list(
    matrix = mm0,
    map = map,
    dropped_terms = dropped_terms,
    rhs = rhs,
    include_intercept = include_intercept
  )
}

## 🟠 Build: Kaplan-Meier step evaluators for IPCW and Cox baseline ===============================
km_surv_eval <- function(eval_time, km_time, km_surv, left_limit = FALSE) {
  eval_time <- as.numeric(eval_time)
  if (length(km_time) == 0L) {
    return(rep(1, length(eval_time)))
  }
  
  tiny <- 1e-10
  t_use <- if (left_limit) pmax(eval_time - tiny, 0) else eval_time
  idx <- findInterval(t_use, km_time, left.open = FALSE, rightmost.closed = TRUE)
  out <- rep(1, length(eval_time))
  ok <- idx > 0L
  out[ok] <- km_surv[idx[ok]]
  pmin(pmax(out, 0), 1)
}

build_ipcw_helper <- function(time, event) {
  censor_event <- as.integer(event == 0L)
  sf <- survival::survfit(survival::Surv(time, censor_event) ~ 1)
  km_time <- sf$time %||% numeric(0)
  km_surv <- sf$surv %||% numeric(0)
  
  list(
    km_time = km_time,
    km_surv = km_surv,
    G_at = function(t) km_surv_eval(t, km_time, km_surv, left_limit = FALSE),
    G_left = function(t) km_surv_eval(t, km_time, km_surv, left_limit = TRUE)
  )
}

make_step_surv_fun <- function(basehaz_df, zero_tail = TRUE) {
  bh <- basehaz_df[order(basehaz_df$time), , drop = FALSE]
  bh <- bh[!duplicated(bh$time, fromLast = TRUE), , drop = FALSE]
  times <- bh$time
  surv <- exp(-bh$hazard)
  last_event_time <- if (length(times) > 0L) max(times) else NA_real_
  
  function(t) {
    t <- as.numeric(t)
    out <- rep(1, length(t))
    
    if (length(times) > 0L) {
      idx <- findInterval(t, vec = times, left.open = FALSE, rightmost.closed = TRUE)
      ok <- idx > 0L
      out[ok] <- surv[pmin(idx[ok], length(surv))]
      if (isTRUE(zero_tail) && is.finite(last_event_time)) {
        out[t > last_event_time] <- 0
      }
    }
    
    pmin(pmax(out, 0), 1)
  }
}

## 🟠 Build: AFT latency survival and density calculators ===============================
aft_surv_density <- function(family, time, eta, log_scale = NULL, tiny = 1e-12) {
  t <- pmax(as.numeric(time), tiny)
  eta <- as.numeric(eta)
  
  if (family == "exp") {
    scale_par <- exp(eta)
    logS <- -t / scale_par
    logf <- -eta - t / scale_par
  } else if (family == "weibull") {
    sigma <- exp(log_scale)
    shape <- 1 / sigma
    log_t_over_scale <- log(t) - eta
    z <- exp(shape * log_t_over_scale)
    logS <- -z
    logf <- log(shape) - shape * eta + (shape - 1) * log(t) - z
  } else if (family == "lnorm") {
    sigma <- exp(log_scale)
    z <- (log(t) - eta) / sigma
    logS <- stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)
    logf <- stats::dnorm(z, log = TRUE) - log(t) - log(sigma)
  } else if (family == "llogis") {
    sigma <- exp(log_scale)
    z <- (log(t) - eta) / sigma
    logS <- stats::plogis(-z, log = TRUE)
    logf <- stats::dlogis(z, log = TRUE) - log(t) - log(sigma)
  } else {
    stop("지원하지 않는 AFT family입니다: ", family)
  }
  
  logS <- pmin(logS, 0)
  S <- exp(pmax(logS, log(.Machine$double.xmin)))
  S <- pmin(pmax(S, 0), 1)
  
  list(logS = logS, S = S, logf = logf)
}

## 🟠 Build: start values for parametric cure optimization ===============================
build_parametric_starts <- function(data,
                                    Z,
                                    X,
                                    inc_term_labels_effective,
                                    lat_term_labels_effective,
                                    lat_formula_rhs,
                                    family) {
  k_gamma <- ncol(Z)
  k_beta <- ncol(X)
  has_scale <- family != "exp"
  
  gamma_base <- rep(0, k_gamma)
  if (k_gamma > 0L) {
    fit_glm0 <- try(stats::glm.fit(
      x = Z,
      y = data$event_transition,
      family = stats::quasibinomial(link = "logit")
    ), silent = TRUE)
    
    if (!inherits(fit_glm0, "try-error") && !is.null(fit_glm0$coefficients)) {
      gamma_base <- fit_glm0$coefficients
      gamma_base[!is.finite(gamma_base)] <- 0
      gamma_base <- as.numeric(gamma_base)
    } else if ("(Intercept)" %in% inc_term_labels_effective) {
      gamma_base[match("(Intercept)", inc_term_labels_effective)] <- stats::qlogis(min(max(mean(data$event_transition), 0.05), 0.95))
    }
  }
  
  beta_base <- rep(0, k_beta)
  log_scale_base <- 0
  
  dist_name <- switch(
    family,
    exp = "exponential",
    weibull = "weibull",
    lnorm = "lognormal",
    llogis = "loglogistic"
  )
  
  survreg_formula <- stats::as.formula(
    paste0("survival::Surv(time_years, event_transition) ~ ", lat_formula_rhs)
  )
  
  fit_survreg0 <- try(capture_warnings(
    survival::survreg(
      formula = survreg_formula,
      data = data,
      dist = dist_name
    )
  ), silent = TRUE)
  
  if (!inherits(fit_survreg0, "try-error")) {
    survreg_fit <- fit_survreg0$value
    co <- stats::coef(survreg_fit)
    tmp_beta <- rep(0, length(lat_term_labels_effective))
    names(tmp_beta) <- lat_term_labels_effective
    common_beta <- intersect(names(co), names(tmp_beta))
    tmp_beta[common_beta] <- co[common_beta]
    beta_base <- as.numeric(tmp_beta)
    
    if (has_scale) {
      sc <- survreg_fit$scale %||% 1
      if (is.finite(sc) && sc > 0) {
        log_scale_base <- log(sc)
      }
    }
  } else {
    unc_idx <- which(data$event_transition == 1L)
    if (length(unc_idx) >= max(5L, k_beta + 1L)) {
      fit_lm0 <- try(stats::lm(
        formula = stats::as.formula(paste0("log(time_years) ~ ", lat_formula_rhs)),
        data = data[unc_idx, , drop = FALSE]
      ), silent = TRUE)
      
      if (!inherits(fit_lm0, "try-error")) {
        co <- stats::coef(fit_lm0)
        tmp_beta <- rep(0, length(lat_term_labels_effective))
        names(tmp_beta) <- lat_term_labels_effective
        common_beta <- intersect(names(co), names(tmp_beta))
        tmp_beta[common_beta] <- co[common_beta]
        beta_base <- as.numeric(tmp_beta)
        
        if (has_scale) {
          log_scale_base <- log(max(stats::sigma(fit_lm0), 0.20))
        }
      }
    }
  }
  
  intercept_idx <- match("(Intercept)", inc_term_labels_effective)
  gamma_starts <- list(gamma_base)
  
  if (is.finite(intercept_idx)) {
    for (shift_val in parametric_gamma_intercept_shifts) {
      g <- gamma_base
      g[intercept_idx] <- g[intercept_idx] + shift_val
      gamma_starts[[length(gamma_starts) + 1L]] <- g
    }
    for (p0 in parametric_gamma_mean_uncured_starts) {
      g <- gamma_base
      g[intercept_idx] <- stats::qlogis(min(max(p0, prob_floor), 1 - prob_floor))
      gamma_starts[[length(gamma_starts) + 1L]] <- g
    }
  }
  
  gamma_keys <- vapply(gamma_starts, function(v) paste(round(v, 6), collapse = "|"), character(1))
  gamma_starts <- gamma_starts[!duplicated(gamma_keys)]
  
  starts <- list()
  for (g in gamma_starts) {
    base_par <- c(g, beta_base)
    starts[[length(starts) + 1L]] <- if (has_scale) c(base_par, log_scale_base) else base_par
    if (has_scale) {
      starts[[length(starts) + 1L]] <- c(base_par, log_scale_base + 0.5)
      starts[[length(starts) + 1L]] <- c(base_par, log_scale_base - 0.5)
    }
  }
  
  start_keys <- vapply(starts, function(v) paste(round(v, 6), collapse = "|"), character(1))
  starts[!duplicated(start_keys)]
}

## 🟠 Build: parametric AFT cure fitting engine ===============================
fit_parametric_mixture_cure <- function(data,
                                        incidence_bundle,
                                        latency_bundle,
                                        family) {
  time <- data$time_years
  event <- data$event_transition
  Z <- incidence_bundle$matrix
  X <- latency_bundle$matrix
  
  k_gamma <- ncol(Z)
  k_beta <- ncol(X)
  has_scale <- family != "exp"
  
  nll_fun <- function(par) {
    gamma <- par[seq_len(k_gamma)]
    beta <- par[k_gamma + seq_len(k_beta)]
    log_scale <- if (has_scale) par[k_gamma + k_beta + 1L] else NULL
    
    if (has_scale && (!is.finite(log_scale) || log_scale < -7 || log_scale > 7)) {
      return(big_nll)
    }
    if (any(!is.finite(par))) return(big_nll)
    
    eta_inc <- as.vector(Z %*% gamma)
    pi_uncured <- stats::plogis(eta_inc)
    pi_uncured <- pmin(pmax(pi_uncured, prob_floor), 1 - prob_floor)
    
    eta_lat <- if (k_beta > 0L) as.vector(X %*% beta) else rep(0, length(time))
    lat_obj <- aft_surv_density(
      family = family,
      time = time,
      eta = eta_lat,
      log_scale = log_scale,
      tiny = prob_floor
    )
    
    log_event <- log(pi_uncured) + lat_obj$logf
    cens_term <- pmax((1 - pi_uncured) + pi_uncured * lat_obj$S, prob_floor)
    log_cens <- log(cens_term)
    
    ll <- sum(event * log_event + (1 - event) * log_cens)
    if (!is.finite(ll)) return(big_nll)
    -ll
  }
  
  start_list <- build_parametric_starts(
    data = data,
    Z = Z,
    X = X,
    inc_term_labels_effective = incidence_bundle$map$term_label,
    lat_term_labels_effective = latency_bundle$map$term_label,
    lat_formula_rhs = latency_bundle$rhs,
    family = family
  )
  
  all_warnings <- character()
  best_fit <- NULL
  best_start_index <- NA_integer_
  n_successful_starts <- 0L
  
  fit_start_time <- Sys.time()
  
  for (s_idx in seq_along(start_list)) {
    st <- start_list[[s_idx]]
    
    fit_try <- try(capture_warnings(
      stats::optim(
        par = st,
        fn = nll_fun,
        method = "BFGS",
        control = list(maxit = parametric_maxit, reltol = parametric_reltol)
      )
    ), silent = TRUE)
    
    if (inherits(fit_try, "try-error")) {
      all_warnings <- unique(c(all_warnings, paste0("optim_error_start_", s_idx, ": ", as.character(fit_try))))
      next
    }
    
    all_warnings <- unique(c(all_warnings, fit_try$warnings))
    opt_obj <- fit_try$value
    if (is.finite(opt_obj$value)) {
      n_successful_starts <- n_successful_starts + 1L
    }
    
    if (is.null(best_fit) || (is.finite(opt_obj$value) && opt_obj$value < best_fit$value)) {
      best_fit <- opt_obj
      best_start_index <- s_idx
    }
  }
  
  fit_time_sec <- as.numeric(difftime(Sys.time(), fit_start_time, units = "secs"))
  
  if (is.null(best_fit) || !is.finite(best_fit$value)) {
    return(list(
      fit_ok = FALSE,
      fit_status = "error",
      warnings = unique(c(all_warnings, "모든 parametric 최적화 시도가 실패했습니다.")),
      fit_time_sec = fit_time_sec,
      n_start_attempts = length(start_list),
      n_successful_starts = n_successful_starts
    ))
  }
  
  gamma_hat <- best_fit$par[seq_len(k_gamma)]
  beta_hat <- best_fit$par[k_gamma + seq_len(k_beta)]
  log_scale_hat <- if (has_scale) best_fit$par[k_gamma + k_beta + 1L] else NA_real_
  
  hess_try <- try(capture_warnings(
    stats::optimHess(best_fit$par, fn = nll_fun)
  ), silent = TRUE)
  
  if (inherits(hess_try, "try-error")) {
    all_warnings <- unique(c(all_warnings, paste0("optimHess_error: ", as.character(hess_try))))
    hessian <- matrix(NA_real_, length(best_fit$par), length(best_fit$par))
  } else {
    all_warnings <- unique(c(all_warnings, hess_try$warnings))
    hessian <- hess_try$value
  }
  
  vcov_obj <- safe_inverse_hessian(hessian)
  se_vec <- if (!vcov_obj$singular) sqrt(diag(vcov_obj$vcov)) else rep(NA_real_, length(best_fit$par))
  
  gamma_se <- se_vec[seq_len(k_gamma)]
  beta_se <- se_vec[k_gamma + seq_len(k_beta)]
  log_scale_se <- if (has_scale) se_vec[k_gamma + k_beta + 1L] else NA_real_
  
  names(gamma_hat) <- incidence_bundle$map$term_label
  names(beta_hat) <- latency_bundle$map$term_label
  names(gamma_se) <- incidence_bundle$map$term_label
  names(beta_se) <- latency_bundle$map$term_label
  
  eta_inc_fit <- as.vector(Z %*% unname(gamma_hat))
  pi_uncured_fit <- stats::plogis(eta_inc_fit)
  pi_uncured_fit <- pmin(pmax(pi_uncured_fit, prob_floor), 1 - prob_floor)
  cure_fraction_fit <- 1 - pi_uncured_fit
  
  eta_lat_fit <- if (k_beta > 0L) as.vector(X %*% unname(beta_hat)) else rep(0, nrow(data))
  lat_fit <- aft_surv_density(
    family = family,
    time = time,
    eta = eta_lat_fit,
    log_scale = if (has_scale) log_scale_hat else NULL,
    tiny = prob_floor
  )
  
  log_event_fit <- log(pi_uncured_fit) + lat_fit$logf
  cens_term_fit <- pmax((1 - pi_uncured_fit) + pi_uncured_fit * lat_fit$S, prob_floor)
  log_cens_fit <- log(cens_term_fit)
  logLik_fit <- sum(event * log_event_fit + (1 - event) * log_cens_fit)
  
  n_par <- length(best_fit$par)
  AIC_fit <- 2 * n_par - 2 * logLik_fit
  BIC_fit <- log(nrow(data)) * n_par - 2 * logLik_fit
  
  boundary_flag <- any(pi_uncured_fit < boundary_prob_threshold) ||
    any(pi_uncured_fit > 1 - boundary_prob_threshold) ||
    any(cure_fraction_fit < boundary_prob_threshold) ||
    any(cure_fraction_fit > 1 - boundary_prob_threshold)
  
  list(
    fit_ok = TRUE,
    fit_status = "ok",
    lane = "AFT_param",
    family = family,
    warnings = unique(all_warnings),
    fit_time_sec = fit_time_sec,
    convergence_code = best_fit$convergence,
    convergence_message = best_fit$message %||% "",
    convergence_flag = best_fit$convergence == 0L,
    singular_hessian_flag = vcov_obj$singular,
    hessian_condition_number = vcov_obj$condition_number,
    hessian = hessian,
    vcov = vcov_obj$vcov,
    vcov_method = if (vcov_obj$singular) "not_available_due_to_hessian" else "inverse_hessian_full",
    incidence_target = "uncured_probability",
    gamma = gamma_hat,
    gamma_se = gamma_se,
    beta = beta_hat,
    beta_se = beta_se,
    log_scale = log_scale_hat,
    log_scale_se = log_scale_se,
    scale_natural = if (has_scale) exp(log_scale_hat) else NA_real_,
    scale_se_natural = if (has_scale && is.finite(log_scale_se)) exp(log_scale_hat) * log_scale_se else NA_real_,
    logLik = logLik_fit,
    AIC = AIC_fit,
    BIC = BIC_fit,
    n_parameters = n_par,
    pi_uncured_train = pi_uncured_fit,
    cure_fraction_train = cure_fraction_fit,
    boundary_flag = boundary_flag,
    max_event_time_years = if (any(event == 1L)) max(time[event == 1L]) else NA_real_,
    max_followup_time_years = max(time),
    incidence_bundle = incidence_bundle,
    latency_bundle = latency_bundle,
    has_scale = has_scale,
    n_start_attempts = length(start_list),
    n_successful_starts = n_successful_starts,
    best_start_index = best_start_index,
    best_objective_value = best_fit$value
  )
}

## 🟠 Build: parametric AFT cure prediction engine ===============================
predict_parametric_cure_fit <- function(fit_obj, Znew, Xnew, years) {
  pi_uncured <- stats::plogis(as.vector(Znew %*% unname(fit_obj$gamma)))
  pi_uncured <- pmin(pmax(pi_uncured, prob_floor), 1 - prob_floor)
  
  eta_lat <- if (ncol(Xnew) > 0L) as.vector(Xnew %*% unname(fit_obj$beta)) else rep(0, nrow(Znew))
  
  S_u <- matrix(NA_real_, nrow = nrow(Znew), ncol = length(years))
  S_pop <- matrix(NA_real_, nrow = nrow(Znew), ncol = length(years))
  
  for (j in seq_along(years)) {
    lat_obj <- aft_surv_density(
      family = fit_obj$family,
      time = rep(years[j], nrow(Znew)),
      eta = eta_lat,
      log_scale = if (fit_obj$has_scale) fit_obj$log_scale else NULL,
      tiny = prob_floor
    )
    S_u[, j] <- lat_obj$S
    S_pop[, j] <- (1 - pi_uncured) + pi_uncured * lat_obj$S
  }
  
  colnames(S_u) <- paste0("y", years)
  colnames(S_pop) <- paste0("y", years)
  
  list(
    years = years,
    pi_uncured = pi_uncured,
    cure_fraction = 1 - pi_uncured,
    S_u = S_u,
    S_pop = S_pop
  )
}

## 🟠 Build: fractional logistic and Cox EM subroutines ===============================
fractional_logit_nll <- function(gamma, Z, w) {
  eta <- as.vector(Z %*% gamma)
  p <- stats::plogis(eta)
  p <- pmin(pmax(p, prob_floor), 1 - prob_floor)
  -sum(w * log(p) + (1 - w) * log1p(-p))
}

fit_fractional_logit_mstep <- function(Z, w, start_gamma) {
  fit_try <- try(capture_warnings(
    stats::optim(
      par = start_gamma,
      fn = fractional_logit_nll,
      Z = Z,
      w = w,
      method = "BFGS",
      control = list(maxit = 500L, reltol = 1e-10)
    )
  ), silent = TRUE)
  
  if (inherits(fit_try, "try-error")) {
    list(
      par = start_gamma,
      warnings = paste0("fractional_logit_optim_error: ", as.character(fit_try)),
      convergence = 1L
    )
  } else {
    list(
      par = fit_try$value$par,
      warnings = fit_try$warnings,
      convergence = fit_try$value$convergence
    )
  }
}

fit_weighted_cox_mstep <- function(time, event, X, w, beta_start = NULL) {
  df <- as.data.frame(X)
  df$time <- time
  df$event <- event
  df$w <- w
  
  rhs <- if (ncol(X) > 0L) {
    paste(colnames(X), collapse = " + ")
  } else {
    "1"
  }
  
  fml <- stats::as.formula(paste0("survival::Surv(time, event) ~ ", rhs))
  
  cox_try <- try(capture_warnings(
    survival::coxph(
      formula = fml,
      data = df,
      ties = "breslow",
      singular.ok = TRUE,
      init = beta_start,
      weights = w,
      control = survival::coxph.control(iter.max = coxph_iter_max, eps = 1e-09),
      x = FALSE,
      y = FALSE,
      model = FALSE
    )
  ), silent = TRUE)
  
  if (inherits(cox_try, "try-error")) {
    return(list(
      ok = FALSE,
      warnings = paste0("coxph_error: ", as.character(cox_try))
    ))
  }
  
  cox_fit <- cox_try$value
  beta_hat <- stats::coef(cox_fit)
  beta_hat <- if (length(beta_hat) == 0L) numeric(0) else beta_hat
  if (anyNA(beta_hat)) beta_hat[is.na(beta_hat)] <- 0
  
  bh_try <- try(capture_warnings(
    survival::basehaz(cox_fit, centered = FALSE)
  ), silent = TRUE)
  
  if (inherits(bh_try, "try-error")) {
    return(list(
      ok = FALSE,
      warnings = c(cox_try$warnings, paste0("basehaz_error: ", as.character(bh_try)))
    ))
  }
  
  vc_try <- try(stats::vcov(cox_fit), silent = TRUE)
  beta_vcov <- if (inherits(vc_try, "try-error")) {
    matrix(NA_real_, length(beta_hat), length(beta_hat))
  } else {
    vc_try
  }
  
  list(
    ok = TRUE,
    beta = beta_hat,
    beta_vcov = beta_vcov,
    cox_fit = cox_fit,
    basehaz_df = bh_try$value,
    warnings = unique(c(cox_try$warnings, bh_try$warnings))
  )
}

## 🟠 Build: semiparametric Cox-latency cure EM engine ===============================
fit_ph_cure_em <- function(data,
                           incidence_bundle,
                           latency_bundle) {
  time <- data$time_years
  event <- data$event_transition
  Z <- incidence_bundle$matrix
  X <- latency_bundle$matrix
  
  n <- nrow(data)
  if (sum(event) == 0L) {
    return(list(
      fit_ok = FALSE,
      fit_status = "error",
      warnings = "PH cure EM 적합 불가: transition event가 0개입니다.",
      fit_time_sec = 0
    ))
  }
  
  all_warnings <- character()
  fit_start_time <- Sys.time()
  
  init_glm <- try(stats::glm.fit(
    x = Z,
    y = event,
    family = stats::quasibinomial(link = "logit")
  ), silent = TRUE)
  
  gamma_curr <- if (!inherits(init_glm, "try-error") && !is.null(init_glm$coefficients)) {
    g <- init_glm$coefficients
    g[!is.finite(g)] <- 0
    as.numeric(g)
  } else {
    g <- rep(0, ncol(Z))
    if ("(Intercept)" %in% incidence_bundle$map$term_label) {
      g[match("(Intercept)", incidence_bundle$map$term_label)] <- stats::qlogis(min(max(mean(event), 0.05), 0.95))
    }
    g
  }
  
  init_cox <- fit_weighted_cox_mstep(
    time = time,
    event = event,
    X = X,
    w = rep(1, n),
    beta_start = NULL
  )
  
  if (!isTRUE(init_cox$ok)) {
    return(list(
      fit_ok = FALSE,
      fit_status = "error",
      warnings = unique(c(all_warnings, init_cox$warnings, "초기 Cox PH 적합 실패")),
      fit_time_sec = as.numeric(difftime(Sys.time(), fit_start_time, units = "secs"))
    ))
  }
  
  beta_curr <- as.numeric(init_cox$beta)
  names(beta_curr) <- latency_bundle$map$term_label
  basehaz_df_curr <- init_cox$basehaz_df
  cox_fit_curr <- init_cox$cox_fit
  beta_vcov_curr <- init_cox$beta_vcov
  all_warnings <- unique(c(all_warnings, init_cox$warnings))
  
  converged <- FALSE
  em_trace <- numeric(em_max_iter)
  pseudo_loglik_trace <- numeric(em_max_iter)
  final_w <- rep(1, n)
  
  for (iter in seq_len(em_max_iter)) {
    S0_fun <- make_step_surv_fun(basehaz_df_curr, zero_tail = TRUE)
    eta_inc <- as.vector(Z %*% gamma_curr)
    pi_uncured <- stats::plogis(eta_inc)
    pi_uncured <- pmin(pmax(pi_uncured, prob_floor), 1 - prob_floor)
    
    lp <- if (ncol(X) > 0L) as.vector(X %*% beta_curr) else rep(0, n)
    S0_at_t <- S0_fun(time)
    S_u <- S0_at_t ^ exp(lp)
    
    denom <- pmax((1 - pi_uncured) + pi_uncured * S_u, prob_floor)
    w <- event + (1 - event) * pi_uncured * S_u / denom
    w <- pmin(pmax(w, prob_floor), 1)
    final_w <- w
    
    gamma_update <- fit_fractional_logit_mstep(
      Z = Z,
      w = w,
      start_gamma = gamma_curr
    )
    gamma_new <- as.numeric(gamma_update$par)
    all_warnings <- unique(c(all_warnings, gamma_update$warnings))
    
    cox_update <- fit_weighted_cox_mstep(
      time = time,
      event = event,
      X = X,
      w = w,
      beta_start = beta_curr
    )
    
    if (!isTRUE(cox_update$ok)) {
      return(list(
        fit_ok = FALSE,
        fit_status = "error",
        warnings = unique(c(all_warnings, cox_update$warnings, paste0("EM iteration ", iter, " 에서 Cox PH 적합 실패"))),
        fit_time_sec = as.numeric(difftime(Sys.time(), fit_start_time, units = "secs"))
      ))
    }
    
    beta_new <- if (length(cox_update$beta) > 0L) as.numeric(cox_update$beta) else numeric(0)
    basehaz_df_new <- cox_update$basehaz_df
    cox_fit_new <- cox_update$cox_fit
    beta_vcov_new <- cox_update$beta_vcov
    all_warnings <- unique(c(all_warnings, cox_update$warnings))
    
    incidence_ll <- -fractional_logit_nll(gamma_new, Z, w)
    partial_ll <- cox_fit_new$loglik[2] %||% NA_real_
    pseudo_loglik_trace[iter] <- incidence_ll + partial_ll
    
    old_par <- c(gamma_curr, beta_curr)
    new_par <- c(gamma_new, beta_new)
    rel_change <- max(abs(new_par - old_par) / pmax(abs(old_par), 1e-6))
    em_trace[iter] <- rel_change
    
    gamma_curr <- gamma_new
    beta_curr <- beta_new
    basehaz_df_curr <- basehaz_df_new
    cox_fit_curr <- cox_fit_new
    beta_vcov_curr <- beta_vcov_new
    
    if (is.finite(rel_change) && rel_change < em_tol) {
      converged <- TRUE
      em_trace <- em_trace[seq_len(iter)]
      pseudo_loglik_trace <- pseudo_loglik_trace[seq_len(iter)]
      break
    }
    
    if (iter == em_max_iter) {
      em_trace <- em_trace[seq_len(iter)]
      pseudo_loglik_trace <- pseudo_loglik_trace[seq_len(iter)]
    }
  }
  
  gamma_hess_try <- try(capture_warnings(
    stats::optimHess(
      par = gamma_curr,
      fn = fractional_logit_nll,
      Z = Z,
      w = final_w
    )
  ), silent = TRUE)
  
  if (inherits(gamma_hess_try, "try-error")) {
    gamma_hessian <- matrix(NA_real_, length(gamma_curr), length(gamma_curr))
    all_warnings <- unique(c(all_warnings, paste0("gamma_optimHess_error: ", as.character(gamma_hess_try))))
  } else {
    gamma_hessian <- gamma_hess_try$value
    all_warnings <- unique(c(all_warnings, gamma_hess_try$warnings))
  }
  
  gamma_vcov_obj <- safe_inverse_hessian(gamma_hessian)
  gamma_vcov <- gamma_vcov_obj$vcov
  gamma_se <- if (!gamma_vcov_obj$singular) sqrt(diag(gamma_vcov)) else rep(NA_real_, length(gamma_curr))
  
  beta_se <- tryCatch(sqrt(diag(beta_vcov_curr)), error = function(e) rep(NA_real_, length(beta_curr)))
  if (length(beta_se) == 0L) beta_se <- numeric(0)
  
  names(gamma_curr) <- incidence_bundle$map$term_label
  names(beta_curr) <- latency_bundle$map$term_label
  names(gamma_se) <- incidence_bundle$map$term_label
  names(beta_se) <- latency_bundle$map$term_label
  
  S0_fun_final <- make_step_surv_fun(basehaz_df_curr, zero_tail = TRUE)
  eta_inc_final <- as.vector(Z %*% unname(gamma_curr))
  pi_uncured_final <- stats::plogis(eta_inc_final)
  pi_uncured_final <- pmin(pmax(pi_uncured_final, prob_floor), 1 - prob_floor)
  cure_fraction_final <- 1 - pi_uncured_final
  
  lp_final <- if (ncol(X) > 0L) as.vector(X %*% unname(beta_curr)) else rep(0, n)
  S_u_final <- S0_fun_final(time) ^ exp(lp_final)
  
  boundary_flag <- any(pi_uncured_final < boundary_prob_threshold) ||
    any(pi_uncured_final > 1 - boundary_prob_threshold) ||
    any(cure_fraction_final < boundary_prob_threshold) ||
    any(cure_fraction_final > 1 - boundary_prob_threshold)
  
  fit_time_sec <- as.numeric(difftime(Sys.time(), fit_start_time, units = "secs"))
  
  list(
    fit_ok = TRUE,
    fit_status = if (converged) "ok" else "maxit_reached",
    lane = "PH_semiparam",
    family = "coxph",
    warnings = unique(all_warnings),
    fit_time_sec = fit_time_sec,
    convergence_code = if (converged) 0L else 1L,
    convergence_message = if (converged) "EM converged" else "EM max iteration reached",
    convergence_flag = converged,
    singular_hessian_flag = gamma_vcov_obj$singular,
    hessian_condition_number = gamma_vcov_obj$condition_number,
    incidence_target = "uncured_probability",
    gamma = gamma_curr,
    gamma_se = gamma_se,
    gamma_vcov = gamma_vcov,
    beta = beta_curr,
    beta_se = beta_se,
    beta_vcov = beta_vcov_curr,
    log_scale = NA_real_,
    log_scale_se = NA_real_,
    scale_natural = NA_real_,
    scale_se_natural = NA_real_,
    logLik = NA_real_,
    partial_logLik = cox_fit_curr$loglik[2] %||% NA_real_,
    pseudo_logLik = tail(pseudo_loglik_trace, 1),
    AIC = NA_real_,
    BIC = NA_real_,
    n_parameters = length(gamma_curr) + length(beta_curr),
    pi_uncured_train = pi_uncured_final,
    cure_fraction_train = cure_fraction_final,
    boundary_flag = boundary_flag,
    max_event_time_years = if (any(event == 1L)) max(time[event == 1L]) else NA_real_,
    max_followup_time_years = max(time),
    incidence_bundle = incidence_bundle,
    latency_bundle = latency_bundle,
    basehaz_df = basehaz_df_curr,
    baseline_zero_tail = TRUE,
    em_iterations = length(em_trace),
    em_max_relative_change = if (length(em_trace) > 0L) tail(em_trace, 1) else NA_real_,
    em_trace = em_trace,
    pseudo_logLik_trace = pseudo_loglik_trace,
    final_posterior_uncured_weights = final_w,
    vcov_method = "conditional_block_diagonal"
  )
}

## 🟠 Build: semiparametric Cox-latency cure prediction engine ===============================
predict_ph_cure_fit <- function(fit_obj, Znew, Xnew, years) {
  pi_uncured <- stats::plogis(as.vector(Znew %*% unname(fit_obj$gamma)))
  pi_uncured <- pmin(pmax(pi_uncured, prob_floor), 1 - prob_floor)
  
  lp <- if (ncol(Xnew) > 0L) as.vector(Xnew %*% unname(fit_obj$beta)) else rep(0, nrow(Znew))
  S0_fun <- make_step_surv_fun(fit_obj$basehaz_df, zero_tail = isTRUE(fit_obj$baseline_zero_tail))
  
  S_u <- matrix(NA_real_, nrow = nrow(Znew), ncol = length(years))
  S_pop <- matrix(NA_real_, nrow = nrow(Znew), ncol = length(years))
  
  for (j in seq_along(years)) {
    S0_t <- S0_fun(rep(years[j], nrow(Znew)))
    S_u[, j] <- S0_t ^ exp(lp)
    S_pop[, j] <- (1 - pi_uncured) + pi_uncured * S_u[, j]
  }
  
  colnames(S_u) <- paste0("y", years)
  colnames(S_pop) <- paste0("y", years)
  
  list(
    years = years,
    pi_uncured = pi_uncured,
    cure_fraction = 1 - pi_uncured,
    S_u = S_u,
    S_pop = S_pop
  )
}

## 🟠 Build: IPCW performance metrics and IBS calculators ===============================
weighted_auc_ipcw <- function(r_case, r_ctrl, w_case, w_ctrl) {
  if (length(r_case) == 0L || length(r_ctrl) == 0L) return(NA_real_)
  if (sum(w_case) <= 0 || sum(w_ctrl) <= 0) return(NA_real_)
  
  ord_ctrl <- order(r_ctrl)
  r_ctrl <- r_ctrl[ord_ctrl]
  w_ctrl <- w_ctrl[ord_ctrl]
  
  ctrl_sum_by_r <- tapply(w_ctrl, r_ctrl, sum)
  u <- as.numeric(names(ctrl_sum_by_r))
  ord_u <- order(u)
  u <- u[ord_u]
  ctrl_sum <- as.numeric(ctrl_sum_by_r[ord_u])
  cum_ctrl <- cumsum(ctrl_sum)
  
  less_idx <- findInterval(r_case - 1e-12, u)
  less_w <- ifelse(less_idx > 0L, cum_ctrl[less_idx], 0)
  
  match_idx <- match(r_case, u)
  equal_w <- ifelse(!is.na(match_idx), ctrl_sum[match_idx], 0)
  
  numerator <- sum(w_case * (less_w + 0.5 * equal_w))
  denominator <- sum(w_case) * sum(w_ctrl)
  
  if (denominator <= 0) return(NA_real_)
  numerator / denominator
}

compute_time_dependent_metrics <- function(data, pred_obj, ipcw_helper, years, meta) {
  rows <- vector("list", length(years))
  
  for (j in seq_along(years)) {
    t0 <- years[j]
    surv_pred <- pred_obj$S_pop[, j]
    risk_pred <- 1 - surv_pred
    
    cases <- data$event_transition == 1L & data$time_years <= t0
    controls <- data$time_years > t0
    
    G_t <- ipcw_helper$G_at(t0)
    G_t_used <- pmax(G_t, prob_floor)
    
    weights_brier <- rep(0, nrow(data))
    if (any(cases)) {
      G_left_case <- ipcw_helper$G_left(data$time_years[cases])
      weights_brier[cases] <- 1 / pmax(G_left_case, prob_floor)
    }
    if (any(controls)) {
      weights_brier[controls] <- 1 / G_t_used
    }
    
    outcome_t <- as.integer(cases)
    brier_ipcw <- if (sum(weights_brier) > 0 && is.finite(G_t)) {
      mean(weights_brier * (outcome_t - risk_pred)^2)
    } else {
      NA_real_
    }
    
    auc_cd <- if (any(cases) && any(controls) && is.finite(G_t) && G_t > prob_floor) {
      w_case <- 1 / pmax(ipcw_helper$G_left(data$time_years[cases]), prob_floor)
      w_ctrl <- rep(1 / G_t_used, sum(controls))
      weighted_auc_ipcw(
        r_case = risk_pred[cases],
        r_ctrl = risk_pred[controls],
        w_case = w_case,
        w_ctrl = w_ctrl
      )
    } else {
      NA_real_
    }
    
    rows[[j]] <- data.frame(
      dataset_branch = meta$dataset_branch,
      data_cut = "full",
      subgroup = "overall",
      model_class = "cure_MLE",
      lane = meta$lane,
      family = meta$family_label,
      family_code = meta$family_code,
      spec_id = meta$spec_id,
      model_id = meta$model_id,
      metric_scope = "apparent_training",
      year = t0,
      auc_cd = auc_cd,
      brier_ipcw = brier_ipcw,
      G_t_censoring_survival = G_t,
      n_cases = sum(cases),
      n_controls = sum(controls),
      auc_valid = is.finite(auc_cd),
      brier_valid = is.finite(brier_ipcw),
      beyond_max_followup = if (is.finite(max(data$time_years))) t0 > max(data$time_years) else NA,
      beyond_last_event = if (any(data$event_transition == 1L)) t0 > max(data$time_years[data$event_transition == 1L]) else NA,
      stringsAsFactors = FALSE
    )
  }
  
  perf_df <- do.call(rbind, rows)
  valid_brier <- is.finite(perf_df$brier_ipcw)
  ibs_val <- if (sum(valid_brier) >= 2L) {
    trapz_base(perf_df$year[valid_brier], perf_df$brier_ipcw[valid_brier]) /
      (max(perf_df$year[valid_brier]) - min(perf_df$year[valid_brier]))
  } else {
    NA_real_
  }
  perf_df$IBS_apparent <- ibs_val
  perf_df
}

## 🟠 Build: output assemblers for predictions and summaries ===============================
assemble_subject_prediction_df <- function(data, pred_obj, meta, last_event_time, max_followup_time) {
  out <- vector("list", length(pred_obj$years))
  
  for (j in seq_along(pred_obj$years)) {
    yr <- pred_obj$years[j]
    out[[j]] <- data.frame(
      dataset_branch = meta$dataset_branch,
      data_cut = "full",
      subgroup = "overall",
      model_class = "cure_MLE",
      lane = meta$lane,
      family = meta$family_label,
      family_code = meta$family_code,
      spec_id = meta$spec_id,
      model_id = meta$model_id,
      id_original = data[[id_var]],
      unique_id = data$unique_id,
      site_std = data$site_std,
      site_num = data$site_num,
      sex_num = data$sex_num,
      age_raw = data$age_raw,
      age_s = data$age_s,
      time_days_original = data$time_days_original,
      time_years = data$time_years,
      status_num_observed = data[[status_var]],
      event_transition = data$event_transition,
      year = yr,
      surv = pred_obj$S_pop[, j],
      risk = 1 - pred_obj$S_pop[, j],
      surv_uncured = pred_obj$S_u[, j],
      pi_uncured = pred_obj$pi_uncured,
      cure_fraction = pred_obj$cure_fraction,
      beyond_last_event = if (is.finite(last_event_time)) yr > last_event_time else NA,
      beyond_max_followup = if (is.finite(max_followup_time)) yr > max_followup_time else NA,
      is_extrapolated = if (is.finite(max_followup_time)) yr > max_followup_time else NA,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, out)
}

assemble_mean_prediction_df <- function(pred_subject_df) {
  years <- sort(unique(pred_subject_df$year))
  out <- vector("list", length(years))
  
  for (j in seq_along(years)) {
    dd <- pred_subject_df[pred_subject_df$year == years[j], , drop = FALSE]
    out[[j]] <- data.frame(
      dataset_branch = dd$dataset_branch[1],
      data_cut = dd$data_cut[1],
      subgroup = dd$subgroup[1],
      model_class = dd$model_class[1],
      lane = dd$lane[1],
      family = dd$family[1],
      family_code = dd$family_code[1],
      spec_id = dd$spec_id[1],
      model_id = dd$model_id[1],
      year = years[j],
      mean_surv = mean(dd$surv),
      sd_surv = stats::sd(dd$surv),
      surv_q25 = stats::quantile(dd$surv, 0.25),
      surv_median = stats::median(dd$surv),
      surv_q75 = stats::quantile(dd$surv, 0.75),
      mean_risk = mean(dd$risk),
      sd_risk = stats::sd(dd$risk),
      risk_q25 = stats::quantile(dd$risk, 0.25),
      risk_median = stats::median(dd$risk),
      risk_q75 = stats::quantile(dd$risk, 0.75),
      mean_surv_uncured = mean(dd$surv_uncured),
      mean_pi_uncured = mean(dd$pi_uncured),
      mean_cure_fraction = mean(dd$cure_fraction),
      median_cure_fraction = stats::median(dd$cure_fraction),
      n_subjects = nrow(dd),
      beyond_last_event = unique(dd$beyond_last_event)[1],
      beyond_max_followup = unique(dd$beyond_max_followup)[1],
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, out)
}

assemble_cure_summary_row <- function(fit_obj, mean_pred_df, perf_df, meta) {
  row <- data.frame(
    dataset_branch = meta$dataset_branch,
    data_cut = "full",
    subgroup = "overall",
    model_class = "cure_MLE",
    lane = meta$lane,
    family = meta$family_label,
    family_code = meta$family_code,
    spec_id = meta$spec_id,
    model_id = meta$model_id,
    incidence_target = fit_obj$incidence_target,
    mean_pi_uncured = mean(fit_obj$pi_uncured_train),
    mean_cure_fraction = mean(fit_obj$cure_fraction_train),
    median_cure_fraction = stats::median(fit_obj$cure_fraction_train),
    cure_fraction_q25 = stats::quantile(fit_obj$cure_fraction_train, 0.25),
    cure_fraction_q75 = stats::quantile(fit_obj$cure_fraction_train, 0.75),
    min_cure_fraction = min(fit_obj$cure_fraction_train),
    max_cure_fraction = max(fit_obj$cure_fraction_train),
    boundary_flag = fit_obj$boundary_flag,
    convergence_flag = fit_obj$convergence_flag,
    IBS_apparent = if (all(is.na(perf_df$IBS_apparent))) NA_real_ else unique(stats::na.omit(perf_df$IBS_apparent))[1],
    mean_auc_valid = if (any(is.finite(perf_df$auc_cd))) mean(perf_df$auc_cd, na.rm = TRUE) else NA_real_,
    mean_brier_valid = if (any(is.finite(perf_df$brier_ipcw))) mean(perf_df$brier_ipcw, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (yr in mean_pred_df$year) {
    row[[paste0("mean_surv_y", yr)]] <- mean_pred_df$mean_surv[mean_pred_df$year == yr]
    row[[paste0("mean_risk_y", yr)]] <- mean_pred_df$mean_risk[mean_pred_df$year == yr]
    row[[paste0("mean_su_y", yr)]] <- mean_pred_df$mean_surv_uncured[mean_pred_df$year == yr]
  }
  
  row
}

assemble_param_coef_df <- function(fit_obj, meta) {
  rows <- list()
  
  if (length(fit_obj$gamma) > 0L) {
    gamma_terms <- names(fit_obj$gamma)
    gamma_est <- unname(fit_obj$gamma)
    gamma_se <- unname(fit_obj$gamma_se)
    gamma_z <- gamma_est / gamma_se
    gamma_p <- norm_p_value(gamma_z)
    gamma_ci <- wald_ci(gamma_est, gamma_se)
    
    rows[[length(rows) + 1L]] <- data.frame(
      dataset_branch = meta$dataset_branch,
      data_cut = "full",
      subgroup = "overall",
      model_class = "cure_MLE",
      lane = meta$lane,
      family = meta$family_label,
      family_code = meta$family_code,
      spec_id = meta$spec_id,
      model_id = meta$model_id,
      component = "incidence",
      coefficient_scale = "logit_uncured_odds",
      exp_scale_label = "odds_ratio_uncured",
      term = gamma_terms,
      estimate = gamma_est,
      std_error = gamma_se,
      z_value = gamma_z,
      p_value = gamma_p,
      wald_lcl = gamma_ci[, "lcl"],
      wald_ucl = gamma_ci[, "ucl"],
      exp_estimate = exp(gamma_est),
      exp_wald_lcl = exp(gamma_ci[, "lcl"]),
      exp_wald_ucl = exp(gamma_ci[, "ucl"]),
      se_method = if (meta$lane == "AFT_param") "inverse_hessian" else "conditional_fractional_logit",
      stringsAsFactors = FALSE
    )
  }
  
  if (length(fit_obj$beta) > 0L) {
    beta_terms <- names(fit_obj$beta)
    beta_est <- unname(fit_obj$beta)
    beta_se <- unname(fit_obj$beta_se)
    beta_z <- beta_est / beta_se
    beta_p <- norm_p_value(beta_z)
    beta_ci <- wald_ci(beta_est, beta_se)
    
    rows[[length(rows) + 1L]] <- data.frame(
      dataset_branch = meta$dataset_branch,
      data_cut = "full",
      subgroup = "overall",
      model_class = "cure_MLE",
      lane = meta$lane,
      family = meta$family_label,
      family_code = meta$family_code,
      spec_id = meta$spec_id,
      model_id = meta$model_id,
      component = "latency",
      coefficient_scale = if (meta$lane == "AFT_param") "log_time_ratio" else "log_hazard_ratio",
      exp_scale_label = if (meta$lane == "AFT_param") "time_ratio" else "hazard_ratio",
      term = beta_terms,
      estimate = beta_est,
      std_error = beta_se,
      z_value = beta_z,
      p_value = beta_p,
      wald_lcl = beta_ci[, "lcl"],
      wald_ucl = beta_ci[, "ucl"],
      exp_estimate = exp(beta_est),
      exp_wald_lcl = exp(beta_ci[, "lcl"]),
      exp_wald_ucl = exp(beta_ci[, "ucl"]),
      se_method = if (meta$lane == "AFT_param") "inverse_hessian" else "conditional_weighted_coxph",
      stringsAsFactors = FALSE
    )
  }
  
  if (meta$lane == "AFT_param" && isTRUE(fit_obj$has_scale)) {
    sc_ci <- wald_ci(fit_obj$log_scale, fit_obj$log_scale_se)
    rows[[length(rows) + 1L]] <- data.frame(
      dataset_branch = meta$dataset_branch,
      data_cut = "full",
      subgroup = "overall",
      model_class = "cure_MLE",
      lane = meta$lane,
      family = meta$family_label,
      family_code = meta$family_code,
      spec_id = meta$spec_id,
      model_id = meta$model_id,
      component = "ancillary",
      coefficient_scale = "log_scale_parameter",
      exp_scale_label = "scale_parameter",
      term = "log_scale",
      estimate = fit_obj$log_scale,
      std_error = fit_obj$log_scale_se,
      z_value = fit_obj$log_scale / fit_obj$log_scale_se,
      p_value = norm_p_value(fit_obj$log_scale / fit_obj$log_scale_se),
      wald_lcl = sc_ci[, "lcl"],
      wald_ucl = sc_ci[, "ucl"],
      exp_estimate = fit_obj$scale_natural,
      exp_wald_lcl = exp(sc_ci[, "lcl"]),
      exp_wald_ucl = exp(sc_ci[, "ucl"]),
      se_method = "inverse_hessian",
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, rows)
}

assemble_vcov_long_df <- function(fit_obj, meta) {
  rows <- list()
  
  if (meta$lane == "AFT_param" && !is.null(fit_obj$vcov) && all(dim(fit_obj$vcov) > 0L)) {
    param_labels <- c(
      paste0("incidence::", names(fit_obj$gamma)),
      paste0("latency::", names(fit_obj$beta)),
      if (isTRUE(fit_obj$has_scale)) "ancillary::log_scale" else character(0)
    )
    
    vc <- fit_obj$vcov
    if (length(param_labels) == nrow(vc)) {
      for (r in seq_along(param_labels)) {
        for (c in seq_len(r)) {
          rows[[length(rows) + 1L]] <- data.frame(
            dataset_branch = meta$dataset_branch,
            data_cut = "full",
            subgroup = "overall",
            model_class = "cure_MLE",
            lane = meta$lane,
            family = meta$family_label,
            family_code = meta$family_code,
            spec_id = meta$spec_id,
            model_id = meta$model_id,
            param_row = param_labels[r],
            param_col = param_labels[c],
            covariance = vc[r, c],
            vcov_method = fit_obj$vcov_method,
            cross_block_available = TRUE,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  if (meta$lane == "PH_semiparam") {
    if (!is.null(fit_obj$gamma_vcov) && all(dim(fit_obj$gamma_vcov) > 0L)) {
      g_labels <- paste0("incidence::", names(fit_obj$gamma))
      for (r in seq_along(g_labels)) {
        for (c in seq_len(r)) {
          rows[[length(rows) + 1L]] <- data.frame(
            dataset_branch = meta$dataset_branch,
            data_cut = "full",
            subgroup = "overall",
            model_class = "cure_MLE",
            lane = meta$lane,
            family = meta$family_label,
            family_code = meta$family_code,
            spec_id = meta$spec_id,
            model_id = meta$model_id,
            param_row = g_labels[r],
            param_col = g_labels[c],
            covariance = fit_obj$gamma_vcov[r, c],
            vcov_method = "conditional_block_diagonal_incidence",
            cross_block_available = FALSE,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    if (!is.null(fit_obj$beta_vcov) && all(dim(fit_obj$beta_vcov) > 0L) && length(fit_obj$beta) > 0L) {
      b_labels <- paste0("latency::", names(fit_obj$beta))
      for (r in seq_along(b_labels)) {
        for (c in seq_len(r)) {
          rows[[length(rows) + 1L]] <- data.frame(
            dataset_branch = meta$dataset_branch,
            data_cut = "full",
            subgroup = "overall",
            model_class = "cure_MLE",
            lane = meta$lane,
            family = meta$family_label,
            family_code = meta$family_code,
            spec_id = meta$spec_id,
            model_id = meta$model_id,
            param_row = b_labels[r],
            param_col = b_labels[c],
            covariance = fit_obj$beta_vcov[r, c],
            vcov_method = "conditional_block_diagonal_latency",
            cross_block_available = FALSE,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  do.call(rbind, rows)
}

assemble_cox_baseline_df <- function(fit_obj, meta) {
  bh <- fit_obj$basehaz_df
  if (is.null(bh) || nrow(bh) == 0L) return(NULL)
  
  data.frame(
    dataset_branch = meta$dataset_branch,
    data_cut = "full",
    subgroup = "overall",
    model_class = "cure_MLE",
    lane = meta$lane,
    family = meta$family_label,
    family_code = meta$family_code,
    spec_id = meta$spec_id,
    model_id = meta$model_id,
    time_years = bh$time,
    baseline_cumhaz = bh$hazard,
    baseline_surv = exp(-bh$hazard),
    zero_tail_applied = isTRUE(fit_obj$baseline_zero_tail),
    stringsAsFactors = FALSE
  )
}

assemble_fit_diagnostics_df <- function(fit_obj, meta) {
  rows <- list()
  
  rows[[length(rows) + 1L]] <- data.frame(
    dataset_branch = meta$dataset_branch,
    data_cut = "full",
    subgroup = "overall",
    model_class = "cure_MLE",
    lane = meta$lane,
    family = meta$family_label,
    family_code = meta$family_code,
    spec_id = meta$spec_id,
    model_id = meta$model_id,
    diagnostic_group = "summary",
    iteration = NA_integer_,
    metric_name = "fit_time_sec",
    metric_value = fit_obj$fit_time_sec %||% NA_real_,
    metric_text = NA_character_,
    stringsAsFactors = FALSE
  )
  
  if (meta$lane == "AFT_param") {
    rows[[length(rows) + 1L]] <- data.frame(
      dataset_branch = meta$dataset_branch,
      data_cut = "full",
      subgroup = "overall",
      model_class = "cure_MLE",
      lane = meta$lane,
      family = meta$family_label,
      family_code = meta$family_code,
      spec_id = meta$spec_id,
      model_id = meta$model_id,
      diagnostic_group = "summary",
      iteration = NA_integer_,
      metric_name = "n_start_attempts",
      metric_value = fit_obj$n_start_attempts %||% NA_real_,
      metric_text = NA_character_,
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1L]] <- data.frame(
      dataset_branch = meta$dataset_branch,
      data_cut = "full",
      subgroup = "overall",
      model_class = "cure_MLE",
      lane = meta$lane,
      family = meta$family_label,
      family_code = meta$family_code,
      spec_id = meta$spec_id,
      model_id = meta$model_id,
      diagnostic_group = "summary",
      iteration = NA_integer_,
      metric_name = "n_successful_starts",
      metric_value = fit_obj$n_successful_starts %||% NA_real_,
      metric_text = NA_character_,
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1L]] <- data.frame(
      dataset_branch = meta$dataset_branch,
      data_cut = "full",
      subgroup = "overall",
      model_class = "cure_MLE",
      lane = meta$lane,
      family = meta$family_label,
      family_code = meta$family_code,
      spec_id = meta$spec_id,
      model_id = meta$model_id,
      diagnostic_group = "summary",
      iteration = NA_integer_,
      metric_name = "best_objective_value",
      metric_value = fit_obj$best_objective_value %||% NA_real_,
      metric_text = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  
  if (meta$lane == "PH_semiparam") {
    rows[[length(rows) + 1L]] <- data.frame(
      dataset_branch = meta$dataset_branch,
      data_cut = "full",
      subgroup = "overall",
      model_class = "cure_MLE",
      lane = meta$lane,
      family = meta$family_label,
      family_code = meta$family_code,
      spec_id = meta$spec_id,
      model_id = meta$model_id,
      diagnostic_group = "summary",
      iteration = NA_integer_,
      metric_name = "em_iterations",
      metric_value = fit_obj$em_iterations %||% NA_real_,
      metric_text = NA_character_,
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1L]] <- data.frame(
      dataset_branch = meta$dataset_branch,
      data_cut = "full",
      subgroup = "overall",
      model_class = "cure_MLE",
      lane = meta$lane,
      family = meta$family_label,
      family_code = meta$family_code,
      spec_id = meta$spec_id,
      model_id = meta$model_id,
      diagnostic_group = "summary",
      iteration = NA_integer_,
      metric_name = "em_max_relative_change",
      metric_value = fit_obj$em_max_relative_change %||% NA_real_,
      metric_text = NA_character_,
      stringsAsFactors = FALSE
    )
    
    if (!is.null(fit_obj$em_trace) && length(fit_obj$em_trace) > 0L) {
      for (k in seq_along(fit_obj$em_trace)) {
        rows[[length(rows) + 1L]] <- data.frame(
          dataset_branch = meta$dataset_branch,
          data_cut = "full",
          subgroup = "overall",
          model_class = "cure_MLE",
          lane = meta$lane,
          family = meta$family_label,
          family_code = meta$family_code,
          spec_id = meta$spec_id,
          model_id = meta$model_id,
          diagnostic_group = "em_trace",
          iteration = k,
          metric_name = "max_relative_change",
          metric_value = fit_obj$em_trace[k],
          metric_text = NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }
    
    if (!is.null(fit_obj$pseudo_logLik_trace) && length(fit_obj$pseudo_logLik_trace) > 0L) {
      for (k in seq_along(fit_obj$pseudo_logLik_trace)) {
        rows[[length(rows) + 1L]] <- data.frame(
          dataset_branch = meta$dataset_branch,
          data_cut = "full",
          subgroup = "overall",
          model_class = "cure_MLE",
          lane = meta$lane,
          family = meta$family_label,
          family_code = meta$family_code,
          spec_id = meta$spec_id,
          model_id = meta$model_id,
          diagnostic_group = "em_trace",
          iteration = k,
          metric_name = "pseudo_logLik",
          metric_value = fit_obj$pseudo_logLik_trace[k],
          metric_text = NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  do.call(rbind, rows)
}

# 🔴 Import: raw merged dataset and validate schema ===============================

## 🟠 Read: CSV input and choose age variable ===============================
raw0 <- utils::read.csv(
  file = data_path,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)

required_cols <- c(id_var, site_var, sex_var, time_var, status_var)
missing_cols <- setdiff(required_cols, names(raw0))
if (length(missing_cols) > 0L) {
  stop("필수 컬럼이 없습니다: ", paste(missing_cols, collapse = ", "))
}

age_var <- if (age_var_preferred %in% names(raw0)) {
  age_var_preferred
} else if ("age_int" %in% names(raw0)) {
  "age_int"
} else {
  stop("age_exact_entry 또는 age_int 컬럼이 필요합니다.")
}

## 🟠 Normalize: core variables and endpoint coding ===============================
raw0[[id_var]] <- as.character(raw0[[id_var]])
raw0[[site_var]] <- toupper(trimws(as.character(raw0[[site_var]])))
raw0[[sex_var]] <- suppressWarnings(as.numeric(as.character(raw0[[sex_var]])))
raw0[[time_var]] <- suppressWarnings(as.numeric(as.character(raw0[[time_var]])))
raw0[[status_var]] <- suppressWarnings(as.integer(as.character(raw0[[status_var]])))
raw0[[age_var]] <- suppressWarnings(as.numeric(as.character(raw0[[age_var]])))

site_levels_observed <- sort(unique(raw0[[site_var]]))
if (length(site_levels_observed) != 2L) {
  stop("현재 Step5 코드는 site가 정확히 2수준인 병합 데이터만 가정합니다. 관측된 수준: ", paste(site_levels_observed, collapse = ", "))
}

pnu_site_effective <- if (toupper(pnu_site_label) %in% site_levels_observed) {
  toupper(pnu_site_label)
} else {
  toupper(site_reference_label)
}

site_reference_effective <- if (toupper(site_reference_label) %in% site_levels_observed) {
  toupper(site_reference_label)
} else {
  pnu_site_effective
}

snu_site_effective <- if (toupper(snu_site_label) %in% site_levels_observed) {
  toupper(snu_site_label)
} else {
  setdiff(site_levels_observed, pnu_site_effective)[1]
}

raw0$unique_id <- paste(raw0[[site_var]], raw0[[id_var]], sep = "::")
raw0$event_transition <- as.integer(raw0[[status_var]] == 1L)
raw0$time_days_original <- raw0[[time_var]]
raw0$time_days_adjusted <- pmax(raw0[[time_var]], time_zero_epsilon_days)
raw0$time_years <- raw0$time_days_adjusted / days_per_year
raw0$age_raw <- raw0[[age_var]]
raw0$site_num <- ifelse(raw0[[site_var]] == site_reference_effective, 0, 1)

if (!all(stats::na.omit(unique(raw0[[sex_var]])) %in% c(0, 1))) {
  stop("sex_num은 0/1 코딩이어야 합니다.")
}

## 🟠 Filter: complete-case analytic rows ===============================
core_complete_vars <- c("unique_id", "time_years", "event_transition", "age_raw", sex_var, site_var, "site_num", status_var, "time_days_original")
raw_complete <- raw0[stats::complete.cases(raw0[, core_complete_vars]), , drop = FALSE]

if (nrow(raw_complete) == 0L) {
  stop("complete-case 적용 후 남은 관측치가 없습니다.")
}

# 🔴 Prepare: branch-specific analytic datasets ===============================

## 🟠 Split: PNU SNU merged branches ===============================
branch_data_list <- list(
  PNU = raw_complete[raw_complete[[site_var]] == pnu_site_effective, , drop = FALSE],
  SNU = raw_complete[raw_complete[[site_var]] == snu_site_effective, , drop = FALSE],
  merged = raw_complete
)

if (nrow(branch_data_list$PNU) == 0L) stop("PNU 브랜치에 데이터가 없습니다. pnu_site_label 또는 site_reference_label을 확인해 주세요.")
if (nrow(branch_data_list$SNU) == 0L) stop("SNU 브랜치에 데이터가 없습니다. snu_site_label을 확인해 주세요.")

## 🟠 Scale: branch-specific age_s and QC summaries ===============================
preproc_spec_rows <- list()
qc_rows <- list()
ipcw_helper_list <- list()

for (branch_name in names(branch_data_list)) {
  branch_dat <- branch_data_list[[branch_name]]
  
  age_center <- mean(branch_dat$age_raw)
  age_scale_2sd <- 2 * stats::sd(branch_dat$age_raw)
  if (!is.finite(age_scale_2sd) || age_scale_2sd <= 0) {
    stop("branch ", branch_name, " 에서 age 2SD 스케일을 계산할 수 없습니다.")
  }
  
  branch_dat$age_s <- (branch_dat$age_raw - age_center) / age_scale_2sd
  branch_dat$sex_num <- branch_dat[[sex_var]]
  branch_dat$site_std <- branch_dat[[site_var]]
  
  branch_data_list[[branch_name]] <- branch_dat
  ipcw_helper_list[[branch_name]] <- build_ipcw_helper(
    time = branch_dat$time_years,
    event = branch_dat$event_transition
  )
  
  preproc_spec_rows[[length(preproc_spec_rows) + 1L]] <- data.frame(
    dataset_branch = branch_name,
    age_var = age_var,
    age_center = age_center,
    age_scale_2sd = age_scale_2sd,
    sex_var = sex_var,
    sex_coding = "0=Female,1=Male",
    site_reference_label = site_reference_effective,
    pnu_site_effective = pnu_site_effective,
    snu_site_effective = snu_site_effective,
    site_levels_observed = paste(site_levels_observed, collapse = " | "),
    time_scale = "years_from_entry",
    time_zero_epsilon_days = time_zero_epsilon_days,
    stringsAsFactors = FALSE
  )
  
  qc_rows[[length(qc_rows) + 1L]] <- data.frame(
    dataset_branch = branch_name,
    n_subjects = nrow(branch_dat),
    n_events_transition = sum(branch_dat$event_transition == 1L),
    n_censored_total = sum(branch_dat$event_transition == 0L),
    n_right_censor_status0 = sum(branch_dat[[status_var]] == 0L),
    n_remission_censor_status2 = sum(branch_dat[[status_var]] == 2L),
    n_zero_or_negative_time_adjusted = sum(branch_dat$time_days_original <= 0),
    min_followup_years = min(branch_dat$time_years),
    max_followup_years = max(branch_dat$time_years),
    max_event_years = if (any(branch_dat$event_transition == 1L)) max(branch_dat$time_years[branch_dat$event_transition == 1L]) else NA_real_,
    stringsAsFactors = FALSE
  )
}

preproc_spec_df <- do.call(rbind, preproc_spec_rows)
qc_summary_df <- do.call(rbind, qc_rows)

# 🔴 Design: Step5 cure grid and cached matrices ===============================

## 🟠 Enumerate: cure specification table ===============================
build_spec_table <- function(dataset_branch) {
  grid0 <- expand.grid(
    incidence_int_flag = c(0L, 1L),
    latency_int_flag = c(0L, 1L),
    stringsAsFactors = FALSE
  )
  
  if (dataset_branch == "merged") {
    out <- list()
    for (sf in c(0L, 1L)) {
      tmp <- grid0
      tmp$site_flag <- sf
      out[[length(out) + 1L]] <- tmp
    }
    grid0 <- do.call(rbind, out)
    grid0$spec_id <- paste0("C", grid0$incidence_int_flag, grid0$latency_int_flag, "S", grid0$site_flag)
  } else {
    grid0$site_flag <- 0L
    grid0$spec_id <- paste0("C", grid0$incidence_int_flag, grid0$latency_int_flag)
  }
  
  rows <- vector("list", nrow(grid0))
  for (i in seq_len(nrow(grid0))) {
    inc_terms <- c("age_s", "sex_num")
    lat_terms <- c("age_s", "sex_num")
    
    if (grid0$incidence_int_flag[i] == 1L) inc_terms <- c(inc_terms, "age_s:sex_num")
    if (grid0$latency_int_flag[i] == 1L) lat_terms <- c(lat_terms, "age_s:sex_num")
    if (grid0$site_flag[i] == 1L) {
      inc_terms <- c(inc_terms, "site_num")
      lat_terms <- c(lat_terms, "site_num")
    }
    
    rows[[i]] <- data.frame(
      dataset_branch = dataset_branch,
      spec_id = grid0$spec_id[i],
      incidence_int_flag = grid0$incidence_int_flag[i],
      latency_int_flag = grid0$latency_int_flag[i],
      site_flag = grid0$site_flag[i],
      incidence_terms_raw = collapse_terms(inc_terms),
      latency_terms_raw = collapse_terms(lat_terms),
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, rows)
}

spec_table_list <- lapply(names(branch_data_list), build_spec_table)
spec_table_df <- do.call(rbind, spec_table_list)

## 🟠 Define: family lanes and model classes ===============================
family_table <- data.frame(
  lane = c("AFT_param", "AFT_param", "AFT_param", "AFT_param", "PH_semiparam"),
  family_code = c("exp", "weibull", "llogis", "lnorm", "coxph"),
  family_label = c("Exponential", "Weibull", "Log-logistic", "Log-normal", "Cox-latency"),
  stringsAsFactors = FALSE
)

## 🟠 Cache: design matrices by branch and cure spec ===============================
design_cache <- list()

for (branch_name in names(branch_data_list)) {
  dat_branch <- branch_data_list[[branch_name]]
  specs_branch <- spec_table_df[spec_table_df$dataset_branch == branch_name, , drop = FALSE]
  
  for (i in seq_len(nrow(specs_branch))) {
    spec_row <- specs_branch[i, , drop = FALSE]
    cache_key <- paste(branch_name, spec_row$spec_id, sep = "__")
    
    inc_terms <- strsplit(spec_row$incidence_terms_raw, " \\+ ")[[1]]
    lat_terms <- strsplit(spec_row$latency_terms_raw, " \\+ ")[[1]]
    
    inc_bundle <- make_model_matrix_with_pruning(
      data = dat_branch,
      terms = inc_terms,
      include_intercept = TRUE,
      prefix = "inc"
    )
    
    lat_param_bundle <- make_model_matrix_with_pruning(
      data = dat_branch,
      terms = lat_terms,
      include_intercept = TRUE,
      prefix = "latp"
    )
    
    lat_ph_bundle <- make_model_matrix_with_pruning(
      data = dat_branch,
      terms = lat_terms,
      include_intercept = FALSE,
      prefix = "latph"
    )
    
    design_cache[[cache_key]] <- list(
      incidence_bundle = inc_bundle,
      latency_param_bundle = lat_param_bundle,
      latency_ph_bundle = lat_ph_bundle
    )
  }
}

## 🟠 Assemble: final Step5 model grid ===============================
model_grid_rows <- list()

for (i in seq_len(nrow(spec_table_df))) {
  spec_row <- spec_table_df[i, , drop = FALSE]
  for (j in seq_len(nrow(family_table))) {
    fam_row <- family_table[j, , drop = FALSE]
    
    model_id <- paste(
      "step5",
      spec_row$dataset_branch,
      fam_row$lane,
      fam_row$family_code,
      spec_row$spec_id,
      sep = "__"
    )
    
    model_grid_rows[[length(model_grid_rows) + 1L]] <- data.frame(
      dataset_branch = spec_row$dataset_branch,
      data_cut = "full",
      subgroup = "overall",
      model_class = "cure_MLE",
      lane = fam_row$lane,
      family_code = fam_row$family_code,
      family_label = fam_row$family_label,
      spec_id = spec_row$spec_id,
      incidence_int_flag = spec_row$incidence_int_flag,
      latency_int_flag = spec_row$latency_int_flag,
      site_flag = spec_row$site_flag,
      incidence_terms_raw = spec_row$incidence_terms_raw,
      latency_terms_raw = spec_row$latency_terms_raw,
      cache_key = paste(spec_row$dataset_branch, spec_row$spec_id, sep = "__"),
      model_id = model_id,
      rds_file = paste0(sanitize_filename(model_id), ".rds"),
      stringsAsFactors = FALSE
    )
  }
}

model_grid_df <- do.call(rbind, model_grid_rows)

# 🔴 Fit: Step5 MLE cure grid over all branches and families ===============================

## 🟠 Initialize: result containers and templates ===============================
fit_registry_rows <- list()
coef_rows <- list()
pred_subject_rows <- list()
pred_mean_rows <- list()
cure_summary_rows <- list()
perf_rows <- list()
warning_rows <- list()
vcov_rows <- list()
cox_baseline_rows <- list()
diag_rows <- list()
manifest_rows <- list()

template_fit_grid <- data.frame(
  dataset_branch = character(), data_cut = character(), subgroup = character(),
  model_class = character(), lane = character(), family = character(), family_code = character(),
  spec_id = character(), model_id = character(), incidence_int_flag = integer(),
  latency_int_flag = integer(), site_flag = integer(), incidence_terms_raw = character(),
  latency_terms_raw = character(), incidence_terms_effective = character(),
  latency_terms_effective = character(), dropped_incidence_terms = character(),
  dropped_latency_terms = character(), incidence_target = character(), n_subjects = integer(),
  n_events = integer(), n_censored = integer(), logLik = numeric(), partial_logLik = numeric(),
  pseudo_logLik = numeric(), AIC = numeric(), BIC = numeric(), n_parameters = integer(),
  IBS_apparent = numeric(), mean_auc_valid = numeric(), mean_brier_valid = numeric(),
  n_auc_valid = integer(), n_brier_valid = integer(), convergence_flag = logical(),
  convergence_code = integer(), convergence_message = character(), boundary_flag = logical(),
  singular_hessian_flag = logical(), hessian_condition_number = numeric(),
  vcov_method = character(), max_event_time_years = numeric(), max_followup_time_years = numeric(),
  mean_pi_uncured = numeric(), mean_cure_fraction = numeric(), em_iterations = integer(),
  em_max_relative_change = numeric(), n_start_attempts = integer(),
  n_successful_starts = integer(), fit_ok = logical(), fit_status = character(),
  fit_time_sec = numeric(), warning_count = integer(), rds_file = character(),
  loaded_existing_rds = logical(), stringsAsFactors = FALSE
)

template_param_est <- data.frame(
  dataset_branch = character(), data_cut = character(), subgroup = character(),
  model_class = character(), lane = character(), family = character(), family_code = character(),
  spec_id = character(), model_id = character(), component = character(),
  coefficient_scale = character(), exp_scale_label = character(), term = character(),
  estimate = numeric(), std_error = numeric(), z_value = numeric(), p_value = numeric(),
  wald_lcl = numeric(), wald_ucl = numeric(), exp_estimate = numeric(),
  exp_wald_lcl = numeric(), exp_wald_ucl = numeric(), se_method = character(),
  stringsAsFactors = FALSE
)

template_subject_pred <- data.frame(
  dataset_branch = character(), data_cut = character(), subgroup = character(),
  model_class = character(), lane = character(), family = character(), family_code = character(),
  spec_id = character(), model_id = character(), id_original = character(),
  unique_id = character(), site_std = character(), site_num = numeric(),
  sex_num = numeric(), age_raw = numeric(), age_s = numeric(),
  time_days_original = numeric(), time_years = numeric(), status_num_observed = integer(),
  event_transition = integer(), year = numeric(), surv = numeric(), risk = numeric(),
  surv_uncured = numeric(), pi_uncured = numeric(), cure_fraction = numeric(),
  beyond_last_event = logical(), beyond_max_followup = logical(), is_extrapolated = logical(),
  stringsAsFactors = FALSE
)

template_mean_pred <- data.frame(
  dataset_branch = character(), data_cut = character(), subgroup = character(),
  model_class = character(), lane = character(), family = character(), family_code = character(),
  spec_id = character(), model_id = character(), year = numeric(), mean_surv = numeric(),
  sd_surv = numeric(), surv_q25 = numeric(), surv_median = numeric(), surv_q75 = numeric(),
  mean_risk = numeric(), sd_risk = numeric(), risk_q25 = numeric(),
  risk_median = numeric(), risk_q75 = numeric(), mean_surv_uncured = numeric(),
  mean_pi_uncured = numeric(), mean_cure_fraction = numeric(),
  median_cure_fraction = numeric(), n_subjects = integer(),
  beyond_last_event = logical(), beyond_max_followup = logical(),
  stringsAsFactors = FALSE
)

template_cure_summary <- data.frame(
  dataset_branch = character(), data_cut = character(), subgroup = character(),
  model_class = character(), lane = character(), family = character(), family_code = character(),
  spec_id = character(), model_id = character(), incidence_target = character(),
  mean_pi_uncured = numeric(), mean_cure_fraction = numeric(), median_cure_fraction = numeric(),
  cure_fraction_q25 = numeric(), cure_fraction_q75 = numeric(),
  min_cure_fraction = numeric(), max_cure_fraction = numeric(),
  boundary_flag = logical(), convergence_flag = logical(),
  IBS_apparent = numeric(), mean_auc_valid = numeric(), mean_brier_valid = numeric(),
  stringsAsFactors = FALSE
)

template_perf <- data.frame(
  dataset_branch = character(), data_cut = character(), subgroup = character(),
  model_class = character(), lane = character(), family = character(), family_code = character(),
  spec_id = character(), model_id = character(), metric_scope = character(),
  year = numeric(), auc_cd = numeric(), brier_ipcw = numeric(),
  G_t_censoring_survival = numeric(), n_cases = integer(), n_controls = integer(),
  auc_valid = logical(), brier_valid = logical(),
  beyond_max_followup = logical(), beyond_last_event = logical(),
  IBS_apparent = numeric(), stringsAsFactors = FALSE
)

template_warning <- data.frame(
  dataset_branch = character(), lane = character(), family = character(),
  spec_id = character(), model_id = character(), warning_type = character(),
  warning_message = character(), stringsAsFactors = FALSE
)

template_vcov <- data.frame(
  dataset_branch = character(), data_cut = character(), subgroup = character(),
  model_class = character(), lane = character(), family = character(), family_code = character(),
  spec_id = character(), model_id = character(), param_row = character(),
  param_col = character(), covariance = numeric(), vcov_method = character(),
  cross_block_available = logical(), stringsAsFactors = FALSE
)

template_baseline <- data.frame(
  dataset_branch = character(), data_cut = character(), subgroup = character(),
  model_class = character(), lane = character(), family = character(), family_code = character(),
  spec_id = character(), model_id = character(), time_years = numeric(),
  baseline_cumhaz = numeric(), baseline_surv = numeric(), zero_tail_applied = logical(),
  stringsAsFactors = FALSE
)

template_diag <- data.frame(
  dataset_branch = character(), data_cut = character(), subgroup = character(),
  model_class = character(), lane = character(), family = character(), family_code = character(),
  spec_id = character(), model_id = character(), diagnostic_group = character(),
  iteration = integer(), metric_name = character(), metric_value = numeric(),
  metric_text = character(), stringsAsFactors = FALSE
)

## 🟠 Iterate: fit or load every Step5 model ===============================
for (i in seq_len(nrow(model_grid_df))) {
  meta <- model_grid_df[i, , drop = FALSE]
  meta_list <- as.list(meta)
  
  message(sprintf("[%s] (%d/%d) %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), i, nrow(model_grid_df), meta$model_id))
  
  dat_branch <- branch_data_list[[meta$dataset_branch]]
  design_obj <- design_cache[[meta$cache_key]]
  ipcw_helper <- ipcw_helper_list[[meta$dataset_branch]]
  
  incidence_bundle <- design_obj$incidence_bundle
  latency_bundle <- if (meta$lane == "AFT_param") design_obj$latency_param_bundle else design_obj$latency_ph_bundle
  
  rds_path <- file.path(export_path, meta$rds_file)
  fit_obj <- NULL
  loaded_existing <- FALSE
  refit_reason <- NA_character_
  
  if (isTRUE(reuse_existing_model_rds) && !isTRUE(overwrite_model_rds) && file.exists(rds_path)) {
    load_try <- try(readRDS(rds_path), silent = TRUE)
    if (!inherits(load_try, "try-error")) {
      temp_fit <- load_try
      version_ok <- identical(temp_fit$script_version %||% "", script_version)
      if (!isTRUE(force_refit_if_script_version_mismatch) || version_ok) {
        fit_obj <- temp_fit
        loaded_existing <- TRUE
      } else {
        refit_reason <- "script_version_mismatch"
      }
    } else {
      refit_reason <- "rds_read_error"
    }
  }
  
  if (is.null(fit_obj)) {
    if (meta$lane == "AFT_param") {
      fit_obj <- fit_parametric_mixture_cure(
        data = dat_branch,
        incidence_bundle = incidence_bundle,
        latency_bundle = latency_bundle,
        family = meta$family_code
      )
    } else {
      fit_obj <- fit_ph_cure_em(
        data = dat_branch,
        incidence_bundle = incidence_bundle,
        latency_bundle = latency_bundle
      )
    }
    
    fit_obj$script_version <- script_version
    fit_obj$meta <- meta_list
    fit_obj$training_data_minimal <- dat_branch[, c(
      "unique_id", id_var, "site_std", "site_num", "sex_num", "age_raw", "age_s",
      "time_days_original", "time_years", status_var, "event_transition"
    ), drop = FALSE]
    fit_obj$loaded_existing_rds <- FALSE
    fit_obj$refit_reason <- refit_reason
    
    saveRDS(fit_obj, rds_path)
  } else {
    fit_obj$loaded_existing_rds <- TRUE
    fit_obj$refit_reason <- NA_character_
  }
  
  manifest_rows <- record_manifest_row(
    manifest_list = manifest_rows,
    file_name = basename(rds_path),
    file_type = "model_rds",
    n_rows = NA_integer_,
    note = if (loaded_existing) "existing model fit reused" else paste("new model fit saved", refit_reason %||% ""),
    model_id = meta$model_id
  )
  
  if (!isTRUE(fit_obj$fit_ok)) {
    fit_registry_rows[[length(fit_registry_rows) + 1L]] <- data.frame(
      dataset_branch = meta$dataset_branch,
      data_cut = "full",
      subgroup = "overall",
      model_class = "cure_MLE",
      lane = meta$lane,
      family = meta$family_label,
      family_code = meta$family_code,
      spec_id = meta$spec_id,
      model_id = meta$model_id,
      incidence_int_flag = meta$incidence_int_flag,
      latency_int_flag = meta$latency_int_flag,
      site_flag = meta$site_flag,
      incidence_terms_raw = meta$incidence_terms_raw,
      latency_terms_raw = meta$latency_terms_raw,
      incidence_terms_effective = collapse_terms(incidence_bundle$map$term_label),
      latency_terms_effective = collapse_terms(latency_bundle$map$term_label),
      dropped_incidence_terms = collapse_terms(incidence_bundle$dropped_terms),
      dropped_latency_terms = collapse_terms(latency_bundle$dropped_terms),
      incidence_target = NA_character_,
      n_subjects = nrow(dat_branch),
      n_events = sum(dat_branch$event_transition == 1L),
      n_censored = sum(dat_branch$event_transition == 0L),
      logLik = NA_real_,
      partial_logLik = NA_real_,
      pseudo_logLik = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      n_parameters = NA_integer_,
      IBS_apparent = NA_real_,
      mean_auc_valid = NA_real_,
      mean_brier_valid = NA_real_,
      n_auc_valid = NA_integer_,
      n_brier_valid = NA_integer_,
      convergence_flag = FALSE,
      convergence_code = NA_integer_,
      convergence_message = paste(fit_obj$warnings %||% "fit failed", collapse = " | "),
      boundary_flag = NA,
      singular_hessian_flag = NA,
      hessian_condition_number = NA_real_,
      vcov_method = NA_character_,
      max_event_time_years = if (any(dat_branch$event_transition == 1L)) max(dat_branch$time_years[dat_branch$event_transition == 1L]) else NA_real_,
      max_followup_time_years = max(dat_branch$time_years),
      mean_pi_uncured = NA_real_,
      mean_cure_fraction = NA_real_,
      em_iterations = fit_obj$em_iterations %||% NA_integer_,
      em_max_relative_change = fit_obj$em_max_relative_change %||% NA_real_,
      n_start_attempts = fit_obj$n_start_attempts %||% NA_integer_,
      n_successful_starts = fit_obj$n_successful_starts %||% NA_integer_,
      fit_ok = FALSE,
      fit_status = fit_obj$fit_status %||% "error",
      fit_time_sec = fit_obj$fit_time_sec %||% NA_real_,
      warning_count = length(fit_obj$warnings %||% character()),
      rds_file = basename(rds_path),
      loaded_existing_rds = loaded_existing,
      stringsAsFactors = FALSE
    )
    
    warning_rows[[length(warning_rows) + 1L]] <- data.frame(
      dataset_branch = meta$dataset_branch,
      lane = meta$lane,
      family = meta$family_label,
      spec_id = meta$spec_id,
      model_id = meta$model_id,
      warning_type = "fit_error",
      warning_message = paste(fit_obj$warnings %||% "unknown error", collapse = " | "),
      stringsAsFactors = FALSE
    )
    
    next
  }
  
  pred_obj <- if (meta$lane == "AFT_param") {
    predict_parametric_cure_fit(
      fit_obj = fit_obj,
      Znew = incidence_bundle$matrix,
      Xnew = latency_bundle$matrix,
      years = prediction_years
    )
  } else {
    predict_ph_cure_fit(
      fit_obj = fit_obj,
      Znew = incidence_bundle$matrix,
      Xnew = latency_bundle$matrix,
      years = prediction_years
    )
  }
  
  pred_subject_df <- assemble_subject_prediction_df(
    data = dat_branch,
    pred_obj = pred_obj,
    meta = meta_list,
    last_event_time = fit_obj$max_event_time_years %||% NA_real_,
    max_followup_time = fit_obj$max_followup_time_years %||% NA_real_
  )
  
  pred_mean_df <- assemble_mean_prediction_df(pred_subject_df)
  perf_df <- compute_time_dependent_metrics(
    data = dat_branch,
    pred_obj = pred_obj,
    ipcw_helper = ipcw_helper,
    years = prediction_years,
    meta = meta_list
  )
  
  if (length(fit_obj$warnings) > 0L) {
    warning_rows[[length(warning_rows) + 1L]] <- data.frame(
      dataset_branch = meta$dataset_branch,
      lane = meta$lane,
      family = meta$family_label,
      spec_id = meta$spec_id,
      model_id = meta$model_id,
      warning_type = "warning",
      warning_message = paste(unique(fit_obj$warnings), collapse = " | "),
      stringsAsFactors = FALSE
    )
  }
  
  if (isTRUE(write_param_estimates)) {
    coef_rows[[length(coef_rows) + 1L]] <- assemble_param_coef_df(
      fit_obj = fit_obj,
      meta = meta_list
    )
  }
  
  if (isTRUE(write_subject_level_predictions)) {
    pred_subject_rows[[length(pred_subject_rows) + 1L]] <- pred_subject_df
  }
  
  if (isTRUE(write_mean_level_predictions)) {
    pred_mean_rows[[length(pred_mean_rows) + 1L]] <- pred_mean_df
  }
  
  if (isTRUE(write_cure_summary)) {
    cure_summary_rows[[length(cure_summary_rows) + 1L]] <- assemble_cure_summary_row(
      fit_obj = fit_obj,
      mean_pred_df = pred_mean_df,
      perf_df = perf_df,
      meta = meta_list
    )
  }
  
  if (isTRUE(write_perf_metrics)) {
    perf_rows[[length(perf_rows) + 1L]] <- perf_df
  }
  
  if (isTRUE(write_vcov_long)) {
    vc_df <- assemble_vcov_long_df(
      fit_obj = fit_obj,
      meta = meta_list
    )
    if (!is.null(vc_df) && nrow(vc_df) > 0L) {
      vcov_rows[[length(vcov_rows) + 1L]] <- vc_df
    }
  }
  
  if (isTRUE(write_cox_baseline_curve) && meta$lane == "PH_semiparam") {
    bh_df <- assemble_cox_baseline_df(
      fit_obj = fit_obj,
      meta = meta_list
    )
    if (!is.null(bh_df) && nrow(bh_df) > 0L) {
      cox_baseline_rows[[length(cox_baseline_rows) + 1L]] <- bh_df
    }
  }
  
  if (isTRUE(write_fit_diagnostics)) {
    diag_rows[[length(diag_rows) + 1L]] <- assemble_fit_diagnostics_df(
      fit_obj = fit_obj,
      meta = meta_list
    )
  }
  
  fit_registry_rows[[length(fit_registry_rows) + 1L]] <- data.frame(
    dataset_branch = meta$dataset_branch,
    data_cut = "full",
    subgroup = "overall",
    model_class = "cure_MLE",
    lane = meta$lane,
    family = meta$family_label,
    family_code = meta$family_code,
    spec_id = meta$spec_id,
    model_id = meta$model_id,
    incidence_int_flag = meta$incidence_int_flag,
    latency_int_flag = meta$latency_int_flag,
    site_flag = meta$site_flag,
    incidence_terms_raw = meta$incidence_terms_raw,
    latency_terms_raw = meta$latency_terms_raw,
    incidence_terms_effective = collapse_terms(fit_obj$incidence_bundle$map$term_label),
    latency_terms_effective = collapse_terms(fit_obj$latency_bundle$map$term_label),
    dropped_incidence_terms = collapse_terms(fit_obj$incidence_bundle$dropped_terms),
    dropped_latency_terms = collapse_terms(fit_obj$latency_bundle$dropped_terms),
    incidence_target = fit_obj$incidence_target %||% NA_character_,
    n_subjects = nrow(dat_branch),
    n_events = sum(dat_branch$event_transition == 1L),
    n_censored = sum(dat_branch$event_transition == 0L),
    logLik = fit_obj$logLik %||% NA_real_,
    partial_logLik = fit_obj$partial_logLik %||% NA_real_,
    pseudo_logLik = fit_obj$pseudo_logLik %||% NA_real_,
    AIC = fit_obj$AIC %||% NA_real_,
    BIC = fit_obj$BIC %||% NA_real_,
    n_parameters = fit_obj$n_parameters %||% NA_integer_,
    IBS_apparent = if (all(is.na(perf_df$IBS_apparent))) NA_real_ else unique(stats::na.omit(perf_df$IBS_apparent))[1],
    mean_auc_valid = if (any(is.finite(perf_df$auc_cd))) mean(perf_df$auc_cd, na.rm = TRUE) else NA_real_,
    mean_brier_valid = if (any(is.finite(perf_df$brier_ipcw))) mean(perf_df$brier_ipcw, na.rm = TRUE) else NA_real_,
    n_auc_valid = sum(is.finite(perf_df$auc_cd)),
    n_brier_valid = sum(is.finite(perf_df$brier_ipcw)),
    convergence_flag = fit_obj$convergence_flag %||% FALSE,
    convergence_code = fit_obj$convergence_code %||% NA_integer_,
    convergence_message = fit_obj$convergence_message %||% "",
    boundary_flag = fit_obj$boundary_flag %||% NA,
    singular_hessian_flag = fit_obj$singular_hessian_flag %||% NA,
    hessian_condition_number = fit_obj$hessian_condition_number %||% NA_real_,
    vcov_method = fit_obj$vcov_method %||% NA_character_,
    max_event_time_years = fit_obj$max_event_time_years %||% NA_real_,
    max_followup_time_years = fit_obj$max_followup_time_years %||% NA_real_,
    mean_pi_uncured = mean(fit_obj$pi_uncured_train),
    mean_cure_fraction = mean(fit_obj$cure_fraction_train),
    em_iterations = fit_obj$em_iterations %||% NA_integer_,
    em_max_relative_change = fit_obj$em_max_relative_change %||% NA_real_,
    n_start_attempts = fit_obj$n_start_attempts %||% NA_integer_,
    n_successful_starts = fit_obj$n_successful_starts %||% NA_integer_,
    fit_ok = fit_obj$fit_ok %||% FALSE,
    fit_status = fit_obj$fit_status %||% "unknown",
    fit_time_sec = fit_obj$fit_time_sec %||% NA_real_,
    warning_count = length(fit_obj$warnings %||% character()),
    rds_file = basename(rds_path),
    loaded_existing_rds = loaded_existing,
    stringsAsFactors = FALSE
  )
}

# 🔴 Combine: final Step5 tables from all fitted models ===============================

## 🟠 Bind: all result containers into data frames ===============================
fit_mc_grid_df <- bind_rows_or_template(fit_registry_rows, template_fit_grid)
param_estimates_mc_grid_df <- bind_rows_or_template(coef_rows, template_param_est)
pred_subject_mc_grid_df <- bind_rows_or_template(pred_subject_rows, template_subject_pred)
pred_mean_mc_grid_df <- bind_rows_or_template(pred_mean_rows, template_mean_pred)
cure_only_mc_grid_df <- bind_rows_or_template(cure_summary_rows, template_cure_summary)
perf_mc_grid_df <- bind_rows_or_template(perf_rows, template_perf)
fit_warnings_df <- bind_rows_or_template(warning_rows, template_warning)
vcov_long_mc_grid_df <- bind_rows_or_template(vcov_rows, template_vcov)
cox_baseline_curve_grid_df <- bind_rows_or_template(cox_baseline_rows, template_baseline)
fit_diagnostics_grid_df <- bind_rows_or_template(diag_rows, template_diag)

## 🟠 Add: wide yearly columns to cure summary when available ===============================
if (nrow(cure_only_mc_grid_df) > 0L && nrow(pred_mean_mc_grid_df) > 0L) {
  years_unique <- sort(unique(pred_mean_mc_grid_df$year))
  for (yr in years_unique) {
    idx_year <- pred_mean_mc_grid_df$year == yr
    tmp <- pred_mean_mc_grid_df[idx_year, c("model_id", "mean_surv", "mean_risk", "mean_surv_uncured"), drop = FALSE]
    names(tmp)[2:4] <- c(
      paste0("mean_surv_y", yr),
      paste0("mean_risk_y", yr),
      paste0("mean_su_y", yr)
    )
    cure_only_mc_grid_df <- merge(cure_only_mc_grid_df, tmp, by = "model_id", all.x = TRUE, sort = FALSE)
  }
}

# 🔴 Export: Step5 CSV outputs and manifest ===============================

## 🟠 Write: Step5 source-of-truth csv files ===============================
file_model_grid_spec <- file.path(export_path, "step5_model_grid_spec.csv")
file_preproc_spec <- file.path(export_path, "step5_preproc_spec.csv")
file_qc_summary <- file.path(export_path, "step5_qc_summary.csv")
file_fit_grid <- file.path(export_path, "step5_fit_mc_grid.csv")
file_cure_only <- file.path(export_path, "step5_cure_only_mc_grid.csv")
file_param_est <- file.path(export_path, "step5_param_estimates_mc_grid.csv")
file_pred_subject <- file.path(export_path, "step5_pred_subject_mc_grid.csv")
file_pred_mean <- file.path(export_path, "step5_pred_mean_mc_grid.csv")
file_perf <- file.path(export_path, "step5_perf_mc_grid.csv")
file_fit_warnings <- file.path(export_path, "step5_fit_warnings.csv")
file_vcov_long <- file.path(export_path, "step5_vcov_long_mc_grid.csv")
file_cox_baseline <- file.path(export_path, "step5_cox_baseline_curve_grid.csv")
file_fit_diag <- file.path(export_path, "step5_fit_diagnostics_grid.csv")
file_manifest <- file.path(export_path, "step5_run_manifest.csv")

write_csv_utf8(model_grid_df, file_model_grid_spec)
manifest_rows <- record_manifest_row(manifest_rows, basename(file_model_grid_spec), "csv", nrow(model_grid_df), "Step5 전체 model grid", NA_character_)

write_csv_utf8(preproc_spec_df, file_preproc_spec)
manifest_rows <- record_manifest_row(manifest_rows, basename(file_preproc_spec), "csv", nrow(preproc_spec_df), "branch별 preprocessing spec", NA_character_)

write_csv_utf8(qc_summary_df, file_qc_summary)
manifest_rows <- record_manifest_row(manifest_rows, basename(file_qc_summary), "csv", nrow(qc_summary_df), "branch별 QC summary", NA_character_)

write_csv_utf8(fit_mc_grid_df, file_fit_grid)
manifest_rows <- record_manifest_row(manifest_rows, basename(file_fit_grid), "csv", nrow(fit_mc_grid_df), "모델 적합 레지스트리와 핵심 요약", NA_character_)

write_csv_utf8(cure_only_mc_grid_df, file_cure_only)
manifest_rows <- record_manifest_row(manifest_rows, basename(file_cure_only), "csv", nrow(cure_only_mc_grid_df), "모델별 cure summary와 1-10년 mean risk/survival", NA_character_)

write_csv_utf8(param_estimates_mc_grid_df, file_param_est)
manifest_rows <- record_manifest_row(manifest_rows, basename(file_param_est), "csv", nrow(param_estimates_mc_grid_df), "incidence/latency/ancillary coefficients", NA_character_)

write_csv_utf8(pred_subject_mc_grid_df, file_pred_subject)
manifest_rows <- record_manifest_row(manifest_rows, basename(file_pred_subject), "csv", nrow(pred_subject_mc_grid_df), "subject-level yearly predictions with covariates", NA_character_)

write_csv_utf8(pred_mean_mc_grid_df, file_pred_mean)
manifest_rows <- record_manifest_row(manifest_rows, basename(file_pred_mean), "csv", nrow(pred_mean_mc_grid_df), "cohort-average yearly predictions and quantiles", NA_character_)

write_csv_utf8(perf_mc_grid_df, file_perf)
manifest_rows <- record_manifest_row(manifest_rows, basename(file_perf), "csv", nrow(perf_mc_grid_df), "time-dependent AUC/Brier/IBS apparent performance", NA_character_)

write_csv_utf8(fit_warnings_df, file_fit_warnings)
manifest_rows <- record_manifest_row(manifest_rows, basename(file_fit_warnings), "csv", nrow(fit_warnings_df), "fit warnings and errors", NA_character_)

write_csv_utf8(vcov_long_mc_grid_df, file_vcov_long)
manifest_rows <- record_manifest_row(manifest_rows, basename(file_vcov_long), "csv", nrow(vcov_long_mc_grid_df), "parameter covariance long table", NA_character_)

write_csv_utf8(cox_baseline_curve_grid_df, file_cox_baseline)
manifest_rows <- record_manifest_row(manifest_rows, basename(file_cox_baseline), "csv", nrow(cox_baseline_curve_grid_df), "Cox-latency baseline cumulative hazard and survival", NA_character_)

write_csv_utf8(fit_diagnostics_grid_df, file_fit_diag)
manifest_rows <- record_manifest_row(manifest_rows, basename(file_fit_diag), "csv", nrow(fit_diagnostics_grid_df), "optimizer and EM diagnostics", NA_character_)

## 🟠 Write: manifest for all exported files ===============================
manifest_df <- do.call(rbind, manifest_rows)
manifest_df$script_version <- script_version
manifest_df$data_path <- data_path
manifest_df$export_path <- export_path
manifest_df$created_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

write_csv_utf8(manifest_df, file_manifest)

message("Step5 완료: ", export_path)