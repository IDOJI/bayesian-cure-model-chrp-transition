# 🔴 Configure: paths and screening options ===============================

data_path <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
export_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦3.Step3_non-zero cure fraction screening/attachments'

config <- list(
  seed = 20260316L,
  time_scale_days_per_year = 365.25,
  time_epsilon_days = 0.5,
  alpha_screen = 0.05,
  
  min_group_n = 30L,
  min_group_events = 5L,
  min_group_censored = 5L,
  
  hsu_primary_z = "age_s",
  hsu_primary_kappa_target = 3,
  hsu_kappa_grid = c(0.5, 1, 3, 5, 7),
  hsu_theta_grid_length = 161L,
  hsu_n_resamples = 499L,
  hsu_primary_baseline_candidates = c("weibull", "loglogistic"),
  hsu_include_sex_in_latency = TRUE,
  hsu_include_site_in_latency_for_merged = TRUE,
  hsu_kappa_plateau_fraction = 0.95,
  
  xie_n_boot = 499L,
  
  receus_threshold_pi_cure = 0.025,
  receus_threshold_r = 0.05,
  receus_candidate_families = c("exponential", "weibull", "gamma", "loglogistic"),
  
  run_practical_followup = FALSE,
  save_bootstrap_vectors = FALSE
)

dir.create(export_path, recursive = TRUE, showWarnings = FALSE)
set.seed(config$seed)
options(stringsAsFactors = FALSE, scipen = 999)

# 🔴 Import: required namespace checks ===============================

## 🟠 Validate: package availability ===============================

required_pkgs <- c("survival")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("다음 패키지가 필요합니다: ", paste(missing_pkgs, collapse = ", "))
}

# 🔴 Define: utility helpers and numerical kernels ===============================

## 🟠 Create: small scalar helpers ===============================

trim_ws <- function(x) {
  trimws(as.character(x), which = "both")
}

sanitize_name <- function(x) {
  x <- trim_ws(x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "NA_LABEL", x)
}

clamp_prob <- function(p, eps = 1e-8) {
  pmin(pmax(p, eps), 1 - eps)
}

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, type = 8))
}

logspace_add <- function(logx, logy) {
  m <- pmax(logx, logy)
  m + log(exp(logx - m) + exp(logy - m))
}

drop_constant_cols <- function(X) {
  if (is.null(X)) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }
  X <- as.matrix(X)
  if (ncol(X) == 0) {
    return(X)
  }
  keep <- apply(X, 2, function(v) {
    stats::sd(v, na.rm = TRUE) > 0 && all(is.finite(v))
  })
  if (!any(keep)) {
    return(matrix(numeric(0), nrow = nrow(X), ncol = 0))
  }
  X[, keep, drop = FALSE]
}

safe_inverse <- function(M, ridge = 1e-8) {
  M <- (M + t(M)) / 2
  ee <- eigen(M, symmetric = TRUE)
  vals <- ee$values
  vals[!is.finite(vals) | vals < ridge] <- ridge
  ee$vectors %*% diag(1 / vals, nrow = length(vals)) %*% t(ee$vectors)
}

rbind_fill_df <- function(x) {
  if (length(x) == 0) {
    return(data.frame())
  }
  nm <- unique(unlist(lapply(x, names), use.names = FALSE))
  out <- lapply(x, function(df) {
    miss <- setdiff(nm, names(df))
    if (length(miss) > 0) {
      for (m in miss) df[[m]] <- NA
    }
    df[nm]
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

calc_tau_admin <- function(time, event) {
  if (any(event == 0L)) {
    return(max(time[event == 0L], na.rm = TRUE))
  }
  max(time, na.rm = TRUE)
}

# ## 🟠 Create: Kaplan-Meier and follow-up summaries ===============================

km_surv_at_times <- function(time, event, times) {
  fit <- survival::survfit(survival::Surv(time, event) ~ 1)
  s <- summary(fit, times = times, extend = TRUE)$surv
  if (length(s) == 0) {
    return(rep(1, length(times)))
  }
  as.numeric(s)
}

reverse_km_median <- function(time, event) {
  fit <- survival::survfit(survival::Surv(time, 1L - event) ~ 1)
  tb <- suppressWarnings(summary(fit)$table)
  if (length(tb) == 0) return(NA_real_)
  med <- unname(tb["median"])
  if (is.null(med) || length(med) == 0) return(NA_real_)
  as.numeric(med)
}

compute_followup_context <- function(df_group) {
  time <- df_group$time_years
  event <- df_group$event
  tau_max <- max(time, na.rm = TRUE)
  tau_admin <- calc_tau_admin(time, event)
  last_event <- if (any(event == 1L)) max(time[event == 1L], na.rm = TRUE) else NA_real_
  km_tail_surv <- km_surv_at_times(time, event, tau_max)
  n_censored_after_last_event <- if (is.finite(last_event)) {
    sum(event == 0L & time > last_event)
  } else {
    sum(event == 0L)
  }
  plateau_length <- if (is.finite(last_event)) max(0, tau_max - last_event) else NA_real_
  
  cci <- if (is.finite(tau_admin) && tau_admin > 0) {
    sum(time) / sum(ifelse(event == 1L, time, tau_admin))
  } else {
    NA_real_
  }
  
  spt <- if (is.finite(tau_admin) && tau_admin > 0) {
    sum(ifelse(event == 1L, tau_admin, pmin(time, tau_admin))) / (nrow(df_group) * tau_admin)
  } else {
    NA_real_
  }
  
  list(
    tau_max = tau_max,
    tau_admin = tau_admin,
    last_event = last_event,
    plateau_length = plateau_length,
    n_censored_after_last_event = n_censored_after_last_event,
    km_tail_surv = km_tail_surv,
    reverse_km_median_followup = reverse_km_median(time, event),
    cci = cci,
    spt = spt
  )
}

# ## 🟠 Create: data validation and preprocessing helpers ===============================

validate_input_columns <- function(df) {
  required_cols <- c("id", "site", "sex_num", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("입력 데이터에 다음 필수 컬럼이 없습니다: ", paste(missing_cols, collapse = ", "))
  }
  if (!("age_exact_entry" %in% names(df) || "age_int" %in% names(df))) {
    stop("입력 데이터에 age_exact_entry 또는 age_int 중 하나가 필요합니다.")
  }
}

choose_age_column <- function(df) {
  if ("age_exact_entry" %in% names(df) && any(!is.na(df$age_exact_entry))) {
    return(as.numeric(df$age_exact_entry))
  }
  as.numeric(df$age_int)
}

# 🔴 Read: merged cohort and preprocess analysis variables ===============================

## 🟠 Load: raw CSV and harmonize primary endpoint ===============================

raw_df <- read.csv(
  file = data_path,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "NaN", "NULL")
)

validate_input_columns(raw_df)

raw_df$id <- as.character(raw_df$id)
raw_df$site <- trim_ws(raw_df$site)
raw_df$sex_num <- as.numeric(raw_df$sex_num)
raw_df$days_followup <- as.numeric(raw_df$days_followup)
raw_df$status_num <- as.integer(raw_df$status_num)

if ("sex_fact" %in% names(raw_df)) {
  raw_df$sex_fact <- trim_ws(raw_df$sex_fact)
}

if (any(is.na(raw_df$sex_num)) && "sex_fact" %in% names(raw_df)) {
  raw_df$sex_num[is.na(raw_df$sex_num) & tolower(raw_df$sex_fact) == "female"] <- 0
  raw_df$sex_num[is.na(raw_df$sex_num) & tolower(raw_df$sex_fact) == "male"] <- 1
}

analysis_df <- raw_df
analysis_df$age_raw <- choose_age_column(analysis_df)
analysis_df$event <- as.integer(analysis_df$status_num == 1L)
analysis_df$censor_type <- ifelse(analysis_df$status_num == 2L, "remission", "right_censoring")
analysis_df$unique_id <- paste0(analysis_df$site, "::", analysis_df$id)

analysis_df <- analysis_df[
  !is.na(analysis_df$unique_id) &
    !is.na(analysis_df$site) &
    !is.na(analysis_df$age_raw) &
    !is.na(analysis_df$sex_num) &
    !is.na(analysis_df$days_followup) &
    !is.na(analysis_df$status_num),
]

analysis_df <- analysis_df[
  analysis_df$sex_num %in% c(0, 1) &
    analysis_df$days_followup >= 0,
]

analysis_df$time_days <- pmax(analysis_df$days_followup, config$time_epsilon_days)
analysis_df$time_years <- analysis_df$time_days / config$time_scale_days_per_year

age_center <- mean(analysis_df$age_raw, na.rm = TRUE)
age_scale_2sd <- 2 * stats::sd(analysis_df$age_raw, na.rm = TRUE)
if (!is.finite(age_scale_2sd) || age_scale_2sd <= 0) {
  age_scale_2sd <- 1
}
analysis_df$age_s <- (analysis_df$age_raw - age_center) / age_scale_2sd

analysis_df <- analysis_df[order(analysis_df$site, analysis_df$unique_id, analysis_df$time_years), ]
rownames(analysis_df) <- NULL

preproc_spec <- data.frame(
  age_center = age_center,
  age_scale_2sd = age_scale_2sd,
  sex_coding = "0=Female,1=Male",
  time_unit = "years_since_entry",
  endpoint = "transition_only",
  remission_as_censor = TRUE,
  stringsAsFactors = FALSE
)

# 🔴 Construct: merged-overall and site-specific screening sets ===============================

## 🟠 Build: screening groups from site labels ===============================

site_levels <- sort(unique(analysis_df$site))

analysis_groups <- list(
  list(
    group_id = "merged_overall",
    group_display = "MERGED",
    group_kind = "merged_overall",
    site_label = NA_character_,
    data = analysis_df
  )
)

for (s in site_levels) {
  analysis_groups[[length(analysis_groups) + 1L]] <- list(
    group_id = paste0("site__", sanitize_name(s)),
    group_display = s,
    group_kind = "site_specific",
    site_label = s,
    data = analysis_df[analysis_df$site == s, , drop = FALSE]
  )
}

# 🔴 Implement: Hsu sup-score screening machinery ===============================

## 🟠 Create: covariate builders for Hsu screening ===============================

build_hsu_latency_matrix <- function(df_group, group_kind, cfg) {
  X_parts <- list()
  
  if ("age_s" %in% names(df_group)) {
    X_parts[["age_s"]] <- matrix(df_group$age_s, ncol = 1, dimnames = list(NULL, "age_s"))
  }
  
  if (isTRUE(cfg$hsu_include_sex_in_latency) && "sex_num" %in% names(df_group)) {
    X_parts[["sex_num"]] <- matrix(df_group$sex_num, ncol = 1, dimnames = list(NULL, "sex_num"))
  }
  
  if (identical(group_kind, "merged_overall") && isTRUE(cfg$hsu_include_site_in_latency_for_merged)) {
    mm <- model.matrix(~ site, data = df_group)
    if (ncol(mm) > 1) {
      X_parts[["site"]] <- mm[, -1, drop = FALSE]
    }
  }
  
  if (length(X_parts) == 0) {
    return(matrix(numeric(0), nrow = nrow(df_group), ncol = 0))
  }
  
  X <- do.call(cbind, X_parts)
  X <- drop_constant_cols(X)
  
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("x", seq_len(ncol(X)))
  }
  X
}

build_hsu_z_vector <- function(df_group, cfg) {
  z_name <- cfg$hsu_primary_z
  if (!z_name %in% names(df_group)) {
    return(NULL)
  }
  z <- as.numeric(df_group[[z_name]])
  if (!all(is.finite(z))) {
    return(NULL)
  }
  if (stats::sd(z) <= 0) {
    return(NULL)
  }
  z
}

# ## 🟠 Create: PH null model components for Hsu ===============================

ph_null_obs_components <- function(par, dist, time, event, X) {
  p_x <- if (is.null(dim(X))) 0L else ncol(X)
  log_t <- log(time)
  
  alpha <- exp(par[1])
  lambda <- exp(par[2])
  
  if (!is.finite(alpha) || !is.finite(lambda)) return(NULL)
  if (alpha <= 0.05 || alpha >= 20) return(NULL)
  if (lambda <= 1e-8 || lambda >= 1e4) return(NULL)
  
  beta <- if (p_x > 0) par[3:(2 + p_x)] else numeric(0)
  if (length(beta) > 0 && any(!is.finite(beta))) return(NULL)
  if (length(beta) > 0 && any(abs(beta) > 20)) return(NULL)
  
  eta <- if (p_x > 0) as.vector(X %*% beta) else rep(0, length(time))
  eta <- pmax(pmin(eta, 700), -700)
  exp_eta <- exp(eta)
  
  if (dist == "weibull") {
    t_alpha <- exp(pmin(alpha * log_t, 700))
    E <- lambda * t_alpha * exp_eta
    
    logS <- -E
    ll_i <- ifelse(
      event == 1L,
      log(lambda) + log(alpha) + (alpha - 1) * log_t + eta + logS,
      logS
    )
    
    a1 <- ifelse(event == 1L, 1 + alpha * log_t - E * alpha * log_t, -E * alpha * log_t)
    a2 <- ifelse(event == 1L, 1 - E, -E)
    dls1 <- -E * alpha * log_t
    dls2 <- -E
    
    if (p_x > 0) {
      a_beta <- X * ifelse(event == 1L, 1 - E, -E)
      dls_beta <- -X * E
      a_mat <- cbind(a1, a2, a_beta)
      dlogS_mat <- cbind(dls1, dls2, dls_beta)
    } else {
      a_mat <- cbind(a1, a2)
      dlogS_mat <- cbind(dls1, dls2)
    }
  } else if (dist == "loglogistic") {
    A <- exp(pmin(log(lambda) + alpha * log_t, 700))
    H0 <- log1p(A)
    
    logS <- -exp_eta * H0
    ll_i <- ifelse(
      event == 1L,
      log(lambda) + log(alpha) + (alpha - 1) * log_t - log1p(A) + eta + logS,
      logS
    )
    
    common_a <- A * alpha * log_t / (1 + A)
    common_l <- A / (1 + A)
    
    a1 <- ifelse(
      event == 1L,
      1 + alpha * log_t - (1 + exp_eta) * common_a,
      -exp_eta * common_a
    )
    a2 <- ifelse(
      event == 1L,
      1 - (1 + exp_eta) * common_l,
      -exp_eta * common_l
    )
    dls1 <- -exp_eta * common_a
    dls2 <- -exp_eta * common_l
    
    if (p_x > 0) {
      a_beta <- X * ifelse(event == 1L, 1 - exp_eta * H0, -exp_eta * H0)
      dls_beta <- -X * (exp_eta * H0)
      a_mat <- cbind(a1, a2, a_beta)
      dlogS_mat <- cbind(dls1, dls2, dls_beta)
    } else {
      a_mat <- cbind(a1, a2)
      dlogS_mat <- cbind(dls1, dls2)
    }
  } else {
    stop("지원하지 않는 Hsu baseline입니다: ", dist)
  }
  
  if (any(!is.finite(ll_i))) return(NULL)
  logS_cap <- pmax(logS, log(.Machine$double.xmin))
  invS <- exp(pmin(-logS_cap, 700))
  
  list(
    ll_i = ll_i,
    logS = logS_cap,
    invS = invS,
    a_mat = a_mat,
    dlogS_mat = dlogS_mat,
    alpha = alpha,
    lambda = lambda,
    beta = beta
  )
}

negloglik_ph_null <- function(par, dist, time, event, X) {
  comps <- ph_null_obs_components(par = par, dist = dist, time = time, event = event, X = X)
  if (is.null(comps)) {
    return(1e50)
  }
  val <- -sum(comps$ll_i)
  if (!is.finite(val)) return(1e50)
  val
}

make_hsu_null_starts <- function(time, event, X, dist) {
  p_x <- if (is.null(dim(X))) 0L else ncol(X)
  event_time <- if (any(event == 1L)) stats::median(time[event == 1L]) else stats::median(time)
  event_time <- max(event_time, 1e-3)
  
  log_alpha_grid <- log(c(0.7, 1, 1.5))
  log_lambda_grid <- log(c(0.5, 1, 2) / event_time)
  
  start_grid <- expand.grid(
    log_alpha = log_alpha_grid,
    log_lambda = log_lambda_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (dist %in% c("weibull", "loglogistic")) {
    start_mat <- as.matrix(start_grid)
  } else {
    start_mat <- as.matrix(start_grid)
  }
  
  if (p_x > 0) {
    beta0 <- rep(0, p_x)
    start_mat <- cbind(start_mat, matrix(rep(beta0, each = nrow(start_mat)), nrow = nrow(start_mat)))
  }
  
  colnames(start_mat) <- c("log_alpha", "log_lambda", if (p_x > 0) colnames(X) else NULL)
  start_mat
}

run_multistart_optim <- function(fn, start_mat) {
  best <- NULL
  for (i in seq_len(nrow(start_mat))) {
    start_par <- as.numeric(start_mat[i, ])
    res <- tryCatch({
      nm <- optim(
        par = start_par,
        fn = fn,
        method = "Nelder-Mead",
        control = list(maxit = 5000, reltol = 1e-10)
      )
      stats::optim(
        par = nm$par,
        fn = fn,
        method = "BFGS",
        hessian = TRUE,
        control = list(maxit = 5000, reltol = 1e-10)
      )
    }, error = function(e) NULL)
    
    if (!is.null(res) && is.finite(res$value)) {
      if (is.null(best) || res$value < best$value) {
        best <- res
      }
    }
  }
  best
}

fit_hsu_null_model <- function(time, event, X, dist) {
  starts <- make_hsu_null_starts(time = time, event = event, X = X, dist = dist)
  fit <- run_multistart_optim(
    fn = function(par) negloglik_ph_null(par = par, dist = dist, time = time, event = event, X = X),
    start_mat = starts
  )
  
  if (is.null(fit)) {
    return(list(
      converged = FALSE,
      dist = dist,
      fit = NULL,
      aic = NA_real_,
      info_inv = NULL,
      components = NULL
    ))
  }
  
  hess <- fit$hessian
  if (is.null(hess) || any(!is.finite(hess))) {
    hess <- optimHess(
      par = fit$par,
      fn = function(par) negloglik_ph_null(par = par, dist = dist, time = time, event = event, X = X)
    )
  }
  
  info_inv <- safe_inverse(hess)
  comps <- ph_null_obs_components(par = fit$par, dist = dist, time = time, event = event, X = X)
  
  list(
    converged = TRUE,
    dist = dist,
    fit = fit,
    aic = 2 * length(fit$par) + 2 * fit$value,
    info_inv = info_inv,
    components = comps
  )
}

compute_hsu_supscore_from_fit <- function(time, event, z, fit_obj, theta_grid, n_boot, alpha_level) {
  n <- length(time)
  comps <- fit_obj$components
  A <- comps$a_mat
  dlogS_mat <- comps$dlogS_mat
  invS <- comps$invS
  
  W <- exp(outer(z, theta_grid, "*"))
  base_term <- ((1 - event) * invS) - 1
  bmat <- W * matrix(base_term, nrow = n, ncol = length(theta_grid))
  u_obs <- colSums(bmat)
  
  M <- W * matrix((1 - event) * invS, nrow = n, ncol = length(theta_grid))
  hmat <- crossprod(M, dlogS_mat) / n
  Q <- hmat %*% fit_obj$info_inv
  correction <- A %*% t(Q)
  cmat <- bmat - correction
  
  denom <- pmax(colSums(cmat^2), .Machine$double.eps)
  stat_grid <- ifelse(u_obs >= 0, (u_obs^2) / denom, 0)
  max_idx <- which.max(stat_grid)[1]
  T_obs <- stat_grid[max_idx]
  theta_at_max <- theta_grid[max_idx]
  
  if (n_boot <= 0) {
    return(list(
      T_obs = T_obs,
      theta_at_max = theta_at_max,
      p_value = NA_real_,
      crit_value = NA_real_,
      stat_grid = stat_grid,
      theta_grid = theta_grid,
      boot_stats = NULL
    ))
  }
  
  Xi <- matrix(stats::rnorm(n * n_boot), nrow = n, ncol = n_boot)
  U_boot <- crossprod(cmat, Xi)
  score_boot <- sweep(U_boot^2, 1, denom, "/")
  score_boot[U_boot < 0] <- 0
  T_boot <- apply(score_boot, 2, max)
  
  p_value <- (1 + sum(T_boot >= T_obs, na.rm = TRUE)) / (length(T_boot) + 1)
  crit_value <- safe_quantile(T_boot, 1 - alpha_level)
  
  list(
    T_obs = T_obs,
    theta_at_max = theta_at_max,
    p_value = p_value,
    crit_value = crit_value,
    stat_grid = stat_grid,
    theta_grid = theta_grid,
    boot_stats = if (isTRUE(config$save_bootstrap_vectors)) T_boot else NULL
  )
}

select_hsu_primary_kappa <- function(hsu_detail_one_baseline, cfg) {
  detail <- hsu_detail_one_baseline[order(hsu_detail_one_baseline$kappa), , drop = FALSE]
  if (nrow(detail) == 0) return(NA_real_)
  Tmax <- max(detail$Tn, na.rm = TRUE)
  if (!is.finite(Tmax)) return(NA_real_)
  eligible <- detail$kappa[detail$Tn >= cfg$hsu_kappa_plateau_fraction * Tmax]
  if (length(eligible) == 0) {
    return(detail$kappa[which.max(detail$Tn)[1]])
  }
  min(eligible)
}

run_hsu_supscore_group <- function(df_group, group_id, group_kind, cfg) {
  out_summary <- data.frame(
    group_id = group_id,
    hsu_primary_baseline = NA_character_,
    hsu_primary_kappa = NA_real_,
    hsu_stat = NA_real_,
    hsu_p = NA_real_,
    hsu_crit = NA_real_,
    hsu_theta_at_max = NA_real_,
    hsu_signal = NA,
    hsu_sensitivity_min_p = NA_real_,
    hsu_note = NA_character_,
    stringsAsFactors = FALSE
  )
  out_detail <- list()
  fit_store <- list()
  
  z <- build_hsu_z_vector(df_group, cfg)
  if (is.null(z)) {
    out_summary$hsu_note <- "Hsu skipped: hsu_primary_z unavailable or variance=0."
    return(list(summary = out_summary, detail = data.frame(), fits = fit_store))
  }
  
  X <- build_hsu_latency_matrix(df_group, group_kind, cfg)
  time <- df_group$time_years
  event <- df_group$event
  
  if (nrow(df_group) < cfg$min_group_n || sum(event == 1L) < cfg$min_group_events || sum(event == 0L) < cfg$min_group_censored) {
    out_summary$hsu_note <- "Hsu skipped: insufficient sample size/event/censor counts."
    return(list(summary = out_summary, detail = data.frame(), fits = fit_store))
  }
  
  baseline_fits <- lapply(cfg$hsu_primary_baseline_candidates, function(dist_name) {
    fit_hsu_null_model(time = time, event = event, X = X, dist = dist_name)
  })
  names(baseline_fits) <- cfg$hsu_primary_baseline_candidates
  
  detail_rows <- list()
  
  for (dist_name in names(baseline_fits)) {
    fit_obj <- baseline_fits[[dist_name]]
    fit_store[[dist_name]] <- fit_obj
    
    if (!isTRUE(fit_obj$converged)) {
      detail_rows[[length(detail_rows) + 1L]] <- data.frame(
        group_id = group_id,
        baseline = dist_name,
        baseline_aic = NA_real_,
        baseline_converged = FALSE,
        kappa = NA_real_,
        Tn = NA_real_,
        p_value = NA_real_,
        crit_value = NA_real_,
        theta_at_max = NA_real_,
        n_theta = NA_integer_,
        n_boot = cfg$hsu_n_resamples,
        stringsAsFactors = FALSE
      )
      next
    }
    
    for (kappa in cfg$hsu_kappa_grid) {
      theta_grid <- seq(-kappa, kappa, length.out = cfg$hsu_theta_grid_length)
      stat_obj <- compute_hsu_supscore_from_fit(
        time = time,
        event = event,
        z = z,
        fit_obj = fit_obj,
        theta_grid = theta_grid,
        n_boot = cfg$hsu_n_resamples,
        alpha_level = cfg$alpha_screen
      )
      
      fit_store[[dist_name]][[paste0("kappa_", sanitize_name(kappa))]] <- stat_obj
      
      detail_rows[[length(detail_rows) + 1L]] <- data.frame(
        group_id = group_id,
        baseline = dist_name,
        baseline_aic = fit_obj$aic,
        baseline_converged = TRUE,
        kappa = kappa,
        Tn = stat_obj$T_obs,
        p_value = stat_obj$p_value,
        crit_value = stat_obj$crit_value,
        theta_at_max = stat_obj$theta_at_max,
        n_theta = length(theta_grid),
        n_boot = cfg$hsu_n_resamples,
        stringsAsFactors = FALSE
      )
    }
  }
  
  detail_df <- rbind_fill_df(detail_rows)
  if (nrow(detail_df) == 0 || all(is.na(detail_df$baseline_aic))) {
    out_summary$hsu_note <- "Hsu failed: no successful baseline fits."
    return(list(summary = out_summary, detail = detail_df, fits = fit_store))
  }
  
  successful_baselines <- unique(detail_df$baseline[is.finite(detail_df$baseline_aic)])
  baseline_aic_df <- unique(detail_df[is.finite(detail_df$baseline_aic), c("baseline", "baseline_aic")])
  baseline_primary <- baseline_aic_df$baseline[which.min(baseline_aic_df$baseline_aic)[1]]
  
  primary_subset <- detail_df[detail_df$baseline == baseline_primary & is.finite(detail_df$Tn), , drop = FALSE]
  kappa_selected <- select_hsu_primary_kappa(primary_subset, cfg)
  primary_row <- primary_subset[which.min(abs(primary_subset$kappa - kappa_selected))[1], , drop = FALSE]
  
  out_summary$hsu_primary_baseline <- baseline_primary
  out_summary$hsu_primary_kappa <- kappa_selected
  out_summary$hsu_stat <- primary_row$Tn
  out_summary$hsu_p <- primary_row$p_value
  out_summary$hsu_crit <- primary_row$crit_value
  out_summary$hsu_theta_at_max <- primary_row$theta_at_max
  out_summary$hsu_signal <- ifelse(is.na(primary_row$p_value), NA, primary_row$p_value < cfg$alpha_screen)
  out_summary$hsu_sensitivity_min_p <- suppressWarnings(min(detail_df$p_value, na.rm = TRUE))
  out_summary$hsu_note <- sprintf(
    "Primary baseline=%s selected by null-model AIC; primary kappa chosen as the smallest kappa reaching %.0f%% of max Tn for that baseline.",
    baseline_primary, 100 * cfg$hsu_kappa_plateau_fraction
  )
  
  list(summary = out_summary, detail = detail_df, fits = fit_store)
}

# 🔴 Implement: strict sufficient follow-up test by Xie-style extremes ===============================

## 🟠 Create: Xie test statistic and bootstrap routine ===============================

run_xie_stat_once <- function(time, event) {
  n <- length(time)
  if (!any(event == 1L)) {
    return(list(
      p_hat_n = NA_real_,
      p_hat_g = NA_real_,
      epsilon = NA_real_,
      Tn = NA_real_,
      q_n = NA_real_,
      alpha_n = NA_real_,
      t_n = max(time, na.rm = TRUE),
      t_k = NA_real_
    ))
  }
  
  t_n <- max(time, na.rm = TRUE)
  t_k <- max(time[event == 1L], na.rm = TRUE)
  
  if (2 * (t_n - t_k) < t_n) {
    epsilon <- (9 / 8) * t_n - (1 / 4) * t_k
  } else {
    epsilon <- t_n
  }
  
  t1 <- max(t_n - epsilon / 2, 0)
  t2 <- max(t_n - epsilon, 0)
  
  surv_vals <- km_surv_at_times(time, event, c(t_n, t1, t2))
  F_tn <- 1 - surv_vals[1]
  F_t1 <- 1 - surv_vals[2]
  F_t2 <- 1 - surv_vals[3]
  
  denom <- 2 * F_t1 - F_t2 - F_tn
  p_hat_n <- F_tn
  
  if (!is.finite(denom) || abs(denom) < 1e-12) {
    p_hat_g <- p_hat_n
  } else {
    p_hat_g <- F_t2 + (F_t1 - F_t2)^2 / denom
  }
  
  if (!is.finite(p_hat_g) || p_hat_g < p_hat_n) {
    p_hat_g <- p_hat_n
  }
  
  p_hat_g <- pmin(pmax(p_hat_g, p_hat_n), 1)
  Tn <- p_hat_g - p_hat_n
  
  interval_lower <- 2 * t_k - t_n
  Nn <- sum(event == 1L & time > interval_lower & time <= t_k)
  q_n <- Nn / n
  alpha_n <- (1 - Nn / n)^n
  
  list(
    p_hat_n = p_hat_n,
    p_hat_g = p_hat_g,
    epsilon = epsilon,
    Tn = Tn,
    q_n = q_n,
    alpha_n = alpha_n,
    t_n = t_n,
    t_k = t_k
  )
}

run_xie_group <- function(df_group, group_id, cfg) {
  time <- df_group$time_years
  event <- df_group$event
  
  out <- data.frame(
    group_id = group_id,
    xie_t_n = NA_real_,
    xie_t_k = NA_real_,
    xie_epsilon = NA_real_,
    xie_p_hat_n = NA_real_,
    xie_p_hat_g = NA_real_,
    xie_stat = NA_real_,
    xie_p = NA_real_,
    xie_crit = NA_real_,
    xie_flag = NA_character_,
    legacy_q_n = NA_real_,
    legacy_alpha_n = NA_real_,
    xie_note = NA_character_,
    stringsAsFactors = FALSE
  )
  
  if (nrow(df_group) < cfg$min_group_n || sum(event == 1L) < cfg$min_group_events || sum(event == 0L) < cfg$min_group_censored) {
    out$xie_note <- "Xie skipped: insufficient sample size/event/censor counts."
    return(list(summary = out, boot_stats = NULL))
  }
  
  obs <- run_xie_stat_once(time = time, event = event)
  if (!is.finite(obs$Tn)) {
    out$xie_note <- "Xie skipped: no uncensored event time available."
    return(list(summary = out, boot_stats = NULL))
  }
  
  T_boot <- rep(NA_real_, cfg$xie_n_boot)
  for (b in seq_len(cfg$xie_n_boot)) {
    idx <- sample.int(n = nrow(df_group), size = nrow(df_group), replace = TRUE)
    boot_stat <- run_xie_stat_once(time = time[idx], event = event[idx])
    T_boot[b] <- boot_stat$Tn
  }
  
  centered_boot <- T_boot - obs$Tn
  crit <- safe_quantile(centered_boot, 1 - cfg$alpha_screen)
  pval <- (1 + sum(centered_boot >= obs$Tn, na.rm = TRUE)) / (sum(is.finite(centered_boot)) + 1)
  
  out$xie_t_n <- obs$t_n
  out$xie_t_k <- obs$t_k
  out$xie_epsilon <- obs$epsilon
  out$xie_p_hat_n <- obs$p_hat_n
  out$xie_p_hat_g <- obs$p_hat_g
  out$xie_stat <- obs$Tn
  out$xie_p <- pval
  out$xie_crit <- crit
  out$xie_flag <- ifelse(
    is.na(pval),
    NA_character_,
    ifelse(pval < cfg$alpha_screen, "evidence_insufficient_followup", "fail_to_reject_sufficient_followup")
  )
  out$legacy_q_n <- obs$q_n
  out$legacy_alpha_n <- obs$alpha_n
  out$xie_note <- "Primary strict sufficient follow-up test using Xie-style extremes; centered naive bootstrap p-value."
  
  list(summary = out, boot_stats = if (isTRUE(cfg$save_bootstrap_vectors)) T_boot else NULL)
}

# 🔴 Implement: RECeUS-AIC screening machinery ===============================

## 🟠 Create: candidate density and survival kernels ===============================

dist_logf_logS <- function(dist, time, par) {
  if (dist == "exponential") {
    scale <- exp(par[1])
    if (!is.finite(scale) || scale <= 1e-8 || scale >= 1e4) return(NULL)
    rate <- 1 / scale
    logf <- stats::dexp(time, rate = rate, log = TRUE)
    logS <- stats::pexp(time, rate = rate, lower.tail = FALSE, log.p = TRUE)
    return(list(
      logf = logf,
      logS = logS,
      shape = NA_real_,
      scale = scale,
      rate = rate
    ))
  }
  
  if (dist == "weibull") {
    shape <- exp(par[1])
    scale <- exp(par[2])
    if (!is.finite(shape) || !is.finite(scale) || shape <= 0.05 || shape >= 30 || scale <= 1e-8 || scale >= 1e4) return(NULL)
    logf <- stats::dweibull(time, shape = shape, scale = scale, log = TRUE)
    logS <- stats::pweibull(time, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
    return(list(
      logf = logf,
      logS = logS,
      shape = shape,
      scale = scale,
      rate = NA_real_
    ))
  }
  
  if (dist == "gamma") {
    shape <- exp(par[1])
    scale <- exp(par[2])
    if (!is.finite(shape) || !is.finite(scale) || shape <= 0.05 || shape >= 30 || scale <= 1e-8 || scale >= 1e4) return(NULL)
    logf <- stats::dgamma(time, shape = shape, scale = scale, log = TRUE)
    logS <- stats::pgamma(time, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
    return(list(
      logf = logf,
      logS = logS,
      shape = shape,
      scale = scale,
      rate = NA_real_
    ))
  }
  
  if (dist == "loglogistic") {
    shape <- exp(par[1])
    scale <- exp(par[2])
    if (!is.finite(shape) || !is.finite(scale) || shape <= 0.05 || shape >= 30 || scale <= 1e-8 || scale >= 1e4) return(NULL)
    A <- (time / scale)^shape
    logS <- -log1p(A)
    logf <- log(shape) - log(scale) + (shape - 1) * (log(time) - log(scale)) - 2 * log1p(A)
    return(list(
      logf = logf,
      logS = logS,
      shape = shape,
      scale = scale,
      rate = NA_real_
    ))
  }
  
  stop("지원하지 않는 분포입니다: ", dist)
}

negloglik_receus_nocure <- function(par, dist, time, event) {
  comps <- dist_logf_logS(dist = dist, time = time, par = par)
  if (is.null(comps)) return(1e50)
  ll_i <- ifelse(event == 1L, comps$logf, comps$logS)
  val <- -sum(ll_i)
  if (!is.finite(val)) return(1e50)
  val
}

negloglik_receus_cure <- function(par, dist, time, event) {
  pi_cure <- inv_logit(par[length(par)])
  pi_cure <- clamp_prob(pi_cure)
  base_par <- par[-length(par)]
  comps <- dist_logf_logS(dist = dist, time = time, par = base_par)
  if (is.null(comps)) return(1e50)
  
  log_pi <- log(pi_cure)
  log1m_pi <- log1p(-pi_cure)
  log_mix_surv <- logspace_add(log_pi, log1m_pi + comps$logS)
  
  ll_i <- ifelse(
    event == 1L,
    log1m_pi + comps$logf,
    log_mix_surv
  )
  
  val <- -sum(ll_i)
  if (!is.finite(val)) return(1e50)
  val
}

make_receus_start_grid <- function(time, event, dist, mixture, km_tail_surv) {
  event_time <- if (any(event == 1L)) stats::median(time[event == 1L]) else stats::median(time)
  event_time <- max(event_time, 1e-3)
  mean_time <- max(mean(time), 1e-3)
  
  if (dist == "exponential") {
    base_scale <- max(sum(time) / max(sum(event), 1), 1e-3)
    base_grid <- matrix(log(base_scale * c(0.5, 1, 2)), ncol = 1)
    colnames(base_grid) <- "log_scale"
  } else {
    log_shape_grid <- log(c(0.7, 1, 1.5))
    log_scale_grid <- log(c(event_time, mean_time, stats::quantile(time, 0.75, na.rm = TRUE)))
    base_grid <- expand.grid(
      log_shape = log_shape_grid,
      log_scale = as.numeric(log_scale_grid),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    base_grid <- as.matrix(unique(base_grid))
    colnames(base_grid) <- c("log_shape", "log_scale")
  }
  
  if (!mixture) {
    return(base_grid)
  }
  
  pi_start <- clamp_prob(km_tail_surv, eps = 0.05)
  cure_grid <- qlogis(clamp_prob(c(pi_start, 0.10, 0.25, 0.40), eps = 0.01))
  start_list <- lapply(seq_len(nrow(base_grid)), function(i) {
    cbind(
      matrix(rep(base_grid[i, ], each = length(cure_grid)), ncol = ncol(base_grid), byrow = FALSE),
      cure_logit = cure_grid
    )
  })
  start_mat <- do.call(rbind, start_list)
  colnames(start_mat) <- c(colnames(base_grid), "cure_logit")
  start_mat
}

fit_receus_candidate <- function(time, event, dist, mixture, km_tail_surv) {
  starts <- make_receus_start_grid(time = time, event = event, dist = dist, mixture = mixture, km_tail_surv = km_tail_surv)
  fn <- if (mixture) {
    function(par) negloglik_receus_cure(par = par, dist = dist, time = time, event = event)
  } else {
    function(par) negloglik_receus_nocure(par = par, dist = dist, time = time, event = event)
  }
  
  fit <- run_multistart_optim(fn = fn, start_mat = starts)
  
  if (is.null(fit)) {
    return(list(
      converged = FALSE,
      model_class = if (mixture) "mixture_cure" else "nocure",
      family = dist
    ))
  }
  
  loglik <- -fit$value
  k <- length(fit$par)
  aic <- 2 * k - 2 * loglik
  bic <- log(length(time)) * k - 2 * loglik
  
  if (mixture) {
    pi_cure <- inv_logit(fit$par[length(fit$par)])
    pi_cure <- clamp_prob(pi_cure)
    base_par <- fit$par[-length(fit$par)]
  } else {
    pi_cure <- 0
    base_par <- fit$par
  }
  
  comps <- dist_logf_logS(dist = dist, time = time, par = base_par)
  tau_admin <- calc_tau_admin(time, event)
  tau_surv <- dist_logf_logS(dist = dist, time = tau_admin, par = base_par)
  
  if (is.null(comps) || is.null(tau_surv)) {
    return(list(
      converged = FALSE,
      model_class = if (mixture) "mixture_cure" else "nocure",
      family = dist
    ))
  }
  
  S_u_tau <- exp(tau_surv$logS)
  S_pop_tau <- if (mixture) {
    pi_cure + (1 - pi_cure) * S_u_tau
  } else {
    S_u_tau
  }
  
  r_hat <- if (mixture) S_u_tau / S_pop_tau else NA_real_
  
  list(
    converged = TRUE,
    fit = fit,
    model_class = if (mixture) "mixture_cure" else "nocure",
    family = dist,
    logLik = loglik,
    AIC = aic,
    BIC = bic,
    pi_cure_hat = pi_cure,
    r_hat = r_hat,
    tau_admin = tau_admin,
    shape_hat = comps$shape,
    scale_hat = comps$scale,
    rate_hat = comps$rate
  )
}

run_receus_group <- function(df_group, group_id, cfg) {
  time <- df_group$time_years
  event <- df_group$event
  
  out_summary <- data.frame(
    group_id = group_id,
    receus_best_model = NA_character_,
    receus_best_model_class = NA_character_,
    receus_aic = NA_real_,
    receus_bic = NA_real_,
    receus_pi_cure = NA_real_,
    receus_r = NA_real_,
    receus_decision = NA,
    receus_reason = NA_character_,
    stringsAsFactors = FALSE
  )
  
  if (nrow(df_group) < cfg$min_group_n || sum(event == 1L) < cfg$min_group_events || sum(event == 0L) < cfg$min_group_censored) {
    out_summary$receus_reason <- "RECeUS skipped: insufficient sample size/event/censor counts."
    return(list(summary = out_summary, detail = data.frame(), fits = list()))
  }
  
  tau_admin <- calc_tau_admin(time, event)
  km_tail_surv <- km_surv_at_times(time, event, tau_admin)
  
  candidate_rows <- list()
  fit_store <- list()
  
  for (fam in cfg$receus_candidate_families) {
    fit_nc <- fit_receus_candidate(time = time, event = event, dist = fam, mixture = FALSE, km_tail_surv = km_tail_surv)
    fit_mc <- fit_receus_candidate(time = time, event = event, dist = fam, mixture = TRUE, km_tail_surv = km_tail_surv)
    
    fit_store[[paste0("nocure_", fam)]] <- fit_nc
    fit_store[[paste0("cure_", fam)]] <- fit_mc
    
    candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
      group_id = group_id,
      model_class = fit_nc$model_class,
      family = fit_nc$family,
      candidate_id = paste(fit_nc$model_class, fit_nc$family, sep = "::"),
      converged = isTRUE(fit_nc$converged),
      logLik = ifelse(isTRUE(fit_nc$converged), fit_nc$logLik, NA_real_),
      AIC = ifelse(isTRUE(fit_nc$converged), fit_nc$AIC, NA_real_),
      BIC = ifelse(isTRUE(fit_nc$converged), fit_nc$BIC, NA_real_),
      pi_cure_hat = ifelse(isTRUE(fit_nc$converged), fit_nc$pi_cure_hat, NA_real_),
      r_hat = ifelse(isTRUE(fit_nc$converged), fit_nc$r_hat, NA_real_),
      tau_admin = ifelse(isTRUE(fit_nc$converged), fit_nc$tau_admin, NA_real_),
      shape_hat = ifelse(isTRUE(fit_nc$converged), fit_nc$shape_hat, NA_real_),
      scale_hat = ifelse(isTRUE(fit_nc$converged), fit_nc$scale_hat, NA_real_),
      rate_hat = ifelse(isTRUE(fit_nc$converged), fit_nc$rate_hat, NA_real_),
      stringsAsFactors = FALSE
    )
    
    candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
      group_id = group_id,
      model_class = fit_mc$model_class,
      family = fit_mc$family,
      candidate_id = paste(fit_mc$model_class, fit_mc$family, sep = "::"),
      converged = isTRUE(fit_mc$converged),
      logLik = ifelse(isTRUE(fit_mc$converged), fit_mc$logLik, NA_real_),
      AIC = ifelse(isTRUE(fit_mc$converged), fit_mc$AIC, NA_real_),
      BIC = ifelse(isTRUE(fit_mc$converged), fit_mc$BIC, NA_real_),
      pi_cure_hat = ifelse(isTRUE(fit_mc$converged), fit_mc$pi_cure_hat, NA_real_),
      r_hat = ifelse(isTRUE(fit_mc$converged), fit_mc$r_hat, NA_real_),
      tau_admin = ifelse(isTRUE(fit_mc$converged), fit_mc$tau_admin, NA_real_),
      shape_hat = ifelse(isTRUE(fit_mc$converged), fit_mc$shape_hat, NA_real_),
      scale_hat = ifelse(isTRUE(fit_mc$converged), fit_mc$scale_hat, NA_real_),
      rate_hat = ifelse(isTRUE(fit_mc$converged), fit_mc$rate_hat, NA_real_),
      stringsAsFactors = FALSE
    )
  }
  
  detail_df <- rbind_fill_df(candidate_rows)
  detail_ok <- detail_df[detail_df$converged & is.finite(detail_df$AIC), , drop = FALSE]
  
  if (nrow(detail_ok) == 0) {
    out_summary$receus_reason <- "RECeUS failed: no candidate model converged."
    return(list(summary = out_summary, detail = detail_df, fits = fit_store))
  }
  
  detail_ok <- detail_ok[order(detail_ok$AIC, detail_ok$BIC, detail_ok$model_class), , drop = FALSE]
  best_row <- detail_ok[1, , drop = FALSE]
  
  out_summary$receus_best_model <- paste(best_row$model_class, best_row$family, sep = "::")
  out_summary$receus_best_model_class <- best_row$model_class
  out_summary$receus_aic <- best_row$AIC
  out_summary$receus_bic <- best_row$BIC
  out_summary$receus_pi_cure <- best_row$pi_cure_hat
  out_summary$receus_r <- best_row$r_hat
  
  if (best_row$model_class == "nocure") {
    out_summary$receus_decision <- FALSE
    out_summary$receus_reason <- "AIC selected a no-cure model."
  } else if (best_row$pi_cure_hat <= cfg$receus_threshold_pi_cure) {
    out_summary$receus_decision <- FALSE
    out_summary$receus_reason <- sprintf("Mixture cure selected, but pi_cure_hat <= %.3f.", cfg$receus_threshold_pi_cure)
  } else if (!is.finite(best_row$r_hat) || best_row$r_hat >= cfg$receus_threshold_r) {
    out_summary$receus_decision <- FALSE
    out_summary$receus_reason <- sprintf("Mixture cure selected, but r_hat >= %.3f.", cfg$receus_threshold_r)
  } else {
    out_summary$receus_decision <- TRUE
    out_summary$receus_reason <- "Mixture cure selected and both RECeUS thresholds satisfied."
  }
  
  list(summary = out_summary, detail = detail_df, fits = fit_store)
}

# 🔴 Execute: group-wise Step3 screening loop ===============================

## 🟠 Iterate: screening methods over merged-overall and site-specific groups ===============================

group_summary_rows <- list()
hsu_detail_rows <- list()
receus_detail_rows <- list()
followup_detail_rows <- list()
group_objects <- list()

for (grp in analysis_groups) {
  gid <- grp$group_id
  gkind <- grp$group_kind
  gdisp <- grp$group_display
  gsite <- grp$site_label
  gdat <- grp$data
  
  message("Running Step 3 screening: ", gid)
  
  basic_fu <- compute_followup_context(gdat)
  
  basic_row <- data.frame(
    group_id = gid,
    group_display = gdisp,
    group_kind = gkind,
    site_label = gsite,
    n = nrow(gdat),
    n_event = sum(gdat$event == 1L),
    n_censored = sum(gdat$event == 0L),
    n_remission_censor = sum(gdat$status_num == 2L),
    n_right_censor = sum(gdat$status_num == 0L),
    censor_rate = mean(gdat$event == 0L),
    tau_max_years = basic_fu$tau_max,
    tau_admin_years = basic_fu$tau_admin,
    last_event_years = basic_fu$last_event,
    plateau_length_years = basic_fu$plateau_length,
    n_censored_after_last_event = basic_fu$n_censored_after_last_event,
    km_tail_survival = basic_fu$km_tail_surv,
    reverseKM_median_followup_years = basic_fu$reverse_km_median_followup,
    cci = basic_fu$cci,
    spt = basic_fu$spt,
    stringsAsFactors = FALSE
  )
  
  hsu_res <- run_hsu_supscore_group(df_group = gdat, group_id = gid, group_kind = gkind, cfg = config)
  xie_res <- run_xie_group(df_group = gdat, group_id = gid, cfg = config)
  receus_res <- run_receus_group(df_group = gdat, group_id = gid, cfg = config)
  
  practical_row <- data.frame(
    group_id = gid,
    practical_run = FALSE,
    practical_stat = NA_real_,
    practical_p = NA_real_,
    practical_flag = "not_run_by_design",
    stringsAsFactors = FALSE
  )
  
  merged_summary <- Reduce(
    function(x, y) merge(x, y, by = "group_id", all = TRUE),
    list(
      basic_row,
      hsu_res$summary,
      xie_res$summary,
      practical_row,
      receus_res$summary
    )
  )
  
  merged_summary$notes <- paste(
    na.omit(unique(c(
      merged_summary$hsu_note,
      merged_summary$xie_note,
      merged_summary$receus_reason
    ))),
    collapse = " | "
  )
  
  group_summary_rows[[length(group_summary_rows) + 1L]] <- merged_summary
  hsu_detail_rows[[length(hsu_detail_rows) + 1L]] <- hsu_res$detail
  receus_detail_rows[[length(receus_detail_rows) + 1L]] <- receus_res$detail
  followup_detail_rows[[length(followup_detail_rows) + 1L]] <- Reduce(
    function(x, y) merge(x, y, by = "group_id", all = TRUE),
    list(
      basic_row[, c(
        "group_id",
        "tau_max_years",
        "tau_admin_years",
        "last_event_years",
        "plateau_length_years",
        "n_censored_after_last_event",
        "km_tail_survival",
        "reverseKM_median_followup_years",
        "cci",
        "spt"
      )],
      xie_res$summary,
      practical_row
    )
  )
  
  group_objects[[gid]] <- list(
    basic = basic_row,
    hsu = hsu_res,
    xie = xie_res,
    practical = list(summary = practical_row, note = "Optional practical sufficient follow-up test not run."),
    receus = receus_res
  )
}

screening_summary_df <- rbind_fill_df(group_summary_rows)
hsu_detail_df <- rbind_fill_df(hsu_detail_rows)
receus_detail_df <- rbind_fill_df(receus_detail_rows)
followup_detail_df <- rbind_fill_df(followup_detail_rows)

# 🔴 Aggregate: grouped adaptation summary for merged interpretation ===============================

## 🟠 Derive: at-least-one-site and pooled-vs-site rules ===============================

site_rows <- screening_summary_df[screening_summary_df$group_kind == "site_specific", , drop = FALSE]
merged_row <- screening_summary_df[screening_summary_df$group_kind == "merged_overall", , drop = FALSE]

grouped_summary_df <- data.frame(
  rule_id = c(
    "receus_at_least_one_site",
    "receus_all_sites",
    "hsu_at_least_one_site",
    "xie_any_site_evidence_insufficient",
    "xie_all_sites_fail_to_reject_insufficient",
    "merged_receus_primary",
    "merged_hsu_primary_signal",
    "merged_xie_primary_flag"
  ),
  value = c(
    if (nrow(site_rows) > 0) any(site_rows$receus_decision %in% TRUE, na.rm = TRUE) else NA,
    if (nrow(site_rows) > 0) all(site_rows$receus_decision %in% TRUE, na.rm = TRUE) else NA,
    if (nrow(site_rows) > 0) any(site_rows$hsu_signal %in% TRUE, na.rm = TRUE) else NA,
    if (nrow(site_rows) > 0) any(site_rows$xie_p < config$alpha_screen, na.rm = TRUE) else NA,
    if (nrow(site_rows) > 0) all(site_rows$xie_p >= config$alpha_screen, na.rm = TRUE) else NA,
    if (nrow(merged_row) > 0) merged_row$receus_decision[1] else NA,
    if (nrow(merged_row) > 0) merged_row$hsu_signal[1] else NA,
    if (nrow(merged_row) > 0) merged_row$xie_flag[1] else NA
  ),
  note = c(
    "Kouadio-style grouped adaptation: at least one site meets RECeUS criteria.",
    "All site-specific groups meet RECeUS criteria.",
    "At least one site-specific group has Hsu primary p < alpha.",
    "At least one site-specific group shows evidence against sufficient follow-up by Xie.",
    "All site-specific groups fail to reject insufficient follow-up by Xie.",
    "Secondary pooled merged-overall RECeUS conclusion.",
    "Secondary pooled merged-overall Hsu conclusion.",
    "Secondary pooled merged-overall Xie conclusion."
  ),
  stringsAsFactors = FALSE
)

if (nrow(merged_row) > 0) {
  screening_summary_df$receus_at_least_one_site <- ifelse(
    screening_summary_df$group_kind == "merged_overall",
    grouped_summary_df$value[grouped_summary_df$rule_id == "receus_at_least_one_site"],
    NA
  )
  screening_summary_df$hsu_at_least_one_site <- ifelse(
    screening_summary_df$group_kind == "merged_overall",
    grouped_summary_df$value[grouped_summary_df$rule_id == "hsu_at_least_one_site"],
    NA
  )
  screening_summary_df$xie_any_site_evidence_insufficient <- ifelse(
    screening_summary_df$group_kind == "merged_overall",
    grouped_summary_df$value[grouped_summary_df$rule_id == "xie_any_site_evidence_insufficient"],
    NA
  )
}

# 🔴 Export: step3 tables, manifest, and reusable objects ===============================

## 🟠 Write: source-of-truth CSV tables ===============================

screening_summary_file <- file.path(export_path, "step3_screening_summary.csv")
hsu_detail_file <- file.path(export_path, "step3_hsu_detail.csv")
receus_detail_file <- file.path(export_path, "step3_receus_candidates.csv")
followup_detail_file <- file.path(export_path, "step3_followup_detail.csv")
grouped_summary_file <- file.path(export_path, "step3_grouped_summary.csv")
objects_rds_file <- file.path(export_path, "step3_objects.rds")
manifest_file <- file.path(export_path, "step3_export_manifest.csv")

utils::write.csv(screening_summary_df, screening_summary_file, row.names = FALSE, na = "")
utils::write.csv(hsu_detail_df, hsu_detail_file, row.names = FALSE, na = "")
utils::write.csv(receus_detail_df, receus_detail_file, row.names = FALSE, na = "")
utils::write.csv(followup_detail_df, followup_detail_file, row.names = FALSE, na = "")
utils::write.csv(grouped_summary_df, grouped_summary_file, row.names = FALSE, na = "")

## 🟠 Save: reusable fitted objects as RDS ===============================

step3_object <- list(
  config = config,
  preproc_spec = preproc_spec,
  site_levels = site_levels,
  analysis_data = analysis_df,
  analysis_groups = lapply(analysis_groups, function(x) x[c("group_id", "group_display", "group_kind", "site_label")]),
  results = group_objects,
  exports = list(
    screening_summary_file = screening_summary_file,
    hsu_detail_file = hsu_detail_file,
    receus_detail_file = receus_detail_file,
    followup_detail_file = followup_detail_file,
    grouped_summary_file = grouped_summary_file
  ),
  session_info = utils::capture.output(sessionInfo())
)

saveRDS(step3_object, objects_rds_file)

## 🟠 Record: manifest for downstream reuse ===============================

export_manifest_df <- data.frame(
  file_role = c(
    "screening_summary",
    "hsu_detail",
    "receus_candidates",
    "followup_detail",
    "grouped_summary",
    "objects_rds"
  ),
  file_name = basename(c(
    screening_summary_file,
    hsu_detail_file,
    receus_detail_file,
    followup_detail_file,
    grouped_summary_file,
    objects_rds_file
  )),
  absolute_path = c(
    screening_summary_file,
    hsu_detail_file,
    receus_detail_file,
    followup_detail_file,
    grouped_summary_file,
    objects_rds_file
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(export_manifest_df, manifest_file, row.names = FALSE, na = "")

## 🟠 Check: exported file count ceiling ===============================

exported_files <- c(
  screening_summary_file,
  hsu_detail_file,
  receus_detail_file,
  followup_detail_file,
  grouped_summary_file,
  objects_rds_file,
  manifest_file
)

if (length(exported_files) > 20) {
  stop("생성 파일 수가 20개를 초과했습니다: ", length(exported_files))
}

message("Step 3 screening complete.")
message("Exported files:")
for (fp in exported_files) {
  message(" - ", fp)
}