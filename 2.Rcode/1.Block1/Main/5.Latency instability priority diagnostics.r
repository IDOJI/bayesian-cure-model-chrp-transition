# Configure -----------------------------------------------------------------
find_repo_root <- function(start_dir) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)

  repeat {
    has_repo_markers <- dir.exists(file.path(current_dir, ".git")) ||
      (
        dir.exists(file.path(current_dir, "0.Data")) &&
          dir.exists(file.path(current_dir, "2.Rcode")) &&
          dir.exists(file.path(current_dir, "3.Results files"))
      )

    if (has_repo_markers) {
      return(current_dir)
    }

    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      break
    }
    current_dir <- parent_dir
  }

  stop("Could not locate the repository root.", call. = FALSE)
}

command_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", command_args, value = TRUE)
search_start_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

repo_root <- find_repo_root(search_start_dir)
dropbox_root_default <- "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model"

analysis_dataset_file <- Sys.getenv(
  "LATENCY_DIAG_ANALYSIS_DATASETS_RDS",
  unset = file.path(dropbox_root_default, "2.Block2", "stage1_analysis_datasets.rds")
)
frequentist_fit_cache_file <- Sys.getenv(
  "LATENCY_DIAG_FREQUENTIST_FIT_CACHE",
  unset = file.path(dropbox_root_default, "1.Block1", "3.MLE Mixture Cure", "mle_mixture_cure_lognormal_fitted_objects.rds")
)
bayesian_model_registry_file <- Sys.getenv(
  "LATENCY_DIAG_BAYESIAN_MODEL_REGISTRY",
  unset = file.path(dropbox_root_default, "1.Block1", "4.Bayesian Mixture Cure", "bayesian_mixture_cure_model_registry.csv")
)
output_dir <- Sys.getenv(
  "LATENCY_DIAG_OUTPUT_DIR",
  unset = file.path(repo_root, "3.Results files", "3.Block3", "latency_instability_priority")
)

tail_n <- as.integer(Sys.getenv("LATENCY_DIAG_TAIL_N", unset = "5"))
joint_draws_max <- as.integer(Sys.getenv("LATENCY_DIAG_JOINT_DRAWS_MAX", unset = "1000"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(rstan)
  library(scales)
  library(survival)
  library(tibble)
  library(tidyr)
})

options(stringsAsFactors = FALSE, scipen = 999)

# Helpers -------------------------------------------------------------------
`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}

assert_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s does not exist: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

trim_rhs <- function(rhs) {
  rhs <- gsub("\\s+", " ", trimws(as.character(rhs)))
  rhs[is.na(rhs)] <- NA_character_
  rhs
}

clamp_prob <- function(x, eps = 1e-12) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

weighted_mean_safe <- function(x, w) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  sum(x[ok] * w[ok]) / sum(w[ok])
}

capture_with_warnings <- function(expr) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  list(
    value = if (inherits(value, "error")) NULL else value,
    warnings = warnings,
    error_message = if (inherits(value, "error")) conditionMessage(value) else NA_character_
  )
}

make_rhs_formula <- function(rhs) {
  rhs <- trim_rhs(rhs)
  if (is.na(rhs) || rhs == "") {
    stats::as.formula("~ 1")
  } else {
    stats::as.formula(paste("~", rhs))
  }
}

make_surv_formula <- function(rhs) {
  rhs <- trim_rhs(rhs)
  if (is.na(rhs) || rhs == "") {
    stats::as.formula("survival::Surv(time_year_model, event_main) ~ 1")
  } else {
    stats::as.formula(paste("survival::Surv(time_year_model, event_main) ~", rhs))
  }
}

model_matrix_from_rhs <- function(df, rhs) {
  stats::model.matrix(make_rhs_formula(rhs), data = df)
}

make_psd_matrix <- function(x, eigen_floor = 1e-10) {
  eig <- eigen(x, symmetric = TRUE)
  eig$values[eig$values < eigen_floor] <- eigen_floor
  eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
}

invert_hessian <- function(hessian_mat, eigen_floor = 1e-8) {
  sym_hessian <- 0.5 * (hessian_mat + t(hessian_mat))
  eig <- eigen(sym_hessian, symmetric = TRUE)
  eig$values[eig$values < eigen_floor] <- eigen_floor
  eig$vectors %*% diag(1 / eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
}

draw_mvn <- function(n_draws, mu, sigma, seed_value) {
  sigma_psd <- make_psd_matrix(sigma)
  eig <- eigen(sigma_psd, symmetric = TRUE)
  transform_mat <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow = length(eig$values))

  set.seed(seed_value)
  z <- matrix(stats::rnorm(length(mu) * n_draws), nrow = length(mu), ncol = n_draws)
  draws <- matrix(mu, nrow = length(mu), ncol = n_draws) + transform_mat %*% z
  t(draws)
}

ensure_draw_matrix <- function(x, ncol_expected) {
  if (is.null(dim(x))) {
    return(matrix(as.numeric(x), ncol = ncol_expected))
  }
  matrix(
    as.numeric(x),
    nrow = dim(x)[1L],
    ncol = ncol_expected,
    dimnames = dimnames(x)
  )
}

safe_quantile <- function(x, prob) {
  x <- as.numeric(x)
  if (!any(is.finite(x))) {
    return(NA_real_)
  }
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 8))
}

summarize_vector <- function(x) {
  tibble(
    mean = mean(x, na.rm = TRUE),
    sd = if (length(x) <= 1L) 0 else stats::sd(x, na.rm = TRUE),
    q025 = safe_quantile(x, 0.025),
    q50 = safe_quantile(x, 0.50),
    q975 = safe_quantile(x, 0.975)
  )
}

format_removed_ids <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) {
    return(NA_character_)
  }
  paste(x, collapse = " | ")
}

# Data preparation -----------------------------------------------------------
prepare_dataset <- function(df, dataset_name) {
  required_cols <- c("site", "id", "sex_num", "age_exact_entry", "age_s", "status_num", "time_year")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("[%s] Dataset is missing required columns: %s", dataset_name, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  out <- tibble::as_tibble(df)

  if (!("unique_person_id" %in% names(out))) {
    out <- out %>%
      mutate(unique_person_id = paste(site, id, sep = "_"))
  }

  out <- out %>%
    mutate(
      unique_person_id = as.character(unique_person_id),
      site = factor(as.character(site)),
      id = as.character(id),
      sex_num = as.integer(as.numeric(sex_num)),
      age_exact_entry = as.numeric(age_exact_entry),
      age_s = as.numeric(age_s),
      status_num = as.integer(as.numeric(status_num)),
      time_year = as.numeric(time_year),
      time_year_model = pmax(as.numeric(time_year), 1e-10),
      event_main = as.integer(status_num == 1L),
      right_censor_flag = as.integer(status_num == 0L),
      remission_flag = as.integer(status_num == 2L),
      censor_main = as.integer(status_num %in% c(0L, 2L)),
      censor_detail = dplyr::case_when(
        status_num == 2L ~ "remission",
        status_num == 0L ~ "right_censoring",
        status_num == 1L ~ "transition",
        TRUE ~ "other"
      )
    )

  if (nrow(out) != dplyr::n_distinct(out$unique_person_id)) {
    stop(sprintf("[%s] `unique_person_id` must be unique.", dataset_name), call. = FALSE)
  }

  out
}

frequentist_spec_registry <- tibble(
  dataset = c("PNU", "SNU", "merged_no_site", "merged_site_adjusted"),
  dataset_label = c("PNU", "SNU", "Merged", "Merged (site-adjusted)"),
  source_dataset = c("PNU", "SNU", "merged", "merged"),
  formula_name = c("base", "base", "base", "site_added"),
  incidence_rhs = c("age_s + sex_num", "age_s + sex_num", "age_s + sex_num", "age_s + sex_num + site"),
  latency_rhs = c("age_s + sex_num", "age_s + sex_num", "age_s + sex_num", "age_s + sex_num + site"),
  site_adjustment_flag = c(FALSE, FALSE, FALSE, TRUE)
)

# Parametric cure fitting ----------------------------------------------------
survreg_family_map <- c(
  lognormal = "lognormal",
  weibull = "weibull",
  loglogistic = "loglogistic",
  exponential = "exponential"
)

family_has_aux_scale <- function(family_name) {
  !identical(as.character(family_name), "exponential")
}

get_family_scale <- function(family_name, log_scale_value = NA_real_) {
  if (family_has_aux_scale(family_name)) {
    exp(as.numeric(log_scale_value))
  } else {
    1
  }
}

aft_surv_density <- function(time, lp, family_name, log_scale_value = NA_real_) {
  dist_name <- unname(survreg_family_map[[family_name]])
  if (is.na(dist_name)) {
    stop(sprintf("Unsupported family: %s", family_name), call. = FALSE)
  }

  time <- pmax(as.numeric(time), 1e-10)
  lp <- as.numeric(lp)
  scale_value <- get_family_scale(family_name, log_scale_value)

  surv <- 1 - survival::psurvreg(
    q = time,
    mean = lp,
    scale = scale_value,
    distribution = dist_name
  )
  dens <- survival::dsurvreg(
    x = time,
    mean = lp,
    scale = scale_value,
    distribution = dist_name
  )

  list(
    surv = pmin(pmax(as.numeric(surv), 0), 1),
    dens = pmax(as.numeric(dens), 1e-300),
    scale = scale_value
  )
}

fit_parametric_mixture_cure <- function(df, incidence_rhs, latency_rhs, family_name, dataset_name) {
  X_inc <- model_matrix_from_rhs(df, incidence_rhs)
  X_lat <- model_matrix_from_rhs(df, latency_rhs)
  y_time <- pmax(as.numeric(df$time_year_model), 1e-10)
  y_event <- as.integer(df$event_main)

  start_gamma <- tryCatch(
    stats::coef(
      stats::glm(
        stats::as.formula(paste("event_main ~", trim_rhs(incidence_rhs))),
        data = df,
        family = stats::binomial(link = "logit")
      )
    ),
    error = function(e) rep(0, ncol(X_inc))
  )
  start_gamma <- start_gamma[colnames(X_inc)]
  start_gamma[is.na(start_gamma)] <- 0

  start_survreg <- tryCatch(
    survival::survreg(
      formula = make_surv_formula(latency_rhs),
      data = df,
      dist = unname(survreg_family_map[[family_name]])
    ),
    error = function(e) NULL
  )

  if (!is.null(start_survreg)) {
    start_beta <- stats::coef(start_survreg)[colnames(X_lat)]
    start_beta[is.na(start_beta)] <- 0
    start_log_scale <- if (family_has_aux_scale(family_name)) log(start_survreg$scale) else numeric()
  } else {
    start_beta <- rep(0, ncol(X_lat))
    names(start_beta) <- colnames(X_lat)
    start_log_scale <- if (family_has_aux_scale(family_name)) 0 else numeric()
  }

  par_start <- c(start_gamma, start_beta, start_log_scale)

  neg_loglik <- function(par) {
    gamma <- par[seq_len(ncol(X_inc))]
    beta <- par[ncol(X_inc) + seq_len(ncol(X_lat))]
    log_scale_value <- if (family_has_aux_scale(family_name)) par[[length(par)]] else NA_real_

    uncured_prob <- clamp_prob(plogis(drop(X_inc %*% gamma)))
    lp_lat <- drop(X_lat %*% beta)
    dens_surv <- aft_surv_density(
      time = y_time,
      lp = lp_lat,
      family_name = family_name,
      log_scale_value = log_scale_value
    )
    surv_uncured <- pmax(dens_surv$surv, 1e-300)
    dens_uncured <- pmax(dens_surv$dens, 1e-300)

    loglik <- sum(
      y_event * (log(uncured_prob) + log(dens_uncured)) +
        (1 - y_event) * log((1 - uncured_prob) + uncured_prob * surv_uncured)
    )

    if (!is.finite(loglik)) {
      return(1e12)
    }

    -loglik
  }

  attempt_settings <- list(
    list(method = "BFGS", control = list(maxit = 3000, reltol = 1e-10)),
    list(method = "Nelder-Mead", control = list(maxit = 6000, reltol = 1e-10))
  )

  fit_attempts <- lapply(attempt_settings, function(setting) {
    out <- capture_with_warnings(
      stats::optim(
        par = par_start,
        fn = neg_loglik,
        method = setting$method,
        control = setting$control
      )
    )
    out$method <- setting$method
    out
  })

  valid_attempts <- Filter(
    function(x) !is.null(x$value) && is.list(x$value) && is.finite(x$value$value),
    fit_attempts
  )

  if (length(valid_attempts) == 0L) {
    error_messages <- unique(na.omit(vapply(fit_attempts, function(x) x$error_message, character(1))))
    return(list(
      fit_ok = FALSE,
      family = family_name,
      dataset_name = dataset_name,
      error_message = paste(error_messages, collapse = " | ")
    ))
  }

  converged_attempts <- Filter(
    function(x) identical(as.integer(x$value$convergence), 0L),
    valid_attempts
  )
  candidate_attempts <- if (length(converged_attempts) > 0L) converged_attempts else valid_attempts
  best_attempt <- candidate_attempts[[which.min(vapply(candidate_attempts, function(x) x$value$value, numeric(1)))]]
  best_fit <- best_attempt$value

  gamma_hat <- best_fit$par[seq_len(ncol(X_inc))]
  beta_hat <- best_fit$par[ncol(X_inc) + seq_len(ncol(X_lat))]
  log_scale_hat <- if (family_has_aux_scale(family_name)) best_fit$par[[length(best_fit$par)]] else NA_real_

  list(
    fit_ok = identical(as.integer(best_fit$convergence), 0L),
    family = family_name,
    dataset_name = dataset_name,
    gamma = gamma_hat,
    beta = beta_hat,
    log_scale = log_scale_hat,
    incidence_rhs = incidence_rhs,
    latency_rhs = latency_rhs,
    incidence_coef_names = colnames(X_inc),
    latency_coef_names = colnames(X_lat),
    site_levels = levels(df$site),
    optimizer_method = best_attempt$method,
    loglik = -best_fit$value,
    convergence_code = as.integer(best_fit$convergence),
    warnings = unique(best_attempt$warnings),
    aic = 2 * length(best_fit$par) - 2 * (-best_fit$value)
  )
}

predict_parametric_mixture_cure <- function(fit, newdata, times) {
  pred_data <- tibble::as_tibble(newdata) %>%
    mutate(site = factor(as.character(site), levels = fit$site_levels))

  X_inc <- model_matrix_from_rhs(pred_data, fit$incidence_rhs)
  X_lat <- model_matrix_from_rhs(pred_data, fit$latency_rhs)

  uncured_prob <- clamp_prob(plogis(drop(X_inc %*% fit$gamma)))
  lp_lat <- drop(X_lat %*% fit$beta)

  surv_mat <- sapply(times, function(tt) {
    if (tt <= 0) {
      rep(1, nrow(pred_data))
    } else {
      aft_surv_density(
        time = rep(tt, nrow(pred_data)),
        lp = lp_lat,
        family_name = fit$family,
        log_scale_value = fit$log_scale
      )$surv
    }
  })

  if (is.null(dim(surv_mat))) {
    surv_mat <- matrix(surv_mat, nrow = nrow(pred_data), ncol = length(times))
  }

  overall_surv_mat <- (1 - uncured_prob) + uncured_prob * surv_mat

  list(
    susceptible_fraction = uncured_prob,
    susceptible_only_survival = surv_mat,
    susceptible_only_risk = 1 - surv_mat,
    overall_survival = overall_surv_mat,
    overall_risk = 1 - overall_surv_mat
  )
}

subject_specific_median_latency <- function(lp_vec, family_name, log_scale_value) {
  dist_name <- unname(survreg_family_map[[family_name]])
  scale_value <- get_family_scale(family_name, log_scale_value)
  survival::qsurvreg(
    p = 0.5,
    mean = as.numeric(lp_vec),
    scale = scale_value,
    distribution = dist_name
  )
}

solve_cohort_uncured_median <- function(lp_vec, family_name, log_scale_value) {
  dist_name <- unname(survreg_family_map[[family_name]])
  scale_value <- get_family_scale(family_name, log_scale_value)

  lower_bound <- suppressWarnings(min(
    survival::qsurvreg(1e-08, mean = as.numeric(lp_vec), scale = scale_value, distribution = dist_name),
    na.rm = TRUE
  ))
  upper_bound <- suppressWarnings(max(
    survival::qsurvreg(1 - 1e-08, mean = as.numeric(lp_vec), scale = scale_value, distribution = dist_name),
    na.rm = TRUE
  ))

  lower_bound <- pmax(lower_bound, 1e-10)
  if (!is.finite(lower_bound) || lower_bound <= 0) {
    lower_bound <- 1e-10
  }
  if (!is.finite(upper_bound) || upper_bound <= lower_bound) {
    upper_bound <- max(exp(mean(lp_vec, na.rm = TRUE)), 1)
  }

  objective <- function(time_value) {
    mean(
      survival::psurvreg(
        q = time_value,
        mean = as.numeric(lp_vec),
        scale = scale_value,
        distribution = dist_name
      ),
      na.rm = TRUE
    ) - 0.5
  }

  f_lower <- objective(lower_bound)
  f_upper <- objective(upper_bound)

  guard <- 0L
  while (is.finite(f_upper) && f_upper < 0 && guard < 12L) {
    upper_bound <- upper_bound * 2
    f_upper <- objective(upper_bound)
    guard <- guard + 1L
  }

  if (!is.finite(f_lower) || !is.finite(f_upper) || f_lower > 0 || f_upper < 0) {
    return(NA_real_)
  }

  stats::uniroot(objective, interval = c(lower_bound, upper_bound), tol = 1e-08)$root
}

summarize_fit_metrics <- function(spec_row, analysis_df, fit, horizons = c(5, 10)) {
  pred <- predict_parametric_mixture_cure(fit = fit, newdata = analysis_df, times = horizons)

  pred_data <- tibble::as_tibble(analysis_df) %>%
    mutate(site = factor(as.character(site), levels = fit$site_levels))
  X_lat <- model_matrix_from_rhs(pred_data, fit$latency_rhs)
  lp_lat <- drop(X_lat %*% fit$beta)
  subject_median <- subject_specific_median_latency(lp_lat, fit$family, fit$log_scale)
  cohort_uncured_median <- solve_cohort_uncured_median(lp_lat, fit$family, fit$log_scale)
  uncured_prob <- pred$susceptible_fraction

  tibble(
    dataset = spec_row$dataset[[1]],
    dataset_label = spec_row$dataset_label[[1]],
    source_dataset = spec_row$source_dataset[[1]],
    family = fit$family,
    formula_name = spec_row$formula_name[[1]],
    site_adjustment_flag = as.logical(spec_row$site_adjustment_flag[[1]]),
    n_subjects = nrow(analysis_df),
    n_transition = sum(analysis_df$event_main == 1L, na.rm = TRUE),
    n_censored = sum(analysis_df$censor_main == 1L, na.rm = TRUE),
    overall_risk_5y = mean(pred$overall_risk[, horizons == 5], na.rm = TRUE),
    overall_risk_10y = mean(pred$overall_risk[, horizons == 10], na.rm = TRUE),
    cure_fraction = mean(1 - uncured_prob, na.rm = TRUE),
    susceptible_fraction = mean(uncured_prob, na.rm = TRUE),
    susceptible_weighted_subject_latency_median_year = weighted_mean_safe(subject_median, uncured_prob),
    cohort_uncured_median_latency_year = cohort_uncured_median,
    latency_scale = get_family_scale(fit$family, fit$log_scale),
    loglik = as.numeric(fit$loglik %||% NA_real_),
    aic = as.numeric(fit$aic %||% NA_real_),
    optimizer_method = as.character(fit$optimizer_method %||% NA_character_),
    convergence_code = as.integer(fit$convergence_code %||% NA_integer_),
    warnings_text = if (length(fit$warnings %||% character()) > 0L) paste(unique(fit$warnings), collapse = " | ") else NA_character_,
    fit_ok = isTRUE(fit$fit_ok)
  )
}

# Influence analysis ---------------------------------------------------------
select_tail_cases <- function(df, max_n) {
  bind_rows(
    df %>%
      filter(event_main == 0L) %>%
      arrange(desc(time_year), desc(remission_flag), desc(right_censor_flag), unique_person_id) %>%
      mutate(selection_group = "late_censored", selection_rank = row_number()) %>%
      slice_head(n = max_n),
    df %>%
      filter(event_main == 1L) %>%
      arrange(desc(time_year), unique_person_id) %>%
      mutate(selection_group = "late_event", selection_rank = row_number()) %>%
      slice_head(n = max_n)
  ) %>%
    select(unique_person_id, time_year, status_num, censor_detail, selection_group, selection_rank)
}

build_influence_scenarios <- function(df, max_n) {
  selected_cases <- select_tail_cases(df, max_n)
  scenario_rows <- list(
    tibble(
      scenario_type = "baseline",
      removal_scope = "none",
      removal_rank = 0L,
      removed_ids = NA_character_
    )
  )

  if (nrow(selected_cases) == 0L) {
    return(bind_rows(scenario_rows))
  }

  loo_rows <- selected_cases %>%
    mutate(
      scenario_type = "leave_one_out",
      removal_scope = selection_group,
      removal_rank = selection_rank,
      removed_ids = unique_person_id
    ) %>%
    select(scenario_type, removal_scope, removal_rank, removed_ids)

  last_k_rows <- bind_rows(lapply(unique(selected_cases$selection_group), function(group_name) {
    group_tbl <- selected_cases %>%
      filter(selection_group == group_name) %>%
      arrange(selection_rank)
    bind_rows(lapply(seq_len(nrow(group_tbl)), function(k) {
      tibble(
        scenario_type = "leave_last_k_out",
        removal_scope = group_name,
        removal_rank = k,
        removed_ids = format_removed_ids(group_tbl$unique_person_id[seq_len(k)])
      )
    }))
  }))

  bind_rows(scenario_rows, loo_rows, last_k_rows)
}

run_influence_analysis <- function(analysis_datasets) {
  message("Running leave-one-out / leave-last-k-out diagnostics...")

  results <- vector("list", nrow(frequentist_spec_registry))
  idx <- 0L

  for (ii in seq_len(nrow(frequentist_spec_registry))) {
    spec_row <- frequentist_spec_registry[ii, , drop = FALSE]
    source_name <- spec_row$source_dataset[[1]]
    dataset_name <- spec_row$dataset[[1]]
    analysis_df <- prepare_dataset(analysis_datasets[[source_name]], dataset_name)
    scenarios <- build_influence_scenarios(analysis_df, tail_n)

    spec_results <- vector("list", nrow(scenarios))

    for (jj in seq_len(nrow(scenarios))) {
      scenario_row <- scenarios[jj, , drop = FALSE]
      removed_ids_chr <- if (is.na(scenario_row$removed_ids[[1]])) character() else strsplit(scenario_row$removed_ids[[1]], " \\| ", fixed = FALSE)[[1]]
      perturbed_df <- analysis_df %>%
        filter(!(unique_person_id %in% removed_ids_chr))

      fit_result <- fit_parametric_mixture_cure(
        df = perturbed_df,
        incidence_rhs = spec_row$incidence_rhs[[1]],
        latency_rhs = spec_row$latency_rhs[[1]],
        family_name = "lognormal",
        dataset_name = dataset_name
      )

      if (isTRUE(fit_result$fit_ok)) {
        metric_row <- summarize_fit_metrics(spec_row, perturbed_df, fit_result)
      } else {
        metric_row <- tibble(
          dataset = dataset_name,
          dataset_label = spec_row$dataset_label[[1]],
          source_dataset = source_name,
          family = "lognormal",
          formula_name = spec_row$formula_name[[1]],
          site_adjustment_flag = as.logical(spec_row$site_adjustment_flag[[1]]),
          n_subjects = nrow(perturbed_df),
          n_transition = sum(perturbed_df$event_main == 1L, na.rm = TRUE),
          n_censored = sum(perturbed_df$censor_main == 1L, na.rm = TRUE),
          overall_risk_5y = NA_real_,
          overall_risk_10y = NA_real_,
          cure_fraction = NA_real_,
          susceptible_fraction = NA_real_,
          susceptible_weighted_subject_latency_median_year = NA_real_,
          cohort_uncured_median_latency_year = NA_real_,
          latency_scale = NA_real_,
          loglik = NA_real_,
          aic = NA_real_,
          optimizer_method = NA_character_,
          convergence_code = NA_integer_,
          warnings_text = fit_result$error_message %||% NA_character_,
          fit_ok = FALSE
        )
      }

      removed_case_tbl <- analysis_df %>%
        filter(unique_person_id %in% removed_ids_chr)

      spec_results[[jj]] <- metric_row %>%
        mutate(
          scenario_type = scenario_row$scenario_type[[1]],
          removal_scope = scenario_row$removal_scope[[1]],
          removal_rank = as.integer(scenario_row$removal_rank[[1]]),
          removed_ids = scenario_row$removed_ids[[1]],
          removed_n = length(removed_ids_chr),
          removed_max_followup_year = if (nrow(removed_case_tbl) > 0L) max(removed_case_tbl$time_year, na.rm = TRUE) else NA_real_,
          removed_status_mix = if (nrow(removed_case_tbl) > 0L) paste(sort(unique(removed_case_tbl$censor_detail)), collapse = " | ") else NA_character_
        )
    }

    idx <- idx + 1L
    results[[idx]] <- bind_rows(spec_results)
  }

  influence_tbl <- bind_rows(results)
  baseline_tbl <- influence_tbl %>%
    filter(scenario_type == "baseline") %>%
    transmute(
      dataset,
      baseline_overall_risk_5y = overall_risk_5y,
      baseline_overall_risk_10y = overall_risk_10y,
      baseline_cure_fraction = cure_fraction,
      baseline_weighted_subject_latency_median = susceptible_weighted_subject_latency_median_year,
      baseline_cohort_uncured_median = cohort_uncured_median_latency_year,
      baseline_latency_scale = latency_scale
    )

  influence_tbl %>%
    left_join(baseline_tbl, by = "dataset") %>%
    mutate(
      delta_overall_risk_5y = overall_risk_5y - baseline_overall_risk_5y,
      delta_overall_risk_10y = overall_risk_10y - baseline_overall_risk_10y,
      delta_cure_fraction = cure_fraction - baseline_cure_fraction,
      delta_weighted_subject_latency_median = susceptible_weighted_subject_latency_median_year - baseline_weighted_subject_latency_median,
      delta_cohort_uncured_median = cohort_uncured_median_latency_year - baseline_cohort_uncured_median,
      delta_latency_scale = latency_scale - baseline_latency_scale
    ) %>%
    arrange(match(dataset, frequentist_spec_registry$dataset), scenario_type, removal_scope, removal_rank)
}

# Family sensitivity ---------------------------------------------------------
run_family_sensitivity <- function(analysis_datasets) {
  message("Running latency family sensitivity diagnostics...")

  family_order <- c("lognormal", "weibull", "loglogistic", "exponential")

  bind_rows(lapply(seq_len(nrow(frequentist_spec_registry)), function(ii) {
    spec_row <- frequentist_spec_registry[ii, , drop = FALSE]
    source_name <- spec_row$source_dataset[[1]]
    dataset_name <- spec_row$dataset[[1]]
    analysis_df <- prepare_dataset(analysis_datasets[[source_name]], dataset_name)

    bind_rows(lapply(family_order, function(family_name) {
      fit_result <- fit_parametric_mixture_cure(
        df = analysis_df,
        incidence_rhs = spec_row$incidence_rhs[[1]],
        latency_rhs = spec_row$latency_rhs[[1]],
        family_name = family_name,
        dataset_name = dataset_name
      )

      if (isTRUE(fit_result$fit_ok)) {
        summarize_fit_metrics(spec_row, analysis_df, fit_result)
      } else {
        tibble(
          dataset = dataset_name,
          dataset_label = spec_row$dataset_label[[1]],
          source_dataset = source_name,
          family = family_name,
          formula_name = spec_row$formula_name[[1]],
          site_adjustment_flag = as.logical(spec_row$site_adjustment_flag[[1]]),
          n_subjects = nrow(analysis_df),
          n_transition = sum(analysis_df$event_main == 1L, na.rm = TRUE),
          n_censored = sum(analysis_df$censor_main == 1L, na.rm = TRUE),
          overall_risk_5y = NA_real_,
          overall_risk_10y = NA_real_,
          cure_fraction = NA_real_,
          susceptible_fraction = NA_real_,
          susceptible_weighted_subject_latency_median_year = NA_real_,
          cohort_uncured_median_latency_year = NA_real_,
          latency_scale = NA_real_,
          loglik = NA_real_,
          aic = NA_real_,
          optimizer_method = NA_character_,
          convergence_code = NA_integer_,
          warnings_text = fit_result$error_message %||% NA_character_,
          fit_ok = FALSE
        )
      }
    }))
  })) %>%
    group_by(dataset) %>%
    mutate(
      best_aic = if (any(is.finite(aic))) min(aic[is.finite(aic)], na.rm = TRUE) else NA_real_,
      delta_aic_vs_best = aic - best_aic,
      best_family_by_aic = if (any(is.finite(aic))) family[which.min(if_else(is.finite(aic), aic, Inf))] else NA_character_
    ) %>%
    ungroup() %>%
    select(-best_aic) %>%
    arrange(match(dataset, frequentist_spec_registry$dataset), match(family, family_order))
}

# Frequentist baseline validation -------------------------------------------
run_baseline_validation <- function(analysis_datasets, frequentist_fit_cache) {
  message("Validating reconstructed lognormal baseline against saved frequentist fits...")

  bind_rows(lapply(seq_len(nrow(frequentist_spec_registry)), function(ii) {
    spec_row <- frequentist_spec_registry[ii, , drop = FALSE]
    dataset_name <- spec_row$dataset[[1]]
    source_name <- spec_row$source_dataset[[1]]
    analysis_df <- prepare_dataset(analysis_datasets[[source_name]], dataset_name)
    cached_fit <- frequentist_fit_cache[[dataset_name]]

    if (is.null(cached_fit)) {
      return(tibble(dataset = dataset_name, validation_status = "missing_cached_fit"))
    }

    spec_metrics <- summarize_fit_metrics(spec_row, analysis_df, c(
      list(
        fit_ok = TRUE,
        family = "lognormal",
        gamma = cached_fit$gamma,
        beta = cached_fit$beta,
        log_scale = cached_fit$log_sigma,
        incidence_rhs = cached_fit$incidence_rhs,
        latency_rhs = cached_fit$latency_rhs,
        site_levels = cached_fit$site_levels,
        loglik = cached_fit$loglik,
        aic = 2 * (length(cached_fit$gamma) + length(cached_fit$beta) + 1L) - 2 * cached_fit$loglik,
        optimizer_method = cached_fit$optimizer_method,
        convergence_code = cached_fit$convergence_code,
        warnings = cached_fit$warnings
      )
    ))

    refit_result <- fit_parametric_mixture_cure(
      df = analysis_df,
      incidence_rhs = spec_row$incidence_rhs[[1]],
      latency_rhs = spec_row$latency_rhs[[1]],
      family_name = "lognormal",
      dataset_name = dataset_name
    )

    if (!isTRUE(refit_result$fit_ok)) {
      return(spec_metrics %>% mutate(validation_status = "refit_failed"))
    }

    refit_metrics <- summarize_fit_metrics(spec_row, analysis_df, refit_result)

    spec_metrics %>%
      transmute(
        dataset,
        dataset_label,
        cached_overall_risk_10y = overall_risk_10y,
        cached_cure_fraction = cure_fraction,
        cached_weighted_subject_latency_median = susceptible_weighted_subject_latency_median_year,
        cached_cohort_uncured_median = cohort_uncured_median_latency_year
      ) %>%
      bind_cols(
        refit_metrics %>%
          transmute(
            refit_overall_risk_10y = overall_risk_10y,
            refit_cure_fraction = cure_fraction,
            refit_weighted_subject_latency_median = susceptible_weighted_subject_latency_median_year,
            refit_cohort_uncured_median = cohort_uncured_median_latency_year
          )
      ) %>%
      mutate(
        delta_overall_risk_10y = refit_overall_risk_10y - cached_overall_risk_10y,
        delta_cure_fraction = refit_cure_fraction - cached_cure_fraction,
        delta_weighted_subject_latency_median = refit_weighted_subject_latency_median - cached_weighted_subject_latency_median,
        delta_cohort_uncured_median = refit_cohort_uncured_median - cached_cohort_uncured_median,
        validation_status = "ok"
      )
  }))
}

# Bayesian joint uncertainty -------------------------------------------------
build_bayesian_design_matrices <- function(df, model_row) {
  z_i <- as.integer(df$sex_num)
  x20_i <- as.integer(df$age_exact_entry >= 20 & df$age_exact_entry < 30)
  x30_i <- as.integer(df$age_exact_entry >= 30)
  s_i <- as.integer(as.character(df$site) == "SNU")

  X_inc <- cbind(
    sex_num = z_i,
    age20_29 = x20_i,
    age30plus = x30_i,
    sex_x_age20_29 = z_i * x20_i,
    sex_x_age30plus = z_i * x30_i
  )

  if (isTRUE(model_row$incidence_site_indicator[[1]])) {
    X_inc <- cbind(X_inc, site_SNU = s_i)
  }

  X_lat <- switch(
    as.character(model_row$structural_model_id[[1]]),
    "PNU-L0" = cbind(age_s = as.numeric(df$age_s), sex_num = z_i),
    "SNU-L0" = cbind(age_s = as.numeric(df$age_s), sex_num = z_i),
    "MERGED-I0-L0S0" = cbind(age_s = as.numeric(df$age_s), sex_num = z_i),
    "MERGED-I1-L0S1" = cbind(age_s = as.numeric(df$age_s), sex_num = z_i, site_SNU = s_i),
    stop("Unknown structural model id: ", model_row$structural_model_id[[1]], call. = FALSE)
  )

  list(
    X_inc = unclass(as.matrix(X_inc)),
    X_lat = unclass(as.matrix(X_lat))
  )
}

select_joint_draw_indices <- function(n_total, max_draws) {
  if (n_total <= max_draws) {
    return(seq_len(n_total))
  }
  unique(round(seq(1, n_total, length.out = max_draws)))
}

solve_bayesian_cohort_uncured_median <- function(mu_vec, sigma_value) {
  objective <- function(log_time_value) {
    mean(stats::pnorm((log_time_value - mu_vec) / sigma_value), na.rm = TRUE) - 0.5
  }

  lower_bound <- min(mu_vec - 12 * sigma_value, na.rm = TRUE)
  upper_bound <- max(mu_vec + 12 * sigma_value, na.rm = TRUE)
  f_lower <- objective(lower_bound)
  f_upper <- objective(upper_bound)

  if (!is.finite(f_lower) || !is.finite(f_upper)) {
    return(NA_real_)
  }
  if (f_lower > 0) {
    lower_bound <- lower_bound - 10 * sigma_value - 5
    f_lower <- objective(lower_bound)
  }
  if (f_upper < 0) {
    upper_bound <- upper_bound + 10 * sigma_value + 5
    f_upper <- objective(upper_bound)
  }
  if (!is.finite(f_lower) || !is.finite(f_upper) || f_lower > 0 || f_upper < 0) {
    return(NA_real_)
  }

  exp(stats::uniroot(objective, interval = c(lower_bound, upper_bound), tol = 1e-08)$root)
}

extract_bayesian_joint_draws <- function(model_row, analysis_df) {
  fit_path <- as.character(model_row$fit_rds_path[[1]])
  assert_file_exists(fit_path, "Bayesian fit RDS")
  fit <- readRDS(fit_path)

  design_bundle <- build_bayesian_design_matrices(analysis_df, model_row)
  ext <- rstan::extract(fit, pars = c("alpha_inc", "beta_inc", "gamma0", "gamma_lat", "log_sigma"), permuted = TRUE)

  beta_inc <- ensure_draw_matrix(ext$beta_inc, ncol(design_bundle$X_inc))
  gamma_lat <- ensure_draw_matrix(ext$gamma_lat, ncol(design_bundle$X_lat))
  draw_keep_idx <- select_joint_draw_indices(nrow(beta_inc), joint_draws_max)

  beta_inc <- beta_inc[draw_keep_idx, , drop = FALSE]
  gamma_lat <- gamma_lat[draw_keep_idx, , drop = FALSE]
  alpha_inc <- as.numeric(ext$alpha_inc)[draw_keep_idx]
  gamma0 <- as.numeric(ext$gamma0)[draw_keep_idx]
  log_sigma <- as.numeric(ext$log_sigma)[draw_keep_idx]

  eta_inc_mat <- design_bundle$X_inc %*% t(beta_inc)
  eta_inc_mat <- sweep(eta_inc_mat, 2, alpha_inc, "+")
  pi_mat <- plogis(eta_inc_mat)
  pi_mat <- pmin(pmax(pi_mat, 1e-12), 1 - 1e-12)

  mu_lat_mat <- design_bundle$X_lat %*% t(gamma_lat)
  mu_lat_mat <- sweep(mu_lat_mat, 2, gamma0, "+")
  sigma_vec <- exp(log_sigma)

  weighted_subject_latency_median <- vapply(seq_along(sigma_vec), function(jj) {
    weighted_mean_safe(exp(mu_lat_mat[, jj]), pi_mat[, jj])
  }, numeric(1))

  cohort_uncured_median <- vapply(seq_along(sigma_vec), function(jj) {
    solve_bayesian_cohort_uncured_median(mu_lat_mat[, jj], sigma_vec[[jj]])
  }, numeric(1))

  tibble(
    dataset = as.character(model_row$dataset[[1]]),
    dataset_label = as.character(model_row$dataset_label[[1]]),
    model_id = as.character(model_row$model_id[[1]]),
    prior_branch = as.character(model_row$prior_branch[[1]]),
    site_adjustment_flag = as.logical(model_row$site_adjustment_flag[[1]]),
    draw_id = seq_along(sigma_vec),
    cure_fraction = 1 - colMeans(pi_mat),
    susceptible_fraction = colMeans(pi_mat),
    latency_sigma = sigma_vec,
    weighted_subject_latency_median_year = weighted_subject_latency_median,
    cohort_uncured_median_latency_year = cohort_uncured_median,
    n_draws_total = nrow(ensure_draw_matrix(ext$beta_inc, ncol(design_bundle$X_inc))),
    n_draws_used = length(draw_keep_idx)
  )
}

run_bayesian_joint_uncertainty <- function(analysis_datasets) {
  message("Running Bayesian joint uncertainty diagnostics...")
  model_registry <- readr::read_csv(bayesian_model_registry_file, show_col_types = FALSE, progress = FALSE)

  draw_tbl <- bind_rows(lapply(seq_len(nrow(model_registry)), function(ii) {
    model_row <- model_registry[ii, , drop = FALSE]
    source_name <- as.character(model_row$source_dataset_key[[1]])
    analysis_df <- prepare_dataset(analysis_datasets[[source_name]], as.character(model_row$dataset[[1]]))
    extract_bayesian_joint_draws(model_row, analysis_df)
  }))

  summary_tbl <- draw_tbl %>%
    group_by(dataset, dataset_label, model_id, prior_branch, site_adjustment_flag) %>%
    summarise(
      n_draws_total = first(n_draws_total),
      n_draws_used = first(n_draws_used),
      cor_cure_vs_sigma = stats::cor(cure_fraction, latency_sigma, use = "pairwise.complete.obs"),
      cor_cure_vs_weighted_subject_median = stats::cor(cure_fraction, weighted_subject_latency_median_year, use = "pairwise.complete.obs"),
      cor_cure_vs_cohort_uncured_median = stats::cor(cure_fraction, cohort_uncured_median_latency_year, use = "pairwise.complete.obs"),
      cor_sigma_vs_cohort_uncured_median = stats::cor(latency_sigma, cohort_uncured_median_latency_year, use = "pairwise.complete.obs"),
      cure_fraction_mean = mean(cure_fraction, na.rm = TRUE),
      cure_fraction_q025 = safe_quantile(cure_fraction, 0.025),
      cure_fraction_q975 = safe_quantile(cure_fraction, 0.975),
      latency_sigma_mean = mean(latency_sigma, na.rm = TRUE),
      latency_sigma_q025 = safe_quantile(latency_sigma, 0.025),
      latency_sigma_q975 = safe_quantile(latency_sigma, 0.975),
      weighted_subject_median_mean = mean(weighted_subject_latency_median_year, na.rm = TRUE),
      weighted_subject_median_q025 = safe_quantile(weighted_subject_latency_median_year, 0.025),
      weighted_subject_median_q975 = safe_quantile(weighted_subject_latency_median_year, 0.975),
      cohort_uncured_median_mean = mean(cohort_uncured_median_latency_year, na.rm = TRUE),
      cohort_uncured_median_q025 = safe_quantile(cohort_uncured_median_latency_year, 0.025),
      cohort_uncured_median_q975 = safe_quantile(cohort_uncured_median_latency_year, 0.975),
      .groups = "drop"
    ) %>%
    arrange(match(dataset, c("PNU", "SNU", "merged_no_site", "merged_site_adjusted")), prior_branch)

  list(draws = draw_tbl, summary = summary_tbl)
}

make_joint_scatter_plot <- function(draw_tbl, y_var, y_label, log10_y = FALSE) {
  plot_obj <- ggplot(draw_tbl, aes(x = cure_fraction, y = .data[[y_var]], color = prior_branch)) +
    geom_point(alpha = 0.22, size = 0.7, stroke = 0) +
    facet_wrap(~ dataset_label, scales = "free_y") +
    scale_color_manual(
      values = c(
        anchor_informed = "#1f77b4",
        neutral_no_external_info = "#d62728"
      ),
      breaks = c("anchor_informed", "neutral_no_external_info"),
      labels = c("Anchor-informed", "Neutral")
    ) +
    labs(
      title = sprintf("Bayesian joint uncertainty: cure fraction vs %s", y_label),
      x = "Posterior cure fraction",
      y = y_label,
      color = "Prior branch"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold")
    )

  if (isTRUE(log10_y)) {
    plot_obj <- plot_obj +
      scale_y_log10(labels = label_number(accuracy = 0.1))
  }

  plot_obj
}

write_summary_markdown <- function(path, influence_tbl, family_tbl, bayes_summary_tbl) {
  top_influence_tbl <- influence_tbl %>%
    filter(scenario_type != "baseline") %>%
    mutate(abs_delta_latency = abs(delta_weighted_subject_latency_median)) %>%
    arrange(desc(abs_delta_latency)) %>%
    slice_head(n = 12)

  family_best_tbl <- family_tbl %>%
    filter(fit_ok) %>%
    arrange(dataset, delta_aic_vs_best, aic) %>%
    group_by(dataset, dataset_label) %>%
    slice_head(n = 4) %>%
    ungroup()

  bayes_focus_tbl <- bayes_summary_tbl %>%
    transmute(
      dataset = dataset_label,
      prior_branch,
      cor_cure_vs_sigma = round(cor_cure_vs_sigma, 3),
      cor_cure_vs_cohort_uncured_median = round(cor_cure_vs_cohort_uncured_median, 3),
      cohort_uncured_median_mean = round(cohort_uncured_median_mean, 2),
      cohort_uncured_median_q975 = round(cohort_uncured_median_q975, 2)
    )

  lines <- c(
    "# Latency instability priority diagnostics",
    "",
    "## Leave-one-out / leave-last-k-out",
    ""
  )

  if (nrow(top_influence_tbl) == 0L) {
    lines <- c(lines, "- No influence rows were generated.", "")
  } else {
    lines <- c(
      lines,
      sprintf(
        "- %s / %s / rank %d: delta weighted latency median = %.2f years; delta cure fraction = %.3f; removed = %s",
        top_influence_tbl$dataset_label,
        top_influence_tbl$removal_scope,
        top_influence_tbl$removal_rank,
        top_influence_tbl$delta_weighted_subject_latency_median,
        top_influence_tbl$delta_cure_fraction,
        top_influence_tbl$removed_ids
      ),
      ""
    )
  }

  lines <- c(lines, "## Latency family sensitivity", "")
  if (nrow(family_best_tbl) == 0L) {
    lines <- c(lines, "- No family-sensitivity rows were generated.", "")
  } else {
    lines <- c(
      lines,
      sprintf(
        "- %s / %s: AIC=%.2f, cure fraction=%.3f, weighted latency median=%.2f years",
        family_best_tbl$dataset_label,
        family_best_tbl$family,
        family_best_tbl$aic,
        family_best_tbl$cure_fraction,
        family_best_tbl$susceptible_weighted_subject_latency_median_year
      ),
      ""
    )
  }

  lines <- c(lines, "## Bayesian joint uncertainty", "")
  if (nrow(bayes_focus_tbl) == 0L) {
    lines <- c(lines, "- No Bayesian joint-uncertainty rows were generated.", "")
  } else {
    lines <- c(
      lines,
      sprintf(
        "- %s / %s: cor(cure, sigma)=%.3f; cor(cure, cohort uncured median)=%.3f; posterior mean uncured median=%.2f y; q975=%.2f y",
        bayes_focus_tbl$dataset,
        bayes_focus_tbl$prior_branch,
        bayes_focus_tbl$cor_cure_vs_sigma,
        bayes_focus_tbl$cor_cure_vs_cohort_uncured_median,
        bayes_focus_tbl$cohort_uncured_median_mean,
        bayes_focus_tbl$cohort_uncured_median_q975
      ),
      ""
    )
  }

  writeLines(lines, con = path)
}

# Run -----------------------------------------------------------------------
assert_file_exists(analysis_dataset_file, "Analysis dataset RDS")
assert_file_exists(frequentist_fit_cache_file, "Frequentist fit cache RDS")
assert_file_exists(bayesian_model_registry_file, "Bayesian model registry CSV")

analysis_datasets <- readRDS(analysis_dataset_file)
frequentist_fit_cache <- readRDS(frequentist_fit_cache_file)

influence_file <- file.path(output_dir, "latency_instability_influence_log_normal.csv")
influence_top_file <- file.path(output_dir, "latency_instability_influence_top.csv")
family_file <- file.path(output_dir, "latency_instability_family_sensitivity.csv")
family_plot_file <- file.path(output_dir, "latency_instability_family_sensitivity_plot.png")
baseline_validation_file <- file.path(output_dir, "latency_instability_baseline_validation.csv")
bayes_draws_file <- file.path(output_dir, "latency_instability_bayesian_joint_draws.csv")
bayes_summary_file <- file.path(output_dir, "latency_instability_bayesian_joint_summary.csv")
bayes_sigma_plot_file <- file.path(output_dir, "latency_instability_bayesian_cure_vs_sigma.png")
bayes_median_plot_file <- file.path(output_dir, "latency_instability_bayesian_cure_vs_uncured_median.png")
summary_md_file <- file.path(output_dir, "latency_instability_priority_summary.md")

baseline_validation_tbl <- run_baseline_validation(analysis_datasets, frequentist_fit_cache)
influence_tbl <- run_influence_analysis(analysis_datasets)
family_tbl <- run_family_sensitivity(analysis_datasets)
bayes_joint <- run_bayesian_joint_uncertainty(analysis_datasets)

influence_top_tbl <- influence_tbl %>%
  filter(scenario_type != "baseline") %>%
  mutate(
    abs_delta_cure_fraction = abs(delta_cure_fraction),
    abs_delta_weighted_subject_latency_median = abs(delta_weighted_subject_latency_median),
    abs_delta_cohort_uncured_median = abs(delta_cohort_uncured_median)
  ) %>%
  arrange(
    desc(abs_delta_weighted_subject_latency_median),
    desc(abs_delta_cure_fraction)
  ) %>%
  group_by(dataset, removal_scope) %>%
  slice_head(n = 5) %>%
  ungroup()

family_plot_tbl <- family_tbl %>%
  filter(fit_ok)

family_plot <- ggplot(
  family_plot_tbl,
  aes(x = family, y = susceptible_weighted_subject_latency_median_year, color = family)
) +
  geom_point(size = 2.6) +
  facet_wrap(~ dataset_label, scales = "free_y") +
  scale_y_log10(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "Latency family sensitivity: weighted uncured-latency median",
    x = "Latency family",
    y = "Weighted subject-specific latency median (years)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

bayes_sigma_plot <- make_joint_scatter_plot(
  draw_tbl = bayes_joint$draws,
  y_var = "latency_sigma",
  y_label = "Latency sigma",
  log10_y = FALSE
)

bayes_median_plot <- make_joint_scatter_plot(
  draw_tbl = bayes_joint$draws,
  y_var = "cohort_uncured_median_latency_year",
  y_label = "Cohort uncured median latency (years)",
  log10_y = TRUE
)

readr::write_csv(baseline_validation_tbl, baseline_validation_file)
readr::write_csv(influence_tbl, influence_file)
readr::write_csv(influence_top_tbl, influence_top_file)
readr::write_csv(family_tbl, family_file)
readr::write_csv(bayes_joint$draws, bayes_draws_file)
readr::write_csv(bayes_joint$summary, bayes_summary_file)

ggplot2::ggsave(family_plot_file, family_plot, width = 10, height = 6, dpi = 320)
ggplot2::ggsave(bayes_sigma_plot_file, bayes_sigma_plot, width = 10, height = 6, dpi = 320)
ggplot2::ggsave(bayes_median_plot_file, bayes_median_plot, width = 10, height = 6, dpi = 320)

write_summary_markdown(summary_md_file, influence_tbl, family_tbl, bayes_joint$summary)

message("Baseline validation CSV: ", normalizePath(baseline_validation_file, winslash = "/", mustWork = FALSE))
message("Influence diagnostics CSV: ", normalizePath(influence_file, winslash = "/", mustWork = FALSE))
message("Influence top CSV: ", normalizePath(influence_top_file, winslash = "/", mustWork = FALSE))
message("Family sensitivity CSV: ", normalizePath(family_file, winslash = "/", mustWork = FALSE))
message("Bayesian joint draws CSV: ", normalizePath(bayes_draws_file, winslash = "/", mustWork = FALSE))
message("Bayesian joint summary CSV: ", normalizePath(bayes_summary_file, winslash = "/", mustWork = FALSE))
message("Summary markdown: ", normalizePath(summary_md_file, winslash = "/", mustWork = FALSE))
