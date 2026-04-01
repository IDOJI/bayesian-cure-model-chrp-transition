# Configure Paths ---------------------------------------------------------
find_repo_root <- function(start_dir) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)

  repeat {
    stage7_reference <- file.path(current_dir, "2.Rcode", "stage7_ Frequentist cure block.r")
    if (file.exists(stage7_reference)) {
      return(current_dir)
    }

    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      break
    }
    current_dir <- parent_dir
  }

  stop(
    "Could not locate the repository root containing `2.Rcode/stage7_ Frequentist cure block.r`.",
    call. = FALSE
  )
}

first_existing_dir <- function(paths) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  existing <- paths[dir.exists(paths)]
  if (length(existing) > 0L) {
    return(existing[[1L]])
  }
  if (length(paths) > 0L) {
    return(paths[[1L]])
  }
  NA_character_
}

command_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", command_args, value = TRUE)
search_start_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

repo_root <- find_repo_root(search_start_dir)
results_root_default <- file.path(repo_root, "3.Results files")

sys_name <- Sys.info()[["sysname"]]
dropbox_results_root_default <- switch(
  sys_name,
  "Darwin" = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  "Windows" = "C:/Users/clair/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  NA_character_
)

stage1_path_default <- first_existing_dir(
  c(
    file.path(results_root_default, "stage1_Backbone lock"),
    file.path(dropbox_results_root_default, "stage1_Backbone lock")
  )
)

stage1_path <- Sys.getenv("STAGE7_SIMPLE_STAGE1_PATH", unset = stage1_path_default)
export_path <- Sys.getenv("STAGE7_SIMPLE_EXPORT_PATH", unset = results_root_default)

stage1_analysis_datasets_file <- file.path(stage1_path, "stage1_analysis_datasets.rds")
stage1_formula_registry_file <- file.path(stage1_path, "stage1_formula_registry.csv")

horizon_years <- 1:10
curve_horizon_max_year <- 10
curve_step_year <- 0.05
time_origin_epsilon_year <- 1e-10
main_risk_scale <- "transition_only_main"
model_id_value <- "frequentist_mixture_cure_lognormal"
plot_width_in <- 10
plot_height_in <- 6
plot_dpi <- 320

horizon_summary_file <- file.path(
  export_path,
  "stage7_simple_lognormal_mixture_cure_horizon_summary.csv"
)
plot_source_file <- file.path(
  export_path,
  "stage7_simple_lognormal_mixture_cure_plot_source.csv"
)
fit_object_rds_file <- file.path(
  export_path,
  "stage7_simple_lognormal_mixture_cure_fitted_objects.rds"
)
plot_rds_file <- file.path(
  export_path,
  "stage7_simple_lognormal_mixture_cure_plot_objects.rds"
)
overall_survival_png_file <- file.path(
  export_path,
  "stage7_simple_lognormal_mixture_cure_overall_survival.png"
)
overall_risk_png_file <- file.path(
  export_path,
  "stage7_simple_lognormal_mixture_cure_overall_risk.png"
)
susceptible_survival_png_file <- file.path(
  export_path,
  "stage7_simple_lognormal_mixture_cure_susceptible_only_survival.png"
)
susceptible_risk_png_file <- file.path(
  export_path,
  "stage7_simple_lognormal_mixture_cure_susceptible_only_risk.png"
)

dataset_model_registry <- tibble::tibble(
  dataset = c("PNU", "SNU", "merged_no_site", "merged_site_adjusted"),
  dataset_label = c("PNU", "SNU", "merged", "merged (site-adjusted)"),
  source_dataset = c("PNU", "SNU", "merged", "merged"),
  formula_name = c("base", "base", "base", "site_added"),
  site_adjustment_flag = c(FALSE, FALSE, FALSE, TRUE)
)

dataset_palette <- c(
  "PNU" = "#1B4332",
  "SNU" = "#2A6F97",
  "merged_no_site" = "#C1666B",
  "merged_site_adjusted" = "#B8860B"
)

# Load Packages -----------------------------------------------------------
required_packages <- c("readr", "dplyr", "tibble", "ggplot2", "scales", "survival")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install required packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(survival)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# Define Helpers ----------------------------------------------------------
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

assert_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

read_csv_checked <- function(path, label) {
  assert_file_exists(path, label)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
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

trim_rhs <- function(rhs) {
  rhs <- gsub("\\s+", " ", trimws(as.character(rhs)))
  rhs[is.na(rhs)] <- NA_character_
  rhs
}

clamp_prob <- function(x, eps = 1e-12) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
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

lognormal_surv_density <- function(time, lp, log_sigma) {
  time <- pmax(as.numeric(time), time_origin_epsilon_year)
  lp <- as.numeric(lp)
  sigma <- exp(as.numeric(log_sigma))
  z <- (log(time) - lp) / sigma
  surv <- 1 - stats::pnorm(z)
  dens <- stats::dnorm(z) / (sigma * time)
  list(
    surv = pmin(pmax(surv, 0), 1),
    dens = pmax(dens, 1e-300)
  )
}

validate_inputs <- function(analysis_datasets, formula_registry) {
  if (!is.list(analysis_datasets) || !all(c("PNU", "SNU", "merged") %in% names(analysis_datasets))) {
    stop("Stage 1 analysis dataset bundle must contain `PNU`, `SNU`, and `merged`.", call. = FALSE)
  }

  required_formula_cols <- c(
    "dataset",
    "formula_id",
    "formula_name",
    "formula_label",
    "formula_rhs",
    "site_branch",
    "interaction_branch"
  )
  missing_formula_cols <- setdiff(required_formula_cols, names(formula_registry))
  if (length(missing_formula_cols) > 0L) {
    stop(
      sprintf(
        "Stage 1 formula registry is missing required columns: %s",
        paste(missing_formula_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

prepare_dataset <- function(df, dataset_name) {
  required_cols <- c("site", "id", "sex_num", "age_exact_entry", "age_s", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("[%s] Dataset is missing required columns: %s", dataset_name, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  out <- tibble::as_tibble(df)

  if (!("unique_person_id" %in% names(out))) {
    out$unique_person_id <- paste(out$site, out$id, sep = "_")
  }

  if (!("time_year" %in% names(out))) {
    if (!("days_followup" %in% names(out))) {
      stop(sprintf("[%s] Dataset must contain `time_year` or `days_followup`.", dataset_name), call. = FALSE)
    }
    out$time_year <- as.numeric(out$days_followup) / 365.25
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
      time_year_model = pmax(as.numeric(time_year), time_origin_epsilon_year),
      event_main = as.integer(status_num == 1L),
      censor_main = as.integer(status_num %in% c(0L, 2L))
    )

  if (nrow(out) != dplyr::n_distinct(out$unique_person_id)) {
    stop(sprintf("[%s] `unique_person_id` must be unique.", dataset_name), call. = FALSE)
  }

  if (anyNA(out[, c("unique_person_id", "site", "id", "sex_num", "age_exact_entry", "age_s", "status_num", "time_year")])) {
    stop(sprintf("[%s] Missing values detected in required analysis columns.", dataset_name), call. = FALSE)
  }

  if (any(out$time_year < 0, na.rm = TRUE)) {
    stop(sprintf("[%s] Negative follow-up times are not allowed.", dataset_name), call. = FALSE)
  }

  if (any(!out$status_num %in% c(0L, 1L, 2L))) {
    stop(sprintf("[%s] `status_num` must be coded as 0/1/2 only.", dataset_name), call. = FALSE)
  }

  out
}

select_model_specs <- function(formula_registry) {
  formula_tbl <- formula_registry %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      formula_name = trimws(as.character(formula_name)),
      formula_rhs = trim_rhs(formula_rhs)
    )

  selected_specs <- dataset_model_registry %>%
    left_join(
      formula_tbl,
      by = c(
        "source_dataset" = "dataset",
        "formula_name" = "formula_name"
      )
    )

  if (anyNA(selected_specs$formula_id) || anyNA(selected_specs$formula_rhs)) {
    missing_keys <- selected_specs$dataset[is.na(selected_specs$formula_id) | is.na(selected_specs$formula_rhs)]
    stop(
      sprintf(
        "Could not resolve formula definitions for dataset versions: %s",
        paste(unique(missing_keys), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Follow the retained Stage 7 simplified rule for the merged site-adjusted
  # mixture cure branch: keep `site` in both the incidence part and the
  # latency part. This matches the earlier `site_in_both` / `base_siteadjusted_both`
  # selection, rather than placing `site` in only one submodel.
  selected_specs %>%
    mutate(
      incidence_rhs = formula_rhs,
      latency_rhs = formula_rhs,
      site_in_incidence = grepl("\\bsite\\b", incidence_rhs),
      site_in_latency = grepl("\\bsite\\b", latency_rhs)
    ) %>%
    select(
      dataset,
      dataset_label,
      source_dataset,
      formula_name,
      formula_id,
      formula_label,
      formula_rhs,
      site_branch,
      interaction_branch,
      site_adjustment_flag,
      incidence_rhs,
      latency_rhs,
      site_in_incidence,
      site_in_latency
    )
}

fit_lognormal_mixture_cure <- function(df, incidence_rhs, latency_rhs, dataset_name) {
  X_inc <- model_matrix_from_rhs(df, incidence_rhs)
  X_lat <- model_matrix_from_rhs(df, latency_rhs)
  y_time <- pmax(as.numeric(df$time_year_model), time_origin_epsilon_year)
  y_event <- as.integer(df$event_main)

  start_gamma <- tryCatch(
    {
      stats::coef(
        stats::glm(
          stats::as.formula(paste("event_main ~", trim_rhs(incidence_rhs))),
          data = df,
          family = stats::binomial(link = "logit")
        )
      )
    },
    error = function(e) rep(0, ncol(X_inc))
  )
  start_gamma <- start_gamma[colnames(X_inc)]
  start_gamma[is.na(start_gamma)] <- 0

  start_survreg <- tryCatch(
    survival::survreg(
      formula = make_surv_formula(latency_rhs),
      data = df,
      dist = "lognormal"
    ),
    error = function(e) NULL
  )

  if (!is.null(start_survreg)) {
    start_beta <- stats::coef(start_survreg)[colnames(X_lat)]
    start_beta[is.na(start_beta)] <- 0
    start_log_sigma <- log(start_survreg$scale)
  } else {
    start_beta <- rep(0, ncol(X_lat))
    names(start_beta) <- colnames(X_lat)
    start_log_sigma <- 0
  }

  par_start <- c(start_gamma, start_beta, start_log_sigma)

  neg_loglik <- function(par) {
    gamma <- par[seq_len(ncol(X_inc))]
    beta <- par[ncol(X_inc) + seq_len(ncol(X_lat))]
    log_sigma <- par[[length(par)]]

    uncured_prob <- clamp_prob(plogis(drop(X_inc %*% gamma)))
    lp_lat <- drop(X_lat %*% beta)
    dens_surv <- lognormal_surv_density(y_time, lp_lat, log_sigma)
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
    list(method = "BFGS", control = list(maxit = 2000, reltol = 1e-10)),
    list(method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-10))
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
    stop(
      sprintf(
        "[%s] Log-normal mixture cure fit failed: %s",
        dataset_name,
        paste(error_messages, collapse = " | ")
      ),
      call. = FALSE
    )
  }

  converged_attempts <- Filter(
    function(x) identical(as.integer(x$value$convergence), 0L),
    valid_attempts
  )
  candidate_attempts <- if (length(converged_attempts) > 0L) converged_attempts else valid_attempts
  best_index <- which.min(vapply(candidate_attempts, function(x) x$value$value, numeric(1)))
  best_attempt <- candidate_attempts[[best_index]]
  best_fit <- best_attempt$value

  if (!identical(as.integer(best_fit$convergence), 0L)) {
    stop(
      sprintf(
        "[%s] Log-normal mixture cure optimizer did not converge. Best attempt used `%s` with convergence code %s.",
        dataset_name,
        best_attempt$method,
        as.character(best_fit$convergence)
      ),
      call. = FALSE
    )
  }

  gamma_hat <- best_fit$par[seq_len(ncol(X_inc))]
  beta_hat <- best_fit$par[ncol(X_inc) + seq_len(ncol(X_lat))]
  log_sigma_hat <- best_fit$par[[length(best_fit$par)]]

  list(
    gamma = gamma_hat,
    beta = beta_hat,
    log_sigma = log_sigma_hat,
    incidence_rhs = incidence_rhs,
    latency_rhs = latency_rhs,
    incidence_coef_names = colnames(X_inc),
    latency_coef_names = colnames(X_lat),
    site_levels = levels(df$site),
    optimizer_method = best_attempt$method,
    loglik = -best_fit$value,
    convergence_code = as.integer(best_fit$convergence),
    warnings = unique(best_attempt$warnings)
  )
}

load_fitted_object_cache <- function(path) {
  if (!file.exists(path)) {
    return(list())
  }

  cache_obj <- readRDS(path)
  if (!is.list(cache_obj)) {
    return(list())
  }

  cache_obj
}

has_usable_cached_fit <- function(fit_obj, spec_row) {
  required_fields <- c(
    "gamma",
    "beta",
    "log_sigma",
    "incidence_rhs",
    "latency_rhs",
    "site_levels"
  )

  if (!is.list(fit_obj) || !all(required_fields %in% names(fit_obj))) {
    return(FALSE)
  }

  identical(as.character(fit_obj$incidence_rhs), as.character(spec_row$incidence_rhs[[1]])) &&
    identical(as.character(fit_obj$latency_rhs), as.character(spec_row$latency_rhs[[1]]))
}

predict_lognormal_mixture_cure <- function(fit, newdata, times) {
  pred_data <- tibble::as_tibble(newdata) %>%
    mutate(site = factor(as.character(site), levels = fit$site_levels))

  X_inc <- model_matrix_from_rhs(pred_data, fit$incidence_rhs)
  X_lat <- model_matrix_from_rhs(pred_data, fit$latency_rhs)

  uncured_prob <- clamp_prob(plogis(drop(X_inc %*% fit$gamma)))
  lp_lat <- drop(X_lat %*% fit$beta)

  susceptible_survival_mat <- sapply(times, function(tt) {
    if (tt <= 0) {
      rep(1, nrow(pred_data))
    } else {
      lognormal_surv_density(rep(tt, nrow(pred_data)), lp_lat, fit$log_sigma)$surv
    }
  })

  if (is.null(dim(susceptible_survival_mat))) {
    susceptible_survival_mat <- matrix(
      susceptible_survival_mat,
      nrow = nrow(pred_data),
      ncol = length(times)
    )
  }

  overall_survival_mat <- (1 - uncured_prob) + uncured_prob * susceptible_survival_mat

  list(
    susceptible_fraction = uncured_prob,
    susceptible_only_survival = susceptible_survival_mat,
    susceptible_only_risk = 1 - susceptible_survival_mat,
    overall_survival = overall_survival_mat,
    overall_risk = 1 - overall_survival_mat
  )
}

build_subject_prediction_df <- function(df, spec_row, fit, times) {
  pred <- predict_lognormal_mixture_cure(fit = fit, newdata = df, times = times)

  bind_rows(lapply(seq_along(times), function(j) {
    tibble(
      unique_person_id = as.character(df$unique_person_id),
      dataset = spec_row$dataset[[1]],
      dataset_label = spec_row$dataset_label[[1]],
      source_dataset = spec_row$source_dataset[[1]],
      model_id = model_id_value,
      site_adjustment_flag = as.logical(spec_row$site_adjustment_flag[[1]]),
      formula_name = spec_row$formula_name[[1]],
      formula_id = spec_row$formula_id[[1]],
      formula_label = spec_row$formula_label[[1]],
      incidence_rhs = spec_row$incidence_rhs[[1]],
      latency_rhs = spec_row$latency_rhs[[1]],
      site_in_incidence = as.logical(spec_row$site_in_incidence[[1]]),
      site_in_latency = as.logical(spec_row$site_in_latency[[1]]),
      risk_scale = main_risk_scale,
      time_year = as.numeric(times[[j]]),
      susceptible_fraction = as.numeric(pred$susceptible_fraction),
      susceptible_only_survival_prob = as.numeric(pred$susceptible_only_survival[, j]),
      susceptible_only_risk_prob = as.numeric(pred$susceptible_only_risk[, j]),
      overall_survival_prob = as.numeric(pred$overall_survival[, j]),
      overall_risk_prob = as.numeric(pred$overall_risk[, j])
    )
  }))
}

make_curve_plot <- function(curve_df, value_col, y_label, title_text) {
  yearly_points <- curve_df %>%
    filter(time_horizon_year %in% horizon_years)

  ggplot(curve_df, aes(x = time_horizon_year, y = .data[[value_col]], color = dataset)) +
    geom_line(linewidth = 1.05) +
    geom_point(
      data = yearly_points,
      aes(x = time_horizon_year, y = .data[[value_col]], color = dataset),
      size = 1.8
    ) +
    scale_color_manual(values = dataset_palette) +
    scale_x_continuous(
      breaks = 0:curve_horizon_max_year,
      limits = c(0, curve_horizon_max_year),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = title_text,
      subtitle = "Cohort-level means of subject-level predictions from frequentist log-normal mixture cure models",
      x = "Years after cohort entry",
      y = y_label,
      color = "Dataset"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
}

save_plot_png <- function(plot_object, output_file) {
  ggplot2::ggsave(
    filename = output_file,
    plot = plot_object,
    width = plot_width_in,
    height = plot_height_in,
    dpi = plot_dpi,
    units = "in"
  )
}

# Read Stage 1 Inputs -----------------------------------------------------
assert_file_exists(stage1_analysis_datasets_file, "Stage 1 analysis datasets RDS")
assert_file_exists(stage1_formula_registry_file, "Stage 1 formula registry CSV")

analysis_datasets <- readRDS(stage1_analysis_datasets_file)
formula_registry <- read_csv_checked(stage1_formula_registry_file, "Stage 1 formula registry CSV")

validate_inputs(analysis_datasets, formula_registry)
model_specs <- select_model_specs(formula_registry)

# Fit Models and Build Subject-Level Predictions --------------------------
curve_time_grid <- sort(unique(c(0, seq(0, curve_horizon_max_year, by = curve_step_year), horizon_years)))

fitted_object_cache <- load_fitted_object_cache(fit_object_rds_file)
subject_prediction_list <- vector("list", nrow(model_specs))

for (ii in seq_len(nrow(model_specs))) {
  spec_row <- model_specs[ii, , drop = FALSE]
  dataset_name <- spec_row$dataset[[1]]
  source_dataset_name <- spec_row$source_dataset[[1]]

  analysis_df <- prepare_dataset(
    df = analysis_datasets[[source_dataset_name]],
    dataset_name = dataset_name
  )

  cached_fit_obj <- fitted_object_cache[[dataset_name]]
  if (has_usable_cached_fit(cached_fit_obj, spec_row)) {
    fit_obj <- cached_fit_obj
    message(sprintf("[%s] Using cached fitted object from %s", dataset_name, basename(fit_object_rds_file)))
  } else {
    message(sprintf("[%s] Fitting log-normal mixture cure model and updating %s", dataset_name, basename(fit_object_rds_file)))
    fit_obj <- fit_lognormal_mixture_cure(
      df = analysis_df,
      incidence_rhs = spec_row$incidence_rhs[[1]],
      latency_rhs = spec_row$latency_rhs[[1]],
      dataset_name = dataset_name
    )
    fitted_object_cache[[dataset_name]] <- fit_obj
    saveRDS(fitted_object_cache, fit_object_rds_file)
  }

  subject_prediction_list[[ii]] <- build_subject_prediction_df(
    df = analysis_df,
    spec_row = spec_row,
    fit = fit_obj,
    times = curve_time_grid
  )
}

subject_prediction_df <- bind_rows(subject_prediction_list)

# Aggregate Cohort-Level Predictions --------------------------------------
cohort_curve_df <- subject_prediction_df %>%
  group_by(
    dataset,
    dataset_label,
    source_dataset,
    model_id,
    site_adjustment_flag,
    formula_name,
    formula_id,
    formula_label,
    incidence_rhs,
    latency_rhs,
    site_in_incidence,
    site_in_latency,
    risk_scale,
    time_year
  ) %>%
  summarise(
    susceptible_fraction = mean(susceptible_fraction, na.rm = TRUE),
    overall_survival_prob = mean(overall_survival_prob, na.rm = TRUE),
    overall_risk_prob = mean(overall_risk_prob, na.rm = TRUE),
    susceptible_only_survival_prob = mean(susceptible_only_survival_prob, na.rm = TRUE),
    susceptible_only_risk_prob = mean(susceptible_only_risk_prob, na.rm = TRUE),
    n_subjects = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(match(dataset, dataset_model_registry$dataset), time_year)

plot_source_df <- cohort_curve_df %>%
  transmute(
    dataset,
    model_formula = if_else(
      incidence_rhs == latency_rhs,
      paste0("incidence = latency = ", incidence_rhs),
      paste0("incidence = ", incidence_rhs, " ; latency = ", latency_rhs)
    ),
    time_horizon_year = time_year,
    susceptible_fraction,
    cure_fraction = 1 - susceptible_fraction,
    overall_survival_prob,
    overall_risk_prob,
    susceptible_only_survival_prob,
    susceptible_only_risk_prob
  ) %>%
  arrange(match(dataset, dataset_model_registry$dataset), time_horizon_year)

horizon_summary_df <- plot_source_df %>%
  filter(time_horizon_year %in% horizon_years)

# Build Plots -------------------------------------------------------------
overall_survival_plot <- make_curve_plot(
  curve_df = plot_source_df,
  value_col = "overall_survival_prob",
  y_label = "Overall survival probability",
  title_text = "Frequentist log-normal mixture cure: overall survival"
)

overall_risk_plot <- make_curve_plot(
  curve_df = plot_source_df,
  value_col = "overall_risk_prob",
  y_label = "Overall cumulative risk",
  title_text = "Frequentist log-normal mixture cure: overall cumulative risk"
)

susceptible_survival_plot <- make_curve_plot(
  curve_df = plot_source_df,
  value_col = "susceptible_only_survival_prob",
  y_label = "Susceptible-only survival probability",
  title_text = "Frequentist log-normal mixture cure: susceptible-only survival"
)

susceptible_risk_plot <- make_curve_plot(
  curve_df = plot_source_df,
  value_col = "susceptible_only_risk_prob",
  y_label = "Susceptible-only cumulative risk",
  title_text = "Frequentist log-normal mixture cure: susceptible-only cumulative risk"
)

plot_registry <- list(
  overall_survival_curve = overall_survival_plot,
  overall_risk_curve = overall_risk_plot,
  susceptible_only_survival_curve = susceptible_survival_plot,
  susceptible_only_risk_curve = susceptible_risk_plot
)

# Write Outputs -----------------------------------------------------------
readr::write_csv(horizon_summary_df, horizon_summary_file)
readr::write_csv(plot_source_df, plot_source_file)
saveRDS(fitted_object_cache, fit_object_rds_file)
saveRDS(plot_registry, plot_rds_file)

legacy_curve_data_file <- file.path(
  export_path,
  "stage7_simple_lognormal_mixture_cure_curve_data.csv"
)
if (file.exists(legacy_curve_data_file)) {
  unlink(legacy_curve_data_file)
}

save_plot_png(overall_survival_plot, overall_survival_png_file)
save_plot_png(overall_risk_plot, overall_risk_png_file)
save_plot_png(susceptible_survival_plot, susceptible_survival_png_file)
save_plot_png(susceptible_risk_plot, susceptible_risk_png_file)

message("Stage 7 simplified log-normal mixture cure export completed.")
message("Horizon summary CSV: ", normalizePath(horizon_summary_file, winslash = "/", mustWork = FALSE))
message("Plot source CSV: ", normalizePath(plot_source_file, winslash = "/", mustWork = FALSE))
message("Fitted object RDS: ", normalizePath(fit_object_rds_file, winslash = "/", mustWork = FALSE))
message("Plot RDS: ", normalizePath(plot_rds_file, winslash = "/", mustWork = FALSE))
