# Configure: raw-data Block 5 paths and constants ================================
dropbox_project_root <- Sys.getenv(
  "CHR_COHORTS_DROPBOX_ROOT",
  unset = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model"
)
merged_file <- Sys.getenv(
  "BLOCK5_RAW_MERGED_FILE",
  unset = file.path(
    dropbox_project_root,
    "0.Data",
    "2.Preprocessed data",
    "Preprocessed_Merged_PNUH_SNUH_Data.csv"
  )
)
legacy_stage6_export_path <- Sys.getenv(
  "BLOCK5_RAW_LEGACY_STAGE6_EXPORT_PATH",
  unset = file.path(dropbox_project_root, "old", "stage6_Cure-appropriateness screening")
)
export_path <- Sys.getenv(
  "BLOCK5_RAW_EXPORT_PATH",
  unset = file.path(dropbox_project_root, "5.Block5", "sub_raw_recompute")
)

main_output_dir <- export_path
supporting_output_dir <- file.path(export_path, "sub_supporting_tables")
metadata_output_dir <- file.path(export_path, "sub_metadata")
published_output_dir <- Sys.getenv(
  "BLOCK5_RAW_PUBLISHED_DIR",
  unset = dirname(export_path)
)
sync_published_outputs <- tolower(
  Sys.getenv("BLOCK5_RAW_SYNC_PUBLISHED_OUTPUTS", unset = "TRUE")
) %in% c("1", "true", "t", "yes", "y")

main_risk_scale <- "transition_only_main"
alpha_screening <- 0.05
receus_pi_threshold <- 0.025
receus_r_threshold <- 0.05
primary_receus_families <- c("lognormal")
sensitivity_receus_families <- c("exponential", "weibull", "gamma", "loglogistic")
import_legacy_xie_if_available <- tolower(
  Sys.getenv("BLOCK5_RAW_IMPORT_LEGACY_XIE", unset = "TRUE")
) %in% c("1", "true", "t", "yes", "y")

# Attach: packages and options ===================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(survival)
  library(tibble)
  library(purrr)
})

options(stringsAsFactors = FALSE, scipen = 999)

# Define: file helpers ===========================================================
assert_exists <- function(path, label = basename(path)) {
  if (!file.exists(path) && !dir.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

make_temp_export_path <- function(path) {
  ext <- tools::file_ext(path)
  tempfile(
    pattern = paste0(tools::file_path_sans_ext(basename(path)), "_"),
    tmpdir = dirname(path),
    fileext = if (nzchar(ext)) paste0(".", ext) else ""
  )
}

replace_file_atomically <- function(temp_path, final_path) {
  if (file.exists(final_path)) {
    unlink(final_path)
  }
  renamed <- file.rename(temp_path, final_path)
  if (isTRUE(renamed)) {
    return(TRUE)
  }
  copied <- file.copy(temp_path, final_path, overwrite = TRUE)
  if (isTRUE(copied)) {
    unlink(temp_path)
    return(TRUE)
  }
  FALSE
}

safe_write_csv_atomic <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_path <- make_temp_export_path(path)
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  readr::write_csv(df, temp_path)
  ok <- replace_file_atomically(temp_path, path)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to write CSV atomically: %s", path), call. = FALSE)
  }
  invisible(path)
}

safe_save_plot <- function(plot_object, path, width, height, dpi = 320, bg = "white") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_path <- make_temp_export_path(path)
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  ggplot2::ggsave(
    filename = temp_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = bg,
    limitsize = FALSE
  )
  ok <- replace_file_atomically(temp_path, path)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to write plot atomically: %s", path), call. = FALSE)
  }
  invisible(path)
}

safe_save_pdf_atomic <- function(path, width, height, plot_fun) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_path <- make_temp_export_path(path)
  pdf_open <- FALSE
  on.exit({
    if (pdf_open) {
      try(grDevices::dev.off(), silent = TRUE)
    }
    if (file.exists(temp_path)) {
      unlink(temp_path)
    }
  }, add = TRUE)

  grDevices::pdf(temp_path, width = width, height = height, onefile = TRUE)
  pdf_open <- TRUE
  plot_fun()
  grDevices::dev.off()
  pdf_open <- FALSE

  ok <- replace_file_atomically(temp_path, path)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to write PDF atomically: %s", path), call. = FALSE)
  }
  invisible(path)
}

standard_error_column_patterns <- c(
  "(^|_)sd$",
  "(^|_)se$",
  "std_error$",
  "std\\.error$",
  "stderr$",
  "posterior_sd$",
  "uncertainty_sd$",
  "robust_se$"
)

identify_standard_error_columns <- function(df) {
  col_names <- names(df)
  if (is.null(col_names) || length(col_names) == 0L) {
    return(character())
  }

  col_names[vapply(
    col_names,
    function(col_name) any(grepl(standard_error_column_patterns, col_name, ignore.case = TRUE)),
    logical(1)
  )]
}

identify_standard_error_id_columns <- function(df) {
  preferred_cols <- c(
    "dataset",
    "dataset_label",
    "source_dataset",
    "analysis_branch",
    "latency_family",
    "family",
    "selected_by_primary_aic",
    "receus_candidate_flag",
    "screening_class",
    "model_family"
  )

  intersect(preferred_cols, names(df))
}

empty_standard_error_registry <- function() {
  tibble(
    source_object = character(),
    n_rows = integer(),
    n_columns = integer(),
    n_standard_error_columns = integer(),
    standard_error_columns = character()
  )
}

empty_standard_error_long <- function() {
  tibble(
    source_object = character(),
    row_id = integer(),
    standard_error_column = character(),
    standard_error_value = numeric(),
    standard_error_value_raw = character()
  )
}

build_standard_error_registry_entry <- function(df, source_name) {
  se_cols <- identify_standard_error_columns(df)

  tibble(
    source_object = source_name,
    n_rows = nrow(df),
    n_columns = ncol(df),
    n_standard_error_columns = length(se_cols),
    standard_error_columns = if (length(se_cols) > 0L) paste(se_cols, collapse = "|") else NA_character_
  )
}

build_standard_error_long_table <- function(df, source_name) {
  se_cols <- identify_standard_error_columns(df)
  if (length(se_cols) == 0L || nrow(df) == 0L) {
    return(empty_standard_error_long())
  }

  base_df <- tibble::as_tibble(df) %>%
    mutate(row_id = dplyr::row_number())
  id_cols <- identify_standard_error_id_columns(base_df)

  bind_rows(lapply(se_cols, function(se_col) {
    out <- tibble(
      source_object = source_name,
      row_id = base_df$row_id,
      standard_error_column = se_col,
      standard_error_value = suppressWarnings(as.numeric(base_df[[se_col]])),
      standard_error_value_raw = as.character(base_df[[se_col]])
    )

    if (length(id_cols) > 0L) {
      out <- bind_cols(out, base_df[, id_cols, drop = FALSE])
    }

    out
  }))
}

build_standard_error_export_bundle <- function(named_tables) {
  named_tables <- named_tables[vapply(named_tables, function(x) inherits(x, "data.frame"), logical(1))]

  registry_tbl <- bind_rows(lapply(names(named_tables), function(table_name) {
    build_standard_error_registry_entry(named_tables[[table_name]], table_name)
  }))
  long_tbl <- bind_rows(lapply(names(named_tables), function(table_name) {
    build_standard_error_long_table(named_tables[[table_name]], table_name)
  }))

  if (nrow(registry_tbl) == 0L) {
    registry_tbl <- empty_standard_error_registry()
  }
  if (nrow(long_tbl) == 0L) {
    long_tbl <- empty_standard_error_long()
  }

  list(
    registry = registry_tbl,
    long = long_tbl
  )
}

# Define: utility helpers ========================================================
coerce_numeric_text <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

bounded_prob <- function(x, eps = 1e-8) {
  pmin(pmax(x, eps), 1 - eps)
}

logspace_add <- function(logx, logy) {
  m <- pmax(logx, logy)
  m + log(exp(logx - m) + exp(logy - m))
}

safe_scalar_character <- function(x, default = NA_character_) {
  if (length(x) == 0L || is.null(x)) {
    return(default)
  }
  x0 <- as.character(x[[1]])
  if (is.na(x0) || !nzchar(x0)) default else x0
}

format_fixed_value <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

# Define: ingestion and dataset preparation ======================================
read_main_dataset <- function(path) {
  assert_exists(path, label = "Merged preprocessed dataset")
  readr::read_csv(
    file = path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
}

prepare_transition_dataset <- function(df, dataset_key, site_mode = c("single", "merged")) {
  site_mode <- match.arg(site_mode)
  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(sprintf("[%s] Missing required columns: %s", dataset_key, paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  out <- df %>%
    mutate(
      id = trimws(as.character(id)),
      site = trimws(as.character(site)),
      sex_num = as.integer(coerce_numeric_text(sex_num)),
      age_exact_entry = coerce_numeric_text(age_exact_entry),
      days_followup = coerce_numeric_text(days_followup),
      status_num = as.integer(coerce_numeric_text(status_num))
    ) %>%
    mutate(
      time_year = pmax(days_followup / 365.25, 1e-8),
      event_transition = as.integer(status_num == 1L),
      dataset_key = dataset_key,
      risk_scale = main_risk_scale
    )

  if (nrow(out) == 0L) {
    stop(sprintf("[%s] Dataset has zero rows.", dataset_key), call. = FALSE)
  }
  if (anyNA(out[, required_cols])) {
    stop(sprintf("[%s] Missing values in required backbone columns.", dataset_key), call. = FALSE)
  }
  if (any(out$days_followup < 0)) {
    stop(sprintf("[%s] Negative follow-up detected.", dataset_key), call. = FALSE)
  }
  if (any(!out$sex_num %in% c(0L, 1L))) {
    stop(sprintf("[%s] sex_num must be coded as 0/1.", dataset_key), call. = FALSE)
  }
  if (any(!out$status_num %in% c(0L, 1L, 2L))) {
    stop(sprintf("[%s] status_num must be coded as 0/1/2.", dataset_key), call. = FALSE)
  }

  n_site_levels <- dplyr::n_distinct(out$site)
  if (site_mode == "single" && n_site_levels != 1L) {
    stop(sprintf("[%s] Expected exactly one site for a single-site dataset.", dataset_key), call. = FALSE)
  }
  if (site_mode == "merged" && n_site_levels < 2L) {
    stop(sprintf("[%s] Expected at least two sites for a merged dataset.", dataset_key), call. = FALSE)
  }

  out
}

# Define: survival kernels =======================================================
family_npar <- function(family) {
  switch(
    family,
    exponential = 1L,
    weibull = 2L,
    gamma = 2L,
    loglogistic = 2L,
    lognormal = 2L,
    stop(sprintf("Unsupported family: %s", family), call. = FALSE)
  )
}

family_log_density <- function(time, family, theta_raw) {
  switch(
    family,
    exponential = {
      rate <- exp(theta_raw[1])
      dexp(time, rate = rate, log = TRUE)
    },
    weibull = {
      shape <- exp(theta_raw[1])
      scale <- exp(theta_raw[2])
      dweibull(time, shape = shape, scale = scale, log = TRUE)
    },
    gamma = {
      shape <- exp(theta_raw[1])
      scale <- exp(theta_raw[2])
      dgamma(time, shape = shape, scale = scale, log = TRUE)
    },
    loglogistic = {
      shape <- exp(theta_raw[1])
      scale <- exp(theta_raw[2])
      z <- pmax(time / scale, 1e-12)
      log(shape) - log(scale) + (shape - 1) * log(z) - 2 * log1p(z^shape)
    },
    lognormal = {
      meanlog <- theta_raw[1]
      sdlog <- exp(theta_raw[2])
      dlnorm(time, meanlog = meanlog, sdlog = sdlog, log = TRUE)
    },
    stop(sprintf("Unsupported family: %s", family), call. = FALSE)
  )
}

family_log_survival <- function(time, family, theta_raw) {
  switch(
    family,
    exponential = {
      rate <- exp(theta_raw[1])
      pexp(time, rate = rate, lower.tail = FALSE, log.p = TRUE)
    },
    weibull = {
      shape <- exp(theta_raw[1])
      scale <- exp(theta_raw[2])
      pweibull(time, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
    },
    gamma = {
      shape <- exp(theta_raw[1])
      scale <- exp(theta_raw[2])
      pgamma(time, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
    },
    loglogistic = {
      shape <- exp(theta_raw[1])
      scale <- exp(theta_raw[2])
      z <- pmax(time / scale, 1e-12)
      -log1p(z^shape)
    },
    lognormal = {
      meanlog <- theta_raw[1]
      sdlog <- exp(theta_raw[2])
      plnorm(time, meanlog = meanlog, sdlog = sdlog, lower.tail = FALSE, log.p = TRUE)
    },
    stop(sprintf("Unsupported family: %s", family), call. = FALSE)
  )
}

family_survival <- function(time, family, theta_raw) {
  exp(family_log_survival(time = time, family = family, theta_raw = theta_raw))
}

# Define: likelihoods and optimizers =============================================
make_family_start_values <- function(time, event, family, cure = FALSE) {
  uncensored_time <- time[event == 1L]
  uncensored_time <- uncensored_time[is.finite(uncensored_time) & uncensored_time > 0]
  if (length(uncensored_time) == 0L) {
    uncensored_time <- time[is.finite(time) & time > 0]
  }

  m_time <- stats::median(uncensored_time)
  if (!is.finite(m_time) || m_time <= 0) {
    m_time <- stats::median(time[time > 0], na.rm = TRUE)
  }
  if (!is.finite(m_time) || m_time <= 0) {
    m_time <- 1
  }

  base_theta <- switch(
    family,
    exponential = c(log(1 / max(mean(uncensored_time), 1e-4))),
    weibull = c(log(1), log(max(m_time, 1e-4))),
    gamma = c(log(1), log(max(m_time, 1e-4))),
    loglogistic = c(log(1), log(max(m_time, 1e-4))),
    lognormal = c(log(max(m_time, 1e-4)), log(1)),
    stop(sprintf("Unsupported family: %s", family), call. = FALSE)
  )

  if (!cure) {
    return(list(base_theta))
  }

  km_fit <- survival::survfit(survival::Surv(time, event) ~ 1)
  km_tail <- if (length(km_fit$surv) == 0L) 0.05 else max(min(tail(km_fit$surv, 1), 0.90), 0.001)
  cure_grid <- unique(c(0.02, 0.05, 0.10, 0.20, 0.35, km_tail))
  cure_grid <- cure_grid[cure_grid > 0 & cure_grid < 0.95]
  starts <- vector("list", length(cure_grid))
  for (i in seq_along(cure_grid)) {
    starts[[i]] <- c(qlogis(bounded_prob(cure_grid[i])), base_theta)
  }
  starts
}

negloglik_noncure <- function(par, time, event, family) {
  logf <- suppressWarnings(family_log_density(time = time, family = family, theta_raw = par))
  logS <- suppressWarnings(family_log_survival(time = time, family = family, theta_raw = par))
  val <- -(sum(event * logf + (1 - event) * logS))
  if (!is.finite(val)) 1e50 else val
}

negloglik_cure <- function(par, time, event, family) {
  eta <- par[1]
  theta <- par[-1]
  p_cure <- bounded_prob(plogis(eta))
  logf <- suppressWarnings(family_log_density(time = time, family = family, theta_raw = theta))
  logS_uncured <- suppressWarnings(family_log_survival(time = time, family = family, theta_raw = theta))
  log_event <- log1p(-p_cure) + logf
  log_censor <- logspace_add(log(p_cure), log1p(-p_cure) + logS_uncured)
  val <- -(sum(event * log_event + (1 - event) * log_censor))
  if (!is.finite(val)) 1e50 else val
}

run_best_optim <- function(start_list, objective_fn, time, event, family) {
  best <- NULL
  for (start in start_list) {
    fit_try <- try(
      stats::optim(
        par = start,
        fn = objective_fn,
        time = time,
        event = event,
        family = family,
        method = "BFGS",
        hessian = FALSE,
        control = list(maxit = 5000, reltol = 1e-10)
      ),
      silent = TRUE
    )
    if (inherits(fit_try, "try-error")) {
      next
    }
    current <- list(
      par = fit_try$par,
      value = fit_try$value,
      convergence = fit_try$convergence,
      message = safe_scalar_character(fit_try$message)
    )
    if (is.null(best) || (is.finite(current$value) && current$value < best$value)) {
      best <- current
    }
  }
  best
}

fit_survival_family <- function(df, family, model_type = c("noncure", "cure")) {
  model_type <- match.arg(model_type)
  time <- df$time_year
  event <- df$event_transition
  starts <- make_family_start_values(time = time, event = event, family = family, cure = identical(model_type, "cure"))
  best <- run_best_optim(
    start_list = starts,
    objective_fn = if (identical(model_type, "cure")) negloglik_cure else negloglik_noncure,
    time = time,
    event = event,
    family = family
  )

  k <- family_npar(family) + if (identical(model_type, "cure")) 1L else 0L
  if (is.null(best)) {
    return(list(
      model_type = model_type,
      family = family,
      converged = FALSE,
      logLik = NA_real_,
      k = k,
      AIC = NA_real_,
      par = NA_real_,
      cure_fraction_hat = if (identical(model_type, "cure")) NA_real_ else 0,
      uncured_survival_fn = function(t) rep(NA_real_, length(t)),
      total_survival_fn = function(t) rep(NA_real_, length(t)),
      detail_message = "Optimization failed for all starting values."
    ))
  }

  logLik <- -best$value
  if (identical(model_type, "noncure")) {
    return(list(
      model_type = model_type,
      family = family,
      converged = isTRUE(best$convergence == 0L) && is.finite(logLik),
      logLik = logLik,
      k = k,
      AIC = 2 * k - 2 * logLik,
      par = best$par,
      cure_fraction_hat = 0,
      uncured_survival_fn = function(t) family_survival(time = pmax(t, 1e-8), family = family, theta_raw = best$par),
      total_survival_fn = function(t) family_survival(time = pmax(t, 1e-8), family = family, theta_raw = best$par),
      detail_message = best$message
    ))
  }

  eta <- best$par[1]
  theta <- best$par[-1]
  p_cure <- bounded_prob(plogis(eta))
  list(
    model_type = model_type,
    family = family,
    converged = isTRUE(best$convergence == 0L) && is.finite(logLik),
    logLik = logLik,
    k = k,
    AIC = 2 * k - 2 * logLik,
    par = best$par,
    cure_fraction_hat = p_cure,
    uncured_survival_fn = function(t) family_survival(time = pmax(t, 1e-8), family = family, theta_raw = theta),
    total_survival_fn = function(t) p_cure + (1 - p_cure) * family_survival(time = pmax(t, 1e-8), family = family, theta_raw = theta),
    detail_message = best$message
  )
}

# Define: screening summaries ====================================================
compute_km_tail_survival <- function(df) {
  fit <- survival::survfit(survival::Surv(time_year, event_transition) ~ 1, data = df)
  if (length(fit$surv) == 0L) {
    return(1)
  }
  as.numeric(tail(fit$surv, 1))
}

compute_tail_summary <- function(df) {
  t_n <- max(df$time_year, na.rm = TRUE)
  event_times <- df$time_year[df$event_transition == 1L]
  t_star <- if (length(event_times) == 0L) NA_real_ else max(event_times, na.rm = TRUE)
  tibble(
    n_total = nrow(df),
    n_event = sum(df$event_transition == 1L),
    n_censor_main = sum(df$event_transition == 0L),
    tau_year = t_n,
    largest_event_time_year = t_star,
    plateau_length_year = ifelse(is.na(t_star), NA_real_, t_n - t_star),
    km_tail_survival = compute_km_tail_survival(df),
    km_tail_noncure = 1 - compute_km_tail_survival(df)
  )
}

classify_receus_primary <- function(model_type, cure_fraction_hat, receus_ratio_hat) {
  if (!identical(model_type, "cure")) {
    return("unsupportive")
  }
  if (is.na(cure_fraction_hat) || is.na(receus_ratio_hat)) {
    return("unsupportive")
  }
  if (cure_fraction_hat > receus_pi_threshold && receus_ratio_hat < receus_r_threshold) {
    return("supportive")
  }
  "unsupportive"
}

compute_maller_zhou_dn <- function(cure_fit, noncure_fit) {
  if (is.null(cure_fit) || is.null(noncure_fit) || !isTRUE(cure_fit$converged) || !isTRUE(noncure_fit$converged)) {
    return(tibble(d_n = NA_real_, p_value = NA_real_, contextual_signal_flag = NA))
  }
  d_n <- 2 * (cure_fit$logLik - noncure_fit$logLik)
  if (!is.finite(d_n) || d_n < 0) {
    d_n <- max(d_n, 0)
  }
  p_value <- if (d_n <= 0) 1 else 0.5 * stats::pchisq(d_n, df = 1, lower.tail = FALSE)
  tibble(
    d_n = d_n,
    p_value = p_value,
    contextual_signal_flag = !is.na(p_value) & p_value < alpha_screening
  )
}

read_legacy_xie_summary <- function(path) {
  if (!isTRUE(import_legacy_xie_if_available)) {
    return(NULL)
  }
  summary_file <- file.path(path, "stage6_screening_summary.csv")
  if (!file.exists(summary_file)) {
    return(NULL)
  }
  xie_raw <- readr::read_csv(summary_file, show_col_types = FALSE, progress = FALSE)
  xie_raw %>%
    filter(dataset_key %in% c("PNU", "SNU", "merged__site_free")) %>%
    transmute(
      dataset = dplyr::case_when(
        dataset_key == "PNU" ~ "PNU",
        dataset_key == "SNU" ~ "SNU",
        TRUE ~ "merged"
      ),
      xie_statistic = as.numeric(xie_T_n),
      xie_bootstrap_p_value = as.numeric(xie_centered_bootstrap_p_value),
      xie_contradiction_flag = !is.na(xie_bootstrap_p_value) & xie_bootstrap_p_value < alpha_screening,
      xie_method_label = "legacy_stage6_centered_bootstrap",
      xie_interpretation_note = dplyr::case_when(
        !is.na(xie_bootstrap_p_value) & xie_bootstrap_p_value < alpha_screening ~ "Legacy Stage 6 Xie bootstrap output indicates a contradiction signal against sufficient follow-up.",
        !is.na(xie_bootstrap_p_value) ~ "Legacy Stage 6 Xie bootstrap output does not indicate a contradiction signal against sufficient follow-up.",
        TRUE ~ "Legacy Stage 6 Xie bootstrap output is unavailable."
      )
    )
}

make_xie_status_row <- function(dataset, legacy_xie_tbl) {
  if (!is.null(legacy_xie_tbl) && dataset %in% legacy_xie_tbl$dataset) {
    return(legacy_xie_tbl %>% filter(dataset == .env$dataset) %>% slice(1))
  }
  tibble(
    dataset = dataset,
    xie_statistic = NA_real_,
    xie_bootstrap_p_value = NA_real_,
    xie_contradiction_flag = NA,
    xie_method_label = "not_implemented_in_raw_recompute_script",
    xie_interpretation_note = "This raw-data recompute script does not approximate Xie with a proxy; supply legacy Xie outputs or implement the native Xie bootstrap separately."
  )
}

# Define: dataset-level workflow =================================================
fit_family_pool <- function(df, families) {
  fit_rows <- vector("list", length(families) * 2L)
  nm <- character(length(fit_rows))
  idx <- 1L
  for (family in families) {
    fit_rows[[idx]] <- fit_survival_family(df = df, family = family, model_type = "noncure")
    nm[[idx]] <- paste0("noncure__", family)
    idx <- idx + 1L
    fit_rows[[idx]] <- fit_survival_family(df = df, family = family, model_type = "cure")
    nm[[idx]] <- paste0("cure__", family)
    idx <- idx + 1L
  }
  names(fit_rows) <- nm
  fit_rows
}

fit_list_to_table <- function(dataset, fit_list, pool_name) {
  purrr::imap_dfr(fit_list, function(fit_obj, fit_name) {
    tibble(
      dataset = dataset,
      selection_pool = pool_name,
      fit_name = fit_name,
      model_type = fit_obj$model_type,
      latency_family = fit_obj$family,
      converged = isTRUE(fit_obj$converged),
      logLik = fit_obj$logLik,
      effective_parameter_count = fit_obj$k,
      AIC = fit_obj$AIC,
      cure_fraction_hat = fit_obj$cure_fraction_hat,
      detail_message = fit_obj$detail_message
    )
  })
}

run_block5_raw_for_dataset <- function(df, dataset, legacy_xie_tbl = NULL) {
  primary_fit_list <- fit_family_pool(df = df, families = primary_receus_families)
  primary_fit_table <- fit_list_to_table(dataset = dataset, fit_list = primary_fit_list, pool_name = "primary") %>%
    mutate(
      selected_aic = is.finite(AIC) & AIC == min(AIC[is.finite(AIC)], na.rm = TRUE)
    )

  selected_row <- primary_fit_table %>% filter(selected_aic) %>% slice(1)
  if (nrow(selected_row) == 0L) {
    stop(sprintf("[%s] No converged primary RECeUS candidate model available.", dataset), call. = FALSE)
  }

  selected_fit <- primary_fit_list[[selected_row$fit_name[[1]]]]
  tau_year <- max(df$time_year, na.rm = TRUE)
  receus_ratio_hat <- if (identical(selected_fit$model_type, "cure")) {
    s_u <- selected_fit$uncured_survival_fn(tau_year)
    s_t <- selected_fit$total_survival_fn(tau_year)
    as.numeric(s_u / s_t)
  } else {
    NA_real_
  }
  cure_fraction_hat <- if (identical(selected_fit$model_type, "cure")) selected_fit$cure_fraction_hat else 0
  receus_primary_flag <- classify_receus_primary(
    model_type = selected_fit$model_type,
    cure_fraction_hat = cure_fraction_hat,
    receus_ratio_hat = receus_ratio_hat
  )

  receus_primary_summary <- tibble(
    dataset = dataset,
    selected_model_type = selected_fit$model_type,
    selected_family = selected_fit$family,
    selected_aic = selected_row$AIC[[1]],
    analysis_time_year = tau_year,
    cure_fraction_hat = cure_fraction_hat,
    receus_ratio_hat = receus_ratio_hat,
    receus_pi_threshold = receus_pi_threshold,
    receus_r_threshold = receus_r_threshold,
    receus_primary_flag = receus_primary_flag,
    screening_note = dplyr::case_when(
      identical(selected_fit$model_type, "noncure") ~ sprintf(
        "Prespecified lognormal comparison selected the non-cure %s model, so the main-model-family RECeUS screen is unsupportive by construction.",
        selected_fit$family
      ),
      receus_primary_flag == "supportive" ~ sprintf(
        "Prespecified lognormal comparison selected a cure %s model with cure_fraction_hat=%.4f and receus_ratio_hat=%.4f, so the main-model-family RECeUS screen is supportive.",
        selected_fit$family,
        cure_fraction_hat,
        receus_ratio_hat
      ),
      TRUE ~ sprintf(
        "Prespecified lognormal comparison selected a cure %s model with cure_fraction_hat=%.4f and receus_ratio_hat=%.4f, but the main-model-family RECeUS screen is unsupportive.",
        selected_fit$family,
        cure_fraction_hat,
        receus_ratio_hat
      )
    )
  )

  sensitivity_fit_table <- NULL
  sensitivity_fit_list <- list()
  if (length(sensitivity_receus_families) > 0L) {
    sensitivity_fit_list <- fit_family_pool(df = df, families = sensitivity_receus_families)
    sensitivity_fit_table <- fit_list_to_table(
      dataset = dataset,
      fit_list = sensitivity_fit_list,
      pool_name = "sensitivity"
    )
  }

  cure_candidate_lists <- c(primary_fit_list, sensitivity_fit_list)
  cure_candidate_summary <- purrr::imap_dfr(cure_candidate_lists, function(fit_obj, fit_name) {
    if (!startsWith(fit_name, "cure__")) {
      return(NULL)
    }
    receus_ratio_hat_candidate <- if (isTRUE(fit_obj$converged)) {
      s_u <- fit_obj$uncured_survival_fn(tau_year)
      s_t <- fit_obj$total_survival_fn(tau_year)
      as.numeric(s_u / s_t)
    } else {
      NA_real_
    }
    cure_fraction_hat_candidate <- if (isTRUE(fit_obj$converged)) fit_obj$cure_fraction_hat else NA_real_
    tibble(
      dataset = dataset,
      fit_name = fit_name,
      selection_pool = ifelse(fit_name %in% names(primary_fit_list), "primary", "sensitivity"),
      latency_family = fit_obj$family,
      converged = isTRUE(fit_obj$converged),
      cure_fraction_hat = cure_fraction_hat_candidate,
      receus_ratio_hat = receus_ratio_hat_candidate,
      receus_candidate_flag = classify_receus_primary(
        model_type = "cure",
        cure_fraction_hat = cure_fraction_hat_candidate,
        receus_ratio_hat = receus_ratio_hat_candidate
      ),
      selected_by_primary_aic = identical(selected_fit$model_type, "cure") && identical(fit_name, selected_row$fit_name[[1]])
    )
  })

  maller_fit_list <- c(primary_fit_list, sensitivity_fit_list)
  maller_tbl <- compute_maller_zhou_dn(
    cure_fit = maller_fit_list[["cure__exponential"]],
    noncure_fit = maller_fit_list[["noncure__exponential"]]
  ) %>%
    mutate(
      dataset = dataset,
      latency_family = "exponential",
      contextual_note = dplyr::case_when(
        contextual_signal_flag %in% TRUE ~ "Exponential working-model Maller-Zhou boundary test is positive, but remains contextual only.",
        contextual_signal_flag %in% FALSE ~ "Exponential working-model Maller-Zhou boundary test is not positive and remains contextual only.",
        TRUE ~ "Exponential working-model Maller-Zhou boundary test could not be evaluated and remains contextual only."
      )
    ) %>%
    select(dataset, latency_family, d_n, p_value, contextual_signal_flag, contextual_note)

  tail_summary <- compute_tail_summary(df) %>%
    mutate(dataset = dataset, .before = 1)

  xie_status <- make_xie_status_row(dataset = dataset, legacy_xie_tbl = legacy_xie_tbl)

  screening_summary <- tail_summary %>%
    left_join(receus_primary_summary, by = "dataset") %>%
    left_join(maller_tbl, by = "dataset") %>%
    left_join(xie_status, by = "dataset") %>%
    mutate(
      block5_screening_class = dplyr::case_when(
        receus_primary_flag == "supportive" & xie_contradiction_flag %in% TRUE ~ "primary_supportive_but_xie_contradicted",
        receus_primary_flag == "supportive" ~ "primary_supportive",
        TRUE ~ "primary_unsupportive"
      ),
      block5_screening_note = dplyr::case_when(
        receus_primary_flag == "supportive" & xie_contradiction_flag %in% TRUE ~ "Main-model-family lognormal RECeUS screen is supportive, but Xie contradicts sufficient follow-up.",
        receus_primary_flag == "supportive" ~ "Main-model-family lognormal RECeUS screen is supportive and Xie does not add a contradiction signal.",
        TRUE ~ "Main-model-family lognormal RECeUS screen is unsupportive; Xie and Maller-Zhou remain secondary/contextual."
      )
    )

  list(
    primary_fit_table = primary_fit_table,
    sensitivity_fit_table = sensitivity_fit_table,
    cure_candidate_summary = cure_candidate_summary,
    receus_primary_summary = receus_primary_summary,
    maller_tbl = maller_tbl,
    xie_status = xie_status,
    tail_summary = tail_summary,
    screening_summary = screening_summary
  )
}

# Execute: data loading ==========================================================
raw_merged <- read_main_dataset(merged_file) %>%
  mutate(site = trimws(as.character(site)))

dataset_registry <- tibble(
  dataset = c("PNU", "SNU", "merged"),
  data_slice = c("single_site", "single_site", "merged_site_free"),
  site_rule = c("site == PNU", "site == SNU", "all rows"),
  risk_scale = main_risk_scale,
  note = c(
    "Single-site Block 5 raw-data recomputation on the transition-only scale.",
    "Single-site Block 5 raw-data recomputation on the transition-only scale.",
    "Merged raw-tail Block 5 recomputation without a site-adjusted screening branch."
  )
)

dataset_inputs <- list(
  PNU = raw_merged %>% filter(toupper(site) == "PNU") %>% prepare_transition_dataset(dataset_key = "PNU", site_mode = "single"),
  SNU = raw_merged %>% filter(toupper(site) == "SNU") %>% prepare_transition_dataset(dataset_key = "SNU", site_mode = "single"),
  merged = raw_merged %>% prepare_transition_dataset(dataset_key = "merged", site_mode = "merged")
)

legacy_xie_tbl <- read_legacy_xie_summary(legacy_stage6_export_path)

# Execute: per-dataset raw recompute =============================================
result_list <- purrr::imap(dataset_inputs, function(df, dataset) {
  run_block5_raw_for_dataset(df = df, dataset = dataset, legacy_xie_tbl = legacy_xie_tbl)
})

primary_fit_summary <- bind_rows(purrr::map(result_list, "primary_fit_table"))
sensitivity_fit_summary <- bind_rows(purrr::compact(purrr::map(result_list, "sensitivity_fit_table")))
cure_candidate_summary <- bind_rows(purrr::map(result_list, "cure_candidate_summary"))
receus_primary_summary <- bind_rows(purrr::map(result_list, "receus_primary_summary"))
maller_zhou_contextual_summary <- bind_rows(purrr::map(result_list, "maller_tbl"))
xie_status_summary <- bind_rows(purrr::map(result_list, "xie_status"))
tail_summary <- bind_rows(purrr::map(result_list, "tail_summary"))
screening_summary <- bind_rows(purrr::map(result_list, "screening_summary"))

# Visualize: raw recompute summaries ============================================
dataset_levels <- c("PNU", "SNU", "merged")
component_levels <- c("Lognormal primary", "Xie contradiction", "Maller-Zhou contextual")

model_fit_plot_data <- primary_fit_summary %>%
  mutate(
    dataset_order = match(dataset, dataset_levels),
    model_label = paste(model_type, latency_family, sep = "::")
  )

model_fit_plot <- model_fit_plot_data %>%
  mutate(dataset = factor(dataset, levels = dataset_levels)) %>%
  ggplot(aes(x = reorder(model_label, AIC), y = AIC, fill = model_type)) +
  geom_col() +
  facet_wrap(~dataset, scales = "free_y") +
  geom_point(
    data = subset(model_fit_plot_data, selected_aic),
    aes(x = model_label, y = AIC),
    color = "black",
    size = 2.5
  ) +
  coord_flip() +
  labs(
    title = "Block 5 raw-data recompute: prespecified lognormal primary model pool",
    x = "Candidate model",
    y = "AIC"
  ) +
  theme_bw(base_size = 11)

receus_primary_map_plot_data <- receus_primary_summary %>%
  mutate(
    dataset_order = match(dataset, dataset_levels),
    receus_ratio_plot = ifelse(is.na(receus_ratio_hat), receus_r_threshold * 1.1, receus_ratio_hat),
    plot_shape = ifelse(selected_model_type == "noncure", "noncure_selected", "cure_selected"),
    plot_label = ifelse(selected_model_type == "noncure", paste0(dataset, " (non-cure)"), dataset)
  )

receus_primary_map <- receus_primary_map_plot_data %>%
  ggplot(aes(x = cure_fraction_hat, y = receus_ratio_plot, color = receus_primary_flag, shape = plot_shape, label = plot_label)) +
  geom_hline(yintercept = receus_r_threshold, linetype = 2) +
  geom_vline(xintercept = receus_pi_threshold, linetype = 2) +
  geom_point(size = 3) +
  geom_text(nudge_y = 0.015, size = 3.2, show.legend = FALSE) +
  scale_shape_manual(values = c(cure_selected = 16, noncure_selected = 17)) +
  labs(
    title = "Block 5 raw-data recompute: prespecified lognormal primary-screen map",
    x = "Estimated cure fraction",
    y = "RECeUS ratio",
    shape = NULL
  ) +
  theme_bw(base_size = 11)

screening_overview_plot_data <- bind_rows(
  screening_summary %>%
    transmute(
      dataset,
      component = "Lognormal primary",
      component_status = paste("Primary", receus_primary_flag),
      component_value_label = sprintf(
        "AIC=%s\npi=%s, r=%s",
        format_fixed_value(selected_aic, digits = 2),
        format_fixed_value(cure_fraction_hat, digits = 3),
        ifelse(is.na(receus_ratio_hat), "NA", format_fixed_value(receus_ratio_hat, digits = 3))
      )
    ),
  screening_summary %>%
    transmute(
      dataset,
      component = "Xie contradiction",
      component_status = ifelse(xie_contradiction_flag %in% TRUE, "Xie contradiction", "Xie no contradiction"),
      component_value_label = sprintf(
        "T=%s\np=%s",
        format_fixed_value(xie_statistic, digits = 3),
        format_fixed_value(xie_bootstrap_p_value, digits = 3)
      )
    ),
  screening_summary %>%
    transmute(
      dataset,
      component = "Maller-Zhou contextual",
      component_status = ifelse(contextual_signal_flag %in% TRUE, "Contextual signal present", "Contextual signal absent"),
      component_value_label = sprintf(
        "d=%s\np=%s",
        format_fixed_value(d_n, digits = 3),
        format_fixed_value(p_value, digits = 3)
      )
    )
) %>%
  mutate(
    dataset_order = match(dataset, dataset_levels),
    component_order = match(component, component_levels)
  )

screening_overview_plot <- screening_overview_plot_data %>%
  mutate(
    dataset = factor(dataset, levels = dataset_levels),
    component = factor(component, levels = component_levels)
  ) %>%
  ggplot(aes(x = component, y = dataset, fill = component_status)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = component_value_label), color = "white", size = 3.1, lineheight = 0.95, fontface = "bold") +
  scale_fill_manual(
    values = c(
      "Primary supportive" = "#54A24B",
      "Primary unsupportive" = "#E45756",
      "Xie contradiction" = "#E45756",
      "Xie no contradiction" = "#4C78A8",
      "Contextual signal present" = "#B279A2",
      "Contextual signal absent" = "#BDBDBD"
    )
  ) +
  labs(
    title = "Block 5 raw-data recompute: screening overview",
    subtitle = "Prespecified lognormal primary screen, legacy Xie contradiction check, and contextual Maller-Zhou signal",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

receus_candidate_map_plot_data <- cure_candidate_summary %>%
  filter(is.finite(cure_fraction_hat), is.finite(receus_ratio_hat)) %>%
  mutate(
    dataset_order = match(dataset, dataset_levels),
    receus_pi_threshold = receus_pi_threshold,
    receus_r_threshold = receus_r_threshold
  )

receus_candidate_map <- receus_candidate_map_plot_data %>%
  mutate(dataset = factor(dataset, levels = dataset_levels)) %>%
  ggplot(aes(x = receus_ratio_hat, y = cure_fraction_hat, color = receus_candidate_flag, shape = selected_by_primary_aic)) +
  geom_vline(xintercept = receus_r_threshold, linewidth = 0.4, linetype = "dashed", color = "#666666") +
  geom_hline(yintercept = receus_pi_threshold, linewidth = 0.4, linetype = "dashed", color = "#666666") +
  geom_point(size = 3) +
  geom_text(aes(label = latency_family), nudge_y = 0.02, size = 3, check_overlap = TRUE, show.legend = FALSE) +
  facet_wrap(~dataset, nrow = 1) +
  scale_color_manual(values = c("supportive" = "#54A24B", "unsupportive" = "#E45756")) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17)) +
  labs(
    title = "Block 5 raw-data recompute: cure-candidate RECeUS map",
    subtitle = "Published primary uses prespecified lognormal; other cure families are shown as exploratory sensitivity candidates",
    x = "RECeUS ratio",
    y = "Estimated cure fraction",
    color = "RECeUS flag",
    shape = "Primary AIC-selected cure fit"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

# Export: outputs ================================================================
output_paths <- list(
  screening_summary = file.path(main_output_dir, "block5_raw_screening_summary.csv"),
  model_fit_plot = file.path(main_output_dir, "block5_raw_model_fit_summary.png"),
  receus_primary_map = file.path(main_output_dir, "block5_raw_receus_primary_map.png"),
  screening_overview_plot = file.path(main_output_dir, "block5_raw_screening_overview.png"),
  receus_candidate_map = file.path(main_output_dir, "block5_raw_receus_candidate_map.png"),
  visual_summary_pdf = file.path(main_output_dir, "block5_raw_visual_summary.pdf"),
  dataset_registry = file.path(supporting_output_dir, "block5_raw_dataset_registry.csv"),
  primary_fit_summary = file.path(supporting_output_dir, "block5_raw_primary_model_fit_summary.csv"),
  model_fit_plot_data = file.path(supporting_output_dir, "block5_raw_model_fit_plot_data.csv"),
  sensitivity_fit_summary = file.path(supporting_output_dir, "block5_raw_sensitivity_model_fit_summary.csv"),
  cure_candidate_summary = file.path(supporting_output_dir, "block5_raw_cure_candidate_summary.csv"),
  receus_candidate_map_plot_data = file.path(supporting_output_dir, "block5_raw_receus_candidate_map_plot_data.csv"),
  receus_primary_summary = file.path(supporting_output_dir, "block5_raw_receus_primary_summary.csv"),
  receus_primary_map_plot_data = file.path(supporting_output_dir, "block5_raw_receus_primary_map_plot_data.csv"),
  maller_zhou_contextual_summary = file.path(supporting_output_dir, "block5_raw_maller_zhou_contextual_summary.csv"),
  xie_status_summary = file.path(supporting_output_dir, "block5_raw_xie_status_summary.csv"),
  tail_summary = file.path(supporting_output_dir, "block5_raw_tail_summary.csv"),
  screening_overview_plot_data = file.path(supporting_output_dir, "block5_raw_screening_overview_plot_data.csv"),
  standard_error_registry = file.path(supporting_output_dir, "block5_raw_standard_error_table_registry.csv"),
  standard_error_long = file.path(supporting_output_dir, "block5_raw_standard_error_long.csv"),
  export_manifest = file.path(metadata_output_dir, "block5_raw_export_manifest.csv")
)

block5_raw_standard_error_bundle <- build_standard_error_export_bundle(list(
  screening_summary = screening_summary,
  dataset_registry = dataset_registry,
  primary_fit_summary = primary_fit_summary,
  model_fit_plot_data = model_fit_plot_data,
  sensitivity_fit_summary = sensitivity_fit_summary,
  cure_candidate_summary = cure_candidate_summary,
  receus_candidate_map_plot_data = receus_candidate_map_plot_data,
  receus_primary_summary = receus_primary_summary,
  receus_primary_map_plot_data = receus_primary_map_plot_data,
  maller_zhou_contextual_summary = maller_zhou_contextual_summary,
  xie_status_summary = xie_status_summary,
  tail_summary = tail_summary,
  screening_overview_plot_data = screening_overview_plot_data
))

safe_write_csv_atomic(screening_summary, output_paths$screening_summary)
safe_save_plot(model_fit_plot, output_paths$model_fit_plot, width = 12, height = 7)
safe_save_plot(receus_primary_map, output_paths$receus_primary_map, width = 9, height = 6)
safe_save_plot(screening_overview_plot, output_paths$screening_overview_plot, width = 10, height = 4.8)
safe_save_plot(receus_candidate_map, output_paths$receus_candidate_map, width = 12, height = 4.8)
safe_save_pdf_atomic(
  output_paths$visual_summary_pdf,
  width = 12,
  height = 5,
  plot_fun = function() {
    print(screening_overview_plot)
    print(receus_candidate_map)
  }
)
safe_write_csv_atomic(dataset_registry, output_paths$dataset_registry)
safe_write_csv_atomic(primary_fit_summary, output_paths$primary_fit_summary)
safe_write_csv_atomic(model_fit_plot_data, output_paths$model_fit_plot_data)
safe_write_csv_atomic(sensitivity_fit_summary, output_paths$sensitivity_fit_summary)
safe_write_csv_atomic(cure_candidate_summary, output_paths$cure_candidate_summary)
safe_write_csv_atomic(receus_candidate_map_plot_data, output_paths$receus_candidate_map_plot_data)
safe_write_csv_atomic(receus_primary_summary, output_paths$receus_primary_summary)
safe_write_csv_atomic(receus_primary_map_plot_data, output_paths$receus_primary_map_plot_data)
safe_write_csv_atomic(maller_zhou_contextual_summary, output_paths$maller_zhou_contextual_summary)
safe_write_csv_atomic(xie_status_summary, output_paths$xie_status_summary)
safe_write_csv_atomic(tail_summary, output_paths$tail_summary)
safe_write_csv_atomic(screening_overview_plot_data, output_paths$screening_overview_plot_data)
safe_write_csv_atomic(block5_raw_standard_error_bundle$registry, output_paths$standard_error_registry)
safe_write_csv_atomic(block5_raw_standard_error_bundle$long, output_paths$standard_error_long)

if (isTRUE(sync_published_outputs)) {
  safe_write_csv_atomic(screening_summary, file.path(published_output_dir, "block5_screening_summary.csv"))
  safe_write_csv_atomic(cure_candidate_summary, file.path(published_output_dir, "block5_raw_cure_candidate_summary.csv"))
  safe_write_csv_atomic(primary_fit_summary, file.path(published_output_dir, "block5_raw_primary_model_fit_summary.csv"))
  safe_write_csv_atomic(model_fit_plot_data, file.path(published_output_dir, "block5_raw_model_fit_plot_data.csv"))
  safe_write_csv_atomic(sensitivity_fit_summary, file.path(published_output_dir, "block5_raw_sensitivity_model_fit_summary.csv"))
  safe_write_csv_atomic(receus_primary_summary, file.path(published_output_dir, "block5_raw_receus_primary_summary.csv"))
  safe_write_csv_atomic(receus_primary_map_plot_data, file.path(published_output_dir, "block5_raw_receus_primary_map_plot_data.csv"))
  safe_write_csv_atomic(maller_zhou_contextual_summary, file.path(published_output_dir, "block5_raw_maller_zhou_contextual_summary.csv"))
  safe_write_csv_atomic(xie_status_summary, file.path(published_output_dir, "block5_raw_xie_status_summary.csv"))
  safe_write_csv_atomic(tail_summary, file.path(published_output_dir, "block5_raw_tail_summary.csv"))
  safe_write_csv_atomic(screening_overview_plot_data, file.path(published_output_dir, "block5_raw_screening_overview_plot_data.csv"))
  safe_write_csv_atomic(receus_candidate_map_plot_data, file.path(published_output_dir, "block5_raw_receus_candidate_map_plot_data.csv"))
  safe_write_csv_atomic(block5_raw_standard_error_bundle$registry, file.path(published_output_dir, "block5_raw_standard_error_table_registry.csv"))
  safe_write_csv_atomic(block5_raw_standard_error_bundle$long, file.path(published_output_dir, "block5_raw_standard_error_long.csv"))
  safe_save_plot(screening_overview_plot, file.path(published_output_dir, "block5_screening_overview.png"), width = 10, height = 4.8)
  safe_save_plot(receus_candidate_map, file.path(published_output_dir, "block5_receus_candidate_map.png"), width = 12, height = 4.8)
  safe_save_pdf_atomic(
    file.path(published_output_dir, "block5_visual_summary.pdf"),
    width = 12,
    height = 5,
    plot_fun = function() {
      print(screening_overview_plot)
      print(receus_candidate_map)
    }
  )
}

export_manifest <- tibble(
  artifact_name = names(output_paths)[names(output_paths) != "export_manifest"],
  file_path = unlist(output_paths[names(output_paths) != "export_manifest"], use.names = FALSE)
) %>%
  mutate(
    file_name = basename(file_path),
    source_mode = dplyr::case_when(
      artifact_name == "xie_status_summary" & any(xie_status_summary$xie_method_label == "legacy_stage6_centered_bootstrap") ~ "raw_recompute_plus_legacy_xie_reuse",
      TRUE ~ "raw_recompute_primary_receus_and_contextual_maller"
    )
  )
safe_write_csv_atomic(export_manifest, output_paths$export_manifest)

message("Block 5 raw-data recompute export completed:")
for (ii in seq_len(nrow(export_manifest))) {
  message("  - ", export_manifest$file_path[[ii]])
}
