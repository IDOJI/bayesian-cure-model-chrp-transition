# 🔴 Configure: input and output paths ===============================
find_repo_root <- function(start_dir) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)

  repeat {
    stage5_reference <- file.path(current_dir, "2.Rcode", "stage5_Individualized no-cure comparator.r")
    if (file.exists(stage5_reference)) {
      return(current_dir)
    }

    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      break
    }
    current_dir <- parent_dir
  }

  stop(
    "Could not locate the repository root containing `2.Rcode/stage5_Individualized no-cure comparator.r`.",
    call. = FALSE
  )
}

command_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", command_args, value = TRUE)
search_start_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

repo_root <- find_repo_root(search_start_dir)
script_dir <- file.path(repo_root, "2.Rcode", "Simplified")

sys_name <- Sys.info()[["sysname"]]
results_root_default <- switch(
  sys_name,
  "Darwin" = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  "Windows" = "C:/Users/clair/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  stop("Unsupported OS: ", sys_name)
)

stage1_path <- Sys.getenv(
  "STAGE5_SIMPLE_STAGE1_PATH",
  unset = file.path(results_root_default, "stage1_Backbone lock")
)
export_path <- Sys.getenv(
  "STAGE5_SIMPLE_EXPORT_PATH",
  unset = file.path(script_dir, "exports", "stage5_individualized_nocure_lognormal_simple")
)

stage1_analysis_datasets_file <- file.path(stage1_path, "stage1_analysis_datasets.rds")
stage1_formula_registry_file <- file.path(stage1_path, "stage1_formula_registry.csv")
stage1_horizon_registry_file <- file.path(stage1_path, "stage1_horizon_registry.csv")

survival_yearly_file <- file.path(export_path, "stage5_simple_lognormal_survival_estimates.csv")
risk_yearly_file <- file.path(export_path, "stage5_simple_lognormal_risk_estimates.csv")
fitted_models_rds_file <- file.path(export_path, "stage5_simple_lognormal_fitted_models.rds")
survival_plot_rds_file <- file.path(export_path, "stage5_simple_lognormal_survival_plot.rds")
risk_plot_rds_file <- file.path(export_path, "stage5_simple_lognormal_risk_plot.rds")
survival_plot_png_file <- file.path(export_path, "stage5_simple_lognormal_survival_plot.png")
risk_plot_png_file <- file.path(export_path, "stage5_simple_lognormal_risk_plot.png")

time_origin_epsilon_year <- 1e-08
curve_step_year <- 0.05
curve_horizon_max_year <- 10

dataset_version_registry <- tibble::tibble(
  dataset_version_key = c("PNU", "SNU", "merged", "merged_site_adjusted"),
  dataset_version_label = c("PNU", "SNU", "merged", "merged (site-adjusted)"),
  source_dataset = c("PNU", "SNU", "merged", "merged"),
  formula_name = c("base", "base", "base", "site_added"),
  expected_site_branch = c("site_free", "site_free", "site_free", "site_adjusted"),
  expected_interaction_branch = rep("no_age_sex_interaction", 4L)
)

dataset_palette <- c(
  "PNU" = "#1B4332",
  "SNU" = "#2A6F97",
  "merged" = "#C1666B",
  "merged (site-adjusted)" = "#B8860B"
)

required_packages <- c("readr", "dplyr", "tibble", "ggplot2", "scales", "survival", "flexsurv")
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
  library(flexsurv)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# 🔴 Define: helpers ===============================
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

parse_flag <- function(x, default = FALSE) {
  if (length(x) == 0L || is.na(x) || trimws(x) == "") {
    return(default)
  }

  normalized_x <- tolower(trimws(as.character(x)))
  if (normalized_x %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (normalized_x %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }

  default
}

clip_prob <- function(x) {
  pmin(pmax(as.numeric(x), 0), 1)
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

validate_inputs <- function(analysis_datasets, formula_registry, horizon_registry) {
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

  required_horizon_cols <- c(
    "dataset",
    "horizon_year",
    "horizon_days",
    "support_tier",
    "horizon_evidence_class",
    "claim_restriction_flag",
    "interpretation_note"
  )
  missing_horizon_cols <- setdiff(required_horizon_cols, names(horizon_registry))
  if (length(missing_horizon_cols) > 0L) {
    stop(
      sprintf(
        "Stage 1 horizon registry is missing required columns: %s",
        paste(missing_horizon_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

prepare_dataset <- function(df, dataset_name) {
  required_cols <- c("unique_person_id", "site", "id", "sex_num", "age_exact_entry", "age_s")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("[%s] Dataset is missing required columns: %s", dataset_name, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!("time_year" %in% names(df))) {
    if (!("days_followup" %in% names(df))) {
      stop(sprintf("[%s] Dataset must contain `time_year` or `days_followup`.", dataset_name), call. = FALSE)
    }
    df$time_year <- as.numeric(df$days_followup) / 365.25
  }

  if (!("event_main" %in% names(df))) {
    if (!("status_num" %in% names(df))) {
      stop(sprintf("[%s] Dataset must contain `event_main` or `status_num`.", dataset_name), call. = FALSE)
    }
    df$event_main <- as.integer(as.numeric(df$status_num) == 1)
  }

  out <- df %>%
    mutate(
      unique_person_id = as.character(unique_person_id),
      site = factor(as.character(site)),
      id = as.character(id),
      sex_num = as.integer(sex_num),
      age_exact_entry = as.numeric(age_exact_entry),
      age_s = as.numeric(age_s),
      time_year = as.numeric(time_year),
      time_year_model = pmax(as.numeric(time_year), time_origin_epsilon_year),
      event_main = as.integer(event_main)
    )

  if (nrow(out) != dplyr::n_distinct(out$unique_person_id)) {
    stop(sprintf("[%s] `unique_person_id` must be unique.", dataset_name), call. = FALSE)
  }

  if (anyNA(out$time_year) || any(out$time_year < 0, na.rm = TRUE)) {
    stop(sprintf("[%s] Invalid `time_year` values detected.", dataset_name), call. = FALSE)
  }

  if (anyNA(out$event_main) || any(!out$event_main %in% c(0L, 1L))) {
    stop(sprintf("[%s] `event_main` must be coded as 0/1.", dataset_name), call. = FALSE)
  }

  if (anyNA(out$sex_num) || any(!out$sex_num %in% c(0L, 1L))) {
    stop(sprintf("[%s] `sex_num` must be coded as 0/1.", dataset_name), call. = FALSE)
  }

  out
}

select_model_specs <- function(formula_registry) {
  formula_tbl <- formula_registry %>%
    mutate(dataset = normalize_dataset_label(dataset))

  selected_specs <- dataset_version_registry %>%
    left_join(
      formula_tbl,
      by = c(
        "source_dataset" = "dataset",
        "formula_name" = "formula_name",
        "expected_site_branch" = "site_branch",
        "expected_interaction_branch" = "interaction_branch"
      )
    )

  if (anyNA(selected_specs$formula_id) || anyNA(selected_specs$formula_rhs)) {
    missing_keys <- selected_specs$dataset_version_key[is.na(selected_specs$formula_id) | is.na(selected_specs$formula_rhs)]
    stop(
      sprintf(
        "Could not resolve simplified Stage 5 formulas for dataset versions: %s",
        paste(unique(missing_keys), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  selected_specs %>%
    select(
      dataset_version_key,
      dataset_version_label,
      source_dataset,
      formula_id,
      formula_name,
      formula_label,
      formula_rhs,
      site_branch = expected_site_branch,
      interaction_branch = expected_interaction_branch
    )
}

get_horizon_table <- function(horizon_registry, dataset_name) {
  horizon_registry %>%
    filter(dataset == dataset_name) %>%
    transmute(
      source_dataset = as.character(dataset),
      horizon_year = as.integer(horizon_year),
      horizon_days = as.numeric(horizon_days),
      support_tier = as.character(support_tier),
      horizon_evidence_class = as.character(horizon_evidence_class),
      claim_restriction_flag = as.character(claim_restriction_flag),
      interpretation_note = as.character(interpretation_note)
    ) %>%
    distinct() %>%
    arrange(horizon_year)
}

extract_flexsurv_matrix <- function(summary_list, target_times) {
  if (!is.list(summary_list)) {
    summary_list <- list(summary_list)
  }

  t(
    vapply(
      summary_list,
      FUN.VALUE = numeric(length(target_times)),
      FUN = function(x) {
        xx <- as.data.frame(x)
        if ("time" %in% names(xx) && nrow(xx) != length(target_times)) {
          matched_idx <- vapply(
            target_times,
            FUN.VALUE = integer(1),
            FUN = function(tt) which.min(abs(xx$time - tt))
          )
          xx <- xx[matched_idx, , drop = FALSE]
        }
        if ("est" %in% names(xx)) {
          return(as.numeric(xx$est))
        }
        as.numeric(xx[[ncol(xx)]])
      }
    )
  )
}

predict_parametric_survival <- function(fit, newdata, target_times) {
  surv_list <- summary(
    fit,
    newdata = newdata,
    t = target_times,
    type = "survival",
    ci = FALSE,
    se = FALSE
  )

  extract_flexsurv_matrix(surv_list, target_times)
}

fit_single_lognormal_model <- function(spec_row, analysis_datasets) {
  dataset_name <- spec_row$source_dataset[[1]]
  formula_rhs <- spec_row$formula_rhs[[1]]

  df <- prepare_dataset(analysis_datasets[[dataset_name]], dataset_name)

  surv_formula <- stats::as.formula(
    paste0("survival::Surv(time_year_model, event_main) ~ ", formula_rhs)
  )

  fit <- flexsurv::flexsurvreg(
    formula = surv_formula,
    data = df,
    dist = "lnorm"
  )
}

build_model_exports_from_fit <- function(spec_row, fit, analysis_datasets, horizon_registry, curve_times) {
  dataset_name <- spec_row$source_dataset[[1]]
  dataset_version_key <- spec_row$dataset_version_key[[1]]
  dataset_version_label <- spec_row$dataset_version_label[[1]]
  formula_rhs <- spec_row$formula_rhs[[1]]

  df <- prepare_dataset(analysis_datasets[[dataset_name]], dataset_name)
  horizon_tbl <- get_horizon_table(horizon_registry, dataset_name)

  yearly_survival_mat <- predict_parametric_survival(
    fit = fit,
    newdata = df,
    target_times = horizon_tbl$horizon_year
  )
  yearly_survival <- clip_prob(colMeans(yearly_survival_mat, na.rm = TRUE))

  curve_survival_mat <- predict_parametric_survival(
    fit = fit,
    newdata = df,
    target_times = curve_times
  )
  curve_survival <- clip_prob(colMeans(curve_survival_mat, na.rm = TRUE))

  yearly_summary <- horizon_tbl %>%
    mutate(
      dataset_version_key = dataset_version_key,
      dataset_version_label = dataset_version_label,
      source_dataset = dataset_name,
      model_family = "lognormal",
      formula_id = spec_row$formula_id[[1]],
      formula_name = spec_row$formula_name[[1]],
      formula_label = spec_row$formula_label[[1]],
      formula_rhs = formula_rhs,
      site_branch = spec_row$site_branch[[1]],
      interaction_branch = spec_row$interaction_branch[[1]],
      estimated_survival_probability = yearly_survival,
      estimated_risk_probability = 1 - yearly_survival,
      n_subjects_averaged = nrow(df),
      n_transition_events_total = sum(df$event_main == 1L, na.rm = TRUE),
      max_followup_years = max(df$time_year, na.rm = TRUE)
    ) %>%
    select(
      dataset_version_key,
      dataset_version_label,
      source_dataset,
      model_family,
      formula_id,
      formula_name,
      formula_label,
      formula_rhs,
      site_branch,
      interaction_branch,
      horizon_year,
      horizon_days,
      estimated_survival_probability,
      estimated_risk_probability,
      n_subjects_averaged,
      n_transition_events_total,
      max_followup_years,
      support_tier,
      horizon_evidence_class,
      claim_restriction_flag,
      interpretation_note
    )

  curve_df <- tibble(
    dataset_version_key = dataset_version_key,
    dataset_version_label = dataset_version_label,
    source_dataset = dataset_name,
    model_family = "lognormal",
    formula_id = spec_row$formula_id[[1]],
    formula_name = spec_row$formula_name[[1]],
    formula_label = spec_row$formula_label[[1]],
    formula_rhs = formula_rhs,
    site_branch = spec_row$site_branch[[1]],
    interaction_branch = spec_row$interaction_branch[[1]],
    time_year = c(0, curve_times),
    estimated_survival_probability = c(1, curve_survival),
    estimated_risk_probability = c(0, 1 - curve_survival),
    n_subjects_averaged = nrow(df),
    n_transition_events_total = sum(df$event_main == 1L, na.rm = TRUE)
  )

  list(
    yearly_summary = yearly_summary,
    curve_df = curve_df
  )
}

validate_fitted_model_cache <- function(fitted_models, model_specs) {
  if (!is.list(fitted_models) || length(fitted_models) != nrow(model_specs)) {
    return(FALSE)
  }

  required_names <- model_specs$dataset_version_key
  if (!identical(sort(names(fitted_models)), sort(required_names))) {
    return(FALSE)
  }

  ordered_models <- fitted_models[required_names]
  all(vapply(ordered_models, inherits, logical(1), what = "flexsurvreg"))
}

make_survival_plot <- function(curve_df, yearly_survival_df) {
  ggplot(curve_df, aes(x = time_year, y = estimated_survival_probability, color = dataset_version_label)) +
    geom_line(linewidth = 1.1) +
    geom_point(
      data = yearly_survival_df,
      aes(x = horizon_year, y = estimated_survival_probability, color = dataset_version_label),
      size = 2
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
      title = "Simplified Stage 5 Log-normal Survival After Cohort Entry",
      subtitle = "Population-averaged predictions from subject-level log-normal no-cure models",
      x = "Years after cohort entry",
      y = "Estimated survival probability",
      color = "Dataset version"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
}

make_risk_plot <- function(curve_df, yearly_risk_df) {
  ggplot(curve_df, aes(x = time_year, y = estimated_risk_probability, color = dataset_version_label)) +
    geom_line(linewidth = 1.1) +
    geom_point(
      data = yearly_risk_df,
      aes(x = horizon_year, y = estimated_risk_probability, color = dataset_version_label),
      size = 2
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
      title = "Simplified Stage 5 Log-normal Cumulative Risk After Cohort Entry",
      subtitle = "Population-averaged predictions from subject-level log-normal no-cure models",
      x = "Years after cohort entry",
      y = "Estimated risk probability",
      color = "Dataset version"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
}

save_plot_png <- function(plot_object, output_file, width = 10, height = 6, dpi = 300) {
  ggplot2::ggsave(
    filename = output_file,
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    units = "in"
  )
}

# 🔴 Read: Stage 1 inputs ===============================
assert_file_exists(stage1_analysis_datasets_file, "Stage 1 analysis datasets RDS")
assert_file_exists(stage1_formula_registry_file, "Stage 1 formula registry CSV")
assert_file_exists(stage1_horizon_registry_file, "Stage 1 horizon registry CSV")

analysis_datasets <- readRDS(stage1_analysis_datasets_file)
formula_registry <- read_csv_checked(stage1_formula_registry_file, "Stage 1 formula registry CSV")
horizon_registry <- read_csv_checked(stage1_horizon_registry_file, "Stage 1 horizon registry CSV")

validate_inputs(analysis_datasets, formula_registry, horizon_registry)
model_specs <- select_model_specs(formula_registry)

curve_times <- seq(curve_step_year, curve_horizon_max_year, by = curve_step_year)
reuse_existing_fitted_models <- parse_flag(
  Sys.getenv("STAGE5_SIMPLE_REUSE_FITTED_MODELS", unset = "true"),
  default = TRUE
)
force_refit_fitted_models <- parse_flag(
  Sys.getenv("STAGE5_SIMPLE_FORCE_REFIT", unset = "false"),
  default = FALSE
)

# 🔴 Fit: simplified Stage 5 log-normal models ===============================
fitted_models <- NULL

if (isTRUE(reuse_existing_fitted_models) && !isTRUE(force_refit_fitted_models) && file.exists(fitted_models_rds_file)) {
  fitted_models <- tryCatch(readRDS(fitted_models_rds_file), error = function(e) NULL)
  if (validate_fitted_model_cache(fitted_models, model_specs)) {
    fitted_models <- fitted_models[model_specs$dataset_version_key]
    message("Loaded cached fitted models from: ", fitted_models_rds_file)
  } else {
    fitted_models <- NULL
    message("Cached fitted model RDS was unavailable or invalid; refitting models.")
  }
}

if (is.null(fitted_models)) {
  fitted_models <- lapply(
    seq_len(nrow(model_specs)),
    function(i) {
      fit_single_lognormal_model(
        spec_row = model_specs[i, , drop = FALSE],
        analysis_datasets = analysis_datasets
      )
    }
  )
  names(fitted_models) <- model_specs$dataset_version_key
  saveRDS(fitted_models, fitted_models_rds_file)
  message("Saved fitted models to: ", fitted_models_rds_file)
}

fit_results <- lapply(
  seq_len(nrow(model_specs)),
  function(i) {
    build_model_exports_from_fit(
      spec_row = model_specs[i, , drop = FALSE],
      fit = fitted_models[[model_specs$dataset_version_key[[i]]]],
      analysis_datasets = analysis_datasets,
      horizon_registry = horizon_registry,
      curve_times = curve_times
    )
  }
)
names(fit_results) <- model_specs$dataset_version_key

yearly_summary <- bind_rows(lapply(fit_results, `[[`, "yearly_summary")) %>%
  mutate(
    dataset_version_label = factor(
      dataset_version_label,
      levels = dataset_version_registry$dataset_version_label
    )
  ) %>%
  arrange(match(dataset_version_key, dataset_version_registry$dataset_version_key), horizon_year) %>%
  mutate(dataset_version_label = as.character(dataset_version_label))

curve_df <- bind_rows(lapply(fit_results, `[[`, "curve_df")) %>%
  mutate(
    dataset_version_label = factor(
      dataset_version_label,
      levels = dataset_version_registry$dataset_version_label
    )
  ) %>%
  arrange(match(dataset_version_key, dataset_version_registry$dataset_version_key), time_year) %>%
  mutate(dataset_version_label = as.character(dataset_version_label))

csv_export_columns <- c(
  "dataset_version_key",
  "model_family",
  "formula_rhs",
  "horizon_year",
  "estimated_survival_probability",
  "estimated_risk_probability"
)

survival_plot_input <- yearly_summary
risk_plot_input <- yearly_summary

survival_yearly <- yearly_summary %>%
  select(all_of(csv_export_columns))

risk_yearly <- yearly_summary %>%
  select(all_of(csv_export_columns))

survival_plot_obj <- make_survival_plot(curve_df, survival_plot_input)
risk_plot_obj <- make_risk_plot(curve_df, risk_plot_input)

# 🔴 Export: simplified outputs only ===============================
readr::write_csv(survival_yearly, survival_yearly_file)
readr::write_csv(risk_yearly, risk_yearly_file)
saveRDS(survival_plot_obj, survival_plot_rds_file)
saveRDS(risk_plot_obj, risk_plot_rds_file)
save_plot_png(survival_plot_obj, survival_plot_png_file)
save_plot_png(risk_plot_obj, risk_plot_png_file)

message("Saved simplified Stage 5 log-normal outputs to: ", export_path)
