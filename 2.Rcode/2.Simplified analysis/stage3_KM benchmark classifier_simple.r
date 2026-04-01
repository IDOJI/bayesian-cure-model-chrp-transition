# 🔴 Configure: input and output paths ===============================
find_repo_root <- function(start_dir) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)

  repeat {
    stage3_reference <- file.path(current_dir, "2.Rcode", "stage3_KM benchmark classifier.r")
    if (file.exists(stage3_reference)) {
      return(current_dir)
    }

    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      break
    }
    current_dir <- parent_dir
  }

  stop(
    "Could not locate the repository root containing `2.Rcode/stage3_KM benchmark classifier.r`.",
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
  "STAGE3_SIMPLE_STAGE1_PATH",
  unset = file.path(results_root_default, "stage1_Backbone lock")
)
export_path <- Sys.getenv(
  "STAGE3_SIMPLE_EXPORT_PATH",
  unset = file.path(repo_root, "3.Results files")
)

stage1_analysis_datasets_file <- file.path(stage1_path, "stage1_analysis_datasets.rds")
stage1_horizon_registry_file <- file.path(stage1_path, "stage1_horizon_registry.csv")

yearly_estimates_file <- file.path(export_path, "stage3_simple_km_yearly_estimates.csv")
km_fit_objects_file <- file.path(export_path, "stage3_simple_km_fit_objects.rds")
survival_plot_rds_file <- file.path(export_path, "stage3_simple_km_survival_plot.rds")
risk_plot_rds_file <- file.path(export_path, "stage3_simple_km_risk_plot.rds")
survival_plot_png_file <- file.path(export_path, "stage3_simple_km_survival_plot.png")
risk_plot_png_file <- file.path(export_path, "stage3_simple_km_risk_plot.png")
refit_km_models <- identical(
  toupper(Sys.getenv("STAGE3_SIMPLE_REFIT_KM", unset = "FALSE")),
  "TRUE"
)

required_datasets <- c("PNU", "SNU", "merged")
plot_dataset_label_lookup <- c(PNU = "PNU", SNU = "SNU", merged = "Merged")
dataset_palette <- c(PNU = "#1B4332", SNU = "#2A6F97", Merged = "#C1666B")

required_packages <- c("readr", "dplyr", "tibble", "survival", "ggplot2", "scales")
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
  library(survival)
  library(ggplot2)
  library(scales)
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

validate_analysis_datasets <- function(analysis_datasets) {
  if (!is.list(analysis_datasets) || !all(required_datasets %in% names(analysis_datasets))) {
    stop("Stage 1 analysis dataset bundle must contain `PNU`, `SNU`, and `merged`.", call. = FALSE)
  }
  invisible(TRUE)
}

validate_horizon_registry <- function(horizon_registry) {
  required_horizon_cols <- c("dataset", "horizon_year", "horizon_days")
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

validate_cached_fit_objects <- function(fit_objects) {
  if (!is.list(fit_objects) || !all(required_datasets %in% names(fit_objects))) {
    stop("Cached KM fit RDS must contain `PNU`, `SNU`, and `merged` fit objects.", call. = FALSE)
  }

  for (dataset_name in required_datasets) {
    if (!inherits(fit_objects[[dataset_name]], "survfit")) {
      stop(sprintf("[%s] Cached fit object is not a `survfit` object.", dataset_name), call. = FALSE)
    }
  }

  invisible(TRUE)
}

prepare_dataset <- function(df, dataset_name) {
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

  df <- df %>%
    mutate(
      time_year = as.numeric(time_year),
      event_main = as.integer(event_main)
    )

  if (anyNA(df$time_year) || anyNA(df$event_main)) {
    stop(sprintf("[%s] Missing values detected in `time_year` or `event_main`.", dataset_name), call. = FALSE)
  }

  if (any(df$time_year < 0, na.rm = TRUE)) {
    stop(sprintf("[%s] Negative follow-up time detected.", dataset_name), call. = FALSE)
  }

  if (any(!df$event_main %in% c(0L, 1L), na.rm = TRUE)) {
    stop(sprintf("[%s] `event_main` must be coded as 0/1.", dataset_name), call. = FALSE)
  }

  df
}

get_dataset_horizons <- function(horizon_registry, dataset_name) {
  horizon_tbl <- horizon_registry %>%
    filter(dataset == dataset_name) %>%
    transmute(
      dataset = as.character(dataset),
      horizon_year = as.integer(horizon_year),
      horizon_days = as.numeric(horizon_days)
    ) %>%
    distinct() %>%
    arrange(horizon_year)

  if (nrow(horizon_tbl) == 0L) {
    stop(sprintf("[%s] No horizon rows found in Stage 1 horizon registry.", dataset_name), call. = FALSE)
  }

  horizon_tbl
}

extract_yearly_km_estimates <- function(fit, horizons_tbl, dataset_name) {
  fit_summary <- summary(fit, times = horizons_tbl$horizon_year, extend = TRUE)
  survival_values <- pmin(pmax(as.numeric(fit_summary$surv), 0), 1)

  tibble(
    dataset = dataset_name,
    horizon_year = horizons_tbl$horizon_year,
    estimated_survival_probability = survival_values,
    estimated_risk_probability = 1 - survival_values
  )
}

extract_full_km_curve <- function(fit, dataset_name) {
  curve_time <- c(0, as.numeric(fit$time))
  curve_survival <- c(1, as.numeric(fit$surv))

  tibble(
    dataset = dataset_name,
    time_year = curve_time,
    estimated_survival_probability = pmin(pmax(curve_survival, 0), 1),
    estimated_risk_probability = 1 - pmin(pmax(curve_survival, 0), 1)
  )
}

fit_single_dataset_km_object <- function(dataset_name, analysis_datasets) {
  df <- prepare_dataset(analysis_datasets[[dataset_name]], dataset_name)
  survival::survfit(
    survival::Surv(time_year, event_main) ~ 1,
    data = df
  )
}

build_single_dataset_outputs <- function(dataset_name, fit_objects, horizon_registry) {
  fit <- fit_objects[[dataset_name]]
  horizons_tbl <- get_dataset_horizons(horizon_registry, dataset_name)

  list(
    yearly_estimates = extract_yearly_km_estimates(fit, horizons_tbl, dataset_name),
    full_curve = extract_full_km_curve(fit, dataset_name)
  )
}

make_survival_plot <- function(curve_df, yearly_df) {
  ggplot(curve_df, aes(x = time_year, y = estimated_survival_probability, color = plot_dataset_label)) +
    geom_step(linewidth = 1.1) +
    geom_point(
      data = yearly_df,
      aes(x = horizon_year, y = estimated_survival_probability, color = plot_dataset_label),
      size = 2
    ) +
    scale_color_manual(values = dataset_palette) +
    scale_x_continuous(
      breaks = sort(unique(yearly_df$horizon_year)),
      expand = expansion(mult = c(0.01, 0.03))
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = "Kaplan-Meier Survival After Cohort Entry",
      x = "Years after cohort entry",
      y = "Estimated survival probability",
      color = "Dataset"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold")
    )
}

make_risk_plot <- function(curve_df, yearly_df) {
  ggplot(curve_df, aes(x = time_year, y = estimated_risk_probability, color = plot_dataset_label)) +
    geom_step(linewidth = 1.1) +
    geom_point(
      data = yearly_df,
      aes(x = horizon_year, y = estimated_risk_probability, color = plot_dataset_label),
      size = 2
    ) +
    scale_color_manual(values = dataset_palette) +
    scale_x_continuous(
      breaks = sort(unique(yearly_df$horizon_year)),
      expand = expansion(mult = c(0.01, 0.03))
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = "Kaplan-Meier Cumulative Risk After Cohort Entry",
      x = "Years after cohort entry",
      y = "Estimated risk probability",
      color = "Dataset"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold")
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
assert_file_exists(stage1_horizon_registry_file, "Stage 1 horizon registry CSV")

horizon_registry <- read_csv_checked(stage1_horizon_registry_file, "Stage 1 horizon registry CSV")
validate_horizon_registry(horizon_registry)

if (file.exists(km_fit_objects_file) && !isTRUE(refit_km_models)) {
  km_fit_objects <- readRDS(km_fit_objects_file)
  validate_cached_fit_objects(km_fit_objects)
  message("Loaded cached KM fit objects from: ", km_fit_objects_file)
} else {
  assert_file_exists(stage1_analysis_datasets_file, "Stage 1 analysis datasets RDS")
  analysis_datasets <- readRDS(stage1_analysis_datasets_file)
  validate_analysis_datasets(analysis_datasets)

  km_fit_objects <- lapply(
    required_datasets,
    fit_single_dataset_km_object,
    analysis_datasets = analysis_datasets
  )
  names(km_fit_objects) <- required_datasets
  saveRDS(km_fit_objects, km_fit_objects_file)
  message("Saved KM fit objects to: ", km_fit_objects_file)
}

# 🔴 Fit: KM models and assemble exports ===============================
dataset_results <- lapply(
  required_datasets,
  build_single_dataset_outputs,
  fit_objects = km_fit_objects,
  horizon_registry = horizon_registry
)
names(dataset_results) <- required_datasets

yearly_estimates <- bind_rows(lapply(dataset_results, `[[`, "yearly_estimates")) %>%
  mutate(
    dataset = factor(dataset, levels = required_datasets)
  ) %>%
  arrange(dataset, horizon_year) %>%
  mutate(
    dataset = as.character(dataset)
  )

full_curve <- bind_rows(lapply(dataset_results, `[[`, "full_curve")) %>%
  mutate(
    dataset = factor(dataset, levels = required_datasets)
  ) %>%
  arrange(dataset, time_year) %>%
  mutate(
    dataset = as.character(dataset)
  )

yearly_plot_data <- yearly_estimates %>%
  mutate(
    plot_dataset_label = factor(
      unname(plot_dataset_label_lookup[dataset]),
      levels = c("PNU", "SNU", "Merged")
    )
  )

full_curve_plot_data <- full_curve %>%
  mutate(
    plot_dataset_label = factor(
      unname(plot_dataset_label_lookup[dataset]),
      levels = c("PNU", "SNU", "Merged")
    )
  )

survival_plot_obj <- make_survival_plot(full_curve_plot_data, yearly_plot_data)
risk_plot_obj <- make_risk_plot(full_curve_plot_data, yearly_plot_data)

# 🔴 Export: simplified outputs only ===============================
readr::write_csv(yearly_estimates, yearly_estimates_file)
saveRDS(survival_plot_obj, survival_plot_rds_file)
saveRDS(risk_plot_obj, risk_plot_rds_file)
save_plot_png(survival_plot_obj, survival_plot_png_file)
save_plot_png(risk_plot_obj, risk_plot_png_file)

message("Saved simplified KM outputs to: ", export_path)
