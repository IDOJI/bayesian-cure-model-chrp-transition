# Configure: direct input and output paths ===============================
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

  stop(
    "Could not locate the repository root.",
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
script_dir <- file.path(repo_root, "2.Rcode", "2.Simplified Analysis")

source_data_file <- Sys.getenv(
  "STAGE3_SIMPLE_DATA_FILE",
  unset = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/0.Data/2.Preprocessed data/Preprocessed_Merged_PNUH_SNUH_Data.csv"
)
export_path <- Sys.getenv(
  "STAGE3_SIMPLE_EXPORT_PATH",
  unset = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/1.Modeling/1.KM"
)

pnu_site_label <- Sys.getenv("STAGE3_SIMPLE_PNU_SITE_LABEL", unset = "PNU")
snu_site_label <- Sys.getenv("STAGE3_SIMPLE_SNU_SITE_LABEL", unset = "SNU")
refit_km_models <- identical(
  toupper(Sys.getenv("STAGE3_SIMPLE_REFIT_KM", unset = "FALSE")),
  "TRUE"
)

yearly_estimates_file <- file.path(export_path, "stage3_simple_km_yearly_estimates.csv")
dataset_summary_file <- file.path(export_path, "stage3_simple_km_dataset_summary.csv")
plot_manifest_file <- file.path(export_path, "stage3_simple_km_plot_manifest.csv")
km_fit_objects_file <- file.path(export_path, "stage3_simple_km_fit_objects.rds")
survival_plot_rds_file <- file.path(export_path, "stage3_simple_km_survival_plot.rds")
risk_plot_rds_file <- file.path(export_path, "stage3_simple_km_risk_plot.rds")
survival_plot_png_file <- file.path(export_path, "stage3_simple_km_survival_plot.png")
risk_plot_png_file <- file.path(export_path, "stage3_simple_km_risk_plot.png")

required_datasets <- c("PNU", "SNU", "merged")
common_horizons_year <- 1:10
plot_dataset_label_lookup <- c(PNU = "PNU", SNU = "SNU", merged = "Merged")
dataset_palette <- c(PNU = "#1B4332", SNU = "#2A6F97", Merged = "#C1666B")
count_metric_palette <- c(
  "Number at risk" = "#4C78A8",
  "Cumulative transitions" = "#E07A5F"
)
plot_scenarios <- list(
  list(key = "all_cohorts", label = "PNU, SNU, and Merged", datasets = required_datasets),
  list(key = "pnu_only", label = "PNU Only", datasets = c("PNU")),
  list(key = "snu_only", label = "SNU Only", datasets = c("SNU")),
  list(key = "merged_only", label = "Merged Only", datasets = c("merged"))
)

required_packages <- c("readr", "dplyr", "tibble", "survival", "ggplot2", "scales", "patchwork")
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
  library(patchwork)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# Define: helpers ===============================
assert_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

read_csv_checked <- function(path, label) {
  assert_file_exists(path, label)
  readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
}

coerce_numeric_text <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

standardize_known_site_labels <- function(df, pnu_label, snu_label) {
  if (!"site" %in% names(df)) {
    stop("Merged input must contain column `site`.", call. = FALSE)
  }

  df %>%
    mutate(
      site = trimws(as.character(site)),
      site = dplyr::case_when(
        toupper(site) == toupper(pnu_label) ~ pnu_label,
        toupper(site) == toupper(snu_label) ~ snu_label,
        TRUE ~ site
      )
    )
}

validate_source_data <- function(df, pnu_label, snu_label) {
  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "Merged preprocessed data is missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  observed_sites <- sort(unique(trimws(as.character(df$site))))
  expected_sites <- c(pnu_label, snu_label)
  unexpected_sites <- setdiff(observed_sites, expected_sites)
  if (length(unexpected_sites) > 0L) {
    stop(
      sprintf(
        "Unexpected site labels found in merged input: %s",
        paste(unexpected_sites, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

split_site_dataset <- function(merged_df, site_label, dataset_name) {
  out <- merged_df %>%
    filter(toupper(trimws(as.character(site))) == toupper(site_label)) %>%
    mutate(site = site_label)

  if (nrow(out) == 0L) {
    stop(
      sprintf("[%s] No rows found in merged input for site label `%s`.", dataset_name, site_label),
      call. = FALSE
    )
  }

  out
}

prepare_backbone_dataset <- function(df, dataset_name, site_mode = c("single", "merged")) {
  site_mode <- match.arg(site_mode)

  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("[%s] Missing required columns: %s", dataset_name, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  out <- df %>%
    mutate(
      id = trimws(as.character(id)),
      site = trimws(as.character(site)),
      sex_num = as.integer(coerce_numeric_text(sex_num)),
      age_exact_entry = coerce_numeric_text(age_exact_entry),
      days_followup = coerce_numeric_text(days_followup),
      status_num = as.integer(coerce_numeric_text(status_num))
    )

  if (nrow(out) == 0L) {
    stop(sprintf("[%s] Dataset has zero rows after loading.", dataset_name), call. = FALSE)
  }

  if (anyNA(out[, required_cols])) {
    stop(sprintf("[%s] Missing values detected in required columns.", dataset_name), call. = FALSE)
  }

  if (any(out$id == "", na.rm = TRUE)) {
    stop(sprintf("[%s] Blank `id` values detected.", dataset_name), call. = FALSE)
  }

  if (any(out$site == "", na.rm = TRUE)) {
    stop(sprintf("[%s] Blank `site` values detected.", dataset_name), call. = FALSE)
  }

  if (any(!out$sex_num %in% c(0L, 1L), na.rm = TRUE)) {
    stop(sprintf("[%s] `sex_num` must be coded as 0/1 only.", dataset_name), call. = FALSE)
  }

  if (any(!out$status_num %in% c(0L, 1L, 2L), na.rm = TRUE)) {
    stop(sprintf("[%s] `status_num` must be coded as 0/1/2 only.", dataset_name), call. = FALSE)
  }

  if (any(out$days_followup < 0, na.rm = TRUE)) {
    stop(sprintf("[%s] Negative `days_followup` values detected.", dataset_name), call. = FALSE)
  }

  n_site_levels <- dplyr::n_distinct(out$site)
  if (site_mode == "single" && n_site_levels != 1L) {
    stop(sprintf("[%s] Single-cohort input must contain exactly one site.", dataset_name), call. = FALSE)
  }
  if (site_mode == "merged" && n_site_levels < 2L) {
    stop(sprintf("[%s] Merged input must contain at least two site levels.", dataset_name), call. = FALSE)
  }

  age_sd <- stats::sd(out$age_exact_entry)
  if (is.na(age_sd) || age_sd <= 0) {
    stop(sprintf("[%s] `age_exact_entry` must have positive standard deviation.", dataset_name), call. = FALSE)
  }

  age_mean <- mean(out$age_exact_entry)

  out <- out %>%
    mutate(
      unique_person_id = paste(site, id, sep = "_"),
      time_year = days_followup / 365.25,
      event_main = as.integer(status_num == 1L),
      right_censor_flag = as.integer(status_num == 0L),
      remission_flag = as.integer(status_num == 2L),
      censor_main = as.integer(status_num %in% c(0L, 2L)),
      sex_label = factor(if_else(sex_num == 0L, "Female", "Male"), levels = c("Female", "Male")),
      age_s = (age_exact_entry - age_mean) / (2 * age_sd)
    )

  if (nrow(out) != dplyr::n_distinct(out$unique_person_id)) {
    stop(sprintf("[%s] `site + id` is not unique within dataset.", dataset_name), call. = FALSE)
  }

  out
}

validate_analysis_datasets <- function(analysis_datasets) {
  if (!is.list(analysis_datasets) || !all(required_datasets %in% names(analysis_datasets))) {
    stop("Analysis dataset bundle must contain `PNU`, `SNU`, and `merged`.", call. = FALSE)
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

build_dataset_summary <- function(analysis_datasets) {
  bind_rows(lapply(required_datasets, function(dataset_name) {
    df <- analysis_datasets[[dataset_name]]
    tibble(
      dataset = dataset_name,
      n_subjects = nrow(df),
      n_events_transition = sum(df$status_num == 1L, na.rm = TRUE),
      n_events_remission = sum(df$status_num == 2L, na.rm = TRUE),
      n_right_censored = sum(df$status_num == 0L, na.rm = TRUE),
      max_followup_days = max(df$days_followup, na.rm = TRUE),
      max_followup_years = max(df$time_year, na.rm = TRUE)
    )
  })) %>%
    arrange(match(dataset, required_datasets))
}

derive_support_tier <- function(dataset_name, horizon_year) {
  h0 <- as.integer(horizon_year)

  if (dataset_name == "PNU") {
    if (h0 == 1L) return("primary_supported")
    if (h0 == 2L) return("sensitivity")
    return("projection")
  }

  if (dataset_name %in% c("SNU", "merged")) {
    if (h0 %in% c(1L, 2L)) return("primary_supported")
    if (h0 <= 5L) return("secondary")
    return("projection")
  }

  stop(sprintf("Unknown dataset `%s` in support-tier derivation.", dataset_name), call. = FALSE)
}

derive_horizon_evidence_class <- function(dataset_name, horizon_year) {
  h0 <- as.integer(horizon_year)

  if (dataset_name == "PNU") {
    if (h0 == 1L) return("directly_observed_data_supported")
    if (h0 == 2L) return("partly_model_dependent")
    return("mostly_extrapolated")
  }

  if (dataset_name %in% c("SNU", "merged")) {
    if (h0 %in% c(1L, 2L)) return("directly_observed_data_supported")
    if (h0 <= 5L) return("partly_model_dependent")
    return("mostly_extrapolated")
  }

  stop(sprintf("Unknown dataset `%s` in horizon-evidence derivation.", dataset_name), call. = FALSE)
}

derive_claim_restriction_flag <- function(horizon_evidence_class) {
  if (is.na(horizon_evidence_class)) {
    return(NA_character_)
  }

  if (horizon_evidence_class == "directly_observed_data_supported") {
    return("primary_claim_allowed")
  }
  if (horizon_evidence_class == "partly_model_dependent") {
    return("secondary_or_sensitivity_only")
  }
  "projection_only"
}

derive_interpretation_note <- function(dataset_name, horizon_year, support_tier, claim_restriction_flag) {
  h0 <- as.integer(horizon_year)

  if (support_tier == "primary_supported") {
    return("Primary supported horizon with comparatively direct follow-up support.")
  }
  if (dataset_name == "PNU" && h0 == 2L) {
    return("Sensitivity horizon for PNU; partly model-dependent and not for primary claims.")
  }
  if (support_tier == "secondary") {
    return("Secondary horizon with growing tail uncertainty; interpret with model-dependence explicitly acknowledged.")
  }
  if (claim_restriction_flag == "projection_only") {
    return("Projection-dominant horizon; mostly extrapolated and not eligible for primary claims.")
  }
  "Common comparison horizon retained for cross-model comparability."
}

build_horizon_registry <- function() {
  tibble::as_tibble(expand.grid(
    dataset = required_datasets,
    horizon_year = common_horizons_year,
    stringsAsFactors = FALSE
  )) %>%
    mutate(
      dataset_key = dataset,
      horizon_id = paste0("year_", horizon_year),
      horizon_days = horizon_year * 365.25,
      support_tier = mapply(derive_support_tier, dataset, horizon_year, USE.NAMES = FALSE),
      horizon_evidence_class = mapply(derive_horizon_evidence_class, dataset, horizon_year, USE.NAMES = FALSE),
      claim_restriction_flag = mapply(derive_claim_restriction_flag, horizon_evidence_class, USE.NAMES = FALSE),
      interpretation_note = mapply(
        derive_interpretation_note,
        dataset,
        horizon_year,
        support_tier,
        claim_restriction_flag,
        USE.NAMES = FALSE
      )
    ) %>%
    arrange(match(dataset, required_datasets), horizon_year)
}

build_analysis_datasets_from_source <- function(source_file, pnu_label, snu_label) {
  source_df <- read_csv_checked(source_file, "Merged preprocessed data CSV")
  source_df <- standardize_known_site_labels(source_df, pnu_label, snu_label)
  validate_source_data(source_df, pnu_label, snu_label)

  analysis_datasets <- list(
    PNU = prepare_backbone_dataset(
      split_site_dataset(source_df, pnu_label, "PNU"),
      dataset_name = "PNU",
      site_mode = "single"
    ),
    SNU = prepare_backbone_dataset(
      split_site_dataset(source_df, snu_label, "SNU"),
      dataset_name = "SNU",
      site_mode = "single"
    ),
    merged = prepare_backbone_dataset(
      source_df,
      dataset_name = "merged",
      site_mode = "merged"
    )
  )

  validate_analysis_datasets(analysis_datasets)
  analysis_datasets
}

get_dataset_horizons <- function(horizon_registry, dataset_name) {
  horizon_registry %>%
    filter(dataset == dataset_name) %>%
    transmute(
      dataset = as.character(dataset),
      horizon_year = as.integer(horizon_year),
      horizon_days = as.numeric(horizon_days)
    ) %>%
    distinct() %>%
    arrange(horizon_year)
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
  clipped_survival <- pmin(pmax(curve_survival, 0), 1)

  tibble(
    dataset = dataset_name,
    time_year = curve_time,
    estimated_survival_probability = clipped_survival,
    estimated_risk_probability = 1 - clipped_survival
  )
}

fit_single_dataset_km_object <- function(dataset_name, analysis_datasets) {
  survival::survfit(
    survival::Surv(time_year, event_main) ~ 1,
    data = analysis_datasets[[dataset_name]]
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

make_plot_data_for_scenario <- function(curve_df, yearly_df, scenario) {
  scenario_curve <- curve_df %>%
    filter(dataset %in% scenario$datasets)

  scenario_yearly <- yearly_df %>%
    filter(dataset %in% scenario$datasets)

  if (nrow(scenario_curve) == 0L || nrow(scenario_yearly) == 0L) {
    stop(sprintf("No plot data available for scenario `%s`.", scenario$key), call. = FALSE)
  }

  list(curve = scenario_curve, yearly = scenario_yearly)
}

make_pnu_followup_line_df <- function(scenario, dataset_summary) {
  if (!"PNU" %in% scenario$datasets) {
    return(tibble())
  }

  dataset_summary %>%
    filter(dataset == "PNU") %>%
    transmute(
      xintercept = max_followup_years,
      line_label = sprintf("PNU max follow-up = %.2f years", max_followup_years)
    )
}

make_km_plot <- function(curve_df, yearly_df, scenario, dataset_summary, y_col, title_prefix, y_label) {
  plot_data <- make_plot_data_for_scenario(curve_df, yearly_df, scenario)
  vertical_line_df <- make_pnu_followup_line_df(scenario, dataset_summary)

  aes_mapping <- ggplot2::aes(
    x = time_year,
    y = .data[[y_col]],
    color = plot_dataset_label
  )

  point_mapping <- ggplot2::aes(
    x = horizon_year,
    y = .data[[y_col]],
    color = plot_dataset_label
  )

  plot_object <- ggplot(plot_data$curve, aes_mapping) +
    geom_step(linewidth = 1.1) +
    geom_point(
      data = plot_data$yearly,
      mapping = point_mapping,
      size = 2
    )

  if (nrow(vertical_line_df) > 0L) {
    plot_object <- plot_object +
      geom_vline(
        data = vertical_line_df,
        aes(xintercept = xintercept),
        linewidth = 0.8,
        linetype = "22",
        color = dataset_palette[["PNU"]]
      )
  }

  subtitle_text <- if (nrow(vertical_line_df) > 0L) {
    vertical_line_df$line_label[[1L]]
  } else {
    NULL
  }

  legend_position <- if (length(scenario$datasets) > 1L) "top" else "none"
  palette_values <- dataset_palette[as.character(unique(plot_data$curve$plot_dataset_label))]

  plot_object +
    scale_color_manual(values = palette_values) +
    scale_x_continuous(
      breaks = 1:10,
      expand = expansion(mult = c(0.01, 0.03))
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(xlim = c(0, 10)) +
    labs(
      title = paste(title_prefix, "-", scenario$label),
      subtitle = subtitle_text,
      x = "Years after cohort entry (k = 1, 2, ..., 10)",
      y = y_label,
      color = "Dataset"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = legend_position,
      plot.title = element_text(face = "bold")
    )
}

build_plot_file_stem <- function(scenario_key, plot_type) {
  paste0("stage3_simple_km_", plot_type, "_plot_", scenario_key)
}

build_plot_manifest_row <- function(scenario_key, plot_type, png_file, rds_file) {
  tibble(
    scenario_key = scenario_key,
    plot_type = plot_type,
    png_file = png_file,
    rds_file = rds_file
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

build_dataset_count_plot_file_stem <- function(dataset_name, plot_type) {
  paste0(
    "stage3_simple_km_",
    plot_type,
    "_plot_",
    tolower(dataset_name),
    "_with_counts"
  )
}

cleanup_existing_with_counts_outputs <- function(output_dir) {
  old_files <- list.files(
    output_dir,
    pattern = "^stage3_simple_km_(survival|risk)_plot_(pnu|snu|merged)_with_counts\\.(png|rds)$",
    full.names = TRUE
  )

  if (length(old_files) == 0L) {
    return(invisible(character()))
  }

  removed_flag <- file.remove(old_files)
  if (any(!removed_flag)) {
    stop(
      sprintf(
        "Failed to remove existing with-counts files: %s",
        paste(basename(old_files[!removed_flag]), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(old_files[removed_flag])
}

format_half_year_label <- function(x) {
  ifelse(
    abs(x - round(x)) < 1e-08,
    as.character(as.integer(round(x))),
    formatC(x, format = "f", digits = 1)
  )
}

build_count_bar_data <- function(
    dataset_name,
    analysis_df,
    interval_year = 0.5,
    max_year = 10
) {
  time_grid <- seq(0, max_year, by = interval_year)

  tibble(
    dataset = dataset_name,
    time_year = time_grid,
    n_at_risk = vapply(
      time_grid,
      function(tt) sum(analysis_df$time_year >= tt, na.rm = TRUE),
      integer(1)
    ),
    n_transition_cumulative = vapply(
      time_grid,
      function(tt) sum(analysis_df$status_num == 1L & analysis_df$time_year <= tt, na.rm = TRUE),
      integer(1)
    )
  )
}

make_dataset_curve_component <- function(
    curve_df,
    yearly_df,
    dataset_name,
    dataset_summary,
    plot_type = c("survival", "risk")
) {
  plot_type <- match.arg(plot_type)
  dataset_label <- unname(plot_dataset_label_lookup[[dataset_name]])
  dataset_color <- unname(dataset_palette[[dataset_label]])
  curve_subset <- curve_df %>% filter(dataset == dataset_name)
  yearly_subset <- yearly_df %>% filter(dataset == dataset_name)

  subtitle_text <- NULL
  if (identical(dataset_name, "PNU")) {
    pnu_followup_years <- dataset_summary %>%
      filter(dataset == "PNU") %>%
      pull(max_followup_years)
    subtitle_text <- sprintf("PNU max follow-up = %.2f years", pnu_followup_years[[1L]])
  }

  y_col <- if (identical(plot_type, "risk")) {
    "estimated_risk_probability"
  } else {
    "estimated_survival_probability"
  }

  curve_plot <- ggplot(
    curve_subset,
    aes(x = time_year, y = .data[[y_col]])
  ) +
    geom_step(linewidth = 1.1, color = dataset_color) +
    geom_point(
      data = yearly_subset,
      aes(x = horizon_year, y = .data[[y_col]]),
      inherit.aes = FALSE,
      size = 2,
      color = dataset_color
    ) +
    scale_x_continuous(
      breaks = 1:10,
      expand = expansion(mult = c(0.01, 0.03))
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(xlim = c(0, 10)) +
    labs(
      title = if (identical(plot_type, "risk")) {
        paste("Kaplan-Meier Estimated Risk Probability with Bar Counts -", dataset_label)
      } else {
        paste("Kaplan-Meier Estimated Survival Probability with Bar Counts -", dataset_label)
      },
      subtitle = subtitle_text,
      x = NULL,
      y = if (identical(plot_type, "risk")) {
        "Estimated risk probability (1 - survival probability)"
      } else {
        "Estimated survival probability"
      }
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(5.5, 5.5, 0, 5.5)
    )

  if (identical(dataset_name, "PNU")) {
    pnu_vertical_line <- dataset_summary %>%
      filter(dataset == "PNU") %>%
      transmute(xintercept = max_followup_years)

    curve_plot <- curve_plot +
      geom_vline(
        data = pnu_vertical_line,
        aes(xintercept = xintercept),
        linewidth = 0.8,
        linetype = "22",
        color = dataset_color
      )
  }

  curve_plot
}

make_count_bar_component <- function(
    count_bar_data,
    dataset_name,
    dataset_summary
) {
  dataset_label <- unname(plot_dataset_label_lookup[[dataset_name]])
  dataset_color <- unname(dataset_palette[[dataset_label]])

  count_long_df <- bind_rows(
    count_bar_data %>%
      transmute(time_year, metric = "Number at risk", value = n_at_risk),
    count_bar_data %>%
      transmute(time_year, metric = "Cumulative transitions", value = n_transition_cumulative)
  ) %>%
    mutate(
      metric = factor(metric, levels = c("Number at risk", "Cumulative transitions"))
    )

  bar_plot <- ggplot(count_long_df, aes(x = time_year, y = value, fill = metric)) +
    geom_col(width = 0.40, alpha = 0.90) +
    geom_text(
      aes(label = value),
      vjust = -0.20,
      size = 2.4
    ) +
    facet_grid(metric ~ ., scales = "free_y", switch = "y") +
    scale_fill_manual(values = count_metric_palette, guide = "none") +
    scale_x_continuous(
      breaks = seq(0, 10, by = 0.5),
      labels = format_half_year_label,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.18))
    ) +
    coord_cartesian(xlim = c(0, 10), clip = "off") +
    labs(
      x = "Years after cohort entry (6-month intervals)",
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.placement = "outside",
      strip.background = element_rect(fill = "grey95"),
      strip.text.y.left = element_text(angle = 0, face = "bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      panel.spacing = unit(0.25, "lines"),
      plot.margin = margin(0, 5.5, 5.5, 5.5)
    )

  if (identical(dataset_name, "PNU")) {
    pnu_vertical_line <- dataset_summary %>%
      filter(dataset == "PNU") %>%
      transmute(xintercept = max_followup_years)

    bar_plot <- bar_plot +
      geom_vline(
        data = pnu_vertical_line,
        aes(xintercept = xintercept),
        linewidth = 0.6,
        linetype = "22",
        color = dataset_color
      )
  }

  bar_plot
}

make_dataset_count_plot <- function(
    curve_df,
    yearly_df,
    dataset_name,
    analysis_df,
    dataset_summary,
    plot_type = c("survival", "risk")
) {
  plot_type <- match.arg(plot_type)

  curve_component <- make_dataset_curve_component(
    curve_df = curve_df,
    yearly_df = yearly_df,
    dataset_name = dataset_name,
    dataset_summary = dataset_summary,
    plot_type = plot_type
  )
  count_bar_data <- build_count_bar_data(
    dataset_name = dataset_name,
    analysis_df = analysis_df
  )
  bar_component <- make_count_bar_component(
    count_bar_data = count_bar_data,
    dataset_name = dataset_name,
    dataset_summary = dataset_summary
  )

  curve_component / bar_component +
    patchwork::plot_layout(heights = c(2.4, 1.9))
}

save_patchwork_png <- function(plot_object, output_file, width = 12, height = 9, dpi = 300) {
  ggplot2::ggsave(
    filename = output_file,
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    units = "in"
  )
}

# Read: direct source data and build analysis inputs ===============================
message("Reading source data directly from: ", source_data_file)
analysis_datasets <- build_analysis_datasets_from_source(
  source_file = source_data_file,
  pnu_label = pnu_site_label,
  snu_label = snu_site_label
)
horizon_registry <- build_horizon_registry()
dataset_summary <- build_dataset_summary(analysis_datasets)

# Fit: KM models directly from source data or cached RDS ===============================
if (file.exists(km_fit_objects_file) && !isTRUE(refit_km_models)) {
  km_fit_objects <- readRDS(km_fit_objects_file)
  validate_cached_fit_objects(km_fit_objects)
  message("Loaded cached KM fit objects from: ", km_fit_objects_file)
} else {
  km_fit_objects <- lapply(
    required_datasets,
    fit_single_dataset_km_object,
    analysis_datasets = analysis_datasets
  )
  names(km_fit_objects) <- required_datasets
  saveRDS(km_fit_objects, km_fit_objects_file)
  message("Saved directly refit KM objects to: ", km_fit_objects_file)
}

dataset_results <- lapply(
  required_datasets,
  build_single_dataset_outputs,
  fit_objects = km_fit_objects,
  horizon_registry = horizon_registry
)
names(dataset_results) <- required_datasets

yearly_estimates <- bind_rows(lapply(dataset_results, `[[`, "yearly_estimates")) %>%
  mutate(dataset = factor(dataset, levels = required_datasets)) %>%
  arrange(dataset, horizon_year) %>%
  mutate(dataset = as.character(dataset))

full_curve <- bind_rows(lapply(dataset_results, `[[`, "full_curve")) %>%
  mutate(dataset = factor(dataset, levels = required_datasets)) %>%
  arrange(dataset, time_year) %>%
  mutate(dataset = as.character(dataset))

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

plot_manifest <- tibble()
plot_objects <- list()

for (scenario in plot_scenarios) {
  survival_plot_obj <- make_km_plot(
    curve_df = full_curve_plot_data,
    yearly_df = yearly_plot_data,
    scenario = scenario,
    dataset_summary = dataset_summary,
    y_col = "estimated_survival_probability",
    title_prefix = "Kaplan-Meier Estimated Survival Probability",
    y_label = "Estimated survival probability"
  )

  risk_plot_obj <- make_km_plot(
    curve_df = full_curve_plot_data,
    yearly_df = yearly_plot_data,
    scenario = scenario,
    dataset_summary = dataset_summary,
    y_col = "estimated_risk_probability",
    title_prefix = "Kaplan-Meier Estimated Risk Probability",
    y_label = "Estimated risk probability (1 - survival probability)"
  )

  survival_file_stem <- build_plot_file_stem(scenario$key, "survival")
  risk_file_stem <- build_plot_file_stem(scenario$key, "risk")

  survival_png_file <- file.path(export_path, paste0(survival_file_stem, ".png"))
  survival_rds_file <- file.path(export_path, paste0(survival_file_stem, ".rds"))
  risk_png_file <- file.path(export_path, paste0(risk_file_stem, ".png"))
  risk_rds_file <- file.path(export_path, paste0(risk_file_stem, ".rds"))

  saveRDS(survival_plot_obj, survival_rds_file)
  saveRDS(risk_plot_obj, risk_rds_file)
  save_plot_png(survival_plot_obj, survival_png_file)
  save_plot_png(risk_plot_obj, risk_png_file)

  plot_objects[[scenario$key]] <- list(
    survival = survival_plot_obj,
    risk = risk_plot_obj
  )

  plot_manifest <- bind_rows(
    plot_manifest,
    build_plot_manifest_row(scenario$key, "survival", survival_png_file, survival_rds_file),
    build_plot_manifest_row(scenario$key, "risk", risk_png_file, risk_rds_file)
  )

  if (identical(scenario$key, "all_cohorts")) {
    saveRDS(survival_plot_obj, survival_plot_rds_file)
    saveRDS(risk_plot_obj, risk_plot_rds_file)
    save_plot_png(survival_plot_obj, survival_plot_png_file)
    save_plot_png(risk_plot_obj, risk_plot_png_file)
  }
}

removed_with_count_files <- cleanup_existing_with_counts_outputs(export_path)
if (length(removed_with_count_files) > 0L) {
  message("Removed existing with-counts files from: ", export_path)
}

for (dataset_name in required_datasets) {
  dataset_survival_count_plot <- make_dataset_count_plot(
    curve_df = full_curve_plot_data,
    yearly_df = yearly_plot_data,
    dataset_name = dataset_name,
    analysis_df = analysis_datasets[[dataset_name]],
    dataset_summary = dataset_summary,
    plot_type = "survival"
  )

  dataset_risk_count_plot <- make_dataset_count_plot(
    curve_df = full_curve_plot_data,
    yearly_df = yearly_plot_data,
    dataset_name = dataset_name,
    analysis_df = analysis_datasets[[dataset_name]],
    dataset_summary = dataset_summary,
    plot_type = "risk"
  )

  survival_count_file_stem <- build_dataset_count_plot_file_stem(dataset_name, "survival")
  risk_count_file_stem <- build_dataset_count_plot_file_stem(dataset_name, "risk")

  survival_count_png_file <- file.path(export_path, paste0(survival_count_file_stem, ".png"))
  survival_count_rds_file <- file.path(export_path, paste0(survival_count_file_stem, ".rds"))
  risk_count_png_file <- file.path(export_path, paste0(risk_count_file_stem, ".png"))
  risk_count_rds_file <- file.path(export_path, paste0(risk_count_file_stem, ".rds"))

  saveRDS(dataset_survival_count_plot, survival_count_rds_file)
  saveRDS(dataset_risk_count_plot, risk_count_rds_file)
  save_patchwork_png(dataset_survival_count_plot, survival_count_png_file)
  save_patchwork_png(dataset_risk_count_plot, risk_count_png_file)

  plot_manifest <- bind_rows(
    plot_manifest,
    build_plot_manifest_row(
      scenario_key = paste0(tolower(dataset_name), "_with_counts"),
      plot_type = "survival_with_counts",
      png_file = survival_count_png_file,
      rds_file = survival_count_rds_file
    ),
    build_plot_manifest_row(
      scenario_key = paste0(tolower(dataset_name), "_with_counts"),
      plot_type = "risk_with_counts",
      png_file = risk_count_png_file,
      rds_file = risk_count_rds_file
    )
  )
}

# Export: simplified outputs ===============================
readr::write_csv(yearly_estimates, yearly_estimates_file)
readr::write_csv(dataset_summary, dataset_summary_file)
readr::write_csv(plot_manifest, plot_manifest_file)

message("Saved simplified KM outputs to: ", export_path)
