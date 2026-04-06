# Configure Paths ---------------------------------------------------------
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
block1_root_default <- "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/1.Block1"

km_input_path <- Sys.getenv(
  "STAGE3_SIMPLE_EXPORT_PATH",
  unset = file.path(block1_root_default, "1.KM")
)
mle_input_path <- Sys.getenv(
  "STAGE7_SIMPLE_EXPORT_PATH",
  unset = file.path(block1_root_default, "3.MLE Mixture Cure")
)
export_path <- Sys.getenv(
  "RESULTS_INTEGRATION_EXPORT_PATH",
  unset = file.path(block1_root_default, "5.Results Integration", "km_vs_mle_cure_integrated_curves")
)

km_dataset_summary_file <- file.path(km_input_path, "simple_km_dataset_summary.csv")
km_survival_plot_file <- file.path(km_input_path, "simple_km_survival_plot.rds")
km_risk_plot_file <- file.path(km_input_path, "simple_km_risk_plot.rds")
mle_horizon_summary_file <- file.path(mle_input_path, "mle_mixture_cure_lognormal_horizon_summary.csv")
mle_plot_source_file <- file.path(mle_input_path, "mle_mixture_cure_lognormal_plot_source.csv")
mle_detail_annotation_file <- file.path(mle_input_path, "mle_mixture_cure_lognormal_detail_annotation_table.csv")

horizon_years <- 1:10
detail_tick_step_year <- 0.5
curve_horizon_max_year <- 10

dataset_order <- c("PNU", "SNU", "merged", "merged_site_adjusted")
dataset_label_lookup <- c(
  "PNU" = "PNU",
  "SNU" = "SNU",
  "merged" = "Merged",
  "merged_site_adjusted" = "Merged (site-adjusted)"
)
model_order <- c("KM", "MLE cure")
model_palette <- c(
  "KM" = "#1B4332",
  "MLE cure" = "#C1666B"
)
count_fill_palette <- c(
  "Number at risk" = "#8FA3BF",
  "Cumulative transitions" = "#D98C6C"
)
estimate_fill_palette <- c(
  "KM" = "#EAF4EC",
  "MLE cure" = "#FBE7EA"
)

plot_width_in <- 14
plot_height_in <- 10
plot_dpi <- 320

# Load Packages -----------------------------------------------------------
required_packages <- c("readr", "dplyr", "tibble", "ggplot2", "scales", "patchwork")
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
  library(patchwork)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# Define Helpers ----------------------------------------------------------
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

read_rds_checked <- function(path, label) {
  assert_file_exists(path, label)
  readRDS(path)
}

clip_probability <- function(x) {
  pmin(pmax(as.numeric(x), 0), 1)
}

cleanup_existing_integrated_outputs <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    return(invisible(character()))
  }

  stale_files <- list.files(output_dir, pattern = ".*\\.(csv|rds|png)$", full.names = TRUE)
  if (length(stale_files) == 0L) {
    return(invisible(character()))
  }

  removed_flag <- file.remove(stale_files)
  if (any(!removed_flag)) {
    stop(
      sprintf(
        "Failed to remove stale integrated KM-vs-cure outputs: %s",
        paste(basename(stale_files[!removed_flag]), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(stale_files[removed_flag])
}

normalize_dataset_name <- function(x) {
  x <- trimws(as.character(x))

  dplyr::case_when(
    toupper(x) == "PNU" ~ "PNU",
    toupper(x) == "SNU" ~ "SNU",
    tolower(x) %in% c("merged", "merged_no_site") ~ "merged",
    tolower(x) == "merged_site_adjusted" ~ "merged_site_adjusted",
    TRUE ~ x
  )
}

format_half_year_label <- function(x) {
  ifelse(
    abs(x - round(x)) < 1e-8,
    as.character(as.integer(round(x))),
    formatC(x, format = "f", digits = 1)
  )
}

extract_point_layer_df <- function(plot_object, dataset_col, time_col, value_col, label) {
  candidate_layers <- lapply(plot_object$layers, function(layer) {
    if (is.null(layer$data) || inherits(layer$data, "waiver")) {
      return(NULL)
    }
    tibble::as_tibble(layer$data)
  })
  candidate_layers <- Filter(Negate(is.null), candidate_layers)
  candidate_layers <- Filter(
    function(df) all(c(dataset_col, time_col, value_col) %in% names(df)),
    candidate_layers
  )

  if (length(candidate_layers) == 0L) {
    stop(label, " does not contain a compatible yearly point layer.", call. = FALSE)
  }

  horizon_layer <- which(vapply(candidate_layers, function(df) {
    values <- suppressWarnings(as.numeric(df[[time_col]]))
    values <- values[is.finite(values)]
    identical(sort(unique(as.integer(round(values)))), as.integer(horizon_years))
  }, logical(1)))

  if (length(horizon_layer) > 0L) {
    return(candidate_layers[[horizon_layer[[1L]]]])
  }

  candidate_layers[[which.min(vapply(candidate_layers, nrow, integer(1)))]]
}

extract_km_curve_data <- function(plot_object, value_col) {
  if (!inherits(plot_object, "ggplot")) {
    stop("KM plot object must be a ggplot object.", call. = FALSE)
  }

  curve_df <- tibble::as_tibble(plot_object$data)
  yearly_df <- extract_point_layer_df(
    plot_object = plot_object,
    dataset_col = "dataset",
    time_col = "horizon_year",
    value_col = value_col,
    label = "KM plot object"
  )

  list(
    curve = curve_df %>%
      transmute(
        dataset = normalize_dataset_name(dataset),
        time_year = as.numeric(time_year),
        estimate_value = clip_probability(.data[[value_col]])
      ) %>%
      filter(dataset %in% dataset_order, time_year >= 0, time_year <= curve_horizon_max_year),
    yearly = yearly_df %>%
      transmute(
        dataset = normalize_dataset_name(dataset),
        horizon_year = as.integer(round(as.numeric(horizon_year))),
        estimate_value = clip_probability(.data[[value_col]])
      ) %>%
      filter(dataset %in% dataset_order, horizon_year %in% horizon_years)
  )
}

prepare_mle_curve_data <- function(plot_source_df, horizon_summary_df, value_col) {
  curve_df <- plot_source_df %>%
    transmute(
      dataset = normalize_dataset_name(dataset),
      time_year = as.numeric(time_horizon_year),
      estimate_value = clip_probability(.data[[value_col]])
    ) %>%
    filter(dataset %in% dataset_order, time_year >= 0, time_year <= curve_horizon_max_year) %>%
    arrange(match(dataset, dataset_order), time_year)

  yearly_df <- horizon_summary_df %>%
    transmute(
      dataset = normalize_dataset_name(dataset),
      horizon_year = as.integer(round(as.numeric(time_horizon_year))),
      estimate_value = clip_probability(.data[[value_col]])
    ) %>%
    filter(dataset %in% dataset_order, horizon_year %in% horizon_years) %>%
    distinct(dataset, horizon_year, .keep_all = TRUE) %>%
    arrange(match(dataset, dataset_order), horizon_year)

  list(curve = curve_df, yearly = yearly_df)
}

build_integrated_plot_source_table <- function(
  survival_curve_df,
  survival_yearly_df,
  risk_curve_df,
  risk_yearly_df,
  annotation_df
) {
  bind_rows(
    survival_curve_df %>%
      transmute(
        metric = "survival",
        component_type = "curve",
        dataset = as.character(dataset),
        dataset_label = unname(dataset_label_lookup[dataset]),
        model = as.character(model),
        time_year = as.numeric(time_year),
        horizon_year = NA_integer_,
        estimate_value = as.numeric(estimate_value),
        n_at_risk = NA_integer_,
        n_transition_cumulative = NA_integer_
      ),
    survival_yearly_df %>%
      transmute(
        metric = "survival",
        component_type = "yearly_point",
        dataset = as.character(dataset),
        dataset_label = unname(dataset_label_lookup[dataset]),
        model = as.character(model),
        time_year = as.numeric(horizon_year),
        horizon_year = as.integer(horizon_year),
        estimate_value = as.numeric(estimate_value),
        n_at_risk = NA_integer_,
        n_transition_cumulative = NA_integer_
      ),
    risk_curve_df %>%
      transmute(
        metric = "risk",
        component_type = "curve",
        dataset = as.character(dataset),
        dataset_label = unname(dataset_label_lookup[dataset]),
        model = as.character(model),
        time_year = as.numeric(time_year),
        horizon_year = NA_integer_,
        estimate_value = as.numeric(estimate_value),
        n_at_risk = NA_integer_,
        n_transition_cumulative = NA_integer_
      ),
    risk_yearly_df %>%
      transmute(
        metric = "risk",
        component_type = "yearly_point",
        dataset = as.character(dataset),
        dataset_label = unname(dataset_label_lookup[dataset]),
        model = as.character(model),
        time_year = as.numeric(horizon_year),
        horizon_year = as.integer(horizon_year),
        estimate_value = as.numeric(estimate_value),
        n_at_risk = NA_integer_,
        n_transition_cumulative = NA_integer_
      ),
    annotation_df %>%
      transmute(
        metric = "followup",
        component_type = "count_panel",
        dataset = as.character(dataset),
        dataset_label = unname(dataset_label_lookup[dataset]),
        model = NA_character_,
        time_year = as.numeric(time_horizon_year),
        horizon_year = NA_integer_,
        estimate_value = NA_real_,
        n_at_risk = as.integer(n_at_risk),
        n_transition_cumulative = as.integer(n_transition_cumulative)
      )
  ) %>%
    mutate(
      .metric_order = match(metric, c("survival", "risk", "followup")),
      .dataset_order = match(dataset, dataset_order),
      .component_order = match(component_type, c("curve", "yearly_point", "count_panel")),
      .model_order = ifelse(is.na(model), length(model_order) + 1L, match(model, model_order)),
      .time_order = dplyr::coalesce(time_year, as.numeric(horizon_year))
    ) %>%
    arrange(.metric_order, .dataset_order, .component_order, .model_order, .time_order) %>%
    select(-starts_with("."))
}

build_count_panels <- function(annotation_df, dataset_name, dataset_summary) {
  x_limits <- c(-detail_tick_step_year * 0.4, curve_horizon_max_year + detail_tick_step_year * 0.4)
  pnu_followup_year <- dataset_summary %>%
    filter(dataset == "PNU") %>%
    pull(max_followup_years)
  pnu_followup_year <- if (length(pnu_followup_year) > 0L) pnu_followup_year[[1L]] else NA_real_

  at_risk_plot <- ggplot(annotation_df, aes(x = time_horizon_year, y = n_at_risk)) +
    geom_col(width = detail_tick_step_year * 0.72, fill = count_fill_palette[["Number at risk"]]) +
    geom_text(aes(label = n_at_risk), vjust = -0.25, size = 2.8) +
    scale_x_continuous(
      breaks = seq(0, curve_horizon_max_year, by = detail_tick_step_year),
      labels = format_half_year_label,
      limits = x_limits,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(x = NULL, y = "n at risk") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 5.5, b = 0, l = 5.5)
    )

  transition_plot <- ggplot(annotation_df, aes(x = time_horizon_year, y = n_transition_cumulative)) +
    geom_col(width = detail_tick_step_year * 0.72, fill = count_fill_palette[["Cumulative transitions"]]) +
    geom_text(aes(label = n_transition_cumulative), vjust = -0.25, size = 2.8) +
    scale_x_continuous(
      breaks = seq(0, curve_horizon_max_year, by = detail_tick_step_year),
      labels = format_half_year_label,
      limits = x_limits,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(
      x = "Years after cohort entry",
      y = "cumulative\ntransitions"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(size = 7),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5)
    )

  if (identical(dataset_name, "PNU") && is.finite(pnu_followup_year)) {
    at_risk_plot <- at_risk_plot +
      geom_vline(
        xintercept = pnu_followup_year,
        linetype = "dashed",
        linewidth = 0.6,
        color = model_palette[["KM"]]
      )

    transition_plot <- transition_plot +
      geom_vline(
        xintercept = pnu_followup_year,
        linetype = "dashed",
        linewidth = 0.6,
        color = model_palette[["KM"]]
      )
  }

  list(at_risk_plot = at_risk_plot, transition_plot = transition_plot)
}

build_estimate_panel <- function(estimate_df, metric_label) {
  plot_df <- estimate_df %>%
    mutate(
      model = factor(model, levels = model_order),
      horizon_year = factor(horizon_year, levels = horizon_years)
    )

  ggplot(plot_df, aes(x = horizon_year, y = model)) +
    geom_tile(aes(fill = model), color = "white", linewidth = 0.6, height = 0.92) +
    geom_text(aes(label = scales::percent(estimate_value, accuracy = 0.1)), size = 3.2, fontface = "bold") +
    scale_fill_manual(values = estimate_fill_palette, guide = "none") +
    labs(
      x = NULL,
      y = NULL,
      title = paste0("Yearly estimated ", metric_label)
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 10),
      plot.margin = margin(t = 2, r = 5.5, b = 2, l = 5.5)
    )
}

build_curve_panel <- function(curve_df, yearly_df, dataset_name, dataset_summary, metric) {
  metric <- match.arg(metric, choices = c("survival", "risk"))
  dataset_label <- unname(dataset_label_lookup[[dataset_name]])
  pnu_followup_year <- dataset_summary %>%
    filter(dataset == "PNU") %>%
    pull(max_followup_years)
  pnu_followup_year <- if (length(pnu_followup_year) > 0L) pnu_followup_year[[1L]] else NA_real_

  y_label <- if (identical(metric, "survival")) {
    "Estimated survival probability"
  } else {
    "Estimated cumulative risk probability"
  }

  title_text <- if (identical(metric, "survival")) {
    paste0("KM vs MLE cure estimated survival curves - ", dataset_label)
  } else {
    paste0("KM vs MLE cure estimated cumulative risk curves - ", dataset_label)
  }

  subtitle_text <- "KM curve extracted from benchmark RDS; MLE cure curve extracted from frequentist log-normal mixture cure exports"
  if (identical(dataset_name, "PNU") && is.finite(pnu_followup_year)) {
    subtitle_text <- paste0(
      subtitle_text,
      sprintf(" | Dashed line = PNU max observed follow-up (%.2f years)", pnu_followup_year)
    )
  }

  plot_df <- curve_df %>%
    filter(dataset == dataset_name) %>%
    mutate(model = factor(model, levels = model_order))

  yearly_points <- yearly_df %>%
    filter(dataset == dataset_name) %>%
    mutate(model = factor(model, levels = model_order))

  curve_plot <- ggplot(plot_df, aes(x = time_year, y = estimate_value, color = model, linetype = model)) +
    geom_line(linewidth = 1.05) +
    geom_point(
      data = yearly_points,
      aes(x = horizon_year, y = estimate_value, color = model),
      inherit.aes = FALSE,
      size = 2.0
    ) +
    scale_color_manual(values = model_palette) +
    scale_linetype_manual(values = c("KM" = "solid", "MLE cure" = "22")) +
    scale_x_continuous(
      breaks = seq(0, curve_horizon_max_year, by = detail_tick_step_year),
      labels = format_half_year_label,
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
      subtitle = subtitle_text,
      x = NULL,
      y = y_label,
      color = "Model",
      linetype = "Model"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5)
    )

  if (identical(dataset_name, "PNU") && is.finite(pnu_followup_year)) {
    curve_plot <- curve_plot +
      geom_vline(
        xintercept = pnu_followup_year,
        linetype = "dashed",
        linewidth = 0.7,
        color = model_palette[["KM"]]
      )
  }

  curve_plot
}

build_combined_plot <- function(
  dataset_name,
  metric,
  curve_df,
  yearly_df,
  annotation_df,
  dataset_summary
) {
  metric <- match.arg(metric, choices = c("survival", "risk"))

  curve_panel <- build_curve_panel(
    curve_df = curve_df,
    yearly_df = yearly_df,
    dataset_name = dataset_name,
    dataset_summary = dataset_summary,
    metric = metric
  )

  estimate_panel <- build_estimate_panel(
    estimate_df = yearly_df %>% filter(dataset == dataset_name),
    metric_label = if (identical(metric, "survival")) "survival probabilities" else "cumulative risk probabilities"
  )

  count_panels <- build_count_panels(
    annotation_df = annotation_df,
    dataset_name = dataset_name,
    dataset_summary = dataset_summary
  )

  curve_panel / estimate_panel / count_panels$at_risk_plot / count_panels$transition_plot +
    patchwork::plot_layout(heights = c(3.6, 1.2, 1.4, 1.4))
}

save_patchwork_png <- function(plot_object, output_file) {
  ggplot2::ggsave(
    filename = output_file,
    plot = plot_object,
    width = plot_width_in,
    height = plot_height_in,
    dpi = plot_dpi,
    units = "in"
  )
}

# Read Inputs -------------------------------------------------------------
km_dataset_summary <- read_csv_checked(km_dataset_summary_file, "KM dataset summary")
km_survival_plot <- read_rds_checked(km_survival_plot_file, "KM survival plot RDS")
km_risk_plot <- read_rds_checked(km_risk_plot_file, "KM risk plot RDS")
mle_horizon_summary <- read_csv_checked(mle_horizon_summary_file, "MLE horizon summary")
mle_plot_source <- read_csv_checked(mle_plot_source_file, "MLE plot source")
mle_detail_annotation <- read_csv_checked(mle_detail_annotation_file, "MLE detail annotation table")

# Prepare Sources ---------------------------------------------------------
km_survival_data <- extract_km_curve_data(km_survival_plot, "estimated_survival_probability")
km_risk_data <- extract_km_curve_data(km_risk_plot, "estimated_risk_probability")
mle_survival_data <- prepare_mle_curve_data(mle_plot_source, mle_horizon_summary, "overall_survival_prob")
mle_risk_data <- prepare_mle_curve_data(mle_plot_source, mle_horizon_summary, "overall_risk_prob")

km_survival_data$curve <- bind_rows(
  km_survival_data$curve,
  km_survival_data$curve %>% filter(dataset == "merged") %>% mutate(dataset = "merged_site_adjusted")
)
km_survival_data$yearly <- bind_rows(
  km_survival_data$yearly,
  km_survival_data$yearly %>% filter(dataset == "merged") %>% mutate(dataset = "merged_site_adjusted")
)
km_risk_data$curve <- bind_rows(
  km_risk_data$curve,
  km_risk_data$curve %>% filter(dataset == "merged") %>% mutate(dataset = "merged_site_adjusted")
)
km_risk_data$yearly <- bind_rows(
  km_risk_data$yearly,
  km_risk_data$yearly %>% filter(dataset == "merged") %>% mutate(dataset = "merged_site_adjusted")
)

survival_curve_df <- bind_rows(
  km_survival_data$curve %>% mutate(model = "KM"),
  mle_survival_data$curve %>% mutate(model = "MLE cure")
) %>%
  mutate(model = factor(model, levels = model_order))

survival_yearly_df <- bind_rows(
  km_survival_data$yearly %>% mutate(model = "KM"),
  mle_survival_data$yearly %>% mutate(model = "MLE cure")
) %>%
  mutate(model = factor(model, levels = model_order))

risk_curve_df <- bind_rows(
  km_risk_data$curve %>% mutate(model = "KM"),
  mle_risk_data$curve %>% mutate(model = "MLE cure")
) %>%
  mutate(model = factor(model, levels = model_order))

risk_yearly_df <- bind_rows(
  km_risk_data$yearly %>% mutate(model = "KM"),
  mle_risk_data$yearly %>% mutate(model = "MLE cure")
) %>%
  mutate(model = factor(model, levels = model_order))

annotation_df <- mle_detail_annotation %>%
  transmute(
    dataset = normalize_dataset_name(source_dataset),
    time_horizon_year = as.numeric(time_horizon_year),
    n_at_risk = as.integer(n_at_risk),
    n_transition_cumulative = as.integer(n_transition_cumulative)
  ) %>%
  filter(dataset %in% dataset_order)

annotation_df <- bind_rows(
  annotation_df,
  annotation_df %>% filter(dataset == "merged") %>% mutate(dataset = "merged_site_adjusted")
)

integrated_estimate_table <- bind_rows(
  survival_yearly_df %>% mutate(metric = "survival"),
  risk_yearly_df %>% mutate(metric = "risk")
) %>%
  mutate(
    dataset_label = unname(dataset_label_lookup[dataset]),
    metric = factor(metric, levels = c("survival", "risk")),
    model = factor(model, levels = model_order)
  ) %>%
  arrange(metric, match(dataset, dataset_order), horizon_year, model)

plot_registry <- list()
plot_manifest <- tibble()

cleanup_existing_integrated_outputs(export_path)
plot_source_table <- build_integrated_plot_source_table(
  survival_curve_df = survival_curve_df,
  survival_yearly_df = survival_yearly_df,
  risk_curve_df = risk_curve_df,
  risk_yearly_df = risk_yearly_df,
  annotation_df = annotation_df
)

# Build Plots -------------------------------------------------------------
for (dataset_name in dataset_order) {
  dataset_annotation_df <- annotation_df %>%
    filter(dataset == dataset_name)

  if (nrow(dataset_annotation_df) == 0L) {
    stop(sprintf("No annotation rows found for dataset `%s`.", dataset_name), call. = FALSE)
  }

  survival_plot <- build_combined_plot(
    dataset_name = dataset_name,
    metric = "survival",
    curve_df = survival_curve_df,
    yearly_df = survival_yearly_df,
    annotation_df = dataset_annotation_df,
    dataset_summary = km_dataset_summary
  )

  risk_plot <- build_combined_plot(
    dataset_name = dataset_name,
    metric = "risk",
    curve_df = risk_curve_df,
    yearly_df = risk_yearly_df,
    annotation_df = dataset_annotation_df,
    dataset_summary = km_dataset_summary
  )

  dataset_slug <- tolower(dataset_name)
  survival_png_file <- file.path(export_path, paste0("km_vs_mle_cure_survival_curve_with_counts__", dataset_slug, ".png"))
  survival_rds_file <- file.path(export_path, paste0("km_vs_mle_cure_survival_curve_with_counts__", dataset_slug, ".rds"))
  risk_png_file <- file.path(export_path, paste0("km_vs_mle_cure_risk_curve_with_counts__", dataset_slug, ".png"))
  risk_rds_file <- file.path(export_path, paste0("km_vs_mle_cure_risk_curve_with_counts__", dataset_slug, ".rds"))

  saveRDS(survival_plot, survival_rds_file)
  saveRDS(risk_plot, risk_rds_file)
  save_patchwork_png(survival_plot, survival_png_file)
  save_patchwork_png(risk_plot, risk_png_file)

  plot_registry[[paste0(dataset_name, "__survival")]] <- survival_plot
  plot_registry[[paste0(dataset_name, "__risk")]] <- risk_plot

  plot_manifest <- bind_rows(
    plot_manifest,
    tibble(
      dataset = dataset_name,
      metric = "survival",
      png_file = survival_png_file,
      rds_file = survival_rds_file
    ),
    tibble(
      dataset = dataset_name,
      metric = "risk",
      png_file = risk_png_file,
      rds_file = risk_rds_file
    )
  )
}

# Write Outputs -----------------------------------------------------------
readr::write_csv(
  plot_source_table,
  file.path(export_path, "km_vs_mle_cure_plot_source.csv")
)
readr::write_csv(
  integrated_estimate_table,
  file.path(export_path, "km_vs_mle_cure_integrated_yearly_estimates.csv")
)
readr::write_csv(
  annotation_df,
  file.path(export_path, "km_vs_mle_cure_count_annotation_source.csv")
)
readr::write_csv(
  plot_manifest,
  file.path(export_path, "km_vs_mle_cure_plot_manifest.csv")
)
saveRDS(
  plot_registry,
  file.path(export_path, "km_vs_mle_cure_plot_registry.rds")
)

message("Integrated KM vs MLE cure plots saved to: ", normalizePath(export_path, winslash = "/", mustWork = FALSE))
