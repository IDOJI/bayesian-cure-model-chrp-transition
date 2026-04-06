# Configure Paths ---------------------------------------------------------
find_repo_root <- function(start_dir) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)

  repeat {
    has_repo_markers <- dir.exists(file.path(current_dir, ".git")) ||
      (
        dir.exists(file.path(current_dir, "0.Data")) &&
          dir.exists(file.path(current_dir, "2.Rcode"))
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
modeling_root <- Sys.getenv(
  "MIXTURE_CURE_MODELING_RESULTS_ROOT",
  unset = block1_root_default
)
results_root <- file.path(repo_root, "3.Results files")
results_integration_root <- Sys.getenv(
  "MIXTURE_CURE_RESULTS_INTEGRATION_ROOT",
  unset = file.path(block1_root_default, "5.Results Integration")
)
export_path <- Sys.getenv(
  "SIMPLIFIED_DATASET_COMPARISON_EXPORT_PATH",
  unset = file.path(results_integration_root, "dataset_level_model_comparisons")
)

stage3_survival_plot_candidates <- c(
  file.path(modeling_root, "1.KM", "simple_km_survival_plot.rds"),
  file.path(results_root, "stage3_simple_km_survival_plot.rds")
)
stage3_risk_plot_candidates <- c(
  file.path(modeling_root, "1.KM", "simple_km_risk_plot.rds"),
  file.path(results_root, "stage3_simple_km_risk_plot.rds")
)
stage5_survival_plot_candidates <- c(
  file.path(modeling_root, "2.MLE No-Cure", "mle_nocure_lognormal_survival_plot.rds"),
  file.path(modeling_root, "2.MLE No-Cure", "stage5_simple_lognormal_survival_plot.rds"),
  file.path(
    results_root,
    "stage5_individualized_nocure_lognormal_simple",
    "stage5_simple_lognormal_survival_plot.rds"
  )
)
stage5_risk_plot_candidates <- c(
  file.path(modeling_root, "2.MLE No-Cure", "mle_nocure_lognormal_risk_plot.rds"),
  file.path(modeling_root, "2.MLE No-Cure", "stage5_simple_lognormal_risk_plot.rds"),
  file.path(
    results_root,
    "stage5_individualized_nocure_lognormal_simple",
    "stage5_simple_lognormal_risk_plot.rds"
  )
)
stage7_curve_candidates <- c(
  file.path(modeling_root, "3.MLE Mixture Cure", "mle_mixture_cure_lognormal_plot_source.csv"),
  file.path(modeling_root, "3.MLE Mixture Cure", "stage7_simple_lognormal_mixture_cure_plot_source.csv"),
  file.path(results_root, "stage7_simple_lognormal_mixture_cure_plot_source.csv")
)
stage8_curve_candidates <- c(
  file.path(modeling_root, "4.Bayesian Mixture Cure", "bayesian_mixture_cure_plot_source.csv"),
  file.path(modeling_root, "4.Bayesian Mixture Cure", "bayesian_mixture_cure_curve_data.csv"),
  file.path(results_root, "stage8A_simple_bayesian_transition_only_curve_data.csv")
)
integrated_risk_dir_candidates <- c(
  file.path(results_integration_root, "integrated_risk_probabilities_by_horizon"),
  file.path(results_root, "integrated_risk_probabilities_by_horizon")
)

dataset_group_order <- c("PNU", "SNU", "Merged")
dataset_display_lookup <- c(
  "PNU" = "PNU",
  "SNU" = "SNU",
  "Merged" = "Merged"
)

pnu_snu_series_order <- c(
  "KM",
  "Freq No-Cure",
  "Freq Cure Overall",
  "Freq Cure Susceptible",
  "Bayes Cure Overall",
  "Bayes Cure Susceptible"
)

merged_series_order <- c(
  "KM",
  "Freq No-Cure",
  "Freq No-Cure\n(site-adj)",
  "Freq Cure Overall",
  "Freq Cure Susceptible",
  "Freq Cure Overall\n(site-adj)",
  "Freq Cure Susceptible\n(site-adj)",
  "Bayes Cure Overall",
  "Bayes Cure Susceptible",
  "Bayes Cure Overall\n(site-adj)",
  "Bayes Cure Susceptible\n(site-adj)"
)

all_series_order <- unique(c(pnu_snu_series_order, merged_series_order))
series_color_map <- c(
  "KM" = "#10B981",
  "Freq No-Cure" = "#F59E0B",
  "Freq No-Cure\n(site-adj)" = "#F59E0B",
  "Freq Cure Overall" = "#FF4D6D",
  "Freq Cure Susceptible" = "#FF4D6D",
  "Freq Cure Overall\n(site-adj)" = "#FF4D6D",
  "Freq Cure Susceptible\n(site-adj)" = "#FF4D6D",
  "Bayes Cure Overall" = "#3B82F6",
  "Bayes Cure Susceptible" = "#3B82F6",
  "Bayes Cure Overall\n(site-adj)" = "#3B82F6",
  "Bayes Cure Susceptible\n(site-adj)" = "#3B82F6"
)
series_linetype_map <- c(
  "KM" = "solid",
  "Freq No-Cure" = "solid",
  "Freq No-Cure\n(site-adj)" = "longdash",
  "Freq Cure Overall" = "solid",
  "Freq Cure Susceptible" = "dotdash",
  "Freq Cure Overall\n(site-adj)" = "longdash",
  "Freq Cure Susceptible\n(site-adj)" = "twodash",
  "Bayes Cure Overall" = "solid",
  "Bayes Cure Susceptible" = "dotdash",
  "Bayes Cure Overall\n(site-adj)" = "longdash",
  "Bayes Cure Susceptible\n(site-adj)" = "twodash"
)
heatmap_palette <- c("#FFF3B0", "#FFD166", "#F77F00", "#E63946")

# Load Packages -----------------------------------------------------------
required_packages <- c("ggplot2", "dplyr", "tibble", "readr", "scales")
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
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(readr)
  library(scales)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(export_path, "survival_curves"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(export_path, "risk_curves"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(export_path, "risk_heatmaps"), recursive = TRUE, showWarnings = FALSE)

# Define Helpers ----------------------------------------------------------
assert_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

resolve_existing_path <- function(paths, label) {
  existing_paths <- paths[file.exists(paths)]

  if (length(existing_paths) == 0L) {
    stop(
      sprintf(
        "%s not found. Checked paths:\n- %s",
        label,
        paste(paths, collapse = "\n- ")
      ),
      call. = FALSE
    )
  }

  normalizePath(existing_paths[[1L]], winslash = "/", mustWork = TRUE)
}

resolve_existing_dir <- function(paths, label) {
  existing_paths <- paths[dir.exists(paths)]

  if (length(existing_paths) == 0L) {
    stop(
      sprintf(
        "%s not found. Checked paths:\n- %s",
        label,
        paste(paths, collapse = "\n- ")
      ),
      call. = FALSE
    )
  }

  normalizePath(existing_paths[[1L]], winslash = "/", mustWork = TRUE)
}

assert_dir_exists <- function(path, label) {
  if (!dir.exists(path)) {
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

cleanup_existing_dataset_comparison_outputs <- function(export_dir) {
  stale_files <- list.files(
    path = export_dir,
    pattern = "\\.png$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(stale_files) > 0L) {
    file.remove(stale_files)
  }

  invisible(stale_files)
}

normalize_dataset_variant <- function(x) {
  x <- trimws(as.character(x))

  dplyr::case_when(
    toupper(x) == "PNU" ~ "PNU",
    toupper(x) == "SNU" ~ "SNU",
    tolower(x) %in% c("merged", "merged_no_site") ~ "merged",
    tolower(x) == "merged_site_adjusted" ~ "merged_site_adjusted",
    TRUE ~ x
  )
}

dataset_group_from_variant <- function(variant_key) {
  variant_key <- normalize_dataset_variant(variant_key)

  dplyr::case_when(
    variant_key == "PNU" ~ "PNU",
    variant_key == "SNU" ~ "SNU",
    variant_key %in% c("merged", "merged_site_adjusted") ~ "Merged",
    TRUE ~ variant_key
  )
}

series_order_for_group <- function(dataset_group) {
  if (identical(dataset_group, "Merged")) {
    return(merged_series_order)
  }
  pnu_snu_series_order
}

round_down_to_step <- function(x, step = 0.05) {
  floor(x / step) * step
}

round_up_to_step <- function(x, step = 0.05) {
  ceiling(x / step) * step
}

compute_heatmap_fill_limits <- function(probabilities) {
  valid_probabilities <- probabilities[is.finite(probabilities)]

  if (length(valid_probabilities) == 0L) {
    return(c(0, 1))
  }

  lower_quantile <- unname(stats::quantile(valid_probabilities, probs = 0.05, na.rm = TRUE))
  upper_quantile <- unname(stats::quantile(valid_probabilities, probs = 0.95, na.rm = TRUE))
  quantile_span <- upper_quantile - lower_quantile
  padding <- max(0.02, quantile_span * 0.08)

  fill_min <- max(0, round_down_to_step(lower_quantile - padding, step = 0.05))
  fill_max <- min(1, round_up_to_step(upper_quantile + padding, step = 0.05))

  if ((fill_max - fill_min) < 0.15) {
    fill_min <- max(0, fill_min - 0.05)
    fill_max <- min(1, fill_max + 0.05)
  }

  if (fill_max <= fill_min) {
    fill_max <- min(1, fill_min + 0.10)
  }

  c(fill_min, fill_max)
}

freq_nocure_label <- function(variant_key) {
  if (identical(variant_key, "merged_site_adjusted")) {
    return("Freq No-Cure\n(site-adj)")
  }
  "Freq No-Cure"
}

freq_cure_label <- function(variant_key, scope) {
  base_label <- if (identical(scope, "overall")) {
    "Freq Cure Overall"
  } else {
    "Freq Cure Susceptible"
  }

  if (identical(variant_key, "merged_site_adjusted")) {
    paste0(base_label, "\n(site-adj)")
  } else {
    base_label
  }
}

bayes_cure_label <- function(variant_key, scope) {
  base_label <- if (identical(scope, "overall")) {
    "Bayes Cure Overall"
  } else {
    "Bayes Cure Susceptible"
  }

  if (identical(variant_key, "merged_site_adjusted")) {
    paste0(base_label, "\n(site-adj)")
  } else {
    base_label
  }
}

build_stage3_curve_df <- function(plot_file, value_col, curve_type) {
  plot_object <- read_rds_checked(plot_file, paste("Stage 3", curve_type, "plot RDS"))
  tibble::as_tibble(plot_object$data) %>%
    transmute(
      dataset_variant = normalize_dataset_variant(dataset),
      dataset_group = dataset_group_from_variant(dataset),
      series_label = "KM",
      curve_type = curve_type,
      time_year = as.numeric(time_year),
      probability = as.numeric(.data[[value_col]])
    )
}

build_stage5_curve_df <- function(plot_file, value_col, curve_type) {
  plot_object <- read_rds_checked(plot_file, paste("Stage 5", curve_type, "plot RDS"))
  tibble::as_tibble(plot_object$data) %>%
    transmute(
      dataset_variant = normalize_dataset_variant(dataset_version_key),
      dataset_group = dataset_group_from_variant(dataset_version_key),
      series_label = vapply(
        normalize_dataset_variant(dataset_version_key),
        freq_nocure_label,
        character(1)
      ),
      curve_type = curve_type,
      time_year = as.numeric(time_year),
      probability = as.numeric(.data[[value_col]])
    )
}

build_stage7_curve_df <- function(path) {
  curve_df <- read_csv_checked(path, "Stage 7 curve source CSV") %>%
    mutate(
      dataset_variant = normalize_dataset_variant(dataset),
      dataset_group = dataset_group_from_variant(dataset)
    )
  time_col <- if ("time_horizon_year" %in% names(curve_df)) "time_horizon_year" else "time_year"

  bind_rows(
    curve_df %>%
      transmute(
        dataset_variant = dataset_variant,
        dataset_group = dataset_group,
        series_label = vapply(dataset_variant, freq_cure_label, character(1), scope = "overall"),
        curve_type = "survival",
        time_year = as.numeric(.data[[time_col]]),
        probability = as.numeric(overall_survival_prob)
      ),
    curve_df %>%
      transmute(
        dataset_variant = dataset_variant,
        dataset_group = dataset_group,
        series_label = vapply(dataset_variant, freq_cure_label, character(1), scope = "susceptible"),
        curve_type = "survival",
        time_year = as.numeric(.data[[time_col]]),
        probability = as.numeric(susceptible_only_survival_prob)
      ),
    curve_df %>%
      transmute(
        dataset_variant = dataset_variant,
        dataset_group = dataset_group,
        series_label = vapply(dataset_variant, freq_cure_label, character(1), scope = "overall"),
        curve_type = "risk",
        time_year = as.numeric(.data[[time_col]]),
        probability = as.numeric(overall_risk_prob)
      ),
    curve_df %>%
      transmute(
        dataset_variant = dataset_variant,
        dataset_group = dataset_group,
        series_label = vapply(dataset_variant, freq_cure_label, character(1), scope = "susceptible"),
        curve_type = "risk",
        time_year = as.numeric(.data[[time_col]]),
        probability = as.numeric(susceptible_only_risk_prob)
      )
  )
}

build_stage8_curve_df <- function(path) {
  curve_df <- read_csv_checked(path, "Stage 8A curve data CSV") %>%
    {
      if ("prior_branch" %in% names(.)) {
        dplyr::filter(., prior_branch == "anchor_informed")
      } else {
        .
      }
    } %>%
    mutate(
      dataset_variant = normalize_dataset_variant(dataset),
      dataset_group = dataset_group_from_variant(dataset)
    )
  time_col <- if ("time_year" %in% names(curve_df)) "time_year" else "time_horizon_year"

  bind_rows(
    curve_df %>%
      transmute(
        dataset_variant = dataset_variant,
        dataset_group = dataset_group,
        series_label = vapply(dataset_variant, bayes_cure_label, character(1), scope = "overall"),
        curve_type = "survival",
        time_year = as.numeric(.data[[time_col]]),
        probability = as.numeric(overall_survival_prob)
      ),
    curve_df %>%
      transmute(
        dataset_variant = dataset_variant,
        dataset_group = dataset_group,
        series_label = vapply(dataset_variant, bayes_cure_label, character(1), scope = "susceptible"),
        curve_type = "survival",
        time_year = as.numeric(.data[[time_col]]),
        probability = as.numeric(susceptible_only_survival_prob)
      ),
    curve_df %>%
      transmute(
        dataset_variant = dataset_variant,
        dataset_group = dataset_group,
        series_label = vapply(dataset_variant, bayes_cure_label, character(1), scope = "overall"),
        curve_type = "risk",
        time_year = as.numeric(.data[[time_col]]),
        probability = as.numeric(overall_risk_prob)
      ),
    curve_df %>%
      transmute(
        dataset_variant = dataset_variant,
        dataset_group = dataset_group,
        series_label = vapply(dataset_variant, bayes_cure_label, character(1), scope = "susceptible"),
        curve_type = "risk",
        time_year = as.numeric(.data[[time_col]]),
        probability = as.numeric(susceptible_only_risk_prob)
      )
  )
}

build_curve_comparison_df <- function() {
  bind_rows(
    build_stage3_curve_df(stage3_survival_plot_file, "estimated_survival_probability", "survival"),
    build_stage3_curve_df(stage3_risk_plot_file, "estimated_risk_probability", "risk"),
    build_stage5_curve_df(stage5_survival_plot_file, "estimated_survival_probability", "survival"),
    build_stage5_curve_df(stage5_risk_plot_file, "estimated_risk_probability", "risk"),
    build_stage7_curve_df(stage7_curve_file),
    build_stage8_curve_df(stage8_curve_file)
  ) %>%
    filter(
      dataset_group %in% dataset_group_order,
      series_label %in% all_series_order,
      !is.na(probability),
      !is.na(time_year),
      time_year >= 0,
      time_year <= 10
    )
}

build_yearly_risk_heatmap_df <- function() {
  assert_dir_exists(integrated_risk_dir, "Integrated risk probability folder")

  csv_files <- list.files(
    path = integrated_risk_dir,
    pattern = "^risk_probabilities_year_[0-9]{2}\\.csv$",
    full.names = TRUE
  )

  if (length(csv_files) == 0L) {
    stop("No integrated risk probability CSV files found.", call. = FALSE)
  }

  raw_df <- bind_rows(lapply(csv_files, function(csv_file) {
    horizon_year <- as.integer(sub("^.*_([0-9]{2})\\.csv$", "\\1", basename(csv_file)))

    readr::read_csv(csv_file, show_col_types = FALSE, progress = FALSE) %>%
      mutate(horizon_year = horizon_year)
  })) %>%
    mutate(
      dataset_variant = normalize_dataset_variant(dataset_name),
      dataset_group = dataset_group_from_variant(dataset_name)
    )

  raw_df %>%
    rowwise() %>%
    mutate(
      series_label = dplyr::case_when(
        dataset_variant %in% c("PNU", "SNU", "merged") && model_name == "KM benchmark" ~ "KM",
        dataset_variant %in% c("PNU", "SNU", "merged") && model_name == "No-cure lognormal" ~ "Freq No-Cure",
        dataset_variant == "merged_site_adjusted" && model_name == "No-cure lognormal" ~ "Freq No-Cure\n(site-adj)",
        dataset_variant %in% c("PNU", "SNU", "merged") && model_name == "Frequentist mixture cure (lognormal)" ~ "Freq Cure Overall",
        dataset_variant %in% c("PNU", "SNU", "merged") && model_name == "Frequentist mixture cure (lognormal, susceptible only)" ~ "Freq Cure Susceptible",
        dataset_variant == "merged_site_adjusted" && model_name == "Frequentist mixture cure (lognormal)" ~ "Freq Cure Overall\n(site-adj)",
        dataset_variant == "merged_site_adjusted" && model_name == "Frequentist mixture cure (lognormal, susceptible only)" ~ "Freq Cure Susceptible\n(site-adj)",
        dataset_variant %in% c("PNU", "SNU", "merged") && model_name == "Bayesian transition-only cure" ~ "Bayes Cure Overall",
        dataset_variant %in% c("PNU", "SNU", "merged") && model_name == "Bayesian transition-only cure (susceptible only)" ~ "Bayes Cure Susceptible",
        dataset_variant == "merged_site_adjusted" && model_name == "Bayesian transition-only cure" ~ "Bayes Cure Overall\n(site-adj)",
        dataset_variant == "merged_site_adjusted" && model_name == "Bayesian transition-only cure (susceptible only)" ~ "Bayes Cure Susceptible\n(site-adj)",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup() %>%
    filter(!is.na(series_label)) %>%
    filter(!(dataset_variant == "merged_site_adjusted" & model_name == "KM benchmark")) %>%
    transmute(
      dataset_group = dataset_group,
      series_label = series_label,
      horizon_year = as.integer(horizon_year),
      probability = as.numeric(estimated_risk_probability)
    ) %>%
    distinct() %>%
    arrange(dataset_group, series_label, horizon_year)
}

make_curve_plot <- function(curve_df, dataset_group_key, curve_type_key) {
  dataset_curve_df <- curve_df %>%
    filter(dataset_group == dataset_group_key, curve_type == curve_type_key) %>%
    mutate(
      series_label = factor(
        series_label,
        levels = series_order_for_group(dataset_group_key)
      )
    ) %>%
    arrange(series_label, time_year)

  if (nrow(dataset_curve_df) == 0L) {
    stop("No curve rows found for dataset group: ", dataset_group_key, call. = FALSE)
  }

  if (identical(curve_type_key, "risk")) {
    y_min <- 0
    y_max <- min(1, ceiling((max(dataset_curve_df$probability, na.rm = TRUE) + 0.05) * 10) / 10)
  } else {
    y_min <- max(0, floor((min(dataset_curve_df$probability, na.rm = TRUE) - 0.05) * 10) / 10)
    y_max <- 1
  }

  legend_columns <- if (identical(dataset_group_key, "Merged")) 3 else 3

  title_prefix <- if (identical(curve_type_key, "risk")) {
    "Risk Probability Curves"
  } else {
    "Survival Probability Curves"
  }

  subtitle_text <- if (identical(dataset_group_key, "Merged")) {
    paste(
      "Merged comparison includes non-adjusted and site-adjusted models,",
      "plus overall and susceptible-only outputs from mixture cure fits."
    )
  } else {
    paste(
      "Comparison includes KM, frequentist no-cure,",
      "and overall/susceptible outputs from frequentist and Bayesian mixture cure fits."
    )
  }

  ggplot(
    dataset_curve_df,
    aes(x = time_year, y = probability, color = series_label, linetype = series_label)
  ) +
    geom_line(linewidth = 1.05) +
    scale_color_manual(
      values = series_color_map,
      breaks = series_order_for_group(dataset_group_key),
      drop = FALSE,
      guide = guide_legend(
        ncol = legend_columns,
        byrow = TRUE,
        override.aes = list(linewidth = 1.25)
      )
    ) +
    scale_linetype_manual(
      values = series_linetype_map,
      breaks = series_order_for_group(dataset_group_key),
      drop = FALSE,
      guide = guide_legend(
        ncol = legend_columns,
        byrow = TRUE
      )
    ) +
    scale_x_continuous(
      breaks = 0:10,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(xlim = c(0, 10), ylim = c(y_min, y_max), expand = FALSE) +
    labs(
      title = paste0(dataset_display_lookup[[dataset_group_key]], ": ", title_prefix),
      subtitle = subtitle_text,
      x = "Years after cohort entry",
      y = if (identical(curve_type_key, "risk")) {
        "Estimated risk probability"
      } else {
        "Estimated survival probability"
      },
      color = "Series",
      linetype = "Series"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.width = grid::unit(2.8, "lines"),
      legend.key.height = grid::unit(0.9, "lines"),
      legend.spacing.x = grid::unit(0.5, "lines"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.margin = margin(t = 8, r = 12, b = 8, l = 8)
    )
}

make_heatmap_plot <- function(heatmap_df, dataset_group_key) {
  dataset_heatmap_df <- heatmap_df %>%
    filter(dataset_group == dataset_group_key)

  fill_limits <- compute_heatmap_fill_limits(dataset_heatmap_df$probability)
  label_threshold <- mean(fill_limits)

  dataset_heatmap_df <- dataset_heatmap_df %>%
    mutate(
      series_label = factor(
        series_label,
        levels = rev(series_order_for_group(dataset_group_key))
      ),
      label_text = scales::percent(probability, accuracy = 0.1),
      label_color = ifelse(probability >= label_threshold, "white", "black")
    )

  subtitle_text <- if (identical(dataset_group_key, "Merged")) {
    "Rows compare the 11 plotted series for Merged, including site-adjusted variants."
  } else {
    "Rows compare all plotted model outputs for this dataset."
  }

  ggplot(
    dataset_heatmap_df,
    aes(x = factor(horizon_year), y = series_label, fill = probability)
  ) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(
      aes(label = label_text, color = label_color),
      size = 3.1,
      fontface = "bold"
    ) +
    scale_color_identity() +
    scale_fill_gradientn(
      colours = heatmap_palette,
      limits = fill_limits,
      oob = scales::squish,
      labels = scales::label_percent(accuracy = 1),
      name = paste0(
        "Risk\n(",
        scales::percent(fill_limits[[1]], accuracy = 1),
        "-",
        scales::percent(fill_limits[[2]], accuracy = 1),
        ")"
      )
    ) +
    labs(
      title = paste0(dataset_display_lookup[[dataset_group_key]], ": Risk Probability Heatmap"),
      subtitle = subtitle_text,
      x = "Years after cohort entry",
      y = NULL,
      caption = paste0(
        "Fill scale is automatically adjusted to ",
        scales::percent(fill_limits[[1]], accuracy = 1),
        "-",
        scales::percent(fill_limits[[2]], accuracy = 1),
        " using the heatmap's central value range; values outside that range are saturated."
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold", lineheight = 0.95),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(hjust = 0),
      legend.position = "right",
      plot.margin = margin(t = 8, r = 12, b = 8, l = 8)
    )
}

save_plot_outputs <- function(plot_object, output_file, width_in, height_in, dpi = 320) {
  ggplot2::ggsave(
    filename = output_file,
    plot = plot_object,
    width = width_in,
    height = height_in,
    dpi = dpi,
    units = "in"
  )
  invisible(output_file)
}

# Execute Plot Generation -------------------------------------------------
stage3_survival_plot_file <- resolve_existing_path(
  stage3_survival_plot_candidates,
  "Stage 3 survival plot RDS"
)
stage3_risk_plot_file <- resolve_existing_path(
  stage3_risk_plot_candidates,
  "Stage 3 risk plot RDS"
)
stage5_survival_plot_file <- resolve_existing_path(
  stage5_survival_plot_candidates,
  "Stage 5 survival plot RDS"
)
stage5_risk_plot_file <- resolve_existing_path(
  stage5_risk_plot_candidates,
  "Stage 5 risk plot RDS"
)
stage7_curve_file <- resolve_existing_path(
  stage7_curve_candidates,
  "Stage 7 curve source CSV"
)
stage8_curve_file <- resolve_existing_path(
  stage8_curve_candidates,
  "Stage 8A curve data CSV"
)
integrated_risk_dir <- resolve_existing_dir(
  integrated_risk_dir_candidates,
  "Integrated risk probability folder"
)

curve_comparison_df <- build_curve_comparison_df()
heatmap_comparison_df <- build_yearly_risk_heatmap_df()
cleanup_existing_dataset_comparison_outputs(export_path)

for (dataset_group in dataset_group_order) {
  survival_plot <- make_curve_plot(
    curve_df = curve_comparison_df,
    dataset_group_key = dataset_group,
    curve_type_key = "survival"
  )
  risk_plot <- make_curve_plot(
    curve_df = curve_comparison_df,
    dataset_group_key = dataset_group,
    curve_type_key = "risk"
  )
  heatmap_plot <- make_heatmap_plot(
    heatmap_df = heatmap_comparison_df,
    dataset_group_key = dataset_group
  )

  curve_width <- if (identical(dataset_group, "Merged")) 14 else 12
  curve_height <- if (identical(dataset_group, "Merged")) 9.8 else 8.2
  heatmap_width <- if (identical(dataset_group, "Merged")) 12.5 else 10.5
  heatmap_height <- if (identical(dataset_group, "Merged")) 9.6 else 6.4

  save_plot_outputs(
    survival_plot,
    file.path(export_path, "survival_curves", paste0(tolower(dataset_group), "_survival_curve.png")),
    width_in = curve_width,
    height_in = curve_height
  )
  save_plot_outputs(
    risk_plot,
    file.path(export_path, "risk_curves", paste0(tolower(dataset_group), "_risk_curve.png")),
    width_in = curve_width,
    height_in = curve_height
  )
  save_plot_outputs(
    heatmap_plot,
    file.path(export_path, "risk_heatmaps", paste0(tolower(dataset_group), "_risk_heatmap.png")),
    width_in = heatmap_width,
    height_in = heatmap_height
  )
}

message("Dataset-level comparison plots written to: ", export_path)
