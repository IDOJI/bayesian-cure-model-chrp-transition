# Configure Paths ---------------------------------------------------------
find_repo_root <- function(start_dir) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)

  repeat {
    simplified_reference <- file.path(
      current_dir,
      "2.Rcode",
      "Simplified",
      "stage3_KM benchmark classifier_simple.r"
    )
    if (file.exists(simplified_reference)) {
      return(current_dir)
    }

    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      break
    }
    current_dir <- parent_dir
  }

  stop(
    paste(
      "Could not locate the repository root containing",
      "`2.Rcode/Simplified/stage3_KM benchmark classifier_simple.r`."
    ),
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
results_root <- file.path(repo_root, "3.Results files")
export_path <- Sys.getenv(
  "SIMPLIFIED_PLOT_INTEGRATION_EXPORT_PATH",
  unset = file.path(results_root, "integrated_model_probability_plots")
)

stage3_survival_plot_file <- file.path(results_root, "stage3_simple_km_survival_plot.rds")
stage3_risk_plot_file <- file.path(results_root, "stage3_simple_km_risk_plot.rds")
stage5_survival_plot_file <- file.path(
  results_root,
  "stage5_individualized_nocure_lognormal_simple",
  "stage5_simple_lognormal_survival_plot.rds"
)
stage5_risk_plot_file <- file.path(
  results_root,
  "stage5_individualized_nocure_lognormal_simple",
  "stage5_simple_lognormal_risk_plot.rds"
)
stage7_plot_bundle_file <- file.path(results_root, "stage7_simple_lognormal_mixture_cure_plot_objects.rds")
stage8_plot_bundle_file <- file.path(results_root, "stage8A_simple_bayesian_transition_only_plot_objects.rds")

plot_width_in <- 12
plot_height_in <- 8
plot_dpi <- 320
horizon_years <- 1:10
dataset_order <- c("PNU", "SNU", "merged", "merged_site_adjusted")
dataset_label_lookup <- c(
  "PNU" = "PNU",
  "SNU" = "SNU",
  "merged" = "merged",
  "merged_site_adjusted" = "merged (site-adjusted)"
)
model_order <- c(
  "KM benchmark",
  "No-cure lognormal",
  "Frequentist mixture cure (lognormal)",
  "Bayesian transition-only cure"
)
model_palette <- c(
  "KM benchmark" = "#1B4332",
  "No-cure lognormal" = "#2A6F97",
  "Frequentist mixture cure (lognormal)" = "#C1666B",
  "Bayesian transition-only cure" = "#B8860B"
)
scope_order <- c("Overall cohort reference", "Susceptible subgroup estimate")
scope_linetypes <- c(
  "Overall cohort reference" = "22",
  "Susceptible subgroup estimate" = "solid"
)

# Load Packages -----------------------------------------------------------
required_packages <- c("ggplot2", "dplyr", "tibble", "scales")
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
  library(scales)
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

read_rds_checked <- function(path, label) {
  assert_file_exists(path, label)
  readRDS(path)
}

assert_is_ggplot <- function(plot_object, label) {
  if (!inherits(plot_object, "ggplot")) {
    stop(sprintf("%s is not a ggplot object.", label), call. = FALSE)
  }
  invisible(plot_object)
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

clip_probability <- function(x) {
  pmin(pmax(as.numeric(x), 0), 1)
}

assert_columns_present <- function(df, required_columns, label) {
  missing_columns <- setdiff(required_columns, names(df))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "%s is missing required columns: %s",
        label,
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(df)
}

extract_plot_curve_and_points <- function(
  plot_object,
  dataset_col,
  curve_time_col,
  point_time_col,
  value_col,
  model_name,
  draw_style = c("line", "step"),
  estimate_scope = "Overall cohort reference",
  label
) {
  draw_style <- match.arg(draw_style)
  assert_is_ggplot(plot_object, label)

  if (length(plot_object$layers) < 2L) {
    stop(sprintf("%s must contain at least two layers.", label), call. = FALSE)
  }

  curve_df <- tibble::as_tibble(plot_object$data)
  point_df <- tibble::as_tibble(plot_object$layers[[2L]]$data)

  assert_columns_present(curve_df, c(dataset_col, curve_time_col, value_col), paste(label, "curve data"))
  assert_columns_present(point_df, c(dataset_col, point_time_col, value_col), paste(label, "point data"))

  curve_out <- curve_df %>%
    transmute(
      dataset_key = normalize_dataset_name(.data[[dataset_col]]),
      dataset_label = unname(dataset_label_lookup[normalize_dataset_name(.data[[dataset_col]])]),
      model_name = model_name,
      estimate_scope = estimate_scope,
      draw_style = draw_style,
      time_year = as.numeric(.data[[curve_time_col]]),
      estimate_value = clip_probability(.data[[value_col]])
    ) %>%
    filter(
      !is.na(dataset_key),
      !is.na(dataset_label),
      !is.na(time_year),
      !is.na(estimate_value),
      time_year >= 0,
      time_year <= max(horizon_years)
    )

  point_out <- point_df %>%
    transmute(
      dataset_key = normalize_dataset_name(.data[[dataset_col]]),
      dataset_label = unname(dataset_label_lookup[normalize_dataset_name(.data[[dataset_col]])]),
      model_name = model_name,
      estimate_scope = estimate_scope,
      horizon_year = as.integer(round(as.numeric(.data[[point_time_col]]))),
      estimate_value = clip_probability(.data[[value_col]])
    ) %>%
    filter(
      !is.na(dataset_key),
      !is.na(dataset_label),
      !is.na(horizon_year),
      !is.na(estimate_value),
      horizon_year %in% horizon_years
    )

  list(curve = curve_out, points = point_out)
}

apply_plot_ordering <- function(df, include_draw_style = FALSE) {
  out <- df %>%
    mutate(
      dataset_key = factor(dataset_key, levels = dataset_order),
      dataset_label = factor(unname(dataset_label_lookup[as.character(dataset_key)]), levels = unname(dataset_label_lookup[dataset_order])),
      model_name = factor(model_name, levels = model_order),
      estimate_scope = factor(estimate_scope, levels = scope_order)
    ) %>%
    arrange(dataset_key, model_name)

  if (!include_draw_style) {
    return(out)
  }

  out %>%
    mutate(draw_style = factor(draw_style, levels = c("step", "line")))
}

duplicate_dataset_rows <- function(
  plot_source,
  source_dataset_key,
  target_dataset_key
) {
  source_dataset_key <- as.character(source_dataset_key)
  target_dataset_key <- as.character(target_dataset_key)
  target_dataset_label <- unname(dataset_label_lookup[[target_dataset_key]])

  curve_dup <- plot_source$curve %>%
    filter(as.character(dataset_key) == source_dataset_key) %>%
    mutate(
      dataset_key = target_dataset_key,
      dataset_label = target_dataset_label
    )

  point_dup <- plot_source$points %>%
    filter(as.character(dataset_key) == source_dataset_key) %>%
    mutate(
      dataset_key = target_dataset_key,
      dataset_label = target_dataset_label
    )

  list(
    curve = bind_rows(plot_source$curve, curve_dup),
    points = bind_rows(plot_source$points, point_dup)
  )
}

validate_plot_sources <- function(curve_df, point_df) {
  curve_required <- c(
    "dataset_key",
    "dataset_label",
    "model_name",
    "estimate_scope",
    "draw_style",
    "time_year",
    "estimate_value"
  )
  point_required <- c(
    "dataset_key",
    "dataset_label",
    "model_name",
    "estimate_scope",
    "horizon_year",
    "estimate_value"
  )

  assert_columns_present(curve_df, curve_required, "Integrated curve source")
  assert_columns_present(point_df, point_required, "Integrated point source")

  if (anyNA(curve_df$estimate_value) || any(curve_df$estimate_value < 0 | curve_df$estimate_value > 1)) {
    stop("Integrated curve source contains invalid probabilities.", call. = FALSE)
  }
  if (anyNA(point_df$estimate_value) || any(point_df$estimate_value < 0 | point_df$estimate_value > 1)) {
    stop("Integrated point source contains invalid probabilities.", call. = FALSE)
  }
  if (!identical(sort(unique(as.integer(point_df$horizon_year))), as.integer(horizon_years))) {
    stop("Integrated point source must contain horizons 1 through 10.", call. = FALSE)
  }

  invisible(TRUE)
}

make_integrated_plot <- function(
  curve_df,
  point_df,
  title_text,
  subtitle_text,
  y_label,
  include_scope_legend = FALSE,
  caption_text = NULL
) {
  curve_df <- apply_plot_ordering(curve_df, include_draw_style = TRUE)
  point_df <- apply_plot_ordering(point_df, include_draw_style = FALSE)
  validate_plot_sources(curve_df, point_df)

  plot_obj <- ggplot()

  step_curve_df <- curve_df %>% filter(draw_style == "step")
  line_curve_df <- curve_df %>% filter(draw_style == "line")

  if (include_scope_legend) {
    if (nrow(step_curve_df) > 0L) {
      plot_obj <- plot_obj +
        geom_step(
          data = step_curve_df,
          aes(
            x = time_year,
            y = estimate_value,
            color = model_name,
            linetype = estimate_scope,
            group = interaction(dataset_key, model_name, estimate_scope)
          ),
          linewidth = 1.05,
          direction = "hv"
        )
    }

    if (nrow(line_curve_df) > 0L) {
      plot_obj <- plot_obj +
        geom_line(
          data = line_curve_df,
          aes(
            x = time_year,
            y = estimate_value,
            color = model_name,
            linetype = estimate_scope,
            group = interaction(dataset_key, model_name, estimate_scope)
          ),
          linewidth = 1.05
        )
    }
  } else {
    if (nrow(step_curve_df) > 0L) {
      plot_obj <- plot_obj +
        geom_step(
          data = step_curve_df,
          aes(
            x = time_year,
            y = estimate_value,
            color = model_name,
            group = interaction(dataset_key, model_name)
          ),
          linewidth = 1.05,
          direction = "hv"
        )
    }

    if (nrow(line_curve_df) > 0L) {
      plot_obj <- plot_obj +
        geom_line(
          data = line_curve_df,
          aes(
            x = time_year,
            y = estimate_value,
            color = model_name,
            group = interaction(dataset_key, model_name)
          ),
          linewidth = 1.05
        )
    }
  }

  plot_obj <- plot_obj +
    geom_point(
      data = point_df,
      aes(
        x = horizon_year,
        y = estimate_value,
        color = model_name
      ),
      size = 1.9
    ) +
    facet_wrap(~dataset_label, ncol = 2) +
    scale_color_manual(values = model_palette, drop = FALSE) +
    scale_x_continuous(
      breaks = 0:10,
      limits = c(0, 10),
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
      x = "Years after cohort entry",
      y = y_label,
      color = "Model",
      caption = caption_text
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(hjust = 0)
    ) +
    guides(
      color = guide_legend(order = 1, nrow = 2, byrow = TRUE)
    )

  if (include_scope_legend) {
    plot_obj <- plot_obj +
      scale_linetype_manual(values = scope_linetypes, drop = FALSE) +
      guides(
        color = guide_legend(order = 1, nrow = 2, byrow = TRUE),
        linetype = guide_legend(order = 2)
      ) +
      labs(linetype = "Curve type")
  }

  plot_obj
}

save_plot_outputs <- function(plot_object, export_dir, file_stub) {
  png_file <- file.path(export_dir, paste0(file_stub, ".png"))
  rds_file <- file.path(export_dir, paste0(file_stub, ".rds"))

  saveRDS(plot_object, rds_file)
  ggplot2::ggsave(
    filename = png_file,
    plot = plot_object,
    width = plot_width_in,
    height = plot_height_in,
    dpi = plot_dpi,
    units = "in"
  )

  invisible(list(rds = rds_file, png = png_file))
}

save_dataset_specific_pngs <- function(
  curve_df,
  point_df,
  title_text,
  subtitle_text,
  y_label,
  include_scope_legend,
  file_stub_prefix,
  export_dir,
  caption_text = NULL
) {
  dataset_export_dir <- file.path(export_dir, "by_dataset")
  dir.create(dataset_export_dir, recursive = TRUE, showWarnings = FALSE)

  for (dataset_key_value in dataset_order) {
    filtered_curve_df <- curve_df %>%
      filter(as.character(dataset_key) == dataset_key_value)
    filtered_point_df <- point_df %>%
      filter(as.character(dataset_key) == dataset_key_value)

    if (nrow(filtered_curve_df) == 0L || nrow(filtered_point_df) == 0L) {
      next
    }

    single_plot <- make_integrated_plot(
      curve_df = filtered_curve_df,
      point_df = filtered_point_df,
      title_text = paste0(title_text, " - ", dataset_label_lookup[[dataset_key_value]]),
      subtitle_text = subtitle_text,
      y_label = y_label,
      include_scope_legend = include_scope_legend,
      caption_text = caption_text
    )

    png_file <- file.path(
      dataset_export_dir,
      paste0(file_stub_prefix, "_", dataset_key_value, ".png")
    )

    ggplot2::ggsave(
      filename = png_file,
      plot = single_plot,
      width = 8,
      height = 6,
      dpi = plot_dpi,
      units = "in"
    )
  }

  invisible(dataset_export_dir)
}

build_plot_sources <- function() {
  stage3_survival_plot <- read_rds_checked(stage3_survival_plot_file, "Stage 3 survival plot RDS")
  stage3_risk_plot <- read_rds_checked(stage3_risk_plot_file, "Stage 3 risk plot RDS")
  stage5_survival_plot <- read_rds_checked(stage5_survival_plot_file, "Stage 5 survival plot RDS")
  stage5_risk_plot <- read_rds_checked(stage5_risk_plot_file, "Stage 5 risk plot RDS")
  stage7_plot_bundle <- read_rds_checked(stage7_plot_bundle_file, "Stage 7 plot bundle RDS")
  stage8_plot_bundle <- read_rds_checked(stage8_plot_bundle_file, "Stage 8A plot bundle RDS")

  if (!is.list(stage7_plot_bundle) || !all(c(
    "overall_survival_curve",
    "overall_risk_curve",
    "susceptible_only_survival_curve",
    "susceptible_only_risk_curve"
  ) %in% names(stage7_plot_bundle))) {
    stop("Stage 7 plot bundle is missing one or more expected plot objects.", call. = FALSE)
  }

  if (!is.list(stage8_plot_bundle) || !all(c(
    "overall_survival_curve",
    "overall_risk_curve",
    "susceptible_only_survival_curve",
    "susceptible_only_risk_curve"
  ) %in% names(stage8_plot_bundle))) {
    stop("Stage 8A plot bundle is missing one or more expected plot objects.", call. = FALSE)
  }

  list(
    overall_survival = list(
      stage3 = duplicate_dataset_rows(
        plot_source = extract_plot_curve_and_points(
          plot_object = stage3_survival_plot,
          dataset_col = "dataset",
          curve_time_col = "time_year",
          point_time_col = "horizon_year",
          value_col = "estimated_survival_probability",
          model_name = "KM benchmark",
          draw_style = "step",
          estimate_scope = "Overall cohort reference",
          label = "Stage 3 survival plot"
        ),
        source_dataset_key = "merged",
        target_dataset_key = "merged_site_adjusted"
      ),
      stage5 = extract_plot_curve_and_points(
        plot_object = stage5_survival_plot,
        dataset_col = "dataset_version_key",
        curve_time_col = "time_year",
        point_time_col = "horizon_year",
        value_col = "estimated_survival_probability",
        model_name = "No-cure lognormal",
        draw_style = "line",
        estimate_scope = "Overall cohort reference",
        label = "Stage 5 survival plot"
      ),
      stage7 = extract_plot_curve_and_points(
        plot_object = stage7_plot_bundle[["overall_survival_curve"]],
        dataset_col = "dataset",
        curve_time_col = "time_horizon_year",
        point_time_col = "time_horizon_year",
        value_col = "overall_survival_prob",
        model_name = "Frequentist mixture cure (lognormal)",
        draw_style = "line",
        estimate_scope = "Overall cohort reference",
        label = "Stage 7 overall survival plot"
      ),
      stage8 = extract_plot_curve_and_points(
        plot_object = stage8_plot_bundle[["overall_survival_curve"]],
        dataset_col = "dataset",
        curve_time_col = "time_year",
        point_time_col = "time_year",
        value_col = "overall_survival_prob",
        model_name = "Bayesian transition-only cure",
        draw_style = "line",
        estimate_scope = "Overall cohort reference",
        label = "Stage 8A overall survival plot"
      )
    ),
    overall_risk = list(
      stage3 = duplicate_dataset_rows(
        plot_source = extract_plot_curve_and_points(
          plot_object = stage3_risk_plot,
          dataset_col = "dataset",
          curve_time_col = "time_year",
          point_time_col = "horizon_year",
          value_col = "estimated_risk_probability",
          model_name = "KM benchmark",
          draw_style = "step",
          estimate_scope = "Overall cohort reference",
          label = "Stage 3 risk plot"
        ),
        source_dataset_key = "merged",
        target_dataset_key = "merged_site_adjusted"
      ),
      stage5 = extract_plot_curve_and_points(
        plot_object = stage5_risk_plot,
        dataset_col = "dataset_version_key",
        curve_time_col = "time_year",
        point_time_col = "horizon_year",
        value_col = "estimated_risk_probability",
        model_name = "No-cure lognormal",
        draw_style = "line",
        estimate_scope = "Overall cohort reference",
        label = "Stage 5 risk plot"
      ),
      stage7 = extract_plot_curve_and_points(
        plot_object = stage7_plot_bundle[["overall_risk_curve"]],
        dataset_col = "dataset",
        curve_time_col = "time_horizon_year",
        point_time_col = "time_horizon_year",
        value_col = "overall_risk_prob",
        model_name = "Frequentist mixture cure (lognormal)",
        draw_style = "line",
        estimate_scope = "Overall cohort reference",
        label = "Stage 7 overall risk plot"
      ),
      stage8 = extract_plot_curve_and_points(
        plot_object = stage8_plot_bundle[["overall_risk_curve"]],
        dataset_col = "dataset",
        curve_time_col = "time_year",
        point_time_col = "time_year",
        value_col = "overall_risk_prob",
        model_name = "Bayesian transition-only cure",
        draw_style = "line",
        estimate_scope = "Overall cohort reference",
        label = "Stage 8A overall risk plot"
      )
    ),
    susceptible_survival = list(
      stage3 = duplicate_dataset_rows(
        plot_source = extract_plot_curve_and_points(
          plot_object = stage3_survival_plot,
          dataset_col = "dataset",
          curve_time_col = "time_year",
          point_time_col = "horizon_year",
          value_col = "estimated_survival_probability",
          model_name = "KM benchmark",
          draw_style = "step",
          estimate_scope = "Overall cohort reference",
          label = "Stage 3 survival plot"
        ),
        source_dataset_key = "merged",
        target_dataset_key = "merged_site_adjusted"
      ),
      stage5 = extract_plot_curve_and_points(
        plot_object = stage5_survival_plot,
        dataset_col = "dataset_version_key",
        curve_time_col = "time_year",
        point_time_col = "horizon_year",
        value_col = "estimated_survival_probability",
        model_name = "No-cure lognormal",
        draw_style = "line",
        estimate_scope = "Overall cohort reference",
        label = "Stage 5 survival plot"
      ),
      stage7 = extract_plot_curve_and_points(
        plot_object = stage7_plot_bundle[["susceptible_only_survival_curve"]],
        dataset_col = "dataset",
        curve_time_col = "time_horizon_year",
        point_time_col = "time_horizon_year",
        value_col = "susceptible_only_survival_prob",
        model_name = "Frequentist mixture cure (lognormal)",
        draw_style = "line",
        estimate_scope = "Susceptible subgroup estimate",
        label = "Stage 7 susceptible-only survival plot"
      ),
      stage8 = extract_plot_curve_and_points(
        plot_object = stage8_plot_bundle[["susceptible_only_survival_curve"]],
        dataset_col = "dataset",
        curve_time_col = "time_year",
        point_time_col = "time_year",
        value_col = "susceptible_only_survival_prob",
        model_name = "Bayesian transition-only cure",
        draw_style = "line",
        estimate_scope = "Susceptible subgroup estimate",
        label = "Stage 8A susceptible-only survival plot"
      )
    ),
    susceptible_risk = list(
      stage3 = duplicate_dataset_rows(
        plot_source = extract_plot_curve_and_points(
          plot_object = stage3_risk_plot,
          dataset_col = "dataset",
          curve_time_col = "time_year",
          point_time_col = "horizon_year",
          value_col = "estimated_risk_probability",
          model_name = "KM benchmark",
          draw_style = "step",
          estimate_scope = "Overall cohort reference",
          label = "Stage 3 risk plot"
        ),
        source_dataset_key = "merged",
        target_dataset_key = "merged_site_adjusted"
      ),
      stage5 = extract_plot_curve_and_points(
        plot_object = stage5_risk_plot,
        dataset_col = "dataset_version_key",
        curve_time_col = "time_year",
        point_time_col = "horizon_year",
        value_col = "estimated_risk_probability",
        model_name = "No-cure lognormal",
        draw_style = "line",
        estimate_scope = "Overall cohort reference",
        label = "Stage 5 risk plot"
      ),
      stage7 = extract_plot_curve_and_points(
        plot_object = stage7_plot_bundle[["susceptible_only_risk_curve"]],
        dataset_col = "dataset",
        curve_time_col = "time_horizon_year",
        point_time_col = "time_horizon_year",
        value_col = "susceptible_only_risk_prob",
        model_name = "Frequentist mixture cure (lognormal)",
        draw_style = "line",
        estimate_scope = "Susceptible subgroup estimate",
        label = "Stage 7 susceptible-only risk plot"
      ),
      stage8 = extract_plot_curve_and_points(
        plot_object = stage8_plot_bundle[["susceptible_only_risk_curve"]],
        dataset_col = "dataset",
        curve_time_col = "time_year",
        point_time_col = "time_year",
        value_col = "susceptible_only_risk_prob",
        model_name = "Bayesian transition-only cure",
        draw_style = "line",
        estimate_scope = "Susceptible subgroup estimate",
        label = "Stage 8A susceptible-only risk plot"
      )
    )
  )
}

collapse_sources <- function(source_list, key) {
  bind_rows(lapply(source_list, `[[`, key))
}

# Execute Plot Integration ------------------------------------------------
plot_sources <- build_plot_sources()

overall_survival_plot <- make_integrated_plot(
  curve_df = collapse_sources(plot_sources$overall_survival, "curve"),
  point_df = collapse_sources(plot_sources$overall_survival, "points"),
  title_text = "Integrated Overall Survival Probability Across Simplified Models",
  subtitle_text = "Facets show dataset versions; lines show cohort-level model predictions and yearly summaries",
  y_label = "Estimated survival probability",
  include_scope_legend = FALSE
)

overall_risk_plot <- make_integrated_plot(
  curve_df = collapse_sources(plot_sources$overall_risk, "curve"),
  point_df = collapse_sources(plot_sources$overall_risk, "points"),
  title_text = "Integrated Overall Risk Probability Across Simplified Models",
  subtitle_text = "Facets show dataset versions; lines show cohort-level model predictions and yearly summaries",
  y_label = "Estimated risk probability",
  include_scope_legend = FALSE
)

susceptible_survival_plot <- make_integrated_plot(
  curve_df = collapse_sources(plot_sources$susceptible_survival, "curve"),
  point_df = collapse_sources(plot_sources$susceptible_survival, "points"),
  title_text = "Integrated Susceptible-Group Survival Comparison Across Simplified Models",
  subtitle_text = paste(
    "Mixture cure models contribute susceptible-only curves;",
    "KM and no-cure models are overlaid as overall-cohort references."
  ),
  y_label = "Estimated survival probability",
  include_scope_legend = TRUE,
  caption_text = paste(
    "Reference lines from KM and no-cure models do not represent a separately estimated susceptible subgroup."
  )
)

susceptible_risk_plot <- make_integrated_plot(
  curve_df = collapse_sources(plot_sources$susceptible_risk, "curve"),
  point_df = collapse_sources(plot_sources$susceptible_risk, "points"),
  title_text = "Integrated Susceptible-Group Risk Comparison Across Simplified Models",
  subtitle_text = paste(
    "Mixture cure models contribute susceptible-only curves;",
    "KM and no-cure models are overlaid as overall-cohort references."
  ),
  y_label = "Estimated risk probability",
  include_scope_legend = TRUE,
  caption_text = paste(
    "Reference lines from KM and no-cure models do not represent a separately estimated susceptible subgroup."
  )
)

save_plot_outputs(
  overall_survival_plot,
  export_dir = export_path,
  file_stub = "integrated_overall_survival_probability_plot"
)
save_dataset_specific_pngs(
  curve_df = collapse_sources(plot_sources$overall_survival, "curve"),
  point_df = collapse_sources(plot_sources$overall_survival, "points"),
  title_text = "Integrated Overall Survival Probability Across Simplified Models",
  subtitle_text = "Dataset-specific export with the same model overlays and yearly summary points",
  y_label = "Estimated survival probability",
  include_scope_legend = FALSE,
  file_stub_prefix = "integrated_overall_survival_probability_plot",
  export_dir = export_path
)
save_plot_outputs(
  overall_risk_plot,
  export_dir = export_path,
  file_stub = "integrated_overall_risk_probability_plot"
)
save_dataset_specific_pngs(
  curve_df = collapse_sources(plot_sources$overall_risk, "curve"),
  point_df = collapse_sources(plot_sources$overall_risk, "points"),
  title_text = "Integrated Overall Risk Probability Across Simplified Models",
  subtitle_text = "Dataset-specific export with the same model overlays and yearly summary points",
  y_label = "Estimated risk probability",
  include_scope_legend = FALSE,
  file_stub_prefix = "integrated_overall_risk_probability_plot",
  export_dir = export_path
)
save_plot_outputs(
  susceptible_survival_plot,
  export_dir = export_path,
  file_stub = "integrated_susceptible_comparison_survival_probability_plot"
)
save_dataset_specific_pngs(
  curve_df = collapse_sources(plot_sources$susceptible_survival, "curve"),
  point_df = collapse_sources(plot_sources$susceptible_survival, "points"),
  title_text = "Integrated Susceptible-Group Survival Comparison Across Simplified Models",
  subtitle_text = "Dataset-specific export; reference lines are overall-cohort curves for KM and no-cure models",
  y_label = "Estimated survival probability",
  include_scope_legend = TRUE,
  file_stub_prefix = "integrated_susceptible_comparison_survival_probability_plot",
  export_dir = export_path,
  caption_text = paste(
    "Reference lines from KM and no-cure models do not represent a separately estimated susceptible subgroup."
  )
)
save_plot_outputs(
  susceptible_risk_plot,
  export_dir = export_path,
  file_stub = "integrated_susceptible_comparison_risk_probability_plot"
)
save_dataset_specific_pngs(
  curve_df = collapse_sources(plot_sources$susceptible_risk, "curve"),
  point_df = collapse_sources(plot_sources$susceptible_risk, "points"),
  title_text = "Integrated Susceptible-Group Risk Comparison Across Simplified Models",
  subtitle_text = "Dataset-specific export; reference lines are overall-cohort curves for KM and no-cure models",
  y_label = "Estimated risk probability",
  include_scope_legend = TRUE,
  file_stub_prefix = "integrated_susceptible_comparison_risk_probability_plot",
  export_dir = export_path,
  caption_text = paste(
    "Reference lines from KM and no-cure models do not represent a separately estimated susceptible subgroup."
  )
)

message("Integrated probability plots written to: ", export_path)
