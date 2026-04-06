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
results_root <- file.path(repo_root, "3.Results files")
results_integration_root <- Sys.getenv(
  "MIXTURE_CURE_RESULTS_INTEGRATION_ROOT",
  unset = file.path(block1_root_default, "5.Results Integration")
)
input_path_candidates <- c(
  file.path(results_integration_root, "integrated_risk_probabilities_by_horizon"),
  file.path(results_root, "integrated_risk_probabilities_by_horizon")
)
input_path <- Sys.getenv(
  "SIMPLIFIED_RISK_VISUAL_INPUT_PATH",
  unset = input_path_candidates[[1L]]
)
export_path <- Sys.getenv(
  "SIMPLIFIED_RISK_VISUAL_EXPORT_PATH",
  unset = file.path(results_integration_root, "integrated_risk_probability_visualizations")
)

horizon_years <- 1:10
dataset_order <- c("PNU", "SNU", "merged", "merged_site_adjusted")
dataset_label_lookup <- c(
  "PNU" = "PNU",
  "SNU" = "SNU",
  "merged" = "Merged",
  "merged_site_adjusted" = "Merged (site-adjusted)"
)

overall_model_order <- c(
  "KM benchmark",
  "No-cure lognormal",
  "Frequentist mixture cure (lognormal)",
  "Bayesian transition-only cure"
)
overall_model_label_lookup <- c(
  "KM benchmark" = "KM",
  "No-cure lognormal" = "No-cure",
  "Frequentist mixture cure (lognormal)" = "Frequentist cure",
  "Bayesian transition-only cure" = "Bayesian cure"
)
overall_model_palette <- c(
  "KM" = "#00C853",
  "No-cure" = "#FFB300",
  "Frequentist cure" = "#FF4D6D",
  "Bayesian cure" = "#2979FF"
)

susceptible_model_order <- c(
  "KM benchmark",
  "No-cure lognormal",
  "Frequentist mixture cure (lognormal, susceptible only)",
  "Bayesian transition-only cure (susceptible only)"
)
susceptible_model_label_lookup <- c(
  "KM benchmark" = "KM",
  "No-cure lognormal" = "No-cure",
  "Frequentist mixture cure (lognormal, susceptible only)" = "Frequentist cure",
  "Bayesian transition-only cure (susceptible only)" = "Bayesian cure"
)
susceptible_model_palette <- c(
  "KM" = "#00C853",
  "No-cure" = "#FFB300",
  "Frequentist cure" = "#FF4D6D",
  "Bayesian cure" = "#2979FF"
)

all_model_order <- c(
  "KM benchmark",
  "No-cure lognormal",
  "Frequentist mixture cure (lognormal)",
  "Frequentist mixture cure (lognormal, susceptible only)",
  "Bayesian transition-only cure",
  "Bayesian transition-only cure (susceptible only)"
)
all_model_label_lookup <- c(
  "KM benchmark" = "KM",
  "No-cure lognormal" = "No-cure",
  "Frequentist mixture cure (lognormal)" = "Freq cure",
  "Frequentist mixture cure (lognormal, susceptible only)" = "Freq cure\n(susc.)",
  "Bayesian transition-only cure" = "Bayes cure",
  "Bayesian transition-only cure (susceptible only)" = "Bayes cure\n(susc.)"
)
all_model_family_lookup <- c(
  "KM benchmark" = "KM",
  "No-cure lognormal" = "No-cure",
  "Frequentist mixture cure (lognormal)" = "Frequentist cure",
  "Frequentist mixture cure (lognormal, susceptible only)" = "Frequentist cure",
  "Bayesian transition-only cure" = "Bayesian cure",
  "Bayesian transition-only cure (susceptible only)" = "Bayesian cure"
)
all_model_scope_lookup <- c(
  "KM benchmark" = "Overall/reference",
  "No-cure lognormal" = "Overall/reference",
  "Frequentist mixture cure (lognormal)" = "Overall/reference",
  "Frequentist mixture cure (lognormal, susceptible only)" = "Susceptible only",
  "Bayesian transition-only cure" = "Overall/reference",
  "Bayesian transition-only cure (susceptible only)" = "Susceptible only"
)
all_model_palette <- c(
  "KM" = "#00C853",
  "No-cure" = "#FFB300",
  "Frequentist cure" = "#FF4D6D",
  "Bayesian cure" = "#2979FF"
)
estimate_scope_alpha <- c(
  "Overall/reference" = 0.98,
  "Susceptible only" = 0.62
)
dataset_focus_order <- c("PNU", "SNU", "merged_pair")
dataset_focus_label_lookup <- c(
  "PNU" = "PNU",
  "SNU" = "SNU",
  "merged_pair" = "Merged vs site-adjusted"
)

heatmap_palette <- c("#FFF4CC", "#FFD166", "#F77F00", "#E63946")

# Load Packages -----------------------------------------------------------
required_packages <- c("ggplot2", "dplyr", "readr", "tibble", "scales", "stringr")
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
  library(readr)
  library(tibble)
  library(scales)
  library(stringr)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(export_path, "barplots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(export_path, "heatmaps"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(export_path, "horizon_specific_combined_barplots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(export_path, "dataset_horizon_barplots"), recursive = TRUE, showWarnings = FALSE)

# Define Helpers ----------------------------------------------------------
assert_dir_exists <- function(path, label) {
  if (!dir.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
  invisible(path)
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

cleanup_existing_visual_outputs <- function(export_dir) {
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

clip_probability <- function(x) {
  pmin(pmax(as.numeric(x), 0), 1)
}

read_integrated_risk_files <- function(path) {
  assert_dir_exists(path, "Integrated risk probability folder")

  csv_files <- list.files(
    path = path,
    pattern = "^risk_probabilities_year_[0-9]{2}\\.csv$",
    full.names = TRUE
  )

  if (length(csv_files) == 0L) {
    stop("No integrated risk probability CSV files were found.", call. = FALSE)
  }

  risk_df <- bind_rows(lapply(csv_files, function(csv_file) {
    horizon_match <- stringr::str_match(basename(csv_file), "risk_probabilities_year_([0-9]{2})\\.csv$")
    horizon_year <- as.integer(horizon_match[, 2])

    readr::read_csv(csv_file, show_col_types = FALSE, progress = FALSE) %>%
      mutate(
        horizon_year = horizon_year,
        estimated_risk_probability = clip_probability(estimated_risk_probability)
      )
  }))

  risk_df %>%
    mutate(
      dataset_name = as.character(dataset_name),
      model_name = as.character(model_name),
      horizon_year = as.integer(horizon_year)
    ) %>%
    filter(horizon_year %in% horizon_years)
}

prepare_plot_df <- function(df, model_order, model_label_lookup) {
  df %>%
    filter(model_name %in% model_order) %>%
    mutate(
      dataset_name = factor(dataset_name, levels = dataset_order),
      dataset_label = factor(
        unname(dataset_label_lookup[as.character(dataset_name)]),
        levels = unname(dataset_label_lookup[dataset_order])
      ),
      model_name = factor(model_name, levels = model_order),
      plot_model_label = factor(
        unname(model_label_lookup[as.character(model_name)]),
        levels = unname(model_label_lookup[model_order])
      ),
      horizon_label = factor(
        paste0(horizon_year, "-year"),
        levels = paste0(horizon_years, "-year")
      )
    ) %>%
    arrange(horizon_year, dataset_name, model_name)
}

prepare_dataset_horizon_plot_df <- function(df) {
  dataset_variant_levels <- unname(dataset_label_lookup[dataset_order])

  df %>%
    filter(model_name %in% all_model_order) %>%
    mutate(
      dataset_name = factor(dataset_name, levels = dataset_order),
      focus_key = dplyr::case_when(
        as.character(dataset_name) == "PNU" ~ "PNU",
        as.character(dataset_name) == "SNU" ~ "SNU",
        as.character(dataset_name) %in% c("merged", "merged_site_adjusted") ~ "merged_pair",
        TRUE ~ NA_character_
      ),
      focus_label = factor(
        unname(dataset_focus_label_lookup[focus_key]),
        levels = unname(dataset_focus_label_lookup[dataset_focus_order])
      ),
      dataset_variant_label = factor(
        unname(dataset_label_lookup[as.character(dataset_name)]),
        levels = dataset_variant_levels
      ),
      model_name = factor(model_name, levels = all_model_order),
      plot_model_label = factor(
        unname(all_model_label_lookup[as.character(model_name)]),
        levels = unname(all_model_label_lookup[all_model_order])
      ),
      model_family = factor(
        unname(all_model_family_lookup[as.character(model_name)]),
        levels = names(all_model_palette)
      ),
      estimate_scope = factor(
        unname(all_model_scope_lookup[as.character(model_name)]),
        levels = names(estimate_scope_alpha)
      ),
      label_text = scales::percent(estimated_risk_probability, accuracy = 0.1),
      label_y = pmin(estimated_risk_probability + 0.03, 1.04)
    ) %>%
    filter(!is.na(focus_key)) %>%
    arrange(horizon_year, focus_key, dataset_name, model_name)
}

make_bar_plot <- function(df, model_palette, title_text, subtitle_text, caption_text = NULL) {
  ggplot(df, aes(x = dataset_label, y = estimated_risk_probability, fill = plot_model_label)) +
    geom_col(
      position = position_dodge2(width = 0.82, preserve = "single"),
      width = 0.72
    ) +
    facet_wrap(~horizon_label, ncol = 5) +
    scale_fill_manual(values = model_palette, drop = FALSE) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = NULL,
      y = "Estimated risk probability",
      fill = "Model",
      caption = caption_text
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 25, hjust = 1),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(hjust = 0)
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}

make_dataset_horizon_bar_plot <- function(df, horizon_value, focus_key) {
  plot_df <- df %>%
    filter(horizon_year == horizon_value, focus_key == !!focus_key)

  if (nrow(plot_df) == 0L) {
    stop(
      sprintf(
        "No rows found for dataset focus `%s` at horizon %s.",
        focus_key,
        horizon_value
      ),
      call. = FALSE
    )
  }

  focus_label <- unique(as.character(plot_df$focus_label))
  focus_label <- focus_label[!is.na(focus_label)][1]
  is_merged_comparison <- identical(focus_key, "merged_pair")

  plot_object <- ggplot(
    plot_df,
    aes(
      x = plot_model_label,
      y = estimated_risk_probability,
      fill = model_family,
      alpha = estimate_scope
    )
  ) +
    geom_col(width = 0.72, color = "white", linewidth = 0.35) +
    geom_text(
      aes(y = label_y, label = label_text),
      size = 3.1,
      fontface = "bold",
      vjust = 0,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = all_model_palette, drop = FALSE) +
    scale_alpha_manual(values = estimate_scope_alpha, drop = FALSE) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0, 0.04))
    ) +
    coord_cartesian(ylim = c(0, 1.06), clip = "off") +
    labs(
      title = paste0(focus_label, ": ", horizon_value, "-year risk probability comparison"),
      subtitle = if (is_merged_comparison) {
        paste(
          "Merged and merged site-adjusted estimates are shown side by side.",
          "Lighter cure-model bars indicate susceptible-only estimates."
        )
      } else {
        paste(
          "Bars compare all model-based risk estimates for this dataset.",
          "Lighter cure-model bars indicate susceptible-only estimates."
        )
      },
      x = NULL,
      y = "Estimated risk probability",
      fill = "Model family",
      alpha = "Estimate scope",
      caption = "KM and no-cure values are overall-cohort references; only cure models add susceptible-only estimates."
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 9, lineheight = 0.9),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(hjust = 0),
      panel.grid.minor = element_blank()
    ) +
    guides(
      fill = guide_legend(order = 1, nrow = 1, byrow = TRUE),
      alpha = guide_legend(order = 2, nrow = 1, byrow = TRUE)
    )

  if (is_merged_comparison) {
    plot_object <- plot_object +
      facet_wrap(~dataset_variant_label, ncol = 2) +
      theme(strip.text = element_text(face = "bold"))
  }

  plot_object
}

build_combined_horizon_bar_df <- function(overall_df, susceptible_df) {
  bind_rows(
    overall_df %>%
      mutate(comparison_panel = "Overall risk"),
    susceptible_df %>%
      mutate(comparison_panel = "Susceptible comparison risk")
  ) %>%
    mutate(
      comparison_panel = factor(
        comparison_panel,
        levels = c("Overall risk", "Susceptible comparison risk")
      ),
      plot_model_label = factor(
        as.character(plot_model_label),
        levels = c("KM", "No-cure", "Frequentist cure", "Bayesian cure")
      ),
      label_text = scales::percent(estimated_risk_probability, accuracy = 0.1),
      label_y = pmin(estimated_risk_probability + 0.025, 1.02)
    )
}

make_horizon_combined_bar_plot <- function(df, horizon_value) {
  plot_df <- df %>%
    filter(horizon_year == horizon_value)

  bar_position <- position_dodge2(width = 0.86, preserve = "single", padding = 0.14)

  ggplot(plot_df, aes(x = dataset_label, y = estimated_risk_probability, fill = plot_model_label)) +
    geom_col(
      position = bar_position,
      width = 0.74,
      color = "white",
      linewidth = 0.3
    ) +
    geom_text(
      aes(y = label_y, label = label_text, group = plot_model_label),
      position = bar_position,
      size = 3,
      fontface = "bold",
      vjust = 0
    ) +
    facet_wrap(~comparison_panel, ncol = 1) +
    scale_fill_manual(values = overall_model_palette, drop = FALSE) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0, 0.04))
    ) +
    coord_cartesian(ylim = c(0, 1.05), clip = "off") +
    labs(
      title = paste0("Risk Comparison at ", horizon_value, "-Year Horizon"),
      subtitle = paste(
        "Overall and susceptible-comparison risks are shown together.",
        "Mixture cure bars in the lower panel use susceptible-only risk."
      ),
      x = NULL,
      y = "Estimated risk probability",
      fill = "Model",
      caption = "KM and no-cure bars in the lower panel are cohort-level references rather than separately estimated susceptible-subgroup risks."
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 18, hjust = 1),
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(hjust = 0),
      panel.spacing.y = grid::unit(1.2, "lines")
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}

make_heatmap_plot <- function(df, dataset_key, title_text, subtitle_text, caption_text = NULL) {
  dataset_df <- df %>%
    filter(as.character(dataset_name) == dataset_key) %>%
    mutate(
      plot_model_label = factor(
        as.character(plot_model_label),
        levels = rev(levels(plot_model_label))
      )
    )

  ggplot(dataset_df, aes(x = factor(horizon_year), y = plot_model_label, fill = estimated_risk_probability)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(
      aes(label = scales::percent(estimated_risk_probability, accuracy = 0.1)),
      size = 3.1,
      fontface = "bold"
    ) +
    scale_fill_gradientn(
      colours = heatmap_palette,
      limits = c(0, 1),
      labels = scales::label_percent(accuracy = 1),
      name = "Risk"
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Years after cohort entry",
      y = NULL,
      caption = caption_text
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(hjust = 0),
      legend.position = "right"
    )
}

save_png <- function(plot_object, output_file, width_in, height_in, dpi = 320) {
  ggplot2::ggsave(
    filename = output_file,
    plot = plot_object,
    width = width_in,
    height = height_in,
    dpi = dpi,
    units = "in"
  )
}

# Execute Visualization ---------------------------------------------------
if (identical(input_path, input_path_candidates[[1L]])) {
  input_path <- resolve_existing_dir(input_path_candidates, "Integrated risk probability folder")
}

risk_df <- read_integrated_risk_files(input_path)
cleanup_existing_visual_outputs(export_path)

overall_df <- prepare_plot_df(
  df = risk_df,
  model_order = overall_model_order,
  model_label_lookup = overall_model_label_lookup
)
susceptible_df <- prepare_plot_df(
  df = risk_df,
  model_order = susceptible_model_order,
  model_label_lookup = susceptible_model_label_lookup
)
combined_horizon_bar_df <- build_combined_horizon_bar_df(
  overall_df = overall_df,
  susceptible_df = susceptible_df
)
dataset_horizon_df <- prepare_dataset_horizon_plot_df(risk_df)

overall_bar_plot <- make_bar_plot(
  df = overall_df,
  model_palette = overall_model_palette,
  title_text = "Overall Risk Comparison by Horizon and Dataset",
  subtitle_text = "Each panel is a time horizon; bars compare the four datasets across the main model estimates."
)
save_png(
  plot_object = overall_bar_plot,
  output_file = file.path(export_path, "barplots", "overall_risk_comparison_barplot.png"),
  width_in = 16,
  height_in = 9
)

susceptible_bar_plot <- make_bar_plot(
  df = susceptible_df,
  model_palette = susceptible_model_palette,
  title_text = "Susceptible-Comparison Risk by Horizon and Dataset",
  subtitle_text = "Each panel is a time horizon; mixture cure bars use susceptible-only risks, while KM and no-cure remain reference values.",
  caption_text = "KM and no-cure bars are cohort-level references rather than separately estimated susceptible-subgroup risks."
)
save_png(
  plot_object = susceptible_bar_plot,
  output_file = file.path(export_path, "barplots", "susceptible_comparison_risk_barplot.png"),
  width_in = 16,
  height_in = 9
)

for (horizon_value in horizon_years) {
  horizon_plot <- make_horizon_combined_bar_plot(
    df = combined_horizon_bar_df,
    horizon_value = horizon_value
  )

  save_png(
    plot_object = horizon_plot,
    output_file = file.path(
      export_path,
      "horizon_specific_combined_barplots",
      paste0("risk_comparison_horizon_", sprintf("%02d", horizon_value), ".png")
    ),
    width_in = 11.5,
    height_in = 8.2
  )
}

for (focus_key in dataset_focus_order) {
  for (horizon_value in horizon_years) {
    dataset_horizon_plot <- make_dataset_horizon_bar_plot(
      df = dataset_horizon_df,
      horizon_value = horizon_value,
      focus_key = focus_key
    )

    save_png(
      plot_object = dataset_horizon_plot,
      output_file = file.path(
        export_path,
        "dataset_horizon_barplots",
        paste0(
          "risk_comparison_",
          focus_key,
          "_horizon_",
          sprintf("%02d", horizon_value),
          ".png"
        )
      ),
      width_in = if (identical(focus_key, "merged_pair")) 11.4 else 8.8,
      height_in = 6.2
    )
  }
}

for (dataset_key in dataset_order) {
  overall_heatmap <- make_heatmap_plot(
    df = overall_df,
    dataset_key = dataset_key,
    title_text = paste0("Overall Risk Heatmap - ", dataset_label_lookup[[dataset_key]]),
    subtitle_text = "Annotated heatmap for fast horizon-by-model comparison."
  )
  save_png(
    plot_object = overall_heatmap,
    output_file = file.path(
      export_path,
      "heatmaps",
      paste0("overall_risk_heatmap_", dataset_key, ".png")
    ),
    width_in = 9,
    height_in = 4.8
  )

  susceptible_heatmap <- make_heatmap_plot(
    df = susceptible_df,
    dataset_key = dataset_key,
    title_text = paste0("Susceptible-Comparison Risk Heatmap - ", dataset_label_lookup[[dataset_key]]),
    subtitle_text = "Mixture cure rows use susceptible-only risk; KM and no-cure remain reference values.",
    caption_text = "KM and no-cure rows are cohort-level references rather than separately estimated susceptible-subgroup risks."
  )
  save_png(
    plot_object = susceptible_heatmap,
    output_file = file.path(
      export_path,
      "heatmaps",
      paste0("susceptible_comparison_risk_heatmap_", dataset_key, ".png")
    ),
    width_in = 9,
    height_in = 4.8
  )
}

message("Integrated risk probability visualizations written to: ", export_path)
