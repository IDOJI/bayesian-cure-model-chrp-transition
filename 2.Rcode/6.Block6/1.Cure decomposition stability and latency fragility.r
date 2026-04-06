# Configure: Block 6 paths and export layout =====================================
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
current_script_file <- file.path(
  repo_root,
  "2.Rcode",
  "6.Block6",
  "1.Cure decomposition stability and latency fragility.r"
)
dropbox_project_root <- Sys.getenv(
  "BLOCK6_DROPBOX_PROJECT_ROOT",
  unset = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model"
)
export_path <- Sys.getenv(
  "BLOCK6_EXPORT_PATH",
  unset = file.path(dropbox_project_root, "6.Block6")
)

main_output_dir <- export_path
supporting_table_dir <- file.path(export_path, "2.Supporting Tables")
technical_dir <- file.path(export_path, "3.Technical")
staging_dir <- Sys.getenv(
  "BLOCK6_STAGING_DIR",
  unset = file.path(tempdir(), "block6_generated_raw")
)
legacy_generated_raw_dir <- file.path(technical_dir, "generated_raw")

analysis_dataset_file <- Sys.getenv(
  "BLOCK6_ANALYSIS_DATASETS_RDS",
  unset = file.path(dropbox_project_root, "2.Block2", "stage1_analysis_datasets.rds")
)
frequentist_fit_cache_file <- Sys.getenv(
  "BLOCK6_FREQUENTIST_FIT_CACHE",
  unset = file.path(dropbox_project_root, "1.Block1", "3.MLE Mixture Cure", "mle_mixture_cure_lognormal_fitted_objects.rds")
)
bayesian_model_registry_file <- Sys.getenv(
  "BLOCK6_BAYESIAN_MODEL_REGISTRY",
  unset = file.path(dropbox_project_root, "1.Block1", "4.Bayesian Mixture Cure", "bayesian_mixture_cure_model_registry.csv")
)
block4_late_tail_summary_file <- Sys.getenv(
  "BLOCK6_BLOCK4_LATE_TAIL_SUMMARY",
  unset = file.path(dropbox_project_root, "4.Block4", "block4_late_tail_narrative_summary.csv")
)

tail_n <- as.integer(Sys.getenv("BLOCK6_TAIL_N", unset = "5"))
joint_draws_max <- as.integer(Sys.getenv("BLOCK6_JOINT_DRAWS_MAX", unset = "1000"))

core_script <- file.path(
  repo_root,
  "2.Rcode",
  "1.Block1",
  "Main",
  "5.Latency instability priority diagnostics.r"
)

suppressPackageStartupMessages({
  library(readr)
  library(tibble)
  library(dplyr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE, scipen = 999)

# Helpers ========================================================================
`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}

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

safe_copy_file <- function(src, dest) {
  assert_exists(src, sprintf("Source file for %s", basename(dest)))
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  temp_path <- make_temp_export_path(dest)
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  ok <- file.copy(src, temp_path, overwrite = TRUE)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to copy file: %s", src), call. = FALSE)
  }
  ok <- replace_file_atomically(temp_path, dest)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to finalize copied file: %s", dest), call. = FALSE)
  }
  invisible(dest)
}

safe_write_lines_atomic <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_path <- make_temp_export_path(path)
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  writeLines(lines, con = temp_path)
  ok <- replace_file_atomically(temp_path, path)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to write text atomically: %s", path), call. = FALSE)
  }
  invisible(path)
}

normalize_existing_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

normalize_output_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

make_manifest_tbl <- function(file_paths, object_names, descriptions, export_root = export_path) {
  tibble(
    object_name = object_names,
    description = descriptions,
    file_path = normalize_output_path(file_paths),
    file_exists = file.exists(file_paths),
    relative_path = dplyr::if_else(
      startsWith(normalize_output_path(file_paths), paste0(normalize_output_path(export_root), "/")),
      sub(paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", normalize_output_path(export_root)), "/?"), "", normalize_output_path(file_paths)),
      normalize_output_path(file_paths)
    )
  )
}

write_block6_readme <- function(path) {
  lines <- c(
    "# Block 6 README",
    "",
    "Block 6 summarizes cure-decomposition stability and latency fragility.",
    "It standardizes three priority diagnostics that directly qualify cure-fraction and uncured-latency interpretation:",
    "",
    "- leave-one-out / leave-last-k-out influence analysis",
    "- latency-family sensitivity analysis",
    "- Bayesian joint uncertainty between cure fraction and latency scale",
    "- supported-horizon overlay against Block 4 primary-supported and latest-stable follow-up windows",
    "",
    "Top-level files are the interpretation-facing bundle.",
    "Top-level plots are rebuilt from exported CSV files so the visual outputs remain reproducible from saved tables.",
    "Detailed long-form tables and raw posterior-draw exports are stored under `2.Supporting Tables`.",
    "Execution manifest and code-reference files are stored under `3.Technical`.",
    "A file-by-file description is provided in `block6_file_guide.md`."
  )
  safe_write_lines_atomic(lines, path)
}

write_block6_file_guide <- function(path, entry_tbl) {
  section_order <- c("Top-level files", "Supporting tables", "Technical files")
  lines <- c(
    "# Block 6 File Guide",
    "",
    "This file explains the exported Block 6 artifacts."
  )

  for (section_name in section_order) {
    section_tbl <- entry_tbl %>%
      filter(section == section_name)

    if (nrow(section_tbl) == 0L) {
      next
    }

    lines <- c(lines, "", paste0("## ", section_name), "")
    for (ii in seq_len(nrow(section_tbl))) {
      lines <- c(
        lines,
        sprintf(
          "- `%s`: %s",
          basename(section_tbl$file_path[[ii]]),
          section_tbl$description[[ii]]
        )
      )
    }
  }

  safe_write_lines_atomic(lines, path)
}

make_block4_support_lookup <- function(block4_late_tail_path) {
  support_tbl <- readr::read_csv(block4_late_tail_path, show_col_types = FALSE, progress = FALSE) %>%
    filter(
      subgroup_kind == "overall",
      dataset %in% c("PNU", "SNU", "merged")
    ) %>%
    transmute(
      source_dataset = as.character(dataset),
      primary_supported_max_year = as.numeric(primary_supported_max_year),
      latest_stable_horizon_year = as.numeric(latest_stable_horizon_year),
      latest_supported_horizon_year = as.numeric(latest_supported_horizon_year),
      first_instability_horizon_year = as.numeric(first_instability_horizon_year),
      late_tail_status = as.character(late_tail_status),
      block4_narrative = as.character(narrative_text)
    )

  bind_rows(
    support_tbl %>%
      filter(source_dataset %in% c("PNU", "SNU")) %>%
      transmute(
        dataset = source_dataset,
        block4_source_dataset = source_dataset,
        primary_supported_max_year,
        latest_stable_horizon_year,
        latest_supported_horizon_year,
        first_instability_horizon_year,
        late_tail_status,
        block4_narrative
      ),
    support_tbl %>%
      filter(source_dataset == "merged") %>%
      transmute(
        dataset = "merged_no_site",
        block4_source_dataset = source_dataset,
        primary_supported_max_year,
        latest_stable_horizon_year,
        latest_supported_horizon_year,
        first_instability_horizon_year,
        late_tail_status,
        block4_narrative
      ),
    support_tbl %>%
      filter(source_dataset == "merged") %>%
      transmute(
        dataset = "merged_site_adjusted",
        block4_source_dataset = source_dataset,
        primary_supported_max_year,
        latest_stable_horizon_year,
        latest_supported_horizon_year,
        first_instability_horizon_year,
        late_tail_status,
        block4_narrative
      )
  )
}

safe_ratio <- function(x, y) {
  ifelse(is.na(x) | is.na(y) | y <= 0, NA_real_, x / y)
}

build_supported_horizon_overlay <- function(family_tbl, bayes_summary_tbl, support_lookup) {
  frequentist_overlay <- bind_rows(
    family_tbl %>%
      transmute(
        analysis_source = "frequentist",
        dataset = as.character(dataset),
        dataset_label = as.character(dataset_label),
        model_label = paste("frequentist", family, sep = "_"),
        latency_family = as.character(family),
        prior_branch = NA_character_,
        summary_metric = "weighted_subject_median",
        point_estimate_year = as.numeric(susceptible_weighted_subject_latency_median_year),
        interval_q025_year = NA_real_,
        interval_q975_year = NA_real_
      ),
    family_tbl %>%
      transmute(
        analysis_source = "frequentist",
        dataset = as.character(dataset),
        dataset_label = as.character(dataset_label),
        model_label = paste("frequentist", family, sep = "_"),
        latency_family = as.character(family),
        prior_branch = NA_character_,
        summary_metric = "cohort_uncured_median",
        point_estimate_year = as.numeric(cohort_uncured_median_latency_year),
        interval_q025_year = NA_real_,
        interval_q975_year = NA_real_
      )
  )

  bayesian_overlay <- bind_rows(
    bayes_summary_tbl %>%
      transmute(
        analysis_source = "bayesian",
        dataset = as.character(dataset),
        dataset_label = as.character(dataset_label),
        model_label = as.character(model_id),
        latency_family = "lognormal",
        prior_branch = as.character(prior_branch),
        summary_metric = "weighted_subject_median",
        point_estimate_year = as.numeric(weighted_subject_median_mean),
        interval_q025_year = as.numeric(weighted_subject_median_q025),
        interval_q975_year = as.numeric(weighted_subject_median_q975)
      ),
    bayes_summary_tbl %>%
      transmute(
        analysis_source = "bayesian",
        dataset = as.character(dataset),
        dataset_label = as.character(dataset_label),
        model_label = as.character(model_id),
        latency_family = "lognormal",
        prior_branch = as.character(prior_branch),
        summary_metric = "cohort_uncured_median",
        point_estimate_year = as.numeric(cohort_uncured_median_mean),
        interval_q025_year = as.numeric(cohort_uncured_median_q025),
        interval_q975_year = as.numeric(cohort_uncured_median_q975)
      )
  )

  bind_rows(frequentist_overlay, bayesian_overlay) %>%
    left_join(support_lookup, by = "dataset") %>%
    mutate(
      point_minus_primary_supported_year = point_estimate_year - primary_supported_max_year,
      point_minus_latest_stable_year = point_estimate_year - latest_stable_horizon_year,
      point_ratio_to_primary_supported = safe_ratio(point_estimate_year, primary_supported_max_year),
      point_ratio_to_latest_stable = safe_ratio(point_estimate_year, latest_stable_horizon_year),
      point_beyond_primary_supported = !is.na(point_estimate_year) & !is.na(primary_supported_max_year) & point_estimate_year > (primary_supported_max_year + 1e-10),
      point_beyond_latest_stable = !is.na(point_estimate_year) & !is.na(latest_stable_horizon_year) & point_estimate_year > (latest_stable_horizon_year + 1e-10),
      q025_beyond_primary_supported = !is.na(interval_q025_year) & !is.na(primary_supported_max_year) & interval_q025_year > (primary_supported_max_year + 1e-10),
      q025_beyond_latest_stable = !is.na(interval_q025_year) & !is.na(latest_stable_horizon_year) & interval_q025_year > (latest_stable_horizon_year + 1e-10)
    ) %>%
    arrange(
      match(dataset, c("PNU", "SNU", "merged_no_site", "merged_site_adjusted")),
      match(analysis_source, c("frequentist", "bayesian")),
      latency_family,
      prior_branch,
      summary_metric
    )
}

describe_support_status <- function(point_estimate_year, primary_supported_max_year, latest_stable_horizon_year) {
  if (is.na(point_estimate_year)) {
    return("support status unavailable")
  }
  if (!is.na(latest_stable_horizon_year) && point_estimate_year > (latest_stable_horizon_year + 1e-10)) {
    return(sprintf(
      "exceeds latest stable horizon %.0f y by %.2f y",
      latest_stable_horizon_year,
      point_estimate_year - latest_stable_horizon_year
    ))
  }
  if (!is.na(primary_supported_max_year) && point_estimate_year > (primary_supported_max_year + 1e-10)) {
    if (!is.na(latest_stable_horizon_year)) {
      return(sprintf(
        "exceeds primary-supported horizon %.0f y by %.2f y but remains below latest stable horizon %.0f y",
        primary_supported_max_year,
        point_estimate_year - primary_supported_max_year,
        latest_stable_horizon_year
      ))
    }
    return(sprintf(
      "exceeds primary-supported horizon %.0f y by %.2f y",
      primary_supported_max_year,
      point_estimate_year - primary_supported_max_year
    ))
  }
  if (!is.na(primary_supported_max_year)) {
    return(sprintf("remains within primary-supported horizon %.0f y", primary_supported_max_year))
  }
  "support status unavailable"
}

append_supported_horizon_section <- function(summary_path, overlay_tbl) {
  existing_lines <- readLines(summary_path, warn = FALSE)

  frequentist_focus <- overlay_tbl %>%
    filter(
      analysis_source == "frequentist",
      latency_family == "lognormal",
      summary_metric == "weighted_subject_median"
    ) %>%
    select(
      dataset,
      dataset_label,
      freq_point = point_estimate_year,
      primary_supported_max_year,
      latest_stable_horizon_year
    )

  bayesian_focus <- overlay_tbl %>%
    filter(
      analysis_source == "bayesian",
      prior_branch == "anchor_informed",
      summary_metric == "cohort_uncured_median"
    ) %>%
    select(
      dataset,
      bayes_point = point_estimate_year,
      bayes_q025 = interval_q025_year,
      bayes_q975 = interval_q975_year,
      primary_supported_max_year,
      latest_stable_horizon_year
    )

  overlay_focus <- frequentist_focus %>%
    left_join(
      bayesian_focus %>%
        select(dataset, bayes_point, bayes_q025, bayes_q975),
      by = "dataset"
    ) %>%
    arrange(match(dataset, c("PNU", "SNU", "merged_no_site", "merged_site_adjusted")))

  section_lines <- c("## Supported-horizon overlay", "")

  for (ii in seq_len(nrow(overlay_focus))) {
    row <- overlay_focus[ii, , drop = FALSE]
    freq_text <- describe_support_status(
      row$freq_point[[1]],
      row$primary_supported_max_year[[1]],
      row$latest_stable_horizon_year[[1]]
    )
    bayes_text <- describe_support_status(
      row$bayes_point[[1]],
      row$primary_supported_max_year[[1]],
      row$latest_stable_horizon_year[[1]]
    )

    section_lines <- c(
      section_lines,
      sprintf(
        "- %s: frequentist lognormal weighted median = %.2f y and %s; Bayesian anchor posterior mean uncured median = %.2f y (95%% interval %.2f to %.2f y) and %s.",
        row$dataset_label[[1]],
        row$freq_point[[1]],
        freq_text,
        row$bayes_point[[1]],
        row$bayes_q025[[1]],
        row$bayes_q975[[1]],
        bayes_text
      )
    )
  }

  safe_write_lines_atomic(c(existing_lines, "", section_lines), summary_path)
}

build_family_plot_from_csv <- function(path) {
  family_plot_tbl <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE) %>%
    filter(fit_ok)

  ggplot(
    family_plot_tbl,
    aes(x = family, y = susceptible_weighted_subject_latency_median_year, color = family)
  ) +
    geom_point(size = 2.6) +
    facet_wrap(~ dataset_label, scales = "free_y") +
    scale_y_log10(labels = scales::label_number(accuracy = 0.1)) +
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
}

build_joint_plot_from_csv <- function(path, y_var, y_label, log10_y = FALSE) {
  draw_tbl <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE) %>%
    mutate(
      prior_branch = factor(
        prior_branch,
        levels = c("anchor_informed", "neutral_no_external_info")
      )
    )

  plot_obj <- ggplot(
    draw_tbl,
    aes(x = cure_fraction, y = .data[[y_var]], color = prior_branch)
  ) +
    geom_point(alpha = 0.28, size = 0.9) +
    facet_wrap(~ dataset_label, scales = "free_y") +
    labs(
      title = paste0("Bayesian joint uncertainty: cure fraction vs ", y_label),
      x = "Cure fraction",
      y = y_label,
      color = "Prior branch"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  if (isTRUE(log10_y)) {
    plot_obj <- plot_obj + scale_y_log10(labels = scales::label_number(accuracy = 0.1))
  }

  plot_obj
}

remove_legacy_generated_raw_export <- function(legacy_dir, active_staging_dir) {
  legacy_norm <- normalize_output_path(legacy_dir)
  staging_norm <- normalize_output_path(active_staging_dir)

  if (dir.exists(legacy_dir) && !identical(legacy_norm, staging_norm)) {
    unlink(legacy_dir, recursive = TRUE, force = TRUE)
  }

  invisible(legacy_dir)
}

# Run core diagnostics ===========================================================
assert_exists(core_script, "Block 6 core diagnostics script")
assert_exists(analysis_dataset_file, "Block 6 analysis dataset RDS")
assert_exists(frequentist_fit_cache_file, "Block 6 frequentist fit cache")
assert_exists(bayesian_model_registry_file, "Block 6 Bayesian model registry")
assert_exists(block4_late_tail_summary_file, "Block 4 late-tail narrative summary")

dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supporting_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(technical_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(
  LATENCY_DIAG_ANALYSIS_DATASETS_RDS = analysis_dataset_file,
  LATENCY_DIAG_FREQUENTIST_FIT_CACHE = frequentist_fit_cache_file,
  LATENCY_DIAG_BAYESIAN_MODEL_REGISTRY = bayesian_model_registry_file,
  LATENCY_DIAG_OUTPUT_DIR = staging_dir,
  LATENCY_DIAG_TAIL_N = as.character(tail_n),
  LATENCY_DIAG_JOINT_DRAWS_MAX = as.character(joint_draws_max)
)

source(core_script, local = new.env(parent = globalenv()))

# Sync exported files into Block 6 layout ========================================
staging_files <- list(
  baseline_validation = file.path(staging_dir, "latency_instability_baseline_validation.csv"),
  influence_long = file.path(staging_dir, "latency_instability_influence_log_normal.csv"),
  influence_top = file.path(staging_dir, "latency_instability_influence_top.csv"),
  family_sensitivity = file.path(staging_dir, "latency_instability_family_sensitivity.csv"),
  family_plot = file.path(staging_dir, "latency_instability_family_sensitivity_plot.png"),
  bayes_draws = file.path(staging_dir, "latency_instability_bayesian_joint_draws.csv"),
  bayes_summary = file.path(staging_dir, "latency_instability_bayesian_joint_summary.csv"),
  bayes_sigma_plot = file.path(staging_dir, "latency_instability_bayesian_cure_vs_sigma.png"),
  bayes_uncured_median_plot = file.path(staging_dir, "latency_instability_bayesian_cure_vs_uncured_median.png"),
  summary_md = file.path(staging_dir, "latency_instability_priority_summary.md")
)

main_output_paths <- list(
  summary_md = file.path(main_output_dir, "block6_summary.md"),
  influence_top = file.path(main_output_dir, "block6_influence_top.csv"),
  family_sensitivity = file.path(main_output_dir, "block6_family_sensitivity.csv"),
  bayes_summary = file.path(main_output_dir, "block6_bayesian_joint_summary.csv"),
  supported_horizon_overlay = file.path(main_output_dir, "block6_supported_horizon_overlay.csv"),
  family_plot = file.path(main_output_dir, "block6_plot_family_sensitivity.png"),
  bayes_sigma_plot = file.path(main_output_dir, "block6_plot_bayesian_cure_vs_sigma.png"),
  bayes_uncured_median_plot = file.path(main_output_dir, "block6_plot_bayesian_cure_vs_uncured_median.png"),
  file_guide = file.path(main_output_dir, "block6_file_guide.md"),
  readme = file.path(main_output_dir, "README.md")
)

supporting_output_paths <- list(
  baseline_validation = file.path(supporting_table_dir, "block6_baseline_validation.csv"),
  influence_long = file.path(supporting_table_dir, "block6_influence_log_normal.csv"),
  bayes_draws = file.path(supporting_table_dir, "block6_bayesian_joint_draws.csv")
)

technical_output_paths <- list(
  export_manifest = file.path(technical_dir, "block6_export_manifest.csv"),
  code_reference = file.path(technical_dir, "block6_code_reference.csv")
)

safe_copy_file(staging_files$summary_md, main_output_paths$summary_md)
safe_copy_file(staging_files$influence_top, main_output_paths$influence_top)
safe_copy_file(staging_files$family_sensitivity, main_output_paths$family_sensitivity)
safe_copy_file(staging_files$bayes_summary, main_output_paths$bayes_summary)

safe_copy_file(staging_files$baseline_validation, supporting_output_paths$baseline_validation)
safe_copy_file(staging_files$influence_long, supporting_output_paths$influence_long)
safe_copy_file(staging_files$bayes_draws, supporting_output_paths$bayes_draws)

write_block6_readme(main_output_paths$readme)

family_tbl <- readr::read_csv(main_output_paths$family_sensitivity, show_col_types = FALSE, progress = FALSE)
bayes_summary_tbl <- readr::read_csv(main_output_paths$bayes_summary, show_col_types = FALSE, progress = FALSE)
support_lookup <- make_block4_support_lookup(block4_late_tail_summary_file)
supported_horizon_overlay_tbl <- build_supported_horizon_overlay(family_tbl, bayes_summary_tbl, support_lookup)
safe_write_csv_atomic(supported_horizon_overlay_tbl, main_output_paths$supported_horizon_overlay)
append_supported_horizon_section(main_output_paths$summary_md, supported_horizon_overlay_tbl)

family_plot <- build_family_plot_from_csv(main_output_paths$family_sensitivity)
bayes_sigma_plot <- build_joint_plot_from_csv(
  supporting_output_paths$bayes_draws,
  y_var = "latency_sigma",
  y_label = "Latency sigma",
  log10_y = FALSE
)
bayes_uncured_median_plot <- build_joint_plot_from_csv(
  supporting_output_paths$bayes_draws,
  y_var = "cohort_uncured_median_latency_year",
  y_label = "Cohort uncured median latency (years)",
  log10_y = TRUE
)

ggplot2::ggsave(main_output_paths$family_plot, family_plot, width = 10, height = 6, dpi = 320)
ggplot2::ggsave(main_output_paths$bayes_sigma_plot, bayes_sigma_plot, width = 10, height = 6, dpi = 320)
ggplot2::ggsave(main_output_paths$bayes_uncured_median_plot, bayes_uncured_median_plot, width = 10, height = 6, dpi = 320)

remove_legacy_generated_raw_export(legacy_generated_raw_dir, staging_dir)

code_reference_tbl <- tibble(
  code_role = c("block6_entry_script", "block6_core_diagnostics_script"),
  file_path = c(normalize_existing_path(current_script_file), normalize_existing_path(core_script)),
  description = c(
    "Main Block 6 entry script that configures Dropbox export layout and syncs interpretation-facing outputs.",
    "Core diagnostics engine reused by Block 6 for influence, family-sensitivity, and Bayesian joint-uncertainty analysis."
  )
)

safe_write_csv_atomic(code_reference_tbl, technical_output_paths$code_reference)

export_entry_tbl <- tibble(
  section = c(
    rep("Top-level files", 9),
    rep("Supporting tables", 3),
    rep("Technical files", 2)
  ),
  object_name = c(
    "summary_md",
    "influence_top",
    "family_sensitivity",
    "bayesian_joint_summary",
    "supported_horizon_overlay",
    "family_sensitivity_plot",
    "bayesian_cure_vs_sigma_plot",
    "bayesian_cure_vs_uncured_median_plot",
    "export_readme",
    "baseline_validation",
    "influence_long",
    "bayesian_joint_draws",
    "code_reference",
    "export_manifest"
  ),
  description = c(
    "Short Markdown interpretation summary for Block 6.",
    "Top influential leave-one-out / leave-last-k-out cases for Block 6.",
    "Top-level latency-family sensitivity comparison table for Block 6.",
    "Top-level Bayesian joint-uncertainty summary table for Block 6.",
    "Top-level overlay comparing model-implied latency medians against Block 4 primary-supported and latest-stable horizons.",
    "Top-level latency-family sensitivity plot rebuilt from the exported family-sensitivity CSV.",
    "Top-level Bayesian cure-vs-sigma scatter plot rebuilt from the exported Bayesian draw CSV.",
    "Top-level Bayesian cure-vs-uncured-median scatter plot rebuilt from the exported Bayesian draw CSV.",
    "Top-level guide to the Block 6 exported files.",
    "Supporting-table validation of reconstructed lognormal baseline metrics against saved Block 1 fits.",
    "Long-form leave-one-out / leave-last-k-out diagnostics for Block 6.",
    "Supporting-table posterior-draw export for Block 6 Bayesian joint uncertainty.",
    "Reference table for the Block 6 entry and core analysis scripts.",
    "Manifest of exported Block 6 files."
  ),
  file_path = c(
    main_output_paths$summary_md,
    main_output_paths$influence_top,
    main_output_paths$family_sensitivity,
    main_output_paths$bayes_summary,
    main_output_paths$supported_horizon_overlay,
    main_output_paths$family_plot,
    main_output_paths$bayes_sigma_plot,
    main_output_paths$bayes_uncured_median_plot,
    main_output_paths$readme,
    supporting_output_paths$baseline_validation,
    supporting_output_paths$influence_long,
    supporting_output_paths$bayes_draws,
    technical_output_paths$code_reference,
    technical_output_paths$export_manifest
  )
)

write_block6_file_guide(main_output_paths$file_guide, export_entry_tbl)

manifest_tbl <- make_manifest_tbl(
  file_paths = c(
    main_output_paths$summary_md,
    main_output_paths$influence_top,
    main_output_paths$family_sensitivity,
    main_output_paths$bayes_summary,
    main_output_paths$supported_horizon_overlay,
    main_output_paths$family_plot,
    main_output_paths$bayes_sigma_plot,
    main_output_paths$bayes_uncured_median_plot,
    main_output_paths$file_guide,
    main_output_paths$readme,
    unlist(supporting_output_paths, use.names = FALSE),
    technical_output_paths$code_reference,
    technical_output_paths$export_manifest
  ),
  object_names = c(
    "summary_md",
    "influence_top",
    "family_sensitivity",
    "bayesian_joint_summary",
    "supported_horizon_overlay",
    "family_sensitivity_plot",
    "bayesian_cure_vs_sigma_plot",
    "bayesian_cure_vs_uncured_median_plot",
    "file_guide_md",
    "export_readme",
    "baseline_validation",
    "influence_long",
    "bayesian_joint_draws",
    "code_reference",
    "export_manifest"
  ),
  descriptions = c(
    "Short Markdown interpretation summary for Block 6.",
    "Top influential leave-one-out / leave-last-k-out cases for Block 6.",
    "Top-level latency-family sensitivity comparison table for Block 6.",
    "Top-level Bayesian joint-uncertainty summary table for Block 6.",
    "Top-level overlay comparing model-implied latency medians against Block 4 primary-supported and latest-stable horizons.",
    "Top-level latency-family sensitivity plot rebuilt from the exported family-sensitivity CSV.",
    "Top-level Bayesian cure-vs-sigma scatter plot rebuilt from the exported Bayesian draw CSV.",
    "Top-level Bayesian cure-vs-uncured-median scatter plot rebuilt from the exported Bayesian draw CSV.",
    "Markdown file-by-file guide for the Block 6 export bundle.",
    "Top-level guide to the Block 6 exported files.",
    "Supporting-table validation of reconstructed lognormal baseline metrics against saved Block 1 fits.",
    "Long-form leave-one-out / leave-last-k-out diagnostics for Block 6.",
    "Supporting-table posterior-draw export for Block 6 Bayesian joint uncertainty.",
    "Reference table for the Block 6 entry and core analysis scripts.",
    "Manifest of exported Block 6 files."
  )
)

safe_write_csv_atomic(manifest_tbl, technical_output_paths$export_manifest)

message("Block 6 export completed:")
message("  Main output directory      : ", normalize_output_path(main_output_dir))
message("  Supporting tables directory: ", normalize_output_path(supporting_table_dir))
message("  Technical directory        : ", normalize_output_path(technical_dir))
