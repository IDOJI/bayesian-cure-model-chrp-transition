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
results_integration_root <- Sys.getenv(
  "BLOCK1_RESULTS_INTEGRATION_ROOT",
  unset = file.path(block1_root_default, "5.Results Integration")
)
legacy_root <- file.path(results_integration_root, "legacy")
r_model_root <- block1_root_default
rscript_bin <- Sys.getenv(
  "BLOCK1_RESULTS_INTEGRATION_RSCRIPT",
  unset = "/Library/Frameworks/R.framework/Resources/bin/Rscript"
)

required_scripts <- c(
  "export_mixture_cure_fraction_tables.r",
  "export_simple_model_probability_comparison_tables.r",
  "integrate_simple_model_probability_plots.r",
  "plot_km_vs_mle_cure_integrated_curves.r",
  "plot_km_vs_mle_nocure_cure_integrated_curves.r"
)

script_dir <- file.path(repo_root, "2.Rcode", "1.Block1", "Sub_Results Integration")
script_paths <- file.path(script_dir, required_scripts)

missing_scripts <- script_paths[!file.exists(script_paths)]
if (length(missing_scripts) > 0L) {
  stop(
    sprintf(
      "Missing required results-integration scripts:\n- %s",
      paste(normalizePath(missing_scripts, winslash = "/", mustWork = FALSE), collapse = "\n- ")
    ),
    call. = FALSE
  )
}

dir.create(results_integration_root, recursive = TRUE, showWarnings = FALSE)

archive_existing_results <- function(output_root, archive_root) {
  dir.create(archive_root, recursive = TRUE, showWarnings = FALSE)

  existing_paths <- list.files(
    output_root,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE
  )
  existing_paths <- existing_paths[basename(existing_paths) != "legacy"]

  if (length(existing_paths) == 0L) {
    return(invisible(NULL))
  }

  archive_name <- format(Sys.time(), "%Y%m%d_%H%M%S")
  archive_dir <- file.path(archive_root, archive_name)
  suffix <- 1L
  while (dir.exists(archive_dir) || file.exists(archive_dir)) {
    suffix <- suffix + 1L
    archive_dir <- file.path(archive_root, sprintf("%s_%02d", archive_name, suffix))
  }

  dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
  move_single_path <- function(source_path, target_path) {
    if (file.rename(source_path, target_path)) {
      return(TRUE)
    }

    copied_flag <- file.copy(
      from = source_path,
      to = target_path,
      recursive = TRUE,
      copy.mode = TRUE,
      copy.date = TRUE
    )

    if (!isTRUE(copied_flag)) {
      return(FALSE)
    }

    unlink(source_path, recursive = TRUE, force = TRUE)
    !file.exists(source_path)
  }

  moved_flag <- vapply(
    seq_along(existing_paths),
    function(ii) move_single_path(existing_paths[[ii]], file.path(archive_dir, basename(existing_paths[[ii]]))),
    logical(1)
  )

  if (any(!moved_flag)) {
    stop(
      sprintf(
        "Failed to archive existing results-integration contents:\n- %s",
        paste(basename(existing_paths[!moved_flag]), collapse = "\n- ")
      ),
      call. = FALSE
    )
  }

  invisible(normalizePath(archive_dir, winslash = "/", mustWork = TRUE))
}

write_results_integration_readme <- function(output_root, archived_dir = NULL) {
  archived_line <- if (is.null(archived_dir)) {
    "- No previous export bundle was present to archive for this run."
  } else {
    sprintf("- Previous contents were archived to `%s` before regeneration.", archived_dir)
  }

  readme_lines <- c(
    "# Block 1 Results Integration",
    "",
    "This folder contains the cleaned Block 1 cross-model integration outputs.",
    "",
    "## Active export folders",
    "",
    "- `mixture_cure_fraction_tables`: cure-fraction summaries and dataset-level comparison tables.",
    "- `simple_model_probability_comparison_tables`: yearly probability comparison tables across the four Block 1 models.",
    "- `integrated_model_probability_plots`: overall and susceptible-only cross-model probability plots.",
    "- `km_vs_mle_cure_integrated_curves`: KM versus frequentist mixture-cure plots with count panels.",
    "- `km_vs_mle_nocure_cure_integrated_curves`: KM versus no-cure versus mixture-cure plots with count panels.",
    "",
    "## Legacy handling",
    "",
    archived_line,
    "- Older result bundles are kept under the `legacy/` subfolder instead of being mixed with the active export set.",
    "",
    "## Regeneration entry point",
    "",
    "Run `run_block1_results_integration_exports.r` to archive the current bundle and rebuild the active outputs."
  )

  writeLines(readme_lines, con = file.path(output_root, "README.md"), useBytes = TRUE)
  invisible(file.path(output_root, "README.md"))
}

run_export_script <- function(script_path, output_root) {
  script_name <- basename(script_path)
  env_vars <- c(
    BLOCK1_RESULTS_INTEGRATION_ROOT = output_root,
    MIXTURE_CURE_RESULTS_INTEGRATION_ROOT = output_root,
    SIMPLE_MODELING_RESULTS_ROOT = r_model_root,
    MIXTURE_CURE_MODELING_RESULTS_ROOT = r_model_root,
    STAGE3_SIMPLE_EXPORT_PATH = file.path(r_model_root, "1.KM"),
    STAGE5_SIMPLE_EXPORT_PATH = file.path(r_model_root, "2.MLE No-Cure"),
    STAGE7_SIMPLE_EXPORT_PATH = file.path(r_model_root, "3.MLE Mixture Cure")
  )

  if (identical(script_name, "export_mixture_cure_fraction_tables.r")) {
    env_vars <- c(
      env_vars,
      MIXTURE_CURE_FRACTION_EXPORT_PATH = file.path(output_root, "mixture_cure_fraction_tables")
    )
  } else if (identical(script_name, "export_simple_model_probability_comparison_tables.r")) {
    env_vars <- c(
      env_vars,
      SIMPLE_MODEL_COMPARISON_TABLE_EXPORT_PATH = file.path(output_root, "simple_model_probability_comparison_tables")
    )
  } else if (identical(script_name, "integrate_simple_model_probability_plots.r")) {
    env_vars <- c(
      env_vars,
      SIMPLIFIED_PLOT_INTEGRATION_EXPORT_PATH = file.path(output_root, "integrated_model_probability_plots")
    )
  } else if (identical(script_name, "plot_km_vs_mle_cure_integrated_curves.r")) {
    env_vars <- c(
      env_vars,
      RESULTS_INTEGRATION_EXPORT_PATH = file.path(output_root, "km_vs_mle_cure_integrated_curves")
    )
  } else if (identical(script_name, "plot_km_vs_mle_nocure_cure_integrated_curves.r")) {
    env_vars <- c(
      env_vars,
      RESULTS_INTEGRATION_EXPORT_PATH = file.path(output_root, "km_vs_mle_nocure_cure_integrated_curves")
    )
  }

  env_names <- names(env_vars)
  old_env <- Sys.getenv(env_names, unset = NA_character_)
  names(old_env) <- env_names
  do.call(Sys.setenv, as.list(env_vars))
  on.exit({
    restore_values <- old_env[!is.na(old_env)]
    if (length(restore_values) > 0L) {
      do.call(Sys.setenv, as.list(restore_values))
    }
    unset_names <- names(old_env)[is.na(old_env)]
    if (length(unset_names) > 0L) {
      Sys.unsetenv(unset_names)
    }
  }, add = TRUE)

  exit_status <- system2(
    command = rscript_bin,
    args = c(shQuote(normalizePath(script_path, winslash = "/", mustWork = TRUE))),
    stdout = "",
    stderr = ""
  )

  if (!identical(exit_status, 0L)) {
    stop(
      sprintf(
        "Results-integration script failed: %s (exit status %s)",
        basename(script_path),
        as.character(exit_status)
      ),
      call. = FALSE
    )
  }

  invisible(exit_status)
}

archived_dir <- archive_existing_results(results_integration_root, legacy_root)

for (script_path in script_paths) {
  run_export_script(script_path, results_integration_root)
}

write_results_integration_readme(results_integration_root, archived_dir)

message("Block 1 results-integration bundle written to: ", normalizePath(results_integration_root, winslash = "/", mustWork = FALSE))
if (!is.null(archived_dir)) {
  message("Archived previous bundle to: ", archived_dir)
}
