# 🔴 Configure: existing run paths and refresh switches ===============================
run_root_dir <- "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/STAGE6"
export_path <- run_root_dir

refresh_existing_outputs <- TRUE
refresh_stage6_export_manifest <- TRUE
update_stage6_screening_bundle <- FALSE

bootstrap_hist_bins <- 35L
bootstrap_plot_width <- 11
bootstrap_plot_height <- 8.5
bootstrap_density_adjust <- 1

# 🔴 Initialize: packages and runtime options ===============================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(stringr)
})

options(stringsAsFactors = FALSE, scipen = 999)

# 🔴 Define: reusable file and parsing helpers ===============================
assert_exists <- function(path, label) {
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
  temp_path <- make_temp_export_path(path)
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  readr::write_csv(df, temp_path)
  ok <- replace_file_atomically(temp_path, path)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to write CSV atomically: %s", path), call. = FALSE)
  }
  invisible(path)
}

safe_save_rds_atomic <- function(object, path) {
  temp_path <- make_temp_export_path(path)
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  saveRDS(object, temp_path)
  ok <- replace_file_atomically(temp_path, path)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to write RDS atomically: %s", path), call. = FALSE)
  }
  invisible(path)
}

safe_save_pdf_atomic <- function(path, width, height, plot_fun) {
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

as_integer_or_na <- function(x) {
  suppressWarnings(as.integer(as.character(x)))
}

as_numeric_or_na <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

safe_divide <- function(numerator, denominator) {
  ifelse(is.na(denominator) | denominator == 0, NA_real_, numerator / denominator)
}

extract_list_scalar <- function(x, name, default = NA) {
  if (is.null(x) || is.null(x[[name]]) || length(x[[name]]) == 0L) {
    return(default)
  }
  x[[name]][[1]]
}

strip_formula_rhs <- function(formula_text) {
  if (is.null(formula_text) || length(formula_text) == 0L) {
    return(NA_character_)
  }
  
  formula_text <- as.character(formula_text[[1]])
  if (is.na(formula_text) || !nzchar(formula_text)) {
    return(NA_character_)
  }
  
  gsub("^\\s*~\\s*", "", formula_text)
}

safe_quantile <- function(x, probs) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(rep(NA_real_, length(probs)))
  }
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, type = 8, names = FALSE))
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  mean(x)
}

safe_sd <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) <= 1L) {
    return(NA_real_)
  }
  stats::sd(x)
}

safe_min <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  min(x)
}

safe_max <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  max(x)
}

# 🔴 Define: bootstrap reconstruction helpers ===============================
collect_completed_bootstrap_stats_long <- function(registry_df, task_type, dataset_key, working_family = NA_character_) {
  subset_df <- registry_df %>%
    dplyr::filter(task_type == !!task_type, dataset_key == !!dataset_key, status == "completed")
  
  if (!is.na(working_family)) {
    subset_df <- subset_df %>% dplyr::filter(working_family == !!working_family)
  }
  
  if (nrow(subset_df) == 0L) {
    return(tibble(rep_id = integer(), statistic = numeric()))
  }
  
  shard_rows <- list()
  for (ii in seq_len(nrow(subset_df))) {
    shard_file <- subset_df$shard_file[[ii]]
    if (!file.exists(shard_file)) {
      next
    }
    
    obj <- tryCatch(readRDS(shard_file), error = function(e) NULL)
    if (is.null(obj)) {
      next
    }
    
    rep_ids <- as_integer_or_na(obj$rep_ids)
    stats <- as_numeric_or_na(obj$statistics)
    if (length(rep_ids) == 0L || length(stats) == 0L) {
      next
    }
    
    n_use <- min(length(rep_ids), length(stats))
    rep_ids <- rep_ids[seq_len(n_use)]
    stats <- stats[seq_len(n_use)]
    keep_idx <- which(!is.na(rep_ids))
    if (length(keep_idx) == 0L) {
      next
    }
    
    shard_rows[[length(shard_rows) + 1L]] <- tibble(
      rep_id = rep_ids[keep_idx],
      statistic = stats[keep_idx]
    )
  }
  
  if (length(shard_rows) == 0L) {
    return(tibble(rep_id = integer(), statistic = numeric()))
  }
  
  bind_rows(shard_rows) %>%
    group_by(rep_id) %>%
    slice_tail(n = 1L) %>%
    ungroup() %>%
    arrange(rep_id)
}

finalize_xie_from_bootstrap_stats <- function(observed_object, bootstrap_stats, alpha_screening) {
  observed <- observed_object$observed
  bootstrap_stats <- as.numeric(bootstrap_stats)
  
  if (is.na(observed$T_n)) {
    return(list(observed = observed, p_centered = NA_real_, p_raw = NA_real_, critical_centered = NA_real_))
  }
  
  if (isTRUE(observed$degenerate) && observed$T_n == 0) {
    return(list(observed = observed, p_centered = 1, p_raw = 1, critical_centered = 0))
  }
  
  centered_stats <- bootstrap_stats - observed$T_n
  n_valid_centered <- sum(is.finite(centered_stats))
  n_valid_raw <- sum(is.finite(bootstrap_stats))
  
  critical_centered <- if (n_valid_centered > 0L) {
    as.numeric(stats::quantile(centered_stats, probs = 1 - alpha_screening, na.rm = TRUE, type = 8, names = FALSE))
  } else {
    NA_real_
  }
  
  p_centered <- if (n_valid_centered > 0L) {
    (sum(centered_stats >= observed$T_n, na.rm = TRUE) + 1) / (n_valid_centered + 1)
  } else {
    NA_real_
  }
  
  p_raw <- if (n_valid_raw > 0L) {
    (sum(bootstrap_stats >= observed$T_n, na.rm = TRUE) + 1) / (n_valid_raw + 1)
  } else {
    NA_real_
  }
  
  list(observed = observed, p_centered = p_centered, p_raw = p_raw, critical_centered = critical_centered)
}

finalize_hsu_family_from_bootstrap_stats <- function(observed_object, bootstrap_stats, alpha_screening) {
  observed <- observed_object$result
  bootstrap_stats <- as.numeric(bootstrap_stats)
  working_family <- observed_object$working_family
  
  if (sum(is.finite(bootstrap_stats)) == 0L || is.null(observed$null_fit) || !isTRUE(observed$null_fit$converged) || is.na(observed$statistic)) {
    return(list(
      statistic = observed$statistic,
      p_value = NA_real_,
      critical_value = NA_real_,
      working_family = working_family,
      selected_formula = NA_character_,
      latency_formula = NA_character_,
      note = paste0("HSU bootstrap export could not finalize for family `", working_family, "`.")
    ))
  }
  
  n_valid <- sum(is.finite(bootstrap_stats))
  p_value <- if (n_valid > 0L) {
    (sum(bootstrap_stats >= observed$statistic, na.rm = TRUE) + 1) / (n_valid + 1)
  } else {
    NA_real_
  }
  
  critical_value <- if (n_valid > 0L) {
    as.numeric(stats::quantile(bootstrap_stats, probs = 1 - alpha_screening, na.rm = TRUE, type = 8, names = FALSE))
  } else {
    NA_real_
  }
  
  list(
    statistic = observed$statistic,
    p_value = p_value,
    critical_value = critical_value,
    working_family = working_family,
    selected_formula = extract_list_scalar(observed, "selected_formula", NA_character_),
    latency_formula = extract_list_scalar(observed, "latency_formula", NA_character_),
    note = paste0(
      "Self-contained pragmatic approximation to the Hsu sup-score test using a ",
      working_family,
      " AFT working latency model and sharded parametric bootstrap under the no-cure null."
    )
  )
}

hsu_method_name <- function(family) {
  switch(
    as.character(family),
    weibull = "hsu_supscore_weibull_approx",
    lognormal = "hsu_supscore_lognormal_approx",
    loglogistic = "hsu_supscore_loglogistic_approx",
    stop(sprintf("Unsupported HSU family for method naming: %s", family), call. = FALSE)
  )
}

source_dataset_from_key <- function(dataset_key) {
  dplyr::case_when(
    dataset_key == "PNU" ~ "PNU",
    dataset_key == "SNU" ~ "SNU",
    dataset_key %in% c("merged__site_free", "merged__site_adjusted") ~ "merged",
    TRUE ~ NA_character_
  )
}

analysis_variant_from_key <- function(dataset_key) {
  dplyr::case_when(
    dataset_key %in% c("PNU", "SNU") ~ "single_cohort",
    dataset_key == "merged__site_free" ~ "site_free",
    dataset_key == "merged__site_adjusted" ~ "site_adjusted",
    TRUE ~ NA_character_
  )
}

build_bootstrap_distribution_table <- function(registry_df, observed_dir, alpha_screening, candidate_hsu_families) {
  dataset_order_full <- c("PNU", "SNU", "merged__site_free", "merged__site_adjusted")
  xie_dataset_keys <- c("PNU", "SNU", "merged__site_free")
  
  out_rows <- list()
  
  for (dataset_key in xie_dataset_keys) {
    observed_file <- file.path(observed_dir, paste0("xie__", dataset_key, "__observed.rds"))
    if (!file.exists(observed_file)) {
      next
    }
    
    observed_object <- readRDS(observed_file)
    shard_tbl <- collect_completed_bootstrap_stats_long(registry_df, task_type = "xie", dataset_key = dataset_key)
    if (nrow(shard_tbl) == 0L) {
      next
    }
    
    finalized <- finalize_xie_from_bootstrap_stats(observed_object, shard_tbl$statistic, alpha_screening = alpha_screening)
    out_rows[[length(out_rows) + 1L]] <- tibble(
      record_type = "distribution",
      task_type = "xie",
      dataset_key = dataset_key,
      source_dataset = source_dataset_from_key(dataset_key),
      analysis_variant = analysis_variant_from_key(dataset_key),
      working_family = NA_character_,
      method_name = "xie_sufficient_followup",
      statistic_name = "T_n",
      rep_id = shard_tbl$rep_id,
      bootstrap_statistic = as.numeric(shard_tbl$statistic),
      centered_statistic = as.numeric(shard_tbl$statistic) - finalized$observed$T_n,
      observed_statistic = finalized$observed$T_n,
      p_value = finalized$p_centered,
      critical_value = finalized$critical_centered,
      critical_scale = "centered",
      plot_metric_name = "centered_bootstrap_statistic",
      plot_statistic_value = as.numeric(shard_tbl$statistic) - finalized$observed$T_n,
      plot_observed_line = finalized$observed$T_n,
      plot_critical_line = finalized$critical_centered,
      valid_statistic = is.finite(as.numeric(shard_tbl$statistic)),
      tau_year = as.numeric(observed_object$tau_year),
      selected_formula_rhs = NA_character_,
      latency_formula_rhs = NA_character_,
      note = "Centered Xie bootstrap distribution: histogram is drawn for T*_n - T_n and compared against the observed T_n and centered critical value."
    )
  }
  
  for (dataset_key in dataset_order_full) {
    for (family_name in candidate_hsu_families) {
      observed_file <- file.path(observed_dir, paste0("hsu__", dataset_key, "__", family_name, "__observed.rds"))
      if (!file.exists(observed_file)) {
        next
      }
      
      observed_object <- readRDS(observed_file)
      shard_tbl <- collect_completed_bootstrap_stats_long(
        registry_df,
        task_type = "hsu",
        dataset_key = dataset_key,
        working_family = family_name
      )
      if (nrow(shard_tbl) == 0L) {
        next
      }
      
      finalized <- finalize_hsu_family_from_bootstrap_stats(observed_object, shard_tbl$statistic, alpha_screening = alpha_screening)
      out_rows[[length(out_rows) + 1L]] <- tibble(
        record_type = "distribution",
        task_type = "hsu",
        dataset_key = dataset_key,
        source_dataset = source_dataset_from_key(dataset_key),
        analysis_variant = analysis_variant_from_key(dataset_key),
        working_family = family_name,
        method_name = hsu_method_name(family_name),
        statistic_name = "sup_lr_boot",
        rep_id = shard_tbl$rep_id,
        bootstrap_statistic = as.numeric(shard_tbl$statistic),
        centered_statistic = as.numeric(shard_tbl$statistic) - finalized$statistic,
        observed_statistic = finalized$statistic,
        p_value = finalized$p_value,
        critical_value = finalized$critical_value,
        critical_scale = "raw",
        plot_metric_name = "bootstrap_statistic",
        plot_statistic_value = as.numeric(shard_tbl$statistic),
        plot_observed_line = finalized$statistic,
        plot_critical_line = finalized$critical_value,
        valid_statistic = is.finite(as.numeric(shard_tbl$statistic)),
        tau_year = as.numeric(observed_object$tau_year),
        selected_formula_rhs = strip_formula_rhs(finalized$selected_formula),
        latency_formula_rhs = strip_formula_rhs(finalized$latency_formula),
        note = extract_list_scalar(finalized, "note", NA_character_)
      )
    }
  }
  
  dist_df <- bind_rows(out_rows)
  if (nrow(dist_df) == 0L) {
    return(tibble())
  }
  
  summary_df <- dist_df %>%
    group_by(
      task_type, dataset_key, source_dataset, analysis_variant, working_family,
      method_name, statistic_name, observed_statistic, p_value, critical_value,
      critical_scale, plot_metric_name, plot_observed_line, plot_critical_line,
      tau_year, selected_formula_rhs, latency_formula_rhs, note
    ) %>%
    summarise(
      n_bootstrap = n(),
      n_valid = sum(valid_statistic),
      valid_fraction = safe_divide(n_valid, n_bootstrap),
      bootstrap_mean = safe_mean(bootstrap_statistic),
      bootstrap_sd = safe_sd(bootstrap_statistic),
      bootstrap_min = safe_min(bootstrap_statistic),
      bootstrap_q025 = safe_quantile(bootstrap_statistic, probs = 0.025),
      bootstrap_q500 = safe_quantile(bootstrap_statistic, probs = 0.500),
      bootstrap_q975 = safe_quantile(bootstrap_statistic, probs = 0.975),
      centered_mean = safe_mean(centered_statistic),
      centered_sd = safe_sd(centered_statistic),
      centered_q025 = safe_quantile(centered_statistic, probs = 0.025),
      centered_q500 = safe_quantile(centered_statistic, probs = 0.500),
      centered_q975 = safe_quantile(centered_statistic, probs = 0.975),
      plot_min = safe_min(plot_statistic_value),
      plot_max = safe_max(plot_statistic_value),
      .groups = "drop"
    )
  
  dist_df %>%
    left_join(
      summary_df,
      by = c(
        "task_type", "dataset_key", "source_dataset", "analysis_variant", "working_family",
        "method_name", "statistic_name", "observed_statistic", "p_value", "critical_value",
        "critical_scale", "plot_metric_name", "plot_observed_line", "plot_critical_line",
        "tau_year", "selected_formula_rhs", "latency_formula_rhs", "note"
      )
    ) %>%
    arrange(
      factor(dataset_key, levels = dataset_order_full),
      task_type,
      working_family,
      rep_id
    )
}

generate_bootstrap_histogram_pdf <- function(bootstrap_df, final_path, bins, width, height, density_adjust) {
  dataset_order_full <- c("PNU", "SNU", "merged__site_free", "merged__site_adjusted")
  family_order <- c("weibull", "lognormal", "loglogistic")
  
  if (nrow(bootstrap_df) == 0L) {
    safe_save_pdf_atomic(final_path, width = width, height = height, plot_fun = function() {
      p <- ggplot() +
        annotate("text", x = 0, y = 0, label = "No completed bootstrap shard results were found.", size = 5) +
        theme_void()
      print(p)
    })
    return(invisible(TRUE))
  }
  
  plot_df <- bootstrap_df %>%
    mutate(
      plot_group_id = case_when(
        task_type == "xie" ~ paste0("xie__", dataset_key),
        TRUE ~ paste0("hsu__", dataset_key, "__", working_family)
      )
    )
  
  group_df <- plot_df %>%
    distinct(
      plot_group_id, task_type, dataset_key, working_family, method_name,
      plot_metric_name, observed_statistic, p_value, plot_observed_line,
      plot_critical_line, n_bootstrap, n_valid, valid_fraction
    ) %>%
    arrange(
      factor(dataset_key, levels = dataset_order_full),
      task_type,
      factor(working_family, levels = family_order)
    )
  
  safe_save_pdf_atomic(final_path, width = width, height = height, plot_fun = function() {
    for (ii in seq_len(nrow(group_df))) {
      current_group <- group_df[ii, , drop = FALSE]
      current_df <- plot_df %>% filter(plot_group_id == current_group$plot_group_id[[1]])
      
      title_text <- if (current_group$task_type[[1]] == "xie") {
        paste0("Stage 6 bootstrap distribution: Xie | ", current_group$dataset_key[[1]])
      } else {
        paste0(
          "Stage 6 bootstrap distribution: HSU-",
          current_group$working_family[[1]],
          " | ",
          current_group$dataset_key[[1]]
        )
      }
      
      subtitle_text <- paste0(
        "Observed = ", formatC(current_group$plot_observed_line[[1]], format = "f", digits = 4),
        "; Critical = ", formatC(current_group$plot_critical_line[[1]], format = "f", digits = 4),
        "; p = ", formatC(current_group$p_value[[1]], format = "f", digits = 4),
        "; valid = ", current_group$n_valid[[1]], "/", current_group$n_bootstrap[[1]]
      )
      
      p <- ggplot(current_df, aes(x = plot_statistic_value)) +
        geom_histogram(bins = bins, fill = "#4C78A8", colour = "white", linewidth = 0.3, na.rm = TRUE) +
        geom_density(adjust = density_adjust, linewidth = 0.6, colour = "#1F1F1F", na.rm = TRUE) +
        geom_vline(aes(xintercept = plot_observed_line), colour = "#D95F02", linewidth = 0.7) +
        geom_vline(aes(xintercept = plot_critical_line), colour = "#7570B3", linewidth = 0.7, linetype = "dashed") +
        labs(
          title = title_text,
          subtitle = subtitle_text,
          x = current_group$plot_metric_name[[1]],
          y = "Count"
        ) +
        theme_minimal(base_size = 11) +
        theme(panel.grid.minor = element_blank())
      
      print(p)
    }
  })
  
  invisible(TRUE)
}

# 🔴 Load: existing stage-six run products ===============================
run_root_dir <- normalizePath(run_root_dir, winslash = "/", mustWork = TRUE)
export_path <- normalizePath(export_path, winslash = "/", mustWork = FALSE)
observed_dir <- file.path(run_root_dir, "observed")
state_dir <- file.path(run_root_dir, "state")

stage6_shard_registry_file <- file.path(run_root_dir, "stage6_shard_registry.csv")
stage6_export_manifest_file <- file.path(export_path, "stage6_export_manifest.csv")
stage6_screening_bundle_file <- file.path(export_path, "stage6_screening_bundle.rds")

assert_exists(observed_dir, "observed_dir")
assert_exists(state_dir, "state_dir")
assert_exists(stage6_shard_registry_file, "stage6_shard_registry.csv")
assert_exists(stage6_screening_bundle_file, "stage6_screening_bundle.rds")

registry_df <- readr::read_csv(stage6_shard_registry_file, show_col_types = FALSE, progress = FALSE) %>%
  mutate(
    task_type = as.character(task_type),
    dataset_key = as.character(dataset_key),
    working_family = as.character(working_family),
    status = as.character(status)
  )

screening_bundle <- readRDS(stage6_screening_bundle_file)
state_bundle_file <- file.path(state_dir, "stage6_state_bundle.rds")
state_bundle <- if (file.exists(state_bundle_file)) readRDS(state_bundle_file) else screening_bundle

alpha_screening <- as.numeric(extract_list_scalar(state_bundle$config, "alpha_screening", 0.05))
candidate_hsu_families <- state_bundle$config$candidate_hsu_families
if (is.null(candidate_hsu_families) || length(candidate_hsu_families) == 0L) {
  candidate_hsu_families <- c("weibull", "lognormal", "loglogistic")
}
candidate_hsu_families <- as.character(candidate_hsu_families)
if (length(candidate_hsu_families) == 1L && grepl("\\|", candidate_hsu_families[[1]])) {
  candidate_hsu_families <- strsplit(candidate_hsu_families[[1]], "\\|")[[1]]
}

# 🔴 Build: bootstrap exports from existing shard results ===============================
bootstrap_export_df <- build_bootstrap_distribution_table(
  registry_df = registry_df,
  observed_dir = observed_dir,
  alpha_screening = alpha_screening,
  candidate_hsu_families = candidate_hsu_families
)

if (nrow(bootstrap_export_df) == 0L) {
  warning("No completed Stage 6 bootstrap shard results were found. The exported CSV and PDF will contain only placeholders.", call. = FALSE)
}

# 🔴 Export: refreshed bootstrap files and manifest entries ===============================
bootstrap_csv_file <- file.path(export_path, "stage6_bootstrap_results.csv")
bootstrap_pdf_file <- file.path(export_path, "stage6_bootstrap_histograms.pdf")

if (!isTRUE(refresh_existing_outputs) && (file.exists(bootstrap_csv_file) || file.exists(bootstrap_pdf_file))) {
  stop("Bootstrap export files already exist and `refresh_existing_outputs` is FALSE.", call. = FALSE)
}

safe_write_csv_atomic(bootstrap_export_df, bootstrap_csv_file)
generate_bootstrap_histogram_pdf(
  bootstrap_df = bootstrap_export_df,
  final_path = bootstrap_pdf_file,
  bins = bootstrap_hist_bins,
  width = bootstrap_plot_width,
  height = bootstrap_plot_height,
  density_adjust = bootstrap_density_adjust
)

if (isTRUE(update_stage6_screening_bundle) && is.list(screening_bundle) && !is.null(screening_bundle$outputs)) {
  screening_bundle$outputs$bootstrap_results <- bootstrap_export_df
  safe_save_rds_atomic(screening_bundle, stage6_screening_bundle_file)
}

if (isTRUE(refresh_stage6_export_manifest)) {
  manifest_updates <- tibble(
    file_name = c("stage6_bootstrap_results.csv", "stage6_bootstrap_histograms.pdf"),
    object_name = c("bootstrap_results", "bootstrap_histogram_pdf"),
    description = c(
      "Long-format bootstrap source-of-truth table reconstructed from completed Stage 6 shard files, including plot-ready columns and repeated group summaries.",
      "Multipage histogram PDF generated directly from stage6_bootstrap_results.csv."
    ),
    file_path = c(bootstrap_csv_file, bootstrap_pdf_file)
  )
  
  if (file.exists(stage6_export_manifest_file)) {
    manifest_df <- readr::read_csv(stage6_export_manifest_file, show_col_types = FALSE, progress = FALSE)
    manifest_df <- manifest_df %>%
      filter(!file_name %in% manifest_updates$file_name) %>%
      bind_rows(manifest_updates)
  } else {
    manifest_df <- manifest_updates
  }
  
  safe_write_csv_atomic(manifest_df, stage6_export_manifest_file)
}

message("Stage 6 bootstrap exports refreshed from existing shard results:")
message("  - ", bootstrap_csv_file)
message("  - ", bootstrap_pdf_file)
message("  - ", stage6_export_manifest_file)
