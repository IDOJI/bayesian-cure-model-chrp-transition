# Configure: Block 5 paths and constants ========================================
dropbox_project_root <- Sys.getenv(
  "CHR_COHORTS_DROPBOX_ROOT",
  unset = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model"
)
stage1_data_path <- Sys.getenv(
  "BLOCK5_STAGE1_DATA_PATH",
  unset = file.path(dropbox_project_root, "old", "stage1_Backbone lock")
)
legacy_stage6_export_path <- Sys.getenv(
  "BLOCK5_LEGACY_STAGE6_EXPORT_PATH",
  unset = file.path(dropbox_project_root, "old", "stage6_Cure-appropriateness screening")
)
export_path <- Sys.getenv(
  "BLOCK5_EXPORT_PATH",
  unset = file.path(dropbox_project_root, "5.Block5", "archive_legacy_reuse")
)
main_output_dir <- export_path
supporting_table_dir <- file.path(export_path, "sub_supporting_tables")
bootstrap_detail_dir <- file.path(export_path, "sub_bootstrap")
metadata_dir <- file.path(export_path, "sub_metadata")

block5_spec_version <- "block5_cure_appropriateness_v1_20260402"
alpha_screening <- 0.05
receus_pi_threshold <- 0.025
receus_r_threshold <- 0.05
write_xie_bootstrap_histograms <- tolower(
  Sys.getenv("BLOCK5_WRITE_XIE_BOOTSTRAP_HISTOGRAMS", unset = "FALSE")
) %in% c("1", "true", "t", "yes", "y")

dataset_order <- c("PNU", "SNU", "merged")
dataset_map <- tibble::tribble(
  ~legacy_dataset_key, ~dataset, ~source_dataset, ~analysis_branch, ~dataset_label, ~source_note,
  "PNU", "PNU", "PNU", "single_cohort", "PNU", "Direct single-cohort Block 5 screening on the PNU subset.",
  "SNU", "SNU", "SNU", "single_cohort", "SNU", "Direct single-cohort Block 5 screening on the SNU subset.",
  "merged__site_free", "merged", "merged", "site_free", "Merged", "Merged raw-tail screening is reported on the site-free branch because RECeUS, Xie, and Maller-Zhou are raw-tail methods."
)

dataset_map_join <- dplyr::rename(
  dataset_map,
  block5_source_dataset = source_dataset,
  block5_analysis_branch = analysis_branch,
  block5_dataset_label = dataset_label,
  block5_source_note = source_note
)

reference_registry <- tibble::tribble(
  ~method_name, ~block5_role, ~reference_file, ~implementation_note,
  "RECeUS / RECeUS-AIC", "Primary cure-appropriateness screening", "R8. Selukar S, Othus M. RECeUS - Ratio estimation of censored uncured subjects, a different approach for assessing cure model appropriateness in studies with long-term survivors. Stat Med. 2023. (PubMed).pdf", "Primary Block 5 gate follows the RECeUS paper threshold logic pi_hat_n > 0.025 and r_hat_n < 0.05, with AIC-based family selection retained as the misspecification safeguard described in the paper.",
  "Xie sufficient follow-up test", "Secondary contradiction-oriented follow-up check", "R7. Xie P, Escobar-Bach M, Van Keilegom I. Testing for sufficient follow-up in censored survival data by using extremes. Biometrical Journal. 2024. (PubMed).pdf", "Block 5 keeps the centered bootstrap p-value as a contradiction check against follow-up adequacy and does not use Xie as a stand-alone positive gate.",
  "Maller-Zhou family", "Older contextual reference only", "R5. Maller RA, Zhou S. Testing for the presence of immune or cured individuals in censored survival data. Biometrics. 1995. (PubMed).pdf", "Block 5 reports the d_n statistic and legacy tail summaries contextually only and does not use them as the primary adjudication rule."
)

# Attach: packages and options ===================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(stringr)
})

options(stringsAsFactors = FALSE, scipen = 999)

# Define: helpers =================================================================
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

safe_save_plot <- function(plot_object, path, width, height, dpi = 320, bg = "white") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp_path <- make_temp_export_path(path)
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  ggplot2::ggsave(
    filename = temp_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = bg,
    limitsize = FALSE
  )
  ok <- replace_file_atomically(temp_path, path)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to write plot atomically: %s", path), call. = FALSE)
  }
  invisible(path)
}

safe_save_pdf_atomic <- function(path, width, height, plot_fun) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
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

safe_read_csv <- function(path) {
  assert_exists(path, label = basename(path))
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

safe_read_optional_csv <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

cleanup_duplicate_root_file <- function(root_dir, target_path) {
  root_path <- file.path(root_dir, basename(target_path))
  if (normalizePath(root_path, winslash = "/", mustWork = FALSE) == normalizePath(target_path, winslash = "/", mustWork = FALSE)) {
    return(invisible(FALSE))
  }
  if (file.exists(target_path) && file.exists(root_path)) {
    unlink(root_path)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

standard_error_column_patterns <- c(
  "(^|_)sd$",
  "(^|_)se$",
  "std_error$",
  "std\\.error$",
  "stderr$",
  "posterior_sd$",
  "uncertainty_sd$",
  "robust_se$"
)

identify_standard_error_columns <- function(df) {
  col_names <- names(df)
  if (is.null(col_names) || length(col_names) == 0L) {
    return(character())
  }

  col_names[vapply(
    col_names,
    function(col_name) any(grepl(standard_error_column_patterns, col_name, ignore.case = TRUE)),
    logical(1)
  )]
}

identify_standard_error_id_columns <- function(df) {
  preferred_cols <- c(
    "dataset",
    "dataset_label",
    "source_dataset",
    "analysis_branch",
    "family",
    "latency_family",
    "selected_by_aic",
    "receus_candidate_flag",
    "screening_class",
    "method_name"
  )

  intersect(preferred_cols, names(df))
}

empty_standard_error_registry <- function() {
  tibble(
    source_object = character(),
    n_rows = integer(),
    n_columns = integer(),
    n_standard_error_columns = integer(),
    standard_error_columns = character()
  )
}

empty_standard_error_long <- function() {
  tibble(
    source_object = character(),
    row_id = integer(),
    standard_error_column = character(),
    standard_error_value = numeric(),
    standard_error_value_raw = character()
  )
}

build_standard_error_registry_entry <- function(df, source_name) {
  se_cols <- identify_standard_error_columns(df)

  tibble(
    source_object = source_name,
    n_rows = nrow(df),
    n_columns = ncol(df),
    n_standard_error_columns = length(se_cols),
    standard_error_columns = if (length(se_cols) > 0L) paste(se_cols, collapse = "|") else NA_character_
  )
}

build_standard_error_long_table <- function(df, source_name) {
  se_cols <- identify_standard_error_columns(df)
  if (length(se_cols) == 0L || nrow(df) == 0L) {
    return(empty_standard_error_long())
  }

  base_df <- tibble::as_tibble(df) %>%
    mutate(row_id = dplyr::row_number())
  id_cols <- identify_standard_error_id_columns(base_df)

  bind_rows(lapply(se_cols, function(se_col) {
    out <- tibble(
      source_object = source_name,
      row_id = base_df$row_id,
      standard_error_column = se_col,
      standard_error_value = suppressWarnings(as.numeric(base_df[[se_col]])),
      standard_error_value_raw = as.character(base_df[[se_col]])
    )

    if (length(id_cols) > 0L) {
      out <- bind_cols(out, base_df[, id_cols, drop = FALSE])
    }

    out
  }))
}

build_standard_error_export_bundle <- function(named_tables) {
  named_tables <- named_tables[vapply(named_tables, function(x) inherits(x, "data.frame"), logical(1))]

  registry_tbl <- bind_rows(lapply(names(named_tables), function(table_name) {
    build_standard_error_registry_entry(named_tables[[table_name]], table_name)
  }))
  long_tbl <- bind_rows(lapply(names(named_tables), function(table_name) {
    build_standard_error_long_table(named_tables[[table_name]], table_name)
  }))

  if (nrow(registry_tbl) == 0L) {
    registry_tbl <- empty_standard_error_registry()
  }
  if (nrow(long_tbl) == 0L) {
    long_tbl <- empty_standard_error_long()
  }

  list(
    registry = registry_tbl,
    long = long_tbl
  )
}

collapse_pipe <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) {
    return(NA_character_)
  }
  paste(sort(x), collapse = "|")
}

make_receus_threshold_note <- function(pi_hat_n, r_hat_n) {
  if (!is.finite(pi_hat_n) || !is.finite(r_hat_n)) {
    return("RECeUS thresholds could not be evaluated because one or both summary quantities are missing.")
  }
  if (pi_hat_n > receus_pi_threshold && r_hat_n < receus_r_threshold) {
    return("RECeUS threshold rule is supportive because pi_hat_n exceeds 0.025 and r_hat_n is below 0.05.")
  }
  if (pi_hat_n <= receus_pi_threshold && r_hat_n < receus_r_threshold) {
    return("RECeUS threshold rule is unsupportive because pi_hat_n does not exceed 0.025.")
  }
  if (pi_hat_n > receus_pi_threshold && r_hat_n >= receus_r_threshold) {
    return("RECeUS threshold rule is unsupportive because r_hat_n is not below 0.05.")
  }
  "RECeUS threshold rule is unsupportive because pi_hat_n does not exceed 0.025 and r_hat_n is not below 0.05."
}

make_block5_screening_class <- function(receus_primary_flag, xie_contradiction_flag, maller_zhou_contextual_signal_flag) {
  if (is.na(receus_primary_flag) || !nzchar(receus_primary_flag)) {
    return("not_classified")
  }
  if (receus_primary_flag == "supportive" && isTRUE(xie_contradiction_flag)) {
    return("primary_supportive_but_xie_contradicted")
  }
  if (receus_primary_flag == "supportive") {
    return("primary_supportive_without_xie_contradiction")
  }
  if (receus_primary_flag == "equivocal" && isTRUE(xie_contradiction_flag)) {
    return("primary_equivocal_with_xie_contradiction")
  }
  if (receus_primary_flag == "equivocal") {
    return("primary_equivocal")
  }
  if (receus_primary_flag == "unsupportive" && isTRUE(maller_zhou_contextual_signal_flag)) {
    return("primary_unsupportive_with_contextual_tail_signal")
  }
  "primary_unsupportive"
}

make_block5_screening_note <- function(receus_primary_flag, xie_contradiction_flag, maller_zhou_contextual_signal_flag) {
  receus_sentence <- dplyr::case_when(
    receus_primary_flag == "supportive" ~ "Primary RECeUS-AIC supports cure-aware modeling on the paper's threshold rule.",
    receus_primary_flag == "equivocal" ~ "Primary RECeUS-AIC is equivocal on the paper's threshold rule.",
    receus_primary_flag == "unsupportive" ~ "Primary RECeUS-AIC does not support cure-aware modeling on the paper's threshold rule.",
    TRUE ~ "Primary RECeUS-AIC could not be classified."
  )
  xie_sentence <- if (isTRUE(xie_contradiction_flag)) {
    "Xie shows a contradiction signal against sufficient follow-up."
  } else {
    "Xie does not show a contradiction signal against sufficient follow-up."
  }
  mz_sentence <- if (isTRUE(maller_zhou_contextual_signal_flag)) {
    "Maller-Zhou shows a contextual tail signal, but it is not used as the primary adjudication rule."
  } else {
    "Maller-Zhou does not add a positive contextual tail signal."
  }
  paste(receus_sentence, xie_sentence, mz_sentence)
}

make_dataset_factor <- function(x) {
  factor(x, levels = dataset_order)
}

# Read: legacy Stage 6 inputs ====================================================
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

stage1_dataset_registry <- safe_read_csv(file.path(stage1_data_path, "stage1_dataset_registry.csv"))
legacy_followup_summary <- safe_read_csv(file.path(legacy_stage6_export_path, "stage6_followup_sufficiency_summary.csv"))
legacy_screening_summary <- safe_read_csv(file.path(legacy_stage6_export_path, "stage6_screening_summary.csv"))
legacy_screening_method_results <- safe_read_csv(file.path(legacy_stage6_export_path, "stage6_screening_method_results.csv"))
legacy_receus_candidate_fits <- safe_read_csv(file.path(legacy_stage6_export_path, "stage6_receus_candidate_fits.csv"))
legacy_bootstrap_results <- safe_read_optional_csv(file.path(legacy_stage6_export_path, "stage6_bootstrap_results.csv"))

# Validate: compatibility with the current Stage 1 backbone ======================
stage1_counts <- stage1_dataset_registry %>%
  filter(dataset %in% c("PNU", "SNU", "merged")) %>%
  transmute(dataset = as.character(dataset), stage1_n_rows = as.integer(n_rows))

legacy_counts <- legacy_followup_summary %>%
  inner_join(dataset_map_join, by = c("dataset_key" = "legacy_dataset_key")) %>%
  transmute(dataset = as.character(dataset), legacy_n_total = as.integer(n_total))

count_check <- stage1_counts %>%
  inner_join(legacy_counts, by = "dataset")

if (!all(count_check$stage1_n_rows == count_check$legacy_n_total)) {
  mismatch_text <- count_check %>%
    filter(stage1_n_rows != legacy_n_total) %>%
    transmute(message = paste0(dataset, ": stage1=", stage1_n_rows, ", legacy_stage6=", legacy_n_total)) %>%
    pull(message) %>%
    paste(collapse = "; ")
  stop(sprintf("Legacy Stage 6 outputs are not compatible with the current Stage 1 backbone. %s", mismatch_text), call. = FALSE)
}

# Transform: Block 5 data tables =================================================
data_backbone_summary <- legacy_followup_summary %>%
  inner_join(dataset_map_join, by = c("dataset_key" = "legacy_dataset_key")) %>%
  transmute(
    dataset,
    source_dataset = block5_source_dataset,
    analysis_branch = block5_analysis_branch,
    dataset_label = block5_dataset_label,
    n_total = as.integer(n_total),
    n_event = as.integer(n_event),
    n_censor_main = as.integer(n_censor_main),
    n_right_censor = as.integer(n_right_censor),
    n_remission = as.integer(n_remission),
    censor_rate_main = as.numeric(censor_rate_main),
    person_time_year = as.numeric(person_time_year),
    tau_year = as.numeric(tau_year),
    largest_event_time_year = as.numeric(largest_event_time_year),
    plateau_length_year = as.numeric(plateau_length_year),
    n_censored_after_last_event = as.integer(n_censored_after_last_event),
    km_tail_survival = as.numeric(km_tail_survival),
    km_tail_noncure = as.numeric(km_tail_noncure),
    xie_epsilon = as.numeric(xie_epsilon),
    xie_p_hat_n = as.numeric(xie_p_hat_n),
    xie_p_hat_G = as.numeric(xie_p_hat_G),
    xie_T_n = as.numeric(xie_T_n),
    maller_zhou_alpha_n = as.numeric(mz_alpha_n),
    maller_zhou_q_n = as.numeric(mz_q_n),
    source_note = block5_source_note
  ) %>%
  arrange(make_dataset_factor(dataset))

receus_candidate_fits <- legacy_receus_candidate_fits %>%
  inner_join(dataset_map_join, by = c("dataset_key" = "legacy_dataset_key")) %>%
  transmute(
    dataset,
    source_dataset = block5_source_dataset,
    analysis_branch = block5_analysis_branch,
    dataset_label = block5_dataset_label,
    family = as.character(family),
    converged = as.logical(converged),
    loglik = as.numeric(loglik),
    aic = as.numeric(aic),
    pi_hat_n = as.numeric(cure_fraction_hat),
    susceptible_survival_tau = as.numeric(susceptible_survival_tau),
    overall_survival_tau = as.numeric(overall_survival_tau),
    r_hat_n = as.numeric(receus_ratio_hat),
    tau_year = as.numeric(tau_year),
    receus_candidate_flag = as.character(receus_flag),
    selected_by_aic = as.logical(selected_by_aic),
    threshold_note = mapply(make_receus_threshold_note, as.numeric(cure_fraction_hat), as.numeric(receus_ratio_hat)),
    source_note = dplyr::coalesce(as.character(note), block5_source_note)
  ) %>%
  arrange(make_dataset_factor(dataset), aic, family)

receus_candidate_overview <- receus_candidate_fits %>%
  group_by(dataset, source_dataset, analysis_branch, dataset_label) %>%
  reframe(
    supportive_candidate_count = sum(receus_candidate_flag == "supportive", na.rm = TRUE),
    supportive_candidate_families = collapse_pipe(family[receus_candidate_flag == "supportive"]),
    best_aic_family = family[which.min(aic)][1]
  )

receus_primary_summary <- legacy_screening_summary %>%
  inner_join(dataset_map_join, by = c("dataset_key" = "legacy_dataset_key")) %>%
  transmute(
    dataset,
    source_dataset = block5_source_dataset,
    analysis_branch = block5_analysis_branch,
    dataset_label = block5_dataset_label,
    selected_family = as.character(receus_aic_selected_family),
    pi_hat_n = as.numeric(receus_aic_cure_fraction),
    r_hat_n = as.numeric(receus_aic_ratio),
    pi_threshold = receus_pi_threshold,
    r_threshold = receus_r_threshold,
    receus_primary_flag = as.character(receus_aic_flag),
    threshold_note = mapply(make_receus_threshold_note, as.numeric(receus_aic_cure_fraction), as.numeric(receus_aic_ratio)),
    source_note = block5_source_note
  ) %>%
  left_join(receus_candidate_overview, by = c("dataset", "source_dataset", "analysis_branch", "dataset_label")) %>%
  mutate(
    any_supportive_candidate_flag = supportive_candidate_count > 0L,
    interpretation_note = dplyr::case_when(
      receus_primary_flag == "supportive" ~ "AIC-selected RECeUS fit supports cure-appropriateness on the paper's threshold rule.",
      receus_primary_flag == "equivocal" ~ "AIC-selected RECeUS fit is equivocal on the paper's threshold rule.",
      receus_primary_flag == "unsupportive" & any_supportive_candidate_flag ~ "AIC-selected RECeUS fit is unsupportive even though at least one nonselected candidate family crosses the paper's support threshold.",
      receus_primary_flag == "unsupportive" ~ "AIC-selected RECeUS fit is unsupportive on the paper's threshold rule.",
      TRUE ~ "RECeUS primary result could not be classified."
    )
  ) %>%
  arrange(make_dataset_factor(dataset))

xie_contradiction_summary <- legacy_screening_summary %>%
  inner_join(dataset_map_join, by = c("dataset_key" = "legacy_dataset_key")) %>%
  left_join(
    legacy_screening_method_results %>%
      filter(method_name == "xie_sufficient_followup") %>%
      transmute(
        dataset_key,
        xie_bootstrap_p_value = as.numeric(p_value),
        xie_method_note = as.character(note)
      ),
    by = "dataset_key"
  ) %>%
  transmute(
    dataset,
    source_dataset = block5_source_dataset,
    analysis_branch = block5_analysis_branch,
    dataset_label = block5_dataset_label,
    xie_epsilon = as.numeric(xie_epsilon),
    p_hat_n = as.numeric(xie_p_hat_n),
    p_hat_G = as.numeric(xie_p_hat_G),
    T_n = as.numeric(xie_T_n),
    bootstrap_p_value = as.numeric(xie_bootstrap_p_value),
    alpha_screening = alpha_screening,
    contradiction_flag = !is.na(bootstrap_p_value) & bootstrap_p_value < alpha_screening,
    interpretation_note = dplyr::case_when(
      !is.na(bootstrap_p_value) & bootstrap_p_value < alpha_screening ~ "Xie provides a contradiction signal against sufficient follow-up at alpha = 0.05.",
      !is.na(bootstrap_p_value) ~ "Xie does not provide a contradiction signal against sufficient follow-up at alpha = 0.05.",
      TRUE ~ "Xie bootstrap p-value is unavailable, so the contradiction check could not be evaluated."
    ),
    source_note = dplyr::coalesce(xie_method_note, block5_source_note)
  ) %>%
  arrange(make_dataset_factor(dataset))

maller_zhou_contextual_summary <- legacy_screening_summary %>%
  inner_join(dataset_map_join, by = c("dataset_key" = "legacy_dataset_key")) %>%
  transmute(
    dataset,
    source_dataset = block5_source_dataset,
    analysis_branch = block5_analysis_branch,
    dataset_label = block5_dataset_label,
    d_n_statistic = as.numeric(maller_zhou_dn_stat),
    d_n_p_value = as.numeric(maller_zhou_dn_p_value),
    d_n_flag = as.character(maller_zhou_dn_flag),
    alpha_n = as.numeric(maller_zhou_alpha_n),
    q_n = as.numeric(maller_zhou_q_n),
    contextual_signal_flag = d_n_flag == "supportive",
    contextual_note = dplyr::case_when(
      d_n_flag == "supportive" ~ "Maller-Zhou shows a positive contextual tail signal, but Block 5 treats it as historical/contextual evidence only.",
      d_n_flag == "unsupportive" ~ "Maller-Zhou does not show a positive contextual tail signal and remains contextual only.",
      TRUE ~ "Maller-Zhou result is unavailable or unclassified and remains contextual only."
    ),
    source_note = block5_source_note
  ) %>%
  arrange(make_dataset_factor(dataset))

block5_screening_summary <- data_backbone_summary %>%
  select(dataset, source_dataset, analysis_branch, dataset_label, n_total, n_event, tau_year, plateau_length_year, km_tail_survival, km_tail_noncure, source_note) %>%
  left_join(
    receus_primary_summary %>%
      select(dataset, selected_family, pi_hat_n, r_hat_n, receus_primary_flag, supportive_candidate_count, supportive_candidate_families),
    by = "dataset"
  ) %>%
  left_join(
    xie_contradiction_summary %>%
      select(dataset, xie_epsilon, p_hat_n, p_hat_G, T_n, bootstrap_p_value, contradiction_flag),
    by = "dataset"
  ) %>%
  left_join(
    maller_zhou_contextual_summary %>%
      select(dataset, d_n_statistic, d_n_p_value, d_n_flag, alpha_n, q_n, contextual_signal_flag),
    by = "dataset"
  ) %>%
  mutate(
    block5_screening_class = mapply(
      make_block5_screening_class,
      receus_primary_flag,
      contradiction_flag,
      contextual_signal_flag
    ),
    block5_screening_note = mapply(
      make_block5_screening_note,
      receus_primary_flag,
      contradiction_flag,
      contextual_signal_flag
    ),
    interpretation_boundary = "Block 5 is a supporting interpretation block only; RECeUS is primary, Xie is contradiction-oriented, and Maller-Zhou is contextual only."
  ) %>%
  arrange(make_dataset_factor(dataset))

screening_overview_plot_data <- bind_rows(
  receus_primary_summary %>%
    transmute(dataset, component = "RECeUS-AIC primary", component_status = paste("Primary", receus_primary_flag)),
  xie_contradiction_summary %>%
    transmute(
      dataset,
      component = "Xie contradiction",
      component_status = ifelse(contradiction_flag, "Xie contradiction", "Xie no contradiction")
    ),
  maller_zhou_contextual_summary %>%
    transmute(
      dataset,
      component = "Maller-Zhou contextual",
      component_status = ifelse(contextual_signal_flag, "Contextual signal present", "Contextual signal absent")
    )
) %>%
  mutate(
    dataset = make_dataset_factor(dataset),
    component = factor(component, levels = c("RECeUS-AIC primary", "Xie contradiction", "Maller-Zhou contextual"))
  )

screening_overview_plot <- ggplot(
  screening_overview_plot_data,
  aes(x = component, y = dataset, fill = component_status)
) +
  geom_tile(color = "white", linewidth = 0.6) +
  scale_fill_manual(
    values = c(
      "Primary supportive" = "#54A24B",
      "Primary equivocal" = "#ECA82C",
      "Primary unsupportive" = "#E45756",
      "Xie contradiction" = "#E45756",
      "Xie no contradiction" = "#4C78A8",
      "Contextual signal present" = "#B279A2",
      "Contextual signal absent" = "#BDBDBD"
    )
  ) +
  labs(
    title = "Block 5 screening overview",
    subtitle = "Primary RECeUS-AIC, Xie contradiction check, and contextual Maller-Zhou signal",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

receus_candidate_plot <- receus_candidate_fits %>%
  mutate(dataset = make_dataset_factor(dataset)) %>%
  ggplot(aes(x = r_hat_n, y = pi_hat_n, color = receus_candidate_flag, shape = selected_by_aic)) +
  geom_vline(xintercept = receus_r_threshold, linewidth = 0.4, linetype = "dashed", color = "#666666") +
  geom_hline(yintercept = receus_pi_threshold, linewidth = 0.4, linetype = "dashed", color = "#666666") +
  geom_point(size = 3) +
  geom_text(aes(label = family), nudge_y = 0.02, size = 3, check_overlap = TRUE, show.legend = FALSE) +
  facet_wrap(~dataset, nrow = 1) +
  scale_color_manual(values = c("supportive" = "#54A24B", "equivocal" = "#ECA82C", "unsupportive" = "#E45756")) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17)) +
  labs(
    title = "RECeUS candidate-family map",
    subtitle = "Dashed lines show the paper's default thresholds pi_hat_n > 0.025 and r_hat_n < 0.05",
    x = "r_hat_n",
    y = "pi_hat_n",
    color = "RECeUS flag",
    shape = "AIC-selected"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

# Transform: Xie bootstrap reuse =================================================
xie_bootstrap_results <- NULL
xie_bootstrap_plot_index <- NULL
xie_bootstrap_pdf_path <- file.path(bootstrap_detail_dir, "block5_xie_bootstrap_histograms.pdf")

if (!is.null(legacy_bootstrap_results)) {
  xie_bootstrap_results <- legacy_bootstrap_results %>%
    inner_join(dataset_map_join, by = c("dataset_key" = "legacy_dataset_key")) %>%
    filter(task_type == "xie", record_type == "distribution") %>%
    transmute(
      dataset,
      source_dataset = block5_source_dataset,
      analysis_branch = block5_analysis_branch,
      dataset_label = block5_dataset_label,
      rep_id = as.integer(rep_id),
      bootstrap_statistic = as.numeric(plot_statistic_value),
      observed_statistic = as.numeric(plot_observed_line),
      critical_value = as.numeric(plot_critical_line),
      bootstrap_p_value = as.numeric(p_value),
      n_bootstrap = as.integer(n_bootstrap),
      n_valid = as.integer(n_valid),
      valid_fraction = as.numeric(valid_fraction),
      tau_year = as.numeric(tau_year),
      source_note = "Rebuilt from the legacy Stage 6 bootstrap distribution table without rerunning the bootstrap."
    ) %>%
    arrange(make_dataset_factor(dataset), rep_id)

  xie_bootstrap_plot_index <- xie_bootstrap_results %>%
    group_by(dataset, source_dataset, analysis_branch, dataset_label) %>%
    reframe(
      observed_statistic = observed_statistic[[1]],
      critical_value = critical_value[[1]],
      bootstrap_p_value = bootstrap_p_value[[1]],
      n_bootstrap = n_bootstrap[[1]],
      n_valid = n_valid[[1]],
      valid_fraction = valid_fraction[[1]],
      tau_year = tau_year[[1]],
      png_file_name = paste0("block5_xie_bootstrap_histogram__", dataset[[1]], ".png")
    ) %>%
    arrange(make_dataset_factor(dataset))
}

if (!is.null(xie_bootstrap_results) && isTRUE(write_xie_bootstrap_histograms)) {
  xie_plot_list <- list()
  for (ii in seq_len(nrow(xie_bootstrap_plot_index))) {
    current_index <- xie_bootstrap_plot_index[ii, , drop = FALSE]
    current_df <- xie_bootstrap_results %>%
      filter(dataset == current_index$dataset[[1]])
    p <- ggplot(current_df, aes(x = bootstrap_statistic)) +
      geom_histogram(bins = 35, fill = "#4C78A8", color = "white", linewidth = 0.3) +
      geom_vline(xintercept = current_index$observed_statistic[[1]], color = "#E45756", linewidth = 0.8) +
      geom_vline(xintercept = current_index$critical_value[[1]], color = "#1F1F1F", linewidth = 0.8, linetype = "dashed") +
      labs(
        title = paste0("Xie centered bootstrap distribution: ", current_index$dataset_label[[1]]),
        subtitle = paste0(
          "Observed T_n = ", signif(current_index$observed_statistic[[1]], 4),
          " | critical value = ", signif(current_index$critical_value[[1]], 4),
          " | p = ", signif(current_index$bootstrap_p_value[[1]], 4)
        ),
        x = "Centered bootstrap statistic",
        y = "Count"
      ) +
      theme_minimal(base_size = 12)
    xie_plot_list[[current_index$dataset[[1]]]] <- p
    safe_save_plot(
      p,
      file.path(bootstrap_detail_dir, current_index$png_file_name[[1]]),
      width = 9,
      height = 6
    )
  }

  safe_save_pdf_atomic(
    xie_bootstrap_pdf_path,
    width = 9,
    height = 6,
    plot_fun = function() {
      for (nm in names(xie_plot_list)) {
        print(xie_plot_list[[nm]])
      }
    }
  )
}

# Export: Block 5 outputs ========================================================
main_output_paths <- list(
  screening_summary = file.path(main_output_dir, "block5_screening_summary.csv"),
  screening_overview_plot = file.path(main_output_dir, "block5_screening_overview.png"),
  receus_candidate_map = file.path(main_output_dir, "block5_receus_candidate_map.png"),
  visual_summary_pdf = file.path(main_output_dir, "block5_visual_summary.pdf")
)
supporting_output_paths <- list(
  data_backbone_summary = file.path(supporting_table_dir, "block5_data_backbone_summary.csv"),
  receus_primary_summary = file.path(supporting_table_dir, "block5_receus_primary_summary.csv"),
  receus_candidate_fits = file.path(supporting_table_dir, "block5_receus_candidate_fits.csv"),
  xie_contradiction_summary = file.path(supporting_table_dir, "block5_xie_contradiction_summary.csv"),
  maller_zhou_contextual_summary = file.path(supporting_table_dir, "block5_maller_zhou_contextual_summary.csv"),
  screening_overview_plot_data = file.path(supporting_table_dir, "block5_screening_overview_plot_data.csv"),
  standard_error_registry = file.path(supporting_table_dir, "block5_standard_error_table_registry.csv"),
  standard_error_long = file.path(supporting_table_dir, "block5_standard_error_long.csv")
)
bootstrap_output_paths <- list(
  xie_bootstrap_results = file.path(bootstrap_detail_dir, "block5_xie_bootstrap_results.csv"),
  xie_bootstrap_plot_index = file.path(bootstrap_detail_dir, "block5_xie_bootstrap_plot_index.csv"),
  xie_bootstrap_histograms_pdf = xie_bootstrap_pdf_path
)
metadata_output_paths <- list(
  reference_registry = file.path(metadata_dir, "block5_reference_registry.csv"),
  export_manifest = file.path(metadata_dir, "block5_export_manifest.csv")
)

block5_standard_error_bundle <- build_standard_error_export_bundle(list(
  data_backbone_summary = data_backbone_summary,
  receus_primary_summary = receus_primary_summary,
  receus_candidate_fits = receus_candidate_fits,
  xie_contradiction_summary = xie_contradiction_summary,
  maller_zhou_contextual_summary = maller_zhou_contextual_summary,
  screening_summary = block5_screening_summary,
  screening_overview_plot_data = screening_overview_plot_data,
  xie_bootstrap_results = xie_bootstrap_results,
  xie_bootstrap_plot_index = xie_bootstrap_plot_index
))

safe_write_csv_atomic(data_backbone_summary, supporting_output_paths$data_backbone_summary)
safe_write_csv_atomic(receus_primary_summary, supporting_output_paths$receus_primary_summary)
safe_write_csv_atomic(receus_candidate_fits, supporting_output_paths$receus_candidate_fits)
safe_write_csv_atomic(xie_contradiction_summary, supporting_output_paths$xie_contradiction_summary)
safe_write_csv_atomic(maller_zhou_contextual_summary, supporting_output_paths$maller_zhou_contextual_summary)
safe_write_csv_atomic(block5_screening_summary, main_output_paths$screening_summary)
safe_write_csv_atomic(reference_registry, metadata_output_paths$reference_registry)
safe_write_csv_atomic(screening_overview_plot_data, supporting_output_paths$screening_overview_plot_data)
safe_write_csv_atomic(block5_standard_error_bundle$registry, supporting_output_paths$standard_error_registry)
safe_write_csv_atomic(block5_standard_error_bundle$long, supporting_output_paths$standard_error_long)

safe_save_plot(screening_overview_plot, main_output_paths$screening_overview_plot, width = 10, height = 4.8)
safe_save_plot(receus_candidate_plot, main_output_paths$receus_candidate_map, width = 12, height = 4.8)
safe_save_pdf_atomic(
  main_output_paths$visual_summary_pdf,
  width = 12,
  height = 5,
  plot_fun = function() {
    print(screening_overview_plot)
    print(receus_candidate_plot)
  }
)

if (!is.null(xie_bootstrap_results)) {
  safe_write_csv_atomic(xie_bootstrap_results, bootstrap_output_paths$xie_bootstrap_results)
}
if (!is.null(xie_bootstrap_plot_index)) {
  safe_write_csv_atomic(xie_bootstrap_plot_index, bootstrap_output_paths$xie_bootstrap_plot_index)
}

export_manifest <- tibble::tribble(
  ~artifact_name, ~description, ~file_path,
  "data_backbone_summary", "Dataset-level Block 5 tail and follow-up summary restricted to the new three-dataset workflow.", supporting_output_paths$data_backbone_summary,
  "receus_primary_summary", "Primary RECeUS-AIC Block 5 summary using the paper's default threshold rule.", supporting_output_paths$receus_primary_summary,
  "receus_candidate_fits", "Candidate-family RECeUS fits for the AIC model class described in the RECeUS paper.", supporting_output_paths$receus_candidate_fits,
  "xie_contradiction_summary", "Secondary contradiction-oriented Xie sufficient-follow-up summary with bootstrap p-values.", supporting_output_paths$xie_contradiction_summary,
  "maller_zhou_contextual_summary", "Older Maller-Zhou contextual reference summary.", supporting_output_paths$maller_zhou_contextual_summary,
  "screening_summary", "Integrated Block 5 summary that preserves the RECeUS primary, Xie contradiction, Maller-Zhou contextual hierarchy.", main_output_paths$screening_summary,
  "reference_registry", "Reference registry linking each Block 5 method to the paper used for the implementation boundary.", metadata_output_paths$reference_registry,
  "screening_overview_plot_data", "Long plot-ready table for the Block 5 overview tile chart.", supporting_output_paths$screening_overview_plot_data,
  "standard_error_registry", "Registry of Block 5 tables that contain standard-error or uncertainty columns.", supporting_output_paths$standard_error_registry,
  "standard_error_long", "Long-form Block 5 standard-error table assembled from the exported screening sources.", supporting_output_paths$standard_error_long,
  "screening_overview_plot", "Overview plot for the spec-aligned Block 5 hierarchy.", main_output_paths$screening_overview_plot,
  "receus_candidate_map", "Scatterplot of RECeUS candidate families against the paper's thresholds.", main_output_paths$receus_candidate_map,
  "visual_summary_pdf", "Multi-page PDF bundling the main Block 5 visual summaries.", main_output_paths$visual_summary_pdf
) %>%
  mutate(
    file_name = basename(file_path),
    source_mode = "rebuilt_from_legacy_stage6_outputs_without_rerunning_heavy_model_fits"
  )

if (!is.null(xie_bootstrap_results)) {
  export_manifest <- bind_rows(
    export_manifest,
    tibble(
      artifact_name = "xie_bootstrap_results",
      description = "Filtered Xie bootstrap distribution rebuilt from the legacy Stage 6 bootstrap table.",
      file_path = bootstrap_output_paths$xie_bootstrap_results,
      file_name = basename(bootstrap_output_paths$xie_bootstrap_results),
      source_mode = "reused_legacy_xie_bootstrap_distribution"
    )
  )
}
if (!is.null(xie_bootstrap_plot_index)) {
  export_manifest <- bind_rows(
    export_manifest,
    tibble(
      artifact_name = "xie_bootstrap_plot_index",
      description = "Index of regenerated Xie bootstrap histogram exports.",
      file_path = bootstrap_output_paths$xie_bootstrap_plot_index,
      file_name = basename(bootstrap_output_paths$xie_bootstrap_plot_index),
      source_mode = "reused_legacy_xie_bootstrap_distribution"
    )
  )
}
if (file.exists(xie_bootstrap_pdf_path)) {
  export_manifest <- bind_rows(
    export_manifest,
    tibble(
      artifact_name = "xie_bootstrap_histograms_pdf",
      description = "Multi-page PDF of regenerated Xie bootstrap histograms only.",
      file_path = xie_bootstrap_pdf_path,
      file_name = basename(xie_bootstrap_pdf_path),
      source_mode = "reused_legacy_xie_bootstrap_distribution"
    )
  )
}

safe_write_csv_atomic(export_manifest, metadata_output_paths$export_manifest)

for (target_path in c(
  unlist(supporting_output_paths, use.names = FALSE),
  unlist(bootstrap_output_paths, use.names = FALSE),
  unlist(metadata_output_paths, use.names = FALSE)
)) {
  cleanup_duplicate_root_file(export_path, target_path)
}

message("Block 5 export completed:")
for (ii in seq_len(nrow(export_manifest))) {
  message("  - ", export_manifest$file_path[[ii]])
}
