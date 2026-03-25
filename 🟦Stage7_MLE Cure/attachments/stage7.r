# 🔴 Configure: runtime paths and rerun-safe switches ===============================
run_root_dir <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/STAGE7/stage7_run__2026-03-24_shard50_auto20_v1__00009fe20000'
export_path <- run_root_dir

stage1_data_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.데이터 분석 New/🟦Stage1/attachments'
stage1_analysis_datasets_file <- file.path(stage1_data_path, "stage1_analysis_datasets.rds")

stage6_screening_file <- NA_character_
# Supported alternatives:
# - '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.데이터 분석 New/🟦Stage6/attachments/stage6_carry_forward_flag_table.csv'
# - '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.데이터 분석 New/🟦Stage6/attachments/stage6_screening_bundle.rds'
# The code will extract outputs$carry_forward_flag_table from the RDS bundle when needed.

refresh_exports_even_when_outputs_exist <- TRUE
run_bootstrap_extraction <- TRUE
write_full_metric_values_long_csv <- FALSE

top_n_variable_metric_cells <- 20L
top_n_failure_patterns <- 30L
top_n_failure_models <- 30L

pnu_threshold_suppress_from_year <- 3L
common_horizons_year <- 1:10
threshold_plot_horizons_year <- c(1L, 2L, 5L, 10L)

# 🔴 Initialize: packages and global options ===============================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(survival)
  library(stringr)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# 🔴 Define: reusable file helpers ===============================
## 🟠 Define: path and I/O guards ===============================
assert_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

normalize_existing_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

safe_read_csv <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

safe_write_csv <- function(df, path) {
  readr::write_csv(df, path)
  normalize_existing_path(path)
}

safe_save_rds <- function(object, path) {
  saveRDS(object, file = path)
  normalize_existing_path(path)
}

make_output_path <- function(file_name) {
  file.path(export_path, file_name)
}

# 🔴 Define: scalar and table helpers ===============================
## 🟠 Define: low-level coercion helpers ===============================
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

coalesce_character <- function(...) {
  xs <- list(...)
  if (length(xs) == 0L) {
    return(character(0))
  }
  out <- as.character(xs[[1]])
  if (length(xs) > 1L) {
    for (ii in 2:length(xs)) {
      cand <- as.character(xs[[ii]])
      idx <- is.na(out) | out == ""
      out[idx] <- cand[idx]
    }
  }
  out
}

col_or_chr <- function(df, nm, default = NA_character_) {
  if (!nm %in% names(df)) return(rep(default, nrow(df)))
  as.character(df[[nm]])
}

col_or_num <- function(df, nm, default = NA_real_) {
  if (!nm %in% names(df)) return(rep(default, nrow(df)))
  suppressWarnings(as.numeric(df[[nm]]))
}

col_or_lgl <- function(df, nm, default = NA) {
  if (!nm %in% names(df)) return(rep(default, nrow(df)))
  as.logical(df[[nm]])
}

safe_quantile <- function(x, prob) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, names = FALSE, na.rm = TRUE))
}

collapse_unique_sorted <- function(x, empty_value = "") {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) return(empty_value)
  paste(sort(x), collapse = "|")
}

is_success_status <- function(x) {
  stringr::str_starts(coalesce(as.character(x), ""), "ok")
}

make_metric_key <- function(df) {
  paste0(
    df$table_name,
    " | ", ifelse(is.na(df$dataset) | df$dataset == "", "<no_dataset>", df$dataset),
    " | ", ifelse(is.na(df$formula_variant) | df$formula_variant == "", "<no_formula>", df$formula_variant),
    " | ", ifelse(is.na(df$model_id) | df$model_id == "", "<no_model>", df$model_id),
    ifelse(!is.na(df$component) & nzchar(df$component), paste0(" | ", df$component), ""),
    ifelse(!is.na(df$term) & nzchar(df$term), paste0(" | ", df$term), ""),
    ifelse(!is.na(df$horizon_year), paste0(" | h=", df$horizon_year), ""),
    ifelse(!is.na(df$threshold), paste0(" | thr=", format(df$threshold, trim = TRUE, scientific = FALSE)), ""),
    " | ", df$metric_name
  )
}

# 🔴 Define: Stage 6 carry-forward readers ===============================
## 🟠 Define: carry-forward file parsers ===============================
make_empty_stage6_carry_forward_tbl <- function() {
  tibble(
    dataset_key = character(),
    source_dataset = character(),
    analysis_variant = character(),
    final_decision_flag = character(),
    carry_forward_stage7 = logical(),
    carry_forward_stage8 = logical(),
    screening_repeated_in_stage8 = logical(),
    presence_evidence_label = character(),
    followup_evidence_label = character(),
    screening_note = character()
  )
}

read_stage6_carry_forward <- function(path) {
  assert_file_exists(path, "Stage 6 screening file")
  ext <- tolower(tools::file_ext(path))
  
  if (ext == "csv") {
    tbl <- safe_read_csv(path)
  } else if (ext == "rds") {
    obj <- readRDS(path)
    if (is.data.frame(obj)) {
      tbl <- tibble::as_tibble(obj)
    } else if (is.list(obj) && !is.null(obj$outputs) && !is.null(obj$outputs$carry_forward_flag_table)) {
      tbl <- tibble::as_tibble(obj$outputs$carry_forward_flag_table)
    } else {
      stop("Unsupported Stage 6 RDS structure. Expected a data frame or screening_bundle with outputs$carry_forward_flag_table.", call. = FALSE)
    }
  } else {
    stop("Stage 6 screening file must be either CSV or RDS.", call. = FALSE)
  }
  
  required_cols <- c(
    "dataset_key", "source_dataset", "analysis_variant", "final_decision_flag",
    "carry_forward_stage7", "carry_forward_stage8", "screening_repeated_in_stage8"
  )
  missing_cols <- setdiff(required_cols, names(tbl))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Stage 6 carry-forward table is missing required columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  
  tbl %>%
    mutate(
      dataset_key = as.character(dataset_key),
      source_dataset = as.character(source_dataset),
      analysis_variant = as.character(analysis_variant)
    )
}

load_stage6_carry_forward_optional <- function(path = NA_character_) {
  if (is.null(path) || length(path) == 0L || is.na(path) || !nzchar(path)) {
    return(make_empty_stage6_carry_forward_tbl())
  }
  
  if (!file.exists(path)) {
    warning(
      sprintf("Stage 6 screening file not found: %s. Proceeding without Stage 6 carry-forward.", path),
      call. = FALSE
    )
    return(make_empty_stage6_carry_forward_tbl())
  }
  
  read_stage6_carry_forward(path)
}

make_stage6_join_key <- function(dataset, formula_variant = NA_character_, site_branch = NA_character_) {
  dataset_chr <- as.character(dataset)
  formula_chr <- as.character(formula_variant)
  site_branch_chr <- as.character(site_branch)
  
  dplyr::case_when(
    dataset_chr %in% c("PNU", "SNU") ~ dataset_chr,
    dataset_chr == "merged" & site_branch_chr == "site_adjusted" ~ "merged__site_adjusted",
    dataset_chr == "merged" ~ "merged__site_free",
    TRUE ~ NA_character_
  )
}

build_stage6_registry <- function(stage6_tbl) {
  stage6_tbl %>%
    mutate(
      stage6_join_key = dataset_key,
      stage6_flag_joined = TRUE
    ) %>%
    rename_with(
      .fn = ~ ifelse(.x %in% c("stage6_join_key", "stage6_flag_joined"), .x, paste0("stage6__", .x))
    )
}

join_stage6_flags <- function(df, stage6_registry) {
  if (nrow(df) == 0L) {
    return(df)
  }
  
  existing_stage6_cols <- grep("^stage6_", names(df), value = TRUE)
  
  join_key <- make_stage6_join_key(
    dataset = col_or_chr(df, "dataset"),
    formula_variant = col_or_chr(df, "formula_variant"),
    site_branch = col_or_chr(df, "site_branch")
  )
  
  out <- df %>%
    select(-any_of(existing_stage6_cols)) %>%
    mutate(stage6_join_key = join_key) %>%
    left_join(stage6_registry, by = "stage6_join_key") %>%
    mutate(
      stage6_flag_joined = dplyr::coalesce(stage6_flag_joined, FALSE),
      stage6__screening_note = dplyr::coalesce(
        stage6__screening_note,
        if_else(stage6_flag_joined, NA_character_, "Stage 6 carry-forward was not joined in this refresh run.")
      )
    )
  
  out
}

# 🔴 Define: main Stage 7 artifact readers ===============================
## 🟠 Define: current output loaders ===============================
load_stage7_outputs <- function(root_dir) {
  required_files <- c(
    "stage7_run_manifest.csv",
    "stage7_fit_registry.csv",
    "stage7_coefficients.csv",
    "stage7_subject_predictions.csv",
    "stage7_risk_summary.csv",
    "stage7_delta_risk.csv",
    "stage7_threshold_metrics.csv"
  )
  
  for (file_name in required_files) {
    assert_file_exists(file.path(root_dir, file_name), file_name)
  }
  
  list(
    run_manifest = safe_read_csv(file.path(root_dir, "stage7_run_manifest.csv")),
    fit_registry = safe_read_csv(file.path(root_dir, "stage7_fit_registry.csv")),
    coefficients = safe_read_csv(file.path(root_dir, "stage7_coefficients.csv")),
    subject_predictions = safe_read_csv(file.path(root_dir, "stage7_subject_predictions.csv")),
    risk_summary = safe_read_csv(file.path(root_dir, "stage7_risk_summary.csv")),
    delta_risk = safe_read_csv(file.path(root_dir, "stage7_delta_risk.csv")),
    threshold_metrics = safe_read_csv(file.path(root_dir, "stage7_threshold_metrics.csv")),
    bootstrap_model_qc = if (file.exists(file.path(root_dir, "stage7_bootstrap_model_qc.csv"))) {
      safe_read_csv(file.path(root_dir, "stage7_bootstrap_model_qc.csv"))
    } else {
      tibble()
    },
    shard_registry = if (file.exists(file.path(root_dir, "stage7_shard_registry.csv"))) {
      safe_read_csv(file.path(root_dir, "stage7_shard_registry.csv"))
    } else {
      tibble()
    },
    merged_bootstrap_rds_file = file.path(root_dir, "stage7_bootstrap_merged_results.rds")
  )
}

# 🔴 Define: IPCW rebuilding logic ===============================
## 🟠 Define: IPCW reconstruction from Stage 1 backbone ===============================
compute_ipcw_frame <- function(df, horizons_year) {
  censor_fit <- survival::survfit(survival::Surv(time_year, censor_main) ~ 1, data = df)
  
  step_value <- function(x, y, at, yleft = 1, yright = NULL, lower = 0, upper = 1) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    at <- as.numeric(at)
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
    if (length(x) == 0L || length(y) == 0L) {
      out <- rep(yleft, length(at))
      return(pmin(pmax(out, lower), upper))
    }
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    if (is.null(yright)) {
      yright <- tail(y, 1)
    }
    out <- approx(
      x = x,
      y = y,
      xout = at,
      method = "constant",
      f = 0,
      yleft = yleft,
      yright = yright,
      ties = "ordered"
    )$y
    pmin(pmax(out, lower), upper)
  }
  
  G_at_t <- step_value(censor_fit$time, censor_fit$surv, horizons_year, yleft = 1, lower = 0, upper = 1)
  G_at_obs <- step_value(censor_fit$time, censor_fit$surv, df$time_year, yleft = 1, lower = 0, upper = 1)
  
  bind_rows(lapply(seq_along(horizons_year), function(idx) {
    hh <- as.numeric(horizons_year[idx])
    G_h <- pmax(G_at_t[idx], 1e-10)
    is_case <- df$event_main == 1L & df$time_year <= hh
    is_control <- df$time_year > hh
    
    tibble(
      dataset = unique(as.character(df$dataset)),
      unique_person_id = df$unique_person_id,
      horizon_year = as.integer(hh),
      event_by_horizon = as.integer(is_case),
      known_case_or_control = as.integer(is_case | is_control),
      weight_case = if_else(is_case, 1 / pmax(G_at_obs, 1e-10), 0),
      weight_control = if_else(is_control, 1 / G_h, 0)
    )
  }))
}

rebuild_ipcw_registry <- function(stage1_analysis_datasets_file, horizons_year) {
  assert_file_exists(stage1_analysis_datasets_file, "stage1_analysis_datasets.rds")
  dataset_map <- readRDS(stage1_analysis_datasets_file)
  required_names <- c("PNU", "SNU", "merged")
  if (!is.list(dataset_map) || !all(required_names %in% names(dataset_map))) {
    stop("stage1_analysis_datasets.rds must contain PNU, SNU, and merged.", call. = FALSE)
  }
  bind_rows(lapply(required_names, function(nm) {
    df <- tibble::as_tibble(dataset_map[[nm]])
    if (!"dataset" %in% names(df)) {
      df$dataset <- nm
    }
    compute_ipcw_frame(df, horizons_year = horizons_year)
  }))
}

# 🔴 Define: threshold suppression rules ===============================
## 🟠 Define: PNU late-horizon masking ===============================
apply_threshold_reporting_rules <- function(threshold_metrics_tbl) {
  metric_cols_to_blank <- c(
    "positive_classification_rate",
    "weighted_tp",
    "weighted_fn",
    "weighted_fp",
    "weighted_tn",
    "sensitivity",
    "specificity",
    "ppv",
    "false_positive_weighted_per_n",
    "false_positive_burden_non_event",
    "false_positive_burden_primary",
    "false_positive_burden_all",
    "unnecessary_high_risk_per_100_population",
    "unnecessary_high_risk_per_100_non_event",
    "net_benefit",
    "false_positive_burden_primary_ci_lower",
    "false_positive_burden_primary_ci_upper",
    "false_positive_weighted_per_n_ci_lower",
    "false_positive_weighted_per_n_ci_upper",
    "net_benefit_ci_lower",
    "net_benefit_ci_upper",
    "false_positive_burden_ci_lower",
    "false_positive_burden_ci_upper"
  )
  
  existing_metric_cols <- intersect(metric_cols_to_blank, names(threshold_metrics_tbl))
  
  suppression_flag <- threshold_metrics_tbl$dataset == "PNU" & as.integer(threshold_metrics_tbl$horizon_year) >= as.integer(pnu_threshold_suppress_from_year)
  
  out <- threshold_metrics_tbl %>%
    mutate(
      threshold_projection_suppressed_flag = suppression_flag,
      threshold_reporting_status = if_else(threshold_projection_suppressed_flag, "suppressed_projection", "reported"),
      threshold_suppression_reason = if_else(
        threshold_projection_suppressed_flag,
        "PNU threshold-based metrics suppressed from year 3 onward because they are projection-only and not primary-supported.",
        NA_character_
      )
    )
  
  if (length(existing_metric_cols) > 0L) {
    out[existing_metric_cols] <- lapply(out[existing_metric_cols], function(col) {
      col[suppression_flag] <- NA
      col
    })
  }
  
  out
}

# 🔴 Define: bootstrap extraction readers ===============================
## 🟠 Define: bootstrap object normalization ===============================
required_bootstrap_tables <- c(
  "fit_registry",
  "coefficients",
  "risk_summary",
  "delta_risk",
  "threshold_metrics"
)

standardize_bootstrap_results <- function(x) {
  if (!is.list(x)) {
    stop("`stage7_bootstrap_merged_results.rds` must contain a list.", call. = FALSE)
  }
  for (nm in required_bootstrap_tables) {
    if (!nm %in% names(x) || is.null(x[[nm]])) {
      x[[nm]] <- tibble()
    }
    x[[nm]] <- tibble::as_tibble(x[[nm]])
  }
  x[required_bootstrap_tables]
}

# 🔴 Define: bootstrap failure summaries ===============================
## 🟠 Define: replicate-level failure extractors ===============================
normalize_fit_registry <- function(fit_registry_raw) {
  if (nrow(fit_registry_raw) == 0L) {
    return(
      tibble(
        rep_id = integer(),
        dataset = character(),
        formula_variant = character(),
        site_branch = character(),
        interaction_branch = character(),
        model_id = character(),
        model_class = character(),
        model_block = character(),
        model_family = character(),
        latency_type = character(),
        matched_nocure_model_id = character(),
        fit_status = character(),
        prediction_status = character(),
        overall_status = character(),
        fit_success = logical(),
        prediction_success = logical(),
        overall_success = logical(),
        failure_type = character(),
        convergence_code = double(),
        warning_count_total = double(),
        primary_error_message = character(),
        primary_warning_message = character(),
        model_key = character(),
        failure_key = character()
      )
    )
  }
  
  fit_success_fallback <- is_success_status(col_or_chr(fit_registry_raw, "fit_status"))
  prediction_success_fallback <- is_success_status(col_or_chr(fit_registry_raw, "prediction_status"))
  
  fit_success_vec <- dplyr::coalesce(
    col_or_lgl(fit_registry_raw, "fit_component_success"),
    fit_success_fallback
  )
  
  prediction_success_vec <- dplyr::coalesce(
    col_or_lgl(fit_registry_raw, "prediction_component_success"),
    prediction_success_fallback
  )
  
  warning_count_total_vec <- dplyr::coalesce(
    col_or_num(fit_registry_raw, "warning_count_total"),
    col_or_num(fit_registry_raw, "fit_warning_count") + col_or_num(fit_registry_raw, "prediction_warning_count"),
    col_or_num(fit_registry_raw, "fit_warning_count"),
    col_or_num(fit_registry_raw, "prediction_warning_count"),
    rep(0, nrow(fit_registry_raw))
  )
  
  tibble(
    rep_id = as.integer(col_or_num(fit_registry_raw, "rep_id")),
    dataset = col_or_chr(fit_registry_raw, "dataset", "<missing_dataset>"),
    formula_variant = col_or_chr(fit_registry_raw, "formula_variant", "<missing_formula>"),
    site_branch = col_or_chr(fit_registry_raw, "site_branch"),
    interaction_branch = col_or_chr(fit_registry_raw, "interaction_branch"),
    model_id = col_or_chr(fit_registry_raw, "model_id", "<missing_model>"),
    model_class = col_or_chr(fit_registry_raw, "model_class"),
    model_block = col_or_chr(fit_registry_raw, "model_block"),
    model_family = col_or_chr(fit_registry_raw, "model_family"),
    latency_type = col_or_chr(fit_registry_raw, "latency_type"),
    matched_nocure_model_id = col_or_chr(fit_registry_raw, "matched_nocure_model_id"),
    fit_status = col_or_chr(fit_registry_raw, "fit_status"),
    prediction_status = col_or_chr(fit_registry_raw, "prediction_status"),
    overall_status = col_or_chr(fit_registry_raw, "overall_status"),
    fit_success = as.logical(fit_success_vec),
    prediction_success = as.logical(prediction_success_vec),
    convergence_code = col_or_num(fit_registry_raw, "convergence_code"),
    warning_count_total = as.numeric(warning_count_total_vec),
    primary_error_message = coalesce_character(
      col_or_chr(fit_registry_raw, "bootstrap_error"),
      col_or_chr(fit_registry_raw, "error_message"),
      col_or_chr(fit_registry_raw, "fit_error_message"),
      col_or_chr(fit_registry_raw, "prediction_error_message"),
      col_or_chr(fit_registry_raw, "fit_message"),
      col_or_chr(fit_registry_raw, "prediction_message"),
      col_or_chr(fit_registry_raw, "message")
    ),
    primary_warning_message = coalesce_character(
      col_or_chr(fit_registry_raw, "warning_message"),
      col_or_chr(fit_registry_raw, "fit_warning_message"),
      col_or_chr(fit_registry_raw, "prediction_warning_message")
    )
  ) %>%
    mutate(
      overall_success = fit_success & prediction_success,
      failure_type = dplyr::case_when(
        !fit_success & !prediction_success ~ "fit_and_prediction_failed",
        !fit_success ~ "fit_failed",
        fit_success & !prediction_success ~ "prediction_failed",
        TRUE ~ "success"
      ),
      model_key = paste(dataset, formula_variant, model_id, sep = "__"),
      failure_key = ifelse(
        failure_type == "success",
        NA_character_,
        paste(dataset, formula_variant, model_id, failure_type, sep = "::")
      )
    )
}

make_failure_events <- function(fit_registry_norm) {
  fit_registry_norm %>%
    filter(!overall_success) %>%
    transmute(
      rep_id,
      dataset,
      formula_variant,
      site_branch,
      interaction_branch,
      model_id,
      model_class,
      model_block,
      model_family,
      latency_type,
      matched_nocure_model_id,
      failure_type,
      fit_status,
      prediction_status,
      overall_status,
      convergence_code,
      warning_count_total,
      primary_error_message,
      primary_warning_message
    ) %>%
    arrange(rep_id, dataset, formula_variant, model_id)
}

make_replicate_overview <- function(fit_registry_norm) {
  fit_registry_norm %>%
    group_by(rep_id) %>%
    summarise(
      n_model_rows = n(),
      n_fail_total = sum(!overall_success, na.rm = TRUE),
      n_fit_fail = sum(!fit_success, na.rm = TRUE),
      n_prediction_fail = sum(fit_success & !prediction_success, na.rm = TRUE),
      n_warning_rows = sum(warning_count_total > 0, na.rm = TRUE),
      fit_success_rate = mean(fit_success, na.rm = TRUE),
      prediction_success_rate = mean(prediction_success, na.rm = TRUE),
      overall_success_rate = mean(overall_success, na.rm = TRUE),
      failed_dataset_set = collapse_unique_sorted(dataset[!overall_success], empty_value = ""),
      failed_model_set = collapse_unique_sorted(model_id[!overall_success], empty_value = ""),
      failure_signature = ifelse(
        sum(!overall_success, na.rm = TRUE) == 0,
        "NO_FAILURE",
        paste(sort(unique(failure_key[!is.na(failure_key)])), collapse = " | ")
      ),
      pnu_fail_n = sum(dataset == "PNU" & !overall_success, na.rm = TRUE),
      snu_fail_n = sum(dataset == "SNU" & !overall_success, na.rm = TRUE),
      merged_fail_n = sum(dataset == "merged" & !overall_success, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(rep_id)
}

make_failure_pattern_summary <- function(replicate_overview) {
  replicate_overview %>%
    group_by(failure_signature) %>%
    summarise(
      n_replicates = n(),
      n_replicates_with_failures = sum(n_fail_total > 0, na.rm = TRUE),
      mean_failures_per_replicate = mean(n_fail_total, na.rm = TRUE),
      min_rep_id = min(rep_id, na.rm = TRUE),
      max_rep_id = max(rep_id, na.rm = TRUE),
      example_rep_ids = paste(head(rep_id, 20), collapse = "|"),
      .groups = "drop"
    ) %>%
    arrange(desc(n_replicates), desc(mean_failures_per_replicate), failure_signature)
}

make_failure_by_model <- function(failure_events) {
  if (nrow(failure_events) == 0L) {
    return(
      tibble(
        dataset = character(),
        formula_variant = character(),
        model_id = character(),
        failure_type = character(),
        n_failure_rows = integer(),
        n_replicates = integer(),
        example_rep_ids = character()
      )
    )
  }
  
  failure_events %>%
    group_by(dataset, formula_variant, model_id, failure_type) %>%
    summarise(
      n_failure_rows = n(),
      n_replicates = n_distinct(rep_id),
      example_rep_ids = paste(head(sort(unique(rep_id)), 20), collapse = "|"),
      .groups = "drop"
    ) %>%
    arrange(desc(n_replicates), desc(n_failure_rows), dataset, formula_variant, model_id)
}

# 🔴 Define: metric summarization without giant long pivot ===============================
## 🟠 Define: memory-safe metric summaries ===============================
metric_id_cols <- c(
  "dataset", "formula_variant", "site_branch", "interaction_branch",
  "model_id", "model_class", "model_block", "model_family", "latency_type",
  "matched_nocure_model_id", "component", "term", "horizon_year", "threshold"
)

numeric_metric_cols <- function(df, structural_cols = character(), extra_exclude = character()) {
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  out <- setdiff(numeric_cols, structural_cols)
  out <- setdiff(out, extra_exclude)
  out
}

complete_metric_id_columns <- function(df) {
  default_values <- list(
    dataset = NA_character_,
    formula_variant = NA_character_,
    site_branch = NA_character_,
    interaction_branch = NA_character_,
    model_id = NA_character_,
    model_class = NA_character_,
    model_block = NA_character_,
    model_family = NA_character_,
    latency_type = NA_character_,
    matched_nocure_model_id = NA_character_,
    component = NA_character_,
    term = NA_character_,
    horizon_year = NA_integer_,
    threshold = NA_real_
  )
  
  for (nm in names(default_values)) {
    if (!nm %in% names(df)) {
      df[[nm]] <- default_values[[nm]]
    }
  }
  
  df
}

get_metric_specs <- function(bootstrap_results) {
  list(
    coefficients = list(
      df = bootstrap_results$coefficients,
      structural_cols = c("rep_id"),
      metric_cols = intersect("estimate", names(bootstrap_results$coefficients))
    ),
    risk_summary = list(
      df = bootstrap_results$risk_summary,
      structural_cols = c("rep_id", "horizon_year"),
      metric_cols = numeric_metric_cols(
        bootstrap_results$risk_summary,
        structural_cols = c("rep_id", "horizon_year"),
        extra_exclude = c("n_subject")
      )
    ),
    delta_risk = list(
      df = bootstrap_results$delta_risk,
      structural_cols = c("rep_id", "horizon_year"),
      metric_cols = numeric_metric_cols(
        bootstrap_results$delta_risk,
        structural_cols = c("rep_id", "horizon_year"),
        extra_exclude = c("n_subject")
      )
    ),
    threshold_metrics = list(
      df = bootstrap_results$threshold_metrics,
      structural_cols = c("rep_id", "horizon_year", "threshold"),
      metric_cols = numeric_metric_cols(
        bootstrap_results$threshold_metrics,
        structural_cols = c("rep_id", "horizon_year", "threshold"),
        extra_exclude = c("n_subject")
      )
    )
  )
}

summarize_one_metric_column <- function(df, table_name, metric_col) {
  if (nrow(df) == 0L || !metric_col %in% names(df) || !"rep_id" %in% names(df)) {
    return(tibble())
  }
  
  group_cols <- intersect(metric_id_cols, names(df))
  
  work_df <- df %>%
    mutate(
      metric_value = suppressWarnings(as.numeric(.data[[metric_col]])),
      rep_id = as.integer(rep_id)
    ) %>%
    select(all_of(c("rep_id", group_cols, "metric_value")))
  
  if (all(is.na(work_df$metric_value))) {
    return(tibble())
  }
  
  work_df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      table_name = table_name,
      metric_name = metric_col,
      n_replicates_nonmissing = n_distinct(rep_id[!is.na(metric_value)]),
      mean_value = mean(metric_value, na.rm = TRUE),
      sd_value = stats::sd(metric_value, na.rm = TRUE),
      min_value = suppressWarnings(min(metric_value, na.rm = TRUE)),
      q25_value = safe_quantile(metric_value, 0.25),
      median_value = safe_quantile(metric_value, 0.50),
      q75_value = safe_quantile(metric_value, 0.75),
      max_value = suppressWarnings(max(metric_value, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      min_value = ifelse(is.infinite(min_value), NA_real_, min_value),
      max_value = ifelse(is.infinite(max_value), NA_real_, max_value)
    )
}

summarize_metric_table <- function(df, table_name, metric_cols) {
  if (nrow(df) == 0L || length(metric_cols) == 0L) {
    return(tibble())
  }
  
  bind_rows(lapply(metric_cols, function(metric_col) {
    summarize_one_metric_column(df, table_name, metric_col)
  }))
}

extract_top_metric_values <- function(top_cells, bootstrap_results) {
  if (nrow(top_cells) == 0L) {
    return(
      tibble(
        metric_cell_id = integer(),
        table_name = character(),
        metric_name = character(),
        rep_id = integer(),
        dataset = character(),
        formula_variant = character(),
        site_branch = character(),
        interaction_branch = character(),
        model_id = character(),
        model_class = character(),
        model_block = character(),
        model_family = character(),
        latency_type = character(),
        matched_nocure_model_id = character(),
        component = character(),
        term = character(),
        horizon_year = integer(),
        threshold = double(),
        metric_value = double(),
        metric_key = character()
      )
    )
  }
  
  out <- list()
  
  for (table_name in unique(top_cells$table_name)) {
    df <- bootstrap_results[[table_name]]
    if (is.null(df) || nrow(df) == 0L) next
    
    table_cells <- top_cells %>% filter(table_name == !!table_name)
    
    for (metric_name in unique(table_cells$metric_name)) {
      if (!metric_name %in% names(df)) next
      
      cells_sub <- table_cells %>% filter(metric_name == !!metric_name)
      join_cols <- intersect(metric_id_cols, names(df))
      
      src <- df %>%
        mutate(
          rep_id = as.integer(rep_id),
          metric_value = suppressWarnings(as.numeric(.data[[metric_name]]))
        ) %>%
        select(all_of(c("rep_id", join_cols, "metric_value")))
      
      cell_map <- cells_sub %>%
        select(metric_cell_id, metric_key, all_of(join_cols)) %>%
        distinct()
      
      joined <- dplyr::inner_join(
        src,
        cell_map,
        by = join_cols,
        na_matches = "na"
      ) %>%
        mutate(
          table_name = table_name,
          metric_name = metric_name
        ) %>%
        complete_metric_id_columns() %>%
        select(
          metric_cell_id, table_name, metric_name, rep_id,
          dataset, formula_variant, site_branch, interaction_branch,
          model_id, model_class, model_block, model_family, latency_type,
          matched_nocure_model_id, component, term, horizon_year, threshold,
          metric_value, metric_key
        )
      
      out[[paste(table_name, metric_name, sep = "__")]] <- joined
    }
  }
  
  bind_rows(out) %>% arrange(metric_cell_id, rep_id)
}

write_metric_values_long_streaming <- function(bootstrap_results, out_file) {
  if (file.exists(out_file)) file.remove(out_file)
  metric_specs <- get_metric_specs(bootstrap_results)
  
  for (table_name in names(metric_specs)) {
    spec <- metric_specs[[table_name]]
    df <- spec$df
    metric_cols <- spec$metric_cols
    
    if (nrow(df) == 0L || length(metric_cols) == 0L || !"rep_id" %in% names(df)) next
    
    for (metric_name in metric_cols) {
      if (!metric_name %in% names(df)) next
      
      chunk_df <- df %>%
        transmute(
          table_name = table_name,
          metric_name = metric_name,
          rep_id = as.integer(rep_id),
          dataset = if ("dataset" %in% names(df)) as.character(dataset) else NA_character_,
          formula_variant = if ("formula_variant" %in% names(df)) as.character(formula_variant) else NA_character_,
          site_branch = if ("site_branch" %in% names(df)) as.character(site_branch) else NA_character_,
          interaction_branch = if ("interaction_branch" %in% names(df)) as.character(interaction_branch) else NA_character_,
          model_id = if ("model_id" %in% names(df)) as.character(model_id) else NA_character_,
          model_class = if ("model_class" %in% names(df)) as.character(model_class) else NA_character_,
          model_block = if ("model_block" %in% names(df)) as.character(model_block) else NA_character_,
          model_family = if ("model_family" %in% names(df)) as.character(model_family) else NA_character_,
          latency_type = if ("latency_type" %in% names(df)) as.character(latency_type) else NA_character_,
          matched_nocure_model_id = if ("matched_nocure_model_id" %in% names(df)) as.character(matched_nocure_model_id) else NA_character_,
          component = if ("component" %in% names(df)) as.character(component) else NA_character_,
          term = if ("term" %in% names(df)) as.character(term) else NA_character_,
          horizon_year = if ("horizon_year" %in% names(df)) as.integer(horizon_year) else NA_integer_,
          threshold = if ("threshold" %in% names(df)) as.numeric(threshold) else NA_real_,
          metric_value = suppressWarnings(as.numeric(.data[[metric_name]]))
        )
      
      readr::write_csv(chunk_df, out_file, append = file.exists(out_file))
    }
  }
  
  normalize_existing_path(out_file)
}

# 🔴 Load: current artifacts and supporting inputs ===============================
## 🟠 Load: Stage 7 and Stage 6 objects ===============================
stage7_objects <- load_stage7_outputs(run_root_dir)
stage6_carry_forward_tbl <- load_stage6_carry_forward_optional(stage6_screening_file)
stage6_registry <- build_stage6_registry(stage6_carry_forward_tbl)

# 🔴 Transform: apply Stage 6 carry-forward flags ===============================
## 🟠 Transform: join interpretation flags into Stage 7 tables ===============================
fit_registry_tbl <- join_stage6_flags(stage7_objects$fit_registry, stage6_registry)
coefficients_tbl <- join_stage6_flags(stage7_objects$coefficients, stage6_registry)
subject_predictions_tbl <- join_stage6_flags(stage7_objects$subject_predictions, stage6_registry)
risk_summary_tbl <- join_stage6_flags(stage7_objects$risk_summary, stage6_registry)
delta_risk_tbl <- join_stage6_flags(stage7_objects$delta_risk, stage6_registry)
threshold_metrics_tbl <- join_stage6_flags(stage7_objects$threshold_metrics, stage6_registry)
bootstrap_model_qc_tbl <- join_stage6_flags(stage7_objects$bootstrap_model_qc, stage6_registry)

# 🔴 Rebuild: IPCW registry with dataset key ===============================
## 🟠 Rebuild: dataset-aware IPCW output ===============================
ipcw_registry_tbl <- rebuild_ipcw_registry(stage1_analysis_datasets_file, horizons_year = common_horizons_year)

# 🔴 Transform: threshold reporting suppression ===============================
## 🟠 Transform: blank unsupported PNU threshold horizons ===============================
threshold_metrics_tbl <- apply_threshold_reporting_rules(threshold_metrics_tbl)

# 🔴 Regenerate: visual summary PDF from patched tables ===============================
## 🟠 Regenerate: plot data and PDF export ===============================
plot_risk_tbl <- risk_summary_tbl %>%
  mutate(plot_group = paste(dataset, formula_variant, sep = " | "))

plot_delta_tbl <- delta_risk_tbl %>%
  mutate(plot_group = paste(dataset, formula_variant, sep = " | "))

plot_threshold_tbl <- threshold_metrics_tbl %>%
  filter(horizon_year %in% threshold_plot_horizons_year) %>%
  filter(!threshold_projection_suppressed_flag) %>%
  mutate(plot_group = paste(dataset, formula_variant, sep = " | "))

visual_summary_pdf <- make_output_path("stage7_visual_summary.pdf")
grDevices::pdf(visual_summary_pdf, width = 12, height = 8)

if (nrow(plot_risk_tbl) > 0L) {
  print(
    ggplot(plot_risk_tbl, aes(x = horizon_year, y = mean_risk_overall, color = model_id, group = model_id)) +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      geom_point(size = 1.2, na.rm = TRUE) +
      geom_ribbon(
        aes(ymin = mean_risk_ci_lower, ymax = mean_risk_ci_upper, fill = model_id),
        alpha = 0.12,
        linewidth = 0,
        inherit.aes = TRUE,
        na.rm = TRUE,
        show.legend = FALSE
      ) +
      facet_wrap(~ plot_group, scales = "free_y") +
      scale_x_continuous(breaks = common_horizons_year) +
      labs(
        title = "Stage 7: overall predicted transition risk by horizon",
        x = "Years after cohort entry",
        y = "Mean predicted overall risk",
        color = "Model"
      ) +
      theme_bw()
  )
}

if (nrow(plot_delta_tbl) > 0L) {
  print(
    ggplot(plot_delta_tbl, aes(x = horizon_year, y = mean_delta_risk_nc_minus_cure, color = model_id, group = model_id)) +
      geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      geom_point(size = 1.2, na.rm = TRUE) +
      geom_ribbon(
        aes(ymin = mean_delta_risk_ci_lower, ymax = mean_delta_risk_ci_upper, fill = model_id),
        alpha = 0.12,
        linewidth = 0,
        inherit.aes = TRUE,
        na.rm = TRUE,
        show.legend = FALSE
      ) +
      facet_wrap(~ plot_group, scales = "free_y") +
      scale_x_continuous(breaks = common_horizons_year) +
      labs(
        title = "Stage 7: DeltaRisk_NC_C(t) = risk(no-cure) - risk(cure)",
        x = "Years after cohort entry",
        y = "Mean delta risk (no-cure minus cure)",
        color = "Cure model"
      ) +
      theme_bw()
  )
}

if (nrow(plot_threshold_tbl) > 0L) {
  print(
    ggplot(plot_threshold_tbl, aes(x = threshold, y = false_positive_burden_primary, color = model_id, group = model_id)) +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      geom_point(size = 1.2, na.rm = TRUE) +
      geom_ribbon(
        aes(ymin = false_positive_burden_primary_ci_lower, ymax = false_positive_burden_primary_ci_upper, fill = model_id),
        alpha = 0.12,
        linewidth = 0,
        inherit.aes = TRUE,
        na.rm = TRUE,
        show.legend = FALSE
      ) +
      facet_grid(horizon_year ~ plot_group, scales = "free_y") +
      labs(
        title = "Stage 7: false-positive burden among non-events",
        x = "Risk threshold",
        y = "False-positive burden among non-events (IPCW-weighted)",
        color = "Model"
      ) +
      theme_bw()
  )
  
  print(
    ggplot(plot_threshold_tbl, aes(x = threshold, y = net_benefit, color = model_id, group = model_id)) +
      geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      geom_point(size = 1.2, na.rm = TRUE) +
      geom_ribbon(
        aes(ymin = net_benefit_ci_lower, ymax = net_benefit_ci_upper, fill = model_id),
        alpha = 0.12,
        linewidth = 0,
        inherit.aes = TRUE,
        na.rm = TRUE,
        show.legend = FALSE
      ) +
      facet_grid(horizon_year ~ plot_group, scales = "free_y") +
      labs(
        title = "Stage 7: net benefit across thresholds",
        x = "Risk threshold",
        y = "Net benefit (IPCW-weighted)",
        color = "Model"
      ) +
      theme_bw()
  )
}

grDevices::dev.off()
visual_summary_pdf <- normalize_existing_path(visual_summary_pdf)

# 🔴 Export: patched Stage 7 tables ===============================
## 🟠 Export: refreshed CSV outputs without refitting ===============================
fit_registry_file <- safe_write_csv(fit_registry_tbl, make_output_path("stage7_fit_registry.csv"))
coefficients_file <- safe_write_csv(coefficients_tbl, make_output_path("stage7_coefficients.csv"))
subject_predictions_file <- safe_write_csv(subject_predictions_tbl, make_output_path("stage7_subject_predictions.csv"))
risk_summary_file <- safe_write_csv(risk_summary_tbl, make_output_path("stage7_risk_summary.csv"))
delta_risk_file <- safe_write_csv(delta_risk_tbl, make_output_path("stage7_delta_risk.csv"))
threshold_metrics_file <- safe_write_csv(threshold_metrics_tbl, make_output_path("stage7_threshold_metrics.csv"))
ipcw_registry_file <- safe_write_csv(ipcw_registry_tbl, make_output_path("stage7_ipcw_registry.csv"))
bootstrap_model_qc_file <- safe_write_csv(bootstrap_model_qc_tbl, make_output_path("stage7_bootstrap_model_qc.csv"))
stage6_registry_file <- safe_write_csv(stage6_carry_forward_tbl, make_output_path("stage7_stage6_carry_forward_registry.csv"))

# 🔴 Extract: memory-safe bootstrap summaries ===============================
## 🟠 Extract: replicate-level outputs from saved bootstrap RDS ===============================
replicate_overview_file <- NA_character_
failure_events_file <- NA_character_
failure_pattern_file <- NA_character_
failure_by_model_file <- NA_character_
metric_distribution_file <- NA_character_
top_metric_cells_file <- NA_character_
top_metric_values_file <- NA_character_
diagnostic_pdf_file <- NA_character_

if (isTRUE(run_bootstrap_extraction) && file.exists(stage7_objects$merged_bootstrap_rds_file)) {
  bootstrap_results <- standardize_bootstrap_results(readRDS(stage7_objects$merged_bootstrap_rds_file))
  
  fit_registry_norm <- normalize_fit_registry(bootstrap_results$fit_registry)
  failure_events <- make_failure_events(fit_registry_norm)
  replicate_overview <- make_replicate_overview(fit_registry_norm)
  failure_pattern_summary <- make_failure_pattern_summary(replicate_overview)
  failure_by_model <- make_failure_by_model(failure_events)
  
  metric_specs <- get_metric_specs(bootstrap_results)
  metric_distribution_summary <- bind_rows(
    lapply(names(metric_specs), function(table_name) {
      spec <- metric_specs[[table_name]]
      summarize_metric_table(spec$df, table_name, spec$metric_cols)
    })
  ) %>%
    mutate(metric_key = make_metric_key(.)) %>%
    arrange(table_name, dataset, formula_variant, model_id, component, term, horizon_year, threshold, metric_name)
  
  top_variable_metric_cells <- metric_distribution_summary %>%
    filter(n_replicates_nonmissing >= 5, is.finite(sd_value), !is.na(sd_value)) %>%
    arrange(desc(sd_value), desc(n_replicates_nonmissing), metric_key) %>%
    slice_head(n = top_n_variable_metric_cells) %>%
    mutate(metric_cell_id = row_number())
  
  top_metric_values_long <- extract_top_metric_values(top_variable_metric_cells, bootstrap_results)
  
  if (isTRUE(write_full_metric_values_long_csv)) {
    write_metric_values_long_streaming(
      bootstrap_results = bootstrap_results,
      out_file = make_output_path("stage7_metric_values_long_full.csv")
    )
  }
  
  replicate_overview_file <- safe_write_csv(replicate_overview, make_output_path("stage7_replicate_overview.csv"))
  failure_events_file <- safe_write_csv(failure_events, make_output_path("stage7_replicate_failure_events.csv"))
  failure_pattern_file <- safe_write_csv(failure_pattern_summary, make_output_path("stage7_replicate_failure_pattern_summary.csv"))
  failure_by_model_file <- safe_write_csv(failure_by_model, make_output_path("stage7_replicate_failure_by_model.csv"))
  metric_distribution_file <- safe_write_csv(metric_distribution_summary, make_output_path("stage7_metric_distribution_summary.csv"))
  top_metric_cells_file <- safe_write_csv(top_variable_metric_cells, make_output_path("stage7_top_variable_metric_cells.csv"))
  top_metric_values_file <- safe_write_csv(top_metric_values_long, make_output_path("stage7_top_variable_metric_values_long.csv"))
  
  diagnostic_pdf_file <- make_output_path("stage7_bootstrap_extraction_plots.pdf")
  grDevices::pdf(diagnostic_pdf_file, width = 12, height = 8)
  
  if (nrow(replicate_overview) > 0L) {
    print(
      ggplot(replicate_overview, aes(x = n_fail_total)) +
        geom_histogram(binwidth = 1, boundary = -0.5, closed = "right") +
        labs(
          title = "Bootstrap replicate failure count distribution",
          x = "Failed model rows per replicate",
          y = "Replicate count"
        ) +
        theme_minimal()
    )
    
    print(
      ggplot(replicate_overview, aes(x = rep_id, y = n_fail_total)) +
        geom_point(alpha = 0.7) +
        geom_line(alpha = 0.35) +
        labs(
          title = "Failure count by bootstrap replicate",
          x = "Replicate ID",
          y = "Failed model rows"
        ) +
        theme_minimal()
    )
  }
  
  if (nrow(failure_pattern_summary) > 0L) {
    plot_df <- failure_pattern_summary %>%
      slice_head(n = top_n_failure_patterns) %>%
      mutate(pattern_label = stringr::str_wrap(failure_signature, width = 80))
    
    print(
      ggplot(plot_df, aes(x = reorder(pattern_label, n_replicates), y = n_replicates)) +
        geom_col() +
        coord_flip() +
        labs(
          title = sprintf("Top %s failure signatures", top_n_failure_patterns),
          x = "Failure signature",
          y = "Replicate count"
        ) +
        theme_minimal()
    )
  }
  
  if (nrow(failure_by_model) > 0L) {
    plot_df <- failure_by_model %>%
      slice_head(n = top_n_failure_models) %>%
      mutate(model_label = stringr::str_wrap(paste(dataset, formula_variant, model_id, failure_type, sep = " | "), width = 70))
    
    print(
      ggplot(plot_df, aes(x = reorder(model_label, n_replicates), y = n_replicates)) +
        geom_col() +
        coord_flip() +
        labs(
          title = sprintf("Top %s failing model patterns", top_n_failure_models),
          x = "Dataset | formula | model | failure type",
          y = "Replicate count"
        ) +
        theme_minimal()
    )
  }
  
  if (nrow(top_metric_values_long) > 0L) {
    plot_df <- top_metric_values_long %>%
      mutate(metric_key_wrapped = stringr::str_wrap(metric_key, width = 70))
    
    print(
      ggplot(plot_df, aes(x = metric_key_wrapped, y = metric_value)) +
        geom_boxplot(outlier.alpha = 0.25) +
        coord_flip() +
        labs(
          title = sprintf("Top %s most variable metric cells", top_n_variable_metric_cells),
          x = "Metric cell",
          y = "Metric value across replicates"
        ) +
        theme_minimal()
    )
  }
  
  grDevices::dev.off()
  diagnostic_pdf_file <- normalize_existing_path(diagnostic_pdf_file)
}

# 🔴 Assemble: refreshed Stage 7 output manifest ===============================
## 🟠 Assemble: one manifest covering refreshed and extracted outputs ===============================
output_manifest <- tibble(
  file_name = c(
    basename(fit_registry_file),
    basename(coefficients_file),
    basename(subject_predictions_file),
    basename(risk_summary_file),
    basename(delta_risk_file),
    basename(threshold_metrics_file),
    basename(ipcw_registry_file),
    basename(bootstrap_model_qc_file),
    basename(stage6_registry_file),
    basename(visual_summary_pdf),
    if (!is.na(replicate_overview_file)) basename(replicate_overview_file) else NA_character_,
    if (!is.na(failure_events_file)) basename(failure_events_file) else NA_character_,
    if (!is.na(failure_pattern_file)) basename(failure_pattern_file) else NA_character_,
    if (!is.na(failure_by_model_file)) basename(failure_by_model_file) else NA_character_,
    if (!is.na(metric_distribution_file)) basename(metric_distribution_file) else NA_character_,
    if (!is.na(top_metric_cells_file)) basename(top_metric_cells_file) else NA_character_,
    if (!is.na(top_metric_values_file)) basename(top_metric_values_file) else NA_character_,
    if (!is.na(diagnostic_pdf_file)) basename(diagnostic_pdf_file) else NA_character_,
    basename(stage7_objects$merged_bootstrap_rds_file)
  ),
  file_path = c(
    fit_registry_file,
    coefficients_file,
    subject_predictions_file,
    risk_summary_file,
    delta_risk_file,
    threshold_metrics_file,
    ipcw_registry_file,
    bootstrap_model_qc_file,
    stage6_registry_file,
    visual_summary_pdf,
    replicate_overview_file,
    failure_events_file,
    failure_pattern_file,
    failure_by_model_file,
    metric_distribution_file,
    top_metric_cells_file,
    top_metric_values_file,
    diagnostic_pdf_file,
    normalize_existing_path(stage7_objects$merged_bootstrap_rds_file)
  ),
  description = c(
    "Patched fit registry with Stage 6 carry-forward flags.",
    "Patched coefficient table with Stage 6 carry-forward flags.",
    "Patched subject-by-horizon predictions with Stage 6 carry-forward flags.",
    "Patched risk summary with Stage 6 carry-forward flags.",
    "Patched delta-risk table with Stage 6 carry-forward flags.",
    "Patched threshold metrics with dataset-aware suppression for PNU years 3-10.",
    "Rebuilt IPCW registry including the dataset key.",
    "Patched bootstrap model QC with Stage 6 carry-forward flags.",
    "Source-of-truth Stage 6 carry-forward table used for the joins.",
    "Regenerated Stage 7 visual summary based on patched tables.",
    "Bootstrap replicate overview.",
    "Bootstrap failure events by replicate and model row.",
    "Bootstrap failure-pattern frequency table.",
    "Bootstrap failure counts aggregated by model.",
    "Across-replicate summary statistics for each metric cell.",
    "Top variable metric cells selected from the summary table.",
    "Replicate-level raw values for the top variable metric cells.",
    "Bootstrap extraction diagnostic PDF.",
    "Saved merged bootstrap results RDS reused without rerunning bootstrap."
  )
) %>%
  filter(!is.na(file_name)) %>%
  mutate(
    exists = file.exists(file_path),
    size_bytes = ifelse(exists, file.info(file_path)$size, NA_real_)
  )

output_manifest_file <- safe_write_csv(output_manifest, make_output_path("stage7_output_manifest.csv"))

message("✅ Stage 7 refresh completed without refitting the main models or rerunning bootstrap.")
print(output_manifest)
