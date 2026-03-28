# 🔴 Configure: old/new Stage-1 paths ===============================
old_stage1_dir <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/stage1_Backbone lock_old'
new_stage1_dir <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/stage1_Backbone lock'
export_path <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/stage1_check output'

numeric_tolerance <- 1e-10
max_example_ids_per_field <- 5L
required_dataset_names <- c("PNU", "SNU", "merged")

core_backbone_compare_cols <- c(
  "unique_person_id",
  "site",
  "id",
  "sex_num",
  "age_exact_entry",
  "age_s",
  "days_followup",
  "status_num",
  "time_year",
  "event_main",
  "censor_main"
)

metadata_fit_input_names <- c(
  "unique_person_id_rule",
  "analysis_time_variable",
  "reporting_time_variable",
  "event_definition",
  "main_censoring_definition",
  "age_variable",
  "age_scaling_rule"
)

# 🔴 Initialize: packages and runtime options ===============================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# 🔴 Define: file and scalar helpers ===============================
## 🟠 Define: path guards and serializers ===============================
assert_dir_exists <- function(path, label) {
  if (!dir.exists(path)) {
    stop(sprintf("%s does not exist: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

assert_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s does not exist: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

normalize_existing_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

make_output_path <- function(file_name) {
  file.path(export_path, file_name)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    return(tibble())
  }
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

safe_write_csv <- function(df, path) {
  readr::write_csv(df, path)
  normalize_existing_path(path)
}

scalar_to_compare_string <- function(x) {
  if (length(x) == 0L) {
    return(NA_character_)
  }
  if (inherits(x, "POSIXt")) {
    return(format(x, tz = "UTC", usetz = TRUE))
  }
  if (inherits(x, "Date")) {
    return(as.character(x))
  }
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (is.logical(x)) {
    return(ifelse(is.na(x), NA_character_, ifelse(x, "TRUE", "FALSE")))
  }
  if (is.numeric(x)) {
    return(ifelse(is.na(x), NA_character_, format(signif(x, 15), scientific = FALSE, trim = TRUE)))
  }
  if (is.integer(x)) {
    return(ifelse(is.na(x), NA_character_, as.character(x)))
  }
  if (is.list(x)) {
    return(ifelse(lengths(x) == 0L, NA_character_, vapply(x, function(xx) paste(as.character(xx), collapse = "|"), character(1))))
  }
  out <- trimws(as.character(x))
  out[is.na(x)] <- NA_character_
  out
}

vector_to_compare_string <- function(x) {
  out <- scalar_to_compare_string(x)
  if (length(out) == 1L && length(x) > 1L) {
    out <- rep(out, length(x))
  }
  out
}

first_nonempty <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) return(NA_character_)
  x[[1]]
}

# 🔴 Define: Stage-1 artifact readers ===============================
## 🟠 Define: optional bundle and registry loaders ===============================
read_stage1_analysis_datasets <- function(stage1_dir) {
  direct_file <- file.path(stage1_dir, "stage1_analysis_datasets.rds")
  bundle_file <- file.path(stage1_dir, "stage1_backbone_bundle.rds")
  
  if (file.exists(direct_file)) {
    return(readRDS(direct_file))
  }
  
  if (file.exists(bundle_file)) {
    bundle <- readRDS(bundle_file)
    if (is.list(bundle) && !is.null(bundle$datasets)) {
      return(bundle$datasets)
    }
  }
  
  stop(
    sprintf(
      "Could not find `stage1_analysis_datasets.rds` or dataset contents inside `stage1_backbone_bundle.rds` under: %s",
      stage1_dir
    ),
    call. = FALSE
  )
}

read_stage1_bundle_optional <- function(stage1_dir) {
  bundle_file <- file.path(stage1_dir, "stage1_backbone_bundle.rds")
  if (!file.exists(bundle_file)) {
    return(NULL)
  }
  readRDS(bundle_file)
}

read_stage1_artifacts <- function(stage1_dir, label) {
  assert_dir_exists(stage1_dir, label)
  
  list(
    label = label,
    stage1_dir = normalize_existing_path(stage1_dir),
    analysis_datasets = read_stage1_analysis_datasets(stage1_dir),
    backbone_bundle = read_stage1_bundle_optional(stage1_dir),
    dataset_registry = safe_read_csv(file.path(stage1_dir, "stage1_dataset_registry.csv")),
    formula_registry = safe_read_csv(file.path(stage1_dir, "stage1_formula_registry.csv")),
    horizon_registry = safe_read_csv(file.path(stage1_dir, "stage1_horizon_registry.csv")),
    threshold_registry = safe_read_csv(file.path(stage1_dir, "stage1_threshold_registry.csv")),
    metadata_registry = safe_read_csv(file.path(stage1_dir, "stage1_metadata_registry.csv")),
    modeling_registry = safe_read_csv(file.path(stage1_dir, "stage1_modeling_registry.csv"))
  )
}

# 🔴 Define: dataset standardization helpers ===============================
## 🟠 Define: backbone derivations and coercions ===============================
coerce_numeric_text <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

coerce_integer_text <- function(x) {
  suppressWarnings(as.integer(as.character(x)))
}

standardize_stage1_dataset <- function(df, dataset_name) {
  if (!is.data.frame(df)) {
    stop(sprintf("[%s] Stage-1 dataset object must be a data frame.", dataset_name), call. = FALSE)
  }
  
  df <- tibble::as_tibble(df)
  required_minimal_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_minimal_cols <- setdiff(required_minimal_cols, names(df))
  
  if (length(missing_minimal_cols) > 0L) {
    stop(
      sprintf("[%s] Missing required minimal columns: %s", dataset_name, paste(missing_minimal_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  
  df <- df %>%
    mutate(
      id = trimws(as.character(id)),
      site = trimws(as.character(site)),
      sex_num = coerce_integer_text(sex_num),
      age_exact_entry = coerce_numeric_text(age_exact_entry),
      days_followup = coerce_numeric_text(days_followup),
      status_num = coerce_integer_text(status_num)
    )
  
  if (!"unique_person_id" %in% names(df)) {
    df <- df %>% mutate(unique_person_id = paste(site, id, sep = "_"))
  } else {
    df <- df %>% mutate(unique_person_id = trimws(as.character(unique_person_id)))
  }
  
  if (!"time_year" %in% names(df)) {
    df <- df %>% mutate(time_year = days_followup / 365.25)
  } else {
    df <- df %>% mutate(time_year = coerce_numeric_text(time_year))
  }
  
  if (!"event_main" %in% names(df)) {
    df <- df %>% mutate(event_main = as.integer(status_num == 1L))
  } else {
    df <- df %>% mutate(event_main = coerce_integer_text(event_main))
  }
  
  if (!"censor_main" %in% names(df)) {
    df <- df %>% mutate(censor_main = as.integer(status_num %in% c(0L, 2L)))
  } else {
    df <- df %>% mutate(censor_main = coerce_integer_text(censor_main))
  }
  
  if (!"age_s" %in% names(df)) {
    age_sd <- stats::sd(df$age_exact_entry, na.rm = TRUE)
    if (!is.finite(age_sd) || age_sd <= 0) {
      stop(sprintf("[%s] Could not derive `age_s` because `age_exact_entry` has non-positive SD.", dataset_name), call. = FALSE)
    }
    df <- df %>% mutate(age_s = (age_exact_entry - mean(age_exact_entry, na.rm = TRUE)) / (2 * age_sd))
  } else {
    df <- df %>% mutate(age_s = coerce_numeric_text(age_s))
  }
  
  required_compare_cols <- core_backbone_compare_cols
  missing_compare_cols <- setdiff(required_compare_cols, names(df))
  if (length(missing_compare_cols) > 0L) {
    stop(
      sprintf("[%s] Missing compare-ready columns after standardization: %s", dataset_name, paste(missing_compare_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  
  if (nrow(df) == 0L) {
    stop(sprintf("[%s] Dataset has zero rows.", dataset_name), call. = FALSE)
  }
  
  if (anyNA(df[, required_compare_cols])) {
    stop(sprintf("[%s] Missing values detected in compare-ready columns.", dataset_name), call. = FALSE)
  }
  
  if (any(duplicated(df$unique_person_id))) {
    stop(sprintf("[%s] `unique_person_id` must be unique for equivalence checking.", dataset_name), call. = FALSE)
  }
  
  df %>%
    select(all_of(required_compare_cols)) %>%
    arrange(unique_person_id)
}

# 🔴 Define: dataset comparison engine ===============================
## 🟠 Define: row-wise and column-wise comparators ===============================
compare_dataset_pair <- function(old_df, new_df, dataset_name, numeric_tol, max_examples) {
  numeric_fields <- c("sex_num", "age_exact_entry", "age_s", "days_followup", "status_num", "time_year", "event_main", "censor_main")
  key_col <- "unique_person_id"
  compare_cols <- setdiff(core_backbone_compare_cols, key_col)
  
  old_ids <- old_df %>% select(unique_person_id)
  new_ids <- new_df %>% select(unique_person_id)
  
  key_only_old <- anti_join(old_ids, new_ids, by = key_col) %>%
    mutate(dataset = dataset_name, presence_side = "old_only")
  key_only_new <- anti_join(new_ids, old_ids, by = key_col) %>%
    mutate(dataset = dataset_name, presence_side = "new_only")
  
  common_old <- semi_join(old_df, new_ids, by = key_col)
  common_new <- semi_join(new_df, old_ids, by = key_col)
  joined <- inner_join(common_old, common_new, by = key_col, suffix = c(".old", ".new"))
  
  col_diff_list <- vector("list", length(compare_cols))
  row_diff_list <- vector("list", length(compare_cols))
  
  for (ii in seq_along(compare_cols)) {
    col_name <- compare_cols[[ii]]
    old_vec <- joined[[paste0(col_name, ".old")]]
    new_vec <- joined[[paste0(col_name, ".new")]]
    
    na_mismatch <- xor(is.na(old_vec), is.na(new_vec))
    
    if (col_name %in% numeric_fields) {
      abs_diff <- ifelse(is.na(old_vec) | is.na(new_vec), NA_real_, abs(as.numeric(old_vec) - as.numeric(new_vec)))
      value_mismatch <- !is.na(abs_diff) & abs_diff > numeric_tol
      max_abs_diff <- suppressWarnings(max(abs_diff, na.rm = TRUE))
      if (!is.finite(max_abs_diff)) max_abs_diff <- NA_real_
    } else {
      abs_diff <- rep(NA_real_, length(old_vec))
      value_mismatch <- !is.na(old_vec) & !is.na(new_vec) & (as.character(old_vec) != as.character(new_vec))
      max_abs_diff <- NA_real_
    }
    
    mismatch_flag <- na_mismatch | value_mismatch
    mismatch_ids <- joined[[key_col]][mismatch_flag]
    example_ids <- paste(head(mismatch_ids, max_examples), collapse = "|")
    if (!nzchar(example_ids)) example_ids <- NA_character_
    
    col_diff_list[[ii]] <- tibble(
      dataset = dataset_name,
      variable = col_name,
      n_compared = nrow(joined),
      n_diff = sum(mismatch_flag, na.rm = TRUE),
      max_abs_diff = max_abs_diff,
      example_unique_person_ids = example_ids
    )
    
    if (sum(mismatch_flag, na.rm = TRUE) > 0L) {
      row_diff_list[[ii]] <- tibble(
        dataset = dataset_name,
        unique_person_id = joined[[key_col]][mismatch_flag],
        variable = col_name,
        old_value = vector_to_compare_string(old_vec[mismatch_flag]),
        new_value = vector_to_compare_string(new_vec[mismatch_flag]),
        abs_diff = abs_diff[mismatch_flag]
      )
    } else {
      row_diff_list[[ii]] <- tibble(
        dataset = character(),
        unique_person_id = character(),
        variable = character(),
        old_value = character(),
        new_value = character(),
        abs_diff = double()
      )
    }
  }
  
  column_differences <- bind_rows(col_diff_list)
  row_differences <- bind_rows(row_diff_list)
  total_value_diff_cells <- sum(column_differences$n_diff, na.rm = TRUE)
  
  summary_tbl <- tibble(
    dataset = dataset_name,
    old_n_rows = nrow(old_df),
    new_n_rows = nrow(new_df),
    old_n_unique_person_id = dplyr::n_distinct(old_df$unique_person_id),
    new_n_unique_person_id = dplyr::n_distinct(new_df$unique_person_id),
    only_old_key_n = nrow(key_only_old),
    only_new_key_n = nrow(key_only_new),
    value_difference_cell_n = total_value_diff_cells,
    exact_key_match_flag = nrow(key_only_old) == 0L && nrow(key_only_new) == 0L,
    exact_value_match_flag = total_value_diff_cells == 0L,
    backbone_equivalent_flag = nrow(key_only_old) == 0L && nrow(key_only_new) == 0L && total_value_diff_cells == 0L
  )
  
  key_differences <- bind_rows(key_only_old, key_only_new) %>%
    select(dataset, presence_side, unique_person_id)
  
  list(
    summary = summary_tbl,
    column_differences = column_differences,
    row_differences = row_differences,
    key_differences = key_differences
  )
}

# 🔴 Define: registry comparison engine ===============================
## 🟠 Define: key-based registry comparators ===============================
compare_table_on_keys <- function(old_tbl, new_tbl, table_name, key_cols, compare_cols) {
  old_tbl <- tibble::as_tibble(old_tbl)
  new_tbl <- tibble::as_tibble(new_tbl)
  
  missing_key_old <- setdiff(key_cols, names(old_tbl))
  missing_key_new <- setdiff(key_cols, names(new_tbl))
  missing_compare_old <- setdiff(compare_cols, names(old_tbl))
  missing_compare_new <- setdiff(compare_cols, names(new_tbl))
  
  summary_tbl <- tibble(
    table_name = table_name,
    old_n_rows = nrow(old_tbl),
    new_n_rows = nrow(new_tbl),
    missing_key_old = if (length(missing_key_old) == 0L) NA_character_ else paste(missing_key_old, collapse = "|"),
    missing_key_new = if (length(missing_key_new) == 0L) NA_character_ else paste(missing_key_new, collapse = "|"),
    missing_compare_old = if (length(missing_compare_old) == 0L) NA_character_ else paste(missing_compare_old, collapse = "|"),
    missing_compare_new = if (length(missing_compare_new) == 0L) NA_character_ else paste(missing_compare_new, collapse = "|"),
    duplicate_key_old_flag = NA,
    duplicate_key_new_flag = NA,
    only_old_key_n = NA_integer_,
    only_new_key_n = NA_integer_,
    value_difference_cell_n = NA_integer_,
    equivalent_flag = FALSE
  )
  
  empty_value_differences <- tibble::as_tibble(setNames(
    replicate(length(c("table_name", key_cols, "variable", "old_value", "new_value")), character(), simplify = FALSE),
    c("table_name", key_cols, "variable", "old_value", "new_value")
  ))
  empty_key_differences <- tibble::as_tibble(setNames(
    replicate(length(c("table_name", "presence_side", key_cols)), character(), simplify = FALSE),
    c("table_name", "presence_side", key_cols)
  ))
  
  if (length(missing_key_old) > 0L || length(missing_key_new) > 0L || length(missing_compare_old) > 0L || length(missing_compare_new) > 0L) {
    return(list(
      summary = summary_tbl,
      value_differences = empty_value_differences,
      key_differences = empty_key_differences
    ))
  }
  
  old_cmp <- old_tbl %>%
    select(all_of(c(key_cols, compare_cols))) %>%
    mutate(across(all_of(c(key_cols, compare_cols)), vector_to_compare_string))
  
  new_cmp <- new_tbl %>%
    select(all_of(c(key_cols, compare_cols))) %>%
    mutate(across(all_of(c(key_cols, compare_cols)), vector_to_compare_string))
  
  duplicate_key_old_flag <- any(duplicated(old_cmp[key_cols]))
  duplicate_key_new_flag <- any(duplicated(new_cmp[key_cols]))
  
  if (duplicate_key_old_flag || duplicate_key_new_flag) {
    summary_tbl <- summary_tbl %>%
      mutate(
        duplicate_key_old_flag = duplicate_key_old_flag,
        duplicate_key_new_flag = duplicate_key_new_flag,
        only_old_key_n = NA_integer_,
        only_new_key_n = NA_integer_,
        value_difference_cell_n = NA_integer_,
        equivalent_flag = FALSE
      )
    
    return(list(
      summary = summary_tbl,
      value_differences = empty_value_differences,
      key_differences = empty_key_differences
    ))
  }
  
  key_only_old <- anti_join(old_cmp %>% select(all_of(key_cols)), new_cmp %>% select(all_of(key_cols)), by = key_cols) %>%
    mutate(table_name = table_name, presence_side = "old_only")
  key_only_new <- anti_join(new_cmp %>% select(all_of(key_cols)), old_cmp %>% select(all_of(key_cols)), by = key_cols) %>%
    mutate(table_name = table_name, presence_side = "new_only")
  
  joined <- inner_join(old_cmp, new_cmp, by = key_cols, suffix = c(".old", ".new"))
  
  diff_tables <- lapply(compare_cols, function(col_name) {
    old_vec <- joined[[paste0(col_name, ".old")]]
    new_vec <- joined[[paste0(col_name, ".new")]]
    mismatch_flag <- xor(is.na(old_vec), is.na(new_vec)) | (!is.na(old_vec) & !is.na(new_vec) & old_vec != new_vec)
    
    if (!any(mismatch_flag)) {
      return(empty_value_differences[0, , drop = FALSE])
    }
    
    tibble(
      table_name = table_name,
      variable = col_name,
      old_value = old_vec[mismatch_flag],
      new_value = new_vec[mismatch_flag]
    ) %>%
      bind_cols(joined[mismatch_flag, key_cols, drop = FALSE])
  })
  
  value_differences <- bind_rows(diff_tables)
  if (ncol(value_differences) == 0L) value_differences <- empty_value_differences
  
  summary_tbl <- summary_tbl %>%
    mutate(
      duplicate_key_old_flag = duplicate_key_old_flag,
      duplicate_key_new_flag = duplicate_key_new_flag,
      only_old_key_n = nrow(key_only_old),
      only_new_key_n = nrow(key_only_new),
      value_difference_cell_n = nrow(value_differences),
      equivalent_flag = nrow(key_only_old) == 0L && nrow(key_only_new) == 0L && nrow(value_differences) == 0L
    )
  
  key_differences <- bind_rows(key_only_old, key_only_new) %>%
    select(table_name, presence_side, all_of(key_cols))
  if (nrow(key_differences) == 0L && ncol(key_differences) == 0L) key_differences <- empty_key_differences
  
  list(
    summary = summary_tbl,
    value_differences = value_differences,
    key_differences = key_differences
  )
}

# 🔴 Define: decision helpers ===============================
## 🟠 Define: reuse classification and narratives ===============================
get_equivalent_flag <- function(summary_tbl, table_name) {
  row <- summary_tbl %>% filter(.data$table_name == .env$table_name)
  if (nrow(row) != 1L) {
    return(FALSE)
  }
  isTRUE(row$equivalent_flag[[1]])
}

make_decision_text <- function(fit_input_backbone_equivalent, reporting_contract_equivalent) {
  if (fit_input_backbone_equivalent && reporting_contract_equivalent) {
    return("existing_stage7_fit_results_can_be_reused_without_refit")
  }
  if (fit_input_backbone_equivalent && !reporting_contract_equivalent) {
    return("existing_stage7_fit_results_can_be_reused_but_stage7_postfit_exports_should_be_refreshed")
  }
  "stage7_refit_required"
}

make_interpretation_text <- function(fit_input_backbone_equivalent, reporting_contract_equivalent) {
  if (fit_input_backbone_equivalent && reporting_contract_equivalent) {
    return("Subject-level Stage-1 modeling backbone and shared Stage-7 postfit contract match; reuse is safe without refit.")
  }
  if (fit_input_backbone_equivalent && !reporting_contract_equivalent) {
    return("Subject-level Stage-1 modeling backbone matches, so model fitting can be reused; however, horizon/threshold reporting contract changed and postfit exports should be regenerated.")
  }
  "At least one fit-critical Stage-1 element changed (subject-level backbone and/or shared fit-input contract), so Stage 7 should be refit against the new backbone."
}

# 🔴 Load: old and new Stage-1 artifacts ===============================
## 🟠 Read: old and new output folders ===============================
old_stage1 <- read_stage1_artifacts(old_stage1_dir, "Old Stage 1")
new_stage1 <- read_stage1_artifacts(new_stage1_dir, "New Stage 1")

# 🔴 Compare: subject-level Stage-1 analysis datasets ===============================
## 🟠 Evaluate: PNU, SNU, and merged backbone equality ===============================
dataset_compare_results <- lapply(required_dataset_names, function(dataset_name) {
  if (!dataset_name %in% names(old_stage1$analysis_datasets)) {
    stop(sprintf("Old Stage-1 outputs are missing dataset `%s`.", dataset_name), call. = FALSE)
  }
  if (!dataset_name %in% names(new_stage1$analysis_datasets)) {
    stop(sprintf("New Stage-1 outputs are missing dataset `%s`.", dataset_name), call. = FALSE)
  }
  
  old_df <- standardize_stage1_dataset(old_stage1$analysis_datasets[[dataset_name]], dataset_name)
  new_df <- standardize_stage1_dataset(new_stage1$analysis_datasets[[dataset_name]], dataset_name)
  
  compare_dataset_pair(
    old_df = old_df,
    new_df = new_df,
    dataset_name = dataset_name,
    numeric_tol = numeric_tolerance,
    max_examples = max_example_ids_per_field
  )
})

dataset_backbone_summary <- bind_rows(lapply(dataset_compare_results, `[[`, "summary"))
dataset_column_differences <- bind_rows(lapply(dataset_compare_results, `[[`, "column_differences"))
dataset_row_value_differences <- bind_rows(lapply(dataset_compare_results, `[[`, "row_differences"))
dataset_key_presence_differences <- bind_rows(lapply(dataset_compare_results, `[[`, "key_differences"))

# 🔴 Compare: shared fit-input and reporting registries ===============================
## 🟠 Evaluate: metadata, formula, horizon, and threshold contracts ===============================
metadata_core_old <- old_stage1$metadata_registry %>%
  filter(.data$metadata_name %in% metadata_fit_input_names) %>%
  select(metadata_name, metadata_value)

metadata_core_new <- new_stage1$metadata_registry %>%
  filter(.data$metadata_name %in% metadata_fit_input_names) %>%
  select(metadata_name, metadata_value)

metadata_core_cmp <- compare_table_on_keys(
  old_tbl = metadata_core_old,
  new_tbl = metadata_core_new,
  table_name = "metadata_core_fit_input",
  key_cols = c("metadata_name"),
  compare_cols = c("metadata_value")
)

formula_common_cols <- c("formula_rhs", "uses_site", "uses_age_sex_interaction", "site_branch", "interaction_branch")
formula_cmp <- compare_table_on_keys(
  old_tbl = old_stage1$formula_registry,
  new_tbl = new_stage1$formula_registry,
  table_name = "formula_registry_fit_input",
  key_cols = c("dataset", "formula_name"),
  compare_cols = formula_common_cols
)

horizon_common_cols <- intersect(c("interpretation_tier", "primary_supported_flag"), intersect(names(old_stage1$horizon_registry), names(new_stage1$horizon_registry)))
horizon_cmp <- compare_table_on_keys(
  old_tbl = old_stage1$horizon_registry,
  new_tbl = new_stage1$horizon_registry,
  table_name = "horizon_registry_reporting_contract",
  key_cols = c("dataset", "horizon_year"),
  compare_cols = horizon_common_cols
)

threshold_compare_cols <- intersect(c("threshold_label", "positive_rule", "editable_in_script_config"), intersect(names(old_stage1$threshold_registry), names(new_stage1$threshold_registry)))
threshold_cmp <- compare_table_on_keys(
  old_tbl = old_stage1$threshold_registry,
  new_tbl = new_stage1$threshold_registry,
  table_name = "threshold_registry_reporting_contract",
  key_cols = c("threshold"),
  compare_cols = threshold_compare_cols
)

registry_comparison_summary <- bind_rows(
  metadata_core_cmp$summary,
  formula_cmp$summary,
  horizon_cmp$summary,
  threshold_cmp$summary
)

registry_value_differences <- bind_rows(
  metadata_core_cmp$value_differences,
  formula_cmp$value_differences,
  horizon_cmp$value_differences,
  threshold_cmp$value_differences
)

registry_key_presence_differences <- bind_rows(
  metadata_core_cmp$key_differences,
  formula_cmp$key_differences,
  horizon_cmp$key_differences,
  threshold_cmp$key_differences
)

# 🔴 Decide: whether Stage-7 refit is required ===============================
## 🟠 Assemble: reuse verdict and rationale ===============================
core_dataset_identical <- all(dataset_backbone_summary$backbone_equivalent_flag)
metadata_fit_input_equivalent <- get_equivalent_flag(registry_comparison_summary, "metadata_core_fit_input")
formula_fit_input_equivalent <- get_equivalent_flag(registry_comparison_summary, "formula_registry_fit_input")
horizon_reporting_equivalent <- get_equivalent_flag(registry_comparison_summary, "horizon_registry_reporting_contract")
threshold_reporting_equivalent <- get_equivalent_flag(registry_comparison_summary, "threshold_registry_reporting_contract")

fit_input_backbone_equivalent <- core_dataset_identical && metadata_fit_input_equivalent && formula_fit_input_equivalent
reporting_contract_equivalent <- horizon_reporting_equivalent && threshold_reporting_equivalent

reuse_decision <- make_decision_text(
  fit_input_backbone_equivalent = fit_input_backbone_equivalent,
  reporting_contract_equivalent = reporting_contract_equivalent
)

reuse_interpretation <- make_interpretation_text(
  fit_input_backbone_equivalent = fit_input_backbone_equivalent,
  reporting_contract_equivalent = reporting_contract_equivalent
)

overall_decision <- tibble(
  old_stage1_dir = normalize_existing_path(old_stage1_dir),
  new_stage1_dir = normalize_existing_path(new_stage1_dir),
  fit_input_backbone_equivalent = fit_input_backbone_equivalent,
  reporting_contract_equivalent = reporting_contract_equivalent,
  core_dataset_identical = core_dataset_identical,
  metadata_fit_input_equivalent = metadata_fit_input_equivalent,
  formula_fit_input_equivalent = formula_fit_input_equivalent,
  horizon_reporting_equivalent = horizon_reporting_equivalent,
  threshold_reporting_equivalent = threshold_reporting_equivalent,
  reuse_decision = reuse_decision,
  interpretation = reuse_interpretation
)

# 🔴 Export: comparison outputs and manifest ===============================
## 🟠 Save: CSV artifacts for review ===============================
overall_decision_file <- safe_write_csv(overall_decision, make_output_path("stage1_equivalence_overall_decision.csv"))
dataset_backbone_summary_file <- safe_write_csv(dataset_backbone_summary, make_output_path("stage1_equivalence_dataset_backbone_summary.csv"))
dataset_column_differences_file <- safe_write_csv(dataset_column_differences, make_output_path("stage1_equivalence_dataset_column_differences.csv"))
dataset_row_value_differences_file <- safe_write_csv(dataset_row_value_differences, make_output_path("stage1_equivalence_dataset_row_value_differences.csv"))
dataset_key_presence_differences_file <- safe_write_csv(dataset_key_presence_differences, make_output_path("stage1_equivalence_dataset_key_presence_differences.csv"))
registry_comparison_summary_file <- safe_write_csv(registry_comparison_summary, make_output_path("stage1_equivalence_registry_comparison_summary.csv"))
registry_value_differences_file <- safe_write_csv(registry_value_differences, make_output_path("stage1_equivalence_registry_value_differences.csv"))
registry_key_presence_differences_file <- safe_write_csv(registry_key_presence_differences, make_output_path("stage1_equivalence_registry_key_presence_differences.csv"))

export_manifest <- tibble(
  file_name = c(
    basename(overall_decision_file),
    basename(dataset_backbone_summary_file),
    basename(dataset_column_differences_file),
    basename(dataset_row_value_differences_file),
    basename(dataset_key_presence_differences_file),
    basename(registry_comparison_summary_file),
    basename(registry_value_differences_file),
    basename(registry_key_presence_differences_file)
  ),
  file_path = c(
    overall_decision_file,
    dataset_backbone_summary_file,
    dataset_column_differences_file,
    dataset_row_value_differences_file,
    dataset_key_presence_differences_file,
    registry_comparison_summary_file,
    registry_value_differences_file,
    registry_key_presence_differences_file
  ),
  description = c(
    "One-row final decision on whether existing Stage-7 results can be reused or whether refit is required.",
    "Dataset-level summary comparing old and new PNU/SNU/merged subject-level Stage-1 backbones.",
    "Column-level difference counts for the subject-level Stage-1 backbone.",
    "Row-level value differences for backbone variables among matched unique_person_id values.",
    "unique_person_id values present only in old or only in new Stage-1 outputs.",
    "Summary of shared fit-input and reporting-contract registry comparisons.",
    "Cell-level value differences for the compared registries.",
    "Key presence differences for the compared registries."
  )
) %>%
  mutate(
    exists = file.exists(file_path),
    size_bytes = ifelse(exists, file.info(file_path)$size, NA_real_)
  )

export_manifest_file <- safe_write_csv(export_manifest, make_output_path("stage1_equivalence_export_manifest.csv"))

message("✅ Stage-1 old/new equivalence check finished.")
message("Reuse decision: ", overall_decision$reuse_decision[[1]])
message("Interpretation: ", overall_decision$interpretation[[1]])
message("Manifest: ", export_manifest_file)
