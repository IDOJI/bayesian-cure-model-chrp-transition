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
  "SIMPLIFIED_RISK_INTEGRATION_EXPORT_PATH",
  unset = file.path(results_root, "integrated_risk_probabilities_by_horizon")
)

stage3_input_file <- file.path(results_root, "stage3_simple_km_yearly_estimates.csv")
stage5_input_file <- file.path(
  results_root,
  "stage5_individualized_nocure_lognormal_simple",
  "stage5_simple_lognormal_risk_estimates.csv"
)
stage7_input_file <- file.path(results_root, "stage7_simple_lognormal_mixture_cure_horizon_summary.csv")
stage8_input_file <- file.path(results_root, "stage8A_simple_bayesian_transition_only_horizon_summary.csv")

required_horizons <- 1:10
dataset_order <- c("PNU", "SNU", "merged", "merged_site_adjusted")
model_order <- c(
  "KM benchmark",
  "No-cure lognormal",
  "Frequentist mixture cure (lognormal)",
  "Frequentist mixture cure (lognormal, susceptible only)",
  "Bayesian transition-only cure",
  "Bayesian transition-only cure (susceptible only)"
)

# Load Packages -----------------------------------------------------------
required_packages <- c("readr", "dplyr", "tibble")
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

duplicate_dataset_rows <- function(df, source_dataset_name, target_dataset_name) {
  source_dataset_name <- as.character(source_dataset_name)
  target_dataset_name <- as.character(target_dataset_name)

  duplicated_rows <- df %>%
    filter(dataset_name == source_dataset_name) %>%
    mutate(dataset_name = target_dataset_name)

  bind_rows(df, duplicated_rows)
}

validate_integrated_table <- function(df) {
  required_columns <- c("dataset_name", "model_name", "horizon_year", "estimated_risk_probability")
  missing_columns <- setdiff(required_columns, names(df))

  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "Integrated table is missing required columns: %s",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (anyNA(df$dataset_name) || any(df$dataset_name == "")) {
    stop("Integrated table contains missing `dataset_name` values.", call. = FALSE)
  }

  if (anyNA(df$model_name) || any(df$model_name == "")) {
    stop("Integrated table contains missing `model_name` values.", call. = FALSE)
  }

  if (
    anyNA(df$horizon_year) ||
    !identical(sort(unique(as.integer(df$horizon_year))), as.integer(required_horizons))
  ) {
    stop("Integrated table must contain horizons 1 through 10.", call. = FALSE)
  }

  if (anyNA(df$estimated_risk_probability)) {
    stop("Integrated table contains missing `estimated_risk_probability` values.", call. = FALSE)
  }

  if (any(df$estimated_risk_probability < 0 | df$estimated_risk_probability > 1)) {
    stop("Integrated table contains risk probabilities outside [0, 1].", call. = FALSE)
  }

  duplicate_rows <- df %>%
    count(dataset_name, model_name, horizon_year, name = "row_count") %>%
    filter(row_count > 1L)

  if (nrow(duplicate_rows) > 0L) {
    stop("Integrated table contains duplicated dataset-model-horizon rows.", call. = FALSE)
  }

  invisible(TRUE)
}

standardize_stage3 <- function(path) {
  stage3_df <- read_csv_checked(path, "Stage 3 yearly estimate file")

  standard_df <- stage3_df %>%
    transmute(
      dataset_name = normalize_dataset_name(dataset),
      model_name = "KM benchmark",
      horizon_year = as.integer(horizon_year),
      estimated_risk_probability = clip_probability(estimated_risk_probability)
    )

  duplicate_dataset_rows(
    df = standard_df,
    source_dataset_name = "merged",
    target_dataset_name = "merged_site_adjusted"
  )
}

standardize_stage5 <- function(path) {
  stage5_df <- read_csv_checked(path, "Stage 5 risk estimate file")

  stage5_df %>%
    transmute(
      dataset_name = normalize_dataset_name(dataset_version_key),
      model_name = "No-cure lognormal",
      horizon_year = as.integer(horizon_year),
      estimated_risk_probability = clip_probability(estimated_risk_probability)
    )
}

standardize_stage7 <- function(path) {
  stage7_df <- read_csv_checked(path, "Stage 7 horizon summary file")

  bind_rows(
    stage7_df %>%
      transmute(
        dataset_name = normalize_dataset_name(dataset),
        model_name = "Frequentist mixture cure (lognormal)",
        horizon_year = as.integer(time_horizon_year),
        estimated_risk_probability = clip_probability(overall_risk_prob)
      ),
    stage7_df %>%
      transmute(
        dataset_name = normalize_dataset_name(dataset),
        model_name = "Frequentist mixture cure (lognormal, susceptible only)",
        horizon_year = as.integer(time_horizon_year),
        estimated_risk_probability = clip_probability(susceptible_only_risk_prob)
      )
  )
}

standardize_stage8 <- function(path) {
  stage8_df <- read_csv_checked(path, "Stage 8A horizon summary file")

  bind_rows(
    stage8_df %>%
      transmute(
        dataset_name = normalize_dataset_name(dataset),
        model_name = "Bayesian transition-only cure",
        horizon_year = as.integer(horizon_year),
        estimated_risk_probability = clip_probability(overall_risk_prob)
      ),
    stage8_df %>%
      transmute(
        dataset_name = normalize_dataset_name(dataset),
        model_name = "Bayesian transition-only cure (susceptible only)",
        horizon_year = as.integer(horizon_year),
        estimated_risk_probability = clip_probability(susceptible_only_risk_prob)
      )
  )
}

write_horizon_csv <- function(df, export_dir, horizon_year) {
  horizon_df <- df %>%
    filter(horizon_year == !!horizon_year) %>%
    mutate(
      dataset_name = factor(dataset_name, levels = dataset_order),
      model_name = factor(model_name, levels = model_order)
    ) %>%
    arrange(dataset_name, model_name) %>%
    mutate(
      dataset_name = as.character(dataset_name),
      model_name = as.character(model_name)
    ) %>%
    select(dataset_name, model_name, estimated_risk_probability)

  output_file <- file.path(
    export_dir,
    sprintf("risk_probabilities_year_%02d.csv", as.integer(horizon_year))
  )

  readr::write_csv(horizon_df, output_file)
  invisible(output_file)
}

# Execute Integration -----------------------------------------------------
integrated_risk_df <- bind_rows(
  standardize_stage3(stage3_input_file),
  standardize_stage5(stage5_input_file),
  standardize_stage7(stage7_input_file),
  standardize_stage8(stage8_input_file)
)

validate_integrated_table(integrated_risk_df)

for (year_value in required_horizons) {
  write_horizon_csv(integrated_risk_df, export_path, year_value)
}

message("Integrated risk probability files written to: ", export_path)
