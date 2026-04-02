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
dropbox_root <- Sys.getenv(
  "CHRPMIXTURE_DROPBOX_ROOT",
  unset = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model"
)
modeling_root <- Sys.getenv(
  "SIMPLE_MODELING_RESULTS_ROOT",
  unset = file.path(dropbox_root, "1.Modeling")
)
export_path <- Sys.getenv(
  "SIMPLE_MODEL_COMPARISON_TABLE_EXPORT_PATH",
  unset = file.path(
    dropbox_root,
    "2.Results Integration",
    "simple_model_probability_comparison_tables"
  )
)

stage3_input_file <- file.path(modeling_root, "1.KM", "stage3_simple_km_yearly_estimates.csv")
stage5_input_file <- file.path(
  modeling_root,
  "2.MLE No-Cure",
  "stage5_simple_lognormal_survival_estimates.csv"
)
stage7_input_file <- file.path(
  modeling_root,
  "3.MLE Mixture Cure",
  "stage7_simple_lognormal_mixture_cure_horizon_summary.csv"
)
stage8_input_candidates <- c(
  file.path(modeling_root, "4.Bayesian Mixture Cure", "bayesian_mixture_cure_horizon_summary.csv"),
  file.path(repo_root, "3.Results files", "stage8A_simple_bayesian_transition_only_horizon_summary.csv")
)

required_horizons <- 1:10
round_digits <- as.integer(Sys.getenv("SIMPLE_MODEL_COMPARISON_ROUND_DIGITS", unset = "6"))

# Load Packages -----------------------------------------------------------
required_packages <- c("readr", "dplyr", "tidyr", "tibble")
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
  library(tidyr)
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

read_csv_checked <- function(path, label) {
  assert_file_exists(path, label)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
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

dataset_group_from_variant <- function(dataset_variant) {
  dataset_variant <- normalize_dataset_variant(dataset_variant)

  dplyr::case_when(
    dataset_variant == "PNU" ~ "PNU",
    dataset_variant == "SNU" ~ "SNU",
    dataset_variant %in% c("merged", "merged_site_adjusted") ~ "Merged",
    TRUE ~ dataset_variant
  )
}

clip_probability <- function(x) {
  pmin(pmax(as.numeric(x), 0), 1)
}

round_probability <- function(x) {
  round(as.numeric(x), digits = round_digits)
}

validate_long_table <- function(df) {
  required_columns <- c(
    "dataset_variant",
    "dataset_group",
    "model_family",
    "model_scope",
    "probability_type",
    "horizon_year",
    "estimate_probability"
  )
  missing_columns <- setdiff(required_columns, names(df))

  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "Integrated long table is missing required columns: %s",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  duplicate_rows <- df %>%
    count(
      dataset_variant,
      model_family,
      model_scope,
      probability_type,
      horizon_year,
      name = "row_count"
    ) %>%
    filter(row_count > 1L)

  if (nrow(duplicate_rows) > 0L) {
    stop("Integrated long table contains duplicated rows.", call. = FALSE)
  }

  missing_horizons <- df %>%
    count(
      dataset_variant,
      dataset_group,
      model_family,
      model_scope,
      probability_type,
      name = "n_horizons"
    ) %>%
    filter(n_horizons != length(required_horizons))

  if (nrow(missing_horizons) > 0L) {
    stop("At least one dataset-model-probability series is missing yearly horizons 1-10.", call. = FALSE)
  }

  unexpected_horizons <- setdiff(sort(unique(df$horizon_year)), required_horizons)
  if (length(unexpected_horizons) > 0L) {
    stop("Unexpected yearly horizons detected in integrated long table.", call. = FALSE)
  }

  if (anyNA(df$estimate_probability)) {
    stop("Integrated long table contains missing probabilities.", call. = FALSE)
  }

  invisible(TRUE)
}

model_row_label <- function(dataset_group, dataset_variant, model_family, model_scope) {
  if (identical(model_family, "KM benchmark")) {
    base_label <- "KM benchmark"
  } else if (identical(model_family, "Frequentist no-cure lognormal")) {
    base_label <- "Frequentist no-cure lognormal"
  } else if (identical(model_family, "Frequentist mixture cure")) {
    base_label <- if (identical(model_scope, "susceptible")) {
      "Frequentist mixture cure (susceptible only)"
    } else {
      "Frequentist mixture cure (overall)"
    }
  } else if (identical(model_family, "Bayesian transition-only cure")) {
    base_label <- if (identical(model_scope, "susceptible")) {
      "Bayesian transition-only cure (susceptible only)"
    } else {
      "Bayesian transition-only cure (overall)"
    }
  } else {
    base_label <- model_family
  }

  if (!identical(dataset_group, "Merged")) {
    return(base_label)
  }

  variant_suffix <- if (identical(dataset_variant, "merged_site_adjusted")) {
    "[merged(site-adj)]"
  } else {
    "[merged]"
  }

  paste(base_label, variant_suffix)
}

model_row_order <- function(dataset_group, dataset_variant, model_family, model_scope) {
  if (!identical(dataset_group, "Merged")) {
    return(dplyr::case_when(
      model_family == "KM benchmark" ~ 1,
      model_family == "Frequentist no-cure lognormal" ~ 2,
      model_family == "Frequentist mixture cure" & model_scope == "overall" ~ 3,
      model_family == "Frequentist mixture cure" & model_scope == "susceptible" ~ 4,
      model_family == "Bayesian transition-only cure" & model_scope == "overall" ~ 5,
      model_family == "Bayesian transition-only cure" & model_scope == "susceptible" ~ 6,
      TRUE ~ 999
    ))
  }

  dplyr::case_when(
    model_family == "KM benchmark" & dataset_variant == "merged" ~ 1,
    model_family == "Frequentist no-cure lognormal" & dataset_variant == "merged" ~ 2,
    model_family == "Frequentist no-cure lognormal" & dataset_variant == "merged_site_adjusted" ~ 3,
    model_family == "Frequentist mixture cure" & model_scope == "overall" & dataset_variant == "merged" ~ 4,
    model_family == "Frequentist mixture cure" & model_scope == "susceptible" & dataset_variant == "merged" ~ 5,
    model_family == "Frequentist mixture cure" & model_scope == "overall" & dataset_variant == "merged_site_adjusted" ~ 6,
    model_family == "Frequentist mixture cure" & model_scope == "susceptible" & dataset_variant == "merged_site_adjusted" ~ 7,
    model_family == "Bayesian transition-only cure" & model_scope == "overall" & dataset_variant == "merged" ~ 8,
    model_family == "Bayesian transition-only cure" & model_scope == "susceptible" & dataset_variant == "merged" ~ 9,
    model_family == "Bayesian transition-only cure" & model_scope == "overall" & dataset_variant == "merged_site_adjusted" ~ 10,
    model_family == "Bayesian transition-only cure" & model_scope == "susceptible" & dataset_variant == "merged_site_adjusted" ~ 11,
    TRUE ~ 999
  )
}

build_comparison_table <- function(long_df, dataset_group, probability_type) {
  comparison_df <- long_df %>%
    filter(
      dataset_group == !!dataset_group,
      probability_type == !!probability_type
    ) %>%
    mutate(
      row_label = mapply(
        model_row_label,
        dataset_group = dataset_group,
        dataset_variant = dataset_variant,
        model_family = model_family,
        model_scope = model_scope,
        USE.NAMES = FALSE
      ),
      row_order = mapply(
        model_row_order,
        dataset_group = dataset_group,
        dataset_variant = dataset_variant,
        model_family = model_family,
        model_scope = model_scope,
        USE.NAMES = FALSE
      )
    )

  duplicate_rows <- comparison_df %>%
    count(row_label, horizon_year, name = "row_count") %>%
    filter(row_count > 1L)

  if (nrow(duplicate_rows) > 0L) {
    stop(
      sprintf(
        "Duplicated rows detected while building %s / %s comparison table.",
        dataset_group,
        probability_type
      ),
      call. = FALSE
    )
  }

  comparison_df %>%
    transmute(
      model = row_label,
      row_order = row_order,
      horizon_year = horizon_year,
      estimate_probability = round_probability(estimate_probability)
    ) %>%
    tidyr::pivot_wider(
      names_from = horizon_year,
      values_from = estimate_probability,
      names_glue = "{horizon_year} year"
    ) %>%
    arrange(row_order, model) %>%
    select(-row_order)
}

write_comparison_table <- function(table_df, export_dir, dataset_group, probability_type) {
  dataset_stub <- dplyr::case_when(
    dataset_group == "PNU" ~ "pnu",
    dataset_group == "SNU" ~ "snu",
    dataset_group == "Merged" ~ "merged",
    TRUE ~ tolower(dataset_group)
  )

  probability_stub <- dplyr::case_when(
    probability_type == "survival" ~ "survival_probability",
    probability_type == "risk" ~ "cumulative_risk_probability",
    TRUE ~ probability_type
  )

  output_file <- file.path(
    export_dir,
    sprintf("%s_%s_comparison_table.csv", dataset_stub, probability_stub)
  )

  readr::write_csv(table_df, output_file)
  invisible(output_file)
}

standardize_stage3 <- function(path) {
  stage3_df <- read_csv_checked(path, "Stage 3 yearly estimate file")

  bind_rows(
    stage3_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "KM benchmark",
        model_scope = "overall",
        probability_type = "survival",
        horizon_year = as.integer(horizon_year),
        estimate_probability = clip_probability(estimated_survival_probability)
      ),
    stage3_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "KM benchmark",
        model_scope = "overall",
        probability_type = "risk",
        horizon_year = as.integer(horizon_year),
        estimate_probability = clip_probability(estimated_risk_probability)
      )
  )
}

standardize_stage5 <- function(path) {
  stage5_df <- read_csv_checked(path, "Stage 5 yearly estimate file")

  bind_rows(
    stage5_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset_version_key),
        dataset_group = dataset_group_from_variant(dataset_version_key),
        model_family = "Frequentist no-cure lognormal",
        model_scope = "overall",
        probability_type = "survival",
        horizon_year = as.integer(horizon_year),
        estimate_probability = clip_probability(estimated_survival_probability)
      ),
    stage5_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset_version_key),
        dataset_group = dataset_group_from_variant(dataset_version_key),
        model_family = "Frequentist no-cure lognormal",
        model_scope = "overall",
        probability_type = "risk",
        horizon_year = as.integer(horizon_year),
        estimate_probability = clip_probability(estimated_risk_probability)
      )
  )
}

standardize_stage7 <- function(path) {
  stage7_df <- read_csv_checked(path, "Stage 7 horizon summary file")

  bind_rows(
    stage7_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "Frequentist mixture cure",
        model_scope = "overall",
        probability_type = "survival",
        horizon_year = as.integer(time_horizon_year),
        estimate_probability = clip_probability(overall_survival_prob)
      ),
    stage7_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "Frequentist mixture cure",
        model_scope = "overall",
        probability_type = "risk",
        horizon_year = as.integer(time_horizon_year),
        estimate_probability = clip_probability(overall_risk_prob)
      ),
    stage7_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "Frequentist mixture cure",
        model_scope = "susceptible",
        probability_type = "survival",
        horizon_year = as.integer(time_horizon_year),
        estimate_probability = clip_probability(susceptible_only_survival_prob)
      ),
    stage7_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "Frequentist mixture cure",
        model_scope = "susceptible",
        probability_type = "risk",
        horizon_year = as.integer(time_horizon_year),
        estimate_probability = clip_probability(susceptible_only_risk_prob)
      )
  )
}

standardize_stage8 <- function(path) {
  stage8_df <- read_csv_checked(path, "Stage 8 horizon summary file")

  bind_rows(
    stage8_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "Bayesian transition-only cure",
        model_scope = "overall",
        probability_type = "survival",
        horizon_year = as.integer(horizon_year),
        estimate_probability = clip_probability(overall_survival_prob)
      ),
    stage8_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "Bayesian transition-only cure",
        model_scope = "overall",
        probability_type = "risk",
        horizon_year = as.integer(horizon_year),
        estimate_probability = clip_probability(overall_risk_prob)
      ),
    stage8_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "Bayesian transition-only cure",
        model_scope = "susceptible",
        probability_type = "survival",
        horizon_year = as.integer(horizon_year),
        estimate_probability = clip_probability(susceptible_only_survival_prob)
      ),
    stage8_df %>%
      transmute(
        dataset_variant = normalize_dataset_variant(dataset),
        dataset_group = dataset_group_from_variant(dataset),
        model_family = "Bayesian transition-only cure",
        model_scope = "susceptible",
        probability_type = "risk",
        horizon_year = as.integer(horizon_year),
        estimate_probability = clip_probability(susceptible_only_risk_prob)
      )
  )
}

# Execute Integration -----------------------------------------------------
resolved_stage8_input_file <- resolve_existing_path(
  stage8_input_candidates,
  "Stage 8 Bayesian horizon summary file"
)

metadata_df <- tibble::tibble(
  input_label = c("stage3", "stage5", "stage7", "stage8"),
  input_path = c(
    normalizePath(stage3_input_file, winslash = "/", mustWork = TRUE),
    normalizePath(stage5_input_file, winslash = "/", mustWork = TRUE),
    normalizePath(stage7_input_file, winslash = "/", mustWork = TRUE),
    resolved_stage8_input_file
  )
)

integrated_long_df <- bind_rows(
  standardize_stage3(stage3_input_file),
  standardize_stage5(stage5_input_file),
  standardize_stage7(stage7_input_file),
  standardize_stage8(resolved_stage8_input_file)
) %>%
  filter(horizon_year %in% required_horizons) %>%
  arrange(
    factor(dataset_group, levels = c("PNU", "SNU", "Merged")),
    factor(dataset_variant, levels = c("PNU", "SNU", "merged", "merged_site_adjusted")),
    model_family,
    model_scope,
    probability_type,
    horizon_year
  )

validate_long_table(integrated_long_df)

readr::write_csv(
  metadata_df,
  file.path(export_path, "input_file_registry.csv")
)

readr::write_csv(
  integrated_long_df %>%
    mutate(estimate_probability = round_probability(estimate_probability)),
  file.path(export_path, "integrated_model_probabilities_long.csv")
)

comparison_specs <- tibble::tibble(
  dataset_group = c("PNU", "PNU", "SNU", "SNU", "Merged", "Merged"),
  probability_type = c("survival", "risk", "survival", "risk", "survival", "risk")
)

output_registry <- comparison_specs %>%
  rowwise() %>%
  mutate(
    output_file = write_comparison_table(
      table_df = build_comparison_table(integrated_long_df, dataset_group, probability_type),
      export_dir = export_path,
      dataset_group = dataset_group,
      probability_type = probability_type
    )
  ) %>%
  ungroup()

readr::write_csv(
  output_registry,
  file.path(export_path, "output_table_registry.csv")
)

message("Comparison tables written to: ", export_path)
