# Configure Paths ---------------------------------------------------------
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

  stop(
    "Could not locate the repository root.",
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
results_root_default <- file.path(repo_root, "3.Results files")
block1_root_default <- "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/1.Block1"

source_data_file <- Sys.getenv(
  "STAGE7_SIMPLE_DATA_FILE",
  unset = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/0.Data/2.Preprocessed data/Preprocessed_Merged_PNUH_SNUH_Data.csv"
)
export_path <- Sys.getenv(
  "STAGE7_SIMPLE_EXPORT_PATH",
  unset = file.path(block1_root_default, "3.MLE Mixture Cure")
)
pnu_site_label <- Sys.getenv("STAGE7_SIMPLE_PNU_SITE_LABEL", unset = "PNU")
snu_site_label <- Sys.getenv("STAGE7_SIMPLE_SNU_SITE_LABEL", unset = "SNU")
reuse_existing_fit_rds <- identical(
  toupper(trimws(Sys.getenv("STAGE7_SIMPLE_REUSE_EXISTING_FIT_RDS", unset = "TRUE"))),
  "TRUE"
)

horizon_years <- 1:10
curve_horizon_max_year <- 10
curve_step_year <- 0.05
detail_tick_step_year <- 0.5
time_origin_epsilon_year <- 1e-10
main_risk_scale <- "transition_only_main"
model_id_value <- "frequentist_mixture_cure_lognormal"
stage7_ci_draws <- as.integer(Sys.getenv("STAGE7_SIMPLE_CI_DRAWS", unset = "4000"))
output_prefix <- Sys.getenv(
  "STAGE7_SIMPLE_OUTPUT_PREFIX",
  unset = "mle_mixture_cure_lognormal"
)
analysis_label <- Sys.getenv(
  "STAGE7_SIMPLE_ANALYSIS_LABEL",
  unset = "Frequentist log-normal mixture cure"
)
plot_width_in <- 10
plot_height_in <- 6
plot_dpi <- 320
detail_plot_width_in <- 14
detail_plot_height_in <- 8

horizon_summary_file <- file.path(
  export_path,
  paste0(output_prefix, "_horizon_summary.csv")
)
fit_summary_file <- file.path(
  export_path,
  paste0(output_prefix, "_fit_summary.csv")
)
plot_source_file <- file.path(
  export_path,
  paste0(output_prefix, "_plot_source.csv")
)
fit_object_rds_file <- file.path(
  export_path,
  paste0(output_prefix, "_fitted_objects.rds")
)
plot_rds_file <- file.path(
  export_path,
  paste0(output_prefix, "_plot_objects.rds")
)
detail_annotation_file <- file.path(
  export_path,
  paste0(output_prefix, "_detail_annotation_table.csv")
)
latency_plot_source_file <- file.path(
  export_path,
  paste0(output_prefix, "_latency_plot_source.csv")
)
latency_median_distribution_png_file <- file.path(
  export_path,
  paste0(output_prefix, "_uncured_latency_median_distribution.png")
)
latency_mean_distribution_png_file <- file.path(
  export_path,
  paste0(output_prefix, "_uncured_latency_mean_distribution.png")
)
latency_summary_png_file <- file.path(
  export_path,
  paste0(output_prefix, "_uncured_latency_summary_comparison.png")
)
latency_aft_effect_png_file <- file.path(
  export_path,
  paste0(output_prefix, "_latency_aft_time_ratio_plot.png")
)
standard_error_registry_file <- file.path(
  export_path,
  paste0(output_prefix, "_standard_error_table_registry.csv")
)
standard_error_long_file <- file.path(
  export_path,
  paste0(output_prefix, "_standard_error_long.csv")
)

dataset_model_registry <- tibble::tibble(
  dataset = c("PNU", "SNU", "merged_no_site", "merged_site_adjusted"),
  dataset_label = c("PNU", "SNU", "merged", "merged (site-adjusted)"),
  source_dataset = c("PNU", "SNU", "merged", "merged"),
  formula_name = c("base", "base", "base", "site_added"),
  site_adjustment_flag = c(FALSE, FALSE, FALSE, TRUE)
)

dataset_palette <- c(
  "PNU" = "#1B4332",
  "SNU" = "#2A6F97",
  "merged_no_site" = "#C1666B",
  "merged_site_adjusted" = "#B8860B"
)

plot_dataset_label_lookup <- c(
  "PNU" = "PNU",
  "SNU" = "SNU",
  "merged_no_site" = "Merged",
  "merged_site_adjusted" = "Merged (site-adjusted)"
)

# Load Packages -----------------------------------------------------------
required_packages <- c("readr", "dplyr", "tibble", "ggplot2", "scales", "survival")
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
  library(ggplot2)
  library(scales)
  library(survival)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# Define Helpers ----------------------------------------------------------
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

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

cleanup_existing_stage7_outputs <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    return(invisible(character()))
  }

  sub_dir <- file.path(output_dir, "sub")
  if (dir.exists(sub_dir)) {
    unlink(sub_dir, recursive = TRUE, force = TRUE)
  }

  stale_files <- list.files(
    output_dir,
    pattern = "^(stage7_simple_lognormal_mixture_cure_|mle_mixture_cure_lognormal_).*\\.(csv|rds|png)$",
    full.names = TRUE
  )

  if (length(stale_files) == 0L) {
    return(invisible(character()))
  }

  removed_flag <- file.remove(stale_files)
  if (any(!removed_flag)) {
    stop(
      sprintf(
        "Failed to remove stale frequentist cure outputs: %s",
        paste(basename(stale_files[!removed_flag]), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(stale_files[removed_flag])
}

organize_stage7_outputs <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    return(invisible(character()))
  }

  sub_dir <- file.path(output_dir, "sub")
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)

  keep_files <- c(
    "readme.md",
    paste0(output_prefix, "_fitted_objects.rds"),
    paste0(output_prefix, "_horizon_summary.csv"),
    paste0(output_prefix, "_fit_summary.csv"),
    paste0(output_prefix, "_plot_source.csv"),
    paste0(output_prefix, "_detail_annotation_table.csv"),
    paste0(output_prefix, "_latency_plot_source.csv"),
    paste0(output_prefix, "_standard_error_table_registry.csv"),
    paste0(output_prefix, "_standard_error_long.csv"),
    paste0(output_prefix, "_uncured_latency_summary_comparison.png"),
    paste0(output_prefix, "_latency_aft_time_ratio_plot.png")
  )

  top_level_files <- list.files(
    output_dir,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE
  )
  top_level_files <- top_level_files[file.info(top_level_files)$isdir %in% FALSE]

  move_files <- top_level_files[!(basename(top_level_files) %in% keep_files)]

  if (length(move_files) == 0L) {
    return(invisible(character()))
  }

  moved_flag <- file.rename(move_files, file.path(sub_dir, basename(move_files)))
  if (any(!moved_flag)) {
    stop(
      sprintf(
        "Failed to move frequentist cure support files into `sub`: %s",
        paste(basename(move_files[!moved_flag]), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(file.path(sub_dir, basename(move_files[moved_flag])))
}

write_mle_mixture_cure_readme <- function(output_dir) {
  readme_lines <- c(
    "이미지와 함께 아래 설명을 주면 ChatGPT가 훨씬 안정적으로 해석할 수 있습니다.",
    "",
    "**공통 설명**",
    "- 이 그림들은 `frequentist lognormal mixture-cure model` 결과입니다.",
    "- `survival`과 `risk`는 개별 예측값의 코호트 평균입니다.",
    "- `uncured latency`는 전체 코호트가 아니라 `susceptible/uncured subgroup` 기준입니다.",
    "- `latency AFT effect`는 HR이 아니라 `AFT time ratio`입니다.",
    "- `mean latency`는 lognormal tail 영향으로 매우 커질 수 있으므로 `median latency`를 우선 해석합니다.",
    "- `Merged`와 `Merged (site-adjusted)`는 서로 다른 모델 버전입니다.",
    "- PNU 그림의 수직선은 PNU의 최대 관측 추적 종료 시점입니다.",
    "",
    "**파일별 설명**",
    sprintf("- `%s_uncured_latency_median_distribution.png`", output_prefix),
    "  - 각 점은 개인별 `uncured subgroup` 기준 `median time to transition`입니다.",
    "  - boxplot은 dataset별 분포를 요약합니다.",
    "  - y축은 로그축입니다.",
    "",
    sprintf("- `%s_uncured_latency_mean_distribution.png`", output_prefix),
    "  - 각 점은 개인별 `uncured subgroup` 기준 `mean time to transition`입니다.",
    "  - lognormal tail 때문에 큰 값이 나올 수 있습니다.",
    "  - 해석은 보조적으로 하고, median plot을 더 중시합니다.",
    "",
    sprintf("- `%s_uncured_latency_summary_comparison.png`", output_prefix),
    "  - dataset별 대표 latency를 비교한 요약 그림입니다.",
    "  - `weighted mean median latency`와 `weighted mean mean latency`를 보여줍니다.",
    "  - 가중치는 각 개인의 `susceptible_fraction`입니다.",
    "",
    sprintf("- `%s_latency_aft_time_ratio_plot.png`", output_prefix),
    "  - uncured subgroup의 latency 부분 회귀효과를 `time ratio`로 표현한 그림입니다.",
    "  - `time ratio > 1`이면 전이까지 시간이 길어짐, `< 1`이면 짧아짐을 뜻합니다.",
    "  - 이것은 hazard ratio가 아닙니다.",
    "",
    sprintf("- `%s_estimated_survival_curve*.png`", output_prefix),
    "  - mixture-cure model 기반 전체 생존확률 곡선입니다.",
    "  - x축은 cohort entry 이후 연수입니다.",
    "",
    sprintf("- `%s_estimated_risk_curve*.png`", output_prefix),
    "  - 위 생존확률의 `1 - survival`인 누적 위험확률 곡선입니다.",
    "",
    sprintf("- `%s_*with_counts*.png`", output_prefix),
    "  - 위 곡선 아래에 6개월 간격 `n at risk`와 `cumulative transitions`를 bar plot으로 함께 표시한 그림입니다.",
    "",
    sprintf("- `%s_standard_error_table_registry.csv`", output_prefix),
    "  - 어떤 export table에 `sd/se` 계열 불확실성 컬럼이 포함됐는지 정리한 registry입니다.",
    "",
    sprintf("- `%s_standard_error_long.csv`", output_prefix),
    "  - `sd/se` 계열 불확실성 값을 long format으로 펼친 table입니다."
  )

  writeLines(readme_lines, con = file.path(output_dir, "readme.md"), useBytes = TRUE)
  invisible(file.path(output_dir, "readme.md"))
}

coerce_numeric_text <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
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
    "unique_person_id",
    "dataset",
    "dataset_label",
    "plot_dataset_label",
    "source_dataset",
    "model_id",
    "formula_name",
    "formula_id",
    "formula_label",
    "incidence_rhs",
    "latency_rhs",
    "site_adjustment_flag",
    "site_in_incidence",
    "site_in_latency",
    "risk_scale",
    "time_year",
    "time_horizon_year",
    "horizon_year",
    "component_type",
    "interval_type",
    "latency_term",
    "metric",
    "metric_label",
    "effect_model"
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

standardize_known_site_labels <- function(df, pnu_label, snu_label) {
  if (!"site" %in% names(df)) {
    stop("Merged input must contain column `site`.", call. = FALSE)
  }

  df %>%
    mutate(
      site = trimws(as.character(site)),
      site = dplyr::case_when(
        toupper(site) == toupper(pnu_label) ~ pnu_label,
        toupper(site) == toupper(snu_label) ~ snu_label,
        TRUE ~ site
      )
    )
}

validate_source_data <- function(df, pnu_label, snu_label) {
  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "Merged preprocessed data is missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  observed_sites <- sort(unique(trimws(as.character(df$site))))
  expected_sites <- c(pnu_label, snu_label)
  unexpected_sites <- setdiff(observed_sites, expected_sites)
  if (length(unexpected_sites) > 0L) {
    stop(
      sprintf(
        "Unexpected site labels found in merged input: %s",
        paste(unexpected_sites, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

split_site_dataset <- function(merged_df, site_label, dataset_name) {
  out <- merged_df %>%
    filter(toupper(trimws(as.character(site))) == toupper(site_label)) %>%
    mutate(site = site_label)

  if (nrow(out) == 0L) {
    stop(
      sprintf("[%s] No rows found in merged input for site label `%s`.", dataset_name, site_label),
      call. = FALSE
    )
  }

  out
}

prepare_backbone_dataset <- function(df, dataset_name, site_mode = c("single", "merged")) {
  site_mode <- match.arg(site_mode)

  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("[%s] Missing required columns: %s", dataset_name, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  out <- df %>%
    mutate(
      id = trimws(as.character(id)),
      site = trimws(as.character(site)),
      sex_num = as.integer(coerce_numeric_text(sex_num)),
      age_exact_entry = coerce_numeric_text(age_exact_entry),
      days_followup = coerce_numeric_text(days_followup),
      status_num = as.integer(coerce_numeric_text(status_num))
    )

  if (nrow(out) == 0L) {
    stop(sprintf("[%s] Dataset has zero rows after loading.", dataset_name), call. = FALSE)
  }

  if (anyNA(out[, required_cols])) {
    stop(sprintf("[%s] Missing values detected in required columns.", dataset_name), call. = FALSE)
  }

  if (any(out$id == "", na.rm = TRUE)) {
    stop(sprintf("[%s] Blank `id` values detected.", dataset_name), call. = FALSE)
  }

  if (any(out$site == "", na.rm = TRUE)) {
    stop(sprintf("[%s] Blank `site` values detected.", dataset_name), call. = FALSE)
  }

  if (any(!out$sex_num %in% c(0L, 1L), na.rm = TRUE)) {
    stop(sprintf("[%s] `sex_num` must be coded as 0/1 only.", dataset_name), call. = FALSE)
  }

  if (any(!out$status_num %in% c(0L, 1L, 2L), na.rm = TRUE)) {
    stop(sprintf("[%s] `status_num` must be coded as 0/1/2 only.", dataset_name), call. = FALSE)
  }

  if (any(out$days_followup < 0, na.rm = TRUE)) {
    stop(sprintf("[%s] Negative `days_followup` values detected.", dataset_name), call. = FALSE)
  }

  n_site_levels <- dplyr::n_distinct(out$site)
  if (site_mode == "single" && n_site_levels != 1L) {
    stop(sprintf("[%s] Single-cohort input must contain exactly one site.", dataset_name), call. = FALSE)
  }
  if (site_mode == "merged" && n_site_levels < 2L) {
    stop(sprintf("[%s] Merged input must contain at least two site levels.", dataset_name), call. = FALSE)
  }

  age_sd <- stats::sd(out$age_exact_entry)
  if (is.na(age_sd) || age_sd <= 0) {
    stop(sprintf("[%s] `age_exact_entry` must have positive standard deviation.", dataset_name), call. = FALSE)
  }

  age_mean <- mean(out$age_exact_entry)

  out %>%
    mutate(
      unique_person_id = paste(site, id, sep = "_"),
      time_year = days_followup / 365.25,
      event_main = as.integer(status_num == 1L),
      right_censor_flag = as.integer(status_num == 0L),
      remission_flag = as.integer(status_num == 2L),
      censor_main = as.integer(status_num %in% c(0L, 2L)),
      sex_label = factor(if_else(sex_num == 0L, "Female", "Male"), levels = c("Female", "Male")),
      age_s = (age_exact_entry - age_mean) / (2 * age_sd)
    )
}

build_analysis_datasets_from_source <- function(source_file, pnu_label, snu_label) {
  assert_file_exists(source_file, "Merged preprocessed data CSV")
  source_df <- readr::read_csv(
    source_file,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
  source_df <- standardize_known_site_labels(source_df, pnu_label, snu_label)
  validate_source_data(source_df, pnu_label, snu_label)

  analysis_datasets <- list(
    PNU = prepare_backbone_dataset(
      split_site_dataset(source_df, pnu_label, "PNU"),
      dataset_name = "PNU",
      site_mode = "single"
    ),
    SNU = prepare_backbone_dataset(
      split_site_dataset(source_df, snu_label, "SNU"),
      dataset_name = "SNU",
      site_mode = "single"
    ),
    merged = prepare_backbone_dataset(
      source_df,
      dataset_name = "merged",
      site_mode = "merged"
    )
  )

  if (!all(c("PNU", "SNU", "merged") %in% names(analysis_datasets))) {
    stop("Analysis dataset bundle must contain `PNU`, `SNU`, and `merged`.", call. = FALSE)
  }

  analysis_datasets
}

summarize_observed_followup <- function(analysis_datasets) {
  bind_rows(lapply(names(analysis_datasets), function(dataset_name) {
    tibble(
      source_dataset = dataset_name,
      max_observed_followup_year = max(as.numeric(analysis_datasets[[dataset_name]]$time_year), na.rm = TRUE)
    )
  }))
}

make_formula_registry <- function() {
  formula_tbl <- bind_rows(
    tibble::tibble(
      dataset = c("PNU", "PNU", "SNU", "SNU"),
      formula_name = c("base", "interaction", "base", "interaction"),
      formula_label = c("Base", "Interaction", "Base", "Interaction"),
      formula_rhs = c(
        "age_s + sex_num",
        "age_s + sex_num + age_s:sex_num",
        "age_s + sex_num",
        "age_s + sex_num + age_s:sex_num"
      )
    ),
    tibble::tibble(
      dataset = c("merged", "merged", "merged", "merged"),
      formula_name = c("base", "interaction", "site_added", "site_interaction"),
      formula_label = c("Base", "Interaction", "Site-added", "Site + interaction"),
      formula_rhs = c(
        "age_s + sex_num",
        "age_s + sex_num + age_s:sex_num",
        "age_s + sex_num + site",
        "age_s + sex_num + age_s:sex_num + site"
      )
    )
  )

  formula_tbl %>%
    mutate(
      formula_id = paste(dataset, formula_name, sep = "__"),
      formula_full = paste("~", formula_rhs),
      uses_site = grepl("\\bsite\\b", formula_rhs),
      uses_age_sex_interaction = grepl("age_s:sex_num", formula_rhs, fixed = TRUE),
      site_branch = if_else(uses_site, "site_adjusted", "site_free"),
      interaction_branch = if_else(uses_age_sex_interaction, "age_sex_interaction", "no_age_sex_interaction")
    ) %>%
    select(
      dataset,
      formula_id,
      formula_name,
      formula_label,
      formula_rhs,
      formula_full,
      uses_site,
      uses_age_sex_interaction,
      site_branch,
      interaction_branch
    )
}

normalize_dataset_label <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    toupper(x) == "PNU" ~ "PNU",
    toupper(x) == "SNU" ~ "SNU",
    tolower(x) == "merged" ~ "merged",
    TRUE ~ x
  )
}

trim_rhs <- function(rhs) {
  rhs <- gsub("\\s+", " ", trimws(as.character(rhs)))
  rhs[is.na(rhs)] <- NA_character_
  rhs
}

clamp_prob <- function(x, eps = 1e-12) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

weighted_mean_safe <- function(x, w) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  sum(x[ok] * w[ok]) / sum(w[ok])
}

capture_with_warnings <- function(expr) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  list(
    value = if (inherits(value, "error")) NULL else value,
    warnings = warnings,
    error_message = if (inherits(value, "error")) conditionMessage(value) else NA_character_
  )
}

make_rhs_formula <- function(rhs) {
  rhs <- trim_rhs(rhs)
  if (is.na(rhs) || rhs == "") {
    stats::as.formula("~ 1")
  } else {
    stats::as.formula(paste("~", rhs))
  }
}

make_surv_formula <- function(rhs) {
  rhs <- trim_rhs(rhs)
  if (is.na(rhs) || rhs == "") {
    stats::as.formula("survival::Surv(time_year_model, event_main) ~ 1")
  } else {
    stats::as.formula(paste("survival::Surv(time_year_model, event_main) ~", rhs))
  }
}

model_matrix_from_rhs <- function(df, rhs) {
  stats::model.matrix(make_rhs_formula(rhs), data = df)
}

lognormal_surv_density <- function(time, lp, log_sigma) {
  time <- pmax(as.numeric(time), time_origin_epsilon_year)
  lp <- as.numeric(lp)
  sigma <- exp(as.numeric(log_sigma))
  z <- (log(time) - lp) / sigma
  surv <- 1 - stats::pnorm(z)
  dens <- stats::dnorm(z) / (sigma * time)
  list(
    surv = pmin(pmax(surv, 0), 1),
    dens = pmax(dens, 1e-300)
  )
}

summarize_draws <- function(x) {
  tibble(
    sd = stats::sd(x, na.rm = TRUE),
    q025 = unname(stats::quantile(x, probs = 0.025, na.rm = TRUE, names = FALSE)),
    q50 = unname(stats::quantile(x, probs = 0.5, na.rm = TRUE, names = FALSE)),
    q975 = unname(stats::quantile(x, probs = 0.975, na.rm = TRUE, names = FALSE))
  )
}

summarize_draws_with_prefix <- function(x, prefix) {
  summarize_draws(x) %>%
    rename_with(~ paste0(prefix, "_", .x), everything())
}

make_psd_matrix <- function(x, eigen_floor = 1e-10) {
  eig <- eigen(x, symmetric = TRUE)
  eig$values[eig$values < eigen_floor] <- eigen_floor
  eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
}

invert_hessian <- function(hessian_mat, eigen_floor = 1e-8) {
  sym_hessian <- 0.5 * (hessian_mat + t(hessian_mat))
  eig <- eigen(sym_hessian, symmetric = TRUE)
  eig$values[eig$values < eigen_floor] <- eigen_floor
  eig$vectors %*% diag(1 / eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
}

draw_mvn <- function(n_draws, mu, sigma, seed_value) {
  sigma_psd <- make_psd_matrix(sigma)
  eig <- eigen(sigma_psd, symmetric = TRUE)
  transform_mat <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow = length(eig$values))

  set.seed(seed_value)
  z <- matrix(stats::rnorm(length(mu) * n_draws), nrow = length(mu), ncol = n_draws)
  draws <- matrix(mu, nrow = length(mu), ncol = n_draws) + transform_mat %*% z
  t(draws)
}

build_prediction_draw_bundle <- function(df, fit, n_draws, seed_value) {
  pred_data <- tibble::as_tibble(df) %>%
    mutate(site = factor(as.character(site), levels = fit$site_levels))

  X_inc <- model_matrix_from_rhs(pred_data, fit$incidence_rhs)
  X_lat <- model_matrix_from_rhs(pred_data, fit$latency_rhs)
  y_time <- pmax(as.numeric(pred_data$time_year_model), time_origin_epsilon_year)
  y_event <- as.integer(pred_data$event_main)

  neg_loglik <- function(par) {
    gamma <- par[seq_len(ncol(X_inc))]
    beta <- par[ncol(X_inc) + seq_len(ncol(X_lat))]
    log_sigma <- par[[length(par)]]

    uncured_prob_inner <- clamp_prob(plogis(drop(X_inc %*% gamma)))
    lp_lat <- drop(X_lat %*% beta)
    dens_surv <- lognormal_surv_density(y_time, lp_lat, log_sigma)
    surv_uncured <- pmax(dens_surv$surv, 1e-300)
    dens_uncured <- pmax(dens_surv$dens, 1e-300)

    loglik <- sum(
      y_event * (log(uncured_prob_inner) + log(dens_uncured)) +
        (1 - y_event) * log((1 - uncured_prob_inner) + uncured_prob_inner * surv_uncured)
    )

    if (!is.finite(loglik)) {
      return(1e12)
    }

    -loglik
  }

  par_hat <- c(fit$gamma, fit$beta, fit$log_sigma)
  hessian_mat <- tryCatch(
    stats::optimHess(par = par_hat, fn = neg_loglik),
    error = function(e) NULL
  )
  if (is.null(hessian_mat) || any(!is.finite(hessian_mat))) {
    stop("Failed to compute a finite Hessian for frequentist mixture-cure CI export.", call. = FALSE)
  }

  covariance_mat <- tryCatch(
    invert_hessian(hessian_mat),
    error = function(e) NULL
  )
  if (is.null(covariance_mat) || any(!is.finite(covariance_mat))) {
    stop("Failed to derive a finite covariance approximation for frequentist mixture-cure CI export.", call. = FALSE)
  }

  parameter_draws <- draw_mvn(
    n_draws = n_draws,
    mu = par_hat,
    sigma = covariance_mat,
    seed_value = seed_value
  )

  gamma_draws <- parameter_draws[, seq_len(ncol(X_inc)), drop = FALSE]
  beta_draws <- parameter_draws[, ncol(X_inc) + seq_len(ncol(X_lat)), drop = FALSE]
  log_sigma_draws <- parameter_draws[, length(par_hat)]

  pi_matrix <- plogis(X_inc %*% t(gamma_draws))
  pi_matrix <- pmin(pmax(pi_matrix, 1e-12), 1 - 1e-12)
  latency_lp_matrix <- X_lat %*% t(beta_draws)
  sigma_matrix <- matrix(exp(log_sigma_draws), nrow = nrow(X_lat), ncol = nrow(parameter_draws), byrow = TRUE)

  list(
    X_inc = X_inc,
    X_lat = X_lat,
    pi_matrix = pi_matrix,
    latency_lp_matrix = latency_lp_matrix,
    sigma_matrix = sigma_matrix,
    beta_draws = beta_draws,
    log_sigma_draws = log_sigma_draws,
    n_draws = nrow(parameter_draws)
  )
}

build_curve_interval_summary_df <- function(spec_row, draw_bundle, times) {
  susceptible_fraction_draws <- clamp_prob(colMeans(draw_bundle$pi_matrix, na.rm = TRUE))
  cure_fraction_draws <- clamp_prob(1 - susceptible_fraction_draws)

  bind_rows(lapply(times, function(tt) {
    time_value <- as.numeric(tt)
    if (time_value <= 0) {
      susceptible_survival_draws <- rep(1, draw_bundle$n_draws)
      overall_survival_draws <- rep(1, draw_bundle$n_draws)
    } else {
      z_matrix <- (log(pmax(time_value, time_origin_epsilon_year)) - draw_bundle$latency_lp_matrix) / draw_bundle$sigma_matrix
      susceptible_survival_draws <- clamp_prob(colMeans(1 - stats::pnorm(z_matrix), na.rm = TRUE))
      overall_survival_draws <- clamp_prob(colMeans((1 - draw_bundle$pi_matrix) + draw_bundle$pi_matrix * (1 - stats::pnorm(z_matrix)), na.rm = TRUE))
    }

    susceptible_risk_draws <- clamp_prob(1 - susceptible_survival_draws)
    overall_risk_draws <- clamp_prob(1 - overall_survival_draws)

    tibble(
      dataset = spec_row$dataset[[1]],
      time_year = time_value
    ) %>%
      bind_cols(
        summarize_draws_with_prefix(susceptible_fraction_draws, "susceptible_fraction"),
        summarize_draws_with_prefix(cure_fraction_draws, "cure_fraction"),
        summarize_draws_with_prefix(overall_survival_draws, "overall_survival_prob"),
        summarize_draws_with_prefix(overall_risk_draws, "overall_risk_prob"),
        summarize_draws_with_prefix(susceptible_survival_draws, "susceptible_only_survival_prob"),
        summarize_draws_with_prefix(susceptible_risk_draws, "susceptible_only_risk_prob")
      )
  }))
}

build_latency_summary_interval_df <- function(spec_row, draw_bundle) {
  pi_matrix <- draw_bundle$pi_matrix
  latency_median_matrix <- exp(draw_bundle$latency_lp_matrix)
  latency_mean_matrix <- exp(draw_bundle$latency_lp_matrix + 0.5 * (draw_bundle$sigma_matrix^2))
  sigma_draws <- exp(draw_bundle$log_sigma_draws)

  weighted_susceptible_fraction_draws <- vapply(
    seq_len(draw_bundle$n_draws),
    function(j) weighted_mean_safe(pi_matrix[, j], pi_matrix[, j]),
    numeric(1)
  )
  median_latency_median_draws <- apply(latency_median_matrix, 2, stats::median, na.rm = TRUE)
  median_latency_mean_draws <- apply(latency_mean_matrix, 2, stats::median, na.rm = TRUE)
  weighted_latency_median_draws <- vapply(
    seq_len(draw_bundle$n_draws),
    function(j) weighted_mean_safe(latency_median_matrix[, j], pi_matrix[, j]),
    numeric(1)
  )
  weighted_latency_mean_draws <- vapply(
    seq_len(draw_bundle$n_draws),
    function(j) weighted_mean_safe(latency_mean_matrix[, j], pi_matrix[, j]),
    numeric(1)
  )

  tibble(dataset = spec_row$dataset[[1]]) %>%
    bind_cols(
      summarize_draws_with_prefix(colMeans(pi_matrix, na.rm = TRUE), "mean_susceptible_fraction"),
      summarize_draws_with_prefix(weighted_susceptible_fraction_draws, "weighted_mean_susceptible_fraction"),
      summarize_draws_with_prefix(colMeans(latency_median_matrix, na.rm = TRUE), "mean_subject_specific_latency_median_year"),
      summarize_draws_with_prefix(median_latency_median_draws, "median_subject_specific_latency_median_year"),
      summarize_draws_with_prefix(weighted_latency_median_draws, "susceptible_weighted_mean_subject_specific_latency_median_year"),
      summarize_draws_with_prefix(colMeans(latency_mean_matrix, na.rm = TRUE), "mean_subject_specific_latency_mean_year"),
      summarize_draws_with_prefix(median_latency_mean_draws, "median_subject_specific_latency_mean_year"),
      summarize_draws_with_prefix(weighted_latency_mean_draws, "susceptible_weighted_mean_subject_specific_latency_mean_year"),
      summarize_draws_with_prefix(sigma_draws, "latency_sigma")
    )
}

build_latency_aft_effect_interval_df <- function(spec_row, fit, draw_bundle) {
  latency_coef_names <- fit$latency_coef_names %||% names(fit$beta)
  if (length(latency_coef_names) != ncol(draw_bundle$beta_draws)) {
    latency_coef_names <- names(fit$beta)
  }

  bind_rows(lapply(seq_len(ncol(draw_bundle$beta_draws)), function(j) {
    beta_draws_j <- draw_bundle$beta_draws[, j]
    time_ratio_draws <- exp(beta_draws_j)
    percent_change_draws <- 100 * (time_ratio_draws - 1)

    tibble(
      dataset = spec_row$dataset[[1]],
      latency_term = as.character(latency_coef_names[[j]])
    ) %>%
      bind_cols(
        summarize_draws_with_prefix(beta_draws_j, "latency_log_time"),
        summarize_draws_with_prefix(time_ratio_draws, "latency_time_ratio"),
        summarize_draws_with_prefix(percent_change_draws, "latency_percent_change_in_time"),
        summarize_draws_with_prefix(exp(draw_bundle$log_sigma_draws), "latency_sigma")
      )
  }))
}

validate_inputs <- function(analysis_datasets, formula_registry) {
  if (!is.list(analysis_datasets) || !all(c("PNU", "SNU", "merged") %in% names(analysis_datasets))) {
    stop("Stage 1 analysis dataset bundle must contain `PNU`, `SNU`, and `merged`.", call. = FALSE)
  }

  required_formula_cols <- c(
    "dataset",
    "formula_id",
    "formula_name",
    "formula_label",
    "formula_rhs",
    "site_branch",
    "interaction_branch"
  )
  missing_formula_cols <- setdiff(required_formula_cols, names(formula_registry))
  if (length(missing_formula_cols) > 0L) {
    stop(
      sprintf(
        "Stage 1 formula registry is missing required columns: %s",
        paste(missing_formula_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

prepare_dataset <- function(df, dataset_name) {
  required_cols <- c("site", "id", "sex_num", "age_exact_entry", "age_s", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("[%s] Dataset is missing required columns: %s", dataset_name, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  out <- tibble::as_tibble(df)

  if (!("unique_person_id" %in% names(out))) {
    out$unique_person_id <- paste(out$site, out$id, sep = "_")
  }

  if (!("time_year" %in% names(out))) {
    if (!("days_followup" %in% names(out))) {
      stop(sprintf("[%s] Dataset must contain `time_year` or `days_followup`.", dataset_name), call. = FALSE)
    }
    out$time_year <- as.numeric(out$days_followup) / 365.25
  }

  out <- out %>%
    mutate(
      unique_person_id = as.character(unique_person_id),
      site = factor(as.character(site)),
      id = as.character(id),
      sex_num = as.integer(as.numeric(sex_num)),
      age_exact_entry = as.numeric(age_exact_entry),
      age_s = as.numeric(age_s),
      status_num = as.integer(as.numeric(status_num)),
      time_year = as.numeric(time_year),
      time_year_model = pmax(as.numeric(time_year), time_origin_epsilon_year),
      event_main = as.integer(status_num == 1L),
      censor_main = as.integer(status_num %in% c(0L, 2L))
    )

  if (nrow(out) != dplyr::n_distinct(out$unique_person_id)) {
    stop(sprintf("[%s] `unique_person_id` must be unique.", dataset_name), call. = FALSE)
  }

  if (anyNA(out[, c("unique_person_id", "site", "id", "sex_num", "age_exact_entry", "age_s", "status_num", "time_year")])) {
    stop(sprintf("[%s] Missing values detected in required analysis columns.", dataset_name), call. = FALSE)
  }

  if (any(out$time_year < 0, na.rm = TRUE)) {
    stop(sprintf("[%s] Negative follow-up times are not allowed.", dataset_name), call. = FALSE)
  }

  if (any(!out$status_num %in% c(0L, 1L, 2L))) {
    stop(sprintf("[%s] `status_num` must be coded as 0/1/2 only.", dataset_name), call. = FALSE)
  }

  out
}

select_model_specs <- function(formula_registry) {
  formula_tbl <- formula_registry %>%
    mutate(
      dataset = normalize_dataset_label(dataset),
      formula_name = trimws(as.character(formula_name)),
      formula_rhs = trim_rhs(formula_rhs)
    )

  selected_specs <- dataset_model_registry %>%
    left_join(
      formula_tbl,
      by = c(
        "source_dataset" = "dataset",
        "formula_name" = "formula_name"
      )
    )

  if (anyNA(selected_specs$formula_id) || anyNA(selected_specs$formula_rhs)) {
    missing_keys <- selected_specs$dataset[is.na(selected_specs$formula_id) | is.na(selected_specs$formula_rhs)]
    stop(
      sprintf(
        "Could not resolve formula definitions for dataset versions: %s",
        paste(unique(missing_keys), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Follow the retained Stage 7 simplified rule for the merged site-adjusted
  # mixture cure branch: keep `site` in both the incidence part and the
  # latency part. This matches the earlier `site_in_both` / `base_siteadjusted_both`
  # selection, rather than placing `site` in only one submodel.
  selected_specs %>%
    mutate(
      incidence_rhs = formula_rhs,
      latency_rhs = formula_rhs,
      site_in_incidence = grepl("\\bsite\\b", incidence_rhs),
      site_in_latency = grepl("\\bsite\\b", latency_rhs)
    ) %>%
    select(
      dataset,
      dataset_label,
      source_dataset,
      formula_name,
      formula_id,
      formula_label,
      formula_rhs,
      site_branch,
      interaction_branch,
      site_adjustment_flag,
      incidence_rhs,
      latency_rhs,
      site_in_incidence,
      site_in_latency
    )
}

fit_lognormal_mixture_cure <- function(df, incidence_rhs, latency_rhs, dataset_name) {
  X_inc <- model_matrix_from_rhs(df, incidence_rhs)
  X_lat <- model_matrix_from_rhs(df, latency_rhs)
  y_time <- pmax(as.numeric(df$time_year_model), time_origin_epsilon_year)
  y_event <- as.integer(df$event_main)

  start_gamma <- tryCatch(
    {
      stats::coef(
        stats::glm(
          stats::as.formula(paste("event_main ~", trim_rhs(incidence_rhs))),
          data = df,
          family = stats::binomial(link = "logit")
        )
      )
    },
    error = function(e) rep(0, ncol(X_inc))
  )
  start_gamma <- start_gamma[colnames(X_inc)]
  start_gamma[is.na(start_gamma)] <- 0

  start_survreg <- tryCatch(
    survival::survreg(
      formula = make_surv_formula(latency_rhs),
      data = df,
      dist = "lognormal"
    ),
    error = function(e) NULL
  )

  if (!is.null(start_survreg)) {
    start_beta <- stats::coef(start_survreg)[colnames(X_lat)]
    start_beta[is.na(start_beta)] <- 0
    start_log_sigma <- log(start_survreg$scale)
  } else {
    start_beta <- rep(0, ncol(X_lat))
    names(start_beta) <- colnames(X_lat)
    start_log_sigma <- 0
  }

  par_start <- c(start_gamma, start_beta, start_log_sigma)

  neg_loglik <- function(par) {
    gamma <- par[seq_len(ncol(X_inc))]
    beta <- par[ncol(X_inc) + seq_len(ncol(X_lat))]
    log_sigma <- par[[length(par)]]

    uncured_prob <- clamp_prob(plogis(drop(X_inc %*% gamma)))
    lp_lat <- drop(X_lat %*% beta)
    dens_surv <- lognormal_surv_density(y_time, lp_lat, log_sigma)
    surv_uncured <- pmax(dens_surv$surv, 1e-300)
    dens_uncured <- pmax(dens_surv$dens, 1e-300)

    loglik <- sum(
      y_event * (log(uncured_prob) + log(dens_uncured)) +
        (1 - y_event) * log((1 - uncured_prob) + uncured_prob * surv_uncured)
    )

    if (!is.finite(loglik)) {
      return(1e12)
    }

    -loglik
  }

  attempt_settings <- list(
    list(method = "BFGS", control = list(maxit = 2000, reltol = 1e-10)),
    list(method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-10))
  )

  fit_attempts <- lapply(attempt_settings, function(setting) {
    out <- capture_with_warnings(
      stats::optim(
        par = par_start,
        fn = neg_loglik,
        method = setting$method,
        control = setting$control
      )
    )
    out$method <- setting$method
    out
  })

  valid_attempts <- Filter(
    function(x) !is.null(x$value) && is.list(x$value) && is.finite(x$value$value),
    fit_attempts
  )

  if (length(valid_attempts) == 0L) {
    error_messages <- unique(na.omit(vapply(fit_attempts, function(x) x$error_message, character(1))))
    stop(
      sprintf(
        "[%s] Log-normal mixture cure fit failed: %s",
        dataset_name,
        paste(error_messages, collapse = " | ")
      ),
      call. = FALSE
    )
  }

  converged_attempts <- Filter(
    function(x) identical(as.integer(x$value$convergence), 0L),
    valid_attempts
  )
  candidate_attempts <- if (length(converged_attempts) > 0L) converged_attempts else valid_attempts
  best_index <- which.min(vapply(candidate_attempts, function(x) x$value$value, numeric(1)))
  best_attempt <- candidate_attempts[[best_index]]
  best_fit <- best_attempt$value

  if (!identical(as.integer(best_fit$convergence), 0L)) {
    stop(
      sprintf(
        "[%s] Log-normal mixture cure optimizer did not converge. Best attempt used `%s` with convergence code %s.",
        dataset_name,
        best_attempt$method,
        as.character(best_fit$convergence)
      ),
      call. = FALSE
    )
  }

  gamma_hat <- best_fit$par[seq_len(ncol(X_inc))]
  beta_hat <- best_fit$par[ncol(X_inc) + seq_len(ncol(X_lat))]
  log_sigma_hat <- best_fit$par[[length(best_fit$par)]]

  list(
    gamma = gamma_hat,
    beta = beta_hat,
    log_sigma = log_sigma_hat,
    incidence_rhs = incidence_rhs,
    latency_rhs = latency_rhs,
    incidence_coef_names = colnames(X_inc),
    latency_coef_names = colnames(X_lat),
    site_levels = levels(df$site),
    optimizer_method = best_attempt$method,
    loglik = -best_fit$value,
    convergence_code = as.integer(best_fit$convergence),
    warnings = unique(best_attempt$warnings)
  )
}

load_fitted_object_cache <- function(path) {
  if (!file.exists(path)) {
    return(list())
  }

  cache_obj <- readRDS(path)
  if (!is.list(cache_obj)) {
    return(list())
  }

  cache_obj
}

has_usable_cached_fit <- function(fit_obj, spec_row) {
  required_fields <- c(
    "gamma",
    "beta",
    "log_sigma",
    "incidence_rhs",
    "latency_rhs",
    "site_levels"
  )

  if (!is.list(fit_obj) || !all(required_fields %in% names(fit_obj))) {
    return(FALSE)
  }

  identical(as.character(fit_obj$incidence_rhs), as.character(spec_row$incidence_rhs[[1]])) &&
    identical(as.character(fit_obj$latency_rhs), as.character(spec_row$latency_rhs[[1]]))
}

format_named_estimates <- function(values) {
  if (length(values) == 0L) {
    return(NA_character_)
  }

  paste(
    paste0(
      names(values),
      "=",
      formatC(as.numeric(values), digits = 6L, format = "fg")
    ),
    collapse = "; "
  )
}

build_fit_summary_df <- function(fitted_object_cache, model_specs, analysis_datasets) {
  bind_rows(lapply(seq_len(nrow(model_specs)), function(ii) {
    spec_row <- model_specs[ii, , drop = FALSE]
    dataset_name <- spec_row$dataset[[1]]
    source_dataset_name <- spec_row$source_dataset[[1]]
    analysis_df <- prepare_dataset(
      df = analysis_datasets[[source_dataset_name]],
      dataset_name = dataset_name
    )
    fit_obj <- fitted_object_cache[[dataset_name]]

    if (is.null(fit_obj)) {
      stop(sprintf("Missing fitted object for dataset `%s`.", dataset_name), call. = FALSE)
    }

    incidence_coef_names <- fit_obj$incidence_coef_names %||% names(fit_obj$gamma)
    latency_coef_names <- fit_obj$latency_coef_names %||% names(fit_obj$beta)
    incidence_coef_values <- as.numeric(fit_obj$gamma)
    latency_coef_values <- as.numeric(fit_obj$beta)

    if (length(incidence_coef_names) != length(incidence_coef_values)) {
      incidence_coef_names <- paste0("incidence_coef_", seq_along(incidence_coef_values))
    }
    if (length(latency_coef_names) != length(latency_coef_values)) {
      latency_coef_names <- paste0("latency_coef_", seq_along(latency_coef_values))
    }

    warning_values <- unique(as.character(fit_obj$warnings %||% character()))
    warning_values <- warning_values[nzchar(warning_values)]

    tibble(
      dataset = dataset_name,
      dataset_label = spec_row$dataset_label[[1]],
      source_dataset = source_dataset_name,
      model_id = model_id_value,
      formula_name = spec_row$formula_name[[1]],
      formula_id = spec_row$formula_id[[1]],
      formula_label = spec_row$formula_label[[1]],
      incidence_rhs = spec_row$incidence_rhs[[1]],
      latency_rhs = spec_row$latency_rhs[[1]],
      site_adjustment_flag = as.logical(spec_row$site_adjustment_flag[[1]]),
      site_in_incidence = as.logical(spec_row$site_in_incidence[[1]]),
      site_in_latency = as.logical(spec_row$site_in_latency[[1]]),
      risk_scale = main_risk_scale,
      n_subjects = nrow(analysis_df),
      n_transition = sum(analysis_df$event_main == 1L, na.rm = TRUE),
      n_censored = sum(analysis_df$censor_main == 1L, na.rm = TRUE),
      optimizer_method = as.character(fit_obj$optimizer_method %||% NA_character_),
      convergence_code = as.integer(fit_obj$convergence_code %||% NA_integer_),
      loglik = as.numeric(fit_obj$loglik %||% NA_real_),
      latency_sigma = exp(as.numeric(fit_obj$log_sigma %||% NA_real_)),
      incidence_n_parameters = length(incidence_coef_values),
      latency_n_parameters = length(latency_coef_values),
      incidence_coefficients = format_named_estimates(stats::setNames(incidence_coef_values, incidence_coef_names)),
      latency_coefficients = format_named_estimates(stats::setNames(latency_coef_values, latency_coef_names)),
      warnings_text = if (length(warning_values) > 0L) paste(warning_values, collapse = " | ") else NA_character_,
      has_warnings = length(warning_values) > 0L
    )
  })) %>%
    arrange(match(dataset, dataset_model_registry$dataset))
}

predict_lognormal_mixture_cure <- function(fit, newdata, times) {
  pred_data <- tibble::as_tibble(newdata) %>%
    mutate(site = factor(as.character(site), levels = fit$site_levels))

  X_inc <- model_matrix_from_rhs(pred_data, fit$incidence_rhs)
  X_lat <- model_matrix_from_rhs(pred_data, fit$latency_rhs)

  uncured_prob <- clamp_prob(plogis(drop(X_inc %*% fit$gamma)))
  lp_lat <- drop(X_lat %*% fit$beta)

  susceptible_survival_mat <- sapply(times, function(tt) {
    if (tt <= 0) {
      rep(1, nrow(pred_data))
    } else {
      lognormal_surv_density(rep(tt, nrow(pred_data)), lp_lat, fit$log_sigma)$surv
    }
  })

  if (is.null(dim(susceptible_survival_mat))) {
    susceptible_survival_mat <- matrix(
      susceptible_survival_mat,
      nrow = nrow(pred_data),
      ncol = length(times)
    )
  }

  overall_survival_mat <- (1 - uncured_prob) + uncured_prob * susceptible_survival_mat

  list(
    susceptible_fraction = uncured_prob,
    susceptible_only_survival = susceptible_survival_mat,
    susceptible_only_risk = 1 - susceptible_survival_mat,
    overall_survival = overall_survival_mat,
    overall_risk = 1 - overall_survival_mat
  )
}

build_subject_prediction_df <- function(df, spec_row, fit, times) {
  pred <- predict_lognormal_mixture_cure(fit = fit, newdata = df, times = times)

  bind_rows(lapply(seq_along(times), function(j) {
    tibble(
      unique_person_id = as.character(df$unique_person_id),
      dataset = spec_row$dataset[[1]],
      dataset_label = spec_row$dataset_label[[1]],
      source_dataset = spec_row$source_dataset[[1]],
      model_id = model_id_value,
      site_adjustment_flag = as.logical(spec_row$site_adjustment_flag[[1]]),
      formula_name = spec_row$formula_name[[1]],
      formula_id = spec_row$formula_id[[1]],
      formula_label = spec_row$formula_label[[1]],
      incidence_rhs = spec_row$incidence_rhs[[1]],
      latency_rhs = spec_row$latency_rhs[[1]],
      site_in_incidence = as.logical(spec_row$site_in_incidence[[1]]),
      site_in_latency = as.logical(spec_row$site_in_latency[[1]]),
      risk_scale = main_risk_scale,
      time_year = as.numeric(times[[j]]),
      susceptible_fraction = as.numeric(pred$susceptible_fraction),
      susceptible_only_survival_prob = as.numeric(pred$susceptible_only_survival[, j]),
      susceptible_only_risk_prob = as.numeric(pred$susceptible_only_risk[, j]),
      overall_survival_prob = as.numeric(pred$overall_survival[, j]),
      overall_risk_prob = as.numeric(pred$overall_risk[, j])
    )
  }))
}

build_subject_latency_summary_df <- function(df, spec_row, fit) {
  pred_data <- tibble::as_tibble(df) %>%
    mutate(site = factor(as.character(site), levels = fit$site_levels))

  X_inc <- model_matrix_from_rhs(pred_data, fit$incidence_rhs)
  X_lat <- model_matrix_from_rhs(pred_data, fit$latency_rhs)

  susceptible_fraction <- clamp_prob(plogis(drop(X_inc %*% fit$gamma)))
  lp_lat <- drop(X_lat %*% fit$beta)
  sigma_lat <- exp(fit$log_sigma)

  tibble(
    unique_person_id = as.character(df$unique_person_id),
    dataset = spec_row$dataset[[1]],
    dataset_label = spec_row$dataset_label[[1]],
    source_dataset = spec_row$source_dataset[[1]],
    model_id = model_id_value,
    site_adjustment_flag = as.logical(spec_row$site_adjustment_flag[[1]]),
    formula_name = spec_row$formula_name[[1]],
    formula_id = spec_row$formula_id[[1]],
    formula_label = spec_row$formula_label[[1]],
    incidence_rhs = spec_row$incidence_rhs[[1]],
    latency_rhs = spec_row$latency_rhs[[1]],
    site_in_incidence = as.logical(spec_row$site_in_incidence[[1]]),
    site_in_latency = as.logical(spec_row$site_in_latency[[1]]),
    risk_scale = main_risk_scale,
    susceptible_fraction = as.numeric(susceptible_fraction),
    latency_linear_predictor = as.numeric(lp_lat),
    latency_sigma = sigma_lat,
    latency_median_year = exp(lp_lat),
    latency_mean_year = exp(lp_lat + 0.5 * sigma_lat^2)
  )
}

summarize_uncured_latency <- function(latency_subject_df) {
  latency_subject_df %>%
    group_by(
      dataset,
      dataset_label,
      source_dataset,
      model_id,
      site_adjustment_flag,
      formula_name,
      formula_id,
      formula_label,
      incidence_rhs,
      latency_rhs,
      site_in_incidence,
      site_in_latency,
      risk_scale
    ) %>%
    summarise(
      n_subjects = dplyr::n(),
      mean_susceptible_fraction = mean(susceptible_fraction, na.rm = TRUE),
      weighted_mean_susceptible_fraction = weighted_mean_safe(susceptible_fraction, susceptible_fraction),
      mean_subject_specific_latency_median_year = mean(latency_median_year, na.rm = TRUE),
      median_subject_specific_latency_median_year = stats::median(latency_median_year, na.rm = TRUE),
      susceptible_weighted_mean_subject_specific_latency_median_year = weighted_mean_safe(latency_median_year, susceptible_fraction),
      mean_subject_specific_latency_mean_year = mean(latency_mean_year, na.rm = TRUE),
      median_subject_specific_latency_mean_year = stats::median(latency_mean_year, na.rm = TRUE),
      susceptible_weighted_mean_subject_specific_latency_mean_year = weighted_mean_safe(latency_mean_year, susceptible_fraction),
      latency_sigma = dplyr::first(latency_sigma),
      .groups = "drop"
    ) %>%
    arrange(match(dataset, dataset_model_registry$dataset))
}

build_latency_aft_effect_df <- function(spec_row, fit) {
  latency_coef_names <- fit$latency_coef_names %||% names(fit$beta)
  beta_hat <- as.numeric(fit$beta)
  if (length(latency_coef_names) != length(beta_hat)) {
    latency_coef_names <- names(fit$beta)
  }

  tibble(
    dataset = spec_row$dataset[[1]],
    dataset_label = spec_row$dataset_label[[1]],
    source_dataset = spec_row$source_dataset[[1]],
    model_id = model_id_value,
    site_adjustment_flag = as.logical(spec_row$site_adjustment_flag[[1]]),
    formula_name = spec_row$formula_name[[1]],
    formula_id = spec_row$formula_id[[1]],
    formula_label = spec_row$formula_label[[1]],
    incidence_rhs = spec_row$incidence_rhs[[1]],
    latency_rhs = spec_row$latency_rhs[[1]],
    site_in_incidence = as.logical(spec_row$site_in_incidence[[1]]),
    site_in_latency = as.logical(spec_row$site_in_latency[[1]]),
    risk_scale = main_risk_scale,
    latency_term = as.character(latency_coef_names),
    latency_log_time_estimate = beta_hat,
    latency_time_ratio = exp(beta_hat),
    latency_percent_change_in_time = 100 * (exp(beta_hat) - 1),
    latency_sigma = exp(fit$log_sigma),
    effect_model = "lognormal_aft_conditional_on_uncured"
  )
}

make_latency_subject_distribution_plot <- function(latency_subject_df, value_col, y_label, title_text) {
  plot_df <- latency_subject_df %>%
    mutate(
      plot_dataset_label = dplyr::coalesce(unname(plot_dataset_label_lookup[dataset]), dataset_label)
    )

  ggplot(plot_df, aes(x = plot_dataset_label, y = .data[[value_col]], color = dataset)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.6, width = 0.55) +
    geom_jitter(width = 0.12, alpha = 0.35, size = 1.1) +
    scale_color_manual(
      values = dataset_palette,
      breaks = dataset_model_registry$dataset,
      labels = unname(plot_dataset_label_lookup[dataset_model_registry$dataset])
    ) +
    scale_y_log10(labels = scales::label_number(accuracy = 0.1)) +
    labs(
      title = title_text,
      subtitle = "Subject-level uncured-latency summaries; y-axis shown on log10 scale",
      x = NULL,
      y = y_label,
      color = "Dataset"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 12, hjust = 1)
    )
}

make_latency_summary_plot <- function(latency_summary_df) {
  summary_long_df <- bind_rows(
    latency_summary_df %>%
      transmute(
        dataset,
        plot_dataset_label = dplyr::coalesce(unname(plot_dataset_label_lookup[dataset]), dataset_label),
        metric = "Weighted mean median latency",
        value_year = susceptible_weighted_mean_subject_specific_latency_median_year
      ),
    latency_summary_df %>%
      transmute(
        dataset,
        plot_dataset_label = dplyr::coalesce(unname(plot_dataset_label_lookup[dataset]), dataset_label),
        metric = "Weighted mean mean latency",
        value_year = susceptible_weighted_mean_subject_specific_latency_mean_year
      )
  )

  ggplot(summary_long_df, aes(x = value_year, y = plot_dataset_label, color = dataset)) +
    geom_point(size = 3) +
    facet_wrap(~ metric, scales = "free_x") +
    scale_color_manual(
      values = dataset_palette,
      breaks = dataset_model_registry$dataset,
      labels = unname(plot_dataset_label_lookup[dataset_model_registry$dataset])
    ) +
    scale_x_log10(labels = scales::label_number(accuracy = 0.1)) +
    labs(
      title = "Uncured-latency summary comparison across datasets",
      subtitle = "Susceptible-fraction weighted latency summaries; x-axis shown on log10 scale",
      x = "Years",
      y = NULL,
      color = "Dataset"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
}

make_latency_aft_effect_plot <- function(latency_aft_effect_df) {
  plot_df <- latency_aft_effect_df %>%
    mutate(
      plot_dataset_label = dplyr::coalesce(unname(plot_dataset_label_lookup[dataset]), dataset_label),
      latency_term_label = dplyr::case_when(
        latency_term == "(Intercept)" ~ "Intercept",
        latency_term == "age_s" ~ "Age (scaled)",
        latency_term == "sex_num" ~ "Male vs Female",
        grepl("^site", latency_term) ~ "Site effect",
        TRUE ~ latency_term
      )
    )
  term_levels <- rev(unique(plot_df$latency_term_label))
  plot_df <- plot_df %>%
    mutate(latency_term_label = factor(latency_term_label, levels = term_levels))

  ggplot(plot_df, aes(x = latency_time_ratio, y = latency_term_label, color = dataset)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_point(size = 2.8, position = position_dodge(width = 0.5)) +
    facet_wrap(~ plot_dataset_label) +
    scale_color_manual(
      values = dataset_palette,
      breaks = dataset_model_registry$dataset,
      labels = unname(plot_dataset_label_lookup[dataset_model_registry$dataset])
    ) +
    scale_x_log10(labels = scales::label_number(accuracy = 0.1)) +
    labs(
      title = "Latency AFT effects among the uncured subgroup",
      subtitle = "Time ratio > 1 indicates longer time to transition; < 1 indicates shorter time",
      x = "Time ratio (log10 scale)",
      y = NULL,
      color = "Dataset"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
}

build_plot_group_registry <- function() {
  tibble::tibble(
    plot_group = c("all_cohorts", "pnu_only", "snu_only", "merged_only"),
    group_title = c(
      "PNU, SNU, and merged cohorts",
      "PNU cohort",
      "SNU cohort",
      "Merged cohort models"
    ),
    dataset_keys = list(
      c("PNU", "SNU", "merged_no_site", "merged_site_adjusted"),
      "PNU",
      "SNU",
      c("merged_no_site", "merged_site_adjusted")
    ),
    include_pnu_reference = c(TRUE, TRUE, FALSE, FALSE)
  )
}

build_dataset_detail_plot_registry <- function() {
  tibble::tibble(
    plot_group = c("pnu_detail", "snu_detail", "merged_detail"),
    group_title = c(
      "PNU cohort",
      "SNU cohort",
      "Merged cohort models"
    ),
    dataset_keys = list(
      "PNU",
      "SNU",
      c("merged_no_site", "merged_site_adjusted")
    ),
    source_dataset = c("PNU", "SNU", "merged"),
    include_pnu_reference = c(TRUE, FALSE, FALSE)
  )
}

make_plot_output_file <- function(export_dir, metric_name, plot_group) {
  file.path(
    export_dir,
    paste0(output_prefix, "_", metric_name, "__", plot_group, ".png")
  )
}

format_half_year_tick_labels <- function(x) {
  ifelse(
    abs(x - round(x)) < 1e-9,
    formatC(x, format = "f", digits = 0),
    formatC(x, format = "f", digits = 1)
  )
}

derive_interval_colnames <- function(value_col) {
  c(
    lower = paste0(value_col, "_q025"),
    upper = paste0(value_col, "_q975")
  )
}

build_detail_annotation_table <- function(analysis_df, times) {
  tibble(
    time_horizon_year = as.numeric(times),
    n_at_risk = vapply(
      times,
      function(tt) sum(as.numeric(analysis_df$time_year) >= tt, na.rm = TRUE),
      integer(1)
    ),
    n_transition_cumulative = vapply(
      times,
      function(tt) sum(as.integer(analysis_df$event_main) == 1L & as.numeric(analysis_df$time_year) <= tt, na.rm = TRUE),
      integer(1)
    )
  )
}

remove_existing_detail_plot_outputs <- function(export_dir, annotation_file) {
  if (!dir.exists(export_dir)) {
    return(invisible(NULL))
  }

  existing_detail_pngs <- list.files(
    export_dir,
    pattern = paste0("^", output_prefix, "_.*with_counts__.*\\.png$"),
    full.names = TRUE
  )
  existing_targets <- unique(c(existing_detail_pngs, annotation_file[file.exists(annotation_file)]))

  if (length(existing_targets) > 0L) {
    unlink(existing_targets)
  }

  invisible(existing_targets)
}

make_curve_plot <- function(
  curve_df,
  value_col,
  y_label,
  title_text,
  pnu_reference_year = NA_real_,
  x_breaks = horizon_years,
  x_labels = scales::label_number(accuracy = 1),
  y_lower_limit = 0,
  caption_text = NULL
) {
  dataset_labels <- curve_df %>%
    distinct(dataset, plot_dataset_label) %>%
    arrange(match(dataset, dataset_model_registry$dataset))
  interval_cols <- derive_interval_colnames(value_col)
  has_interval_cols <- all(interval_cols %in% names(curve_df))
  yearly_points <- curve_df %>%
    filter(time_horizon_year %in% horizon_years)
  show_legend <- nrow(dataset_labels) > 1L
  subtitle_text <- "Cohort-level means of subject-level predictions from frequentist log-normal mixture cure models"

  if (is.finite(pnu_reference_year) && any(curve_df$dataset == "PNU")) {
    subtitle_text <- paste0(
      subtitle_text,
      sprintf(" | Dashed line = PNU max observed follow-up (%.2f years)", pnu_reference_year)
    )
  }

  plot_object <- ggplot(curve_df, aes(x = time_horizon_year, y = .data[[value_col]], color = dataset))

  if (has_interval_cols) {
    plot_object <- plot_object +
      geom_ribbon(
        aes(
          x = time_horizon_year,
          ymin = .data[[interval_cols[["lower"]]]],
          ymax = .data[[interval_cols[["upper"]]]],
          fill = dataset,
          group = dataset
        ),
        inherit.aes = FALSE,
        alpha = 0.14,
        color = NA
      ) +
      geom_linerange(
        data = yearly_points,
        aes(
          x = time_horizon_year,
          ymin = .data[[interval_cols[["lower"]]]],
          ymax = .data[[interval_cols[["upper"]]]],
          color = dataset
        ),
        inherit.aes = FALSE,
        linewidth = 0.45,
        alpha = 0.80
      )
  }

  plot_object <- plot_object +
    geom_line(linewidth = 1.05) +
    geom_point(
      data = yearly_points,
      aes(x = time_horizon_year, y = .data[[value_col]], color = dataset),
      size = 1.8
    ) +
    scale_fill_manual(
      values = dataset_palette,
      breaks = dataset_labels$dataset,
      labels = dataset_labels$plot_dataset_label,
      guide = "none"
    ) +
    scale_color_manual(
      values = dataset_palette,
      breaks = dataset_labels$dataset,
      labels = dataset_labels$plot_dataset_label
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      limits = c(0, curve_horizon_max_year),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      breaks = seq(0, 1, by = 0.2),
      limits = c(y_lower_limit, 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Years after cohort entry (k)",
      y = y_label,
      color = "Dataset",
      caption = caption_text
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = if (show_legend) "top" else "none",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 9, hjust = 0),
      plot.caption.position = "plot",
      plot.margin = margin(
        t = 5.5,
        r = 5.5,
        b = if (y_lower_limit < 0) 28 else 5.5,
        l = 5.5
      ),
      axis.text.x = element_text(size = if (length(x_breaks) > length(horizon_years)) 7 else 9)
    ) +
    coord_cartesian(clip = "off")

  if (is.finite(pnu_reference_year) && any(curve_df$dataset == "PNU")) {
    plot_object <- plot_object +
      geom_vline(
        xintercept = pnu_reference_year,
        linetype = "dashed",
        linewidth = 0.7,
        color = dataset_palette[["PNU"]]
      )
  }

  plot_object
}

make_detail_curve_plot <- function(
  curve_df,
  value_col,
  y_label,
  title_text,
  analysis_df,
  pnu_reference_year = NA_real_
) {
  detail_times <- seq(0, curve_horizon_max_year, by = detail_tick_step_year)
  detail_bar_x_limits <- c(
    -detail_tick_step_year * 0.4,
    curve_horizon_max_year + detail_tick_step_year * 0.4
  )
  curve_plot <- make_curve_plot(
    curve_df = curve_df,
    value_col = value_col,
    y_label = y_label,
    title_text = title_text,
    pnu_reference_year = pnu_reference_year,
    x_breaks = detail_times,
    x_labels = format_half_year_tick_labels,
    y_lower_limit = 0,
    caption_text = NULL
  ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5)
    )

  annotation_tbl <- build_detail_annotation_table(analysis_df, detail_times)

  at_risk_plot <- ggplot(annotation_tbl, aes(x = time_horizon_year, y = n_at_risk)) +
    geom_col(width = detail_tick_step_year * 0.72, fill = "#8FA3BF") +
    geom_text(aes(label = n_at_risk), vjust = -0.25, size = 2.8) +
    scale_x_continuous(
      breaks = detail_times,
      labels = format_half_year_tick_labels,
      limits = detail_bar_x_limits,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(x = NULL, y = "n at risk") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 5.5, b = 0, l = 5.5)
    )

  transition_plot <- ggplot(annotation_tbl, aes(x = time_horizon_year, y = n_transition_cumulative)) +
    geom_col(width = detail_tick_step_year * 0.72, fill = "#D98C6C") +
    geom_text(aes(label = n_transition_cumulative), vjust = -0.25, size = 2.8) +
    scale_x_continuous(
      breaks = detail_times,
      labels = format_half_year_tick_labels,
      limits = detail_bar_x_limits,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(
      x = "Years after cohort entry (k)",
      y = "cumulative\ntransitions"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(size = 7),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5)
    )

  list(
    curve_plot = curve_plot,
    at_risk_plot = at_risk_plot,
    transition_plot = transition_plot,
    annotation_tbl = annotation_tbl
  )
}

save_plot_png <- function(plot_object, output_file, width = plot_width_in, height = plot_height_in) {
  ggplot2::ggsave(
    filename = output_file,
    plot = plot_object,
    width = width,
    height = height,
    dpi = plot_dpi,
    units = "in"
  )
}

save_stacked_plot_png <- function(plot_bundle, output_file, width = detail_plot_width_in, height = detail_plot_height_in) {
  grDevices::png(filename = output_file, width = width, height = height, units = "in", res = plot_dpi)
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()
  grid::pushViewport(
    grid::viewport(
      layout = grid::grid.layout(
        nrow = 3,
        ncol = 1,
        heights = grid::unit(c(3.8, 1.4, 1.4), "null")
      )
    )
  )

  print(plot_bundle$curve_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(plot_bundle$at_risk_plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(plot_bundle$transition_plot, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
}

# Build Direct Analysis Inputs --------------------------------------------
analysis_datasets <- build_analysis_datasets_from_source(
  source_file = source_data_file,
  pnu_label = pnu_site_label,
  snu_label = snu_site_label
)
formula_registry <- make_formula_registry()
observed_followup_summary <- summarize_observed_followup(analysis_datasets)
pnu_max_observed_followup_year <- observed_followup_summary %>%
  filter(source_dataset == "PNU") %>%
  pull(max_observed_followup_year) %>%
  .[[1L]]

validate_inputs(analysis_datasets, formula_registry)
model_specs <- select_model_specs(formula_registry)

# Fit Models and Build Subject-Level Predictions --------------------------
curve_time_grid <- sort(unique(c(0, seq(0, curve_horizon_max_year, by = curve_step_year), horizon_years)))

fitted_object_cache <- if (reuse_existing_fit_rds) load_fitted_object_cache(fit_object_rds_file) else list()
subject_prediction_list <- vector("list", nrow(model_specs))
latency_subject_list <- vector("list", nrow(model_specs))
latency_aft_effect_list <- vector("list", nrow(model_specs))
curve_interval_summary_list <- vector("list", nrow(model_specs))
latency_summary_interval_list <- vector("list", nrow(model_specs))
latency_aft_interval_list <- vector("list", nrow(model_specs))

for (ii in seq_len(nrow(model_specs))) {
  spec_row <- model_specs[ii, , drop = FALSE]
  dataset_name <- spec_row$dataset[[1]]
  source_dataset_name <- spec_row$source_dataset[[1]]

  analysis_df <- prepare_dataset(
    df = analysis_datasets[[source_dataset_name]],
    dataset_name = dataset_name
  )

  cached_fit_obj <- fitted_object_cache[[dataset_name]]
  if (reuse_existing_fit_rds && has_usable_cached_fit(cached_fit_obj, spec_row)) {
    fit_obj <- cached_fit_obj
    message(sprintf("[%s] Using cached fitted object from %s", dataset_name, basename(fit_object_rds_file)))
  } else {
    message(sprintf("[%s] Fitting log-normal mixture cure model and updating %s", dataset_name, basename(fit_object_rds_file)))
    fit_obj <- fit_lognormal_mixture_cure(
      df = analysis_df,
      incidence_rhs = spec_row$incidence_rhs[[1]],
      latency_rhs = spec_row$latency_rhs[[1]],
      dataset_name = dataset_name
    )
    fitted_object_cache[[dataset_name]] <- fit_obj
    saveRDS(fitted_object_cache, fit_object_rds_file)
  }

  subject_prediction_list[[ii]] <- build_subject_prediction_df(
    df = analysis_df,
    spec_row = spec_row,
    fit = fit_obj,
    times = curve_time_grid
  )

  latency_subject_list[[ii]] <- build_subject_latency_summary_df(
    df = analysis_df,
    spec_row = spec_row,
    fit = fit_obj
  )

  latency_aft_effect_list[[ii]] <- build_latency_aft_effect_df(
    spec_row = spec_row,
    fit = fit_obj
  )

  draw_bundle <- build_prediction_draw_bundle(
    df = analysis_df,
    fit = fit_obj,
    n_draws = stage7_ci_draws,
    seed_value = 930000L + ii
  )
  curve_interval_summary_list[[ii]] <- build_curve_interval_summary_df(
    spec_row = spec_row,
    draw_bundle = draw_bundle,
    times = curve_time_grid
  )
  latency_summary_interval_list[[ii]] <- build_latency_summary_interval_df(
    spec_row = spec_row,
    draw_bundle = draw_bundle
  )
  latency_aft_interval_list[[ii]] <- build_latency_aft_effect_interval_df(
    spec_row = spec_row,
    fit = fit_obj,
    draw_bundle = draw_bundle
  )
}

subject_prediction_df <- bind_rows(subject_prediction_list)
latency_subject_df <- bind_rows(latency_subject_list)
latency_summary_df <- summarize_uncured_latency(latency_subject_df)
latency_aft_effect_df <- bind_rows(latency_aft_effect_list)
curve_interval_summary_df <- bind_rows(curve_interval_summary_list)
latency_summary_interval_df <- bind_rows(latency_summary_interval_list)
latency_aft_interval_df <- bind_rows(latency_aft_interval_list)
fit_summary_df <- build_fit_summary_df(
  fitted_object_cache = fitted_object_cache,
  model_specs = model_specs,
  analysis_datasets = analysis_datasets
)

latency_summary_df <- latency_summary_df %>%
  left_join(latency_summary_interval_df, by = "dataset") %>%
  mutate(interval_type = "95% CI")

latency_aft_effect_df <- latency_aft_effect_df %>%
  left_join(latency_aft_interval_df, by = c("dataset", "latency_term")) %>%
  mutate(interval_type = "95% CI")

# Aggregate Cohort-Level Predictions --------------------------------------
cohort_curve_df <- subject_prediction_df %>%
  group_by(
    dataset,
    dataset_label,
    source_dataset,
    model_id,
    site_adjustment_flag,
    formula_name,
    formula_id,
    formula_label,
    incidence_rhs,
    latency_rhs,
    site_in_incidence,
    site_in_latency,
    risk_scale,
    time_year
  ) %>%
  summarise(
    susceptible_fraction = mean(susceptible_fraction, na.rm = TRUE),
    overall_survival_prob = mean(overall_survival_prob, na.rm = TRUE),
    overall_risk_prob = mean(overall_risk_prob, na.rm = TRUE),
    susceptible_only_survival_prob = mean(susceptible_only_survival_prob, na.rm = TRUE),
    susceptible_only_risk_prob = mean(susceptible_only_risk_prob, na.rm = TRUE),
    n_subjects = dplyr::n(),
    .groups = "drop"
  ) %>%
  left_join(
    curve_interval_summary_df,
    by = c("dataset", "time_year")
  ) %>%
  arrange(match(dataset, dataset_model_registry$dataset), time_year)

plot_source_df <- cohort_curve_df %>%
  transmute(
    dataset,
    plot_dataset_label = dplyr::coalesce(
      unname(plot_dataset_label_lookup[dataset]),
      dataset_label
    ),
    model_formula = if_else(
      incidence_rhs == latency_rhs,
      paste0("incidence = latency = ", incidence_rhs),
      paste0("incidence = ", incidence_rhs, " ; latency = ", latency_rhs)
    ),
    interval_type = "95% CI",
    time_horizon_year = time_year,
    susceptible_fraction,
    susceptible_fraction_sd,
    susceptible_fraction_q025,
    susceptible_fraction_q50,
    susceptible_fraction_q975,
    cure_fraction = 1 - susceptible_fraction,
    cure_fraction_sd,
    cure_fraction_q025,
    cure_fraction_q50,
    cure_fraction_q975,
    overall_survival_prob,
    overall_survival_prob_sd,
    overall_survival_prob_q025,
    overall_survival_prob_q50,
    overall_survival_prob_q975,
    overall_risk_prob,
    overall_risk_prob_sd,
    overall_risk_prob_q025,
    overall_risk_prob_q50,
    overall_risk_prob_q975,
    susceptible_only_survival_prob,
    susceptible_only_survival_prob_sd,
    susceptible_only_survival_prob_q025,
    susceptible_only_survival_prob_q50,
    susceptible_only_survival_prob_q975,
    susceptible_only_risk_prob,
    susceptible_only_risk_prob_sd,
    susceptible_only_risk_prob_q025,
    susceptible_only_risk_prob_q50,
    susceptible_only_risk_prob_q975
  ) %>%
  arrange(match(dataset, dataset_model_registry$dataset), time_horizon_year)

horizon_summary_df <- plot_source_df %>%
  filter(time_horizon_year %in% horizon_years)

# Build Plots -------------------------------------------------------------
plot_group_registry <- build_plot_group_registry()
detail_plot_group_registry <- build_dataset_detail_plot_registry()
plot_registry <- list()
detail_annotation_list <- list()

for (ii in seq_len(nrow(plot_group_registry))) {
  plot_group <- plot_group_registry$plot_group[[ii]]
  group_title <- plot_group_registry$group_title[[ii]]
  group_dataset_keys <- plot_group_registry$dataset_keys[[ii]]
  include_pnu_reference <- isTRUE(plot_group_registry$include_pnu_reference[[ii]])
  group_curve_df <- plot_source_df %>%
    filter(dataset %in% group_dataset_keys)

  if (nrow(group_curve_df) == 0L) {
    stop(sprintf("No plotting rows found for plot group `%s`.", plot_group), call. = FALSE)
  }

  pnu_reference_year <- if (include_pnu_reference) pnu_max_observed_followup_year else NA_real_

  survival_plot <- make_curve_plot(
    curve_df = group_curve_df,
    value_col = "overall_survival_prob",
    y_label = "Estimated survival probability",
    title_text = paste0("Estimated survival probability: ", group_title),
    pnu_reference_year = pnu_reference_year
  )

  risk_plot <- make_curve_plot(
    curve_df = group_curve_df,
    value_col = "overall_risk_prob",
    y_label = "Estimated risk probability (1 - survival)",
    title_text = paste0("Estimated risk probability: ", group_title),
    pnu_reference_year = pnu_reference_year
  )

  plot_registry[[paste0(plot_group, "__survival_curve")]] <- survival_plot
  plot_registry[[paste0(plot_group, "__risk_curve")]] <- risk_plot

  save_plot_png(
    survival_plot,
    make_plot_output_file(export_path, "estimated_survival_curve", plot_group)
  )
  save_plot_png(
    risk_plot,
    make_plot_output_file(export_path, "estimated_risk_curve", plot_group)
  )
}

remove_existing_detail_plot_outputs(
  export_dir = export_path,
  annotation_file = detail_annotation_file
)

for (ii in seq_len(nrow(detail_plot_group_registry))) {
  plot_group <- detail_plot_group_registry$plot_group[[ii]]
  group_title <- detail_plot_group_registry$group_title[[ii]]
  group_dataset_keys <- detail_plot_group_registry$dataset_keys[[ii]]
  source_dataset_name <- detail_plot_group_registry$source_dataset[[ii]]
  include_pnu_reference <- isTRUE(detail_plot_group_registry$include_pnu_reference[[ii]])
  analysis_df <- analysis_datasets[[source_dataset_name]]
  group_curve_df <- plot_source_df %>%
    filter(dataset %in% group_dataset_keys)

  if (nrow(group_curve_df) == 0L) {
    stop(sprintf("No plotting rows found for detail plot group `%s`.", plot_group), call. = FALSE)
  }

  pnu_reference_year <- if (include_pnu_reference) pnu_max_observed_followup_year else NA_real_

  detail_annotation_list[[ii]] <- build_detail_annotation_table(
    analysis_df = analysis_df,
    times = seq(0, curve_horizon_max_year, by = detail_tick_step_year)
  ) %>%
    mutate(
      plot_group = plot_group,
      source_dataset = source_dataset_name
    )

  survival_detail_plot_bundle <- make_detail_curve_plot(
    curve_df = group_curve_df,
    value_col = "overall_survival_prob",
    y_label = "Estimated survival probability",
    title_text = paste0("Estimated survival probability with risk counts: ", group_title),
    analysis_df = analysis_df,
    pnu_reference_year = pnu_reference_year
  )

  risk_detail_plot_bundle <- make_detail_curve_plot(
    curve_df = group_curve_df,
    value_col = "overall_risk_prob",
    y_label = "Estimated risk probability (1 - survival)",
    title_text = paste0("Estimated risk probability with risk counts: ", group_title),
    analysis_df = analysis_df,
    pnu_reference_year = pnu_reference_year
  )

  plot_registry[[paste0(plot_group, "__survival_curve_with_counts")]] <- survival_detail_plot_bundle
  plot_registry[[paste0(plot_group, "__risk_curve_with_counts")]] <- risk_detail_plot_bundle

  save_stacked_plot_png(
    survival_detail_plot_bundle,
    make_plot_output_file(export_path, "estimated_survival_curve_with_counts", plot_group),
    width = detail_plot_width_in,
    height = detail_plot_height_in
  )
  save_stacked_plot_png(
    risk_detail_plot_bundle,
    make_plot_output_file(export_path, "estimated_risk_curve_with_counts", plot_group),
    width = detail_plot_width_in,
    height = detail_plot_height_in
  )
}

detail_annotation_df <- bind_rows(detail_annotation_list)

latency_plot_source_df <- bind_rows(
  latency_subject_df %>%
    mutate(component_type = "subject_distribution"),
  latency_summary_df %>%
    mutate(component_type = "dataset_summary"),
  latency_aft_effect_df %>%
    mutate(component_type = "aft_effect")
) %>%
  mutate(
    component_type = factor(
      component_type,
      levels = c("subject_distribution", "dataset_summary", "aft_effect")
    )
  ) %>%
  arrange(match(dataset, dataset_model_registry$dataset), component_type) %>%
  mutate(component_type = as.character(component_type))

latency_median_distribution_plot <- make_latency_subject_distribution_plot(
  latency_subject_df = latency_subject_df,
  value_col = "latency_median_year",
  y_label = "Subject-specific median latency (years)",
  title_text = "Uncured-latency distribution: median time to transition"
)

latency_mean_distribution_plot <- make_latency_subject_distribution_plot(
  latency_subject_df = latency_subject_df,
  value_col = "latency_mean_year",
  y_label = "Subject-specific mean latency (years)",
  title_text = "Uncured-latency distribution: mean time to transition"
)

latency_summary_plot <- make_latency_summary_plot(latency_summary_df)
latency_aft_effect_plot <- make_latency_aft_effect_plot(latency_aft_effect_df)

plot_registry[["uncured_latency_median_distribution"]] <- latency_median_distribution_plot
plot_registry[["uncured_latency_mean_distribution"]] <- latency_mean_distribution_plot
plot_registry[["uncured_latency_summary_comparison"]] <- latency_summary_plot
plot_registry[["latency_aft_time_ratio_plot"]] <- latency_aft_effect_plot

standard_error_export_bundle <- build_standard_error_export_bundle(list(
  horizon_summary = horizon_summary_df,
  fit_summary = fit_summary_df,
  plot_source = plot_source_df,
  detail_annotation = detail_annotation_df,
  latency_plot_source = latency_plot_source_df
))

# Write Outputs -----------------------------------------------------------
cleanup_existing_stage7_outputs(export_path)
readr::write_csv(horizon_summary_df, horizon_summary_file)
readr::write_csv(fit_summary_df, fit_summary_file)
readr::write_csv(plot_source_df, plot_source_file)
readr::write_csv(detail_annotation_df, detail_annotation_file)
readr::write_csv(latency_plot_source_df, latency_plot_source_file)
readr::write_csv(standard_error_export_bundle$registry, standard_error_registry_file)
readr::write_csv(standard_error_export_bundle$long, standard_error_long_file)
saveRDS(fitted_object_cache, fit_object_rds_file)
saveRDS(plot_registry, plot_rds_file)
save_plot_png(latency_median_distribution_plot, latency_median_distribution_png_file)
save_plot_png(latency_mean_distribution_plot, latency_mean_distribution_png_file)
save_plot_png(latency_summary_plot, latency_summary_png_file, width = 12, height = 6)
save_plot_png(latency_aft_effect_plot, latency_aft_effect_png_file, width = 12, height = 6)

legacy_curve_data_file <- file.path(
  export_path,
  paste0(output_prefix, "_curve_data.csv")
)
if (file.exists(legacy_curve_data_file)) {
  unlink(legacy_curve_data_file)
}

organize_stage7_outputs(export_path)
write_mle_mixture_cure_readme(export_path)

message(analysis_label, " export completed.")
message("Horizon summary CSV: ", normalizePath(horizon_summary_file, winslash = "/", mustWork = FALSE))
message("Fit summary CSV: ", normalizePath(fit_summary_file, winslash = "/", mustWork = FALSE))
message("Plot source CSV: ", normalizePath(plot_source_file, winslash = "/", mustWork = FALSE))
message("Detail annotation CSV: ", normalizePath(detail_annotation_file, winslash = "/", mustWork = FALSE))
message("Latency plot source CSV: ", normalizePath(latency_plot_source_file, winslash = "/", mustWork = FALSE))
message("Standard-error registry CSV: ", normalizePath(standard_error_registry_file, winslash = "/", mustWork = FALSE))
message("Standard-error long CSV: ", normalizePath(standard_error_long_file, winslash = "/", mustWork = FALSE))
message("Fitted object RDS: ", normalizePath(fit_object_rds_file, winslash = "/", mustWork = FALSE))
message("Plot RDS: ", normalizePath(plot_rds_file, winslash = "/", mustWork = FALSE))
