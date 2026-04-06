# Configure paths ---------------------------------------------------------------
source_file <- Sys.getenv(
  "CHR_P_DATA_SUMMARY_SOURCE_FILE",
  unset = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/0.Data/2.Preprocessed data/Preprocessed_Merged_PNUH_SNUH_Data.csv"
)

dictionary_file <- Sys.getenv(
  "CHR_P_DATA_SUMMARY_DICTIONARY_FILE",
  unset = "/Users/ido/Documents/GitHub/bayesian-cure-model-chrp-transition/0.Data/3.Dictionary/Merged Data Dictionary_\U0001f1ec\U0001f1e7ENG.md"
)

export_path <- Sys.getenv(
  "CHR_P_DATA_SUMMARY_EXPORT_PATH",
  unset = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/0.Data/3.Data Summary"
)

pnu_site_label <- "PNU"
snu_site_label <- "SNU"
days_per_year <- 365.25

dataset_order <- c("PNU", "SNU", "Merged")
scenario_order <- c("observed_status", "remission_censored")
section_order <- c(
  "cohort_overview",
  "core_missingness",
  "data_flags",
  "sex_distribution",
  "original_status_distribution",
  "analysis_status_distribution",
  "age_summary",
  "followup_summary"
)

sex_report_levels <- c("Female", "Male", "Missing", "Unexpected_code")
original_status_levels <- c("right_censoring", "transition", "remission", "Missing", "Unexpected_code")
analysis_status_levels <- list(
  observed_status = c("right_censoring", "transition", "remission", "Missing", "Unexpected_code"),
  remission_censored = c("censoring_including_remission", "transition", "Missing", "Unexpected_code")
)

scenario_labels <- c(
  observed_status = "Observed status: remission retained as a separate outcome",
  remission_censored = "Transition-only status: remission treated as censoring"
)

followup_dataset_palette <- c(
  "PNU" = "#1B4332",
  "SNU" = "#2A6F97",
  "Merged" = "#C1666B"
)

followup_panel_order <- c(
  "Observed status\nRight censoring",
  "Observed status\nTransition",
  "Observed status\nRemission",
  "Remission as censoring\nCensoring + remission",
  "Remission as censoring\nTransition"
)

followup_panel_title_lookup <- c(
  "Observed status\nRight censoring" = "Observed status: right censoring",
  "Observed status\nTransition" = "Observed status: transition",
  "Observed status\nRemission" = "Observed status: remission",
  "Remission as censoring\nCensoring + remission" = "Remission treated as censoring: censoring + remission",
  "Remission as censoring\nTransition" = "Remission treated as censoring: transition"
)

followup_panel_file_stub_lookup <- c(
  "Observed status\nRight censoring" = "observed_status__right_censoring",
  "Observed status\nTransition" = "observed_status__transition",
  "Observed status\nRemission" = "observed_status__remission",
  "Remission as censoring\nCensoring + remission" = "remission_as_censoring__censoring_plus_remission",
  "Remission as censoring\nTransition" = "remission_as_censoring__transition"
)

master_output_file <- "DESC__data_summary__master.csv"
compact_output_file <- "DESC__data_summary__compact.csv"
compact_observed_output_file <- "DESC__data_summary__compact__observed_status.csv"
compact_remission_censored_output_file <- "DESC__data_summary__compact__remission_censored.csv"
workbook_output_file <- "DESC__data_summary__tables.xlsx"
followup_boxplot_png_file <- "PLOT__data_summary__followup_time_boxplot.png"
followup_boxplot_pdf_file <- "PLOT__data_summary__followup_time_boxplot.pdf"
followup_boxplot_panel_dir_name <- "PLOT__data_summary__followup_time_boxplot_panels"
followup_jitter_png_file <- "PLOT__data_summary__followup_time_jitter.png"
followup_jitter_pdf_file <- "PLOT__data_summary__followup_time_jitter.pdf"

# Initialize packages and runtime options --------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# Define helpers ----------------------------------------------------------------
read_input_dataset <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Input file does not exist: %s", path), call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))
  if (ext != "csv") {
    stop(sprintf("Only CSV input is supported for this script: %s", path), call. = FALSE)
  }

  readr::read_csv(
    file = path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
}

coerce_numeric_text <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

count_missing_like <- function(x) {
  if (is.character(x)) {
    return(sum(is.na(x) | trimws(x) == ""))
  }
  sum(is.na(x))
}

safe_mean <- function(x) {
  if (sum(!is.na(x)) == 0L) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (sum(!is.na(x)) <= 1L) return(NA_real_)
  stats::sd(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (sum(!is.na(x)) == 0L) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

safe_quantile <- function(x, prob) {
  if (sum(!is.na(x)) == 0L) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 7))
}

safe_min <- function(x) {
  if (sum(!is.na(x)) == 0L) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (sum(!is.na(x)) == 0L) return(NA_real_)
  max(x, na.rm = TRUE)
}

safe_sum <- function(x) {
  if (sum(!is.na(x)) == 0L) return(NA_real_)
  sum(x, na.rm = TRUE)
}

classify_sex <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "Missing",
    x == 0L ~ "Female",
    x == 1L ~ "Male",
    TRUE ~ "Unexpected_code"
  )
}

classify_original_status <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "Missing",
    x == 0L ~ "right_censoring",
    x == 1L ~ "transition",
    x == 2L ~ "remission",
    TRUE ~ "Unexpected_code"
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

prepare_summary_dataset <- function(df, pnu_label, snu_label, days_per_year_value) {
  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  has_age_exact_followup <- "age_exact_followup" %in% names(df)

  out <- standardize_known_site_labels(df, pnu_label = pnu_label, snu_label = snu_label) %>%
    mutate(
      id = trimws(as.character(id)),
      sex_num = as.integer(coerce_numeric_text(sex_num)),
      age_exact_entry = coerce_numeric_text(age_exact_entry),
      days_followup = coerce_numeric_text(days_followup),
      status_num = as.integer(coerce_numeric_text(status_num))
    )

  if (has_age_exact_followup) {
    out <- out %>% mutate(age_exact_followup = coerce_numeric_text(age_exact_followup))
  } else {
    out <- out %>% mutate(age_exact_followup = NA_real_)
  }

  unexpected_sites <- setdiff(
    unique(out$site[!is.na(out$site) & out$site != ""]),
    c(pnu_label, snu_label)
  )

  if (length(unexpected_sites) > 0L) {
    stop(
      sprintf(
        "Unexpected site labels detected after standardization: %s",
        paste(sort(unexpected_sites), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!any(out$site == pnu_label, na.rm = TRUE)) {
    stop(sprintf("No rows found for site label `%s`.", pnu_label), call. = FALSE)
  }

  if (!any(out$site == snu_label, na.rm = TRUE)) {
    stop(sprintf("No rows found for site label `%s`.", snu_label), call. = FALSE)
  }

  out %>%
    mutate(
      site_id = dplyr::if_else(
        is.na(site) | site == "" | is.na(id) | id == "",
        NA_character_,
        paste(site, id, sep = "__")
      ),
      sex_label = classify_sex(sex_num),
      original_status_label = classify_original_status(status_num),
      age_exact_entry_clean = dplyr::if_else(!is.na(age_exact_entry) & age_exact_entry >= 0, age_exact_entry, NA_real_),
      age_exact_followup_clean = dplyr::if_else(!is.na(age_exact_followup) & age_exact_followup >= 0, age_exact_followup, NA_real_),
      days_followup_clean = dplyr::if_else(!is.na(days_followup) & days_followup >= 0, days_followup, NA_real_),
      years_followup_clean = days_followup_clean / days_per_year_value
    )
}

split_dataset_variants <- function(df, pnu_label, snu_label) {
  out <- list(
    PNU = df %>% filter(site == pnu_label),
    SNU = df %>% filter(site == snu_label),
    Merged = df
  )

  zero_row_names <- names(out)[vapply(out, nrow, integer(1)) == 0L]
  if (length(zero_row_names) > 0L) {
    stop(sprintf("Zero-row dataset variants detected: %s", paste(zero_row_names, collapse = ", ")), call. = FALSE)
  }

  out
}

apply_summary_scenario <- function(df, scenario_name) {
  scenario_name <- match.arg(scenario_name, choices = scenario_order)

  if (scenario_name == "observed_status") {
    return(
      df %>%
        mutate(
          analysis_scenario = scenario_name,
          scenario_label = unname(scenario_labels[[scenario_name]]),
          analysis_status_label = original_status_label,
          event_main = as.integer(status_num == 1L),
          censor_main = as.integer(status_num == 0L),
          remission_recoded_to_censoring = 0L
        )
    )
  }

  df %>%
    mutate(
      analysis_scenario = scenario_name,
      scenario_label = unname(scenario_labels[[scenario_name]]),
      analysis_status_label = dplyr::case_when(
        is.na(status_num) ~ "Missing",
        status_num == 1L ~ "transition",
        status_num %in% c(0L, 2L) ~ "censoring_including_remission",
        TRUE ~ "Unexpected_code"
      ),
      event_main = as.integer(status_num == 1L),
      censor_main = as.integer(status_num %in% c(0L, 2L)),
      remission_recoded_to_censoring = as.integer(status_num == 2L)
    )
}

initialize_master_row <- function(scenario_name, scenario_label, section, dataset, variable, stratum_type, stratum_value, summary_type, unit, definition) {
  tibble(
    analysis_scenario = scenario_name,
    scenario_label = scenario_label,
    section = section,
    dataset = dataset,
    variable = variable,
    stratum_type = stratum_type,
    stratum_value = stratum_value,
    summary_type = summary_type,
    unit = unit,
    definition = definition,
    n_total = NA_integer_,
    n_nonmissing = NA_integer_,
    n_missing = NA_integer_,
    count = NA_real_,
    denominator = NA_real_,
    percent = NA_real_,
    mean = NA_real_,
    sd = NA_real_,
    median = NA_real_,
    q1 = NA_real_,
    q3 = NA_real_,
    iqr = NA_real_,
    min = NA_real_,
    max = NA_real_,
    sum_value = NA_real_
  )
}

build_single_value_row <- function(scenario_name, scenario_label, section, dataset, variable, value, unit, definition, value_field = c("count", "sum_value"), stratum_type = "overall", stratum_value = "overall") {
  value_field <- match.arg(value_field)
  out <- initialize_master_row(
    scenario_name = scenario_name,
    scenario_label = scenario_label,
    section = section,
    dataset = dataset,
    variable = variable,
    stratum_type = stratum_type,
    stratum_value = stratum_value,
    summary_type = "single_value",
    unit = unit,
    definition = definition
  )
  out[[value_field]] <- as.numeric(value)
  out
}

build_count_percent_row <- function(scenario_name, scenario_label, section, dataset, variable, stratum_type, stratum_value, count, denominator, definition, summary_type = "count_summary", unit = "count") {
  out <- initialize_master_row(
    scenario_name = scenario_name,
    scenario_label = scenario_label,
    section = section,
    dataset = dataset,
    variable = variable,
    stratum_type = stratum_type,
    stratum_value = stratum_value,
    summary_type = summary_type,
    unit = unit,
    definition = definition
  )
  out$n_total <- as.integer(denominator)
  out$n_nonmissing <- as.integer(denominator)
  out$count <- as.numeric(count)
  out$denominator <- as.numeric(denominator)
  out$percent <- if (is.na(denominator) || denominator <= 0) NA_real_ else 100 * as.numeric(count) / as.numeric(denominator)
  out
}

build_missingness_row <- function(scenario_name, scenario_label, dataset, variable, missing_count, denominator, definition) {
  out <- initialize_master_row(
    scenario_name = scenario_name,
    scenario_label = scenario_label,
    section = "core_missingness",
    dataset = dataset,
    variable = variable,
    stratum_type = "overall",
    stratum_value = "overall",
    summary_type = "missingness",
    unit = "count",
    definition = definition
  )
  out$n_total <- as.integer(denominator)
  out$n_nonmissing <- as.integer(denominator - missing_count)
  out$n_missing <- as.integer(missing_count)
  out$count <- as.numeric(missing_count)
  out$denominator <- as.numeric(denominator)
  out$percent <- if (denominator <= 0) NA_real_ else 100 * as.numeric(missing_count) / as.numeric(denominator)
  out
}

build_numeric_summary_row <- function(scenario_name, scenario_label, section, dataset, variable, x, unit, definition, stratum_type = "overall", stratum_value = "overall") {
  x <- as.numeric(x)
  n_total_value <- length(x)
  n_nonmissing_value <- sum(!is.na(x))
  q1_value <- safe_quantile(x, prob = 0.25)
  q3_value <- safe_quantile(x, prob = 0.75)

  out <- initialize_master_row(
    scenario_name = scenario_name,
    scenario_label = scenario_label,
    section = section,
    dataset = dataset,
    variable = variable,
    stratum_type = stratum_type,
    stratum_value = stratum_value,
    summary_type = "numeric_summary",
    unit = unit,
    definition = definition
  )

  out$n_total <- as.integer(n_total_value)
  out$n_nonmissing <- as.integer(n_nonmissing_value)
  out$n_missing <- as.integer(n_total_value - n_nonmissing_value)
  out$mean <- safe_mean(x)
  out$sd <- safe_sd(x)
  out$median <- safe_median(x)
  out$q1 <- q1_value
  out$q3 <- q3_value
  out$iqr <- if (is.na(q1_value) || is.na(q3_value)) NA_real_ else q3_value - q1_value
  out$min <- safe_min(x)
  out$max <- safe_max(x)
  out$sum_value <- safe_sum(x)
  out
}

build_distribution_rows <- function(scenario_name, scenario_label, dataset, section, variable, values, report_levels, denominator, definition_prefix) {
  bind_rows(lapply(report_levels, function(level_name) {
    build_count_percent_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label,
      section = section,
      dataset = dataset,
      variable = variable,
      stratum_type = variable,
      stratum_value = level_name,
      count = sum(values == level_name, na.rm = TRUE),
      denominator = denominator,
      definition = sprintf("%s: %s", definition_prefix, level_name)
    )
  }))
}

create_dataset_master_rows <- function(df, dataset_name, scenario_name) {
  scenario_label_value <- unique(df$scenario_label)
  if (length(scenario_label_value) != 1L) {
    stop(sprintf("[%s][%s] Scenario label must be unique within dataset rows.", scenario_name, dataset_name), call. = FALSE)
  }
  scenario_label_value <- scenario_label_value[[1]]

  n_total <- nrow(df)
  nonmissing_site_id <- df$site_id[!is.na(df$site_id)]
  unique_site_id_n <- dplyr::n_distinct(nonmissing_site_id)
  duplicated_site_id_extra_n <- sum(duplicated(nonmissing_site_id))

  missingness_variables <- c("id", "site", "sex_num", "age_exact_entry", "age_exact_followup", "days_followup", "status_num")

  missingness_rows <- bind_rows(lapply(missingness_variables, function(variable_name) {
    build_missingness_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label_value,
      dataset = dataset_name,
      variable = variable_name,
      missing_count = count_missing_like(df[[variable_name]]),
      denominator = n_total,
      definition = sprintf("Missing values for `%s` in the source data.", variable_name)
    )
  }))

  cohort_rows <- bind_rows(
    build_single_value_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label_value,
      section = "cohort_overview",
      dataset = dataset_name,
      variable = "n_subjects",
      value = n_total,
      unit = "subjects",
      definition = "Number of subjects (rows) in the dataset variant.",
      value_field = "count"
    ),
    build_single_value_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label_value,
      section = "cohort_overview",
      dataset = dataset_name,
      variable = "n_unique_site_id",
      value = unique_site_id_n,
      unit = "subjects",
      definition = "Number of unique subjects using `site + id` as the subject key.",
      value_field = "count"
    ),
    build_single_value_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label_value,
      section = "cohort_overview",
      dataset = dataset_name,
      variable = "n_duplicated_site_id_extra_rows",
      value = duplicated_site_id_extra_n,
      unit = "rows",
      definition = "Number of extra duplicated rows beyond the first row within duplicated `site + id`.",
      value_field = "count"
    ),
    build_single_value_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label_value,
      section = "cohort_overview",
      dataset = dataset_name,
      variable = "total_observed_followup_days",
      value = safe_sum(df$days_followup_clean),
      unit = "days",
      definition = "Total observed follow-up across all rows, using non-negative follow-up only.",
      value_field = "sum_value"
    ),
    build_single_value_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label_value,
      section = "cohort_overview",
      dataset = dataset_name,
      variable = "total_observed_followup_years",
      value = safe_sum(df$years_followup_clean),
      unit = "person_years",
      definition = "Total observed follow-up across all rows, expressed in person-years.",
      value_field = "sum_value"
    )
  )

  flag_rows <- bind_rows(
    build_count_percent_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label_value,
      section = "data_flags",
      dataset = dataset_name,
      variable = "sex_num_invalid_code",
      stratum_type = "flag",
      stratum_value = "sex_num_invalid_code",
      count = sum(!is.na(df$sex_num) & !df$sex_num %in% c(0L, 1L)),
      denominator = n_total,
      definition = "Rows with `sex_num` not coded as 0 or 1.",
      summary_type = "flag_count"
    ),
    build_count_percent_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label_value,
      section = "data_flags",
      dataset = dataset_name,
      variable = "status_num_invalid_code",
      stratum_type = "flag",
      stratum_value = "status_num_invalid_code",
      count = sum(!is.na(df$status_num) & !df$status_num %in% c(0L, 1L, 2L)),
      denominator = n_total,
      definition = "Rows with `status_num` not coded as 0, 1, or 2.",
      summary_type = "flag_count"
    ),
    build_count_percent_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label_value,
      section = "data_flags",
      dataset = dataset_name,
      variable = "age_exact_entry_negative",
      stratum_type = "flag",
      stratum_value = "age_exact_entry_negative",
      count = sum(!is.na(df$age_exact_entry) & df$age_exact_entry < 0),
      denominator = n_total,
      definition = "Rows with negative `age_exact_entry`.",
      summary_type = "flag_count"
    ),
    build_count_percent_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label_value,
      section = "data_flags",
      dataset = dataset_name,
      variable = "age_exact_followup_negative",
      stratum_type = "flag",
      stratum_value = "age_exact_followup_negative",
      count = sum(!is.na(df$age_exact_followup) & df$age_exact_followup < 0),
      denominator = n_total,
      definition = "Rows with negative `age_exact_followup`.",
      summary_type = "flag_count"
    ),
    build_count_percent_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label_value,
      section = "data_flags",
      dataset = dataset_name,
      variable = "days_followup_negative",
      stratum_type = "flag",
      stratum_value = "days_followup_negative",
      count = sum(!is.na(df$days_followup) & df$days_followup < 0),
      denominator = n_total,
      definition = "Rows with negative `days_followup`.",
      summary_type = "flag_count"
    )
  )

  sex_rows <- build_distribution_rows(
    scenario_name = scenario_name,
    scenario_label = scenario_label_value,
    dataset = dataset_name,
    section = "sex_distribution",
    variable = "sex_label",
    values = df$sex_label,
    report_levels = sex_report_levels,
    denominator = n_total,
    definition_prefix = "Sex distribution"
  )

  original_status_rows <- build_distribution_rows(
    scenario_name = scenario_name,
    scenario_label = scenario_label_value,
    dataset = dataset_name,
    section = "original_status_distribution",
    variable = "original_status_label",
    values = df$original_status_label,
    report_levels = original_status_levels,
    denominator = n_total,
    definition_prefix = "Original observed outcome status distribution"
  )

  scenario_analysis_levels <- analysis_status_levels[[scenario_name]]
  analysis_status_rows <- build_distribution_rows(
    scenario_name = scenario_name,
    scenario_label = scenario_label_value,
    dataset = dataset_name,
    section = "analysis_status_distribution",
    variable = "analysis_status_label",
    values = df$analysis_status_label,
    report_levels = scenario_analysis_levels,
    denominator = n_total,
    definition_prefix = "Scenario-specific outcome status distribution"
  )

  if (scenario_name == "remission_censored") {
    analysis_status_rows <- bind_rows(
      analysis_status_rows,
      build_count_percent_row(
        scenario_name = scenario_name,
        scenario_label = scenario_label_value,
        section = "analysis_status_distribution",
        dataset = dataset_name,
        variable = "remission_recoded_to_censoring",
        stratum_type = "analysis_note",
        stratum_value = "remission_recoded_to_censoring",
        count = sum(df$remission_recoded_to_censoring == 1L, na.rm = TRUE),
        denominator = n_total,
        definition = "Rows with original remission recoded to censoring in the transition-only scenario.",
        summary_type = "flag_count"
      )
    )
  }

  age_rows <- bind_rows(
    build_numeric_summary_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label_value,
      section = "age_summary",
      dataset = dataset_name,
      variable = "age_exact_entry_clean",
      x = df$age_exact_entry_clean,
      unit = "years",
      definition = "Distribution of age at entry, using non-negative values only."
    ),
    build_numeric_summary_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label_value,
      section = "age_summary",
      dataset = dataset_name,
      variable = "age_exact_followup_clean",
      x = df$age_exact_followup_clean,
      unit = "years",
      definition = "Distribution of age at the end of follow-up, using non-negative values only."
    )
  )

  followup_rows <- bind_rows(
    build_numeric_summary_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label_value,
      section = "followup_summary",
      dataset = dataset_name,
      variable = "days_followup_clean",
      x = df$days_followup_clean,
      unit = "days",
      definition = "Distribution of observed follow-up time in days, using non-negative values only."
    ),
    build_numeric_summary_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label_value,
      section = "followup_summary",
      dataset = dataset_name,
      variable = "years_followup_clean",
      x = df$years_followup_clean,
      unit = "years",
      definition = "Distribution of observed follow-up time in years, using non-negative values only."
    )
  )

  bind_rows(
    cohort_rows,
    missingness_rows,
    flag_rows,
    sex_rows,
    original_status_rows,
    analysis_status_rows,
    age_rows,
    followup_rows
  )
}

fmt_number <- function(x, digits = 2) {
  if (length(x) == 0L || is.na(x)) return("NA")
  formatC(x, format = "f", digits = digits, big.mark = ",")
}

fmt_integer <- function(x) {
  if (length(x) == 0L || is.na(x)) return("NA")
  formatC(round(x), format = "f", digits = 0, big.mark = ",")
}

format_count_percent <- function(count, percent, digits_percent = 1) {
  if (is.na(count)) return("NA")
  paste0(fmt_integer(count), " (", fmt_number(percent, digits = digits_percent), "%)")
}

format_mean_sd <- function(mean_value, sd_value, digits = 2) {
  if (is.na(mean_value) || is.na(sd_value)) return("NA")
  paste0(fmt_number(mean_value, digits = digits), " +- ", fmt_number(sd_value, digits = digits))
}

format_median_iqr <- function(median_value, q1_value, q3_value, digits = 2) {
  if (is.na(median_value) || is.na(q1_value) || is.na(q3_value)) return("NA")
  paste0(fmt_number(median_value, digits = digits), " [", fmt_number(q1_value, digits = digits), ", ", fmt_number(q3_value, digits = digits), "]")
}

format_range <- function(min_value, max_value, digits = 2) {
  if (is.na(min_value) || is.na(max_value)) return("NA")
  paste0(fmt_number(min_value, digits = digits), " to ", fmt_number(max_value, digits = digits))
}

make_followup_boxplot_status_label <- function(x) {
  dplyr::case_when(
    x == "right_censoring" ~ "Right censoring",
    x == "transition" ~ "Transition",
    x == "remission" ~ "Remission",
    x == "censoring_including_remission" ~ "Censoring + remission",
    x == "Missing" ~ "Missing",
    x == "Unexpected_code" ~ "Unexpected code",
    TRUE ~ as.character(x)
  )
}

build_followup_boxplot_data <- function(analysis_ready_df, scenario_name, pnu_label, snu_label) {
  scenario_df <- apply_summary_scenario(analysis_ready_df, scenario_name)
  dataset_variants <- split_dataset_variants(
    df = scenario_df,
    pnu_label = pnu_label,
    snu_label = snu_label
  )

  bind_rows(lapply(names(dataset_variants), function(dataset_name) {
    dataset_variants[[dataset_name]] %>%
      transmute(
        analysis_scenario = scenario_name,
        scenario_label = unname(scenario_labels[[scenario_name]]),
        dataset = dataset_name,
        analysis_status_label,
        plot_status_label = make_followup_boxplot_status_label(analysis_status_label),
        years_followup_clean
      )
  }))
}

prepare_followup_boxplot_display_df <- function(plot_df) {
  plot_df %>%
    filter(
      !is.na(years_followup_clean),
      plot_status_label %in% c("Right censoring", "Transition", "Remission", "Censoring + remission")
    ) %>%
    mutate(
      analysis_scenario = factor(analysis_scenario, levels = scenario_order),
      dataset = factor(dataset, levels = dataset_order),
      plot_dataset_label = factor(as.character(dataset), levels = dataset_order),
      plot_status_label = factor(
        plot_status_label,
        levels = c("Right censoring", "Transition", "Remission", "Censoring + remission")
      ),
      scenario_short_label = dplyr::case_when(
        analysis_scenario == "observed_status" ~ "Observed status",
        analysis_scenario == "remission_censored" ~ "Remission as censoring",
        TRUE ~ as.character(analysis_scenario)
      ),
      plot_panel_label = paste0(scenario_short_label, "\n", as.character(plot_status_label))
    ) %>%
    mutate(
      plot_panel_label = factor(plot_panel_label, levels = followup_panel_order)
    )
}

build_followup_boxplot_base <- function(plot_df) {
  ggplot(plot_df, aes(x = plot_dataset_label, y = years_followup_clean, color = dataset)) +
    geom_boxplot(
      linewidth = 0.6,
      width = 0.55,
      outlier.shape = NA
    ) +
    geom_jitter(width = 0.12, alpha = 0.35, size = 1.1) +
    scale_color_manual(
      values = followup_dataset_palette,
      breaks = dataset_order,
      labels = dataset_order
    ) +
    labs(
      title = "Follow-up time distribution by dataset",
      x = NULL,
      y = "Follow-up time (years)",
      color = "Dataset"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "#F2F2F2", colour = "#D0D0D0"),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 12, hjust = 1)
    )
}

build_followup_boxplot <- function(plot_df) {
  plot_df <- prepare_followup_boxplot_display_df(plot_df)

  build_followup_boxplot_base(plot_df) +
    facet_wrap(~ plot_panel_label, ncol = 3, scales = "free_x") +
    labs(
      subtitle = "Styled after the stage7 mixture-cure plots; boxplot with overlaid subject-level jitter"
    )
}

build_followup_single_panel_boxplot <- function(plot_df, panel_label) {
  plot_df <- plot_df %>%
    prepare_followup_boxplot_display_df() %>%
    filter(as.character(plot_panel_label) == panel_label)

  build_followup_boxplot_base(plot_df) +
    labs(
      subtitle = unname(followup_panel_title_lookup[[panel_label]])
    )
}

build_followup_boxplot_panel_registry <- function(plot_df) {
  prepared_df <- prepare_followup_boxplot_display_df(plot_df)

  tibble(
    plot_panel_label = followup_panel_order,
    panel_file_stub = unname(followup_panel_file_stub_lookup[followup_panel_order]),
    panel_title = unname(followup_panel_title_lookup[followup_panel_order])
  ) %>%
    filter(.data$plot_panel_label %in% unique(as.character(prepared_df$plot_panel_label)))
}

build_followup_jitter_plot <- function(plot_df) {
  plot_fill_values <- c(
    "Right censoring" = "#B8C0C7",
    "Transition" = "#E07A5F",
    "Remission" = "#5C9E8F",
    "Censoring + remission" = "#4F6D7A"
  )

  plot_df <- plot_df %>%
    filter(
      !is.na(years_followup_clean),
      plot_status_label %in% names(plot_fill_values)
    ) %>%
    mutate(
      analysis_scenario = factor(analysis_scenario, levels = scenario_order),
      dataset = factor(dataset, levels = dataset_order),
      plot_status_label = factor(plot_status_label, levels = names(plot_fill_values))
    )

  ggplot(
    plot_df,
    aes(x = dataset, y = years_followup_clean, fill = plot_status_label)
  ) +
    geom_point(
      shape = 21,
      size = 2.2,
      alpha = 0.7,
      stroke = 0.2,
      position = position_jitterdodge(
        jitter.width = 0.16,
        jitter.height = 0,
        dodge.width = 0.82,
        seed = 42
      )
    ) +
    facet_wrap(~ scenario_label, ncol = 1, scales = "free_x") +
    scale_fill_manual(
      values = plot_fill_values,
      drop = TRUE,
      name = "Status coding"
    ) +
    labs(
      title = "Follow-up time distribution by dataset",
      subtitle = "Each point is one subject; descriptive jitter plot only",
      x = NULL,
      y = "Follow-up time (years)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "#F2F2F2", colour = "#D0D0D0"),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(face = "bold")
    )
}

get_master_row <- function(master_table, scenario_name, dataset_name, section_name, variable_name, stratum_type_name = "overall", stratum_value_name = "overall") {
  out <- master_table %>%
    filter(
      .data$analysis_scenario == scenario_name,
      .data$dataset == dataset_name,
      .data$section == section_name,
      .data$variable == variable_name,
      .data$stratum_type == stratum_type_name,
      .data$stratum_value == stratum_value_name
    )

  if (nrow(out) == 0L) {
    return(NULL)
  }

  out %>% slice(1)
}

get_distribution_total <- function(master_table, scenario_name, dataset_name, section_name, variable_name, stratum_values) {
  out <- master_table %>%
    filter(
      .data$analysis_scenario == scenario_name,
      .data$dataset == dataset_name,
      .data$section == section_name,
      .data$variable == variable_name,
      .data$stratum_value %in% stratum_values
    )

  if (nrow(out) == 0L) {
    return(tibble(count = NA_real_, percent = NA_real_, denominator = NA_real_))
  }

  denominator_value <- out$denominator[match(TRUE, !is.na(out$denominator))]
  if (length(denominator_value) == 0L || is.na(denominator_value)) {
    denominator_value <- NA_real_
  }

  count_value <- sum(out$count, na.rm = TRUE)
  percent_value <- if (is.na(denominator_value) || denominator_value <= 0) NA_real_ else 100 * count_value / denominator_value

  tibble(
    count = count_value,
    percent = percent_value,
    denominator = denominator_value
  )
}

make_compact_row <- function(scenario_name, scenario_label, panel, row_label, pnu_value, snu_value, merged_value) {
  tibble(
    analysis_scenario = scenario_name,
    scenario_label = scenario_label,
    panel = panel,
    row_label = row_label,
    PNU = pnu_value,
    SNU = snu_value,
    Merged = merged_value
  )
}

build_compact_table <- function(master_table, scenario_name) {
  scenario_label_value <- unname(scenario_labels[[scenario_name]])
  compact_rows <- list()
  row_index <- 0L

  add_compact_row <- function(panel, row_label, value_builder) {
    row_index <<- row_index + 1L
    values <- lapply(dataset_order, value_builder)
    compact_rows[[row_index]] <<- make_compact_row(
      scenario_name = scenario_name,
      scenario_label = scenario_label_value,
      panel = panel,
      row_label = row_label,
      pnu_value = values[[1]],
      snu_value = values[[2]],
      merged_value = values[[3]]
    )
  }

  has_age_followup_summary <- any(
    master_table$analysis_scenario == scenario_name &
      master_table$section == "age_summary" &
      master_table$variable == "age_exact_followup_clean" &
      !is.na(master_table$n_nonmissing) &
      master_table$n_nonmissing > 0L
  )

  has_sex_missing_or_invalid <- any(
    master_table$analysis_scenario == scenario_name &
      master_table$section == "sex_distribution" &
      master_table$variable == "sex_label" &
      master_table$stratum_value %in% c("Missing", "Unexpected_code") &
      !is.na(master_table$count) &
      master_table$count > 0
  )

  has_status_missing_or_invalid <- any(
    master_table$analysis_scenario == scenario_name &
      master_table$section == "analysis_status_distribution" &
      master_table$variable == "analysis_status_label" &
      master_table$stratum_value %in% c("Missing", "Unexpected_code") &
      !is.na(master_table$count) &
      master_table$count > 0
  )

  add_compact_row("Cohort", "N subjects", function(dataset_name) {
    out <- get_master_row(master_table, scenario_name, dataset_name, "cohort_overview", "n_subjects")
    if (is.null(out)) return("NA")
    fmt_integer(out$count)
  })

  add_compact_row("Sex", "Female, n (%)", function(dataset_name) {
    out <- get_master_row(master_table, scenario_name, dataset_name, "sex_distribution", "sex_label", "sex_label", "Female")
    if (is.null(out)) return("NA")
    format_count_percent(out$count, out$percent, digits_percent = 1)
  })

  add_compact_row("Sex", "Male, n (%)", function(dataset_name) {
    out <- get_master_row(master_table, scenario_name, dataset_name, "sex_distribution", "sex_label", "sex_label", "Male")
    if (is.null(out)) return("NA")
    format_count_percent(out$count, out$percent, digits_percent = 1)
  })

  if (has_sex_missing_or_invalid) {
    add_compact_row("Sex", "Missing/invalid sex, n (%)", function(dataset_name) {
      out <- get_distribution_total(master_table, scenario_name, dataset_name, "sex_distribution", "sex_label", c("Missing", "Unexpected_code"))
      format_count_percent(out$count, out$percent, digits_percent = 1)
    })
  }

  add_compact_row("Age at entry (years)", "Mean +- SD", function(dataset_name) {
    out <- get_master_row(master_table, scenario_name, dataset_name, "age_summary", "age_exact_entry_clean")
    if (is.null(out)) return("NA")
    format_mean_sd(out$mean, out$sd, digits = 2)
  })

  add_compact_row("Age at entry (years)", "Median [Q1, Q3]", function(dataset_name) {
    out <- get_master_row(master_table, scenario_name, dataset_name, "age_summary", "age_exact_entry_clean")
    if (is.null(out)) return("NA")
    format_median_iqr(out$median, out$q1, out$q3, digits = 2)
  })

  add_compact_row("Age at entry (years)", "Min to max", function(dataset_name) {
    out <- get_master_row(master_table, scenario_name, dataset_name, "age_summary", "age_exact_entry_clean")
    if (is.null(out)) return("NA")
    format_range(out$min, out$max, digits = 2)
  })

  if (has_age_followup_summary) {
    add_compact_row("Age at end of follow-up (years)", "Mean +- SD", function(dataset_name) {
      out <- get_master_row(master_table, scenario_name, dataset_name, "age_summary", "age_exact_followup_clean")
      if (is.null(out)) return("NA")
      format_mean_sd(out$mean, out$sd, digits = 2)
    })

    add_compact_row("Age at end of follow-up (years)", "Median [Q1, Q3]", function(dataset_name) {
      out <- get_master_row(master_table, scenario_name, dataset_name, "age_summary", "age_exact_followup_clean")
      if (is.null(out)) return("NA")
      format_median_iqr(out$median, out$q1, out$q3, digits = 2)
    })

    add_compact_row("Age at end of follow-up (years)", "Min to max", function(dataset_name) {
      out <- get_master_row(master_table, scenario_name, dataset_name, "age_summary", "age_exact_followup_clean")
      if (is.null(out)) return("NA")
      format_range(out$min, out$max, digits = 2)
    })
  }

  add_compact_row("Follow-up duration (days)", "Mean +- SD", function(dataset_name) {
    out <- get_master_row(master_table, scenario_name, dataset_name, "followup_summary", "days_followup_clean")
    if (is.null(out)) return("NA")
    format_mean_sd(out$mean, out$sd, digits = 1)
  })

  add_compact_row("Follow-up duration (days)", "Median [Q1, Q3]", function(dataset_name) {
    out <- get_master_row(master_table, scenario_name, dataset_name, "followup_summary", "days_followup_clean")
    if (is.null(out)) return("NA")
    format_median_iqr(out$median, out$q1, out$q3, digits = 1)
  })

  add_compact_row("Follow-up duration (days)", "Min to max", function(dataset_name) {
    out <- get_master_row(master_table, scenario_name, dataset_name, "followup_summary", "days_followup_clean")
    if (is.null(out)) return("NA")
    format_range(out$min, out$max, digits = 1)
  })

  add_compact_row("Follow-up duration (years)", "Mean +- SD", function(dataset_name) {
    out <- get_master_row(master_table, scenario_name, dataset_name, "followup_summary", "years_followup_clean")
    if (is.null(out)) return("NA")
    format_mean_sd(out$mean, out$sd, digits = 2)
  })

  add_compact_row("Follow-up duration (years)", "Median [Q1, Q3]", function(dataset_name) {
    out <- get_master_row(master_table, scenario_name, dataset_name, "followup_summary", "years_followup_clean")
    if (is.null(out)) return("NA")
    format_median_iqr(out$median, out$q1, out$q3, digits = 2)
  })

  add_compact_row("Follow-up duration (years)", "Min to max", function(dataset_name) {
    out <- get_master_row(master_table, scenario_name, dataset_name, "followup_summary", "years_followup_clean")
    if (is.null(out)) return("NA")
    format_range(out$min, out$max, digits = 2)
  })

  if (scenario_name == "observed_status") {
    add_compact_row("Outcome coding used in this table", "Right censoring, n (%)", function(dataset_name) {
      out <- get_master_row(master_table, scenario_name, dataset_name, "analysis_status_distribution", "analysis_status_label", "analysis_status_label", "right_censoring")
      if (is.null(out)) return("NA")
      format_count_percent(out$count, out$percent, digits_percent = 1)
    })

    add_compact_row("Outcome coding used in this table", "Transition, n (%)", function(dataset_name) {
      out <- get_master_row(master_table, scenario_name, dataset_name, "analysis_status_distribution", "analysis_status_label", "analysis_status_label", "transition")
      if (is.null(out)) return("NA")
      format_count_percent(out$count, out$percent, digits_percent = 1)
    })

    add_compact_row("Outcome coding used in this table", "Remission, n (%)", function(dataset_name) {
      out <- get_master_row(master_table, scenario_name, dataset_name, "analysis_status_distribution", "analysis_status_label", "analysis_status_label", "remission")
      if (is.null(out)) return("NA")
      format_count_percent(out$count, out$percent, digits_percent = 1)
    })
  } else {
    add_compact_row("Outcome coding used in this table", "Censoring incl. remission, n (%)", function(dataset_name) {
      out <- get_master_row(master_table, scenario_name, dataset_name, "analysis_status_distribution", "analysis_status_label", "analysis_status_label", "censoring_including_remission")
      if (is.null(out)) return("NA")
      format_count_percent(out$count, out$percent, digits_percent = 1)
    })

    add_compact_row("Outcome coding used in this table", "Transition, n (%)", function(dataset_name) {
      out <- get_master_row(master_table, scenario_name, dataset_name, "analysis_status_distribution", "analysis_status_label", "analysis_status_label", "transition")
      if (is.null(out)) return("NA")
      format_count_percent(out$count, out$percent, digits_percent = 1)
    })

    add_compact_row("Outcome coding used in this table", "Original remission recoded, n (%)", function(dataset_name) {
      out <- get_master_row(master_table, scenario_name, dataset_name, "analysis_status_distribution", "remission_recoded_to_censoring", "analysis_note", "remission_recoded_to_censoring")
      if (is.null(out)) return("NA")
      format_count_percent(out$count, out$percent, digits_percent = 1)
    })
  }

  if (has_status_missing_or_invalid) {
    add_compact_row("Outcome coding used in this table", "Missing/invalid status, n (%)", function(dataset_name) {
      out <- get_distribution_total(master_table, scenario_name, dataset_name, "analysis_status_distribution", "analysis_status_label", c("Missing", "Unexpected_code"))
      format_count_percent(out$count, out$percent, digits_percent = 1)
    })
  }

  add_compact_row("Observed follow-up", "Total observed person-years", function(dataset_name) {
    out <- get_master_row(master_table, scenario_name, dataset_name, "cohort_overview", "total_observed_followup_years")
    if (is.null(out)) return("NA")
    fmt_number(out$sum_value, digits = 2)
  })

  bind_rows(compact_rows)
}

# Read and transform source data ------------------------------------------------
raw_df <- read_input_dataset(source_file)
analysis_ready_df <- prepare_summary_dataset(
  df = raw_df,
  pnu_label = pnu_site_label,
  snu_label = snu_site_label,
  days_per_year_value = days_per_year
)

scenario_master_tables <- lapply(scenario_order, function(scenario_name) {
  scenario_df <- apply_summary_scenario(analysis_ready_df, scenario_name)
  dataset_variants <- split_dataset_variants(
    df = scenario_df,
    pnu_label = pnu_site_label,
    snu_label = snu_site_label
  )

  bind_rows(lapply(names(dataset_variants), function(dataset_name) {
    create_dataset_master_rows(
      df = dataset_variants[[dataset_name]],
      dataset_name = dataset_name,
      scenario_name = scenario_name
    )
  }))
})

master_table <- bind_rows(scenario_master_tables) %>%
  mutate(
    analysis_scenario = factor(analysis_scenario, levels = scenario_order),
    dataset = factor(dataset, levels = dataset_order),
    section = factor(section, levels = section_order)
  ) %>%
  arrange(analysis_scenario, section, dataset, variable, stratum_type, stratum_value) %>%
  mutate(
    analysis_scenario = as.character(analysis_scenario),
    dataset = as.character(dataset),
    section = as.character(section)
  )

compact_table <- bind_rows(lapply(scenario_order, function(scenario_name) {
  build_compact_table(master_table, scenario_name = scenario_name)
}))

compact_observed_table <- compact_table %>% filter(analysis_scenario == "observed_status")
compact_remission_censored_table <- compact_table %>% filter(analysis_scenario == "remission_censored")

followup_boxplot_data <- bind_rows(lapply(scenario_order, function(scenario_name) {
  build_followup_boxplot_data(
    analysis_ready_df = analysis_ready_df,
    scenario_name = scenario_name,
    pnu_label = pnu_site_label,
    snu_label = snu_site_label
  )
}))

followup_boxplot <- build_followup_boxplot(followup_boxplot_data)
followup_boxplot_panel_registry <- build_followup_boxplot_panel_registry(followup_boxplot_data)
followup_jitter_plot <- build_followup_jitter_plot(followup_boxplot_data)

metadata_table <- tibble(
  item = c(
    "source_file",
    "dictionary_file",
    "export_path",
    "created_at",
    "days_per_year",
    "scenario_observed_status",
    "scenario_remission_censored"
  ),
  value = c(
    source_file,
    dictionary_file,
    export_path,
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    as.character(days_per_year),
    scenario_labels[["observed_status"]],
    scenario_labels[["remission_censored"]]
  )
)

# Export tables -----------------------------------------------------------------
readr::write_csv(master_table, file.path(export_path, master_output_file))
readr::write_csv(compact_table, file.path(export_path, compact_output_file))
readr::write_csv(compact_observed_table, file.path(export_path, compact_observed_output_file))
readr::write_csv(compact_remission_censored_table, file.path(export_path, compact_remission_censored_output_file))

ggplot2::ggsave(
  filename = file.path(export_path, followup_boxplot_png_file),
  plot = followup_boxplot,
  width = 12,
  height = 9,
  dpi = 320
)

ggplot2::ggsave(
  filename = file.path(export_path, followup_boxplot_pdf_file),
  plot = followup_boxplot,
  width = 12,
  height = 9
)

followup_boxplot_panel_export_dir <- file.path(export_path, followup_boxplot_panel_dir_name)
dir.create(followup_boxplot_panel_export_dir, recursive = TRUE, showWarnings = FALSE)

invisible(lapply(seq_len(nrow(followup_boxplot_panel_registry)), function(ii) {
  panel_label <- followup_boxplot_panel_registry$plot_panel_label[[ii]]
  panel_file_stub <- followup_boxplot_panel_registry$panel_file_stub[[ii]]
  panel_plot <- build_followup_single_panel_boxplot(
    plot_df = followup_boxplot_data,
    panel_label = panel_label
  )

  ggplot2::ggsave(
    filename = file.path(
      followup_boxplot_panel_export_dir,
      paste0("PLOT__data_summary__followup_time_boxplot__", panel_file_stub, ".png")
    ),
    plot = panel_plot,
    width = 7,
    height = 5,
    dpi = 320
  )

  ggplot2::ggsave(
    filename = file.path(
      followup_boxplot_panel_export_dir,
      paste0("PLOT__data_summary__followup_time_boxplot__", panel_file_stub, ".pdf")
    ),
    plot = panel_plot,
    width = 7,
    height = 5
  )

  invisible(NULL)
}))

ggplot2::ggsave(
  filename = file.path(export_path, followup_jitter_png_file),
  plot = followup_jitter_plot,
  width = 12,
  height = 9,
  dpi = 320
)

ggplot2::ggsave(
  filename = file.path(export_path, followup_jitter_pdf_file),
  plot = followup_jitter_plot,
  width = 12,
  height = 9
)

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    x = list(
      metadata = metadata_table,
      master = master_table,
      compact = compact_table,
      compact_observed = compact_observed_table,
      compact_remission_censored = compact_remission_censored_table
    ),
    path = file.path(export_path, workbook_output_file)
  )
}

message("Data summary tables written to: ", export_path)
