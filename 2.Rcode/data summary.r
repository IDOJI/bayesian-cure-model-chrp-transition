# 🔴 Configure: paths and pre-analysis data-inspection settings ===============================
merged_file <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/data/MERGED_dataset3_pnu_snu.csv'
export_path <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results/data summary'

pnu_site_label <- "PNU"
snu_site_label <- "SNU"

main_risk_scale <- "transition_only_main"
supplementary_risk_scale <- "transition_cif_competing"

days_per_year <- 365.25
output_stem <- "DESC__preanalysis_clinical_survival_data_inspection__PNU_SNU_merged"

# 🔴 Initialize: packages and runtime options ===============================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(lubridate)
  library(tibble)
  library(purrr)
  library(stringr)
  library(survival)
})

options(stringsAsFactors = FALSE, scipen = 999)

dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

dataset_order <- c("PNU", "SNU", "merged")
section_order <- c(
  "cohort_overview",
  "sex_distribution",
  "baseline_age_summary",
  "age_by_sex",
  "event_distribution",
  "event_rate_summary",
  "followup_numeric_summary",
  "followup_reverse_km_summary",
  "status_by_sex",
  "calendar_distribution"
)

# 🔴 Define: descriptive helper functions ===============================
## 🟠 Define: input and coercion utilities ===============================
read_input_dataset <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Input file does not exist: %s", path), call. = FALSE)
  }
  
  ext <- tolower(tools::file_ext(path))
  
  if (ext %in% c("csv", "txt")) {
    return(
      readr::read_csv(
        file = path,
        col_types = readr::cols(.default = readr::col_character()),
        show_col_types = FALSE,
        progress = FALSE
      )
    )
  }
  
  if (ext == "rds") {
    return(readRDS(path))
  }
  
  stop(sprintf("Unsupported input extension for `%s`.", path), call. = FALSE)
}

coerce_numeric_text <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

parse_date_ymd_safe <- function(x) {
  y <- suppressWarnings(lubridate::ymd(x, quiet = TRUE))
  if (all(is.na(y))) {
    y <- suppressWarnings(as.Date(x))
  }
  y
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

## 🟠 Define: backbone validation and derivation ===============================
prepare_backbone_dataset <- function(df, pnu_label, snu_label) {
  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  
  df <- standardize_known_site_labels(df, pnu_label = pnu_label, snu_label = snu_label) %>%
    mutate(
      id = trimws(as.character(id)),
      sex_num = as.integer(coerce_numeric_text(sex_num)),
      age_exact_entry = coerce_numeric_text(age_exact_entry),
      days_followup = coerce_numeric_text(days_followup),
      status_num = as.integer(coerce_numeric_text(status_num))
    )
  
  if (nrow(df) == 0) {
    stop("Merged input has zero rows.", call. = FALSE)
  }
  
  if (anyNA(df[c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")])) {
    stop("Missing values detected in required backbone columns.", call. = FALSE)
  }
  
  if (any(df$id == "")) {
    stop("Blank `id` values detected.", call. = FALSE)
  }
  
  if (any(df$site == "")) {
    stop("Blank `site` values detected.", call. = FALSE)
  }
  
  if (any(!df$site %in% c(pnu_label, snu_label))) {
    stop("Unexpected `site` labels detected after standardization.", call. = FALSE)
  }
  
  if (any(!df$sex_num %in% c(0L, 1L))) {
    stop("`sex_num` must be coded as 0/1 only.", call. = FALSE)
  }
  
  if (any(!df$status_num %in% c(0L, 1L, 2L))) {
    stop("`status_num` must be coded as 0/1/2 only.", call. = FALSE)
  }
  
  if (any(df$days_followup < 0)) {
    stop("Negative `days_followup` values detected.", call. = FALSE)
  }
  
  if ("date_entry" %in% names(df)) {
    df <- df %>% mutate(date_entry = parse_date_ymd_safe(date_entry))
    if (anyNA(df$date_entry)) {
      stop("`date_entry` could not be parsed completely as Date.", call. = FALSE)
    }
  }
  
  df <- df %>%
    mutate(
      site_id = paste(site, id, sep = "__"),
      sex_label = case_when(
        sex_num == 0L ~ "Female",
        sex_num == 1L ~ "Male",
        TRUE ~ NA_character_
      ),
      status_label = case_when(
        status_num == 0L ~ "right_censoring",
        status_num == 1L ~ "transition",
        status_num == 2L ~ "remission",
        TRUE ~ NA_character_
      ),
      years_followup = days_followup / days_per_year,
      is_transition = status_num == 1L,
      is_remission = status_num == 2L,
      is_right_censoring = status_num == 0L,
      is_main_transition_censor = status_num != 1L
    )

  if ("date_entry" %in% names(df)) {
    df <- df %>%
      mutate(
        date_followup = date_entry + lubridate::days(as.integer(round(days_followup))),
        entry_year = lubridate::year(date_entry),
        end_year = lubridate::year(date_followup)
      )
  }

  df
}

split_dataset_variants <- function(df, pnu_label, snu_label) {
  out <- list(
    PNU = df %>% filter(site == pnu_label),
    SNU = df %>% filter(site == snu_label),
    merged = df
  )
  
  zero_names <- names(out)[vapply(out, nrow, integer(1)) == 0L]
  if (length(zero_names) > 0) {
    stop(sprintf("Zero-row dataset variants detected: %s", paste(zero_names, collapse = ", ")), call. = FALSE)
  }
  
  out
}

## 🟠 Define: row builders for source-of-truth export ===============================
empty_numeric_fields <- function() {
  tibble(
    n_total = NA_integer_,
    n_nonmissing = NA_integer_,
    count = NA_real_,
    denominator = NA_real_,
    proportion = NA_real_,
    percent = NA_real_,
    person_years = NA_real_,
    rate_per100py = NA_real_,
    rate_lcl95_per100py = NA_real_,
    rate_ucl95_per100py = NA_real_,
    mean = NA_real_,
    sd = NA_real_,
    median = NA_real_,
    q1 = NA_real_,
    q3 = NA_real_,
    iqr = NA_real_,
    min = NA_real_,
    max = NA_real_,
    lcl95 = NA_real_,
    ucl95 = NA_real_
  )
}

make_base_row <- function(section, dataset, measure, definition, unit = "not_applicable", stratum_type = "overall", stratum_value = "overall", risk_scale = "not_applicable") {
  tibble(
    section = section,
    dataset = factor(dataset, levels = dataset_order),
    measure = measure,
    definition = definition,
    unit = unit,
    stratum_type = stratum_type,
    stratum_value = stratum_value,
    risk_scale = risk_scale
  ) %>% bind_cols(empty_numeric_fields())
}

build_single_value_row <- function(section, dataset, measure, value, definition, unit = "not_applicable", stratum_type = "overall", stratum_value = "overall", risk_scale = "not_applicable") {
  out <- make_base_row(section, dataset, measure, definition, unit, stratum_type, stratum_value, risk_scale)
  if (identical(unit, "person_years")) {
    out$person_years <- as.numeric(value)
  } else {
    out$count <- as.numeric(value)
  }
  out
}

build_count_row <- function(section, dataset, measure, count, denominator, definition, stratum_type = "overall", stratum_value = "overall", risk_scale = "not_applicable") {
  out <- make_base_row(section, dataset, measure, definition, unit = "count", stratum_type, stratum_value, risk_scale)
  out$n_total <- as.integer(denominator)
  out$n_nonmissing <- as.integer(denominator)
  out$count <- as.numeric(count)
  out$denominator <- as.numeric(denominator)
  out$proportion <- ifelse(is.na(denominator) || denominator <= 0, NA_real_, as.numeric(count) / as.numeric(denominator))
  out$percent <- 100 * out$proportion
  out
}

build_numeric_summary_row <- function(section, dataset, measure, x, definition, unit, stratum_type = "overall", stratum_value = "overall", risk_scale = "not_applicable") {
  x <- x[!is.na(x)]
  out <- make_base_row(section, dataset, measure, definition, unit, stratum_type, stratum_value, risk_scale)
  out$n_total <- length(x)
  out$n_nonmissing <- length(x)
  
  if (length(x) == 0L) {
    return(out)
  }
  
  out$mean <- mean(x)
  out$sd <- if (length(x) >= 2L) stats::sd(x) else 0
  out$median <- stats::median(x)
  out$q1 <- unname(stats::quantile(x, probs = 0.25, names = FALSE, type = 2))
  out$q3 <- unname(stats::quantile(x, probs = 0.75, names = FALSE, type = 2))
  out$iqr <- out$q3 - out$q1
  out$min <- min(x)
  out$max <- max(x)
  out
}

scale_numeric_summary_row <- function(row_df, new_measure, new_unit, multiplier) {
  scale_cols <- c("mean", "sd", "median", "q1", "q3", "iqr", "min", "max", "lcl95", "ucl95")
  out <- row_df
  out$measure <- new_measure
  out$unit <- new_unit
  for (col_name in intersect(scale_cols, names(out))) {
    if (is.numeric(out[[col_name]])) {
      out[[col_name]] <- out[[col_name]] * multiplier
    }
  }
  out
}

build_poisson_rate_row <- function(section, dataset, measure, count, person_years, definition, stratum_type = "overall", stratum_value = "overall", risk_scale = "not_applicable") {
  out <- make_base_row(section, dataset, measure, definition, unit = "per100_person_years", stratum_type, stratum_value, risk_scale)
  out$count <- as.numeric(count)
  out$person_years <- as.numeric(person_years)
  
  if (is.na(person_years) || person_years <= 0) {
    return(out)
  }
  
  pt <- stats::poisson.test(x = count, T = person_years)
  out$rate_per100py <- unname(pt$estimate) * 100
  out$rate_lcl95_per100py <- pt$conf.int[1] * 100
  out$rate_ucl95_per100py <- pt$conf.int[2] * 100
  out
}

build_reverse_km_row <- function(section, dataset, measure, time, censor_event, definition, unit, risk_scale) {
  out <- make_base_row(section, dataset, measure, definition, unit, stratum_type = "overall", stratum_value = "overall", risk_scale = risk_scale)
  
  time <- as.numeric(time)
  censor_event <- as.integer(censor_event)
  keep <- !is.na(time) & !is.na(censor_event)
  time <- time[keep]
  censor_event <- censor_event[keep]
  
  out$n_total <- length(time)
  out$n_nonmissing <- length(time)
  out$count <- sum(censor_event == 1L)
  out$denominator <- length(time)
  out$proportion <- ifelse(length(time) == 0L, NA_real_, sum(censor_event == 1L) / length(time))
  out$percent <- 100 * out$proportion
  
  if (length(time) == 0L || sum(censor_event == 1L) == 0L) {
    return(out)
  }
  
  fit <- survival::survfit(survival::Surv(time, censor_event) ~ 1)
  tab <- fit$table
  
  out$median <- unname(tab["median"])
  out$lcl95 <- unname(tab["0.95LCL"])
  out$ucl95 <- unname(tab["0.95UCL"])
  out$min <- min(time)
  out$max <- max(time)
  out
}

bind_section_rows <- function(...) {
  dplyr::bind_rows(...) %>%
    mutate(
      dataset = factor(as.character(dataset), levels = dataset_order),
      section = factor(section, levels = section_order)
    ) %>%
    arrange(section, dataset, stratum_type, stratum_value, measure)
}

## 🟠 Define: dataset-specific summary engines ===============================
make_cohort_overview_rows <- function(df, dataset_name) {
  dplyr::bind_rows(
    build_single_value_row(
      section = "cohort_overview",
      dataset = dataset_name,
      measure = "n_rows",
      value = nrow(df),
      definition = "Number of rows/subjects in the analysis dataset.",
      unit = "subjects"
    ),
    build_single_value_row(
      section = "cohort_overview",
      dataset = dataset_name,
      measure = "n_unique_site_id",
      value = dplyr::n_distinct(df$site_id),
      definition = "Number of unique site+id identifiers.",
      unit = "subjects"
    ),
    build_single_value_row(
      section = "cohort_overview",
      dataset = dataset_name,
      measure = "n_duplicate_site_id_rows",
      value = nrow(df) - dplyr::n_distinct(df$site_id),
      definition = "Number of rows exceeding the unique site+id count; zero is expected for one-subject-per-row data.",
      unit = "rows"
    ),
    build_single_value_row(
      section = "cohort_overview",
      dataset = dataset_name,
      measure = "total_person_years_observed",
      value = sum(df$years_followup, na.rm = TRUE),
      definition = "Total observed person-time computed as days_followup / 365.25 across all subjects.",
      unit = "person_years"
    )
  )
}

make_sex_distribution_rows <- function(df, dataset_name) {
  total_n <- nrow(df)
  sex_counts <- df %>% count(sex_label, name = "count")
  
  purrr::pmap_dfr(
    list(sex_counts$sex_label, sex_counts$count),
    function(sex_value, count_value) {
      build_count_row(
        section = "sex_distribution",
        dataset = dataset_name,
        measure = paste0("sex_", str_to_lower(sex_value)),
        count = count_value,
        denominator = total_n,
        definition = "Sex distribution based on sex_num coding (0 = Female, 1 = Male).",
        stratum_type = "sex",
        stratum_value = sex_value
      )
    }
  )
}

make_baseline_age_rows <- function(df, dataset_name) {
  build_numeric_summary_row(
    section = "baseline_age_summary",
    dataset = dataset_name,
    measure = "age_exact_entry",
    x = df$age_exact_entry,
    definition = "Exact age at cohort entry (years) using the merged analysis backbone field age_exact_entry.",
    unit = "years"
  )
}

make_age_by_sex_rows <- function(df, dataset_name) {
  df %>%
    group_by(sex_label) %>%
    group_modify(~ build_numeric_summary_row(
      section = "age_by_sex",
      dataset = dataset_name,
      measure = "age_exact_entry_by_sex",
      x = .x$age_exact_entry,
      definition = "Exact age at cohort entry summarized within each sex stratum.",
      unit = "years",
      stratum_type = "sex",
      stratum_value = unique(.x$sex_label)
    )) %>%
    ungroup()
}

make_event_distribution_rows <- function(df, dataset_name) {
  total_n <- nrow(df)
  status_counts <- df %>% count(status_label, name = "count")
  
  purrr::pmap_dfr(
    list(status_counts$status_label, status_counts$count),
    function(status_value, count_value) {
      build_count_row(
        section = "event_distribution",
        dataset = dataset_name,
        measure = paste0("status_", status_value),
        count = count_value,
        denominator = total_n,
        definition = "Observed outcome distribution using status_num coding: 0 = right_censoring, 1 = transition, 2 = remission.",
        stratum_type = "status",
        stratum_value = status_value
      )
    }
  )
}

make_event_rate_rows <- function(df, dataset_name) {
  total_py <- sum(df$years_followup, na.rm = TRUE)
  dplyr::bind_rows(
    build_poisson_rate_row(
      section = "event_rate_summary",
      dataset = dataset_name,
      measure = "transition_incidence_rate",
      count = sum(df$is_transition, na.rm = TRUE),
      person_years = total_py,
      definition = "Observed transition incidence rate per 100 person-years based on total observed person-time.",
      stratum_type = "event_type",
      stratum_value = "transition",
      risk_scale = main_risk_scale
    ),
    build_poisson_rate_row(
      section = "event_rate_summary",
      dataset = dataset_name,
      measure = "remission_incidence_rate",
      count = sum(df$is_remission, na.rm = TRUE),
      person_years = total_py,
      definition = "Observed remission incidence rate per 100 person-years based on total observed person-time.",
      stratum_type = "event_type",
      stratum_value = "remission",
      risk_scale = supplementary_risk_scale
    )
  )
}

make_followup_numeric_rows <- function(df, dataset_name) {
  observed_days <- build_numeric_summary_row(
    section = "followup_numeric_summary",
    dataset = dataset_name,
    measure = "followup_observed_all_days",
    x = df$days_followup,
    definition = "Observed analysis time T = days_followup for all subjects, regardless of event or censoring status.",
    unit = "days"
  )
  
  non_transition_days <- build_numeric_summary_row(
    section = "followup_numeric_summary",
    dataset = dataset_name,
    measure = "followup_non_transition_subjects_days",
    x = df$days_followup[df$status_num != 1L],
    definition = "Observed follow-up among subjects without observed transition (status 0 or 2); useful on the main transition-only scale where these observations act as censoring for transition.",
    unit = "days",
    risk_scale = main_risk_scale
  )
  
  right_censor_days <- build_numeric_summary_row(
    section = "followup_numeric_summary",
    dataset = dataset_name,
    measure = "followup_right_censored_only_days",
    x = df$days_followup[df$status_num == 0L],
    definition = "Observed follow-up among right-censored subjects only (status 0).",
    unit = "days"
  )
  
  transition_days <- build_numeric_summary_row(
    section = "followup_numeric_summary",
    dataset = dataset_name,
    measure = "time_to_transition_days",
    x = df$days_followup[df$status_num == 1L],
    definition = "Observed time to transition among transitioned subjects only.",
    unit = "days",
    risk_scale = main_risk_scale
  )
  
  remission_days <- build_numeric_summary_row(
    section = "followup_numeric_summary",
    dataset = dataset_name,
    measure = "time_to_remission_days",
    x = df$days_followup[df$status_num == 2L],
    definition = "Observed time to remission among remission subjects only.",
    unit = "days",
    risk_scale = supplementary_risk_scale
  )
  
  rows_days <- bind_rows(
    observed_days,
    non_transition_days,
    right_censor_days,
    transition_days,
    remission_days
  )
  
  rows_years <- bind_rows(
    scale_numeric_summary_row(observed_days, "followup_observed_all_years", "years", 1 / days_per_year),
    scale_numeric_summary_row(non_transition_days, "followup_non_transition_subjects_years", "years", 1 / days_per_year),
    scale_numeric_summary_row(right_censor_days, "followup_right_censored_only_years", "years", 1 / days_per_year),
    scale_numeric_summary_row(transition_days, "time_to_transition_years", "years", 1 / days_per_year),
    scale_numeric_summary_row(remission_days, "time_to_remission_years", "years", 1 / days_per_year)
  )
  
  bind_rows(rows_days, rows_years)
}

make_followup_reverse_km_rows <- function(df, dataset_name) {
  transition_scale_days <- build_reverse_km_row(
    section = "followup_reverse_km_summary",
    dataset = dataset_name,
    measure = "reverse_km_followup_transition_scale_days",
    time = df$days_followup,
    censor_event = as.integer(df$status_num != 1L),
    definition = "Reverse Kaplan-Meier follow-up for the main transition-only analysis, treating non-transition observations (status 0 or 2) as censoring events in the censoring-time distribution.",
    unit = "days",
    risk_scale = main_risk_scale
  )
  
  competing_scale_days <- build_reverse_km_row(
    section = "followup_reverse_km_summary",
    dataset = dataset_name,
    measure = "reverse_km_followup_competing_scale_days",
    time = df$days_followup,
    censor_event = as.integer(df$status_num == 0L),
    definition = "Reverse Kaplan-Meier follow-up for the remission-sensitive competing-risk scale, treating only right-censoring (status 0) as censoring events in the censoring-time distribution.",
    unit = "days",
    risk_scale = supplementary_risk_scale
  )
  
  bind_rows(
    transition_scale_days,
    scale_numeric_summary_row(transition_scale_days, "reverse_km_followup_transition_scale_years", "years", 1 / days_per_year),
    competing_scale_days,
    scale_numeric_summary_row(competing_scale_days, "reverse_km_followup_competing_scale_years", "years", 1 / days_per_year)
  )
}

make_status_by_sex_rows <- function(df, dataset_name) {
  df %>%
    count(sex_label, status_label, name = "count") %>%
    group_by(sex_label) %>%
    mutate(denominator = sum(count)) %>%
    ungroup() %>%
    mutate(
      proportion = count / denominator,
      percent = 100 * proportion,
      measure = paste0("status_within_", str_to_lower(sex_label), "_", status_label),
      definition = "Status composition within each sex stratum; denominator is the sex-specific total.",
      section = "status_by_sex",
      dataset = dataset_name,
      unit = "count",
      stratum_type = "sex",
      stratum_value = sex_label,
      risk_scale = "not_applicable"
    ) %>%
    transmute(
      section,
      dataset = factor(dataset, levels = dataset_order),
      measure,
      definition,
      unit,
      stratum_type,
      stratum_value,
      risk_scale,
      n_total = as.integer(denominator),
      n_nonmissing = as.integer(denominator),
      count = as.numeric(count),
      denominator = as.numeric(denominator),
      proportion,
      percent,
      person_years = NA_real_,
      rate_per100py = NA_real_,
      rate_lcl95_per100py = NA_real_,
      rate_ucl95_per100py = NA_real_,
      mean = NA_real_,
      sd = NA_real_,
      median = NA_real_,
      q1 = NA_real_,
      q3 = NA_real_,
      iqr = NA_real_,
      min = NA_real_,
      max = NA_real_,
      lcl95 = NA_real_,
      ucl95 = NA_real_
    )
}

make_calendar_rows <- function(df, dataset_name) {
  if (!"entry_year" %in% names(df) || all(is.na(df$entry_year))) {
    return(tibble())
  }
  
  entry_rows <- df %>%
    count(entry_year, name = "count") %>%
    mutate(
      section = "calendar_distribution",
      dataset = dataset_name,
      measure = "entry_year_distribution",
      definition = "Calendar distribution of cohort entry years.",
      unit = "count",
      stratum_type = "entry_year",
      stratum_value = as.character(entry_year),
      risk_scale = "not_applicable",
      denominator = nrow(df),
      proportion = count / denominator,
      percent = 100 * proportion
    )
  
  end_rows <- df %>%
    count(end_year, name = "count") %>%
    mutate(
      section = "calendar_distribution",
      dataset = dataset_name,
      measure = "followup_end_year_distribution",
      definition = "Calendar distribution of follow-up end years derived from date_entry + days_followup.",
      unit = "count",
      stratum_type = "end_year",
      stratum_value = as.character(end_year),
      risk_scale = "not_applicable",
      denominator = nrow(df),
      proportion = count / denominator,
      percent = 100 * proportion
    )
  
  bind_rows(entry_rows, end_rows) %>%
    transmute(
      section = factor(section, levels = section_order),
      dataset = factor(dataset, levels = dataset_order),
      measure,
      definition,
      unit,
      stratum_type,
      stratum_value,
      risk_scale,
      n_total = as.integer(denominator),
      n_nonmissing = as.integer(denominator),
      count = as.numeric(count),
      denominator = as.numeric(denominator),
      proportion,
      percent,
      person_years = NA_real_,
      rate_per100py = NA_real_,
      rate_lcl95_per100py = NA_real_,
      rate_ucl95_per100py = NA_real_,
      mean = NA_real_,
      sd = NA_real_,
      median = NA_real_,
      q1 = NA_real_,
      q3 = NA_real_,
      iqr = NA_real_,
      min = NA_real_,
      max = NA_real_,
      lcl95 = NA_real_,
      ucl95 = NA_real_
    )
}

make_dataset_summary_rows <- function(df, dataset_name) {
  bind_section_rows(
    make_cohort_overview_rows(df, dataset_name),
    make_sex_distribution_rows(df, dataset_name),
    make_baseline_age_rows(df, dataset_name),
    make_age_by_sex_rows(df, dataset_name),
    make_event_distribution_rows(df, dataset_name),
    make_event_rate_rows(df, dataset_name),
    make_followup_numeric_rows(df, dataset_name),
    make_followup_reverse_km_rows(df, dataset_name),
    make_status_by_sex_rows(df, dataset_name),
    make_calendar_rows(df, dataset_name)
  )
}

## 🟠 Define: workbook export helpers ===============================
write_summary_workbook <- function(summary_rows, workbook_path) {
  section_split <- summary_rows %>%
    mutate(section = as.character(section), dataset = as.character(dataset)) %>%
    group_split(section, .keep = TRUE)
  
  section_names <- summary_rows %>%
    mutate(section = as.character(section)) %>%
    distinct(section) %>%
    pull(section)
  
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "all_rows")
    openxlsx::writeData(wb, sheet = "all_rows", x = summary_rows)
    
    for (i in seq_along(section_split)) {
      sheet_name <- stringr::str_sub(section_names[i], 1, 31)
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet = sheet_name, x = section_split[[i]])
    }
    
    openxlsx::saveWorkbook(wb, workbook_path, overwrite = TRUE)
    return(invisible(TRUE))
  }
  
  if (requireNamespace("writexl", quietly = TRUE)) {
    x_list <- c(list(all_rows = summary_rows), base::setNames(section_split, nm = stringr::str_sub(section_names, 1, 31)))
    writexl::write_xlsx(x = x_list, path = workbook_path)
    return(invisible(TRUE))
  }
  
  warning("Neither `openxlsx` nor `writexl` is available; workbook export skipped.", call. = FALSE)
  invisible(FALSE)
}

# 🔴 Load: merged descriptive input ===============================
raw_merged <- read_input_dataset(merged_file)

# 🔴 Validate: backbone fields and derived descriptors ===============================
analysis_df <- prepare_backbone_dataset(raw_merged, pnu_label = pnu_site_label, snu_label = snu_site_label)
datasets <- split_dataset_variants(analysis_df, pnu_label = pnu_site_label, snu_label = snu_site_label)

# 🔴 Summarize: pre-analysis inspection rows ===============================
summary_rows <- purrr::imap_dfr(datasets, make_dataset_summary_rows) %>%
  mutate(
    dataset = factor(as.character(dataset), levels = dataset_order),
    section = factor(as.character(section), levels = section_order)
  ) %>%
  arrange(section, dataset, stratum_type, stratum_value, measure)

if (nrow(summary_rows) == 0) {
  stop("No descriptive summary rows were generated.", call. = FALSE)
}

# 🔴 Export: consolidated descriptive artifacts ===============================
summary_csv_file <- file.path(export_path, paste0(output_stem, ".csv"))
summary_xlsx_file <- file.path(export_path, paste0(output_stem, ".xlsx"))

readr::write_csv(summary_rows, summary_csv_file, na = "")
write_summary_workbook(summary_rows, summary_xlsx_file)

export_registry <- tibble(
  artifact = c("summary_csv", "summary_xlsx"),
  file_path = c(summary_csv_file, summary_xlsx_file),
  exists = c(file.exists(summary_csv_file), file.exists(summary_xlsx_file))
)

print(export_registry)

cat("\n[Completed] Pre-analysis clinical survival data-inspection export finished.\n")
cat("- CSV source of truth: ", summary_csv_file, "\n", sep = "")
if (file.exists(summary_xlsx_file)) {
  cat("- Workbook convenience file: ", summary_xlsx_file, "\n", sep = "")
} else {
  cat("- Workbook convenience file: skipped because no xlsx writer package was available.\n", sep = "")
}
