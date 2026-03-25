# 🔴 Configure: paths and fixed backbone ===============================
merged_file <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
export_path <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.데이터 분석 New/🟦Stage1/attachments"

data_path <- dirname(merged_file)

pnu_site_label <- "PNU"
snu_site_label <- "SNU"

common_horizons_year <- 1:10
risk_thresholds <- c(0.05, 0.10, 0.15, 0.20)  # user-editable

# 🔴 Initialize: packages and options ===============================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

options(stringsAsFactors = FALSE, scipen = 999)

dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

dataset_order <- c("PNU", "SNU", "merged")

common_horizons_year <- as.integer(common_horizons_year)
risk_thresholds <- sort(unique(as.numeric(risk_thresholds)))

if (!identical(common_horizons_year, 1:10)) {
  stop("Stage 1 requires `common_horizons_year <- 1:10`.", call. = FALSE)
}

if (length(risk_thresholds) == 0 || anyNA(risk_thresholds) || any(risk_thresholds <= 0 | risk_thresholds >= 1)) {
  stop("`risk_thresholds` must be unique non-missing probabilities strictly between 0 and 1.", call. = FALSE)
}

if (toupper(pnu_site_label) == toupper(snu_site_label)) {
  stop("`pnu_site_label` and `snu_site_label` must be different.", call. = FALSE)
}

horizon_string <- paste(common_horizons_year, collapse = ",")
threshold_string <- paste(format(risk_thresholds, trim = TRUE, scientific = FALSE), collapse = ",")

# 🔴 Define: backbone helpers ===============================
## 🟠 Define: input readers ===============================
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

normalize_existing_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

## 🟠 Define: site harmonizers ===============================
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

split_site_dataset <- function(merged_df, site_label, dataset_name) {
  out <- merged_df %>%
    filter(toupper(trimws(as.character(site))) == toupper(site_label)) %>%
    mutate(site = site_label)
  
  if (nrow(out) == 0) {
    stop(
      sprintf("[%s] No rows found in merged input for site label `%s`.", dataset_name, site_label),
      call. = FALSE
    )
  }
  
  out
}

## 🟠 Define: backbone validators ===============================
prepare_backbone_dataset <- function(df, dataset_name, site_mode = c("single", "merged")) {
  site_mode <- match.arg(site_mode)
  
  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(
      sprintf("[%s] Missing required columns: %s", dataset_name, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  
  df <- df %>%
    mutate(
      id = trimws(as.character(id)),
      site = trimws(as.character(site)),
      sex_num = as.integer(coerce_numeric_text(sex_num)),
      age_exact_entry = coerce_numeric_text(age_exact_entry),
      days_followup = coerce_numeric_text(days_followup),
      status_num = as.integer(coerce_numeric_text(status_num))
    )
  
  if (nrow(df) == 0) {
    stop(sprintf("[%s] Dataset has zero rows after loading.", dataset_name), call. = FALSE)
  }
  
  if (anyNA(df[required_cols])) {
    stop(sprintf("[%s] Missing values detected in required backbone columns.", dataset_name), call. = FALSE)
  }
  
  if (any(df$id == "", na.rm = TRUE)) {
    stop(sprintf("[%s] Blank `id` values detected.", dataset_name), call. = FALSE)
  }
  
  if (any(df$site == "", na.rm = TRUE)) {
    stop(sprintf("[%s] Blank `site` values detected.", dataset_name), call. = FALSE)
  }
  
  if (any(!df$sex_num %in% c(0L, 1L))) {
    stop(sprintf("[%s] `sex_num` must be coded as 0/1 only.", dataset_name), call. = FALSE)
  }
  
  if (any(!df$status_num %in% c(0L, 1L, 2L))) {
    stop(sprintf("[%s] `status_num` must be coded as 0/1/2 only.", dataset_name), call. = FALSE)
  }
  
  if (any(df$days_followup < 0)) {
    stop(sprintf("[%s] Negative `days_followup` values detected.", dataset_name), call. = FALSE)
  }
  
  n_site_levels <- dplyr::n_distinct(df$site)
  
  if (site_mode == "single" && n_site_levels != 1L) {
    stop(sprintf("[%s] Single-cohort input must contain exactly one site.", dataset_name), call. = FALSE)
  }
  
  if (site_mode == "merged" && n_site_levels < 2L) {
    stop(sprintf("[%s] Merged input must contain at least two site levels.", dataset_name), call. = FALSE)
  }
  
  age_mean <- mean(df$age_exact_entry)
  age_sd <- stats::sd(df$age_exact_entry)
  
  if (is.na(age_sd) || age_sd <= 0) {
    stop(sprintf("[%s] `age_exact_entry` must have positive standard deviation.", dataset_name), call. = FALSE)
  }
  
  df <- df %>%
    mutate(
      unique_person_id = paste(site, id, sep = "_"),
      time_year = days_followup / 365.25,
      event_main = as.integer(status_num == 1L),
      right_censor_flag = as.integer(status_num == 0L),
      remission_flag = as.integer(status_num == 2L),
      censor_main = as.integer(status_num %in% c(0L, 2L)),
      sex_label = factor(if_else(sex_num == 0L, "Female", "Male"), levels = c("Female", "Male")),
      status_label = factor(
        case_when(
          status_num == 0L ~ "right_censoring",
          status_num == 2L ~ "remission",
          TRUE ~ "transition"
        ),
        levels = c("right_censoring", "remission", "transition")
      ),
      age_s = (age_exact_entry - age_mean) / (2 * age_sd)
    )
  
  if (nrow(df) != dplyr::n_distinct(df$unique_person_id)) {
    stop(sprintf("[%s] `site + id` is not unique within dataset.", dataset_name), call. = FALSE)
  }
  
  scaling_row <- tibble::tibble(
    dataset = dataset_name,
    variable = "age_exact_entry",
    scaled_variable = "age_s",
    center_mean = age_mean,
    scale_sd = age_sd,
    scale_two_sd = 2 * age_sd,
    scaling_rule = "(age_exact_entry - mean(age_exact_entry)) / (2 * sd(age_exact_entry))"
  )
  
  list(data = df, scaling = scaling_row)
}

## 🟠 Define: registry builders ===============================
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
      dataset, formula_id, formula_name, formula_label,
      formula_rhs, formula_full,
      uses_site, uses_age_sex_interaction,
      site_branch, interaction_branch
    )
}

make_dataset_registry <- function(analysis_datasets, source_lookup) {
  bind_rows(lapply(names(analysis_datasets), function(dataset_name) {
    df <- analysis_datasets[[dataset_name]]
    source_desc <- unname(source_lookup[[dataset_name]])
    
    tibble::tibble(
      dataset = dataset_name,
      source_description = source_desc,
      source_type = ifelse(grepl("^derived_from", source_desc), "derived", "file_input"),
      n_rows = nrow(df),
      n_unique_person_id = dplyr::n_distinct(df$unique_person_id),
      n_site_levels = dplyr::n_distinct(df$site),
      site_values = paste(sort(unique(df$site)), collapse = "|"),
      unique_person_id_variable = "unique_person_id",
      unique_person_id_rule = "site + id",
      time_variable = "days_followup",
      reporting_time_variable = "time_year",
      event_definition = "status_num == 1",
      censoring_definition = "status_num %in% c(0, 2)",
      age_variable = "age_exact_entry",
      scaled_age_variable = "age_s",
      sex_variable = "sex_num",
      keep_site_free_branch = TRUE,
      keep_site_adjusted_branch = dataset_name == "merged",
      keep_age_sex_interaction_branch = TRUE
    )
  }))
}

make_dataset_summary <- function(analysis_datasets) {
  bind_rows(lapply(names(analysis_datasets), function(dataset_name) {
    df <- analysis_datasets[[dataset_name]]
    
    tibble::tibble(
      dataset = dataset_name,
      n_rows = nrow(df),
      n_unique_person_id = dplyr::n_distinct(df$unique_person_id),
      n_status1_transition = sum(df$event_main),
      n_status2_remission = sum(df$remission_flag),
      n_status0_right_censoring = sum(df$right_censor_flag),
      n_main_censoring_0_or_2 = sum(df$censor_main),
      person_time_days = sum(df$days_followup),
      person_time_years = sum(df$time_year),
      median_followup_days = stats::median(df$days_followup),
      median_followup_years = stats::median(df$time_year),
      max_followup_days = max(df$days_followup),
      max_followup_years = max(df$time_year)
    )
  }))
}

make_age_summary_block <- function(df, dataset_name) {
  overall <- tibble::tibble(
    dataset = dataset_name,
    summary_domain = "age_overall",
    group_value = "overall",
    n = nrow(df),
    proportion = 1,
    age_mean = mean(df$age_exact_entry),
    age_sd = stats::sd(df$age_exact_entry),
    age_median = stats::median(df$age_exact_entry),
    age_q1 = stats::quantile(df$age_exact_entry, probs = 0.25, names = FALSE),
    age_q3 = stats::quantile(df$age_exact_entry, probs = 0.75, names = FALSE),
    age_min = min(df$age_exact_entry),
    age_max = max(df$age_exact_entry)
  )
  
  by_sex <- df %>%
    group_by(sex_label) %>%
    summarise(
      n = n(),
      age_mean = mean(age_exact_entry),
      age_sd = stats::sd(age_exact_entry),
      age_median = stats::median(age_exact_entry),
      age_q1 = stats::quantile(age_exact_entry, probs = 0.25, names = FALSE),
      age_q3 = stats::quantile(age_exact_entry, probs = 0.75, names = FALSE),
      age_min = min(age_exact_entry),
      age_max = max(age_exact_entry),
      .groups = "drop"
    ) %>%
    mutate(
      dataset = dataset_name,
      summary_domain = "age_by_sex",
      group_value = as.character(sex_label),
      proportion = n / sum(n)
    ) %>%
    select(dataset, summary_domain, group_value, n, proportion, age_mean, age_sd, age_median, age_q1, age_q3, age_min, age_max)
  
  by_site <- df %>%
    group_by(site) %>%
    summarise(
      n = n(),
      age_mean = mean(age_exact_entry),
      age_sd = stats::sd(age_exact_entry),
      age_median = stats::median(age_exact_entry),
      age_q1 = stats::quantile(age_exact_entry, probs = 0.25, names = FALSE),
      age_q3 = stats::quantile(age_exact_entry, probs = 0.75, names = FALSE),
      age_min = min(age_exact_entry),
      age_max = max(age_exact_entry),
      .groups = "drop"
    ) %>%
    mutate(
      dataset = dataset_name,
      summary_domain = "age_by_site",
      group_value = as.character(site),
      proportion = n / sum(n)
    ) %>%
    select(dataset, summary_domain, group_value, n, proportion, age_mean, age_sd, age_median, age_q1, age_q3, age_min, age_max)
  
  bind_rows(overall, by_sex, by_site)
}

make_distribution_block <- function(df, dataset_name) {
  sex_block <- df %>%
    group_by(sex_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      dataset = dataset_name,
      summary_domain = "sex_distribution",
      group_value = as.character(sex_label),
      proportion = n / sum(n),
      age_mean = NA_real_,
      age_sd = NA_real_,
      age_median = NA_real_,
      age_q1 = NA_real_,
      age_q3 = NA_real_,
      age_min = NA_real_,
      age_max = NA_real_
    ) %>%
    select(dataset, summary_domain, group_value, n, proportion, age_mean, age_sd, age_median, age_q1, age_q3, age_min, age_max)
  
  site_block <- df %>%
    group_by(site) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      dataset = dataset_name,
      summary_domain = "site_distribution",
      group_value = as.character(site),
      proportion = n / sum(n),
      age_mean = NA_real_,
      age_sd = NA_real_,
      age_median = NA_real_,
      age_q1 = NA_real_,
      age_q3 = NA_real_,
      age_min = NA_real_,
      age_max = NA_real_
    ) %>%
    select(dataset, summary_domain, group_value, n, proportion, age_mean, age_sd, age_median, age_q1, age_q3, age_min, age_max)
  
  bind_rows(sex_block, site_block)
}

make_age_sex_site_summary <- function(analysis_datasets) {
  bind_rows(lapply(names(analysis_datasets), function(dataset_name) {
    df <- analysis_datasets[[dataset_name]]
    
    bind_rows(
      make_age_summary_block(df, dataset_name),
      make_distribution_block(df, dataset_name)
    )
  })) %>%
    arrange(match(dataset, dataset_order), summary_domain, group_value)
}

make_horizon_registry <- function() {
  tibble::as_tibble(expand.grid(
    dataset = dataset_order,
    horizon_year = common_horizons_year,
    stringsAsFactors = FALSE
  )) %>%
    mutate(
      horizon_id = paste0("year_", horizon_year),
      horizon_days = horizon_year * 365.25,
      interpretation_tier = case_when(
        dataset == "PNU" & horizon_year == 1 ~ "primary-supported",
        dataset == "PNU" & horizon_year == 2 ~ "sensitivity",
        dataset == "PNU" & horizon_year >= 3 ~ "projection",
        dataset %in% c("SNU", "merged") & horizon_year <= 2 ~ "primary-supported",
        dataset %in% c("SNU", "merged") & horizon_year <= 5 ~ "secondary",
        TRUE ~ "projection"
      ),
      primary_supported_flag = interpretation_tier == "primary-supported",
      interpretation_note = case_when(
        dataset == "PNU" & horizon_year >= 2 ~ "Interpret cautiously; later stages must add empirical follow-up maturity and instability flags.",
        dataset %in% c("SNU", "merged") & horizon_year >= 6 ~ "Projection zone; later stages must report late-horizon instability explicitly.",
        TRUE ~ "Common comparison horizon retained for cross-model comparability."
      )
    ) %>%
    arrange(match(dataset, dataset_order), horizon_year)
}

make_threshold_registry <- function() {
  tibble::tibble(
    threshold_id = paste0("threshold_", seq_along(risk_thresholds)),
    threshold = risk_thresholds,
    threshold_label = paste0(format(100 * risk_thresholds, trim = TRUE, scientific = FALSE), "%"),
    positive_rule = "Classify as high risk when predicted risk >= threshold",
    applies_to = "all_datasets_all_model_classes_all_horizons",
    editable_in_script_config = TRUE
  )
}

make_metadata_registry <- function(source_lookup) {
  tibble::tibble(
    metadata_group = c(
      "stage", "stage", "stage",
      "documents", "documents",
      "identification",
      "time", "time",
      "event", "event",
      "covariates", "covariates",
      "horizons", "horizons",
      "thresholds", "thresholds",
      "interpretation", "interpretation",
      "outputs",
      "inputs", "inputs", "inputs", "inputs", "inputs"
    ),
    metadata_name = c(
      "stage_name", "stage_role", "model_fitting_allowed",
      "canonical_common_rules", "canonical_framework",
      "unique_person_id_rule",
      "analysis_time_variable", "reporting_time_variable",
      "event_definition", "main_censoring_definition",
      "age_variable", "age_scaling_rule",
      "horizon_vector", "horizon_rule",
      "threshold_vector", "threshold_edit_rule",
      "site_effect_interpretation_rule", "remission_reanalysis_note",
      "save_folder_rule",
      "data_path", "export_path", "merged_source", "pnu_source", "snu_source"
    ),
    metadata_value = c(
      "Stage 1 backbone lock",
      "Freeze one common comparison backbone only",
      "FALSE",
      "2.General Model Specifications_🇬🇧ENG.md",
      "5.Model Specification Framework_🇬🇧ENG.md",
      "site + id",
      "days_followup",
      "time_year = days_followup / 365.25",
      "status_num == 1",
      "status_num %in% c(0, 2)",
      "age_exact_entry",
      "(age_exact_entry - mean(age_exact_entry)) / (2 * sd(age_exact_entry))",
      horizon_string,
      "Fixed to 1:10 in Stage 1",
      threshold_string,
      "Edit `risk_thresholds` near the top of the script if needed",
      "Interpret site as a proxy for broader treatment context, care pathway, selection, or follow-up structure.",
      "Dedicated remission-sensitive reanalysis belongs to a later sensitivity stage, not Stage 1.",
      "Write all outputs into one `export_path` folder without creating extra subfolders.",
      data_path,
      export_path,
      unname(source_lookup[["merged"]]),
      unname(source_lookup[["PNU"]]),
      unname(source_lookup[["SNU"]])
    )
  )
}

make_modeling_registry <- function(formula_registry, dataset_registry) {
  formula_registry %>%
    left_join(dataset_registry %>% select(dataset, n_site_levels, site_values), by = "dataset") %>%
    mutate(
      stage = "Stage 1",
      stage_role = "backbone_only_no_model_fit",
      analysis_time_variable = "days_followup",
      reporting_time_variable = "time_year",
      event_definition = "status_num == 1",
      censoring_definition = "status_num %in% c(0, 2)",
      horizon_vector = horizon_string,
      threshold_vector = threshold_string,
      incidence_formula_rhs = formula_rhs,
      latency_formula_rhs = formula_rhs,
      intended_model_classes = "KM benchmark|non-cure|frequentist cure|Bayesian cure"
    ) %>%
    select(
      dataset, formula_id, formula_name, formula_label,
      formula_rhs, formula_full,
      incidence_formula_rhs, latency_formula_rhs,
      uses_site, uses_age_sex_interaction,
      site_branch, interaction_branch,
      n_site_levels, site_values,
      analysis_time_variable, reporting_time_variable,
      event_definition, censoring_definition,
      horizon_vector, threshold_vector,
      intended_model_classes, stage, stage_role
    )
}

# 🔴 Load: merged source and split cohorts ===============================
## 🟠 Read: merged input ===============================
merged_raw <- read_input_dataset(merged_file)
merged_raw <- standardize_known_site_labels(merged_raw, pnu_site_label, snu_site_label)

observed_site_values <- sort(unique(trimws(as.character(merged_raw$site))))
allowed_site_values <- c(pnu_site_label, snu_site_label)
unexpected_sites <- setdiff(observed_site_values, allowed_site_values)

if (length(unexpected_sites) > 0) {
  stop(
    sprintf(
      "Merged input contains unexpected site labels: %s",
      paste(unexpected_sites, collapse = ", ")
    ),
    call. = FALSE
  )
}

if (!all(allowed_site_values %in% observed_site_values)) {
  stop(
    sprintf(
      "Merged input must contain both `%s` and `%s` in `site`.",
      pnu_site_label, snu_site_label
    ),
    call. = FALSE
  )
}

## 🟠 Derive: PNU and SNU from merged ===============================
pnu_raw <- split_site_dataset(merged_raw, pnu_site_label, "PNU")
snu_raw <- split_site_dataset(merged_raw, snu_site_label, "SNU")

source_lookup <- c(
  PNU = paste0("derived_from_merged:", pnu_site_label),
  SNU = paste0("derived_from_merged:", snu_site_label),
  merged = normalize_existing_path(merged_file)
)

# 🔴 Construct: Stage-1 analysis datasets ===============================
## 🟠 Prepare: processed cohorts ===============================
pnu_prepped <- prepare_backbone_dataset(pnu_raw, "PNU", site_mode = "single")
snu_prepped <- prepare_backbone_dataset(snu_raw, "SNU", site_mode = "single")
merged_prepped <- prepare_backbone_dataset(merged_raw, "merged", site_mode = "merged")

analysis_datasets <- list(
  PNU = pnu_prepped$data,
  SNU = snu_prepped$data,
  merged = merged_prepped$data
)

scaling_registry <- bind_rows(
  pnu_prepped$scaling,
  snu_prepped$scaling,
  merged_prepped$scaling
) %>%
  arrange(match(dataset, dataset_order))

# 🔴 Assemble: registry and summary tables ===============================
## 🟠 Create: stage-one tables ===============================
dataset_registry <- make_dataset_registry(analysis_datasets, source_lookup) %>%
  arrange(match(dataset, dataset_order))

dataset_summary <- make_dataset_summary(analysis_datasets) %>%
  arrange(match(dataset, dataset_order))

age_sex_site_summary <- make_age_sex_site_summary(analysis_datasets)

formula_registry <- make_formula_registry()

horizon_registry <- make_horizon_registry()

threshold_registry <- make_threshold_registry()

metadata_registry <- make_metadata_registry(source_lookup)

modeling_registry <- make_modeling_registry(formula_registry, dataset_registry)

# 🔴 Export: reusable outputs ===============================
## 🟠 Save: CSV tables ===============================
readr::write_csv(dataset_registry, file.path(export_path, "stage1_dataset_registry.csv"))
readr::write_csv(dataset_summary, file.path(export_path, "stage1_dataset_summary.csv"))
readr::write_csv(age_sex_site_summary, file.path(export_path, "stage1_age_sex_site_summary.csv"))
readr::write_csv(scaling_registry, file.path(export_path, "stage1_scaling_registry.csv"))
readr::write_csv(formula_registry, file.path(export_path, "stage1_formula_registry.csv"))
readr::write_csv(modeling_registry, file.path(export_path, "stage1_modeling_registry.csv"))
readr::write_csv(horizon_registry, file.path(export_path, "stage1_horizon_registry.csv"))
readr::write_csv(threshold_registry, file.path(export_path, "stage1_threshold_registry.csv"))
readr::write_csv(metadata_registry, file.path(export_path, "stage1_metadata_registry.csv"))

## 🟠 Save: reusable RDS bundles ===============================
analysis_datasets_rds <- list(
  PNU = analysis_datasets[["PNU"]],
  SNU = analysis_datasets[["SNU"]],
  merged = analysis_datasets[["merged"]]
)

saveRDS(analysis_datasets_rds, file.path(export_path, "stage1_analysis_datasets.rds"))

backbone_bundle <- list(
  stage = "Stage 1",
  created_at = as.character(Sys.time()),
  session_info = utils::sessionInfo(),
  config = list(
    data_path = data_path,
    export_path = export_path,
    merged_file = merged_file,
    pnu_site_label = pnu_site_label,
    snu_site_label = snu_site_label,
    common_horizons_year = common_horizons_year,
    risk_thresholds = risk_thresholds
  ),
  source_lookup = as.list(source_lookup),
  datasets = analysis_datasets_rds,
  registries = list(
    dataset_registry = dataset_registry,
    dataset_summary = dataset_summary,
    age_sex_site_summary = age_sex_site_summary,
    scaling_registry = scaling_registry,
    formula_registry = formula_registry,
    modeling_registry = modeling_registry,
    horizon_registry = horizon_registry,
    threshold_registry = threshold_registry,
    metadata_registry = metadata_registry
  )
)

saveRDS(backbone_bundle, file.path(export_path, "stage1_backbone_bundle.rds"))

export_file_names <- c(
  "stage1_dataset_registry.csv",
  "stage1_dataset_summary.csv",
  "stage1_age_sex_site_summary.csv",
  "stage1_scaling_registry.csv",
  "stage1_formula_registry.csv",
  "stage1_modeling_registry.csv",
  "stage1_horizon_registry.csv",
  "stage1_threshold_registry.csv",
  "stage1_metadata_registry.csv",
  "stage1_analysis_datasets.rds",
  "stage1_backbone_bundle.rds",
  "stage1_export_manifest.csv"
)

export_object_names <- c(
  "dataset_registry",
  "dataset_summary",
  "age_sex_site_summary",
  "scaling_registry",
  "formula_registry",
  "modeling_registry",
  "horizon_registry",
  "threshold_registry",
  "metadata_registry",
  "analysis_datasets_rds",
  "backbone_bundle",
  "export_manifest"
)

export_descriptions <- c(
  "Dataset-level Stage 1 registry",
  "Dataset counts, event/remission/censoring counts, and person-time summary",
  "Age, sex, and site summary table",
  "Age scaling constants for age_s by dataset",
  "Common formula registry",
  "Common modeling registry table",
  "Dataset-by-horizon interpretation registry",
  "User-editable threshold registry",
  "Stage 1 backbone metadata",
  "Reusable processed datasets for later stages",
  "Reusable Stage 1 bundle with datasets and registries",
  "Manifest of all Stage 1 exported files"
)

export_manifest <- tibble::tibble(
  file_name = export_file_names,
  object_name = export_object_names,
  description = export_descriptions,
  file_path = file.path(export_path, export_file_names)
)

readr::write_csv(export_manifest, file.path(export_path, "stage1_export_manifest.csv"))