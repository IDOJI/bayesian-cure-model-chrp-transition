# 🔴 Configure: simplified transition-only project paths ===============================
sys_name <- Sys.info()[["sysname"]]

results_root_default <- switch(
  sys_name,
  "Darwin" = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  "Windows" = "C:/Users/clair/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  stop("Unsupported OS: ", sys_name)
)

results_root <- Sys.getenv("STAGE11_RESULTS_ROOT", unset = results_root_default)
export_path <- Sys.getenv(
  "STAGE11_EXPORT_PATH",
  unset = file.path(results_root, "stage11_Simplified transition-only cure project")
)

stage_paths <- list(
  stage3 = file.path(results_root, "stage3_KM benchmark classifier"),
  stage5 = file.path(results_root, "stage5_Individualized no-cure comparator"),
  stage7 = file.path(results_root, "stage7_Frequentist cure block"),
  stage8 = file.path(results_root, "stage8A_Bayesian transition-only cure")
)

stage11_spec_paths <- c(
  integrated_spec = normalizePath(
    "1.Model specifciation/Integrated_Modeling_Master_Specification_Transition_Only_Simplified.md",
    winslash = "/",
    mustWork = FALSE
  ),
  bayesian_spec = normalizePath(
    "1.Model specifciation/Bayesian_Modeling_Specification_Transition_Only_Simplified.md",
    winslash = "/",
    mustWork = FALSE
  )
)

required_packages <- c("readr", "dplyr", "tibble", "tidyr", "ggplot2", "scales")
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
  library(tidyr)
  library(ggplot2)
  library(scales)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# 🔴 Define: helpers ===============================
`%||%` <- function(x, y) if (is.null(x)) y else x

normalize_existing_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
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

make_temp_output_path <- function(path, tag = "tmp") {
  stamp <- paste0(format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%08d", sample.int(99999999, 1)))
  dir <- dirname(path)
  base <- basename(path)

  if (grepl("\\.csv\\.gz$", base, ignore.case = TRUE)) {
    base <- sub("\\.csv\\.gz$", paste0("_", tag, "_", stamp, ".csv.gz"), base, ignore.case = TRUE)
  } else if (grepl("\\.csv$", base, ignore.case = TRUE)) {
    base <- sub("\\.csv$", paste0("_", tag, "_", stamp, ".csv"), base, ignore.case = TRUE)
  } else if (grepl("\\.pdf$", base, ignore.case = TRUE)) {
    base <- sub("\\.pdf$", paste0("_", tag, "_", stamp, ".pdf"), base, ignore.case = TRUE)
  } else if (grepl("\\.md$", base, ignore.case = TRUE)) {
    base <- sub("\\.md$", paste0("_", tag, "_", stamp, ".md"), base, ignore.case = TRUE)
  } else {
    base <- paste0(base, "_", tag, "_", stamp)
  }

  file.path(dir, base)
}

safe_promote_file <- function(tmp_path, final_path) {
  dir.create(dirname(final_path), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(final_path)) {
    unlink(final_path)
  }

  ok <- file.rename(tmp_path, final_path)
  if (!isTRUE(ok)) {
    ok_copy <- file.copy(tmp_path, final_path, overwrite = TRUE)
    unlink(tmp_path)
    if (!isTRUE(ok_copy)) {
      stop("Failed to replace file: ", final_path, call. = FALSE)
    }
  }

  invisible(TRUE)
}

write_csv_preserve_schema <- function(df, path) {
  tmp <- make_temp_output_path(path, tag = "tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  readr::write_csv(df, tmp)
  safe_promote_file(tmp, path)
}

write_text_file <- function(lines, path) {
  tmp <- make_temp_output_path(path, tag = "tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writeLines(lines, tmp, useBytes = TRUE)
  safe_promote_file(tmp, path)
}

copy_file_atomic <- function(from, to) {
  assert_file_exists(from, basename(from))
  tmp <- make_temp_output_path(to, tag = "tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  ok <- file.copy(from, tmp, overwrite = TRUE)
  if (!isTRUE(ok)) {
    stop("Failed to copy file from `", from, "` to temporary path.", call. = FALSE)
  }
  safe_promote_file(tmp, to)
}

first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA)
  }
  x[[1L]]
}

dataset_version_registry <- tibble(
  dataset_version_key = c("PNU", "SNU", "merged", "merged_site_adjusted"),
  dataset_version_label = c("PNU", "SNU", "merged", "merged (site-adjusted)"),
  dataset_order = c(1L, 2L, 3L, 4L),
  site_control_status = c("not_applicable_single_site", "not_applicable_single_site", "site_unadjusted", "site_adjusted")
)

model_registry_display <- tibble(
  model_key = c(
    "kaplan_meier",
    "non_cure_lognormal",
    "frequentist_mixture_cure_lognormal",
    "bayesian_mixture_cure_lognormal"
  ),
  model_label = c(
    "Kaplan-Meier",
    "Non-cure log-normal",
    "Frequentist mixture cure (log-normal)",
    "Bayesian mixture cure (log-normal, anchor-informed)"
  ),
  model_order = c(1L, 2L, 3L, 4L)
)

model_palette <- c(
  "Kaplan-Meier" = "#1B4332",
  "Non-cure log-normal" = "#2A9D8F",
  "Frequentist mixture cure (log-normal)" = "#C96A23",
  "Bayesian mixture cure (log-normal, anchor-informed)" = "#1D4ED8"
)

projection_regions_from_table <- function(df) {
  df %>%
    filter(claim_restriction_flag %in% c("projection_only", "projection_plus_prior_sensitive")) %>%
    group_by(dataset_version_key, dataset_version_label) %>%
    summarise(
      projection_start = min(as.integer(horizon), na.rm = TRUE),
      .groups = "drop"
    )
}

augment_common_fields <- function(df) {
  df %>%
    left_join(dataset_version_registry, by = "dataset_version_key") %>%
    left_join(model_registry_display, by = "model_key") %>%
    mutate(
      projection_or_extrapolation_flag = claim_restriction_flag %in% c("projection_only", "projection_plus_prior_sensitive"),
      survival_percent = 100 * survival_probability,
      cumulative_risk_percent = 100 * cumulative_risk
    ) %>%
    arrange(dataset_order, horizon, model_order)
}

augment_susceptible_fields <- function(df) {
  df %>%
    left_join(dataset_version_registry, by = "dataset_version_key") %>%
    left_join(model_registry_display, by = "model_key") %>%
    mutate(
      projection_or_extrapolation_flag = claim_restriction_flag %in% c("projection_only", "projection_plus_prior_sensitive"),
      susceptible_survival_percent = 100 * susceptible_survival_probability,
      susceptible_cumulative_risk_percent = 100 * susceptible_cumulative_risk
    ) %>%
    arrange(dataset_order, horizon, model_order)
}

# 🔴 Read: upstream source tables ===============================
stage3_subject_predictions <- read_csv_checked(
  file.path(stage_paths$stage3, "stage3_km_subject_horizon_predictions.csv"),
  "Stage 3 KM subject prediction table"
)
stage5_subject_risk <- read_csv_checked(
  file.path(stage_paths$stage5, "stage5_subject_horizon_risk_long.csv"),
  "Stage 5 subject-by-horizon risk table"
)
stage7_risk_summary <- read_csv_checked(
  file.path(stage_paths$stage7, "stage7_risk_summary.csv"),
  "Stage 7 risk summary table"
)
stage8_performance <- read_csv_checked(
  file.path(stage_paths$stage8, "bayes_stage8a_performance_classification_long.csv"),
  "Stage 8A performance/classification table"
)
stage8_uncured <- read_csv_checked(
  file.path(stage_paths$stage8, "bayes_stage8a_uncured_only_decomposition_panel.csv"),
  "Stage 8A uncured decomposition table"
)
stage8_coefficient_summary <- read_csv_checked(
  file.path(stage_paths$stage8, "bayes_stage8a_coefficient_summary.csv"),
  "Stage 8A coefficient summary"
)
stage8_diagnostic_pdf <- file.path(stage_paths$stage8, "bayes_stage8a_diagnostic_plots.pdf")
assert_file_exists(stage8_diagnostic_pdf, "Stage 8A diagnostic PDF")

# 🔴 Build: Table 1 overall population comparison ===============================
km_overall <- stage3_subject_predictions %>%
  filter(
    (dataset_key == "PNU" & benchmark_id == "km_overall") |
      (dataset_key == "SNU" & benchmark_id == "km_overall") |
      (dataset_key == "merged__site_free" & benchmark_id == "km_overall") |
      (dataset_key == "merged__site_adjusted" & benchmark_id == "km_site_adjusted")
  ) %>%
  mutate(
    dataset_version_key = case_when(
      dataset_key == "PNU" ~ "PNU",
      dataset_key == "SNU" ~ "SNU",
      dataset_key == "merged__site_free" ~ "merged",
      dataset_key == "merged__site_adjusted" ~ "merged_site_adjusted",
      TRUE ~ NA_character_
    ),
    model_key = "kaplan_meier",
    source_stage = "Stage 3",
    source_dataset_key = dataset_key,
    source_model_id = benchmark_id,
    horizon = as.integer(horizon_year)
  ) %>%
  group_by(dataset_version_key, horizon) %>%
  summarise(
    model_key = first_non_missing(model_key),
    source_stage = first_non_missing(source_stage),
    source_dataset_key = first_non_missing(source_dataset_key),
    source_model_id = first_non_missing(source_model_id),
    survival_probability = mean(as.numeric(km_survival), na.rm = TRUE),
    cumulative_risk = mean(as.numeric(km_risk), na.rm = TRUE),
    support_tier = first_non_missing(as.character(support_tier)),
    horizon_evidence_class = first_non_missing(as.character(horizon_evidence_class)),
    claim_restriction_flag = first_non_missing(as.character(claim_restriction_flag)),
    interpretation_note = first_non_missing(as.character(interpretation_note)),
    n_subjects_averaged = dplyr::n(),
    .groups = "drop"
  )

km_selection_registry <- km_overall %>%
  distinct(dataset_version_key, model_key, source_stage, source_dataset_key, source_model_id) %>%
  mutate(selection_rule = "Stage 3 subject-level KM predictions averaged to overall population level; merged site-adjusted uses km_site_adjusted.")

noncure_overall <- stage5_subject_risk %>%
  filter(model_family == "lognormal") %>%
  mutate(
    dataset_version_key = case_when(
      dataset_key == "PNU" &
        formula_label == "Base" &
        site_branch == "site_free" &
        interaction_branch == "no_age_sex_interaction" ~ "PNU",
      dataset_key == "SNU" &
        formula_label == "Base" &
        site_branch == "site_free" &
        interaction_branch == "no_age_sex_interaction" ~ "SNU",
      dataset_key == "merged" &
        formula_label == "Base" &
        site_branch == "site_free" &
        interaction_branch == "no_age_sex_interaction" ~ "merged",
      dataset_key == "merged" &
        formula_label == "Site-added" &
        site_branch == "site_adjusted" &
        interaction_branch == "no_age_sex_interaction" ~ "merged_site_adjusted",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(dataset_version_key)) %>%
  group_by(dataset_version_key, horizon = as.integer(horizon_year)) %>%
  summarise(
    model_key = first_non_missing("non_cure_lognormal"),
    source_stage = first_non_missing("Stage 5"),
    source_dataset_key = first_non_missing(dataset_key),
    source_model_id = first_non_missing(model_id),
    survival_probability = mean(as.numeric(survival_pred), na.rm = TRUE),
    cumulative_risk = mean(as.numeric(risk_pred), na.rm = TRUE),
    support_tier = first_non_missing(as.character(support_tier)),
    horizon_evidence_class = first_non_missing(as.character(horizon_evidence_class)),
    claim_restriction_flag = first_non_missing(as.character(claim_restriction_flag)),
    interpretation_note = first_non_missing(as.character(interpretation_note)),
    n_subjects_averaged = dplyr::n(),
    .groups = "drop"
  )

noncure_selection_registry <- stage5_subject_risk %>%
  filter(model_family == "lognormal") %>%
  mutate(
    dataset_version_key = case_when(
      dataset_key == "PNU" & formula_label == "Base" & site_branch == "site_free" & interaction_branch == "no_age_sex_interaction" ~ "PNU",
      dataset_key == "SNU" & formula_label == "Base" & site_branch == "site_free" & interaction_branch == "no_age_sex_interaction" ~ "SNU",
      dataset_key == "merged" & formula_label == "Base" & site_branch == "site_free" & interaction_branch == "no_age_sex_interaction" ~ "merged",
      dataset_key == "merged" & formula_label == "Site-added" & site_branch == "site_adjusted" & interaction_branch == "no_age_sex_interaction" ~ "merged_site_adjusted",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(dataset_version_key)) %>%
  distinct(dataset_version_key, dataset_key, model_id, formula_label, site_branch, interaction_branch, formula_scope) %>%
  transmute(
    dataset_version_key,
    model_key = "non_cure_lognormal",
    source_stage = "Stage 5",
    source_dataset_key = dataset_key,
    source_model_id = model_id,
    selection_rule = paste0(
      "Stage 5 lognormal no-cure model using formula_label=",
      formula_label,
      ", site_branch=",
      site_branch,
      ", interaction_branch=",
      interaction_branch,
      "."
    )
  )

freq_cure_selected <- stage7_risk_summary %>%
  filter(model_family == "lnorm", model_class == "parametric_mixture_cure") %>%
  mutate(
    dataset_version_key = case_when(
      dataset_key == "PNU" &
        formula_variant == "base" &
        site_branch == "site_free" &
        interaction_branch == "no_age_sex_interaction" ~ "PNU",
      dataset_key == "SNU" &
        formula_variant == "base" &
        site_branch == "site_free" &
        interaction_branch == "no_age_sex_interaction" ~ "SNU",
      dataset_key == "merged__site_free" &
        formula_variant == "base_sitefree" &
        site_branch == "site_free" &
        site_placement_label == "site_in_neither" ~ "merged",
      dataset_key == "merged__site_adjusted" &
        formula_variant == "base_siteadjusted_both" &
        site_branch == "site_adjusted" ~ "merged_site_adjusted",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(dataset_version_key))

freq_cure_overall <- freq_cure_selected %>%
  transmute(
    dataset_version_key,
    model_key = "frequentist_mixture_cure_lognormal",
    source_stage = "Stage 7",
    source_dataset_key = dataset_key,
    source_model_id = model_id,
    horizon = as.integer(horizon_year),
    survival_probability = as.numeric(mean_survival_overall),
    cumulative_risk = as.numeric(mean_risk_overall),
    support_tier = as.character(support_tier),
    horizon_evidence_class = as.character(horizon_evidence_class),
    claim_restriction_flag = as.character(claim_restriction_flag),
    interpretation_note = as.character(interpretation_note)
  )

freq_cure_susceptible <- freq_cure_selected %>%
  transmute(
    dataset_version_key,
    model_key = "frequentist_mixture_cure_lognormal",
    source_stage = "Stage 7",
    source_dataset_key = dataset_key,
    source_model_id = model_id,
    horizon = as.integer(horizon_year),
    susceptible_survival_probability = as.numeric(mean_survival_susceptible),
    susceptible_cumulative_risk = as.numeric(mean_risk_susceptible),
    cure_fraction = as.numeric(mean_cure_fraction),
    susceptible_fraction = if_else(is.na(mean_cure_fraction), NA_real_, 1 - as.numeric(mean_cure_fraction)),
    support_tier = as.character(support_tier),
    horizon_evidence_class = as.character(horizon_evidence_class),
    claim_restriction_flag = as.character(claim_restriction_flag),
    interpretation_note = as.character(interpretation_note)
  )

freq_cure_selection_registry <- freq_cure_selected %>%
  distinct(dataset_version_key, dataset_key, model_id, formula_variant, site_branch, site_placement_label, incidence_rhs, latency_rhs) %>%
  transmute(
    dataset_version_key,
    model_key = "frequentist_mixture_cure_lognormal",
    source_stage = "Stage 7",
    source_dataset_key = dataset_key,
    source_model_id = model_id,
    selection_rule = paste0(
      "Stage 7 lognormal parametric mixture cure model using formula_variant=",
      formula_variant,
      ", site_placement=",
      site_placement_label,
      ", incidence_rhs=",
      incidence_rhs,
      ", latency_rhs=",
      latency_rhs,
      "."
    )
  )

stage8_cohort_only <- stage8_performance %>%
  filter(is.na(threshold), family_code == "LN", prior_branch == "anchor_informed", admissibility_flag %in% TRUE) %>%
  mutate(
    dataset_version_key = case_when(
      dataset_key == "PNU" &
        formula_anchor == "base" &
        site_placement == "site_in_neither" &
        site_prior_family == "not_applicable" ~ "PNU",
      dataset_key == "SNU" &
        formula_anchor == "base" &
        site_placement == "site_in_neither" &
        site_prior_family == "not_applicable" ~ "SNU",
      dataset_key == "merged" &
        formula_anchor == "base" &
        site_placement == "site_in_neither" &
        site_prior_family == "not_applicable" ~ "merged",
      dataset_key == "merged" &
        formula_anchor == "site_added" &
        site_placement == "site_in_both" &
        site_prior_family == "normal_0_1_main" ~ "merged_site_adjusted",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(dataset_version_key))

bayes_overall <- stage8_cohort_only %>%
  transmute(
    dataset_version_key,
    model_key = "bayesian_mixture_cure_lognormal",
    source_stage = "Stage 8A",
    source_dataset_key = dataset_key,
    source_model_id = model_id,
    horizon = as.integer(horizon),
    survival_probability = as.numeric(survival_mean),
    cumulative_risk = as.numeric(risk_mean),
    support_tier = as.character(support_tier),
    horizon_evidence_class = as.character(horizon_evidence_class),
    claim_restriction_flag = as.character(claim_restriction_flag),
    interpretation_note = as.character(interpretation_note)
  )

stage8_uncured_selected <- stage8_uncured %>%
  filter(family_code == "LN", prior_branch == "anchor_informed") %>%
  mutate(
    dataset_version_key = case_when(
      dataset_key == "PNU" &
        formula_anchor == "base" &
        site_placement == "site_in_neither" &
        site_prior_family == "not_applicable" ~ "PNU",
      dataset_key == "SNU" &
        formula_anchor == "base" &
        site_placement == "site_in_neither" &
        site_prior_family == "not_applicable" ~ "SNU",
      dataset_key == "merged" &
        formula_anchor == "base" &
        site_placement == "site_in_neither" &
        site_prior_family == "not_applicable" ~ "merged",
      dataset_key == "merged" &
        formula_anchor == "site_added" &
        site_placement == "site_in_both" &
        site_prior_family == "normal_0_1_main" ~ "merged_site_adjusted",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(dataset_version_key))

bayes_susceptible <- stage8_uncured_selected %>%
  transmute(
    dataset_version_key,
    model_key = "bayesian_mixture_cure_lognormal",
    source_stage = "Stage 8A",
    source_dataset_key = dataset_key,
    source_model_id = model_id,
    horizon = as.integer(horizon),
    susceptible_survival_probability = as.numeric(uncured_survival),
    susceptible_cumulative_risk = as.numeric(uncured_risk),
    cure_fraction = as.numeric(cure_fraction),
    susceptible_fraction = as.numeric(susceptible_fraction),
    support_tier = as.character(support_tier),
    horizon_evidence_class = as.character(horizon_evidence_class),
    claim_restriction_flag = as.character(claim_restriction_flag),
    interpretation_note = NA_character_
  )

bayes_selection_registry <- stage8_cohort_only %>%
  distinct(dataset_version_key, dataset_key, model_id, structural_model_id, formula_anchor, site_placement, site_prior_family) %>%
  transmute(
    dataset_version_key,
    model_key = "bayesian_mixture_cure_lognormal",
    source_stage = "Stage 8A",
    source_dataset_key = dataset_key,
    source_model_id = model_id,
    selection_rule = paste0(
      "Stage 8A anchor-informed LN model using formula_anchor=",
      formula_anchor,
      ", site_placement=",
      site_placement,
      ", site_prior_family=",
      site_prior_family,
      ", structural_model_id=",
      structural_model_id,
      "."
    )
  )

bayes_posterior_index <- stage8_cohort_only %>%
  distinct(
    dataset_version_key,
    dataset_key,
    model_id,
    structural_model_id,
    formula_anchor,
    site_placement,
    site_prior_family
  )

table1_overall <- bind_rows(
  km_overall %>% select(-n_subjects_averaged),
  noncure_overall %>% select(-n_subjects_averaged),
  freq_cure_overall,
  bayes_overall
) %>%
  augment_common_fields() %>%
  mutate(table_name = "Table 1. Overall survival and risk comparison")

table2_susceptible <- bind_rows(
  freq_cure_susceptible,
  bayes_susceptible
) %>%
  augment_susceptible_fields() %>%
  mutate(
    table_name = "Table 2. Susceptible-only quantities from mixture cure models",
    risk_scale = "transition_only_main",
    event_definition = "status_num == 1",
    remission_handling = "status_num == 2 treated as censoring"
  )

model_selection_registry <- bind_rows(
  km_selection_registry,
  noncure_selection_registry,
  freq_cure_selection_registry,
  bayes_selection_registry
) %>%
  left_join(dataset_version_registry, by = "dataset_version_key") %>%
  left_join(model_registry_display, by = "model_key") %>%
  arrange(dataset_order, model_order)

# 🔴 Build: Bayesian posterior summary and notes ===============================
bayes_posterior_summary <- stage8_coefficient_summary %>%
  filter(family_code == "LN", prior_branch == "anchor_informed") %>%
  inner_join(
    bayes_posterior_index,
    by = c("dataset_key", "model_id", "structural_model_id", "formula_anchor", "site_prior_family")
  ) %>%
  filter(!is.na(dataset_version_key)) %>%
  left_join(dataset_version_registry, by = "dataset_version_key") %>%
  select(
    dataset_version_key,
    dataset_version_label,
    model_id,
    structural_model_id,
    formula_anchor,
    site_placement,
    site_prior_family,
    parameter,
    mean,
    sd,
    q025,
    q50,
    q975
  ) %>%
  mutate(dataset_order = match(dataset_version_key, dataset_version_registry$dataset_version_key)) %>%
  arrange(dataset_order, model_id, parameter) %>%
  select(-dataset_order)

bayesian_note_lines <- c(
  "# Stage 11 Bayesian prior and diagnostics note",
  "",
  "This simplified project keeps the Stage 8A transition-only Bayesian mixture cure structure only.",
  "",
  "- Outcome rule: event = `status_num == 1`; remission (`status_num == 2`) is handled as censoring.",
  "- Latency family: log-normal only.",
  "- Main Bayesian branch: anchor-informed prior only.",
  "- Excluded from the simplified project: neutral/no-external-information branch, Stage 8B remission-aware branch, and all sensitivity analyses.",
  "- For merged site-adjusted reporting, the retained Bayesian model is the `site_added` / `site_in_both` / `normal_0_1_main` branch.",
  "- The external meta-analytic information is retained as an incidence-side prior/calibration anchor rather than as a remission model.",
  "- Posterior summaries are exported in `stage11_bayesian_posterior_summary.csv`.",
  "- The chain-level trace plots are carried forward by copying the original Stage 8A diagnostic PDF into this output folder as the Bayesian trace-plot reference file.",
  "",
  "Important interpretation note:",
  "This Stage 11 script is a lightweight reporting layer. It does not refit the Bayesian model; it reuses the locked Stage 8A outputs."
)

analysis_assumption_lines <- c(
  "# Stage 11 simplified project assumptions",
  "",
  "- This sub-project reuses upstream Stage 3, Stage 5, Stage 7, and Stage 8A outputs rather than refitting them.",
  "- The simplified comparison keeps only four models: KM, non-cure log-normal, frequentist mixture cure log-normal, and Bayesian mixture cure log-normal.",
  "- The merged site-adjusted cure model is operationally defined as the branch that includes site in both incidence and latency submodels.",
  "- The Bayesian main analysis retains the `anchor_informed` branch and the `normal_0_1_main` site prior family for the merged site-adjusted model.",
  "- Susceptible-only quantities are exported only for cure models because KM and non-cure models do not define a cure-fraction decomposition.",
  "- Horizons marked `projection_only` or `projection_plus_prior_sensitive` should be interpreted as late-horizon extrapolation/projection."
)

reuse_report_intro_lines <- c(
  "# Stage 11 reuse-of-existing-results report",
  "",
  "## Summary",
  "",
  "- Stage 11 did not refit the retained Stage 3, Stage 5, Stage 7, or Stage 8A models.",
  "- Instead, Stage 11 reused locked upstream result files and reorganized them into the simplified transition-only cure project outputs.",
  "- All Stage 11 outputs were written under the Stage 11 Dropbox export folder."
)

source_reuse_lines <- c(
  "",
  "## Reused source files and where they were used",
  "",
  paste0(
    "- `",
    normalize_existing_path(file.path(stage_paths$stage3, "stage3_km_subject_horizon_predictions.csv")),
    "` -> used to build the Kaplan-Meier rows in Table 1 and the KM series in Plot Set 1. For `merged (site-adjusted)`, `km_site_adjusted` subject-level predictions were averaged to the overall population level."
  ),
  paste0(
    "- `",
    normalize_existing_path(file.path(stage_paths$stage5, "stage5_subject_horizon_risk_long.csv")),
    "` -> used to build the non-cure log-normal rows in Table 1 and the non-cure log-normal series in Plot Set 1."
  ),
  paste0(
    "- `",
    normalize_existing_path(file.path(stage_paths$stage7, "stage7_risk_summary.csv")),
    "` -> used to build the frequentist mixture cure log-normal rows in Table 1, the frequentist susceptible-only rows in Table 2, the frequentist series in Plot Set 1, and the frequentist susceptible-only series in Plot Set 2."
  ),
  paste0(
    "- `",
    normalize_existing_path(file.path(stage_paths$stage8, "bayes_stage8a_performance_classification_long.csv")),
    "` -> used to build the Bayesian mixture cure log-normal overall rows in Table 1 and the Bayesian overall series in Plot Set 1."
  ),
  paste0(
    "- `",
    normalize_existing_path(file.path(stage_paths$stage8, "bayes_stage8a_uncured_only_decomposition_panel.csv")),
    "` -> used to build the Bayesian susceptible-only rows in Table 2 and the Bayesian susceptible-only series in Plot Set 2."
  ),
  paste0(
    "- `",
    normalize_existing_path(file.path(stage_paths$stage8, "bayes_stage8a_coefficient_summary.csv")),
    "` -> used to build `stage11_bayesian_posterior_summary.csv`."
  ),
  paste0(
    "- `",
    normalize_existing_path(stage8_diagnostic_pdf),
    "` -> copied forward as `stage11_bayesian_trace_plots_reference_from_stage8a.pdf` to document the Stage 8A trace-plot/diagnostic reference used by the simplified project."
  )
)

selection_report_lines <- c(
  "",
  "## Retained upstream model/result selections",
  ""
)

for (ii in seq_len(nrow(model_selection_registry))) {
  row_now <- model_selection_registry[ii, ]
  selection_report_lines <- c(
    selection_report_lines,
    paste0(
      "- `",
      row_now$dataset_version_label[[1L]],
      "` / `",
      row_now$model_label[[1L]],
      "` <- source stage `",
      row_now$source_stage[[1L]],
      "` / dataset key `",
      row_now$source_dataset_key[[1L]],
      "` / source model id `",
      row_now$source_model_id[[1L]],
      "`. ",
      row_now$selection_rule[[1L]]
    )
  )
}

output_report_lines <- c(
  "",
  "## Stage 11 outputs produced from those reused results",
  "",
  "- `stage11_table1_overall_survival_and_risk_comparison.csv`",
  "- `stage11_table2_susceptible_only_comparison.csv`",
  "- `stage11_plot_set1_overall_survival_and_risk_curves.pdf`",
  "- `stage11_plot_set2_susceptible_only_survival_and_risk_curves.pdf`",
  "- `stage11_bayesian_posterior_summary.csv`",
  "- `stage11_bayesian_trace_plots_reference_from_stage8a.pdf`",
  "- `stage11_model_selection_registry.csv`",
  "- `stage11_export_manifest.csv`"
)

reuse_report_lines <- c(
  reuse_report_intro_lines,
  source_reuse_lines,
  selection_report_lines,
  output_report_lines
)

# 🔴 Draw: plot sets ===============================
projection_regions <- projection_regions_from_table(table1_overall)

plot_overall_curves <- function(df, value_col, y_label, title_text, subtitle_text) {
  ggplot(df, aes(x = horizon, y = .data[[value_col]], color = model_label, group = model_label)) +
    geom_rect(
      data = projection_regions,
      aes(
        xmin = projection_start - 0.5,
        xmax = 10.5,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "grey92",
      alpha = 0.8
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    facet_wrap(~ dataset_version_label, ncol = 2) +
    scale_color_manual(values = model_palette, breaks = model_registry_display$model_label) +
    scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Years since cohort entry",
      y = y_label,
      color = "Model"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

plot_susceptible_curves <- function(df, value_col, y_label, title_text, subtitle_text) {
  ggplot(df, aes(x = horizon, y = .data[[value_col]], color = model_label, group = model_label)) +
    geom_rect(
      data = projection_regions,
      aes(
        xmin = projection_start - 0.5,
        xmax = 10.5,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "grey92",
      alpha = 0.8
    ) +
    geom_line(linewidth = 0.95) +
    geom_point(size = 1.6) +
    facet_wrap(~ dataset_version_label, ncol = 2) +
    scale_color_manual(
      values = model_palette[names(model_palette) %in% unique(df$model_label)],
      breaks = unique(df$model_label)
    ) +
    scale_x_continuous(breaks = 1:10, limits = c(1, 10)) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "Years since cohort entry",
      y = y_label,
      color = "Model"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

plot_set1_pdf <- file.path(export_path, "stage11_plot_set1_overall_survival_and_risk_curves.pdf")
plot_set2_pdf <- file.path(export_path, "stage11_plot_set2_susceptible_only_survival_and_risk_curves.pdf")

plot1_survival <- plot_overall_curves(
  table1_overall,
  value_col = "survival_probability",
  y_label = "Overall survival probability",
  title_text = "Stage 11 overall survival curves",
  subtitle_text = "Shaded region marks projection-dominant late horizons."
)

plot1_risk <- plot_overall_curves(
  table1_overall,
  value_col = "cumulative_risk",
  y_label = "Overall cumulative risk",
  title_text = "Stage 11 overall cumulative risk curves",
  subtitle_text = "Shaded region marks projection-dominant late horizons."
)

plot2_survival <- plot_susceptible_curves(
  table2_susceptible,
  value_col = "susceptible_survival_probability",
  y_label = "Susceptible-only survival probability",
  title_text = "Stage 11 susceptible-only survival curves",
  subtitle_text = "Only cure models are shown because susceptible-only quantities are not defined for KM or non-cure models."
)

plot2_risk <- plot_susceptible_curves(
  table2_susceptible,
  value_col = "susceptible_cumulative_risk",
  y_label = "Susceptible-only cumulative risk",
  title_text = "Stage 11 susceptible-only cumulative risk curves",
  subtitle_text = "Only cure models are shown because susceptible-only quantities are not defined for KM or non-cure models."
)

grDevices::pdf(plot_set1_pdf, width = 12, height = 8, onefile = TRUE)
print(plot1_survival)
print(plot1_risk)
grDevices::dev.off()

grDevices::pdf(plot_set2_pdf, width = 12, height = 8, onefile = TRUE)
print(plot2_survival)
print(plot2_risk)
grDevices::dev.off()

# 🔴 Export: tables, notes, and trace reference ===============================
table1_path <- file.path(export_path, "stage11_table1_overall_survival_and_risk_comparison.csv")
table2_path <- file.path(export_path, "stage11_table2_susceptible_only_comparison.csv")
posterior_summary_path <- file.path(export_path, "stage11_bayesian_posterior_summary.csv")
selection_registry_path <- file.path(export_path, "stage11_model_selection_registry.csv")
overall_plot_source_path <- file.path(export_path, "stage11_plot_source_overall.csv")
susceptible_plot_source_path <- file.path(export_path, "stage11_plot_source_susceptible_only.csv")
bayesian_note_path <- file.path(export_path, "stage11_bayesian_prior_note.md")
assumption_note_path <- file.path(export_path, "stage11_analysis_assumptions.md")
trace_reference_path <- file.path(export_path, "stage11_bayesian_trace_plots_reference_from_stage8a.pdf")
reuse_report_path <- file.path(export_path, "stage11_reused_existing_results_report.md")
spec_registry_path <- file.path(export_path, "stage11_spec_registry.csv")
manifest_path <- file.path(export_path, "stage11_export_manifest.csv")

write_csv_preserve_schema(table1_overall, table1_path)
write_csv_preserve_schema(table2_susceptible, table2_path)
write_csv_preserve_schema(bayes_posterior_summary, posterior_summary_path)
write_csv_preserve_schema(model_selection_registry, selection_registry_path)
write_csv_preserve_schema(table1_overall, overall_plot_source_path)
write_csv_preserve_schema(table2_susceptible, susceptible_plot_source_path)
write_text_file(bayesian_note_lines, bayesian_note_path)
write_text_file(analysis_assumption_lines, assumption_note_path)
write_text_file(reuse_report_lines, reuse_report_path)
copy_file_atomic(stage8_diagnostic_pdf, trace_reference_path)

spec_registry <- tibble(
  spec_type = c("integrated_master_spec", "bayesian_spec"),
  spec_path = unname(stage11_spec_paths)
)
write_csv_preserve_schema(spec_registry, spec_registry_path)

export_manifest <- tibble(
  file_name = c(
    basename(table1_path),
    basename(table2_path),
    basename(posterior_summary_path),
    basename(selection_registry_path),
    basename(overall_plot_source_path),
    basename(susceptible_plot_source_path),
    basename(plot_set1_pdf),
    basename(plot_set2_pdf),
    basename(trace_reference_path),
    basename(reuse_report_path),
    basename(bayesian_note_path),
    basename(assumption_note_path),
    basename(spec_registry_path),
    basename(manifest_path)
  ),
  description = c(
    "Table 1 overall survival probability and cumulative risk comparison across the four retained models.",
    "Table 2 susceptible-only survival probability and cumulative risk comparison for the retained cure models.",
    "Posterior summary for the retained Bayesian log-normal mixture cure models.",
    "Registry documenting which upstream model IDs were retained for each dataset version.",
    "Figure-ready source table for overall survival/risk plots.",
    "Figure-ready source table for susceptible-only plots.",
    "Two-page PDF containing overall survival and overall cumulative risk curves.",
    "Two-page PDF containing susceptible-only survival and susceptible-only cumulative risk curves.",
    "Copied Stage 8A diagnostic PDF retained as the trace-plot reference file.",
    "Markdown report documenting which upstream result files were reused and where they were used in Stage 11.",
    "Short note describing the retained Bayesian prior rationale and diagnostic carry-forward rule.",
    "Short note recording the simplifying assumptions used by Stage 11.",
    "Absolute-path registry for the new simplified specification documents.",
    "Manifest of Stage 11 exported outputs."
  ),
  file_path = c(
    table1_path,
    table2_path,
    posterior_summary_path,
    selection_registry_path,
    overall_plot_source_path,
    susceptible_plot_source_path,
    plot_set1_pdf,
    plot_set2_pdf,
    trace_reference_path,
    reuse_report_path,
    bayesian_note_path,
    assumption_note_path,
    spec_registry_path,
    manifest_path
  )
)

write_csv_preserve_schema(export_manifest, manifest_path)

message("Stage 11 simplified transition-only cure project completed.")
message("Export path: ", normalize_existing_path(export_path))
