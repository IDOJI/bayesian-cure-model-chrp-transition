# 🔴 Configure: Stage 8B paths and locked options ===============================
sys_name <- Sys.info()[["sysname"]]

project_root <- switch(
  sys_name,
  "Darwin" = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  "Windows" = "C:/Users/clair/Dropbox/Data Analysis/Survival Analysis On CHR-P_Results",
  "Linux" = "/mnt/data/chrp_results",
  stop("Unsupported OS: ", sys_name)
)

stage1_export_path <- file.path(project_root, "stage1_Backbone lock")
stage6_export_path <- file.path(project_root, "stage6_Cure-appropriateness screening")
stage8a_export_path <- file.path(project_root, "stage8A_Bayesian transition-only cure")
export_path <- file.path(project_root, "stage8B_Bayesian remission-sensitive competing-risk extension")
data_path <- stage1_export_path

pnu_site_label <- "PNU"
snu_site_label <- "SNU"

run_model_ids <- NULL
include_merged_incidence_site_supplementary <- TRUE
fit_prior_branches <- c("anchor_informed", "neutral_no_external_info")
fit_site_prior_families <- c("normal_0_1_main", "student_t3_0_1_sensitivity")

reuse_existing_stage8b_rds <- TRUE
save_full_stanfit_rds <- FALSE

stan_chains <- 4L
stan_iter <- 2000L
stan_warmup <- 1000L
stan_thin <- 1L
stan_seed <- 20260328L
stan_adapt_delta <- 0.95
stan_max_treedepth <- 12L
stan_refresh <- 0L

prior_predictive_draws <- 250L
posterior_prediction_draws <- 300L
remission_cut_years <- c(0, 1, 2, 5, 10, 100)
integration_step_year <- 0.05

ppc_tolerance_abs <- 0.12
ess_min_threshold <- 400
rhat_max_threshold <- 1.01
degenerate_draw_fraction_threshold <- 0.90
degenerate_subject_fraction_threshold <- 0.95
tiny_susceptible_prob <- 0.01
huge_susceptible_prob <- 0.99
tiny_median_years <- 0.05
huge_median_years <- 50

prior_predictive_horizons <- c(1L, 2L, 5L, 10L)
prior_tail_warning_mean_years <- 100
prior_tail_warning_q975_years <- 200

prior_materiality_risk <- 0.03
prior_materiality_false_positive_burden <- 0.03
prior_materiality_nb <- 0.01

stage1_bundle_file <- file.path(stage1_export_path, "stage1_backbone_bundle.rds")
stage1_analysis_datasets_file <- file.path(stage1_export_path, "stage1_analysis_datasets.rds")
stage1_horizon_registry_file <- file.path(stage1_export_path, "stage1_horizon_registry.csv")
stage1_threshold_registry_file <- file.path(stage1_export_path, "stage1_threshold_registry.csv")
stage1_scaling_registry_file <- file.path(stage1_export_path, "stage1_scaling_registry.csv")
stage1_dataset_registry_file <- file.path(stage1_export_path, "stage1_dataset_registry.csv")
stage1_export_manifest_file <- file.path(stage1_export_path, "stage1_export_manifest.csv")

stage6_screening_flag_csv <- file.path(stage6_export_path, "stage6_carry_forward_flag_table.csv")

stage8a_model_registry_csv <- file.path(stage8a_export_path, "bayes_stage8a_model_registry.csv")
stage8a_performance_csv <- file.path(stage8a_export_path, "bayes_stage8a_performance_classification_long.csv")

stage8b_diagnostic_pdf_file <- file.path(export_path, "bayes_stage8b_diagnostic_plots.pdf")
stage8b_plot_transition_png_file <- file.path(export_path, "bayes_stage8b_plot_transition_cif.png")
stage8b_plot_remission_png_file <- file.path(export_path, "bayes_stage8b_plot_remission_cif.png")
stage8b_plot_nb_png_file <- file.path(export_path, "bayes_stage8b_plot_net_benefit.png")
stage8b_plot_ppc_png_file <- file.path(export_path, "bayes_stage8b_plot_ppc.png")

# 🔴 Initialize: packages and runtime options ===============================
core_packages <- c(
  "readr", "dplyr", "tibble", "tidyr", "ggplot2", "survival",
  "matrixStats", "purrr", "tools"
)

missing_core_packages <- core_packages[
  !vapply(core_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_core_packages) > 0L) {
  stop(
    "Install required packages before running this script: ",
    paste(missing_core_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(survival)
  library(matrixStats)
  library(purrr)
})

options(stringsAsFactors = FALSE, scipen = 999)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

# 🔴 Define: utility helpers and I/O tools ===============================
`%||%` <- function(x, y) if (is.null(x)) y else x

assert_existing_file <- function(path, label = "File") {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

read_delimited_or_rds <- function(path) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  if (!file.exists(path)) {
    stop(sprintf("Input file does not exist: %s", path), call. = FALSE)
  }
  ext <- tolower(tools::file_ext(path))
  path_low <- tolower(path)
  if (grepl("\\.(csv|txt)(\\.gz)?$", path_low)) {
    return(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  }
  if (ext == "rds") return(readRDS(path))
  stop(sprintf("Unsupported file extension for `%s`.", path), call. = FALSE)
}

make_temp_output_path <- function(path, tag = "tmp") {
  dir_now <- dirname(path)
  base_now <- basename(path)
  stamp <- paste0(format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%08d", sample.int(99999999, 1)))
  if (grepl("\\.csv\\.gz$", base_now, ignore.case = TRUE)) {
    base_now <- sub("\\.csv\\.gz$", paste0("_", tag, "_", stamp, ".csv.gz"), base_now, ignore.case = TRUE)
  } else if (grepl("\\.csv$", base_now, ignore.case = TRUE)) {
    base_now <- sub("\\.csv$", paste0("_", tag, "_", stamp, ".csv"), base_now, ignore.case = TRUE)
  } else if (grepl("\\.pdf$", base_now, ignore.case = TRUE)) {
    base_now <- sub("\\.pdf$", paste0("_", tag, "_", stamp, ".pdf"), base_now, ignore.case = TRUE)
  } else if (grepl("\\.png$", base_now, ignore.case = TRUE)) {
    base_now <- sub("\\.png$", paste0("_", tag, "_", stamp, ".png"), base_now, ignore.case = TRUE)
  } else if (grepl("\\.rds$", base_now, ignore.case = TRUE)) {
    base_now <- sub("\\.rds$", paste0("_", tag, "_", stamp, ".rds"), base_now, ignore.case = TRUE)
  } else {
    base_now <- paste0(base_now, "_", tag, "_", stamp)
  }
  file.path(dir_now, base_now)
}

safe_promote_file <- function(tmp_path, final_path) {
  if (!file.exists(tmp_path)) {
    stop(sprintf("Temporary file does not exist: %s", tmp_path), call. = FALSE)
  }
  dir.create(dirname(final_path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(final_path)) unlink(final_path)
  ok <- file.rename(tmp_path, final_path)
  if (!isTRUE(ok)) {
    ok_copy <- file.copy(tmp_path, final_path, overwrite = TRUE)
    unlink(tmp_path)
    if (!isTRUE(ok_copy)) {
      stop(sprintf("Failed to replace file: %s", final_path), call. = FALSE)
    }
  }
  invisible(TRUE)
}

write_csv_preserve_schema <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- make_temp_output_path(path, tag = "tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  readr::write_csv(df, tmp)
  safe_promote_file(tmp, path)
}

pdf_file_is_usable <- function(path) {
  if (!file.exists(path)) return(FALSE)
  info <- file.info(path)
  isTRUE(!is.na(info$size[[1]]) && info$size[[1]] > 0L)
}

png_file_is_usable <- function(path) {
  if (!file.exists(path)) return(FALSE)
  info <- file.info(path)
  isTRUE(!is.na(info$size[[1]]) && info$size[[1]] > 0L)
}

nrow_or_zero <- function(df) {
  if (is.null(df) || !inherits(df, "data.frame")) return(0L)
  nrow(df)
}

bind_rows_safe <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) == 0L) return(tibble())
  bind_rows(x)
}

drop_existing_columns <- function(df, cols) {
  if (is.null(df) || !inherits(df, "data.frame") || length(cols) == 0L) return(df)
  dplyr::select(df, -any_of(cols))
}

left_join_replacing_columns <- function(x, y, by) {
  if (is.null(x) || !inherits(x, "data.frame")) return(x)
  if (is.null(y) || !inherits(y, "data.frame") || nrow(y) == 0L) return(x)
  by_cols <- unname(by)
  if (any(!by_cols %in% names(x)) || any(!by_cols %in% names(y))) return(x)
  join_cols <- setdiff(names(y), by_cols)
  x %>% drop_existing_columns(join_cols) %>% left_join(y, by = by)
}

first_existing_name <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) return(NULL)
  hit[[1L]]
}

format_number <- function(x, digits = 3L) {
  if (length(x) == 0L || is.na(x) || !is.finite(x)) return("NA")
  formatC(as.numeric(x), format = "f", digits = digits)
}

elapsed_seconds <- function(start_time) {
  as.numeric(difftime(Sys.time(), start_time, units = "secs"))
}

emit_progress <- function(done, total, model_id, detail) {
  pct <- if (total <= 0L) 100 else 100 * done / total
  message(sprintf("[%d/%d | %5.1f%%] %s %s", as.integer(done), as.integer(total), pct, model_id, detail))
}

is_localhost_connection_warning <- function(w) {
  grepl("closing unused connection .*<-localhost:", conditionMessage(w))
}

safe_quantile <- function(x, p) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

summary_scalar <- function(x) {
  tibble(
    mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    q025 = safe_quantile(x, 0.025),
    q50 = safe_quantile(x, 0.500),
    q975 = safe_quantile(x, 0.975)
  )
}

summarize_cols_matrix <- function(mat) {
  if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1L)
  q <- matrixStats::colQuantiles(mat, probs = c(0.025, 0.500, 0.975), na.rm = TRUE)
  tibble(
    mean = matrixStats::colMeans2(mat, na.rm = TRUE),
    q025 = q[, 1],
    q50 = q[, 2],
    q975 = q[, 3]
  )
}

# 🔴 Define: Stage 1 and Stage 6 input loaders ===============================
empty_stage1_bundle <- function() {
  list(
    datasets = NULL,
    registries = list(
      dataset_registry = tibble(),
      horizon_registry = tibble(),
      threshold_registry = tibble(),
      scaling_registry = tibble(),
      metadata_registry = tibble()
    ),
    config = list()
  )
}

load_stage1_backbone <- function(bundle_file, datasets_file, dataset_registry_file, horizon_registry_file, threshold_registry_file, scaling_registry_file) {
  out <- empty_stage1_bundle()

  if (file.exists(bundle_file)) {
    obj <- readRDS(bundle_file)
    if (is.list(obj)) {
      out$config <- obj$config %||% list()
      out$datasets <- obj$datasets %||% NULL
      regs <- obj$registries %||% list()
      out$registries$dataset_registry <- tibble::as_tibble(regs$dataset_registry %||% tibble())
      out$registries$horizon_registry <- tibble::as_tibble(regs$horizon_registry %||% tibble())
      out$registries$threshold_registry <- tibble::as_tibble(regs$threshold_registry %||% tibble())
      out$registries$scaling_registry <- tibble::as_tibble(regs$scaling_registry %||% tibble())
      out$registries$metadata_registry <- tibble::as_tibble(regs$metadata_registry %||% tibble())
    }
  }

  if (is.null(out$datasets) && file.exists(datasets_file)) out$datasets <- readRDS(datasets_file)
  if (nrow_or_zero(out$registries$dataset_registry) == 0L && file.exists(dataset_registry_file)) out$registries$dataset_registry <- read_delimited_or_rds(dataset_registry_file)
  if (nrow_or_zero(out$registries$horizon_registry) == 0L && file.exists(horizon_registry_file)) out$registries$horizon_registry <- read_delimited_or_rds(horizon_registry_file)
  if (nrow_or_zero(out$registries$threshold_registry) == 0L && file.exists(threshold_registry_file)) out$registries$threshold_registry <- read_delimited_or_rds(threshold_registry_file)
  if (nrow_or_zero(out$registries$scaling_registry) == 0L && file.exists(scaling_registry_file)) out$registries$scaling_registry <- read_delimited_or_rds(scaling_registry_file)

  if (is.null(out$datasets) || !is.list(out$datasets) || length(out$datasets) == 0L) {
    stop("Stage 1 datasets could not be loaded from stage1_backbone_bundle.rds or stage1_analysis_datasets.rds.", call. = FALSE)
  }

  req_datasets <- c("PNU", "SNU", "merged")
  if (any(!req_datasets %in% names(out$datasets))) {
    stop("Stage 1 datasets must contain PNU, SNU, and merged.", call. = FALSE)
  }

  out$registries$horizon_registry <- tibble::as_tibble(out$registries$horizon_registry)
  out$registries$threshold_registry <- tibble::as_tibble(out$registries$threshold_registry)
  out$registries$dataset_registry <- tibble::as_tibble(out$registries$dataset_registry)
  out$registries$scaling_registry <- tibble::as_tibble(out$registries$scaling_registry)

  out
}

build_threshold_vector_from_stage1 <- function(threshold_registry, fallback = c(0.05, 0.10, 0.15, 0.20)) {
  out <- tibble::as_tibble(threshold_registry)
  if (nrow(out) == 0L || !("threshold" %in% names(out))) {
    return(as.numeric(fallback))
  }
  sort(unique(safe_numeric(out$threshold)))
}

normalize_dataset_time <- function(df) {
  df <- tibble::as_tibble(df)
  if (!("time_year" %in% names(df))) {
    if (!("days_followup" %in% names(df))) {
      stop("Dataset must contain time_year or days_followup.", call. = FALSE)
    }
    df$time_year <- as.numeric(df$days_followup) / 365.25
  }
  df
}

normalize_unique_person_id <- function(df) {
  df <- tibble::as_tibble(df)
  if (!("unique_person_id" %in% names(df))) {
    if (!all(c("site", "id") %in% names(df))) {
      stop("Dataset must contain unique_person_id or both site and id.", call. = FALSE)
    }
    df$unique_person_id <- paste0(as.character(df$site), "::", as.character(df$id))
  }
  df
}

empty_stage6_lookup <- function() {
  tibble(
    dataset = character(),
    formula_anchor = character(),
    cure_model_eligibility_flag = character(),
    primary_gate_method = character(),
    primary_gate_flag = character(),
    receus_primary_class = character(),
    presence_modifier_flag = character(),
    cure_presence_support_flag = character(),
    followup_contradiction_flag = character(),
    followup_not_contradicted_flag = character(),
    screening_note = character()
  )
}

normalize_stage6_formula_anchor <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    x %in% c("base", "Base") ~ "base",
    x %in% c("interaction", "Interaction") ~ "interaction",
    x %in% c("site_added", "Site-added", "site-added") ~ "site_added",
    x %in% c("site_interaction", "Site + interaction", "site+interaction") ~ "site_interaction",
    TRUE ~ x
  )
}

expand_stage6_dataset_key <- function(dataset_value, formula_anchor_value = NA_character_) {
  dataset_value <- trimws(as.character(dataset_value %||% NA_character_))
  formula_anchor_value <- trimws(as.character(formula_anchor_value %||% NA_character_))

  if (identical(dataset_value, "merged__site_free")) {
    return(tibble(dataset = "merged", formula_anchor = c("base", "interaction")))
  }
  if (identical(dataset_value, "merged__site_adjusted")) {
    return(tibble(dataset = "merged", formula_anchor = c("site_added", "site_interaction")))
  }
  if (dataset_value %in% c("PNU", "SNU", "merged")) {
    return(tibble(
      dataset = dataset_value,
      formula_anchor = if (!is.na(formula_anchor_value) && nzchar(formula_anchor_value)) formula_anchor_value else "ALL"
    ))
  }

  tibble(
    dataset = dataset_value,
    formula_anchor = if (!is.na(formula_anchor_value) && nzchar(formula_anchor_value)) formula_anchor_value else "ALL"
  )
}

read_stage6_screening_lookup <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(empty_stage6_lookup())
  }

  df <- tibble::as_tibble(read_delimited_or_rds(path))

  dataset_col <- first_existing_name(df, c("dataset_key", "dataset", "source_dataset", "cohort", "screening_dataset_key"))
  formula_col <- first_existing_name(df, c("formula_anchor", "formula_name", "formula_type"))
  eligibility_col <- first_existing_name(df, c("cure_model_eligibility_flag", "stage6_final_class", "final_decision_flag", "screening_flag"))
  primary_method_col <- first_existing_name(df, c("primary_gate_method"))
  primary_flag_col <- first_existing_name(df, c("primary_gate_flag", "receus_primary_flag"))
  receus_col <- first_existing_name(df, c("receus_primary_class", "primary_gate_class"))
  presence_modifier_col <- first_existing_name(df, c("presence_modifier_flag"))
  presence_col <- first_existing_name(df, c("cure_presence_support_flag", "presence_support_flag", "presence_modifier_support"))
  followup_contradiction_col <- first_existing_name(df, c("followup_contradiction_flag"))
  followup_col <- first_existing_name(df, c("followup_not_contradicted_flag", "followup_contradiction_resolved_flag"))
  note_col <- first_existing_name(df, c("screening_note", "screening_detail", "screening_context", "primary_gate_method"))

  if (is.null(dataset_col) || is.null(eligibility_col)) {
    return(empty_stage6_lookup())
  }

  raw_out <- tibble(
    dataset_value = as.character(df[[dataset_col]]),
    formula_anchor_value = if (!is.null(formula_col)) normalize_stage6_formula_anchor(df[[formula_col]]) else NA_character_,
    cure_model_eligibility_flag = as.character(df[[eligibility_col]]),
    primary_gate_method = if (!is.null(primary_method_col)) as.character(df[[primary_method_col]]) else NA_character_,
    primary_gate_flag = if (!is.null(primary_flag_col)) as.character(df[[primary_flag_col]]) else NA_character_,
    receus_primary_class = if (!is.null(receus_col)) as.character(df[[receus_col]]) else NA_character_,
    presence_modifier_flag = if (!is.null(presence_modifier_col)) as.character(df[[presence_modifier_col]]) else NA_character_,
    cure_presence_support_flag = if (!is.null(presence_col)) as.character(df[[presence_col]]) else NA_character_,
    followup_contradiction_flag = if (!is.null(followup_contradiction_col)) as.character(df[[followup_contradiction_col]]) else NA_character_,
    followup_not_contradicted_flag = if (!is.null(followup_col)) as.character(df[[followup_col]]) else NA_character_,
    screening_note = if (!is.null(note_col)) as.character(df[[note_col]]) else NA_character_
  )

  out <- bind_rows(lapply(seq_len(nrow(raw_out)), function(i) {
    one <- raw_out[i, , drop = FALSE]
    expanded <- expand_stage6_dataset_key(one$dataset_value[[1]], one$formula_anchor_value[[1]])
    carried <- dplyr::select(one, -dataset_value, -formula_anchor_value)
    carried <- carried[rep(1L, nrow(expanded)), , drop = FALSE]
    dplyr::bind_cols(expanded, carried)
  })) %>%
    filter(dataset %in% c("PNU", "SNU", "merged")) %>%
    mutate(formula_anchor = ifelse(is.na(formula_anchor) | formula_anchor == "", "ALL", formula_anchor)) %>%
    group_by(dataset, formula_anchor) %>%
    summarise(
      cure_model_eligibility_flag = paste(unique(na.omit(cure_model_eligibility_flag)), collapse = "|"),
      primary_gate_method = paste(unique(na.omit(primary_gate_method)), collapse = "|"),
      primary_gate_flag = paste(unique(na.omit(primary_gate_flag)), collapse = "|"),
      receus_primary_class = paste(unique(na.omit(receus_primary_class)), collapse = "|"),
      presence_modifier_flag = paste(unique(na.omit(presence_modifier_flag)), collapse = "|"),
      cure_presence_support_flag = paste(unique(na.omit(cure_presence_support_flag)), collapse = "|"),
      followup_contradiction_flag = paste(unique(na.omit(followup_contradiction_flag)), collapse = "|"),
      followup_not_contradicted_flag = paste(unique(na.omit(followup_not_contradicted_flag)), collapse = "|"),
      screening_note = paste(unique(na.omit(screening_note)), collapse = "; "),
      .groups = "drop"
    )

  if (nrow(out) == 0L) {
    return(empty_stage6_lookup())
  }
  out
}

build_stage6_model_lookup <- function(screening_lookup, model_grid) {
  if (nrow_or_zero(screening_lookup) == 0L) {
    return(model_grid %>% transmute(
      model_id = model_id,
      cure_model_eligibility_flag = NA_character_,
      primary_gate_method = NA_character_,
      primary_gate_flag = NA_character_,
      receus_primary_class = NA_character_,
      presence_modifier_flag = NA_character_,
      cure_presence_support_flag = NA_character_,
      followup_contradiction_flag = NA_character_,
      followup_not_contradicted_flag = NA_character_,
      screening_note = NA_character_
    ))
  }

  model_base <- model_grid %>% distinct(model_id, dataset, formula_anchor)

  exact_lookup <- model_base %>%
    left_join(
      screening_lookup %>% filter(formula_anchor != "ALL"),
      by = c("dataset", "formula_anchor")
    )

  all_lookup <- model_base %>%
    left_join(
      screening_lookup %>%
        filter(formula_anchor == "ALL") %>%
        select(dataset, cure_model_eligibility_flag, primary_gate_method, primary_gate_flag,
               receus_primary_class, presence_modifier_flag, cure_presence_support_flag,
               followup_contradiction_flag, followup_not_contradicted_flag, screening_note) %>%
        rename_with(~ paste0(.x, "_all"), -dataset),
      by = "dataset"
    )

  exact_lookup %>%
    left_join(all_lookup, by = c("model_id", "dataset", "formula_anchor")) %>%
    transmute(
      model_id = model_id,
      cure_model_eligibility_flag = coalesce(cure_model_eligibility_flag, cure_model_eligibility_flag_all),
      primary_gate_method = coalesce(primary_gate_method, primary_gate_method_all),
      primary_gate_flag = coalesce(primary_gate_flag, primary_gate_flag_all),
      receus_primary_class = coalesce(receus_primary_class, receus_primary_class_all),
      presence_modifier_flag = coalesce(presence_modifier_flag, presence_modifier_flag_all),
      cure_presence_support_flag = coalesce(cure_presence_support_flag, cure_presence_support_flag_all),
      followup_contradiction_flag = coalesce(followup_contradiction_flag, followup_contradiction_flag_all),
      followup_not_contradicted_flag = coalesce(followup_not_contradicted_flag, followup_not_contradicted_flag_all),
      screening_note = coalesce(screening_note, screening_note_all)
    )
}

# 🔴 Define: Stage 8A loaders and shared horizon governance ===============================
load_stage8a_outputs <- function(model_registry_csv, performance_csv) {
  out <- list(
    model_registry = tibble(),
    performance_classification_long = tibble()
  )
  if (file.exists(model_registry_csv)) out$model_registry <- tibble::as_tibble(read_delimited_or_rds(model_registry_csv))
  if (file.exists(performance_csv)) out$performance_classification_long <- tibble::as_tibble(read_delimited_or_rds(performance_csv))
  out
}

map_stage8a_to_stage8b_structural_id <- function(structural_model_id, dataset_key = NA_character_) {
  structural_model_id <- as.character(structural_model_id %||% NA_character_)
  dataset_key <- as.character(dataset_key %||% NA_character_)

  vapply(seq_along(structural_model_id), function(ii) {
    id_now <- structural_model_id[[ii]]
    ds_now <- dataset_key[[ii]]

    if (!nzchar(id_now) || is.na(id_now)) return(id_now)

    if (grepl("^(PNU|SNU)-L[01]$", id_now)) {
      latency_branch <- sub("^(PNU|SNU)-(L[01])$", "\\2", id_now)
      remission_branch <- sub("^L", "R", latency_branch)
      return(paste0(sub("-(L[01])$", "", id_now), "-", latency_branch, remission_branch))
    }

    if (grepl("^MERGED-I[01]-L[01]S[01]$", id_now)) {
      incidence_tag <- sub("^MERGED-(I[01])-L[01]S[01]$", "\\1", id_now)
      latency_branch <- sub("^MERGED-I[01]-(L[01]S[01])$", "\\1", id_now)
      remission_branch <- sub("^L", "R", latency_branch)
      prefix <- if (identical(incidence_tag, "I1")) "MERGED-I1" else "MERGED"
      return(paste(prefix, latency_branch, remission_branch, sep = "-"))
    }

    if (identical(ds_now, "merged") && grepl("^MERGED-L[01]S[01]$", id_now)) {
      latency_branch <- sub("^MERGED-(L[01]S[01])$", "\\1", id_now)
      remission_branch <- sub("^L", "R", latency_branch)
      return(paste("MERGED", latency_branch, remission_branch, sep = "-"))
    }

    id_now
  }, character(1))
}

build_stage8a_delta_keys <- function(stage8a_outputs) {
  perf <- tibble::as_tibble(stage8a_outputs$performance_classification_long)
  out <- list(
    cohort_key = tibble(),
    class_key = tibble()
  )

  if (nrow_or_zero(perf) == 0L) {
    return(out)
  }

  perf_dataset_col <- first_existing_name(perf, c("dataset", "dataset_key"))
  perf_horizon_col <- first_existing_name(perf, c("horizon_year", "horizon"))
  perf_threshold_col <- first_existing_name(perf, c("threshold", "risk_threshold"))
  perf_transition_col <- first_existing_name(perf, c("risk_mean", "meanRisk_Bayes_mean", "transition_risk_mean", "transition_only_risk_mean", "transition_cif_mean"))
  perf_cure_col <- first_existing_name(perf, c("cure_fraction_mean", "cohort_mean_cure_fraction_mean"))

  has_model_keys <- all(c("structural_model_id", "formula_anchor", "family_code") %in% names(perf))
  has_branch_keys <- all(c("prior_branch", "site_prior_family") %in% names(perf))
  if (!has_model_keys || !has_branch_keys || is.null(perf_dataset_col) || is.null(perf_horizon_col)) {
    return(out)
  }

  perf_base <- perf %>%
    mutate(
      dataset_key = as.character(.data[[perf_dataset_col]]),
      horizon = as.integer(safe_numeric(.data[[perf_horizon_col]])),
      structural_model_id = map_stage8a_to_stage8b_structural_id(
        structural_model_id = as.character(structural_model_id),
        dataset_key = dataset_key
      )
    )

  if (!is.null(perf_transition_col) || !is.null(perf_cure_col)) {
    out$cohort_key <- perf_base %>%
      transmute(
        dataset_key = dataset_key,
        structural_model_id = as.character(structural_model_id),
        formula_anchor = as.character(formula_anchor),
        family_code = as.character(family_code),
        prior_branch = as.character(prior_branch),
        site_prior_family = as.character(site_prior_family),
        horizon = horizon,
        stage8a_transition_risk = if (!is.null(perf_transition_col)) safe_numeric(.data[[perf_transition_col]]) else NA_real_,
        stage8a_cure_fraction = if (!is.null(perf_cure_col)) safe_numeric(.data[[perf_cure_col]]) else NA_real_
      ) %>%
      distinct()
  }

  class_metric_cols <- c("false_positive_burden_mean", "FP100_mean", "NB_mean", "PPV_mean", "TPR_mean")
  if (!is.null(perf_threshold_col) && all(class_metric_cols %in% names(perf_base))) {
    out$class_key <- perf_base %>%
      transmute(
        dataset_key = dataset_key,
        structural_model_id = as.character(structural_model_id),
        formula_anchor = as.character(formula_anchor),
        family_code = as.character(family_code),
        prior_branch = as.character(prior_branch),
        site_prior_family = as.character(site_prior_family),
        horizon = horizon,
        threshold = as.numeric(safe_numeric(.data[[perf_threshold_col]])),
        stage8a_false_positive_burden = safe_numeric(false_positive_burden_mean),
        stage8a_FP100 = safe_numeric(FP100_mean),
        stage8a_NB = safe_numeric(NB_mean),
        stage8a_PPV = safe_numeric(PPV_mean),
        stage8a_TPR = safe_numeric(TPR_mean)
      ) %>%
      distinct()
  }

  out
}

stage8b_support_tier_from_grid <- function(dataset, horizon_year) {
  dataset <- as.character(dataset)
  horizon_year <- as.integer(safe_numeric(horizon_year))
  dplyr::case_when(
    dataset == "PNU" & horizon_year == 1L ~ "primary_supported",
    dataset == "PNU" & horizon_year == 2L ~ "sensitivity",
    dataset == "PNU" & horizon_year >= 3L ~ "projection",
    dataset %in% c("SNU", "merged") & horizon_year %in% c(1L, 2L) ~ "primary_supported",
    dataset %in% c("SNU", "merged") & horizon_year <= 5L ~ "secondary",
    TRUE ~ "projection"
  )
}

stage8b_horizon_evidence_from_grid <- function(dataset, horizon_year) {
  dataset <- as.character(dataset)
  horizon_year <- as.integer(safe_numeric(horizon_year))
  dplyr::case_when(
    dataset == "PNU" & horizon_year == 1L ~ "directly_observed_data_supported",
    dataset == "PNU" & horizon_year == 2L ~ "partly_model_dependent",
    dataset == "PNU" & horizon_year >= 3L ~ "mostly_extrapolated",
    dataset %in% c("SNU", "merged") & horizon_year %in% c(1L, 2L) ~ "directly_observed_data_supported",
    dataset %in% c("SNU", "merged") & horizon_year <= 5L ~ "partly_model_dependent",
    TRUE ~ "mostly_extrapolated"
  )
}

stage8b_claim_restriction_from_evidence <- function(horizon_evidence_class) {
  dplyr::case_when(
    horizon_evidence_class == "directly_observed_data_supported" ~ "primary_claim_allowed",
    horizon_evidence_class == "partly_model_dependent" ~ "secondary_or_sensitivity_only",
    TRUE ~ "projection_only"
  )
}

stage8b_interpretation_note_from_grid <- function(dataset, horizon_year, support_tier, horizon_evidence_class, claim_restriction_flag) {
  dataset <- as.character(dataset)
  horizon_year <- as.integer(safe_numeric(horizon_year))
  dplyr::case_when(
    claim_restriction_flag == "primary_claim_allowed" ~ "Primary supported horizon with comparatively direct follow-up support.",
    dataset == "PNU" & horizon_year == 2L ~ "Sensitivity horizon for PNU; partly model-dependent and not for primary claims.",
    claim_restriction_flag == "secondary_or_sensitivity_only" ~ "Secondary or sensitivity horizon with explicit model dependence.",
    TRUE ~ "Projection-dominant horizon."
  )
}

build_horizon_annotation_from_stage1 <- function(horizon_registry, datasets, horizons) {
  dataset_levels <- if (is.list(datasets)) names(datasets) else as.character(datasets)
  out <- tibble::as_tibble(horizon_registry)

  if (nrow(out) == 0L) {
    out <- tidyr::crossing(dataset = dataset_levels, horizon_year = as.integer(horizons))
  }

  if (!("dataset" %in% names(out))) {
    if ("dataset_key" %in% names(out)) {
      out$dataset <- out$dataset_key
    } else {
      stop("Stage 1 horizon registry must contain dataset or dataset_key.", call. = FALSE)
    }
  }
  if (!("dataset_key" %in% names(out))) out$dataset_key <- out$dataset
  if (!("horizon_year" %in% names(out))) {
    horizon_col <- first_existing_name(out, c("horizon", "year"))
    if (is.null(horizon_col)) {
      stop("Stage1 horizon registry must contain horizon_year (or horizon / year).", call. = FALSE)
    }
    out$horizon_year <- out[[horizon_col]]
  }
  if (!("support_tier" %in% names(out))) out$support_tier <- NA_character_
  if (!("support_tier_standard" %in% names(out))) out$support_tier_standard <- NA_character_
  if (!("horizon_evidence_class" %in% names(out))) out$horizon_evidence_class <- NA_character_
  if (!("claim_restriction_flag" %in% names(out))) out$claim_restriction_flag <- NA_character_
  if (!("interpretation_note" %in% names(out))) out$interpretation_note <- NA_character_

  out %>%
    mutate(
      dataset = as.character(dataset),
      dataset_key = coalesce(as.character(dataset_key), as.character(dataset)),
      horizon_year = as.integer(safe_numeric(horizon_year)),
      support_tier = coalesce(
        dplyr::na_if(as.character(support_tier), ""),
        dplyr::na_if(as.character(support_tier_standard), ""),
        stage8b_support_tier_from_grid(dataset, horizon_year)
      ),
      horizon_evidence_class = coalesce(
        dplyr::na_if(as.character(horizon_evidence_class), ""),
        stage8b_horizon_evidence_from_grid(dataset, horizon_year)
      ),
      claim_restriction_flag = coalesce(
        dplyr::na_if(as.character(claim_restriction_flag), ""),
        stage8b_claim_restriction_from_evidence(horizon_evidence_class)
      ),
      interpretation_note = coalesce(
        dplyr::na_if(as.character(interpretation_note), ""),
        stage8b_interpretation_note_from_grid(dataset, horizon_year, support_tier, horizon_evidence_class, claim_restriction_flag)
      )
    ) %>%
    transmute(
      dataset = as.character(dataset),
      dataset_key = as.character(dataset_key),
      horizon_year = as.integer(horizon_year),
      support_tier = as.character(support_tier),
      horizon_evidence_class = as.character(horizon_evidence_class),
      claim_restriction_flag = as.character(claim_restriction_flag),
      interpretation_note = as.character(interpretation_note)
    ) %>%
    filter(dataset %in% c("PNU", "SNU", "merged"), horizon_year %in% as.integer(horizons)) %>%
    distinct(dataset, horizon_year, .keep_all = TRUE) %>%
    arrange(factor(dataset, levels = c("PNU", "SNU", "merged")), horizon_year)
}

make_stage8b_support_registry <- function(horizon_registry_stage1, datasets = c("PNU", "SNU", "merged"), horizons = 1:10) {
  out <- build_horizon_annotation_from_stage1(horizon_registry_stage1, datasets, horizons)

  if (nrow(out) == 0L) {
    return(tidyr::expand_grid(dataset = datasets, horizon_year = as.integer(horizons)) %>%
             mutate(
               dataset_key = dataset,
               support_tier = NA_character_,
               support_tier_standard = NA_character_,
               horizon_evidence_class = NA_character_,
               claim_restriction_flag = NA_character_,
               interpretation_note = NA_character_
             ))
  }

  out %>%
    transmute(
      dataset = as.character(dataset),
      dataset_key = as.character(dataset_key),
      horizon_year = as.integer(safe_numeric(horizon_year)),
      support_tier = as.character(support_tier),
      support_tier_standard = as.character(support_tier),
      horizon_evidence_class = as.character(horizon_evidence_class),
      claim_restriction_flag = as.character(claim_restriction_flag),
      interpretation_note = as.character(interpretation_note)
    ) %>%
    filter(dataset %in% datasets, horizon_year %in% as.integer(horizons)) %>%
    distinct(dataset, horizon_year, .keep_all = TRUE) %>%
    arrange(factor(dataset, levels = datasets), horizon_year)
}

build_risk_set_support_table <- function(datasets, horizons) {
  bind_rows(lapply(names(datasets), function(ds) {
    dat <- normalize_unique_person_id(normalize_dataset_time(datasets[[ds]]))
    bind_rows(lapply(as.integer(horizons), function(hh) {
      tibble(
        dataset_key = ds,
        horizon = as.integer(hh),
        n_total = nrow(dat),
        n_transition_by_horizon = sum(dat$status_num == 1L & dat$time_year <= hh, na.rm = TRUE),
        n_remission_by_horizon = sum(dat$status_num == 2L & dat$time_year <= hh, na.rm = TRUE),
        n_right_censoring_by_horizon = sum(dat$status_num == 0L & dat$time_year <= hh, na.rm = TRUE),
        n_observed_at_or_beyond_horizon = sum(dat$time_year >= hh, na.rm = TRUE),
        risk_set_fraction = mean(dat$time_year >= hh, na.rm = TRUE)
      )
    }))
  }))
}

# 🔴 Define: Stage 8B model menu and prior specifications ===============================
build_stage8b_model_grid <- function(
  include_merged_incidence_site_supplementary = TRUE,
  fit_prior_branches = c("anchor_informed", "neutral_no_external_info"),
  fit_site_prior_families = c("normal_0_1_main", "student_t3_0_1_sensitivity")
) {
  single_rows <- tibble::tribble(
    ~dataset, ~structural_model_id, ~formula_anchor, ~transition_latency_branch, ~remission_branch, ~incidence_site_indicator, ~latency_site_indicator, ~remission_site_indicator, ~interaction_indicator, ~site_placement_label, ~is_supplementary_branch,
    "PNU", "PNU-L0R0", "base", "L0", "R0", FALSE, FALSE, FALSE, FALSE, "site_in_neither", FALSE,
    "PNU", "PNU-L1R1", "interaction", "L1", "R1", FALSE, FALSE, FALSE, TRUE, "site_in_neither", FALSE,
    "SNU", "SNU-L0R0", "base", "L0", "R0", FALSE, FALSE, FALSE, FALSE, "site_in_neither", FALSE,
    "SNU", "SNU-L1R1", "interaction", "L1", "R1", FALSE, FALSE, FALSE, TRUE, "site_in_neither", FALSE
  )

  merged_rows <- tibble::tribble(
    ~dataset, ~structural_model_id, ~formula_anchor, ~transition_latency_branch, ~remission_branch, ~incidence_site_indicator, ~latency_site_indicator, ~remission_site_indicator, ~interaction_indicator, ~site_placement_label, ~is_supplementary_branch,
    "merged", "MERGED-L0S0-R0S0", "base", "L0S0", "R0S0", FALSE, FALSE, FALSE, FALSE, "site_in_neither", FALSE,
    "merged", "MERGED-L1S0-R1S0", "interaction", "L1S0", "R1S0", FALSE, FALSE, FALSE, TRUE, "site_in_neither", FALSE,
    "merged", "MERGED-L0S1-R0S1", "site_added", "L0S1", "R0S1", FALSE, TRUE, TRUE, FALSE, "site_in_latency_only", FALSE,
    "merged", "MERGED-L1S1-R1S1", "site_interaction", "L1S1", "R1S1", FALSE, TRUE, TRUE, TRUE, "site_in_latency_only", FALSE
  )

  if (isTRUE(include_merged_incidence_site_supplementary)) {
    merged_rows <- bind_rows(
      merged_rows,
      tibble::tribble(
        ~dataset, ~structural_model_id, ~formula_anchor, ~transition_latency_branch, ~remission_branch, ~incidence_site_indicator, ~latency_site_indicator, ~remission_site_indicator, ~interaction_indicator, ~site_placement_label, ~is_supplementary_branch,
        "merged", "MERGED-I1-L0S0-R0S0", "base", "L0S0", "R0S0", TRUE, FALSE, FALSE, FALSE, "site_in_incidence_only", TRUE,
        "merged", "MERGED-I1-L1S0-R1S0", "interaction", "L1S0", "R1S0", TRUE, FALSE, FALSE, TRUE, "site_in_incidence_only", TRUE,
        "merged", "MERGED-I1-L0S1-R0S1", "site_added", "L0S1", "R0S1", TRUE, TRUE, TRUE, FALSE, "site_in_both", TRUE,
        "merged", "MERGED-I1-L1S1-R1S1", "site_interaction", "L1S1", "R1S1", TRUE, TRUE, TRUE, TRUE, "site_in_both", TRUE
      )
    )
  }

  family_rows <- tibble::tribble(
    ~family_code, ~latency_family, ~family_id,
    "E", "exponential", 1L,
    "W", "weibull", 2L,
    "LN", "lognormal", 3L,
    "LL", "loglogistic", 4L
  )

  base_grid <- tidyr::crossing(bind_rows(single_rows, merged_rows), family_rows) %>%
    mutate(has_any_site_term = incidence_site_indicator | latency_site_indicator | remission_site_indicator)

  rows_out <- list()
  for (ii in seq_len(nrow(base_grid))) {
    row_i <- base_grid[ii, , drop = FALSE]
    site_prior_candidates <- if (isTRUE(row_i$has_any_site_term[[1]])) fit_site_prior_families else "not_applicable"
    rows_out[[ii]] <- tidyr::crossing(
      row_i,
      tibble(prior_branch = fit_prior_branches, site_prior_family = site_prior_candidates)
    )
  }

  bind_rows(rows_out) %>%
    mutate(
      dataset_key = dataset,
      branch = "Stage8B",
      risk_scale = "transition_cif_competing",
      retained_fit_id = paste(structural_model_id, family_code, site_prior_family, sep = "__"),
      model_id = paste(structural_model_id, family_code, prior_branch, site_prior_family, sep = "__")
    ) %>%
    arrange(factor(dataset, levels = c("PNU", "SNU", "merged")), structural_model_id, family_code, prior_branch, site_prior_family)
}

build_stage8b_prior_specs <- function() {
  list(
    anchor_informed = list(
      prior_branch = "anchor_informed",
      alpha_prior_center = -9.581369553169,
      mu_beta_inc_base = c(0.419871845822, 0.907608052926, 0.586202561451, 0.466865123863, 0.037997248763),
      sd_beta_inc_base = c(0.132789397422, 0.173731076538, 0.191221553945, 0.270393197518, 0.302838606651),
      sd_delta = 3.5,
      sd_gamma0 = 2.5,
      sd_gamma = 1.0,
      sd_rho_W = 0.35,
      sd_log_sigma_LN = 0.50,
      sd_psi_LL = 0.50,
      sd_xi0 = 2.5,
      sd_rem = 1.5,
      sd_log_rho_rem = 1.5,
      sd_site = 1.0
    ),
    neutral_no_external_info = list(
      prior_branch = "neutral_no_external_info",
      alpha_prior_center = 0.0,
      mu_beta_inc_base = rep(0, 5),
      sd_beta_inc_base = rep(2.0, 5),
      sd_delta = 3.5,
      sd_gamma0 = 2.5,
      sd_gamma = 1.0,
      sd_rho_W = 0.35,
      sd_log_sigma_LN = 0.50,
      sd_psi_LL = 0.50,
      sd_xi0 = 2.5,
      sd_rem = 1.5,
      sd_log_rho_rem = 1.5,
      sd_site = 1.0
    )
  )
}

make_stage8b_design_bundle <- function(df, model_row, prior_spec, snu_label, remission_cut_years) {
  df <- normalize_unique_person_id(normalize_dataset_time(df))
  z_i <- as.integer(df$sex_num)
  x20_i <- as.integer(df$age_exact_entry >= 20 & df$age_exact_entry < 30)
  x30_i <- as.integer(df$age_exact_entry >= 30)
  s_i <- as.integer(df$site == snu_label)

  X_inc <- cbind(
    sex_num = z_i,
    age20_29 = x20_i,
    age30plus = x30_i,
    sex_x_age20_29 = z_i * x20_i,
    sex_x_age30plus = z_i * x30_i
  )
  mu_beta_inc <- prior_spec$mu_beta_inc_base
  sd_beta_inc <- prior_spec$sd_beta_inc_base

  if (isTRUE(model_row$incidence_site_indicator[[1]])) {
    X_inc <- cbind(X_inc, site_SNU = s_i)
    mu_beta_inc <- c(mu_beta_inc, 0)
    sd_beta_inc <- c(sd_beta_inc, prior_spec$sd_site)
  }

  a_i <- as.numeric(df$age_s)
  az_i <- a_i * z_i

  X_lat <- switch(
    model_row$transition_latency_branch,
    L0 = cbind(age_s = a_i, sex_num = z_i),
    L1 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i),
    L0S0 = cbind(age_s = a_i, sex_num = z_i),
    L1S0 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i),
    L0S1 = cbind(age_s = a_i, sex_num = z_i, site_SNU = s_i),
    L1S1 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i, site_SNU = s_i),
    stop(sprintf("Unknown transition latency branch `%s`.", model_row$transition_latency_branch), call. = FALSE)
  )

  X_rem <- switch(
    model_row$remission_branch,
    R0 = cbind(age_s = a_i, sex_num = z_i),
    R1 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i),
    R0S0 = cbind(age_s = a_i, sex_num = z_i),
    R1S0 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i),
    R0S1 = cbind(age_s = a_i, sex_num = z_i, site_SNU = s_i),
    R1S1 = cbind(age_s = a_i, sex_num = z_i, age_s_x_sex = az_i, site_SNU = s_i),
    stop(sprintf("Unknown remission branch `%s`.", model_row$remission_branch), call. = FALSE)
  )

  sd_gamma_lat <- rep(prior_spec$sd_gamma, ncol(X_lat))
  sd_beta_rem <- rep(prior_spec$sd_rem, ncol(X_rem))
  if (grepl("S1$", model_row$transition_latency_branch) && isTRUE(model_row$latency_site_indicator)) {
    sd_gamma_lat[ncol(X_lat)] <- prior_spec$sd_site
  }
  if (grepl("S1$", model_row$remission_branch) && isTRUE(model_row$remission_site_indicator)) {
    sd_beta_rem[ncol(X_rem)] <- prior_spec$sd_site
  }

  rem_cuts <- as.numeric(remission_cut_years)
  if (length(rem_cuts) < 2L || any(diff(rem_cuts) <= 0)) {
    stop("remission_cut_years must be a strictly increasing numeric vector.", call. = FALSE)
  }

  list(
    X_inc = unclass(as.matrix(X_inc)),
    X_lat = unclass(as.matrix(X_lat)),
    X_rem = unclass(as.matrix(X_rem)),
    mu_beta_inc = as.numeric(mu_beta_inc),
    sd_beta_inc = as.numeric(sd_beta_inc),
    alpha_prior_center = as.numeric(prior_spec$alpha_prior_center),
    sd_gamma_lat = as.numeric(sd_gamma_lat),
    sd_beta_rem = as.numeric(sd_beta_rem),
    rem_cuts = rem_cuts,
    time = as.numeric(df$time_year),
    status = as.integer(df$status_num),
    id_df = df %>% select(unique_person_id, id, site, sex_num, age_exact_entry, age_s)
  )
}

# 🔴 Define: competing-risk references and posterior trajectory helpers ===============================
km_step_eval <- function(survfit_obj, times) {
  if (is.null(survfit_obj) || length(survfit_obj$time) == 0L) return(rep(1, length(times)))
  base_times <- survfit_obj$time
  base_surv <- survfit_obj$surv
  vapply(times, function(tt) {
    idx <- max(c(0L, which(base_times <= tt)))
    if (idx == 0L) 1 else base_surv[[idx]]
  }, numeric(1))
}

build_aalen_johansen_reference <- function(df, horizons) {
  df <- normalize_dataset_time(df)
  event_times <- sort(unique(df$time_year[df$status_num %in% c(1L, 2L)]))
  if (length(event_times) == 0L) {
    return(tibble(
      horizon_year = as.integer(horizons),
      observed_transition_cif = 0,
      observed_remission_cif = 0,
      observed_all_event_free = 1
    ))
  }

  n_total <- nrow(df)
  S_all <- 1
  F1 <- 0
  F2 <- 0

  aj_tbl <- tibble(
    time_year = event_times,
    observed_transition_cif = NA_real_,
    observed_remission_cif = NA_real_,
    observed_all_event_free = NA_real_
  )

  n_risk <- n_total
  for (ii in seq_along(event_times)) {
    tt <- event_times[[ii]]
    d1 <- sum(df$status_num == 1L & df$time_year == tt)
    d2 <- sum(df$status_num == 2L & df$time_year == tt)
    c0 <- sum(df$status_num == 0L & df$time_year == tt)

    if (n_risk <= 0L) {
      aj_tbl$observed_transition_cif[[ii]] <- F1
      aj_tbl$observed_remission_cif[[ii]] <- F2
      aj_tbl$observed_all_event_free[[ii]] <- S_all
      next
    }

    F1 <- F1 + S_all * d1 / n_risk
    F2 <- F2 + S_all * d2 / n_risk
    S_all <- S_all * (1 - (d1 + d2) / n_risk)

    aj_tbl$observed_transition_cif[[ii]] <- F1
    aj_tbl$observed_remission_cif[[ii]] <- F2
    aj_tbl$observed_all_event_free[[ii]] <- S_all

    n_risk <- n_risk - d1 - d2 - c0
  }

  eval_step <- function(values, times, horizons_eval, default_value) {
    vapply(horizons_eval, function(hh) {
      idx <- max(c(0L, which(times <= hh)))
      if (idx == 0L) return(default_value)
      values[[idx]]
    }, numeric(1))
  }

  tibble(
    horizon_year = as.integer(horizons),
    observed_transition_cif = eval_step(aj_tbl$observed_transition_cif, aj_tbl$time_year, horizons, 0),
    observed_remission_cif = eval_step(aj_tbl$observed_remission_cif, aj_tbl$time_year, horizons, 0),
    observed_all_event_free = eval_step(aj_tbl$observed_all_event_free, aj_tbl$time_year, horizons, 1)
  )
}

build_stage8b_ipcw_reference <- function(df, horizons) {
  df <- normalize_dataset_time(df)
  censor_fit <- survival::survfit(survival::Surv(time_year, status_num == 0L) ~ 1, data = df)
  aj_ref <- build_aalen_johansen_reference(df, horizons)

  out_rows <- lapply(horizons, function(hh) {
    G_h <- pmax(km_step_eval(censor_fit, hh), 1e-8)
    G_tminus <- pmax(km_step_eval(censor_fit, pmax(df$time_year - 1e-10, 0)), 1e-8)
    G_t <- pmax(km_step_eval(censor_fit, df$time_year), 1e-8)

    w_case <- ifelse(df$status_num == 1L & df$time_year <= hh, 1 / G_tminus, 0)
    w_control <- ifelse(
      df$status_num == 2L & df$time_year <= hh,
      1 / G_t,
      ifelse(df$time_year > hh, 1 / G_h, 0)
    )

    aj_row <- aj_ref %>% filter(horizon_year == hh)

    tibble(
      horizon_year = as.integer(hh),
      observed_transition_cif = aj_row$observed_transition_cif[[1]],
      observed_remission_cif = aj_row$observed_remission_cif[[1]],
      observed_all_event_free = aj_row$observed_all_event_free[[1]],
      denom_case = sum(w_case),
      denom_control = sum(w_control),
      G_h = G_h,
      w_case = list(w_case),
      w_control = list(w_control)
    )
  })

  bind_rows(out_rows)
}

interval_exposure_vector <- function(t, cuts) {
  lower <- cuts[-length(cuts)]
  upper <- cuts[-1]
  pmax(pmin(t, upper) - lower, 0)
}

interval_index_for_time <- function(t, cuts) {
  lower <- cuts[-length(cuts)]
  upper <- cuts[-1]
  idx <- which(t >= lower & t < upper)
  if (length(idx) == 0L) return(length(lower))
  idx[[1L]]
}

compute_stage8b_linear_terms <- function(draws, X_inc, X_lat, X_rem, alpha_prior_center) {
  eta_inc <- draws$beta_inc %*% t(X_inc)
  eta_inc <- sweep(eta_inc, 1, alpha_prior_center + draws$delta0, FUN = "+")
  pi_mat <- plogis(eta_inc)

  mu_lat <- draws$gamma_lat %*% t(X_lat)
  mu_lat <- sweep(mu_lat, 1, draws$gamma0, FUN = "+")
  median_mat <- exp(mu_lat)

  eta_rem <- draws$beta_rem %*% t(X_rem)
  eta_rem <- sweep(eta_rem, 1, draws$xi0, FUN = "+")

  list(
    pi_mat = pmin(pmax(pi_mat, 1e-12), 1 - 1e-12),
    cure_prob_mat = pmin(pmax(1 - pi_mat, 1e-12), 1 - 1e-12),
    mu_lat_mat = mu_lat,
    median_mat = pmin(pmax(median_mat, 1e-8), 1e8),
    eta_rem_mat = eta_rem
  )
}

transition_family_components <- function(t, family_code, mu_lat_mat, median_mat, draws) {
  if (family_code == "E") {
    lambda <- median_mat / log(2)
    Su <- exp(-t / lambda)
    haz <- 1 / lambda
    fu <- haz * Su
  } else if (family_code == "W") {
    k <- exp(draws$rho_W)
    lambda <- median_mat / ((log(2))^(1 / k))
    ratio <- t / lambda
    Su <- exp(-(ratio^k))
    haz <- (k / lambda) * (ratio^(k - 1))
    fu <- haz * Su
  } else if (family_code == "LN") {
    sigma <- exp(draws$log_sigma_LN)
    z <- (log(t) - mu_lat_mat) / sigma
    Su <- 1 - pnorm(z)
    log_pdf <- dnorm(z, log = TRUE) - log(t) - log(sigma)
    fu <- exp(log_pdf)
    haz <- fu / pmax(Su, 1e-12)
  } else if (family_code == "LL") {
    k <- exp(-draws$psi_LL)
    lambda <- median_mat
    ratio <- t / lambda
    Su <- 1 / (1 + ratio^k)
    haz <- (k / lambda) * (ratio^(k - 1)) / (1 + ratio^k)
    fu <- haz * Su
  } else {
    stop(sprintf("Unknown family code `%s`.", family_code), call. = FALSE)
  }

  list(
    Su = pmin(pmax(Su, 1e-12), 1 - 1e-12),
    haz = pmax(haz, 1e-12),
    fu = pmax(fu, 1e-12)
  )
}

remission_components_piecewise <- function(t, eta_rem_mat, log_rho_rem_mat, rem_cuts) {
  expo <- interval_exposure_vector(t, rem_cuts)
  base_draw <- as.vector(exp(log_rho_rem_mat) %*% expo)
  H2 <- sweep(exp(eta_rem_mat), 1, base_draw, FUN = "*")

  interval_idx <- interval_index_for_time(t, rem_cuts)
  rho_t_draw <- exp(log_rho_rem_mat[, interval_idx])
  h2 <- sweep(exp(eta_rem_mat), 1, rho_t_draw, FUN = "*")

  list(
    H2 = pmax(H2, 0),
    G2 = pmin(pmax(exp(-H2), 1e-12), 1),
    h2 = pmax(h2, 1e-12)
  )
}

build_stage8b_prediction_trajectories <- function(state, draws, model_row, horizons, integration_step_year, rem_cuts) {
  max_h <- max(horizons)
  bounds <- seq(0, max_h, by = integration_step_year)
  if (tail(bounds, 1) < max_h) bounds <- c(bounds, max_h)
  mids <- (bounds[-1] + bounds[-length(bounds)]) / 2
  widths <- diff(bounds)

  n_draw <- nrow(state$pi_mat)
  n_subj <- ncol(state$pi_mat)
  F1_cum <- matrix(0, nrow = n_draw, ncol = n_subj)
  F2_cum <- matrix(0, nrow = n_draw, ncol = n_subj)

  horizon_results <- vector("list", length(horizons))
  names(horizon_results) <- as.character(horizons)
  step_to_horizon <- findInterval(horizons, vec = bounds[-1], left.open = FALSE)
  horizon_counter <- 1L

  for (jj in seq_along(mids)) {
    tt <- max(mids[[jj]], 1e-8)
    comp1 <- transition_family_components(tt, model_row$family_code[[1]], state$mu_lat_mat, state$median_mat, draws)
    comp2 <- remission_components_piecewise(tt, state$eta_rem_mat, draws$log_rho_rem, rem_cuts)

    M1 <- (1 - state$pi_mat) + state$pi_mat * comp1$Su
    trans_inc <- state$pi_mat * comp1$fu * comp2$G2 * widths[[jj]]
    rem_inc <- comp2$h2 * comp2$G2 * M1 * widths[[jj]]

    F1_cum <- F1_cum + trans_inc
    F2_cum <- F2_cum + rem_inc

    while (horizon_counter <= length(horizons) && step_to_horizon[[horizon_counter]] == jj) {
      hh <- horizons[[horizon_counter]]
      comp1_h <- transition_family_components(max(hh, 1e-8), model_row$family_code[[1]], state$mu_lat_mat, state$median_mat, draws)
      comp2_h <- remission_components_piecewise(max(hh, 1e-8), state$eta_rem_mat, draws$log_rho_rem, rem_cuts)
      M1_h <- (1 - state$pi_mat) + state$pi_mat * comp1_h$Su
      S_all_h <- comp2_h$G2 * M1_h
      hazard_trans_pop <- (state$pi_mat * comp1_h$fu) / pmax(M1_h, 1e-12)

      horizon_results[[as.character(hh)]] <- list(
        transition_cif = pmin(pmax(F1_cum, 0), 1),
        remission_cif = pmin(pmax(F2_cum, 0), 1),
        all_event_free = pmin(pmax(S_all_h, 0), 1),
        uncured_survival = comp1_h$Su,
        transition_population_hazard = pmax(hazard_trans_pop, 0),
        remission_hazard = comp2_h$h2
      )
      horizon_counter <- horizon_counter + 1L
    }
  }

  horizon_results
}

ppc_horizons_for_dataset <- function(dataset_name) {
  if (dataset_name == "PNU") return(c(1L, 2L))
  c(1L, 2L, 5L)
}

compute_stage8b_degeneracy <- function(pi_mat, median_mat, transition_cif_supported) {
  near_zero_pi <- rowMeans(pi_mat < tiny_susceptible_prob) > degenerate_subject_fraction_threshold
  near_one_pi <- rowMeans(pi_mat > huge_susceptible_prob) > degenerate_subject_fraction_threshold
  near_zero_median <- rowMeans(median_mat < tiny_median_years) > degenerate_subject_fraction_threshold
  huge_median_flat <- if (length(transition_cif_supported) > 0L) {
    mean_risk_supported <- Reduce("+", transition_cif_supported) / length(transition_cif_supported)
    (rowMeans(median_mat > huge_median_years) > degenerate_subject_fraction_threshold) &
      (rowMeans(mean_risk_supported < 0.01) > degenerate_subject_fraction_threshold)
  } else {
    rep(FALSE, nrow(pi_mat))
  }
  any_problem <- near_zero_pi | near_one_pi | near_zero_median | huge_median_flat

  tibble(
    near_zero_pi_rate = mean(near_zero_pi),
    near_one_pi_rate = mean(near_one_pi),
    near_zero_median_rate = mean(near_zero_median),
    huge_median_flat_risk_rate = mean(huge_median_flat),
    degenerate_flag = mean(any_problem) > degenerate_draw_fraction_threshold
  )
}

compute_stage8b_classification_summary <- function(risk_draws, horizon_row, thresholds) {
  prevalence <- as.numeric(horizon_row$observed_transition_cif)
  w_case <- unlist(horizon_row$w_case)
  w_control <- unlist(horizon_row$w_control)
  denom_case <- as.numeric(horizon_row$denom_case)
  denom_control <- as.numeric(horizon_row$denom_control)
  n_subject <- ncol(risk_draws)

  out_list <- vector("list", length(thresholds))
  for (jj in seq_along(thresholds)) {
    thr <- thresholds[[jj]]
    H_mat <- (risk_draws >= thr) * 1

    if (denom_case > 0) tpr_draw <- as.vector(H_mat %*% w_case) / denom_case else tpr_draw <- rep(NA_real_, nrow(risk_draws))
    if (denom_control > 0) fpr_draw <- as.vector(H_mat %*% w_control) / denom_control else fpr_draw <- rep(NA_real_, nrow(risk_draws))

    pos_rate_draw <- prevalence * tpr_draw + (1 - prevalence) * fpr_draw
    ppv_draw <- ifelse(pos_rate_draw > 0, prevalence * tpr_draw / pos_rate_draw, NA_real_)
    fp_burden_draw <- (1 - prevalence) * fpr_draw
    fp100_draw <- 100 * fp_burden_draw
    fp_count_draw <- n_subject * fp_burden_draw
    nb_draw <- prevalence * tpr_draw - (1 - prevalence) * fpr_draw * (thr / (1 - thr))

    out_list[[jj]] <- tibble(
      threshold = thr,
      positive_rate_mean = mean(pos_rate_draw, na.rm = TRUE),
      positive_rate_q025 = safe_quantile(pos_rate_draw, 0.025),
      positive_rate_q50 = safe_quantile(pos_rate_draw, 0.500),
      positive_rate_q975 = safe_quantile(pos_rate_draw, 0.975),
      FPR_mean = mean(fpr_draw, na.rm = TRUE),
      FPR_q025 = safe_quantile(fpr_draw, 0.025),
      FPR_q50 = safe_quantile(fpr_draw, 0.500),
      FPR_q975 = safe_quantile(fpr_draw, 0.975),
      false_positive_burden_mean = mean(fp_burden_draw, na.rm = TRUE),
      false_positive_burden_q025 = safe_quantile(fp_burden_draw, 0.025),
      false_positive_burden_q50 = safe_quantile(fp_burden_draw, 0.500),
      false_positive_burden_q975 = safe_quantile(fp_burden_draw, 0.975),
      false_positive_count_mean = mean(fp_count_draw, na.rm = TRUE),
      false_positive_count_q025 = safe_quantile(fp_count_draw, 0.025),
      false_positive_count_q50 = safe_quantile(fp_count_draw, 0.500),
      false_positive_count_q975 = safe_quantile(fp_count_draw, 0.975),
      FP100_mean = mean(fp100_draw, na.rm = TRUE),
      FP100_q025 = safe_quantile(fp100_draw, 0.025),
      FP100_q50 = safe_quantile(fp100_draw, 0.500),
      FP100_q975 = safe_quantile(fp100_draw, 0.975),
      PPV_mean = mean(ppv_draw, na.rm = TRUE),
      PPV_q025 = safe_quantile(ppv_draw, 0.025),
      PPV_q50 = safe_quantile(ppv_draw, 0.500),
      PPV_q975 = safe_quantile(ppv_draw, 0.975),
      TPR_mean = mean(tpr_draw, na.rm = TRUE),
      TPR_q025 = safe_quantile(tpr_draw, 0.025),
      TPR_q50 = safe_quantile(tpr_draw, 0.500),
      TPR_q975 = safe_quantile(tpr_draw, 0.975),
      NB_mean = mean(nb_draw, na.rm = TRUE),
      NB_q025 = safe_quantile(nb_draw, 0.025),
      NB_q50 = safe_quantile(nb_draw, 0.500),
      NB_q975 = safe_quantile(nb_draw, 0.975)
    )
  }

  bind_rows(out_list)
}

classify_hazard_shape <- function(hazard_values, tol = 1e-6) {
  hv <- as.numeric(hazard_values)
  if (length(hv) < 2L || all(is.na(hv))) return("undetermined")
  diffs <- diff(hv)
  diffs[abs(diffs) < tol] <- 0
  if (all(diffs >= 0)) {
    if (all(diffs == 0)) return("flat")
    return("monotone_increasing")
  }
  if (all(diffs <= 0)) {
    if (all(diffs == 0)) return("flat")
    return("monotone_decreasing")
  }
  peak_idx <- which.max(hv)
  if (peak_idx > 1L && peak_idx < length(hv)) {
    left_ok <- all(diff(hv[1:peak_idx]) >= -tol)
    right_ok <- all(diff(hv[peak_idx:length(hv)]) <= tol)
    if (left_ok && right_ok) return("unimodal")
  }
  "irregular"
}

# 🔴 Define: prior predictive, anchor update, Stan model, diagnostics ===============================
get_stage8b_site_prior_indices <- function(model_row, design_bundle) {
  list(
    inc_idx = if (isTRUE(model_row$incidence_site_indicator[[1]])) ncol(design_bundle$X_inc) else NA_integer_,
    lat_idx = if (isTRUE(model_row$latency_site_indicator[[1]])) ncol(design_bundle$X_lat) else NA_integer_,
    rem_idx = if (isTRUE(model_row$remission_site_indicator[[1]])) ncol(design_bundle$X_rem) else NA_integer_
  )
}

apply_stage8b_site_prior_draws <- function(mat, idx, site_prior_family, scale) {
  if (is.na(idx) || idx <= 0L || idx > ncol(mat)) return(mat)
  if (identical(site_prior_family, "student_t3_0_1_sensitivity")) {
    mat[, idx] <- stats::rt(nrow(mat), df = 3) * scale
  }
  mat
}

empty_stage8b_prior_predictive_summary <- function() {
  tibble(
    dataset_key = character(),
    branch = character(),
    model_id = character(),
    retained_fit_id = character(),
    risk_scale = character(),
    prior_branch = character(),
    site_prior_family = character(),
    metric = character(),
    horizon = integer(),
    mean = double(),
    sd = double(),
    q025 = double(),
    q50 = double(),
    q975 = double(),
    prior_degenerate_flag = logical()
  )
}

simulate_stage8b_prior_predictive <- function(dataset_df, model_row, design_bundle, prior_spec, n_draws, horizons_eval) {
  idx_info <- get_stage8b_site_prior_indices(model_row, design_bundle)
  K_inc <- ncol(design_bundle$X_inc)
  K_lat <- ncol(design_bundle$X_lat)
  K_rem <- ncol(design_bundle$X_rem)
  J_rem <- length(design_bundle$rem_cuts) - 1L

  sim_draws <- list(
    delta0 = rnorm(n_draws, 0, prior_spec$sd_delta),
    beta_inc = matrix(
      rnorm(
        n_draws * K_inc,
        mean = rep(design_bundle$mu_beta_inc, each = n_draws),
        sd = rep(design_bundle$sd_beta_inc, each = n_draws)
      ),
      nrow = n_draws,
      byrow = FALSE
    ),
    gamma0 = rnorm(n_draws, 0, prior_spec$sd_gamma0),
    gamma_lat = matrix(
      rnorm(n_draws * K_lat, mean = 0, sd = rep(design_bundle$sd_gamma_lat, each = n_draws)),
      nrow = n_draws,
      byrow = FALSE
    ),
    xi0 = rnorm(n_draws, 0, prior_spec$sd_xi0),
    beta_rem = matrix(
      rnorm(n_draws * K_rem, mean = 0, sd = rep(design_bundle$sd_beta_rem, each = n_draws)),
      nrow = n_draws,
      byrow = FALSE
    ),
    log_rho_rem = matrix(
      rnorm(n_draws * J_rem, mean = 0, sd = prior_spec$sd_log_rho_rem),
      nrow = n_draws,
      byrow = FALSE
    ),
    rho_W = rnorm(n_draws, 0, prior_spec$sd_rho_W),
    log_sigma_LN = rnorm(n_draws, 0, prior_spec$sd_log_sigma_LN),
    psi_LL = rnorm(n_draws, 0, prior_spec$sd_psi_LL)
  )

  sim_draws$beta_inc <- apply_stage8b_site_prior_draws(sim_draws$beta_inc, idx_info$inc_idx, model_row$site_prior_family[[1]], prior_spec$sd_site)
  sim_draws$gamma_lat <- apply_stage8b_site_prior_draws(sim_draws$gamma_lat, idx_info$lat_idx, model_row$site_prior_family[[1]], prior_spec$sd_site)
  sim_draws$beta_rem <- apply_stage8b_site_prior_draws(sim_draws$beta_rem, idx_info$rem_idx, model_row$site_prior_family[[1]], prior_spec$sd_site)

  state <- compute_stage8b_linear_terms(
    draws = sim_draws,
    X_inc = design_bundle$X_inc,
    X_lat = design_bundle$X_lat,
    X_rem = design_bundle$X_rem,
    alpha_prior_center = design_bundle$alpha_prior_center
  )

  pred <- build_stage8b_prediction_trajectories(
    state = state,
    draws = sim_draws,
    model_row = model_row,
    horizons = sort(unique(as.integer(horizons_eval))),
    integration_step_year = integration_step_year,
    rem_cuts = design_bundle$rem_cuts
  )

  supported_horizons <- intersect(as.character(ppc_horizons_for_dataset(model_row$dataset[[1]])), names(pred))
  transition_supported <- lapply(supported_horizons, function(hh) pred[[hh]]$transition_cif)
  degeneracy <- compute_stage8b_degeneracy(state$pi_mat, state$median_mat, transition_supported)

  out_rows <- list(
    tibble(dataset_key = model_row$dataset[[1]], branch = "Stage8B", model_id = model_row$model_id[[1]], retained_fit_id = model_row$retained_fit_id[[1]], risk_scale = model_row$risk_scale[[1]], prior_branch = model_row$prior_branch[[1]], site_prior_family = model_row$site_prior_family[[1]], metric = "susceptible_fraction", horizon = NA_integer_, summary_scalar(rowMeans(state$pi_mat)), prior_degenerate_flag = degeneracy$degenerate_flag[[1]]),
    tibble(dataset_key = model_row$dataset[[1]], branch = "Stage8B", model_id = model_row$model_id[[1]], retained_fit_id = model_row$retained_fit_id[[1]], risk_scale = model_row$risk_scale[[1]], prior_branch = model_row$prior_branch[[1]], site_prior_family = model_row$site_prior_family[[1]], metric = "cure_fraction", horizon = NA_integer_, summary_scalar(rowMeans(1 - state$pi_mat)), prior_degenerate_flag = degeneracy$degenerate_flag[[1]]),
    tibble(dataset_key = model_row$dataset[[1]], branch = "Stage8B", model_id = model_row$model_id[[1]], retained_fit_id = model_row$retained_fit_id[[1]], risk_scale = model_row$risk_scale[[1]], prior_branch = model_row$prior_branch[[1]], site_prior_family = model_row$site_prior_family[[1]], metric = "median_susceptible_time", horizon = NA_integer_, summary_scalar(apply(state$median_mat, 1, stats::median)), prior_degenerate_flag = degeneracy$degenerate_flag[[1]])
  )

  for (hh in sort(unique(as.integer(horizons_eval)))) {
    pred_h <- pred[[as.character(hh)]]
    out_rows[[length(out_rows) + 1L]] <- tibble(dataset_key = model_row$dataset[[1]], branch = "Stage8B", model_id = model_row$model_id[[1]], retained_fit_id = model_row$retained_fit_id[[1]], risk_scale = model_row$risk_scale[[1]], prior_branch = model_row$prior_branch[[1]], site_prior_family = model_row$site_prior_family[[1]], metric = "transition_cif", horizon = hh, summary_scalar(rowMeans(pred_h$transition_cif)), prior_degenerate_flag = degeneracy$degenerate_flag[[1]])
    out_rows[[length(out_rows) + 1L]] <- tibble(dataset_key = model_row$dataset[[1]], branch = "Stage8B", model_id = model_row$model_id[[1]], retained_fit_id = model_row$retained_fit_id[[1]], risk_scale = model_row$risk_scale[[1]], prior_branch = model_row$prior_branch[[1]], site_prior_family = model_row$site_prior_family[[1]], metric = "remission_cif", horizon = hh, summary_scalar(rowMeans(pred_h$remission_cif)), prior_degenerate_flag = degeneracy$degenerate_flag[[1]])
    out_rows[[length(out_rows) + 1L]] <- tibble(dataset_key = model_row$dataset[[1]], branch = "Stage8B", model_id = model_row$model_id[[1]], retained_fit_id = model_row$retained_fit_id[[1]], risk_scale = model_row$risk_scale[[1]], prior_branch = model_row$prior_branch[[1]], site_prior_family = model_row$site_prior_family[[1]], metric = "all_event_free", horizon = hh, summary_scalar(rowMeans(pred_h$all_event_free)), prior_degenerate_flag = degeneracy$degenerate_flag[[1]])
  }

  bind_rows(out_rows)
}

annotate_stage8b_prior_predictive <- function(df) {
  if (nrow_or_zero(df) == 0L) {
    return(empty_stage8b_prior_predictive_summary() %>% mutate(prior_tail_warning_flag = logical(), prior_tail_warning_detail = character()))
  }

  tibble::as_tibble(df) %>%
    mutate(
      horizon = as.integer(safe_numeric(horizon)),
      mean = safe_numeric(mean),
      q025 = safe_numeric(q025),
      q50 = safe_numeric(q50),
      q975 = safe_numeric(q975),
      prior_degenerate_flag = as.logical(prior_degenerate_flag),
      prior_tail_warning_flag = metric == "median_susceptible_time" & ((mean > prior_tail_warning_mean_years) | (q975 > prior_tail_warning_q975_years)),
      prior_tail_warning_detail = case_when(
        metric == "median_susceptible_time" & mean > prior_tail_warning_mean_years & q975 > prior_tail_warning_q975_years ~ paste0("Prior median susceptible time mean > ", prior_tail_warning_mean_years, " years and q975 > ", prior_tail_warning_q975_years, " years."),
        metric == "median_susceptible_time" & mean > prior_tail_warning_mean_years ~ paste0("Prior median susceptible time mean > ", prior_tail_warning_mean_years, " years."),
        metric == "median_susceptible_time" & q975 > prior_tail_warning_q975_years ~ paste0("Prior median susceptible time q975 > ", prior_tail_warning_q975_years, " years."),
        TRUE ~ NA_character_
      )
    )
}

make_stage8b_anchor_cells <- function(model_row) {
  age_band_df <- tibble::tribble(
    ~sex_num, ~age20_29, ~age30plus, ~cell_label,
    0L, 0L, 0L, "Female_lt20",
    1L, 0L, 0L, "Male_lt20",
    0L, 1L, 0L, "Female_20_29",
    1L, 1L, 0L, "Male_20_29",
    0L, 0L, 1L, "Female_30plus",
    1L, 0L, 1L, "Male_30plus"
  )

  site_levels <- if (model_row$dataset[[1]] == "merged" && isTRUE(model_row$incidence_site_indicator[[1]])) {
    c("PNU", "SNU")
  } else if (model_row$dataset[[1]] == "merged") {
    "pooled"
  } else {
    model_row$dataset[[1]]
  }

  tidyr::crossing(age_band_df, tibble(site_level = site_levels)) %>%
    mutate(
      sex_x_age20_29 = sex_num * age20_29,
      sex_x_age30plus = sex_num * age30plus,
      site_SNU = as.integer(site_level == "SNU"),
      age_sex_anchor_cell = paste0(cell_label, "__", site_level)
    )
}

make_stage8b_incidence_anchor_update <- function(model_row, design_bundle, prior_spec, draws_compact) {
  if (!identical(model_row$prior_branch[[1]], "anchor_informed")) return(tibble())

  cells <- make_stage8b_anchor_cells(model_row)
  out_rows <- vector("list", nrow(cells))

  for (ii in seq_len(nrow(cells))) {
    one <- cells[ii, , drop = FALSE]
    x_vec <- c(
      sex_num = one$sex_num[[1]],
      age20_29 = one$age20_29[[1]],
      age30plus = one$age30plus[[1]],
      sex_x_age20_29 = one$sex_x_age20_29[[1]],
      sex_x_age30plus = one$sex_x_age30plus[[1]]
    )
    if (isTRUE(model_row$incidence_site_indicator[[1]])) {
      x_vec <- c(x_vec, site_SNU = one$site_SNU[[1]])
    }

    prior_center_logit <- design_bundle$alpha_prior_center + sum(design_bundle$mu_beta_inc * x_vec)
    external_one_year_risk <- plogis(prior_center_logit)
    external_incidence_rate_per10k <- -10000 * log(pmax(1 - external_one_year_risk, 1e-12))

    post_logit <- design_bundle$alpha_prior_center + draws_compact$delta0 + as.vector(draws_compact$beta_inc %*% x_vec)
    post_risk <- plogis(post_logit)

    out_rows[[ii]] <- tibble(
      dataset_key = model_row$dataset[[1]],
      branch = "Stage8B",
      risk_scale = model_row$risk_scale[[1]],
      prior_branch = model_row$prior_branch[[1]],
      retained_fit_id = model_row$retained_fit_id[[1]],
      model_id = model_row$model_id[[1]],
      site_prior_family = model_row$site_prior_family[[1]],
      age_sex_anchor_cell = one$age_sex_anchor_cell[[1]],
      site_level = one$site_level[[1]],
      external_incidence_rate_per10k = external_incidence_rate_per10k,
      external_one_year_risk = external_one_year_risk,
      prior_center_logit = prior_center_logit,
      posterior_mean_logit = mean(post_logit),
      posterior_lower_logit = safe_quantile(post_logit, 0.025),
      posterior_upper_logit = safe_quantile(post_logit, 0.975),
      posterior_mean_one_year_risk = mean(post_risk),
      posterior_lower_one_year_risk = safe_quantile(post_risk, 0.025),
      posterior_upper_one_year_risk = safe_quantile(post_risk, 0.975),
      posterior_minus_prior_logit = mean(post_logit) - prior_center_logit,
      posterior_minus_prior_risk = mean(post_risk) - external_one_year_risk
    )
  }

  bind_rows(out_rows)
}

compile_stage8b_stan_model <- function() {
  stan_code <- r"(
functions {
  real loglogistic_lpdf_custom(real t, real lambda, real k) {
    real ratio;
    ratio = t / lambda;
    return log(k) - log(lambda) + (k - 1) * log(ratio) - 2 * log1p(pow(ratio, k));
  }
  real loglogistic_lccdf_custom(real t, real lambda, real k) {
    return -log1p(pow(t / lambda, k));
  }
  real rem_base_hazard(real t, vector cuts, vector log_rho_rem) {
    int J;
    real out;
    J = num_elements(log_rho_rem);
    out = 0;
    for (j in 1:J) {
      real lower;
      real upper;
      real exposure;
      lower = cuts[j];
      upper = cuts[j + 1];
      exposure = fmax(fmin(t, upper) - lower, 0);
      out += exp(log_rho_rem[j]) * exposure;
    }
    return out;
  }
  int rem_interval_index(real t, vector cuts) {
    int J;
    J = num_elements(cuts) - 1;
    for (j in 1:J) {
      if (t >= cuts[j] && t < cuts[j + 1]) {
        return j;
      }
    }
    return J;
  }
}
data {
  int<lower=1> N;
  vector<lower=1e-8>[N] time;
  int<lower=0, upper=2> status[N];
  int<lower=1> K_inc;
  matrix[N, K_inc] X_inc;
  int<lower=1> K_lat;
  matrix[N, K_lat] X_lat;
  int<lower=1> K_rem;
  matrix[N, K_rem] X_rem;
  int<lower=1, upper=4> family_id;
  real alpha_prior_center;
  vector[K_inc] mu_beta_inc;
  vector<lower=0>[K_inc] sd_beta_inc;
  real<lower=0> sd_delta;
  real<lower=0> sd_gamma0;
  vector<lower=0>[K_lat] sd_gamma_lat;
  real<lower=0> sd_xi0;
  vector<lower=0>[K_rem] sd_beta_rem;
  real<lower=0> sd_log_rho_rem;
  real<lower=0> sd_rho_W;
  real<lower=0> sd_log_sigma_LN;
  real<lower=0> sd_psi_LL;
  int<lower=1> J_rem;
  vector[J_rem + 1] rem_cuts;
  real<lower=0> sd_site_prior;
  int<lower=0> inc_site_index;
  int<lower=0> lat_site_index;
  int<lower=0> rem_site_index;
  int<lower=0, upper=1> use_t_prior_inc_site;
  int<lower=0, upper=1> use_t_prior_lat_site;
  int<lower=0, upper=1> use_t_prior_rem_site;
}
parameters {
  real delta0;
  vector[K_inc] beta_inc;
  real gamma0;
  vector[K_lat] gamma_lat;
  real xi0;
  vector[K_rem] beta_rem;
  vector[J_rem] log_rho_rem;
  real rho_W;
  real log_sigma_LN;
  real psi_LL;
}
model {
  delta0 ~ normal(0, sd_delta);
  gamma0 ~ normal(0, sd_gamma0);
  xi0 ~ normal(0, sd_xi0);
  log_rho_rem ~ normal(0, sd_log_rho_rem);
  rho_W ~ normal(0, sd_rho_W);
  log_sigma_LN ~ normal(0, sd_log_sigma_LN);
  psi_LL ~ normal(0, sd_psi_LL);

  for (j in 1:K_inc) {
    if (use_t_prior_inc_site == 1 && inc_site_index > 0 && j == inc_site_index) {
      target += student_t_lpdf(beta_inc[j] | 3, 0, sd_site_prior);
    } else {
      target += normal_lpdf(beta_inc[j] | mu_beta_inc[j], sd_beta_inc[j]);
    }
  }

  for (j in 1:K_lat) {
    if (use_t_prior_lat_site == 1 && lat_site_index > 0 && j == lat_site_index) {
      target += student_t_lpdf(gamma_lat[j] | 3, 0, sd_site_prior);
    } else {
      target += normal_lpdf(gamma_lat[j] | 0, sd_gamma_lat[j]);
    }
  }

  for (j in 1:K_rem) {
    if (use_t_prior_rem_site == 1 && rem_site_index > 0 && j == rem_site_index) {
      target += student_t_lpdf(beta_rem[j] | 3, 0, sd_site_prior);
    } else {
      target += normal_lpdf(beta_rem[j] | 0, sd_beta_rem[j]);
    }
  }

  for (i in 1:N) {
    real eta_inc;
    real pi_i;
    real mu_lat;
    real m_i;
    real logS1;
    real logf1;
    real eta_rem;
    real H2_base;
    real logG2;
    real h2;
    real logM1;

    eta_inc = alpha_prior_center + delta0 + dot_product(row(X_inc, i), beta_inc);
    pi_i = inv_logit(eta_inc);
    mu_lat = gamma0 + dot_product(row(X_lat, i), gamma_lat);
    m_i = exp(mu_lat);
    eta_rem = xi0 + dot_product(row(X_rem, i), beta_rem);
    H2_base = rem_base_hazard(time[i], rem_cuts, log_rho_rem);
    logG2 = -exp(eta_rem) * H2_base;
    h2 = exp(eta_rem + log_rho_rem[rem_interval_index(time[i], rem_cuts)]);

    if (family_id == 1) {
      real lambda;
      lambda = m_i / log(2);
      logS1 = -time[i] / lambda;
      logf1 = exponential_lpdf(time[i] | 1 / lambda);
    } else if (family_id == 2) {
      real kW;
      real lambda;
      kW = exp(rho_W);
      lambda = m_i / pow(log(2), 1 / kW);
      logS1 = -pow(time[i] / lambda, kW);
      logf1 = weibull_lpdf(time[i] | kW, lambda);
    } else if (family_id == 3) {
      real sigmaLN;
      sigmaLN = exp(log_sigma_LN);
      logS1 = normal_lccdf(log(time[i]) | mu_lat, sigmaLN);
      logf1 = lognormal_lpdf(time[i] | mu_lat, sigmaLN);
    } else {
      real kLL;
      real lambda;
      kLL = exp(-psi_LL);
      lambda = m_i;
      logS1 = loglogistic_lccdf_custom(time[i], lambda, kLL);
      logf1 = loglogistic_lpdf_custom(time[i], lambda, kLL);
    }

    logM1 = log_sum_exp(log1m(pi_i), log(pi_i) + logS1);

    if (status[i] == 1) {
      target += log(pi_i) + logf1 + logG2;
    } else if (status[i] == 2) {
      target += log(h2) + logG2 + logM1;
    } else {
      target += logG2 + logM1;
    }
  }
}
generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    real eta_inc;
    real pi_i;
    real mu_lat;
    real m_i;
    real logS1;
    real logf1;
    real eta_rem;
    real H2_base;
    real logG2;
    real h2;
    real logM1;

    eta_inc = alpha_prior_center + delta0 + dot_product(row(X_inc, i), beta_inc);
    pi_i = inv_logit(eta_inc);
    mu_lat = gamma0 + dot_product(row(X_lat, i), gamma_lat);
    m_i = exp(mu_lat);
    eta_rem = xi0 + dot_product(row(X_rem, i), beta_rem);
    H2_base = rem_base_hazard(time[i], rem_cuts, log_rho_rem);
    logG2 = -exp(eta_rem) * H2_base;
    h2 = exp(eta_rem + log_rho_rem[rem_interval_index(time[i], rem_cuts)]);

    if (family_id == 1) {
      real lambda;
      lambda = m_i / log(2);
      logS1 = -time[i] / lambda;
      logf1 = exponential_lpdf(time[i] | 1 / lambda);
    } else if (family_id == 2) {
      real kW;
      real lambda;
      kW = exp(rho_W);
      lambda = m_i / pow(log(2), 1 / kW);
      logS1 = -pow(time[i] / lambda, kW);
      logf1 = weibull_lpdf(time[i] | kW, lambda);
    } else if (family_id == 3) {
      real sigmaLN;
      sigmaLN = exp(log_sigma_LN);
      logS1 = normal_lccdf(log(time[i]) | mu_lat, sigmaLN);
      logf1 = lognormal_lpdf(time[i] | mu_lat, sigmaLN);
    } else {
      real kLL;
      real lambda;
      kLL = exp(-psi_LL);
      lambda = m_i;
      logS1 = loglogistic_lccdf_custom(time[i], lambda, kLL);
      logf1 = loglogistic_lpdf_custom(time[i], lambda, kLL);
    }

    logM1 = log_sum_exp(log1m(pi_i), log(pi_i) + logS1);

    if (status[i] == 1) {
      log_lik[i] = log(pi_i) + logf1 + logG2;
    } else if (status[i] == 2) {
      log_lik[i] = log(h2) + logG2 + logM1;
    } else {
      log_lik[i] = logG2 + logM1;
    }
  }
}
)"
  rstan::stan_model(model_code = stan_code, model_name = "stage8b_bayesian_competing_risk_cure")
}

extract_stage8b_draws_compact <- function(fit, K_inc, K_lat, K_rem, J_rem) {
  ext <- rstan::extract(
    fit,
    pars = c("delta0", "beta_inc", "gamma0", "gamma_lat", "xi0", "beta_rem", "log_rho_rem", "rho_W", "log_sigma_LN", "psi_LL", "log_lik"),
    permuted = TRUE,
    inc_warmup = FALSE
  )

  beta_inc <- ext$beta_inc
  gamma_lat <- ext$gamma_lat
  beta_rem <- ext$beta_rem
  log_rho_rem <- ext$log_rho_rem
  log_lik <- ext$log_lik

  if (is.null(dim(beta_inc))) beta_inc <- matrix(beta_inc, ncol = K_inc)
  if (is.null(dim(gamma_lat))) gamma_lat <- matrix(gamma_lat, ncol = K_lat)
  if (is.null(dim(beta_rem))) beta_rem <- matrix(beta_rem, ncol = K_rem)
  if (is.null(dim(log_rho_rem))) log_rho_rem <- matrix(log_rho_rem, ncol = J_rem)
  if (is.null(dim(log_lik))) log_lik <- matrix(log_lik, nrow = length(ext$delta0))

  list(
    delta0 = as.numeric(ext$delta0),
    beta_inc = beta_inc,
    gamma0 = as.numeric(ext$gamma0),
    gamma_lat = gamma_lat,
    xi0 = as.numeric(ext$xi0),
    beta_rem = beta_rem,
    log_rho_rem = log_rho_rem,
    rho_W = as.numeric(ext$rho_W),
    log_sigma_LN = as.numeric(ext$log_sigma_LN),
    psi_LL = as.numeric(ext$psi_LL),
    log_lik = log_lik
  )
}

compute_information_criteria <- function(log_lik) {
  out <- list(
    waic = NA_real_,
    p_waic = NA_real_,
    p_waic_high_n = NA_integer_,
    p_waic_high_pct = NA_real_,
    looic = NA_real_,
    p_loo = NA_real_,
    pareto_k_max = NA_real_,
    pareto_k_bad_n = NA_integer_,
    pareto_k_bad_pct = NA_real_,
    pareto_k_very_bad_n = NA_integer_,
    waic_warning_flag = FALSE,
    loo_warning_flag = FALSE,
    info_criteria_warning_detail = NA_character_
  )

  if (is.null(log_lik) || !requireNamespace("loo", quietly = TRUE)) return(out)

  warning_texts <- character()
  collect_warning <- function(w) {
    warning_texts <<- c(warning_texts, conditionMessage(w))
    tryInvokeRestart("muffleWarning")
  }

  waic_obj <- withCallingHandlers(tryCatch(loo::waic(log_lik), error = function(e) e), warning = collect_warning)
  if (!inherits(waic_obj, "error")) {
    if ("waic" %in% rownames(waic_obj$estimates)) out$waic <- as.numeric(waic_obj$estimates["waic", "Estimate"])
    if ("p_waic" %in% rownames(waic_obj$estimates)) out$p_waic <- as.numeric(waic_obj$estimates["p_waic", "Estimate"])
    if (!is.null(waic_obj$pointwise) && "p_waic" %in% colnames(waic_obj$pointwise)) {
      p_waic_pointwise <- as.numeric(waic_obj$pointwise[, "p_waic"])
      out$p_waic_high_n <- as.integer(sum(p_waic_pointwise > 0.4, na.rm = TRUE))
      out$p_waic_high_pct <- if (length(p_waic_pointwise) > 0L) 100 * out$p_waic_high_n / length(p_waic_pointwise) else NA_real_
    }
  } else {
    warning_texts <- c(warning_texts, paste0("waic_error: ", conditionMessage(waic_obj)))
  }

  loo_obj <- withCallingHandlers(tryCatch(loo::loo(log_lik), error = function(e) e), warning = collect_warning)
  if (!inherits(loo_obj, "error")) {
    if ("looic" %in% rownames(loo_obj$estimates)) out$looic <- as.numeric(loo_obj$estimates["looic", "Estimate"])
    if ("p_loo" %in% rownames(loo_obj$estimates)) out$p_loo <- as.numeric(loo_obj$estimates["p_loo", "Estimate"])
    if (!is.null(loo_obj$diagnostics$pareto_k)) {
      pareto_k <- as.numeric(loo_obj$diagnostics$pareto_k)
      out$pareto_k_max <- if (length(pareto_k) > 0L && any(is.finite(pareto_k))) max(pareto_k, na.rm = TRUE) else NA_real_
      out$pareto_k_bad_n <- as.integer(sum(pareto_k > 0.7, na.rm = TRUE))
      out$pareto_k_very_bad_n <- as.integer(sum(pareto_k > 1.0, na.rm = TRUE))
      out$pareto_k_bad_pct <- if (length(pareto_k) > 0L) 100 * out$pareto_k_bad_n / length(pareto_k) else NA_real_
    }
  } else {
    warning_texts <- c(warning_texts, paste0("loo_error: ", conditionMessage(loo_obj)))
  }

  out$waic_warning_flag <- any(grepl("p_waic", warning_texts, fixed = TRUE))
  out$loo_warning_flag <- any(grepl("Pareto k", warning_texts, fixed = TRUE))
  if (length(warning_texts) > 0L) out$info_criteria_warning_detail <- paste(unique(warning_texts), collapse = " | ")
  out
}

bundle_get <- function(bundle, name) {
  if (is.null(bundle) || !is.list(bundle)) return(NULL)
  if (!is.null(bundle[[name]])) return(bundle[[name]])
  outputs <- bundle$outputs %||% list()
  outputs[[name]] %||% NULL
}

load_stage8b_reuse_bundle <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

is_stage8b_bundle_reusable <- function(bundle, model_row, dataset_df) {
  if (is.null(bundle) || !is.list(bundle) || is.null(bundle$model_registry_row)) return(FALSE)
  reg <- tibble::as_tibble(bundle$model_registry_row)
  if (nrow(reg) != 1L) return(FALSE)
  identical(as.character(reg$model_id[[1]]), as.character(model_row$model_id[[1]])) &&
    identical(as.character(reg$dataset_key[[1]]), as.character(model_row$dataset[[1]])) &&
    identical(as.character(reg$fit_status[[1]]), "ok") &&
    identical(as.integer(reg$n[[1]]), nrow(dataset_df))
}

# 🔴 Define: performance summaries and figure-source helpers ===============================
binary_auc_complete_case <- function(score, label) {
  ok <- is.finite(score) & !is.na(label)
  score <- score[ok]
  label <- label[ok]
  if (length(unique(label)) < 2L) return(NA_real_)
  n1 <- sum(label == 1)
  n0 <- sum(label == 0)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

compute_brier_ipcw_simple <- function(df_h) {
  h <- unique(df_h$horizon)
  if (length(h) != 1L) stop("df_h must contain one horizon only")
  h <- as.numeric(h[[1]])

  censor_fit <- survival::survfit(survival::Surv(time_year, status_num == 0L) ~ 1, data = df_h)
  G_h <- pmax(km_step_eval(censor_fit, h), 1e-8)
  G_tminus <- pmax(km_step_eval(censor_fit, pmax(df_h$time_year - 1e-10, 0)), 1e-8)

  y_h <- as.integer(df_h$status_num == 1L & df_h$time_year <= h)
  w <- ifelse(
    df_h$status_num == 1L & df_h$time_year <= h, 1 / G_tminus,
    ifelse(df_h$status_num == 2L & df_h$time_year <= h, 1 / G_tminus,
           ifelse(df_h$time_year > h, 1 / G_h, 0))
  )

  p <- pmin(pmax(df_h$transition_cif_mean, 0), 1)
  tibble(
    brier_ipcw_mean = weighted.mean((y_h - p)^2, w = w, na.rm = TRUE),
    brier_ipcw_weight_sum = sum(w, na.rm = TRUE)
  )
}

compute_discrimination_simple <- function(df_h) {
  h <- unique(df_h$horizon)
  if (length(h) != 1L) stop("df_h must contain one horizon only")
  h <- as.numeric(h[[1]])

  label <- ifelse(df_h$status_num == 1L & df_h$time_year <= h, 1L,
                  ifelse((df_h$status_num == 2L & df_h$time_year <= h) | (df_h$time_year > h), 0L, NA_integer_))

  tibble(
    auc_transition_horizon = binary_auc_complete_case(df_h$transition_cif_mean, label),
    discrimination_n_complete = sum(!is.na(label))
  )
}

compute_calibration_simple <- function(df_h, observed_transition_cif) {
  pred_mean <- mean(df_h$transition_cif_mean, na.rm = TRUE)
  pred_median <- stats::median(df_h$transition_cif_mean, na.rm = TRUE)
  cal_diff <- pred_mean - observed_transition_cif
  tibble(
    calibration_mean_pred = pred_mean,
    calibration_median_pred = pred_median,
    calibration_diff_mean_minus_observed = cal_diff,
    calibration_abs_diff = abs(cal_diff),
    calibration_oe_ratio = ifelse(is.finite(observed_transition_cif) && observed_transition_cif > 0, pred_mean / observed_transition_cif, NA_real_)
  )
}

pivot_prediction_table <- function(df, id_cols, metric_cols) {
  if (nrow(df) == 0L || length(metric_cols) == 0L) return(tibble())
  df %>%
    select(any_of(id_cols), any_of(metric_cols)) %>%
    tidyr::pivot_longer(
      cols = any_of(metric_cols),
      names_to = c("metric", "summary_stat"),
      names_pattern = "^(.*)_(mean|q025|q50|q975)$",
      values_to = "value"
    )
}

metric_to_anchor_delta_field <- function(metric) {
  dplyr::case_when(
    metric == "transition_cif" ~ "delta_risk_anchor_minus_neutral",
    metric == "cure_fraction" ~ "delta_cure_fraction_anchor_minus_neutral",
    metric == "false_positive_burden" ~ "delta_false_positive_burden_anchor_minus_neutral",
    metric == "FP100" ~ "delta_FP100_anchor_minus_neutral",
    metric == "NB" ~ "delta_NB_anchor_minus_neutral",
    metric == "PPV" ~ "delta_PPV_anchor_minus_neutral",
    metric == "TPR" ~ "delta_TPR_anchor_minus_neutral",
    TRUE ~ NA_character_
  )
}

metric_to_stage8a_delta_field <- function(metric) {
  dplyr::case_when(
    metric == "transition_cif" ~ "delta_risk_8B_minus_8A",
    metric == "cure_fraction" ~ "delta_cure_fraction_8B_minus_8A",
    metric == "false_positive_burden" ~ "delta_false_positive_burden_8B_minus_8A",
    metric == "FP100" ~ "delta_FP100_8B_minus_8A",
    metric == "NB" ~ "delta_NB_8B_minus_8A",
    metric == "PPV" ~ "delta_PPV_8B_minus_8A",
    metric == "TPR" ~ "delta_TPR_8B_minus_8A",
    TRUE ~ NA_character_
  )
}

safe_save_png <- function(plot_obj, path, width = 9, height = 6, dpi = 300) {
  if (is.null(plot_obj)) return(invisible(FALSE))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- make_temp_output_path(path, tag = "tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  ggplot2::ggsave(filename = tmp, plot = plot_obj, width = width, height = height, dpi = dpi)
  safe_promote_file(tmp, path)
  invisible(TRUE)
}

safe_generate_stage8b_visuals <- function(posterior_cohort_yearly, posterior_classification, ppc_summary, final_pdf_path,
                                          transition_png_path, remission_png_path, nb_png_path, ppc_png_path) {
  dir.create(dirname(final_pdf_path), recursive = TRUE, showWarnings = FALSE)
  tmp_pdf <- make_temp_output_path(final_pdf_path, tag = "tmp")
  pdf_open <- FALSE

  close_pdf <- function() {
    if (isTRUE(pdf_open)) {
      try(grDevices::dev.off(), silent = TRUE)
      pdf_open <<- FALSE
    }
  }

  on.exit({
    close_pdf()
    if (file.exists(tmp_pdf)) unlink(tmp_pdf)
  }, add = TRUE)

  transition_plot <- NULL
  remission_plot <- NULL
  nb_plot <- NULL
  ppc_plot <- NULL

  if (nrow_or_zero(posterior_cohort_yearly) > 0L) {
    transition_plot <- posterior_cohort_yearly %>%
      ggplot(aes(x = horizon, y = transition_cif_mean, color = model_id, fill = model_id)) +
      geom_ribbon(aes(ymin = transition_cif_q025, ymax = transition_cif_q975), alpha = 0.15, linewidth = 0) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~ dataset_key, scales = "free_y") +
      labs(title = "Stage 8B posterior transition CIF trajectories", x = "Horizon (years)", y = "Transition CIF") +
      theme_bw() + theme(legend.position = "none")

    remission_plot <- posterior_cohort_yearly %>%
      ggplot(aes(x = horizon, y = remission_cif_mean, color = model_id, fill = model_id)) +
      geom_ribbon(aes(ymin = remission_cif_q025, ymax = remission_cif_q975), alpha = 0.15, linewidth = 0) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~ dataset_key, scales = "free_y") +
      labs(title = "Stage 8B posterior remission CIF trajectories", x = "Horizon (years)", y = "Remission CIF") +
      theme_bw() + theme(legend.position = "none")
  }

  if (nrow_or_zero(posterior_classification) > 0L) {
    nb_df <- posterior_classification %>% filter(horizon %in% c(1L, 2L, 5L))
    if (nrow(nb_df) > 0L) {
      nb_plot <- nb_df %>%
        ggplot(aes(x = threshold, y = NB_mean, color = model_id, fill = model_id)) +
        geom_ribbon(aes(ymin = NB_q025, ymax = NB_q975), alpha = 0.15, linewidth = 0) +
        geom_line(linewidth = 0.7) +
        facet_grid(dataset_key ~ horizon, scales = "free_y") +
        labs(title = "Stage 8B net benefit by threshold", x = "Threshold", y = "Net benefit") +
        theme_bw() + theme(legend.position = "none")
    }
  }

  if (nrow_or_zero(ppc_summary) > 0L) {
    ppc_plot <- ppc_summary %>%
      ggplot(aes(x = horizon, y = posterior_mean_transition_cif, color = model_id)) +
      geom_errorbar(aes(ymin = posterior_q025_transition_cif, ymax = posterior_q975_transition_cif), width = 0.12, alpha = 0.6) +
      geom_line(linewidth = 0.6) +
      geom_point(size = 1.6) +
      geom_point(aes(y = observed_transition_cif), shape = 4, size = 2.0, stroke = 0.9, color = "black") +
      facet_wrap(~ dataset_key, scales = "free_y") +
      labs(title = "Stage 8B posterior predictive checks vs observed transition CIF", x = "Horizon (years)", y = "Transition CIF") +
      theme_bw() + theme(legend.position = "none")
  }

  grDevices::pdf(tmp_pdf, width = 11, height = 8.5, onefile = TRUE)
  pdf_open <- TRUE
  if (!is.null(transition_plot)) print(transition_plot)
  if (!is.null(remission_plot)) print(remission_plot)
  if (!is.null(nb_plot)) print(nb_plot)
  if (!is.null(ppc_plot)) print(ppc_plot)
  if (is.null(transition_plot) && is.null(remission_plot) && is.null(nb_plot) && is.null(ppc_plot)) {
    plot.new()
    text(0.5, 0.5, "No Stage 8B diagnostic summary plots were available.")
  }
  close_pdf()
  if (!pdf_file_is_usable(tmp_pdf)) stop("Temporary Stage8B diagnostic PDF was not created correctly.", call. = FALSE)
  safe_promote_file(tmp_pdf, final_pdf_path)

  if (!is.null(transition_plot)) safe_save_png(transition_plot, transition_png_path)
  if (!is.null(remission_plot)) safe_save_png(remission_plot, remission_png_path)
  if (!is.null(nb_plot)) safe_save_png(nb_plot, nb_png_path)
  if (!is.null(ppc_plot)) safe_save_png(ppc_plot, ppc_png_path)

  invisible(TRUE)
}

# 🔴 Override: per-model reuse rule for the complete Stage 8B script ===============================
is_stage8b_bundle_reusable <- function(bundle, model_row, dataset_df) {
  if (is.null(bundle) || !is.list(bundle) || is.null(bundle$model_registry_row)) return(FALSE)
  reg <- tibble::as_tibble(bundle$model_registry_row)
  if (nrow(reg) != 1L) return(FALSE)

  required_names <- c(
    "coefficient_summary",
    "diagnostics_parameter_level",
    "ppc_summary",
    "posterior_subject_profile",
    "posterior_subject_yearly",
    "posterior_cohort_yearly",
    "posterior_classification",
    "prior_predictive_summary",
    "hazard_shape_plausibility",
    "uncured_supporting_decomposition",
    "incidence_anchor_update"
  )

  has_required <- all(vapply(required_names, function(nm) {
    !is.null(bundle_get(bundle, nm))
  }, logical(1)))

  identical(bundle$version %||% "", "stage8b_v2_complete") &&
    has_required &&
    identical(as.character(reg$model_id[[1]]), as.character(model_row$model_id[[1]])) &&
    identical(as.character(reg$dataset_key[[1]]), as.character(model_row$dataset[[1]])) &&
    identical(as.character(reg$fit_status[[1]]), "ok") &&
    identical(as.integer(reg$n[[1]]), nrow(dataset_df))
}

# 🔴 Load: Stage 1 backbone, Stage 6 flags, and Stage 8A reference outputs ===============================
stage1_backbone <- load_stage1_backbone(
  bundle_file = stage1_bundle_file,
  datasets_file = stage1_analysis_datasets_file,
  dataset_registry_file = stage1_dataset_registry_file,
  horizon_registry_file = stage1_horizon_registry_file,
  threshold_registry_file = stage1_threshold_registry_file,
  scaling_registry_file = stage1_scaling_registry_file
)

analysis_datasets <- lapply(stage1_backbone$datasets, function(x) {
  x %>%
    normalize_dataset_time() %>%
    normalize_unique_person_id()
})

dataset_registry_stage1 <- tibble::as_tibble(stage1_backbone$registries$dataset_registry)
scaling_registry <- tibble::as_tibble(stage1_backbone$registries$scaling_registry)
horizon_registry_stage1 <- build_horizon_annotation_from_stage1(
  horizon_registry = stage1_backbone$registries$horizon_registry,
  datasets = analysis_datasets,
  horizons = 1:10
)
thresholds_from_stage1 <- build_threshold_vector_from_stage1(stage1_backbone$registries$threshold_registry)

horizons_year <- sort(unique(as.integer(horizon_registry_stage1$horizon_year)))
risk_thresholds <- sort(unique(as.numeric(thresholds_from_stage1)))

if (!identical(horizons_year, 1:10)) {
  stop("Stage 8B requires Stage 1 horizon grid 1:10.", call. = FALSE)
}
if (length(risk_thresholds) == 0L || anyNA(risk_thresholds) || any(risk_thresholds <= 0 | risk_thresholds >= 1)) {
  stop("Stage 8B risk thresholds must be probabilities strictly between 0 and 1.", call. = FALSE)
}

screening_lookup_stage6 <- read_stage6_screening_lookup(stage6_screening_flag_csv)

stage8a_outputs_raw <- load_stage8a_outputs(
  model_registry_csv = stage8a_model_registry_csv,
  performance_csv = stage8a_performance_csv
)
stage8a_delta_keys <- build_stage8a_delta_keys(stage8a_outputs_raw)

# 🔴 Prepare: reference tables, model grid, and Stan runtime ===============================
dataset_registry_stage8b <- bind_rows(lapply(names(analysis_datasets), function(ds) {
  dat <- analysis_datasets[[ds]]
  tibble(
    dataset_key = ds,
    n = nrow(dat),
    n_transition = sum(dat$status_num == 1L, na.rm = TRUE),
    n_remission = sum(dat$status_num == 2L, na.rm = TRUE),
    n_right_censoring = sum(dat$status_num == 0L, na.rm = TRUE),
    person_time_years = sum(dat$time_year, na.rm = TRUE),
    median_followup_years = stats::median(dat$time_year, na.rm = TRUE),
    max_followup_years = max(dat$time_year, na.rm = TRUE)
  )
}))

aj_registry <- lapply(names(analysis_datasets), function(ds) {
  build_aalen_johansen_reference(analysis_datasets[[ds]], horizons_year) %>%
    mutate(dataset_key = ds)
})
names(aj_registry) <- names(analysis_datasets)

ipcw_registry <- lapply(names(analysis_datasets), function(ds) {
  build_stage8b_ipcw_reference(analysis_datasets[[ds]], horizons_year) %>%
    mutate(dataset_key = ds)
})
names(ipcw_registry) <- names(analysis_datasets)

support_registry <- make_stage8b_support_registry(
  horizon_registry_stage1,
  datasets = names(analysis_datasets),
  horizons = horizons_year
)

risk_set_support_tbl <- build_risk_set_support_table(analysis_datasets, horizons_year)
prior_specs <- build_stage8b_prior_specs()

model_grid <- build_stage8b_model_grid(
  include_merged_incidence_site_supplementary = include_merged_incidence_site_supplementary,
  fit_prior_branches = fit_prior_branches,
  fit_site_prior_families = fit_site_prior_families
)

if (!is.null(run_model_ids)) {
  model_grid <- model_grid %>% filter(model_id %in% run_model_ids)
  if (nrow(model_grid) == 0L) stop("run_model_ids filtered out all Stage 8B models.", call. = FALSE)
}

stage6_model_lookup <- build_stage6_model_lookup(screening_lookup_stage6, model_grid)
model_grid <- model_grid %>% left_join(stage6_model_lookup, by = "model_id")

model_grid$stage8b_rds_path <- file.path(export_path, paste0(model_grid$model_id, "__bayes_stage8b_fit.rds"))
model_grid$reuse_existing_fit <- FALSE

if (isTRUE(reuse_existing_stage8b_rds)) {
  for (ii in seq_len(nrow(model_grid))) {
    bundle_i <- load_stage8b_reuse_bundle(model_grid$stage8b_rds_path[[ii]])
    model_grid$reuse_existing_fit[[ii]] <- is_stage8b_bundle_reusable(
      bundle = bundle_i,
      model_row = model_grid[ii, , drop = FALSE],
      dataset_df = analysis_datasets[[model_grid$dataset[[ii]]]]
    )
  }
}

n_models_reused <- sum(model_grid$reuse_existing_fit)
n_models_to_fit <- sum(!model_grid$reuse_existing_fit)
message("Stage8B reuse plan: ", n_models_reused, " model(s) reused; ", n_models_to_fit, " model(s) require fitting.")

if (n_models_to_fit > 0L) {
  fit_packages <- c("rstan", "posterior", "loo")
  missing_fit_packages <- fit_packages[!vapply(fit_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_fit_packages) > 0L) {
    stop("Install required fitting packages before running this script: ", paste(missing_fit_packages, collapse = ", "), call. = FALSE)
  }

  suppressPackageStartupMessages({
    library(rstan)
    library(posterior)
  })

  rstan_options(auto_write = TRUE)
  options(mc.cores = max(1L, min(stan_chains, parallel::detectCores(logical = TRUE))))
  stan_model_compiled <- compile_stage8b_stan_model()
} else {
  stan_model_compiled <- NULL
}

# 🔴 Execute: per-model Stage 8B fitting or bundle reuse ===============================
registry_rows <- list()
coef_rows <- list()
diag_param_rows <- list()
ppc_rows <- list()
subject_profile_rows <- list()
subject_yearly_rows <- list()
cohort_rows <- list()
classification_rows <- list()
prior_predictive_rows <- list()
hazard_shape_rows <- list()
uncured_support_rows <- list()
anchor_update_rows <- list()

withCallingHandlers({
  for (ii in seq_len(nrow(model_grid))) {
    model_row <- model_grid[ii, , drop = FALSE]
    model_id_now <- model_row$model_id[[1]]
    dataset_name <- model_row$dataset[[1]]
    dataset_df <- analysis_datasets[[dataset_name]]
    model_started_at <- Sys.time()

    emit_progress(
      ii - 1L,
      nrow(model_grid),
      model_id_now,
      paste0(
        "starting Stage8B (dataset=", dataset_name,
        ", family=", model_row$family_code[[1]],
        ", prior=", model_row$prior_branch[[1]],
        ", site-prior=", model_row$site_prior_family[[1]],
        ", reuse=", isTRUE(model_row$reuse_existing_fit[[1]]), ")"
      )
    )

    prior_spec <- prior_specs[[model_row$prior_branch[[1]]]]
    design_bundle <- make_stage8b_design_bundle(dataset_df, model_row, prior_spec, snu_site_label, remission_cut_years)

    set.seed(stan_seed + ii)
    prior_predictive_tbl <- simulate_stage8b_prior_predictive(
      dataset_df = dataset_df,
      model_row = model_row,
      design_bundle = design_bundle,
      prior_spec = prior_spec,
      n_draws = prior_predictive_draws,
      horizons_eval = prior_predictive_horizons
    )
    prior_predictive_rows[[length(prior_predictive_rows) + 1L]] <- prior_predictive_tbl

    if (isTRUE(model_row$reuse_existing_fit[[1]])) {
      bundle_reuse <- load_stage8b_reuse_bundle(model_row$stage8b_rds_path[[1]])

      reg_reuse <- tibble::as_tibble(bundle_get(bundle_reuse, "model_registry_row"))
      if (nrow_or_zero(reg_reuse) > 0L) {
        reg_reuse <- reg_reuse %>% mutate(fit_reused_flag = TRUE)
      }

      registry_rows[[length(registry_rows) + 1L]] <- reg_reuse
      coef_rows[[length(coef_rows) + 1L]] <- tibble::as_tibble(bundle_get(bundle_reuse, "coefficient_summary") %||% tibble())
      diag_param_rows[[length(diag_param_rows) + 1L]] <- tibble::as_tibble(bundle_get(bundle_reuse, "diagnostics_parameter_level") %||% tibble())
      ppc_rows[[length(ppc_rows) + 1L]] <- tibble::as_tibble(bundle_get(bundle_reuse, "ppc_summary") %||% tibble())
      subject_profile_rows[[length(subject_profile_rows) + 1L]] <- tibble::as_tibble(bundle_get(bundle_reuse, "posterior_subject_profile") %||% tibble())
      subject_yearly_rows[[length(subject_yearly_rows) + 1L]] <- tibble::as_tibble(bundle_get(bundle_reuse, "posterior_subject_yearly") %||% tibble())
      cohort_rows[[length(cohort_rows) + 1L]] <- tibble::as_tibble(bundle_get(bundle_reuse, "posterior_cohort_yearly") %||% tibble())
      classification_rows[[length(classification_rows) + 1L]] <- tibble::as_tibble(bundle_get(bundle_reuse, "posterior_classification") %||% tibble())
      hazard_shape_rows[[length(hazard_shape_rows) + 1L]] <- tibble::as_tibble(bundle_get(bundle_reuse, "hazard_shape_plausibility") %||% tibble())
      uncured_support_rows[[length(uncured_support_rows) + 1L]] <- tibble::as_tibble(bundle_get(bundle_reuse, "uncured_supporting_decomposition") %||% tibble())
      anchor_update_rows[[length(anchor_update_rows) + 1L]] <- tibble::as_tibble(bundle_get(bundle_reuse, "incidence_anchor_update") %||% tibble())

      emit_progress(ii, nrow(model_grid), model_id_now, paste0("reused existing Stage8B fit; elapsed=", format_number(elapsed_seconds(model_started_at), digits = 1L), "s"))
      next
    }

    site_idx <- get_stage8b_site_prior_indices(model_row, design_bundle)

    stan_data <- list(
      N = nrow(dataset_df),
      time = as.numeric(design_bundle$time),
      status = as.integer(design_bundle$status),
      K_inc = ncol(design_bundle$X_inc),
      X_inc = design_bundle$X_inc,
      K_lat = ncol(design_bundle$X_lat),
      X_lat = design_bundle$X_lat,
      K_rem = ncol(design_bundle$X_rem),
      X_rem = design_bundle$X_rem,
      family_id = as.integer(model_row$family_id[[1]]),
      alpha_prior_center = as.numeric(design_bundle$alpha_prior_center),
      mu_beta_inc = as.numeric(design_bundle$mu_beta_inc),
      sd_beta_inc = as.numeric(design_bundle$sd_beta_inc),
      sd_delta = as.numeric(prior_spec$sd_delta),
      sd_gamma0 = as.numeric(prior_spec$sd_gamma0),
      sd_gamma_lat = as.numeric(design_bundle$sd_gamma_lat),
      sd_xi0 = as.numeric(prior_spec$sd_xi0),
      sd_beta_rem = as.numeric(design_bundle$sd_beta_rem),
      sd_log_rho_rem = as.numeric(prior_spec$sd_log_rho_rem),
      sd_rho_W = as.numeric(prior_spec$sd_rho_W),
      sd_log_sigma_LN = as.numeric(prior_spec$sd_log_sigma_LN),
      sd_psi_LL = as.numeric(prior_spec$sd_psi_LL),
      J_rem = length(design_bundle$rem_cuts) - 1L,
      rem_cuts = as.numeric(design_bundle$rem_cuts),
      sd_site_prior = as.numeric(prior_spec$sd_site),
      inc_site_index = ifelse(is.na(site_idx$inc_idx), 0L, as.integer(site_idx$inc_idx)),
      lat_site_index = ifelse(is.na(site_idx$lat_idx), 0L, as.integer(site_idx$lat_idx)),
      rem_site_index = ifelse(is.na(site_idx$rem_idx), 0L, as.integer(site_idx$rem_idx)),
      use_t_prior_inc_site = as.integer(!is.na(site_idx$inc_idx) && identical(model_row$site_prior_family[[1]], "student_t3_0_1_sensitivity")),
      use_t_prior_lat_site = as.integer(!is.na(site_idx$lat_idx) && identical(model_row$site_prior_family[[1]], "student_t3_0_1_sensitivity")),
      use_t_prior_rem_site = as.integer(!is.na(site_idx$rem_idx) && identical(model_row$site_prior_family[[1]], "student_t3_0_1_sensitivity"))
    )

    fit_status <- "ok"
    fit_error_message <- NA_character_

    fit <- tryCatch(
      rstan::sampling(
        object = stan_model_compiled,
        data = stan_data,
        chains = stan_chains,
        iter = stan_iter,
        warmup = stan_warmup,
        thin = stan_thin,
        seed = stan_seed + ii,
        refresh = stan_refresh,
        control = list(adapt_delta = stan_adapt_delta, max_treedepth = stan_max_treedepth)
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      fit_status <- "sampling_error"
      fit_error_message <- conditionMessage(fit)

      registry_row <- tibble(
        dataset_key = dataset_name,
        model_id = model_id_now,
        retained_fit_id = model_row$retained_fit_id[[1]],
        structural_model_id = model_row$structural_model_id[[1]],
        formula_anchor = model_row$formula_anchor[[1]],
        transition_latency_branch = model_row$transition_latency_branch[[1]],
        remission_branch = model_row$remission_branch[[1]],
        site_placement_label = model_row$site_placement_label[[1]],
        branch = "Stage8B",
        risk_scale = model_row$risk_scale[[1]],
        prior_branch = model_row$prior_branch[[1]],
        site_prior_family = model_row$site_prior_family[[1]],
        latency_family = model_row$latency_family[[1]],
        family_code = model_row$family_code[[1]],
        fit_status = fit_status,
        fit_error_message = fit_error_message,
        admissibility_flag = FALSE,
        admissibility_reasons = "sampling_error",
        prior_degenerate_flag = any(prior_predictive_tbl$prior_degenerate_flag %in% TRUE),
        posterior_degenerate_flag = NA,
        ppc_gross_contradiction_flag = NA,
        coherence_violation_flag = NA,
        divergences = NA_integer_,
        max_rhat = NA_real_,
        min_bulk_ess = NA_real_,
        min_tail_ess = NA_real_,
        treedepth_exceeded = NA_integer_,
        waic = NA_real_,
        looic = NA_real_,
        p_waic = NA_real_,
        p_waic_high_n = NA_integer_,
        p_waic_high_pct = NA_real_,
        p_loo = NA_real_,
        pareto_k_max = NA_real_,
        pareto_k_bad_n = NA_integer_,
        pareto_k_bad_pct = NA_real_,
        pareto_k_very_bad_n = NA_integer_,
        waic_warning_flag = NA,
        loo_warning_flag = NA,
        info_criteria_warning_detail = NA_character_,
        n = nrow(dataset_df),
        n_transition = sum(dataset_df$status_num == 1L, na.rm = TRUE),
        n_remission = sum(dataset_df$status_num == 2L, na.rm = TRUE),
        n_right_censoring = sum(dataset_df$status_num == 0L, na.rm = TRUE),
        cure_model_eligibility_flag = model_row$cure_model_eligibility_flag[[1]],
        primary_gate_method = model_row$primary_gate_method[[1]],
        primary_gate_flag = model_row$primary_gate_flag[[1]],
        receus_primary_class = model_row$receus_primary_class[[1]],
        presence_modifier_flag = model_row$presence_modifier_flag[[1]],
        cure_presence_support_flag = model_row$cure_presence_support_flag[[1]],
        followup_contradiction_flag = model_row$followup_contradiction_flag[[1]],
        followup_not_contradicted_flag = model_row$followup_not_contradicted_flag[[1]],
        screening_note = model_row$screening_note[[1]],
        rds_path = model_row$stage8b_rds_path[[1]],
        fit_reused_flag = FALSE
      )

      saveRDS(
        list(version = "stage8b_v2_complete", model_registry_row = registry_row),
        model_row$stage8b_rds_path[[1]]
      )

      registry_rows[[length(registry_rows) + 1L]] <- registry_row
      emit_progress(ii, nrow(model_grid), model_id_now, paste0("sampling error after ", format_number(elapsed_seconds(model_started_at), digits = 1L), "s: ", fit_error_message))
      next
    }

    param_names <- c(
      "delta0", "gamma0", "xi0",
      paste0("beta_inc[", seq_len(ncol(design_bundle$X_inc)), "]"),
      paste0("gamma_lat[", seq_len(ncol(design_bundle$X_lat)), "]"),
      paste0("beta_rem[", seq_len(ncol(design_bundle$X_rem)), "]"),
      paste0("log_rho_rem[", seq_len(length(design_bundle$rem_cuts) - 1L), "]")
    )
    if (model_row$family_code[[1]] == "W") param_names <- c(param_names, "rho_W")
    if (model_row$family_code[[1]] == "LN") param_names <- c(param_names, "log_sigma_LN")
    if (model_row$family_code[[1]] == "LL") param_names <- c(param_names, "psi_LL")

    param_array <- posterior::as_draws_array(as.array(fit, pars = param_names))
    param_diag_tbl <- posterior::summarise_draws(
      param_array,
      mean = base::mean,
      sd = stats::sd,
      rhat = posterior::rhat,
      ess_bulk = posterior::ess_bulk,
      ess_tail = posterior::ess_tail
    )
    param_draws_mat <- posterior::as_draws_matrix(param_array)

    coef_tbl <- tibble(
      dataset_key = dataset_name,
      model_id = model_id_now,
      retained_fit_id = model_row$retained_fit_id[[1]],
      structural_model_id = model_row$structural_model_id[[1]],
      formula_anchor = model_row$formula_anchor[[1]],
      family_code = model_row$family_code[[1]],
      latency_family = model_row$latency_family[[1]],
      parameter = colnames(param_draws_mat),
      mean = apply(param_draws_mat, 2, mean),
      sd = apply(param_draws_mat, 2, stats::sd),
      q025 = apply(param_draws_mat, 2, stats::quantile, probs = 0.025, names = FALSE),
      q50 = apply(param_draws_mat, 2, stats::quantile, probs = 0.500, names = FALSE),
      q975 = apply(param_draws_mat, 2, stats::quantile, probs = 0.975, names = FALSE)
    )

    diag_param_tbl_model <- tibble(
      dataset_key = dataset_name,
      model_id = model_id_now,
      retained_fit_id = model_row$retained_fit_id[[1]],
      structural_model_id = model_row$structural_model_id[[1]],
      formula_anchor = model_row$formula_anchor[[1]],
      family_code = model_row$family_code[[1]],
      latency_family = model_row$latency_family[[1]],
      parameter = param_diag_tbl$variable,
      mean = param_diag_tbl$mean,
      sd = param_diag_tbl$sd,
      rhat = param_diag_tbl$rhat,
      ess_bulk = param_diag_tbl$ess_bulk,
      ess_tail = param_diag_tbl$ess_tail
    )

    sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    divergences <- sum(vapply(sampler_params, function(x) sum(x[, "divergent__"]), numeric(1)))
    treedepth_exceeded <- sum(vapply(sampler_params, function(x) sum(x[, "treedepth__"] >= stan_max_treedepth), numeric(1)))

    draws_compact <- extract_stage8b_draws_compact(
      fit = fit,
      K_inc = ncol(design_bundle$X_inc),
      K_lat = ncol(design_bundle$X_lat),
      K_rem = ncol(design_bundle$X_rem),
      J_rem = length(design_bundle$rem_cuts) - 1L
    )

    info_criteria <- compute_information_criteria(draws_compact$log_lik)

    total_draws <- length(draws_compact$delta0)
    set.seed(stan_seed + 200000L + ii)
    keep_draw_idx <- if (total_draws <= posterior_prediction_draws) {
      seq_len(total_draws)
    } else {
      sort(sample(seq_len(total_draws), size = posterior_prediction_draws, replace = FALSE))
    }

    draws_pred <- list(
      delta0 = draws_compact$delta0[keep_draw_idx],
      beta_inc = draws_compact$beta_inc[keep_draw_idx, , drop = FALSE],
      gamma0 = draws_compact$gamma0[keep_draw_idx],
      gamma_lat = draws_compact$gamma_lat[keep_draw_idx, , drop = FALSE],
      xi0 = draws_compact$xi0[keep_draw_idx],
      beta_rem = draws_compact$beta_rem[keep_draw_idx, , drop = FALSE],
      log_rho_rem = draws_compact$log_rho_rem[keep_draw_idx, , drop = FALSE],
      rho_W = draws_compact$rho_W[keep_draw_idx],
      log_sigma_LN = draws_compact$log_sigma_LN[keep_draw_idx],
      psi_LL = draws_compact$psi_LL[keep_draw_idx]
    )

    state <- compute_stage8b_linear_terms(
      draws = draws_pred,
      X_inc = design_bundle$X_inc,
      X_lat = design_bundle$X_lat,
      X_rem = design_bundle$X_rem,
      alpha_prior_center = design_bundle$alpha_prior_center
    )

    pred_trajectories <- build_stage8b_prediction_trajectories(
      state = state,
      draws = draws_pred,
      model_row = model_row,
      horizons = horizons_year,
      integration_step_year = integration_step_year,
      rem_cuts = design_bundle$rem_cuts
    )

    supported_horizons <- ppc_horizons_for_dataset(dataset_name)
    transition_supported <- lapply(as.character(supported_horizons), function(hh) pred_trajectories[[hh]]$transition_cif)
    posterior_degeneracy <- compute_stage8b_degeneracy(state$pi_mat, state$median_mat, transition_supported)

    cure_subj_sum <- summarize_cols_matrix(state$cure_prob_mat)
    susc_subj_sum <- summarize_cols_matrix(state$pi_mat)
    median_subj_sum <- summarize_cols_matrix(state$median_mat)

    subject_profile_tbl <- bind_cols(
      design_bundle$id_df,
      tibble(
        cure_fraction_mean = cure_subj_sum$mean,
        cure_fraction_q025 = cure_subj_sum$q025,
        cure_fraction_q50 = cure_subj_sum$q50,
        cure_fraction_q975 = cure_subj_sum$q975,
        susceptible_fraction_mean = susc_subj_sum$mean,
        susceptible_fraction_q025 = susc_subj_sum$q025,
        susceptible_fraction_q50 = susc_subj_sum$q50,
        susceptible_fraction_q975 = susc_subj_sum$q975,
        median_susceptible_time_mean = median_subj_sum$mean,
        median_susceptible_time_q025 = median_subj_sum$q025,
        median_susceptible_time_q50 = median_subj_sum$q50,
        median_susceptible_time_q975 = median_subj_sum$q975
      )
    ) %>%
      mutate(
        dataset_key = dataset_name,
        model_id = model_id_now,
        retained_fit_id = model_row$retained_fit_id[[1]],
        structural_model_id = model_row$structural_model_id[[1]],
        formula_anchor = model_row$formula_anchor[[1]],
        family_code = model_row$family_code[[1]],
        latency_family = model_row$latency_family[[1]],
        site_placement_label = model_row$site_placement_label[[1]],
        branch = "Stage8B",
        risk_scale = model_row$risk_scale[[1]],
        prior_branch = model_row$prior_branch[[1]],
        site_prior_family = model_row$site_prior_family[[1]]
      ) %>%
      relocate(dataset_key, model_id, retained_fit_id, structural_model_id, formula_anchor, family_code, latency_family, site_placement_label, branch, risk_scale, prior_branch, site_prior_family)

    horizon_refs <- ipcw_registry[[dataset_name]]
    subject_year_rows_model <- list()
    cohort_rows_model <- list()
    class_rows_model <- list()
    ppc_rows_model <- list()

    coherence_error_max <- 0
    hazard_means <- rep(NA_real_, length(horizons_year))
    names(hazard_means) <- paste0("hazard_", horizons_year, "y")

    for (hh in horizons_year) {
      pred_h <- pred_trajectories[[as.character(hh)]]
      horizon_ref <- horizon_refs %>% filter(horizon_year == hh)

      trans_sum <- summarize_cols_matrix(pred_h$transition_cif)
      rem_sum <- summarize_cols_matrix(pred_h$remission_cif)
      free_sum <- summarize_cols_matrix(pred_h$all_event_free)
      uncured_surv_sum <- summarize_cols_matrix(pred_h$uncured_survival)
      uncured_risk_sum <- summarize_cols_matrix(1 - pred_h$uncured_survival)

      coherence_error_max <- max(
        coherence_error_max,
        max(abs(pred_h$transition_cif + pred_h$remission_cif + pred_h$all_event_free - 1), na.rm = TRUE)
      )

      subject_year_rows_model[[length(subject_year_rows_model) + 1L]] <- bind_cols(
        design_bundle$id_df,
        tibble(
          transition_cif_mean = trans_sum$mean,
          transition_cif_q025 = trans_sum$q025,
          transition_cif_q50 = trans_sum$q50,
          transition_cif_q975 = trans_sum$q975,
          remission_cif_mean = rem_sum$mean,
          remission_cif_q025 = rem_sum$q025,
          remission_cif_q50 = rem_sum$q50,
          remission_cif_q975 = rem_sum$q975,
          all_event_free_mean = free_sum$mean,
          all_event_free_q025 = free_sum$q025,
          all_event_free_q50 = free_sum$q50,
          all_event_free_q975 = free_sum$q975,
          uncured_survival_mean = uncured_surv_sum$mean,
          uncured_survival_q025 = uncured_surv_sum$q025,
          uncured_survival_q50 = uncured_surv_sum$q50,
          uncured_survival_q975 = uncured_surv_sum$q975,
          uncured_risk_mean = uncured_risk_sum$mean,
          uncured_risk_q025 = uncured_risk_sum$q025,
          uncured_risk_q50 = uncured_risk_sum$q50,
          uncured_risk_q975 = uncured_risk_sum$q975
        )
      ) %>%
        mutate(
          dataset_key = dataset_name,
          model_id = model_id_now,
          retained_fit_id = model_row$retained_fit_id[[1]],
          structural_model_id = model_row$structural_model_id[[1]],
          formula_anchor = model_row$formula_anchor[[1]],
          family_code = model_row$family_code[[1]],
          latency_family = model_row$latency_family[[1]],
          site_placement_label = model_row$site_placement_label[[1]],
          branch = "Stage8B",
          risk_scale = model_row$risk_scale[[1]],
          prior_branch = model_row$prior_branch[[1]],
          site_prior_family = model_row$site_prior_family[[1]],
          horizon = hh
        ) %>%
        relocate(dataset_key, model_id, retained_fit_id, structural_model_id, formula_anchor, family_code, latency_family, site_placement_label, branch, risk_scale, prior_branch, site_prior_family, horizon)

      transition_cif_draw <- rowMeans(pred_h$transition_cif)
      remission_cif_draw <- rowMeans(pred_h$remission_cif)
      free_draw <- rowMeans(pred_h$all_event_free)
      hazard_draw <- rowMeans(pred_h$transition_population_hazard)
      cure_draw <- rowMeans(state$cure_prob_mat)
      susc_draw <- rowMeans(state$pi_mat)
      uncured_surv_draw <- rowMeans(pred_h$uncured_survival)
      uncured_risk_draw <- rowMeans(1 - pred_h$uncured_survival)

      hazard_means[[paste0("hazard_", hh, "y")]] <- mean(hazard_draw, na.rm = TRUE)

      cohort_rows_model[[length(cohort_rows_model) + 1L]] <- tibble(
        dataset_key = dataset_name,
        model_id = model_id_now,
        retained_fit_id = model_row$retained_fit_id[[1]],
        structural_model_id = model_row$structural_model_id[[1]],
        formula_anchor = model_row$formula_anchor[[1]],
        family_code = model_row$family_code[[1]],
        latency_family = model_row$latency_family[[1]],
        site_placement_label = model_row$site_placement_label[[1]],
        branch = "Stage8B",
        risk_scale = model_row$risk_scale[[1]],
        prior_branch = model_row$prior_branch[[1]],
        site_prior_family = model_row$site_prior_family[[1]],
        horizon = hh,
        transition_cif_mean = mean(transition_cif_draw),
        transition_cif_q025 = safe_quantile(transition_cif_draw, 0.025),
        transition_cif_q50 = safe_quantile(transition_cif_draw, 0.500),
        transition_cif_q975 = safe_quantile(transition_cif_draw, 0.975),
        remission_cif_mean = mean(remission_cif_draw),
        remission_cif_q025 = safe_quantile(remission_cif_draw, 0.025),
        remission_cif_q50 = safe_quantile(remission_cif_draw, 0.500),
        remission_cif_q975 = safe_quantile(remission_cif_draw, 0.975),
        all_event_free_mean = mean(free_draw),
        all_event_free_q025 = safe_quantile(free_draw, 0.025),
        all_event_free_q50 = safe_quantile(free_draw, 0.500),
        all_event_free_q975 = safe_quantile(free_draw, 0.975),
        cohort_mean_cure_fraction_mean = mean(cure_draw),
        cohort_mean_cure_fraction_q025 = safe_quantile(cure_draw, 0.025),
        cohort_mean_cure_fraction_q50 = safe_quantile(cure_draw, 0.500),
        cohort_mean_cure_fraction_q975 = safe_quantile(cure_draw, 0.975),
        cohort_mean_susceptible_fraction_mean = mean(susc_draw),
        cohort_mean_susceptible_fraction_q025 = safe_quantile(susc_draw, 0.025),
        cohort_mean_susceptible_fraction_q50 = safe_quantile(susc_draw, 0.500),
        cohort_mean_susceptible_fraction_q975 = safe_quantile(susc_draw, 0.975),
        mean_uncured_survival_mean = mean(uncured_surv_draw),
        mean_uncured_survival_q025 = safe_quantile(uncured_surv_draw, 0.025),
        mean_uncured_survival_q50 = safe_quantile(uncured_surv_draw, 0.500),
        mean_uncured_survival_q975 = safe_quantile(uncured_surv_draw, 0.975),
        mean_uncured_risk_mean = mean(uncured_risk_draw),
        mean_uncured_risk_q025 = safe_quantile(uncured_risk_draw, 0.025),
        mean_uncured_risk_q50 = safe_quantile(uncured_risk_draw, 0.500),
        mean_uncured_risk_q975 = safe_quantile(uncured_risk_draw, 0.975),
        mean_transition_hazard_mean = mean(hazard_draw),
        mean_transition_hazard_q025 = safe_quantile(hazard_draw, 0.025),
        mean_transition_hazard_q50 = safe_quantile(hazard_draw, 0.500),
        mean_transition_hazard_q975 = safe_quantile(hazard_draw, 0.975)
      )

      ppc_rows_model[[length(ppc_rows_model) + 1L]] <- tibble(
        dataset_key = dataset_name,
        model_id = model_id_now,
        retained_fit_id = model_row$retained_fit_id[[1]],
        structural_model_id = model_row$structural_model_id[[1]],
        formula_anchor = model_row$formula_anchor[[1]],
        family_code = model_row$family_code[[1]],
        branch = "Stage8B",
        risk_scale = model_row$risk_scale[[1]],
        prior_branch = model_row$prior_branch[[1]],
        site_prior_family = model_row$site_prior_family[[1]],
        horizon = hh,
        observed_transition_cif = horizon_ref$observed_transition_cif[[1]],
        observed_remission_cif = horizon_ref$observed_remission_cif[[1]],
        observed_all_event_free = horizon_ref$observed_all_event_free[[1]],
        posterior_mean_transition_cif = mean(transition_cif_draw),
        posterior_q025_transition_cif = safe_quantile(transition_cif_draw, 0.025),
        posterior_q975_transition_cif = safe_quantile(transition_cif_draw, 0.975),
        posterior_mean_remission_cif = mean(remission_cif_draw),
        posterior_q025_remission_cif = safe_quantile(remission_cif_draw, 0.025),
        posterior_q975_remission_cif = safe_quantile(remission_cif_draw, 0.975),
        absolute_difference_transition_cif = abs(mean(transition_cif_draw) - horizon_ref$observed_transition_cif[[1]]),
        gross_contradiction_flag = (
          (hh %in% ppc_horizons_for_dataset(dataset_name)) &&
            (
              horizon_ref$observed_transition_cif[[1]] < safe_quantile(transition_cif_draw, 0.025) ||
                horizon_ref$observed_transition_cif[[1]] > safe_quantile(transition_cif_draw, 0.975)
            ) &&
            abs(mean(transition_cif_draw) - horizon_ref$observed_transition_cif[[1]]) > ppc_tolerance_abs
        )
      )

      class_tbl_h <- compute_stage8b_classification_summary(
        risk_draws = pred_h$transition_cif,
        horizon_row = horizon_ref,
        thresholds = risk_thresholds
      ) %>%
        mutate(
          dataset_key = dataset_name,
          model_id = model_id_now,
          retained_fit_id = model_row$retained_fit_id[[1]],
          structural_model_id = model_row$structural_model_id[[1]],
          formula_anchor = model_row$formula_anchor[[1]],
          family_code = model_row$family_code[[1]],
          latency_family = model_row$latency_family[[1]],
          branch = "Stage8B",
          risk_scale = model_row$risk_scale[[1]],
          prior_branch = model_row$prior_branch[[1]],
          site_prior_family = model_row$site_prior_family[[1]],
          horizon = hh,
          observed_transition_cif = horizon_ref$observed_transition_cif[[1]],
          observed_remission_cif = horizon_ref$observed_remission_cif[[1]],
          observed_all_event_free = horizon_ref$observed_all_event_free[[1]]
        ) %>%
        relocate(dataset_key, model_id, retained_fit_id, structural_model_id, formula_anchor, family_code, latency_family, branch, risk_scale, prior_branch, site_prior_family, horizon, threshold)

      class_rows_model[[length(class_rows_model) + 1L]] <- class_tbl_h
    }

    ppc_model_tbl <- bind_rows(ppc_rows_model)
    subject_year_tbl <- bind_rows(subject_year_rows_model)
    cohort_model_tbl <- bind_rows(cohort_rows_model)
    class_model_tbl <- bind_rows(class_rows_model)

    ppc_gross_contradiction_flag <- any(ppc_model_tbl$gross_contradiction_flag, na.rm = TRUE)
    coherence_violation_flag <- is.finite(coherence_error_max) && coherence_error_max > 0.02

    max_rhat <- if (nrow(param_diag_tbl) > 0L) max(param_diag_tbl$rhat, na.rm = TRUE) else NA_real_
    min_bulk_ess <- if (nrow(param_diag_tbl) > 0L) min(param_diag_tbl$ess_bulk, na.rm = TRUE) else NA_real_
    min_tail_ess <- if (nrow(param_diag_tbl) > 0L) min(param_diag_tbl$ess_tail, na.rm = TRUE) else NA_real_

    admissibility_reasons <- character()
    if (any(prior_predictive_tbl$prior_degenerate_flag %in% TRUE)) admissibility_reasons <- c(admissibility_reasons, "prior_degenerate")
    if (posterior_degeneracy$degenerate_flag[[1]]) admissibility_reasons <- c(admissibility_reasons, "posterior_degenerate")
    if (!is.finite(max_rhat) || max_rhat >= rhat_max_threshold) admissibility_reasons <- c(admissibility_reasons, "rhat")
    if (!is.finite(min_bulk_ess) || min_bulk_ess < ess_min_threshold) admissibility_reasons <- c(admissibility_reasons, "bulk_ess")
    if (!is.finite(min_tail_ess) || min_tail_ess < ess_min_threshold) admissibility_reasons <- c(admissibility_reasons, "tail_ess")
    if (divergences > 0) admissibility_reasons <- c(admissibility_reasons, "divergences")
    if (treedepth_exceeded > 0) admissibility_reasons <- c(admissibility_reasons, "treedepth")
    if (isTRUE(ppc_gross_contradiction_flag)) admissibility_reasons <- c(admissibility_reasons, "ppc")
    if (isTRUE(coherence_violation_flag)) admissibility_reasons <- c(admissibility_reasons, "coherence")

    admissibility_flag <- length(admissibility_reasons) == 0L

    hazard_shape_tbl <- tibble(
      dataset_key = dataset_name,
      model_id = model_id_now,
      retained_fit_id = model_row$retained_fit_id[[1]],
      structural_model_id = model_row$structural_model_id[[1]],
      formula_anchor = model_row$formula_anchor[[1]],
      family_code = model_row$family_code[[1]],
      latency_family = model_row$latency_family[[1]],
      site_placement_label = model_row$site_placement_label[[1]],
      branch = "Stage8B",
      risk_scale = model_row$risk_scale[[1]],
      prior_branch = model_row$prior_branch[[1]],
      site_prior_family = model_row$site_prior_family[[1]],
      hazard_target = "transition_population_cause_specific_hazard",
      hazard_1y = hazard_means[["hazard_1y"]],
      hazard_2y = hazard_means[["hazard_2y"]],
      hazard_3y = hazard_means[["hazard_3y"]],
      hazard_4y = hazard_means[["hazard_4y"]],
      hazard_5y = hazard_means[["hazard_5y"]],
      hazard_6y = hazard_means[["hazard_6y"]],
      hazard_7y = hazard_means[["hazard_7y"]],
      hazard_8y = hazard_means[["hazard_8y"]],
      hazard_9y = hazard_means[["hazard_9y"]],
      hazard_10y = hazard_means[["hazard_10y"]],
      hazard_ratio_10y_vs_1y = safe_numeric(hazard_means[["hazard_10y"]]) / pmax(safe_numeric(hazard_means[["hazard_1y"]]), 1e-12),
      shape_class = classify_hazard_shape(unname(hazard_means))
    )

    cure_draw <- rowMeans(state$cure_prob_mat)
    susc_draw <- rowMeans(state$pi_mat)
    uncured_support_tbl <- bind_rows(lapply(horizons_year, function(hh) {
      pred_h <- pred_trajectories[[as.character(hh)]]
      uncured_surv_draw <- rowMeans(pred_h$uncured_survival)
      uncured_risk_draw <- rowMeans(1 - pred_h$uncured_survival)
      tibble(
        dataset_key = dataset_name,
        model_id = model_id_now,
        retained_fit_id = model_row$retained_fit_id[[1]],
        structural_model_id = model_row$structural_model_id[[1]],
        formula_anchor = model_row$formula_anchor[[1]],
        family_code = model_row$family_code[[1]],
        latency_family = model_row$latency_family[[1]],
        branch = "Stage8B",
        risk_scale = model_row$risk_scale[[1]],
        prior_branch = model_row$prior_branch[[1]],
        site_prior_family = model_row$site_prior_family[[1]],
        horizon = hh,
        cure_fraction_mean = mean(cure_draw),
        cure_fraction_q025 = safe_quantile(cure_draw, 0.025),
        cure_fraction_q50 = safe_quantile(cure_draw, 0.500),
        cure_fraction_q975 = safe_quantile(cure_draw, 0.975),
        susceptible_fraction_mean = mean(susc_draw),
        susceptible_fraction_q025 = safe_quantile(susc_draw, 0.025),
        susceptible_fraction_q50 = safe_quantile(susc_draw, 0.500),
        susceptible_fraction_q975 = safe_quantile(susc_draw, 0.975),
        uncured_survival_mean = mean(uncured_surv_draw),
        uncured_survival_q025 = safe_quantile(uncured_surv_draw, 0.025),
        uncured_survival_q50 = safe_quantile(uncured_surv_draw, 0.500),
        uncured_survival_q975 = safe_quantile(uncured_surv_draw, 0.975),
        uncured_risk_mean = mean(uncured_risk_draw),
        uncured_risk_q025 = safe_quantile(uncured_risk_draw, 0.025),
        uncured_risk_q50 = safe_quantile(uncured_risk_draw, 0.500),
        uncured_risk_q975 = safe_quantile(uncured_risk_draw, 0.975),
        MSTu_mean = NA_real_,
        MSTu_q025 = NA_real_,
        MSTu_q50 = NA_real_,
        MSTu_q975 = NA_real_,
        uncured_mean_support_flag = FALSE
      )
    }))

    anchor_update_tbl <- make_stage8b_incidence_anchor_update(model_row, design_bundle, prior_spec, draws_compact)

    registry_row <- tibble(
      dataset_key = dataset_name,
      model_id = model_id_now,
      retained_fit_id = model_row$retained_fit_id[[1]],
      structural_model_id = model_row$structural_model_id[[1]],
      formula_anchor = model_row$formula_anchor[[1]],
      transition_latency_branch = model_row$transition_latency_branch[[1]],
      remission_branch = model_row$remission_branch[[1]],
      site_placement_label = model_row$site_placement_label[[1]],
      branch = "Stage8B",
      risk_scale = model_row$risk_scale[[1]],
      prior_branch = model_row$prior_branch[[1]],
      site_prior_family = model_row$site_prior_family[[1]],
      latency_family = model_row$latency_family[[1]],
      family_code = model_row$family_code[[1]],
      fit_status = fit_status,
      fit_error_message = fit_error_message,
      admissibility_flag = admissibility_flag,
      admissibility_reasons = if (length(admissibility_reasons) == 0L) "" else paste(admissibility_reasons, collapse = "|"),
      prior_degenerate_flag = any(prior_predictive_tbl$prior_degenerate_flag %in% TRUE),
      posterior_degenerate_flag = posterior_degeneracy$degenerate_flag[[1]],
      ppc_gross_contradiction_flag = ppc_gross_contradiction_flag,
      coherence_violation_flag = coherence_violation_flag,
      divergences = divergences,
      max_rhat = max_rhat,
      min_bulk_ess = min_bulk_ess,
      min_tail_ess = min_tail_ess,
      treedepth_exceeded = treedepth_exceeded,
      waic = info_criteria$waic,
      looic = info_criteria$looic,
      p_waic = info_criteria$p_waic,
      p_waic_high_n = info_criteria$p_waic_high_n,
      p_waic_high_pct = info_criteria$p_waic_high_pct,
      p_loo = info_criteria$p_loo,
      pareto_k_max = info_criteria$pareto_k_max,
      pareto_k_bad_n = info_criteria$pareto_k_bad_n,
      pareto_k_bad_pct = info_criteria$pareto_k_bad_pct,
      pareto_k_very_bad_n = info_criteria$pareto_k_very_bad_n,
      waic_warning_flag = info_criteria$waic_warning_flag,
      loo_warning_flag = info_criteria$loo_warning_flag,
      info_criteria_warning_detail = info_criteria$info_criteria_warning_detail,
      n = nrow(dataset_df),
      n_transition = sum(dataset_df$status_num == 1L, na.rm = TRUE),
      n_remission = sum(dataset_df$status_num == 2L, na.rm = TRUE),
      n_right_censoring = sum(dataset_df$status_num == 0L, na.rm = TRUE),
      cure_model_eligibility_flag = model_row$cure_model_eligibility_flag[[1]],
      primary_gate_method = model_row$primary_gate_method[[1]],
      primary_gate_flag = model_row$primary_gate_flag[[1]],
      receus_primary_class = model_row$receus_primary_class[[1]],
      presence_modifier_flag = model_row$presence_modifier_flag[[1]],
      cure_presence_support_flag = model_row$cure_presence_support_flag[[1]],
      followup_contradiction_flag = model_row$followup_contradiction_flag[[1]],
      followup_not_contradicted_flag = model_row$followup_not_contradicted_flag[[1]],
      screening_note = model_row$screening_note[[1]],
      cohort_mean_cure_fraction_mean = mean(cure_draw),
      cohort_mean_cure_fraction_q025 = safe_quantile(cure_draw, 0.025),
      cohort_mean_cure_fraction_q50 = safe_quantile(cure_draw, 0.500),
      cohort_mean_cure_fraction_q975 = safe_quantile(cure_draw, 0.975),
      rds_path = model_row$stage8b_rds_path[[1]],
      fit_reused_flag = FALSE
    )

    save_obj <- list(
      version = "stage8b_v2_complete",
      model_registry_row = registry_row,
      coefficient_summary = coef_tbl,
      diagnostics_parameter_level = diag_param_tbl_model,
      ppc_summary = ppc_model_tbl,
      posterior_subject_profile = if (isTRUE(admissibility_flag)) subject_profile_tbl else tibble(),
      posterior_subject_yearly = if (isTRUE(admissibility_flag)) subject_year_tbl else tibble(),
      posterior_cohort_yearly = if (isTRUE(admissibility_flag)) cohort_model_tbl else tibble(),
      posterior_classification = if (isTRUE(admissibility_flag)) class_model_tbl else tibble(),
      prior_predictive_summary = prior_predictive_tbl,
      hazard_shape_plausibility = if (isTRUE(admissibility_flag)) hazard_shape_tbl else tibble(),
      uncured_supporting_decomposition = if (isTRUE(admissibility_flag)) uncured_support_tbl else tibble(),
      incidence_anchor_update = if (isTRUE(admissibility_flag)) anchor_update_tbl else tibble(),
      outputs = list(
        model_registry_row = registry_row,
        coefficient_summary = coef_tbl,
        diagnostics_parameter_level = diag_param_tbl_model,
        ppc_summary = ppc_model_tbl,
        posterior_subject_profile = if (isTRUE(admissibility_flag)) subject_profile_tbl else tibble(),
        posterior_subject_yearly = if (isTRUE(admissibility_flag)) subject_year_tbl else tibble(),
        posterior_cohort_yearly = if (isTRUE(admissibility_flag)) cohort_model_tbl else tibble(),
        posterior_classification = if (isTRUE(admissibility_flag)) class_model_tbl else tibble(),
        prior_predictive_summary = prior_predictive_tbl,
        hazard_shape_plausibility = if (isTRUE(admissibility_flag)) hazard_shape_tbl else tibble(),
        uncured_supporting_decomposition = if (isTRUE(admissibility_flag)) uncured_support_tbl else tibble(),
        incidence_anchor_update = if (isTRUE(admissibility_flag)) anchor_update_tbl else tibble()
      )
    )
    if (isTRUE(save_full_stanfit_rds)) save_obj$fit <- fit
    saveRDS(save_obj, model_row$stage8b_rds_path[[1]])

    registry_rows[[length(registry_rows) + 1L]] <- registry_row
    coef_rows[[length(coef_rows) + 1L]] <- coef_tbl
    diag_param_rows[[length(diag_param_rows) + 1L]] <- diag_param_tbl_model
    ppc_rows[[length(ppc_rows) + 1L]] <- ppc_model_tbl
    if (isTRUE(admissibility_flag)) {
      subject_profile_rows[[length(subject_profile_rows) + 1L]] <- subject_profile_tbl
      subject_yearly_rows[[length(subject_yearly_rows) + 1L]] <- subject_year_tbl
      cohort_rows[[length(cohort_rows) + 1L]] <- cohort_model_tbl
      classification_rows[[length(classification_rows) + 1L]] <- class_model_tbl
      hazard_shape_rows[[length(hazard_shape_rows) + 1L]] <- hazard_shape_tbl
      uncured_support_rows[[length(uncured_support_rows) + 1L]] <- uncured_support_tbl
      anchor_update_rows[[length(anchor_update_rows) + 1L]] <- anchor_update_tbl
    }

    rm(fit, param_array, param_draws_mat, draws_compact, draws_pred, state, pred_trajectories)
    gc(verbose = FALSE)

    emit_progress(
      ii,
      nrow(model_grid),
      model_id_now,
      paste0(
        "completed; elapsed=", format_number(elapsed_seconds(model_started_at), digits = 1L),
        "s; WAIC=", format_number(info_criteria$waic, digits = 2L),
        "; LOOIC=", format_number(info_criteria$looic, digits = 2L),
        "; Pareto k max=", format_number(info_criteria$pareto_k_max, digits = 3L)
      )
    )
  }
}, warning = function(w) {
  if (is_localhost_connection_warning(w)) tryInvokeRestart("muffleWarning")
})

# 🔴 Assemble: run-level tables, annotations, and prior-tail governance ===============================
model_order <- model_grid$model_id

model_registry <- bind_rows_safe(registry_rows) %>%
  arrange(factor(model_id, levels = model_order))

coefficient_summary <- bind_rows_safe(coef_rows) %>%
  arrange(factor(model_id, levels = model_order), parameter)

diagnostics_parameter_level <- bind_rows_safe(diag_param_rows) %>%
  arrange(factor(model_id, levels = model_order), parameter)

ppc_summary <- bind_rows_safe(ppc_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon)

posterior_subject_profile <- bind_rows_safe(subject_profile_rows) %>%
  arrange(factor(model_id, levels = model_order), unique_person_id)

posterior_subject_yearly <- bind_rows_safe(subject_yearly_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon, unique_person_id)

posterior_cohort_yearly <- bind_rows_safe(cohort_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon)

posterior_classification_raw <- bind_rows_safe(classification_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon, threshold)

hazard_shape_plausibility <- bind_rows_safe(hazard_shape_rows) %>%
  arrange(factor(model_id, levels = model_order))

uncured_supporting_decomposition <- bind_rows_safe(uncured_support_rows) %>%
  arrange(factor(model_id, levels = model_order), horizon)

incidence_anchor_update <- bind_rows_safe(anchor_update_rows) %>%
  arrange(factor(model_id, levels = model_order), age_sex_anchor_cell)

prior_predictive_summary <- annotate_stage8b_prior_predictive(bind_rows_safe(prior_predictive_rows)) %>%
  arrange(factor(model_id, levels = model_order), metric, prior_branch, horizon)

model_annotation <- model_registry %>%
  select(
    model_id,
    retained_fit_id,
    admissibility_flag,
    admissibility_reasons,
    cure_model_eligibility_flag,
    primary_gate_method,
    primary_gate_flag,
    receus_primary_class,
    presence_modifier_flag,
    cure_presence_support_flag,
    followup_contradiction_flag,
    followup_not_contradicted_flag,
    screening_note,
    fit_reused_flag
  )

horizon_annotation <- support_registry %>%
  transmute(
    dataset_key = dataset,
    horizon = horizon_year,
    support_tier = support_tier,
    horizon_evidence_class = horizon_evidence_class,
    claim_restriction_flag = claim_restriction_flag,
    interpretation_note = interpretation_note
  )

add_horizon_and_model_annotations <- function(df, has_horizon = TRUE) {
  out <- tibble::as_tibble(df)
  if (nrow_or_zero(out) == 0L) return(out)
  if (isTRUE(has_horizon) && all(c("dataset_key", "horizon") %in% names(out))) {
    out <- out %>%
      left_join_replacing_columns(horizon_annotation, by = c("dataset_key", "horizon")) %>%
      left_join_replacing_columns(risk_set_support_tbl, by = c("dataset_key", "horizon"))
  }
  if (all(c("model_id", "retained_fit_id") %in% names(out))) {
    out <- out %>%
      left_join_replacing_columns(model_annotation, by = c("model_id", "retained_fit_id"))
  }
  out
}

posterior_subject_profile <- add_horizon_and_model_annotations(posterior_subject_profile, has_horizon = FALSE)
posterior_subject_yearly <- add_horizon_and_model_annotations(posterior_subject_yearly, has_horizon = TRUE)
posterior_cohort_yearly <- add_horizon_and_model_annotations(posterior_cohort_yearly, has_horizon = TRUE)
posterior_classification_raw <- add_horizon_and_model_annotations(posterior_classification_raw, has_horizon = TRUE)
ppc_summary <- add_horizon_and_model_annotations(ppc_summary, has_horizon = TRUE)
uncured_supporting_decomposition <- add_horizon_and_model_annotations(uncured_supporting_decomposition, has_horizon = TRUE)

prior_warning_by_model <- prior_predictive_summary %>%
  group_by(model_id, retained_fit_id) %>%
  summarise(
    prior_tail_warning_flag = any(prior_tail_warning_flag %in% TRUE),
    prior_tail_warning_detail = paste(unique(na.omit(prior_tail_warning_detail[prior_tail_warning_flag %in% TRUE])), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(prior_tail_warning_detail = dplyr::na_if(prior_tail_warning_detail, ""))

model_registry <- model_registry %>%
  left_join(prior_warning_by_model, by = c("model_id", "retained_fit_id")) %>%
  mutate(prior_tail_warning_flag = coalesce(prior_tail_warning_flag, FALSE))

posterior_cohort_yearly <- posterior_cohort_yearly %>%
  mutate(DeltaRisk_NC_C_mean = mean_uncured_risk_mean - transition_cif_mean)

cohort_delta_long <- bind_rows(
  posterior_cohort_yearly %>%
    transmute(
      dataset_key, branch, risk_scale, retained_fit_id, structural_model_id, formula_anchor, family_code, site_prior_family,
      model_id, prior_branch, horizon, threshold = NA_real_, metric = "transition_cif", estimate = transition_cif_mean
    ),
  posterior_cohort_yearly %>%
    transmute(
      dataset_key, branch, risk_scale, retained_fit_id, structural_model_id, formula_anchor, family_code, site_prior_family,
      model_id, prior_branch, horizon, threshold = NA_real_, metric = "cure_fraction", estimate = cohort_mean_cure_fraction_mean
    )
)

class_delta_long <- bind_rows(
  posterior_classification_raw %>% transmute(dataset_key, branch, risk_scale, retained_fit_id, structural_model_id, formula_anchor, family_code, site_prior_family, model_id, prior_branch, horizon, threshold, metric = "false_positive_burden", estimate = false_positive_burden_mean),
  posterior_classification_raw %>% transmute(dataset_key, branch, risk_scale, retained_fit_id, structural_model_id, formula_anchor, family_code, site_prior_family, model_id, prior_branch, horizon, threshold, metric = "FP100", estimate = FP100_mean),
  posterior_classification_raw %>% transmute(dataset_key, branch, risk_scale, retained_fit_id, structural_model_id, formula_anchor, family_code, site_prior_family, model_id, prior_branch, horizon, threshold, metric = "NB", estimate = NB_mean),
  posterior_classification_raw %>% transmute(dataset_key, branch, risk_scale, retained_fit_id, structural_model_id, formula_anchor, family_code, site_prior_family, model_id, prior_branch, horizon, threshold, metric = "PPV", estimate = PPV_mean),
  posterior_classification_raw %>% transmute(dataset_key, branch, risk_scale, retained_fit_id, structural_model_id, formula_anchor, family_code, site_prior_family, model_id, prior_branch, horizon, threshold, metric = "TPR", estimate = TPR_mean)
)

all_prior_delta_long <- bind_rows(cohort_delta_long, class_delta_long)
anchor_tbl <- all_prior_delta_long %>%
  filter(prior_branch == "anchor_informed") %>%
  rename(model_id_anchor = model_id, anchor_estimate = estimate)
neutral_tbl <- all_prior_delta_long %>%
  filter(prior_branch == "neutral_no_external_info") %>%
  rename(model_id_neutral = model_id, neutral_estimate = estimate)

anchor_vs_neutral_delta <- full_join(
  anchor_tbl %>% select(-prior_branch),
  neutral_tbl %>% select(-prior_branch),
  by = c("dataset_key", "branch", "risk_scale", "retained_fit_id", "structural_model_id", "formula_anchor", "family_code", "site_prior_family", "horizon", "threshold", "metric")
) %>%
  transmute(
    dataset_key = dataset_key,
    branch = branch,
    risk_scale = risk_scale,
    retained_fit_id = retained_fit_id,
    structural_model_id = structural_model_id,
    formula_anchor = formula_anchor,
    family_code = family_code,
    site_prior_family = site_prior_family,
    horizon = horizon,
    threshold = threshold,
    metric = metric,
    delta_field = metric_to_anchor_delta_field(metric),
    model_id_anchor = model_id_anchor,
    model_id_neutral = model_id_neutral,
    anchor_estimate = anchor_estimate,
    neutral_estimate = neutral_estimate,
    delta_anchor_minus_neutral = anchor_estimate - neutral_estimate
  ) %>%
  left_join_replacing_columns(horizon_annotation, by = c("dataset_key", "horizon")) %>%
  arrange(retained_fit_id, horizon, threshold, metric)

prior_tail_sensitive_lookup <- anchor_vs_neutral_delta %>%
  mutate(
    threshold_key = ifelse(is.na(threshold), "__NA__", sprintf("%.10f", threshold)),
    material_flag = case_when(
      metric %in% c("transition_cif", "cure_fraction", "PPV", "TPR") ~ abs(delta_anchor_minus_neutral) > prior_materiality_risk,
      metric == "false_positive_burden" ~ abs(delta_anchor_minus_neutral) > prior_materiality_false_positive_burden,
      metric == "FP100" ~ abs(delta_anchor_minus_neutral) > 100 * prior_materiality_false_positive_burden,
      metric == "NB" ~ abs(delta_anchor_minus_neutral) > prior_materiality_nb,
      TRUE ~ FALSE
    )
  ) %>%
  group_by(dataset_key, retained_fit_id, horizon, threshold_key) %>%
  summarise(prior_tail_sensitive = any(material_flag %in% TRUE, na.rm = TRUE), .groups = "drop")

apply_prior_tail_sensitive <- function(df, has_threshold = FALSE) {
  if (nrow_or_zero(df) == 0L) return(df)
  df %>%
    mutate(threshold_key = if (has_threshold) ifelse(is.na(threshold), "__NA__", sprintf("%.10f", threshold)) else "__NA__") %>%
    left_join(prior_tail_sensitive_lookup, by = c("dataset_key", "retained_fit_id", "horizon", "threshold_key")) %>%
    left_join(prior_warning_by_model %>% select(model_id, retained_fit_id, prior_tail_warning_flag), by = c("model_id", "retained_fit_id")) %>%
    mutate(
      prior_tail_sensitive = coalesce(prior_tail_sensitive, FALSE) | coalesce(prior_tail_warning_flag, FALSE),
      claim_restriction_flag = ifelse(
        claim_restriction_flag == "projection_only" & prior_tail_sensitive,
        "projection_plus_prior_sensitive",
        claim_restriction_flag
      )
    ) %>%
    select(-threshold_key, -prior_tail_warning_flag)
}

posterior_subject_yearly <- apply_prior_tail_sensitive(posterior_subject_yearly, has_threshold = FALSE)
posterior_cohort_yearly <- apply_prior_tail_sensitive(posterior_cohort_yearly, has_threshold = FALSE)
posterior_classification_raw <- apply_prior_tail_sensitive(posterior_classification_raw, has_threshold = TRUE)
ppc_summary <- apply_prior_tail_sensitive(ppc_summary, has_threshold = FALSE)
uncured_supporting_decomposition <- apply_prior_tail_sensitive(uncured_supporting_decomposition, has_threshold = FALSE)

# 🔴 Compare: Stage 8B against Stage 8A and compute horizon-level performance ===============================
stage8a_cohort_key <- tibble::as_tibble(stage8a_delta_keys$cohort_key)
stage8a_class_key <- tibble::as_tibble(stage8a_delta_keys$class_key)

stage8a_has_cohort_key <- nrow_or_zero(stage8a_cohort_key) > 0L &&
  all(c("dataset_key", "structural_model_id", "formula_anchor", "family_code", "prior_branch", "site_prior_family", "horizon", "stage8a_transition_risk", "stage8a_cure_fraction") %in% names(stage8a_cohort_key))

stage8a_has_class_key <- nrow_or_zero(stage8a_class_key) > 0L &&
  all(c("dataset_key", "structural_model_id", "formula_anchor", "family_code", "prior_branch", "site_prior_family", "horizon", "threshold", "stage8a_false_positive_burden", "stage8a_FP100", "stage8a_NB", "stage8a_PPV", "stage8a_TPR") %in% names(stage8a_class_key))

delta_vs_stage8a_parts <- list(
  if (stage8a_has_cohort_key) {
    posterior_cohort_yearly %>%
      left_join(stage8a_cohort_key, by = c("dataset_key", "structural_model_id", "formula_anchor", "family_code", "prior_branch", "site_prior_family", "horizon")) %>%
      transmute(dataset_key, branch = "Stage8B", risk_scale = "transition_cif_competing", model_id, retained_fit_id, structural_model_id, formula_anchor, family_code, prior_branch, site_prior_family, horizon, threshold = NA_real_, metric = "transition_cif", delta_field = metric_to_stage8a_delta_field("transition_cif"), stage8b_estimate = transition_cif_mean, stage8a_estimate = stage8a_transition_risk, delta_8B_minus_8A = transition_cif_mean - stage8a_transition_risk)
  },
  if (stage8a_has_cohort_key) {
    posterior_cohort_yearly %>%
      left_join(stage8a_cohort_key, by = c("dataset_key", "structural_model_id", "formula_anchor", "family_code", "prior_branch", "site_prior_family", "horizon")) %>%
      transmute(dataset_key, branch = "Stage8B", risk_scale = "transition_cif_competing", model_id, retained_fit_id, structural_model_id, formula_anchor, family_code, prior_branch, site_prior_family, horizon, threshold = NA_real_, metric = "cure_fraction", delta_field = metric_to_stage8a_delta_field("cure_fraction"), stage8b_estimate = cohort_mean_cure_fraction_mean, stage8a_estimate = stage8a_cure_fraction, delta_8B_minus_8A = cohort_mean_cure_fraction_mean - stage8a_cure_fraction)
  },
  if (stage8a_has_class_key) {
    posterior_classification_raw %>%
      left_join(stage8a_class_key, by = c("dataset_key", "structural_model_id", "formula_anchor", "family_code", "prior_branch", "site_prior_family", "horizon", "threshold")) %>%
      transmute(dataset_key, branch = "Stage8B", risk_scale = "transition_cif_competing", model_id, retained_fit_id, structural_model_id, formula_anchor, family_code, prior_branch, site_prior_family, horizon, threshold, metric = "false_positive_burden", delta_field = metric_to_stage8a_delta_field("false_positive_burden"), stage8b_estimate = false_positive_burden_mean, stage8a_estimate = stage8a_false_positive_burden, delta_8B_minus_8A = false_positive_burden_mean - stage8a_false_positive_burden)
  },
  if (stage8a_has_class_key) {
    posterior_classification_raw %>%
      left_join(stage8a_class_key, by = c("dataset_key", "structural_model_id", "formula_anchor", "family_code", "prior_branch", "site_prior_family", "horizon", "threshold")) %>%
      transmute(dataset_key, branch = "Stage8B", risk_scale = "transition_cif_competing", model_id, retained_fit_id, structural_model_id, formula_anchor, family_code, prior_branch, site_prior_family, horizon, threshold, metric = "FP100", delta_field = metric_to_stage8a_delta_field("FP100"), stage8b_estimate = FP100_mean, stage8a_estimate = stage8a_FP100, delta_8B_minus_8A = FP100_mean - stage8a_FP100)
  },
  if (stage8a_has_class_key) {
    posterior_classification_raw %>%
      left_join(stage8a_class_key, by = c("dataset_key", "structural_model_id", "formula_anchor", "family_code", "prior_branch", "site_prior_family", "horizon", "threshold")) %>%
      transmute(dataset_key, branch = "Stage8B", risk_scale = "transition_cif_competing", model_id, retained_fit_id, structural_model_id, formula_anchor, family_code, prior_branch, site_prior_family, horizon, threshold, metric = "NB", delta_field = metric_to_stage8a_delta_field("NB"), stage8b_estimate = NB_mean, stage8a_estimate = stage8a_NB, delta_8B_minus_8A = NB_mean - stage8a_NB)
  },
  if (stage8a_has_class_key) {
    posterior_classification_raw %>%
      left_join(stage8a_class_key, by = c("dataset_key", "structural_model_id", "formula_anchor", "family_code", "prior_branch", "site_prior_family", "horizon", "threshold")) %>%
      transmute(dataset_key, branch = "Stage8B", risk_scale = "transition_cif_competing", model_id, retained_fit_id, structural_model_id, formula_anchor, family_code, prior_branch, site_prior_family, horizon, threshold, metric = "PPV", delta_field = metric_to_stage8a_delta_field("PPV"), stage8b_estimate = PPV_mean, stage8a_estimate = stage8a_PPV, delta_8B_minus_8A = PPV_mean - stage8a_PPV)
  },
  if (stage8a_has_class_key) {
    posterior_classification_raw %>%
      left_join(stage8a_class_key, by = c("dataset_key", "structural_model_id", "formula_anchor", "family_code", "prior_branch", "site_prior_family", "horizon", "threshold")) %>%
      transmute(dataset_key, branch = "Stage8B", risk_scale = "transition_cif_competing", model_id, retained_fit_id, structural_model_id, formula_anchor, family_code, prior_branch, site_prior_family, horizon, threshold, metric = "TPR", delta_field = metric_to_stage8a_delta_field("TPR"), stage8b_estimate = TPR_mean, stage8a_estimate = stage8a_TPR, delta_8B_minus_8A = TPR_mean - stage8a_TPR)
  }
)

delta_vs_stage8a <- bind_rows_safe(delta_vs_stage8a_parts)
if (nrow_or_zero(delta_vs_stage8a) > 0L) {
  delta_vs_stage8a <- delta_vs_stage8a %>%
    left_join_replacing_columns(horizon_annotation, by = c("dataset_key", "horizon")) %>%
    arrange(factor(model_id, levels = model_order), horizon, threshold, metric)
}

truth_lookup <- bind_rows(lapply(names(analysis_datasets), function(ds) {
  analysis_datasets[[ds]] %>%
    transmute(dataset_key = ds, unique_person_id, time_year, status_num)
}))

posterior_subject_yearly_eval <- posterior_subject_yearly %>%
  left_join(truth_lookup, by = c("dataset_key", "unique_person_id"))

horizon_performance_tbl <- posterior_subject_yearly_eval %>%
  group_by(
    dataset_key, model_id, retained_fit_id, structural_model_id, formula_anchor, family_code,
    latency_family, branch, risk_scale, prior_branch, site_prior_family, horizon
  ) %>%
  group_modify(~ {
    dataset_now <- as.character(.y$dataset_key[[1]])
    horizon_now <- as.integer(.y$horizon[[1]])
    horizon_ref <- ipcw_registry[[dataset_now]] %>% filter(horizon_year == horizon_now)
    x_eval <- .x %>% mutate(horizon = horizon_now)
    cal_tbl <- compute_calibration_simple(x_eval, observed_transition_cif = horizon_ref$observed_transition_cif[[1]])
    disc_tbl <- compute_discrimination_simple(x_eval)
    brier_tbl <- compute_brier_ipcw_simple(x_eval)
    bind_cols(
      tibble(
        observed_transition_cif = horizon_ref$observed_transition_cif[[1]],
        observed_remission_cif = horizon_ref$observed_remission_cif[[1]],
        observed_all_event_free = horizon_ref$observed_all_event_free[[1]]
      ),
      disc_tbl,
      cal_tbl,
      brier_tbl
    )
  }) %>%
  ungroup() %>%
  group_by(
    dataset_key, model_id, retained_fit_id, structural_model_id, formula_anchor, family_code,
    latency_family, branch, risk_scale, prior_branch, site_prior_family
  ) %>%
  arrange(horizon, .by_group = TRUE) %>%
  mutate(
    ibs_discrete_upto_horizon = vapply(seq_len(n()), function(k) mean(brier_ipcw_mean[seq_len(k)], na.rm = TRUE), numeric(1))
  ) %>%
  ungroup() %>%
  left_join_replacing_columns(horizon_annotation, by = c("dataset_key", "horizon")) %>%
  left_join_replacing_columns(risk_set_support_tbl, by = c("dataset_key", "horizon")) %>%
  left_join_replacing_columns(model_annotation, by = c("model_id", "retained_fit_id")) %>%
  apply_prior_tail_sensitive(has_threshold = FALSE)

cohort_metrics_for_join <- posterior_cohort_yearly %>%
  select(
    dataset_key, model_id, retained_fit_id, structural_model_id, formula_anchor, family_code,
    latency_family, branch, risk_scale, prior_branch, site_prior_family, horizon,
    transition_cif_mean, transition_cif_q025, transition_cif_q50, transition_cif_q975,
    remission_cif_mean, remission_cif_q025, remission_cif_q50, remission_cif_q975,
    all_event_free_mean, all_event_free_q025, all_event_free_q50, all_event_free_q975,
    cohort_mean_cure_fraction_mean, cohort_mean_cure_fraction_q025, cohort_mean_cure_fraction_q50, cohort_mean_cure_fraction_q975,
    cohort_mean_susceptible_fraction_mean, cohort_mean_susceptible_fraction_q025, cohort_mean_susceptible_fraction_q50, cohort_mean_susceptible_fraction_q975,
    mean_uncured_survival_mean, mean_uncured_survival_q025, mean_uncured_survival_q50, mean_uncured_survival_q975,
    mean_uncured_risk_mean, mean_uncured_risk_q025, mean_uncured_risk_q50, mean_uncured_risk_q975,
    mean_transition_hazard_mean, mean_transition_hazard_q025, mean_transition_hazard_q50, mean_transition_hazard_q975,
    DeltaRisk_NC_C_mean
  )

posterior_classification <- posterior_classification_raw %>%
  left_join_replacing_columns(
    cohort_metrics_for_join,
    by = c("dataset_key", "model_id", "retained_fit_id", "structural_model_id", "formula_anchor", "family_code", "latency_family", "branch", "risk_scale", "prior_branch", "site_prior_family", "horizon")
  ) %>%
  left_join_replacing_columns(
    horizon_performance_tbl,
    by = c("dataset_key", "model_id", "retained_fit_id", "structural_model_id", "formula_anchor", "family_code", "latency_family", "branch", "risk_scale", "prior_branch", "site_prior_family", "horizon")
  ) %>%
  arrange(factor(model_id, levels = model_order), horizon, threshold)

# 🔴 Build: long-format prediction, performance, coherence, and figure-source tables ===============================
profile_metric_cols <- grep("^(cure_fraction|susceptible_fraction|median_susceptible_time)_(mean|q025|q50|q975)$", names(posterior_subject_profile), value = TRUE)
yearly_metric_cols <- grep("^(transition_cif|remission_cif|all_event_free|uncured_survival|uncured_risk)_(mean|q025|q50|q975)$", names(posterior_subject_yearly_eval), value = TRUE)
cohort_metric_cols <- grep("^(transition_cif|remission_cif|all_event_free|cohort_mean_cure_fraction|cohort_mean_susceptible_fraction|mean_uncured_survival|mean_uncured_risk|mean_transition_hazard)_(mean|q025|q50|q975)$", names(posterior_cohort_yearly), value = TRUE)
class_metric_cols <- grep("^(positive_rate|FPR|false_positive_burden|false_positive_count|FP100|PPV|TPR|NB)_(mean|q025|q50|q975)$", names(posterior_classification), value = TRUE)

prediction_long <- bind_rows(
  if (nrow_or_zero(posterior_subject_profile) > 0L) {
    pivot_prediction_table(
      posterior_subject_profile %>% mutate(horizon = NA_integer_),
      id_cols = c(
        "dataset_key", "model_id", "retained_fit_id", "structural_model_id", "formula_anchor", "family_code", "latency_family", "site_placement_label",
        "branch", "risk_scale", "prior_branch", "site_prior_family", "horizon", "unique_person_id", "id", "site", "sex_num", "age_exact_entry", "age_s",
        "admissibility_flag", "cure_model_eligibility_flag", "receus_primary_class", "cure_presence_support_flag", "followup_not_contradicted_flag"
      ),
      metric_cols = profile_metric_cols
    ) %>% mutate(prediction_level = "subject_profile")
  },
  if (nrow_or_zero(posterior_subject_yearly_eval) > 0L) {
    pivot_prediction_table(
      posterior_subject_yearly_eval,
      id_cols = c(
        "dataset_key", "model_id", "retained_fit_id", "structural_model_id", "formula_anchor", "family_code", "latency_family", "site_placement_label",
        "branch", "risk_scale", "prior_branch", "site_prior_family", "horizon", "unique_person_id", "id", "site", "sex_num", "age_exact_entry", "age_s",
        "time_year", "status_num", "support_tier", "horizon_evidence_class", "claim_restriction_flag", "prior_tail_sensitive", "admissibility_flag"
      ),
      metric_cols = yearly_metric_cols
    ) %>% mutate(prediction_level = "subject_horizon")
  }
) %>%
  arrange(factor(model_id, levels = model_order), horizon, unique_person_id, metric, summary_stat)

performance_long <- bind_rows(
  if (nrow_or_zero(posterior_cohort_yearly) > 0L) {
    pivot_prediction_table(
      posterior_cohort_yearly %>% mutate(threshold = NA_real_),
      id_cols = c(
        "dataset_key", "model_id", "retained_fit_id", "structural_model_id", "formula_anchor", "family_code", "latency_family",
        "branch", "risk_scale", "prior_branch", "site_prior_family", "horizon", "threshold",
        "support_tier", "horizon_evidence_class", "claim_restriction_flag", "prior_tail_sensitive", "admissibility_flag"
      ),
      metric_cols = cohort_metric_cols
    ) %>% mutate(metric_group = "cohort")
  },
  if (nrow_or_zero(posterior_classification) > 0L) {
    pivot_prediction_table(
      posterior_classification,
      id_cols = c(
        "dataset_key", "model_id", "retained_fit_id", "structural_model_id", "formula_anchor", "family_code", "latency_family",
        "branch", "risk_scale", "prior_branch", "site_prior_family", "horizon", "threshold",
        "support_tier", "horizon_evidence_class", "claim_restriction_flag", "prior_tail_sensitive", "admissibility_flag"
      ),
      metric_cols = class_metric_cols
    ) %>% mutate(metric_group = "classification")
  },
  if (nrow_or_zero(horizon_performance_tbl) > 0L) {
    horizon_performance_tbl %>%
      mutate(threshold = NA_real_) %>%
      select(
        dataset_key, model_id, retained_fit_id, structural_model_id, formula_anchor, family_code, latency_family,
        branch, risk_scale, prior_branch, site_prior_family, horizon, threshold,
        support_tier, horizon_evidence_class, claim_restriction_flag, prior_tail_sensitive, admissibility_flag,
        auc_transition_horizon, discrimination_n_complete,
        calibration_mean_pred, calibration_median_pred, calibration_diff_mean_minus_observed, calibration_abs_diff, calibration_oe_ratio,
        brier_ipcw_mean, brier_ipcw_weight_sum, ibs_discrete_upto_horizon
      ) %>%
      pivot_longer(
        cols = c(
          auc_transition_horizon, discrimination_n_complete,
          calibration_mean_pred, calibration_median_pred, calibration_diff_mean_minus_observed, calibration_abs_diff, calibration_oe_ratio,
          brier_ipcw_mean, brier_ipcw_weight_sum, ibs_discrete_upto_horizon
        ),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(summary_stat = "value", metric_group = "performance")
  }
) %>%
  arrange(factor(model_id, levels = model_order), horizon, threshold, metric, summary_stat)

stage8b_coherence_check <- posterior_cohort_yearly %>%
  transmute(
    dataset_key,
    model_id,
    retained_fit_id,
    structural_model_id,
    formula_anchor,
    family_code,
    latency_family,
    branch,
    risk_scale,
    prior_branch,
    site_prior_family,
    horizon,
    transition_cif_mean,
    remission_cif_mean,
    all_event_free_mean,
    decomposition_sum_mean = transition_cif_mean + remission_cif_mean + all_event_free_mean,
    decomposition_abs_error_mean = abs(transition_cif_mean + remission_cif_mean + all_event_free_mean - 1),
    coherence_flag = abs(transition_cif_mean + remission_cif_mean + all_event_free_mean - 1) <= 0.02,
    support_tier,
    horizon_evidence_class,
    claim_restriction_flag,
    prior_tail_sensitive,
    admissibility_flag
  )

horizon_support_panel <- posterior_cohort_yearly %>%
  select(
    dataset_key, model_id, retained_fit_id, structural_model_id, formula_anchor, family_code,
    prior_branch, site_prior_family, horizon,
    support_tier, horizon_evidence_class, claim_restriction_flag, prior_tail_sensitive, admissibility_flag,
    n_total, n_transition_by_horizon, n_remission_by_horizon, n_right_censoring_by_horizon,
    n_observed_at_or_beyond_horizon, risk_set_fraction,
    transition_cif_mean, cohort_mean_cure_fraction_mean
  )

anchor_vs_neutral_delta_panel <- anchor_vs_neutral_delta %>%
  mutate(
    threshold_key = ifelse(is.na(threshold), "__NA__", sprintf("%.10f", threshold))
  ) %>%
  left_join(prior_tail_sensitive_lookup, by = c("dataset_key", "retained_fit_id", "horizon", "threshold_key")) %>%
  mutate(prior_tail_sensitive = coalesce(prior_tail_sensitive, FALSE)) %>%
  select(-threshold_key) %>%
  select(
    dataset_key, branch, risk_scale, retained_fit_id, structural_model_id, formula_anchor, family_code, site_prior_family,
    horizon, threshold, support_tier, claim_restriction_flag, prior_tail_sensitive, delta_field, delta_anchor_minus_neutral
  ) %>%
  tidyr::pivot_wider(names_from = delta_field, values_from = delta_anchor_minus_neutral)

stage8a_vs_stage8b_delta_panel <- delta_vs_stage8a %>%
  select(
    dataset_key, model_id, branch, risk_scale, retained_fit_id, structural_model_id, formula_anchor, family_code, prior_branch, site_prior_family,
    horizon, threshold, support_tier, horizon_evidence_class, claim_restriction_flag, interpretation_note, delta_field, delta_8B_minus_8A
  ) %>%
  tidyr::pivot_wider(names_from = delta_field, values_from = delta_8B_minus_8A)

incidence_anchor_update_panel <- incidence_anchor_update
uncured_only_decomposition_panel <- uncured_supporting_decomposition

figure_source_tables <- bind_rows(
  horizon_support_panel %>% mutate(panel_name = "horizon_support_panel"),
  anchor_vs_neutral_delta_panel %>% mutate(panel_name = "anchor_vs_neutral_delta_panel"),
  stage8a_vs_stage8b_delta_panel %>% mutate(panel_name = "8A_vs_8B_delta_panel"),
  incidence_anchor_update_panel %>% mutate(panel_name = "incidence_anchor_update_panel"),
  uncured_only_decomposition_panel %>% mutate(panel_name = "uncured_only_decomposition_panel")
) %>%
  relocate(panel_name)

admissible_models_export <- model_registry %>%
  filter(admissibility_flag %in% TRUE) %>%
  arrange(factor(dataset_key, levels = c("PNU", "SNU", "merged")), structural_model_id, family_code, prior_branch, site_prior_family)

n_total_models_stage8b <- nrow(model_registry)
n_admissible_models_stage8b <- nrow(admissible_models_export)
n_inadmissible_models_stage8b <- n_total_models_stage8b - n_admissible_models_stage8b

admissibility_summary_export <- bind_rows(
  tibble(
    summary_dimension = "overall",
    summary_level = "admissible_models",
    n_models = n_admissible_models_stage8b
  ),
  admissible_models_export %>%
    dplyr::count(dataset_key, name = "n_models") %>%
    transmute(summary_dimension = "dataset_key", summary_level = dataset_key, n_models = n_models),
  admissible_models_export %>%
    dplyr::count(family_code, name = "n_models") %>%
    transmute(summary_dimension = "family_code", summary_level = family_code, n_models = n_models),
  admissible_models_export %>%
    dplyr::count(prior_branch, name = "n_models") %>%
    transmute(summary_dimension = "prior_branch", summary_level = prior_branch, n_models = n_models),
  admissible_models_export %>%
    dplyr::count(site_prior_family, name = "n_models") %>%
    transmute(summary_dimension = "site_prior_family", summary_level = site_prior_family, n_models = n_models),
  admissible_models_export %>%
    dplyr::count(formula_anchor, name = "n_models") %>%
    transmute(summary_dimension = "formula_anchor", summary_level = formula_anchor, n_models = n_models),
  admissible_models_export %>%
    dplyr::count(structural_model_id, name = "n_models") %>%
    transmute(summary_dimension = "structural_model_id", summary_level = structural_model_id, n_models = n_models),
  admissible_models_export %>%
    dplyr::count(site_placement_label, name = "n_models") %>%
    transmute(summary_dimension = "site_placement_label", summary_level = site_placement_label, n_models = n_models),
  admissible_models_export %>%
    dplyr::count(dataset_key, family_code, name = "n_models") %>%
    transmute(summary_dimension = "dataset_key__family_code", summary_level = paste(dataset_key, family_code, sep = "/"), n_models = n_models),
  admissible_models_export %>%
    dplyr::count(dataset_key, prior_branch, name = "n_models") %>%
    transmute(summary_dimension = "dataset_key__prior_branch", summary_level = paste(dataset_key, prior_branch, sep = "/"), n_models = n_models),
  admissible_models_export %>%
    dplyr::count(family_code, prior_branch, name = "n_models") %>%
    transmute(summary_dimension = "family_code__prior_branch", summary_level = paste(family_code, prior_branch, sep = "/"), n_models = n_models)
) %>%
  mutate(
    n_total_models = n_total_models_stage8b,
    n_admissible_models = n_admissible_models_stage8b,
    share_of_all_models = ifelse(n_total_models > 0L, n_models / n_total_models, NA_real_),
    share_of_admissible_models = ifelse(n_admissible_models > 0L, n_models / n_admissible_models, NA_real_)
  ) %>%
  arrange(summary_dimension, desc(n_models), summary_level)

inadmissible_models_export <- model_registry %>%
  filter(!(admissibility_flag %in% TRUE)) %>%
  mutate(
    inadmissibility_reason = ifelse(trimws(admissibility_reasons) == "", "unspecified", as.character(admissibility_reasons))
  ) %>%
  arrange(desc(inadmissibility_reason), factor(dataset_key, levels = c("PNU", "SNU", "merged")), structural_model_id, family_code, prior_branch, site_prior_family)

inadmissibility_summary_export <- bind_rows(
  tibble(
    summary_dimension = "overall",
    summary_level = "inadmissible_models",
    n_models = n_inadmissible_models_stage8b
  ),
  inadmissible_models_export %>%
    dplyr::count(inadmissibility_reason, name = "n_models") %>%
    transmute(summary_dimension = "inadmissibility_reason", summary_level = inadmissibility_reason, n_models = n_models),
  inadmissible_models_export %>%
    dplyr::count(dataset_key, name = "n_models") %>%
    transmute(summary_dimension = "dataset_key", summary_level = dataset_key, n_models = n_models),
  inadmissible_models_export %>%
    dplyr::count(family_code, name = "n_models") %>%
    transmute(summary_dimension = "family_code", summary_level = family_code, n_models = n_models),
  inadmissible_models_export %>%
    dplyr::count(prior_branch, name = "n_models") %>%
    transmute(summary_dimension = "prior_branch", summary_level = prior_branch, n_models = n_models),
  inadmissible_models_export %>%
    dplyr::count(site_placement_label, name = "n_models") %>%
    transmute(summary_dimension = "site_placement_label", summary_level = site_placement_label, n_models = n_models),
  inadmissible_models_export %>%
    dplyr::count(dataset_key, inadmissibility_reason, name = "n_models") %>%
    transmute(summary_dimension = "dataset_key__inadmissibility_reason", summary_level = paste(dataset_key, inadmissibility_reason, sep = "/"), n_models = n_models),
  inadmissible_models_export %>%
    dplyr::count(family_code, inadmissibility_reason, name = "n_models") %>%
    transmute(summary_dimension = "family_code__inadmissibility_reason", summary_level = paste(family_code, inadmissibility_reason, sep = "/"), n_models = n_models),
  inadmissible_models_export %>%
    dplyr::count(prior_branch, inadmissibility_reason, name = "n_models") %>%
    transmute(summary_dimension = "prior_branch__inadmissibility_reason", summary_level = paste(prior_branch, inadmissibility_reason, sep = "/"), n_models = n_models)
) %>%
  mutate(
    n_total_models = n_total_models_stage8b,
    n_admissible_models = n_admissible_models_stage8b,
    n_inadmissible_models = n_inadmissible_models_stage8b,
    share_of_all_models = ifelse(n_total_models > 0L, n_models / n_total_models, NA_real_),
    share_of_inadmissible_models = ifelse(n_inadmissible_models > 0L, n_models / n_inadmissible_models, NA_real_)
  ) %>%
  arrange(summary_dimension, desc(n_models), summary_level)

# 🔴 Check: output audit and visualization products ===============================
output_audit <- bind_rows(
  tibble(check_name = "model_registry_row_count", status = ifelse(nrow(model_registry) == nrow(model_grid), "pass", "fail"), observed_value = as.character(nrow(model_registry)), expected_value = as.character(nrow(model_grid)), detail = "One model-registry row should exist per Stage 8B model."),
  tibble(check_name = "posterior_cohort_yearly_row_count", status = ifelse(nrow(posterior_cohort_yearly) == sum(model_registry$admissibility_flag %in% TRUE) * length(horizons_year), "pass", "fail"), observed_value = as.character(nrow(posterior_cohort_yearly)), expected_value = as.character(sum(model_registry$admissibility_flag %in% TRUE) * length(horizons_year)), detail = "Admissible fits should contribute one cohort row per horizon."),
  tibble(check_name = "posterior_classification_row_count", status = ifelse(nrow(posterior_classification) == sum(model_registry$admissibility_flag %in% TRUE) * length(horizons_year) * length(risk_thresholds), "pass", "fail"), observed_value = as.character(nrow(posterior_classification)), expected_value = as.character(sum(model_registry$admissibility_flag %in% TRUE) * length(horizons_year) * length(risk_thresholds)), detail = "Admissible fits should contribute one classification row per horizon-threshold pair."),
  tibble(check_name = "prediction_long_nonempty", status = ifelse(nrow(prediction_long) > 0L, "pass", "fail"), observed_value = as.character(nrow(prediction_long)), expected_value = ">0", detail = "Long-format posterior prediction table should be populated."),
  tibble(check_name = "performance_long_nonempty", status = ifelse(nrow(performance_long) > 0L, "pass", "fail"), observed_value = as.character(nrow(performance_long)), expected_value = ">0", detail = "Long-format performance table should be populated."),
  tibble(check_name = "figure_source_tables_nonempty", status = ifelse(nrow(figure_source_tables) > 0L, "pass", "fail"), observed_value = as.character(nrow(figure_source_tables)), expected_value = ">0", detail = "Combined figure-source tables should be populated."),
  tibble(check_name = "coherence_check_nonempty", status = ifelse(nrow(stage8b_coherence_check) > 0L, "pass", "fail"), observed_value = as.character(nrow(stage8b_coherence_check)), expected_value = ">0", detail = "Stage 8B coherence-check table should be populated."),
  tibble(check_name = "admissible_models_export_nonempty", status = ifelse(nrow(admissible_models_export) > 0L, "pass", "warn"), observed_value = as.character(nrow(admissible_models_export)), expected_value = ">0", detail = "Admissible-model export should be populated whenever any Stage 8B model passes admissibility."),
  tibble(check_name = "admissibility_summary_export_nonempty", status = ifelse(nrow(admissibility_summary_export) > 0L, "pass", "warn"), observed_value = as.character(nrow(admissibility_summary_export)), expected_value = ">0", detail = "Admissibility summary export should be populated whenever Stage 8B model registry exists."),
  tibble(check_name = "inadmissibility_summary_export_nonempty", status = ifelse(nrow(inadmissibility_summary_export) > 0L, "pass", "warn"), observed_value = as.character(nrow(inadmissibility_summary_export)), expected_value = ">0", detail = "Inadmissibility summary export should be populated whenever any Stage 8B model fails admissibility."),
  tibble(check_name = "anchor_vs_neutral_delta_nonempty", status = ifelse(nrow(anchor_vs_neutral_delta) > 0L, "pass", "warn"), observed_value = as.character(nrow(anchor_vs_neutral_delta)), expected_value = ">0", detail = "Anchor-versus-neutral delta table should be populated when both prior branches were fitted."),
  tibble(check_name = "delta_vs_stage8a_nonempty", status = ifelse(nrow(delta_vs_stage8a) > 0L, "pass", "warn"), observed_value = as.character(nrow(delta_vs_stage8a)), expected_value = ">0", detail = "Stage 8B-versus-Stage 8A delta table is populated only when matching Stage 8A outputs were found.")
)

pdf_ok <- tryCatch(
  {
    safe_generate_stage8b_visuals(
      posterior_cohort_yearly = posterior_cohort_yearly %>% filter(admissibility_flag %in% TRUE),
      posterior_classification = posterior_classification %>% filter(admissibility_flag %in% TRUE),
      ppc_summary = ppc_summary %>% filter(admissibility_flag %in% TRUE),
      final_pdf_path = stage8b_diagnostic_pdf_file,
      transition_png_path = stage8b_plot_transition_png_file,
      remission_png_path = stage8b_plot_remission_png_file,
      nb_png_path = stage8b_plot_nb_png_file,
      ppc_png_path = stage8b_plot_ppc_png_file
    )
    TRUE
  },
  error = function(e) {
    warning(paste0("Stage 8B diagnostic PDF generation failed: ", conditionMessage(e)), call. = FALSE)
    FALSE
  }
)

output_audit <- bind_rows(
  output_audit,
  tibble(check_name = "diagnostic_pdf_exists", status = ifelse(pdf_file_is_usable(stage8b_diagnostic_pdf_file), "pass", "fail"), observed_value = ifelse(pdf_file_is_usable(stage8b_diagnostic_pdf_file), "TRUE", "FALSE"), expected_value = "TRUE", detail = "Diagnostic PDF should exist after Stage 8B export.")
)

# 🔴 Record: metadata registry and export manifest ===============================
metadata_registry <- tibble::tribble(
  ~metadata_group, ~metadata_name, ~metadata_value,
  "stage", "stage_name", "Stage 8B remission-sensitive Bayesian competing-risk cure extension",
  "stage", "branch", "Stage8B",
  "stage", "risk_scale", "transition_cif_competing",
  "stage", "main_stage8_branch_reference", "Stage8A transition-only cure branch remains primary",
  "inputs", "stage1_export_path", stage1_export_path,
  "inputs", "stage6_screening_flag_csv", stage6_screening_flag_csv,
  "inputs", "stage8a_export_path", stage8a_export_path,
  "thresholds", "common_horizons_year", paste(horizons_year, collapse = ","),
  "thresholds", "risk_thresholds", paste(format(risk_thresholds, trim = TRUE, scientific = FALSE), collapse = ","),
  "model", "transition_families", paste(c("E", "W", "LN", "LL"), collapse = "|"),
  "model", "prior_branches", paste(unique(model_grid$prior_branch), collapse = "|"),
  "model", "site_prior_families", paste(sort(unique(model_grid$site_prior_family)), collapse = "|"),
  "model", "remission_cut_years", paste(remission_cut_years, collapse = "|"),
  "model", "remission_model", "piecewise_exponential",
  "stan", "stan_chains", as.character(stan_chains),
  "stan", "stan_iter", as.character(stan_iter),
  "stan", "stan_warmup", as.character(stan_warmup),
  "stan", "posterior_prediction_draws", as.character(posterior_prediction_draws),
  "prior", "prior_materiality_risk", as.character(prior_materiality_risk),
  "prior", "prior_materiality_false_positive_burden", as.character(prior_materiality_false_positive_burden),
  "prior", "prior_materiality_nb", as.character(prior_materiality_nb),
  "performance", "discrimination_method", "complete_case_binary_auc_on_horizon-usable_subjects",
  "performance", "calibration_method", "mean_predicted_transition_CIF_vs_IPCW_observed_transition_CIF",
  "performance", "brier_method", "simple_IPCW_Brier_for_transition_with_remission_as_competing_event",
  "performance", "ibs_method", "discrete_mean_of_annual_Brier_values_up_to_horizon",
  "interpretation", "cure_fraction_meaning", "non_susceptibility_to_transition",
  "interpretation", "site_effect_rule", "site_terms_are_structural_context_proxies_not_causal_treatment_effects",
  "interpretation", "cross_scale_rule", "Stage8B does not replace Stage8A; cross-scale deltas are remission-aware change summaries only",
  "implementation", "script_version", "stage8b_v2_complete",
  "implementation", "figure_source_contract", "single_combined_csv_with_panel_name_field"
)

export_manifest <- tibble(
  file_name = c(
    "bayes_stage8b_model_registry.csv",
    "bayes_stage8b_admissible_models.csv",
    "bayes_stage8b_admissibility_summary.csv",
    "bayes_stage8b_inadmissibility_summary.csv",
    "bayes_stage8b_coefficient_summary.csv",
    "bayes_stage8b_diagnostics_parameter_level.csv",
    "bayes_stage8b_posterior_prediction_long.csv.gz",
    "bayes_stage8b_posterior_cohort_yearly.csv",
    "bayes_stage8b_posterior_classification.csv",
    "bayes_stage8b_performance_long.csv.gz",
    "bayes_stage8b_prior_predictive_summary.csv",
    "bayes_stage8b_hazard_plausibility.csv",
    "bayes_stage8b_uncured_decomposition.csv",
    "bayes_stage8b_anchor_vs_neutral_delta.csv",
    "bayes_stage8b_incidence_anchor_update.csv",
    "bayes_stage8b_stage8a_vs_stage8b_delta.csv",
    "bayes_stage8b_figure_source_tables.csv.gz",
    "bayes_stage8b_coherence_check.csv",
    "bayes_stage8b_output_audit.csv",
    "bayes_stage8b_metadata_registry.csv",
    "bayes_stage8b_diagnostic_plots.pdf",
    "bayes_stage8b_bundle.rds",
    "bayes_stage8b_export_manifest.csv"
  ),
  object_name = c(
    "model_registry",
    "admissible_models_export",
    "admissibility_summary_export",
    "inadmissibility_summary_export",
    "coefficient_summary",
    "diagnostics_parameter_level",
    "prediction_long",
    "posterior_cohort_yearly",
    "posterior_classification",
    "performance_long",
    "prior_predictive_summary",
    "hazard_shape_plausibility",
    "uncured_supporting_decomposition",
    "anchor_vs_neutral_delta",
    "incidence_anchor_update",
    "delta_vs_stage8a",
    "figure_source_tables",
    "stage8b_coherence_check",
    "output_audit",
    "metadata_registry",
    "diagnostic_pdf",
    "stage8b_bundle",
    "export_manifest"
  ),
  description = c(
    "Stage 8B model-level registry with admissibility, diagnostics, and Stage 6 carry-forward fields",
    "Row-level export of Stage 8B models that passed admissibility",
    "Admissible-model count summary across dataset, family, prior, formula, and structure dimensions",
    "Inadmissible-model count summary across failure reasons and key design dimensions",
    "Parameter posterior summaries for all Stage 8B fits",
    "Parameter-level convergence diagnostics for Stage 8B fits",
    "Mandatory long-format posterior prediction table as source of truth",
    "Cohort-by-horizon posterior Stage 8B summaries on the transition CIF scale",
    "Horizon-level performance and threshold-classification table enriched with discrimination, calibration, Brier, and IBS summaries",
    "Long-format performance and clinical-usefulness table",
    "Prior predictive summaries for Stage 8B fits",
    "Hazard-shape plausibility table on the 1-10 year annual grid",
    "Cure-model-only supporting decomposition table with uncured-only annual summaries",
    "Mandatory anchor-informed versus neutral prior delta table",
    "Mandatory prior-to-posterior incidence-shape update table",
    "Mandatory remission-aware Stage 8B versus Stage 8A delta table",
    "Combined figure-ready source tables identified by panel_name",
    "Stage 8B coherence-check table for transition CIF + remission CIF + all-event-free survival",
    "Automated Stage 8B output audit",
    "Stage 8B metadata registry aligned to the revised integrated master specification",
    "Diagnostic PDF generated from exported Stage 8B summary tables",
    "Reusable Stage 8B bundle with key outputs and configuration",
    "Manifest of Stage 8B exported files"
  ),
  file_path = file.path(export_path, c(
    "bayes_stage8b_model_registry.csv",
    "bayes_stage8b_admissible_models.csv",
    "bayes_stage8b_admissibility_summary.csv",
    "bayes_stage8b_inadmissibility_summary.csv",
    "bayes_stage8b_coefficient_summary.csv",
    "bayes_stage8b_diagnostics_parameter_level.csv",
    "bayes_stage8b_posterior_prediction_long.csv.gz",
    "bayes_stage8b_posterior_cohort_yearly.csv",
    "bayes_stage8b_posterior_classification.csv",
    "bayes_stage8b_performance_long.csv.gz",
    "bayes_stage8b_prior_predictive_summary.csv",
    "bayes_stage8b_hazard_plausibility.csv",
    "bayes_stage8b_uncured_decomposition.csv",
    "bayes_stage8b_anchor_vs_neutral_delta.csv",
    "bayes_stage8b_incidence_anchor_update.csv",
    "bayes_stage8b_stage8a_vs_stage8b_delta.csv",
    "bayes_stage8b_figure_source_tables.csv.gz",
    "bayes_stage8b_coherence_check.csv",
    "bayes_stage8b_output_audit.csv",
    "bayes_stage8b_metadata_registry.csv",
    "bayes_stage8b_diagnostic_plots.pdf",
    "bayes_stage8b_bundle.rds",
    "bayes_stage8b_export_manifest.csv"
  ))
)

# 🔴 Export: final Stage 8B deliverables ===============================
write_csv_preserve_schema(model_registry, file.path(export_path, "bayes_stage8b_model_registry.csv"))
write_csv_preserve_schema(admissible_models_export, file.path(export_path, "bayes_stage8b_admissible_models.csv"))
write_csv_preserve_schema(admissibility_summary_export, file.path(export_path, "bayes_stage8b_admissibility_summary.csv"))
write_csv_preserve_schema(inadmissibility_summary_export, file.path(export_path, "bayes_stage8b_inadmissibility_summary.csv"))
write_csv_preserve_schema(coefficient_summary, file.path(export_path, "bayes_stage8b_coefficient_summary.csv"))
write_csv_preserve_schema(diagnostics_parameter_level, file.path(export_path, "bayes_stage8b_diagnostics_parameter_level.csv"))
write_csv_preserve_schema(prediction_long, file.path(export_path, "bayes_stage8b_posterior_prediction_long.csv.gz"))
write_csv_preserve_schema(posterior_cohort_yearly, file.path(export_path, "bayes_stage8b_posterior_cohort_yearly.csv"))
write_csv_preserve_schema(posterior_classification, file.path(export_path, "bayes_stage8b_posterior_classification.csv"))
write_csv_preserve_schema(performance_long, file.path(export_path, "bayes_stage8b_performance_long.csv.gz"))
write_csv_preserve_schema(prior_predictive_summary, file.path(export_path, "bayes_stage8b_prior_predictive_summary.csv"))
write_csv_preserve_schema(hazard_shape_plausibility, file.path(export_path, "bayes_stage8b_hazard_plausibility.csv"))
write_csv_preserve_schema(uncured_supporting_decomposition, file.path(export_path, "bayes_stage8b_uncured_decomposition.csv"))
write_csv_preserve_schema(anchor_vs_neutral_delta, file.path(export_path, "bayes_stage8b_anchor_vs_neutral_delta.csv"))
write_csv_preserve_schema(incidence_anchor_update, file.path(export_path, "bayes_stage8b_incidence_anchor_update.csv"))
write_csv_preserve_schema(delta_vs_stage8a, file.path(export_path, "bayes_stage8b_stage8a_vs_stage8b_delta.csv"))
write_csv_preserve_schema(figure_source_tables, file.path(export_path, "bayes_stage8b_figure_source_tables.csv.gz"))
write_csv_preserve_schema(stage8b_coherence_check, file.path(export_path, "bayes_stage8b_coherence_check.csv"))
write_csv_preserve_schema(output_audit, file.path(export_path, "bayes_stage8b_output_audit.csv"))
write_csv_preserve_schema(metadata_registry, file.path(export_path, "bayes_stage8b_metadata_registry.csv"))

saveRDS(
  list(
    stage = "Stage8B",
    version = "stage8b_v2_complete",
    created_at = as.character(Sys.time()),
    session_info = utils::sessionInfo(),
    config = list(
      stage1_export_path = stage1_export_path,
      stage6_export_path = stage6_export_path,
      stage8a_export_path = stage8a_export_path,
      export_path = export_path,
      pnu_site_label = pnu_site_label,
      snu_site_label = snu_site_label,
      horizons_year = horizons_year,
      risk_thresholds = risk_thresholds,
      remission_cut_years = remission_cut_years,
      fit_prior_branches = fit_prior_branches,
      fit_site_prior_families = fit_site_prior_families,
      stan_chains = stan_chains,
      stan_iter = stan_iter,
      stan_warmup = stan_warmup,
      posterior_prediction_draws = posterior_prediction_draws
    ),
    outputs = list(
      model_registry = model_registry,
      admissible_models_export = admissible_models_export,
      admissibility_summary_export = admissibility_summary_export,
      inadmissibility_summary_export = inadmissibility_summary_export,
      coefficient_summary = coefficient_summary,
      diagnostics_parameter_level = diagnostics_parameter_level,
      prediction_long = prediction_long,
      posterior_subject_profile = posterior_subject_profile,
      posterior_subject_yearly = posterior_subject_yearly_eval,
      posterior_cohort_yearly = posterior_cohort_yearly,
      posterior_classification = posterior_classification,
      performance_long = performance_long,
      ppc_summary = ppc_summary,
      prior_predictive_summary = prior_predictive_summary,
      hazard_shape_plausibility = hazard_shape_plausibility,
      uncured_supporting_decomposition = uncured_supporting_decomposition,
      anchor_vs_neutral_delta = anchor_vs_neutral_delta,
      incidence_anchor_update = incidence_anchor_update,
      stage8a_vs_stage8b_delta = delta_vs_stage8a,
      figure_source_tables = figure_source_tables,
      coherence_check = stage8b_coherence_check,
      output_audit = output_audit,
      metadata_registry = metadata_registry
    )
  ),
  file.path(export_path, "bayes_stage8b_bundle.rds")
)

write_csv_preserve_schema(export_manifest, file.path(export_path, "bayes_stage8b_export_manifest.csv"))

message(
  "Stage8B completed. Reused models: ",
  n_models_reused,
  "; newly fitted models: ",
  n_models_to_fit,
  "; diagnostic PDF ok: ",
  pdf_ok,
  "."
)
