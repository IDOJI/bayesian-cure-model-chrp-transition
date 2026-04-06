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
block1_root_default <- "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/1.Block1"
dropbox_root <- Sys.getenv(
  "CHRPMIXTURE_DROPBOX_ROOT",
  unset = "/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model"
)
source_data_file <- Sys.getenv(
  "MIXTURE_CURE_SOURCE_DATA_FILE",
  unset = file.path(
    dropbox_root,
    "0.Data",
    "2.Preprocessed data",
    "Preprocessed_Merged_PNUH_SNUH_Data.csv"
  )
)
modeling_root <- Sys.getenv(
  "MIXTURE_CURE_MODELING_RESULTS_ROOT",
  unset = block1_root_default
)
results_root <- file.path(repo_root, "3.Results files")
results_integration_root <- Sys.getenv(
  "MIXTURE_CURE_RESULTS_INTEGRATION_ROOT",
  unset = file.path(block1_root_default, "5.Results Integration")
)
export_path <- Sys.getenv(
  "MIXTURE_CURE_FRACTION_EXPORT_PATH",
  unset = file.path(results_integration_root, "mixture_cure_fraction_tables")
)

stage7_fit_candidates <- c(
  file.path(modeling_root, "3.MLE Mixture Cure", "mle_mixture_cure_lognormal_fitted_objects.rds"),
  file.path(
    modeling_root,
    "3.MLE Mixture Cure",
    "stage7_simple_lognormal_mixture_cure_fitted_objects.rds"
  ),
  file.path(results_root, "stage7_simple_lognormal_mixture_cure_fitted_objects.rds")
)
stage7_summary_candidates <- c(
  file.path(modeling_root, "3.MLE Mixture Cure", "mle_mixture_cure_lognormal_horizon_summary.csv"),
  file.path(
    modeling_root,
    "3.MLE Mixture Cure",
    "stage7_simple_lognormal_mixture_cure_horizon_summary.csv"
  ),
  file.path(results_root, "stage7_simple_lognormal_mixture_cure_horizon_summary.csv")
)
stage8_registry_candidates <- c(
  file.path(
    modeling_root,
    "4.Bayesian Mixture Cure",
    "bayesian_mixture_cure_model_registry.csv"
  ),
  file.path(results_root, "stage8A_simple_bayesian_transition_only_model_registry.csv")
)
stage8_summary_candidates <- c(
  file.path(
    modeling_root,
    "4.Bayesian Mixture Cure",
    "bayesian_mixture_cure_horizon_summary.csv"
  ),
  file.path(results_root, "stage8A_simple_bayesian_transition_only_horizon_summary.csv")
)

pnu_site_label <- Sys.getenv("MIXTURE_CURE_PNU_SITE_LABEL", unset = "PNU")
snu_site_label <- Sys.getenv("MIXTURE_CURE_SNU_SITE_LABEL", unset = "SNU")
time_origin_epsilon_year <- 1e-08
round_digits <- as.integer(Sys.getenv("MIXTURE_CURE_FRACTION_ROUND_DIGITS", unset = "6"))
stage7_ci_draws <- as.integer(Sys.getenv("MIXTURE_CURE_FRACTION_STAGE7_CI_DRAWS", unset = "4000"))

# Load Packages -----------------------------------------------------------
required_packages <- c("readr", "dplyr", "tidyr", "tibble", "rstan")
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
  library(rstan)
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

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

clip_prob <- function(x, eps = 1e-12) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

round_value <- function(x) {
  ifelse(is.na(x), NA_real_, round(as.numeric(x), digits = round_digits))
}

cleanup_existing_fraction_outputs <- function(export_dir) {
  stale_files <- list.files(
    path = export_dir,
    pattern = "\\.csv$",
    full.names = TRUE
  )

  if (length(stale_files) > 0L) {
    file.remove(stale_files)
  }

  invisible(stale_files)
}

write_fraction_tables_readme <- function(export_dir) {
  readme_lines <- c(
    sprintf("`%s` 폴더의 파일들은 아래 의미입니다.", normalizePath(export_dir, winslash = "/", mustWork = FALSE)),
    "",
    "`input_file_registry.csv`",
    "- 이 통합 작업에 실제로 사용된 입력 파일 경로 목록입니다.",
    "- `source_data`, `mle_mixture_cure_fit`, `mle_mixture_cure_summary`, `bayesian_registry`, `bayesian_summary`가 기록됩니다.",
    "",
    "`frequentist_mixture_cure_fraction_summary.csv`",
    "- Frequentist mixture cure 모델만 따로 모은 요약표입니다.",
    "- 각 dataset variant별 `susceptible_fraction_estimate`, `cure_fraction_estimate`, 그리고 가능한 경우 95% CI가 들어 있습니다.",
    "",
    "`bayesian_mixture_cure_fraction_summary.csv`",
    "- Bayesian mixture cure 모델만 따로 모은 요약표입니다.",
    "- 각 dataset variant별 point estimate와 함께 95% CrI가 들어 있습니다.",
    "",
    "`integrated_mixture_cure_fraction_summary.csv`",
    "- Frequentist와 Bayesian 결과를 한 long-format 테이블로 합친 통합 요약입니다.",
    "- dataset, model family, estimation method, interval type, fit source, point estimate, 95% interval을 한 번에 비교할 때 쓰기 좋습니다.",
    "",
    "`pnu_cure_fraction_comparison_table.csv`",
    "- PNU 데이터셋에 대한 비교표입니다.",
    "- 행은 모델, 열은 추정치와 95% CI/CrI 요약으로 구성됩니다.",
    "",
    "`snu_cure_fraction_comparison_table.csv`",
    "- SNU 데이터셋에 대한 비교표입니다.",
    "- PNU 표와 같은 형식으로 frequentist와 Bayesian cure fraction을 바로 비교할 수 있습니다.",
    "",
    "`merged_cure_fraction_comparison_table.csv`",
    "- merged 계열 비교표입니다.",
    "- `merged`와 `merged(site-adj)` 모델이 함께 들어가 있어서 site adjustment 유무까지 같이 비교할 수 있습니다.",
    "",
    "`output_table_registry.csv`",
    "- 어떤 comparison table이 생성됐는지 파일 경로를 기록한 목록입니다.",
    "- 후속 스크립트에서 자동으로 결과 파일을 읽어갈 때 유용합니다.",
    "",
    "참고로 비교표들에서",
    "- `susceptible_fraction_estimate`는 uncured fraction",
    "- `cure_fraction_estimate`는 cured fraction",
    "- `interval_type`은 `95% CI` 또는 `95% CrI`",
    "를 뜻합니다."
  )

  writeLines(readme_lines, con = file.path(export_dir, "readme.md"), useBytes = TRUE)
  invisible(file.path(export_dir, "readme.md"))
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

standardize_known_site_labels <- function(df, pnu_label, snu_label) {
  if (!("site" %in% names(df))) {
    stop("Input data must contain a `site` column.", call. = FALSE)
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

prepare_analysis_dataset <- function(df, dataset_name, site_mode = c("single", "merged")) {
  site_mode <- match.arg(site_mode)

  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("[%s] Missing required columns: %s", dataset_name, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  out <- tibble::as_tibble(df) %>%
    mutate(
      id = trimws(as.character(id)),
      site = trimws(as.character(site)),
      sex_num = as.integer(safe_numeric(sex_num)),
      age_exact_entry = safe_numeric(age_exact_entry),
      days_followup = safe_numeric(days_followup),
      status_num = as.integer(safe_numeric(status_num))
    )

  if (nrow(out) == 0L) {
    stop(sprintf("[%s] Dataset has zero rows.", dataset_name), call. = FALSE)
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
    stop(sprintf("[%s] `sex_num` must be coded 0/1.", dataset_name), call. = FALSE)
  }
  if (any(!out$status_num %in% c(0L, 1L, 2L), na.rm = TRUE)) {
    stop(sprintf("[%s] `status_num` must be coded 0/1/2.", dataset_name), call. = FALSE)
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

  age_mean <- mean(out$age_exact_entry)
  age_sd <- stats::sd(out$age_exact_entry)
  if (is.na(age_sd) || age_sd <= 0) {
    stop(sprintf("[%s] `age_exact_entry` must have positive SD.", dataset_name), call. = FALSE)
  }

  out %>%
    mutate(
      unique_person_id = paste(site, id, sep = "_"),
      time_year = pmax(days_followup / 365.25, time_origin_epsilon_year),
      time_year_model = pmax(days_followup / 365.25, time_origin_epsilon_year),
      age_s = (age_exact_entry - age_mean) / (2 * age_sd),
      event_main = as.integer(status_num == 1L),
      censor_main = as.integer(status_num %in% c(0L, 2L)),
      site = factor(as.character(site))
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

  list(
    PNU = prepare_analysis_dataset(
      split_site_dataset(source_df, pnu_label, "PNU"),
      dataset_name = "PNU",
      site_mode = "single"
    ),
    SNU = prepare_analysis_dataset(
      split_site_dataset(source_df, snu_label, "SNU"),
      dataset_name = "SNU",
      site_mode = "single"
    ),
    merged = prepare_analysis_dataset(
      source_df,
      dataset_name = "merged",
      site_mode = "merged"
    )
  )
}

trim_rhs <- function(rhs) {
  rhs <- gsub("\\s+", " ", trimws(as.character(rhs)))
  rhs[is.na(rhs)] <- NA_character_
  rhs
}

make_rhs_formula <- function(rhs) {
  rhs <- trim_rhs(rhs)
  if (is.na(rhs) || rhs == "") {
    stats::as.formula("~ 1")
  } else {
    stats::as.formula(paste("~", rhs))
  }
}

model_matrix_from_rhs <- function(df, rhs) {
  stats::model.matrix(make_rhs_formula(rhs), data = df)
}

ensure_draw_matrix <- function(x, ncol_expected) {
  if (is.null(dim(x))) {
    return(matrix(as.numeric(x), ncol = ncol_expected))
  }

  if (length(dim(x)) != 2L) {
    stop("Posterior draws must be a vector or a 2D matrix.", call. = FALSE)
  }

  matrix(
    as.numeric(x),
    nrow = dim(x)[1L],
    ncol = ncol_expected,
    dimnames = dimnames(x)
  )
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

extract_prediction_draws <- function(fit, K_inc, K_lat, max_draws, seed_value) {
  ext <- rstan::extract(
    fit,
    pars = c("alpha_inc", "beta_inc", "gamma0", "gamma_lat", "log_sigma"),
    permuted = TRUE,
    inc_warmup = FALSE
  )

  total_draws <- length(ext$alpha_inc)
  set.seed(seed_value)
  keep_draw_idx <- if (total_draws <= max_draws) {
    seq_len(total_draws)
  } else {
    sort(sample(seq_len(total_draws), size = max_draws, replace = FALSE))
  }

  list(
    alpha_inc = as.numeric(ext$alpha_inc[keep_draw_idx]),
    beta_inc = ensure_draw_matrix(ext$beta_inc, ncol_expected = K_inc)[keep_draw_idx, , drop = FALSE],
    gamma0 = as.numeric(ext$gamma0[keep_draw_idx]),
    gamma_lat = ensure_draw_matrix(ext$gamma_lat, ncol_expected = K_lat)[keep_draw_idx, , drop = FALSE],
    log_sigma = as.numeric(ext$log_sigma[keep_draw_idx]),
    n_draws = length(keep_draw_idx)
  )
}

build_prior_artifacts <- function() {
  list(
    alpha_gp_vdw = -9.581369553169,
    mu_beta_inc_anchor = c(
      0.419871845822,
      0.907608052926,
      0.586202561451,
      0.466865123863,
      0.037997248763
    ),
    sd_beta_inc_anchor = c(
      0.132789397422,
      0.173731076538,
      0.191221553945,
      0.270393197518,
      0.302838606651
    ),
    alpha_inc_sd = 3.5,
    sd_gamma0 = 2.5,
    sd_gamma_lat = 1.0,
    sd_log_sigma = 0.50
  )
}

make_design_bundle <- function(df, model_row, prior_artifacts, snu_label) {
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
  inc_is_site <- rep(0L, ncol(X_inc))

  mu_beta_inc <- prior_artifacts$mu_beta_inc_anchor
  sd_beta_inc <- prior_artifacts$sd_beta_inc_anchor

  if (isTRUE(model_row$incidence_site_indicator[[1L]])) {
    X_inc <- cbind(X_inc, site_SNU = s_i)
    mu_beta_inc <- c(mu_beta_inc, 0)
    sd_beta_inc <- c(sd_beta_inc, 1)
    inc_is_site <- c(inc_is_site, 1L)
  }

  a_i <- as.numeric(df$age_s)

  X_lat <- switch(
    as.character(model_row$structural_model_id[[1L]]),
    "PNU-L0" = cbind(age_s = a_i, sex_num = z_i),
    "SNU-L0" = cbind(age_s = a_i, sex_num = z_i),
    "MERGED-I0-L0S0" = cbind(age_s = a_i, sex_num = z_i),
    "MERGED-I1-L0S1" = cbind(age_s = a_i, sex_num = z_i, site_SNU = s_i),
    stop("Unknown structural model id: ", model_row$structural_model_id[[1L]], call. = FALSE)
  )

  list(
    X_inc = unclass(as.matrix(X_inc)),
    inc_is_site = as.integer(inc_is_site),
    mu_beta_inc = as.numeric(mu_beta_inc),
    sd_beta_inc = as.numeric(sd_beta_inc),
    alpha_inc_prior_mean = as.numeric(prior_artifacts$alpha_gp_vdw),
    alpha_inc_prior_sd = as.numeric(prior_artifacts$alpha_inc_sd),
    X_lat = unclass(as.matrix(X_lat)),
    lat_is_site = as.integer(colnames(X_lat) == "site_SNU"),
    time = as.numeric(df$time_year),
    event = as.integer(df$event_main)
  )
}

summarize_draws <- function(x) {
  tibble(
    mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    q025 = unname(stats::quantile(x, probs = 0.025, na.rm = TRUE)),
    q50 = unname(stats::quantile(x, probs = 0.50, na.rm = TRUE)),
    q975 = unname(stats::quantile(x, probs = 0.975, na.rm = TRUE))
  )
}

model_row_label <- function(dataset_group, dataset_variant, model_family) {
  base_label <- dplyr::case_when(
    model_family == "Frequentist mixture cure" ~ "Frequentist mixture cure",
    model_family == "Bayesian transition-only cure" ~ "Bayesian transition-only cure",
    TRUE ~ model_family
  )

  if (!identical(dataset_group, "Merged")) {
    return(base_label)
  }

  suffix <- if (identical(dataset_variant, "merged_site_adjusted")) {
    "[merged(site-adj)]"
  } else {
    "[merged]"
  }

  paste(base_label, suffix)
}

model_row_order <- function(dataset_group, dataset_variant, model_family) {
  if (!identical(dataset_group, "Merged")) {
    return(dplyr::case_when(
      model_family == "Frequentist mixture cure" ~ 1,
      model_family == "Bayesian transition-only cure" ~ 2,
      TRUE ~ 999
    ))
  }

  dplyr::case_when(
    model_family == "Frequentist mixture cure" && dataset_variant == "merged" ~ 1,
    model_family == "Frequentist mixture cure" && dataset_variant == "merged_site_adjusted" ~ 2,
    model_family == "Bayesian transition-only cure" && dataset_variant == "merged" ~ 3,
    model_family == "Bayesian transition-only cure" && dataset_variant == "merged_site_adjusted" ~ 4,
    TRUE ~ 999
  )
}

build_comparison_table <- function(summary_df, dataset_group) {
  summary_df %>%
    filter(dataset_group == !!dataset_group) %>%
    mutate(
      model = mapply(
        model_row_label,
        dataset_group = dataset_group,
        dataset_variant = dataset_variant,
        model_family = model_family,
        USE.NAMES = FALSE
      ),
      row_order = mapply(
        model_row_order,
        dataset_group = dataset_group,
        dataset_variant = dataset_variant,
        model_family = model_family,
        USE.NAMES = FALSE
      )
    ) %>%
    arrange(row_order, model) %>%
    transmute(
      model = model,
      estimation_method = estimation_method,
      interval_type = interval_type,
      n_subjects = n_subjects,
      prediction_draws = prediction_draws,
      susceptible_fraction_estimate = round_value(susceptible_fraction_estimate),
      cure_fraction_estimate = round_value(cure_fraction_estimate),
      uncertainty_sd = round_value(cure_fraction_sd),
      interval_lower_95 = round_value(cure_fraction_q025),
      interval_median_50 = round_value(cure_fraction_q50),
      interval_upper_95 = round_value(cure_fraction_q975)
    )
}

write_comparison_table <- function(table_df, export_dir, dataset_group) {
  dataset_stub <- dplyr::case_when(
    dataset_group == "PNU" ~ "pnu",
    dataset_group == "SNU" ~ "snu",
    dataset_group == "Merged" ~ "merged",
    TRUE ~ tolower(dataset_group)
  )

  output_file <- file.path(
    export_dir,
    sprintf("%s_cure_fraction_comparison_table.csv", dataset_stub)
  )

  readr::write_csv(table_df, output_file)
  invisible(output_file)
}

format_summary_for_export <- function(df) {
  df %>%
    mutate(
      susceptible_fraction_estimate = round_value(susceptible_fraction_estimate),
      cure_fraction_estimate = round_value(cure_fraction_estimate),
      uncertainty_sd = round_value(cure_fraction_sd),
      interval_lower_95 = round_value(cure_fraction_q025),
      interval_median_50 = round_value(cure_fraction_q50),
      interval_upper_95 = round_value(cure_fraction_q975)
    ) %>%
    select(
      dataset_variant,
      dataset_group,
      model_family,
      estimation_method,
      interval_type,
      source_fit_file,
      fit_object_key,
      site_adjustment_flag,
      n_subjects,
      prediction_draws,
      susceptible_fraction_estimate,
      cure_fraction_estimate,
      uncertainty_sd,
      interval_lower_95,
      interval_median_50,
      interval_upper_95
    )
}

# Compute Frequentist Cure Fractions -------------------------------------
compute_stage7_cure_fraction_summary <- function(analysis_datasets, fit_path, summary_path = NULL) {
  fit_objects <- readRDS(fit_path)

  if (!is.list(fit_objects)) {
    stop("Stage 7 fitted object file must contain a list.", call. = FALSE)
  }

  stage7_registry <- tibble(
    dataset_variant = c("PNU", "SNU", "merged", "merged_site_adjusted"),
    source_dataset = c("PNU", "SNU", "merged", "merged"),
    fit_key = c("PNU", "SNU", "merged_no_site", "merged_site_adjusted"),
    site_adjustment_flag = c(FALSE, FALSE, FALSE, TRUE)
  )

  output_df <- bind_rows(lapply(seq_len(nrow(stage7_registry)), function(ii) {
    row <- stage7_registry[ii, ]
    fit_key <- as.character(row$fit_key[[1L]])
    fit_obj <- fit_objects[[fit_key]]

    if (is.null(fit_obj)) {
      stop(sprintf("Stage 7 fit object missing for `%s`.", fit_key), call. = FALSE)
    }

    analysis_df <- analysis_datasets[[as.character(row$source_dataset[[1L]])]]
    pred_data <- analysis_df %>%
      mutate(site = factor(as.character(site), levels = fit_obj$site_levels))

    X_inc <- model_matrix_from_rhs(pred_data, fit_obj$incidence_rhs)
    uncured_prob <- clip_prob(plogis(drop(X_inc %*% fit_obj$gamma)))
    cure_prob <- 1 - uncured_prob

    X_lat <- model_matrix_from_rhs(pred_data, fit_obj$latency_rhs)
    y_time <- pmax(as.numeric(pred_data$time_year_model), time_origin_epsilon_year)
    y_event <- as.integer(pred_data$event_main)

    neg_loglik <- function(par) {
      gamma <- par[seq_len(ncol(X_inc))]
      beta <- par[ncol(X_inc) + seq_len(ncol(X_lat))]
      log_sigma <- par[[length(par)]]

      uncured_prob_inner <- clip_prob(plogis(drop(X_inc %*% gamma)))
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

    par_hat <- c(fit_obj$gamma, fit_obj$beta, fit_obj$log_sigma)
    hessian_mat <- tryCatch(
      stats::optimHess(par = par_hat, fn = neg_loglik),
      error = function(e) NULL
    )

    ci_summary <- tibble(
      cure_fraction_sd = NA_real_,
      cure_fraction_q025 = NA_real_,
      cure_fraction_q50 = NA_real_,
      cure_fraction_q975 = NA_real_
    )

    if (!is.null(hessian_mat) && all(is.finite(hessian_mat))) {
      covariance_mat <- tryCatch(
        invert_hessian(hessian_mat),
        error = function(e) NULL
      )

      if (!is.null(covariance_mat) && all(is.finite(covariance_mat))) {
        parameter_draws <- draw_mvn(
          n_draws = stage7_ci_draws,
          mu = par_hat,
          sigma = covariance_mat,
          seed_value = 730000L + ii
        )

        gamma_draws <- parameter_draws[, seq_len(ncol(X_inc)), drop = FALSE]
        eta_draws <- X_inc %*% t(gamma_draws)
        uncured_draws <- clip_prob(plogis(eta_draws))
        if (is.null(dim(uncured_draws))) {
          uncured_draws <- matrix(uncured_draws, nrow = nrow(X_inc), ncol = nrow(gamma_draws))
        }
        cure_fraction_draws <- colMeans(1 - uncured_draws)
        ci_summary <- summarize_draws(cure_fraction_draws) %>%
          transmute(
            cure_fraction_sd = sd,
            cure_fraction_q025 = q025,
            cure_fraction_q50 = q50,
            cure_fraction_q975 = q975
          )
      }
    }

    tibble(
      dataset_variant = as.character(row$dataset_variant[[1L]]),
      dataset_group = dataset_group_from_variant(row$dataset_variant[[1L]]),
      model_family = "Frequentist mixture cure",
      estimation_method = "maximum likelihood",
      interval_type = "95% CI",
      source_fit_file = normalizePath(fit_path, winslash = "/", mustWork = TRUE),
      fit_object_key = fit_key,
      site_adjustment_flag = as.logical(row$site_adjustment_flag[[1L]]),
      n_subjects = nrow(pred_data),
      prediction_draws = as.integer(stage7_ci_draws),
      susceptible_fraction_estimate = mean(uncured_prob, na.rm = TRUE),
      cure_fraction_estimate = mean(cure_prob, na.rm = TRUE),
      cure_fraction_sd = ci_summary$cure_fraction_sd[[1L]],
      cure_fraction_q025 = ci_summary$cure_fraction_q025[[1L]],
      cure_fraction_q50 = ci_summary$cure_fraction_q50[[1L]],
      cure_fraction_q975 = ci_summary$cure_fraction_q975[[1L]]
    ) %>%
      mutate(
        prediction_draws = dplyr::if_else(
          is.na(cure_fraction_q025),
          NA_integer_,
          prediction_draws
        )
      )
  }))

  if (!is.null(summary_path) && file.exists(summary_path)) {
    summary_df <- read_csv_checked(summary_path, "Stage 7 horizon summary") %>%
      mutate(dataset_variant = normalize_dataset_variant(dataset)) %>%
      group_by(dataset_variant) %>%
      summarise(reference_cure_fraction = mean(cure_fraction, na.rm = TRUE), .groups = "drop")

    validation_df <- output_df %>%
      left_join(summary_df, by = "dataset_variant") %>%
      mutate(abs_diff = abs(cure_fraction_estimate - reference_cure_fraction))

    if (any(validation_df$abs_diff > 1e-8, na.rm = TRUE)) {
      stop("Stage 7 cure-fraction recomputation does not match the exported horizon summary.", call. = FALSE)
    }
  }

  output_df
}

# Compute Bayesian Cure Fractions ----------------------------------------
compute_stage8_cure_fraction_summary <- function(
  analysis_datasets,
  registry_path,
  summary_path = NULL,
  snu_label = "SNU"
) {
  if (is.null(summary_path) || !file.exists(summary_path)) {
    stop("Bayesian cure-fraction export requires the horizon summary CSV.", call. = FALSE)
  }

  summary_df <- read_csv_checked(summary_path, "Bayesian horizon summary") %>%
    mutate(dataset_variant = normalize_dataset_variant(dataset))

  if ("prior_branch" %in% names(summary_df)) {
    summary_df <- summary_df %>%
      filter(prior_branch == "anchor_informed")
  }

  if (!("cure_fraction" %in% names(summary_df))) {
    summary_df <- summary_df %>%
      mutate(cure_fraction = 1 - as.numeric(susceptible_fraction))
  }

  registry_df <- if (!is.null(registry_path) && file.exists(registry_path)) {
    read_csv_checked(registry_path, "Bayesian model registry") %>%
      mutate(dataset_variant = normalize_dataset_variant(dataset))
  } else {
    tibble()
  }

  if (nrow(registry_df) > 0L && "prior_branch" %in% names(registry_df)) {
    registry_df <- registry_df %>%
      filter(prior_branch == "anchor_informed")
  }

  registry_lookup_df <- if (nrow(registry_df) > 0L) {
    registry_df %>%
      group_by(dataset_variant) %>%
      summarise(
        source_fit_file = dplyr::first(as.character(fit_rds_path)),
        fit_object_key = dplyr::first(as.character(model_id)),
        prediction_draws = dplyr::first(as.integer(posterior_prediction_draws)),
        .groups = "drop"
      )
  } else {
    tibble(
      dataset_variant = character(),
      source_fit_file = character(),
      fit_object_key = character(),
      prediction_draws = integer()
    )
  }

  summary_df %>%
    group_by(dataset_variant) %>%
    summarise(
      dataset_group = dataset_group_from_variant(dplyr::first(dataset_variant)),
      model_family = "Bayesian transition-only cure",
      estimation_method = "Bayesian posterior",
      interval_type = "95% CrI",
      site_adjustment_flag = dplyr::first(as.logical(site_adjustment_flag)),
      n_subjects = dplyr::first(as.integer(n_subjects)),
      susceptible_fraction_estimate = mean(as.numeric(susceptible_fraction), na.rm = TRUE),
      cure_fraction_estimate = mean(as.numeric(cure_fraction), na.rm = TRUE),
      cure_fraction_sd = mean(as.numeric(cure_fraction_sd), na.rm = TRUE),
      cure_fraction_q025 = mean(as.numeric(cure_fraction_q025), na.rm = TRUE),
      cure_fraction_q50 = mean(as.numeric(cure_fraction_q50), na.rm = TRUE),
      cure_fraction_q975 = mean(as.numeric(cure_fraction_q975), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(registry_lookup_df, by = "dataset_variant") %>%
    mutate(
      source_fit_file = dplyr::coalesce(source_fit_file, summary_path),
      fit_object_key = dplyr::coalesce(fit_object_key, paste0("bayesian_", dataset_variant)),
      prediction_draws = suppressWarnings(as.integer(prediction_draws))
    ) %>%
    select(
      dataset_variant,
      dataset_group,
      model_family,
      estimation_method,
      interval_type,
      source_fit_file,
      fit_object_key,
      site_adjustment_flag,
      n_subjects,
      prediction_draws,
      susceptible_fraction_estimate,
      cure_fraction_estimate,
      cure_fraction_sd,
      cure_fraction_q025,
      cure_fraction_q50,
      cure_fraction_q975
    )
}

# Execute Export ----------------------------------------------------------
resolved_stage7_fit_path <- resolve_existing_path(
  stage7_fit_candidates,
  "Stage 7 fitted object file"
)
resolved_stage7_summary_path <- resolve_existing_path(
  stage7_summary_candidates,
  "Stage 7 horizon summary file"
)
resolved_stage8_registry_path <- resolve_existing_path(
  stage8_registry_candidates,
  "Stage 8 model registry file"
)
resolved_stage8_summary_path <- resolve_existing_path(
  stage8_summary_candidates,
  "Stage 8 horizon summary file"
)

analysis_datasets <- build_analysis_datasets_from_source(
  source_file = source_data_file,
  pnu_label = pnu_site_label,
  snu_label = snu_site_label
)

metadata_df <- tibble(
  input_label = c("source_data", "mle_mixture_cure_fit", "mle_mixture_cure_summary", "bayesian_registry", "bayesian_summary"),
  input_path = c(
    normalizePath(source_data_file, winslash = "/", mustWork = TRUE),
    resolved_stage7_fit_path,
    resolved_stage7_summary_path,
    resolved_stage8_registry_path,
    resolved_stage8_summary_path
  )
)

stage7_summary_df <- compute_stage7_cure_fraction_summary(
  analysis_datasets = analysis_datasets,
  fit_path = resolved_stage7_fit_path,
  summary_path = resolved_stage7_summary_path
)

stage8_summary_df <- compute_stage8_cure_fraction_summary(
  analysis_datasets = analysis_datasets,
  registry_path = resolved_stage8_registry_path,
  summary_path = resolved_stage8_summary_path,
  snu_label = snu_site_label
)

integrated_summary_df <- bind_rows(stage7_summary_df, stage8_summary_df) %>%
  mutate(
    dataset_variant = factor(dataset_variant, levels = c("PNU", "SNU", "merged", "merged_site_adjusted")),
    dataset_group = factor(dataset_group, levels = c("PNU", "SNU", "Merged"))
  ) %>%
  arrange(dataset_group, dataset_variant, model_family) %>%
  mutate(
    dataset_variant = as.character(dataset_variant),
    dataset_group = as.character(dataset_group)
  )

cleanup_existing_fraction_outputs(export_path)

readr::write_csv(metadata_df, file.path(export_path, "input_file_registry.csv"))
readr::write_csv(
  format_summary_for_export(stage7_summary_df),
  file.path(export_path, "frequentist_mixture_cure_fraction_summary.csv")
)
readr::write_csv(
  format_summary_for_export(stage8_summary_df),
  file.path(export_path, "bayesian_mixture_cure_fraction_summary.csv")
)
readr::write_csv(
  format_summary_for_export(integrated_summary_df),
  file.path(export_path, "integrated_mixture_cure_fraction_summary.csv")
)

comparison_specs <- tibble(dataset_group = c("PNU", "SNU", "Merged"))
output_registry <- comparison_specs %>%
  rowwise() %>%
  mutate(
    output_file = write_comparison_table(
      table_df = build_comparison_table(integrated_summary_df, dataset_group),
      export_dir = export_path,
      dataset_group = dataset_group
    )
  ) %>%
  ungroup()

readr::write_csv(output_registry, file.path(export_path, "output_table_registry.csv"))
write_fraction_tables_readme(export_path)

message("Mixture cure fraction tables written to: ", export_path)
