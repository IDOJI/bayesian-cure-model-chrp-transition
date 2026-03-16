# 🔴 Configure: user-editable-paths-and-controls ===============================
DATA_PATH <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
EXPORT_PATH <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦6.Step6_formal cure fraction test/attachments'

STEP4_BUNDLE_PATH <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦4.Step4_Non-cure benchmark models/attachments/step4_noncure_bundle.rds"
STEP5_BUNDLE_ALL_PATH <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦5.Step5_MLE MCM/attachments/step5_fit_bundle_all.rds"
STEP5_BUNDLE_BEST_PATH <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦5.Step5_MLE MCM/attachments/step5_fit_bundle_best.rds"

INSTALL_MISSING_PACKAGES <- FALSE
RANDOM_SEED <- 20260301L

RUN_DATASETS <- c("merged", "PNU", "SNU")
SUBGROUP_VAR <- NULL
MAX_SUBGROUP_LEVELS <- 10L

PREDICTION_HORIZONS_YEARS <- 1:10
DAYS_PER_YEAR <- 365.25
KM_CONF_LEVEL <- 0.95
KM_CONF_TYPE <- "log-log"
WRITE_PNG <- TRUE

# 🔴 Bootstrap: package-runtime-and-session-options ===============================
## 🟠 Ensure: required-packages-available ===============================
required_pkgs <- c(
  "survival", "dplyr", "tidyr", "purrr", "readr",
  "tibble", "stringr", "ggplot2"
)

ensure_packages <- function(pkgs, install_missing = FALSE) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0 && isTRUE(install_missing)) {
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
  }
  still_missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing) > 0) {
    stop(
      "다음 패키지가 설치되어 있지 않습니다: ",
      paste(still_missing, collapse = ", "),
      "\nINSTALL_MISSING_PACKAGES <- TRUE 로 바꾸거나 직접 설치 후 다시 실행하세요."
    )
  }
  invisible(lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

options(stringsAsFactors = FALSE)
set.seed(RANDOM_SEED)
if (!dir.exists(EXPORT_PATH)) dir.create(EXPORT_PATH, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(DATA_PATH)) stop("DATA_PATH에 지정한 파일이 존재하지 않습니다: ", DATA_PATH)
ensure_packages(required_pkgs, install_missing = INSTALL_MISSING_PACKAGES)

# 🔴 Declare: helper-functions-for-step1 ===============================
## 🟠 Build: generic-helpers-for-io-and-qc ===============================
safe_file_stub <- function(x) {
  x |>
    stringr::str_replace_all("[^[:alnum:]_\\-]+", "_") |>
    stringr::str_replace_all("_+", "_") |>
    stringr::str_replace_all("^_|_$", "")
}

bind_or_empty <- function(x, template_df) {
  if (length(x) == 0L) return(template_df)
  dplyr::bind_rows(x)
}

stop_if_missing_columns <- function(data, cols) {
  missing_cols <- setdiff(cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("필수 컬럼이 없습니다: ", paste(missing_cols, collapse = ", "))
  }
}

standardize_status_num <- function(df) {
  if ("status_num" %in% names(df)) {
    out <- suppressWarnings(as.integer(as.character(df$status_num)))
    out[!(out %in% c(0L, 1L, 2L))] <- NA_integer_
    return(out)
  }
  if ("status" %in% names(df)) {
    status_chr <- tolower(trimws(as.character(df$status)))
    out <- dplyr::case_when(
      status_chr == "right_censoring" ~ 0L,
      status_chr == "transition" ~ 1L,
      status_chr == "remission" ~ 2L,
      TRUE ~ NA_integer_
    )
    return(out)
  }
  stop("status_num 또는 status 컬럼이 필요합니다.")
}

derive_support_flag <- function(time_vec, max_event_time, max_followup_time) {
  dplyr::case_when(
    is.na(max_event_time) ~ "no_events",
    time_vec <= max_event_time ~ "observed_support",
    time_vec <= max_followup_time ~ "beyond_last_event_censor_only",
    TRUE ~ "beyond_followup"
  )
}

make_analysis_id <- function(dataset_label, subgroup_label) {
  subgroup_stub <- if (identical(subgroup_label, "overall")) {
    "overall"
  } else {
    paste0("subgroup__", safe_file_stub(subgroup_label))
  }
  paste(dataset_label, subgroup_stub, sep = "__")
}

extract_survfit_table_value <- function(fit_obj, field_name) {
  tbl <- fit_obj$table
  if (is.null(tbl) || length(tbl) == 0L) return(NA_real_)
  if (is.matrix(tbl)) tbl <- tbl[1, , drop = TRUE]
  if (!field_name %in% names(tbl)) return(NA_real_)
  val <- suppressWarnings(as.numeric(tbl[[field_name]]))
  if (!is.finite(val)) return(NA_real_)
  val
}

## 🟠 Build: step1-data-preparation-helpers ===============================
prepare_analysis_data <- function(raw_df) {
  stop_if_missing_columns(raw_df, c("id", "site", "days_followup"))
  if (!("status_num" %in% names(raw_df) || "status" %in% names(raw_df))) {
    stop("status_num 또는 status 컬럼이 필요합니다.")
  }
  
  dat <- raw_df |>
    dplyr::mutate(
      id = as.character(id),
      site = toupper(trimws(as.character(site))),
      site_id = paste(site, id, sep = "::"),
      days_followup = suppressWarnings(as.numeric(days_followup)),
      status_num_clean = standardize_status_num(raw_df),
      status = factor(
        dplyr::case_when(
          status_num_clean == 0L ~ "right_censoring",
          status_num_clean == 1L ~ "transition",
          status_num_clean == 2L ~ "remission",
          TRUE ~ NA_character_
        ),
        levels = c("right_censoring", "remission", "transition")
      ),
      event_transition = dplyr::case_when(
        status_num_clean == 1L ~ 1L,
        status_num_clean %in% c(0L, 2L) ~ 0L,
        TRUE ~ NA_integer_
      ),
      time_years_km = days_followup / DAYS_PER_YEAR,
      sex_fact = dplyr::case_when(
        "sex_fact" %in% names(raw_df) ~ as.character(sex_fact),
        "sex_num" %in% names(raw_df) & suppressWarnings(as.numeric(sex_num)) == 0 ~ "Female",
        "sex_num" %in% names(raw_df) & suppressWarnings(as.numeric(sex_num)) == 1 ~ "Male",
        TRUE ~ NA_character_
      ),
      sex_fact = factor(sex_fact, levels = c("Female", "Male"))
    )
  
  exclusion_df <- tibble::tibble(
    row_index = seq_len(nrow(dat)),
    site_id = dat$site_id,
    exclude_reason = NA_character_
  )
  
  exclusion_df$exclude_reason[is.na(dat$id) | trimws(dat$id) == ""] <- "missing_id"
  exclusion_df$exclude_reason[is.na(exclusion_df$exclude_reason) & (is.na(dat$site) | trimws(dat$site) == "")] <- "missing_site"
  exclusion_df$exclude_reason[is.na(exclusion_df$exclude_reason) & is.na(dat$days_followup)] <- "missing_days_followup"
  exclusion_df$exclude_reason[is.na(exclusion_df$exclude_reason) & dat$days_followup < 0] <- "negative_days_followup"
  exclusion_df$exclude_reason[is.na(exclusion_df$exclude_reason) & is.na(dat$status_num_clean)] <- "missing_or_invalid_status"
  
  dat_valid <- dat[is.na(exclusion_df$exclude_reason), , drop = FALSE]
  exclusion_only <- exclusion_df[!is.na(exclusion_df$exclude_reason), , drop = FALSE]
  
  if (nrow(dat_valid) == 0L) {
    stop("전처리 후 남은 분석 데이터가 없습니다.")
  }
  
  list(
    analysis_data = dat_valid,
    initial_exclusions = exclusion_only
  )
}

make_analysis_splits <- function(data, dataset_names, subgroup_var = NULL, max_levels = 10L) {
  split_list <- list()
  message_rows <- list()
  idx <- 1L
  
  for (ds in dataset_names) {
    ds_upper <- toupper(ds)
    dataset_label <- if (identical(ds_upper, "MERGED")) "merged" else ds_upper
    
    branch_df <- if (identical(dataset_label, "merged")) {
      data
    } else {
      dplyr::filter(data, site == dataset_label)
    }
    
    if (nrow(branch_df) == 0L) {
      message_rows[[length(message_rows) + 1L]] <- tibble::tibble(
        analysis_id = paste(dataset_label, "overall", sep = "__"),
        dataset = dataset_label,
        subgroup = "overall",
        message_type = "skip",
        message_text = "해당 dataset branch에 해당하는 행이 없어 건너뜀"
      )
      next
    }
    
    split_list[[idx]] <- list(
      dataset = dataset_label,
      subgroup = "overall",
      data = branch_df
    )
    idx <- idx + 1L
    
    if (!is.null(subgroup_var)) {
      if (!subgroup_var %in% names(branch_df)) {
        stop("SUBGROUP_VAR가 데이터에 존재하지 않습니다: ", subgroup_var)
      }
      
      subgroup_values <- unique(branch_df[[subgroup_var]])
      subgroup_values <- subgroup_values[!is.na(subgroup_values)]
      subgroup_values <- as.character(subgroup_values)
      
      if (length(subgroup_values) > max_levels) {
        stop(
          "SUBGROUP_VAR의 level 수가 너무 많습니다 (", length(subgroup_values), "). ",
          "MAX_SUBGROUP_LEVELS를 늘리거나 subgroup 변수를 바꿔주세요."
        )
      }
      
      subgroup_values <- subgroup_values[order(subgroup_values)]
      
      for (sv in subgroup_values) {
        sub_df <- branch_df |>
          dplyr::filter(as.character(.data[[subgroup_var]]) == sv)
        
        if (nrow(sub_df) == 0L) next
        
        split_list[[idx]] <- list(
          dataset = dataset_label,
          subgroup = paste0(subgroup_var, "=", sv),
          data = sub_df
        )
        idx <- idx + 1L
      }
    }
  }
  
  list(
    splits = split_list,
    messages = bind_or_empty(
      message_rows,
      tibble::tibble(
        analysis_id = character(),
        dataset = character(),
        subgroup = character(),
        message_type = character(),
        message_text = character()
      )
    )
  )
}

## 🟠 Build: step1-km-summary-helpers ===============================
km_at_times <- function(fit_obj, times_vec, max_event_time, max_followup_time) {
  sum_obj <- summary(fit_obj, times = times_vec, extend = TRUE)
  
  lower_vec <- if (is.null(sum_obj$lower)) rep(NA_real_, length(times_vec)) else as.numeric(sum_obj$lower)
  upper_vec <- if (is.null(sum_obj$upper)) rep(NA_real_, length(times_vec)) else as.numeric(sum_obj$upper)
  surv_vec <- as.numeric(sum_obj$surv)
  n_risk_vec <- if (is.null(sum_obj$n.risk)) rep(NA_real_, length(times_vec)) else as.numeric(sum_obj$n.risk)
  
  tibble::tibble(
    year = as.numeric(times_vec),
    km_surv = surv_vec,
    km_risk = 1 - surv_vec,
    km_lcl = lower_vec,
    km_ucl = upper_vec,
    n_risk = n_risk_vec,
    max_event_time = max_event_time,
    max_followup_time = max_followup_time,
    support_flag = derive_support_flag(as.numeric(times_vec), max_event_time, max_followup_time),
    is_beyond_last_event = if (is.na(max_event_time)) NA else as.numeric(times_vec) > max_event_time,
    is_beyond_followup = as.numeric(times_vec) > max_followup_time,
    prediction_source = "cohort_curve"
  )
}

build_curve_dataframe <- function(fit_obj, dataset_label, subgroup_label, analysis_id, max_event_time, max_followup_time) {
  curve_times <- sort(unique(c(0, fit_obj$time, max_followup_time)))
  curve_times <- curve_times[is.finite(curve_times)]
  
  if (length(curve_times) == 0L) {
    curve_times <- c(0)
  }
  
  sum_obj <- summary(fit_obj, times = curve_times, extend = TRUE)
  
  surv_vec <- as.numeric(sum_obj$surv)
  lower_vec <- if (is.null(sum_obj$lower)) surv_vec else as.numeric(sum_obj$lower)
  upper_vec <- if (is.null(sum_obj$upper)) surv_vec else as.numeric(sum_obj$upper)
  n_risk_vec <- if (is.null(sum_obj$n.risk)) rep(NA_real_, length(curve_times)) else as.numeric(sum_obj$n.risk)
  n_event_vec <- if (is.null(sum_obj$n.event)) rep(NA_real_, length(curve_times)) else as.numeric(sum_obj$n.event)
  n_censor_vec <- if (is.null(sum_obj$n.censor)) rep(NA_real_, length(curve_times)) else as.numeric(sum_obj$n.censor)
  
  tibble::tibble(
    analysis_id = analysis_id,
    dataset = dataset_label,
    subgroup = subgroup_label,
    curve_source = "plain_KM",
    model_class = "plain_KM",
    dist = NA_character_,
    time_years = as.numeric(curve_times),
    surv_value = surv_vec,
    risk_value = 1 - surv_vec,
    lower_value = lower_vec,
    upper_value = upper_vec,
    risk_lower_value = 1 - upper_vec,
    risk_upper_value = 1 - lower_vec,
    n_risk = n_risk_vec,
    n_event = n_event_vec,
    n_censor = n_censor_vec,
    max_event_time = max_event_time,
    max_followup_time = max_followup_time,
    support_flag = derive_support_flag(as.numeric(curve_times), max_event_time, max_followup_time)
  )
}

build_registry_row <- function(fit_obj, sub_df, dataset_label, subgroup_label, analysis_id, bundle_rds_path) {
  n_total <- nrow(sub_df)
  n_event_transition <- sum(sub_df$event_transition == 1L, na.rm = TRUE)
  n_right_censoring <- sum(sub_df$status_num_clean == 0L, na.rm = TRUE)
  n_remission <- sum(sub_df$status_num_clean == 2L, na.rm = TRUE)
  n_nontransition <- sum(sub_df$event_transition == 0L, na.rm = TRUE)
  
  max_followup_years <- max(sub_df$time_years_km, na.rm = TRUE)
  last_event_years <- if (n_event_transition > 0L) {
    max(sub_df$time_years_km[sub_df$event_transition == 1L], na.rm = TRUE)
  } else {
    NA_real_
  }
  
  plateau_length_years <- if (is.finite(last_event_years)) {
    max_followup_years - last_event_years
  } else {
    NA_real_
  }
  
  n_risk_at_last_event <- if (is.finite(last_event_years)) {
    tmp <- summary(fit_obj, times = last_event_years, extend = TRUE)
    if (length(tmp$n.risk) == 0L) NA_real_ else as.numeric(tmp$n.risk[1])
  } else {
    NA_real_
  }
  
  censored_after_last_event <- if (is.finite(last_event_years)) {
    sum(sub_df$time_years_km > last_event_years & sub_df$event_transition == 0L, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  km_tail_survival <- if (length(fit_obj$surv) == 0L) 1 else as.numeric(tail(fit_obj$surv, 1))
  km_tail_risk <- 1 - km_tail_survival
  km_median_survival_years <- extract_survfit_table_value(fit_obj, "median")
  
  tibble::tibble(
    analysis_id = analysis_id,
    dataset = dataset_label,
    subgroup = subgroup_label,
    model_class = "plain_KM",
    bundle_component = "km_fit_list",
    bundle_object_name = analysis_id,
    bundle_rds_file = bundle_rds_path,
    n = n_total,
    n_event_transition = n_event_transition,
    n_nontransition = n_nontransition,
    n_right_censoring = n_right_censoring,
    n_remission = n_remission,
    max_followup_years = max_followup_years,
    last_event_years = last_event_years,
    plateau_length_years = plateau_length_years,
    n_risk_at_last_event = n_risk_at_last_event,
    censored_after_last_event = censored_after_last_event,
    km_tail_survival = km_tail_survival,
    km_tail_risk = km_tail_risk,
    km_median_survival_years = km_median_survival_years,
    fit_status = "success",
    error_message = NA_character_,
    warning_message = NA_character_
  )
}

# 🔴 Read: source-data-and-prepare-transition-endpoint ===============================
## 🟠 Import: merged-dataset3-from-csv ===============================
raw_data <- readr::read_csv(
  file = DATA_PATH,
  show_col_types = FALSE,
  guess_max = 100000
)

## 🟠 Transform: step1-analysis-ready-dataset ===============================
prepared <- prepare_analysis_data(raw_data)
dat_analysis <- prepared$analysis_data
initial_exclusions_df <- prepared$initial_exclusions

# 🔴 Construct: dataset-and-subgroup-km-plan ===============================
## 🟠 Create: branch-specific-analysis-splits ===============================
split_info <- make_analysis_splits(
  data = dat_analysis,
  dataset_names = RUN_DATASETS,
  subgroup_var = SUBGROUP_VAR,
  max_levels = MAX_SUBGROUP_LEVELS
)

analysis_splits <- split_info$splits
fit_message_rows <- list()
if (nrow(split_info$messages) > 0L) {
  fit_message_rows[[length(fit_message_rows) + 1L]] <- split_info$messages
}

if (length(analysis_splits) == 0L) {
  stop("KM 적합을 수행할 수 있는 dataset/subgroup split이 없습니다.")
}

# 🔴 Fit: plain-km-across-branches-and-subgroups ===============================
## 🟠 Estimate: overall-km-and-yearly-summaries ===============================
bundle_rds_path <- file.path(EXPORT_PATH, "step1_km_bundle.rds")

km_fit_list <- list()
registry_rows <- list()
yearly_rows <- list()
curve_rows <- list()

for (split_obj in analysis_splits) {
  dataset_label <- split_obj$dataset
  subgroup_label <- split_obj$subgroup
  sub_df <- split_obj$data
  analysis_id <- make_analysis_id(dataset_label, subgroup_label)
  
  if (nrow(sub_df) == 0L) {
    fit_message_rows[[length(fit_message_rows) + 1L]] <- tibble::tibble(
      analysis_id = analysis_id,
      dataset = dataset_label,
      subgroup = subgroup_label,
      message_type = "skip",
      message_text = "분석 대상 행이 0개라서 건너뜀"
    )
    next
  }
  
  km_warnings <- character(0)
  
  fit_obj <- tryCatch(
    withCallingHandlers(
      survival::survfit(
        survival::Surv(time_years_km, event_transition) ~ 1,
        data = sub_df,
        conf.int = KM_CONF_LEVEL,
        conf.type = KM_CONF_TYPE
      ),
      warning = function(w) {
        km_warnings <<- c(km_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  
  if (inherits(fit_obj, "error")) {
    registry_rows[[length(registry_rows) + 1L]] <- tibble::tibble(
      analysis_id = analysis_id,
      dataset = dataset_label,
      subgroup = subgroup_label,
      model_class = "plain_KM",
      bundle_component = "km_fit_list",
      bundle_object_name = analysis_id,
      bundle_rds_file = bundle_rds_path,
      n = nrow(sub_df),
      n_event_transition = sum(sub_df$event_transition == 1L, na.rm = TRUE),
      n_nontransition = sum(sub_df$event_transition == 0L, na.rm = TRUE),
      n_right_censoring = sum(sub_df$status_num_clean == 0L, na.rm = TRUE),
      n_remission = sum(sub_df$status_num_clean == 2L, na.rm = TRUE),
      max_followup_years = max(sub_df$time_years_km, na.rm = TRUE),
      last_event_years = NA_real_,
      plateau_length_years = NA_real_,
      n_risk_at_last_event = NA_real_,
      censored_after_last_event = NA_real_,
      km_tail_survival = NA_real_,
      km_tail_risk = NA_real_,
      km_median_survival_years = NA_real_,
      fit_status = "fit_failed",
      error_message = conditionMessage(fit_obj),
      warning_message = paste(unique(km_warnings), collapse = " | ")
    )
    
    fit_message_rows[[length(fit_message_rows) + 1L]] <- tibble::tibble(
      analysis_id = analysis_id,
      dataset = dataset_label,
      subgroup = subgroup_label,
      message_type = "fit_error",
      message_text = conditionMessage(fit_obj)
    )
    next
  }
  
  km_fit_list[[analysis_id]] <- fit_obj
  
  last_event_years <- if (any(sub_df$event_transition == 1L)) {
    max(sub_df$time_years_km[sub_df$event_transition == 1L], na.rm = TRUE)
  } else {
    NA_real_
  }
  max_followup_years <- max(sub_df$time_years_km, na.rm = TRUE)
  
  yearly_df <- km_at_times(
    fit_obj = fit_obj,
    times_vec = PREDICTION_HORIZONS_YEARS,
    max_event_time = last_event_years,
    max_followup_time = max_followup_years
  ) |>
    dplyr::mutate(
      analysis_id = analysis_id,
      dataset = dataset_label,
      subgroup = subgroup_label,
      model_class = "plain_KM",
      dist = NA_character_
    ) |>
    dplyr::select(
      analysis_id, dataset, subgroup, model_class, dist, year,
      km_surv, km_risk, km_lcl, km_ucl, n_risk,
      max_event_time, max_followup_time,
      support_flag, is_beyond_last_event, is_beyond_followup,
      prediction_source
    )
  
  curve_df <- build_curve_dataframe(
    fit_obj = fit_obj,
    dataset_label = dataset_label,
    subgroup_label = subgroup_label,
    analysis_id = analysis_id,
    max_event_time = last_event_years,
    max_followup_time = max_followup_years
  )
  
  registry_df <- build_registry_row(
    fit_obj = fit_obj,
    sub_df = sub_df,
    dataset_label = dataset_label,
    subgroup_label = subgroup_label,
    analysis_id = analysis_id,
    bundle_rds_path = bundle_rds_path
  ) |>
    dplyr::mutate(
      warning_message = if (length(km_warnings) == 0L) NA_character_ else paste(unique(km_warnings), collapse = " | ")
    )
  
  yearly_rows[[length(yearly_rows) + 1L]] <- yearly_df
  curve_rows[[length(curve_rows) + 1L]] <- curve_df
  registry_rows[[length(registry_rows) + 1L]] <- registry_df
}

km_registry_df <- bind_or_empty(
  registry_rows,
  tibble::tibble(
    analysis_id = character(),
    dataset = character(),
    subgroup = character(),
    model_class = character(),
    bundle_component = character(),
    bundle_object_name = character(),
    bundle_rds_file = character(),
    n = numeric(),
    n_event_transition = numeric(),
    n_nontransition = numeric(),
    n_right_censoring = numeric(),
    n_remission = numeric(),
    max_followup_years = numeric(),
    last_event_years = numeric(),
    plateau_length_years = numeric(),
    n_risk_at_last_event = numeric(),
    censored_after_last_event = numeric(),
    km_tail_survival = numeric(),
    km_tail_risk = numeric(),
    km_median_survival_years = numeric(),
    fit_status = character(),
    error_message = character(),
    warning_message = character()
  )
) |>
  dplyr::arrange(dataset, subgroup)

km_yearly_all <- bind_or_empty(
  yearly_rows,
  tibble::tibble(
    analysis_id = character(),
    dataset = character(),
    subgroup = character(),
    model_class = character(),
    dist = character(),
    year = numeric(),
    km_surv = numeric(),
    km_risk = numeric(),
    km_lcl = numeric(),
    km_ucl = numeric(),
    n_risk = numeric(),
    max_event_time = numeric(),
    max_followup_time = numeric(),
    support_flag = character(),
    is_beyond_last_event = logical(),
    is_beyond_followup = logical(),
    prediction_source = character()
  )
) |>
  dplyr::arrange(dataset, subgroup, year)

km_curve_data <- bind_or_empty(
  curve_rows,
  tibble::tibble(
    analysis_id = character(),
    dataset = character(),
    subgroup = character(),
    curve_source = character(),
    model_class = character(),
    dist = character(),
    time_years = numeric(),
    surv_value = numeric(),
    risk_value = numeric(),
    lower_value = numeric(),
    upper_value = numeric(),
    risk_lower_value = numeric(),
    risk_upper_value = numeric(),
    n_risk = numeric(),
    n_event = numeric(),
    n_censor = numeric(),
    max_event_time = numeric(),
    max_followup_time = numeric(),
    support_flag = character()
  )
) |>
  dplyr::arrange(dataset, subgroup, time_years)

fit_messages_df <- bind_or_empty(
  fit_message_rows,
  tibble::tibble(
    analysis_id = character(),
    dataset = character(),
    subgroup = character(),
    message_type = character(),
    message_text = character()
  )
) |>
  dplyr::arrange(dataset, subgroup, message_type)

overall_keys <- km_registry_df |>
  dplyr::filter(subgroup == "overall", fit_status == "success") |>
  dplyr::pull(analysis_id)

km_overall <- km_fit_list[overall_keys]
km_by_group <- km_fit_list[setdiff(names(km_fit_list), overall_keys)]

# 🔴 Assemble: plot-ready-dataframes-from-source-of-truth ===============================
## 🟠 Draw: survival-and-risk-curves-from-exported-dataframes ===============================
if (WRITE_PNG && nrow(km_curve_data) > 0L) {
  survival_plot <- ggplot2::ggplot(
    km_curve_data,
    ggplot2::aes(x = time_years, y = surv_value)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower_value, ymax = upper_value),
      alpha = 0.15,
      na.rm = TRUE
    ) +
    ggplot2::geom_step(linewidth = 0.7, na.rm = TRUE) +
    ggplot2::facet_grid(subgroup ~ dataset, scales = "fixed") +
    ggplot2::scale_x_continuous(breaks = PREDICTION_HORIZONS_YEARS) +
    ggplot2::labs(
      x = "Years since cohort entry",
      y = "Kaplan-Meier survival probability",
      title = "Step1: plain KM survival curves"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  risk_plot <- ggplot2::ggplot(
    km_curve_data,
    ggplot2::aes(x = time_years, y = risk_value)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = risk_lower_value, ymax = risk_upper_value),
      alpha = 0.15,
      na.rm = TRUE
    ) +
    ggplot2::geom_step(linewidth = 0.7, na.rm = TRUE) +
    ggplot2::facet_grid(subgroup ~ dataset, scales = "fixed") +
    ggplot2::scale_x_continuous(breaks = PREDICTION_HORIZONS_YEARS) +
    ggplot2::labs(
      x = "Years since cohort entry",
      y = "Kaplan-Meier risk = 1 - survival",
      title = "Step1: plain KM risk curves"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  ggplot2::ggsave(
    filename = file.path(EXPORT_PATH, "step1_km_survival_curves.png"),
    plot = survival_plot,
    width = 12,
    height = max(4, 2 + 2 * max(1, dplyr::n_distinct(km_curve_data$subgroup))),
    dpi = 300
  )
  
  ggplot2::ggsave(
    filename = file.path(EXPORT_PATH, "step1_km_risk_curves.png"),
    plot = risk_plot,
    width = 12,
    height = max(4, 2 + 2 * max(1, dplyr::n_distinct(km_curve_data$subgroup))),
    dpi = 300
  )
}

# 🔴 Export: step1-artifacts-and-reusable-bundle ===============================
## 🟠 Persist: km-fits-and-source-of-truth-files ===============================
step1_bundle <- list(
  meta = list(
    step = "Step1_plain_KM",
    data_path = DATA_PATH,
    export_path = EXPORT_PATH,
    run_datasets = RUN_DATASETS,
    subgroup_var = SUBGROUP_VAR,
    max_subgroup_levels = MAX_SUBGROUP_LEVELS,
    prediction_horizons_years = PREDICTION_HORIZONS_YEARS,
    days_per_year = DAYS_PER_YEAR,
    km_conf_level = KM_CONF_LEVEL,
    km_conf_type = KM_CONF_TYPE,
    remission_handling = "remission_as_censoring_in_transition_only_endpoint",
    comparison_note = "KM-vs-model pairwise comparison should be performed after Step7 by loading this bundle.",
    reference_paths = list(
      step4_bundle = STEP4_BUNDLE_PATH,
      step5_bundle_all = STEP5_BUNDLE_ALL_PATH,
      step5_bundle_best = STEP5_BUNDLE_BEST_PATH
    ),
    random_seed = RANDOM_SEED
  ),
  analysis_data_prepared = dat_analysis,
  initial_exclusions = initial_exclusions_df,
  km_fit_list = km_fit_list,
  km_overall = km_overall,
  km_by_group = km_by_group,
  km_registry = km_registry_df,
  km_yearly_all = km_yearly_all,
  km_curve_data = km_curve_data,
  fit_messages = fit_messages_df
)

saveRDS(
  object = step1_bundle,
  file = bundle_rds_path
)

readr::write_csv(
  km_registry_df,
  file.path(EXPORT_PATH, "step1_km_fit_registry.csv"),
  na = ""
)

readr::write_csv(
  km_yearly_all,
  file.path(EXPORT_PATH, "step1_km_yearly_summary.csv"),
  na = ""
)

readr::write_csv(
  km_curve_data,
  file.path(EXPORT_PATH, "step1_km_curve_data.csv"),
  na = ""
)

readr::write_csv(
  initial_exclusions_df,
  file.path(EXPORT_PATH, "step1_km_initial_exclusions.csv"),
  na = ""
)

readr::write_csv(
  fit_messages_df,
  file.path(EXPORT_PATH, "step1_km_fit_messages.csv"),
  na = ""
)