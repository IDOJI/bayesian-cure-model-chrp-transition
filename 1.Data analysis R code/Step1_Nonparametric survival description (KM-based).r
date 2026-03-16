# 🔴 Configure: 로컬 경로와 분석 스위치 ===============================

DATA_PATH <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
EXPORT_PATH <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦1.Step1_KM/attachments'

FILE_PREFIX <- "step1_plainkm"
YEAR_GRID <- 1:10
DATASET_BRANCHES <- c("merged", "PNU", "SNU")
SUBGROUP_VARS <- c("sex_fact")
KM_CONF_TYPE <- "log-log"
MAX_UNIQUE_LEVELS_FOR_GROUPED_KM <- 10L
SAVE_PNG <- TRUE
INSTALL_MISSING_PACKAGES <- FALSE
SEED <- 20260315

# 🔴 Prepare: 패키지와 전역 옵션 초기화 ===============================

## 🟠 Load: 필요한 패키지 ===============================

load_required_packages <- function(pkgs, install_missing = FALSE) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0L) {
    if (!install_missing) {
      stop(
        sprintf(
          "다음 패키지가 설치되어 있지 않습니다: %s\nINSTALL_MISSING_PACKAGES <- TRUE 로 바꾸거나 수동 설치 후 다시 실행하세요.",
          paste(missing_pkgs, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
  }
  
  invisible(lapply(pkgs, library, character.only = TRUE))
}

required_pkgs <- c(
  "readr",
  "dplyr",
  "tidyr",
  "tibble",
  "survival",
  "ggplot2"
)

load_required_packages(required_pkgs, install_missing = INSTALL_MISSING_PACKAGES)

options(stringsAsFactors = FALSE)
options(dplyr.summarise.inform = FALSE)
options(scipen = 999)

set.seed(SEED)

if (!dir.exists(EXPORT_PATH)) {
  dir.create(EXPORT_PATH, recursive = TRUE, showWarnings = FALSE)
}

EXPORT_PATH <- normalizePath(EXPORT_PATH, winslash = "/", mustWork = FALSE)

# 🔴 Define: Step1용 보조 함수들 ===============================

## 🟠 Create: 경로와 파일명 도우미 ===============================

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nchar(x) == 0, "unnamed", x)
}

make_output_path <- function(stem, ext) {
  file.path(EXPORT_PATH, paste0(FILE_PREFIX, "_", stem, ".", ext))
}

analysis_csv_path <- make_output_path("analysis_dataset", "csv")
qc_overview_csv_path <- make_output_path("qc_dataset_overview", "csv")
qc_status_csv_path <- make_output_path("qc_status_distribution", "csv")
qc_missing_csv_path <- make_output_path("qc_missingness", "csv")
km_yearly_csv_path <- make_output_path("km_yearly_summary", "csv")
km_curve_csv_path <- make_output_path("km_curve_data", "csv")
km_registry_csv_path <- make_output_path("km_registry", "csv")
km_bundle_rds_path <- make_output_path("km_bundle", "rds")
overall_png_path <- make_output_path("km_overall_curves", "png")
subgroup_png_path <- make_output_path("km_subgroup_curves", "png")
manifest_csv_path <- make_output_path("output_manifest", "csv")

## 🟠 Validate: 입력 데이터와 subgroup 설정 ===============================

assert_required_columns <- function(df, required_cols) {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("필수 컬럼이 없습니다: %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
}

resolve_dataset_branches <- function(df, requested_branches) {
  available_sites <- sort(unique(as.character(df$site)))
  requested_branches <- unique(requested_branches)
  
  keep_merged <- any(tolower(requested_branches) == "merged")
  requested_sites <- requested_branches[tolower(requested_branches) != "merged"]
  
  missing_sites <- setdiff(requested_sites, available_sites)
  if (length(missing_sites) > 0L) {
    warning(
      sprintf("다음 DATASET_BRANCHES 값은 site에서 찾지 못해 제외합니다: %s", paste(missing_sites, collapse = ", ")),
      call. = FALSE
    )
  }
  
  resolved <- c()
  if (keep_merged) {
    resolved <- c(resolved, "merged")
  }
  resolved <- c(resolved, intersect(requested_sites, available_sites))
  resolved <- unique(resolved)
  
  if (length(resolved) == 0L) {
    stop("유효한 DATASET_BRANCHES가 없습니다. DATASET_BRANCHES 설정을 확인하세요.", call. = FALSE)
  }
  
  resolved
}

resolve_subgroup_vars <- function(df, subgroup_vars, max_unique_levels = 10L) {
  resolved <- character(0)
  
  if (length(subgroup_vars) == 0L) {
    return(resolved)
  }
  
  for (var in subgroup_vars) {
    if (!var %in% names(df)) {
      warning(sprintf("SUBGROUP_VARS의 '%s' 컬럼이 없어 제외합니다.", var), call. = FALSE)
      next
    }
    
    non_missing <- df[[var]][!is.na(df[[var]])]
    
    if (length(non_missing) == 0L) {
      warning(sprintf("SUBGROUP_VARS의 '%s'는 전부 NA라 제외합니다.", var), call. = FALSE)
      next
    }
    
    n_unique <- dplyr::n_distinct(non_missing)
    
    if (is.numeric(df[[var]]) && n_unique > max_unique_levels) {
      warning(
        sprintf(
          "SUBGROUP_VARS의 '%s'는 numeric이고 고유값이 %d개라 plain grouped KM에 부적절하여 제외합니다.",
          var, n_unique
        ),
        call. = FALSE
      )
      next
    }
    
    resolved <- c(resolved, var)
  }
  
  unique(resolved)
}

## 🟠 Transform: Step1 분석용 survival 데이터셋 ===============================

build_analysis_dataset <- function(raw_df) {
  required_cols <- c("id", "site", "days_followup", "status_num")
  assert_required_columns(raw_df, required_cols)
  
  dat <- raw_df %>%
    dplyr::mutate(
      id = as.character(.data$id),
      site = as.character(.data$site),
      subject_key = paste(.data$site, .data$id, sep = "::"),
      days_followup = as.numeric(.data$days_followup),
      status_num = as.integer(.data$status_num),
      status_transition_only = dplyr::case_when(
        .data$status_num == 1L ~ "transition",
        .data$status_num == 2L ~ "remission_as_censor",
        .data$status_num == 0L ~ "right_censoring",
        TRUE ~ NA_character_
      ),
      event_transition = dplyr::if_else(.data$status_num == 1L, 1L, 0L, missing = NA_integer_),
      time_years = .data$days_followup / 365.25
    )
  
  if (anyNA(dat$subject_key)) {
    stop("site + id 기반 subject_key에 NA가 있습니다. id/site 컬럼을 확인하세요.", call. = FALSE)
  }
  
  dup_keys <- dat %>%
    dplyr::count(.data$subject_key, name = "n_dup") %>%
    dplyr::filter(.data$n_dup > 1L)
  
  if (nrow(dup_keys) > 0L) {
    stop(
      sprintf(
        "site + id 기준 중복 키가 발견되었습니다. 예시: %s",
        paste(utils::head(dup_keys$subject_key, 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  if (anyNA(dat$days_followup)) {
    stop("days_followup에 NA가 있습니다. Step1 전에 결측을 처리하세요.", call. = FALSE)
  }
  
  if (any(dat$days_followup < 0, na.rm = TRUE)) {
    stop("days_followup에 음수가 있습니다. 데이터 QC를 확인하세요.", call. = FALSE)
  }
  
  if (anyNA(dat$status_num)) {
    stop("status_num에 NA가 있습니다. Step1 전에 결측을 처리하세요.", call. = FALSE)
  }
  
  allowed_status <- c(0L, 1L, 2L)
  if (any(!dat$status_num %in% allowed_status, na.rm = TRUE)) {
    bad_vals <- sort(unique(dat$status_num[!dat$status_num %in% allowed_status]))
    stop(
      sprintf("status_num에 허용되지 않은 값이 있습니다: %s", paste(bad_vals, collapse = ", ")),
      call. = FALSE
    )
  }
  
  dat
}

subset_dataset_branch <- function(df, dataset_label) {
  if (tolower(dataset_label) == "merged") {
    return(df)
  }
  df %>% dplyr::filter(.data$site == dataset_label)
}

prepare_grouped_data <- function(df, group_var) {
  if (identical(group_var, "overall")) {
    return(df)
  }
  
  df %>%
    dplyr::filter(!is.na(.data[[group_var]])) %>%
    dplyr::mutate(`.__group__` = droplevels(as.factor(.data[[group_var]])))
}

## 🟠 Fit: plain KM 적합 ===============================

fit_plain_km <- function(df, group_var = "overall", conf_type = "log-log") {
  if (identical(group_var, "overall")) {
    return(
      survival::survfit(
        survival::Surv(time_years, event_transition) ~ 1,
        data = df,
        conf.type = conf_type
      )
    )
  }
  
  df_fit <- prepare_grouped_data(df, group_var)
  
  survival::survfit(
    survival::Surv(time_years, event_transition) ~ `.__group__`,
    data = df_fit,
    conf.type = conf_type
  )
}

## 🟠 Extract: yearly summary와 curve table ===============================

vec_or_default <- function(x, n, default = NA_real_) {
  if (is.null(x) || length(x) == 0L) {
    return(rep(default, n))
  }
  x
}

parse_strata_to_level <- function(strata_raw, group_var) {
  if (identical(group_var, "overall")) {
    return(rep("all", length(strata_raw)))
  }
  sub("^[^=]+=", "", strata_raw)
}

build_support_lookup <- function(df_used, group_var) {
  if (identical(group_var, "overall")) {
    return(
      tibble::tibble(
        subgroup_level = "all",
        n_subjects = nrow(df_used),
        n_transition = sum(df_used$event_transition == 1L, na.rm = TRUE),
        n_right_censor = sum(df_used$status_num == 0L, na.rm = TRUE),
        n_remission = sum(df_used$status_num == 2L, na.rm = TRUE),
        n_any_censor = sum(df_used$event_transition == 0L, na.rm = TRUE),
        max_followup_years = max(df_used$time_years, na.rm = TRUE),
        max_event_years = if (any(df_used$event_transition == 1L, na.rm = TRUE)) {
          max(df_used$time_years[df_used$event_transition == 1L], na.rm = TRUE)
        } else {
          NA_real_
        }
      )
    )
  }
  
  df_used %>%
    dplyr::mutate(subgroup_level = as.character(`.__group__`)) %>%
    dplyr::group_by(.data$subgroup_level) %>%
    dplyr::summarise(
      n_subjects = dplyr::n(),
      n_transition = sum(.data$event_transition == 1L, na.rm = TRUE),
      n_right_censor = sum(.data$status_num == 0L, na.rm = TRUE),
      n_remission = sum(.data$status_num == 2L, na.rm = TRUE),
      n_any_censor = sum(.data$event_transition == 0L, na.rm = TRUE),
      max_followup_years = max(.data$time_years, na.rm = TRUE),
      max_event_years = if (any(.data$event_transition == 1L, na.rm = TRUE)) {
        max(.data$time_years[.data$event_transition == 1L], na.rm = TRUE)
      } else {
        NA_real_
      }
    )
}

extract_km_timepoint_summary <- function(fit, df_used, dataset_label, group_var, time_grid, fit_key) {
  support_lookup <- build_support_lookup(df_used, group_var)
  s <- summary(fit, times = time_grid, extend = TRUE)
  
  n_rows <- length(s$time)
  
  out <- tibble::tibble(
    strata_raw = if (!is.null(s$strata)) as.character(s$strata) else rep("overall=all", n_rows),
    time_years = s$time,
    survival_prob = s$surv,
    lcl = vec_or_default(s$lower, n_rows, default = NA_real_),
    ucl = vec_or_default(s$upper, n_rows, default = NA_real_),
    n_risk = vec_or_default(s$n.risk, n_rows, default = NA_real_),
    n_event = vec_or_default(s$n.event, n_rows, default = NA_real_),
    n_censor = vec_or_default(s$n.censor, n_rows, default = NA_real_)
  ) %>%
    dplyr::mutate(
      dataset = dataset_label,
      data_cut = "full",
      group_var = group_var,
      subgroup_level = parse_strata_to_level(.data$strata_raw, group_var),
      fit_key = fit_key,
      model_class = "plain_KM",
      source_type = "cohort_curve"
    ) %>%
    dplyr::left_join(support_lookup, by = "subgroup_level") %>%
    dplyr::mutate(
      risk_prob = 1 - .data$survival_prob,
      risk_lcl = ifelse(is.na(.data$ucl), NA_real_, pmax(0, 1 - .data$ucl)),
      risk_ucl = ifelse(is.na(.data$lcl), NA_real_, pmin(1, 1 - .data$lcl)),
      is_extended_past_last_followup = .data$time_years > .data$max_followup_years,
      is_extended_past_last_event = ifelse(
        is.na(.data$max_event_years),
        .data$time_years > 0,
        .data$time_years > .data$max_event_years
      ),
      stage1_support_flag = dplyr::case_when(
        .data$is_extended_past_last_followup ~ "beyond_last_followup",
        .data$is_extended_past_last_event ~ "beyond_last_event",
        TRUE ~ "within_observed_event_range"
      ),
      panel_label = if (identical(group_var, "overall")) {
        dataset_label
      } else {
        paste(dataset_label, group_var, subgroup_level, sep = " | ")
      }
    ) %>%
    dplyr::arrange(.data$subgroup_level, .data$time_years)
  
  out
}

extract_km_curve_data <- function(fit, df_used, dataset_label, group_var, fit_key) {
  support_lookup <- build_support_lookup(df_used, group_var)
  s <- summary(fit, censored = TRUE)
  
  n_rows <- length(s$time)
  
  curve_df <- tibble::tibble(
    strata_raw = if (!is.null(s$strata)) as.character(s$strata) else rep("overall=all", n_rows),
    time_years = vec_or_default(s$time, n_rows, default = numeric(0)),
    survival_prob = vec_or_default(s$surv, n_rows, default = numeric(0)),
    lcl = vec_or_default(s$lower, n_rows, default = numeric(0)),
    ucl = vec_or_default(s$upper, n_rows, default = numeric(0)),
    n_risk = vec_or_default(s$n.risk, n_rows, default = numeric(0)),
    n_event = vec_or_default(s$n.event, n_rows, default = numeric(0)),
    n_censor = vec_or_default(s$n.censor, n_rows, default = numeric(0))
  )
  
  if (nrow(curve_df) > 0L) {
    curve_df <- curve_df %>%
      dplyr::mutate(
        subgroup_level = parse_strata_to_level(.data$strata_raw, group_var),
        row_type = "observed"
      )
  } else {
    curve_df <- tibble::tibble(
      strata_raw = character(0),
      time_years = numeric(0),
      survival_prob = numeric(0),
      lcl = numeric(0),
      ucl = numeric(0),
      n_risk = numeric(0),
      n_event = numeric(0),
      n_censor = numeric(0),
      subgroup_level = character(0),
      row_type = character(0)
    )
  }
  
  initial_df <- support_lookup %>%
    dplyr::mutate(
      strata_raw = if (identical(group_var, "overall")) {
        "overall=all"
      } else {
        paste0(group_var, "=", .data$subgroup_level)
      },
      time_years = 0,
      survival_prob = 1,
      lcl = 1,
      ucl = 1,
      n_risk = .data$n_subjects,
      n_event = 0,
      n_censor = 0,
      row_type = "initial"
    ) %>%
    dplyr::select(
      .data$strata_raw, .data$time_years, .data$survival_prob, .data$lcl, .data$ucl,
      .data$n_risk, .data$n_event, .data$n_censor, .data$subgroup_level, .data$row_type
    )
  
  out <- dplyr::bind_rows(initial_df, curve_df) %>%
    dplyr::mutate(
      dataset = dataset_label,
      data_cut = "full",
      group_var = group_var,
      fit_key = fit_key,
      model_class = "plain_KM",
      source_type = "cohort_curve",
      risk_prob = 1 - .data$survival_prob,
      panel_label = if (identical(group_var, "overall")) {
        dataset_label
      } else {
        paste(dataset_label, group_var, subgroup_level, sep = " | ")
      },
      row_type = factor(.data$row_type, levels = c("initial", "observed"))
    ) %>%
    dplyr::arrange(.data$subgroup_level, .data$time_years, .data$row_type)
  
  out
}

build_registry_row <- function(df_dataset, df_fit, dataset_label, group_var, fit_key, bundle_slot, subgroup_levels) {
  tibble::tibble(
    fit_key = fit_key,
    dataset = dataset_label,
    data_cut = "full",
    group_var = group_var,
    bundle_slot = bundle_slot,
    bundle_name = fit_key,
    bundle_rds_path = km_bundle_rds_path,
    n_subjects_dataset = nrow(df_dataset),
    n_subjects_fit = nrow(df_fit),
    n_missing_group_var = if (identical(group_var, "overall")) 0L else sum(is.na(df_dataset[[group_var]])),
    n_transition_fit = sum(df_fit$event_transition == 1L, na.rm = TRUE),
    n_right_censor_fit = sum(df_fit$status_num == 0L, na.rm = TRUE),
    n_remission_fit = sum(df_fit$status_num == 2L, na.rm = TRUE),
    n_any_censor_fit = sum(df_fit$event_transition == 0L, na.rm = TRUE),
    n_strata = if (identical(group_var, "overall")) 1L else dplyr::n_distinct(df_fit$.__group__),
    subgroup_levels = subgroup_levels,
    km_conf_type = KM_CONF_TYPE,
    year_grid = paste(YEAR_GRID, collapse = "|"),
    created_in_step = "Step1"
  )
}

## 🟠 Plot: curve table 기반 KM 그림 ===============================

make_overall_km_plot <- function(curve_df) {
  ggplot2::ggplot(
    curve_df,
    ggplot2::aes(
      x = .data$time_years,
      y = .data$survival_prob,
      color = .data$dataset,
      fill = .data$dataset,
      group = interaction(.data$dataset, .data$subgroup_level)
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lcl, ymax = .data$ucl),
      alpha = 0.18,
      linewidth = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::facet_wrap(~dataset, scales = "free_y") +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Step1 Plain KM: 전체 cohort별 Kaplan-Meier 곡선",
      x = "추적시간 (년)",
      y = "추정 생존확률"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey95")
    )
}

make_subgroup_km_plot <- function(curve_df) {
  ggplot2::ggplot(
    curve_df,
    ggplot2::aes(
      x = .data$time_years,
      y = .data$survival_prob,
      color = .data$subgroup_level,
      fill = .data$subgroup_level,
      group = interaction(.data$panel_label, .data$subgroup_level)
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lcl, ymax = .data$ucl),
      alpha = 0.16,
      linewidth = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::facet_wrap(~panel_label, scales = "free_y") +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Step1 Plain KM: subgroup별 Kaplan-Meier 곡선",
      x = "추적시간 (년)",
      y = "추정 생존확률",
      color = "Subgroup"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey95")
    )
}

# 🔴 Read: 원본 CSV를 불러와 Step1 분석 데이터셋 생성 ===============================

## 🟠 Import: merged dataset 로딩 ===============================

if (!file.exists(DATA_PATH)) {
  stop(sprintf("DATA_PATH에 파일이 없습니다: %s", DATA_PATH), call. = FALSE)
}

raw_df <- readr::read_csv(
  file = DATA_PATH,
  show_col_types = FALSE,
  progress = FALSE,
  locale = readr::locale(encoding = "UTF-8")
)

## 🟠 Check: 필수 컬럼과 키 무결성 ===============================

dat_analysis <- build_analysis_dataset(raw_df)
dataset_branches_resolved <- resolve_dataset_branches(dat_analysis, DATASET_BRANCHES)
subgroup_vars_resolved <- resolve_subgroup_vars(
  dat_analysis,
  SUBGROUP_VARS,
  max_unique_levels = MAX_UNIQUE_LEVELS_FOR_GROUPED_KM
)

## 🟠 Export: analysis-ready 데이터와 QC 테이블 ===============================

readr::write_csv(dat_analysis, analysis_csv_path, na = "")

branch_stack <- dplyr::bind_rows(
  lapply(dataset_branches_resolved, function(branch) {
    subset_dataset_branch(dat_analysis, branch) %>%
      dplyr::mutate(dataset = branch)
  })
)

qc_dataset_overview <- branch_stack %>%
  dplyr::group_by(.data$dataset) %>%
  dplyr::summarise(
    n_subjects = dplyr::n(),
    n_unique_subject_key = dplyr::n_distinct(.data$subject_key),
    n_transition = sum(.data$event_transition == 1L, na.rm = TRUE),
    n_right_censor = sum(.data$status_num == 0L, na.rm = TRUE),
    n_remission = sum(.data$status_num == 2L, na.rm = TRUE),
    max_followup_years = max(.data$time_years, na.rm = TRUE),
    median_followup_years = stats::median(.data$time_years, na.rm = TRUE)
  ) %>%
  dplyr::arrange(match(.data$dataset, dataset_branches_resolved))

qc_status_distribution <- branch_stack %>%
  dplyr::mutate(
    status_label = dplyr::case_when(
      .data$status_num == 0L ~ "right_censoring",
      .data$status_num == 1L ~ "transition",
      .data$status_num == 2L ~ "remission",
      TRUE ~ "unknown"
    )
  ) %>%
  dplyr::count(.data$dataset, .data$status_num, .data$status_label, name = "n") %>%
  dplyr::arrange(.data$dataset, .data$status_num)

qc_missingness <- tibble::tibble(
  variable = names(dat_analysis),
  n_missing = vapply(dat_analysis, function(x) sum(is.na(x)), numeric(1))
) %>%
  dplyr::arrange(dplyr::desc(.data$n_missing), .data$variable)

readr::write_csv(qc_dataset_overview, qc_overview_csv_path, na = "")
readr::write_csv(qc_status_distribution, qc_status_csv_path, na = "")
readr::write_csv(qc_missingness, qc_missing_csv_path, na = "")

# 🔴 Run: Step1 plain KM 적합과 요약 추출 ===============================

## 🟠 Iterate: merged와 site branch별 KM 적합 ===============================

km_overall <- list()
km_by_group <- list()
km_registry_rows <- list()
km_yearly_rows <- list()
km_curve_rows <- list()

for (dataset_label in dataset_branches_resolved) {
  dat_dataset <- subset_dataset_branch(dat_analysis, dataset_label)
  
  if (nrow(dat_dataset) == 0L) {
    next
  }
  
  overall_key <- paste(sanitize_filename(dataset_label), "overall", sep = "__")
  overall_fit <- fit_plain_km(
    df = dat_dataset,
    group_var = "overall",
    conf_type = KM_CONF_TYPE
  )
  
  km_overall[[overall_key]] <- overall_fit
  
  km_registry_rows[[overall_key]] <- build_registry_row(
    df_dataset = dat_dataset,
    df_fit = dat_dataset,
    dataset_label = dataset_label,
    group_var = "overall",
    fit_key = overall_key,
    bundle_slot = "km_overall",
    subgroup_levels = "all"
  )
  
  km_yearly_rows[[overall_key]] <- extract_km_timepoint_summary(
    fit = overall_fit,
    df_used = dat_dataset,
    dataset_label = dataset_label,
    group_var = "overall",
    time_grid = YEAR_GRID,
    fit_key = overall_key
  )
  
  km_curve_rows[[overall_key]] <- extract_km_curve_data(
    fit = overall_fit,
    df_used = dat_dataset,
    dataset_label = dataset_label,
    group_var = "overall",
    fit_key = overall_key
  )
  
  if (length(subgroup_vars_resolved) > 0L) {
    for (group_var in subgroup_vars_resolved) {
      dat_group <- prepare_grouped_data(dat_dataset, group_var)
      
      if (nrow(dat_group) == 0L) {
        next
      }
      
      subgroup_key <- paste(sanitize_filename(dataset_label), sanitize_filename(group_var), sep = "__")
      subgroup_fit <- fit_plain_km(
        df = dat_dataset,
        group_var = group_var,
        conf_type = KM_CONF_TYPE
      )
      
      km_by_group[[subgroup_key]] <- subgroup_fit
      
      subgroup_levels <- paste(sort(unique(as.character(dat_group$.__group__))), collapse = "|")
      
      km_registry_rows[[subgroup_key]] <- build_registry_row(
        df_dataset = dat_dataset,
        df_fit = dat_group,
        dataset_label = dataset_label,
        group_var = group_var,
        fit_key = subgroup_key,
        bundle_slot = "km_by_group",
        subgroup_levels = subgroup_levels
      )
      
      km_yearly_rows[[subgroup_key]] <- extract_km_timepoint_summary(
        fit = subgroup_fit,
        df_used = dat_group,
        dataset_label = dataset_label,
        group_var = group_var,
        time_grid = YEAR_GRID,
        fit_key = subgroup_key
      )
      
      km_curve_rows[[subgroup_key]] <- extract_km_curve_data(
        fit = subgroup_fit,
        df_used = dat_group,
        dataset_label = dataset_label,
        group_var = group_var,
        fit_key = subgroup_key
      )
    }
  }
}

km_registry <- dplyr::bind_rows(km_registry_rows) %>%
  dplyr::mutate(
    analysis_csv_path = analysis_csv_path,
    yearly_summary_csv_path = km_yearly_csv_path,
    curve_data_csv_path = km_curve_csv_path,
    registry_csv_path = km_registry_csv_path
  ) %>%
  dplyr::arrange(.data$dataset, .data$group_var)

km_yearly_summary <- dplyr::bind_rows(km_yearly_rows) %>%
  dplyr::arrange(.data$dataset, .data$group_var, .data$subgroup_level, .data$time_years)

km_curve_data <- dplyr::bind_rows(km_curve_rows) %>%
  dplyr::arrange(.data$dataset, .data$group_var, .data$subgroup_level, .data$time_years, .data$row_type)

# 🔴 Save: Step1 bundle과 export 파일 생성 ===============================

## 🟠 Assemble: later load용 bundle 오브젝트 ===============================

km_bundle <- list(
  meta = list(
    step = "Step1",
    purpose = "plain Kaplan-Meier fitting and storage for later model comparison",
    created_at = as.character(Sys.time()),
    seed = SEED,
    data_path = DATA_PATH,
    export_path = EXPORT_PATH,
    file_prefix = FILE_PREFIX,
    year_grid = YEAR_GRID,
    dataset_branches_requested = DATASET_BRANCHES,
    dataset_branches_resolved = dataset_branches_resolved,
    subgroup_vars_requested = SUBGROUP_VARS,
    subgroup_vars_resolved = subgroup_vars_resolved,
    km_conf_type = KM_CONF_TYPE,
    endpoint_definition = "transition only",
    censoring_definition = "right_censoring + remission treated as censoring",
    time_unit = "years"
  ),
  analysis_data = dat_analysis,
  km_overall = km_overall,
  km_by_group = km_by_group,
  km_registry = km_registry,
  km_yearly_summary = km_yearly_summary,
  km_curve_data = km_curve_data,
  qc_dataset_overview = qc_dataset_overview,
  qc_status_distribution = qc_status_distribution,
  qc_missingness = qc_missingness,
  session_info = utils::capture.output(sessionInfo())
)

## 🟠 Export: CSV, RDS, PNG, manifest 파일 ===============================

readr::write_csv(km_registry, km_registry_csv_path, na = "")
readr::write_csv(km_yearly_summary, km_yearly_csv_path, na = "")
readr::write_csv(km_curve_data, km_curve_csv_path, na = "")
saveRDS(km_bundle, km_bundle_rds_path, version = 3)

if (SAVE_PNG) {
  overall_curve_df <- km_curve_data %>%
    dplyr::filter(.data$group_var == "overall")
  
  if (nrow(overall_curve_df) > 0L) {
    p_overall <- make_overall_km_plot(overall_curve_df)
    ggplot2::ggsave(
      filename = overall_png_path,
      plot = p_overall,
      width = 11,
      height = 7,
      units = "in",
      dpi = 300
    )
  }
  
  subgroup_curve_df <- km_curve_data %>%
    dplyr::filter(.data$group_var != "overall")
  
  if (nrow(subgroup_curve_df) > 0L) {
    p_subgroup <- make_subgroup_km_plot(subgroup_curve_df)
    ggplot2::ggsave(
      filename = subgroup_png_path,
      plot = p_subgroup,
      width = 12,
      height = 8,
      units = "in",
      dpi = 300
    )
  }
}

output_manifest <- tibble::tibble(
  file_path = c(
    analysis_csv_path,
    qc_overview_csv_path,
    qc_status_csv_path,
    qc_missing_csv_path,
    km_yearly_csv_path,
    km_curve_csv_path,
    km_registry_csv_path,
    km_bundle_rds_path,
    if (file.exists(overall_png_path)) overall_png_path else NULL,
    if (file.exists(subgroup_png_path)) subgroup_png_path else NULL
  ),
  file_type = c(
    "csv",
    "csv",
    "csv",
    "csv",
    "csv",
    "csv",
    "csv",
    "rds",
    if (file.exists(overall_png_path)) "png" else NULL,
    if (file.exists(subgroup_png_path)) "png" else NULL
  ),
  description = c(
    "Step1 analysis-ready dataset with transition-only endpoint coding",
    "QC overview by dataset branch",
    "QC status distribution by dataset branch",
    "QC missingness table",
    "Plain KM yearly summary at 1-10 years",
    "Plain KM curve data used for PNG generation",
    "Plain KM registry for later bundle loading",
    "Plain KM bundle RDS with survfit objects and exported tables",
    if (file.exists(overall_png_path)) "Overall cohort KM curves generated from km_curve_data" else NULL,
    if (file.exists(subgroup_png_path)) "Subgroup KM curves generated from km_curve_data" else NULL
  )
) %>%
  dplyr::mutate(file_path = normalizePath(.data$file_path, winslash = "/", mustWork = FALSE))

readr::write_csv(output_manifest, manifest_csv_path, na = "")