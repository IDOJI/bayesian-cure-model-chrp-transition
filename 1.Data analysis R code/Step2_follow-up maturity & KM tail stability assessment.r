# 🔴 설정: 경로와 전역 파라미터 지정 ===============================
data_path <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
export_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦2.Step2_follow-up maturity, KM tail stability/attachments'

analysis_horizons_years <- 1:10
days_per_year <- 365.25
km_conf_type <- "log-log"
stable_gap_threshold <- 0.10
time_epsilon_days <- 1e-6

# 기본값은 overall only
# 예: subgroup_vars <- c("sex_fact")
subgroup_vars <- character(0)

# site 표기가 조금 다를 가능성을 대비한 매핑
site_label_map <- list(
  PNU = c("PNU"),
  SNU = c("SNU", "CHR-P", "CHR_P", "CHRP")
)

# exact PTFR/FPT는 반복 방문 스케줄/interval 정보가 필요하므로
# 현재 단일 follow-up 스키마에서는 CCI/SPT를 PTFR의 practical proxy로 사용
compute_ptfr_exact <- FALSE

# 🔴 준비: 패키지와 옵션 초기화 ===============================
options(stringsAsFactors = FALSE, scipen = 999)

required_pkgs <- c("survival", "dplyr", "purrr", "tibble", "readr", "ggplot2")

if (!dir.exists(export_path)) {
  dir.create(export_path, recursive = TRUE, showWarnings = FALSE)
}

# 🟠 함수: 패키지 검사와 파일 경로 유틸리티 ===============================
check_required_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "다음 패키지가 설치되어 있지 않습니다: ",
      paste(missing_pkgs, collapse = ", "),
      "\ninstall.packages(c(",
      paste(sprintf('\"%s\"', missing_pkgs), collapse = ", "),
      ")) 후 다시 실행하세요."
    )
  }
}

check_required_packages(required_pkgs)

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(ggplot2)
})

make_export_file <- function(filename) {
  file.path(export_path, filename)
}

days_to_years <- function(x) {
  as.numeric(x) / days_per_year
}

vec_or_default <- function(x, n, default = NA_real_) {
  if (is.null(x) || length(x) == 0) {
    rep(default, n)
  } else {
    x
  }
}

normalize_site <- function(site_chr) {
  site_chr <- trimws(as.character(site_chr))
  dplyr::case_when(
    site_chr %in% site_label_map$PNU ~ "PNU",
    site_chr %in% site_label_map$SNU ~ "SNU",
    TRUE ~ site_chr
  )
}

safe_median <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    NA_real_
  } else {
    stats::median(x, na.rm = TRUE)
  }
}

empty_survival_curve_df <- function() {
  tibble(
    panel_id = character(),
    panel_label = character(),
    dataset_group = character(),
    subgroup_var = character(),
    subgroup_level = character(),
    curve_family = character(),
    curve_name = character(),
    time_days = numeric(),
    time_years = numeric(),
    surv = numeric(),
    lower_ci = numeric(),
    upper_ci = numeric(),
    n_risk = integer(),
    n_event = integer(),
    n_censor = integer()
  )
}

empty_gap_curve_df <- function() {
  tibble(
    time_days = numeric(),
    time_years = numeric(),
    upper_minus_lower = numeric(),
    upper_minus_km = numeric(),
    km_minus_lower = numeric()
  )
}

# 🟠 함수: 생존객체 정리와 시점 평가 ===============================
quantile_from_survfit <- function(fit, probs = c(0.25, 0.5, 0.75)) {
  if (is.null(fit) || length(fit$time) == 0) {
    out <- rep(NA_real_, length(probs))
    names(out) <- paste0("q", formatC(probs * 100, format = "fg"))
    return(out)
  }
  
  surv_vals <- fit$surv
  time_vals <- fit$time
  
  out <- vapply(
    probs,
    FUN.VALUE = numeric(1),
    FUN = function(p) {
      target <- 1 - p
      idx <- which(surv_vals <= target)[1]
      if (length(idx) == 0 || is.na(idx)) {
        NA_real_
      } else {
        as.numeric(time_vals[idx])
      }
    }
  )
  
  names(out) <- paste0("q", formatC(probs * 100, format = "fg"))
  out
}

eval_survfit_at_times <- function(fit, times_days) {
  if (is.null(fit)) {
    return(
      tibble(
        time_days = times_days,
        surv = NA_real_,
        lower_ci = NA_real_,
        upper_ci = NA_real_,
        n_risk = NA_integer_,
        n_event = NA_integer_,
        n_censor = NA_integer_
      )
    )
  }
  
  s <- summary(fit, times = times_days, extend = TRUE)
  
  tibble(
    time_days = times_days,
    surv = vec_or_default(s$surv, length(times_days), NA_real_),
    lower_ci = vec_or_default(s$lower, length(times_days), NA_real_),
    upper_ci = vec_or_default(s$upper, length(times_days), NA_real_),
    n_risk = as.integer(vec_or_default(s$n.risk, length(times_days), NA_real_)),
    n_event = as.integer(vec_or_default(s$n.event, length(times_days), NA_real_)),
    n_censor = as.integer(vec_or_default(s$n.censor, length(times_days), NA_real_))
  )
}

tidy_survfit_curve <- function(fit, panel_info, curve_family, curve_name) {
  if (is.null(fit)) {
    return(empty_survival_curve_df())
  }
  
  n_subjects <- panel_info$n_subjects
  n_points <- length(fit$time)
  
  base_row <- tibble(
    panel_id = panel_info$panel_id,
    panel_label = panel_info$panel_label,
    dataset_group = panel_info$dataset_group,
    subgroup_var = panel_info$subgroup_var,
    subgroup_level = panel_info$subgroup_level,
    curve_family = curve_family,
    curve_name = curve_name,
    time_days = 0,
    time_years = 0,
    surv = 1,
    lower_ci = 1,
    upper_ci = 1,
    n_risk = as.integer(n_subjects),
    n_event = 0L,
    n_censor = 0L
  )
  
  if (n_points == 0) {
    return(base_row)
  }
  
  curve_rows <- tibble(
    panel_id = panel_info$panel_id,
    panel_label = panel_info$panel_label,
    dataset_group = panel_info$dataset_group,
    subgroup_var = panel_info$subgroup_var,
    subgroup_level = panel_info$subgroup_level,
    curve_family = curve_family,
    curve_name = curve_name,
    time_days = as.numeric(fit$time),
    time_years = days_to_years(fit$time),
    surv = as.numeric(fit$surv),
    lower_ci = vec_or_default(fit$lower, n_points, NA_real_),
    upper_ci = vec_or_default(fit$upper, n_points, NA_real_),
    n_risk = as.integer(vec_or_default(fit$n.risk, n_points, NA_real_)),
    n_event = as.integer(vec_or_default(fit$n.event, n_points, NA_real_)),
    n_censor = as.integer(vec_or_default(fit$n.censor, n_points, NA_real_))
  )
  
  bind_rows(base_row, curve_rows)
}

next_event_time_after_censor <- function(censor_time, event_times, epsilon = 1e-6) {
  if (length(event_times) == 0) {
    return(censor_time + epsilon)
  }
  next_time <- event_times[event_times >= censor_time][1]
  if (is.na(next_time)) censor_time + epsilon else next_time
}

# 🟠 함수: Step2 척도와 stability 요약 ===============================
compute_percentage_method <- function(time_days, event_primary, tau_days) {
  mean((time_days >= tau_days) | (event_primary == 1L & time_days <= tau_days))
}

compute_cci <- function(time_days, event_primary, tau_days) {
  observed_pt <- sum(pmin(time_days, tau_days))
  potential_pt <- sum(ifelse(event_primary == 0L & time_days < tau_days, tau_days, pmin(time_days, tau_days)))
  observed_pt / potential_pt
}

compute_spt <- function(time_days, event_primary, tau_days) {
  n <- length(time_days)
  numerator <- sum(ifelse(event_primary == 0L & time_days < tau_days, time_days, tau_days))
  numerator / (n * tau_days)
}

compute_gap_area <- function(km_fit, upper_fit, lower_fit, last_event_time_days) {
  if (is.na(last_event_time_days) || last_event_time_days <= 0) {
    return(
      list(
        instability_area = NA_real_,
        instability_area_norm = NA_real_,
        gap_curve = empty_gap_curve_df()
      )
    )
  }
  
  candidate_grid <- sort(unique(c(
    0,
    km_fit$time[km_fit$time <= last_event_time_days],
    upper_fit$time[upper_fit$time <= last_event_time_days],
    lower_fit$time[lower_fit$time <= last_event_time_days],
    last_event_time_days
  )))
  
  candidate_grid <- candidate_grid[!is.na(candidate_grid)]
  
  if (length(candidate_grid) < 2) {
    gap_df <- tibble(
      time_days = candidate_grid,
      time_years = days_to_years(candidate_grid),
      upper_minus_lower = 0,
      upper_minus_km = 0,
      km_minus_lower = 0
    )
    
    return(
      list(
        instability_area = 0,
        instability_area_norm = 0,
        gap_curve = gap_df
      )
    )
  }
  
  upper_eval <- eval_survfit_at_times(upper_fit, candidate_grid)$surv
  lower_eval <- eval_survfit_at_times(lower_fit, candidate_grid)$surv
  km_eval <- eval_survfit_at_times(km_fit, candidate_grid)$surv
  
  gap_total <- pmax(upper_eval - lower_eval, 0)
  gap_upper_km <- pmax(upper_eval - km_eval, 0)
  gap_km_lower <- pmax(km_eval - lower_eval, 0)
  
  instability_area <- sum(gap_total[-length(gap_total)] * diff(candidate_grid))
  instability_area_norm <- instability_area / last_event_time_days
  
  gap_df <- tibble(
    time_days = candidate_grid,
    time_years = days_to_years(candidate_grid),
    upper_minus_lower = gap_total,
    upper_minus_km = gap_upper_km,
    km_minus_lower = gap_km_lower
  )
  
  list(
    instability_area = instability_area,
    instability_area_norm = instability_area_norm,
    gap_curve = gap_df
  )
}

analyze_panel <- function(df_panel, panel_info) {
  if (nrow(df_panel) == 0) {
    return(NULL)
  }
  
  n_subjects <- nrow(df_panel)
  n_events <- sum(df_panel$event_primary == 1L, na.rm = TRUE)
  n_censored_total <- sum(df_panel$event_primary == 0L, na.rm = TRUE)
  n_right_censoring <- sum(df_panel$status_num == 0L, na.rm = TRUE)
  n_remission <- sum(df_panel$status_num == 2L, na.rm = TRUE)
  
  km_fit <- survival::survfit(
    survival::Surv(time_days, event_primary) ~ 1,
    data = df_panel,
    conf.type = km_conf_type
  )
  
  reverse_km_fit <- survival::survfit(
    survival::Surv(time_days, 1L - event_primary) ~ 1,
    data = df_panel,
    conf.type = km_conf_type
  )
  
  observed_time_fit <- survival::survfit(
    survival::Surv(df_panel$time_days, rep(1L, n_subjects)) ~ 1,
    conf.type = km_conf_type
  )
  
  eventfree_df <- df_panel %>% filter(event_primary == 0L)
  
  eventfree_followup_fit <- if (nrow(eventfree_df) > 0) {
    survival::survfit(
      survival::Surv(eventfree_df$time_days, rep(1L, nrow(eventfree_df))) ~ 1,
      conf.type = km_conf_type
    )
  } else {
    NULL
  }
  
  last_followup_time_days <- max(df_panel$time_days, na.rm = TRUE)
  last_censor_time_days <- if (n_censored_total > 0) max(df_panel$time_days[df_panel$event_primary == 0L], na.rm = TRUE) else NA_real_
  last_event_time_days <- if (n_events > 0) max(df_panel$time_days[df_panel$event_primary == 1L], na.rm = TRUE) else NA_real_
  plateau_length_days <- if (!is.na(last_event_time_days)) last_followup_time_days - last_event_time_days else NA_real_
  
  n_risk_at_last_event <- if (!is.na(last_event_time_days)) {
    as.integer(summary(km_fit, times = last_event_time_days, extend = TRUE)$n.risk[1])
  } else {
    NA_integer_
  }
  
  censored_after_last_event <- if (!is.na(last_event_time_days)) {
    sum(df_panel$event_primary == 0L & df_panel$time_days > last_event_time_days, na.rm = TRUE)
  } else {
    NA_integer_
  }
  
  if (n_events > 0) {
    event_times <- sort(unique(df_panel$time_days[df_panel$event_primary == 1L]))
    
    upper_time_days <- ifelse(
      df_panel$event_primary == 0L,
      pmax(df_panel$time_days, last_event_time_days + time_epsilon_days),
      df_panel$time_days
    )
    
    upper_fit <- survival::survfit(
      survival::Surv(upper_time_days, df_panel$event_primary) ~ 1,
      conf.type = km_conf_type
    )
    
    lower_time_days <- df_panel$time_days
    lower_time_days[df_panel$event_primary == 0L] <- vapply(
      df_panel$time_days[df_panel$event_primary == 0L],
      FUN.VALUE = numeric(1),
      FUN = next_event_time_after_censor,
      event_times = event_times,
      epsilon = time_epsilon_days
    )
    
    lower_fit <- survival::survfit(
      survival::Surv(lower_time_days, rep(1L, n_subjects)) ~ 1,
      conf.type = km_conf_type
    )
  } else {
    upper_fit <- km_fit
    lower_fit <- km_fit
  }
  
  km_quant <- quantile_from_survfit(km_fit, c(0.25, 0.50, 0.75))
  upper_quant <- quantile_from_survfit(upper_fit, c(0.25, 0.50, 0.75))
  lower_quant <- quantile_from_survfit(lower_fit, c(0.25, 0.50, 0.75))
  reverse_quant <- quantile_from_survfit(reverse_km_fit, c(0.25, 0.50, 0.75))
  
  gap_obj <- compute_gap_area(
    km_fit = km_fit,
    upper_fit = upper_fit,
    lower_fit = lower_fit,
    last_event_time_days = last_event_time_days
  )
  
  yearly_grid_days <- analysis_horizons_years * days_per_year
  
  km_eval <- eval_survfit_at_times(km_fit, yearly_grid_days)
  upper_eval <- eval_survfit_at_times(upper_fit, yearly_grid_days)
  lower_eval <- eval_survfit_at_times(lower_fit, yearly_grid_days)
  reverse_eval <- eval_survfit_at_times(reverse_km_fit, yearly_grid_days)
  observed_eval <- eval_survfit_at_times(observed_time_fit, yearly_grid_days)
  eventfree_eval <- eval_survfit_at_times(eventfree_followup_fit, yearly_grid_days)
  
  percentage_method_vec <- purrr::map_dbl(
    yearly_grid_days,
    ~ compute_percentage_method(df_panel$time_days, df_panel$event_primary, .x)
  )
  cci_vec <- purrr::map_dbl(
    yearly_grid_days,
    ~ compute_cci(df_panel$time_days, df_panel$event_primary, .x)
  )
  spt_vec <- purrr::map_dbl(
    yearly_grid_days,
    ~ compute_spt(df_panel$time_days, df_panel$event_primary, .x)
  )
  
  support_flags <- tibble(
    panel_id = panel_info$panel_id,
    panel_label = panel_info$panel_label,
    dataset_group = panel_info$dataset_group,
    subgroup_var = panel_info$subgroup_var,
    subgroup_level = panel_info$subgroup_level,
    horizon_year = analysis_horizons_years,
    horizon_days = yearly_grid_days,
    km_surv = km_eval$surv,
    km_risk = 1 - km_eval$surv,
    km_lower_ci = km_eval$lower_ci,
    km_upper_ci = km_eval$upper_ci,
    n_risk = km_eval$n_risk,
    n_event_at_horizon = km_eval$n_event,
    n_censor_at_horizon = km_eval$n_censor,
    km_upper_stability = upper_eval$surv,
    km_lower_stability = lower_eval$surv,
    stability_gap = pmax(upper_eval$surv - lower_eval$surv, 0),
    upper_minus_km = pmax(upper_eval$surv - km_eval$surv, 0),
    km_minus_lower = pmax(km_eval$surv - lower_eval$surv, 0),
    reverse_km_surv = reverse_eval$surv,
    observed_time_surv = observed_eval$surv,
    eventfree_followup_surv = eventfree_eval$surv,
    percentage_method = percentage_method_vec,
    cci = cci_vec,
    spt = spt_vec,
    last_event_time_days = last_event_time_days,
    last_followup_time_days = last_followup_time_days,
    plateau_length_days = plateau_length_days
  ) %>%
    mutate(
      ptfr_lower_proxy = cci,
      ptfr_upper_proxy = spt,
      ptfr_exact_estimable = compute_ptfr_exact,
      support_reason = case_when(
        is.na(last_event_time_days) ~ "no_observed_event",
        horizon_days > last_followup_time_days ~ "beyond_last_followup",
        horizon_days > last_event_time_days ~ "after_last_observed_event",
        stability_gap <= stable_gap_threshold ~ "betensky_gap_within_threshold",
        TRUE ~ "betensky_gap_exceeds_threshold"
      ),
      support_flag = case_when(
        support_reason == "betensky_gap_within_threshold" ~ "observed_stable",
        support_reason == "betensky_gap_exceeds_threshold" ~ "observed_unstable",
        TRUE ~ "beyond_reliable_tail"
      )
    )
  
  km_yearly_summary <- support_flags %>%
    transmute(
      panel_id,
      panel_label,
      dataset_group,
      subgroup_var,
      subgroup_level,
      horizon_year,
      horizon_days,
      km_surv,
      km_risk,
      km_lower_ci,
      km_upper_ci,
      n_risk,
      support_flag
    )
  
  survival_curves <- bind_rows(
    tidy_survfit_curve(km_fit, panel_info, "km_stability", "KM"),
    tidy_survfit_curve(upper_fit, panel_info, "km_stability", "Betensky upper"),
    tidy_survfit_curve(lower_fit, panel_info, "km_stability", "Betensky lower"),
    tidy_survfit_curve(reverse_km_fit, panel_info, "followup_distribution", "Reverse KM: C"),
    tidy_survfit_curve(observed_time_fit, panel_info, "followup_distribution", "Observed T"),
    tidy_survfit_curve(eventfree_followup_fit, panel_info, "followup_distribution", "C | C < X")
  )
  
  stability_gap_curves <- gap_obj$gap_curve %>%
    mutate(
      panel_id = panel_info$panel_id,
      panel_label = panel_info$panel_label,
      dataset_group = panel_info$dataset_group,
      subgroup_var = panel_info$subgroup_var,
      subgroup_level = panel_info$subgroup_level
    ) %>%
    relocate(panel_id:subgroup_level)
  
  followup_summary <- tibble(
    panel_id = panel_info$panel_id,
    panel_label = panel_info$panel_label,
    dataset_group = panel_info$dataset_group,
    subgroup_var = panel_info$subgroup_var,
    subgroup_level = panel_info$subgroup_level,
    n_subjects = n_subjects,
    n_events = n_events,
    n_censored_total = n_censored_total,
    n_right_censoring = n_right_censoring,
    n_remission = n_remission,
    last_event_time_days = last_event_time_days,
    last_event_time_years = days_to_years(last_event_time_days),
    last_censor_time_days = last_censor_time_days,
    last_censor_time_years = days_to_years(last_censor_time_days),
    last_followup_time_days = last_followup_time_days,
    last_followup_time_years = days_to_years(last_followup_time_days),
    plateau_length_days = plateau_length_days,
    plateau_length_years = days_to_years(plateau_length_days),
    n_risk_at_last_event = n_risk_at_last_event,
    censored_after_last_event = censored_after_last_event,
    reverse_km_q25_days = reverse_quant["q25"],
    reverse_km_q50_days = reverse_quant["q50"],
    reverse_km_q75_days = reverse_quant["q75"],
    reverse_km_q25_years = days_to_years(reverse_quant["q25"]),
    reverse_km_q50_years = days_to_years(reverse_quant["q50"]),
    reverse_km_q75_years = days_to_years(reverse_quant["q75"]),
    observation_median_days = safe_median(df_panel$time_days),
    observation_median_years = days_to_years(safe_median(df_panel$time_days)),
    eventfree_followup_median_days = safe_median(eventfree_df$time_days),
    eventfree_followup_median_years = days_to_years(safe_median(eventfree_df$time_days)),
    km_q25_days = km_quant["q25"],
    km_q50_days = km_quant["q50"],
    km_q75_days = km_quant["q75"],
    km_upper_q25_days = upper_quant["q25"],
    km_upper_q50_days = upper_quant["q50"],
    km_upper_q75_days = upper_quant["q75"],
    km_lower_q25_days = lower_quant["q25"],
    km_lower_q50_days = lower_quant["q50"],
    km_lower_q75_days = lower_quant["q75"],
    instability_area = gap_obj$instability_area,
    instability_area_norm = gap_obj$instability_area_norm
  )
  
  panel_registry <- tibble(
    panel_id = panel_info$panel_id,
    panel_label = panel_info$panel_label,
    dataset_group = panel_info$dataset_group,
    subgroup_var = panel_info$subgroup_var,
    subgroup_level = panel_info$subgroup_level,
    n_subjects = n_subjects,
    n_events = n_events,
    n_censored_total = n_censored_total,
    n_right_censoring = n_right_censoring,
    n_remission = n_remission
  )
  
  fits <- list(
    km_fit = km_fit,
    reverse_km_fit = reverse_km_fit,
    observed_time_fit = observed_time_fit,
    eventfree_followup_fit = eventfree_followup_fit,
    upper_fit = upper_fit,
    lower_fit = lower_fit
  )
  
  list(
    panel_registry = panel_registry,
    followup_summary = followup_summary,
    km_yearly_summary = km_yearly_summary,
    support_flags = support_flags,
    survival_curves = survival_curves,
    stability_gap_curves = stability_gap_curves,
    fits = fits
  )
}

# 🔴 입력: 원시 데이터 읽기와 컬럼 검증 ===============================
if (!file.exists(data_path)) {
  stop("data_path에 지정된 파일을 찾을 수 없습니다: ", data_path)
}

raw_df <- readr::read_csv(
  file = data_path,
  show_col_types = FALSE,
  progress = FALSE
)

if (!"status_num" %in% names(raw_df)) {
  if ("status" %in% names(raw_df)) {
    raw_df <- raw_df %>%
      mutate(
        status_num = case_when(
          as.character(status) == "right_censoring" ~ 0L,
          as.character(status) == "transition" ~ 1L,
          as.character(status) == "remission" ~ 2L,
          TRUE ~ NA_integer_
        )
      )
  } else {
    stop("status_num 또는 status 컬럼이 필요합니다.")
  }
}

required_cols <- c("id", "site", "days_followup", "status_num")
missing_cols <- setdiff(required_cols, names(raw_df))
if (length(missing_cols) > 0) {
  stop("필수 컬럼이 누락되었습니다: ", paste(missing_cols, collapse = ", "))
}

# 🔴 전처리: Step2용 분석 데이터 구성 ===============================
analysis_df <- raw_df %>%
  mutate(
    id = as.character(id),
    site_raw = trimws(as.character(site)),
    site_std = normalize_site(site_raw),
    subject_uid = paste(site_std, id, sep = "::"),
    time_days = as.numeric(as.character(days_followup)),
    time_years = days_to_years(as.numeric(as.character(days_followup))),
    status_num = as.integer(as.character(status_num)),
    event_primary = as.integer(status_num == 1L),
    status_primary = case_when(
      status_num == 1L ~ "transition",
      status_num == 2L ~ "remission_censored",
      status_num == 0L ~ "right_censored",
      TRUE ~ NA_character_
    )
  )

if (any(is.na(analysis_df$id))) {
  stop("id에 NA가 존재합니다.")
}

if (any(is.na(analysis_df$site_std))) {
  stop("site 정규화 결과에 NA가 존재합니다.")
}

if (any(is.na(analysis_df$time_days))) {
  stop("days_followup을 numeric으로 변환할 수 없는 값이 있습니다.")
}

if (any(analysis_df$time_days < 0, na.rm = TRUE)) {
  stop("days_followup/time_days에 음수가 존재합니다. 원시 데이터를 확인하세요.")
}

if (any(is.na(analysis_df$status_num))) {
  stop("status_num을 integer로 변환할 수 없는 값이 있습니다.")
}

if (!all(na.omit(unique(analysis_df$status_num)) %in% c(0, 1, 2))) {
  stop("status_num은 0, 1, 2만 허용됩니다.")
}

if (anyDuplicated(analysis_df$subject_uid) > 0) {
  dup_ids <- unique(analysis_df$subject_uid[duplicated(analysis_df$subject_uid)])
  stop(
    "site + id 기준 중복이 존재합니다. 예시: ",
    paste(head(dup_ids, 5), collapse = ", "),
    if (length(dup_ids) > 5) " ..." else ""
  )
}

available_sites <- sort(unique(analysis_df$site_std))
if (!all(c("PNU", "SNU") %in% available_sites)) {
  stop(
    "site_std 정규화 후 PNU/SNU를 모두 찾지 못했습니다. 현재 site 값: ",
    paste(available_sites, collapse = ", "),
    "\nsite_label_map을 데이터에 맞게 수정하세요."
  )
}

analysis_df <- analysis_df %>%
  mutate(
    dataset_row_flag = case_when(
      site_std == "PNU" ~ "PNU",
      site_std == "SNU" ~ "SNU",
      TRUE ~ "OTHER"
    )
  )

# 🔴 패널: 데이터셋별·서브그룹별 분석 단위 생성 ===============================
analysis_sets <- list(
  PNU = analysis_df %>% filter(site_std == "PNU"),
  SNU = analysis_df %>% filter(site_std == "SNU"),
  MERGED = analysis_df
)

panel_definitions_list <- list()

for (dataset_name in names(analysis_sets)) {
  panel_definitions_list[[length(panel_definitions_list) + 1]] <- tibble(
    dataset_group = dataset_name,
    subgroup_var = "overall",
    subgroup_level = "overall"
  )
  
  if (length(subgroup_vars) > 0) {
    for (sg in subgroup_vars) {
      if (!sg %in% names(analysis_sets[[dataset_name]])) {
        warning("subgroup_vars에 지정한 컬럼이 존재하지 않아 건너뜁니다: ", sg)
        next
      }
      
      subgroup_levels <- analysis_sets[[dataset_name]] %>%
        filter(!is.na(.data[[sg]])) %>%
        mutate(.sg = as.character(.data[[sg]])) %>%
        distinct(.sg) %>%
        arrange(.sg) %>%
        pull(.sg)
      
      if (length(subgroup_levels) == 0) {
        next
      }
      
      panel_definitions_list[[length(panel_definitions_list) + 1]] <- tibble(
        dataset_group = dataset_name,
        subgroup_var = sg,
        subgroup_level = subgroup_levels
      )
    }
  }
}

panel_definitions <- bind_rows(panel_definitions_list) %>%
  mutate(
    panel_id = make.unique(make.names(paste(dataset_group, subgroup_var, subgroup_level, sep = "__"))),
    panel_label = ifelse(
      subgroup_var == "overall",
      dataset_group,
      paste0(dataset_group, " | ", subgroup_var, "=", subgroup_level)
    )
  )

# 🔴 계산: 패널별 follow-up maturity와 KM stability 산출 ===============================
panel_results <- vector("list", nrow(panel_definitions))

for (i in seq_len(nrow(panel_definitions))) {
  info <- panel_definitions[i, ]
  dataset_group_i <- info$dataset_group[[1]]
  subgroup_var_i <- info$subgroup_var[[1]]
  subgroup_level_i <- info$subgroup_level[[1]]
  
  panel_df <- analysis_sets[[dataset_group_i]]
  
  if (subgroup_var_i != "overall") {
    panel_df <- panel_df %>%
      filter(as.character(.data[[subgroup_var_i]]) == subgroup_level_i)
  }
  
  panel_info <- list(
    panel_id = info$panel_id[[1]],
    panel_label = info$panel_label[[1]],
    dataset_group = dataset_group_i,
    subgroup_var = subgroup_var_i,
    subgroup_level = subgroup_level_i,
    n_subjects = nrow(panel_df)
  )
  
  panel_results[[i]] <- analyze_panel(
    df_panel = panel_df,
    panel_info = panel_info
  )
}

panel_results <- purrr::compact(panel_results)

if (length(panel_results) == 0) {
  stop("분석 가능한 패널이 없습니다. subgroup_vars 또는 site 매핑을 확인하세요.")
}

# 🔴 결합: 표와 객체를 하나의 결과 묶음으로 정리 ===============================
panel_registry <- bind_rows(purrr::map(panel_results, "panel_registry")) %>%
  arrange(factor(dataset_group, levels = c("PNU", "SNU", "MERGED")), subgroup_var, subgroup_level)

followup_summary <- bind_rows(purrr::map(panel_results, "followup_summary")) %>%
  arrange(factor(dataset_group, levels = c("PNU", "SNU", "MERGED")), subgroup_var, subgroup_level)

km_yearly_summary <- bind_rows(purrr::map(panel_results, "km_yearly_summary")) %>%
  arrange(factor(dataset_group, levels = c("PNU", "SNU", "MERGED")), subgroup_var, subgroup_level, horizon_year)

timepoint_support_flags <- bind_rows(purrr::map(panel_results, "support_flags")) %>%
  arrange(factor(dataset_group, levels = c("PNU", "SNU", "MERGED")), subgroup_var, subgroup_level, horizon_year)

survival_curves <- bind_rows(purrr::map(panel_results, "survival_curves")) %>%
  arrange(factor(dataset_group, levels = c("PNU", "SNU", "MERGED")), subgroup_var, subgroup_level, curve_family, curve_name, time_days)

stability_gap_curves <- bind_rows(purrr::map(panel_results, "stability_gap_curves")) %>%
  arrange(factor(dataset_group, levels = c("PNU", "SNU", "MERGED")), subgroup_var, subgroup_level, time_days)

fit_objects <- purrr::map(panel_results, "fits")
fit_object_names <- purrr::map_chr(panel_results, ~ .x$panel_registry$panel_id[[1]])

if (length(fit_objects) != length(fit_object_names)) {
  stop("fit_objects와 fit_object_names의 길이가 다릅니다.")
}

if (anyDuplicated(fit_object_names) > 0) {
  dup_names <- unique(fit_object_names[duplicated(fit_object_names)])
  stop("fit_objects 이름이 중복됩니다: ", paste(dup_names, collapse = ", "))
}

fit_objects <- stats::setNames(fit_objects, fit_object_names)

analysis_bundle <- list(
  config = list(
    data_path = data_path,
    export_path = export_path,
    analysis_horizons_years = analysis_horizons_years,
    days_per_year = days_per_year,
    km_conf_type = km_conf_type,
    stable_gap_threshold = stable_gap_threshold,
    time_epsilon_days = time_epsilon_days,
    subgroup_vars = subgroup_vars,
    site_label_map = site_label_map,
    compute_ptfr_exact = compute_ptfr_exact
  ),
  analysis_data = analysis_df,
  panel_definitions = panel_definitions,
  panel_registry = panel_registry,
  followup_summary = followup_summary,
  km_yearly_summary = km_yearly_summary,
  timepoint_support_flags = timepoint_support_flags,
  survival_curves = survival_curves,
  stability_gap_curves = stability_gap_curves,
  fit_objects = fit_objects,
  session_info = utils::sessionInfo()
)

# 🔴 시각화: CSV 기반 그래프 데이터로 PNG 생성 ===============================
panel_order <- panel_registry$panel_label

survival_curves_plot <- survival_curves %>%
  mutate(panel_label = factor(panel_label, levels = panel_order))

stability_gap_curves_plot <- stability_gap_curves %>%
  mutate(panel_label = factor(panel_label, levels = panel_order))

timepoint_support_flags_plot <- timepoint_support_flags %>%
  mutate(
    panel_label = factor(panel_label, levels = panel_order),
    support_flag = factor(
      support_flag,
      levels = c("observed_stable", "observed_unstable", "beyond_reliable_tail")
    )
  )

plot_km_overlay <- ggplot(
  data = survival_curves_plot %>% filter(curve_family == "km_stability"),
  aes(x = time_years, y = surv, colour = curve_name, linetype = curve_name)
) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~ panel_label, scales = "free_x") +
  labs(
    title = "Step2: KM와 Betensky upper/lower stability overlay",
    x = "Time since cohort entry (years)",
    y = "Survival probability",
    colour = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

plot_followup_distributions <- ggplot(
  data = survival_curves_plot %>% filter(curve_family == "followup_distribution"),
  aes(x = time_years, y = surv, colour = curve_name, linetype = curve_name)
) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~ panel_label, scales = "free_x") +
  labs(
    title = "Step2: follow-up 관련 분포곡선",
    subtitle = "Reverse KM(C), observed T, C | C < X",
    x = "Time since cohort entry (years)",
    y = "Survival / retention probability",
    colour = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

plot_instability_gap <- ggplot(
  data = stability_gap_curves_plot,
  aes(x = time_years, y = upper_minus_lower)
) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~ panel_label, scales = "free_x") +
  labs(
    title = "Step2: Betensky instability gap (upper - lower)",
    x = "Time since cohort entry (years)",
    y = "Gap width"
  ) +
  theme_bw(base_size = 11)

plot_yearly_support_heatmap <- ggplot(
  data = timepoint_support_flags_plot,
  aes(x = factor(horizon_year), y = panel_label, fill = support_flag)
) +
  geom_tile() +
  labs(
    title = "Step2: 연도별 support flag",
    x = "Horizon (years)",
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

km_overlay_png <- make_export_file("step2_km_stability_overlay.png")
followup_dist_png <- make_export_file("step2_followup_distributions.png")
instability_gap_png <- make_export_file("step2_instability_gap.png")
support_heatmap_png <- make_export_file("step2_yearly_support_heatmap.png")

ggsave(filename = km_overlay_png, plot = plot_km_overlay, width = 12, height = 8, dpi = 300)
ggsave(filename = followup_dist_png, plot = plot_followup_distributions, width = 12, height = 8, dpi = 300)
ggsave(filename = instability_gap_png, plot = plot_instability_gap, width = 12, height = 8, dpi = 300)
ggsave(filename = support_heatmap_png, plot = plot_yearly_support_heatmap, width = 10, height = 5, dpi = 300)

# 🔴 저장: CSV와 RDS 및 파일 매니페스트 내보내기 ===============================
panel_registry_csv <- make_export_file("step2_panel_registry.csv")
followup_summary_csv <- make_export_file("step2_followup_summary.csv")
km_yearly_summary_csv <- make_export_file("step2_km_yearly_summary.csv")
support_flags_csv <- make_export_file("step2_timepoint_support_flags.csv")
survival_curves_csv <- make_export_file("step2_survival_curves.csv")
stability_gap_curves_csv <- make_export_file("step2_stability_gap_curves.csv")
fit_objects_rds <- make_export_file("step2_fit_objects.rds")
analysis_bundle_rds <- make_export_file("step2_analysis_bundle.rds")
file_manifest_csv <- make_export_file("step2_file_manifest.csv")

readr::write_csv(panel_registry, panel_registry_csv)
readr::write_csv(followup_summary, followup_summary_csv)
readr::write_csv(km_yearly_summary, km_yearly_summary_csv)
readr::write_csv(timepoint_support_flags, support_flags_csv)
readr::write_csv(survival_curves, survival_curves_csv)
readr::write_csv(stability_gap_curves, stability_gap_curves_csv)

saveRDS(fit_objects, fit_objects_rds)
saveRDS(analysis_bundle, analysis_bundle_rds)

file_manifest <- tibble(
  file_name = c(
    basename(panel_registry_csv),
    basename(followup_summary_csv),
    basename(km_yearly_summary_csv),
    basename(support_flags_csv),
    basename(survival_curves_csv),
    basename(stability_gap_curves_csv),
    basename(fit_objects_rds),
    basename(analysis_bundle_rds),
    basename(km_overlay_png),
    basename(followup_dist_png),
    basename(instability_gap_png),
    basename(support_heatmap_png),
    basename(file_manifest_csv)
  ),
  file_path = c(
    panel_registry_csv,
    followup_summary_csv,
    km_yearly_summary_csv,
    support_flags_csv,
    survival_curves_csv,
    stability_gap_curves_csv,
    fit_objects_rds,
    analysis_bundle_rds,
    km_overlay_png,
    followup_dist_png,
    instability_gap_png,
    support_heatmap_png,
    file_manifest_csv
  ),
  file_type = c(
    "csv", "csv", "csv", "csv", "csv", "csv",
    "rds", "rds",
    "png", "png", "png", "png",
    "csv"
  ),
  description = c(
    "패널별 표본 수와 사건/검열 개요",
    "reverse KM, plateau, Betensky instability area 등 panel-level 요약",
    "1-10년 KM survival/risk/n.risk 요약",
    "연도별 support flag와 completeness proxy(percentage/CCI/SPT)",
    "KM/Betensky/reverse KM/observation/event-free follow-up 곡선 원자료",
    "Betensky gap 및 partial gap 곡선 원자료",
    "KM/reverse KM/Betensky upper-lower survfit 객체",
    "Step2 전체 결과 묶음(설정, 표, 곡선, survfit 객체 포함)",
    "KM와 Betensky upper/lower overlay 그림",
    "follow-up 관련 분포곡선 그림",
    "Betensky instability gap 그림",
    "연도별 support flag heatmap",
    "생성 파일 목록"
  )
)

readr::write_csv(file_manifest, file_manifest_csv)

file_manifest <- file_manifest %>%
  mutate(
    exists = file.exists(file_path),
    size_bytes = ifelse(exists, file.info(file_path)$size, NA_real_)
  )

readr::write_csv(file_manifest, file_manifest_csv)

message("Step2 산출물이 저장되었습니다: ", normalizePath(export_path, winslash = "/"))