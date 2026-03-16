# 🔴 설정: 경로와 전역옵션 ===============================
## 🟠 불러오기: 패키지와 시드 ===============================
suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(stringr)
})

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260301)

## 🟠 지정: 입력경로와 출력경로 ===============================
data_path <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
export_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦2.Step2_follow-up maturity, KM tail stability/attachments'

dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

## 🟠 지정: Step2 분석 파라미터 ===============================
year_days <- 365.25
pred_years <- 1:10
subgroup_vars <- character(0)   # 예: c("sex_fact")
stable_gap_threshold <- 0.05    # Betensky upper-lower limit 차이가 이 값 이하이면 상대적으로 안정적이라고 표시
min_n_risk_stable <- 10L        # 안정적이라고 표시하기 위한 최소 risk set
include_site_specific_branches <- TRUE
include_merged_branch <- TRUE

## 🟠 지정: 출력 파일명 ===============================
file_fu_summary_csv            <- file.path(export_path, "step2_fu_summary.csv")
file_yearly_summary_csv        <- file.path(export_path, "step2_yearly_timepoint_summary.csv")
file_curves_limits_csv         <- file.path(export_path, "step2_curves_observed_limits.csv")
file_curves_reverse_csv        <- file.path(export_path, "step2_curves_reversekm.csv")
file_analysis_index_csv        <- file.path(export_path, "step2_analysis_index.csv")
file_plot_limits_png           <- file.path(export_path, "step2_plot_observed_limits.png")
file_plot_reverse_png          <- file.path(export_path, "step2_plot_reversekm.png")
file_fit_objects_rds           <- file.path(export_path, "step2_fit_objects.rds")

# 🔴 정의: 보조함수 ===============================
## 🟠 확인: 필수컬럼 검증 ===============================
check_required_columns <- function(dat, required_cols) {
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("필수 컬럼이 없습니다: ", paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}

## 🟠 변환: transition-only endpoint 구성 ===============================
prepare_analysis_data <- function(dat, year_days) {
  dat %>%
    mutate(
      site = as.character(site),
      id = as.character(id),
      days_followup = as.numeric(days_followup),
      status_num = as.integer(status_num),
      subject_key = paste(site, id, sep = "::"),
      event_transition = if_else(status_num == 1L, 1L, 0L),
      time_years = days_followup / year_days
    )
}

## 🟠 생성: 데이터셋 브랜치 목록 ===============================
make_branch_list <- function(dat, include_merged = TRUE, include_site_specific = TRUE) {
  out <- list()
  
  if (include_merged) {
    out[["MERGED"]] <- dat
  }
  
  if (include_site_specific) {
    site_values <- sort(unique(dat$site))
    for (sv in site_values) {
      out[[sv]] <- dat %>% filter(site == sv)
    }
  }
  
  out
}

## 🟠 생성: 서브그룹 스트라타 목록 ===============================
make_strata_table <- function(dat, subgroup_vars = character(0)) {
  base_tbl <- tibble(
    subgroup_label = "Overall",
    subgroup_key   = "overall",
    data           = list(dat)
  )
  
  if (length(subgroup_vars) == 0) {
    return(base_tbl)
  }
  
  missing_subgroups <- setdiff(subgroup_vars, names(dat))
  if (length(missing_subgroups) > 0) {
    stop("subgroup_vars에 지정된 컬럼이 데이터에 없습니다: ",
         paste(missing_subgroups, collapse = ", "))
  }
  
  dat2 <- dat %>%
    mutate(across(all_of(subgroup_vars), ~ ifelse(is.na(.x), "NA", as.character(.x))))
  
  split_list <- split(dat2, interaction(dat2[subgroup_vars], drop = TRUE, lex.order = TRUE))
  
  subgroup_tbl <- purrr::imap_dfr(split_list, function(df_i, nm_i) {
    vals <- as.list(df_i[1, subgroup_vars, drop = FALSE])
    subgroup_label <- paste(
      paste0(names(vals), "=", unlist(vals, use.names = FALSE)),
      collapse = " | "
    )
    
    tibble(
      subgroup_label = subgroup_label,
      subgroup_key   = make.names(subgroup_label),
      data           = list(df_i)
    )
  })
  
  bind_rows(base_tbl, subgroup_tbl)
}

## 🟠 계산: KM stability limit용 시간보정 ===============================
calc_small_shift <- function(x) {
  ux <- sort(unique(as.numeric(x)))
  diffs <- diff(ux)
  diffs <- diffs[is.finite(diffs) & diffs > 0]
  
  if (length(diffs) == 0) {
    return(1 / year_days / 1000)
  }
  
  min(diffs) / 1000
}

build_upper_limit_data <- function(dat) {
  times  <- dat$time_years
  events <- dat$event_transition
  ev_times <- sort(unique(times[events == 1]))
  
  shift <- calc_small_shift(times)
  
  out <- dat
  out$time_upper  <- times
  out$event_upper <- events
  
  if (length(ev_times) > 0) {
    max_event <- max(ev_times)
    cens_idx  <- which(events == 0)
    
    if (length(cens_idx) > 0) {
      out$time_upper[cens_idx] <- pmax(times[cens_idx], max_event) + shift
    }
  }
  
  out
}

build_lower_limit_data <- function(dat) {
  times  <- dat$time_years
  events <- dat$event_transition
  ev_times <- sort(unique(times[events == 1]))
  
  shift <- calc_small_shift(times)
  
  out <- dat
  out$time_lower  <- times
  out$event_lower <- events
  
  if (length(ev_times) == 0) {
    return(out)
  }
  
  cens_idx <- which(events == 0)
  
  if (length(cens_idx) > 0) {
    next_event_time <- vapply(times[cens_idx], function(tt) {
      cand <- ev_times[ev_times > tt]
      if (length(cand) > 0) {
        cand[1]
      } else {
        tt + shift
      }
    }, numeric(1))
    
    next_event_time <- pmax(next_event_time, times[cens_idx] + shift)
    
    out$time_lower[cens_idx]  <- next_event_time
    out$event_lower[cens_idx] <- 1L
  }
  
  out
}

## 🟠 적합: survfit 객체 생성 ===============================
fit_surv_objects <- function(dat) {
  observed_fit <- survival::survfit(
    survival::Surv(time_years, event_transition) ~ 1,
    data = dat,
    conf.type = "log-log"
  )
  
  reverse_fit <- survival::survfit(
    survival::Surv(time_years, 1L - event_transition) ~ 1,
    data = dat,
    conf.type = "log-log"
  )
  
  upper_dat <- build_upper_limit_data(dat)
  upper_fit <- survival::survfit(
    survival::Surv(time_upper, event_upper) ~ 1,
    data = upper_dat,
    conf.type = "log-log"
  )
  
  lower_dat <- build_lower_limit_data(dat)
  lower_fit <- survival::survfit(
    survival::Surv(time_lower, event_lower) ~ 1,
    data = lower_dat,
    conf.type = "log-log"
  )
  
  list(
    observed_fit = observed_fit,
    reverse_fit  = reverse_fit,
    upper_fit    = upper_fit,
    lower_fit    = lower_fit,
    upper_data   = upper_dat,
    lower_data   = lower_dat
  )
}

## 🟠 정리: survfit 결과 tidy 변환 ===============================
vec_or_default <- function(x, n, default = NA_real_) {
  if (length(x) == n) return(as.numeric(x))
  if (length(x) == 0) return(rep(default, n))
  rep(as.numeric(x[1]), n)
}

extract_surv_curve <- function(fit, n_start) {
  ss <- summary(fit)
  
  if (length(ss$time) == 0) {
    return(
      tibble(
        time = 0,
        surv = 1,
        std_err = 0,
        lower = 1,
        upper = 1,
        n_risk = n_start,
        n_event = 0,
        n_censor = 0
      )
    )
  }
  
  n_time <- length(ss$time)
  
  tibble(
    time     = c(0, ss$time),
    surv     = c(1, ss$surv),
    std_err  = c(0, vec_or_default(ss$std.err, n_time, default = NA_real_)),
    lower    = c(1, vec_or_default(ss$lower,   n_time, default = NA_real_)),
    upper    = c(1, vec_or_default(ss$upper,   n_time, default = NA_real_)),
    n_risk   = c(n_start, vec_or_default(ss$n.risk,   n_time, default = NA_real_)),
    n_event  = c(0, vec_or_default(ss$n.event,  n_time, default = 0)),
    n_censor = c(0, vec_or_default(ss$n.censor, n_time, default = 0))
  )
}

evaluate_surv_at_times <- function(fit, times) {
  ss <- summary(fit, times = times, extend = TRUE)
  n <- length(times)
  
  tibble(
    time     = times,
    surv     = vec_or_default(ss$surv, n, default = 1),
    std_err  = vec_or_default(ss$std.err, n, default = NA_real_),
    lower    = vec_or_default(ss$lower, n, default = NA_real_),
    upper    = vec_or_default(ss$upper, n, default = NA_real_),
    n_risk   = vec_or_default(ss$n.risk, n, default = NA_real_),
    n_event  = vec_or_default(ss$n.event, n, default = 0),
    n_censor = vec_or_default(ss$n.censor, n, default = 0)
  )
}

extract_survfit_median <- function(fit) {
  tab <- fit$table
  if (is.null(tab)) return(NA_real_)
  
  med <- if (is.matrix(tab)) {
    tab[1, "median"]
  } else {
    unname(tab["median"])
  }
  
  med <- as.numeric(med)
  if (length(med) == 0 || !is.finite(med)) return(NA_real_)
  med
}

calc_step_area <- function(x, y) {
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  
  keep <- !duplicated(x)
  x <- x[keep]
  y <- y[keep]
  
  if (length(x) <= 1) return(0)
  
  sum(diff(x) * head(y, -1), na.rm = TRUE)
}

## 🟠 계산: 연도별 요약과 안정성 플래그 ===============================
make_yearly_summary <- function(obs_fit,
                                rev_fit,
                                upper_fit,
                                lower_fit,
                                pred_years,
                                last_event_time,
                                last_followup_time,
                                censored_after_last_event,
                                stable_gap_threshold,
                                min_n_risk_stable) {
  obs_yearly <- evaluate_surv_at_times(obs_fit,   pred_years) %>%
    rename_with(~ paste0(.x, "_obs"), -time)
  
  rev_yearly <- evaluate_surv_at_times(rev_fit,   pred_years) %>%
    rename_with(~ paste0(.x, "_rev"), -time)
  
  up_yearly <- evaluate_surv_at_times(upper_fit,  pred_years) %>%
    rename_with(~ paste0(.x, "_upper"), -time)
  
  low_yearly <- evaluate_surv_at_times(lower_fit, pred_years) %>%
    rename_with(~ paste0(.x, "_lower"), -time)
  
  tibble(time = pred_years) %>%
    left_join(obs_yearly, by = "time") %>%
    left_join(rev_yearly, by = "time") %>%
    left_join(up_yearly,  by = "time") %>%
    left_join(low_yearly, by = "time") %>%
    mutate(
      S_KM        = surv_obs,
      Risk_KM     = 1 - surv_obs,
      S_upper     = surv_upper,
      S_lower     = surv_lower,
      stability_gap = pmax(S_upper - S_lower, 0),
      within_observed_followup = time <= last_followup_time,
      beyond_last_event        = ifelse(is.na(last_event_time), NA, time > last_event_time),
      support_flag = case_when(
        time > last_followup_time ~ "beyond_reliable_tail",
        is.na(last_event_time) ~ "observed_unstable",
        time > last_event_time & censored_after_last_event > 0L ~ "observed_unstable",
        stability_gap <= stable_gap_threshold & n_risk_obs >= min_n_risk_stable ~ "observed_stable",
        TRUE ~ "observed_unstable"
      )
    )
}

# 🔴 입력: 원자료와 사전검증 ===============================
## 🟠 읽기: MERGED dataset3 ===============================
raw_dat <- readr::read_csv(
  file = data_path,
  show_col_types = FALSE,
  progress = FALSE
)

## 🟠 검사: 분석 필수조건 ===============================
required_cols <- c("id", "site", "days_followup", "status_num")
check_required_columns(raw_dat, required_cols)

# 🔴 변환: Step2 분석데이터 구성 ===============================
## 🟠 파생: transition-only endpoint 데이터 ===============================
dat_analysis <- prepare_analysis_data(raw_dat, year_days = year_days)

if (any(is.na(dat_analysis$days_followup))) {
  stop("days_followup에 NA가 있습니다. Step2 이전에 정리해야 합니다.")
}

if (any(dat_analysis$days_followup < 0, na.rm = TRUE)) {
  stop("days_followup에 음수가 있습니다. Step2 이전에 정리해야 합니다.")
}

if (any(duplicated(dat_analysis$subject_key))) {
  dup_keys <- dat_analysis$subject_key[duplicated(dat_analysis$subject_key)]
  stop("site + id 조합이 유일하지 않습니다. 중복 key 예시: ", paste(head(unique(dup_keys), 5), collapse = ", "))
}

## 🟠 구성: 브랜치와 스트라타 인덱스 ===============================
branch_list <- make_branch_list(
  dat = dat_analysis,
  include_merged = include_merged_branch,
  include_site_specific = include_site_specific_branches
)

analysis_index <- purrr::imap_dfr(branch_list, function(branch_df, branch_name) {
  strata_tbl <- make_strata_table(branch_df, subgroup_vars = subgroup_vars)
  strata_tbl %>%
    mutate(dataset = branch_name) %>%
    select(dataset, subgroup_label, subgroup_key, data)
}) %>%
  mutate(
    analysis_id = make.unique(paste0(dataset, "__", subgroup_key)),
    analysis_label = if_else(subgroup_label == "Overall",
                             dataset,
                             paste0(dataset, " | ", subgroup_label))
  ) %>%
  select(analysis_id, analysis_label, dataset, subgroup_label, subgroup_key, data)

analysis_index_export <- analysis_index %>%
  select(-data)

# 🔴 수행: Step2 follow-up maturity 분석 ===============================
## 🟠 반복: 브랜치별 적합과 요약 ===============================
result_list <- vector("list", length = nrow(analysis_index))
names(result_list) <- analysis_index$analysis_id

for (ii in seq_len(nrow(analysis_index))) {
  ai <- analysis_index[ii, ]
  dat_i <- ai$data[[1]]
  
  n_subject <- nrow(dat_i)
  n_transition <- sum(dat_i$event_transition == 1L, na.rm = TRUE)
  n_censor_transition <- sum(dat_i$event_transition == 0L, na.rm = TRUE)
  n_remission_as_censor <- sum(dat_i$status_num == 2L, na.rm = TRUE)
  
  fits_i <- fit_surv_objects(dat_i)
  
  last_event_time <- if (n_transition > 0) {
    max(dat_i$time_years[dat_i$event_transition == 1L], na.rm = TRUE)
  } else {
    NA_real_
  }
  
  last_followup_time <- max(dat_i$time_years, na.rm = TRUE)
  
  plateau_length <- if (is.na(last_event_time)) {
    NA_real_
  } else {
    last_followup_time - last_event_time
  }
  
  n_risk_at_last_event <- if (is.na(last_event_time)) {
    NA_real_
  } else {
    tmp <- suppressWarnings(summary(fits_i$observed_fit, times = last_event_time, extend = TRUE)$n.risk)
    if (length(tmp) == 0) NA_real_ else as.numeric(tmp[1])
  }
  
  censored_after_last_event <- if (is.na(last_event_time)) {
    NA_integer_
  } else {
    sum(dat_i$event_transition == 0L & dat_i$time_years > last_event_time, na.rm = TRUE)
  }
  
  reverse_km_median_followup <- extract_survfit_median(fits_i$reverse_fit)
  
  obs_curve_i <- extract_surv_curve(fits_i$observed_fit, n_start = n_subject) %>%
    mutate(curve_type = "Observed_KM")
  
  upper_curve_i <- extract_surv_curve(fits_i$upper_fit, n_start = n_subject) %>%
    mutate(curve_type = "Upper_Stability_Limit")
  
  lower_curve_i <- extract_surv_curve(fits_i$lower_fit, n_start = n_subject) %>%
    mutate(curve_type = "Lower_Stability_Limit")
  
  reverse_curve_i <- extract_surv_curve(fits_i$reverse_fit, n_start = n_subject) %>%
    mutate(curve_type = "Reverse_KM")
  
  curves_obs_limits_i <- bind_rows(obs_curve_i, upper_curve_i, lower_curve_i) %>%
    mutate(
      analysis_id    = ai$analysis_id,
      analysis_label = ai$analysis_label,
      dataset        = ai$dataset,
      subgroup_label = ai$subgroup_label,
      subgroup_key   = ai$subgroup_key
    ) %>%
    select(analysis_id, analysis_label, dataset, subgroup_label, subgroup_key,
           curve_type, everything())
  
  curves_reverse_i <- reverse_curve_i %>%
    mutate(
      analysis_id    = ai$analysis_id,
      analysis_label = ai$analysis_label,
      dataset        = ai$dataset,
      subgroup_label = ai$subgroup_label,
      subgroup_key   = ai$subgroup_key
    ) %>%
    select(analysis_id, analysis_label, dataset, subgroup_label, subgroup_key,
           curve_type, everything())
  
  grid_auc <- sort(unique(c(
    curves_obs_limits_i$time[curves_obs_limits_i$curve_type == "Upper_Stability_Limit"],
    curves_obs_limits_i$time[curves_obs_limits_i$curve_type == "Lower_Stability_Limit"],
    0, last_followup_time
  )))
  
  upper_eval_grid <- evaluate_surv_at_times(fits_i$upper_fit, grid_auc)
  lower_eval_grid <- evaluate_surv_at_times(fits_i$lower_fit, grid_auc)
  
  betensky_diff_curve <- pmax(upper_eval_grid$surv - lower_eval_grid$surv, 0)
  norm_denom <- if (!is.na(last_event_time) && is.finite(last_event_time) && last_event_time > 0) {
    last_event_time
  } else if (is.finite(last_followup_time) && last_followup_time > 0) {
    last_followup_time
  } else {
    1
  }
  
  betensky_area_norm <- calc_step_area(grid_auc, betensky_diff_curve) / norm_denom
  
  yearly_i <- make_yearly_summary(
    obs_fit = fits_i$observed_fit,
    rev_fit = fits_i$reverse_fit,
    upper_fit = fits_i$upper_fit,
    lower_fit = fits_i$lower_fit,
    pred_years = pred_years,
    last_event_time = last_event_time,
    last_followup_time = last_followup_time,
    censored_after_last_event = ifelse(is.na(censored_after_last_event), 0L, censored_after_last_event),
    stable_gap_threshold = stable_gap_threshold,
    min_n_risk_stable = min_n_risk_stable
  ) %>%
    mutate(
      analysis_id    = ai$analysis_id,
      analysis_label = ai$analysis_label,
      dataset        = ai$dataset,
      subgroup_label = ai$subgroup_label,
      subgroup_key   = ai$subgroup_key
    ) %>%
    select(
      analysis_id, analysis_label, dataset, subgroup_label, subgroup_key,
      year = time,
      S_KM, Risk_KM, lower_obs, upper_obs, n_risk_obs, n_event_obs, n_censor_obs,
      surv_rev, lower_rev, upper_rev, n_risk_rev,
      S_upper, S_lower, stability_gap, within_observed_followup, beyond_last_event,
      support_flag
    )
  
  max_stable_year <- suppressWarnings(max(yearly_i$year[yearly_i$support_flag == "observed_stable"], na.rm = TRUE))
  if (!is.finite(max_stable_year)) max_stable_year <- NA_real_
  
  max_observed_year <- suppressWarnings(max(yearly_i$year[yearly_i$support_flag != "beyond_reliable_tail"], na.rm = TRUE))
  if (!is.finite(max_observed_year)) max_observed_year <- NA_real_
  
  fu_summary_i <- tibble(
    analysis_id                 = ai$analysis_id,
    analysis_label              = ai$analysis_label,
    dataset                     = ai$dataset,
    subgroup_label              = ai$subgroup_label,
    subgroup_key                = ai$subgroup_key,
    n_subject                   = n_subject,
    n_transition                = n_transition,
    n_censor_transition         = n_censor_transition,
    n_remission_as_censor       = n_remission_as_censor,
    transition_rate             = n_transition / n_subject,
    reverseKM_median_followup   = reverse_km_median_followup,
    last_event_time_years       = last_event_time,
    last_followup_time_years    = last_followup_time,
    plateau_length_years        = plateau_length,
    n_risk_at_last_event        = n_risk_at_last_event,
    censored_after_last_event   = censored_after_last_event,
    betensky_area_norm          = betensky_area_norm,
    stable_gap_threshold        = stable_gap_threshold,
    min_n_risk_stable           = min_n_risk_stable,
    max_observed_year           = max_observed_year,
    max_stable_year             = max_stable_year
  )
  
  result_list[[ii]] <- list(
    fu_summary            = fu_summary_i,
    yearly_summary        = yearly_i,
    curves_obs_limits     = curves_obs_limits_i,
    curves_reverse        = curves_reverse_i,
    fit_objects           = list(
      observed_fit = fits_i$observed_fit,
      reverse_fit  = fits_i$reverse_fit,
      upper_fit    = fits_i$upper_fit,
      lower_fit    = fits_i$lower_fit
    ),
    fit_metadata          = list(
      last_event_time = last_event_time,
      last_followup_time = last_followup_time,
      plateau_length = plateau_length,
      n_risk_at_last_event = n_risk_at_last_event,
      censored_after_last_event = censored_after_last_event,
      reverseKM_median_followup = reverse_km_median_followup,
      betensky_area_norm = betensky_area_norm
    )
  )
}

# 🔴 정리: Step2 산출물 데이터프레임 ===============================
## 🟠 결합: 요약표와 곡선표 ===============================
fu_summary_df <- bind_rows(lapply(result_list, function(x) x$fu_summary))
yearly_summary_df <- bind_rows(lapply(result_list, function(x) x$yearly_summary))
curves_obs_limits_df <- bind_rows(lapply(result_list, function(x) x$curves_obs_limits))
curves_reverse_df <- bind_rows(lapply(result_list, function(x) x$curves_reverse))

## 🟠 생성: 시각화용 데이터 준비 ===============================
curves_obs_limits_plot_df <- curves_obs_limits_df %>%
  mutate(
    curve_type = factor(
      curve_type,
      levels = c("Observed_KM", "Upper_Stability_Limit", "Lower_Stability_Limit")
    )
  )

curves_reverse_plot_df <- curves_reverse_df %>%
  mutate(
    curve_type = factor(curve_type, levels = c("Reverse_KM"))
  )

# 🔴 저장: Step2 결과물 export ===============================
## 🟠 쓰기: CSV와 RDS ===============================
readr::write_csv(analysis_index_export, file_analysis_index_csv, na = "")
readr::write_csv(fu_summary_df,        file_fu_summary_csv,     na = "")
readr::write_csv(yearly_summary_df,    file_yearly_summary_csv, na = "")
readr::write_csv(curves_obs_limits_df, file_curves_limits_csv,  na = "")
readr::write_csv(curves_reverse_df,    file_curves_reverse_csv, na = "")

step2_bundle <- list(
  config = list(
    data_path = data_path,
    export_path = export_path,
    year_days = year_days,
    pred_years = pred_years,
    subgroup_vars = subgroup_vars,
    stable_gap_threshold = stable_gap_threshold,
    min_n_risk_stable = min_n_risk_stable,
    include_site_specific_branches = include_site_specific_branches,
    include_merged_branch = include_merged_branch
  ),
  analysis_input_data = dat_analysis,
  analysis_index = analysis_index_export,
  fit_objects = lapply(result_list, function(x) x$fit_objects),
  fit_metadata = lapply(result_list, function(x) x$fit_metadata),
  outputs = list(
    fu_summary = fu_summary_df,
    yearly_summary = yearly_summary_df,
    curves_obs_limits = curves_obs_limits_df,
    curves_reverse = curves_reverse_df
  )
)

saveRDS(step2_bundle, file = file_fit_objects_rds)

## 🟠 그리기: PNG 시각화 ===============================
n_facets <- length(unique(curves_obs_limits_plot_df$analysis_label))
plot_width_limits <- max(12, 4 * ceiling(sqrt(n_facets)))
plot_height_limits <- max(8,  3 * ceiling(n_facets / ceiling(sqrt(n_facets))))

p_limits <- ggplot(
  curves_obs_limits_plot_df,
  aes(x = time, y = surv, color = curve_type, linetype = curve_type)
) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~ analysis_label, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(name = "Time since cohort entry (years)") +
  labs(
    title = "Step2: Observed KM and Betensky-style stability limits",
    y = "Survival probability",
    color = "Curve",
    linetype = "Curve"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file_plot_limits_png,
  plot = p_limits,
  width = plot_width_limits,
  height = plot_height_limits,
  units = "in",
  dpi = 320
)

n_facets_rev <- length(unique(curves_reverse_plot_df$analysis_label))
plot_width_rev <- max(12, 4 * ceiling(sqrt(n_facets_rev)))
plot_height_rev <- max(8,  3 * ceiling(n_facets_rev / ceiling(sqrt(n_facets_rev))))

p_reverse <- ggplot(
  curves_reverse_plot_df,
  aes(x = time, y = surv, color = curve_type, linetype = curve_type)
) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~ analysis_label, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(name = "Time since cohort entry (years)") +
  labs(
    title = "Step2: Reverse Kaplan-Meier for follow-up maturity",
    y = "Survival probability of censoring time",
    color = "Curve",
    linetype = "Curve"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file_plot_reverse_png,
  plot = p_reverse,
  width = plot_width_rev,
  height = plot_height_rev,
  units = "in",
  dpi = 320
)