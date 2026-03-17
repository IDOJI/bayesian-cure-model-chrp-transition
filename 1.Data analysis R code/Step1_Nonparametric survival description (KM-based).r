# 🔴 설정: 경로와 핵심 옵션 ===============================
script_name <- "step1_plain_km"
script_version <- "v1.0.0"

data_path <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
export_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦1.Step1_KM/attachments'

site_pnu_values <- c("PNU")
site_snu_values <- c("SNU")  # 예: 실제 라벨이 CHR-P라면 c("SNU", "CHR-P")로 수정
time_grid_years <- 1:10
days_per_year <- 365.25

transition_code <- 1L
right_censor_code <- 0L
remission_code <- 2L
allowed_status_codes <- c(right_censor_code, transition_code, remission_code)

conf_type_km <- "log-log"
missing_subgroup_label <- "(Missing)"
subgroup_specs <- c(sex = "sex_fact")  # overall만 원하면 character(0)로 수정
plot_dpi <- 300

# 🟠 준비: 패키지와 실행 옵션 ===============================
required_packages <- c("survival", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
  stop(
    "다음 패키지를 먼저 설치하세요: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

invisible(lapply(required_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

options(stringsAsFactors = FALSE)
if (!file.exists(data_path)) {
  stop("data_path에 해당 파일이 없습니다: ", data_path, call. = FALSE)
}
if (!dir.exists(export_path)) {
  dir.create(export_path, recursive = TRUE, showWarnings = FALSE)
}

# 🔴 함수: 입력 검증과 KM 추출 ===============================
## 🟠 함수: 공통 유틸리티 정의 ===============================
collapse_value <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(x) == 0L) return("")
  paste(x, collapse = "; ")
}

clamp01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

write_csv_utf8 <- function(x, path) {
  utils::write.csv(
    x,
    file = path,
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
}

expand_component <- function(x, default, n_out) {
  if (is.null(x) || length(x) == 0L) {
    return(rep(default, n_out))
  }
  rep_len(x, n_out)
}

extract_survfit_table <- function(fit, times_vec) {
  if (length(times_vec) == 0L) {
    return(
      data.frame(
        time_year = numeric(0),
        survival = numeric(0),
        lower = numeric(0),
        upper = numeric(0),
        n_risk = integer(0),
        n_event = integer(0),
        n_censor = integer(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  s <- summary(fit, times = times_vec, extend = TRUE)
  n_out <- length(times_vec)
  
  surv_vals <- clamp01(expand_component(s$surv, 1, n_out))
  lower_vals <- expand_component(s$lower, NA_real_, n_out)
  upper_vals <- expand_component(s$upper, NA_real_, n_out)
  
  lower_vals[is.na(lower_vals)] <- surv_vals[is.na(lower_vals)]
  upper_vals[is.na(upper_vals)] <- surv_vals[is.na(upper_vals)]
  
  out <- data.frame(
    time_year = times_vec,
    survival = clamp01(surv_vals),
    lower = clamp01(lower_vals),
    upper = clamp01(upper_vals),
    n_risk = as.integer(round(expand_component(s$n.risk, 0, n_out))),
    n_event = as.integer(round(expand_component(s$n.event, 0, n_out))),
    n_censor = as.integer(round(expand_component(s$n.censor, 0, n_out))),
    stringsAsFactors = FALSE
  )
  
  out
}

## 🟠 함수: 품질점검 행 생성 ===============================
make_qc_row <- function(df, qc_scope, dataset_name) {
  data.frame(
    qc_scope = qc_scope,
    dataset = dataset_name,
    n_rows = nrow(df),
    n_unique_subject_uid = length(unique(df$subject_uid)),
    n_duplicate_subject_uid = sum(duplicated(df$subject_uid)),
    n_transition_raw = sum(df$status_num == transition_code, na.rm = TRUE),
    n_right_censor_raw = sum(df$status_num == right_censor_code, na.rm = TRUE),
    n_remission_raw = sum(df$status_num == remission_code, na.rm = TRUE),
    n_transition_analysis = sum(df$event_transition == 1L, na.rm = TRUE),
    n_censor_analysis = sum(df$event_transition == 0L, na.rm = TRUE),
    followup_min_days = min(df$days_followup, na.rm = TRUE),
    followup_median_days = stats::median(df$days_followup, na.rm = TRUE),
    followup_max_days = max(df$days_followup, na.rm = TRUE),
    followup_min_years = min(df$time_years, na.rm = TRUE),
    followup_median_years = stats::median(df$time_years, na.rm = TRUE),
    followup_max_years = max(df$time_years, na.rm = TRUE),
    site_levels = collapse_value(sort(unique(df$site))),
    stringsAsFactors = FALSE
  )
}

## 🟠 함수: 단일 곡선 KM 적합 ===============================
fit_km_curve <- function(df_curve, curve_id, dataset_name, subgroup_name, subgroup_var, subgroup_level) {
  n_total <- nrow(df_curve)
  n_event <- sum(df_curve$event_transition == 1L, na.rm = TRUE)
  n_censor <- sum(df_curve$event_transition == 0L, na.rm = TRUE)
  n_right_censor <- sum(df_curve$status_num == right_censor_code, na.rm = TRUE)
  n_remission <- sum(df_curve$status_num == remission_code, na.rm = TRUE)
  
  max_followup_years <- max(df_curve$time_years, na.rm = TRUE)
  min_followup_years <- min(df_curve$time_years, na.rm = TRUE)
  median_followup_years <- stats::median(df_curve$time_years, na.rm = TRUE)
  
  last_event_time_years <- if (n_event > 0L) {
    max(df_curve$time_years[df_curve$event_transition == 1L], na.rm = TRUE)
  } else {
    NA_real_
  }
  
  censored_after_last_event <- if (is.na(last_event_time_years)) {
    n_censor
  } else {
    sum(df_curve$event_transition == 0L & df_curve$time_years > last_event_time_years, na.rm = TRUE)
  }
  
  fit_obj <- survival::survfit(
    survival::Surv(time_years, event_transition) ~ 1,
    data = df_curve,
    conf.type = conf_type_km
  )
  
  n_risk_at_last_event <- if (is.na(last_event_time_years)) {
    NA_integer_
  } else {
    as.integer(round(summary(fit_obj, times = last_event_time_years, extend = TRUE)$n.risk[1]))
  }
  
  observed_times <- sort(unique(df_curve$time_years))
  observed_times_positive <- observed_times[observed_times > 0]
  
  curve_points <- extract_survfit_table(fit_obj, observed_times_positive)
  curve_data <- rbind(
    data.frame(
      time_year = 0,
      survival = 1,
      lower = 1,
      upper = 1,
      n_risk = as.integer(n_total),
      n_event = 0L,
      n_censor = 0L,
      stringsAsFactors = FALSE
    ),
    curve_points
  )
  
  curve_data$curve_id <- curve_id
  curve_data$dataset <- dataset_name
  curve_data$subgroup_name <- subgroup_name
  curve_data$subgroup_var <- subgroup_var
  curve_data$subgroup_level <- subgroup_level
  curve_data$risk <- 1 - curve_data$survival
  curve_data$survival_lcl <- curve_data$lower
  curve_data$survival_ucl <- curve_data$upper
  curve_data$risk_lcl <- 1 - curve_data$upper
  curve_data$risk_ucl <- 1 - curve_data$lower
  curve_data$cum_event <- cumsum(curve_data$n_event)
  curve_data$cum_censor <- cumsum(curve_data$n_censor)
  curve_data$curve_label <- if (subgroup_name == "overall") {
    paste0(dataset_name, " | overall")
  } else {
    paste0(dataset_name, " | ", subgroup_name, " | ", subgroup_level)
  }
  
  yearly_points <- extract_survfit_table(fit_obj, time_grid_years)
  yearly_points$curve_id <- curve_id
  yearly_points$dataset <- dataset_name
  yearly_points$subgroup_name <- subgroup_name
  yearly_points$subgroup_var <- subgroup_var
  yearly_points$subgroup_level <- subgroup_level
  yearly_points$time_year_label <- paste0(yearly_points$time_year, "y")
  yearly_points$S_KM <- yearly_points$survival
  yearly_points$Risk_KM <- 1 - yearly_points$survival
  yearly_points$S_KM_LCL <- yearly_points$lower
  yearly_points$S_KM_UCL <- yearly_points$upper
  yearly_points$Risk_KM_LCL <- 1 - yearly_points$upper
  yearly_points$Risk_KM_UCL <- 1 - yearly_points$lower
  yearly_points$time_supported_by_followup <- yearly_points$time_year <= max_followup_years
  yearly_points$time_beyond_last_event <- if (is.na(last_event_time_years)) {
    FALSE
  } else {
    yearly_points$time_year > last_event_time_years
  }
  yearly_points$n_total <- n_total
  yearly_points$n_event_total <- n_event
  yearly_points$n_censor_total <- n_censor
  yearly_points$followup_max_years <- max_followup_years
  yearly_points$last_event_time_years <- last_event_time_years
  yearly_points <- yearly_points[
    ,
    c(
      "curve_id", "dataset", "subgroup_name", "subgroup_var", "subgroup_level",
      "time_year", "time_year_label",
      "S_KM", "Risk_KM", "S_KM_LCL", "S_KM_UCL", "Risk_KM_LCL", "Risk_KM_UCL",
      "n_risk", "n_event", "n_censor",
      "time_supported_by_followup", "time_beyond_last_event",
      "n_total", "n_event_total", "n_censor_total",
      "followup_max_years", "last_event_time_years"
    )
  ]
  
  registry_row <- data.frame(
    curve_id = curve_id,
    dataset = dataset_name,
    subgroup_name = subgroup_name,
    subgroup_var = subgroup_var,
    subgroup_level = subgroup_level,
    n_total = n_total,
    n_unique_subject_uid = length(unique(df_curve$subject_uid)),
    n_event_transition = n_event,
    n_censor_transition = n_censor,
    n_right_censor_raw = n_right_censor,
    n_remission_raw = n_remission,
    followup_min_years = min_followup_years,
    followup_median_years = median_followup_years,
    followup_max_years = max_followup_years,
    last_event_time_years = last_event_time_years,
    n_risk_at_last_event = n_risk_at_last_event,
    censored_after_last_event = censored_after_last_event,
    event_fraction_transition = n_event / n_total,
    censor_fraction_transition = n_censor / n_total,
    stringsAsFactors = FALSE
  )
  
  list(
    fit = fit_obj,
    registry = registry_row,
    yearly = yearly_points,
    curve_data = curve_data
  )
}

# 🔴 입력: 생존분석용 데이터 구성 ===============================
## 🟠 읽기: 병합 코호트 불러오기 ===============================
raw_data <- utils::read.csv(
  file = data_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

if (nrow(raw_data) == 0L) {
  stop("입력 CSV가 비어 있습니다.", call. = FALSE)
}

## 🟠 점검: 핵심 컬럼과 식별키 검증 ===============================
required_columns <- c("id", "site", "days_followup", "status_num")
missing_required_columns <- setdiff(required_columns, names(raw_data))
if (length(missing_required_columns) > 0L) {
  stop(
    "필수 컬럼이 없습니다: ",
    paste(missing_required_columns, collapse = ", "),
    call. = FALSE
  )
}

if (length(subgroup_specs) > 0L && is.null(names(subgroup_specs))) {
  stop("subgroup_specs는 이름이 있는 character 벡터여야 합니다.", call. = FALSE)
}

analysis_df <- raw_data
analysis_df$id <- as.character(analysis_df$id)
analysis_df$site <- trimws(as.character(analysis_df$site))
analysis_df$site[analysis_df$site == ""] <- NA_character_

analysis_df$days_followup <- suppressWarnings(as.numeric(analysis_df$days_followup))
analysis_df$status_num <- suppressWarnings(as.integer(as.numeric(analysis_df$status_num)))

if ("sex_fact" %in% names(analysis_df)) {
  analysis_df$sex_fact <- trimws(as.character(analysis_df$sex_fact))
  analysis_df$sex_fact[analysis_df$sex_fact == ""] <- NA_character_
}

if (!("sex_fact" %in% names(analysis_df)) && ("sex_num" %in% names(analysis_df))) {
  sex_num_tmp <- suppressWarnings(as.integer(as.numeric(analysis_df$sex_num)))
  analysis_df$sex_fact <- ifelse(
    is.na(sex_num_tmp), NA_character_,
    ifelse(sex_num_tmp == 0L, "Female",
           ifelse(sex_num_tmp == 1L, "Male", NA_character_))
  )
}

essential_missing <- is.na(analysis_df$id) |
  is.na(analysis_df$site) |
  is.na(analysis_df$days_followup) |
  is.na(analysis_df$status_num)

if (any(essential_missing)) {
  stop(
    "핵심 분석 변수(id/site/days_followup/status_num)에 결측이 있습니다. 결측 행 수: ",
    sum(essential_missing),
    call. = FALSE
  )
}

if (any(!is.finite(analysis_df$days_followup))) {
  stop("days_followup에 비유한(infinite/non-finite) 값이 있습니다.", call. = FALSE)
}

if (any(analysis_df$days_followup < 0)) {
  stop("days_followup에 음수 값이 있습니다.", call. = FALSE)
}

unexpected_status_codes <- sort(setdiff(unique(analysis_df$status_num), allowed_status_codes))
if (length(unexpected_status_codes) > 0L) {
  stop(
    "status_num에 허용되지 않은 코드가 있습니다: ",
    paste(unexpected_status_codes, collapse = ", "),
    call. = FALSE
  )
}

analysis_df$subject_uid <- paste(analysis_df$site, analysis_df$id, sep = "||")
if (any(duplicated(analysis_df$subject_uid))) {
  stop(
    "site + id 조합이 유일키가 아닙니다. 중복 subject_uid 수: ",
    sum(duplicated(analysis_df$subject_uid)),
    call. = FALSE
  )
}

## 🟠 파생: 전이-전용 분석 변수 생성 ===============================
analysis_df$time_years <- analysis_df$days_followup / days_per_year
analysis_df$event_transition <- ifelse(analysis_df$status_num == transition_code, 1L, 0L)
analysis_df$censor_reason_transition <- ifelse(
  analysis_df$status_num == remission_code, "remission",
  ifelse(analysis_df$status_num == right_censor_code, "right_censoring",
         ifelse(analysis_df$status_num == transition_code, "event_transition", "unknown"))
)

analysis_df$site_std <- toupper(trimws(analysis_df$site))
site_pnu_values_std <- toupper(site_pnu_values)
site_snu_values_std <- toupper(site_snu_values)

resolved_subgroup_specs <- subgroup_specs
if (length(resolved_subgroup_specs) > 0L) {
  keep_index <- logical(length(resolved_subgroup_specs))
  names(keep_index) <- names(resolved_subgroup_specs)
  
  for (sg_name in names(resolved_subgroup_specs)) {
    sg_var <- resolved_subgroup_specs[[sg_name]]
    if (!sg_var %in% names(analysis_df)) {
      warning("subgroup 변수 '", sg_var, "'가 없어 '", sg_name, "' 분석을 건너뜁니다.")
      keep_index[sg_name] <- FALSE
    } else {
      analysis_df[[sg_var]] <- as.character(analysis_df[[sg_var]])
      analysis_df[[sg_var]][is.na(analysis_df[[sg_var]]) | trimws(analysis_df[[sg_var]]) == ""] <- missing_subgroup_label
      keep_index[sg_name] <- TRUE
    }
  }
  
  resolved_subgroup_specs <- resolved_subgroup_specs[keep_index]
}

dataset_list <- list(
  PNU = analysis_df[analysis_df$site_std %in% site_pnu_values_std, , drop = FALSE],
  SNU = analysis_df[analysis_df$site_std %in% site_snu_values_std, , drop = FALSE],
  merged = analysis_df
)

if (nrow(dataset_list$PNU) == 0L) {
  stop(
    "PNU 분기가 비어 있습니다. site_pnu_values 설정을 확인하세요. 현재 site 값: ",
    collapse_value(sort(unique(analysis_df$site))),
    call. = FALSE
  )
}
if (nrow(dataset_list$SNU) == 0L) {
  stop(
    "SNU 분기가 비어 있습니다. site_snu_values 설정을 확인하세요. 현재 site 값: ",
    collapse_value(sort(unique(analysis_df$site))),
    call. = FALSE
  )
}

raw_qc <- make_qc_row(analysis_df, qc_scope = "raw_input", dataset_name = "all_rows")
dataset_qc <- do.call(
  rbind,
  c(
    list(raw_qc),
    lapply(names(dataset_list), function(ds_name) {
      make_qc_row(dataset_list[[ds_name]], qc_scope = "analysis_branch", dataset_name = ds_name)
    })
  )
)

# 🔴 연산: Plain KM와 시점 요약 생성 ===============================
## 🟠 적합: 데이터셋별 overall 및 subgroup KM ===============================
curve_counter <- 0L
km_fits <- list()
registry_list <- list()
yearly_list <- list()
curve_data_list <- list()

for (dataset_name in names(dataset_list)) {
  ds_df <- dataset_list[[dataset_name]]
  
  curve_counter <- curve_counter + 1L
  curve_id <- sprintf("km_%03d", curve_counter)
  overall_result <- fit_km_curve(
    df_curve = ds_df,
    curve_id = curve_id,
    dataset_name = dataset_name,
    subgroup_name = "overall",
    subgroup_var = NA_character_,
    subgroup_level = "overall"
  )
  
  km_fits[[curve_id]] <- overall_result$fit
  registry_list[[length(registry_list) + 1L]] <- overall_result$registry
  yearly_list[[length(yearly_list) + 1L]] <- overall_result$yearly
  curve_data_list[[length(curve_data_list) + 1L]] <- overall_result$curve_data
  
  if (length(resolved_subgroup_specs) > 0L) {
    for (sg_name in names(resolved_subgroup_specs)) {
      sg_var <- resolved_subgroup_specs[[sg_name]]
      subgroup_levels <- sort(unique(ds_df[[sg_var]]), na.last = TRUE)
      
      for (sg_level in subgroup_levels) {
        sg_df <- ds_df[ds_df[[sg_var]] == sg_level, , drop = FALSE]
        if (nrow(sg_df) == 0L) next
        
        curve_counter <- curve_counter + 1L
        curve_id <- sprintf("km_%03d", curve_counter)
        
        subgroup_result <- fit_km_curve(
          df_curve = sg_df,
          curve_id = curve_id,
          dataset_name = dataset_name,
          subgroup_name = sg_name,
          subgroup_var = sg_var,
          subgroup_level = as.character(sg_level)
        )
        
        km_fits[[curve_id]] <- subgroup_result$fit
        registry_list[[length(registry_list) + 1L]] <- subgroup_result$registry
        yearly_list[[length(yearly_list) + 1L]] <- subgroup_result$yearly
        curve_data_list[[length(curve_data_list) + 1L]] <- subgroup_result$curve_data
      }
    }
  }
}

## 🟠 결합: registry, yearly, curve data 생성 ===============================
km_registry <- do.call(rbind, registry_list)
km_yearly_summary <- do.call(rbind, yearly_list)
km_curve_data <- do.call(rbind, curve_data_list)

km_registry <- km_registry[order(km_registry$dataset, km_registry$subgroup_name, km_registry$subgroup_level), ]
km_yearly_summary <- km_yearly_summary[order(km_yearly_summary$dataset, km_yearly_summary$subgroup_name, km_yearly_summary$subgroup_level, km_yearly_summary$time_year), ]
km_curve_data <- km_curve_data[order(km_curve_data$dataset, km_curve_data$subgroup_name, km_curve_data$subgroup_level, km_curve_data$time_year), ]

overall_curve_ids <- km_registry$curve_id[km_registry$subgroup_name == "overall"]
group_curve_ids <- km_registry$curve_id[km_registry$subgroup_name != "overall"]

analysis_spec <- data.frame(
  item = c(
    "script_name",
    "script_version",
    "created_at",
    "data_path",
    "export_path",
    "dataset_branches",
    "site_pnu_values",
    "site_snu_values",
    "endpoint_name",
    "event_definition",
    "censoring_definition",
    "time_origin",
    "raw_time_unit",
    "analysis_time_unit",
    "days_per_year",
    "time_grid_years",
    "km_conf_type",
    "subgroup_specs",
    "missing_subgroup_label",
    "raw_site_levels",
    "n_rows_raw_input",
    "n_rows_PNU",
    "n_rows_SNU",
    "n_rows_merged"
  ),
  value = c(
    script_name,
    script_version,
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    data_path,
    export_path,
    collapse_value(names(dataset_list)),
    collapse_value(site_pnu_values),
    collapse_value(site_snu_values),
    "Transition to schizophrenia (transition-only Step1 KM)",
    "event_transition = 1 if status_num == 1",
    "status_num in c(0, 2) treated as censoring; remission is censored",
    "cohort entry",
    "days_followup in days",
    "years since cohort entry",
    as.character(days_per_year),
    paste(time_grid_years, collapse = ", "),
    conf_type_km,
    if (length(resolved_subgroup_specs) == 0L) {
      "overall only"
    } else {
      paste(paste(names(resolved_subgroup_specs), resolved_subgroup_specs, sep = ":"), collapse = "; ")
    },
    missing_subgroup_label,
    collapse_value(sort(unique(analysis_df$site))),
    as.character(nrow(analysis_df)),
    as.character(nrow(dataset_list$PNU)),
    as.character(nrow(dataset_list$SNU)),
    as.character(nrow(dataset_list$merged))
  ),
  stringsAsFactors = FALSE
)

km_rds_object <- list(
  step = "Step1_plain_KM",
  script_name = script_name,
  script_version = script_version,
  created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  config = list(
    data_path = data_path,
    export_path = export_path,
    site_pnu_values = site_pnu_values,
    site_snu_values = site_snu_values,
    time_grid_years = time_grid_years,
    days_per_year = days_per_year,
    transition_code = transition_code,
    right_censor_code = right_censor_code,
    remission_code = remission_code,
    conf_type_km = conf_type_km,
    missing_subgroup_label = missing_subgroup_label,
    subgroup_specs = resolved_subgroup_specs
  ),
  analysis_spec = analysis_spec,
  raw_data_qc = raw_qc,
  dataset_qc = dataset_qc,
  analysis_data = analysis_df,
  km_registry = km_registry,
  km_yearly_summary = km_yearly_summary,
  km_curve_data = km_curve_data,
  km_overall = km_fits[overall_curve_ids],
  km_by_group = km_fits[group_curve_ids],
  km_fits = km_fits,
  session_info = utils::capture.output(sessionInfo())
)

# 🔴 그림: KM 패널 렌더링 ===============================
## 🟠 작성: 곡선 데이터 기반 PNG 저장 ===============================
plot_file <- file.path(export_path, "step1_km_panels.png")

plot_data <- km_curve_data
plot_data$dataset <- factor(plot_data$dataset, levels = names(dataset_list))
subgroup_order <- unique(c("overall", names(resolved_subgroup_specs)))
subgroup_order <- subgroup_order[subgroup_order %in% unique(plot_data$subgroup_name)]
plot_data$subgroup_name <- factor(plot_data$subgroup_name, levels = subgroup_order)

km_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = time_year,
    y = survival,
    color = subgroup_level,
    group = curve_id
  )
) +
  ggplot2::geom_step(linewidth = 0.7, na.rm = TRUE) +
  ggplot2::facet_grid(subgroup_name ~ dataset, scales = "free_x", labeller = ggplot2::label_both) +
  ggplot2::labs(
    title = "Step 1: Plain KM curves",
    x = "Years since cohort entry",
    y = "KM survival probability",
    color = "Curve"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank()
  )

plot_width <- max(10, 3.5 * length(dataset_list))
plot_height <- max(6, 2.6 * length(unique(plot_data$subgroup_name)))

ggplot2::ggsave(
  filename = plot_file,
  plot = km_plot,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi
)

# 🔴 저장: Step1 산출물 내보내기 ===============================
## 🟠 기록: CSV, RDS, manifest 저장 ===============================
analysis_spec_file <- file.path(export_path, "step1_analysis_spec.csv")
data_qc_file <- file.path(export_path, "step1_data_qc.csv")
km_registry_file <- file.path(export_path, "step1_km_registry.csv")
km_yearly_file <- file.path(export_path, "step1_km_yearly_summary.csv")
km_curve_data_file <- file.path(export_path, "step1_km_curve_data.csv")
km_rds_file <- file.path(export_path, "step1_km_objects.rds")
manifest_file <- file.path(export_path, "step1_export_manifest.csv")

write_csv_utf8(analysis_spec, analysis_spec_file)
write_csv_utf8(dataset_qc, data_qc_file)
write_csv_utf8(km_registry, km_registry_file)
write_csv_utf8(km_yearly_summary, km_yearly_file)
write_csv_utf8(km_curve_data, km_curve_data_file)
saveRDS(km_rds_object, file = km_rds_file)

export_manifest <- data.frame(
  file_name = c(
    basename(analysis_spec_file),
    basename(data_qc_file),
    basename(km_registry_file),
    basename(km_yearly_file),
    basename(km_curve_data_file),
    basename(km_rds_file),
    basename(plot_file)
  ),
  file_path = c(
    analysis_spec_file,
    data_qc_file,
    km_registry_file,
    km_yearly_file,
    km_curve_data_file,
    km_rds_file,
    plot_file
  ),
  description = c(
    "Step1 분석 명세",
    "원자료 및 데이터셋 분기 QC",
    "KM curve registry",
    "1-10년 KM 요약표",
    "KM step curve plotting/source data",
    "KM survfit objects and reusable registry",
    "KM panel figure generated from step1_km_curve_data.csv"
  ),
  n_rows = c(
    nrow(analysis_spec),
    nrow(dataset_qc),
    nrow(km_registry),
    nrow(km_yearly_summary),
    nrow(km_curve_data),
    length(km_fits),
    NA_integer_
  ),
  created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  stringsAsFactors = FALSE
)

write_csv_utf8(export_manifest, manifest_file)